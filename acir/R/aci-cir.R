################################################################################
## aci-cir.R - linear-memory compiled forward-CIR summaries
################################################################################


#' Precompute the fixed-lag primitives used by forward CIR rows (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible explicit single-step filter.
#' @param from First interval index needed. An anchor row only ever reads
#'   primitives at or after its own index, so a windowed reduction skips the
#'   intervals before its earliest anchor.
#' @param rec Covariance-policy recorder from `.aci_reg_new()`.
#' @returns Per-interval affine updates and conservative suffix estimates.
#' @noRd
.compiled_forward_primitives <- function(bundle, filter, from = 1L,
                                         rec = .aci_reg_for(NULL, bundle$t)) {
  N <- bundle$N
  l <- bundle$l
  from <- max(1L, min(as.integer(from), N))
  E <- vector("list", N)
  dmu <- matrix(NA_real_, N, l)
  dR <- vector("list", N)
  one_mu <- matrix(NA_real_, N, l)
  one_R <- vector("list", N)
  s_n <- r_n <- e_n <- numeric(N)

  for (n in seq.int(from, N)) {
    rec$j <- n
    co <- .compiled_co(bundle, n)
    aux <- .thmD1_aux_compiled(
      bundle, n, filter$cov[, , n], co = co, rec = rec
    )
    one <- .onelag_stats(
      co, aux, filter$mean[n, ], filter$cov[, , n],
      filter$mean[n + 1L, ], filter$cov[, , n + 1L],
      bundle$x[n + 1L, ] - bundle$x[n, ], bundle$dt, l, rec
    )
    E[[n]] <- aux$E
    one_mu[n, ] <- one$mu
    one_R[[n]] <- one$R
    dmu[n, ] <- one$mu - filter$mean[n, ]
    dR[[n]] <- one$R - filter$cov[, , n]
    s_n[n] <- sqrt(sum(dmu[n, ]^2))
    r_n[n] <- sqrt(sum(dR[[n]]^2))
    e_n[n] <- if (l == 1L) abs(aux$E[1L, 1L]) else norm(aux$E, "2")
  }

  T2 <- Ub <- numeric(N + 1L)
  for (n in N:from) {
    eh2 <- max(1, e_n[n])^2
    T2[n] <- if (n < N) s_n[n + 1L]^2 + eh2 * T2[n + 1L] else 0
    Ub[n] <- if (n < N) r_n[n + 1L] + eh2 * Ub[n + 1L] else 0
    if (!is.finite(T2[n]) || T2[n] > 1e12) T2[n] <- Inf
    if (!is.finite(Ub[n]) || Ub[n] > 1e12) Ub[n] <- Inf
  }

  list(
    E = E, dmu = dmu, dR = dR, one_mu = one_mu, one_R = one_R,
    T2 = T2, Ub = Ub
  )
}


#' Reduce compiled forward lag rows without retaining the triangle (internal)
#'
#' Each anchor row is formed, reduced to its forward-CIR summaries, and then
#' discarded.  Arithmetic, adaptive stopping, zero padding, tail estimates and
#' quadrature are the same as `lag_table(mode = "forward")` followed by
#' [aci_range()].  Peak retained row storage is therefore linear in record
#' length instead of quadratic.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible explicit single-step filter.
#' @param smoother The complete Theorem-3 reference smoother.
#' @param method Forward-CIR functional.
#' @param epsilon Optional subjective thresholds.
#' @param tol,window,max_lag Lag-row retention controls.
#' @param quadrature Objective quadrature convention.
#' @param simpson_close Even-grid Simpson closing convention.
#' @param epsilon_grid Threshold nodes for the MATLAB-compatibility quadrature.
#' @param anchors Integer anchor indices to compute, in reporting order.
#' @param floor_M Strength floor driving the status vocabulary.
#' @param off Subjective read-out offset: 0 counts cells, 1 reads lag time.
#' @param rec Covariance-policy recorder from `.aci_reg_new()`.
#' @returns Raw range, strength, subjective, status and truncation summaries.
#' @noRd
.compiled_forward_reduce <- function(
    bundle, filter, smoother, method, epsilon, tol, window, max_lag,
    quadrature, simpson_close, epsilon_grid, anchors, floor_M, off,
    rec = .aci_reg_for(NULL, bundle$t)) {
  N1 <- bundle$N1
  N <- bundle$N
  l <- bundle$l
  na <- length(anchors)
  prim <- .compiled_forward_primitives(bundle, filter, from = min(anchors),
                                       rec = rec)

  sch <- vector("list", na)
  lam <- diagv <- Rs <- numeric(na)
  for (i in seq_len(na)) {
    j <- anchors[i]
    rec$j <- j
    cs <- .cov_guard_chol(smoother$cov[, , j], rec, "metric_reference")
    sch[[i]] <- cs$ch
    Rs[i] <- cs$R[1L]
    lam[i] <- min(eigen(cs$R, symmetric = TRUE, only.values = TRUE)$values)
    diagv[i] <- .kl_fast(
      smoother$mean[j, ], cs$ch,
      filter$mean[j, ], filter$cov[, , j], rec
    )
  }

  tau <- strength <- rep(NA_real_, na)
  tail <- numeric(na)
  retained_lag <- integer(na)
  status <- character(na)
  subjective <- if (is.null(epsilon)) NULL else
    matrix(0, nrow = na, ncol = length(epsilon))
  structurally_truncated <- FALSE
  ## A scalar hidden state forms each row as one vector expression.
  sp <- if (l == 1L) .cir_scalar_primitives(prim$E, prim$dmu, prim$dR)

  for (i in seq_len(na)) {
    j <- anchors[i]
    rec$j <- j
    row <- diagv[i]
    lag_kept <- 0L
    row_tail <- 0

    if (j <= N && l == 1L) {
      r1 <- .cir_scalar_row(
        j, sp, prim$T2, prim$Ub, N, filter$mean[j, 1L], filter$cov[1L, 1L, j],
        smoother$mean[j, 1L], Rs[i], diagv[i], lam[i], tol, window, max_lag,
        rec
      )
      row <- r1$row
      lag_kept <- r1$lag_kept
      row_tail <- r1$tail
    } else if (j <= N) {
      mu <- as.numeric(prim$one_mu[j, ])
      R <- matrix(prim$one_R[[j]], l, l)
      D <- diag(l)
      row <- c(row, .kl_fast(
        smoother$mean[j, ], sch[[i]], mu, R, rec
      ))
      lag_kept <- 1L
      below_count <- 0L

      if (j < N) for (n in seq.int(j + 1L, N)) {
        lag_length <- (n + 1L) - j
        if (is.finite(max_lag) && lag_length > max_lag) {
          lag_kept <- as.integer(max_lag)
          row_tail <- utils::tail(row, 1L)
          break
        }

        D <- D %*% prim$E[[n - 1L]]
        mu <- mu + drop(D %*% prim$dmu[n, ])
        R <- sym(R + D %*% prim$dR[[n]] %*% t(D))
        Pval <- .kl_fast(smoother$mean[j, ], sch[[i]], mu, R, rec)
        row <- c(row, Pval)
        lag_kept <- lag_length

        if (tol > 0) {
          Dn <- if (l == 1L) abs(D[1L, 1L]) else sqrt(sum(D^2))
          future <- 1.5 * (
            Dn^2 * prim$T2[n] / (2 * lam[i]) +
              Dn^2 * prim$Ub[n] * sqrt(l) / (2 * lam[i])
          )
          ok <- is.finite(future) && (Pval + future) < tol
          below_count <- if (ok) below_count + 1L else 0L
          if (below_count >= window || lag_length >= max_lag) {
            row_tail <- if (is.finite(future)) Pval + future else Pval
            break
          }
        } else if (lag_length >= max_lag) {
          row_tail <- Pval
          break
        }
      }
    }

    retained_lag[i] <- lag_kept
    tail[i] <- row_tail
    full_length <- N1 - j + 1L
    if (length(row) < full_length) {
      structurally_truncated <- TRUE
      row <- c(row, numeric(full_length - length(row)))
    }
    reduced <- .fwd_lengths(
      row, bundle$dt, method, quadrature, simpson_close, epsilon_grid
    )
    tau[i] <- reduced[["tau"]]
    strength[i] <- reduced[["M"]]
    status[i] <- .cir_status(row, strength[i], floor_M)
    if (!is.null(subjective)) {
      for (q in seq_along(epsilon)) {
        hit <- which(row > epsilon[q])
        subjective[i, q] <- if (!length(hit)) 0 else
          bundle$dt * (max(hit) - off)
      }
    }
  }

  list(
    tau = tau, strength = strength, subjective = subjective,
    tail_bound = tail, retained_lag = retained_lag, status = status,
    structurally_truncated = structurally_truncated
  )
}


#' Forward CIR directly from a compiled CGNS run (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter Optional compatible filter path.
#' @param init Optional Gaussian prior.
#' @param method,epsilon,min_M,masked_value,quadrature,simpson_close,anchors,convention,epsilon_grid
#'   Arguments matching [aci_range.lag_table()].
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @param ... Rejected unused arguments.
#' @returns A `cir_result` without a retained lag table.
#' @noRd
.forward_cir_compiled <- function(
    bundle, filter = NULL, init = NULL,
    method = c("exact", "l1_linf"), epsilon = NULL, min_M = "auto",
    masked_value = c("na", "zero"),
    quadrature = c("simpson", "sum", "matlab_eps_grid"),
    simpson_close = c("quadratic", "trapezoid"),
    anchors = NULL, convention = c("count", "lag_time"),
    epsilon_grid = NULL, regularize = NULL, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims",
              "Unused arguments were supplied to aci_range().")
  rec <- .aci_reg_for(regularize, bundle$t)
  method <- match.arg(method)
  masked_value <- match.arg(masked_value)
  quadrature <- match.arg(quadrature)
  simpson_close <- match.arg(simpson_close)
  convention <- match.arg(convention)
  grid <- .cir_compat_grid(quadrature, method, epsilon_grid)
  .validate_compiled_cgns(
    bundle, conditional = bundle$conditional, scalar = FALSE
  )
  if (!is.null(epsilon) &&
      (!is.numeric(epsilon) || !length(epsilon) || any(!is.finite(epsilon)) ||
       any(epsilon < 0)))
    aci_abort(
      "aci_error_dims", "epsilon must contain finite non-negative thresholds."
    )
  idx <- .cir_anchors(anchors, bundle$N1)
  floor_M <- .cir_min_strength(if (is.null(min_M)) "auto" else min_M)

  tol <- getOption("aci.default_tol", 1e-8)
  window <- 3L
  max_lag <- Inf
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    aci_abort("aci_error_dims", "tol must be one finite non-negative number.")

  filt <- filter
  if (!is.null(filt))
    .validate_gaussian_path(
      filt, bundle$obs, bundle$l, "filter", bundle$conditional,
      model = bundle$model, source_model = bundle$source_model
    )
  if (!is.null(filt) && !is.null(init) &&
      !.same_gaussian_init(init, filt$meta$init, bundle$l))
    aci_abort(
      "aci_error_dims",
      "init conflicts with the prior stored on the supplied filter."
    )
  if (!is.null(filt) &&
      (!identical(filt$meta$stepper %||% "explicit", "explicit") ||
       (filt$meta$nsub %||% 1L) != 1L)) {
    aci_warn(
      "aci_warn_stepper",
      paste(
        "lag_table requires the explicit single-step filter/smoother (the",
        "Theorem 3 recursions are exact for that discretization);",
        "recomputing both internally."
      )
    )
    init <- filt$meta$init %||% init
    filt <- NULL
  }
  filt <- filt %||% .cgns_filter_compiled(
    bundle, init = init, stepper = "explicit", nsub = 1L,
    validate = FALSE, regularize = rec
  )
  smoother <- .smoother_thmD1_compiled(
    bundle, filt, validate = FALSE, warn_cost = FALSE, regularize = rec
  )
  reduced <- .compiled_forward_reduce(
    bundle, filt, smoother, method = method, epsilon = epsilon,
    tol = tol, window = window, max_lag = max_lag,
    quadrature = quadrature, simpson_close = simpson_close,
    epsilon_grid = grid, anchors = idx, floor_M = floor_M,
    off = if (identical(convention, "count")) 0 else 1,
    rec = rec
  )

  tau <- reduced$tau
  if (!is.null(min_M)) {
    threshold <- .cir_min_strength(min_M)
    bad <- which(reduced$strength < threshold)
    if (length(bad)) {
      tau[bad] <- if (masked_value == "zero") 0 else NA_real_
      aci_warn("aci_warn_low_signal", sprintf(
        paste0(
          "%d forward CIR values masked (M < %.3g); interpret CIRs jointly ",
          "with the ACI metric (Andreou & Chen 2026, Remark B.4)."
        ),
        length(bad), threshold
      ))
    }
  }

  tt <- bundle$t[idx]
  new_cir_result(
    tt, tau, reduced$strength, "forward", method, bundle$dt,
    interval = cbind(tt, tt + tau, deparse.level = 0L),
    tail_bound = reduced$tail_bound,
    subjective = reduced$subjective,
    status = factor(reduced$status, levels = .cir_status_levels),
    meta = list(
      quadrature = quadrature,
      simpson_close = simpson_close,
      convention = convention,
      epsilon_grid = grid,
      table_conditional = bundle$conditional,
      anchors = idx,
      structurally_truncated = reduced$structurally_truncated,
      regularization = .aci_reg_freeze(rec)
    ),
    bound = .cir_bound(method, quadrature, reduced$structurally_truncated)
  )
}
