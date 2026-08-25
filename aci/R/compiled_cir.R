################################################################################
## compiled_cir.R - linear-memory compiled forward-CIR summaries
################################################################################


#' Precompute the fixed-lag primitives used by forward CIR rows (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible explicit single-step filter.
#' @returns Per-interval affine updates and conservative suffix estimates.
#' @noRd
.compiled_forward_primitives <- function(bundle, filter) {
  N <- bundle$N
  l <- bundle$l
  E <- vector("list", N)
  dmu <- matrix(NA_real_, N, l)
  dR <- vector("list", N)
  one_mu <- matrix(NA_real_, N, l)
  one_R <- vector("list", N)
  s_n <- r_n <- e_n <- numeric(N)

  for (n in seq_len(N)) {
    co <- .compiled_co(bundle, n)
    aux <- .thmD1_aux_compiled(
      bundle, n, filter$cov[, , n], co = co
    )
    one <- .onelag_stats(
      co, aux, filter$mean[n, ], filter$cov[, , n],
      filter$mean[n + 1L, ], filter$cov[, , n + 1L],
      bundle$x[n + 1L, ] - bundle$x[n, ], bundle$dt, l
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
  for (n in N:1L) {
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
#' [forward_cir()].  Peak retained row storage is therefore linear in record
#' length instead of quadratic.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible explicit single-step filter.
#' @param smoother The complete Theorem-3 reference smoother.
#' @param method Forward-CIR functional.
#' @param eps Optional subjective thresholds.
#' @param tol,window,max_lag Lag-row retention controls.
#' @param quadrature Ratio quadrature convention.
#' @param simpson_close Even-grid Simpson closing convention.
#' @returns Raw range, strength, subjective and truncation summaries.
#' @noRd
.compiled_forward_reduce <- function(
    bundle, filter, smoother, method, eps, tol, window, max_lag,
    quadrature, simpson_close) {
  N1 <- bundle$N1
  N <- bundle$N
  l <- bundle$l
  prim <- .compiled_forward_primitives(bundle, filter)

  sch <- vector("list", N1)
  lam <- diagv <- numeric(N1)
  for (j in seq_len(N1)) {
    cs <- .chol_or_floor(smoother$cov[, , j])
    sch[[j]] <- cs$ch
    lam[j] <- min(eigen(cs$R, symmetric = TRUE, only.values = TRUE)$values)
    diagv[j] <- .kl_fast(
      smoother$mean[j, ], cs$ch,
      filter$mean[j, ], filter$cov[, , j]
    )
  }

  tau <- strength <- rep(NA_real_, N1)
  tail <- numeric(N1)
  retained_lag <- integer(N1)
  subjective <- if (is.null(eps)) NULL else
    matrix(0, nrow = N1, ncol = length(eps))
  structurally_truncated <- FALSE

  for (j in seq_len(N1)) {
    row <- diagv[j]
    lag_kept <- 0L
    row_tail <- 0

    if (j <= N) {
      mu <- as.numeric(prim$one_mu[j, ])
      R <- matrix(prim$one_R[[j]], l, l)
      D <- diag(l)
      row <- c(row, .kl_fast(
        smoother$mean[j, ], sch[[j]], mu, R
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
        Pval <- .kl_fast(smoother$mean[j, ], sch[[j]], mu, R)
        row <- c(row, Pval)
        lag_kept <- lag_length

        if (tol > 0) {
          Dn <- if (l == 1L) abs(D[1L, 1L]) else sqrt(sum(D^2))
          future <- 1.5 * (
            Dn^2 * prim$T2[n] / (2 * lam[j]) +
              Dn^2 * prim$Ub[n] * sqrt(l) / (2 * lam[j])
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

    retained_lag[j] <- lag_kept
    tail[j] <- row_tail
    full_length <- N1 - j + 1L
    if (length(row) < full_length) {
      structurally_truncated <- TRUE
      row <- c(row, numeric(full_length - length(row)))
    }
    reduced <- .fwd_lengths(
      row, bundle$dt, method, quadrature, simpson_close
    )
    tau[j] <- reduced[["tau"]]
    strength[j] <- reduced[["M"]]
    if (!is.null(subjective)) {
      for (q in seq_along(eps)) {
        hit <- which(row > eps[q])
        subjective[j, q] <- if (!length(hit)) 0 else
          bundle$dt * (max(hit) - 1L)
      }
    }
  }

  list(
    tau = tau, strength = strength, subjective = subjective,
    tail_bound = tail, retained_lag = retained_lag,
    structurally_truncated = structurally_truncated
  )
}


#' Forward CIR directly from a compiled CGNS run (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter Optional compatible filter path.
#' @param init Optional Gaussian prior.
#' @param method,eps,min_M,masked_value,quadrature,simpson_close Arguments
#'   matching [forward_cir.lag_table()].
#' @param ... Rejected unused arguments.
#' @returns A `cir_result` without a retained lag table.
#' @noRd
.forward_cir_compiled <- function(
    bundle, filter = NULL, init = NULL,
    method = c("exact", "l1_linf"), eps = NULL, min_M = "auto",
    masked_value = c("na", "zero"),
    quadrature = c("simpson", "sum"),
    simpson_close = c("quadratic", "trapezoid"), ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to forward_cir().")
  method <- match.arg(method)
  masked_value <- match.arg(masked_value)
  quadrature <- match.arg(quadrature)
  simpson_close <- match.arg(simpson_close)
  .validate_compiled_cgns(
    bundle, nontarget = bundle$nontarget, scalar = FALSE
  )
  if (!is.null(eps) &&
      (!is.numeric(eps) || !length(eps) || any(!is.finite(eps)) ||
       any(eps < 0)))
    aci_abort(
      "aci_error_dims", "eps must contain finite non-negative thresholds."
    )

  tol <- getOption("aci.default_tol", 1e-8)
  window <- 3L
  max_lag <- Inf
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    aci_abort("aci_error_dims", "tol must be one finite non-negative number.")

  filt <- filter
  if (!is.null(filt))
    .validate_gaussian_path(
      filt, bundle$obs, bundle$l, "filter", bundle$nontarget,
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
    validate = FALSE
  )
  smoother <- .smoother_thmD1_compiled(
    bundle, filt, validate = FALSE, warn_cost = FALSE
  )
  reduced <- .compiled_forward_reduce(
    bundle, filt, smoother, method = method, eps = eps,
    tol = tol, window = window, max_lag = max_lag,
    quadrature = quadrature, simpson_close = simpson_close
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

  new_cir_result(
    bundle$t, tau, reduced$strength, "forward", method, bundle$dt,
    interval = cbind(bundle$t, bundle$t + tau),
    tail_bound = reduced$tail_bound,
    subjective = reduced$subjective,
    meta = list(
      quadrature = quadrature,
      simpson_close = simpson_close,
      table_nontarget = bundle$nontarget,
      structurally_truncated = reduced$structurally_truncated
    ),
    bound = if (reduced$structurally_truncated) {
      if (method == "exact") "objective_on_truncated_table"
      else "lower_ratio_on_truncated_table_only"
    } else NULL
  )
}


#' Slice a compiled CGNS bundle at a prefix boundary (internal)
#'
#' Coefficients are immutable and filtering is causal, so backward-CIR reference
#' times can share one maximal-prefix compilation and filter. This constructor
#' keeps the ordinary compiled contract and provenance on a shorter prefix.
#'
#' @param bundle A compiled CGNS bundle.
#' @param N1 Number of grid points retained, at least two.
#' @returns A prefix `compiled_cgns` bundle.
#' @noRd
.slice_compiled_cgns <- function(bundle, N1) {
  if (!inherits(bundle, "compiled_cgns") || length(N1) != 1L ||
      !is.numeric(N1) || !is.finite(N1) || N1 != floor(N1) ||
      N1 < 2L || N1 > bundle$N1)
    aci_abort(
      "aci_error_compiled_contract",
      "A compiled CGNS bundle and a valid prefix length are required."
    )
  N1 <- as.integer(N1)
  N <- N1 - 1L
  k <- bundle$k
  l <- bundle$l
  source_obs <- observed_trajectory(
    bundle$source_obs$t[seq_len(N1)],
    bundle$source_obs$x[seq_len(N1), , drop = FALSE]
  )
  resolved_obs <- observed_trajectory(
    bundle$t[seq_len(N1)], bundle$x[seq_len(N1), , drop = FALSE]
  )
  coefficients <- list(
    Lx = array(bundle$coefficients$Lx[, , seq_len(N1), drop = FALSE],
               c(k, l, N1)),
    fx = matrix(bundle$coefficients$fx[seq_len(N1), , drop = FALSE], N1, k),
    Ly = array(bundle$coefficients$Ly[, , seq_len(N1), drop = FALSE],
               c(l, l, N1)),
    fy = matrix(bundle$coefficients$fy[seq_len(N1), , drop = FALSE], N1, l),
    gxx = array(bundle$coefficients$gxx[, , seq_len(N1), drop = FALSE],
                c(k, k, N1)),
    gyy = array(bundle$coefficients$gyy[, , seq_len(N1), drop = FALSE],
                c(l, l, N1)),
    gyx = array(bundle$coefficients$gyx[, , seq_len(N1), drop = FALSE],
                c(l, k, N1)),
    gxx_weight = array(
      bundle$coefficients$gxx_weight[, , seq_len(N), drop = FALSE],
      c(k, k, N)
    )
  )
  rs <- list(
    model = bundle$model,
    obs = resolved_obs,
    likelihood_idx = bundle$likelihood_idx,
    tag = bundle$nontarget
  )
  .new_compiled_cgns(
    rs, bundle$source_model, source_obs, coefficients,
    correlated_noise = bundle$correlated_noise,
    realization = paste0(bundle$realization, "_prefix")
  )
}


#' Slice a Gaussian filter path at a prefix boundary (internal)
#'
#' @param filter A filter path generated on the full bundle.
#' @param bundle The matching prefix bundle.
#' @returns A compatible prefix filter.
#' @noRd
.slice_compiled_filter <- function(filter, bundle) {
  N1 <- bundle$N1
  path <- new_da_path(
    bundle$t,
    filter$mean[seq_len(N1), , drop = FALSE],
    filter$cov[, , seq_len(N1), drop = FALSE],
    "filter"
  )
  path$meta <- filter$meta
  path$meta$obs_x <- bundle$x
  path$meta$model <- bundle$model
  path$meta$nontarget <- bundle$nontarget
  path$meta$source_model <- bundle$source_model
  path
}
