################################################################################
## compiled_lag.R - compiled Theorem-3 smoother and lag-table execution
################################################################################


#' Theorem-3 auxiliaries from one compiled coefficient record (internal)
#'
#' @param bundle A validated `compiled_cgns` object.
#' @param j One-based interval index.
#' @param Rf Filtered covariance at the interval start.
#' @param co Optional coefficient record already read with `.compiled_co()`.
#' @returns A list with the auxiliary matrices `E` and `F`.
#' @noRd
.thmD1_aux_compiled <- function(bundle, j, Rf, co = NULL) {
  co <- co %||% .compiled_co(bundle, j)
  Gi <- .compiled_ginv(bundle, j)
  l <- bundle$l
  dt <- bundle$dt
  Rfi <- chol_solve(Rf, diag(l), "Rf")
  gxy <- t(co$gyx)
  Gx <- co$Lx + gxy %*% Rfi
  Gy <- co$Ly + co$gyy %*% Rfi
  K <- Gi %*% Gx
  H <- Rfi %*% (co$Ly %*% Rf + Rf %*% t(co$Ly) + co$gyy)
  E <- diag(l) + (co$gyx %*% Gi %*% Gx - Gy) * dt
  KR <- K %*% Rf
  F_ <- -Rf %*% (
    t(K) +
      (t(Gx) %*% KR %*% t(K) - Rfi %*% t(H) %*% Rf %*% t(K) +
        t(co$Ly) %*% t(K)) * dt -
      t(co$Lx) %*% (Gi + KR %*% t(K) * dt)
  )
  list(E = E, F = F_)
}


#' Complete compiled online Theorem-3 smoother (internal)
#'
#' Preserves the explicit single-step restriction and backward affine recursion
#' while reading every coefficient and Gram inverse from the realised bundle.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filt A compatible explicit single-step filter path.
#' @param validate Validate the bundle and filter provenance before execution.
#' @param warn_cost Reserved for the established cost-warning contract.
#' @returns A `da_path_gaussian` smoother with route `"thmD1"`.
#' @noRd
.smoother_thmD1_compiled <- function(bundle, filt, validate = TRUE,
                                     warn_cost = TRUE) {
  if (isTRUE(validate)) {
    .validate_compiled_cgns(
      bundle, nontarget = bundle$nontarget, scalar = FALSE
    )
    .validate_gaussian_path(
      filt, bundle$obs, bundle$l, "filter", bundle$nontarget,
      model = bundle$model, source_model = bundle$source_model
    )
  }
  if (!identical(filt$meta$stepper %||% "explicit", "explicit") ||
      (filt$meta$nsub %||% 1L) != 1L)
    aci_abort(
      "aci_error_stepper",
      "Theorem 3 smoothing requires an explicit single-step filter."
    )

  N1 <- bundle$N1
  l <- bundle$l
  MU <- matrix(NA_real_, N1, l)
  CV <- array(NA_real_, c(l, l, N1))
  MU[N1, ] <- filt$mean[N1, ]
  CV[, , N1] <- filt$cov[, , N1]
  for (j in (N1 - 1L):1L) {
    co <- .compiled_co(bundle, j)
    aux <- .thmD1_aux_compiled(bundle, j, filt$cov[, , j], co = co)
    st <- .onelag_stats(
      co, aux, filt$mean[j, ], filt$cov[, , j],
      MU[j + 1L, ], CV[, , j + 1L],
      bundle$x[j + 1L, ] - bundle$x[j, ], bundle$dt, l
    )
    MU[j, ] <- st$mu
    CV[, , j] <- st$R
  }
  p <- new_da_path(bundle$t, MU, CV, "smoother")
  p$meta$route <- "thmD1"
  p$meta$nsub <- 1L
  p$meta$init <- filt$meta$init
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$nontarget <- bundle$nontarget
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p
}


#' Compiled row recursion behind a lag table (internal)
#'
#' Its indexing, KL orientation, adaptive freeze rule, retained-cell
#' representation and covariance flooring preserve the established route.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filt The forward filter path.
#' @param smoo The complete Theorem-3 reference smoother.
#' @param mode One of `"forward"`, `"one_lag"`, `"full"` or the internal
#'   `"smoother_only"` mode.
#' @param tol Non-negative adaptive storage tolerance.
#' @param window Consecutive steps below `tol` required before freezing.
#' @param max_lag Maximum retained positive lag.
#' @returns The internal retained-row and storage-diagnostic result list.
#' @noRd
.lagtable_core_compiled <- function(bundle, filt, smoo, mode, tol, window,
                                    max_lag) {
  N1 <- bundle$N1
  N <- bundle$N
  dt <- bundle$dt
  l <- bundle$l
  MUf <- filt$mean
  CVf <- filt$cov
  diagv <- rep(NA_real_, N1)
  if (!is.null(smoo)) for (j in seq_len(N1))
    diagv[j] <- unname(gaussian_kl(
      smoo$mean[j, ], smoo$cov[, , j],
      MUf[j, ], CVf[, , j], decompose = FALSE
    ))

  if (mode == "one_lag") {
    coN1 <- .compiled_co(bundle, N)
    auxN1 <- .thmD1_aux_compiled(bundle, N, CVf[, , N], co = coN1)
    ol <- .onelag_stats(
      coN1, auxN1, MUf[N, ], CVf[, , N], MUf[N1, ], CVf[, , N1],
      bundle$x[N1, ] - bundle$x[N, ], dt, l
    )
    dmu <- ol$mu - MUf[N, ]
    dR <- ol$R - CVf[, , N]
    P <- rep(0, N1)
    D <- diag(l)
    for (jj in N:1) {
      if (jj < N) {
        co <- .compiled_co(bundle, jj)
        E <- .thmD1_aux_compiled(bundle, jj, CVf[, , jj], co = co)$E
        D <- E %*% D
      }
      RsjN <- smoo$cov[, , jj]
      A <- sym(D %*% dR %*% t(D))
      Rlag <- spd_floor(RsjN - A)
      v <- drop(D %*% dmu)
      sig <- 0.5 * sum(v * chol_solve(Rlag, v, "Rlag"))
      trA <- sum(diag(chol_solve(Rlag, A, "Rlag")))
      ld <- logdet_chol(RsjN) - logdet_chol(Rlag)
      P[jj] <- max(sig + 0.5 * (trA - ld), 0)
    }
    return(list(
      diag = diagv, rows = NULL, L = NULL, tailbnd = NULL,
      onelag = P, stop_index = NA_integer_
    ))
  }

  rows <- if (mode != "smoother_only") vector("list", N1) else NULL
  if (!is.null(rows)) for (j in seq_len(N1)) rows[[j]] <- diagv[j]
  L <- rep(NA_integer_, N1)
  tailb <- rep(0, N1)
  act_D <- vector("list", N1)
  act_mu <- vector("list", N1)
  act_R <- vector("list", N1)
  act_cnt <- integer(N1)
  active <- frozen <- logical(N1)
  Ehist <- vector("list", N1)
  smu <- matrix(NA_real_, N1, l)
  scov <- array(NA_real_, c(l, l, N1))
  smu[N1, ] <- MUf[N1, ]
  scov[, , N1] <- CVf[, , N1]

  DMU <- matrix(NA_real_, N, l)
  DRl <- vector("list", N)
  OLmu <- matrix(NA_real_, N, l)
  OLR <- vector("list", N)
  s_n <- r_n <- e_n <- numeric(N)
  for (n in seq_len(N)) {
    co <- .compiled_co(bundle, n)
    aux <- .thmD1_aux_compiled(bundle, n, CVf[, , n], co = co)
    Ehist[[n]] <- aux$E
    ol <- .onelag_stats(
      co, aux, MUf[n, ], CVf[, , n], MUf[n + 1L, ], CVf[, , n + 1L],
      bundle$x[n + 1L, ] - bundle$x[n, ], dt, l
    )
    OLmu[n, ] <- ol$mu
    OLR[[n]] <- ol$R
    DMU[n, ] <- ol$mu - MUf[n, ]
    DRl[[n]] <- ol$R - CVf[, , n]
    s_n[n] <- sqrt(sum(DMU[n, ]^2))
    r_n[n] <- sqrt(sum(DRl[[n]]^2))
    e_n[n] <- if (l == 1L) abs(aux$E[1, 1]) else norm(aux$E, "2")
  }
  T2 <- Ub <- numeric(N + 1L)
  for (n in N:1) {
    eh2 <- max(1, e_n[n])^2
    T2[n] <- if (n < N) s_n[n + 1L]^2 + eh2 * T2[n + 1L] else 0
    Ub[n] <- if (n < N) r_n[n + 1L] + eh2 * Ub[n + 1L] else 0
    if (!is.finite(T2[n]) || T2[n] > 1e12) T2[n] <- Inf
    if (!is.finite(Ub[n]) || Ub[n] > 1e12) Ub[n] <- Inf
  }
  lam_j <- rep(NA_real_, N1)
  sch <- vector("list", N1)
  if (!is.null(smoo)) for (j in seq_len(N1)) {
    cs <- .chol_or_floor(smoo$cov[, , j])
    sch[[j]] <- cs$ch
    lam_j[j] <- min(eigen(cs$R, symmetric = TRUE, only.values = TRUE)$values)
  }

  for (n in seq_len(N)) {
    dmu <- DMU[n, ]
    dR <- DRl[[n]]
    ol <- list(mu = OLmu[n, ], R = OLR[[n]])
    if (!is.null(rows) && !is.null(smoo))
      rows[[n]] <- c(
        rows[[n]], .kl_fast(smoo$mean[n, ], sch[[n]], ol$mu, ol$R)
      )
    for (jj in which(active)) {
      laglen <- (n + 1L) - jj
      if (mode != "smoother_only" && is.finite(max_lag) &&
          laglen > max_lag) {
        frozen[jj] <- TRUE
        L[jj] <- as.integer(max_lag)
        tailb[jj] <- utils::tail(rows[[jj]], 1L)
        active[jj] <- FALSE
        act_D[jj] <- list(NULL)
        act_mu[jj] <- list(NULL)
        act_R[jj] <- list(NULL)
        next
      }
      act_D[[jj]] <- act_D[[jj]] %*% Ehist[[n - 1L]]
      new_mu <- act_mu[[jj]] + drop(act_D[[jj]] %*% dmu)
      new_R <- sym(act_R[[jj]] + act_D[[jj]] %*% dR %*% t(act_D[[jj]]))
      act_mu[[jj]] <- new_mu
      act_R[[jj]] <- new_R
      if (frozen[jj]) next
      Pval <- if (!is.null(smoo)) {
        .kl_fast(smoo$mean[jj, ], sch[[jj]], new_mu, new_R)
      } else {
        NA_real_
      }
      if (!is.null(rows) && !is.na(Pval)) rows[[jj]] <- c(rows[[jj]], Pval)
      if (!is.na(Pval) && tol > 0) {
        Dn <- if (l == 1L) {
          abs(act_D[[jj]][1, 1])
        } else {
          sqrt(sum(act_D[[jj]]^2))
        }
        fut <- 1.5 * (
          Dn^2 * T2[n] / (2 * lam_j[jj]) +
            Dn^2 * Ub[n] * sqrt(l) / (2 * lam_j[jj])
        )
        ok <- is.finite(fut) && (Pval + fut) < tol
        act_cnt[jj] <- if (ok) act_cnt[jj] + 1L else 0L
        if (act_cnt[jj] >= window || laglen >= max_lag) {
          frozen[jj] <- TRUE
          L[jj] <- laglen
          tailb[jj] <- if (is.finite(fut)) Pval + fut else Pval
          if (mode != "smoother_only") {
            active[jj] <- FALSE
            act_D[jj] <- list(NULL)
            act_mu[jj] <- list(NULL)
            act_R[jj] <- list(NULL)
          }
        }
      } else if (laglen >= max_lag) {
        frozen[jj] <- TRUE
        L[jj] <- laglen
        tailb[jj] <- Pval
        if (mode != "smoother_only") {
          active[jj] <- FALSE
          act_D[jj] <- list(NULL)
          act_mu[jj] <- list(NULL)
          act_R[jj] <- list(NULL)
        }
      }
    }
    active[n] <- TRUE
    act_D[[n]] <- diag(l)
    act_mu[[n]] <- ol$mu
    act_R[[n]] <- ol$R
    act_cnt[n] <- 0L
  }
  for (jj in which(active & !frozen)) {
    L[jj] <- N1 - jj
    tailb[jj] <- 0
  }
  if (mode == "smoother_only") for (jj in which(active)) {
    smu[jj, ] <- act_mu[[jj]]
    scov[, , jj] <- act_R[[jj]]
  }
  L[N1] <- 0L
  list(
    diag = diagv, rows = rows, L = L, tailbnd = tailb,
    onelag = NULL, stop_index = NA_integer_, smu = smu, scov = scov
  )
}


#' Assemble a lag table entirely from a compiled CGNS bundle (internal)
#'
#' Public routing remains unchanged. This function provides the complete private
#' integration seam, including supplied-path validation, Theorem-3 reference
#' policy, warnings and the established `lag_table` result structure.
#'
#' @param bundle A compiled CGNS bundle.
#' @param mode One of `"forward"`, `"one_lag"` or `"full"`.
#' @param tol,window,max_lag Adaptive storage controls matching [lag_table()].
#' @param filter Optional compatible filter path.
#' @param smoother Optional compatible Theorem-3 smoother path.
#' @param init Optional Gaussian initialization.
#' @param validate Validate the bundle and tuning arguments.
#' @returns A `lag_table` object.
#' @noRd
.lag_table_compiled <- function(
    bundle, mode = c("forward", "one_lag", "full"),
    tol = getOption("aci.default_tol", 1e-8), window = 3L,
    max_lag = Inf, filter = NULL, smoother = NULL, init = NULL,
    validate = TRUE) {
  mode <- match.arg(mode)
  if (isTRUE(validate))
    .validate_compiled_cgns(
      bundle, nontarget = bundle$nontarget, scalar = FALSE
    )
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    aci_abort("aci_error_dims", "tol must be one finite non-negative number.")
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window < 1L || window != floor(window))
    aci_abort("aci_error_dims", "window must be a positive integer.")
  if (!is.numeric(max_lag) || length(max_lag) != 1L || is.na(max_lag) ||
      max_lag < 1 || (!is.infinite(max_lag) && max_lag != floor(max_lag)))
    aci_abort("aci_error_dims", "max_lag must be a positive integer or Inf.")
  if (mode == "full" && is.finite(max_lag))
    aci_abort("aci_error_dims", "mode = 'full' requires max_lag = Inf.")

  filt <- filter
  smoo <- smoother
  if (!is.null(filt))
    .validate_gaussian_path(
      filt, bundle$obs, bundle$l, "filter", bundle$nontarget,
      model = bundle$model, source_model = bundle$source_model
    )
  if (!is.null(smoo))
    .validate_gaussian_path(
      smoo, bundle$obs, bundle$l, "smoother", bundle$nontarget,
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
    smoo <- NULL
  }
  filt <- filt %||% .cgns_filter_compiled(
    bundle, init = init, stepper = "explicit", nsub = 1L, validate = FALSE
  )
  filt$meta$source_model <- bundle$source_model

  supplied_smoo <- smoo
  if (!is.null(supplied_smoo) &&
      !identical(supplied_smoo$meta$route %||% NULL, "thmD1")) {
    aci_warn(
      "aci_warn_stepper",
      paste(
        "lag_table uses the complete Theorem 3 online smoother as its",
        "finite-step reference; the supplied backward-ODE smoother is",
        "incompatible and is being recomputed."
      )
    )
  }
  smoo <- .smoother_thmD1_compiled(
    bundle, filt, validate = FALSE, warn_cost = FALSE
  )
  if (!is.null(supplied_smoo) &&
      identical(supplied_smoo$meta$route %||% NULL, "thmD1")) {
    scale <- max(
      1,
      max(abs(c(
        smoo$mean, smoo$cov, supplied_smoo$mean, supplied_smoo$cov
      )))
    )
    if (max(abs(smoo$mean - supplied_smoo$mean)) > 1e-10 * scale ||
        max(abs(smoo$cov - supplied_smoo$cov)) > 1e-10 * scale)
      aci_abort(
        "aci_error_model_contract",
        paste(
          "The supplied Theorem 3 smoother was not generated from the",
          "same filter/prior as this lag table."
        )
      )
  }
  smoo$meta$source_model <- bundle$source_model
  eff_tol <- if (mode == "full") 0 else tol
  eff_win <- if (mode == "full") Inf else as.integer(window)
  eff_max_lag <- if (mode == "full") Inf else max_lag
  res <- .lagtable_core_compiled(
    bundle, filt, smoo, mode = mode, tol = eff_tol,
    window = eff_win, max_lag = eff_max_lag
  )
  dec <- .gaussian_kl_path_compiled(
    bundle, smoo, filt, decompose = TRUE, validate = FALSE
  )
  structure(
    list(
      t = bundle$t,
      dt = bundle$dt,
      mode = mode,
      diag = res$diag,
      rows = res$rows,
      L = res$L,
      diag_signal = dec$signal,
      diag_dispersion = dec$dispersion,
      tailbnd = res$tailbnd,
      onelag = res$onelag,
      meta = list(
        nontarget = bundle$nontarget,
        tol = eff_tol,
        window = eff_win,
        max_lag = eff_max_lag,
        init = filt$meta$init,
        source_model = bundle$source_model,
        source_obs_x = bundle$source_obs$x,
        reference_smoother = "thmD1_online_complete",
        stop_index = res$stop_index
      )
    ),
    class = "lag_table"
  )
}
