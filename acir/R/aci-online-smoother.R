################################################################################
## aci-online-smoother.R - compiled Theorem-3 smoother and lag-table execution
################################################################################


#' Theorem-3 auxiliaries from one compiled coefficient record (internal)
#'
#' @param bundle A validated `compiled_cgns` object.
#' @param j One-based interval index.
#' @param Rf Filtered covariance at the interval start.
#' @param co Optional coefficient record already read with `.compiled_co()`.
#' @param rec Covariance-policy recorder; its `j` names the interval index.
#' @returns A list with the auxiliary matrices `E` and `F`.
#' @noRd
.thmD1_aux_compiled <- function(bundle, j, Rf, co = NULL,
                                rec = .aci_reg_for(NULL, bundle$t)) {
  co <- co %||% .compiled_co(bundle, j)
  Gi <- .compiled_ginv(bundle, j)
  l <- bundle$l
  dt <- bundle$dt
  Rfi <- chol_solve(Rf, diag(l), "Rf", rec, "onelag_filter_cov")
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
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A `da_path_gaussian` smoother with route `"thmD1"`.
#' @noRd
.smoother_thmD1_compiled <- function(bundle, filt, validate = TRUE,
                                     warn_cost = TRUE, regularize = NULL) {
  rec <- .aci_reg_for(regularize, bundle$t)
  if (isTRUE(validate)) {
    .validate_compiled_cgns(
      bundle, conditional = bundle$conditional, scalar = FALSE
    )
    .validate_gaussian_path(
      filt, bundle$obs, bundle$l, "filter", bundle$conditional,
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
    rec$j <- j
    co <- .compiled_co(bundle, j)
    aux <- .thmD1_aux_compiled(bundle, j, filt$cov[, , j], co = co, rec = rec)
    st <- .onelag_stats(
      co, aux, filt$mean[j, ], filt$cov[, , j],
      MU[j + 1L, ], CV[, , j + 1L],
      bundle$x[j + 1L, ] - bundle$x[j, ], bundle$dt, l, rec
    )
    MU[j, ] <- st$mu
    CV[, , j] <- st$R
  }
  p <- new_da_path(bundle$t, MU, CV, "smoother")
  p$meta$route <- "thmD1"
  p$meta$scheme <- "theorem3_discrete"
  p$meta$nsub <- 1L
  p$meta$init <- filt$meta$init
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$conditional <- bundle$conditional
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p$meta$regularization <- .aci_reg_freeze(rec)
  p
}


#' Per-interval online-smoother auxiliaries from one forward pass (internal)
#'
#' One sweep over the intervals collecting, for each `n`, the Theorem 3 update
#' matrix `E_n` and the increment that admitting the observation at `n + 1`
#' adds to a step already informed to `n`. Every fixed-lag posterior on the
#' record is an affine function of these, so they are collected once and reused
#' for every requested lag.
#'
#' This is the loop `.lagtable_core_compiled()` also runs to fill its `Ehist`,
#' `DMU` and `DRl`; de-duplicating the two is a post-0.1.0 candidate, recorded
#' in the adoption ledger.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filt A compatible explicit single-step filter path.
#' @param rec Covariance-policy recorder; its `j` names the interval index.
#' @returns A list with `E` (`l` by `l` by `N`), `dmu` (`N` by `l`) and `dR`
#'   (`l` by `l` by `N`).
#' @noRd
.online_aux_compiled <- function(bundle, filt,
                                 rec = .aci_reg_for(NULL, bundle$t)) {
  N <- bundle$N
  l <- bundle$l
  dt <- bundle$dt
  MUf <- filt$mean
  CVf <- filt$cov
  E <- array(NA_real_, c(l, l, N))
  dmu <- matrix(NA_real_, N, l)
  dR <- array(NA_real_, c(l, l, N))
  for (n in seq_len(N)) {
    rec$j <- n
    co <- .compiled_co(bundle, n)
    aux <- .thmD1_aux_compiled(bundle, n, CVf[, , n], co = co, rec = rec)
    ol <- .onelag_stats(
      co, aux, MUf[n, ], CVf[, , n], MUf[n + 1L, ], CVf[, , n + 1L],
      bundle$x[n + 1L, ] - bundle$x[n, ], dt, l, rec
    )
    E[, , n] <- aux$E
    dmu[n, ] <- ol$mu - MUf[n, ]
    dR[, , n] <- ol$R - CVf[, , n]
  }
  list(E = E, dmu = dmu, dR = dR)
}


#' Compose two adjacent online-update ranges (internal)
#'
#' The contribution of a contiguous range of intervals is the triple
#' `(P, s, S)`: the ordered product of `E` over the range, the mean offset it
#' contributes, and the covariance offset. Composition of adjacent ranges is
#' associative, which is what lets a fixed-length window be aggregated in
#' amortized constant work per step instead of re-multiplying the whole window
#' at every anchor.
#'
#' @param a,b Range triples, `a` earlier than `b`.
#' @returns The triple for the concatenated range.
#' @noRd
.online_join <- function(a, b) {
  list(
    P = a$P %*% b$P,
    s = a$s + drop(a$P %*% b$s),
    S = a$S + a$P %*% b$S %*% t(a$P)
  )
}


#' Fixed-lag online moments at every anchor (internal)
#'
#' Slides a window of exactly `L` intervals across the record with a two-stack
#' aggregation queue, so the cost is O(N) range compositions whatever the lag.
#'
#' The accumulated covariance is symmetrized and not otherwise regularised, the
#' convention `.lagtable_core_compiled()` already uses for the same affine
#' accumulation. Every increment it sums has been through the covariance policy
#' inside `.onelag_stats()`.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filt The forward filter path.
#' @param aux Auxiliaries from `.online_aux_compiled()`; unused when `L <= 0`.
#' @param L Whole lag, already capped at `bundle$N`.
#' @returns A list with `mean`, `cov` and the per-anchor `lag_effective`.
#' @noRd
.online_window_compiled <- function(bundle, filt, aux, L) {
  N <- bundle$N
  l <- bundle$l
  MU <- filt$mean
  CV <- filt$cov
  eff <- integer(bundle$N1)
  if (L <= 0L) return(list(mean = MU, cov = CV, lag_effective = eff))
  ## front_agg[[i]] is the aggregate of front-queue indices i, i - 1, ..., 1 in
  ## queue order, with the oldest element on top at index nf; back holds newly
  ## pushed intervals with a running aggregate. Only the front's aggregates are
  ## ever read back, so the front keeps no element list of its own.
  front_agg <- vector("list", N)
  back_el <- vector("list", N)
  nf <- 0L
  nb <- 0L
  back_agg <- NULL
  m_next <- 1L
  for (j in seq_len(N)) {
    hi <- min(j + L - 1L, N)
    while (m_next <= hi) {
      el <- list(
        P = matrix(aux$E[, , m_next], l, l),
        s = aux$dmu[m_next, ],
        S = matrix(aux$dR[, , m_next], l, l)
      )
      nb <- nb + 1L
      back_el[[nb]] <- el
      back_agg <- if (nb == 1L) el else .online_join(back_agg, el)
      m_next <- m_next + 1L
    }
    agg <- if (nf > 0L && nb > 0L) {
      .online_join(front_agg[[nf]], back_agg)
    } else if (nf > 0L) {
      front_agg[[nf]]
    } else {
      back_agg
    }
    MU[j, ] <- MU[j, ] + agg$s
    CV[, , j] <- sym(CV[, , j] + agg$S)
    eff[j] <- hi - j + 1L
    if (nf == 0L) {
      suffix <- NULL
      for (i in seq.int(nb, 1L)) {
        el <- back_el[[i]]
        suffix <- if (i == nb) el else .online_join(el, suffix)
        nf <- nf + 1L
        front_agg[[nf]] <- suffix
      }
      nb <- 0L
      back_agg <- NULL
    }
    nf <- nf - 1L
  }
  list(mean = MU, cov = CV, lag_effective = eff)
}


#' Online moments at requested anchor and lag pairs (internal)
#'
#' Accumulates each anchor's affine update directly, snapshotting the requested
#' lags on the way past. Cost is proportional to the summed requested lags, so
#' this is the cheap route for a sparse grid of pairs; `.online_window_compiled()`
#' is the route for one lag at every anchor.
#'
#' @param bundle A compiled CGNS bundle.
#' @param filt The forward filter path.
#' @param aux Auxiliaries from `.online_aux_compiled()`.
#' @param j Integer vector of anchor indices.
#' @param lag Numeric vector of requested lags, recycled from length one.
#' @returns A list with `mean` (one row per pair), `cov` (one slice per pair)
#'   and the realised `lag_effective`.
#' @noRd
.online_at_compiled <- function(bundle, filt, aux, j, lag) {
  N1 <- bundle$N1
  l <- bundle$l
  j <- as.integer(j)
  np <- length(j)
  if (length(lag) == 1L) lag <- rep(lag, np)
  cap <- pmin(lag, N1 - j)
  mu_out <- matrix(NA_real_, np, l)
  cov_out <- array(NA_real_, c(l, l, np))
  eff <- integer(np)
  Il <- diag(l)
  ord <- order(j, cap)
  i <- 1L
  while (i <= np) {
    jj <- j[ord[i]]
    grp <- i
    while (grp < np && j[ord[grp + 1L]] == jj) grp <- grp + 1L
    idx <- ord[i:grp]
    want <- cap[idx]
    D <- Il
    mu <- filt$mean[jj, ]
    R <- matrix(filt$cov[, , jj], l, l)
    cnt <- 0L
    for (h in which(want == 0L)) {
      mu_out[idx[h], ] <- mu
      cov_out[, , idx[h]] <- R
    }
    Lmax <- max(want)
    if (Lmax > 0L) for (n in seq.int(jj, jj + Lmax - 1L)) {
      mu <- mu + drop(D %*% aux$dmu[n, ])
      R <- sym(R + D %*% matrix(aux$dR[, , n], l, l) %*% t(D))
      cnt <- cnt + 1L
      for (h in which(want == cnt)) {
        mu_out[idx[h], ] <- mu
        cov_out[, , idx[h]] <- R
        eff[idx[h]] <- cnt
      }
      D <- D %*% matrix(aux$E[, , n], l, l)
    }
    i <- grp + 1L
  }
  list(mean = mu_out, cov = cov_out, lag_effective = eff)
}


#' Assemble a fixed-lag online path from a compiled bundle (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filt A compatible explicit single-step filter path.
#' @param lag Non-negative whole number, or `Inf`.
#' @param method `"auto"`, or one of the two routes for testing.
#' @param aux Optional auxiliaries, to reuse one forward pass across lags.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A `da_path_gaussian` of kind `"online"`.
#' @noRd
.da_online_compiled <- function(bundle, filt, lag,
                                method = c("auto", "window", "backward"),
                                aux = NULL, regularize = NULL) {
  method <- match.arg(method)
  rec <- .aci_reg_for(regularize, bundle$t)
  if (!identical(filt$meta$stepper %||% "explicit", "explicit") ||
      (filt$meta$nsub %||% 1L) != 1L)
    aci_abort(
      "aci_error_stepper",
      "The online Theorem 3 smoother requires an explicit single-step filter."
    )
  N1 <- bundle$N1
  N <- bundle$N
  saturating <- is.infinite(lag) || lag >= N
  if (method == "auto")
    method <- if (saturating) "backward" else "window"
  if (identical(method, "backward")) {
    p <- .smoother_thmD1_compiled(bundle, filt, validate = FALSE,
                                  warn_cost = FALSE, regularize = rec)
    MU <- p$mean
    CV <- p$cov
    eff <- as.integer(pmin(if (is.infinite(lag)) N else lag,
                           N1 - seq_len(N1)))
    route <- "thmD1_backward"
  } else {
    L <- if (is.infinite(lag)) N else as.integer(lag)
    ## A lag of zero is the filter itself, and the window route returns its
    ## moments untouched, so the O(N) auxiliary pass is not run for it.
    if (L > 0L) aux <- aux %||% .online_aux_compiled(bundle, filt, rec)
    res <- .online_window_compiled(bundle, filt, aux, L)
    MU <- res$mean
    CV <- res$cov
    eff <- res$lag_effective
    route <- "thmD1_online_window"
  }
  p <- new_da_path(bundle$t, MU, CV, "online")
  p$meta$route <- route
  p$meta$scheme <- "theorem3_discrete"
  p$meta$lag <- lag
  p$meta$lag_effective <- eff
  p$meta$saturated <- all(eff >= (N1 - seq_len(N1)))
  p$meta$stepper <- "explicit"
  p$meta$nsub <- 1L
  p$meta$init <- filt$meta$init
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$conditional <- bundle$conditional
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p$meta$regularization <- .aci_reg_freeze(rec)
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
#' @param mode Either `"forward"` or `"full"`.
#' @param tol Non-negative adaptive storage tolerance.
#' @param window Consecutive steps below `tol` required before freezing.
#' @param max_lag Maximum retained positive lag.
#' @param rec Covariance-policy recorder from `.aci_reg_new()`.
#' @returns The internal retained-row and storage-diagnostic result list.
#' @noRd
.lagtable_core_compiled <- function(bundle, filt, smoo, mode, tol, window,
                                    max_lag,
                                    rec = .aci_reg_for(NULL, bundle$t)) {
  N1 <- bundle$N1
  N <- bundle$N
  dt <- bundle$dt
  l <- bundle$l
  MUf <- filt$mean
  CVf <- filt$cov
  diagv <- rep(NA_real_, N1)
  if (!is.null(smoo)) for (j in seq_len(N1))
    diagv[j] <- unname(aci_metric_pair(
      smoo$mean[j, ], smoo$cov[, , j],
      MUf[j, ], CVf[, , j], decompose = FALSE
    ))

  rows <- vector("list", N1)
  for (j in seq_len(N1)) rows[[j]] <- diagv[j]
  L <- rep(NA_integer_, N1)
  tailb <- rep(0, N1)
  act_D <- vector("list", N1)
  act_mu <- vector("list", N1)
  act_R <- vector("list", N1)
  act_cnt <- integer(N1)
  active <- frozen <- logical(N1)
  Ehist <- vector("list", N1)

  DMU <- matrix(NA_real_, N, l)
  DRl <- vector("list", N)
  OLmu <- matrix(NA_real_, N, l)
  OLR <- vector("list", N)
  s_n <- r_n <- e_n <- numeric(N)
  for (n in seq_len(N)) {
    rec$j <- n
    co <- .compiled_co(bundle, n)
    aux <- .thmD1_aux_compiled(bundle, n, CVf[, , n], co = co, rec = rec)
    Ehist[[n]] <- aux$E
    ol <- .onelag_stats(
      co, aux, MUf[n, ], CVf[, , n], MUf[n + 1L, ], CVf[, , n + 1L],
      bundle$x[n + 1L, ] - bundle$x[n, ], dt, l, rec
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
  lam_j <- Rs_j <- rep(NA_real_, N1)
  sch <- vector("list", N1)
  if (!is.null(smoo)) for (j in seq_len(N1)) {
    rec$j <- j
    cs <- .cov_guard_chol(smoo$cov[, , j], rec, "metric_reference")
    sch[[j]] <- cs$ch
    Rs_j[j] <- cs$R[1L]
    lam_j[j] <- min(eigen(cs$R, symmetric = TRUE, only.values = TRUE)$values)
  }

  if (l == 1L && !is.null(smoo)) {
    ## A scalar hidden state forms each anchor's row as one vector expression
    ## (`.cir_scalar_row()`) instead of advancing every active anchor one cell
    ## per step; the retained cells, freeze index and tail estimate are the
    ## ones the recursion below records.
    sp <- .cir_scalar_primitives(Ehist, DMU, DRl)
    for (j in seq_len(N)) {
      rec$j <- j
      r1 <- .cir_scalar_row(
        j, sp, T2, Ub, N, MUf[j, 1L], CVf[1L, 1L, j], smoo$mean[j, 1L],
        Rs_j[j], diagv[j], lam_j[j], tol, window, max_lag, rec
      )
      rows[[j]] <- r1$row
      L[j] <- r1$lag_kept
      tailb[j] <- r1$tail
    }
    L[N1] <- 0L
    return(list(
      diag = diagv, rows = rows, L = L, tailbnd = tailb,
      onelag = NULL, stop_index = NA_integer_
    ))
  }

  for (n in seq_len(N)) {
    dmu <- DMU[n, ]
    dR <- DRl[[n]]
    ol <- list(mu = OLmu[n, ], R = OLR[[n]])
    if (!is.null(rows) && !is.null(smoo)) {
      rec$j <- n
      rows[[n]] <- c(
        rows[[n]], .kl_fast(smoo$mean[n, ], sch[[n]], ol$mu, ol$R, rec)
      )
    }
    for (jj in which(active)) {
      laglen <- (n + 1L) - jj
      if (is.finite(max_lag) && laglen > max_lag) {
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
        rec$j <- jj
        .kl_fast(smoo$mean[jj, ], sch[[jj]], new_mu, new_R, rec)
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
          active[jj] <- FALSE
          act_D[jj] <- list(NULL)
          act_mu[jj] <- list(NULL)
          act_R[jj] <- list(NULL)
        }
      } else if (laglen >= max_lag) {
        frozen[jj] <- TRUE
        L[jj] <- laglen
        tailb[jj] <- Pval
        active[jj] <- FALSE
        act_D[jj] <- list(NULL)
        act_mu[jj] <- list(NULL)
        act_R[jj] <- list(NULL)
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
  L[N1] <- 0L
  list(
    diag = diagv, rows = rows, L = L, tailbnd = tailb,
    ## Always NULL on the mainline: the reserved fbcir one-lag mode is the
    ## field's producer (its lt_onelag() accessor reads it), so the slot
    ## stays as the family's staging point. See dev/reserve/fbcir/.
    onelag = NULL, stop_index = NA_integer_
  )
}


#' Assemble a lag table entirely from a compiled CGNS bundle (internal)
#'
#' Public routing remains unchanged. This function provides the complete private
#' integration seam, including supplied-path validation, Theorem-3 reference
#' policy, warnings and the established `lag_table` result structure.
#'
#' @param bundle A compiled CGNS bundle.
#' @param mode Either `"forward"` or `"full"`.
#' @param tol,window,max_lag Adaptive storage controls matching [lag_table()].
#' @param filter Optional compatible filter path.
#' @param smoother Optional compatible Theorem-3 smoother path.
#' @param init Optional Gaussian initialization.
#' @param validate Validate the bundle and tuning arguments.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A `lag_table` object.
#' @noRd
.lag_table_compiled <- function(
    bundle, mode = c("forward", "full"),
    tol = getOption("aci.default_tol", 1e-8), window = 3L,
    max_lag = Inf, filter = NULL, smoother = NULL, init = NULL,
    validate = TRUE, regularize = NULL) {
  mode <- match.arg(mode)
  rec <- .aci_reg_for(regularize, bundle$t)
  if (isTRUE(validate))
    .validate_compiled_cgns(
      bundle, conditional = bundle$conditional, scalar = FALSE
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
      filt, bundle$obs, bundle$l, "filter", bundle$conditional,
      model = bundle$model, source_model = bundle$source_model
    )
  if (!is.null(smoo))
    .validate_gaussian_path(
      smoo, bundle$obs, bundle$l, "smoother", bundle$conditional,
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
    bundle, init = init, stepper = "explicit", nsub = 1L, validate = FALSE,
    regularize = rec
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
    bundle, filt, validate = FALSE, warn_cost = FALSE, regularize = rec
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
    window = eff_win, max_lag = eff_max_lag, rec = rec
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
        conditional = bundle$conditional,
        tol = eff_tol,
        window = eff_win,
        max_lag = eff_max_lag,
        init = filt$meta$init,
        source_model = bundle$source_model,
        source_obs_x = bundle$source_obs$x,
        reference_smoother = "thmD1_online_complete",
        scheme = "theorem3_discrete",
        stop_index = res$stop_index,
        regularization = .aci_reg_freeze(rec)
      )
    ),
    class = "lag_table"
  )
}
