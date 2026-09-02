################################################################################
## aci-kernels-matrix.R - complete compiled CGNS matrix execution
################################################################################


#' Read one coefficient record from a compiled CGNS bundle (internal)
#'
#' @param bundle A validated `compiled_cgns` object.
#' @param j One-based grid index.
#' @returns A coefficient list with the same shapes as `eval_coefs()`.
#' @noRd
.compiled_co <- function(bundle, j) {
  co <- bundle$coefficients
  list(
    Lx = matrix(co$Lx[, , j], bundle$k, bundle$l),
    fx = as.numeric(co$fx[j, ]),
    Ly = matrix(co$Ly[, , j], bundle$l, bundle$l),
    fy = as.numeric(co$fy[j, ]),
    gxx = matrix(co$gxx[, , j], bundle$k, bundle$k),
    gyy = matrix(co$gyy[, , j], bundle$l, bundle$l),
    gyx = matrix(co$gyx[, , j], bundle$l, bundle$k)
  )
}


#' Read one realised observation-Gram inverse (internal)
#'
#' @param bundle A validated `compiled_cgns` object.
#' @param j One-based interval index.
#' @returns A `k` by `k` numeric matrix.
#' @noRd
.compiled_ginv <- function(bundle, j) {
  matrix(bundle$coefficients$gxx_weight[, , j], bundle$k, bundle$k)
}


#' Prepare the initial state for a compiled matrix filter (internal)
#'
#' @param bundle A validated `compiled_cgns` object.
#' @param init Optional Gaussian prior.
#' @returns Initial mean/covariance and their provenance.
#' @noRd
.compiled_filter_init <- function(bundle, init) {
  co0 <- .compiled_co(bundle, 1L)
  used_default_cov <- is.null(init) || is.null(init$cov)
  if (used_default_cov) {
    di <- .default_init(bundle$model, co0)
    if (is.null(init)) init <- di else init$cov <- di$cov
    aci_warn(
      "aci_warn_diffuse_init",
      paste(
        "No init$cov supplied; using a diffuse prior. Discard an initial",
        "burn-in window when interpreting results."
      )
    )
  }
  mu <- as.numeric(init$mean %||% rep(0, bundle$l))
  R <- as.matrix(init$cov)
  if (length(mu) != bundle$l ||
      !identical(dim(R), c(bundle$l, bundle$l)) ||
      !is.numeric(R) || any(!is.finite(c(mu, R))))
    aci_abort(
      "aci_error_dims",
      "init$mean/init$cov have incompatible dimensions or values."
    )
  if (max(abs(R - t(R))) > 1e-12 * max(1, max(abs(R))))
    aci_abort("aci_error_spd", "init$cov must be symmetric.")
  R <- sym(R)
  if (is.null(tryCatch(chol(R), error = function(e) NULL)))
    aci_abort(
      "aci_error_spd",
      paste(
        "init$cov must be positive definite; it is not regularized",
        "implicitly."
      )
    )
  list(
    mu = mu,
    R = R,
    actual = list(
      mean = mu,
      cov = R,
      source = if (used_default_cov) "default_diffuse" else "user"
    )
  )
}


#' Compiled matrix conditional-Gaussian filter (internal)
#'
#' This is the matrix/implicit counterpart of the specialised explicit scalar
#' kernel. It preserves the established matrix equations while consuming
#' coefficient arrays and Gram inverses that have already been realised.
#'
#' @param bundle A compiled CGNS bundle.
#' @param init Optional Gaussian prior.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of substeps.
#' @param validate Validate the compiled bundle before execution.
#' @param loglik Accumulate the predictive log-likelihood. `FALSE` skips the
#'   per-step innovation Cholesky and triangular solve entirely and leaves
#'   `meta$loglik` unset; the state recursion is untouched.
#' @param regularize Covariance policy, or a recorder shared with the caller;
#'   see `.aci_reg_for()`.
#' @returns A `da_path_gaussian` filter.
#' @noRd
.cgns_filter_matrix_compiled <- function(
    bundle, init = NULL, stepper = c("explicit", "implicit"), nsub = 1L,
    validate = TRUE, loglik = TRUE, regularize = NULL) {
  stepper <- match.arg(stepper)
  rec <- .aci_reg_for(regularize, bundle$t)
  if (isTRUE(validate))
    .validate_compiled_cgns(
      bundle, conditional = bundle$conditional, scalar = FALSE
    )
  if (length(nsub) != 1L || !is.finite(nsub) || nsub < 1 ||
      nsub != as.integer(nsub))
    aci_abort("aci_error_dims", "nsub must be a positive integer.")
  nsub <- as.integer(nsub)

  N1 <- bundle$N1
  N <- bundle$N
  dt <- bundle$dt
  l <- bundle$l
  h <- dt / nsub
  ini <- .compiled_filter_init(bundle, init)
  mu <- ini$mu
  R <- ini$R
  co0 <- .compiled_co(bundle, 1L)
  if (dt * max(abs(as.matrix(co0$Ly))) > 0.1 && nsub == 1L)
    aci_warn(
      "aci_warn_dt_stability",
      "dt * ||Ly|| > 0.1 at t0; Euler updates may be inaccurate (consider nsub > 1)."
    )

  ia <- bundle$likelihood_idx
  if (!is.numeric(ia) || !length(ia) || any(!is.finite(ia)) ||
      any(ia != floor(ia)) || any(ia < 1L) || any(ia > bundle$k) ||
      anyDuplicated(ia))
    aci_abort(
      "aci_error_dims",
      "likelihood_idx contains invalid observation indices."
    )
  ia <- as.integer(ia)

  MU <- matrix(NA_real_, N1, l)
  CV <- array(NA_real_, c(l, l, N1))
  MU[1L, ] <- mu
  CV[, , 1L] <- R
  Il <- diag(l)
  do_ll <- isTRUE(loglik)
  ll <- if (do_ll) 0 else NULL
  kll <- length(ia)
  stab <- 0

  ## Hot-loop bindings. The coefficient arrays, the rate matrix and the
  ## realised Gram inverses are sliced directly at every index: `.compiled_co()`
  ## builds a seven-element list of `matrix()` copies per step, which the
  ## profile puts at about a quarter of this kernel, and the smoother never
  ## needs `gxx` at all. Setting `dim` on the freshly sliced temporary restores
  ## the shape that `drop` collapses when `k` or `l` is 1, so the slices are the
  ## same matrices `.compiled_co()` returned.
  ##
  ## Inside the recursions `t.default()` and `chol.default()` stand in for the
  ## generics. Every argument here is a plain numeric matrix with no class, so
  ## the methods dispatched to are these, and the results are identical; on this
  ## machine the `UseMethod()` round trip is about 0.5 us per transpose and
  ## 0.6 us per factorisation, which four transposes and a factorisation per
  ## step turn into a measurable share of the kernel.
  k <- bundle$k
  cf <- bundle$coefficients
  a_Lx <- cf$Lx; a_Ly <- cf$Ly; a_gyy <- cf$gyy; a_gyx <- cf$gyx
  a_gxx <- cf$gxx; a_gw <- cf$gxx_weight
  m_fx <- cf$fx; m_fy <- cf$fy; m_rate <- bundle$rate
  d_kl <- c(k, l); d_ll <- c(l, l); d_lk <- c(l, k); d_kk <- c(k, k)
  d_diag <- seq.int(1L, l * l, by = l + 1L)
  d_ia <- seq.int(1L, kll * kll, by = kll + 1L)
  ia_all <- kll == k && identical(ia, seq_len(k))
  is_imp <- stepper == "implicit"

  for (j in seq_len(N)) {
    Lxj <- a_Lx[, , j]; dim(Lxj) <- d_kl
    Lyj <- a_Ly[, , j]; dim(Lyj) <- d_ll
    gyyj <- a_gyy[, , j]; dim(gyyj) <- d_ll
    gyxj <- a_gyx[, , j]; dim(gyxj) <- d_lk
    Gi <- a_gw[, , j]; dim(Gi) <- d_kk
    fxj <- as.numeric(m_fx[j, ])
    fyj <- as.numeric(m_fy[j, ])
    rate <- as.numeric(m_rate[j, ])
    tLx <- t.default(Lxj)
    tgyx <- t.default(gyxj)

    if (do_ll) {
      gxxj <- a_gxx[, , j]; dim(gxxj) <- d_kk
      Lxa <- if (ia_all) Lxj else Lxj[ia, , drop = FALSE]
      gxa <- if (ia_all) gxxj else gxxj[ia, ia, drop = FALSE]
      iota0 <- (rate[ia] - (fxj[ia] + drop(Lxa %*% mu))) * dt
      Sd <- sym(Lxa %*% R %*% t.default(Lxa) * dt + gxa) * dt
      ## The likelihood scores observation row j; the state update below
      ## writes the covariance stored at grid index j + 1.
      rec$j <- j
      cS <- .cov_guard_chol(Sd, rec, "likelihood_innov")
      w <- forwardsolve(t.default(cS$ch), iota0)
      ll <- ll - 0.5 * sum(w * w) - sum(log(cS$ch[d_ia])) -
        0.5 * kll * log(2 * pi)
    }

    if (is_imp) {
      GiLx <- Gi %*% Lxj
      Lt <- Lyj - gyxj %*% GiLx
      gt <- sym(gyyj - gyxj %*% Gi %*% tgyx)
      Ss <- tLx %*% GiLx
    } else {
      tLy <- t.default(Lyj)
      ## Only the diagonal of the Riccati stiffness matrix is scored, and it is
      ## a matrix already, so index it rather than routing through diag() and
      ## as.matrix().
      stab <- max(stab, max(abs((tLx %*% Gi %*% Lxj %*% R)[d_diag])) * h)
    }

    rec$j <- j + 1L
    for (ss in seq_len(nsub)) {
      Kg <- (R %*% tLx + gyxj) %*% Gi
      if (is_imp) {
        Am <- Lyj - Kg %*% Lxj
        bm <- fyj + drop(Kg %*% (rate - fxj))
        BEm <- Il - h * Am
        BEr <- Il - h * Lt
        if (rcond(BEm) < 1e-12 || rcond(BEr) < 1e-12)
          aci_abort(
            "aci_error_stepper",
            "Implicit filter step is singular; increase nsub or reduce dt."
          )
        mu <- drop(solve(BEm, mu + h * bm))
        Q <- .cov_guard(R + h * gt, rec, "filter_implicit_q")
        left <- solve(BEr, Q)
        Rp <- .cov_guard(t.default(solve(BEr, t.default(left))), rec,
                         "filter_implicit_p")
        info <- sym(
          chol_solve(Rp, Il, "implicit predicted covariance", rec,
                     "filter_implicit_p") + h * Ss
        )
        R <- .cov_guard(
          chol_solve(info, Il, "implicit information matrix", rec,
                     "filter_implicit_info"),
          rec, "filter_implicit_r"
        )
      } else {
        mu <- mu + (fyj + drop(Lyj %*% mu)) * h +
          drop(Kg %*% ((rate - (fxj + drop(Lxj %*% mu))) * h))
        R <- .cov_guard(
          R + (Lyj %*% R + R %*% tLy + gyyj - Kg %*% (Lxj %*% R + tgyx)) * h,
          rec, "filter_explicit"
        )
      }
    }
    MU[j + 1L, ] <- mu
    CV[, , j + 1L] <- R
  }

  if (stepper == "explicit" && stab > 1)
    aci_warn(
      "aci_warn_riccati_stiff",
      sprintf(
        paste(
          "Explicit Riccati step is unstable (max ||Lx' gxx^-1 Lx R|| dt = %.3g > 1):",
          "the covariance can overshoot, leave the positive-definite cone, and oscillate.",
          "Use the positivity-preserving implicit stepper, or reduce dt / increase nsub."
        ),
        stab
      )
    )

  p <- new_da_path(bundle$t, MU, CV, "filter")
  p$meta$stepper <- stepper
  p$meta$nsub <- nsub
  ## loglik = FALSE leaves meta$loglik unset (NULL); likelihood_idx still
  ## records which observation rows the likelihood contract would score.
  if (do_ll) p$meta$loglik <- ll
  p$meta$likelihood_idx <- ia
  p$meta$init <- ini$actual
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$conditional <- bundle$conditional
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p$meta$regularization <- .aci_reg_freeze(rec)
  p
}


#' Compiled CGNS filter dispatcher (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param init Optional Gaussian prior.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of substeps.
#' @param validate Validate the bundle before execution.
#' @param loglik Accumulate the predictive log-likelihood; `FALSE` skips that
#'   work on both the scalar and the matrix path.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A `da_path_gaussian` filter.
#' @noRd
.cgns_filter_compiled <- function(
    bundle, init = NULL, stepper = c("explicit", "implicit"), nsub = 1L,
    validate = TRUE, loglik = TRUE, regularize = NULL) {
  stepper <- match.arg(stepper)
  if (length(loglik) != 1L || is.na(loglik) || !is.logical(loglik))
    aci_abort("aci_error_dims", "loglik must be TRUE or FALSE.")
  if (stepper == "explicit" && bundle$k == 1L && bundle$l == 1L)
    return(.cgns_filter_scalar(
      bundle, init, nsub, validate = validate, loglik = loglik,
      regularize = regularize
    ))
  .cgns_filter_matrix_compiled(
    bundle, init, stepper = stepper, nsub = nsub, validate = validate,
    loglik = loglik, regularize = regularize
  )
}


#' Compiled matrix backward-ODE smoother (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible filter path.
#' @param validate Validate the bundle and externally supplied filter.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A `da_path_gaussian` smoother.
#' @noRd
.cgns_smoother_matrix_compiled <- function(bundle, filter, validate = TRUE,
                                           regularize = NULL) {
  rec <- .aci_reg_for(regularize, bundle$t)
  if (isTRUE(validate)) {
    .validate_compiled_cgns(
      bundle, conditional = bundle$conditional, scalar = FALSE
    )
    .validate_gaussian_path(
      filter, bundle$obs, bundle$l, "filter", bundle$conditional,
      model = bundle$model, source_model = bundle$source_model
    )
  }
  N1 <- bundle$N1
  l <- bundle$l
  h <- bundle$dt / max(1L, as.integer(filter$meta$nsub %||% 1L))
  nsub <- max(1L, as.integer(filter$meta$nsub %||% 1L))
  MU <- matrix(NA_real_, N1, l)
  CV <- array(NA_real_, c(l, l, N1))
  MU[N1, ] <- filter$mean[N1, ]
  CV[, , N1] <- filter$cov[, , N1]
  mus <- MU[N1, ]
  Rs <- CV[, , N1]
  Il <- diag(l)

  ## Hot-loop bindings, as in the filter. `gxx` is not part of the backward
  ## recursion and is not read. The filter covariance is factorised inline:
  ## `chol_solve()` re-validates a covariance this kernel has already accepted,
  ## and its jitter ladder is only reachable when the plain factorisation fails,
  ## which is where the full `chol_solve()` call is kept - now under the
  ## covariance policy, so the ladder is entered only when the caller asked for
  ## `regularize = "floor"`.  The inverse itself
  ## is taken with `chol2inv()` rather than two triangular solves against an
  ## identity: the same inverse to the last bit or two, exactly symmetric where
  ## the triangular route is not, and the recursion is budgeted for it.  It is
  ## bit-identical whenever `l == 1`, where the two routes are one division.
  k <- bundle$k
  cf <- bundle$coefficients
  a_Lx <- cf$Lx; a_Ly <- cf$Ly; a_gyy <- cf$gyy; a_gyx <- cf$gyx
  a_gw <- cf$gxx_weight
  m_fx <- cf$fx; m_fy <- cf$fy; m_rate <- bundle$rate
  d_kl <- c(k, l); d_ll <- c(l, l); d_lk <- c(l, k); d_kk <- c(k, k)
  f_mean <- filter$mean; f_cov <- filter$cov

  for (j in (N1 - 1L):1L) {
    rec$j <- j
    Lxj <- a_Lx[, , j]; dim(Lxj) <- d_kl
    Lyj <- a_Ly[, , j]; dim(Lyj) <- d_ll
    gyyj <- a_gyy[, , j]; dim(gyyj) <- d_ll
    gyxj <- a_gyx[, , j]; dim(gyxj) <- d_lk
    Gi <- a_gw[, , j]; dim(Gi) <- d_kk
    fxj <- as.numeric(m_fx[j, ])
    fyj <- as.numeric(m_fy[j, ])
    rate <- as.numeric(m_rate[j, ])
    tgyx <- t.default(gyxj)
    Cgi <- gyxj %*% Gi
    A0 <- Lyj - Cgi %*% Lxj
    B <- sym(gyyj - Cgi %*% tgyx)
    Rf <- f_cov[, , j]; dim(Rf) <- d_ll
    ## finiteness screened explicitly; LAPACK may complete on non-finite
    ## input (see .cov_guard), and a poisoned factor must fall through to the
    ## guarded route below
    chf <- if (all(is.finite(Rf)))
      tryCatch(chol.default((Rf + t.default(Rf)) / 2),
               error = function(e) NULL) else NULL
    Rfi <- if (is.null(chf))
      chol_solve(Rf, Il, "Rf", rec, "smoother_filter_cov") else chol2inv(chf)
    H <- A0 + B %*% Rfi
    tH <- t.default(H)
    fmj <- f_mean[j, ]
    for (ss in seq_len(nsub)) {
      d <- drop(Rfi %*% (fmj - mus))
      mus <- mus + h * (
        -fyj - drop(Lyj %*% mus) + drop(B %*% d) +
          drop(Cgi %*% (-rate + fxj + drop(Lxj %*% mus)))
      )
      Rs <- .cov_guard(Rs + h * (-(H %*% Rs) - Rs %*% tH + B), rec,
                       "smoother_backward")
    }
    MU[j, ] <- mus
    CV[, , j] <- Rs
  }

  p <- new_da_path(bundle$t, MU, CV, "smoother")
  p$meta$route <- if (bundle$correlated_noise)
    "backward_ode_correlated" else "backward_ode"
  p$meta$scheme <- "backward_ode_euler"
  p$meta$nsub <- nsub
  p$meta$init <- filter$meta$init
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$conditional <- bundle$conditional
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p$meta$regularization <- .aci_reg_freeze(rec)
  stopifnot(
    max(abs(p$mean[N1, ] - filter$mean[N1, ])) < 1e-12
  )
  p
}


#' Compiled CGNS smoother dispatcher (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible filter path.
#' @param validate Validate a user-supplied path.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A `da_path_gaussian` smoother.
#' @noRd
.cgns_smoother_compiled <- function(bundle, filter, validate = TRUE,
                                    regularize = NULL) {
  if (bundle$k == 1L && bundle$l == 1L)
    return(.cgns_smoother_scalar(bundle, filter, validate = validate,
                                 regularize = regularize))
  .cgns_smoother_matrix_compiled(bundle, filter, validate = validate,
                                 regularize = regularize)
}


#' Compiled CGNS Gaussian-KL path dispatcher (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param p,q Compatible Gaussian paths.
#' @param decompose Return signal and dispersion components.
#' @param validate Validate supplied paths.
#' @returns A Gaussian KL path data frame.
#' @noRd
.gaussian_kl_path_compiled <- function(
    bundle, p, q, decompose = TRUE, validate = TRUE) {
  if (bundle$k == 1L && bundle$l == 1L)
    return(.gaussian_kl_path_scalar(
      bundle, p, q, decompose = decompose, validate = validate
    ))
  if (isTRUE(validate)) {
    .validate_compiled_cgns(
      bundle, conditional = bundle$conditional, scalar = FALSE
    )
  }
  aci_metric(p, q, decompose = decompose)
}


#' Complete compiled closed-form ACI execution (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param init Optional Gaussian prior.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of substeps.
#' @param keep Retain paths or discard them after metric construction.
#' @param decompose Return signal and dispersion components.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns An `aci_result` without a lag table.
#' @noRd
.aci_cgns_compiled <- function(
    bundle, init = NULL, stepper = c("explicit", "implicit"), nsub = 1L,
    keep = c("paths", "none"), decompose = TRUE, regularize = NULL) {
  stepper <- match.arg(stepper)
  keep <- match.arg(keep)
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  rec <- .aci_reg_for(regularize, bundle$t)
  .validate_compiled_cgns(
    bundle, conditional = bundle$conditional, scalar = FALSE
  )
  filt <- .cgns_filter_compiled(
    bundle, init = init, stepper = stepper, nsub = nsub, validate = FALSE,
    regularize = rec
  )
  smoo <- .cgns_smoother_compiled(bundle, filt, validate = FALSE,
                                  regularize = rec)
  klp <- .gaussian_kl_path_compiled(
    bundle, smoo, filt, decompose = decompose, validate = FALSE
  )
  structure(
    list(
      t = klp$t,
      aci = klp$total,
      signal = if (decompose) klp$signal else NULL,
      dispersion = if (decompose) klp$dispersion else NULL,
      paths = if (keep == "paths")
        list(filter = filt, smoother = smoo) else NULL,
      table = NULL,
      handles = list(
        model = bundle$source_model,
        obs = bundle$source_obs,
        conditional = bundle$conditional,
        init = filt$meta$init
      ),
      meta = list(
        engine = "cgns",
        conditional = bundle$conditional,
        m = NULL,
        smoother_scheme = smoo$meta$scheme,
        table_reference = NULL,
        regularization = .aci_reg_freeze(rec)
      )
    ),
    class = "aci_result"
  )
}
