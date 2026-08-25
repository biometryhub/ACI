################################################################################
## compiled_matrix.R - complete compiled CGNS matrix execution
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
#' @returns A `da_path_gaussian` filter.
#' @noRd
.cgns_filter_matrix_compiled <- function(
    bundle, init = NULL, stepper = c("explicit", "implicit"), nsub = 1L,
    validate = TRUE) {
  stepper <- match.arg(stepper)
  if (isTRUE(validate))
    .validate_compiled_cgns(
      bundle, nontarget = bundle$nontarget, scalar = FALSE
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
  ll <- 0
  kll <- length(ia)
  stab <- 0

  for (j in seq_len(N)) {
    co <- .compiled_co(bundle, j)
    Gi <- .compiled_ginv(bundle, j)
    rate <- as.numeric(bundle$rate[j, ])

    iota0 <- (rate[ia] - (co$fx[ia] +
      drop(co$Lx[ia, , drop = FALSE] %*% mu))) * dt
    Lxa <- co$Lx[ia, , drop = FALSE]
    Sd <- sym(
      Lxa %*% R %*% t(Lxa) * dt +
        co$gxx[ia, ia, drop = FALSE]
    ) * dt
    cS <- .chol_or_floor(Sd)
    w <- forwardsolve(t(cS$ch), iota0)
    ll <- ll - 0.5 * sum(w * w) - sum(log(diag(cS$ch))) -
      0.5 * kll * log(2 * pi)

    if (stepper == "implicit") {
      GiLx <- Gi %*% co$Lx
      Lt <- co$Ly - co$gyx %*% GiLx
      gt <- sym(co$gyy - co$gyx %*% Gi %*% t(co$gyx))
      Ss <- t(co$Lx) %*% GiLx
    } else {
      stab <- max(
        stab,
        max(abs(diag(as.matrix(t(co$Lx) %*% Gi %*% co$Lx %*% R)))) * h
      )
    }

    for (ss in seq_len(nsub)) {
      Kg <- (R %*% t(co$Lx) + co$gyx) %*% Gi
      if (stepper == "implicit") {
        Am <- co$Ly - Kg %*% co$Lx
        bm <- co$fy + drop(Kg %*% (rate - co$fx))
        BEm <- Il - h * Am
        BEr <- Il - h * Lt
        if (rcond(BEm) < 1e-12 || rcond(BEr) < 1e-12)
          aci_abort(
            "aci_error_stepper",
            "Implicit filter step is singular; increase nsub or reduce dt."
          )
        mu <- drop(solve(BEm, mu + h * bm))
        Q <- spd_floor(sym(R + h * gt))
        left <- solve(BEr, Q)
        Rp <- spd_floor(sym(t(solve(BEr, t(left)))))
        info <- sym(
          chol_solve(Rp, Il, "implicit predicted covariance") + h * Ss
        )
        R <- spd_floor(
          sym(chol_solve(info, Il, "implicit information matrix"))
        )
      } else {
        mu <- mu + (co$fy + drop(co$Ly %*% mu)) * h +
          drop(Kg %*% ((rate - (co$fx + drop(co$Lx %*% mu))) * h))
        R <- R + (
          co$Ly %*% R + R %*% t(co$Ly) + co$gyy -
            Kg %*% (co$Lx %*% R + t(co$gyx))
        ) * h
        R <- spd_floor(sym(R))
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
          "the covariance can overshoot, be floored by spd_floor(), and oscillate.",
          "Use the positivity-preserving implicit stepper, or reduce dt / increase nsub."
        ),
        stab
      )
    )

  p <- new_da_path(bundle$t, MU, CV, "filter")
  p$meta$stepper <- stepper
  p$meta$nsub <- nsub
  p$meta$loglik <- ll
  p$meta$likelihood_idx <- ia
  p$meta$init <- ini$actual
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$nontarget <- bundle$nontarget
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p
}


#' Compiled CGNS filter dispatcher (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param init Optional Gaussian prior.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of substeps.
#' @param validate Validate the bundle before execution.
#' @returns A `da_path_gaussian` filter.
#' @noRd
.cgns_filter_compiled <- function(
    bundle, init = NULL, stepper = c("explicit", "implicit"), nsub = 1L,
    validate = TRUE) {
  stepper <- match.arg(stepper)
  if (stepper == "explicit" && bundle$k == 1L && bundle$l == 1L)
    return(.cgns_filter_scalar(bundle, init, nsub, validate = validate))
  .cgns_filter_matrix_compiled(
    bundle, init, stepper = stepper, nsub = nsub, validate = validate
  )
}


#' Compiled matrix backward-ODE smoother (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param filter A compatible filter path.
#' @param validate Validate the bundle and externally supplied filter.
#' @returns A `da_path_gaussian` smoother.
#' @noRd
.cgns_smoother_matrix_compiled <- function(bundle, filter, validate = TRUE) {
  if (isTRUE(validate)) {
    .validate_compiled_cgns(
      bundle, nontarget = bundle$nontarget, scalar = FALSE
    )
    .validate_gaussian_path(
      filter, bundle$obs, bundle$l, "filter", bundle$nontarget,
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

  for (j in (N1 - 1L):1L) {
    co <- .compiled_co(bundle, j)
    Gi <- .compiled_ginv(bundle, j)
    Cgi <- co$gyx %*% Gi
    A0 <- co$Ly - Cgi %*% co$Lx
    B <- sym(co$gyy - Cgi %*% t(co$gyx))
    Rfi <- chol_solve(filter$cov[, , j], Il, "Rf")
    H <- A0 + B %*% Rfi
    rate <- as.numeric(bundle$rate[j, ])
    for (ss in seq_len(nsub)) {
      d <- drop(Rfi %*% (filter$mean[j, ] - mus))
      mus <- mus + h * (
        -co$fy - drop(co$Ly %*% mus) + drop(B %*% d) +
          drop(Cgi %*% (-rate + co$fx + drop(co$Lx %*% mus)))
      )
      Rs <- Rs + h * (-(H %*% Rs) - Rs %*% t(H) + B)
      Rs <- spd_floor(sym(Rs))
    }
    MU[j, ] <- mus
    CV[, , j] <- Rs
  }

  p <- new_da_path(bundle$t, MU, CV, "smoother")
  p$meta$route <- if (bundle$correlated_noise)
    "backward_ode_correlated" else "backward_ode"
  p$meta$nsub <- nsub
  p$meta$init <- filter$meta$init
  p$meta$obs_x <- bundle$x
  p$meta$model <- bundle$model
  p$meta$nontarget <- bundle$nontarget
  p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
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
#' @returns A `da_path_gaussian` smoother.
#' @noRd
.cgns_smoother_compiled <- function(bundle, filter, validate = TRUE) {
  if (bundle$k == 1L && bundle$l == 1L)
    return(.cgns_smoother_scalar(bundle, filter, validate = validate))
  .cgns_smoother_matrix_compiled(bundle, filter, validate = validate)
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
      bundle, nontarget = bundle$nontarget, scalar = FALSE
    )
  }
  gaussian_kl_path(p, q, decompose = decompose)
}


#' Complete compiled closed-form ACI execution (internal)
#'
#' @param bundle A compiled CGNS bundle.
#' @param init Optional Gaussian prior.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of substeps.
#' @param keep Retain paths or discard them after metric construction.
#' @param decompose Return signal and dispersion components.
#' @returns An `aci_result` without a lag table.
#' @noRd
.aci_cgns_compiled <- function(
    bundle, init = NULL, stepper = c("explicit", "implicit"), nsub = 1L,
    keep = c("paths", "none"), decompose = TRUE) {
  stepper <- match.arg(stepper)
  keep <- match.arg(keep)
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  .validate_compiled_cgns(
    bundle, nontarget = bundle$nontarget, scalar = FALSE
  )
  filt <- .cgns_filter_compiled(
    bundle, init = init, stepper = stepper, nsub = nsub, validate = FALSE
  )
  smoo <- .cgns_smoother_compiled(bundle, filt, validate = FALSE)
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
        nontarget = bundle$nontarget,
        init = filt$meta$init
      ),
      meta = list(
        engine = "cgns",
        nontarget = bundle$nontarget,
        m = NULL,
        smoother_scheme = smoo$meta$route,
        table_reference = NULL
      )
    ),
    class = "aci_result"
  )
}
