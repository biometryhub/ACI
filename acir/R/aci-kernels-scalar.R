################################################################################
## aci-kernels-scalar.R - compiled-CGNS contract and scalar proof kernels
################################################################################


#' Refuse a recursion whose result is not finite (internal)
#'
#' The covariance policy owns positivity and finiteness of the second moment.
#' This is the corresponding invariant on the first moment and on the
#' likelihood accumulator, which no policy covers: a mean or a log-likelihood
#' that has overflowed to `Inf` or `NaN` is invalid arithmetic, not a
#' covariance that left the positive-definite cone, so the class is
#' deliberately distinct from `aci_error_covariance_not_spd`. Used by both the
#' scalar and the matrix kernels.
#'
#' @param quantity Name of the quantity, as it appears in the message.
#' @param index One-based grid index at which it is first not finite.
#' @param time Observation time at that index.
#' @param value The offending value.
#' @returns Never returns; raises `aci_error_nonfinite`.
#' @noRd
.aci_stop_nonfinite <- function(quantity, index, time, value) {
  aci_abort("aci_error_nonfinite", sprintf(paste(
    "The %s is not finite: it reached %s at index %d (time %g). This is a",
    "floating-point overflow in the recursion, not a covariance that left the",
    "positive-definite cone, and no regularisation policy recovers it; reduce",
    "the drift magnitude or the horizon."),
    quantity, format(value), index, time),
    quantity = quantity, index = index, time = time, value = value)
}


#' Construct the private compiled-CGNS representation (internal)
#'
#' Both generic and directed realisers terminate at this constructor so kernels
#' see one neutral contract rather than model-specific inputs.
#'
#' @param rs Output from `.resolve_nontarget()`.
#' @param source_model Original model supplied before non-target resolution.
#' @param source_obs Original observation path supplied before resolution.
#' @param coefficients Realised coefficient arrays and interval Gram weights.
#' @param correlated_noise Logical realised-path route flag.
#' @param realization Name of the private realisation route.
#' @returns A private `compiled_cgns` object.
#' @noRd
.new_compiled_cgns <- function(rs, source_model, source_obs, coefficients,
                               correlated_noise, realization) {
  m <- rs$model; ob <- rs$obs
  N1 <- length(ob$t); N <- N1 - 1L; k <- m$k; l <- m$l
  rate <- (ob$x[-1L, , drop = FALSE] -
           ob$x[-N1, , drop = FALSE]) / ob$dt
  out <- list(
    contract_version = 1L,
    source_model = source_model,
    source_obs = source_obs,
    model = m,
    obs = ob,
    conditional = rs$tag,
    likelihood_idx = as.integer(rs$likelihood_idx),
    k = k,
    l = l,
    N1 = N1,
    N = N,
    dt = ob$dt,
    t = ob$t,
    x = ob$x,
    rate = rate,
    coefficients = coefficients,
    correlated_noise = isTRUE(correlated_noise),
    realization = realization,
    provenance = list(source_model = source_model, source_obs = source_obs,
                      model = m, obs = ob, conditional = rs$tag)
  )
  if (k == 1L && l == 1L) {
    out$scalar <- list(
      Lx = as.numeric(coefficients$Lx),
      fx = as.numeric(coefficients$fx),
      Ly = as.numeric(coefficients$Ly),
      fy = as.numeric(coefficients$fy),
      gxx = as.numeric(coefficients$gxx),
      gyy = as.numeric(coefficients$gyy),
      gyx = as.numeric(coefficients$gyx),
      gxx_weight = as.numeric(coefficients$gxx_weight),
      rate = as.numeric(rate)
    )
  }
  structure(out, class = c("compiled_cgns", "aci_compiled"))
}


#' Compile a CGNS run using an authenticated batch realiser when available
#' (internal)
#'
#' This is the execution front door. Unauthenticated or modified models use
#' the deterministic one-pass generic compiler. It is also where a model that
#' declares its own observation set has that specification applied, so every
#' public entry point assimilates the estimand the constructor named.
#'
#' @param model A `cgns_model`.
#' @param obs Anything accepted by [as_obs()].
#' @param conditional Optional non-target specification.
#' @returns A compiled CGNS bundle.
#' @noRd
.compile_cgns_run <- function(model, obs, conditional = NULL) {
  conditional <- .model_estimand_spec(model, conditional)
  .check_prescribed_grid(model, obs)
  descriptor <- if (inherits(model, "cgns_model"))
    .cgns_realizer_descriptor(model) else NULL
  if (is.null(conditional) && !is.null(descriptor) &&
      identical(descriptor$id, "dyad_observed_x_v1"))
    return(.compile_dyad_cgns(model, obs, conditional = NULL))
  if (!is.null(descriptor) &&
      (identical(descriptor$id, "affine_model_v1") ||
       identical(descriptor$id, "enso6_aci_code_v1")))
    return(.compile_affine_cgns(model, obs, conditional = conditional))
  .compile_cgns_complete(model, obs, conditional = conditional)
}


#' Directed realiser for a built-in observed-x dyad (internal)
#'
#' A constructor-created model carries a locked descriptor containing pure
#' constructor data and exact references to all canonical coefficient
#' functions.  If that descriptor cannot be authenticated, generic compilation
#' is used instead; mutable names and metadata are never routing keys.
#'
#' @param model An unmodified `aci_dyad_model(observe = "x")` model.
#' @param obs Anything accepted by `as_obs()`.
#' @param conditional Must be `NULL` for this one-observation directed realiser.
#' @returns A private `compiled_cgns` object with the same contract as the generic
#'   closure fallback.
#' @noRd
.compile_dyad_cgns <- function(model, obs, conditional = NULL) {
  descriptor <- .cgns_realizer_descriptor(model)
  if (!inherits(model, "cgns_model") || model$k != 1L || model$l != 1L ||
      !is.null(conditional) || is.null(descriptor) ||
      !identical(descriptor$id, "dyad_observed_x_v1") ||
      !is.list(descriptor$spec$params))
    aci_abort("aci_error_compiled_contract",
              paste("The directed dyad realiser requires an unconditioned",
                    "observed-x aci_dyad_model()."))
  p <- descriptor$spec$params
  if (any(!is.finite(unlist(p, use.names = FALSE))))
    aci_abort("aci_error_compiled_contract",
              "The directed dyad realiser requires finite scalar parameters.")

  source_obs <- as_obs(obs)
  if (source_obs$k != 1L)
    aci_abort("aci_error_dims",
              "Observation dimension does not match the scalar dyad model.")
  rs <- .resolve_nontarget(model, source_obs, NULL)
  N1 <- length(source_obs$t); N <- N1 - 1L; x <- source_obs$x[, 1L]
  gxx_value <- p$s_x * p$s_x
  gyy_value <- p$s_y * p$s_y
  ## This realiser is the one route that never reaches
  ## `.compiled_precision_path()`, so it carries the same model-contract test
  ## on its own Gram.  The Gram is constant over the record, so one test
  ## covers every interval start and the index it names is the first.
  gxx_matrix <- matrix(gxx_value, 1L, 1L)
  if (all(is.finite(gxx_matrix))) {
    gxx_rcond <- rcond(gxx_matrix)
    if (gxx_rcond < .ACI_GRAM_RCOND_MIN)
      .aci_stop_gram_path(1L, source_obs$t[1L], gxx_rcond)
  }
  Gi <- as.matrix(rs$ginv(gxx_matrix))
  if (!is.numeric(Gi) || !identical(dim(Gi), c(1L, 1L)) ||
      any(!is.finite(Gi)))
    aci_abort("aci_error_gram",
              "Resolved observation-Gram operator is invalid for the directed dyad.")
  coefficients <- list(
    Lx = array(p$gamma * x, c(1L, 1L, N1)),
    fx = matrix(-p$d_x * x + p$f_x, N1, 1L),
    Ly = array(rep(-p$d_y, N1), c(1L, 1L, N1)),
    fy = matrix(-p$gamma * x^2 + p$f_y, N1, 1L),
    gxx = array(rep(gxx_value, N1), c(1L, 1L, N1)),
    gyy = array(rep(gyy_value, N1), c(1L, 1L, N1)),
    gyx = array(0, c(1L, 1L, N1)),
    gxx_weight = array(rep(Gi[1L, 1L], N), c(1L, 1L, N)))
  .new_compiled_cgns(
    rs, model, source_obs, coefficients,
    correlated_noise = isTRUE(model$meta$correlated_noise),
    realization = "dyad_directed")
}


#' Validate compiled-CGNS identity and structure (internal)
#'
#' @param bundle A private `compiled_cgns` object.
#' @param model Optional original model expected by the bundle.
#' @param obs Optional original observation path expected by the bundle.
#' @param conditional Optional expected non-target specification.
#' @param scalar Require a resolved one-observed/one-hidden system.
#' @returns `bundle`, invisibly.
#' @noRd
.validate_compiled_cgns <- function(bundle, model = NULL, obs = NULL,
                                    conditional = NULL, scalar = FALSE) {
  if (!inherits(bundle, "compiled_cgns") || !is.list(bundle))
    aci_abort("aci_error_compiled_contract",
              "Invalid or incomplete compiled-CGNS object.")
  co <- bundle$coefficients
  if (!identical(bundle$contract_version, 1L) ||
      !inherits(bundle$model, "cgns_model") ||
      !inherits(bundle$source_model, "cgns_model") ||
      !inherits(bundle$obs, "obs_traj") ||
      !inherits(bundle$source_obs, "obs_traj") ||
      !is.list(bundle$provenance) ||
      !identical(bundle$source_model, bundle$provenance$source_model) ||
      !identical(bundle$source_obs, bundle$provenance$source_obs) ||
      !identical(bundle$model, bundle$provenance$model) ||
      !identical(bundle$obs, bundle$provenance$obs) ||
      !identical(bundle$conditional %||% NULL,
                 bundle$provenance$conditional %||% NULL) ||
      length(bundle$N1) != 1L || bundle$N1 != length(bundle$t) ||
      bundle$N != bundle$N1 - 1L || bundle$N < 1L ||
      bundle$k != bundle$model$k || bundle$l != bundle$model$l ||
      bundle$source_obs$k != bundle$source_model$k ||
      !identical(dim(bundle$x), c(bundle$N1, bundle$k)) ||
      !identical(dim(bundle$rate), c(bundle$N, bundle$k)) ||
      !identical(bundle$t, bundle$obs$t) ||
      !identical(bundle$x, bundle$obs$x) ||
      !identical(bundle$dt, bundle$obs$dt) ||
      !is.list(co) ||
      !identical(dim(co$Lx), c(bundle$k, bundle$l, bundle$N1)) ||
      !identical(dim(co$fx), c(bundle$N1, bundle$k)) ||
      !identical(dim(co$Ly), c(bundle$l, bundle$l, bundle$N1)) ||
      !identical(dim(co$fy), c(bundle$N1, bundle$l)) ||
      !identical(dim(co$gxx), c(bundle$k, bundle$k, bundle$N1)) ||
      !identical(dim(co$gyy), c(bundle$l, bundle$l, bundle$N1)) ||
      !identical(dim(co$gyx), c(bundle$l, bundle$k, bundle$N1)) ||
      !identical(dim(co$gxx_weight), c(bundle$k, bundle$k, bundle$N)) ||
      any(!vapply(co, function(z) all(is.finite(z)), logical(1))) ||
      length(bundle$correlated_noise) != 1L ||
      !is.logical(bundle$correlated_noise) || is.na(bundle$correlated_noise) ||
      any(!is.finite(c(bundle$t, bundle$x, bundle$rate, bundle$dt))))
    aci_abort("aci_error_compiled_contract",
              "Invalid or incomplete compiled-CGNS object.")
  if (!is.null(model) && !identical(bundle$source_model, model))
    aci_abort("aci_error_model_contract",
              "The compiled coefficients belong to a different model object.")
  if (!is.null(obs)) {
    ob <- as_obs(obs)
    if (!identical(dim(bundle$source_obs$x), dim(ob$x)) ||
        length(bundle$source_obs$t) != length(ob$t) ||
        max(abs(bundle$source_obs$t - ob$t)) >
          1e-10 * max(1, max(abs(ob$t))) ||
        any(abs(bundle$source_obs$x - ob$x) >
              1e-12 * pmax(1, abs(ob$x))))
      aci_abort("aci_error_obs_contract",
                "The compiled coefficients belong to different observations.")
  }
  if (!identical(bundle$conditional %||% NULL, conditional %||% NULL))
    aci_abort("aci_error_nontarget",
              "The compiled coefficients use a different non-target specification.")
  if (isTRUE(scalar)) {
    s <- bundle$scalar
    fields <- c("Lx", "fx", "Ly", "fy", "gxx", "gyy", "gyx")
    if (bundle$k != 1L || bundle$l != 1L || !is.list(s) ||
        any(vapply(s[fields], length, integer(1)) != bundle$N1) ||
        length(s$gxx_weight) != bundle$N || length(s$rate) != bundle$N ||
        any(!is.finite(unlist(s, use.names = FALSE))) ||
        !identical(s$Lx, as.numeric(co$Lx)) ||
        !identical(s$fx, as.numeric(co$fx)) ||
        !identical(s$Ly, as.numeric(co$Ly)) ||
        !identical(s$fy, as.numeric(co$fy)) ||
        !identical(s$gxx, as.numeric(co$gxx)) ||
        !identical(s$gyy, as.numeric(co$gyy)) ||
        !identical(s$gyx, as.numeric(co$gyx)) ||
        !identical(s$gxx_weight, as.numeric(co$gxx_weight)) ||
        !identical(s$rate, as.numeric(bundle$rate)) ||
        !identical(bundle$likelihood_idx, 1L))
      aci_abort("aci_error_compiled_contract",
                "The compiled-CGNS object is not eligible for the scalar kernel.")
  }
  invisible(bundle)
}


#' Scalar equivalent of `spd_floor()` (internal)
#'
#' @param value One covariance value.
#' @returns A positive scalar under the current eigenvalue-floor policy.
#' @noRd
.scalar_spd_floor <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
    aci_abort("aci_error_spd",
              "spd_floor() needs a finite non-empty square matrix and eps > 0.")
  if (value > 0) return(value)
  max(1e-12 * max(abs(value), 1e-300), 1e-300)
}


#' Prepare the initial state for a compiled scalar filter (internal)
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param init Optional Gaussian prior.
#' @returns Initial values and provenance.
#' @noRd
.scalar_filter_init <- function(bundle, init) {
  .check_init_names(init)
  s <- bundle$scalar
  co0 <- list(gyy = matrix(s$gyy[1L], 1L, 1L),
              Ly = matrix(s$Ly[1L], 1L, 1L))
  used_default_cov <- is.null(init) || is.null(init$cov)
  if (used_default_cov) {
    di <- .default_init(bundle$model, co0)
    if (is.null(init)) init <- di else init$cov <- di$cov
    aci_warn("aci_warn_diffuse_init",
      paste("No init$cov supplied; using a diffuse prior. Its opening steps",
            "are prior-dominated; a prior far wider than the hidden state's",
            "own scale can also destabilise the explicit step, which is a",
            "refusal rather than a window to discard."))
  }
  mu <- as.numeric(init$mean %||% 0)
  R <- as.matrix(init$cov)
  if (length(mu) != 1L || !identical(dim(R), c(1L, 1L)) ||
      !is.numeric(R) || any(!is.finite(c(mu, R))))
    aci_abort("aci_error_dims",
              "init$mean/init$cov have incompatible dimensions or values.")
  if (max(abs(R - t(R))) > 1e-12 * max(1, max(abs(R))))
    aci_abort("aci_error_spd", "init$cov must be symmetric.")
  R <- sym(R)
  if (is.null(tryCatch(chol(R), error = function(e) NULL)))
    aci_abort("aci_error_spd",
              "init$cov must be positive definite; it is not regularized implicitly.")
  list(mu = unname(mu), R = unname(R[1L, 1L]),
       actual = list(mean = unname(mu), cov = R,
                     source = if (used_default_cov)
                       "default_diffuse" else "user"))
}


#' Explicit scalar filter state recursion (internal)
#'
#' This kernel assumes its bundle and initial values have already been checked.
#' Predictive likelihood is deliberately separate so state-only work can be
#' benchmarked against a numerically matched comparison kernel.
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param mu Initial hidden mean.
#' @param R Initial hidden variance.
#' @param nsub Positive number of explicit substeps.
#' @param rec Covariance-policy recorder from `.aci_reg_new()`.
#' @returns Numeric filter moments plus the explicit-stability diagnostic.
#' @noRd
.cgns_filter_scalar_kernel <- function(bundle, mu, R, nsub = 1L,
                                       rec = .aci_reg_for(NULL, bundle$t)) {
  s <- bundle$scalar; N1 <- bundle$N1; N <- bundle$N
  h <- bundle$dt / nsub
  Lx <- s$Lx; fx <- s$fx; Ly <- s$Ly; fy <- s$fy
  gyy <- s$gyy; gyx <- s$gyx; Gi <- s$gxx_weight; rate <- s$rate
  MU <- numeric(N1); CV <- numeric(N1)
  MU[1L] <- mu; CV[1L] <- R; stab <- 0
  if (nsub == 1L) {
    for (j in seq_len(N)) {
      Kg <- (R * Lx[j] + gyx[j]) * Gi[j]
      now_stab <- abs(Lx[j] * Gi[j] * Lx[j] * R) * h
      if (now_stab > stab) stab <- now_stab
      mu <- mu + (fy[j] + Ly[j] * mu) * h +
        Kg * ((rate[j] - (fx[j] + Lx[j] * mu)) * h)
      R <- R + (Ly[j] * R + R * Ly[j] + gyy[j] -
        Kg * (Lx[j] * R + gyx[j])) * h
      ## Two tests, in this order, as before: the first is reached only by a
      ## divergent variance and aborts inside the guard, the second is the
      ## policy. Keeping them apart rather than joining them with `||` is
      ## worth about 1% of this kernel.
      if (!is.finite(R)) { rec$j <- j + 1L
        .cov_guard_scalar(R, rec, "filter_explicit") }
      if (R <= 0) { rec$j <- j + 1L
        R <- .cov_guard_scalar(R, rec, "filter_explicit") }
      MU[j + 1L] <- mu; CV[j + 1L] <- R
    }
  } else {
    for (j in seq_len(N)) {
      now_stab <- abs(Lx[j] * Gi[j] * Lx[j] * R) * h
      if (now_stab > stab) stab <- now_stab
      for (ss in seq_len(nsub)) {
        Kg <- (R * Lx[j] + gyx[j]) * Gi[j]
        mu <- mu + (fy[j] + Ly[j] * mu) * h +
          Kg * ((rate[j] - (fx[j] + Lx[j] * mu)) * h)
        R <- R + (Ly[j] * R + R * Ly[j] + gyy[j] -
          Kg * (Lx[j] * R + gyx[j])) * h
        if (!is.finite(R)) { rec$j <- j + 1L
          .cov_guard_scalar(R, rec, "filter_explicit") }
        if (R <= 0) { rec$j <- j + 1L
          R <- .cov_guard_scalar(R, rec, "filter_explicit") }
      }
      MU[j + 1L] <- mu; CV[j + 1L] <- R
    }
  }
  ## One vector test after the recursion rather than one scalar test per step:
  ## this is the shortest kernel in the package, and the whole-path form costs
  ## about a fortieth of the in-loop form while naming the same first index.
  if (!all(is.finite(MU))) {
    bad <- which(!is.finite(MU))[1L]
    .aci_stop_nonfinite("filter mean", bad, bundle$t[bad], MU[bad])
  }
  list(mean = MU, cov = CV, stability = stab)
}


#' Scalar predictive log-likelihood (internal)
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param filter_mean Filter mean at every grid point.
#' @param filter_cov Filter variance at every grid point.
#' @param rec Covariance-policy recorder from `.aci_reg_new()`.
#' @returns Accumulated predictive log-likelihood.
#' @noRd
.cgns_likelihood_scalar_kernel <- function(
    bundle, filter_mean, filter_cov, rec = .aci_reg_for(NULL, bundle$t)) {
  s <- bundle$scalar; dt <- bundle$dt; ll <- 0
  ## Hot-loop bindings, as in the filter. The four vectors are looked up once
  ## rather than once per grid point; the arithmetic and its order are
  ## untouched, so the accumulation is the same additions in the same sequence.
  rate <- s$rate; fx <- s$fx; Lx <- s$Lx; gxx <- s$gxx
  for (j in seq_len(bundle$N)) {
    iota0 <- (rate[j] - (fx[j] + Lx[j] * filter_mean[j])) * dt
    Sd <- (Lx[j]^2 * filter_cov[j] * dt + gxx[j]) * dt
    if (!is.finite(Sd)) { rec$j <- j
      .cov_guard_scalar(Sd, rec, "likelihood_innov") }
    if (Sd <= 0) { rec$j <- j
      Sd <- .cov_guard_scalar(Sd, rec, "likelihood_innov") }
    ch <- sqrt(Sd)
    w <- iota0 / ch
    ll <- ll - 0.5 * w * w - log(ch) - 0.5 * log(2 * pi)
  }
  ## Defence in depth: a finite mean and a positive finite innovation variance
  ## give a finite term, so with the filter guard in place this is a
  ## post-condition.  The index reported is the first scored interval whose
  ## input moments are not finite, and the last scored interval when the
  ## accumulation itself overflowed on finite inputs.
  if (!is.finite(ll)) {
    n <- bundle$N
    bad <- which(!is.finite(filter_mean[seq_len(n)]) |
                   !is.finite(filter_cov[seq_len(n)]))[1L]
    if (is.na(bad)) bad <- n
    .aci_stop_nonfinite("predictive log-likelihood", bad, bundle$t[bad], ll)
  }
  ll
}


#' Construct a filter path from the explicit scalar kernels (internal)
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param init Optional Gaussian prior.
#' @param nsub Positive number of explicit substeps.
#' @param validate Validate the compiled bundle before execution.
#' @param loglik Accumulate the predictive log-likelihood. `FALSE` skips the
#'   likelihood kernel entirely and leaves `meta$loglik` unset; the state
#'   recursion is untouched.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A contract-shaped `da_path_gaussian` filter.
#' @noRd
.cgns_filter_scalar <- function(bundle, init = NULL, nsub = 1L,
                                validate = TRUE, loglik = TRUE,
                                regularize = NULL) {
  rec <- .aci_reg_for(regularize, bundle$t)
  if (isTRUE(validate))
    .validate_compiled_cgns(bundle, conditional = bundle$conditional,
                            scalar = TRUE)
  ## The range test precedes the integrality test, and the integrality test is
  ## `%% 1` rather than `as.integer()`: coercing above the integer range
  ## returns NA with a base warning, which is the coercion this check exists
  ## to protect.  `%%` also leaves nsub = TRUE accepted, which trunc() and
  ## floor() would not.
  if (length(nsub) != 1L || !is.finite(nsub) || nsub < 1 ||
      nsub > .Machine$integer.max || nsub %% 1 != 0)
    aci_abort("aci_error_dims", "nsub must be a positive integer.")
  nsub <- as.integer(nsub)
  ini <- .scalar_filter_init(bundle, init)
  if (bundle$dt * abs(bundle$scalar$Ly[1L]) > 0.1 && nsub == 1L)
    aci_warn("aci_warn_dt_stability",
             "dt * ||Ly|| > 0.1 at t0; Euler updates may be inaccurate (consider nsub > 1).")
  out <- .cgns_filter_scalar_kernel(bundle, ini$mu, ini$R, nsub, rec)
  do_ll <- isTRUE(loglik)
  ll <- if (do_ll)
    .cgns_likelihood_scalar_kernel(bundle, out$mean, out$cov, rec) else NULL
  if (out$stability > 1)
    aci_warn("aci_warn_riccati_stiff", sprintf(paste(
      "Explicit Riccati step is unstable (max ||Lx' gxx^-1 Lx R|| dt = %.3g > 1):",
      "the covariance can overshoot, leave the positive-definite cone, and oscillate.",
      "Use the positivity-preserving implicit stepper, or reduce dt / increase nsub."),
      out$stability))

  MU <- matrix(out$mean, bundle$N1, 1L)
  CV <- array(out$cov, c(1L, 1L, bundle$N1))
  p <- new_da_path(bundle$t, MU, CV, "filter")
  p$meta$stepper <- "explicit"; p$meta$nsub <- nsub
  ## loglik = FALSE leaves meta$loglik unset (NULL); likelihood_idx still
  ## records which observation rows the likelihood contract would score.
  if (do_ll) p$meta$loglik <- ll
  p$meta$likelihood_idx <- bundle$likelihood_idx; p$meta$init <- ini$actual
  p$meta$obs_x <- bundle$x; p$meta$model <- bundle$model
  p$meta$conditional <- bundle$conditional; p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p$meta$regularization <- .aci_reg_report(rec)
  p
}


#' Scalar backward-ODE smoother recursion (internal)
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param filter_mean Filter means.
#' @param filter_cov Filter variances.
#' @param nsub Positive number of substeps inherited from the filter.
#' @param rec Covariance-policy recorder from `.aci_reg_new()`.
#' @returns Numeric smoother moments.
#' @noRd
.cgns_smoother_scalar_kernel <- function(bundle, filter_mean, filter_cov,
                                         nsub = 1L,
                                         rec = .aci_reg_for(NULL, bundle$t)) {
  s <- bundle$scalar; N1 <- bundle$N1; h <- bundle$dt / nsub
  Lx <- s$Lx; fx <- s$fx; Ly <- s$Ly; fy <- s$fy
  gyy <- s$gyy; gyx <- s$gyx; Gi <- s$gxx_weight; rate <- s$rate
  MU <- numeric(N1); CV <- numeric(N1)
  MU[N1] <- filter_mean[N1]; CV[N1] <- filter_cov[N1]
  mus <- MU[N1]; Rs <- CV[N1]
  Cgi <- gyx[seq_len(bundle$N)] * Gi
  B <- gyy[seq_len(bundle$N)] - Cgi * gyx[seq_len(bundle$N)]
  Rfi <- 1 / filter_cov[seq_len(bundle$N)]
  H <- Ly[seq_len(bundle$N)] - Cgi * Lx[seq_len(bundle$N)] + B * Rfi
  if (nsub == 1L) {
    for (j in (N1 - 1L):1L) {
      d <- Rfi[j] * (filter_mean[j] - mus)
      mus <- mus + h * (-fy[j] - Ly[j] * mus + B[j] * d +
        Cgi[j] * (-rate[j] + fx[j] + Lx[j] * mus))
      Rs <- Rs + h * (-H[j] * Rs - Rs * H[j] + B[j])
      if (!is.finite(Rs)) { rec$j <- j
        .cov_guard_scalar(Rs, rec, "smoother_backward") }
      if (Rs <= 0) { rec$j <- j
        Rs <- .cov_guard_scalar(Rs, rec, "smoother_backward") }
      MU[j] <- mus; CV[j] <- Rs
    }
  } else {
    for (j in (N1 - 1L):1L) {
      for (ss in seq_len(nsub)) {
        d <- Rfi[j] * (filter_mean[j] - mus)
        mus <- mus + h * (-fy[j] - Ly[j] * mus + B[j] * d +
          Cgi[j] * (-rate[j] + fx[j] + Lx[j] * mus))
        Rs <- Rs + h * (-H[j] * Rs - Rs * H[j] + B[j])
        if (!is.finite(Rs)) { rec$j <- j
          .cov_guard_scalar(Rs, rec, "smoother_backward") }
        if (Rs <= 0) { rec$j <- j
          Rs <- .cov_guard_scalar(Rs, rec, "smoother_backward") }
      }
      MU[j] <- mus; CV[j] <- Rs
    }
  }
  if (!all(is.finite(MU))) {
    bad <- which(!is.finite(MU))[1L]
    .aci_stop_nonfinite("smoother mean", bad, bundle$t[bad], MU[bad])
  }
  list(mean = MU, cov = CV)
}


#' Construct a smoother path from the scalar kernel (internal)
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param filter A compatible filter path.
#' @param validate Validate a user-supplied path; trusted same-run paths may skip.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns A contract-shaped `da_path_gaussian` smoother.
#' @noRd
.cgns_smoother_scalar <- function(bundle, filter, validate = TRUE,
                                  regularize = NULL) {
  rec <- .aci_reg_for(regularize, bundle$t)
  if (isTRUE(validate)) {
    .validate_compiled_cgns(bundle, conditional = bundle$conditional,
                            scalar = TRUE)
    .validate_gaussian_path(filter, bundle$obs, 1L, "filter",
                            bundle$conditional, model = bundle$model,
                            source_model = bundle$source_model)
  }
  nsub <- max(1L, as.integer(filter$meta$nsub %||% 1L))
  fm <- as.numeric(filter$mean[, 1L]); fc <- as.numeric(filter$cov[1L, 1L, ])
  out <- .cgns_smoother_scalar_kernel(bundle, fm, fc, nsub, rec)
  MU <- matrix(out$mean, bundle$N1, 1L)
  CV <- array(out$cov, c(1L, 1L, bundle$N1))
  p <- new_da_path(bundle$t, MU, CV, "smoother")
  p$meta$route <- if (bundle$correlated_noise)
    "backward_ode_correlated" else "backward_ode"
  p$meta$scheme <- "backward_ode_euler"
  p$meta$nsub <- nsub; p$meta$init <- filter$meta$init
  p$meta$obs_x <- bundle$x; p$meta$model <- bundle$model
  p$meta$conditional <- bundle$conditional; p$meta$engine <- "cgns"
  p$meta$source_model <- bundle$source_model
  p$meta$regularization <- .aci_reg_report(rec)
  stopifnot(max(abs(p$mean[bundle$N1, ] -
                        filter$mean[bundle$N1, ])) < 1e-12)
  p
}


#' Error-free sum of two doubles (internal)
#'
#' Knuth's TwoSum (Ogita, Rump & Oishi, 2005, Algorithm 3.1,
#' doi:10.1137/030601818): `s` is the rounded sum and `e` its rounding error,
#' so that `s + e` equals `a + b` exactly.
#'
#' @param a,b Numeric vectors of equal length whose sums stay finite.
#' @returns A list with elements `s` and `e`.
#' @noRd
.two_sum <- function(a, b) {
  s <- a + b
  t <- s - a
  list(s = s, e = (a - (s - t)) + (b - t))
}


#' Error-free product of two doubles (internal)
#'
#' Dekker's TwoProduct with Veltkamp's splitting (Ogita, Rump & Oishi, 2005,
#' Algorithms 3.2 and 3.3): `p` is the rounded product and `e` its rounding
#' error, so that `p + e` equals `a * b` exactly unless `e` underflows.
#'
#' @param a,b Numeric vectors of equal length, small enough that their
#'   products, and their products with `2^27 + 1`, stay finite.
#' @returns A list with elements `p` and `e`.
#' @noRd
.two_product <- function(a, b) {
  p <- a * b
  split_a <- 134217729 * a # 2^27 + 1, the binary64 splitting factor.
  a_hi <- split_a - (split_a - a)
  a_lo <- a - a_hi
  split_b <- 134217729 * b
  b_hi <- split_b - (split_b - b)
  b_lo <- b - b_hi
  list(p = p,
       e = a_lo * b_lo - (((p - a_hi * b_hi) - a_lo * b_hi) - a_hi * b_lo))
}


#' Sign of the exact sum of doubles (internal)
#'
#' Adds the terms one at a time by Shewchuk's Grow-Expansion, which keeps
#' each rounding error as a further component, so the components always sum
#' exactly to the terms added so far. The components do not overlap and
#' increase in magnitude, so the largest nonzero component has the sign of
#' the sum (Shewchuk, 1997, doi:10.1007/PL00009321).
#'
#' @param terms List of numeric vectors of equal length whose partial sums
#'   stay finite.
#' @returns Numeric vector of -1, 0 and 1, the sign of each exact sum.
#' @noRd
.exact_sum_sign <- function(terms) {
  parts <- list()
  for (term in terms) {
    carry <- term
    for (k in seq_along(parts)) {
      step <- .two_sum(carry, parts[[k]])
      parts[[k]] <- step$e
      carry <- step$s
    }
    parts[[length(parts) + 1L]] <- carry
  }
  out <- numeric(length(terms[[1L]]))
  for (part in parts) out[part != 0] <- sign(part[part != 0])
  out
}


#' Scalar Gaussian KL path arithmetic (internal)
#'
#' @param mu_p,mu_q Numeric mean paths.
#' @param R_p,R_q Numeric variance paths.
#' @param decompose Return signal and dispersion vectors as well as total.
#' @returns A list of numeric vectors.
#' @noRd
.gaussian_kl_scalar_kernel <- function(mu_p, R_p, mu_q, R_q,
                                       decompose = TRUE) {
  if (any(R_q <= 0))
    aci_abort("aci_error_spd", "Matrix (R_q) must be positive definite.")
  if (any(R_p <= 0))
    aci_abort("aci_error_spd", "Matrix (R_p) must be positive definite.")
  n <- length(mu_p)
  if (length(R_p) != n || length(mu_q) != n || length(R_q) != n) {
    ## The branches below index single elements, so recycle the arguments to
    ## the longest length first, or to length zero if any is empty.
    n <- max(n, length(R_p), length(mu_q), length(R_q))
    if (min(length(mu_p), length(R_p), length(mu_q), length(R_q)) == 0L)
      n <- 0L
    if (length(mu_p) != n) mu_p <- rep_len(mu_p, n)
    if (length(R_p) != n) R_p <- rep_len(R_p, n)
    if (length(mu_q) != n) mu_q <- rep_len(mu_q, n)
    if (length(R_q) != n) R_q <- rep_len(R_q, n)
  }
  Lq <- sqrt(R_q)
  signal <- 0.5 * ((mu_q - mu_p) / Lq)^2
  if (is.infinite(max(signal, 0, na.rm = TRUE))) {
    j <- which(is.infinite(signal))
    ## A finite standardised difference above 2^514 gives a signal above
    ## 2^1027, which no rounding brings back into the double range, so only
    ## the other overflows need the work below.
    z <- (mu_q[j] - mu_p[j]) / Lq[j]
    j <- j[!(is.finite(z) & abs(z) > 2^514)]
    half_difference <- 0.5 * (mu_q[j] - mu_p[j])
    large_difference <- is.infinite(half_difference)
    half_difference[large_difference] <-
      0.5 * mu_q[j[large_difference]] - 0.5 * mu_p[j[large_difference]]
    # Halve before squaring. Divide first only when that square would overflow.
    square_first <- abs(half_difference) <= sqrt(.Machine$double.xmax)
    h <- half_difference[square_first]
    signal[j[square_first]] <- 2 * ((h * h) / R_q[j[square_first]])
    h <- half_difference[!square_first]
    signal[j[!square_first]] <- 2 * ((h / R_q[j[!square_first]]) * h)

    # Rounding near the finite limit can still cross the overflow boundary,
    # so compare the scaled squared difference with that boundary exactly,
    # with the error-free sums and products above. Ordinary finite signals
    # never enter this branch.
    near <- !is.finite(signal[j]) |
      signal[j] > .Machine$double.xmax * (1 - 4 * .Machine$double.eps)
    if (any(near)) {
      j <- j[near]
      scale <- 2^floor(floor(log2(R_q[j])) / 2)
      a <- (0.5 * mu_q[j] / scale) / 2^512
      b <- (-0.5 * mu_p[j] / scale) / 2^512
      x <- a + b
      keep <- is.finite(x) & abs(x) >= 0.25 & abs(x) <= 2
      if (any(keep)) {
        j <- j[keep]
        scale <- scale[keep]
        a <- a[keep]
        b <- b[keep]
        x <- x[keep]
        x_error <- .two_sum(a, b)$e
        v <- (R_q[j] / scale) / scale
        # The signal overflows when 2 * (x + x_error)^2 / v reaches
        # 1 - 2^-54, the round-to-nearest overflow midpoint divided by 2^1024.
        # Expand the square into exact products and take the sign of the whole
        # difference 2 * (x + x_error)^2 - v * (1 - 2^-54).
        square <- .two_product(x, x)
        cross <- .two_product(4 * x, x_error)
        error_square <- .two_product(2 * x_error, x_error)
        gap <- .exact_sum_sign(list(
          -v, v * 2^-54, 2 * square$p, 2 * square$e, cross$p, cross$e,
          error_square$p, error_square$e
        ))
        # A half that underflowed to zero decides an exact tie by its sign.
        lost_a <- gap == 0 & a == 0 & mu_q[j] != 0
        lost_b <- gap == 0 & b == 0 & mu_p[j] != 0
        gap[lost_a] <- sign(x[lost_a]) * sign(mu_q[j[lost_a]])
        gap[lost_b] <- -sign(x[lost_b]) * sign(mu_p[j[lost_b]])
        overflow <- gap >= 0
        signal[j[overflow]] <- Inf
        recover <- j[!overflow & !is.finite(signal[j])]
        signal[recover] <- .Machine$double.xmax
      }
    }
  }
  ## Take the logarithm before subtracting one so small ratios survive.
  ## Separate variance logs retain precision if the ratio is subnormal
  ## or underflows, without slowing the ordinary path with a second log.
  ratio <- R_p / R_q
  ratio_min <- min(ratio, 1)
  ratio_max <- max(ratio, 0)
  if (is.na(ratio_min)) {
    # Two infinite variances give a NaN ratio, which stays NaN below.
    ratio_min <- min(ratio, 1, na.rm = TRUE)
    ratio_max <- max(ratio, 0, na.rm = TRUE)
  }
  if (ratio_min == 1 && ratio_max == 1) {
    # Every ratio is exactly 1 or NaN, so ratio - 1 is the dispersion.
    dispersion <- ratio - 1
  } else {
    log_ratio <- log(ratio)
    if (ratio_min < .Machine$double.xmin) {
      subnormal <- which(ratio < .Machine$double.xmin)
      log_ratio[subnormal] <- log(R_p[subnormal]) - log(R_q[subnormal])
    }
    dispersion <- 0.5 * ((ratio - 1) - log_ratio)
    if (is.infinite(ratio_max)) {
      overflow <- is.infinite(ratio)
      # The half-ratio can be finite even when the ratio overflows.
      dispersion[overflow] <- (0.5 * R_p[overflow]) / R_q[overflow] -
        0.5 * (1 + log(R_p[overflow]) - log(R_q[overflow]))
    }
  }
  signal <- pmax(signal, 0); dispersion <- pmax(dispersion, 0)
  total <- signal + dispersion
  if (isTRUE(getOption("aci.debug_assert", FALSE)) && any(total < -1e-10))
    aci_abort("aci_error_internal", "Negative KL encountered.")
  if (isTRUE(decompose))
    list(total = total, signal = signal, dispersion = dispersion) else
      list(total = total)
}


#' Construct a scalar Gaussian KL path (internal)
#'
#' @param bundle A scalar `compiled_cgns` object supplying the trusted grid.
#' @param p,q Compatible scalar Gaussian paths.
#' @param decompose Return the current three-component result when `TRUE`.
#' @returns A data frame shaped like `aci_metric()`.
#' @noRd
.gaussian_kl_path_scalar <- function(bundle, p, q, decompose = TRUE,
                                     validate = TRUE) {
  if (isTRUE(validate))
    .validate_compiled_cgns(bundle, conditional = bundle$conditional,
                            scalar = TRUE)
  if (!inherits(p, "da_path_gaussian") || !inherits(q, "da_path_gaussian"))
    aci_abort("aci_error_dims", "p and q must be Gaussian assimilation paths.")
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  n <- bundle$N1
  if (length(p$t) != n || length(q$t) != n ||
      max(abs(p$t - q$t)) > 1e-10 ||
      max(abs(p$t - bundle$t)) > 1e-10 * max(1, max(abs(bundle$t))) ||
      max(abs(q$t - bundle$t)) > 1e-10 * max(1, max(abs(bundle$t))) ||
      !identical(dim(p$mean), c(n, 1L)) ||
      !identical(dim(q$mean), c(n, 1L)) ||
      !identical(dim(p$cov), c(1L, 1L, n)) ||
      !identical(dim(q$cov), c(1L, 1L, n)) ||
      any(!is.finite(c(p$t, q$t, p$mean, q$mean, p$cov, q$cov))))
    aci_abort("aci_error_dims",
              "aci_metric: incompatible or non-finite moments.")
  out <- .gaussian_kl_scalar_kernel(
    as.numeric(p$mean), as.numeric(p$cov),
    as.numeric(q$mean), as.numeric(q$cov), decompose)
  data.frame(t = p$t, out, check.names = FALSE)
}


#' Complete private scalar ACI proof path (internal)
#'
#' @param bundle A scalar `compiled_cgns` object.
#' @param init Optional Gaussian prior.
#' @param nsub Positive number of explicit substeps.
#' @param keep Retain paths or discard them after metric construction.
#' @param decompose Return signal and dispersion components.
#' @param regularize Covariance policy, or a recorder shared with the caller.
#' @returns An `aci_result` produced by the explicit scalar kernels.
#' @noRd
.aci_scalar_compiled <- function(bundle, init = NULL, nsub = 1L,
                                 keep = c("paths", "none"),
                                 decompose = TRUE, regularize = NULL) {
  keep <- match.arg(keep)
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  rec <- .aci_reg_for(regularize, bundle$t)
  .validate_compiled_cgns(bundle, conditional = bundle$conditional,
                          scalar = TRUE)
  filt <- .cgns_filter_scalar(bundle, init = init, nsub = nsub,
                              validate = FALSE, regularize = rec)
  smoo <- .cgns_smoother_scalar(bundle, filt, validate = FALSE,
                                regularize = rec)
  klp <- .gaussian_kl_path_scalar(bundle, smoo, filt,
                                  decompose = decompose, validate = FALSE)
  structure(list(
    t = klp$t,
    aci = klp$total,
    signal = if (decompose) klp$signal else NULL,
    dispersion = if (decompose) klp$dispersion else NULL,
    paths = if (keep == "paths") list(filter = filt, smoother = smoo) else NULL,
    table = NULL,
    handles = list(model = bundle$source_model, obs = bundle$source_obs,
                   conditional = bundle$conditional, init = filt$meta$init),
    meta = list(engine = "cgns", conditional = bundle$conditional, m = NULL,
                smoother_scheme = smoo$meta$scheme,
                table_reference = NULL,
                regularization = .aci_reg_report(rec))),
    class = "aci_result")
}
