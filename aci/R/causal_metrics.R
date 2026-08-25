################################################################################
## causal_metrics.R - causal layer: ACI, CIR, KL metric
## ########################################################################## ##
##
## Contents:
##   * Gaussian relative entropy, signal/dispersion decomposition, paths:
##       - gaussian_kl, kl_increment, gaussian_kl_path, projected_kl, empirical_kl
##
##   * the ACI metric front-end aci():
##       - aci, print.aci_result, as.data.frame.aci_result
##
##   * forward/backward causal influence ranges (exact + l1_linf):
##       - .fwd_lengths, .bwd_lengths, new_cir_result, print.cir_result, forward_cir,
##         forward_cir.lag_table, forward_cir.aci_result, backward_cir, backward_cir.lag_table,
##         backward_cir.aci_result, cir_pair, backward_cir.cgns_model
##
##   * package-level tidy read-out of the core CIR estimands:
##       - cir_table
##
################################################################################


################################################################################
# Gaussian relative entropy, signal/dispersion decomposition, paths
################################################################################

#' Gaussian relative entropy
#'
#' Relative entropy of one multivariate Gaussian from another, optionally split
#' into its signal and dispersion parts.
#'
#' Public KL values never apply a covariance ridge; callers wanting
#' regularisation opt in explicitly with [spd_floor()]. andreou2026cir
#' (Section 2.2, closing paragraph) states that regularization is expected
#' only in the degenerate limit, which is the published basis for keeping the
#' public KL strict while the internal recursions floor.
#'
#' @param mu_p Numeric vector, mean of the first distribution.
#' @param R_p Covariance matrix of the first distribution.
#' @param mu_q Numeric vector, mean of the second distribution.
#' @param R_q Covariance matrix of the second distribution.
#' @param decompose `TRUE` to return the signal and dispersion parts alongside
#'   the total.
#' @returns A named numeric vector with the `total` and, when `decompose` is
#'   `TRUE`, the `signal` and `dispersion` parts.
#'
#' @seealso [gaussian_kl_path()], [aci()]
#'
#' @examples
#' gaussian_kl(mu_p = 0, R_p = matrix(1), mu_q = 1, R_q = matrix(2))
#'
#' @export
gaussian_kl <- function(mu_p, R_p, mu_q, R_q, decompose = TRUE) {
  mu_p <- as.numeric(mu_p); mu_q <- as.numeric(mu_q)
  R_p <- as.matrix(R_p); R_q <- as.matrix(R_q)
  l <- length(mu_p)
  if (!is.numeric(R_p) || !is.numeric(R_q) || l < 1L || length(mu_q) != l ||
      !identical(dim(R_p), c(l, l)) || !identical(dim(R_q), c(l, l)) ||
      any(!is.finite(c(mu_p, mu_q, R_p, R_q))))
    aci_abort("aci_error_dims", "gaussian_kl: dimension mismatch.")
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  # Public KL values must not depend on an undocumented covariance ridge.
  # Callers wanting regularisation can opt in explicitly with spd_floor().
  Lq <- .strict_chol(R_q, "R_q"); Lp <- .strict_chol(R_p, "R_p")
  d <- mu_q - mu_p
  w <- forwardsolve(t(Lq), d)                     # Lq' \ d
  signal <- 0.5 * sum(w * w)
  A <- forwardsolve(t(Lq), t(Lp))                 # tr(Rq^{-1}Rp) = ||Lq'^{-1} Lp'||_F^2
  tr <- sum(A * A)
  dispersion <- 0.5 * (tr - l + 2 * sum(log(diag(Lq))) - 2 * sum(log(diag(Lp))))
  signal <- max(signal, 0); dispersion <- max(dispersion, 0)
  tot <- signal + dispersion
  if (isTRUE(getOption("aci.debug_assert", FALSE)) && tot < -1e-10)
    aci_abort("aci_error_internal", "Negative KL encountered.")
  if (decompose) c(total = tot, signal = signal, dispersion = dispersion) else tot
}


#' Gaussian KL arithmetic for already dimension-validated values (internal)
#'
#' Covariances remain strict: the reference is checked before the integrating
#' covariance, matching [gaussian_kl()].
#'
#' @param mu_p,R_p Integrating Gaussian moments.
#' @param mu_q,R_q Reference Gaussian moments.
#' @param decompose Return the three-component result.
#' @returns A named numeric vector or scalar total.
#' @noRd
.gaussian_kl_validated <- function(mu_p, R_p, mu_q, R_q,
                                   decompose = TRUE) {
  Lq <- .strict_chol(R_q, "R_q")
  Lp <- .strict_chol(R_p, "R_p")
  d <- mu_q - mu_p
  w <- forwardsolve(t(Lq), d)
  signal <- 0.5 * sum(w * w)
  A <- forwardsolve(t(Lq), t(Lp))
  dispersion <- 0.5 * (
    sum(A * A) - length(mu_p) +
      2 * sum(log(diag(Lq))) - 2 * sum(log(diag(Lp)))
  )
  signal <- max(signal, 0)
  dispersion <- max(dispersion, 0)
  total <- signal + dispersion
  if (isTRUE(getOption("aci.debug_assert", FALSE)) && total < -1e-10)
    aci_abort("aci_error_internal", "Negative KL encountered.")
  if (isTRUE(decompose))
    c(total = total, signal = signal, dispersion = dispersion) else total
}


#' Relative entropy gained by an update
#'
#' Relative entropy of an updated Gaussian from the one it replaced, oriented so
#' that the new distribution is the integrating density.
#'
#' @param mu_old Numeric vector, mean before the update.
#' @param R_old Covariance matrix before the update.
#' @param mu_new Numeric vector, mean after the update.
#' @param R_new Covariance matrix after the update.
#' @returns 1-length numeric, the relative entropy gained.
#'
#' @seealso [gaussian_kl()]
#'
#' @examples
#' kl_increment(mu_old = 0, R_old = matrix(2),
#'              mu_new = 0.5, R_new = matrix(1))
#'
#' @export
kl_increment <- function(mu_old, R_old, mu_new, R_new)
  unname(gaussian_kl(mu_new, R_new, mu_old, R_old, decompose = FALSE))


#' Gaussian relative entropy along a pair of paths
#'
#' Evaluates the relative entropy at every time of two assimilation paths on a
#' common grid. Gaussian relative entropy is oriented as smoother relative to
#' filter.
#'
#' @param p A `da_path_gaussian` object, the integrating distribution.
#' @param q A `da_path_gaussian` object on the same time grid.
#' @param decompose `TRUE` to return the signal and dispersion parts alongside
#'   the total.
#' @returns A data frame with the time column `t` and either `total` alone or
#'   `total`, `signal` and `dispersion`.
#'
#' @seealso [gaussian_kl()], [aci()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- da_filter(m, ob)
#' s <- da_smooth(m, ob, filter = f)
#' head(gaussian_kl_path(s, f))
#'
#' @export
gaussian_kl_path <- function(p, q, decompose = TRUE) {
  if (!inherits(p, "da_path_gaussian") || !inherits(q, "da_path_gaussian"))
    aci_abort("aci_error_dims", "p and q must be Gaussian assimilation paths.")
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  if (length(p$t) != length(q$t) || max(abs(p$t - q$t)) > 1e-10)
    aci_abort("aci_error_dims", "gaussian_kl_path: mismatched time grids.")
  n <- length(p$t)
  if (!is.matrix(p$mean) || !is.matrix(q$mean))
    aci_abort("aci_error_dims", "gaussian_kl_path: path means must be matrices.")
  lp <- ncol(p$mean); lq <- ncol(q$mean)
  if (n < 1L || lp < 1L || lq != lp ||
      !identical(dim(p$mean), c(n, lp)) ||
      !identical(dim(q$mean), c(n, lp)) ||
      !identical(dim(p$cov), c(lp, lp, n)) ||
      !identical(dim(q$cov), c(lp, lp, n)) ||
      any(!is.finite(c(p$t, q$t, p$mean, q$mean, p$cov, q$cov))))
    aci_abort("aci_error_dims", "gaussian_kl_path: incompatible or non-finite moments.")
  if (lp == 1L) {
    R_p <- as.numeric(p$cov)
    R_q <- as.numeric(q$cov)
    # Preserve gaussian_kl()'s per-index validation order: reference covariance
    # first, integrating covariance second, before vectorised arithmetic.
    for (j in seq_len(n)) {
      if (R_q[j] <= 0)
        aci_abort("aci_error_spd", "Matrix (R_q) must be positive definite.")
      if (R_p[j] <= 0)
        aci_abort("aci_error_spd", "Matrix (R_p) must be positive definite.")
    }
    scalar <- .gaussian_kl_scalar_kernel(
      as.numeric(p$mean), R_p, as.numeric(q$mean), R_q,
      decompose = decompose
    )
    return(data.frame(t = p$t, scalar, check.names = FALSE))
  }
  if (!isTRUE(decompose)) {
    total <- vapply(seq_len(n), function(j)
      .gaussian_kl_validated(
        p$mean[j, ], p$cov[, , j], q$mean[j, ], q$cov[, , j],
        decompose = FALSE
      ), numeric(1))
    return(data.frame(t = p$t, total = total))
  }
  out <- matrix(NA_real_, n, 3, dimnames = list(NULL, c("total", "signal", "dispersion")))
  for (j in seq_len(n))
    out[j, ] <- .gaussian_kl_validated(
      p$mean[j, ], p$cov[, , j], q$mean[j, ], q$cov[, , j]
    )
  data.frame(t = p$t, out)
}


#' Relative entropy of two Gaussians projected onto a direction
#'
#' Projects both distributions onto a unit direction and returns the relative
#' entropy of the resulting one-dimensional Gaussians.
#'
#' @param v Numeric direction vector, normalised internally; must be non-zero.
#' @param mu0 Numeric vector, mean of the reference distribution.
#' @param R0 Covariance matrix of the reference distribution.
#' @param muE Numeric vector, mean of the compared distribution.
#' @param RE Covariance matrix of the compared distribution.
#' @returns 1-length numeric, the projected relative entropy.
#'
#' @seealso [sensitive_directions()]
#'
#' @examples
#' projected_kl(v = c(1, 0), mu0 = c(0, 0), R0 = diag(2),
#'              muE = c(1, 0), RE = diag(c(0.5, 1)))
#'
#' @export
projected_kl <- function(v, mu0, R0, muE, RE) {
  v <- as.numeric(v); nv <- sqrt(sum(v^2))
  if (!length(v) || !is.finite(nv) || nv <= 0)
    aci_abort("aci_error_dims", "v must be a finite non-zero direction vector.")
  mu0 <- as.numeric(mu0); muE <- as.numeric(muE)
  R0 <- as.matrix(R0); RE <- as.matrix(RE); l <- length(v)
  if (length(mu0) != l || length(muE) != l || !is.numeric(R0) ||
      !is.numeric(RE) || !identical(dim(R0), c(l, l)) ||
      !identical(dim(RE), c(l, l)) || any(!is.finite(c(mu0, muE, R0, RE))))
    aci_abort("aci_error_dims", "Projected distributions must match the direction dimension.")
  .strict_chol(R0, "R0")
  .strict_chol(RE, "RE")
  v <- v / nv
  s0 <- drop(t(v) %*% R0 %*% v); sE <- drop(t(v) %*% RE %*% v)
  if (!is.finite(s0) || !is.finite(sE) || s0 <= 0 || sE <= 0)
    aci_abort("aci_error_spd", "Projected variances must be positive.")
  dm <- drop(t(v) %*% (muE - mu0))
  0.5 * (sE / s0 + dm^2 / s0 - 1 + log(s0 / sE))
}


#' Relative entropy between two samples
#'
#' Estimates the relative entropy of one sample from another, either by fitting
#' Gaussian moments or by a nearest-neighbour estimator.
#'
#' @param A Numeric matrix of draws from the first distribution, one row per
#'   draw.
#' @param B Numeric matrix of draws from the second distribution.
#' @param estimator Either `"gaussian"` or `"knn"`.
#' @returns 1-length numeric, the estimated relative entropy.
#'
#' @seealso [gaussian_kl()]
#'
#' @examples
#' set.seed(1)
#' empirical_kl(matrix(rnorm(200), ncol = 2), matrix(rnorm(200), ncol = 2))
#'
#' @export
empirical_kl <- function(A, B, estimator = c("gaussian", "knn")) {
  estimator <- match.arg(estimator)
  if (estimator == "knn")
    aci_abort("aci_error_not_implemented", "knn estimator is a v0.2 stub (SPEC-02).")
  A <- as.matrix(A); B <- as.matrix(B)
  if (nrow(A) < 2L || nrow(B) < 2L || ncol(A) != ncol(B) ||
      ncol(A) < 1L || any(!is.finite(c(A, B))))
    aci_abort("aci_error_dims",
              "A and B must be finite matrices with matching columns and at least two rows each.")
  d <- ncol(A)
  if (nrow(A) <= d || nrow(B) <= d)
    aci_abort("aci_error_ensemble_rank",
              "Gaussian empirical KL requires more samples than dimensions in each sample.")
  RA <- stats::cov(A); RB <- stats::cov(B)
  full_rank <- function(R) !is.null(tryCatch(chol(sym(R)), error = function(e) NULL))
  if (!full_rank(RA) || !full_rank(RB))
    aci_abort("aci_error_ensemble_rank",
              "Gaussian empirical KL requires full-rank sample covariance in both samples.")
  gaussian_kl(colMeans(A), RA, colMeans(B), RB, decompose = FALSE)
}


################################################################################
# ACI / aci()
################################################################################

#' Assimilative causal inference
#'
#' Runs the filter and smoother for a model and observed record and scores each
#' time by the relative entropy of the smoother from the filter. A normal
#' `aci()` call uses the supplied-code backward-ODE headline smoother, including
#' its correlated-noise correction, independently of `keep`. [lag_table()] and
#' `aci(table = ...)` instead use the complete online Theorem 3 smoother;
#' their finite-grid diagonal can therefore differ from headline ACI.
#'
#' @param model A `cgns_model` or `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param engine One of `"auto"`, `"cgns"` or `"ensemble"`; `"auto"` selects the
#'   closed-form engine for a conditional-Gaussian model.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param table Optional precomputed `lag_table`, whose online diagonal is used
#'   in place of the headline smoother.
#' @param keep One of `"paths"`, `"table"` or `"none"`, selecting which objects
#'   are retained on the result. It does not select a smoother.
#' @param decompose `TRUE` to retain the signal and dispersion parts.
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param m Ensemble size, used by the ensemble engine only.
#' @param seed Optional non-negative whole number seeding the ensemble.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param inflation Positive multiplicative variance inflation factor.
#' @param ic_sampler Optional function of the ensemble size drawing the
#'   initial members, or a list with the initial `mean` and `cov`; `NULL`
#'   uses the model's default.
#' @param ... Must be empty; unused arguments are an error.
#' @returns An object of class `aci_result`.
#'
#' @seealso [forward_cir()], [backward_cir()], [lag_table()], [gaussian_kl_path()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' a <- aci(m, ob)
#' a
#'
#' @export
aci <- function(model, obs, engine = c("auto", "cgns", "ensemble"),
                nontarget = NULL, table = NULL,
                keep = c("paths", "table", "none"), decompose = TRUE,
                init = NULL, m = 100, seed = NULL,
                stepper = c("explicit", "implicit"), nsub = 1L,
                localization = NULL, inflation = 1, ic_sampler = NULL, ...) {
  dots <- list(...)
  stepper <- match.arg(stepper); keep <- match.arg(keep)
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  engine <- match.arg(engine)
  obs <- as_obs(obs)
  if (!inherits(model, "stochastic_model") || obs$k != model$k)
    aci_abort("aci_error_model_contract",
              "model must be a stochastic_model matching the observation dimension.")
  if (engine == "auto") engine <- if (inherits(model, "cgns_model")) "cgns" else "ensemble"
  if (!is.null(table)) {
    if (length(dots))
      aci_abort("aci_error_dims", "Unused arguments were supplied with a reused lag table.")
    .validate_lag_table_source(table, model, obs, nontarget, init)
    table_init <- table$meta$init %||% init
    return(structure(list(t = table$t, aci = table$diag,
                          signal = if (decompose) table$diag_signal else NULL,
                          dispersion = if (decompose) table$diag_dispersion else NULL,
                          paths = NULL, table = table,
                          handles = list(model = model, obs = obs,
                                         nontarget = table$meta$nontarget,
                                         init = table_init),
                          meta = list(engine = "reused_table",
                                      nontarget = table$meta$nontarget,
                                      smoother_scheme =
                                        table$meta$reference_smoother %||%
                                        "unspecified_lag_table_reference")),
                     class = "aci_result"))
  }
  bundle <- NULL
  if (engine == "cgns") {
    if (length(dots))
      aci_abort("aci_error_dims", "Unused arguments were supplied to the closed-form ACI engine.")
    if (!inherits(model, "cgns_model"))
      aci_abort("aci_error_model_contract", "engine='cgns' requires a cgns_model.")
    bundle <- .compile_cgns_run(model, obs, nontarget)
    filt <- .cgns_filter_compiled(
      bundle, init = init, stepper = stepper, nsub = nsub, validate = FALSE
    )
    smoo <- .cgns_smoother_compiled(bundle, filt, validate = FALSE)
  } else {
    if (!is.numeric(m) || length(m) != 1L || !is.finite(m) ||
        m != floor(m) || m < 2L)
      aci_abort("aci_error_ensemble_rank", "m must be an integer of at least 2.")
    if (m <= model$l)
      aci_abort("aci_error_ensemble_rank", paste(
        "Full-dimensional Gaussian ACI from ensemble moments requires m > l;",
        "localized EnKBF/EnKBS state estimation may use m <= l, but its raw",
        "sample covariance is singular and needs a separately specified",
        "projection or covariance regularizer for KL."))
    fr <- enkbf(model, obs, m = m, seed = seed, nontarget = nontarget,
                localization = localization, inflation = inflation,
                ic_sampler = ic_sampler, ...)
    sm <- enkbs(fr$model, fr$path, fr$noise, localization = localization)
    filt <- as_gaussian(fr$path); smoo <- as_gaussian(sm)
    aci_warn("aci_warn_ensemble_kl", sprintf(
      "ACI from ensemble moments (m = %d): jiang2026enkbs s3.2 finds m ~ 50 faithful; m = 10 preserves timing/sign but distorts magnitudes.", m))
  }
  klp <- if (engine == "cgns")
    .gaussian_kl_path_compiled(
      bundle, smoo, filt, decompose = decompose, validate = FALSE
    ) else gaussian_kl_path(smoo, filt, decompose = decompose)
  tab <- NULL
  if (keep == "table") {
    tab <- if (engine == "cgns")
      .lag_table_compiled(
        bundle, mode = "forward", filter = filt, init = filt$meta$init,
        validate = FALSE
      ) else
      .ensemble_lag_table_from_run(model, obs, fr, sm,
                                   nontarget = nontarget,
                                   localization = localization)
  }
  actual_init <- if (engine == "cgns") filt$meta$init else init
  structure(list(t = klp$t, aci = klp$total,
                 signal = if (decompose) klp$signal else NULL,
                 dispersion = if (decompose) klp$dispersion else NULL,
                 paths = if (keep != "none") list(filter = filt, smoother = smoo) else NULL,
                 table = tab,
                 handles = list(model = model, obs = obs, nontarget = nontarget,
                                init = actual_init),
                 meta = list(engine = engine, nontarget = nontarget,
                             m = if (engine == "ensemble") m else NULL,
                             smoother_scheme = smoo$meta$route %||%
                               if (engine == "ensemble") "enkbs" else "unspecified",
                             table_reference = if (!is.null(tab))
                               tab$meta$reference_smoother else NULL)),
            class = "aci_result")
}


#' Print aci() result
#'
#' @param x An `aci_result` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.aci_result <- function(x, ...) {
  pk <- which.max(x$aci)
  cat(sprintf("<aci_result> engine = %s%s | peak ACI = %.4g at t = %.4g\n",
              x$meta$engine,
              if (!is.null(x$meta$nontarget)) " (conditional)" else "",
              x$aci[pk], x$t[pk]))
  invisible(x)
}


#' aci() result to data.frame
#'
#' @param x An `aci_result` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns A data frame with one row per time, carrying the metric and, when
#'   retained, its signal and dispersion parts.
#' @export
as.data.frame.aci_result <- function(x, ...) {
  d <- data.frame(t = x$t, aci = x$aci)
  if (!is.null(x$signal)) { d$signal <- x$signal; d$dispersion <- x$dispersion }
  d
}


################################################################################
# forward/backward causal influence ranges (exact + l1_linf)
################################################################################

#' Composite Simpson's rule on a unit-spaced grid (internal)
#'
#' @param y Numeric vector of samples on a uniformly spaced grid.
#' @returns 1-length numeric, the integral in units of the sample spacing.
#' @noRd
.simpson <- function(y, close = c("quadratic", "trapezoid")) {
  close <- match.arg(close)
  n <- length(y)
  if (n < 3L) {
    return(if (n == 2L) 0.5 * sum(y) else 0)
  }
  m <- if (n %% 2L == 1L) n else n - 1L
  i <- seq(1L, m - 2L, by = 2L)
  s <- sum(y[i] + 4 * y[i + 1L] + y[i + 2L]) / 3
  ## Composite Simpson needs an odd point count; an even-length grid leaves
  ## one interval over.  ACI_code-main/simps.m:106-115 closes it by fitting a
  ## quadratic through the last three points and integrating that over the
  ## final interval, which reduces on a uniform grid to (5 y_n + 8 y_{n-1} -
  ## y_{n-2}) / 12.  The "trapezoid" close is the package's pre-0.0.21 rule,
  ## retained so earlier results stay reproducible; the two differ by
  ## -(y_n - 2 y_{n-1} + y_{n-2}) / 12, zero on a flat or linear tail.
  if (n %% 2L == 0L)
    s <- s + if (identical(close, "trapezoid")) 0.5 * (y[n - 1L] + y[n]) else
      (5 * y[n] + 8 * y[n - 1L] - y[n - 2L]) / 12
  return(s)
}


#' Calculate tau using Simpson's rule (internal)
#'
#' Taus is a numeric vector containing the causal influence range.
#'
#' @param dt Positive 1-length numeric step.
#' @param p Numeric vector of finite, non-negative one-lag divergences.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param M Max p or a bias-corrected metric.
#' @param direction Character, "fwd" or "bwd".
#' @returns 1-length numeric `tau`
#' @noRd
.calc_tau <- function(dt, p, method, M, direction, quadrature,
                      simpson_close = "quadratic"){
  if (method != "l1_linf") {
    if(direction == "fwd"){
      p <- rev(cummax(rev(p)))
    } else{ # direction = "bwd"
      p <- rev(cummin(rev(p)))
    }
  }
  ## The two published conventions for the ratio: the plain andreou2026cir
  ## Appendix G L1 grid-function sum (eqs. G.8/G.14, the FBCIR scripts'
  ## active lines) and composite Simpson (the ACI and EnKBS codebases'
  ## active lines). The exact layer-cake form always integrates with
  ## Simpson: it is quadrature of a continuum integral, not a grid norm.
  if (method == "l1_linf" && quadrature == "sum") {
    tau <- dt * sum(p) / M
  } else {
    tau <- dt * .simpson(p, close = simpson_close) / M
  }
  return(tau)
}


#' Forward influence duration from one table row (internal)
#'
#' @param p Numeric vector of finite, non-negative forward divergences.
#' @param dt Positive 1-length numeric step.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @returns Named numeric vector with the range `tau` and the strength `M`.
#' @noRd
.fwd_lengths <- function(p, dt, method, quadrature = "simpson",
                        simpson_close = "quadratic") {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p)) || any(p < 0) || length(dt) != 1L || !is.finite(dt) || dt <= 0){
    aci_abort("aci_error_dims", "Forward discrepancy values must be finite and non-negative, with dt > 0.")
  }
  M <- max(p)
  if (!is.finite(M) || M < 1e-14){ return(c(tau = NA_real_, M = M)) }
  # integrated with composite Simpson, following the ACI reference code; the
  # FBCIR scripts' active ratio lines use the plain andreou2026cir Appendix G
  # L1 sum, with Simpson retained there as a commented alternative
  tau <- .calc_tau(p = p, dt = dt, method = method, M = M, direction = "fwd",
                   quadrature = quadrature, simpson_close = simpson_close)
  return(c(tau = tau, M = M))
}


#' Backward influence duration from the one-lag sequence (internal)
#'
#' @param p Numeric vector of finite, non-negative one-lag divergences.
#' @param dt Positive 1-length numeric step.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @returns Named numeric vector with the range `tau` and the strength `M`.
#' @noRd
.bwd_lengths <- function(p, dt, method, quadrature = "sum",
                        simpson_close = "quadratic") {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p)) || any(p < 0) || length(dt) != 1L || !is.finite(dt) || dt <= 0){
    aci_abort("aci_error_dims", "Backward discrepancy values must be finite and non-negative, with dt > 0.")
  }
  p_bc <- abs(p - p[1])  # bias-corrected backward CIR metric (andreou2026cir eq. 13)
  M <- max(p_bc)
  if (!is.finite(M) || M < 1e-14){ return(c(tau = NA_real_, M = M)) }
  # lagged-grid convention from FBCIR is an O(dt) approximation to T' -> T^- in the paper
  tau <- .calc_tau(p = p_bc, dt = dt, method = method, M = M,
                   direction = "bwd", quadrature = quadrature,
                   simpson_close = simpson_close)
  return(c(tau = tau, M = M))
}


#' Resolve the minimum influence strength (internal)
#'
#' @param min_M Either `"auto"` (reads `aci.cir_min_strength` option) or a finite non-negative number.
#' @returns 1-length numeric strength threshold.
#' @noRd
.cir_min_strength <- function(min_M) {
  if (is.character(min_M) && length(min_M) == 1L &&
      !is.na(min_M) && identical(min_M, "auto")) {
    min_M <- getOption("aci.cir_min_strength", 1e-5)
  }
  if (!is.numeric(min_M) || length(min_M) != 1L ||
      !is.finite(min_M) || min_M < 0)
    aci_abort("aci_error_dims", "min_M must be 'auto' or one finite non-negative number.")
  as.numeric(min_M)
}


#' Construct a causal influence range (internal)
#'
#' @param t Numeric vector of anchor times.
#' @param tau Numeric vector of ranges, in the model's time units.
#' @param M Numeric vector of influence strengths.
#' @param direction Either `"forward"` or `"backward"`.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param dt Positive 1-length numeric step.
#' @param interval Two-column matrix of range endpoints.
#' @param tail_bound Optional heuristic tail estimates carried from the table.
#' @param subjective Optional subjective ranges at the requested thresholds.
#' @param meta Optional named list of metadata carried on the object.
#' @param bound Optional label recording the estimator's discrete semantics.
#' @returns An object of class `cir_result`.
#' @noRd
new_cir_result <- function(t, tau, M, direction, method, dt, interval,
                           tail_bound = NULL, subjective = NULL, meta = list(),
                           bound = NULL) {
  structure(list(t = t, tau = tau, M = M, interval = interval,
                 direction = direction, method = method,
                 bound = bound %||% switch(method, l1_linf =
                   if (direction == "forward") "underestimate" else "overestimate",
                   exact = "layer_cake_objective"),
                 tail_bound = tail_bound, subjective = subjective,
                 dt = dt, meta = meta),
            class = "cir_result")
}


#' Print the causal influence range
#'
#' @param x A `cir_result` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.cir_result <- function(x, ...) {
  cat(sprintf("<cir_result> %s | method = %s (%s) | masked/NA: %d of %d\n",
              x$direction, x$method, x$bound, sum(is.na(x$tau)), length(x$tau)))
  ab <- x$meta$above_baseline
  if (is.null(ab) && !is.null(x$meta$per_reference))
    ab <- vapply(x$meta$per_reference,
                 function(mm) isTRUE(mm$above_baseline), logical(1))
  if (!is.null(ab) && length(ab))
    cat(sprintf("  above initial-time baseline: %d of %d anchors\n",
                sum(ab), length(ab)))
  invisible(x)
}


#' Forward causal influence range
#'
#' Summarizes the duration of influence on the discrete time grid, forward from
#' each anchor time. A finite adaptive table is labelled
#' `objective_on_truncated_table`; its `tail_bound` field is a heuristic tail
#' estimate and must not be interpreted as a certified error bound. The
#' `l1_linf` estimators are ratios, integrated with composite Simpson over the
#' whole span, following the ACI reference code.
#'
#' @param x A `lag_table` or `aci_result` object.
#' @param ... Arguments passed to methods.
#' @returns An object of class `cir_result` with `direction` `"forward"`.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' @seealso [backward_cir()], [lag_table()], [cir_pair()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' forward_cir(tb)
#'
#' @export
forward_cir <- function(x, ...) UseMethod("forward_cir")


#' @describeIn forward_cir Forward range from a precomputed lag table.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param eps Optional numeric vector of finite non-negative thresholds at which
#'   the subjective range is also reported. The subjective length follows
#'   andreou2026cir eq. G.7, the lag time of the last exceedance; the ACI
#'   reference script instead multiplies the 1-based cell index by `dt`, so
#'   its read-out sits one grid step above this one.
#' @param min_M Either `"auto"` or one finite non-negative number; ranges whose
#'   strength falls below it are masked.
#' @param masked_value What a masked range reports: `"na"` (the default) keeps
#'   the masking visible as `NA`; `"zero"` follows andreou2026cir Remarks B.4
#'   and C.4 and the FBCIR scripts, which set the length to 0. The paper's
#'   forward figures are computed on an untruncated table; the adaptive
#'   truncation of the default `lag_table()` is a package storage device whose
#'   use is recorded in the result's `bound` label, and `tol = 0` reproduces
#'   the untruncated convention.
#' @param quadrature Quadrature for the `l1_linf` ratio: `"simpson"` (the
#'   default) follows the ACI and EnKBS codebases' active lines, while
#'   `"sum"` is the literal andreou2026cir eq. G.8 L1 grid-function sum,
#'   which the FBCIR scripts' active lines use. The choice is recorded in the
#'   result's `meta$quadrature`; the `exact` form always integrates with
#'   Simpson.
#' @param simpson_close Closing rule for the leftover interval when a grid has
#'   an even number of points: `"quadratic"` (the default) fits a quadratic
#'   through the last three points, following `simps.m` in the ACI reference
#'   code; `"trapezoid"` is the package's pre-0.0.21 rule and reproduces
#'   results reported by earlier versions. Odd-length grids are unaffected.
#' @export
forward_cir.lag_table <- function(x, method = c("exact", "l1_linf"),
                                  eps = NULL, min_M = "auto",
                                  masked_value = c("na", "zero"),
                                  quadrature = c("simpson", "sum"),
                                  simpson_close = c("quadratic", "trapezoid"),
                                  ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to forward_cir().")
  masked_value <- match.arg(masked_value)
  quadrature <- match.arg(quadrature)
  simpson_close <- match.arg(simpson_close)
  method <- match.arg(method)
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  if (!x$mode %in% c("forward", "full"))
    aci_abort("aci_error_dims", "forward_cir needs a 'forward' or 'full' mode table.")
  N1 <- length(x$t)
  structurally_truncated <- any(vapply(seq_len(N1), function(j)
    length(x$rows[[j]]) < (N1 - j + 1L), logical(1)))
  tau <- M <- rep(NA_real_, N1); subj <- NULL
  for (j in seq_len(N1)) {
    p <- lt_row(x, j, pad = "zero")
    r <- .fwd_lengths(p, x$dt, method, quadrature, simpson_close)
    tau[j] <- r["tau"]; M[j] <- r["M"]
  }
  if (!is.null(eps)) {
    if (!is.numeric(eps) || !length(eps) || any(!is.finite(eps)) || any(eps < 0))
      aci_abort("aci_error_dims", "eps must contain finite non-negative thresholds.")
    ## (idx - 1) * dt is the lag time of the last exceedance, andreou2026cir
    ## eq. G.7 with the first cell at lag 0.  The reference script's
    ## subj_CIR_idx*dt counts cells and reads one grid step higher.
    subj <- sapply(eps, function(e) sapply(seq_len(N1), function(j) {
      p <- lt_row(x, j, pad = "zero"); idx <- which(p > e)
      if (!length(idx)) 0 else x$dt * (max(idx) - 1) }))
  }
  if (!is.null(min_M)) {
    mm <- .cir_min_strength(min_M)
    bad <- which(M < mm)
    ## The published convention (andreou2026cir Rmks B.4 and C.4, and the
    ## FBCIR scripts) sets a low-strength length to 0; the package default
    ## reports NA so that masking stays visible in downstream summaries.
    if (length(bad)) { tau[bad] <- if (masked_value == "zero") 0 else NA_real_
      aci_warn("aci_warn_low_signal", sprintf(
        "%d forward CIR values masked (M < %.3g); interpret CIRs jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).",
        length(bad), mm)) }
  }
  new_cir_result(x$t, tau, M, "forward", method, x$dt,
                 interval = cbind(x$t, x$t + tau),
                 tail_bound = x$tailbnd, subjective = subj,
                 meta = list(quadrature = quadrature,
                             simpson_close = simpson_close,
                             table_nontarget = x$meta$nontarget,
                             structurally_truncated = structurally_truncated),
                 bound = if (structurally_truncated) {
                   if (method == "exact") "objective_on_truncated_table"
                   else "lower_ratio_on_truncated_table_only"
                 } else NULL)
}


#' @describeIn forward_cir Forward range from an ACI result, building the table
#'   when it was not retained.
#' @export
forward_cir.aci_result <- function(x, ...) {
  if (is.null(x$table) && identical(x$meta$engine, "ensemble"))
    aci_abort("aci_error_not_implemented", paste(
      "The forward ensemble CIR needs the lagged EnKBS family; recompute",
      "aci(..., engine = 'ensemble', keep = 'table')."))
  if (!is.null(x$table)) return(forward_cir(x$table, ...))
  bundle <- .compile_cgns_run(
    x$handles$model, x$handles$obs, x$handles$nontarget
  )
  .forward_cir_compiled(
    bundle, filter = x$paths$filter, init = x$handles$init, ...
  )
}


#' Backward causal influence range
#'
#' Summarizes the duration of influence on the discrete time grid, backward from
#' a reference time. A backward result is anchored at `T`, but its computable
#' lagged grid and physical interval end at `T - dt`; its bound label records
#' this convention. The backward `l1_linf` ratio is evaluated with the plain
#' Appendix G L1 sum, following the FBCIR code's active line; the exact form
#' integrates its suffix minima with composite Simpson.
#'
#' A backward result also reports the andreou2026cir Section 2.3.4 validity
#' gate in its `meta`: `baseline` (the initial-time information deficit, the
#' one-lag metric at t = 0), `terminal` (the metric at the `T - dt` end of the
#' grid), and `above_baseline` (their comparison, eq. 21 and Remark C.2).
#' A large `M` with `above_baseline = FALSE` means the attribution evidence
#' never rises above the initial-time baseline; treat such a range with
#' caution rather than as a long attribution window.
#'
#' @param x A `lag_table`, `aci_result` or `cgns_model` object.
#' @param ... Arguments passed to methods.
#' @returns An object of class `cir_result` with `direction` `"backward"`.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' @seealso [forward_cir()], [lag_table()], [cir_pair()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "one_lag")
#' backward_cir(tb)
#'
#' @export
backward_cir <- function(x, ...) UseMethod("backward_cir")


#' @describeIn backward_cir Backward range from a one-lag table.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param eps Optional numeric vector of finite non-negative thresholds at which
#'   the subjective range is also reported.
#' @param min_M Either `"auto"` or one finite non-negative number; a range whose
#'   strength falls below it is masked.
#' @param masked_value What a masked range reports: `"na"` (the default) keeps
#'   the masking visible as `NA`; `"zero"` follows andreou2026cir Remarks B.4
#'   and C.4 and the FBCIR scripts, which set the length to 0.
#' @param quadrature Quadrature for the `l1_linf` ratio: `"sum"` (the default)
#'   is the literal andreou2026cir eq. G.14 L1 grid-function sum, the FBCIR
#'   scripts' active line; `"simpson"` is the smoother alternative those
#'   scripts retain in comments. The choice is recorded in the result's
#'   `meta$quadrature`; the `exact` form always integrates with Simpson.
#' @param simpson_close Closing rule for the leftover interval when a grid has
#'   an even number of points: `"quadratic"` (the default) fits a quadratic
#'   through the last three points, following `simps.m` in the ACI reference
#'   code; `"trapezoid"` is the package's pre-0.0.21 rule and reproduces
#'   results reported by earlier versions. Odd-length grids are unaffected.
#' @export
backward_cir.lag_table <- function(x, method = c("exact", "l1_linf"),
                                   eps = NULL, min_M = "auto",
                                   masked_value = c("na", "zero"),
                                   quadrature = c("sum", "simpson"),
                                   simpson_close = c("quadratic", "trapezoid"),
                                   ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to backward_cir().")
  masked_value <- match.arg(masked_value)
  quadrature <- match.arg(quadrature)
  simpson_close <- match.arg(simpson_close)
  method <- match.arg(method)
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  P <- lt_onelag(x)
  r <- .bwd_lengths(P, x$dt, method, quadrature, simpson_close)
  ## andreou2026cir Section 2.3.4 validity gate: causal intensity at the
  ## T - dt end of the grid versus the initial-time information deficit
  ## P^0_T (eq. 21; Remark C.2).  The gate is a diagnostic beside tau and M,
  ## not a change to the centred metric.
  baseline <- unname(P[1L])
  terminal <- unname(P[length(P)])
  Tend <- x$t[length(x$t)]
  lag_end <- Tend - x$dt
  subj <- NULL
  if (!is.null(eps)) {
    if (!is.numeric(eps) || !length(eps) || any(!is.finite(eps)) || any(eps < 0))
      aci_abort("aci_error_dims", "eps must contain finite non-negative thresholds.")
    cc <- abs(P - P[1]); pt <- x$t[seq_along(P)]
    subj <- vapply(eps, function(e) {
      hit <- which(cc <= e)
      if (!length(hit)) NA_real_ else lag_end - max(pt[hit])
    }, numeric(1))
    names(subj) <- format(eps)
  }
  if (!is.null(min_M) && is.finite(r["M"])) {
    mm <- .cir_min_strength(min_M)
    if (r["M"] < mm) {
      r["tau"] <- if (masked_value == "zero") 0 else NA_real_
      aci_warn("aci_warn_low_signal", sprintf(
        "Backward CIR masked (M < %.3g); interpret it jointly with ACI.", mm))
    }
  }
  new_cir_result(Tend, unname(r["tau"]), unname(r["M"]), "backward", method,
                 x$dt, interval = cbind(lag_end - unname(r["tau"]), lag_end),
                 subjective = subj,
                 meta = list(quadrature = quadrature,
                             simpson_close = simpson_close,
                             reference_time = Tend,
                             lagged_grid_end = lag_end,
                             endpoint_convention = "T_minus_dt",
                             dt_caveat = "one-lag limit approximated to O(dt) (andreou2026cir eq. G.12)",
                             baseline = baseline,
                             terminal = terminal,
                             above_baseline = terminal > baseline,
                             table_nontarget = x$meta$nontarget),
                 bound = if (method == "exact")
                   "layer_cake_on_O(dt)_T_minus_dt_grid" else
                   "upper_ratio_on_O(dt)_T_minus_dt_grid")
}


#' @describeIn backward_cir Backward range from an ACI result, at one or more
#'   reference times.
#' @param T Either `"end"` or numeric reference times within the observed grid.
#' @export
backward_cir.aci_result <- function(x, T = "end", method = "exact", ...) {
  if (identical(x$meta$engine, "ensemble") ||
      identical(x$table$meta$engine %||% NULL, "ensemble"))
    aci_abort("aci_error_not_implemented", paste(
      "jiang2026enkbs implements the EnKBS forward CIR but identifies ensemble backward",
      "CIR as future work; this package does not invent that extension."))
  mdl <- x$handles$model; ob <- x$handles$obs; nt <- x$handles$nontarget
  Ts <- if (identical(T, "end")) ob$t[length(ob$t)] else as.numeric(T)
  if (!length(Ts) || any(!is.finite(Ts)) ||
      any(Ts < ob$t[1] - 1e-12) || any(Ts > ob$t[length(ob$t)] + 1e-12))
    aci_abort("aci_error_dims",
              "T must contain finite reference times within the observation grid span.")
  if (length(Ts) > 20)
    aci_warn("aci_warn_truncation",
             "backward_cir over >20 reference times is O(|T| N l^3) in v0.")
  prefix_lengths <- vapply(Ts, function(Tv)
    sum(ob$t <= Tv + 1e-12), integer(1))
  if (any(prefix_lengths < 5L))
    aci_abort("aci_error_dims", "Reference time too early.")
  # Coefficients after the latest requested reference time are not part of any
  # backward-CIR estimand. Compile the maximal required prefix once so multiple
  # references share work without evaluating irrelevant future coefficients.
  max_N1 <- max(prefix_lengths)
  prefix_obs <- if (max_N1 == length(ob$t)) ob else observed_trajectory(
    ob$t[seq_len(max_N1)], ob$x[seq_len(max_N1), , drop = FALSE],
    names = colnames(ob$x)
  )
  bundle <- .compile_cgns_run(mdl, prefix_obs, nt)
  full_filter <- x$paths$filter %||% NULL
  can_reuse <- !is.null(full_filter) &&
    nrow(full_filter$mean) >= max_N1 &&
    identical(full_filter$meta$stepper %||% "explicit", "explicit") &&
    (full_filter$meta$nsub %||% 1L) == 1L
  if (can_reuse) {
    full_filter <- .slice_compiled_filter(full_filter, bundle)
    compatible <- tryCatch({
      .validate_gaussian_path(
        full_filter, bundle$obs, bundle$l, "filter", bundle$nontarget,
        model = bundle$model, source_model = bundle$source_model
      )
      TRUE
    }, error = function(e) FALSE)
    can_reuse <- isTRUE(compatible)
  }
  if (!can_reuse)
    full_filter <- .cgns_filter_compiled(
      bundle, init = x$handles$init, stepper = "explicit", nsub = 1L,
      validate = FALSE
    )
  res <- lapply(prefix_lengths, function(N1) {
    prefix <- .slice_compiled_cgns(bundle, N1)
    prefix_filter <- .slice_compiled_filter(full_filter, prefix)
    tb <- .lag_table_compiled(
      prefix, mode = "one_lag", filter = prefix_filter,
      init = prefix_filter$meta$init, validate = FALSE
    )
    backward_cir(tb, method = method, ...)
  })
  actual_T <- vapply(res, function(z) as.numeric(z$t)[1], numeric(1))
  if (any(abs(actual_T - Ts) > 1e-10))
    aci_warn("aci_warn_grid_snap",
             "Some requested backward-CIR reference times were snapped down to the observation grid.")
  if (length(res) == 1) return(res[[1]])
  subj <- lapply(res, `[[`, "subjective")
  if (all(vapply(subj, is.null, logical(1)))) subj <- NULL
  common_bound <- unique(vapply(res, function(z) z$bound, character(1)))
  new_cir_result(actual_T, sapply(res, `[[`, "tau"), sapply(res, `[[`, "M"),
                 "backward", method, ob$dt,
                 interval = do.call(rbind, lapply(res, `[[`, "interval")),
                 subjective = subj,
                 meta = list(per_reference = lapply(res, `[[`, "meta")),
                 bound = if (length(common_bound) == 1L) common_bound else
                   "mixed_reference_semantics")
}


#' Forward and backward causal influence ranges together
#'
#' Runs the metric once and returns the forward and backward ranges alongside
#' it, so that both directions share one filter pass and one prior.
#'
#' @param model A `cgns_model` object; the closed-form engine is required.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param engine One of `"auto"` or `"cgns"`.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param ... Passed to [aci()].
#' @returns An object of class `cir_pair`, a list with the `aci` result and the
#'   `forward` and `backward` ranges.
#'
#' @seealso [forward_cir()], [backward_cir()]
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' cir_pair(m, ob)
#' }
#'
#' @export
cir_pair <- function(model, obs, engine = "auto", nontarget = NULL, ...) {
  resolved_engine <- if (identical(engine, "auto")) {
    if (inherits(model, "cgns_model")) "cgns" else "ensemble"
  } else engine
  if (!identical(resolved_engine, "cgns"))
    aci_abort("aci_error_not_implemented",
              "cir_pair() requires the closed-form CGNS lag table engine.")
  a <- aci(model, obs, engine = engine, nontarget = nontarget,
           keep = "table", ...)
  fwd_tab <- a$table
  bwd_tab <- lag_table(model, as_obs(obs), mode = "one_lag", nontarget = nontarget,
                       filter = a$paths$filter,
                       init = a$handles$init)
  structure(list(aci = a, forward = forward_cir(fwd_tab),
                 backward = backward_cir(bwd_tab)), class = "cir_pair")
}


# Convenience method: backward CIR straight from (model, obs) without first
# assembling an aci_result (the front door for fitted models).
#' @describeIn backward_cir Backward range straight from a model and its
#'   observations, without first assembling an ACI result.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @export
backward_cir.cgns_model <- function(x, obs, T, nontarget = NULL,
                                    init = NULL, ...) {
  ar <- structure(list(handles = list(model = x, obs = as_obs(obs),
                                      nontarget = nontarget, init = init)),
                  class = "aci_result")
  backward_cir(ar, T = T, ...)
}


################################################################################
# tidy CIR read-out
################################################################################

#' Forward and backward causal influence ranges as one tidy frame
#'
#' Wraps the lag table / CIR machinery and returns a data.frame with one row per
#' anchor: `t`, `tau` (the range, in the model's time units), `strength` (the
#' influence strength there) and `direction`. `nontarget` gives the conditional
#' range. Both directions receive `init`, so the backward pass runs from the same
#' prior as the forward one.
#'
#' @param model A CGNS model.
#' @param obs Observed trajectory.
#' @param init Optional Gaussian initialization.
#' @param at Optional anchor times (used to select forward rows and/or backward
#'   reference times).
#' @param direction Direction or directions to calculate.
#' @param nontarget Optional conditional ACI masking specification.
#' @param tol Adaptive-lag tolerance.
#' @param min_M Minimum discrepancy strength; `"auto"` uses the package
#'   threshold option and masks low-signal ranges.
#' @param method Exact layer-cake objective or the `l1_linf` ratio.
#' @param max_lag Maximum number of positive-lag cells retained per forward row.
#' @param window Consecutive small-tail window used by adaptive storage.
#' @param ... Unused arguments are rejected.
#' @returns A data frame with one row per anchor and direction, carrying `t`,
#'   `tau`, `strength` and `direction`; backward rows also carry
#'   `above_baseline`, the [backward_cir()] validity gate (`NA` on forward
#'   rows).
#'
#' @seealso [forward_cir()], [backward_cir()], [cir_pair()]
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' head(cir_table(m, ob, direction = "forward"))
#' }
#'
#' @export
cir_table <- function(model, obs, init = NULL, at = NULL,
                      direction = c("forward", "backward", "both"),
                      nontarget = NULL, tol = 1e-4, min_M = "auto",
                      method = c("exact", "l1_linf"), max_lag = Inf,
                      window = 3L, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to cir_table().")
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "cir_table() needs a cgns_model.")
  direction <- match.arg(direction)
  method <- match.arg(method)
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    aci_abort("aci_error_dims", "tol must be one finite non-negative number.")
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window < 1L || window != floor(window))
    aci_abort("aci_error_dims", "window must be a positive integer.")
  if (!is.numeric(max_lag) || length(max_lag) != 1L || is.na(max_lag) ||
      max_lag < 1 || (!is.infinite(max_lag) && max_lag != floor(max_lag)))
    aci_abort("aci_error_dims", "max_lag must be a positive integer or Inf.")
  # Validate even on a direction that happens not to use a given tuning
  # argument, so malformed inputs are never silently ignored.
  .cir_min_strength(min_M)
  if (!is.null(at) && (!is.numeric(at) || !length(at) || any(!is.finite(at)) ||
      any(at < min(obs$t) | at > max(obs$t))))
    aci_abort("aci_error_dims", "at must contain finite times within the observation span.")
  out <- list(); forward_table <- NULL
  if (direction %in% c("forward", "both")) {
    tab <- lag_table(model, obs, mode = "forward", tol = tol, init = init,
                     nontarget = nontarget, max_lag = max_lag,
                     window = window)
    forward_table <- tab
    f   <- suppressWarnings(forward_cir(tab, method = method, min_M = min_M))
    idx <- if (is.null(at)) seq_along(f$t) else
      vapply(at, function(a) which.min(abs(f$t - a)), integer(1))
    out$forward <- data.frame(t = f$t[idx], tau = f$tau[idx], strength = f$M[idx],
                              direction = "forward", bound = f$bound,
                              tail_estimate = f$tail_bound[idx],
                              above_baseline = NA)
    attr(out$forward, "lag_table") <- tab
  }
  if (direction %in% c("backward", "both")) {
    Ts <- if (is.null(at)) as.numeric(stats::quantile(obs$t, c(.25, .40, .55, .70, .85))) else at
    b  <- suppressWarnings(backward_cir(model, obs = obs, T = Ts, init = init,
                                        nontarget = nontarget, method = method,
                                        min_M = min_M))
    ab <- if (!is.null(b$meta$per_reference))
      vapply(b$meta$per_reference,
             function(mm) isTRUE(mm$above_baseline), logical(1))
    else rep(isTRUE(b$meta$above_baseline), length(b$t))
    out$backward <- data.frame(t = b$t, tau = b$tau, strength = b$M,
                               direction = "backward", bound = b$bound,
                               tail_estimate = NA_real_,
                               above_baseline = ab)
  }
  res <- do.call(rbind, unname(out))
  rownames(res) <- NULL
  prov <- list(
    source_model = model,
    source_obs_t = obs$t,
    source_obs_x = obs$x,
    nontarget = nontarget,
    actual_init = if (!is.null(forward_table)) forward_table$meta$init else init,
    direction = direction,
    method = method,
    tol = tol,
    min_M = min_M,
    max_lag = max_lag,
    window = as.integer(window),
    source_status = paste(
      "Experimental tidy/application wrapper; the underlying CIR calculations",
      "are provided by the ACI engine, but this wrapper is outside the supplied",
      "papers and MATLAB reference code."))
  attr(res, "aci_provenance") <- prov
  if (!is.null(forward_table)) attr(res, "lag_table") <- forward_table
  class(res) <- c("aci_cir_table", "data.frame")
  res
}
