# Filter, smoother and metric for a conditional Gaussian system ----------------
#
# The three functions below are the numerical core of assimilative causal
# inference. They operate on a general CGNS "components" list rather than on a
# specific model, so any conditional Gaussian nonlinear system can be filtered,
# smoothed and scored by supplying the appropriate per-step coefficients.
#
# The symbols follow the governing equations of the method paper rather than R
# naming convention. They are L_x, f_x, L_y and f_y for the drift coefficients,
# S_xoS_x and its relatives for the noise Grammians, and A_j and B_j for the
# smoother's per-step terms. The correspondence to the equations is the point,
# and the package's lint configuration admits these names deliberately.

#' Conditional Gaussian components list
#'
#' Several functions in aciR share a `comp` argument, a named list holding the
#' per-step coefficients of a conditional Gaussian nonlinear system (CGNS). For
#' an observed process `x` and an unobserved (conditionally Gaussian) process
#' `y`, the coefficients are those of the pair of stochastic differential
#' equations governing `x` and `y`. The list has the following entries.
#'
#' \describe{
#'   \item{`L_x`}{Numeric vector, length `n`. The coupling of the unobserved
#'     component into the drift of the observed process at each time step.}
#'   \item{`f_x`}{Numeric vector, length `n`. The remaining drift of the
#'     observed process at each time step.}
#'   \item{`L_y`}{Numeric scalar, or numeric vector of length `n`. The linear
#'     self-drift of the unobserved component. A scalar is a self-drift
#'     constant in time, as in the supplied dyad model; a vector carries one
#'     value per observation, as the supplied predator-prey model needs. It may
#'     depend on the observed signal but never on the unobserved component,
#'     which is what keeps the system conditionally Gaussian.}
#'   \item{`f_y`}{Numeric vector, length `n`. The remaining drift of the
#'     unobserved process at each time step.}
#'   \item{`S_xoS_x`}{Numeric scalar. The observation-noise covariance,
#'     the product of the observed-process noise coefficient with its
#'     transpose. Must be strictly positive.}
#'   \item{`S_yoS_y`}{Numeric scalar. The latent-noise covariance of the
#'     unobserved process. Must be non-negative.}
#'   \item{`S_yoS_x`}{Numeric scalar. The latent-to-observation noise
#'     cross-covariance.}
#'   \item{`S_xoS_y`}{Numeric scalar. The observation-to-latent noise
#'     cross-covariance, the transpose of `S_yoS_x`. For the scalar systems
#'     this package integrates the two are equal, and a components list in
#'     which they disagree is rejected.}
#' }
#'
#' The joint noise covariance must be positive semidefinite, which for a
#' scalar system means `S_xoS_x * S_yoS_y - S_yoS_x^2` is non-negative. A
#' components list that violates any of these conditions is rejected before
#' the recursion starts.
#'
#' This schema is the package's extension surface for advanced use. Build a
#' components list directly to run the core on a conditional Gaussian system
#' for which aciR supplies no constructor. See [aci_dyad_components()] for a
#' worked example, and [aci_cgns_model()] for the higher-level alternative.
#'
#' @name aci_components
#' @keywords internal
NULL

#' Forward conditional Gaussian filter
#'
#' Runs the forward filter of a conditional Gaussian nonlinear system. Given an
#' observed signal and the system components, it returns the filtered mean and
#' covariance of the unobserved component at each time step, that is the mean
#' and covariance of the unobserved state conditional on the observed path up
#' to and including the current step.
#'
#' The recursion is the closed-form conditional Gaussian filter and is
#' integrated with a first-order (Euler) step of width `dt`. The observed
#' signal is assumed to be complete and sampled on a regular grid of spacing
#' `dt`.
#'
#' The filtered covariance is checked at every step. An explicit Euler scheme
#' can drive the covariance non-positive when `dt` is too large for the
#' system, even for a perfectly admissible model, and the relative entropy
#' that scores the result has no meaning in that state. The recursion
#' therefore stops at the first offending step and names it, rather than
#' returning a trajectory that looks like a result.
#'
#' @param x Numeric vector. The observed signal, one value per time step; at
#'   least two complete, finite observations.
#' @param comp A conditional Gaussian components list; see [aci_components].
#' @param dt Numeric scalar. The integration time step; must be positive.
#' @param mu0 Numeric scalar. The initial filtered mean of the unobserved
#'   component.
#' @param R0 Numeric scalar. The initial filtered covariance of the unobserved
#'   component; must be positive.
#'
#' @returns A list with two numeric vectors, `mean` and `cov`, the filtered
#'   mean and covariance of the unobserved component at each time step.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 2000, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
#' str(filt)
#'
#' @seealso [aci_smoother()], [aci_metric()], [aci()]
#' @export
aci_filter <- function(x, comp, dt, mu0, R0) {
  if (.aci_is_mv(comp)) {
    x <- .aci_check_signal_mv(x)
    .aci_check_positive(dt, "dt")
    comp <- .aci_check_components_mv(comp, x)
    mu0 <- as.numeric(mu0)
    R0 <- as.matrix(R0)
    if (length(mu0) != nrow(R0) || nrow(R0) != ncol(R0)) {
      stop(
        "`mu0` and `R0` must share the unobserved dimension, with `R0` square.",
        call. = FALSE
      )
    }
    if (is.null(.aci_chol(R0))) {
      stop("`R0` must be symmetric positive definite.", call. = FALSE)
    }
    return(.aci_filter_mv(x, comp, dt, mu0, R0))
  }
  .aci_check_signal(x)
  n <- length(x)
  .aci_check_components(comp, n)
  .aci_check_positive(dt, "dt")
  .aci_check_scalar(mu0, "mu0")
  .aci_check_positive(R0, "R0")

  # ---- Forward recursion ----------------------------------------------------
  m <- numeric(n)
  v <- numeric(n)
  m[1L] <- mu0
  v[1L] <- R0
  inv <- 1 / comp$S_xoS_x
  L_y <- .aci_expand(comp$L_y, n)
  mu <- mu0
  R <- R0
  for (j in seq_len(n - 1L) + 1L) {
    dx <- x[j] - x[j - 1L]
    aux <- comp$S_yoS_x + R * comp$L_x[j - 1L]
    mu <- mu + (L_y[j - 1L] * mu + comp$f_y[j - 1L]) * dt +
      aux * inv * (dx - (comp$L_x[j - 1L] * mu + comp$f_x[j - 1L]) * dt)
    R <- R + (2 * L_y[j - 1L] * R + comp$S_yoS_y - aux * inv * aux) * dt
    if (!is.finite(R) || R <= 0) {
      .aci_stop_covariance("filter", j, (j - 1L) * dt, R)
    }
    m[j] <- mu
    v[j] <- R
  }
  list(mean = m, cov = v)
}

#' Backward conditional Gaussian smoother
#'
#' Runs the backward smoother of a conditional Gaussian nonlinear system.
#' Starting from the final filtered estimate, it sweeps backward through time
#' and returns the smoothed mean and covariance of the unobserved component,
#' that is the mean and covariance of the unobserved state conditional on the
#' whole observed path.
#'
#' The recursion is the closed-form conditional Gaussian smoother and is
#' integrated with a first-order (Euler) step of width `dt`. It consumes the
#' filtered trajectory returned by [aci_filter()], and, like the filter,
#' stops at the first step at which the smoothed covariance leaves its domain.
#'
#' At the final index the smoother is the filter by construction, because
#' conditioning on the whole observed path and conditioning on the path up to
#' the final step are the same conditioning. The returned trajectory reproduces
#' that identity exactly.
#'
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param comp A conditional Gaussian components list; see [aci_components].
#' @param dt Numeric scalar. The integration time step; must be positive.
#' @param filt A list with numeric vectors `mean` and `cov`, as returned by
#'   [aci_filter()].
#'
#' @returns A list with two numeric vectors, `mean` and `cov`, the smoothed
#'   mean and covariance of the unobserved component at each time step.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 2000, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
#' smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
#' str(smooth)
#'
#' @seealso [aci_filter()], [aci_metric()], [aci()]
#' @export
aci_smoother <- function(x, comp, dt, filt) {
  if (.aci_is_mv(comp)) {
    x <- .aci_check_signal_mv(x)
    .aci_check_positive(dt, "dt")
    comp <- .aci_check_components_mv(comp, x)
    if (!is.list(filt) || !is.matrix(filt$mean) ||
          length(dim(filt$cov)) != 3L) {
      stop(
        paste0(
          "`filt` must be a vector-valued posterior, with a matrix `mean` ",
          "and a three-dimensional array `cov`, as returned by ",
          "`aci_filter()` on a vector system."
        ),
        call. = FALSE
      )
    }
    return(.aci_smoother_mv(x, comp, dt, filt))
  }
  .aci_check_signal(x)
  n <- length(x)
  .aci_check_components(comp, n)
  .aci_check_positive(dt, "dt")
  .aci_check_posterior(filt, "filt", n)

  # ---- Backward recursion ---------------------------------------------------
  m <- numeric(n)
  v <- numeric(n)
  m[n] <- filt$mean[n]
  v[n] <- filt$cov[n]
  inv <- 1 / comp$S_xoS_x
  L_y <- .aci_expand(comp$L_y, n)
  muT <- m[n]
  RT <- v[n]
  B_j <- comp$S_yoS_y - comp$S_yoS_x * inv * comp$S_xoS_y
  for (j in rev(seq_len(n - 1L))) {
    dx <- x[j + 1L] - x[j]
    A_j <- L_y[j] - comp$S_yoS_x * inv * comp$L_x[j]
    # The backward mean carries two contributions, the reversed prior drift
    # corrected toward the filtered estimate, and the transport term through
    # which the noise cross-covariance enters.
    drift <- L_y[j] * muT + comp$f_y[j] -
      B_j / filt$cov[j] * (filt$mean[j] - muT)
    transport <- comp$S_yoS_x * inv *
      (-dx + (comp$L_x[j] * muT + comp$f_x[j]) * dt)
    muT <- muT - drift * dt + transport
    RT <- RT - (2 * (A_j + B_j / filt$cov[j]) * RT - B_j) * dt
    if (!is.finite(RT) || RT <= 0) {
      .aci_stop_covariance("smoother", j, (j - 1L) * dt, RT)
    }
    m[j] <- muT
    v[j] <- RT
  }
  list(mean = m, cov = v)
}

#' Assimilative causal-information metric
#'
#' Computes the assimilative causal-information metric at each time step. The
#' metric is the relative entropy (Kullback-Leibler divergence) of the smoother
#' posterior of the unobserved component from its filter posterior, with the
#' smoother posterior as the integrating density. A larger value marks a step
#' at which conditioning on the future observations sharpens the estimate of
#' the unobserved component more strongly, and the metric is non-negative.
#'
#' For scalar conditional Gaussian posteriors the relative entropy separates
#' into a signal part, driven by the shift between the smoothed and filtered
#' means relative to the filtered covariance, and a dispersion part, driven by
#' the ratio of the smoothed to the filtered covariance.
#'
#' The dispersion part is evaluated in a form that stays accurate when the two
#' covariances nearly agree. Writing the covariance ratio as `1 + d`, the
#' direct expression subtracts two nearly equal quantities and loses precision
#' exactly where the metric is smallest; the form used here does not. Values
#' that round to a small negative number are clamped to zero, and anything
#' more negative than round-off is an error rather than a result. The number
#' of values clamped this way is recorded by [aci()] and reported by
#' [summary.aci()]; a value that is exactly zero because the posteriors agree,
#' as at the final step, needs no clamp and is not counted.
#'
#' @param filt A list with numeric vectors `mean` and `cov`, the filtered mean
#'   and covariance, as returned by [aci_filter()].
#' @param smooth A list with numeric vectors `mean` and `cov`, the smoothed
#'   mean and covariance, as returned by [aci_smoother()].
#'
#' @returns A numeric vector of the causal-information metric at each time
#'   step, non-negative throughout.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 2000, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
#' smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
#' summary(aci_metric(filt, smooth))
#'
#' @seealso [aci_filter()], [aci_smoother()], [aci()]
#' @export
aci_metric <- function(filt, smooth) {
  .aci_metric_pair(filt, smooth)$value
}

#' Metric values and clamp count for a filter/smoother pair
#'
#' The full metric computation behind [aci_metric()], returning the clamp
#' count alongside the values so that [aci()] can record it. The public
#' function returns only the values; keeping a single implementation of the
#' relative entropy is what rules out drift between the two callers.
#'
#' @param filt A list with numeric vectors `mean` and `cov`, as returned by
#'   [aci_filter()].
#' @param smooth A list with numeric vectors `mean` and `cov`, as returned by
#'   [aci_smoother()].
#'
#' @returns A list with the numeric metric vector `value` and the integer
#'   clamp count `n_clamped`.
#'
#' @noRd
#' @keywords internal
.aci_metric_pair <- function(filt, smooth) {
  # A vector-valued posterior carries a matrix mean; the scalar checks below
  # would reject it, and the relative entropy it needs is the multivariate one.
  if (is.list(filt) && is.matrix(filt$mean)) {
    return(.aci_metric_mv(filt, smooth))
  }
  .aci_check_posterior(filt, "filt")
  .aci_check_posterior(smooth, "smooth", length(filt$mean))

  # ---- Relative entropy of the smoother from the filter ---------------------
  signal <- 0.5 * (smooth$mean - filt$mean)^2 / filt$cov
  ratio_delta <- smooth$cov / filt$cov - 1
  dispersion <- 0.5 * (ratio_delta - log1p(ratio_delta))
  .aci_metric_finish(signal + dispersion)
}

#' Clamp round-off negatives of raw metric values and count the clamps
#'
#' Applies the domain boundary of the causal-information metric to a raw
#' vector of relative-entropy values. A value that is exactly zero needs no
#' clamp and is not counted, which is what keeps the terminal identity out of
#' the count.
#'
#' @param value Numeric vector. Raw metric values, signal plus dispersion.
#'
#' @returns A list with the clamped numeric vector `value` and the integer
#'   count `n_clamped` of values that were negative before the clamp.
#'
#' @noRd
#' @keywords internal
.aci_metric_finish <- function(value) {
  # A relative entropy is finite and non-negative for a genuine posterior pair,
  # so anything else is either floating-point noise where the two posteriors
  # agree, or a defect. The noise is clamped within a documented tolerance; the
  # defect is reported rather than returned.
  #
  # Non-finiteness is tested first and separately. A covariance ratio that
  # overflows yields Inf - Inf = NaN, and NaN would slip past a comparison
  # guard unnoticed, since `NaN <= x` is NA, not TRUE.
  round_off <- 1e-10
  offending <- which(!is.finite(value) | value <= -round_off)
  if (length(offending) > 0L) {
    index <- offending[1L]
    stop(
      sprintf(
        paste0(
          "The causal-information metric is a relative entropy, so it must be ",
          "finite and non-negative, but it is %s at index %d. The posteriors ",
          "supplied are not a filter/smoother pair of one system: this arises ",
          "when their covariances differ by so many orders of magnitude that ",
          "their ratio overflows. If they are a genuine pair, please report ",
          "this as a bug."
        ),
        format(value[index]), index
      ),
      call. = FALSE
    )
  }
  n_clamped <- sum(value < 0)
  value[value < 0] <- 0
  list(value = value, n_clamped = n_clamped)
}
