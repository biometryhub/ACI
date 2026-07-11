# -- conditional Gaussian nonlinear system: filter / smoother / metric ---------
#
# The three functions below are the numerical core of assimilative causal
# inference. They operate on a general CGNS "components" list rather than on a
# specific model, so any conditional Gaussian nonlinear system can be filtered,
# smoothed and scored by supplying the appropriate per-step coefficients.

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
#'   \item{`L_y`}{Numeric scalar. The linear self-drift of the unobserved
#'     component (constant in time for the supplied dyad model).}
#'   \item{`f_y`}{Numeric vector, length `n`. The remaining drift of the
#'     unobserved process at each time step.}
#'   \item{`S_xoS_x`}{Numeric scalar. The observation-noise covariance,
#'     the product of the observed-process noise coefficient with its
#'     transpose.}
#'   \item{`S_yoS_y`}{Numeric scalar. The latent-noise covariance of the
#'     unobserved process.}
#'   \item{`S_yoS_x`}{Numeric scalar. The latent-to-observation noise
#'     cross-covariance.}
#'   \item{`S_xoS_y`}{Numeric scalar. The observation-to-latent noise
#'     cross-covariance, the transpose of `S_yoS_x`.}
#' }
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
#' integrated with a first-order (Euler) step of width `dt`.
#'
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param comp A conditional Gaussian components list; see [aci_components].
#' @param dt Numeric scalar. The integration time step.
#' @param mu0 Numeric scalar. The initial filtered mean of the unobserved
#'   component.
#' @param R0 Numeric scalar. The initial filtered covariance of the unobserved
#'   component.
#'
#' @return A list with two numeric vectors, `mean` and `cov`, the filtered mean
#'   and covariance of the unobserved component at each time step.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' p <- list(
#'   d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
#'   sigma_x = 0.5, sigma_y = 1
#' )
#' set.seed(1)
#' sim <- aci_simulate_dyad(n = 2000, p = p)
#' comp <- aci_dyad_components(sim$x, p)
#' filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
#' str(filt)
#'
#' @seealso [aci_smoother()], [aci_metric()]
#' @export
aci_filter <- function(x, comp, dt, mu0, R0) {
  n <- length(x)
  m <- numeric(n)
  v <- numeric(n)
  m[1] <- mu0
  v[1] <- R0
  inv <- 1 / comp$S_xoS_x
  mu <- mu0
  R <- R0
  for (j in 2:n) {
    dx <- x[j] - x[j - 1]
    aux <- comp$S_yoS_x + R * comp$L_x[j - 1]
    mu <- mu + (comp$L_y * mu + comp$f_y[j - 1]) * dt +
      aux * inv * (dx - (comp$L_x[j - 1] * mu + comp$f_x[j - 1]) * dt)
    R <- R + (2 * comp$L_y * R + comp$S_yoS_y - aux * inv * aux) * dt
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
#' filtered trajectory returned by [aci_filter()].
#'
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param comp A conditional Gaussian components list; see [aci_components].
#' @param dt Numeric scalar. The integration time step.
#' @param filt A list with numeric vectors `mean` and `cov`, as returned by
#'   [aci_filter()].
#'
#' @return A list with two numeric vectors, `mean` and `cov`, the smoothed mean
#'   and covariance of the unobserved component at each time step.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' p <- list(
#'   d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
#'   sigma_x = 0.5, sigma_y = 1
#' )
#' set.seed(1)
#' sim <- aci_simulate_dyad(n = 2000, p = p)
#' comp <- aci_dyad_components(sim$x, p)
#' filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
#' smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
#' str(smooth)
#'
#' @seealso [aci_filter()], [aci_metric()]
#' @export
aci_smoother <- function(x, comp, dt, filt) {
  n <- length(x)
  m <- numeric(n)
  v <- numeric(n)
  m[n] <- filt$mean[n]
  v[n] <- filt$cov[n]
  inv <- 1 / comp$S_xoS_x
  muT <- m[n]
  RT <- v[n]
  for (j in (n - 1):1) {
    dx <- x[j + 1] - x[j]
    A_j <- comp$L_y - comp$S_yoS_x * inv * comp$L_x[j]
    B_j <- comp$S_yoS_y - comp$S_yoS_x * inv * comp$S_xoS_y
    muT <- muT - (comp$L_y * muT + comp$f_y[j] -
      B_j / filt$cov[j] * (filt$mean[j] - muT)) * dt +
      comp$S_yoS_x * inv * (-dx + (comp$L_x[j] * muT + comp$f_x[j]) * dt)
    RT <- RT - (2 * (A_j + B_j / filt$cov[j]) * RT - B_j) * dt
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
#' @param filt A list with numeric vectors `mean` and `cov`, the filtered mean
#'   and covariance, as returned by [aci_filter()].
#' @param smooth A list with numeric vectors `mean` and `cov`, the smoothed
#'   mean and covariance, as returned by [aci_smoother()].
#'
#' @return A numeric vector of the causal-information metric at each time step.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' p <- list(
#'   d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
#'   sigma_x = 0.5, sigma_y = 1
#' )
#' set.seed(1)
#' sim <- aci_simulate_dyad(n = 2000, p = p)
#' comp <- aci_dyad_components(sim$x, p)
#' filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
#' smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
#' aci <- aci_metric(filt, smooth)
#' summary(aci)
#'
#' @seealso [aci_filter()], [aci_smoother()]
#' @export
aci_metric <- function(filt, smooth) {
  signal <- 0.5 * (smooth$mean - filt$mean)^2 / filt$cov
  cov_ratio <- smooth$cov / filt$cov
  dispersion <- 0.5 * (-log(cov_ratio) + cov_ratio - 1)
  signal + dispersion
}
