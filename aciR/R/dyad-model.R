# -- nonlinear dyad model with intermittent extreme events ---------------------
#
# The dyad model pairs an observed process x with an unobserved process y whose
# feedback on x makes x intermittently extreme. It is the conditional Gaussian
# nonlinear system used as the worked example throughout the package: the two
# functions below build its components list and simulate a realisation.

#' Conditional Gaussian components of the nonlinear dyad model
#'
#' Builds the conditional Gaussian components list for the nonlinear dyad model
#' with intermittent extreme events, given an observed signal and the model
#' parameters. The returned list is the `comp` argument consumed by
#' [aci_filter()] and [aci_smoother()].
#'
#' In the dyad model the observed process couples to the unobserved process
#' through the state-dependent term `gamma * x`, and the unobserved process has
#' constant linear self-drift `-d_y`. The two noise sources are independent, so
#' the noise cross-covariances are zero.
#'
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param p A named list of dyad parameters with elements `d_x`, `d_y`,
#'   `gamma`, `F_x`, `F_y`, `sigma_x` and `sigma_y`.
#'
#' @return A conditional Gaussian components list; see [aci_components].
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
#' str(comp)
#'
#' @seealso [aci_simulate_dyad()], [aci_filter()]
#' @export
aci_dyad_components <- function(x, p) {
  list(
    L_x = p$gamma * x,
    f_x = p$F_x - p$d_x * x,
    L_y = -p$d_y,
    f_y = p$F_y - p$gamma * x^2,
    S_xoS_x = p$sigma_x^2,
    S_yoS_y = p$sigma_y^2,
    S_yoS_x = 0,
    S_xoS_y = 0
  )
}

#' Simulate the nonlinear dyad model
#'
#' Simulates a realisation of the nonlinear dyad model with intermittent
#' extreme events by an Euler-Maruyama integration of its pair of stochastic
#' differential equations. The observed process `x` responds to the unobserved
#' process `y` through a state-dependent coupling, which produces intermittent
#' bursts in `x`.
#'
#' The simulation draws standard normal increments and is therefore stochastic;
#' set a seed with [set.seed()] before calling for a reproducible path. The
#' initial conditions are the deterministic equilibria `F_x / d_x` and
#' `F_y / d_y`.
#'
#' @param n Integer scalar. The number of time steps to simulate, including the
#'   initial step.
#' @param p A named list of dyad parameters with elements `d_x`, `d_y`,
#'   `gamma`, `F_x`, `F_y`, `sigma_x` and `sigma_y`. The default is the
#'   parameter set of the reference model.
#' @param dt Numeric scalar. The integration time step.
#'
#' @return A list with three numeric vectors, `t` (time), `x` (the observed
#'   process) and `y` (the unobserved process), each of length `n`.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' set.seed(1)
#' sim <- aci_simulate_dyad(n = 2000)
#' str(sim)
#'
#' @seealso [aci_dyad_components()]
#' @export
aci_simulate_dyad <- function(n,
                              p = list(
                                d_x = 0.5, d_y = 0.5, gamma = 2,
                                F_x = 0.5, F_y = 1,
                                sigma_x = 0.5, sigma_y = 1
                              ),
                              dt = 0.001) {
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- p$F_x / p$d_x
  y[1] <- p$F_y / p$d_y
  for (j in 2:n) {
    dw_x <- rnorm(1)
    dw_y <- rnorm(1)
    l_x <- p$gamma * x[j - 1]
    f_x <- p$F_x - p$d_x * x[j - 1]
    f_y <- p$F_y - p$gamma * x[j - 1]^2
    x[j] <- x[j - 1] + (l_x * y[j - 1] + f_x) * dt +
      p$sigma_x * sqrt(dt) * dw_x
    y[j] <- y[j - 1] + (-p$d_y * y[j - 1] + f_y) * dt +
      p$sigma_y * sqrt(dt) * dw_y
  }
  list(t = seq(0, by = dt, length.out = n), x = x, y = y)
}
