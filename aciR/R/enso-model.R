# -- stochastic ENSO model -----------------------------------------------------
#
# The six-dimensional conditional Gaussian system of the paper's El Nino
# Southern Oscillation case study. Three observables -- central- and
# eastern-Pacific sea-surface temperature anomalies and the Walker circulation
# strength -- and three unobserved variables: the central-Pacific zonal
# current, the western-Pacific thermocline depth, and the intraseasonal
# westerly wind burst.
#
# It is the largest system the package expresses, and it exercises three things
# nothing else does. The coefficients are matrix-valued and vary at every step.
# They are seasonally modulated, on a six-unit period. And BOTH noise
# covariances vary in time: the Walker circulation's is multiplicative in its
# own state, and the latent one of the wind burst depends on the observed central-Pacific
# temperature through a hyperbolic tangent.
#
# That last property is why this constructor could not have been written
# earlier without shipping a defect: the vector validator inverted only the
# first slice of a time-varying observation-noise covariance until that was
# found and fixed, so every ENSO result would have been computed from a
# Walker-circulation noise frozen at its initial value.

#' Parameters of the stochastic ENSO model
#'
#' Builds the parameter list of the ENSO model, with the defaults of the
#' reference implementation. Most parameters are derived from a smaller set of
#' physical ones, so they are computed here rather than left for a caller to
#' keep consistent.
#'
#' @param factor Numeric scalar. The diversity-modulating parameter that
#'   governs how many extreme eastern-Pacific events the model produces.
#' @param b_0 Numeric scalar. High-end estimate of the thermocline tilt.
#' @param mu Numeric scalar. Relative coupling coefficient.
#' @param d_tau Numeric scalar. Damping of the wind burst.
#' @param lambda Numeric scalar. Damping of the Walker circulation.
#' @param m Numeric scalar. Target equilibrium mean of the Walker circulation.
#'
#' @returns A named list of the model's parameters, derived and primitive.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' str(aci_enso_parameters())
#'
#' @seealso [aci_enso_components()], [aci_enso_model()]
#' @export
aci_enso_parameters <- function(factor = 0.65, b_0 = 2.5, mu = 0.5,
                                d_tau = 2, lambda = 2 / 60, m = 2) {
  for (nm in c("factor", "b_0", "mu", "d_tau", "lambda", "m")) {
    .aci_check_scalar(get(nm), nm)
  }
  if (factor <= 0) {
    stop("`factor` must be positive.", call. = FALSE)
  }

  alpha_2 <- 0.125 * factor
  alpha_1 <- alpha_2 / 2 * factor
  gamma_C <- 0.75 * factor
  gamma_E <- 0.75 * factor

  list(
    factor = factor, b_0 = b_0, mu = mu, d_tau = d_tau,
    lambda = lambda, m = m,
    alpha_1 = alpha_1, alpha_2 = alpha_2,
    # Sea-surface-temperature feedbacks into the two ocean variables.
    delta_u = alpha_1 * b_0 * mu,
    delta_h = alpha_2 * b_0 * mu,
    # Collective damping of the ocean adjustment.
    r = 0.25 * factor,
    gamma_C = gamma_C, gamma_E = gamma_E,
    r_C = gamma_C * b_0 * mu / 2,
    r_E = 3 * gamma_E * b_0 * mu / 2,
    zeta_C = gamma_C * b_0 * mu / 2,
    zeta_E = gamma_E * b_0 * mu / 2,
    # Enforces zero climatology in the anomaly model.
    C_u = 0.03 * factor,
    sigma_u = 0.04 * sqrt(factor),
    sigma_h = 0.02 * sqrt(factor),
    sigma_C = 0.04 * sqrt(factor),
    # The eastern-Pacific temperature is noiseless in the original model. The
    # conditional construction needs an invertible observation-noise
    # covariance, so the reference adds the smallest noise that keeps the
    # filter and smoother stable. That is a modelling decision, carried here
    # rather than hidden.
    sigma_E = sqrt(5) * 1e-2 * sqrt(factor)
  )
}

#' Conditional Gaussian components of the stochastic ENSO model
#'
#' Builds the vector-valued components list of the ENSO model from the observed
#' paths, for the configuration in which the central- and eastern-Pacific
#' temperature anomalies and the Walker circulation are observed and the zonal
#' current, thermocline depth and wind burst are not.
#'
#' The observed state is `x = (T_C, T_E, I)` and the unobserved state is
#' `y = (u, h_W, tau)`, both three-dimensional. Every coefficient varies with
#' the observed state, the season, or both.
#'
#' Two noise covariances vary in time and neither is optional. The Walker
#' circulation's observation noise is multiplicative in its own state, with
#' variance `lambda * (4 - I) * I`, which vanishes at the ends of the interval
#' `[0, 4]` and so keeps the process inside it. The latent noise of the wind burst
#' depends on the observed central-Pacific temperature and on the season. A
#' filter that read either at its first step alone would be integrating a
#' different system.
#'
#' Because the Walker-circulation noise vanishes at `I = 0` and `I = 4`, an
#' observed path that reaches either endpoint gives a singular observation-noise
#' covariance. Such a path is rejected rather than regularised: the filter
#' inverts that covariance, and quietly nudging it would be inventing an
#' observation the data does not support.
#'
#' @param T_C,T_E,I Numeric vectors of equal length. The observed
#'   central-Pacific temperature anomaly, eastern-Pacific temperature anomaly,
#'   and Walker circulation strength.
#' @param p A parameter list, as returned by [aci_enso_parameters()].
#' @param time Numeric vector or `NULL`. The observation times, needed because
#'   the coefficients are seasonally modulated. When `NULL`, a regular grid of
#'   spacing `dt` starting at zero is used.
#' @param dt Numeric scalar. The step used to build the default `time` grid.
#'
#' @returns A vector-valued conditional Gaussian components list; see
#'   [aci_components].
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' p <- aci_enso_parameters()
#' n <- 200
#' tt <- seq(0, by = 0.005, length.out = n)
#' comp <- aci_enso_components(
#'   T_C = 0.05 * sin(tt), T_E = 0.04 * cos(tt),
#'   I = 1.6 + 0.2 * sin(0.5 * tt), p = p, time = tt
#' )
#' dim(comp$L_x)
#'
#' @seealso [aci_enso_parameters()], [aci_enso_model()], [aci_conditional()]
#' @export
aci_enso_components <- function(T_C, T_E, I, p = aci_enso_parameters(),
                                time = NULL, dt = 0.005) {
  for (nm in c("T_C", "T_E", "I")) {
    value <- get(nm)
    if (!is.numeric(value) || !is.null(dim(value)) || anyNA(value) ||
          any(!is.finite(value))) {
      stop(
        sprintf("`%s` must be a complete, finite numeric vector.", nm),
        call. = FALSE
      )
    }
  }
  n <- length(T_C)
  if (length(T_E) != n || length(I) != n) {
    stop(
      "`T_C`, `T_E` and `I` must have the same length.", call. = FALSE
    )
  }
  if (n < 2L) {
    stop("At least two observations are needed.", call. = FALSE)
  }
  .aci_check_positive(dt, "dt")
  if (is.null(time)) {
    time <- (seq_len(n) - 1L) * dt
  } else if (!is.numeric(time) || length(time) != n || anyNA(time)) {
    stop(
      "`time` must be a numeric vector as long as the observed paths.",
      call. = FALSE
    )
  }

  # The Walker circulation's observation noise vanishes at the ends of its
  # domain, and the filter inverts it.
  outside <- which(I <= 0 | I >= 4)
  if (length(outside) > 0L) {
    stop(
      sprintf(
        paste0(
          "The Walker circulation `I` must lie strictly inside (0, 4): its ",
          "observation noise has variance lambda * (4 - I) * I, which ",
          "vanishes at both ends and leaves the observation-noise covariance ",
          "singular. It is %g at index %d."
        ),
        I[outside[1L]], outside[1L]
      ),
      call. = FALSE
    )
  }

  # ---- Seasonal modulation --------------------------------------------------
  # A six-unit period; the model's time unit is two months, so this is the
  # annual cycle.
  season <- 2 * pi / 6
  c_1 <- (25 * (T_C + 0.75 / 7.5)^2 + 0.9) *
    (1 + 0.3 * sin(time * season - pi / 6)) * p$factor
  c_2 <- 1.4 * p$factor *
    (1 + 0.3 * sin(time * season + 2 * pi / 6) +
       0.25 * sin(2 * time * season + 2 * pi / 6))

  # ---- Wind-burst feedbacks, modulated by the Walker circulation ------------
  burst <- (1 + (1 - I / 5)) * 0.15 * sqrt(p$factor)
  beta_u <- -0.2 * burst
  beta_h <- -0.4 * burst
  beta_C <- 0.8 * burst
  beta_E <- 1 * burst

  L_x <- array(0, c(3L, 3L, n))
  L_y <- array(0, c(3L, 3L, n))
  S_xoS_x <- array(0, c(3L, 3L, n))
  S_yoS_y <- array(0, c(3L, 3L, n))
  f_x <- matrix(0, 3L, n)
  f_y <- matrix(0, 3L, n)

  # Coupling of the unobserved state into the observed drift. The third row is
  # zero: the Walker circulation is autonomous of the latent variables.
  L_x[1L, 1L, ] <- I / 5 * p$factor
  L_x[1L, 2L, ] <- p$gamma_C
  L_x[1L, 3L, ] <- beta_C
  L_x[2L, 2L, ] <- p$gamma_E
  L_x[2L, 3L, ] <- beta_E

  f_x[1L, ] <- (p$r_C - c_1) * T_C + p$zeta_C * T_E + p$C_u
  f_x[2L, ] <- (p$r_E - c_2) * T_E - p$zeta_E * T_C
  f_x[3L, ] <- -p$lambda * (I - p$m)

  # Self-drift of the unobserved state.
  L_y[1L, 1L, ] <- -p$r
  L_y[1L, 3L, ] <- beta_u
  L_y[2L, 2L, ] <- -p$r
  L_y[2L, 3L, ] <- beta_h
  L_y[3L, 3L, ] <- -p$d_tau

  f_y[1L, ] <- -p$delta_u * (T_C + T_E) / 2
  f_y[2L, ] <- -p$delta_h * (T_C + T_E) / 2

  # ---- Noise ----------------------------------------------------------------
  # Diagonal throughout: the six variables are driven by independent Wiener
  # processes, so the cross-covariance is zero, as in every model of the
  # reference implementation.
  S_xoS_x[1L, 1L, ] <- p$sigma_C^2
  S_xoS_x[2L, 2L, ] <- p$sigma_E^2
  S_xoS_x[3L, 3L, ] <- p$lambda * (4 - I) * I

  sigma_tau <- 0.9 * (tanh(7.5 * T_C) + 1) *
    (1 + 0.3 * cos(time * season + 2 * pi / 6))
  S_yoS_y[1L, 1L, ] <- p$sigma_u^2
  S_yoS_y[2L, 2L, ] <- p$sigma_h^2
  S_yoS_y[3L, 3L, ] <- sigma_tau^2

  list(
    L_x = L_x, f_x = f_x, L_y = L_y, f_y = f_y,
    S_xoS_x = S_xoS_x, S_yoS_y = S_yoS_y,
    S_yoS_x = matrix(0, 3L, 3L)
  )
}

#' Stochastic ENSO model
#'
#' Returns the ENSO model as a labelled description of its structure and
#' parameters. Unlike the scalar models, this one carries no simulator: the
#' six-dimensional system is integrated with a mixed scheme in the reference
#' implementation, and its wind-burst update contains a correction this package
#' does not reproduce. See the design note on that anomaly before simulating
#' the system by any route.
#'
#' Use [aci_enso_components()] with observed paths to obtain the components
#' consumed by [aci_filter()], [aci_smoother()] and [aci_metric()], and
#' [aci_conditional()] to ask the causal question of one observable given the
#' others.
#'
#' @param p A parameter list, as returned by [aci_enso_parameters()].
#'
#' @returns A list describing the model: its `label`, `parameters`, the names
#'   of its `observed` and `unobserved` components, and the `components`
#'   function that builds a components list from observed paths.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_enso_model()
#' model$observed
#' model$unobserved
#'
#' @seealso [aci_enso_components()], [aci_conditional()]
#' @export
aci_enso_model <- function(p = aci_enso_parameters()) {
  if (!is.list(p) || is.null(p$factor)) {
    stop(
      "`p` must be a parameter list from `aci_enso_parameters()`.",
      call. = FALSE
    )
  }
  list(
    label = "stochastic ENSO model with diversity",
    parameters = p,
    observed = c("T_C", "T_E", "I"),
    unobserved = c("u", "h_W", "tau"),
    components = function(T_C, T_E, I, time = NULL, dt = 0.005) {
      aci_enso_components(T_C, T_E, I, p = p, time = time, dt = dt)
    }
  )
}
