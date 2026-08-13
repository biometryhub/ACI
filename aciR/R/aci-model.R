# -- conditional Gaussian nonlinear system as a model object ------------------
#
# The functions below wrap a conditional Gaussian nonlinear system (CGNS) in a
# small `aci_model` object and provide the two high-level entry points that a
# user reaches for first: simulate a realisation of a model, and run the whole
# assimilative causal inference workflow on an observed signal. The object holds
# the model coefficients as functions of the observed signal, and the entry
# points assemble the per-step components consumed by the numerical core in
# aci-core.R. The core is unchanged; these functions are a convenience layer
# over it.

# -- constructors -------------------------------------------------------------

#' Conditional Gaussian nonlinear system model
#'
#' Builds an `aci_model` object for a conditional Gaussian nonlinear system
#' (CGNS). A CGNS pairs an observed process `x` with an unobserved process `y`
#' whose statistics, conditional on the observed path, are Gaussian in closed
#' form. The model is described by the coefficients of the pair of governing
#' stochastic differential equations, and this constructor is the general entry
#' point for supplying them.
#'
#' The observed process obeys `dx = (L_x(x) y + f_x(x)) dt + dW_x` and the
#' unobserved process obeys `dy = (L_y(x) y + f_y(x)) dt + dW_y`, where the
#' drift of the unobserved component is linear in `y`. All four coefficients may
#' vary with the observed signal and are supplied as vectorised functions of
#' `x`; a numeric constant is accepted and treated as constant in time.
#'
#' A self-drift `L_y` that varies with the observed signal keeps the system
#' conditionally Gaussian, because the coefficient is measurable with respect to
#' the observed path. This is what the noisy predator-prey model needs, where
#' the latent population's growth rate is set by the population being watched.
#' What would leave the class is a coefficient depending on the *unobserved*
#' component, and no such system can be built here.
#'
#' The noise is described by its Grammians, the entries of the noise covariance
#' of the pair of processes. `S_xoS_x` and `S_yoS_y` are the observation-noise
#' and latent-noise covariances, and `S_yoS_x` is the noise cross-covariance,
#' zero for independent noise. For the scalar systems this package integrates
#' the cross-covariance is symmetric, so the transpose `S_xoS_y` of the
#' components schema is derived rather than supplied.
#'
#' The noise covariance must be mathematically admissible: `S_xoS_x` strictly
#' positive (the filter inverts it), `S_yoS_y` non-negative, and the joint
#' covariance positive semidefinite, that is
#' `S_xoS_x * S_yoS_y - S_yoS_x^2 >= 0`. A model that violates any of these
#' cannot be constructed. A singular system, whose determinant is zero,
#' describes perfectly correlated noise and is admissible, but it can drive the
#' filtered covariance toward zero; [aci_filter()] reports the step at which
#' that happens rather than returning an uninterpretable result.
#'
#' @param L_x Coupling of the unobserved component into the drift of the
#'   observed process. A vectorised function of the observed signal returning
#'   one value per observation, or a numeric scalar for a constant coupling.
#' @param f_x Remaining drift of the observed process. A vectorised function of
#'   the observed signal, or a numeric scalar.
#' @param L_y Linear self-drift of the unobserved component. A vectorised
#'   function of the observed signal, or a numeric scalar for a self-drift
#'   constant in time.
#' @param f_y Remaining drift of the unobserved process. A vectorised function
#'   of the observed signal, or a numeric scalar.
#' @param S_xoS_x Numeric scalar. The observation-noise covariance; must be
#'   positive.
#' @param S_yoS_y Numeric scalar. The latent-noise covariance of the unobserved
#'   process; must be non-negative.
#' @param S_yoS_x Numeric scalar. The latent-to-observation noise
#'   cross-covariance. The default `0` is independent noise.
#' @param x0 Numeric scalar. The initial value of the observed process, used
#'   when simulating.
#' @param y0 Numeric scalar. The initial value of the unobserved process, used
#'   when simulating and as the default initial filtered mean.
#' @param label Character string. A human-readable name for the model.
#' @param parameters Optional named list of the finite numeric scalars that
#'   define the model, retained for printing and reproducibility.
#'
#' @returns An `aci_model` object: a list of the model coefficients, noise
#'   Grammians and initial state, with class `"aci_model"`.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' # A linear Ornstein-Uhlenbeck pair with constant coupling.
#' model <- aci_cgns_model(
#'   L_x = 1, f_x = function(x) -0.5 * x, L_y = -0.5, f_y = 0,
#'   S_xoS_x = 0.25, S_yoS_y = 1, y0 = 0
#' )
#' sim <- aci_simulate(model, n = 2000, dt = 0.01, seed = 1)
#' fit <- aci(sim$x, model, dt = 0.01)
#' summary(fit)
#'
#' # An inadmissible noise covariance is rejected at construction.
#' try(aci_cgns_model(
#'   L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
#'   S_xoS_x = 1, S_yoS_y = -1
#' ))
#'
#' @seealso [aci_dyad_model()], [aci_simulate()], [aci()], [aci_components]
#' @export
aci_cgns_model <- function(L_x, f_x, L_y, f_y,
                           S_xoS_x, S_yoS_y,
                           S_yoS_x = 0,
                           x0 = 0, y0 = 0,
                           label = "conditional Gaussian nonlinear system",
                           parameters = NULL) {
  # The unobserved component's self-drift is admitted on the same terms as the
  # other coefficients: a constant, or a vectorised function of the observed
  # signal. A system whose latent damping is set by the observed state -- a
  # predator whose growth rate depends on the prey it can see -- is still
  # conditionally Gaussian, because the coefficient is measurable with respect
  # to the observed path.
  l_y_is_constant <- !is.function(L_y)
  .aci_check_noise_covariance(S_xoS_x, S_yoS_y, S_yoS_x)
  .aci_check_scalar(x0, "x0")
  .aci_check_scalar(y0, "y0")
  .aci_check_label(label)
  .aci_check_parameters(parameters)

  model <- list(
    L_x = .aci_as_coef(L_x, "L_x"),
    f_x = .aci_as_coef(f_x, "f_x"),
    L_y = .aci_as_coef(L_y, "L_y"),
    # Retained so the printed form can show the value when it is one, rather
    # than reporting every model's self-drift as state-dependent.
    L_y_constant = if (l_y_is_constant) L_y else NA_real_,
    f_y = .aci_as_coef(f_y, "f_y"),
    S_xoS_x = S_xoS_x,
    S_yoS_y = S_yoS_y,
    S_yoS_x = S_yoS_x,
    # The cross-covariance of a scalar system is its own transpose. The entry
    # is retained so a model widens directly into the components schema.
    S_xoS_y = S_yoS_x,
    x0 = x0,
    y0 = y0,
    label = label,
    parameters = parameters
  )
  class(model) <- "aci_model"
  model
}

#' Nonlinear dyad model
#'
#' Builds the `aci_model` object for the nonlinear dyad model with intermittent
#' extreme events, the flagship conditional Gaussian nonlinear system of the
#' package. The observed process `x` couples to the unobserved process `y`
#' through the state-dependent term `gamma x`, which makes `x` intermittently
#' extreme, and the unobserved process has constant self-drift `-d_y`.
#'
#' The two governing equations are
#' `dx = (- d_x x + gamma x y + F_x) dt + sigma_x dW_x` and
#' `dy = (- d_y y - gamma x^2 + F_y) dt + sigma_y dW_y`, with independent noise.
#' The default parameters are those of the reference implementation.
#'
#' @param d_x Numeric scalar. Damping of the observed process.
#' @param d_y Numeric scalar. Damping of the unobserved process.
#' @param gamma Numeric scalar. Strength of the quadratic coupling.
#' @param F_x Numeric scalar. Constant forcing of the observed process.
#' @param F_y Numeric scalar. Constant forcing of the unobserved process.
#' @param sigma_x Numeric scalar. Noise coefficient of the observed process;
#'   must be non-zero, since its square is the observation-noise covariance.
#' @param sigma_y Numeric scalar. Noise coefficient of the unobserved process.
#'
#' @returns An `aci_model` object; see [aci_cgns_model()].
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' model
#' sim <- aci_simulate(model, n = 5000, seed = 333)
#' fit <- aci(sim$x, model)
#' summary(fit)
#'
#' @seealso [aci_cgns_model()], [aci_simulate()], [aci()]
#' @export
aci_dyad_model <- function(d_x = 0.5, d_y = 0.5, gamma = 2,
                           F_x = 0.5, F_y = 1,
                           sigma_x = 0.5, sigma_y = 1) {
  .aci_check_scalar(d_x, "d_x")
  .aci_check_scalar(d_y, "d_y")
  .aci_check_scalar(gamma, "gamma")
  .aci_check_scalar(F_x, "F_x")
  .aci_check_scalar(F_y, "F_y")
  .aci_check_scalar(sigma_x, "sigma_x")
  .aci_check_scalar(sigma_y, "sigma_y")
  if (d_x == 0 || d_y == 0) {
    stop(
      "`d_x` and `d_y` must be non-zero; the initial state is their quotient.",
      call. = FALSE
    )
  }

  aci_cgns_model(
    L_x = function(x) gamma * x,
    f_x = function(x) F_x - d_x * x,
    L_y = -d_y,
    f_y = function(x) F_y - gamma * x^2,
    S_xoS_x = sigma_x^2,
    S_yoS_y = sigma_y^2,
    S_yoS_x = 0,
    x0 = F_x / d_x,
    y0 = F_y / d_y,
    label = "nonlinear dyad model with intermittent extreme events",
    parameters = list(
      d_x = d_x, d_y = d_y, gamma = gamma, F_x = F_x, F_y = F_y,
      sigma_x = sigma_x, sigma_y = sigma_y
    )
  )
}

# -- simulation ---------------------------------------------------------------

#' Simulate a conditional Gaussian nonlinear system
#'
#' Simulates a realisation of the observed and unobserved processes of an
#' `aci_model` by an Euler-Maruyama integration of its pair of stochastic
#' differential equations. The simulation starts from the model's initial state
#' and draws independent standard normal increments at each step.
#'
#' The random increments come from R's generator, which differs from that of
#' the reference MATLAB implementation, so a simulated path here does not
#' reproduce the reference path even at a matching seed. This is expected: the
#' independent-oracle grade of the numerical core is run on the MATLAB signal
#' itself, not on a fresh simulation.
#'
#' Supplying `seed` makes the path reproducible without disturbing the calling
#' session: the generator state is saved before the draw and restored when the
#' function exits, so a seeded call has no effect on the sequence a caller
#' would otherwise have seen. Leaving `seed` as `NULL` consumes the global
#' stream in the ordinary way.
#'
#' Simulation supports independent noise, that is a model with zero noise
#' cross-covariance.
#'
#' @param model An `aci_model` object; see [aci_cgns_model()].
#' @param n Integer scalar. The number of time steps to simulate, including the
#'   initial step; at least two.
#' @param dt Numeric scalar. The integration time step; must be positive.
#' @param seed Optional integer. If supplied, seeds a reproducible path and
#'   restores the caller's generator state on exit.
#' @param scheme Character scalar. The integration scheme,
#'   `"euler_maruyama"` (the default) or `"milstein"`. The two coincide unless
#'   the diffusion varies with the state it integrates, so `"milstein"`
#'   requires `sigma_x` and `d_sigma_x`.
#' @param sigma_x Optional function of the observed state. Supplying it makes
#'   the observation noise multiplicative, overriding the model's constant
#'   `S_xoS_x` **for the simulation only**. This is a simulation capability:
#'   the filter still requires a constant observation-noise covariance, so a
#'   path generated this way cannot yet be assimilated by this package.
#' @param d_sigma_x Optional function of the observed state, the derivative of
#'   `sigma_x`. Required by the Milstein scheme, and never estimated
#'   numerically on the caller's behalf.
#'
#' @returns A data frame with `n` rows and columns `t` (time), `x` (the
#'   observed process) and `y` (the unobserved process).
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 5000, seed = 333)
#' head(sim)
#'
#' # A seeded call leaves the caller's random stream untouched.
#' set.seed(1)
#' before <- runif(1)
#' set.seed(1)
#' invisible(aci_simulate(model, n = 10, seed = 99))
#' identical(before, runif(1))
#'
#' @seealso [aci_dyad_model()], [aci()]
#' @export
aci_simulate <- function(model, n, dt = 0.001, seed = NULL,
                         scheme = c("euler_maruyama", "milstein"),
                         sigma_x = NULL, d_sigma_x = NULL) {
  if (!inherits(model, "aci_model")) {
    stop("`model` must be an `aci_model`; see `aci_cgns_model()`.",
      call. = FALSE
    )
  }
  scheme <- match.arg(scheme)
  whole_number <- is.numeric(n) && length(n) == 1L && is.finite(n) &&
    n == round(n)
  if (!whole_number || n < 2) {
    stop("`n` must be a single integer of at least two.", call. = FALSE)
  }
  .aci_check_positive(dt, "dt")
  if (model$S_yoS_x != 0 || model$S_xoS_y != 0) {
    stop(
      "`aci_simulate()` supports independent noise only; the model has a ",
      "non-zero noise cross-covariance.",
      call. = FALSE
    )
  }
  if (!is.null(seed)) {
    .aci_check_scalar(seed, "seed")
    # Seeding is contained: the caller's generator state is restored on exit,
    # so a reproducible path costs the caller nothing.
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(
      assign(".Random.seed", old_seed, envir = globalenv()),
      add = TRUE
    )
    set.seed(seed)
  }

  # ---- Noise coefficients ---------------------------------------------------
  #
  # The observation noise may be supplied as a function of the observed state,
  # which makes it multiplicative. That is a SIMULATION capability: the filter
  # still requires a constant observation-noise covariance, so a path
  # simulated this way cannot yet be assimilated by this package.
  #
  # A conditional Gaussian system's noise may depend on the observed state but
  # never on the unobserved one. The Milstein correction is therefore non-zero
  # for the observed process and identically zero for the unobserved one,
  # whose diffusion does not vary with the variable it integrates.
  state_dependent <- !is.null(sigma_x)
  if (state_dependent && !is.function(sigma_x)) {
    stop("`sigma_x` must be a function of the observed state.", call. = FALSE)
  }
  if (identical(scheme, "milstein")) {
    if (!state_dependent) {
      stop(
        paste0(
          "The Milstein scheme differs from Euler-Maruyama only when the ",
          "diffusion varies with the state it integrates. Supply `sigma_x` ",
          "as a function of the observed state, or leave the scheme at ",
          "\"euler_maruyama\"."
        ),
        call. = FALSE
      )
    }
    if (!is.function(d_sigma_x)) {
      stop(
        paste0(
          "`d_sigma_x` must be supplied as a function for the Milstein ",
          "scheme: the correction term is the derivative of the diffusion ",
          "with respect to the state, and this package will not estimate it ",
          "numerically on your behalf."
        ),
        call. = FALSE
      )
    }
  }

  # ---- Integration ----------------------------------------------------------
  n <- as.integer(n)
  sigma_x_const <- sqrt(model$S_xoS_x)
  sigma_y <- sqrt(model$S_yoS_y)
  root_dt <- sqrt(dt)
  x <- numeric(n)
  y <- numeric(n)
  x[1L] <- model$x0
  y[1L] <- model$y0
  # The increments are exogenous, so they are drawn once rather than two at a
  # time inside the recursion.
  dw_x <- stats::rnorm(n - 1L)
  dw_y <- stats::rnorm(n - 1L)
  for (j in seq_len(n - 1L) + 1L) {
    l_x <- model$L_x(x[j - 1L])
    f_x <- model$f_x(x[j - 1L])
    f_y <- model$f_y(x[j - 1L])
    s_x <- if (state_dependent) sigma_x(x[j - 1L]) else sigma_x_const
    x[j] <- x[j - 1L] + (l_x * y[j - 1L] + f_x) * dt +
      s_x * root_dt * dw_x[j - 1L]
    if (identical(scheme, "milstein")) {
      # The Milstein correction: half the diffusion times its own derivative,
      # against the centred square of the increment. With dW = sqrt(dt) * z,
      # the bracket is dt * z^2 - dt.
      x[j] <- x[j] + 0.5 * s_x * d_sigma_x(x[j - 1L]) *
        (dt * dw_x[j - 1L]^2 - dt)
    }
    l_y <- model$L_y(x[j - 1L])
    y[j] <- y[j - 1L] + (l_y * y[j - 1L] + f_y) * dt +
      sigma_y * root_dt * dw_y[j - 1L]
  }
  data.frame(t = seq(0, by = dt, length.out = n), x = x, y = y)
}

# -- inference entry point ----------------------------------------------------

#' Run assimilative causal inference
#'
#' Runs the whole assimilative causal inference workflow on an observed signal
#' for a given `aci_model`. It assembles the conditional Gaussian components
#' from the model and the signal, runs the forward filter and backward smoother,
#' and scores each step with the causal-information metric. This is the single
#' user-facing entry point that ties the numerical core together.
#'
#' The workflow assumes the observed signal is complete and sampled on a
#' regular grid: the closed-form recursions integrate a fixed step, and there
#' is no contract for missing observations or for an irregular grid. By default
#' the time vector is constructed from `dt`; supply `time` instead to have the
#' step derived from an observed grid and the regularity checked.
#'
#' The metric this returns is a statement about the supplied model, not about
#' the world. See the *Assumptions and interpretation* vignette for what an ACI
#' peak does and does not support.
#'
#' @param x Numeric vector. The observed signal, one complete finite value per
#'   time step; at least two.
#' @param model An `aci_model` object; see [aci_cgns_model()].
#' @param dt Numeric scalar. The integration time step; must be positive.
#'   Ignored with a check when `time` is supplied.
#' @param mu0 Numeric scalar. The initial filtered mean of the unobserved
#'   component. Defaults to the model's initial unobserved state `y0`.
#' @param R0 Numeric scalar. The initial filtered covariance of the unobserved
#'   component; must be positive.
#' @param time Optional numeric vector. The observed time grid, one value per
#'   observation, strictly increasing and equally spaced. When supplied, `dt`
#'   is derived from it.
#'
#' @returns An `aci` object: a list with the `model`, the time vector `t`, the
#'   observed signal `x`, the `filter` and `smoother` statistics (each a list of
#'   `mean` and `cov`), the causal-information metric `aci`, the count
#'   `n_clamped` of metric values clamped from a round-off negative to zero,
#'   and the step `dt`. See [summary.aci()], [plot.aci()] and
#'   [as.data.frame.aci()].
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 5000, seed = 333)
#' fit <- aci(sim$x, model)
#' fit
#'
#' # An observed time grid may be supplied instead of a step.
#' fit_time <- aci(sim$x, model, time = sim$t)
#' identical(fit$aci, fit_time$aci)
#'
#' @seealso [aci_filter()], [aci_smoother()], [aci_metric()], [summary.aci()]
#' @export
aci <- function(x, model, dt = 0.001, mu0 = NULL, R0 = 0.1, time = NULL) {
  .aci_check_signal(x)
  if (!inherits(model, "aci_model")) {
    stop("`model` must be an `aci_model`; see `aci_cgns_model()`.",
      call. = FALSE
    )
  }
  .aci_check_positive(dt, "dt")
  if (!is.null(time)) {
    derived <- .aci_check_time(time, length(x))
    if (!missing(dt) && abs(derived - dt) > .aci_tol(derived, dt)) {
      stop(
        sprintf(
          paste0(
            "`dt` and `time` disagree: `time` is spaced %g apart but `dt` is ",
            "%g. Supply one or the other."
          ),
          derived, dt
        ),
        call. = FALSE
      )
    }
    dt <- derived
  }
  if (is.null(mu0)) {
    mu0 <- model$y0
  }
  .aci_check_scalar(mu0, "mu0")
  .aci_check_positive(R0, "R0")

  # ---- Components, then the core --------------------------------------------
  comp <- list(
    L_x = .aci_eval_coef(model$L_x, "L_x", x),
    f_x = .aci_eval_coef(model$f_x, "f_x", x),
    L_y = .aci_eval_coef(model$L_y, "L_y", x),
    f_y = .aci_eval_coef(model$f_y, "f_y", x),
    S_xoS_x = model$S_xoS_x,
    S_yoS_y = model$S_yoS_y,
    S_yoS_x = model$S_yoS_x,
    S_xoS_y = model$S_xoS_y
  )
  filt <- aci_filter(x, comp, dt, mu0, R0)
  smooth <- aci_smoother(x, comp, dt, filt)
  metric <- .aci_metric_pair(filt, smooth)

  result <- list(
    model = model,
    t = if (is.null(time)) seq(0, by = dt, length.out = length(x)) else time,
    x = x,
    filter = filt,
    smoother = smooth,
    aci = metric$value,
    n_clamped = metric$n_clamped,
    dt = dt
  )
  class(result) <- "aci"
  result
}
