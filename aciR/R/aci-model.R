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

# -- internal helpers ---------------------------------------------------------

.aci_check_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(
      sprintf("`%s` must be a single finite numeric value.", name),
      call. = FALSE
    )
  }
  invisible(value)
}

.aci_as_coef <- function(value, name) {
  if (is.function(value)) {
    return(value)
  }
  if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
    force(value)
    return(function(x) rep_len(value, length(x)))
  }
  stop(
    sprintf(
      "`%s` must be a vectorised function of the observed signal or a %s",
      name, "finite numeric scalar."
    ),
    call. = FALSE
  )
}

# -- constructors -------------------------------------------------------------

#' Conditional Gaussian nonlinear system model
#'
#' Builds an `aci_model` object for a conditional Gaussian nonlinear system
#' (CGNS). A CGNS pairs an observed process `x` with an unobserved process `y`
#' whose statistics, conditional on the observed path, are Gaussian in closed
#' form. The model is described by the coefficients of the pair of governing
#' stochastic differential equations, and this constructor is the elegant
#' general entry point for supplying them.
#'
#' The observed process obeys `dx = (L_x(x) y + f_x(x)) dt + dW_x` and the
#' unobserved process obeys `dy = (L_y y + f_y(x)) dt + dW_y`, where the drift
#' of the unobserved component is linear in `y`. The coefficients `L_x`, `f_x`
#' and `f_y` may vary with the observed signal and are supplied as vectorised
#' functions of `x` (a numeric constant is accepted and treated as constant in
#' time). The self-drift `L_y` of the unobserved component is a numeric scalar:
#' the numerical core in this package integrates a time-invariant `L_y`, which
#' covers the reference dyad model and the general class it belongs to. Systems
#' with a time-varying `L_y`, such as the noisy predator-prey model, are
#' roadmapped; see the package `NEWS`.
#'
#' The noise is described by its Grammians, the entries of the noise covariance
#' of the pair of processes. `S_xoS_x` and `S_yoS_y` are the observation-noise
#' and latent-noise covariances, and `S_yoS_x` (with its transpose `S_xoS_y`)
#' is the noise cross-covariance, zero for independent noise.
#'
#' @param L_x Coupling of the unobserved component into the drift of the
#'   observed process. A vectorised function of the observed signal, or a
#'   numeric scalar for a constant coupling.
#' @param f_x Remaining drift of the observed process. A vectorised function of
#'   the observed signal, or a numeric scalar.
#' @param L_y Numeric scalar. The linear self-drift of the unobserved
#'   component, constant in time.
#' @param f_y Remaining drift of the unobserved process. A vectorised function
#'   of the observed signal, or a numeric scalar.
#' @param S_xoS_x Numeric scalar. The observation-noise covariance; must be
#'   positive.
#' @param S_yoS_y Numeric scalar. The latent-noise covariance of the unobserved
#'   process.
#' @param S_yoS_x Numeric scalar. The latent-to-observation noise
#'   cross-covariance. The default `0` is independent noise.
#' @param S_xoS_y Numeric scalar. The observation-to-latent noise
#'   cross-covariance, the transpose of `S_yoS_x`. Defaults to `S_yoS_x`.
#' @param x0 Numeric scalar. The initial value of the observed process, used
#'   when simulating.
#' @param y0 Numeric scalar. The initial value of the unobserved process, used
#'   when simulating and as the default initial filtered mean.
#' @param label Character string. A human-readable name for the model.
#' @param parameters Optional named list of the scalar parameters that define
#'   the model, retained for printing and reproducibility.
#'
#' @return An `aci_model` object: a list of the model coefficients, noise
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
#' summary(fit$aci)
#'
#' @seealso [aci_dyad_model()], [aci_simulate()], [aci()]
#' @export
aci_cgns_model <- function(L_x, f_x, L_y, f_y,
                           S_xoS_x, S_yoS_y,
                           S_yoS_x = 0, S_xoS_y = S_yoS_x,
                           x0 = 0, y0 = 0,
                           label = "conditional Gaussian nonlinear system",
                           parameters = NULL) {
  .aci_check_scalar(L_y, "L_y")
  .aci_check_scalar(S_xoS_x, "S_xoS_x")
  .aci_check_scalar(S_yoS_y, "S_yoS_y")
  .aci_check_scalar(S_yoS_x, "S_yoS_x")
  .aci_check_scalar(S_xoS_y, "S_xoS_y")
  .aci_check_scalar(x0, "x0")
  .aci_check_scalar(y0, "y0")
  if (S_xoS_x <= 0) {
    stop(
      "`S_xoS_x` must be positive; it is the observation-noise covariance.",
      call. = FALSE
    )
  }
  model <- list(
    L_x = .aci_as_coef(L_x, "L_x"),
    f_x = .aci_as_coef(f_x, "f_x"),
    L_y = L_y,
    f_y = .aci_as_coef(f_y, "f_y"),
    S_xoS_x = S_xoS_x,
    S_yoS_y = S_yoS_y,
    S_yoS_x = S_yoS_x,
    S_xoS_y = S_xoS_y,
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
#' @param sigma_x Numeric scalar. Noise coefficient of the observed process.
#' @param sigma_y Numeric scalar. Noise coefficient of the unobserved process.
#'
#' @return An `aci_model` object; see [aci_cgns_model()].
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
#' summary(fit$aci)
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
#' itself, not on a fresh simulation. Set `seed` for a path that is reproducible
#' within R.
#'
#' Simulation supports independent noise, that is a model with zero noise
#' cross-covariance.
#'
#' @param model An `aci_model` object; see [aci_cgns_model()].
#' @param n Integer scalar. The number of time steps to simulate, including the
#'   initial step; at least two.
#' @param dt Numeric scalar. The integration time step.
#' @param seed Optional integer. If supplied, passed to [set.seed()] for a
#'   reproducible path.
#'
#' @return A data frame with `n` rows and columns `t` (time), `x` (the observed
#'   process) and `y` (the unobserved process).
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
#' @seealso [aci_dyad_model()], [aci()]
#' @export
aci_simulate <- function(model, n, dt = 0.001, seed = NULL) {
  if (!inherits(model, "aci_model")) {
    stop("`model` must be an `aci_model`; see `aci_cgns_model()`.",
      call. = FALSE
    )
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
    n < 2 || n != round(n)) {
    stop("`n` must be a single integer of at least two.", call. = FALSE)
  }
  .aci_check_scalar(dt, "dt")
  if (dt <= 0) {
    stop("`dt` must be positive.", call. = FALSE)
  }
  if (model$S_yoS_x != 0 || model$S_xoS_y != 0) {
    stop(
      "`aci_simulate()` supports independent noise only; the model has a ",
      "non-zero noise cross-covariance.",
      call. = FALSE
    )
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- as.integer(n)
  sigma_x <- sqrt(model$S_xoS_x)
  sigma_y <- sqrt(model$S_yoS_y)
  root_dt <- sqrt(dt)
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- model$x0
  y[1] <- model$y0
  for (j in 2:n) {
    dw_x <- rnorm(1)
    dw_y <- rnorm(1)
    l_x <- model$L_x(x[j - 1])
    f_x <- model$f_x(x[j - 1])
    f_y <- model$f_y(x[j - 1])
    x[j] <- x[j - 1] + (l_x * y[j - 1] + f_x) * dt + sigma_x * root_dt * dw_x
    y[j] <- y[j - 1] + (model$L_y * y[j - 1] + f_y) * dt +
      sigma_y * root_dt * dw_y
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
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param model An `aci_model` object; see [aci_cgns_model()].
#' @param dt Numeric scalar. The integration time step.
#' @param mu0 Numeric scalar. The initial filtered mean of the unobserved
#'   component. Defaults to the model's initial unobserved state `y0`.
#' @param R0 Numeric scalar. The initial filtered covariance of the unobserved
#'   component; must be positive.
#'
#' @return An `aci` object: a list with the `model`, the time vector `t`, the
#'   observed signal `x`, the `filter` and `smoother` statistics (each a list of
#'   `mean` and `cov`), the causal-information metric `aci`, and the step `dt`.
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
#' str(fit$filter)
#'
#' @seealso [aci_filter()], [aci_smoother()], [aci_metric()]
#' @export
aci <- function(x, model, dt = 0.001, mu0 = NULL, R0 = 0.1) {
  if (!is.numeric(x) || length(x) < 2L) {
    stop("`x` must be a numeric vector of at least two observations.",
      call. = FALSE
    )
  }
  if (!inherits(model, "aci_model")) {
    stop("`model` must be an `aci_model`; see `aci_cgns_model()`.",
      call. = FALSE
    )
  }
  .aci_check_scalar(dt, "dt")
  if (dt <= 0) {
    stop("`dt` must be positive.", call. = FALSE)
  }
  if (is.null(mu0)) {
    mu0 <- model$y0
  }
  .aci_check_scalar(mu0, "mu0")
  .aci_check_scalar(R0, "R0")
  if (R0 <= 0) {
    stop("`R0` must be positive.", call. = FALSE)
  }
  comp <- list(
    L_x = model$L_x(x),
    f_x = model$f_x(x),
    L_y = model$L_y,
    f_y = model$f_y(x),
    S_xoS_x = model$S_xoS_x,
    S_yoS_y = model$S_yoS_y,
    S_yoS_x = model$S_yoS_x,
    S_xoS_y = model$S_xoS_y
  )
  filt <- aci_filter(x, comp, dt, mu0, R0)
  smooth <- aci_smoother(x, comp, dt, filt)
  metric <- aci_metric(filt, smooth)
  result <- list(
    model = model,
    t = seq(0, by = dt, length.out = length(x)),
    x = x,
    filter = filt,
    smoother = smooth,
    aci = metric,
    dt = dt
  )
  class(result) <- "aci"
  result
}

# -- print methods ------------------------------------------------------------

#' Print a conditional Gaussian nonlinear system model
#'
#' Prints a compact, one-field-per-line summary of an `aci_model` object: its
#' label, the self-drift of the unobserved component, the noise Grammians, the
#' initial state and, when present, the named scalar parameters that define the
#' model.
#'
#' @param x An `aci_model` object; see [aci_cgns_model()].
#' @param ... Ignored, for compatibility with [print()].
#'
#' @return The model `x`, invisibly.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' print(aci_dyad_model())
#'
#' @seealso [aci_cgns_model()], [aci_dyad_model()]
#' @export
print.aci_model <- function(x, ...) {
  cat("<aci_model> ", x$label, "\n", sep = "")
  cat(sprintf("  unobserved self-drift L_y: %g\n", x$L_y))
  cat(sprintf(
    "  noise Grammians: S_xoS_x = %g, S_yoS_y = %g, S_yoS_x = %g\n",
    x$S_xoS_x, x$S_yoS_y, x$S_yoS_x
  ))
  cat(sprintf("  initial state: x0 = %g, y0 = %g\n", x$x0, x$y0))
  if (!is.null(x$parameters)) {
    pars <- paste(names(x$parameters), unlist(x$parameters),
      sep = " = ", collapse = ", "
    )
    cat("  parameters: ", pars, "\n", sep = "")
  }
  invisible(x)
}

#' Print an assimilative causal inference result
#'
#' Prints a compact summary of an `aci` object: the model it was run for, the
#' number of steps, the integration step and time span, and a five-number
#' summary of the causal-information metric.
#'
#' @param x An `aci` object, as returned by [aci()].
#' @param ... Ignored, for compatibility with [print()].
#'
#' @return The result `x`, invisibly.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' fit <- aci(aci_simulate(model, n = 2000, seed = 1)$x, model)
#' print(fit)
#'
#' @seealso [aci()]
#' @export
print.aci <- function(x, ...) {
  n <- length(x$x)
  cat("<aci> assimilative causal inference\n")
  cat("  model: ", x$model$label, "\n", sep = "")
  cat(sprintf(
    "  steps: %d, dt: %g, time span: [0, %g]\n", n, x$dt, x$t[n]
  ))
  cat("  causal-information metric:\n")
  print(summary(x$aci))
  invisible(x)
}
