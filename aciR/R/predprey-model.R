# Noisy predator-prey model ----------------------------------------------------
#
# A stochastic Lotka-Volterra pair. Unlike the dyad, neither variable is
# privileged. The reference implementation studies the system twice, once with
# the prey observed and the predator latent, and once the other way round. Both
# directions are supported here, because the causal question the method asks is
# directional and the two answers are not each other's mirror image.
#
# This is the model that requires a time-varying self-drift. In both directions
# the latent variable's own damping depends on the OBSERVED state. The
# predator's growth rate is set by how much prey it can see, and the prey's by
# how many predators. That is still a conditional Gaussian system, since the
# coefficient is measurable with respect to the observed path, but it is not a
# constant.

#' Validate a predator-prey parameter list
#'
#' @param p The parameter list.
#' @param name The argument name, used in the error message.
#'
#' @returns `p`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_predprey_parameters <- function(p, name = "p") {
  required <- c("alpha", "beta", "gamma", "delta", "sigma_x", "sigma_y")
  if (!is.list(p)) {
    stop(
      sprintf(
        "`%s` must be a named list of predator-prey parameters; it is %s.",
        name, .aci_describe(p)
      ),
      call. = FALSE
    )
  }
  absent <- setdiff(required, names(p))
  if (length(absent) > 0L) {
    stop(
      sprintf(
        "`%s` is missing the parameter%s %s.",
        name,
        if (length(absent) > 1L) "s" else "",
        paste(sprintf("`%s`", absent), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  for (entry in required) {
    .aci_check_scalar(p[[entry]], sprintf("%s$%s", name, entry))
  }
  invisible(p)
}

#' Validate a predator-prey direction
#'
#' @param direction The direction string.
#'
#' @returns The matched direction.
#'
#' @noRd
#' @keywords internal
.aci_check_direction <- function(direction) {
  if (!is.character(direction) || length(direction) != 1L ||
        !direction %in% c("predator_to_prey", "prey_to_predator")) {
    stop(
      paste0(
        "`direction` must be either \"predator_to_prey\", in which the prey ",
        "is observed and the predator is latent, or \"prey_to_predator\", in ",
        "which the predator is observed and the prey is latent."
      ),
      call. = FALSE
    )
  }
  direction
}

#' Conditional Gaussian components of the noisy predator-prey model
#'
#' Builds the conditional Gaussian components list for a stochastic
#' Lotka-Volterra predator-prey pair, given an observed signal, the model
#' parameters, and which of the two causal directions is being studied.
#'
#' The system is
#' \deqn{\mathrm{d}x = (\beta x y - \alpha x)\,\mathrm{d}t +
#'   \sigma_x\,\mathrm{d}W_x}
#' \deqn{\mathrm{d}y = (\gamma y - \delta x y)\,\mathrm{d}t +
#'   \sigma_y\,\mathrm{d}W_y}
#' with `x` the predator and `y` the prey. Either variable may be taken as the
#' observed one, and the choice determines which causal question is asked.
#'
#' In both directions the latent variable's self-drift depends on the observed
#' state, so `L_y` is a vector rather than a scalar. This is what distinguishes
#' the model from the dyad, whose latent self-drift is the constant `-d_y`.
#'
#' @param x Numeric vector. The observed signal, one value per time step. Which
#'   population this is depends on `direction`.
#' @param p A named list with elements `alpha`, `beta`, `gamma`, `delta`,
#'   `sigma_x` and `sigma_y`, each a finite numeric scalar. The `parameters`
#'   entry of an [aci_predprey_model()] is exactly this list.
#' @param direction Character scalar. `"predator_to_prey"` observes the prey and
#'   treats the predator as latent, asking whether the predator influences the
#'   future of the prey. `"prey_to_predator"` is the converse.
#'
#' @returns A conditional Gaussian components list; see [aci_components].
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_predprey_model("prey_to_predator")
#' sim <- aci_simulate(model, n = 2000, dt = 0.005, seed = 1)
#' comp <- aci_predprey_components(sim$x, model$parameters, model$direction)
#' str(comp)
#'
#' @seealso [aci_components], [aci_predprey_model()], [aci_dyad_components()]
#' @export
aci_predprey_components <- function(x, p, direction) {
  .aci_check_signal(x)
  .aci_check_predprey_parameters(p)
  direction <- .aci_check_direction(direction)

  if (identical(direction, "predator_to_prey")) {
    # The prey is observed and the predator is latent. The prey equation
    # supplies the observed drift, and the predator's own growth rate is set by
    # the prey it sees.
    comp <- list(
      L_x = -p$delta * x,
      f_x = p$gamma * x,
      L_y = p$beta * x - p$alpha,
      f_y = rep(0, length(x)),
      S_xoS_x = p$sigma_y^2,
      S_yoS_y = p$sigma_x^2
    )
  } else {
    # The predator is observed and the prey is latent.
    comp <- list(
      L_x = p$beta * x,
      f_x = -p$alpha * x,
      L_y = p$gamma - p$delta * x,
      f_y = rep(0, length(x)),
      S_xoS_x = p$sigma_x^2,
      S_yoS_y = p$sigma_y^2
    )
  }
  comp$S_yoS_x <- 0
  comp$S_xoS_y <- 0
  comp
}

#' Noisy predator-prey model
#'
#' Builds the stochastic Lotka-Volterra predator-prey model as an `aci_model`
#' object, ready for [aci_simulate()] and [aci()]. The parameter defaults are
#' those of the reference implementation.
#'
#' The model is studied in two directions, and `direction` selects which. The
#' two are genuinely different questions rather than one question and its
#' mirror. The observed process differs, the latent process differs, and the
#' causal-information metric is not symmetric between them.
#'
#' @param direction Character scalar, `"prey_to_predator"` (the default) or
#'   `"predator_to_prey"`. See [aci_predprey_components()].
#' @param alpha Numeric scalar. The predator's natural death rate.
#' @param beta Numeric scalar. The effect of prey availability on predator
#'   growth.
#' @param gamma Numeric scalar. The prey's natural growth rate.
#' @param delta Numeric scalar. The effect of predator presence on prey
#'   decline.
#' @param sigma_x,sigma_y Numeric scalars. The noise amplitudes of the predator
#'   and prey processes.
#' @param x0,y0 Numeric scalars. Initial predator and prey populations.
#'
#' @returns An object of class `aci_model`.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_predprey_model()
#' model
#'
#' @seealso [aci_dyad_model()], [aci_predprey_components()], [aci()]
#' @export
aci_predprey_model <- function(direction = "prey_to_predator",
                               alpha = 0.4, beta = 0.1, gamma = 1.1,
                               delta = 0.4, sigma_x = 0.3, sigma_y = 0.3,
                               x0 = 4, y0 = 4) {
  direction <- .aci_check_direction(direction)
  p <- list(
    alpha = alpha, beta = beta, gamma = gamma, delta = delta,
    sigma_x = sigma_x, sigma_y = sigma_y
  )
  .aci_check_predprey_parameters(p)
  .aci_check_scalar(x0, "x0")
  .aci_check_scalar(y0, "y0")

  # The direction decides which population is watched, so it decides which
  # noise coefficient the filter has to invert and which initial state is the
  # observed one.
  if (identical(direction, "predator_to_prey")) {
    if (sigma_y == 0) {
      stop(
        paste0(
          "`sigma_y` must be non-zero for direction \"predator_to_prey\": ",
          "the prey is then the observed process, and its noise covariance ",
          "is what the filter inverts."
        ),
        call. = FALSE
      )
    }
    model <- aci_cgns_model(
      L_x = function(x) -delta * x,
      f_x = function(x) gamma * x,
      L_y = function(x) beta * x - alpha,
      f_y = 0,
      S_xoS_x = sigma_y^2,
      S_yoS_y = sigma_x^2,
      S_yoS_x = 0,
      x0 = y0,
      y0 = x0,
      label = "noisy predator-prey model (predator to prey)",
      parameters = p
    )
  } else {
    if (sigma_x == 0) {
      stop(
        paste0(
          "`sigma_x` must be non-zero for direction \"prey_to_predator\": ",
          "the predator is then the observed process, and its noise ",
          "covariance is what the filter inverts."
        ),
        call. = FALSE
      )
    }
    model <- aci_cgns_model(
      L_x = function(x) beta * x,
      f_x = function(x) -alpha * x,
      L_y = function(x) gamma - delta * x,
      f_y = 0,
      S_xoS_x = sigma_x^2,
      S_yoS_y = sigma_y^2,
      S_yoS_x = 0,
      x0 = x0,
      y0 = y0,
      label = "noisy predator-prey model (prey to predator)",
      parameters = p
    )
  }
  model$direction <- direction
  model
}
