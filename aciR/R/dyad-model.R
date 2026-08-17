# Nonlinear dyad model with intermittent extreme events ------------------------
#
# The dyad model pairs an observed process x with an unobserved process y whose
# feedback on x makes x intermittently extreme. It is the conditional Gaussian
# nonlinear system used as the worked example throughout the package, and the
# system the packaged independent-oracle fixtures were generated from.
#
# The function below builds the model's components list directly from an
# observed signal and a parameter list. It is the worked example of the
# components schema documented at ?aci_components. The arithmetic mirrors the
# reference implementation's, term for term, which is what makes it a useful
# template for a system the package supplies no constructor for.

#' Validate a dyad parameter list
#'
#' @param p The parameter list.
#' @param name The argument name, used in the error message.
#'
#' @returns `p`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_dyad_parameters <- function(p, name = "p") {
  required <- c("d_x", "d_y", "gamma", "F_x", "F_y", "sigma_x", "sigma_y")
  if (!is.list(p)) {
    stop(
      sprintf(
        "`%s` must be a named list of dyad parameters; it is %s.",
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
  if (p$sigma_x == 0) {
    stop(
      sprintf(
        paste0(
          "`%s$sigma_x` must be non-zero: its square is the ",
          "observation-noise covariance, which the filter inverts."
        ),
        name
      ),
      call. = FALSE
    )
  }
  invisible(p)
}

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
#' Most users should reach for [aci_dyad_model()] and [aci()] instead. This
#' function is the worked example of the components schema, and it serves as
#' the template for a conditional Gaussian system the package supplies no
#' constructor for. See [aci_components] for the schema itself.
#'
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param p A named list of dyad parameters with elements `d_x`, `d_y`,
#'   `gamma`, `F_x`, `F_y`, `sigma_x` and `sigma_y`, each a finite numeric
#'   scalar. The `parameters` entry of an [aci_dyad_model()] is exactly this
#'   list.
#'
#' @returns A conditional Gaussian components list; see [aci_components].
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 2000, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' str(comp)
#'
#' @seealso [aci_components], [aci_dyad_model()], [aci_filter()]
#' @export
aci_dyad_components <- function(x, p) {
  .aci_check_signal(x)
  .aci_check_dyad_parameters(p)
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
