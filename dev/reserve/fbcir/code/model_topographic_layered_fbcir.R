## acir reserve file
## Origin: aci/R/benchmark_models.R:1052-1151
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: FBCIR_code-main/topographic.m repository case study; includes .fbcir_layered_topographic_parameters.
## Verbatim copy from the aci 0.0.30 sources; not modified.

#' Parameters of the seven-state layered topographic model (internal)
#'
#' @param params Optional named list of parameter overrides.
#' @returns The completed named parameter list.
#' @noRd
.fbcir_layered_topographic_parameters <- function(params) {
  defaults <- list(beta = 2, d = 1 / 40, p = 1,
                   topographic_strength = 2)
  if (is.null(params)) params <- list()
  if (!is.list(params))
    aci_abort("aci_error_model_contract", "params must be a list or NULL.")
  unknown <- setdiff(names(params), names(defaults))
  if (length(unknown))
    aci_abort("aci_error_model_contract",
              sprintf("Unknown layered-topographic parameter(s): %s.",
                      paste(unknown, collapse = ", ")))
  p <- utils::modifyList(defaults, params)
  vals <- unlist(p, use.names = FALSE)
  if (length(vals) != length(defaults) || !is.numeric(vals) || any(!is.finite(vals)) ||
      p$d <= 0 || p$topographic_strength <= 0)
    aci_abort("aci_error_model_contract",
              "Layered-topographic parameters must be finite, with d and topographic_strength positive.")
  p
}

#' Seven-dimensional layered topographic FBCIR repository example
#'
#' This is the additional model in `FBCIR_code-main/topographic.m`, not the
#' spectral topographic model from the ACI extreme-events paper.  The observed
#' state is `(v1,...,v6)` and the hidden state is the mean flow `u`.
#'
#' @param target_mode Streamfunction mode (1, 2, or 3) treated as the target
#'   pair by conditional ACI; this changes metadata, not the physical model.
#' @param params Optional named parameter overrides.
#' @return A `cgns_model`.
#'
#' @examples
#' model_topographic_layered_fbcir()
#'
#' @export
model_topographic_layered_fbcir <- function(target_mode = 1L, params = NULL) {
  if (!is.numeric(target_mode) || length(target_mode) != 1L ||
      !is.finite(target_mode) || target_mode != as.integer(target_mode) ||
      !(target_mode %in% 1:3))
    aci_abort("aci_error_model_contract", "target_mode must be one of 1, 2, or 3.")
  target_mode <- as.integer(target_mode)
  p <- .fbcir_layered_topographic_parameters(params)
  modes <- 3L
  H <- p$topographic_strength * sqrt(2) / 2^(seq_len(modes))
  omega <- H / sqrt(2)
  sigma_v <- c(rep(1 / (20 * sqrt(2)), 4),
               rep(1 / (20 * 3^p$p * sqrt(2)), 2))
  sigma_u <- 1 / (20 * sqrt(2))

  m <- cgns_model(
    Lx = function(t, x) matrix(c(
      x[2] - omega[1],
      -x[1],
      2 * x[4] - omega[2] / 2,
      -2 * x[3],
      3 * x[6] - omega[3] / (sqrt(2) * 3^(p$p + 1)),
      -3 * x[5]
    ), 6, 1),
    fx = function(t, x) c(
      -p$d * x[1] - p$beta * x[2],
      -p$d * x[2] + p$beta * x[1],
      -p$d * x[3] - p$beta * x[4] / 2,
      -p$d * x[4] + p$beta * x[3] / 2,
      -p$d * x[5] - p$beta * x[6] / 3,
      -p$d * x[6] + p$beta * x[5] / 3
    ),
    Ly = function(t, x) matrix(-p$d, 1, 1),
    fy = function(t, x)
      omega[1] * x[1] + 2 * omega[2] * x[3] +
      omega[3] * x[5] / (sqrt(2) * 3^(p$p - 1)),
    # topographic.m first defines sigma_v and then sets S_x=diag(sigma_v)/sqrt(2)
    Sx1 = function(t, x) diag(sigma_v / sqrt(2), 6),
    Sy2 = function(t, x) matrix(sigma_u, 1, 1),
    k = 6, l = 1,
    name = sprintf("topographic_layered_fbcir[u->psi%d|other_modes]", target_mode)
  )
  target <- (2L * target_mode - 1L):(2L * target_mode)
  m$meta$params <- c(p, list(H = H, omega = omega,
                             sigma_v = sigma_v, sigma_u = sigma_u))
  m$meta$vars <- list(observed = paste0("v", 1:6), hidden = "u")
  m$meta$target_mode <- target_mode
  m$meta$target_obs_idx <- target
  m$meta$conditioning_obs_idx <- setdiff(1:6, target)
  m$meta$causal_link <- sprintf("u -> psi%d | other streamfunction modes", target_mode)
  m$meta$ic_default <- list(x0 = stats::setNames(rep(0, 6), paste0("v", 1:6)),
                            y0 = c(u = 0))
  m$meta$provenance <- "repository_case_study_named_in_andreou2026cir_code_availability"
  m$meta$source_files <- "FBCIR_code-main/topographic.m"
  m$meta$paper_scope <- FALSE
  m$meta$matlab_pathwise_rng_parity <- FALSE
  m$meta$simulation_convention <- paste(
    "Euler-Maruyama coefficients match the source. R and MATLAB draw observed",
    "and hidden Wiener channels in different RNG orders, so equal seeds are not pathwise comparable.")
  m
}
