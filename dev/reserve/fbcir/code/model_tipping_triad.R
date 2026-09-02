## acir reserve file
## Origin: aci/R/benchmark_models.R:179-282
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: andreou2026cir Sec. 3.1 eqs. (24a)-(24c); FBCIR_code-main/climate_tipping_*.m.
## Verbatim copy from the aci 0.0.30 sources; not modified.

#' Tipping triad (andreou2026cir Section 3.1, eq. 24a-c)
#'
#' @param eps Fast/slow time-scale ratio. The paper studies both regimes: the
#'   default `0.01` gives the bistable slow-transition regime, while the
#'   supplied FBCIR scripts run the intermittent regime at `0.1`.
#' @param params Optional complete parameter list.
#' @param partition Joint or conditional causal partition.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#' @examples
#' model_tipping_triad()
#'
#' @export
model_tipping_triad <- function(eps = 0.01, params = NULL,
                                partition = c("joint", "y_given_gamma",
                                              "gamma_given_y")) {
  partition <- match.arg(partition)
  defaults <- list(d_x = 1/3, alpha = 4, sigma_x = 0.2, d_y = 0.2,
                   beta = -0.8, sigma_y = 0.3, d_gamma = 0.5,
                   gamma_bar = 1, sigma_gamma = 2)
  p <- .complete_scalar_params(params, defaults, "tipping-triad",
    c("d_x", "alpha", "sigma_x", "d_y", "sigma_y", "d_gamma",
      "sigma_gamma"))
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    aci_abort("aci_error_model_contract", "eps must be one positive finite number.")

  if (partition == "joint") {
    # Observed x; hidden state ordered as (y, gamma).
    m <- cgns_model(
      Lx = function(t, x) matrix(c(-p$alpha, 0), 1, 2),
      fx = function(t, x) x - p$d_x * x^3,
      Ly = function(t, x) matrix(c(-p$d_y / eps, 0, x, -p$d_gamma), 2, 2),
      fy = function(t, x) c(p$beta, p$d_gamma * p$gamma_bar),
      Sx1 = function(t, x) matrix(p$sigma_x, 1, 1),
      Sy2 = function(t, x) diag(c(p$sigma_y / sqrt(eps), p$sigma_gamma)),
      k = 1, l = 2, name = sprintf("tipping_triad[joint,eps=%g]", eps))
    m$meta$vars <- list(observed = "x", hidden = c("y", "gamma"))
    m$meta$conditioning <- NULL
    m$meta$target_obs_idx <- 1L
    m$meta$conditioning_obs_idx <- integer()
    m$meta$causal_link <- "(y,gamma) -> x"
    m$meta$ic_default <- list(x0 = c(x = 0), y0 = c(y = 0, gamma = 0))
  } else if (partition == "y_given_gamma") {
    # Conditional FBCIR question y(t) -> x | gamma.  The observed state is
    # ordered as (x, gamma), so gamma remains in the dynamics and can be masked
    # from the likelihood (or supplied as a prescribed forcing) by the caller.
    m <- cgns_model(
      Lx = function(t, x) matrix(c(-p$alpha, 0), 2, 1),
      fx = function(t, x) c(x[1] - p$d_x * x[1]^3,
                             -p$d_gamma * (x[2] - p$gamma_bar)),
      Ly = function(t, x) matrix(-p$d_y / eps, 1, 1),
      fy = function(t, x) p$beta + x[2] * x[1],
      Sx1 = function(t, x) diag(c(p$sigma_x, p$sigma_gamma)),
      Sy2 = function(t, x) matrix(p$sigma_y / sqrt(eps), 1, 1),
      k = 2, l = 1, name = sprintf("tipping_triad[y->x|gamma,eps=%g]", eps))
    m$meta$vars <- list(observed = c("x", "gamma"), hidden = "y")
    m$meta$conditioning <- list(target = "x", conditioned_on = "gamma")
    m$meta$target_obs_idx <- 1L
    m$meta$conditioning_obs_idx <- 2L
    m$meta$causal_link <- "y -> x | gamma"
    m$meta$ic_default <- list(x0 = c(x = 0, gamma = 0), y0 = c(y = 0))
  } else {
    # Conditional FBCIR question gamma(t) -> y | x.  The observed state is
    # ordered as (x, y); x is the resolved conditioning channel.
    m <- cgns_model(
      Lx = function(t, x) matrix(c(0, x[1]), 2, 1),
      fx = function(t, x) c(x[1] - p$d_x * x[1]^3 - p$alpha * x[2],
                             p$beta - (p$d_y / eps) * x[2]),
      Ly = function(t, x) matrix(-p$d_gamma, 1, 1),
      fy = function(t, x) p$d_gamma * p$gamma_bar,
      Sx1 = function(t, x) diag(c(p$sigma_x, p$sigma_y / sqrt(eps))),
      Sy2 = function(t, x) matrix(p$sigma_gamma, 1, 1),
      k = 2, l = 1, name = sprintf("tipping_triad[gamma->y|x,eps=%g]", eps))
    m$meta$vars <- list(observed = c("x", "y"), hidden = "gamma")
    m$meta$conditioning <- list(target = "y", conditioned_on = "x")
    m$meta$target_obs_idx <- 2L
    m$meta$conditioning_obs_idx <- 1L
    m$meta$causal_link <- "gamma -> y | x"
    m$meta$ic_default <- list(x0 = c(x = 0, y = 0), y0 = c(gamma = 0))
  }
  m$meta$params <- c(p, eps = eps)
  m$meta$partition <- partition
  m$meta$provenance <- "andreou2026cir Section 3.1, equations (24a)-(24c)"
  m$meta$source_files <- switch(partition,
    joint = character(),
    y_given_gamma = "FBCIR_code-main/climate_tipping_y_bifurcation_driven.m",
    gamma_given_y = "FBCIR_code-main/climate_tipping_gamma_noise_induced.m")
  m$meta$source_status <- if (partition == "joint")
    "paper checked; joint partition is not a supplied MATLAB causal query" else
    "paper + MATLAB checked"
  m$meta$partition_caveat <- if (partition == "joint") {
    paste("The joint posterior (y,gamma) | x is the paper's CGNS form, but",
          "the supplied scripts analyse the two conditional causal links.")
  } else {
    paste("The conditioning channel remains in the observed state so its",
          "dynamics are retained; conditional ACI must mask its likelihood",
          "channel or prescribe its realised path, as in the supplied script.")
  }
  m
}
