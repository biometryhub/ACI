## acir reserve file
## Origin: aci/R/benchmark_models.R:285-318
## Source package: aci 0.0.30, git tree 97f6b124
## Category: paper-extremes
## Intended release: after the three MATLAB-backed families
## Reason: moser2026extremes eqs. (4.3)-(4.4); paper checked, no corresponding MATLAB.
## Verbatim copy from the aci 0.0.30 sources; not modified.

#' Hidden damping + forcing pathways model (moser2026extremes eq. 4.3-4.4)
#'
#' @param params Optional complete parameter list.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' model_pathways()
#'
#' @export
model_pathways <- function(params = NULL) {
  defaults <- list(d_u = 0.8, d_g = 0.8, d_b = 1, c = 1.2,
                   F_u = 1, F_g = 1, s_u = 0.5, s_g = 2, s_b = 2.5)
  p <- .complete_scalar_params(params, defaults, "pathways",
    c("d_u", "d_g", "d_b", "c", "s_u", "s_g", "s_b"))
  m <- cgns_model(
    Lx = function(t, x) matrix(c(p$c * x, 1), 1, 2),
    fx = function(t, x) -p$d_u * x + p$F_u,
    Ly = function(t, x) diag(c(-p$d_g, -p$d_b)),
    fy = function(t, x) c(-p$c * x^2 + p$F_g, 0),
    Sx1 = function(t, x) matrix(p$s_u, 1, 1),
    Sy2 = function(t, x) diag(c(p$s_g, p$s_b)),
    k = 1, l = 2, name = "pathways (u | gamma, b)")
  m$meta$params <- p
  m$meta$vars <- list(observed = "u", hidden = c("gamma", "b"))
  m$meta$provenance <- "moser2026extremes equations (4.3)-(4.4)"
  m$meta$source_status <- "paper checked; no corresponding MATLAB supplied"
  m$meta$anti_damping_threshold <- p$d_u / p$c
  m$meta$onset_defaults <- list(kappa = 0.8, T_pre = 1.5)
  m$meta$ic_default <- list(x0 = 1, y0 = c(0, 0))
  m
}
