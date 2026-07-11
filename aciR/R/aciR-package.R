#' aciR: Assimilative Causal Inference for Conditional Gaussian Nonlinear
#' Systems
#'
#' The aciR package implements assimilative causal inference (ACI) for
#' conditional Gaussian nonlinear systems (CGNS). A CGNS pairs an observed
#' signal with an unobserved component whose conditional statistics, given the
#' observed path, are Gaussian and available in closed form.
#'
#' The package provides the forward filter ([aci_filter()]), the backward
#' smoother ([aci_smoother()]) and the relative-entropy causal-information
#' metric ([aci_metric()]) that contrasts the smoother posterior against the
#' filter posterior. These three functions form the numerical core and are
#' written against a general CGNS components list, so they apply beyond any
#' single model. The nonlinear dyad model with intermittent extreme events is
#' supplied as a worked example through [aci_dyad_components()] and
#' [aci_simulate_dyad()].
#'
#' Above the core sits a model layer. [aci_cgns_model()] wraps a CGNS in an
#' `aci_model` object, [aci_dyad_model()] builds the flagship dyad model,
#' [aci_simulate()] draws a realisation of a model, and [aci()] is the single
#' entry point that runs the filter, the smoother and the metric together.
#'
#' The method reimplemented here is due to Andreou, Chen and Bollt (2026).
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @importFrom stats rnorm
#' @keywords internal
"_PACKAGE"
