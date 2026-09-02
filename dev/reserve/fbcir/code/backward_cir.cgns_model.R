## acir reserve file
## Origin: aci/R/causal_metrics.R:1024-1038
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Convenience backward_cir method; dies with the backward family.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Convenience backward_cir method dispatching straight from a model.

# Convenience method: backward CIR straight from (model, obs) without first
# assembling an aci_result (the front door for fitted models).
#' @describeIn backward_cir Backward range straight from a model and its
#'   observations, without first assembling an ACI result.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @export
backward_cir.cgns_model <- function(x, obs, T, nontarget = NULL,
                                    init = NULL, ...) {
  ar <- structure(list(handles = list(model = x, obs = as_obs(obs),
                                      nontarget = nontarget, init = init)),
                  class = "aci_result")
  backward_cir(ar, T = T, ...)
}
