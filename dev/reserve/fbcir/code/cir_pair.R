## acir reserve file
## Origin: aci/R/causal_metrics.R:982-1021
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Runs forward and backward together; meaningless without a public backward CIR.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Runs forward and backward together; senseless without a public backward CIR (plan section 3.3).

#' Forward and backward causal influence ranges together
#'
#' Runs the metric once and returns the forward and backward ranges alongside
#' it, so that both directions share one filter pass and one prior.
#'
#' @param model A `cgns_model` object; the closed-form engine is required.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param engine One of `"auto"` or `"cgns"`.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param ... Passed to [aci()].
#' @returns An object of class `cir_pair`, a list with the `aci` result and the
#'   `forward` and `backward` ranges.
#'
#' @seealso [forward_cir()], [backward_cir()]
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' cir_pair(m, ob)
#' }
#'
#' @export
cir_pair <- function(model, obs, engine = "auto", nontarget = NULL, ...) {
  resolved_engine <- if (identical(engine, "auto")) {
    if (inherits(model, "cgns_model")) "cgns" else "ensemble"
  } else engine
  if (!identical(resolved_engine, "cgns"))
    aci_abort("aci_error_not_implemented",
              "cir_pair() requires the closed-form CGNS lag table engine.")
  a <- aci(model, obs, engine = engine, nontarget = nontarget,
           keep = "table", ...)
  fwd_tab <- a$table
  bwd_tab <- lag_table(model, as_obs(obs), mode = "one_lag", nontarget = nontarget,
                       filter = a$paths$filter,
                       init = a$handles$init)
  structure(list(aci = a, forward = forward_cir(fwd_tab),
                 backward = backward_cir(bwd_tab)), class = "cir_pair")
}
