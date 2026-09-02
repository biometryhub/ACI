## acir reserve file
## Origin: aci/R/compiled_cir.R:293-351
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Prefix-slicing helper existing solely for backward_cir.aci_result.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Prefix-slicing helper existing solely for backward_cir.aci_result's multi-reference prefix sharing.

#' Slice a compiled CGNS bundle at a prefix boundary (internal)
#'
#' Coefficients are immutable and filtering is causal, so backward-CIR reference
#' times can share one maximal-prefix compilation and filter. This constructor
#' keeps the ordinary compiled contract and provenance on a shorter prefix.
#'
#' @param bundle A compiled CGNS bundle.
#' @param N1 Number of grid points retained, at least two.
#' @returns A prefix `compiled_cgns` bundle.
#' @noRd
.slice_compiled_cgns <- function(bundle, N1) {
  if (!inherits(bundle, "compiled_cgns") || length(N1) != 1L ||
      !is.numeric(N1) || !is.finite(N1) || N1 != floor(N1) ||
      N1 < 2L || N1 > bundle$N1)
    aci_abort(
      "aci_error_compiled_contract",
      "A compiled CGNS bundle and a valid prefix length are required."
    )
  N1 <- as.integer(N1)
  N <- N1 - 1L
  k <- bundle$k
  l <- bundle$l
  source_obs <- observed_trajectory(
    bundle$source_obs$t[seq_len(N1)],
    bundle$source_obs$x[seq_len(N1), , drop = FALSE]
  )
  resolved_obs <- observed_trajectory(
    bundle$t[seq_len(N1)], bundle$x[seq_len(N1), , drop = FALSE]
  )
  coefficients <- list(
    Lx = array(bundle$coefficients$Lx[, , seq_len(N1), drop = FALSE],
               c(k, l, N1)),
    fx = matrix(bundle$coefficients$fx[seq_len(N1), , drop = FALSE], N1, k),
    Ly = array(bundle$coefficients$Ly[, , seq_len(N1), drop = FALSE],
               c(l, l, N1)),
    fy = matrix(bundle$coefficients$fy[seq_len(N1), , drop = FALSE], N1, l),
    gxx = array(bundle$coefficients$gxx[, , seq_len(N1), drop = FALSE],
                c(k, k, N1)),
    gyy = array(bundle$coefficients$gyy[, , seq_len(N1), drop = FALSE],
                c(l, l, N1)),
    gyx = array(bundle$coefficients$gyx[, , seq_len(N1), drop = FALSE],
                c(l, k, N1)),
    gxx_weight = array(
      bundle$coefficients$gxx_weight[, , seq_len(N), drop = FALSE],
      c(k, k, N)
    )
  )
  rs <- list(
    model = bundle$model,
    obs = resolved_obs,
    likelihood_idx = bundle$likelihood_idx,
    tag = bundle$nontarget
  )
  .new_compiled_cgns(
    rs, bundle$source_model, source_obs, coefficients,
    correlated_noise = bundle$correlated_noise,
    realization = paste0(bundle$realization, "_prefix")
  )
}
