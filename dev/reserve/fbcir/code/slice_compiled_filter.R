## acir reserve file
## Origin: aci/R/compiled_cir.R:354-374
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Prefix-slicing helper existing solely for backward_cir.aci_result.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Prefix-slicing helper for the filter path, same sole consumer.

#' Slice a Gaussian filter path at a prefix boundary (internal)
#'
#' @param filter A filter path generated on the full bundle.
#' @param bundle The matching prefix bundle.
#' @returns A compatible prefix filter.
#' @noRd
.slice_compiled_filter <- function(filter, bundle) {
  N1 <- bundle$N1
  path <- new_da_path(
    bundle$t,
    filter$mean[seq_len(N1), , drop = FALSE],
    filter$cov[, , seq_len(N1), drop = FALSE],
    "filter"
  )
  path$meta <- filter$meta
  path$meta$obs_x <- bundle$x
  path$meta$model <- bundle$model
  path$meta$nontarget <- bundle$nontarget
  path$meta$source_model <- bundle$source_model
  path
}
