## acir reserve file
## Origin: aci/R/assimilation.R:486-513
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: One-lag accessor; its only consumers are backward CIR and the one_lag arm of as.data.frame.lag_table.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Accessor for the one-lag sequence; only consumers are backward CIR and the one_lag arm of as.data.frame.lag_table.

#' One-lag sequence of a lag table
#'
#' Accesses the single-observation-lag entries used by the backward causal
#' influence range.
#'
#' @param x A `lag_table` object.
#' @returns Numeric vector of the one-lag divergences.
#'
#' @seealso [backward_cir()], [lt_row()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "one_lag")
#' head(lt_onelag(tb))
#'
#' @export
lt_onelag <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  out <- x$onelag %||% aci_abort("aci_error_dims", "Table was not built in 'one_lag' mode.")
  # A one-lag discrepancy exists for the N intervals only. The internal N+1
  # allocation has a terminal zero sentinel used by the sweep; exposing that
  # sentinel makes every suffix minimum zero and corrupts exact backward CIR.
  if (length(out) == length(x$t)) out <- out[-length(out)]
  out
}
