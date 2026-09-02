## acir reserve file
## Origin: aci/R/assimilation.R:585-610
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: Lag-table storage diagnostic with zero callers in aci R/ or tests/.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Lag-table storage diagnostic with zero callers in aci R/ or tests/; can return with adaptive-lag work (plan section 3.3).

#' Retained-cell profile of a lag table
#'
#' Reports how many lag cells adaptive storage retained at each anchor time.
#'
#' @param x A `lag_table` object.
#' @returns A data frame with one row per anchor time, giving the retained cell
#'   count and the corresponding lag reach.
#'
#' @seealso [lt_tail_bound()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' head(truncation_profile(tb))
#'
#' @export
truncation_profile <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  if (is.null(x$L) || is.null(x$tailbnd))
    aci_abort("aci_error_dims", "This lag table mode has no truncation profile.")
  data.frame(j = seq_along(x$t), t = x$t, L = x$L,
             tail_estimate = x$tailbnd)
}
