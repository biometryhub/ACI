## acir reserve file
## Origin: aci/R/causal_metrics.R:1041-1162
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: Self-disclaimed at aci causal_metrics.R:1154-1157 as outside the supplied papers and MATLAB reference code.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Self-labelled at aci causal_metrics.R:1154-1157 as outside the supplied papers and MATLAB reference code; includes its section banner.

################################################################################
# tidy CIR read-out
################################################################################

#' Forward and backward causal influence ranges as one tidy frame
#'
#' Wraps the lag table / CIR machinery and returns a data.frame with one row per
#' anchor: `t`, `tau` (the range, in the model's time units), `strength` (the
#' influence strength there) and `direction`. `nontarget` gives the conditional
#' range. Both directions receive `init`, so the backward pass runs from the same
#' prior as the forward one.
#'
#' @param model A CGNS model.
#' @param obs Observed trajectory.
#' @param init Optional Gaussian initialization.
#' @param at Optional anchor times (used to select forward rows and/or backward
#'   reference times).
#' @param direction Direction or directions to calculate.
#' @param nontarget Optional conditional ACI masking specification.
#' @param tol Adaptive-lag tolerance.
#' @param min_M Minimum discrepancy strength; `"auto"` uses the package
#'   threshold option and masks low-signal ranges.
#' @param method Exact layer-cake objective or the `l1_linf` ratio.
#' @param max_lag Maximum number of positive-lag cells retained per forward row.
#' @param window Consecutive small-tail window used by adaptive storage.
#' @param ... Unused arguments are rejected.
#' @returns A data frame with one row per anchor and direction, carrying `t`,
#'   `tau`, `strength` and `direction`; backward rows also carry
#'   `above_baseline`, the [backward_cir()] validity gate (`NA` on forward
#'   rows).
#'
#' @seealso [forward_cir()], [backward_cir()], [cir_pair()]
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' head(cir_table(m, ob, direction = "forward"))
#' }
#'
#' @export
cir_table <- function(model, obs, init = NULL, at = NULL,
                      direction = c("forward", "backward", "both"),
                      nontarget = NULL, tol = 1e-4, min_M = "auto",
                      method = c("exact", "l1_linf"), max_lag = Inf,
                      window = 3L, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to cir_table().")
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "cir_table() needs a cgns_model.")
  direction <- match.arg(direction)
  method <- match.arg(method)
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    aci_abort("aci_error_dims", "tol must be one finite non-negative number.")
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window < 1L || window != floor(window))
    aci_abort("aci_error_dims", "window must be a positive integer.")
  if (!is.numeric(max_lag) || length(max_lag) != 1L || is.na(max_lag) ||
      max_lag < 1 || (!is.infinite(max_lag) && max_lag != floor(max_lag)))
    aci_abort("aci_error_dims", "max_lag must be a positive integer or Inf.")
  # Validate even on a direction that happens not to use a given tuning
  # argument, so malformed inputs are never silently ignored.
  .cir_min_strength(min_M)
  if (!is.null(at) && (!is.numeric(at) || !length(at) || any(!is.finite(at)) ||
      any(at < min(obs$t) | at > max(obs$t))))
    aci_abort("aci_error_dims", "at must contain finite times within the observation span.")
  out <- list(); forward_table <- NULL
  if (direction %in% c("forward", "both")) {
    tab <- lag_table(model, obs, mode = "forward", tol = tol, init = init,
                     nontarget = nontarget, max_lag = max_lag,
                     window = window)
    forward_table <- tab
    f   <- suppressWarnings(forward_cir(tab, method = method, min_M = min_M))
    idx <- if (is.null(at)) seq_along(f$t) else
      vapply(at, function(a) which.min(abs(f$t - a)), integer(1))
    out$forward <- data.frame(t = f$t[idx], tau = f$tau[idx], strength = f$M[idx],
                              direction = "forward", bound = f$bound,
                              tail_estimate = f$tail_bound[idx],
                              above_baseline = NA)
    attr(out$forward, "lag_table") <- tab
  }
  if (direction %in% c("backward", "both")) {
    Ts <- if (is.null(at)) as.numeric(stats::quantile(obs$t, c(.25, .40, .55, .70, .85))) else at
    b  <- suppressWarnings(backward_cir(model, obs = obs, T = Ts, init = init,
                                        nontarget = nontarget, method = method,
                                        min_M = min_M))
    ab <- if (!is.null(b$meta$per_reference))
      vapply(b$meta$per_reference,
             function(mm) isTRUE(mm$above_baseline), logical(1))
    else rep(isTRUE(b$meta$above_baseline), length(b$t))
    out$backward <- data.frame(t = b$t, tau = b$tau, strength = b$M,
                               direction = "backward", bound = b$bound,
                               tail_estimate = NA_real_,
                               above_baseline = ab)
  }
  res <- do.call(rbind, unname(out))
  rownames(res) <- NULL
  prov <- list(
    source_model = model,
    source_obs_t = obs$t,
    source_obs_x = obs$x,
    nontarget = nontarget,
    actual_init = if (!is.null(forward_table)) forward_table$meta$init else init,
    direction = direction,
    method = method,
    tol = tol,
    min_M = min_M,
    max_lag = max_lag,
    window = as.integer(window),
    source_status = paste(
      "Experimental tidy/application wrapper; the underlying CIR calculations",
      "are provided by the ACI engine, but this wrapper is outside the supplied",
      "papers and MATLAB reference code."))
  attr(res, "aci_provenance") <- prov
  if (!is.null(forward_table)) attr(res, "lag_table") <- forward_table
  class(res) <- c("aci_cir_table", "data.frame")
  res
}
