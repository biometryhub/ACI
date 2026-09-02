## acir reserve file
## Origin: aci/R/causal_metrics.R:553-573
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Backward-only reducer, andreou2026cir eq. 13.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Backward-only reducer .bwd_lengths (bias-corrected backward CIR metric, andreou2026cir eq. 13).

#' Backward influence duration from the one-lag sequence (internal)
#'
#' @param p Numeric vector of finite, non-negative one-lag divergences.
#' @param dt Positive 1-length numeric step.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @returns Named numeric vector with the range `tau` and the strength `M`.
#' @noRd
.bwd_lengths <- function(p, dt, method, quadrature = "sum",
                        simpson_close = "quadratic") {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p)) || any(p < 0) || length(dt) != 1L || !is.finite(dt) || dt <= 0){
    aci_abort("aci_error_dims", "Backward discrepancy values must be finite and non-negative, with dt > 0.")
  }
  p_bc <- abs(p - p[1])  # bias-corrected backward CIR metric (andreou2026cir eq. 13)
  M <- max(p_bc)
  if (!is.finite(M) || M < 1e-14){ return(c(tau = NA_real_, M = M)) }
  # lagged-grid convention from FBCIR is an O(dt) approximation to T' -> T^- in the paper
  tau <- .calc_tau(p = p_bc, dt = dt, method = method, M = M,
                   direction = "bwd", quadrature = quadrature,
                   simpson_close = simpson_close)
  return(c(tau = tau, M = M))
}
