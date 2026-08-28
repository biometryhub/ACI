## acir reserve file
## Origin: aci/R/causal_metrics.R:210-247
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: Only in-package consumer is the excluded extremes family (aci extremes.R:498,503,518,581).
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Only in-package consumer is the excluded extremes family (R/extremes.R:498,503,518,581); moser2026extremes scope.

#' Relative entropy of two Gaussians projected onto a direction
#'
#' Projects both distributions onto a unit direction and returns the relative
#' entropy of the resulting one-dimensional Gaussians.
#'
#' @param v Numeric direction vector, normalised internally; must be non-zero.
#' @param mu0 Numeric vector, mean of the reference distribution.
#' @param R0 Covariance matrix of the reference distribution.
#' @param muE Numeric vector, mean of the compared distribution.
#' @param RE Covariance matrix of the compared distribution.
#' @returns 1-length numeric, the projected relative entropy.
#'
#' @seealso [sensitive_directions()]
#'
#' @examples
#' projected_kl(v = c(1, 0), mu0 = c(0, 0), R0 = diag(2),
#'              muE = c(1, 0), RE = diag(c(0.5, 1)))
#'
#' @export
projected_kl <- function(v, mu0, R0, muE, RE) {
  v <- as.numeric(v); nv <- sqrt(sum(v^2))
  if (!length(v) || !is.finite(nv) || nv <= 0)
    aci_abort("aci_error_dims", "v must be a finite non-zero direction vector.")
  mu0 <- as.numeric(mu0); muE <- as.numeric(muE)
  R0 <- as.matrix(R0); RE <- as.matrix(RE); l <- length(v)
  if (length(mu0) != l || length(muE) != l || !is.numeric(R0) ||
      !is.numeric(RE) || !identical(dim(R0), c(l, l)) ||
      !identical(dim(RE), c(l, l)) || any(!is.finite(c(mu0, muE, R0, RE))))
    aci_abort("aci_error_dims", "Projected distributions must match the direction dimension.")
  .strict_chol(R0, "R0")
  .strict_chol(RE, "RE")
  v <- v / nv
  s0 <- drop(t(v) %*% R0 %*% v); sE <- drop(t(v) %*% RE %*% v)
  if (!is.finite(s0) || !is.finite(sE) || s0 <= 0 || sE <= 0)
    aci_abort("aci_error_spd", "Projected variances must be positive.")
  dm <- drop(t(v) %*% (muE - mu0))
  0.5 * (sE / s0 + dm^2 / s0 - 1 + log(s0 / sE))
}
