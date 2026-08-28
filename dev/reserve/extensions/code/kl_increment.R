## acir reserve file
## Origin: aci/R/causal_metrics.R:112-131
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: DEAD in aci 0.0.30: zero callers in R/, tests/ or vignettes/.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Zero callers in R/, tests/ or vignettes at 0.0.30; dropped from the 0.1.0 public surface per acir-0.1.0-plan.md section 3.3.

#' Relative entropy gained by an update
#'
#' Relative entropy of an updated Gaussian from the one it replaced, oriented so
#' that the new distribution is the integrating density.
#'
#' @param mu_old Numeric vector, mean before the update.
#' @param R_old Covariance matrix before the update.
#' @param mu_new Numeric vector, mean after the update.
#' @param R_new Covariance matrix after the update.
#' @returns 1-length numeric, the relative entropy gained.
#'
#' @seealso [gaussian_kl()]
#'
#' @examples
#' kl_increment(mu_old = 0, R_old = matrix(2),
#'              mu_new = 0.5, R_new = matrix(1))
#'
#' @export
kl_increment <- function(mu_old, R_old, mu_new, R_new)
  unname(gaussian_kl(mu_new, R_new, mu_old, R_old, decompose = FALSE))
