## acir reserve file
## Origin: aci/R/causal_metrics.R:250-287
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: Sample-based KL with an explicit not-implemented knn stub (aci causal_metrics.R:270-271); no in-package consumer.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Sample/ensemble-flavoured KL; its knn estimator is an explicit not-implemented stub (aci causal_metrics.R:270-271).

#' Relative entropy between two samples
#'
#' Estimates the relative entropy of one sample from another, either by fitting
#' Gaussian moments or by a nearest-neighbour estimator.
#'
#' @param A Numeric matrix of draws from the first distribution, one row per
#'   draw.
#' @param B Numeric matrix of draws from the second distribution.
#' @param estimator Either `"gaussian"` or `"knn"`.
#' @returns 1-length numeric, the estimated relative entropy.
#'
#' @seealso [gaussian_kl()]
#'
#' @examples
#' set.seed(1)
#' empirical_kl(matrix(rnorm(200), ncol = 2), matrix(rnorm(200), ncol = 2))
#'
#' @export
empirical_kl <- function(A, B, estimator = c("gaussian", "knn")) {
  estimator <- match.arg(estimator)
  if (estimator == "knn")
    aci_abort("aci_error_not_implemented", "knn estimator is a v0.2 stub (SPEC-02).")
  A <- as.matrix(A); B <- as.matrix(B)
  if (nrow(A) < 2L || nrow(B) < 2L || ncol(A) != ncol(B) ||
      ncol(A) < 1L || any(!is.finite(c(A, B))))
    aci_abort("aci_error_dims",
              "A and B must be finite matrices with matching columns and at least two rows each.")
  d <- ncol(A)
  if (nrow(A) <= d || nrow(B) <= d)
    aci_abort("aci_error_ensemble_rank",
              "Gaussian empirical KL requires more samples than dimensions in each sample.")
  RA <- stats::cov(A); RB <- stats::cov(B)
  full_rank <- function(R) !is.null(tryCatch(chol(sym(R)), error = function(e) NULL))
  if (!full_rank(RA) || !full_rank(RB))
    aci_abort("aci_error_ensemble_rank",
              "Gaussian empirical KL requires full-rank sample covariance in both samples.")
  gaussian_kl(colMeans(A), RA, colMeans(B), RB, decompose = FALSE)
}
