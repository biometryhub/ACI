## acir reserve file
## Origin: aci/R/assimilation.R:902-1025
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Gaspari-Cohn localization and inflation; sole consumers are ensemble.R and test-06.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The whole Gaspari-Cohn localization and inflation section (gaspari_cohn, localization_spec, apply_inflation). Sole consumers are R/ensemble.R and test-06-ensemble.R.

################################################################################
# Gaspari-Cohn localization and inflation
################################################################################

#' Gaspari-Cohn localization weight
#'
#' Evaluates the compactly supported fifth-order piecewise-rational correlation
#' function used to taper ensemble covariances with distance.
#'
#' @param r Numeric vector of non-negative distances, scaled by the
#'   localization radius.
#' @returns Numeric vector of weights in the unit interval, zero at scaled
#'   distances beyond two.
#'
#' @references
#' Gaspari, G. and Cohn, S. E. (1999). Construction of correlation functions
#' in two and three dimensions. *Quarterly Journal of the Royal Meteorological
#' Society* **125**, 723-757. \doi{10.1002/qj.49712555417}
#'
#' @seealso [localization_spec()]
#'
#' @examples
#' gaspari_cohn(seq(0, 2.5, by = 0.5))
#'
#' @export
gaspari_cohn <- function(r) {
  if (!is.numeric(r) || any(!is.finite(r)))
    aci_abort("aci_error_dims", "r must contain finite numeric distances.")
  attrs <- attributes(r)
  r <- abs(r); out <- numeric(length(r))
  i1 <- r < 1; i2 <- r >= 1 & r < 2
  z <- r[i1]
  out[i1] <- 1 - (5/3) * z^2 + (5/8) * z^3 + (1/2) * z^4 - (1/4) * z^5
  z <- r[i2]
  out[i2] <- 4 - 5 * z + (5/3) * z^2 + (5/8) * z^3 - (1/2) * z^4 +
             (1/12) * z^5 - 2 / (3 * z)
  out <- pmax(out, 0)
  attributes(out) <- attrs
  out
}


#' Localization specification for the ensemble engine
#'
#' Builds the taper matrix applied to ensemble cross-covariances, from state
#' and observation coordinates and a localization radius.
#'
#' @param coords Numeric vector or matrix of hidden-state coordinates.
#' @param radius Positive 1-length numeric localization radius.
#' @param coords_obs Optional coordinates of the observed channels; `NULL`
#'   reuses `coords`.
#' @param distance Either `"euclidean"` or `"cyclic"`.
#' @param period Positive 1-length numeric domain period, required when
#'   `distance` is `"cyclic"`.
#' @returns A list carrying the taper matrix and the settings used to build it.
#'
#' @seealso [gaspari_cohn()], [enkbf()]
#'
#' @examples
#' localization_spec(coords = 1:5, radius = 2)
#'
#' @export
localization_spec <- function(coords, radius, coords_obs = NULL,
                              distance = c("euclidean", "cyclic"), period = NULL) {
  distance <- match.arg(distance)
  if (length(radius) != 1L || !is.finite(radius) || radius <= 0)
    aci_abort("aci_error_dims", "radius must be a finite positive scalar.")
  if (distance == "cyclic" && is.null(period))
    aci_abort("aci_error_dims", "cyclic distance needs `period`.")
  co <- as.matrix(coords); cobs <- if (is.null(coords_obs)) co else as.matrix(coords_obs)
  if (!nrow(co) || !ncol(co) || ncol(cobs) != ncol(co) ||
      any(!is.finite(c(co, cobs))))
    aci_abort("aci_error_dims", "coords and coords_obs must be finite matrices with matching columns.")
  if (distance == "cyclic" &&
      (length(period) != 1L || !is.finite(period) || period <= 0))
    aci_abort("aci_error_dims", "period must be a finite positive scalar.")
  if (distance == "cyclic" && ncol(co) != 1L)
    aci_abort("aci_error_dims", "cyclic distance currently requires one-dimensional coordinates.")
  dfun <- function(a, b) {
    d <- outer(seq_len(nrow(a)), seq_len(nrow(b)),
               Vectorize(function(i, j) sqrt(sum((a[i, ] - b[j, ])^2))))
    if (distance == "cyclic") {
      d <- outer(a[, 1], b[, 1], function(x, y) {
        dd <- abs((x - y) %% period)
        pmin(dd, period - dd)
      })
    }
    d
  }
  structure(list(radius = radius,
                 C2 = gaspari_cohn(dfun(co, co) / radius),        # l x l (hidden-hidden)
                 C1 = gaspari_cohn(dfun(co, cobs) / radius)),     # l x k (hidden-obs)
            class = "localization_spec")
}


#' Apply multiplicative inflation to ensemble members
#'
#' Scales the anomalies of an ensemble about its mean. The inflation parameter
#' is named for the variance factor, while the anomalies are multiplied by its
#' square root.
#'
#' @param members Numeric matrix of ensemble members.
#' @param delta2 Positive 1-length numeric variance inflation factor.
#' @returns A matrix of the same dimension as `members`, with inflated
#'   anomalies about the unchanged ensemble mean.
#'
#' @seealso [enkbf()]
#'
#' @examples
#' apply_inflation(matrix(rnorm(20), nrow = 2), delta2 = 1.1)
#'
#' @export
apply_inflation <- function(members, delta2) {
  members <- as.matrix(members)
  if (ncol(members) < 2L || any(!is.finite(members)))
    aci_abort("aci_error_dims", "members must be a finite matrix with at least two columns.")
  if (length(delta2) != 1L || !is.finite(delta2) || delta2 < 1)
    aci_abort("aci_error_model_contract", "inflation delta^2 must be a finite scalar >= 1.")
  if (delta2 == 1) return(members)
  yb <- rowMeans(members)
  # jiang2026enkbs eq. 15 anomaly rescaling
  yb + sqrt(delta2) * (members - yb)
}
