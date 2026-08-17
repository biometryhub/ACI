# Vector-valued conditional Gaussian core --------------------------------------
#
# The scalar core in aci-core.R integrates one observed and one unobserved
# component. This file integrates the general case, an observed process of
# dimension n_x and an unobserved process of dimension n_y, with matrix-valued
# coefficients that may vary in time.
#
# The scalar path is NOT routed through here. At one dimension the two agree
# bit for bit, and that agreement is graded, but the scalar recursion is some
# thirty times faster and is the code the package's oldest oracle grades.
# Replacing a validated, fast implementation with a general, slower one that
# computes the same numbers would gain nothing, so the public functions
# dispatch on the shape of the components and the scalar path is left as it is.
#
# Three things need more care here than they did in one dimension.
#
#   * The covariance must stay symmetric as well as positive definite. An
#     explicit Euler step breaks symmetry at round-off, and the asymmetry
#     compounds, so each step re-symmetrises rather than assuming it holds.
#   * Positive-definiteness is checked by attempting a Cholesky factorisation,
#     which is the definition rather than a proxy for it.
#   * The relative entropy is evaluated through Cholesky factors, never an
#     explicit inverse and never a determinant computed as a product.

#' Attempt a Cholesky factorisation, returning NULL rather than raising
#'
#' @param m A numeric matrix.
#'
#' @returns The upper-triangular Cholesky factor, or `NULL` when `m` is not
#'   symmetric positive definite.
#'
#' @noRd
#' @keywords internal
.aci_chol <- function(m) {
  if (!all(is.finite(m))) {
    return(NULL)
  }
  tryCatch(chol(m), error = function(e) NULL)
}

#' Coefficient slice at one time step
#'
#' Coefficients are held either as a fixed matrix, constant in time, or as an
#' array whose last margin indexes time. This returns the slice either way, so
#' the recursions read the same whichever was supplied.
#'
#' @param a A matrix, or a three-dimensional array with time last.
#' @param j Integer scalar. The time index.
#'
#' @returns A matrix.
#'
#' @noRd
#' @keywords internal
.aci_slice <- function(a, j) {
  if (length(dim(a)) == 3L) a[, , j, drop = FALSE][, , 1L] else a
}

#' Column of a possibly-constant vector coefficient
#'
#' @param m A numeric vector, constant in time, or a matrix with time in
#'   columns.
#' @param j Integer scalar. The time index.
#'
#' @returns A numeric vector.
#'
#' @noRd
#' @keywords internal
.aci_column <- function(m, j) {
  if (is.matrix(m)) m[, j] else m
}

#' Forward conditional Gaussian filter, vector case
#'
#' @param x Numeric matrix. The observed signal, one row per observed
#'   component and one column per time step.
#' @param comp A vector-valued conditional Gaussian components list.
#' @param dt Numeric scalar. The integration time step.
#' @param mu0 Numeric vector. The initial filtered mean.
#' @param R0 Numeric matrix. The initial filtered covariance.
#'
#' @returns A list with the numeric matrix `mean` and the numeric array `cov`.
#'
#' @noRd
#' @keywords internal
.aci_filter_mv <- function(x, comp, dt, mu0, R0) {
  n <- ncol(x)
  n_y <- length(mu0)
  mean <- matrix(0, n_y, n)
  cov <- array(0, c(n_y, n_y, n))
  mean[, 1L] <- mu0
  cov[, , 1L] <- R0
  mu <- mu0
  R <- R0

  for (j in seq_len(n - 1L) + 1L) {
    L_x <- .aci_slice(comp$L_x, j - 1L)
    L_y <- .aci_slice(comp$L_y, j - 1L)
    S_yoS_y <- .aci_slice(comp$S_yoS_y, j - 1L)
    S_yoS_x <- .aci_slice(comp$S_yoS_x, j - 1L)
    inv <- .aci_slice(comp$S_xoS_x_inv, j - 1L)
    dx <- x[, j] - x[, j - 1L]

    # The gain carries the noise cross-covariance and the covariance-weighted
    # coupling together; at zero cross-covariance it reduces to the familiar
    # R L_x' form.
    aux <- S_yoS_x + R %*% t(L_x)
    mu <- mu + (L_y %*% mu + .aci_column(comp$f_y, j - 1L)) * dt +
      aux %*% inv %*%
        (dx - (L_x %*% mu + .aci_column(comp$f_x, j - 1L)) * dt)
    R <- R + (L_y %*% R + R %*% t(L_y) + S_yoS_y -
                aux %*% inv %*% t(aux)) * dt
    R <- (R + t(R)) / 2

    if (is.null(.aci_chol(R))) {
      .aci_stop_covariance("filter", j, (j - 1L) * dt, NA_real_)
    }
    mean[, j] <- mu
    cov[, , j] <- R
  }
  list(mean = mean, cov = cov)
}

#' Backward conditional Gaussian smoother, vector case
#'
#' @param x Numeric matrix. The observed signal.
#' @param comp A vector-valued conditional Gaussian components list.
#' @param dt Numeric scalar. The integration time step.
#' @param filt A list with `mean` and `cov`, as returned by
#'   [.aci_filter_mv()].
#'
#' @returns A list with the numeric matrix `mean` and the numeric array `cov`.
#'
#' @noRd
#' @keywords internal
.aci_smoother_mv <- function(x, comp, dt, filt) {
  n <- ncol(x)
  n_y <- nrow(filt$mean)
  mean <- matrix(0, n_y, n)
  cov <- array(0, c(n_y, n_y, n))
  mean[, n] <- filt$mean[, n]
  cov[, , n] <- filt$cov[, , n]
  muT <- filt$mean[, n]
  RT <- .aci_slice(filt$cov, n)

  for (j in rev(seq_len(n - 1L))) {
    L_x <- .aci_slice(comp$L_x, j)
    L_y <- .aci_slice(comp$L_y, j)
    S_yoS_y <- .aci_slice(comp$S_yoS_y, j)
    S_yoS_x <- .aci_slice(comp$S_yoS_x, j)
    S_xoS_y <- t(S_yoS_x)
    inv <- .aci_slice(comp$S_xoS_x_inv, j)
    dx <- x[, j + 1L] - x[, j]

    A_j <- L_y - S_yoS_x %*% inv %*% L_x
    B_j <- S_yoS_y - S_yoS_x %*% inv %*% S_xoS_y
    R_f <- .aci_slice(filt$cov, j)
    # Solve against the filtered covariance rather than inverting it. The
    # quantity wanted is B R_f^{-1}, and a solve is both faster and better
    # conditioned than forming the inverse to multiply by it.
    b_over_r <- t(solve(R_f, t(B_j)))

    drift <- L_y %*% muT + .aci_column(comp$f_y, j) -
      b_over_r %*% (filt$mean[, j] - muT)
    transport <- S_yoS_x %*% inv %*%
      (-dx + (L_x %*% muT + .aci_column(comp$f_x, j)) * dt)
    muT <- muT - drift * dt + transport

    damping <- A_j + b_over_r
    RT <- RT - (damping %*% RT + RT %*% t(damping) - B_j) * dt
    RT <- (RT + t(RT)) / 2

    if (is.null(.aci_chol(RT))) {
      .aci_stop_covariance("smoother", j, (j - 1L) * dt, NA_real_)
    }
    mean[, j] <- muT
    cov[, , j] <- RT
  }
  list(mean = mean, cov = cov)
}

#' Causal-information metric, vector case
#'
#' The relative entropy of the smoother posterior from the filter posterior for
#' multivariate Gaussians. Evaluated through Cholesky factors throughout, with
#' the quadratic form by triangular solve, the trace as a sum of squares, and
#' the log-determinant ratio as a difference of sums of logged diagonals.
#' Forming the inverse, or the determinant as a product, would lose accuracy in
#' exactly the regime where the metric is smallest.
#'
#' @param filt A list with `mean` and `cov`, the filtered posterior.
#' @param smooth A list with `mean` and `cov`, the smoothed posterior.
#'
#' @returns A list with the numeric metric vector `value` and the integer
#'   clamp count `n_clamped`.
#'
#' @noRd
#' @keywords internal
.aci_metric_mv <- function(filt, smooth) {
  n <- ncol(filt$mean)
  n_y <- nrow(filt$mean)
  value <- numeric(n)

  for (j in seq_len(n)) {
    R_f <- .aci_slice(filt$cov, j)
    R_s <- .aci_slice(smooth$cov, j)
    c_f <- .aci_chol(R_f)
    c_s <- .aci_chol(R_s)
    if (is.null(c_f) || is.null(c_s)) {
      stop(
        sprintf(
          paste0(
            "The posterior covariance at index %d is not symmetric positive ",
            "definite, so the relative entropy at that step is undefined."
          ),
          j
        ),
        call. = FALSE
      )
    }
    d <- smooth$mean[, j] - filt$mean[, j]
    # For the quadratic form, with R_f = c_f' c_f, solving c_f' z = d gives
    # z = R_f^{-1/2} d, so the form is its sum of squares.
    q <- backsolve(c_f, d, transpose = TRUE)
    # For the trace, tr(R_f^-1 R_s) = || c_s c_f^-1 ||_F^2, obtained by
    # solving rather than inverting. Writing it as a sum of squares also makes
    # the term manifestly non-negative, which the equivalent-looking
    # tr(c_f^-1 R_s c_f^-1) is not. Those two agree only when the factors are
    # diagonal, so a system without off-diagonal structure cannot tell them
    # apart.
    m <- t(backsolve(c_f, t(c_s), transpose = TRUE))
    log_det_ratio <- 2 * (sum(log(diag(c_s))) - sum(log(diag(c_f))))
    value[j] <- 0.5 * (sum(q^2) + sum(m^2) - n_y - log_det_ratio)
  }
  .aci_metric_finish(value)
}
