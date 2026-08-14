# -- composite Simpson quadrature ---------------------------------------------
#
# The objective causal influence range is an integral of the causal-information
# metric over the lagged observational time, and of the subjective range over
# the threshold. The reference implementation evaluates both with a composite
# Simpson rule in preference to the trapezoidal rule it also carries, so this
# package uses Simpson too: on a smooth integrand the two differ by more than
# the tolerance the rest of the package is graded to, and the choice of rule is
# therefore part of reproducing the published quantity rather than an
# implementation detail.
#
# Simpson's rule carries its own oracle. It is exact for polynomials up to
# cubic, which is a statement about the rule rather than about any
# implementation of it, so the tests grade this code against integrals whose
# values are known in closed form and which no implementation here produced.

#' Composite Simpson quadrature
#'
#' Integrates sampled values by the composite Simpson 1/3 rule. With an even
#' number of intervals the rule applies throughout. With an odd number the
#' rule is applied to all but the final interval, which is closed with the
#' Simpson 3/8 rule over the last three, so the result stays third-order
#' accurate rather than dropping to the trapezoidal order at the tail.
#'
#' @param y Numeric vector. The sampled integrand; at least two values.
#' @param x Numeric vector or `NULL`. The abscissae, the same length as `y`
#'   and strictly increasing. When `NULL`, unit spacing is assumed and the
#'   result is the integral with respect to the sample index, which the caller
#'   scales by the step width.
#'
#' @returns A numeric scalar, the estimated integral.
#'
#' @noRd
#' @keywords internal
.aci_simpson <- function(y, x = NULL) {
  n <- length(y)
  if (!is.numeric(y) || n < 2L || anyNA(y)) {
    stop(
      "`y` must be a numeric vector of at least two non-missing values.",
      call. = FALSE
    )
  }
  if (is.null(x)) {
    x <- seq_len(n)
  } else if (length(x) != n || !is.numeric(x) || anyNA(x)) {
    stop(
      "`x` must be numeric, non-missing and the same length as `y`.",
      call. = FALSE
    )
  } else if (any(diff(x) <= 0)) {
    stop("`x` must be strictly increasing.", call. = FALSE)
  }

  if (n == 2L) {
    # A single interval admits no Simpson panel; the trapezoid is the only
    # rule available and is exact for the linear integrand it can resolve.
    return((x[2L] - x[1L]) * (y[1L] + y[2L]) / 2)
  }

  intervals <- n - 1L
  # A Simpson panel spans two intervals, so an odd interval count leaves one
  # over. The rule runs over the longest even-interval prefix and the final
  # interval is closed separately.
  if (intervals %% 2L == 0L) {
    .aci_simpson_13(y, x, 1L, n)
  } else {
    .aci_simpson_13(y, x, 1L, n - 1L) + .aci_simpson_closure(y, x, n)
  }
}

#' Composite Simpson 1/3 rule over an even number of intervals
#'
#' @param y Numeric vector. The sampled integrand.
#' @param x Numeric vector. The abscissae.
#' @param from,to Integer scalars. First and last index of the range; the
#'   number of intervals between them must be even and at least two.
#'
#' @returns A numeric scalar.
#'
#' @noRd
#' @keywords internal
.aci_simpson_13 <- function(y, x, from, to) {
  if (to - from < 2L) {
    return(0)
  }
  idx <- seq.int(from, to)
  m <- length(idx)
  # Panel by panel, so that a non-uniform grid is integrated correctly rather
  # than by assuming the spacing of the first panel throughout.
  left <- seq.int(1L, m - 2L, by = 2L)
  h1 <- x[idx[left + 1L]] - x[idx[left]]
  h2 <- x[idx[left + 2L]] - x[idx[left + 1L]]
  h <- h1 + h2
  # Simpson on an unequal pair of intervals; this reduces to the familiar
  # h/6 * (f0 + 4 f1 + f2) when the two are equal.
  sum(
    h / 6 * (
      (2 - h2 / h1) * y[idx[left]] +
        h^2 / (h1 * h2) * y[idx[left + 1L]] +
        (2 - h1 / h2) * y[idx[left + 2L]]
    )
  )
}

#' Close an odd interval count over the final interval
#'
#' Integrates the quadratic through the last three samples over the last
#' interval alone. This is the closure the reference implementation uses, and
#' the choice is deliberate rather than incidental: it is what makes the
#' objective causal influence range reproduce the published number.
#'
#' An earlier version closed with a Simpson 3/8 panel over the last three
#' intervals. That is a slightly more accurate rule on an equally spaced grid
#' -- measurably so -- but it is a different number, and on a truncated
#' comparison horizon the difference reaches 1e-7 on the quantity the method
#' leads with. It also assumed equal spacing, which is wrong for the
#' logarithmic threshold grid the exact objective range integrates over. The
#' form below is exact for unequal spacing, so both problems close together.
#'
#' Derived rather than transcribed. With the last three abscissae at
#' `t = -h1, 0, h2` relative to the penultimate point, integrating the Lagrange
#' basis over `[0, h2]` gives the weights below; at `h1 == h2` they reduce to
#' the familiar `h/12 * (-y0 + 8*y1 + 5*y2)`, which the tests assert.
#'
#' @param y Numeric vector. The sampled integrand.
#' @param x Numeric vector. The abscissae.
#' @param n Integer scalar. Index of the final sample.
#'
#' @returns A numeric scalar.
#'
#' @noRd
#' @keywords internal
.aci_simpson_closure <- function(y, x, n) {
  h1 <- x[n - 1L] - x[n - 2L]
  h2 <- x[n] - x[n - 1L]
  span <- h1 + h2
  -h2^3 / (6 * h1 * span) * y[n - 2L] +
    h2 * (h2 + 3 * h1) / (6 * h1) * y[n - 1L] +
    h2 * (2 * h2 + 3 * h1) / (6 * span) * y[n]
}
