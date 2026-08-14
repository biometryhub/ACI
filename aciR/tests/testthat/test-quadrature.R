# -- grading the composite Simpson rule ----------------------------------------
#
# Simpson's rule is exact for polynomials up to cubic. That is a property of
# the rule, established independently of any implementation of it, so the
# closed-form integrals below are an oracle this package did not author. A
# transcription error that still produced a plausible-looking quadrature would
# fail exactness on a cubic almost surely.

test_that("Simpson is exact for cubics when the interval count is even", {
  # With an even interval count the composite 1/3 rule applies throughout and
  # the classical exactness holds. The closed-form integrals are an oracle this
  # package did not author.
  for (n in c(5L, 21L, 33L)) {
    x <- seq(-1.3, 2.7, length.out = n)
    for (cf in list(c(2, 0, 0, 0), c(1, -3, 0, 0), c(0.5, 1, -2, 0),
                    c(-1, 2, 0.5, -3))) {
      y <- cf[1L] + cf[2L] * x + cf[3L] * x^2 + cf[4L] * x^3
      anti <- function(t) {
        cf[1L] * t + cf[2L] * t^2 / 2 + cf[3L] * t^3 / 3 + cf[4L] * t^4 / 4
      }
      exact <- anti(max(x)) - anti(min(x))
      expect_equal(aciR:::.aci_simpson(y, x), exact, tolerance = 1e-10)
    }
  }
})

test_that("an odd interval count is exact for quadratics, on any spacing", {
  # The odd case closes the final interval with the quadratic through the last
  # three samples, so exactness drops from cubic to quadratic there. That is
  # the reference implementation's rule and is a deliberate choice, recorded in
  # .aci_simpson_closure(); asserting cubic exactness here would be asserting a
  # rule this package does not use.
  #
  # Unequal spacing is the case that matters: the exact objective range
  # integrates over a LOGARITHMIC threshold grid, and the closure this replaced
  # assumed equal spacing and was silently wrong there.
  grids <- list(
    seq(-1.3, 2.7, length.out = 6L),
    seq(-1.3, 2.7, length.out = 22L),
    c(0, 0.1, 0.4, 0.45, 1.0, 2.0),
    10^seq(-6, 0.5, length.out = 8L)
  )
  for (x in grids) {
    expect_true(length(x) %% 2L == 0L)          # odd interval count
    for (cf in list(c(2, 0, 0), c(1, -3, 0), c(0.5, 1, -2))) {
      y <- cf[1L] + cf[2L] * x + cf[3L] * x^2
      anti <- function(t) cf[1L] * t + cf[2L] * t^2 / 2 + cf[3L] * t^3 / 3
      exact <- anti(max(x)) - anti(min(x))
      expect_equal(aciR:::.aci_simpson(y, x), exact,
                   tolerance = 1e-10 * max(1, abs(exact)))
    }
  }
})

test_that("the closure reduces to the equal-spacing rule it generalises", {
  # Derived, not transcribed: integrating the Lagrange basis through the last
  # three abscissae over the final interval must collapse to h/12 * (-y0 + 8 y1
  # + 5 y2) when the two spacings agree. A regression here means the general
  # weights have drifted from the rule they are supposed to generalise.
  for (h in c(0.25, 1, 3.5)) {
    for (y in list(c(2.1, -0.7, 3.3), c(0, 1, 0), c(-4, -4, -4))) {
      x <- c(0, h, 2 * h)
      expect_equal(
        aciR:::.aci_simpson_closure(y, x, 3L),
        h / 12 * (-y[1L] + 8 * y[2L] + 5 * y[3L]),
        tolerance = 1e-12
      )
    }
  }
})

test_that("Simpson converges at fourth order on a non-polynomial", {
  # Exactness alone would also hold for a rule that happened to be right on
  # cubics by accident; the convergence order pins the rule itself.
  f <- function(t) exp(sin(t))
  # stats::integrate() is adaptive Gauss-Kronrod from base R, independent of
  # anything in this package, and converges far tighter than the tolerance
  # the order estimate below needs.
  exact <- stats::integrate(f, 0, 2, rel.tol = 1e-12)$value
  err <- vapply(c(9L, 17L, 33L, 65L), function(n) {
    x <- seq(0, 2, length.out = n)
    abs(aciR:::.aci_simpson(f(x), x) - exact)
  }, numeric(1))
  order <- log2(err[-length(err)] / err[-1])
  expect_true(all(order > 3.7 & order < 4.3))
})

test_that("unit spacing integrates with respect to the sample index", {
  y <- c(1, 4, 9, 16, 25)
  # With unit abscissae the result must equal the explicit-abscissae call.
  expect_equal(aciR:::.aci_simpson(y), aciR:::.aci_simpson(y, seq_along(y)))
  # And a two-point call is the trapezoid.
  expect_equal(aciR:::.aci_simpson(c(2, 6)), 4)
  expect_equal(aciR:::.aci_simpson(c(2, 6), c(0, 3)), 12)
})

test_that("Simpson rejects malformed input", {
  expect_error(aciR:::.aci_simpson(1), "at least two")
  expect_error(aciR:::.aci_simpson(c(1, NA, 3)), "non-missing")
  expect_error(aciR:::.aci_simpson(c(1, 2, 3), c(1, 2)), "same length")
  expect_error(
    aciR:::.aci_simpson(c(1, 2, 3), c(1, 1, 2)), "strictly increasing"
  )
})
