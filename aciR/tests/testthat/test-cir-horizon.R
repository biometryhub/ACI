# -- causal influence range: comparison horizon and the exact objective range --
#
# Both quantities exist because the package was graded against the reference
# implementation rather than against a transcription of it, and the two were
# found to answer slightly different questions. `horizon` makes the reference's
# convention reachable; `objective_exact` supplies the definition the source
# paper gives, of which the existing `objective` is the efficient underestimate.
#
# The tests below avoid restating the implementation. Where a value has to be
# predicted it is predicted from the DEFINITION -- the length of the comparison
# sequence, the quadrature over the threshold grid -- not from the code.

.cir_fixture <- function(n = 1200L, seed = 4L) {
  set.seed(seed)
  dt <- 0.002
  t <- seq.int(0L, n) * dt
  x <- cumsum(c(0.7, stats::rnorm(n, sd = 0.05)))
  comp <- list(
    L_x = 0.8 + 0.35 * sin(3 * t),
    L_y = -1.3,
    f_x = 0.2 * cos(2 * t) - 0.1 * x,
    f_y = 0.45 - 0.15 * x^2,
    S_xoS_x = 0.34, S_yoS_y = 1.16, S_yoS_x = 0.5, S_xoS_y = 0.5
  )
  list(x = x, comp = comp, dt = dt,
       filt = aci_filter(x, comp, dt, mu0 = 0.3, R0 = 0.1))
}

test_that("the default horizon leaves every reported quantity unchanged", {
  f <- .cir_fixture()
  window <- seq.int(60L, 400L, by = 20L)
  base <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = window)
  full <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = window,
                  horizon = length(f$x))

  # `NULL` and "the whole record" must be the same request, not merely similar:
  # anything else would mean the default silently changed when the argument was
  # introduced.
  expect_identical(base$peak, full$peak)
  expect_identical(base$objective, full$objective)
  expect_identical(base$subjective, full$subjective)
  expect_identical(base$saturated, full$saturated)
})

test_that("the horizon sets the length of the comparison sequence", {
  f <- .cir_fixture()
  aux <- .aci_online_aux(f$x, f$comp, f$dt, f$filt)
  n <- length(f$x)

  # The row at time j compares against observations j .. horizon, so it holds
  # horizon - j + 1 values. Predicted from the definition, not from the code.
  for (h in c(300L, 700L, n)) {
    for (j in c(60L, 200L, 299L)) {
      row <- .aci_cir_row(aux, f$filt, j, n, horizon = h)
      expect_length(row, min(n - j + 1L, h - j + 1L))
    }
  }
})

test_that("a shorter horizon cannot lengthen a subjective range", {
  f <- .cir_fixture()
  window <- seq.int(60L, 300L, by = 30L)
  long <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = window,
                  margin = 0.001)
  short <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = window,
                   margin = 0.001, horizon = 700L)

  # A range is the last position at which the divergence still exceeds a
  # threshold. Looking at less of the record can only move that position
  # earlier, never later -- so this is a property of the quantity, and would
  # fail if the horizon were applied to the fully informed posterior as well.
  both <- is.finite(long$subjective) & is.finite(short$subjective)
  expect_true(any(both))
  expect_true(all(short$subjective[both] <= long$subjective[both] + 1e-12))
})

test_that("the fully informed posterior ignores the horizon", {
  f <- .cir_fixture()
  aux <- .aci_online_aux(f$x, f$comp, f$dt, f$filt)
  n <- length(f$x)
  j <- 80L

  # The truncated row must be the prefix of the untruncated one. If the horizon
  # were also applied to the reference posterior each row would be compared
  # against a different target and the prefix property would break -- which is
  # exactly the mistake this argument is easiest to get wrong in.
  full <- .aci_cir_row(aux, f$filt, j, n, horizon = n)
  cut <- .aci_cir_row(aux, f$filt, j, n, horizon = 500L)
  expect_equal(cut, full[seq_along(cut)], tolerance = 0)
})

test_that("a horizon too close to the reported time marks the time unresolved", {
  f <- .cir_fixture()
  # Two later observations is not a range. The reference records zero here,
  # which cannot be told apart from "no detectable influence"; this must be
  # marked instead.
  r <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = c(100L, 101L),
               horizon = 102L, margin = 0.001)
  expect_true(all(r$saturated))
  expect_true(all(is.na(r$objective)))
  expect_true(all(is.na(r$objective_exact)))
})

test_that("`horizon` rejects values that are not whole and at least one", {
  f <- .cir_fixture(n = 300L)
  args <- list(f$x, f$comp, f$dt, filt = f$filt, window = 40L:60L)
  expect_error(do.call(aci_cir, c(args, list(horizon = 2.5))), "whole number")
  expect_error(do.call(aci_cir, c(args, list(horizon = 0L))), "whole number")
  expect_error(do.call(aci_cir, c(args, list(horizon = -3L))), "whole number")
})

test_that("the exact objective range is the subjective grid integrated", {
  f <- .cir_fixture()
  window <- seq.int(60L, 300L, by = 30L)
  epsilon <- 10^seq(-6, 0.5, length.out = 65L)
  r <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = window,
               epsilon = epsilon, margin = 0.001)

  # Recomputed from the returned subjective ranges and peak by the definition
  # in the source paper -- the subjective ranges averaged over the threshold
  # grid -- rather than by re-running the internal path.
  ord <- order(epsilon)
  for (i in seq_along(window)) {
    if (!is.finite(r$objective_exact[i])) {
      next
    }
    expected <- .aci_simpson(r$subjective[ord, i], epsilon[ord]) / r$peak[i]
    expect_equal(r$objective_exact[i], expected, tolerance = 1e-10)
  }
  expect_true(any(is.finite(r$objective_exact)))
})

test_that("the exact objective range needs the unmasked subjective ranges", {
  f <- .cir_fixture()
  window <- seq.int(60L, 300L, by = 30L)
  # A quadrature over the threshold grid needs every node. If the masked
  # (NA-bearing) subjective matrix were integrated instead, a tightened margin
  # would blank the exact range wherever any single threshold ran long -- so a
  # margin that leaves some subjective entries NA must still yield finite exact
  # ranges at times that are not themselves saturated.
  r <- aci_cir(f$x, f$comp, f$dt, filt = f$filt, window = window, margin = 0.5)
  masked <- apply(r$subjective, 2L, function(col) any(is.na(col)))
  live <- masked & !r$saturated
  if (any(live)) {
    expect_true(all(is.finite(r$objective_exact[live])))
  } else {
    succeed("no reported time had a masked threshold at this margin")
  }
})

test_that("the vector path accepts a horizon and truncates the same way", {
  # Components built directly rather than through a constructor: this test is
  # about the horizon, and routing it through model construction would make a
  # failure ambiguous between the two.
  set.seed(19L)
  n <- 500L
  dt <- 0.01
  x <- rbind(cumsum(stats::rnorm(n + 1L, sd = 0.05)),
             cumsum(stats::rnorm(n + 1L, sd = 0.05)))
  comp <- list(
    L_x = matrix(c(-0.6, 0.2, 0.1, -0.5), 2L, 2L),
    L_y = matrix(c(-0.9, 0.15, 0.05, -0.7), 2L, 2L),
    f_x = matrix(0.05, 2L, n + 1L),
    f_y = matrix(0.10, 2L, n + 1L),
    S_xoS_x = diag(c(0.16, 0.09)),
    S_yoS_y = diag(c(0.25, 0.12)),
    S_yoS_x = matrix(0, 2L, 2L)
  )
  filt <- aci_filter(x, comp, dt, mu0 = c(0, 0), R0 = diag(2L) * 0.1)
  window <- seq.int(30L, 150L, by = 30L)
  long <- aci_cir(x, comp, dt, filt = filt, window = window, margin = 0.001)
  short <- aci_cir(x, comp, dt, filt = filt, window = window, margin = 0.001,
                   horizon = 300L)
  both <- is.finite(long$subjective) & is.finite(short$subjective)
  expect_true(any(both))
  expect_true(all(short$subjective[both] <= long$subjective[both] + 1e-12))
})
