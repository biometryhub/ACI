# -- grading the fixed-lag online smoother -------------------------------------
#
# The online smoother is graded three ways, none of which needs an external
# reference implementation:
#
#   1. At lag zero it conditions on no future observation, so it must BE the
#      filter. This is an exact identity and is asserted as one.
#
#   2. The implementation reconstructs the update products from cumulative
#      logarithms instead of forming the full O(n^2) triangle of the reference
#      algorithm. A literal transcription of that triangle is carried in this
#      file and the two must agree to machine precision. This is what makes
#      the redesign safe: the fast path is graded against the slow one.
#
#   3. At full lag the estimator is the DISCRETE smoother, while aci_smoother()
#      integrates the CONTINUOUS backward equation. These are different
#      objects that agree only as dt goes to zero, so the identity asserted is
#      the convergence RATE, not equality. A first-order rate is a sharp test:
#      a mis-stated recursion does not converge at a clean order.
#
# Property 3 is the one that is easy to get wrong in the other direction --
# asserting equality here would fail, and asserting nothing would let a wrong
# recursion through.

.aci_test_dyad_params <- function() {
  list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
}

# A deterministic signal keeps every resolution in the convergence study a
# refinement of one underlying function rather than a fresh random path.
.aci_test_signal <- function(dt, t_end = 0.5) {
  tt <- seq(0, t_end, by = dt)
  1 + 0.5 * sin(6 * tt) + 0.2 * cos(11 * tt)
}

.aci_test_pieces <- function(dt, t_end = 0.5) {
  x <- .aci_test_signal(dt, t_end)
  comp <- aci_dyad_components(x, .aci_test_dyad_params())
  filt <- aci_filter(x, comp, dt, mu0 = 2, R0 = 0.1)
  list(x = x, comp = comp, filt = filt, dt = dt)
}

# Literal transcription of the reference algorithm's nested triangle: every
# earlier step is updated explicitly by an ordered product formed factor by
# factor. Deliberately the slow, obvious implementation.
.aci_online_literal <- function(x, comp, dt, filt) {
  n <- length(x)
  aux <- aciR:::.aci_online_aux(x, comp, dt, filt)
  e_j <- aux$E
  prev <- filt$mean[1L]
  for (nn in seq_len(n - 1L) + 1L) {
    cur <- numeric(nn)
    cur[nn] <- filt$mean[nn]
    cur[nn - 1L] <- aux$innov_mean[nn - 1L] + filt$mean[nn - 1L]
    innov <- cur[nn - 1L] - filt$mean[nn - 1L]
    if (nn >= 3L) {
      for (j in rev(seq_len(nn - 2L))) {
        cur[j] <- prev[j] + prod(e_j[j:(nn - 2L)]) * innov
      }
    }
    prev <- cur
  }
  prev
}

test_that("a lag of zero is exactly the filter", {
  p <- .aci_test_pieces(1e-3)
  online <- aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = 0)

  # Exact, not approximate: with no future observation admitted there is
  # nothing for the recursion to add.
  expect_identical(online$mean, p$filt$mean)
  expect_identical(online$cov, p$filt$cov)
  expect_identical(online$lag_effective, 0L)
})

test_that("the logarithmic reconstruction matches the literal triangle", {
  for (dt in c(2e-3, 1e-3)) {
    p <- .aci_test_pieces(dt, t_end = 0.4)
    fast <- aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = Inf)
    slow <- .aci_online_literal(p$x, p$comp, p$dt, p$filt)
    expect_lt(max(abs(fast$mean - slow)), 1e-12)
  }
})

test_that("the update products are non-degenerate", {
  # The whole recursion beyond lag zero is carried by the ordered products of
  # E. If those were ever unity or zero the tests above would still pass while
  # grading nothing, which is exactly the failure mode this package exists to
  # rule out. Assert the products are genuinely damping before trusting them.
  p <- .aci_test_pieces(1e-3, t_end = 0.2)
  aux <- aciR:::.aci_online_aux(p$x, p$comp, p$dt, p$filt)
  d1 <- aciR:::.aci_online_product(aux$cum_log, aux$cum_sign, 1L, 1L)
  d50 <- aciR:::.aci_online_product(aux$cum_log, aux$cum_sign, 1L, 50L)

  expect_true(all(is.finite(c(d1, d50))))
  expect_lt(abs(d1), 1)
  expect_gt(abs(d1), 0)
  # Geometric decay: a fifty-step product must be far smaller than a one-step
  # product, which is the spectral property the truncation relies on.
  expect_lt(abs(d50), abs(d1))
  expect_gt(max(abs(aux$innov_mean)), 0)
})

test_that("full lag converges to the backward smoother at first order", {
  # These are different objects -- the discrete smoother against the Euler
  # discretisation of the continuous backward equation -- so equality is the
  # wrong assertion and the rate is the right one.
  steps <- c(4e-3, 2e-3, 1e-3, 5e-4)
  gaps <- vapply(steps, function(dt) {
    p <- .aci_test_pieces(dt)
    smooth <- aci_smoother(p$x, p$comp, p$dt, p$filt)
    online <- aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = Inf)
    max(abs(online$mean - smooth$mean))
  }, numeric(1))

  # Halving dt must halve the gap: ratios cluster at two, so the estimated
  # order clusters at one.
  order <- log2(gaps[-length(gaps)] / gaps[-1])
  expect_true(all(order > 0.9 & order < 1.1))
  expect_true(all(diff(gaps) < 0))
})

test_that("a longer lag moves the estimate away from the filter", {
  p <- .aci_test_pieces(1e-3)
  distance <- vapply(c(0, 5, 25, 125), function(lag) {
    online <- aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = lag)
    max(abs(online$mean - p$filt$mean))
  }, numeric(1))

  # Monotone in the lag: each additional admitted observation can only add
  # information, never remove it.
  expect_true(all(diff(distance) > 0))
})

test_that("the online smoother rejects an inadmissible lag", {
  p <- .aci_test_pieces(2e-3, t_end = 0.1)
  for (bad in list(-1, 2.5, NA_real_, c(1, 2), "5")) {
    expect_error(
      aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = bad),
      "whole number of time steps"
    )
  }
  expect_error(
    aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = 1, tol = 0),
    "`tol` must be positive"
  )
})
