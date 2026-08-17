# Mathematical identities of the conditional Gaussian core, checked against
# oracles the package did not produce: the algebraic fixed points of the
# Kalman-Bucy equations, derived here from the governing ODEs and solved in
# closed form. These are independent of both the R recursion under test and
# the MATLAB fixtures, and they reach parts of the core the dyad fixture
# cannot: the dyad has independent noise, so its oracle never exercises the
# terms carrying the noise cross-covariance.

# Analytic fixed points --------------------------------------------------------

# The filtered covariance of a constant-coefficient CGNS obeys
#   dR/dt = 2 L_y R + S_yoS_y - (S_yoS_x + R L_x)^2 / S_xoS_x
# Setting the derivative to zero and clearing the denominator gives the
# algebraic Riccati equation
#   L_x^2 R^2 - 2 (L_y S_xoS_x - S_yoS_x L_x) R
#     - (S_yoS_y S_xoS_x - S_yoS_x^2) = 0
# whose positive root is the stationary filtered covariance.
.riccati_cov <- function(L_x, L_y, S_xoS_x, S_yoS_y, S_yoS_x) {
  k <- L_y * S_xoS_x - S_yoS_x * L_x
  (k + sqrt(k^2 + L_x^2 * (S_yoS_y * S_xoS_x - S_yoS_x^2))) / L_x^2
}

# Driving the filter with an exactly linear signal, dx = slope * dt, makes the
# mean recursion a linear ODE with constant coefficients; its stationary point
# solves 0 = L_y mu + f_y + aux (slope - L_x mu - f_x) / S_xoS_x.
.stationary_mean <- function(aux, L_x, L_y, f_x, f_y, S_xoS_x, slope) {
  inv <- 1 / S_xoS_x
  -(f_y + aux * inv * (slope - f_x)) / (L_y - aux * inv * L_x)
}

test_that("the filter matches the analytic Riccati fixed point", {
  # A system with a non-zero noise cross-covariance: S_yoS_x enters the filter
  # through aux = S_yoS_x + R L_x, a term the packaged dyad fixture pins at
  # zero and therefore never grades.
  L_x <- 0.8
  L_y <- -0.5
  f_x <- 0.3
  f_y <- 0.4
  S_xoS_x <- 0.5
  S_yoS_y <- 0.9
  S_yoS_x <- 0.4
  slope <- 0.25
  dt <- 0.001
  n <- 20001L
  expect_gt(S_xoS_x * S_yoS_y - S_yoS_x^2, 0)

  R_star <- .riccati_cov(L_x, L_y, S_xoS_x, S_yoS_y, S_yoS_x)
  # The oracle checks itself: R* must annihilate the Riccati residual.
  residual <- 2 * L_y * R_star + S_yoS_y -
    (S_yoS_x + R_star * L_x)^2 / S_xoS_x
  expect_lt(abs(residual), 1e-12)
  expect_gt(R_star, 0)

  aux <- S_yoS_x + R_star * L_x
  mu_star <- .stationary_mean(aux, L_x, L_y, f_x, f_y, S_xoS_x, slope)

  x <- 1 + slope * (seq_len(n) - 1L) * dt
  comp <- list(
    L_x = rep(L_x, n), f_x = rep(f_x, n), L_y = L_y, f_y = rep(f_y, n),
    S_xoS_x = S_xoS_x, S_yoS_y = S_yoS_y,
    S_yoS_x = S_yoS_x, S_xoS_y = S_yoS_x
  )
  filt <- aci_filter(x, comp, dt, mu0 = mu_star, R0 = R_star)

  # Started at the fixed point, the recursion must stay there: the fixed point
  # of the ODE is exactly the fixed point of its Euler map, so this is a
  # machine-precision statement, not a discretisation-limited one.
  expect_lt(max(abs(filt$cov - R_star)), 1e-12)
  expect_lt(max(abs(filt$mean - mu_star)), 1e-10)
})

test_that("the smoother matches its analytic fixed point", {
  L_x <- 0.8
  L_y <- -0.5
  f_x <- 0.3
  f_y <- 0.4
  S_xoS_x <- 0.5
  S_yoS_y <- 0.9
  S_yoS_x <- 0.4
  slope <- 0.25
  dt <- 0.001
  n <- 20001L

  R_star <- .riccati_cov(L_x, L_y, S_xoS_x, S_yoS_y, S_yoS_x)
  inv <- 1 / S_xoS_x
  aux <- S_yoS_x + R_star * L_x
  mu_star <- .stationary_mean(aux, L_x, L_y, f_x, f_y, S_xoS_x, slope)

  # The smoother's per-step terms, and the fixed points of its backward
  # recursion once the filter sits at R* and mu*.
  A_j <- L_y - S_yoS_x * inv * L_x
  B_j <- S_yoS_y - S_yoS_x * inv * S_yoS_x
  RT_star <- B_j / (2 * (A_j + B_j / R_star))
  muT_numerator <- -f_y + (B_j / R_star) * mu_star -
    S_yoS_x * inv * (slope - f_x)
  muT_denominator <- L_y + B_j / R_star - S_yoS_x * inv * L_x
  muT_star <- muT_numerator / muT_denominator
  # The backward recursion only relaxes onto its fixed point when the map
  # contracts, which this rate asserts.
  expect_gt(2 * (A_j + B_j / R_star), 0)

  x <- 1 + slope * (seq_len(n) - 1L) * dt
  comp <- list(
    L_x = rep(L_x, n), f_x = rep(f_x, n), L_y = L_y, f_y = rep(f_y, n),
    S_xoS_x = S_xoS_x, S_yoS_y = S_yoS_y,
    S_yoS_x = S_yoS_x, S_xoS_y = S_yoS_x
  )
  filt <- aci_filter(x, comp, dt, mu0 = mu_star, R0 = R_star)
  smooth <- aci_smoother(x, comp, dt, filt)

  # Deep in the interior the backward sweep has relaxed onto its fixed point,
  # away from the terminal boundary condition.
  interior <- 1000L
  expect_lt(abs(smooth$cov[interior] - RT_star), 1e-10)
  expect_lt(abs(smooth$mean[interior] - muT_star), 1e-10)

  # Smoothing cannot be less certain than filtering in the stationary interior.
  expect_lt(RT_star, R_star)
})

# Terminal identity ------------------------------------------------------------

test_that("the smoother equals the filter at the final step", {
  # Conditioning on the whole observed path and on the path up to the final
  # step are the same conditioning, so this identity is exact by construction
  # rather than approached.
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 5000L, seed = 333)
  fit <- aci(sim$x, model)
  n <- length(sim$x)

  expect_identical(fit$smoother$mean[n], fit$filter$mean[n])
  expect_identical(fit$smoother$cov[n], fit$filter$cov[n])
  expect_identical(fit$aci[n], 0)
})

# Zero-information identity ----------------------------------------------------

.zero_information_metric <- function(dt) {
  n <- as.integer(4 / dt) + 1L
  x <- 1 + 0.1 * (seq_len(n) - 1L) * dt
  comp <- list(
    L_x = rep(0, n), f_x = rep(0.2, n), L_y = -0.5, f_y = rep(0.3, n),
    S_xoS_x = 0.5, S_yoS_y = 0.9, S_yoS_x = 0, S_xoS_y = 0
  )
  filt <- aci_filter(x, comp, dt, mu0 = 0.6, R0 = 0.4)
  smooth <- aci_smoother(x, comp, dt, filt)
  max(aci_metric(filt, smooth))
}

test_that("a signal carrying no information scores zero", {
  # With no coupling of the unobserved component into the observed drift and
  # no noise cross-covariance, the observed process is independent of the
  # unobserved one: the future of the signal cannot sharpen the estimate, so
  # the smoother is the filter and the metric is zero.
  #
  # The identity is exact in continuous time only. Under an explicit Euler
  # scheme the forward and backward sweeps do not invert one another exactly,
  # so the metric vanishes with the step rather than at it. Asserting equality
  # would be asserting a property of the mathematics against an approximation
  # of it.
  expect_lt(.zero_information_metric(0.001), 1e-5)
})

test_that("the zero-information residual is second order in the step", {
  coarse <- .zero_information_metric(0.002)
  fine <- .zero_information_metric(0.001)
  finer <- .zero_information_metric(0.0005)

  expect_lt(fine, coarse)
  expect_lt(finer, fine)
  # The metric is quadratic in a posterior discrepancy that is itself first
  # order in the step, so halving the step quarters the residual.
  expect_equal(coarse / fine, 4, tolerance = 0.15)
  expect_equal(fine / finer, 4, tolerance = 0.15)
})

# Metric domain ----------------------------------------------------------------

test_that("the metric is exactly zero for identical posteriors", {
  post <- list(mean = c(1, 2, 3), cov = c(0.5, 0.5, 0.5))
  expect_identical(aci_metric(post, post), c(0, 0, 0))
})

test_that("exact zeros are not counted as round-off clamps", {
  # Identical posteriors give exact zeros with nothing to clamp. Whether a
  # genuine clamp fires on a given trajectory depends on the platform's
  # rounding, so the counting itself is pinned here on constructed values
  # rather than on a run that may or may not produce one.
  post <- list(mean = c(1, 2, 3), cov = c(0.5, 0.5, 0.5))
  pair <- .aci_metric_pair(post, post)
  expect_identical(pair$value, c(0, 0, 0))
  expect_identical(pair$n_clamped, 0L)

  finished <- .aci_metric_finish(c(0, -1e-12, 2, 1e-3))
  expect_identical(finished$value, c(0, 0, 2, 1e-3))
  expect_identical(finished$n_clamped, 1L)

  # Anything more negative than round-off is a defect, not a clamp.
  expect_error(.aci_metric_finish(c(1, -1e-3)), "finite and non-negative")
})

test_that("the metric is accurate near a covariance ratio of one", {
  # The direct form of the dispersion term, 0.5 * (-log(r) + r - 1), subtracts
  # two nearly equal quantities and collapses to exactly zero here, losing the
  # value entirely. The cancellation-resistant form recovers it.
  post <- list(mean = c(1, 2, 3), cov = c(0.5, 0.5, 0.5))
  near <- list(mean = c(1, 2, 3), cov = c(0.5, 0.5, 0.5) * (1 + 1e-12))
  value <- aci_metric(post, near)

  expect_true(all(value >= 0))
  # Analytically the dispersion is delta^2 / 4 for a ratio of 1 + delta.
  expect_equal(value, rep(1e-24 / 4, 3L), tolerance = 1e-3)
})

test_that("the metric is non-negative across a full dyad trajectory", {
  model <- aci_dyad_model()
  fit <- aci(aci_simulate(model, n = 5000L, seed = 333)$x, model)
  expect_true(all(is.finite(fit$aci)))
  expect_true(all(fit$aci >= 0))
})
