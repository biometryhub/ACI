# Observation-noise scale and time resolution --------------------------------
# These tests describe the current limitation; they do not grade its inaccurate
# outputs as reference results. A future diagnostic must retain the valid cases.

.observation_scale_case <- function(noise = 1e-8, coupling = 1, k = 1L,
                                    drift = -0.5, process = 0.5,
                                    constant = FALSE, prior = 1,
                                    observed_unit = 1, hidden_unit = 1,
                                    time_unit = 1) {
  force(noise); force(coupling); force(k); force(drift); force(process)
  force(constant); force(observed_unit); force(hidden_unit); force(time_unit)
  m <- aci_model(
    Lx = function(t, x) diag(observed_unit * coupling /
                              (hidden_unit * time_unit), k),
    fx = function(t, x) rep(0, k),
    Ly = function(t, x) diag(drift / time_unit, k),
    fy = function(t, x) rep(0, k),
    Sx1 = function(t, x) diag(observed_unit / sqrt(time_unit) *
      if (constant || abs(t / time_unit - 0.01) < 1e-12) noise else 0.5, k),
    Sy2 = function(t, x) diag(hidden_unit * process / sqrt(time_unit), k),
    k = k, l = k
  )
  t <- c(0, 0.01, 0.02, 0.03)
  list(model = m,
       obs = observed_trajectory(time_unit * t,
         observed_unit * matrix(rep(t, k), ncol = k)),
       init = list(mean = rep(0, k), cov = diag(hidden_unit^2 * prior, k)))
}

.observation_scale_filter <- function(case, stepper = "implicit", nsub = 1L) {
  aci_filter(case$model, case$obs, init = case$init,
             stepper = stepper, nsub = nsub, regularize = "none")
}

test_that("near-zero noise can return a small positive covariance silently", {
  for (k in c(1L, 2L)) {
    case <- .observation_scale_case(k = k)
    f <- expect_no_condition(.observation_scale_filter(case))
    expect_true(all(is.finite(f$mean)))
    expect_true(all(is.finite(f$cov)))
    expect_true(is.finite(f$meta$loglik))
    expect_true(all(diag(matrix(f$cov[, , 3], nrow = k)) > 0))
    expect_lt(f$cov[1, 1, 3] / f$cov[1, 1, 2], 1e-12)
    expect_identical(f$meta$regularization$n_events, 0L)
    expect_error(.observation_scale_filter(case, "explicit"),
                 class = "aci_error_covariance_not_spd")
  }
})

test_that("near-zero noise need not make explicit integration fail", {
  case <- .observation_scale_case(coupling = 1e-10)
  for (stepper in c("explicit", "implicit")) {
    f <- expect_no_condition(.observation_scale_filter(case, stepper))
    expect_true(all(is.finite(f$mean)))
    expect_gt(min(f$cov), 0.9)
    expect_identical(f$meta$regularization$n_events, 0L)
  }
})

test_that("a precise observation of a static hidden state is valid", {
  case <- .observation_scale_case(constant = TRUE, drift = 0, process = 0)
  f <- expect_no_condition(.observation_scale_filter(case))
  # With no hidden dynamics, R' = -R^2/G and R(0) = 1, so
  # R(t) = 1/(1+t/G). The large first reduction is exact, not a defect.
  exact <- 1 / (1 + case$obs$t / 1e-16)
  expect_lt(max(abs(as.numeric(f$cov) / exact - 1)), 1e-12)
  expect_gt(f$cov[1, 1, 1] / f$cov[1, 1, 2], 1e12)
  expect_identical(f$meta$regularization$n_events, 0L)
})

test_that("consistent state and time units preserve the hidden posterior", {
  base <- .observation_scale_filter(
    .observation_scale_case(noise = 1, constant = TRUE))
  for (units in list(c(1e-9, 1, 1), c(1e-9, 1e3, 60))) {
    changed <- expect_no_condition(.observation_scale_filter(
      .observation_scale_case(noise = 1, constant = TRUE,
        observed_unit = units[1], hidden_unit = units[2], time_unit = units[3])))
    expect_equal(changed$mean / units[2], base$mean, tolerance = 1e-12)
    expect_equal(changed$cov / units[2]^2, base$cov, tolerance = 1e-12)
    expect_identical(changed$meta$regularization$n_events, 0L)
  }
})

test_that("a finite positive return does not certify time resolution", {
  # R' = -R + 0.25 - R^2/1e-6. Starting at the positive stationary
  # root makes the exact covariance constant, independently of the record.
  exact <- 0.25 / (sqrt(0.25 + 0.25 / 1e-6) + 0.5)
  case <- .observation_scale_case(noise = 1e-3, constant = TRUE, prior = exact)
  coarse <- expect_no_condition(.observation_scale_filter(case))
  fine <- expect_no_condition(.observation_scale_filter(case, nsub = 1000L))
  coarse_error <- abs(coarse$cov[1, 1, 3] / exact - 1)
  fine_error <- abs(fine$cov[1, 1, 3] / exact - 1)
  expect_gt(coarse_error, 0.5)
  expect_lt(fine_error, 0.01)
  expect_lt(fine_error, coarse_error)
  expect_identical(coarse$meta$regularization$n_events, 0L)
  expect_identical(fine$meta$regularization$n_events, 0L)
})
