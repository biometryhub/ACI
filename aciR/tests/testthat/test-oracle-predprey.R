# -- binding independent-oracle grade: the time-varying self-drift -------------
#
# The dyad and cross fixtures both pin the latent self-drift at a constant, so
# neither can grade a system whose latent damping moves. The predator-prey
# model is the reference implementation's own such system: in both causal
# directions the latent variable's damping is set by the observed state.
#
# Both directions are graded, and they are genuinely different questions rather
# than one question and its mirror -- the observed process differs, the latent
# process differs, and the metric is not symmetric between them.
#
# Runs from the installed fixtures and never skips.

.aci_pp_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

.aci_pp_parameters <- function() {
  list(
    alpha = 0.4, beta = 0.1, gamma = 1.1, delta = 0.4,
    sigma_x = 0.3, sigma_y = 0.3
  )
}

.aci_pp_signal <- function() {
  read.csv(.aci_pp_fixture("predprey_signal.csv"))
}

test_that("the predator-prey fixture exercises a moving self-drift", {
  # The point of this fixture is the coefficient it varies. If the latent
  # self-drift were ever constant here, the grade below would be a second copy
  # of the dyad grade and would not reach the path it exists to reach.
  signal <- .aci_pp_signal()
  p <- .aci_pp_parameters()

  for (direction in c("predator_to_prey", "prey_to_predator")) {
    observed <- if (direction == "predator_to_prey") {
      signal$prey
    } else {
      signal$predator
    }
    comp <- aci_predprey_components(observed, p, direction)
    expect_length(comp$L_y, length(observed))
    expect_gt(diff(range(comp$L_y)), 0.5)
    # The drift changes sign over the record, so the latent process is damped
    # at some times and driven at others -- a regime a constant cannot express.
    expect_lt(min(comp$L_y), 0)
    expect_gt(max(comp$L_y), 0)
  }
})

test_that("both predator-prey directions reproduce the MATLAB oracle to 1e-6", {
  signal <- .aci_pp_signal()
  p <- .aci_pp_parameters()
  dt <- 0.005

  for (direction in c("predator_to_prey", "prey_to_predator")) {
    observed <- if (direction == "predator_to_prey") {
      signal$prey
    } else {
      signal$predator
    }
    ref <- read.csv(
      .aci_pp_fixture(sprintf("predprey_reference_%s.csv", direction))
    )
    comp <- aci_predprey_components(observed, p, direction)
    filt <- aci_filter(observed, comp, dt, mu0 = 4, R0 = 0.1)
    smooth <- aci_smoother(observed, comp, dt, filt)
    metric <- aci_metric(filt, smooth)
    idx <- ref$index

    expect_equal(filt$mean[idx], ref$filter_mean, tolerance = 1e-6)
    expect_equal(filt$cov[idx], ref$filter_cov, tolerance = 1e-6)
    expect_equal(smooth$mean[idx], ref$smoother_mean, tolerance = 1e-6)
    expect_equal(smooth$cov[idx], ref$smoother_cov, tolerance = 1e-6)
    expect_equal(metric[idx], ref$ACI_metric, tolerance = 1e-6)

    max_abs_error <- max(
      abs(filt$mean[idx] - ref$filter_mean),
      abs(filt$cov[idx] - ref$filter_cov),
      abs(smooth$mean[idx] - ref$smoother_mean),
      abs(smooth$cov[idx] - ref$smoother_cov),
      abs(metric[idx] - ref$ACI_metric)
    )
    expect_lt(max_abs_error, 1e-6)
  }
})

test_that("the two directions are different questions", {
  # A symmetric answer would mean the metric was measuring the coupling rather
  # than its direction, which is the misreading the method most invites.
  signal <- .aci_pp_signal()
  p <- .aci_pp_parameters()
  dt <- 0.005

  metrics <- lapply(c("predator_to_prey", "prey_to_predator"), function(d) {
    observed <- if (d == "predator_to_prey") signal$prey else signal$predator
    comp <- aci_predprey_components(observed, p, d)
    filt <- aci_filter(observed, comp, dt, mu0 = 4, R0 = 0.1)
    aci_metric(filt, aci_smoother(observed, comp, dt, filt))
  })
  expect_gt(max(abs(metrics[[1L]] - metrics[[2L]])), 1e-3)
})

test_that("the model layer agrees with the hand-built components", {
  signal <- .aci_pp_signal()
  for (direction in c("predator_to_prey", "prey_to_predator")) {
    observed <- if (direction == "predator_to_prey") {
      signal$prey
    } else {
      signal$predator
    }
    model <- aci_predprey_model(direction)
    built <- aci_predprey_components(observed, model$parameters, direction)
    expect_equal(model$L_x(observed), built$L_x)
    expect_equal(model$L_y(observed), built$L_y)
    expect_equal(model$f_x(observed), built$f_x)
    expect_identical(model$S_xoS_x, built$S_xoS_x)
    expect_identical(model$S_yoS_y, built$S_yoS_y)
    # A state-dependent self-drift has no constant to report.
    expect_true(is.na(model$L_y_constant))
  }
})
