test_that("the CGNS core runs end to end and returns well-formed output", {
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  set.seed(333)
  sim <- aci_simulate_dyad(n = 3000, p = p)
  comp <- aci_dyad_components(sim$x, p)
  filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
  smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
  aci <- aci_metric(filt, smooth)

  expect_length(filt$mean, length(sim$x))
  expect_length(smooth$mean, length(sim$x))
  expect_length(aci, length(sim$x))
  expect_true(all(is.finite(filt$mean)))
  expect_true(all(is.finite(smooth$mean)))
  expect_true(all(filt$cov > 0))
  expect_true(all(smooth$cov > 0))
})

test_that("the causal-information metric is non-negative", {
  # The metric is a Kullback-Leibler divergence, hence non-negative up to
  # floating-point rounding at steps where the smoother and filter agree.
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  set.seed(333)
  sim <- aci_simulate_dyad(n = 3000, p = p)
  comp <- aci_dyad_components(sim$x, p)
  filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
  smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
  aci <- aci_metric(filt, smooth)

  expect_true(all(is.finite(aci)))
  expect_true(all(aci >= -1e-10))
})

test_that("aci_dyad_components has the expected shape", {
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  x <- c(1, 2, 3)
  comp <- aci_dyad_components(x, p)

  expect_named(
    comp,
    c(
      "L_x", "f_x", "L_y", "f_y",
      "S_xoS_x", "S_yoS_y", "S_yoS_x", "S_xoS_y"
    )
  )
  expect_length(comp$L_x, length(x))
  expect_length(comp$L_y, 1L)
  expect_identical(comp$S_yoS_x, 0)
  expect_identical(comp$S_xoS_y, 0)
})
