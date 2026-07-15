test_that("the CGNS core runs end to end and returns well-formed output", {
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 3000L, seed = 333)
  comp <- aci_dyad_components(sim$x, model$parameters)
  filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
  smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
  metric <- aci_metric(filt, smooth)

  expect_length(filt$mean, length(sim$x))
  expect_length(smooth$mean, length(sim$x))
  expect_length(metric, length(sim$x))
  expect_true(all(is.finite(filt$mean)))
  expect_true(all(is.finite(smooth$mean)))
  expect_true(all(filt$cov > 0))
  expect_true(all(smooth$cov > 0))
})

test_that("aci_dyad_components has the expected shape", {
  model <- aci_dyad_model()
  x <- c(1, 2, 3)
  comp <- aci_dyad_components(x, model$parameters)

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

test_that("aci_dyad_components validates its parameter list", {
  x <- c(1, 2, 3)
  p <- aci_dyad_model()$parameters
  expect_error(aci_dyad_components(x, p[-1L]), "missing the parameter")
  expect_error(aci_dyad_components(x, "not a list"), "named list")
  p$sigma_x <- 0
  expect_error(aci_dyad_components(x, p), "sigma_x.*must be non-zero")
  p$sigma_x <- NA_real_
  expect_error(aci_dyad_components(x, p), "single finite numeric")
})

test_that("the components list of the dyad model matches its constructor", {
  # The two routes to the same system must agree exactly: one implementation
  # of the equations, reachable two ways.
  model <- aci_dyad_model()
  x <- c(0.3, 1.1, -0.4, 2.0)
  comp <- aci_dyad_components(x, model$parameters)

  expect_equal(model$L_x(x), comp$L_x)
  expect_equal(model$f_x(x), comp$f_x)
  expect_equal(model$L_y, comp$L_y)
  expect_equal(model$f_y(x), comp$f_y)
  expect_equal(model$S_xoS_x, comp$S_xoS_x)
  expect_equal(model$S_yoS_y, comp$S_yoS_y)
  expect_equal(model$S_yoS_x, comp$S_yoS_x)
  expect_equal(model$S_xoS_y, comp$S_xoS_y)
})
