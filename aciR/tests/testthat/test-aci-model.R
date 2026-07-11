test_that("aci_cgns_model builds a well-formed model object", {
  model <- aci_cgns_model(
    L_x = function(x) 2 * x,
    f_x = function(x) 0.5 - 0.5 * x,
    L_y = -0.5,
    f_y = function(x) 1 - 2 * x^2,
    S_xoS_x = 0.25,
    S_yoS_y = 1
  )

  expect_s3_class(model, "aci_model")
  expect_true(is.function(model$L_x))
  expect_true(is.function(model$f_x))
  expect_true(is.function(model$f_y))
  expect_identical(model$L_y, -0.5)
  expect_identical(model$S_yoS_x, 0)
  expect_identical(model$S_xoS_y, 0)
})

test_that("aci_cgns_model accepts a numeric constant as a coefficient", {
  model <- aci_cgns_model(
    L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
    S_xoS_x = 1, S_yoS_y = 1
  )
  expect_equal(model$L_x(c(1, 2, 3)), c(1, 1, 1))
  expect_equal(model$f_x(numeric(0)), numeric(0))
})

test_that("aci_cgns_model rejects malformed input", {
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = c(1, 2), f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1
    ),
    "L_y"
  )
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 0, S_yoS_y = 1
    ),
    "positive"
  )
  expect_error(
    aci_cgns_model(
      L_x = "not a coefficient", f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1
    ),
    "L_x"
  )
})

test_that("aci_dyad_model matches the hand-built dyad components", {
  model <- aci_dyad_model()
  expect_s3_class(model, "aci_model")

  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  x <- c(0.3, 1.1, -0.4, 2.0)
  comp_ref <- aci_dyad_components(x, p)

  expect_equal(model$L_x(x), comp_ref$L_x)
  expect_equal(model$f_x(x), comp_ref$f_x)
  expect_equal(model$L_y, comp_ref$L_y)
  expect_equal(model$f_y(x), comp_ref$f_y)
  expect_equal(model$S_xoS_x, comp_ref$S_xoS_x)
  expect_equal(model$S_yoS_y, comp_ref$S_yoS_y)
  expect_equal(model$x0, p$F_x / p$d_x)
  expect_equal(model$y0, p$F_y / p$d_y)
})

test_that("aci_simulate returns a well-formed, reproducible path", {
  model <- aci_dyad_model()
  sim_a <- aci_simulate(model, n = 1000, seed = 333)
  sim_b <- aci_simulate(model, n = 1000, seed = 333)

  expect_s3_class(sim_a, "data.frame")
  expect_named(sim_a, c("t", "x", "y"))
  expect_equal(nrow(sim_a), 1000)
  expect_equal(sim_a$x[1], model$x0)
  expect_equal(sim_a$y[1], model$y0)
  expect_true(all(is.finite(sim_a$x)))
  expect_identical(sim_a, sim_b)
})

test_that("aci_simulate rejects malformed input and cross-noise", {
  model <- aci_dyad_model()
  expect_error(aci_simulate(model, n = 1), "at least two")
  expect_error(aci_simulate(model, n = 100, dt = -1), "positive")
  expect_error(aci_simulate(list(), n = 100), "aci_model")

  cross <- aci_cgns_model(
    L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
    S_xoS_x = 1, S_yoS_y = 1, S_yoS_x = 0.2
  )
  expect_error(aci_simulate(cross, n = 100), "independent noise")
})

test_that("aci runs end to end and returns a structured result", {
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 2000, seed = 333)
  fit <- aci(sim$x, model)

  expect_s3_class(fit, "aci")
  expect_named(
    fit,
    c("model", "t", "x", "filter", "smoother", "aci", "dt")
  )
  expect_length(fit$aci, length(sim$x))
  expect_length(fit$t, length(sim$x))
  expect_true(all(is.finite(fit$aci)))
  expect_true(all(fit$aci >= -1e-10))
  expect_true(all(fit$filter$cov > 0))
})

test_that("aci matches the low-level pipeline exactly", {
  model <- aci_dyad_model()
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  sim <- aci_simulate(model, n = 2000, seed = 1)
  x <- sim$x

  comp <- aci_dyad_components(x, p)
  filt <- aci_filter(x, comp, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
  smooth <- aci_smoother(x, comp, dt = 0.001, filt)
  metric <- aci_metric(filt, smooth)

  fit <- aci(x, model, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
  expect_identical(fit$filter$mean, filt$mean)
  expect_identical(fit$smoother$mean, smooth$mean)
  expect_identical(fit$aci, metric)
})

test_that("aci defaults the initial filtered mean to the model state", {
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 500, seed = 7)
  default_fit <- aci(sim$x, model)
  explicit_fit <- aci(sim$x, model, mu0 = model$y0)
  expect_identical(default_fit$filter$mean, explicit_fit$filter$mean)
})

test_that("aci rejects malformed input", {
  model <- aci_dyad_model()
  expect_error(aci(1, model), "at least two")
  expect_error(aci(c(1, 2, 3), list()), "aci_model")
  expect_error(aci(c(1, 2, 3), model, R0 = 0), "positive")
})

test_that("print methods return their object invisibly", {
  model <- aci_dyad_model()
  fit <- aci(aci_simulate(model, n = 500, seed = 1)$x, model)

  expect_output(print(model), "aci_model")
  expect_output(print(fit), "assimilative causal inference")
  expect_invisible(print(model))
  expect_invisible(print(fit))
})
