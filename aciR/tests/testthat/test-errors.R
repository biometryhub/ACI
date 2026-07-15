# The error messages are part of the public interface: a user reads them far
# more often than they read the reference manual, and a package that rejects
# an input without saying why has only half-enforced its contract. These
# snapshots lock the wording so it cannot drift silently.

test_that("the model constructor's error messages are stable", {
  expect_snapshot(error = TRUE, {
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = -1
    )
  })
  expect_snapshot(error = TRUE, {
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, S_yoS_x = 1.5
    )
  })
  expect_snapshot(error = TRUE, {
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 0, S_yoS_y = 1
    )
  })
  expect_snapshot(error = TRUE, {
    aci_cgns_model(
      L_x = "nope", f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1
    )
  })
})

test_that("the observed-signal error messages are stable", {
  model <- aci_dyad_model()
  expect_snapshot(error = TRUE, aci(c(1, NA, 3), model))
  expect_snapshot(error = TRUE, aci(c(1, Inf, 3), model))
  expect_snapshot(error = TRUE, aci(1, model))
  expect_snapshot(error = TRUE, aci(matrix(1:10, nrow = 2L), model))
})

test_that("the coefficient-contract error messages are stable", {
  scalar_coef <- aci_cgns_model(
    L_x = function(x) 1, f_x = 0, L_y = -0.5, f_y = 0,
    S_xoS_x = 1, S_yoS_y = 1
  )
  expect_snapshot(error = TRUE, aci(c(1, 2, 3, 4), scalar_coef))
})

test_that("the components-schema error messages are stable", {
  model <- aci_dyad_model()
  x <- c(1, 2, 3, 4)
  comp <- aci_dyad_components(x, model$parameters)
  expect_snapshot(error = TRUE, aci_filter(x, comp[-1L], 0.001, 2, 0.1))

  asymmetric <- comp
  asymmetric$S_yoS_x <- 0.1
  asymmetric$S_xoS_y <- 0.2
  expect_snapshot(error = TRUE, aci_filter(x, asymmetric, 0.001, 2, 0.1))
})

test_that("the posterior-contract error messages are stable", {
  filt <- list(mean = c(1, 2, 3), cov = c(1, 1, 1))
  expect_snapshot(
    error = TRUE,
    aci_metric(filt, list(mean = c(1, 2), cov = c(1, 1)))
  )
  expect_snapshot(
    error = TRUE,
    aci_metric(filt, list(mean = c(1, 2, 3), cov = c(1, 0, 1)))
  )
})

test_that("the runtime covariance-guard message is stable", {
  model <- aci_dyad_model()
  x <- aci_simulate(model, n = 200L, seed = 1)$x
  expect_snapshot(error = TRUE, aci(x, model, dt = 1))
})

test_that("the time-grid error messages are stable", {
  model <- aci_dyad_model()
  x <- aci_simulate(model, n = 5L, seed = 1)$x
  expect_snapshot(error = TRUE, aci(x, model, time = c(0, 1, 2, 4, 5)))
  expect_snapshot(error = TRUE, aci(x, model, time = c(0, 1, 2, 2, 3)))
})
