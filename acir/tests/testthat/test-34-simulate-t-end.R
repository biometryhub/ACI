# The horizon argument is `t_end`; the retired name `T` is accepted with a
# typed warning until 0.2.0, gives the same draws, and is refused beside
# `t_end`.

test_that("t_end and the retired T give identical draws, T with a warning", {
  m <- aci_dyad_model()
  new <- simulate(m, seed = 3, t_end = 0.5, dt = 0.01)
  expect_warning(old <- simulate(m, seed = 3, T = 0.5, dt = 0.01),
                 class = "aci_warning_deprecated")
  expect_identical(old, new)
})

test_that("passing both names is an error, and an unknown name still is", {
  m <- aci_dyad_model()
  expect_error(simulate(m, seed = 3, t_end = 0.5, T = 0.5, dt = 0.01),
               class = "aci_error_dims")
  expect_error(simulate(m, seed = 3, t_end = 0.5, dt = 0.01, horizon = 1),
               class = "aci_error_dims")
})

test_that("the named verb passes t_end through", {
  m <- aci_dyad_model()
  expect_identical(aci_simulate(m, seed = 5, t_end = 0.2, dt = 0.01),
                   simulate(m, seed = 5, t_end = 0.2, dt = 0.01))
})
