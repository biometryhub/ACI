# Documented ENSO startup example --------------------------------------------
#
# The @examples block of aci_enso_model() is the only worked configuration in
# the package for a multi-state hidden partition, so its numbers are pinned
# here.  The record is aci_enso_model(hidden = c("u", "hW", "tau")) simulated
# at seed 12, t_end 4, dt 0.005, burn_in 0, and the prior is the componentwise
# Ornstein-Uhlenbeck stationary scale gyy_ii / (2 |Ly_ii|) read off the model's
# own coefficients at the first observation.  These values are tied to that
# seed/horizon/step and to the aci_code coefficient preset; if the preset is
# revised this file is the tripwire and the help text has to be re-measured.

.enso_example_setup <- function() {
  m <- aci_enso_model(hidden = c("u", "hW", "tau"))
  ob <- as_obs(simulate(m, seed = 12, t_end = 4, dt = 0.005, burn_in = 0))
  d0 <- function(f) diag(as.matrix(f(ob$t[1], ob$x[1, ])))
  gyy <- d0(function(t, x) tcrossprod(m$Sy1(t, x)) + tcrossprod(m$Sy2(t, x)))
  Ly <- d0(m$Ly)
  ini <- list(mean = m$meta$ic_default$y0, cov = diag(gyy / (2 * abs(Ly))))
  list(model = m, obs = ob, init = ini)
}


test_that("the bare constructor default is the two-state partition", {
  expect_identical(aci_enso_model()$meta$vars$hidden, c("hW", "tau"))
})


test_that("the documented ENSO configuration runs clean", {
  r <- .enso_example_setup()
  expect_equal(
    diag(r$init$cov), c(0.0032, 0.0008, 0.283697443577),
    tolerance = 1e-12, ignore_attr = TRUE
  )
  a <- expect_no_condition(aci(r$model, r$obs, init = r$init))
  expect_identical(a$meta$regularization$policy, "none")
  expect_false(a$meta$regularization$fired)
  expect_identical(a$meta$regularization$n_events, 0L)
  expect_equal(max(a$aci), 1.83086547999, tolerance = 1e-8)
})


test_that("the automatic prior is refused at the first update", {
  r <- .enso_example_setup()
  e <- tryCatch(
    suppressWarnings(aci(r$model, r$obs)),
    aci_error_covariance_not_spd = function(e) e
  )
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_identical(e$site, "filter_explicit")
  expect_identical(e$index, 2L)
  expect_equal(e$time, 0.005, tolerance = 1e-12)
  expect_equal(e$value, -7.149360359, tolerance = 1e-6)
})


test_that("refining the step runs the automatic-prior record instead", {
  r <- .enso_example_setup()
  a20 <- suppressWarnings(aci(r$model, r$obs, nsub = 20))
  expect_identical(a20$meta$regularization$n_events, 0L)
  expect_equal(max(a20$aci), 4.821106, tolerance = 1e-5)
  imp <- suppressWarnings(aci(r$model, r$obs, stepper = "implicit"))
  expect_identical(imp$meta$regularization$n_events, 0L)
  expect_equal(max(imp$aci), 4.821235, tolerance = 1e-5)
})
