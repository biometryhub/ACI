# Non-finite recursions are refusals, not results -----------------------------
#
# Overflow in a mean recursion or a likelihood accumulator is a classed
# refusal naming the quantity, the grid index and the time, not a returned
# Gaussian path with non-finite means. Finite records are untouched.

quietly_nf <- function(expr) withCallingHandlers(
  expr, warning = function(w) invokeRestart("muffleWarning"))

overflow_scalar_model <- function(fy_value) aci_model(
  Lx = function(t, x) matrix(0, 1, 1),
  fx = function(t, x) 0,
  Ly = function(t, x) matrix(0, 1, 1),
  fy = function(t, x) fy_value,
  Sx1 = function(t, x) matrix(1, 1, 1),
  Sy2 = function(t, x) matrix(sqrt(0.5), 1, 1),
  k = 1, l = 1)

overflow_record <- function() observed_trajectory(0:5, matrix(0, 6, 1))
overflow_init <- list(mean = 0, cov = matrix(1, 1, 1))


test_that("a diverging scalar filter mean is refused with its index and time", {
  m <- overflow_scalar_model(6e307)
  ob <- overflow_record()
  for (stepper in c("explicit", "implicit"))
    for (nsub in list(1L, 4L))
      for (loglik in c(TRUE, FALSE))
        for (regularize in c("none", "floor")) {
          info <- paste(stepper, nsub, loglik, regularize)
          e <- expect_error(
            quietly_nf(aci_filter(m, ob, init = overflow_init,
                                  stepper = stepper, nsub = nsub,
                                  loglik = loglik, regularize = regularize)),
            class = "aci_error_nonfinite", info = info)
          expect_identical(e$quantity, "filter mean", info = info)
          expect_identical(e$index, 4L, info = info)
          expect_identical(e$time, 3, info = info)
          ## Invalid arithmetic is distinguished from a covariance that left
          ## the positive-definite cone.
          expect_false(inherits(e, "aci_error_covariance_not_spd"),
                       info = info)
          expect_false(inherits(e, "aci_error_spd"), info = info)
          expect_match(conditionMessage(e), "not finite")
        }
})


test_that("the matrix analogue is refused at the same index", {
  m <- aci_model(
    Lx = function(t, x) matrix(0, 2, 2),
    fx = function(t, x) c(0, 0),
    Ly = function(t, x) matrix(0, 2, 2),
    fy = function(t, x) c(6e307, 0),
    Sx1 = function(t, x) diag(2),
    Sy2 = function(t, x) diag(2),
    k = 2, l = 2)
  ob <- observed_trajectory(0:5, matrix(0, 6, 2))
  ini <- list(mean = c(0, 0), cov = diag(2))
  for (stepper in c("explicit", "implicit"))
    for (nsub in list(1L, 4L)) {
      info <- paste(stepper, nsub)
      e <- expect_error(
        quietly_nf(aci_filter(m, ob, init = ini, stepper = stepper,
                              nsub = nsub)),
        class = "aci_error_nonfinite", info = info)
      expect_identical(e$quantity, "filter mean", info = info)
      expect_identical(e$index, 4L, info = info)
      expect_identical(e$time, 3, info = info)
    }
})


test_that("higher-level consumers refuse rather than return", {
  m <- overflow_scalar_model(6e307)
  ob <- overflow_record()
  ini <- overflow_init
  calls <- list(
    aci = function() aci(m, ob, init = ini),
    aci_smoother = function() aci_smoother(m, ob, init = ini),
    online_lag2 = function() aci_online(m, ob, lag = 2, init = ini),
    online_full = function() aci_online(m, ob, lag = Inf, init = ini),
    lag_table = function() lag_table(m, ob, mode = "forward", init = ini)
  )
  for (nm in names(calls)) {
    e <- expect_error(quietly_nf(calls[[nm]]()), class = "aci_error_nonfinite",
                      info = nm)
    expect_identical(e$index, 4L, info = nm)
  }
})


test_that("the smoother mean carries the same invariant", {
  ## Reached through the kernel with a supplied filter, because on the public
  ## route the filter guard fires first on any model whose drift is large
  ## enough to overflow the backward recursion. The supplied filter variance
  ## keeps the backward Riccati step inside the positive-definite cone, so the
  ## only invariant this can violate is the one on the mean.
  m <- overflow_scalar_model(6e307)
  ob <- observed_trajectory(0:8, matrix(0, 9, 1))
  bundle <- .compile_cgns_run(m, ob)
  n1 <- bundle$N1
  e <- expect_error(
    .cgns_smoother_scalar_kernel(bundle, rep(0, n1), rep(2, n1)),
    class = "aci_error_nonfinite")
  expect_identical(e$quantity, "smoother mean")
  expect_false(inherits(e, "aci_error_covariance_not_spd"))
  expect_true(e$index >= 1L && e$index <= n1)
})


test_that("a non-finite likelihood accumulator is refused", {
  ## The mean guard makes this a post-condition on the public routes, so it is
  ## reached through the kernel with moments the filter would never store.
  m <- aci_dyad_model()
  ob <- as_obs(simulate(m, seed = 1, t_end = 0.05, dt = 0.01))
  bundle <- .compile_cgns_run(m, ob)
  n1 <- bundle$N1
  e <- expect_error(
    .cgns_likelihood_scalar_kernel(bundle, rep(Inf, n1), rep(1, n1)),
    class = "aci_error_nonfinite")
  expect_identical(e$quantity, "predictive log-likelihood")
  expect_identical(e$index, 1L)
  expect_identical(e$time, bundle$t[1L])
})


test_that(paste0("a finite record keeps its moments exactly, with and ",
                 "without the likelihood"), {
  m <- overflow_scalar_model(1e307)
  ob <- overflow_record()
  expected_mean <- c(0, 1e307, 2e307, 3e307, 4e307, 5e307)
  expected_cov <- c(1, 1.5, 2, 2.5, 3, 3.5)
  for (loglik in c(TRUE, FALSE)) {
    f <- quietly_nf(aci_filter(m, ob, init = overflow_init, loglik = loglik))
    expect_identical(as.numeric(f$mean), expected_mean, info = loglik)
    expect_identical(as.numeric(f$cov), expected_cov, info = loglik)
    if (loglik) expect_true(is.finite(f$meta$loglik)) else
      expect_null(f$meta$loglik)
    expect_identical(f$meta$regularization$n_events, 0L, info = loglik)
  }
  ## The sum of those means is 1.5e308, within a factor of 1.2 of DBL_MAX, so
  ## a guard written as is.finite(sum(.)) would refuse this valid record.
  expect_false(is.finite(sum(expected_mean * 2)))
  expect_true(all(is.finite(expected_mean)))
})


test_that("the library baselines are unchanged", {
  m <- aci_dyad_model()
  ob <- as_obs(simulate(m, seed = 1, t_end = 2, dt = 0.01))
  a <- quietly_nf(aci(m, ob, init = list(mean = 0, cov = matrix(1, 1, 1)),
                      keep = "paths"))
  expect_equal(max(a$aci), 3.1797417591261574, tolerance = 1e-12)
  expect_equal(a$paths$filter$meta$loglik, 323.32851067503259,
               tolerance = 1e-12)
  e <- aci_enso_model(hidden = c("u", "hW", "tau"))
  obe <- as_obs(simulate(e, seed = 12, t_end = 4, dt = 0.005, burn_in = 0))
  fe <- quietly_nf(aci_filter(e, obe, stepper = "implicit",
                              init = list(mean = rep(0, 3),
                                          cov = diag(0.5, 3))))
  expect_equal(fe$meta$loglik, 9742.8037719546464, tolerance = 1e-12)
  expect_true(all(is.finite(fe$mean)))
})
