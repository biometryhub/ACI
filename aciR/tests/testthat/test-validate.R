# The contract boundary: every malformed input must fail before the numerical
# recursion starts, and must say what was wrong. These tests exercise the
# validators through the public functions rather than directly, because the
# contract they protect is a public one.

.aci_test_signal <- function(n = 50L) {
  sim <- aci_simulate(aci_dyad_model(), n = n, seed = 1)
  sim$x
}

.aci_test_comp <- function(x) {
  aci_dyad_components(x, aci_dyad_model()$parameters)
}

# -- noise covariance admissibility (T1, T2, T3) ------------------------------

test_that("the constructor rejects a negative latent-noise covariance", {
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = -1
    ),
    "S_yoS_y.*non-negative"
  )
})

test_that("the constructor rejects a non-positive observation covariance", {
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 0, S_yoS_y = 1
    ),
    "positive"
  )
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = -1, S_yoS_y = 1
    ),
    "positive"
  )
})

test_that("the constructor rejects an indefinite joint noise covariance", {
  # |S_yoS_x| > sqrt(S_xoS_x * S_yoS_y) = 1 is not a covariance matrix.
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, S_yoS_x = 1.5
    ),
    "positive semidefinite"
  )
})

test_that("a singular joint noise covariance is admissible", {
  # Perfectly correlated noise: the determinant is exactly zero. The model is
  # mathematically admissible, so construction must succeed; whether a given
  # run stays numerically well posed is the per-step guard's business.
  expect_s3_class(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, S_yoS_x = 1
    ),
    "aci_model"
  )
  # A zero latent-noise covariance is the other singular boundary.
  expect_s3_class(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 0
    ),
    "aci_model"
  )
})

# -- label and parameters (T4) ------------------------------------------------

test_that("the constructor validates label and parameters", {
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, label = NA_character_
    ),
    "non-missing character"
  )
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, label = c("a", "b")
    ),
    "single non-missing character"
  )
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, parameters = list(1, 2)
    ),
    "named list of finite numeric scalars"
  )
  expect_error(
    aci_cgns_model(
      L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1, parameters = list(a = "x")
    ),
    "named list of finite numeric scalars"
  )
})

# -- coefficient contract (T5, T6) --------------------------------------------

test_that("a coefficient returning a scalar is rejected, never recycled", {
  # The old failure: this constructed, then produced NA filter entries.
  model <- aci_cgns_model(
    L_x = function(x) 1, f_x = 0, L_y = -0.5, f_y = 0,
    S_xoS_x = 1, S_yoS_y = 1
  )
  expect_error(aci(.aci_test_signal(), model), "L_x.*expected length 50, got 1")
  expect_error(aci(.aci_test_signal(), model), "not recycled")
})

test_that("a coefficient returning the wrong type is rejected", {
  model <- aci_cgns_model(
    L_x = function(x) rep("a", length(x)), f_x = 0, L_y = -0.5, f_y = 0,
    S_xoS_x = 1, S_yoS_y = 1
  )
  expect_error(aci(.aci_test_signal(), model), "L_x.*must return a numeric")
})

test_that("a coefficient returning a non-finite value is rejected", {
  model_na <- aci_cgns_model(
    L_x = function(x) rep(NA_real_, length(x)), f_x = 0, L_y = -0.5, f_y = 0,
    S_xoS_x = 1, S_yoS_y = 1
  )
  expect_error(aci(.aci_test_signal(), model_na), "L_x.*NA at index 1")

  model_inf <- aci_cgns_model(
    L_x = 1, f_x = 0, L_y = -0.5,
    f_y = function(x) c(Inf, rep(0, length(x) - 1L)),
    S_xoS_x = 1, S_yoS_y = 1
  )
  expect_error(
    aci(.aci_test_signal(), model_inf),
    "f_y.*infinite value at index 1"
  )
})

test_that("a coefficient that is neither function nor scalar is rejected", {
  expect_error(
    aci_cgns_model(
      L_x = "not a coefficient", f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1
    ),
    "L_x"
  )
  expect_error(
    aci_cgns_model(
      L_x = c(1, 2), f_x = 0, L_y = -0.5, f_y = 0,
      S_xoS_x = 1, S_yoS_y = 1
    ),
    "L_x"
  )
})

# -- observed signal contract (T7) --------------------------------------------

test_that("a signal with missing or infinite values is rejected", {
  model <- aci_dyad_model()
  # The old failure: this returned an aci object full of NA metric values.
  expect_error(aci(c(1, NA, 3), model), "complete and finite")
  expect_error(aci(c(1, NA, 3), model), "index 2")
  expect_error(aci(c(1, NaN, 3), model), "NaN at index 2")
  expect_error(aci(c(1, Inf, 3), model), "infinite value at index 2")
})

test_that("a signal of the wrong type or shape is rejected", {
  model <- aci_dyad_model()
  expect_error(aci(matrix(1:10, nrow = 2L), model), "plain numeric vector")
  expect_error(aci(letters[1:5], model), "plain numeric vector")
  expect_error(aci(1, model), "at least two")
  expect_error(aci(numeric(0), model), "at least two")
})

test_that("the low-level entry points enforce the same signal contract", {
  # The old failure: aci_filter(1, ...) reached an opaque replacement error
  # from a descending loop range rather than failing at the boundary.
  comp <- .aci_test_comp(.aci_test_signal())
  expect_error(aci_filter(1, comp, 0.001, 2, 0.1), "at least two")
  expect_error(aci_smoother(1, comp, 0.001, list(mean = 1, cov = 1)),
    "at least two"
  )
  expect_error(aci_filter(c(1, NA, 3), comp, 0.001, 2, 0.1),
    "complete and finite"
  )
})

test_that("dt, mu0 and R0 are validated at every entry point", {
  model <- aci_dyad_model()
  x <- .aci_test_signal()
  comp <- .aci_test_comp(x)
  expect_error(aci(x, model, dt = -1), "positive")
  expect_error(aci(x, model, dt = 0), "positive")
  expect_error(aci(x, model, R0 = 0), "positive")
  expect_error(aci(x, model, mu0 = NA_real_), "single finite numeric")
  expect_error(aci_filter(x, comp, dt = -1, 2, 0.1), "positive")
  expect_error(aci_filter(x, comp, dt = 0.001, mu0 = 2, R0 = -1), "positive")
})

# -- components schema (T8) ---------------------------------------------------

test_that("a components list missing an entry is rejected", {
  x <- .aci_test_signal()
  comp <- .aci_test_comp(x)
  expect_error(aci_filter(x, comp[-1L], 0.001, 2, 0.1), "missing the component")
  expect_error(aci_filter(x, comp[-1L], 0.001, 2, 0.1), "aci_components")
  expect_error(aci_filter(x, "not a list", 0.001, 2, 0.1), "components list")
})

test_that("a components list with a misshapen coefficient is rejected", {
  x <- .aci_test_signal()
  comp <- .aci_test_comp(x)
  short <- comp
  short$L_x <- short$L_x[1:3]
  expect_error(aci_filter(x, short, 0.001, 2, 0.1), "one value per observation")

  bad <- comp
  bad$f_x[5L] <- NA_real_
  expect_error(aci_filter(x, bad, 0.001, 2, 0.1), "f_x.*NA at index 5")
})

test_that("an asymmetric noise cross-covariance is rejected", {
  x <- .aci_test_signal()
  comp <- .aci_test_comp(x)
  comp$S_yoS_x <- 0.1
  comp$S_xoS_y <- 0.2
  expect_error(aci_filter(x, comp, 0.001, 2, 0.1), "symmetric")
})

test_that("an inadmissible components covariance is rejected", {
  x <- .aci_test_signal()
  comp <- .aci_test_comp(x)
  comp$S_yoS_y <- -1
  expect_error(aci_filter(x, comp, 0.001, 2, 0.1), "non-negative")
})

# -- posterior contract (T9) --------------------------------------------------

test_that("a malformed posterior is rejected before the metric", {
  filt <- list(mean = c(1, 2, 3), cov = c(1, 1, 1))
  expect_error(
    aci_metric(filt, list(mean = c(1, 2), cov = c(1, 1))),
    "must cover 3 time steps; it covers 2"
  )
  expect_error(
    aci_metric(filt, list(mean = c(1, 2, 3), cov = c(1, 0, 1))),
    "finite and strictly positive"
  )
  expect_error(
    aci_metric(list(mean = c(1, 2, 3), cov = c(1, -1, 1)), filt),
    "finite and strictly positive"
  )
  expect_error(aci_metric(filt, list(mu = 1)), "list with numeric")
  expect_error(
    aci_metric(filt, list(mean = c(1, 2, 3), cov = c(1, 1))),
    "equal length"
  )
})

test_that("a posterior with a non-finite mean is rejected", {
  filt <- list(mean = c(1, NA, 3), cov = c(1, 1, 1))
  smooth <- list(mean = c(1, 2, 3), cov = c(1, 1, 1))
  expect_error(aci_metric(filt, smooth), "filt\\$mean.*NA at index 2")
  expect_error(
    aci_metric(list(mean = c(1, Inf, 3), cov = c(1, 1, 1)), smooth),
    "infinite value at index 2"
  )
})

test_that("a non-finite time value is rejected", {
  model <- aci_dyad_model()
  x <- .aci_test_signal(5L)
  expect_error(aci(x, model, time = c(0, 1, NA, 3, 4)), "NA at index 3")
  expect_error(aci(x, model, time = c(0, 1, Inf, 3, 4)), "infinite value")
})

test_that("an error names a function it was handed by mistake", {
  # The describing helper has a branch for functions, because handing a
  # coefficient function where a components list belongs is an easy slip and
  # "a function" is more use than "closure of length 1".
  x <- .aci_test_signal()
  expect_error(
    aci_filter(x, function(z) z, 0.001, 2, 0.1),
    "it is a function"
  )
})

test_that("the smoother rejects a filter posterior of the wrong length", {
  x <- .aci_test_signal()
  comp <- .aci_test_comp(x)
  expect_error(
    aci_smoother(x, comp, 0.001, list(mean = c(1, 2), cov = c(1, 1))),
    "must cover 50 time steps; it covers 2"
  )
})

# -- runtime covariance guard (T10) -------------------------------------------

test_that("the metric refuses a posterior pair whose ratio overflows", {
  # Both covariances are positive and finite and both means are finite, so the
  # posterior contract is satisfied; but the ratio overflows to Inf and the
  # dispersion term evaluates Inf - Inf. The result would be NaN, and a
  # comparison guard alone would not catch it, because NaN <= x is NA rather
  # than TRUE. The metric must refuse rather than return a missing value.
  filt <- list(mean = 0, cov = 1e-300)
  smooth <- list(mean = 0, cov = 1e300)
  expect_error(aci_metric(filt, smooth), "finite and non-negative")
  expect_error(aci_metric(filt, smooth), "overflows")
})

test_that("the smoother covariance guard fires on its own", {
  # The filter is well posed throughout -- its covariance climbs from a tiny
  # initial value toward its stationary point -- but the smoother divides by
  # that tiny initial covariance on its final backward step, which drives its
  # own covariance negative. This reaches the smoother's guard without the
  # filter's guard firing first.
  n <- 200L
  x <- seq(1, 1.2, length.out = n)
  comp <- list(
    L_x = rep(0, n), f_x = rep(0.2, n), L_y = -0.5, f_y = rep(0.3, n),
    S_xoS_x = 0.5, S_yoS_y = 1, S_yoS_x = 0, S_xoS_y = 0
  )
  filt <- aci_filter(x, comp, dt = 0.001, mu0 = 0.6, R0 = 1e-6)
  expect_true(all(filt$cov > 0))

  expect_error(aci_smoother(x, comp, dt = 0.001, filt), "smoother covariance")
  expect_error(aci_smoother(x, comp, dt = 0.001, filt), "index 1 ")
})

test_that("an unstable discretisation fails at the first bad step", {
  # A perfectly admissible model, integrated with a step far too large for it:
  # the explicit Euler update drives the filtered covariance non-positive.
  # The old behaviour returned NaN metric values; the contract is to stop and
  # name the step.
  model <- aci_dyad_model()
  x <- .aci_test_signal(200L)
  expect_error(aci(x, model, dt = 1), "filter covariance")
  expect_error(aci(x, model, dt = 1), "index \\d+")
  expect_error(aci(x, model, dt = 1), "Reduce `dt`")
})

# -- time grid (T17) ----------------------------------------------------------

test_that("an irregular time grid is rejected", {
  model <- aci_dyad_model()
  x <- .aci_test_signal(10L)
  irregular <- c(0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008,
    0.009, 0.01
  )
  expect_error(aci(x, model, time = irregular), "equally spaced")
  expect_error(aci(x, model, time = rev(seq(0, by = 0.001, length.out = 10L))),
    "strictly increasing"
  )
  expect_error(aci(x, model, time = seq(0, by = 0.001, length.out = 5L)),
    "one value per observation"
  )
  expect_error(aci(x, model, time = "a"), "plain numeric vector")
})

test_that("a time grid that contradicts dt is rejected", {
  model <- aci_dyad_model()
  x <- .aci_test_signal(10L)
  expect_error(
    aci(x, model, dt = 0.5, time = seq(0, by = 0.001, length.out = 10L)),
    "disagree"
  )
})

test_that("a regular time grid derives the step and agrees with dt", {
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 500L, seed = 3)
  by_dt <- aci(sim$x, model, dt = 0.001)
  by_time <- aci(sim$x, model, time = sim$t)
  expect_equal(by_time$dt, 0.001)
  expect_identical(by_time$aci, by_dt$aci)
})
