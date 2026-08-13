# -- independent-oracle grade: the vector core with correlated noise -----------
#
# Read the scope of this grade before citing it.
#
# oracle/aci_oracle_mv.m has NO upstream counterpart. Every scalar model in the
# reference implementation sets the noise cross-covariance to zero, and all
# five ENSO scripts state that the cross-interaction terms are absent. So the
# matrix-valued cross-noise terms are graded here against a SECOND INDEPENDENT
# IMPLEMENTATION of the published equations, in a different language and
# runtime. Agreement refutes an R-side transcription error; it is not an
# authors'-reference grounding.
#
# The primary grounding for those terms is analytic and lives in
# test-identities-mv.R: block-diagonal algebraic Riccati fixed points solved in
# closed form, which depend on neither implementation. Neither test alone is
# sufficient; together they cover the path.

.aci_mv_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

.aci_mv_noise <- function() {
  s <- matrix(
    c(0.60, 0.10, 0.25, 0.05,
      0.20, 0.50, 0.10, 0.30,
      0.15, 0.05, 0.70, 0.10,
      0.05, 0.10, 0.15, 0.55),
    nrow = 4L, ncol = 4L
  )
  s_x <- s[1:2, ]
  s_y <- s[3:4, ]
  list(
    S_xoS_x = s_x %*% t(s_x),
    S_yoS_y = s_y %*% t(s_y),
    S_yoS_x = s_y %*% t(s_x)
  )
}

.aci_mv_components <- function(signal) {
  n <- nrow(signal)
  x <- rbind(signal$x1, signal$x2)
  comp <- .aci_mv_noise()
  L_x <- array(0, c(2L, 2L, n))
  L_y <- array(0, c(2L, 2L, n))
  f_x <- matrix(0, 2L, n)
  f_y <- matrix(0, 2L, n)
  for (j in seq_len(n)) {
    x1 <- x[1L, j]
    x2 <- x[2L, j]
    tt <- signal$t[j]
    L_x[, , j] <- matrix(
      c(0.8 + 0.3 * x1, 0.1 * x2, 0.2 * sin(tt), 0.6 - 0.2 * x1), 2L, 2L
    )
    L_y[, , j] <- matrix(
      c(-1.2 + 0.1 * x1, 0.2, 0.3, -0.9 - 0.1 * x2), 2L, 2L
    )
    f_x[, j] <- c(0.4 - 0.5 * x1, -0.3 * x2 + 0.2)
    f_y[, j] <- c(0.5 - 0.2 * x1^2, 0.1 - 0.15 * x2)
  }
  comp$L_x <- L_x
  comp$L_y <- L_y
  comp$f_x <- f_x
  comp$f_y <- f_y
  list(x = x, comp = comp)
}

test_that("the fixture actually carries a non-zero matrix cross-covariance", {
  # The terms this fixture exists to grade are exactly the ones every upstream
  # model annihilates. If the cross-block were ever zero here, every comparison
  # below would still pass while grading nothing -- which is the failure this
  # package was built to make impossible.
  noise <- .aci_mv_noise()
  expect_gt(max(abs(noise$S_yoS_x)), 0.1)
  expect_gt(sum(noise$S_yoS_x^2), 0.1)
  # Off-diagonal structure too, or the Cholesky paths reduce to scalar ones.
  expect_gt(abs(noise$S_yoS_y[1L, 2L]), 0.05)
  expect_gt(abs(noise$S_xoS_x[1L, 2L]), 0.05)
})

test_that("the vector core reproduces the MATLAB oracle to 1e-6", {
  signal <- read.csv(.aci_mv_fixture("mv_signal.csv"))
  ref <- read.csv(.aci_mv_fixture("mv_reference.csv"))
  pieces <- .aci_mv_components(signal)
  dt <- 0.002

  filt <- aci_filter(
    pieces$x, pieces$comp, dt, mu0 = c(0.8, 0.2), R0 = 0.2 * diag(2L)
  )
  smooth <- aci_smoother(pieces$x, pieces$comp, dt, filt)
  metric <- aci_metric(filt, smooth)
  idx <- ref$index

  expect_lt(max(abs(filt$mean[1L, idx] - ref$fm1)), 1e-6)
  expect_lt(max(abs(filt$mean[2L, idx] - ref$fm2)), 1e-6)
  expect_lt(max(abs(filt$cov[1L, 1L, idx] - ref$fc11)), 1e-6)
  expect_lt(max(abs(filt$cov[1L, 2L, idx] - ref$fc12)), 1e-6)
  expect_lt(max(abs(filt$cov[2L, 2L, idx] - ref$fc22)), 1e-6)
  expect_lt(max(abs(smooth$mean[1L, idx] - ref$sm1)), 1e-6)
  expect_lt(max(abs(smooth$mean[2L, idx] - ref$sm2)), 1e-6)
  expect_lt(max(abs(smooth$cov[1L, 1L, idx] - ref$sc11)), 1e-6)
  expect_lt(max(abs(smooth$cov[1L, 2L, idx] - ref$sc12)), 1e-6)
  expect_lt(max(abs(smooth$cov[2L, 2L, idx] - ref$sc22)), 1e-6)
  expect_lt(max(abs(metric[idx] - ref$ACI_metric)), 1e-6)

  # Non-degeneracy: the covariances must carry genuine off-diagonal structure,
  # or the matrix path was never distinguished from two scalar ones.
  expect_gt(max(abs(ref$fc12)), 1e-3)
  expect_gt(max(ref$ACI_metric), 1e-3)
})

test_that("a one-dimensional vector system is the scalar system", {
  # The matrix core, instantiated at one dimension, must reproduce the grade the
  # scalar core carries against the authors' own reference. This is what ties
  # the new code to the package's oldest oracle.
  signal <- read.csv(.aci_mv_fixture("dyad_signal_x.csv"), header = FALSE)
  ref <- read.csv(.aci_mv_fixture("dyad_reference.csv"))
  x <- signal$V2
  p <- list(gamma = 2, d_x = 0.5, F_x = 0.5, F_y = 1, d_y = 0.5)
  n <- length(x)

  comp <- list(
    L_x = array(p$gamma * x, c(1L, 1L, n)),
    f_x = matrix(p$F_x - p$d_x * x, 1L, n),
    L_y = matrix(-p$d_y, 1L, 1L),
    f_y = matrix(p$F_y - p$gamma * x^2, 1L, n),
    S_xoS_x = matrix(0.25, 1L, 1L),
    S_yoS_y = matrix(1, 1L, 1L),
    S_yoS_x = matrix(0, 1L, 1L)
  )
  expect_true(aciR:::.aci_is_mv(comp))

  filt <- aci_filter(matrix(x, 1L, n), comp, 0.001, mu0 = 2, R0 = matrix(0.1))
  smooth <- aci_smoother(matrix(x, 1L, n), comp, 0.001, filt)
  metric <- aci_metric(filt, smooth)
  idx <- seq(1L, n, by = 100L)

  expect_lt(max(abs(filt$mean[1L, idx] - ref$filter_mean)), 1e-6)
  expect_lt(max(abs(smooth$cov[1L, 1L, idx] - ref$smoother_cov)), 1e-6)
  expect_lt(max(abs(metric[idx] - ref$ACI_metric)), 1e-6)
})

test_that("conditioning removes a component's weight without removing it", {
  signal <- read.csv(.aci_mv_fixture("mv_signal.csv"))
  pieces <- .aci_mv_components(signal)
  dt <- 0.002

  conditioned <- aci_conditional(pieces$comp, target = 1L)
  # The weight of the non-target component is exactly zero, and the target's
  # own weight is the inverse of its own noise variance.
  expect_identical(conditioned$S_xoS_x_inv[2L, 2L], 0)
  expect_identical(conditioned$S_xoS_x_inv[1L, 2L], 0)
  expect_equal(
    conditioned$S_xoS_x_inv[1L, 1L],
    1 / pieces$comp$S_xoS_x[1L, 1L],
    tolerance = 1e-12
  )

  # The system still runs, and gives a different answer from the
  # unconditional one -- conditioning changes the estimand.
  filt <- aci_filter(
    pieces$x, conditioned, dt, mu0 = c(0.8, 0.2), R0 = 0.2 * diag(2L)
  )
  smooth <- aci_smoother(pieces$x, conditioned, dt, filt)
  metric <- aci_metric(filt, smooth)
  expect_true(all(is.finite(metric)))
  expect_true(all(metric >= 0))

  base_filt <- aci_filter(
    pieces$x, pieces$comp, dt, mu0 = c(0.8, 0.2), R0 = 0.2 * diag(2L)
  )
  base_metric <- aci_metric(
    base_filt, aci_smoother(pieces$x, pieces$comp, dt, base_filt)
  )
  expect_gt(max(abs(metric - base_metric)), 1e-3)
})

test_that("aci_conditional validates its target", {
  signal <- read.csv(.aci_mv_fixture("mv_signal.csv"))
  comp <- .aci_mv_components(signal)$comp

  expect_error(aci_conditional(comp, target = 0L), "distinct whole numbers")
  expect_error(aci_conditional(comp, target = 5L), "distinct whole numbers")
  expect_error(aci_conditional(comp, target = c(1L, 1L)), "distinct whole")
  expect_error(
    aci_conditional(comp, target = c(1L, 2L)), "nothing to condition"
  )
  expect_error(
    aci_conditional(
      list(L_x = 1, f_x = 1, L_y = 1, f_y = 1, S_xoS_x = 1, S_yoS_y = 1,
           S_yoS_x = 0, S_xoS_y = 0),
      target = 1L
    ),
    "vector-valued components list"
  )
})
