# Independent-oracle grade of the cross-noise path -----------------------------
#
# The dyad oracle grades the core against the authors' reference, but that
# reference (like every scalar model in the reference implementation) sets
# the noise cross-covariance to zero. The terms carrying it therefore execute on
# every dyad run with the term annihilated, and the fixture never grades them,
# while aci_cgns_model(S_yoS_x = ..) exposes those same terms publicly.
#
# This test closes that gap from the transient side: the fixtures come from a
# second implementation of the published CGNS equations, in a different language
# and runtime, on the reference dyad with correlated noise switched on. The
# stationary side is closed analytically in test-identities.R, against the
# algebraic Kalman-Bucy fixed points, which are independent of both
# implementations. Neither test alone is sufficient; together they cover the
# path.

.aci_cross_params <- function() {
  # The Grammians of the harness's noise decomposition:
  #   Sx_1 = 0.6, Sx_2 = 0.3, Sy_1 = 0.5, Sy_2 = 0.8
  list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    S_xoS_x = 0.6^2 + 0.3^2,
    S_yoS_y = 0.5^2 + 0.8^2,
    S_yoS_x = 0.5 * 0.6 + 0.8 * 0.3
  )
}

.aci_cross_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

.aci_cross_components <- function(x) {
  p <- .aci_cross_params()
  list(
    L_x = p$gamma * x,
    f_x = p$F_x - p$d_x * x,
    L_y = -p$d_y,
    f_y = p$F_y - p$gamma * x^2,
    S_xoS_x = p$S_xoS_x,
    S_yoS_y = p$S_yoS_y,
    S_yoS_x = p$S_yoS_x,
    S_xoS_y = p$S_yoS_x
  )
}

test_that("the fixture actually exercises a non-zero cross-covariance", {
  # The point of the fixture is the term it switches on. If this ever reads
  # zero, the test below is grading nothing and must not be allowed to pass
  # quietly, which is exactly how the gap it closes came about.
  p <- .aci_cross_params()
  expect_gt(abs(p$S_yoS_x), 0.5)
  expect_gt(p$S_xoS_x * p$S_yoS_y - p$S_yoS_x^2, 0)
})

test_that("the cross-noise core reproduces the MATLAB oracle to 1e-6", {
  signal_csv <- .aci_cross_fixture("cross_signal_x.csv")
  reference_csv <- .aci_cross_fixture("cross_reference.csv")

  sig <- read.csv(signal_csv, header = FALSE)
  x <- sig$V2

  p <- .aci_cross_params()
  dt <- 0.001
  comp <- .aci_cross_components(x)
  filt <- aci_filter(x, comp, dt, mu0 = p$F_y / p$d_y, R0 = 0.1)
  smooth <- aci_smoother(x, comp, dt, filt)
  metric <- aci_metric(filt, smooth)

  ref <- read.csv(reference_csv)
  idx <- seq(1, length(x), by = 100)
  expect_equal(length(idx), nrow(ref))

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
})

test_that("the high-level entry reproduces the cross-noise oracle to 1e-6", {
  # The same grade through aci_cgns_model(), which is how a user reaches the
  # cross-noise path.
  signal_csv <- .aci_cross_fixture("cross_signal_x.csv")
  reference_csv <- .aci_cross_fixture("cross_reference.csv")

  sig <- read.csv(signal_csv, header = FALSE)
  x <- sig$V2
  p <- .aci_cross_params()

  model <- aci_cgns_model(
    L_x = function(x) p$gamma * x,
    f_x = function(x) p$F_x - p$d_x * x,
    L_y = -p$d_y,
    f_y = function(x) p$F_y - p$gamma * x^2,
    S_xoS_x = p$S_xoS_x,
    S_yoS_y = p$S_yoS_y,
    S_yoS_x = p$S_yoS_x,
    y0 = p$F_y / p$d_y,
    label = "dyad model with correlated noise"
  )
  fit <- aci(x, model, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)

  ref <- read.csv(reference_csv)
  idx <- seq(1, length(x), by = 100)

  expect_equal(fit$filter$mean[idx], ref$filter_mean, tolerance = 1e-6)
  expect_equal(fit$filter$cov[idx], ref$filter_cov, tolerance = 1e-6)
  expect_equal(fit$smoother$mean[idx], ref$smoother_mean, tolerance = 1e-6)
  expect_equal(fit$smoother$cov[idx], ref$smoother_cov, tolerance = 1e-6)
  expect_equal(fit$aci[idx], ref$ACI_metric, tolerance = 1e-6)
})
