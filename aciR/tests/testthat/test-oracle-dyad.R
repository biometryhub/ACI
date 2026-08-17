# Binding independent-oracle grade of the packaged core vs the MATLAB oracle ---
#
# This is the package's independent-oracle test. It runs aciR's own conditional
# Gaussian core (aci_filter / aci_smoother / aci_metric) on the MATLAB input
# signal and asserts the outputs reproduce the MATLAB reference at the 301
# sampled indices. The two fixtures were produced by an independent MATLAB
# implementation and are shipped in the package's inst/extdata, so this test is
# self-contained: it runs in the development tree, under R CMD check and on an
# installed copy alike, and never skips. The dyad parameters, initial filtered
# mean mu0 = F_y / d_y = 2, initial covariance R0 = 0.1 and step dt = 0.001 are
# those of the reference run.

.aci_oracle_dyad_params <- function() {
  list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
}

.aci_oracle_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

test_that("the packaged core reproduces the MATLAB oracle to 1e-6", {
  signal_csv <- .aci_oracle_fixture("dyad_signal_x.csv")
  reference_csv <- .aci_oracle_fixture("dyad_reference.csv")

  # The signal fixture is headerless: column one is time, column two the
  # observed signal x that the reference run was driven with.
  sig <- read.csv(signal_csv, header = FALSE)
  x <- sig$V2

  p <- .aci_oracle_dyad_params()
  dt <- 0.001
  comp <- aci_dyad_components(x, p)
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

  # A single hard absolute-error gate that records the worst deviation across
  # all five series in the reporter, independent of the relative tolerance
  # semantics used by expect_equal above.
  max_abs_error <- max(
    abs(filt$mean[idx] - ref$filter_mean),
    abs(filt$cov[idx] - ref$filter_cov),
    abs(smooth$mean[idx] - ref$smoother_mean),
    abs(smooth$cov[idx] - ref$smoother_cov),
    abs(metric[idx] - ref$ACI_metric)
  )
  expect_lt(max_abs_error, 1e-6)
})

test_that("the high-level aci() entry reproduces the MATLAB oracle to 1e-6", {
  signal_csv <- .aci_oracle_fixture("dyad_signal_x.csv")
  reference_csv <- .aci_oracle_fixture("dyad_reference.csv")

  sig <- read.csv(signal_csv, header = FALSE)
  x <- sig$V2

  p <- .aci_oracle_dyad_params()
  model <- aci_dyad_model()
  fit <- aci(x, model, dt = 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)

  ref <- read.csv(reference_csv)
  idx <- seq(1, length(x), by = 100)

  expect_equal(fit$filter$mean[idx], ref$filter_mean, tolerance = 1e-6)
  expect_equal(fit$filter$cov[idx], ref$filter_cov, tolerance = 1e-6)
  expect_equal(fit$smoother$mean[idx], ref$smoother_mean, tolerance = 1e-6)
  expect_equal(fit$smoother$cov[idx], ref$smoother_cov, tolerance = 1e-6)
  expect_equal(fit$aci[idx], ref$ACI_metric, tolerance = 1e-6)
})
