# Independent-oracle grade: the packaged CGNS core, run on the MATLAB input
# signal, must reproduce the MATLAB reference outputs to numerical tolerance.
# The fixtures live in the graduate workspace's oracle/ directory, outside the
# package, and were produced by an independent MATLAB implementation. The test
# skips gracefully when the fixtures are not present (for example on an
# installed copy without the development tree).

test_that("packaged core matches the MATLAB oracle to < 1e-6", {
  oracle_dir <- normalizePath(
    file.path(testthat::test_path(), "..", "..", "..", "oracle"),
    mustWork = FALSE
  )
  signal_csv <- file.path(oracle_dir, "dyad_signal_x.csv")
  reference_csv <- file.path(oracle_dir, "dyad_reference.csv")
  testthat::skip_if_not(
    file.exists(signal_csv) && file.exists(reference_csv),
    "oracle fixtures not available"
  )

  sig <- read.csv(signal_csv, header = FALSE)
  x <- sig$V2
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  dt <- 0.001
  comp <- aci_dyad_components(x, p)
  filt <- aci_filter(x, comp, dt, mu0 = p$F_y / p$d_y, R0 = 0.1)
  smooth <- aci_smoother(x, comp, dt, filt)
  aci <- aci_metric(filt, smooth)

  ref <- read.csv(reference_csv)
  idx <- seq(1, length(x), by = 100)
  expect_equal(length(idx), nrow(ref))

  expect_lt(max(abs(filt$mean[idx] - ref$filter_mean)), 1e-6)
  expect_lt(max(abs(filt$cov[idx] - ref$filter_cov)), 1e-6)
  expect_lt(max(abs(smooth$mean[idx] - ref$smoother_mean)), 1e-6)
  expect_lt(max(abs(smooth$cov[idx] - ref$smoother_cov)), 1e-6)
  expect_lt(max(abs(aci[idx] - ref$ACI_metric)), 1e-6)
})

test_that("the high-level aci() entry matches the MATLAB oracle to < 1e-6", {
  oracle_dir <- normalizePath(
    file.path(testthat::test_path(), "..", "..", "..", "oracle"),
    mustWork = FALSE
  )
  signal_csv <- file.path(oracle_dir, "dyad_signal_x.csv")
  reference_csv <- file.path(oracle_dir, "dyad_reference.csv")
  testthat::skip_if_not(
    file.exists(signal_csv) && file.exists(reference_csv),
    "oracle fixtures not available"
  )

  sig <- read.csv(signal_csv, header = FALSE)
  x <- sig$V2
  model <- aci_dyad_model()
  fit <- aci(x, model, dt = 0.001, mu0 = 2, R0 = 0.1)

  ref <- read.csv(reference_csv)
  idx <- seq(1, length(x), by = 100)

  expect_lt(max(abs(fit$filter$mean[idx] - ref$filter_mean)), 1e-6)
  expect_lt(max(abs(fit$filter$cov[idx] - ref$filter_cov)), 1e-6)
  expect_lt(max(abs(fit$smoother$mean[idx] - ref$smoother_mean)), 1e-6)
  expect_lt(max(abs(fit$smoother$cov[idx] - ref$smoother_cov)), 1e-6)
  expect_lt(max(abs(fit$aci[idx] - ref$ACI_metric)), 1e-6)
})
