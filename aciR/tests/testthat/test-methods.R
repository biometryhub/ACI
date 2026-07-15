.aci_test_fit <- function(n = 2000L) {
  model <- aci_dyad_model()
  aci(aci_simulate(model, n = n, seed = 1)$x, model)
}

test_that("print methods return their object invisibly", {
  model <- aci_dyad_model()
  fit <- .aci_test_fit(500L)

  expect_output(print(model), "aci_model")
  expect_output(print(fit), "assimilative causal inference")
  expect_invisible(print(model))
  expect_invisible(print(fit))
})

test_that("summary.aci reports the metric and the stability diagnostics", {
  fit <- .aci_test_fit()
  s <- summary(fit)

  expect_s3_class(s, "summary.aci")
  expect_named(
    s,
    c(
      "n", "dt", "span", "label", "metric", "peak", "peak_time",
      "min_filter_cov", "min_smoother_cov", "terminal_residual", "n_clamped"
    )
  )
  expect_identical(s$n, length(fit$x))
  expect_identical(s$dt, fit$dt)
  expect_identical(s$peak, max(fit$aci))
  expect_identical(s$peak_time, fit$t[which.max(fit$aci)])
  expect_identical(s$min_filter_cov, min(fit$filter$cov))

  # The terminal identity is exact, so its residual is a true zero rather than
  # a small number: the diagnostic earns its place by being able to detect
  # otherwise.
  expect_identical(s$terminal_residual, 0)
  expect_gt(s$min_filter_cov, 0)
})

test_that("print.summary.aci returns its object invisibly", {
  s <- summary(.aci_test_fit(500L))
  expect_output(print(s), "diagnostics")
  expect_output(print(s), "terminal identity residual")
  expect_invisible(print(s))
})

test_that("as.data.frame.aci flattens the result to one row per step", {
  fit <- .aci_test_fit(500L)
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_named(
    df,
    c(
      "t", "x", "filter_mean", "filter_cov",
      "smoother_mean", "smoother_cov", "aci"
    )
  )
  expect_identical(nrow(df), length(fit$x))
  expect_identical(df$aci, fit$aci)
  expect_identical(df$t, fit$t)
  expect_identical(df$smoother_cov, fit$smoother$cov)
})

test_that("plot.aci draws without error and returns its object invisibly", {
  fit <- .aci_test_fit(500L)
  path <- tempfile(fileext = ".png")
  on.exit(unlink(path), add = TRUE)

  grDevices::png(path)
  expect_invisible(plot(fit))
  # The device must close before the file exists on disk, so the assertion
  # that anything was drawn belongs after dev.off(), not inside on.exit().
  grDevices::dev.off()

  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)
})

test_that("plot.aci restores the caller's graphical parameters", {
  # The method splits the device into two panels; a user's session must not be
  # left in that state.
  fit <- .aci_test_fit(200L)
  path <- tempfile(fileext = ".png")

  grDevices::png(path)
  before <- graphics::par("mfrow")
  invisible(plot(fit))
  after <- graphics::par("mfrow")
  grDevices::dev.off()
  unlink(path)

  expect_identical(after, before)
  expect_identical(after, c(1L, 1L))
})
