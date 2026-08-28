## acir reserve file
## Origin: aci/tests/testthat/test-09-interface.R
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: Whole test file for the formula front-end.
## Verbatim copy from the aci 0.0.30 sources; not modified.

# Identifiable, well-conditioned dyad regime: zero-mean x decorrelates the
# polynomial design; a fluctuating hidden state keeps the posterior trackable.
pp_fit <- list(d_x = 1, gamma = 2, f_x = 1.5, s_x = 0.7,
               d_y = 1, f_y = 0.5, s_y = 2)
make_fit_data <- function(seed = 5, T = 8, keep_latent = TRUE) {
  s <- simulate(model_dyad(params = pp_fit), seed = seed, T = T,
                dt = 5e-3, burn_in = 1)
  d <- data.frame(t = s$obs$t, x1 = s$obs$x[, 1])
  if (keep_latent) d$y <- s$hidden[, 1]
  list(d = d, sim = s)
}

test_that("F1: full-mode aci_fit recovers the dyad and speaks the generics", {
  md <- make_fit_data()
  fit <- suppressWarnings(aci_fit(y ~ x1, data = md$d))
  expect_s3_class(fit, "aci_fit"); expect_equal(fit$mode, "full")
  eq <- fit$report$equations
  expect_true(all(c("x1", "x1*y1") %in% eq$dx1$kept))
  expect_lt(abs(eq$dx1$coef[["x1*y1"]] - 2) / 2, 0.35)
  # plumbing identity: the pipeline reproduces direct LS on the same design
  Z <- as.matrix(md$d[, c("x1", "y")]); n <- nrow(Z)
  dtt <- stats::median(diff(md$d$t))
  Th <- eval_library(cgns_library(1, 1, 2), Z[-n, , drop = FALSE], md$d$t[-n])
  kd <- eq$dx1$kept
  bd <- solve(crossprod(Th[, kd]), crossprod(Th[, kd], diff(Z[, 1]) / dtt))
  expect_lt(max(abs(drop(bd) - eq$dx1$coef)), 1e-6)
  # standard generics
  expect_output(print(fit), "peak ACI")
  expect_type(coef(fit), "list"); expect_equal(nobs(fit), nrow(md$d))
  p <- predict(fit); expect_equal(names(p), c("t", "mean", "sd"))
  expect_equal(length(fitted(fit)), nrow(md$d))
  r <- residuals(fit); expect_equal(dim(r), c(nrow(md$d) - 1L, 1L))
  ss <- simulate(fit, seed = 1); expect_s3_class(ss, "aci_sim")
  grDevices::pdf(NULL)
  plot(fit, "latent"); plot(fit, "metric")
  plot(suppressWarnings(aci(fit$model, fit$obs, init = fit$init)))
  grDevices::dev.off()
  # the reconstruction tracks the withheld truth
  expect_gt(cor(fitted(fit), md$sim$hidden[, 1]), 0.65)
})

test_that("F2: cir() extractor and summary() with surrogate verdict", {
  md <- make_fit_data(seed = 8, T = 5)
  fit <- suppressWarnings(aci_fit(y ~ x1, data = md$d))
  # min_M = NULL: the extractor is under test here, not low-signal masking.
  fc <- cir(fit, min_M = NULL)
  expect_equal(fc$direction, rep("forward", nrow(fc)))
  expect_true(any(is.finite(fc$tau)))
  bc <- cir(fit, direction = "backward", at = c(2, 3.5))
  expect_equal(nrow(bc), 2L)
  expect_true(all(bc$tau >= -1e-8, na.rm = TRUE))
  sm <- suppressWarnings(summary(fit, verdict = TRUE, B = 9))
  expect_output(print(sm), "Learned equations")
  expect_true(sm$verdict$p_value >= 0 && sm$verdict$p_value <= 1)
  expect_lt(sm$verdict$p_value, 0.2)          # coupling detected
  expect_true(sm$verdict$null_type %in% c("decoupled", "linearized"))
})

test_that("F3: energy pairs in aci_fit; partial mode runs (experimental)", {
  md <- make_fit_data(seed = 5)
  fit <- suppressWarnings(aci_fit(y ~ x1, data = md$d,
                                  energy_pairs = list(c(obs = 1, hid = 1))))
  eq <- fit$report$equations
  expect_lt(abs(eq$dx1$coef[["x1*y1"]] + eq$dy$coef[["x1^2"]]), 1e-8)
  mp <- make_fit_data(seed = 5, T = 4, keep_latent = FALSE)
  fitp <- suppressWarnings(aci_fit(y ~ x1, data = mp$d, n_samples = 8))
  expect_match(fitp$mode, "partial")
  expect_match(fitp$report$mode, "experimental")
  expect_equal(nrow(predict(fitp)), nrow(mp$d))
  expect_output(print(fitp), "peak ACI")
})

test_that("F4: marshalling -- groups, Date time, and clear failures", {
  md <- make_fit_data(seed = 3, T = 3)
  d2 <- md$d
  d2$g <- rep(c("a", "b"), each = ceiling(nrow(d2) / 2))[seq_len(nrow(d2))]
  fit <- suppressWarnings(aci_fit(y ~ x1, data = d2, group = "g"))
  expect_equal(fit$n_segments, 2L)
  dd <- md$d
  dd$date <- as.Date("2026-01-01") + seq_len(nrow(dd)) - 1
  dd$t <- NULL
  fitd <- suppressWarnings(aci_fit(y ~ x1, data = dd))
  expect_s3_class(fitd, "aci_fit")            # Date time auto-detected
  expect_error(suppressWarnings(aci_fit(y ~ x1 + missing_col, data = md$d)),
               class = "aci_error_obs_contract")
  expect_error(aci_fit(~x1, data = md$d), class = "aci_error_model_contract")
})

test_that("cir cache stores tables and cannot mix estimation methods", {
  md <- make_fit_data(seed = 4, T = 2)
  fit <- suppressWarnings(aci_fit(y ~ x1, data = md$d))
  l1 <- suppressWarnings(cir(fit, full = TRUE, method = "l1_linf", tol = 0))
  exact <- suppressWarnings(cir(fit, full = TRUE, method = "exact", tol = 0))
  # Masking left at its default so this matches `exact` term for term
  direct <- suppressWarnings(
    forward_cir(lag_table(fit$model, fit$obs, mode = "forward",
                          tol = 0, init = fit$init), method = "exact"))
  expect_equal(exact$tau, direct$tau)
  expect_false(identical(l1$tau, exact$tau))
})

test_that("aci_fit rejects irregular, mixed-latent, and mixed-dt segments", {
  md <- make_fit_data(seed = 10, T = 1)
  irregular <- md$d
  irregular$t[3] <- irregular$t[3] + 0.001
  expect_error(aci_fit(y ~ x1, irregular), class = "aci_error_obs_contract")
  a <- md$d; b <- md$d; b$y <- NULL
  expect_error(aci_fit(y ~ x1, list(a, b)), class = "aci_error_obs_contract")
  b <- md$d[seq(1, nrow(md$d), by = 2), ]
  expect_error(aci_fit(y ~ x1, list(md$d, b)), class = "aci_error_obs_contract")
})


test_that("plot.aci_fit latent view includes truth and the 2 sd band in its range", {
  sim <- simulate(model_dyad(), seed = 8, T = 2, dt = 0.01, burn_in = 0.5)
  dat <- data.frame(t = sim$obs$t, x1 = sim$obs$x[, 1], y = sim$hidden[, 1])
  fit <- suppressWarnings(aci_fit(y ~ x1, data = dat))
  p <- predict(fit)
  grDevices::pdf(NULL)
  plot(fit, "latent", truth = dat$y)
  usr <- graphics::par("usr")
  grDevices::dev.off()
  # par("usr") extends the requested ylim by 4%, so containment is strict
  expect_lte(usr[3], min(dat$y, p$mean - 2 * p$sd))
  expect_gte(usr[4], max(dat$y, p$mean + 2 * p$sd))
})
