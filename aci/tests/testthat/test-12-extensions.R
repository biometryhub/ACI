
test_that("nil check scores the requested cause and effect rather than all channels", {
  m <- cgns_model(
    Lx = function(t, x) diag(c(0, 2), 2),
    fx = function(t, x) -x,
    Ly = function(t, x) diag(c(-1, -1), 2),
    fy = function(t, x) c(0, 0),
    Sx1 = function(t, x) diag(c(0.4, 0.4), 2),
    Sy2 = function(t, x) diag(c(0.7, 0.7), 2),
    k = 2, l = 2,
    meta = list(ic_default = list(x0 = c(0, 0), y0 = c(0, 0))))
  sim <- simulate(m, seed = 91, T = 1, dt = 0.005)
  full <- suppressWarnings(aci(m, sim$obs))
  expect_gt(max(full$aci), 1e-4)

  selected_null <- suppressWarnings(nil_causality_check(
    m, sim$obs, direction = list(cause = 1, effect = 1), floor = 1e-4))
  expect_true(selected_null$structural_pass)
  expect_true(selected_null$empirical_pass)
  expect_equal(selected_null$direction, list(cause = 1L, effect = 1L))

  selected_link <- suppressWarnings(nil_causality_check(
    m, sim$obs, direction = list(cause = 2, effect = 2), floor = 1e-4))
  expect_false(selected_link$structural_pass)
  expect_false(selected_link$empirical_pass)
  expect_error(nil_causality_check(
    m, sim$obs, direction = list(cause = 3, effect = 1)),
    class = "aci_error_direction")
})

test_that("ensemble/oracle cross-validation shares one initial Gaussian prior", {
  m <- conditionally_linear_model(
    lambda_x = 0.8, lambda_y = -1, fx = function(t, x) -x,
    fy = 0, sigma_x = 0.4, sigma_y = 0.6)
  s <- simulate(m, seed = 4, T = 0.2, dt = 0.01,
                ic = list(x0 = 0, y0 = 0))
  ini <- list(mean = 2.5, cov = matrix(0.3, 1, 1))
  cv <- suppressWarnings(cross_validate(
    m, s$obs, m_grid = 6, n_rep = 1, seed = 7, init = ini))
  expect_equal(attr(cv, "initial_prior")$mean, ini$mean)
  expect_equal(attr(cv, "initial_prior")$cov, ini$cov)
  expect_match(attr(cv, "source_status"), "same initial Gaussian prior")
  expect_error(suppressWarnings(cross_validate(
    m, s$obs, m_grid = 6, n_rep = 1, init = ini,
    ic_sampler = ini)), class = "aci_error_model_contract")
  expect_error(cross_validate(m, s$obs, m_grid = 1, n_rep = 1),
               class = "aci_error_ensemble_rank")
  expect_error(cross_validate(m, s$obs, m_grid = 6, n_rep = 0),
               class = "aci_error_dims")
})

test_that("surrogate calibration preserves an absolute model clock", {
  m <- cgns_model(
    Lx = function(t, x) matrix(0.3, 1, 1),
    fx = function(t, x) -x + 0.1 * sin(t),
    Ly = function(t, x) matrix(-1, 1, 1),
    fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.5, 1, 1),
    Sy2 = function(t, x) matrix(0.6, 1, 1), k = 1, l = 1,
    meta = list(ic_default = list(x0 = 0, y0 = 0)))
  tt <- 10 + 0.01 * 0:30
  ob <- observed_trajectory(tt, matrix(0.1 * sin(tt), ncol = 1))
  st <- suppressWarnings(nil_surrogate_test(
    m, ob, B = 1, burn_frac = 0.2, seed = 2,
    init = list(mean = 0, cov = matrix(1, 1, 1))))
  expect_true(st$calibration$time_origin_preserved)
  expect_equal(st$calibration$time_origin, 10)
  expect_equal(st$requested_B, 1L)
  expect_match(st$calibration$source_status, "Experimental")
  expect_silent(capture.output(print(st)))
})

test_that("calibration returns fitted multichannel noise and RMSE", {
  build <- function(par) cgns_model(
    Lx = function(t, x) matrix(0, 2, 1),
    fx = function(t, x) c(-par$d * x[1], -2 * par$d * x[2]),
    Ly = function(t, x) matrix(-1, 1, 1),
    fy = function(t, x) 0,
    Sx1 = function(t, x) diag(c(0.1, 0.2), 2),
    Sy2 = function(t, x) matrix(0.3, 1, 1),
    k = 2, l = 1, name = "two-channel calibration fixture")
  tt <- seq(0, 1, by = 0.05)
  det <- simulate_deterministic(build(list(d = 0.4)), tt, c(1, 2), 0)
  Y <- det$x + cbind(0.02 * sin(7 * tt), 0.05 * cos(5 * tt))
  colnames(Y) <- c("canopy", "moisture")
  ob <- observed_trajectory(tt, Y, names = colnames(Y))
  cal <- calibrate_cgns(
    build, ob, start = list(d = 0.5), bounds = list(d = c(0.1, 1)),
    y0 = 0, restarts = 1, sigma_floor = 0.005, maxit = 100)
  fitted_sd <- sqrt(diag(cgns_grams(cal$model, tt[1], Y[1, ])$gxx))
  expect_equal(unname(fitted_sd), unname(cal$sigma), tolerance = 1e-10)
  expect_equal(names(cal$rmse_by_channel), colnames(Y))
  expect_equal(cal$rmse, sqrt(mean((cal$det$x - Y)^2)))
  expect_match(cal$meta$noise_calibration, "hidden-process noise were not fitted")
  expect_silent(calibrate_cgns(
    build, ob, start = list(d = 0.5), bounds = list(d = c(0.1, 1)),
    y0 = 0, weights = c(moisture = 2, canopy = 1),
    restarts = 1, sigma_floor = 0.005, maxit = 30))
})

test_that("cir_table forwards method, truncation, and strength controls with provenance", {
  m <- conditionally_linear_model(
    lambda_x = 0.7, lambda_y = -1,
    fx = function(t, x) -0.4 * x, fy = 0,
    sigma_x = 0.45, sigma_y = 0.6)
  sim <- simulate(m, seed = 112, T = 0.12, dt = 0.01,
                  ic = list(x0 = 0, y0 = 0))
  ini <- list(mean = 0, cov = matrix(0.5, 1, 1))

  ct <- suppressWarnings(cir_table(
    m, sim$obs, init = ini, at = 0.10, direction = "both",
    method = "l1_linf", tol = 0, min_M = 0,
    max_lag = 1, window = 2))
  expect_s3_class(ct, "aci_cir_table")
  expect_equal(ct$direction, c("forward", "backward"))
  expect_match(ct$bound[ct$direction == "forward"],
               "lower_ratio_on_truncated_table_only")
  expect_match(ct$bound[ct$direction == "backward"],
               "upper_ratio_on_O\\(dt\\)_T_minus_dt_grid")
  tab <- attr(ct, "lag_table", exact = TRUE)
  expect_s3_class(tab, "lag_table")
  expect_equal(tab$meta$max_lag, 1)
  expect_equal(tab$meta$window, 2L)
  p <- attr(ct, "aci_provenance", exact = TRUE)
  expect_identical(p$source_model, m)
  expect_equal(p$source_obs_x, sim$obs$x)
  expect_equal(p$method, "l1_linf")
  expect_equal(p$min_M, 0)
  expect_equal(p$max_lag, 1)
  expect_match(p$source_status, "Experimental")

  masked <- suppressWarnings(cir_table(
    m, sim$obs, init = ini, at = 0.10, direction = "both",
    method = "exact", min_M = 1e9, max_lag = 1, window = 1))
  expect_true(all(is.na(masked$tau)))
  expect_error(cir_table(m, sim$obs, nonsense = 1), class = "aci_error_dims")
  expect_error(cir_table(m, sim$obs, window = 1.5), class = "aci_error_dims")
  expect_error(cir_table(m, sim$obs, max_lag = 0), class = "aci_error_dims")
  expect_error(cir_table(m, sim$obs, min_M = -1), class = "aci_error_dims")
})

test_that("aci_check handles CIR reuse and checks full covariance PSD ordering", {
  m <- cgns_model(
    Lx = function(t, x) diag(c(0.7, 0.35), 2),
    fx = function(t, x) -0.4 * x,
    Ly = function(t, x) diag(c(-0.8, -1.2), 2),
    fy = function(t, x) c(0, 0),
    Sx1 = function(t, x) diag(c(0.45, 0.55), 2),
    Sy2 = function(t, x) diag(c(0.5, 0.75), 2),
    k = 2, l = 2, name = "two-latent diagnostics",
    meta = list(ic_default = list(x0 = c(0, 0), y0 = c(0.2, -0.3))))
  sim <- simulate(m, seed = 113, T = 0.10, dt = 0.01)
  ini <- list(mean = c(0, 0), cov = diag(c(0.6, 0.9)))

  no_cir <- NULL
  expect_output(
    no_cir <- aci_check(m, sim$obs, init = ini, test = FALSE,
                        plot = FALSE, cir = FALSE),
    "full covariance PSD order")
  expect_null(no_cir$cir)
  manual_min <- vapply(seq_along(sim$obs$t), function(j)
    min(eigen((no_cir$filter$cov[, , j] - no_cir$smoother$cov[, , j] +
                 t(no_cir$filter$cov[, , j] - no_cir$smoother$cov[, , j])) / 2,
              symmetric = TRUE, only.values = TRUE)$values), numeric(1))
  expect_equal(no_cir$consistency_min_eigen, manual_min, tolerance = 1e-14)
  expect_equal(no_cir$consistency_by_step,
               manual_min >= -no_cir$consistency_tolerance)
  expect_match(no_cir$consistency_definition, "positive semidefinite")

  ct <- suppressWarnings(cir_table(
    m, sim$obs, init = ini, direction = "forward",
    min_M = 0, max_lag = 2, window = 1))
  reused <- NULL
  expect_output(
    reused <- aci_check(m, sim$obs, init = ini, test = FALSE,
                        plot = FALSE, cir = ct),
    "forward CIR")
  expect_identical(reused$cir, ct)

  computed <- NULL
  expect_output(
    computed <- aci_check(m, sim$obs, init = ini, test = FALSE,
                          plot = FALSE, cir = TRUE),
    "forward CIR")
  expect_s3_class(computed$cir, "aci_cir_table")

  bare <- ct
  attr(bare, "aci_provenance") <- NULL
  expect_error(aci_check(m, sim$obs, init = ini, test = FALSE,
                         plot = FALSE, cir = bare),
               class = "aci_error_dims")
  changed_obs <- observed_trajectory(sim$obs$t, sim$obs$x + 1e-4)
  expect_error(aci_check(m, changed_obs, init = ini, test = FALSE,
                         plot = FALSE, cir = ct),
               class = "aci_error_dims")
  expect_error(aci_check(m, sim$obs, init = ini, test = FALSE,
                         plot = FALSE, cir = FALSE, mystery = 1),
               class = "aci_error_dims")
})

test_that("osse_twin aggregates RMSE and coverage across every hidden component", {
  m <- cgns_model(
    Lx = function(t, x) diag(c(0.8, 0.25), 2),
    fx = function(t, x) c(-0.3 * x[1], -0.6 * x[2]),
    Ly = function(t, x) diag(c(-0.7, -1.4), 2),
    fy = function(t, x) c(0, 0),
    Sx1 = function(t, x) diag(c(0.4, 0.65), 2),
    Sy2 = function(t, x) diag(c(0.45, 0.8), 2),
    k = 2, l = 2, name = "two-latent OSSE fixture",
    meta = list(ic_default = list(x0 = c(0, 0), y0 = c(0.3, -0.4))))
  ini <- list(mean = c(0, 0), cov = diag(c(0.7, 1.1)))
  tw <- suppressWarnings(osse_twin(
    m, init = ini, T = 0.12, dt = 0.01, nrep = 1, seed = 31,
    stepper = "implicit", nsub = 2))
  expect_s3_class(tw, "osse_twin")
  expect_equal(tw$l, 2L)
  expect_equal(tw$filter,
               sqrt(mean(tw$per_component$filter^2)), tolerance = 1e-12)
  expect_equal(tw$smoother,
               sqrt(mean(tw$per_component$smoother^2)), tolerance = 1e-12)
  expect_equal(tw$coverage, mean(tw$per_component$coverage), tolerance = 1e-12)
  expect_gt(abs(tw$filter - tw$per_component$filter[1]), 1e-8)
  expect_equal(dim(tw$per_twin_by_component$filter), c(1L, 2L))
  expect_equal(tw$meta$stepper, "implicit")
  expect_equal(tw$meta$nsub, 2L)
  expect_match(tw$meta$aggregation, "all hidden components")
  expect_silent(capture.output(print(tw)))

  expect_error(osse_twin(m, ini, T = 0.1, dt = 0.01, nrep = 0),
               class = "aci_error_dims")
  expect_error(osse_twin(m, ini, T = 0.1, dt = 0.03),
               class = "aci_error_dims")
  expect_error(osse_twin(m, ini, T = 0.1, dt = 0.01, seed = 1.5),
               class = "aci_error_dims")
  expect_error(osse_twin(m, ini, T = 0.1, dt = 0.01, nsub = 0),
               class = "aci_error_dims")
  expect_error(osse_twin(m, ini, T = 0.1, dt = 0.01, y0 = 0),
               class = "aci_error_dims")
  bad_init <- list(mean = c(0, 0), cov = matrix(c(1, 2, 2, 1), 2))
  expect_error(osse_twin(m, bad_init, T = 0.1, dt = 0.01),
               class = "aci_error_spd")
})

test_that("agronomy observation and deterministic on-ramps reject malformed inputs", {
  tt <- 0:4
  vv <- c(0.1, 0.2, 0.15, 0.3, 0.4)
  expect_error(screen_outliers(tt, vv, iter = -1), class = "aci_error_dims")
  expect_error(screen_outliers(tt, vv, iter = 1.5), class = "aci_error_dims")
  zero_iter <- screen_outliers(tt, vv, iter = 0)
  expect_true(all(zero_iter))
  expect_equal(attr(zero_iter, "n_dropped"), 0L)

  expect_error(
    as_uniform_trajectory(tt, list(canopy = vv[-1]), dt = 1),
    class = "aci_error_dims")
  expect_error(
    as_uniform_trajectory(tt, list(canopy = vv), dt = 1,
                          keep = c(TRUE, TRUE)),
    class = "aci_error_dims")
  expect_error(
    as_uniform_trajectory(c(0, 1, 1, 2, 3), list(canopy = vv), dt = 1),
    class = "aci_error_obs_contract")
  expect_error(
    as_uniform_trajectory(tt, list(canopy = vv), dt = 10),
    class = "aci_error_obs_contract")
  with_gap <- vv
  with_gap[2] <- NA_real_
  expect_s3_class(
    as_uniform_trajectory(tt, list(canopy = with_gap), dt = 1,
                          keep = c(TRUE, FALSE, TRUE, TRUE, TRUE)),
    "obs_traj")

  m <- conditionally_linear_model(
    lambda_x = 0.5, lambda_y = -1,
    fx = function(t, x) -x, fy = 0,
    sigma_x = 0.4, sigma_y = 0.6)
  expect_error(simulate_deterministic(m, numeric(0), 0, 0),
               class = "aci_error_dims")
  expect_error(simulate_deterministic(m, c(0, 0.1, 0.1), 0, 0),
               class = "aci_error_dims")
  expect_error(simulate_deterministic(m, c(0, Inf), 0, 0),
               class = "aci_error_dims")
  expect_error(simulate_deterministic(m, 0:1, numeric(0), 0),
               class = "aci_error_dims")
  expect_error(simulate_deterministic(m, 0:1, 0, NA_real_),
               class = "aci_error_dims")
  one <- simulate_deterministic(m, 2, x0 = 0.25, y0 = -0.5)
  expect_equal(one$t, 2)
  expect_equal(one$x, matrix(0.25, 1, 1))
  expect_equal(one$y, matrix(-0.5, 1, 1))
})


test_that("print.nilcheck reports the sensitivities and both flags", {
  nc <- suppressWarnings(nil_causality_check(
    model_dyad(), direction = list(cause = 1, effect = 1), T = 1))
  expect_output(print(nc), "direct sensitivity")
  expect_output(print(nc), "structural_pass: FALSE")
  expect_output(print(nc), "TRUE flags the direction as nil")
})
