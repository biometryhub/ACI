test_that("E1: EnKBF/EnKBS moments approach the CGNS oracle", {
  m <- model_dyad()
  s <- simulate(m, seed = 55, T = 1.5, dt = 5e-3)
  init <- list(mean = 0, cov = diag(1, 1))
  filt <- suppressWarnings(da_filter(m, s$obs, init = init))
  smoo <- da_smooth(m, s$obs, filter = filt)
  fr <- enkbf(m, s$obs, m = 800, seed = 8, ic_sampler = init)
  sm <- enkbs(m, fr$path, fr$noise)
  fg <- as_gaussian(fr$path); sg <- as_gaussian(sm)
  burn <- seq_len(60)
  sc <- stats::sd(smoo$mean[-burn, 1])
  expect_lt(sqrt(mean((fg$mean[-burn, 1] - filt$mean[-burn, 1])^2)) / sc, 0.15)
  expect_lt(sqrt(mean((sg$mean[-burn, 1] - smoo$mean[-burn, 1])^2)) / sc, 0.15)
})
test_that("E5/E6: noise reuse enforced; no inflation pathway backward", {
  m <- model_dyad()
  s <- simulate(m, seed = 56, T = 0.5, dt = 5e-3)
  fr <- enkbf(m, s$obs, m = 30, seed = 1)
  bad <- fr$noise; bad$B[1, 2, 1] <- bad$B[1, 2, 1] + 1
  expect_error(enkbs(m, fr$path, bad), class = "aci_error_noise_mismatch")
  expect_false("inflation" %in% names(formals(enkbs)))
})

test_that("P2/P3 forward ensemble CIR uses a complete lagged EnKBS family", {
  m <- model_dyad()
  s <- simulate(m, seed = 91, T = 0.06, dt = 0.01)
  a <- suppressWarnings(aci(m, s$obs, engine = "ensemble", m = 30,
                            seed = 7, keep = "table"))
  tab <- a$table
  expect_s3_class(tab, "lag_table")
  expect_identical(tab$meta$engine, "ensemble")
  expect_identical(tab$meta$reference_smoother, "enkbs_full_horizon")
  expect_equal(tab$diag, a$aci, tolerance = 0)
  expect_equal(length(tab$rows), length(tab$t))
  expect_true(all(vapply(tab$rows, function(z)
    all(is.finite(z)) && all(z >= 0) && tail(z, 1) == 0, logical(1))))
  expect_identical(tab$L, length(tab$t) - seq_along(tab$t))

  f <- suppressWarnings(forward_cir(a, min_M = 0))
  expect_s3_class(f, "cir_result")
  expect_identical(f$direction, "forward")
  expect_error(backward_cir(a), class = "aci_error_not_implemented")

  no_table <- suppressWarnings(aci(m, s$obs, engine = "ensemble", m = 30,
                                   seed = 7, keep = "paths"))
  expect_error(forward_cir(no_table), class = "aci_error_not_implemented")

  direct <- suppressWarnings(ensemble_lag_table(m, s$obs, m = 30, seed = 7))
  expect_equal(direct$diag, tab$diag, tolerance = 0)
  expect_equal(direct$rows, tab$rows, tolerance = 0)
  reused <- aci(m, s$obs, table = direct)
  expect_equal(reused$aci, direct$diag, tolerance = 0)
  expect_identical(reused$meta$smoother_scheme, "enkbs_full_horizon")
})

test_that("ensemble lagged posteriors approach the CGNS finite-lag oracle", {
  m <- model_dyad()
  s <- simulate(m, seed = 12, T = 0.05, dt = 0.005)
  init <- list(mean = 0, cov = matrix(1, 1, 1))
  exact <- lag_table(m, s$obs, mode = "full", init = init)
  ens <- suppressWarnings(ensemble_lag_table(
    m, s$obs, m = 800, seed = 3, ic_sampler = init))
  x <- unlist(exact$rows, use.names = FALSE)
  y <- unlist(ens$rows, use.names = FALSE)
  expect_lt(sqrt(mean((x - y)^2)) / max(x), 0.10)
})

test_that("ensemble engine rejects unsupported noise contracts and bad m", {
  mc <- cgns_model(
    Lx = function(t, x) matrix(0, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy1 = function(t, x) matrix(0.5, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
  s <- simulate(mc, seed = 4, T = 0.1, dt = 0.01,
                ic = list(x0 = 0, y0 = 0))
  expect_error(enkbf(mc, s$obs, m = 20),
               class = "aci_error_ensemble_noise_contract")
  expect_error(enkbf(model_dyad(), s$obs, m = 1),
               class = "aci_error_ensemble_rank")
  expect_error(enkbf(model_dyad(), s$obs, m = 20, inflation = 0.9),
               class = "aci_error_dims")
  expect_error(enkbf(model_dyad(), s$obs, m = 20, seed = -1),
               class = "aci_error_dims")
  expect_error(enkbf(model_dyad(), s$obs, m = 20, seed = 1.5),
               class = "aci_error_dims")
  expect_error(enkbf(model_dyad(), s$obs, m = 20, unknown = TRUE),
               class = "aci_error_dims")
  expect_error(aci(model_dyad(), s$obs, engine = "ensemble", m = NA_real_),
               class = "aci_error_ensemble_rank")
  expect_error(aci(model_dyad(), s$obs, engine = "ensemble", m = 4.5),
               class = "aci_error_ensemble_rank")
  expect_error(ensemble_lag_table(model_dyad(), s$obs, m = 1),
               class = "aci_error_ensemble_rank")
})
test_that("E7: Gaspari-Cohn taper properties", {
  expect_equal(gaspari_cohn(0), 1)
  expect_equal(gaspari_cohn(2), 0)
  expect_true(all(diff(gaspari_cohn(seq(0, 2.5, by = 0.01))) <= 1e-12))
  set.seed(2); l <- 8
  P <- crossprod(matrix(rnorm(l * l), l)) + diag(l)
  C <- gaspari_cohn(abs(outer(1:l, 1:l, "-")) / 3)
  ev <- eigen(C * P, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), -1e-10)                 # Schur product preserves PSD
})

test_that("localization preserves matrices and supports m below l", {
  m <- model_l96(n = 12)
  s <- simulate(m, seed = 14, T = 0.04, dt = 0.01)
  loc <- localization_spec(m$meta$coords$hidden, radius = 3,
                           coords_obs = m$meta$coords$obs,
                           distance = "cyclic", period = 12)
  expect_identical(dim(loc$C1), c(m$l, m$k))
  expect_identical(dim(loc$C2), c(m$l, m$l))
  fr <- enkbf(m, s$obs, m = 5, seed = 2, localization = loc)
  sm <- enkbs(m, fr$path, fr$noise)
  expect_s3_class(sm, "da_path_ensemble")
  expect_error(enkbs(m, fr$path, fr$noise,
                     localization = localization_spec(
                       m$meta$coords$hidden, radius = 2,
                       coords_obs = m$meta$coords$obs,
                       distance = "cyclic", period = 12)),
               class = "aci_error_dims")
  fr0 <- enkbf(m, s$obs, m = 5, seed = 2)
  expect_error(enkbs(m, fr0$path, fr0$noise),
               class = "aci_error_ensemble_rank")
  expect_error(aci(m, s$obs, engine = "ensemble", m = 5,
                   seed = 2, localization = loc),
               class = "aci_error_ensemble_rank")

  # An all-ones taper is rank one and must not obtain an undocumented ridge
  # from the generic safe-Cholesky helper.
  bad_loc <- list(C1 = matrix(1, m$l, m$k), C2 = matrix(1, m$l, m$l))
  bad_fr <- enkbf(m, s$obs, m = 5, seed = 2, localization = bad_loc)
  expect_error(enkbs(m, bad_fr$path, bad_fr$noise),
               class = "aci_error_ensemble_rank")
})

test_that("ensemble smoother verifies its source observations and model", {
  m <- model_dyad(); s <- simulate(m, seed = 7, T = 0.1, dt = 0.01)
  f <- da_filter.stochastic_model(m, s$obs, m = 10, seed = 3)
  bad_obs <- observed_trajectory(s$obs$t, s$obs$x + 0.01)
  expect_error(da_smooth.stochastic_model(m, bad_obs, filter = f),
               class = "aci_error_dims")
  other <- model_dyad(params = list(d_x = 0.6, gamma = 2, f_x = 0.5,
                                    s_x = 0.5, d_y = 0.5, f_y = 1, s_y = 1))
  expect_error(da_smooth.stochastic_model(other, s$obs, filter = f),
               class = "aci_error_model_contract")
  fr <- enkbf(m, s$obs, m = 10, seed = 3)
  expect_error(enkbs(other, fr$path, fr$noise),
               class = "aci_error_model_contract")
})

test_that("ensemble conversion and smoothing never hide rank or rerun filters", {
  bad <- new_ens_path(0:1, array(1:8, c(2, 2, 2)), "filter")
  expect_error(as_gaussian(bad), class = "aci_error_ensemble_rank")
  wrong_grid <- new_ens_path(0, array(rnorm(16), c(1, 2, 8)), "filter")
  expect_error(as_gaussian(wrong_grid), class = "aci_error_dims")

  m <- model_dyad(); s <- simulate(m, seed = 8, T = 0.1, dt = 0.01)
  fr <- enkbf(m, s$obs, m = 8, seed = 12)
  wrong_kind <- fr$path; wrong_kind$kind <- "smoother"
  expect_error(enkbs(m, wrong_kind, fr$noise),
               class = "aci_error_noise_mismatch")
  expected <- da_smooth.stochastic_model(m, s$obs, filter = fr)
  got <- sample_paths(m, s$obs, n_samples = 99, method = "enkbs",
                      filter = fr, seed = 999)
  expect_equal(got, expected$members, tolerance = 0)
  expect_error(sample_paths(m, s$obs, n_samples = 8, method = "enkbs",
                            filter = fr, init = list(mean = 0, cov = 1)),
               class = "aci_error_dims")

  sampler <- function(members) matrix(seq(-1, 1, length.out = members), 1, members)
  fr_fixed <- enkbf(m, s$obs, m = 8, ic_sampler = sampler, seed = 4)
  sm_fixed <- da_smooth.stochastic_model(
    m, s$obs, m = 8, ic_sampler = sampler, noise = fr_fixed$noise, seed = 999)
  expect_equal(sm_fixed$members,
               enkbs(fr_fixed$model, fr_fixed$path, fr_fixed$noise)$members,
               tolerance = 0)
  sampled <- sample_paths(m, s$obs, n_samples = 8, method = "enkbs",
                          init = list(mean = 0, cov = matrix(0.2, 1, 1)), seed = 6)
  expect_identical(dim(sampled), c(1L, length(s$obs$t), 8L))
})
test_that("nil_causality_check passes on a decoupled system, fails when coupled", {
  m0 <- cgns_model(Lx = function(t, x) matrix(0, 1, 1),
                   fx = function(t, x) -x, Ly = function(t, x) matrix(-1, 1, 1),
                   fy = function(t, x) 0, Sx1 = function(t, x) matrix(0.4, 1, 1),
                   Sy2 = function(t, x) matrix(0.6, 1, 1), k = 1, l = 1)
  nc <- suppressWarnings(nil_causality_check(m0, direction = list(cause = 1, effect = 1),
                                             T = 3, dt = 2e-3))
  expect_true(nc$structural_pass); expect_true(nc$empirical_pass)
  m1 <- model_dyad()
  nc1 <- suppressWarnings(nil_causality_check(m1, direction = list(cause = 1, effect = 1),
                                              T = 3, dt = 2e-3))
  expect_false(nc1$structural_pass); expect_false(nc1$empirical_pass)
})

test_that("surrogate-null verdict separates coupled from decoupled systems", {
  # A y-identifiable regime: fluctuating hidden state, lower observation
  # noise. (The default-parameter dyad can sit in a weakly-identifiable
  # window where near-zero LR is the *correct* answer -- see T1c.)
  pp <- list(d_x = 0.5, gamma = 2, f_x = 0.5, s_x = 0.3,
             d_y = 1, f_y = 0, s_y = 1.5)
  m <- model_dyad(params = pp)
  s <- simulate(m, seed = 21, T = 3, dt = 5e-3, burn_in = 0.5)
  ini <- list(mean = 0, cov = diag(1, 1))
  st <- nil_surrogate_test(m, s$obs, B = 19, seed = 2, init = ini)
  expect_lt(st$p_value, 0.1)                       # real coupling detected
  expect_gt(st$observed_lr, max(st$null_lr))
  m0 <- cgns_model(Lx = function(t, x) matrix(0, 1, 1),
                   fx = m$fx, Ly = m$Ly, fy = m$fy,
                   Sx1 = m$Sx1, Sy2 = m$Sy2, k = 1, l = 1)
  s0 <- simulate(m0, seed = 22, T = 3, dt = 5e-3, burn_in = 0.5)
  st0 <- nil_surrogate_test(m, s0$obs, B = 19, seed = 3, init = ini)
  expect_gt(st0$p_value, 0.2)                      # nil data not flagged
})
