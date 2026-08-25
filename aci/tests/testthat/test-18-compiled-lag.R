.expect_compiled_lag_table_equal <- function(actual, expected,
                                             tolerance = 1e-12) {
  expect_s3_class(actual, "lag_table")
  expect_identical(names(actual), names(expected))
  expect_equal(actual$t, expected$t, tolerance = 0)
  expect_equal(actual$dt, expected$dt, tolerance = 0)
  expect_identical(actual$mode, expected$mode)
  expect_equal(actual$diag, expected$diag, tolerance = tolerance)
  expect_equal(actual$rows, expected$rows, tolerance = tolerance)
  expect_identical(actual$L, expected$L)
  expect_equal(actual$diag_signal, expected$diag_signal,
               tolerance = tolerance)
  expect_equal(actual$diag_dispersion, expected$diag_dispersion,
               tolerance = tolerance)
  expect_equal(actual$tailbnd, expected$tailbnd, tolerance = tolerance)
  expect_equal(actual$onelag, expected$onelag, tolerance = tolerance)
  expect_identical(names(actual$meta), names(expected$meta))
  expect_identical(actual$meta$nontarget, expected$meta$nontarget)
  expect_equal(actual$meta$tol, expected$meta$tol, tolerance = 0)
  expect_equal(actual$meta$window, expected$meta$window, tolerance = 0)
  expect_equal(actual$meta$max_lag, expected$meta$max_lag, tolerance = 0)
  expect_equal(actual$meta$init, expected$meta$init, tolerance = tolerance)
  expect_identical(actual$meta$source_model, expected$meta$source_model)
  expect_equal(actual$meta$source_obs_x, expected$meta$source_obs_x,
               tolerance = 0)
  expect_identical(actual$meta$reference_smoother,
                   expected$meta$reference_smoother)
  expect_identical(actual$meta$stop_index, expected$meta$stop_index)
}


.expected_thmD1_aux <- function(co, Rf, dt) {
  Rf <- as.matrix(Rf)
  l <- nrow(Rf)
  k <- nrow(co$gxx)
  Rfi <- chol_solve(Rf, diag(l), "Rf")
  Gi <- chol_solve(co$gxx, diag(k), "gxx")
  Gx <- co$Lx + t(co$gyx) %*% Rfi
  Gy <- co$Ly + co$gyy %*% Rfi
  K <- Gi %*% Gx
  H <- Rfi %*% (co$Ly %*% Rf + Rf %*% t(co$Ly) + co$gyy)
  E <- diag(l) + (co$gyx %*% Gi %*% Gx - Gy) * dt
  KR <- K %*% Rf
  F <- -Rf %*% (
    t(K) +
      (t(Gx) %*% KR %*% t(K) - Rfi %*% t(H) %*% Rf %*% t(K) +
         t(co$Ly) %*% t(K)) * dt -
      t(co$Lx) %*% (Gi + KR %*% t(K) * dt)
  )
  list(E = E, F = F)
}


.compiled_lag_matrix_setup <- function() {
  model <- cgns_model(
    Lx = function(t, x) matrix(c(
      0.65 + 0.04 * x[1], -0.12,
      0.08, 0.45 + 0.03 * x[2]
    ), 2, 2, byrow = TRUE),
    fx = function(t, x) c(-0.25 * x[1] + 0.08,
                           -0.18 * x[2] - 0.03),
    Ly = function(t, x) matrix(c(-0.7, -0.04, 0.09, -0.55),
                                2, 2, byrow = TRUE),
    fy = function(t, x) c(0.1 * sin(t) - 0.06 * x[1],
                           0.12 * cos(t) + 0.04 * x[2]),
    Sx1 = function(t, x) matrix(c(0.65, 0.05, 0, 0.55), 2, 2),
    Sy1 = function(t, x) matrix(c(0.05, 0.01, 0, 0.04), 2, 2),
    Sy2 = function(t, x) matrix(c(0.72, 0.03, 0, 0.64), 2, 2),
    k = 2, l = 2, name = "compiled-lag-matrix"
  )
  t <- seq(0, 0.12, by = 0.01)
  x <- cbind(0.25 + 0.08 * sin(4 * t),
             -0.15 + 0.06 * cos(3 * t))
  obs <- observed_trajectory(t, x)
  init <- list(
    mean = c(0.15, -0.08),
    cov = matrix(c(0.42, 0.025, 0.025, 0.31), 2, 2)
  )
  list(
    model = model, obs = obs, init = init,
    bundle = .compile_cgns_complete(model, obs)
  )
}


test_that("compiled Theorem-3 auxiliaries match the published matrix equation", {
  model <- model_dyad()
  sim <- simulate(model, seed = 31, T = 0.12, dt = 0.01, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 1.8, cov = matrix(0.2, 1, 1))
  bundle <- .compile_cgns_complete(model, obs)
  compiled_filter <- .cgns_filter_compiled(bundle, init)

  for (j in seq_len(bundle$N)) {
    expected_aux <- .expected_thmD1_aux(
      eval_coefs(model, obs$t[j], obs$x[j, ]),
      compiled_filter$cov[, , j], obs$dt
    )
    compiled_aux <- .thmD1_aux_compiled(
      bundle, j, compiled_filter$cov[, , j]
    )
    expect_equal(compiled_aux, expected_aux, tolerance = 1e-14)
  }

  compiled <- .smoother_thmD1_compiled(bundle, compiled_filter)
  expected_mean <- matrix(NA_real_, bundle$N1, bundle$l)
  expected_cov <- array(NA_real_, c(bundle$l, bundle$l, bundle$N1))
  expected_mean[bundle$N1, ] <- compiled_filter$mean[bundle$N1, ]
  expected_cov[, , bundle$N1] <- compiled_filter$cov[, , bundle$N1]
  for (j in bundle$N:1L) {
    co <- eval_coefs(model, obs$t[j], obs$x[j, ])
    aux <- .expected_thmD1_aux(co, compiled_filter$cov[, , j], obs$dt)
    one <- .onelag_stats(
      co, aux, compiled_filter$mean[j, ], compiled_filter$cov[, , j],
      expected_mean[j + 1L, ], expected_cov[, , j + 1L],
      obs$x[j + 1L, ] - obs$x[j, ], obs$dt, bundle$l
    )
    expected_mean[j, ] <- one$mu
    expected_cov[, , j] <- one$R
  }
  expect_equal(compiled$mean, expected_mean, tolerance = 1e-13)
  expect_equal(compiled$cov, expected_cov, tolerance = 1e-13)
  expect_identical(compiled$meta$route, "thmD1")
  expect_identical(compiled$meta$nsub, 1L)
})


test_that("compiled scalar lag tables match all established storage modes", {
  model <- model_dyad()
  sim <- simulate(model, seed = 91, T = 0.14, dt = 0.01, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  bundle <- .compile_cgns_complete(model, obs)

  settings <- list(
    one_lag = list(tol = 1e-7, window = 2L, max_lag = Inf),
    forward = list(tol = 1e-7, window = 2L, max_lag = 6L),
    full = list(tol = 1e-7, window = 2L, max_lag = Inf)
  )
  for (mode in names(settings)) {
    z <- settings[[mode]]
    public <- lag_table(
      model, obs, mode = mode, tol = z$tol, window = z$window,
      max_lag = z$max_lag, init = init
    )
    compiled <- .lag_table_compiled(
      bundle, mode = mode, tol = z$tol, window = z$window,
      max_lag = z$max_lag, init = init
    )
    .expect_compiled_lag_table_equal(compiled, public, tolerance = 1e-12)
  }
})


test_that("compiled correlated matrix lag execution matches public assembly", {
  ds <- .compiled_lag_matrix_setup()
  expect_true(ds$bundle$correlated_noise)
  compiled_filter <- .cgns_filter_compiled(ds$bundle, ds$init)
  compiled_smoother <- .smoother_thmD1_compiled(
    ds$bundle, compiled_filter
  )
  expect_identical(compiled_smoother$meta$route, "thmD1")

  for (mode in c("one_lag", "forward", "full")) {
    public <- lag_table(
      ds$model, ds$obs, mode = mode, tol = 1e-8, window = 2L,
      init = ds$init
    )
    compiled <- .lag_table_compiled(
      ds$bundle, mode = mode, tol = 1e-8, window = 2L,
      init = ds$init
    )
    .expect_compiled_lag_table_equal(compiled, public, tolerance = 1e-11)
  }
})


test_that("compiled lag core retains internal smoother-only state", {
  ds <- .compiled_lag_matrix_setup()
  compiled_filter <- .cgns_filter_compiled(ds$bundle, ds$init)
  compiled_smoother <- .smoother_thmD1_compiled(
    ds$bundle, compiled_filter
  )
  compiled <- .lagtable_core_compiled(
    ds$bundle, compiled_filter, compiled_smoother,
    mode = "smoother_only", tol = 0, window = Inf, max_lag = Inf
  )
  expect_null(compiled$rows)
  expect_equal(compiled$diag,
               gaussian_kl_path(compiled_smoother, compiled_filter)$total,
               tolerance = 1e-11)
  expect_equal(compiled$smu, compiled_smoother$mean, tolerance = 1e-11)
  expect_equal(compiled$scov, compiled_smoother$cov, tolerance = 1e-11)
  expect_true(all(compiled$tailbnd == 0))
})


test_that("compiled Theorem-3 stepper and lag-table warning policies remain", {
  model <- model_dyad()
  sim <- simulate(model, seed = 7, T = 0.08, dt = 0.01, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  bundle <- .compile_cgns_complete(model, obs)
  implicit <- .cgns_filter_compiled(
    bundle, init, stepper = "implicit", nsub = 1L
  )
  expect_error(
    .smoother_thmD1_compiled(bundle, implicit),
    class = "aci_error_stepper"
  )
  expect_warning(
    rebuilt <- .lag_table_compiled(
      bundle, mode = "one_lag", filter = implicit
    ),
    class = "aci_warn_stepper"
  )
  expect_s3_class(rebuilt, "lag_table")
  expect_identical(rebuilt$meta$init, implicit$meta$init)

  explicit <- .cgns_filter_compiled(bundle, init)
  backward <- .cgns_smoother_compiled(bundle, explicit)
  expect_warning(
    recomputed <- .lag_table_compiled(
      bundle, mode = "one_lag", filter = explicit, smoother = backward
    ),
    class = "aci_warn_stepper"
  )
  expect_s3_class(recomputed, "lag_table")
  expect_identical(recomputed$meta$reference_smoother,
                   "thmD1_online_complete")
})
