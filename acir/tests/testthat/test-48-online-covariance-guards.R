# Covariance policy on the fixed-lag online smoother --------------------------

.ocg_linear <- function(l = 1L) {
  ob <- observed_trajectory((0:12) * 0.09, rep(0, 13))
  if (l == 1L) {
    m <- aci_linear_model(lambda_x = 2, lambda_y = 1, fx = 0, fy = 0,
                          sigma_x = 0.5, sigma_y = 0.05)
    init <- list(mean = 0, cov = matrix(0.5))
  } else {
    m <- aci_linear_model(lambda_x = matrix(c(2, 0), 1, 2),
                          lambda_y = diag(c(1, -0.3)), fx = 0, fy = c(0, 0),
                          sigma_x = 0.5, sigma_y = c(0.05, 0.7), k = 1, l = 2)
    init <- list(mean = c(0, 0), cov = diag(c(0.5, 0.5)))
  }
  list(model = m, obs = ob, init = init)
}

.ocg_closure <- function() {
  r <- .ocg_linear(1L)
  r$model <- aci_model(
    Lx = function(t, x) matrix(2), fx = function(t, x) 0,
    Ly = function(t, x) matrix(1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.5), Sy2 = function(t, x) matrix(0.05),
    k = 1L, l = 1L, name = "closure-linear"
  )
  r
}

.ocg_dyad <- function() {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 1, t_end = 1.5, dt = 0.005)
  idx <- seq.int(1L, length(s$obs$t), by = 20L)
  list(
    model = m,
    obs = observed_trajectory(s$obs$t[idx], s$obs$x[idx, , drop = FALSE]),
    init = list(mean = s$hidden[1, 1], cov = matrix(0.1))
  )
}

.ocg_err <- function(expr) tryCatch(expr, error = function(e) e)

## The dyad record sits at the explicit stepper's stability limit, and its
## filter warns about that; only that warning is muffled.
.ocg_quiet_stiff <- function(expr) withCallingHandlers(
  expr,
  aci_warn_riccati_stiff = function(w) invokeRestart("muffleWarning")
)

.ocg_min_eig <- function(p) apply(p$cov, 3L, function(R)
  min(eigen(R, symmetric = TRUE, only.values = TRUE)$values))


test_that("a finite lag stops on a covariance that is not positive definite", {
  cases <- list(
    list(r = .ocg_linear(1L), lag = 8, index = 1L, time = 0),
    list(r = .ocg_linear(2L), lag = 8, index = 1L, time = 0),
    list(r = .ocg_closure(), lag = 8, index = 1L, time = 0),
    list(r = .ocg_dyad(), lag = 2, index = 3L, time = 0.2)
  )
  for (cs in cases) {
    r <- cs$r
    f <- .ocg_quiet_stiff(aci_filter(r$model, r$obs, init = r$init))
    ## a model and filter saved together and read back
    rt <- unserialize(serialize(list(model = r$model, filter = f), NULL))
    runs <- list(list(r$model, NULL), list(r$model, f),
                 list(rt$model, rt$filter))
    for (run in runs) {
      e <- .ocg_err(.ocg_quiet_stiff(aci_online(
        run[[1L]], r$obs, lag = cs$lag, filter = run[[2L]], init = r$init
      )))
      expect_s3_class(e, "aci_error_covariance_not_spd")
      expect_s3_class(e, "aci_error_spd")
      expect_identical(e$site, "smoother_window")
      expect_identical(e$role, "fixed-lag smoothed covariance")
      expect_identical(e$index, cs$index)
      expect_equal(e$time, cs$time)
      expect_lt(e$value, 0)
      expect_match(conditionMessage(e),
                   sprintf("index %d (time %g)", cs$index, cs$time),
                   fixed = TRUE)
    }
  }
})


test_that("the stopping value is the covariance the window route forms", {
  r <- .ocg_linear(1L)
  e <- .ocg_err(aci_online(r$model, r$obs, lag = 8, init = r$init))
  expect_equal(e$value, -0.00787249, tolerance = 1e-5)
  r2 <- .ocg_linear(2L)
  e2 <- .ocg_err(aci_online(r2$model, r2$obs, lag = 8, init = r2$init))
  expect_equal(e2$value, e$value, tolerance = 1e-12)
})


test_that("the anchor and lag pair route applies the same policy", {
  r <- .ocg_linear(1L)
  b <- .compile_cgns_run(r$model, as_obs(r$obs), NULL)
  f <- aci_filter(r$model, r$obs, init = r$init)
  aux <- .online_aux_compiled(b, f)
  e <- .ocg_err(.online_at_compiled(b, f, aux, c(2L, 1L), c(8, 8)))
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_identical(e$site, "smoother_window")
  expect_identical(e$index, 1L)

  rec <- .aci_reg_new("floor", b$t)
  got <- .online_at_compiled(b, f, aux, c(2L, 1L), c(8, 8), rec)
  win <- .online_window_compiled(b, f, aux, 8L, .aci_reg_new("floor", b$t))
  expect_identical(rec$n, 1L)
  expect_equal(got$cov[, , 2L], win$cov[, , 1L], tolerance = 1e-12)
  expect_equal(got$cov[, , 1L], win$cov[, , 2L], tolerance = 1e-12)
  expect_gt(got$cov[1L, 1L, 2L], 0)
})


test_that("regularize = 'floor' floors the covariance and records the event", {
  for (r in list(.ocg_linear(1L), .ocg_linear(2L), .ocg_closure())) {
    ws <- character(0)
    on <- withCallingHandlers(
      aci_online(r$model, r$obs, lag = 8, init = r$init, regularize = "floor"),
      warning = function(w) {
        ws <<- c(ws, class(w)[1L])
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(ws, "aci_warn_regularized")
    reg <- on$meta$regularization
    expect_identical(reg$policy, "floor")
    expect_true(reg$fired)
    expect_identical(reg$n_events, 1L)
    expect_identical(reg$sites$site, "smoother_window")
    expect_identical(reg$sites$role, "fixed-lag smoothed covariance")
    expect_identical(reg$sites$first_index, 1L)
    expect_equal(reg$sites$first_time, 0)
    expect_lt(reg$sites$worst_value, 0)
    expect_true(all(.ocg_min_eig(on) > 0))
    ## the other anchors keep the values the default policy accepts for them
    b <- .compile_cgns_run(r$model, as_obs(r$obs), NULL)
    f <- aci_filter(r$model, r$obs, init = r$init)
    rest <- .online_at_compiled(b, f, .online_aux_compiled(b, f), 2:13, 8)
    expect_equal(on$cov[, , -1L, drop = FALSE], rest$cov, tolerance = 1e-12)
    expect_equal(on$mean[-1L, , drop = FALSE], rest$mean, tolerance = 1e-12)
  }
})


test_that("the option aci.regularize reaches the fixed-lag check", {
  r <- .ocg_linear(1L)
  old <- options(aci.regularize = "floor")
  on.exit(options(old), add = TRUE)
  on <- suppressWarnings(aci_online(r$model, r$obs, lag = 8, init = r$init))
  expect_identical(on$meta$regularization$sites$site, "smoother_window")
})


test_that("lags next to the failing one run and report nothing", {
  ## lag 6 keeps every variance positive on the record where lag 7 does not
  r <- .ocg_linear(1L)
  f <- aci_filter(r$model, r$obs, init = r$init)
  expect_s3_class(
    .ocg_err(aci_online(r$model, r$obs, lag = 7, filter = f)),
    "aci_error_covariance_not_spd"
  )
  on6 <- expect_no_warning(aci_online(r$model, r$obs, lag = 6, filter = f))
  expect_false(on6$meta$regularization$fired)
  expect_identical(on6$meta$regularization$n_events, 0L)
  expect_true(all(on6$cov > 0))
  fl6 <- expect_no_warning(
    aci_online(r$model, r$obs, lag = 6, filter = f, regularize = "floor")
  )
  expect_identical(fl6$cov, on6$cov)
  expect_identical(fl6$mean, on6$mean)

  ## the dyad record: lag 1 runs where lag 2 stops
  d <- .ocg_dyad()
  fd <- .ocg_quiet_stiff(aci_filter(d$model, d$obs, init = d$init))
  on1 <- expect_no_warning(aci_online(d$model, d$obs, lag = 1, filter = fd))
  expect_false(on1$meta$regularization$fired)
  expect_true(all(on1$cov > 0))
})


test_that("the check returns the symmetric part unchanged when it passes", {
  rec <- .aci_reg_new("none", c(0, 1))
  rec$j <- 1L
  set.seed(48)
  for (l in 1:4) {
    A <- matrix(rnorm(l * l), l, l)
    S <- crossprod(A) + diag(0.1, l) + matrix(rnorm(l * l, sd = 1e-3), l, l)
    expect_identical(.cov_guard(S, rec, "smoother_window"), sym(S))
  }
  expect_identical(rec$n, 0L)
})
