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

## That record's filter variance also falls below one step of hidden noise at
## the fourth interval; the tests of the stop muffle that warning by class.
.ocg_quiet_range <- function(expr) withCallingHandlers(
  expr,
  aci_warn_leading_order = function(w) invokeRestart("muffleWarning")
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
      e <- .ocg_err(.ocg_quiet_range(.ocg_quiet_stiff(aci_online(
        run[[1L]], r$obs, lag = cs$lag, filter = run[[2L]], init = r$init
      ))))
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

  ## the dyad record: lag 1 runs where lag 2 stops, and reports the interval
  ## at which the filter variance is below one step of hidden noise
  d <- .ocg_dyad()
  fd <- .ocg_quiet_stiff(aci_filter(d$model, d$obs, init = d$init))
  ws <- list()
  on1 <- withCallingHandlers(
    aci_online(d$model, d$obs, lag = 1, filter = fd),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_length(ws, 1L)
  expect_s3_class(ws[[1L]], "aci_warn_leading_order")
  expect_identical(ws[[1L]]$index, 4L)
  expect_equal(ws[[1L]]$time, 0.3)
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


# A Theorem 3 update outside its range ----------------------------------------

## One step of the dyad's hidden noise is 0.01 at this step; `ratio` is that
## over the prior variance.
.ocg_tight <- function(ratio) {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 1, t_end = 0.5, dt = 0.01)
  list(model = m, obs = s$obs,
       init = list(mean = s$hidden[1, 1], cov = matrix(0.01 / ratio)))
}

## A linear record that `aci()` accepts with a step ratio of 1.02.
.ocg_accepted <- function() {
  m <- aci_linear_model(lambda_x = 5, lambda_y = -0.5, fx = 0, fy = 0,
                        sigma_x = 0.1, sigma_y = 1)
  s <- simulate(m, seed = 3, t_end = 1, dt = 0.01, ic = list(x0 = 0, y0 = 0.3))
  list(model = m, obs = s$obs, init = list(mean = 0.3, cov = matrix(0.01 / 1.02)))
}

## Two hidden components with correlated hidden and observed noise.
.ocg_correlated <- function(prior) {
  m <- .compiled_oracle_mv_model()
  s <- simulate(m, seed = 1, t_end = 0.2, dt = 0.01)
  list(model = m, obs = as_obs(s),
       init = list(mean = c(0.8, 0.2), cov = diag(prior, 2)))
}

.ocg_n_range <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, "aci_warn_leading_order")) n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  n
}


test_that("a step ratio above 1 warns once on every Theorem 3 route", {
  r <- .ocg_tight(10)
  f <- aci_filter(r$model, r$obs, init = r$init)
  rt <- unserialize(serialize(list(model = r$model, filter = f), NULL))
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = Inf, init = r$init)), 1L)
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = Inf, filter = f)), 1L)
  expect_identical(.ocg_n_range(
    aci_online(rt$model, r$obs, lag = Inf, filter = rt$filter)), 1L)
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = 5, init = r$init)), 1L)
  expect_identical(.ocg_n_range(
    lag_table(r$model, r$obs, init = r$init)), 1L)
  expect_identical(.ocg_n_range(
    lag_table(r$model, r$obs, init = r$init, mode = "full")), 1L)

  ## a closure model takes the same route
  cl <- .ocg_closure()
  cl$init <- list(mean = 0, cov = matrix(0.05^2 * 0.09 / 4))
  expect_identical(.ocg_n_range(
    aci_online(cl$model, cl$obs, lag = Inf, init = cl$init)), 1L)

  ## a record aci() accepts: the table it builds, and each range route
  a <- .ocg_accepted()
  expect_identical(.ocg_n_range(
    aci(a$model, a$obs, init = a$init, keep = "table")), 1L)
  for (keep in c("none", "paths")) {
    res <- aci(a$model, a$obs, init = a$init, keep = keep)
    expect_identical(.ocg_n_range(aci_range(res)), 1L)
  }
})


test_that("the warning names the interval and carries the ratio", {
  r <- .ocg_tight(10)
  w <- tryCatch(aci_online(r$model, r$obs, lag = Inf, init = r$init),
                warning = function(w) w)
  expect_s3_class(w, "aci_warn_leading_order")
  expect_s3_class(w, "aci_warning")
  expect_equal(w$ratio, 10, tolerance = 1e-12)
  expect_identical(w$index, 1L)
  expect_equal(w$time, 0)
  expect_match(conditionMessage(w), "index 1 (time 0)", fixed = TRUE)

  ## the call still returns, and the first anchor shows what the warning is
  ## about: a smoothed variance above the filter variance
  f <- aci_filter(r$model, r$obs, init = r$init)
  on <- suppressWarnings(aci_online(r$model, r$obs, lag = Inf, filter = f))
  expect_s3_class(on, "da_path_gaussian")
  expect_gt(on$cov[1L, 1L, 1L], 100 * f$cov[1L, 1L, 1L])
})


test_that("a step ratio below 1 raises nothing", {
  for (r in list(.ocg_tight(0.99), .ocg_tight(0.1))) {
    expect_identical(.ocg_n_range(
      aci_online(r$model, r$obs, lag = Inf, init = r$init)), 0L)
    expect_identical(.ocg_n_range(
      aci_online(r$model, r$obs, lag = 5, init = r$init)), 0L)
    expect_identical(.ocg_n_range(
      lag_table(r$model, r$obs, init = r$init)), 0L)
  }
  ## the boundary is between 0.99 and 1.01
  r <- .ocg_tight(1.01)
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = Inf, init = r$init)), 1L)
  ## the automatic prior
  r <- .ocg_tight(1)
  expect_identical(.ocg_n_range(aci_online(r$model, r$obs, lag = Inf)), 0L)
  expect_identical(.ocg_n_range(lag_table(r$model, r$obs)), 0L)
})


test_that("with correlated noise the ratio uses the conditional hidden noise", {
  ## `gyy Rf^-1 dt` alone is 1.18 here, and nothing is wrong
  below <- .ocg_correlated(0.006)
  expect_identical(.ocg_n_range(
    aci_online(below$model, below$obs, lag = Inf, init = below$init)), 0L)
  expect_identical(.ocg_n_range(
    aci_online(below$model, below$obs, lag = 5, init = below$init)), 0L)
  expect_identical(.ocg_n_range(
    lag_table(below$model, below$obs, init = below$init)), 0L)

  above <- .ocg_correlated(0.004)
  expect_identical(.ocg_n_range(
    aci_online(above$model, above$obs, lag = 5, init = above$init)), 1L)
  expect_identical(.ocg_n_range(
    lag_table(above$model, above$obs, init = above$init)), 1L)
  w <- tryCatch(
    aci_online(above$model, above$obs, lag = Inf, init = above$init),
    warning = function(w) w
  )
  expect_s3_class(w, "aci_warn_leading_order")
  expect_identical(w$index, 1L)
  b <- .compile_cgns_run(above$model, above$obs, NULL)
  co <- .compiled_co(b, 1L)
  Q <- co$gyy - co$gyx %*% .compiled_ginv(b, 1L) %*% t(co$gyx)
  want <- max(eigen(Q, symmetric = TRUE, only.values = TRUE)$values) *
    b$dt / 0.004
  expect_gt(want, 1)
  expect_equal(w$ratio, want, tolerance = 1e-12)
})


test_that("the scalar and per-interval kernels report the same ratio once", {
  r <- .ocg_tight(10)
  b <- .compile_cgns_run(r$model, as_obs(r$obs), NULL)
  f <- aci_filter(r$model, r$obs, init = r$init)
  rec <- .aci_reg_new("none", b$t)
  w1 <- tryCatch(.online_aux_scalar(b, f, rec = rec), warning = function(w) w)
  expect_s3_class(w1, "aci_warn_leading_order")
  expect_true(rec$step_warned)
  ## the latch: the same recorder is silent afterwards, on either kernel
  expect_no_warning(.online_aux_scalar(b, f, rec = rec))
  expect_no_warning(.thmD1_aux_compiled(b, 1L, f$cov[, , 1L], rec = rec))
  w2 <- tryCatch(
    .thmD1_aux_compiled(b, 1L, f$cov[, , 1L], rec = .aci_reg_new("none", b$t)),
    warning = function(w) w
  )
  expect_s3_class(w2, "aci_warn_leading_order")
  expect_equal(w2$ratio, w1$ratio, tolerance = 1e-12)
  expect_identical(w2$index, w1$index)
  ## an interval below 1 is silent on a fresh recorder
  expect_no_warning(
    .thmD1_aux_compiled(b, 2L, f$cov[, , 2L], rec = .aci_reg_new("none", b$t))
  )
})


test_that("a floored filter covariance is left out of the step ratio", {
  r <- .ocg_tight(0.1)
  b <- .compile_cgns_run(r$model, as_obs(r$obs), NULL)
  bad <- aci_filter(r$model, r$obs, init = r$init)
  bad$cov[1L, 1L, 10L] <- -1e-9
  ws <- character(0)
  sm <- withCallingHandlers(
    .smoother_thmD1_compiled(b, bad, validate = FALSE, warn_cost = FALSE,
                             regularize = "floor"),
    warning = function(w) {
      ws <<- c(ws, class(w)[1L])
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(ws, "aci_warn_regularized")
  expect_identical(sm$meta$regularization$n_events, 1L)
})


## The dyad with a wide prior, whose explicit Riccati step overshoots at the
## second point; the filter floors that variance and records one event there.
.ocg_floored <- function() {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 1, t_end = 1, dt = 0.01)
  init <- list(mean = s$hidden[1, 1], cov = matrix(20))
  f <- suppressWarnings(aci_filter(m, s$obs, init = init, regularize = "floor"))
  list(model = m, obs = s$obs, init = init, filter = f)
}

## The first warning of a class, and the call's error if it stops after it.
.ocg_first <- function(expr, class = "aci_warn_leading_order") {
  got <- NULL
  tryCatch(withCallingHandlers(expr, warning = function(w) {
    if (is.null(got) && inherits(w, class)) got <<- w
    invokeRestart("muffleWarning")
  }), error = function(e) NULL)
  got
}


test_that("a filter floor ends the step ratio on every route", {
  r <- .ocg_floored()
  sites <- r$filter$meta$regularization$sites
  expect_identical(sites$site, "filter_explicit")
  expect_identical(as.integer(sites$first_index), 2L)
  rt <- unserialize(serialize(list(model = r$model, filter = r$filter), NULL))
  ## a filter supplied after a floor
  expect_identical(.ocg_n_range(
    lag_table(r$model, r$obs, filter = r$filter)), 0L)
  expect_identical(.ocg_n_range(
    lag_table(rt$model, r$obs, filter = rt$filter)), 0L)
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = Inf, filter = r$filter)), 0L)
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = 5, filter = r$filter)), 0L)
  ## a filter floored inside the call
  expect_identical(.ocg_n_range(
    lag_table(r$model, r$obs, init = r$init, regularize = "floor")), 0L)
  expect_identical(.ocg_n_range(
    aci(r$model, r$obs, init = r$init, regularize = "floor", keep = "table")),
    0L)
  expect_identical(.ocg_n_range(
    aci_online(r$model, r$obs, lag = Inf, init = r$init,
               regularize = "floor")), 0L)
})


test_that("both kernels leave out the intervals from a filter floor on", {
  r <- .ocg_floored()
  b <- .compile_cgns_run(r$model, as_obs(r$obs), NULL)
  rec <- .aci_reg_new("none", b$t)
  expect_identical(.aci_ratio_end(rec, r$filter), 2)
  expect_identical(.aci_ratio_end(rec, list(meta = list())), Inf)
  ## the per-interval kernel: the floored interval is far above 1
  w <- tryCatch(
    .thmD1_aux_compiled(b, 2L, r$filter$cov[, , 2L],
                        rec = .aci_reg_new("none", b$t)),
    warning = function(w) w
  )
  expect_s3_class(w, "aci_warn_leading_order")
  expect_gt(w$ratio, 1e8)
  expect_no_warning(
    .thmD1_aux_compiled(b, 2L, r$filter$cov[, , 2L],
                        rec = .aci_reg_new("none", b$t), ratio_end = 2)
  )
  ## the scalar kernel reads the same index off the filter
  expect_no_warning(
    .online_aux_scalar(b, r$filter, rec = .aci_reg_new("none", b$t))
  )
  ## a filter floor formed in the same call is read off the recorder
  rec <- .aci_reg_new("floor", b$t)
  rec$j <- 3L
  .aci_reg_note(rec, "filter_explicit", 3L, b$t[3L], -1)
  expect_identical(.aci_ratio_end(rec, list(meta = list())), 3)
  ## a floor at a smoother site does not end the ratio
  rec <- .aci_reg_new("floor", b$t)
  .aci_reg_note(rec, "smoother_onelag", 3L, b$t[3L], -1)
  expect_identical(.aci_ratio_end(rec, list(meta = list())), Inf)
})


test_that("a ratio above 1 before the first floor still warns", {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 7, t_end = 8, dt = 0.01, burn_in = 1)
  ob <- observed_trajectory(s$obs$t, s$obs$x, names = "x")
  idx <- seq.int(1L, length(ob$t), by = 20L)
  ob <- observed_trajectory(ob$t[idx], ob$x[idx, , drop = FALSE], names = "x")
  init <- list(mean = 0, cov = matrix(1, 1, 1))
  f <- suppressWarnings(aci_filter(m, ob, init = init, regularize = "floor"))
  expect_identical(as.integer(f$meta$regularization$sites$first_index), 37L)
  expect_identical(.ocg_n_range(
    lag_table(m, ob, filter = f, regularize = "floor")), 1L)
  w <- .ocg_first(lag_table(m, ob, filter = f, regularize = "floor"))
  expect_identical(w$index, 35L)
  expect_gt(w$ratio, 1)
})


test_that("the message names the prior only at the first interval", {
  w1 <- .ocg_first(with(.ocg_tight(10),
                        aci_online(model, obs, lag = Inf, init = init)))
  expect_identical(w1$index, 1L)
  expect_match(conditionMessage(w1), "prior variance", fixed = TRUE)

  ## the coarse dyad: its filter variance falls below one step of hidden
  ## noise at the fourth interval, and the call stops after the warning
  w4 <- .ocg_first(with(.ocg_dyad(),
                        aci_online(model, obs, lag = Inf, init = init)))
  expect_identical(w4$index, 4L)
  expect_gt(w4$ratio, 1)
  expect_match(conditionMessage(w4), "Reduce dt.", fixed = TRUE)
  expect_no_match(conditionMessage(w4), "prior", fixed = TRUE)
})
