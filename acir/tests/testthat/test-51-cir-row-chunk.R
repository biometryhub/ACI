# Scalar rows that the covariance policy splits into many passes --------------

## Primitives whose row at anchor `j` has the variances `target` at its first
## cells; later cells keep the sign of the last one.
.crc_record <- function(target, N, j = 20L, R_f = 1, seed = 20261125L) {
  set.seed(seed)
  E <- stats::runif(N, 0.99, 0.9995)
  dmu <- stats::rnorm(N, sd = 0.05)
  dR <- stats::runif(N, -1e-4, 1e-4)
  m <- length(target)
  d <- c(1, cumprod(E[seq.int(j, length.out = m - 1L)]))
  dR[seq.int(j, length.out = m)] <- diff(c(R_f, target)) / (d * d)
  list(
    E = E, dmu = dmu, dR = dR, j = j, R_f = R_f,
    sp = .cir_scalar_primitives(
      lapply(E, function(e) matrix(e, 1L, 1L)), matrix(dmu, ncol = 1L),
      lapply(dR, function(r) matrix(r, 1L, 1L))
    )
  )
}

## The row cell by cell. A variance that is not positive is floored for its
## own divergence, by the rule of `spd_floor()`, and carried on unfloored.
.crc_recursion <- function(rec, mu_f = 0.3, max_lag = Inf) {
  E <- rec$E
  N <- length(E)
  j <- rec$j
  d <- c(1, cumprod(E[j:(N - 1L)]))
  mu_s <- mu_f + sum(d * rec$dmu[j:N])
  R_s <- rec$R_f + sum(d * d * rec$dR[j:N])
  chp <- sqrt(R_s)
  kl <- function(mu, R) {
    if (!(R > 0)) R <- max(1e-12 * max(abs(R), 1e-300), 1e-300)
    chq <- sqrt(R)
    w <- (mu - mu_s) / chq
    A <- chp / chq
    max(0.5 * (w * w + A * A - 1) + log(chq) - log(chp), 0)
  }
  cells <- as.integer(min(N - j + 1L, max_lag))
  row <- numeric(cells + 1L)
  raw <- numeric(cells)
  row[1L] <- kl(mu_f, rec$R_f)
  mu <- mu_f
  R <- rec$R_f
  D <- 1
  for (i in seq_len(cells)) {
    n <- j + i - 1L
    if (i > 1L) D <- D * E[n - 1L]
    mu <- mu + D * rec$dmu[n]
    R <- R + D * rec$dR[n] * D
    raw[i] <- R
    row[i + 1L] <- kl(mu, R)
  }
  list(row = row, raw = raw, mu_s = mu_s, R_s = R_s, mu_f = mu_f)
}

.crc_row <- function(rec, ref, reg, max_lag = Inf) {
  N <- length(rec$E)
  .cir_scalar_row(
    rec$j, rec$sp, numeric(N + 1L), numeric(N + 1L), N, ref$mu_f, rec$R_f,
    ref$mu_s, ref$R_s, ref$row[1L], ref$R_s, 0, Inf, max_lag, reg
  )
}

## The value of a call with every warning it raised, muffled.
.crc_run <- function(expr) {
  w <- list()
  value <- withCallingHandlers(expr, warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

.crc_count <- function(w, class) {
  sum(vapply(w, inherits, logical(1L), what = class))
}

.crc_gap <- function(got, ref) max(abs(got - ref) / pmax(abs(ref), 1))

.crc_gate <- 1e-12

## The dyad sampled every 0.1: the third anchor's row has a variance below
## zero at 12 of its 13 cells.
.crc_dyad <- function() {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 1, t_end = 1.5, dt = 0.005)
  idx <- seq.int(1L, length(s$obs$t), by = 20L)
  list(
    model = m,
    obs = observed_trajectory(s$obs$t[idx], s$obs$x[idx, , drop = FALSE]),
    init = list(mean = s$hidden[1, 1], cov = matrix(0.1))
  )
}

.crc_closure <- function(n, dt) {
  list(
    model = aci_model(
      Lx = function(t, x) matrix(2), fx = function(t, x) 0,
      Ly = function(t, x) matrix(1), fy = function(t, x) 0,
      Sx1 = function(t, x) matrix(0.5), Sy2 = function(t, x) matrix(0.05),
      k = 1L, l = 1L, name = "closure-linear"
    ),
    obs = observed_trajectory((0:(n - 1L)) * dt, rep(0, n)),
    init = list(mean = 0, cov = matrix(0.5))
  )
}

.crc_metric_events <- function(x) {
  s <- x$meta$regularization$sites
  s$n[s$site == "metric_reference"]
}


test_that("a run of 40 floored cells gives the recursion's row", {
  rec <- .crc_record(c(rep(1, 5L), rep(-2, 40L), 1), N = 120L)
  ref <- .crc_recursion(rec)
  expect_identical(which(!(ref$raw > 0)), 6:45)
  reg <- .aci_reg_for("floor", NA_real_)
  r <- .crc_run(.crc_row(rec, ref, reg))
  expect_length(r$warnings, 0L)
  got <- r$value
  expect_length(got$row, length(ref$row))
  expect_true(all(is.finite(got$row)))
  expect_lt(.crc_gap(got$row, ref$row), .crc_gate)
  expect_identical(.aci_reg_freeze(reg)$n_events, 40L)
  expect_identical(got$lag_kept, 101L)
  expect_false(got$frozen)
})


test_that("the row is the recursion's when max_lag cuts inside the run", {
  rec <- .crc_record(c(rep(1, 5L), rep(-2, 40L), 1), N = 120L)
  ref <- .crc_recursion(rec, max_lag = 30)
  reg <- .aci_reg_for("floor", NA_real_)
  r <- .crc_run(.crc_row(rec, ref, reg, max_lag = 30))
  expect_length(r$warnings, 0L)
  got <- r$value
  expect_length(got$row, 31L)
  expect_lt(.crc_gap(got$row, ref$row), .crc_gate)
  expect_identical(.aci_reg_freeze(reg)$n_events, 25L)
  expect_identical(got$lag_kept, 30L)
  expect_identical(got$tail, got$row[31L])
})


test_that("five floored cells between kept ones take eleven passes", {
  rec <- .crc_record(rep(c(1, -2), length.out = 11L), N = 120L)
  ref <- .crc_recursion(rec)
  expect_identical(which(!(ref$raw > 0)), c(2L, 4L, 6L, 8L, 10L))
  reg <- .aci_reg_for("floor", NA_real_)
  r <- .crc_run(.crc_row(rec, ref, reg))
  expect_length(r$warnings, 0L)
  expect_lt(.crc_gap(r$value$row, ref$row), .crc_gate)
  expect_identical(.aci_reg_freeze(reg)$n_events, 5L)
})


test_that("a row of exactly ten passes raises no warning", {
  ## one kept cell, eight floored ones, one kept cell that ends the record
  rec <- .crc_record(c(1, rep(-2, 8L), 1), N = 29L)
  ref <- .crc_recursion(rec)
  expect_identical(which(!(ref$raw > 0)), 2:9)
  expect_length(ref$raw, 10L)
  reg <- .aci_reg_for("floor", NA_real_)
  r <- .crc_run(.crc_row(rec, ref, reg))
  expect_length(r$warnings, 0L)
  expect_lt(.crc_gap(r$value$row, ref$row), .crc_gate)
  expect_identical(.aci_reg_freeze(reg)$n_events, 8L)
})


test_that("the strict policy still stops at the first such cell", {
  rec <- .crc_record(c(rep(1, 5L), rep(-2, 40L), 1), N = 120L)
  ref <- .crc_recursion(rec)
  expect_error(
    .crc_row(rec, ref, .aci_reg_for("none", NA_real_)),
    class = "aci_error_spd"
  )
  d <- .crc_dyad()
  expect_error(
    .crc_run(lag_table(d$model, d$obs, init = d$init, regularize = "none")),
    class = "aci_error_spd"
  )
})


test_that("lag_table() returns the floored dyad table", {
  d <- .crc_dyad()
  r <- .crc_run(
    lag_table(d$model, d$obs, init = d$init, regularize = "floor")
  )
  lt <- r$value
  expect_s3_class(lt, "lag_table")
  expect_identical(.crc_count(r$warnings, "aci_warning"), length(r$warnings))
  expect_identical(.crc_count(r$warnings, "aci_warn_regularized"), 1L)
  expect_true(all(is.finite(unlist(lt$rows))))
  expect_identical(lengths(lt$rows)[3L], 14L)
  expect_true(lt$meta$regularization$fired)
  expect_identical(.crc_metric_events(lt), 12L)

  cut <- .crc_run(lag_table(
    d$model, d$obs, init = d$init, regularize = "floor", max_lag = 12
  ))
  expect_identical(.crc_count(cut$warnings, "aci_warning"),
                   length(cut$warnings))
  expect_identical(lengths(cut$value$rows)[3L], 13L)
  expect_identical(cut$value$rows[[3L]], lt$rows[[3L]][1:13])
})


test_that("aci_range() on an aci() result returns for every keep", {
  d <- .crc_dyad()
  tab <- .crc_run(aci_range(.crc_run(
    lag_table(d$model, d$obs, init = d$init, regularize = "floor")
  )$value))$value
  for (keep in c("none", "paths", "table")) {
    made <- .crc_run(
      aci(d$model, d$obs, init = d$init, regularize = "floor", keep = keep)
    )
    r <- .crc_run(aci_range(made$value))
    w <- c(made$warnings, r$warnings)
    expect_s3_class(r$value, "cir_result")
    expect_identical(.crc_count(w, "aci_warning"), length(w), info = keep)
    expect_identical(.crc_count(w, "aci_warn_regularized"), 1L, info = keep)
    expect_identical(.crc_metric_events(r$value), 12L, info = keep)
    expect_equal(r$value$tau, tab$tau, tolerance = .crc_gate, info = keep)
    expect_identical(r$value$status, tab$status, info = keep)
  }
})


test_that("a model built from closures returns on both routes", {
  d <- .crc_closure(20L, 0.09)
  r <- .crc_run(
    lag_table(d$model, d$obs, init = d$init, regularize = "floor")
  )
  expect_s3_class(r$value, "lag_table")
  expect_identical(.crc_count(r$warnings, "aci_warning"), length(r$warnings))
  expect_identical(.crc_count(r$warnings, "aci_warn_regularized"), 1L)
  expect_true(all(is.finite(unlist(r$value$rows))))
  expect_identical(.crc_metric_events(r$value), 26L)

  made <- .crc_run(
    aci(d$model, d$obs, init = d$init, regularize = "floor", keep = "none")
  )
  s <- .crc_run(aci_range(made$value))
  w <- c(made$warnings, s$warnings)
  expect_identical(.crc_count(w, "aci_warning"), length(w))
  expect_identical(.crc_metric_events(s$value), 26L)
  expect_equal(s$value$tau, .crc_run(aci_range(r$value))$value$tau,
               tolerance = .crc_gate)
})


test_that("a table whose row takes exactly ten passes raises no base warning", {
  d <- .crc_closure(13L, 0.12)
  r <- .crc_run(
    lag_table(d$model, d$obs, init = d$init, regularize = "floor")
  )
  expect_s3_class(r$value, "lag_table")
  expect_identical(.crc_count(r$warnings, "aci_warning"), length(r$warnings))
  expect_identical(.crc_count(r$warnings, "aci_warn_regularized"), 1L)
  expect_identical(.crc_metric_events(r$value), 10L)
})
