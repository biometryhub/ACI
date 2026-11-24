max_lag_setup <- function() {
  m <- aci_dyad_model()
  ob <- as_obs(simulate(m, seed = 1, t_end = 3, dt = 0.01))
  list(m = m, ob = ob, init = list(mean = 0, cov = matrix(0.5)))
}

# The range of a cut table with its max_lag record removed: every number as
# the table gives it, and the status from the row alone.
uncapped_range <- function(tb, ...) {
  tb$meta$max_lag <- Inf
  suppressWarnings(aci_range(tb, ...))
}

cut_rows <- function(tb) {
  kept <- lengths(tb$rows)
  kept < length(tb$t) - seq_along(tb$t) + 1L & kept > tb$meta$max_lag
}

last_kept <- function(tb) vapply(tb$rows, function(r) r[length(r)], numeric(1))

expect_cut_status <- function(tb, floor_M = 1e-5, ...) {
  got <- suppressWarnings(aci_range(tb, ...))
  ref <- uncapped_range(tb, ...)
  expected <- as.character(ref$status)
  up <- cut_rows(tb) & last_kept(tb) > floor_M
  expected[up] <- "censored"
  expect_identical(as.character(got$status), expected)
  got$status <- ref$status <- NULL
  expect_identical(got, ref)
  invisible(up)
}

test_that("a row cut by max_lag above the strength floor is censored", {
  s <- max_lag_setup()
  tb <- lag_table(s$m, s$ob, init = s$init, max_lag = 5)
  r <- suppressWarnings(aci_range(tb))
  expect_identical(as.vector(table(r$status)), c(0L, 299L, 0L, 2L))
  expect_identical(as.character(r$status[100]), "censored")
  expect_gt(last_kept(tb)[100], 1e-5)
  up <- expect_cut_status(tb)
  expect_identical(sum(up), 295L)
  expect_cut_status(lag_table(s$m, s$ob, init = s$init, max_lag = 5, tol = 0))
  expect_cut_status(tb, min_M = 1e-3, epsilon = c(1e-3, 1e-5),
                    floor_M = 1e-3)
})

test_that("a row cut by max_lag at or below the floor keeps its status", {
  s <- max_lag_setup()
  tb <- lag_table(s$m, s$ob, init = s$init, max_lag = 200)
  full <- suppressWarnings(
    aci_range(lag_table(s$m, s$ob, init = s$init, tol = 0))
  )
  low <- cut_rows(tb) & last_kept(tb) <= 1e-5
  expect_identical(sum(cut_rows(tb)), 100L)
  r <- suppressWarnings(aci_range(tb))
  expect_identical(r$status, full$status)
  expect_true(any(r$status[low] == "resolved"))
  expect_cut_status(tb)
  # The comparison with the floor is strict, as for any other exceedance.
  p <- c(1, 0.5, 2e-5, 0, 0, 0, 0, 0)
  expect_identical(.cir_status(p, 1, 2e-5, kept = 3L), "resolved")
  expect_identical(.cir_status(p, 1, 1e-5, kept = 3L), "censored")
  expect_identical(.cir_status(p, 1, 1e-5), "resolved")
})

test_that("full and tol-truncated tables keep their status", {
  s <- max_lag_setup()
  full <- lag_table(s$m, s$ob, init = s$init, tol = 0)
  r_full <- suppressWarnings(aci_range(full))
  expect_identical(as.vector(table(r_full$status)), c(28L, 271L, 0L, 2L))
  expect_identical(
    suppressWarnings(aci_range(lag_table(s$m, s$ob, init = s$init,
                                         mode = "full")))$status,
    r_full$status
  )
  for (tb in list(lag_table(s$m, s$ob, init = s$init),
                  lag_table(s$m, s$ob, init = s$init, tol = 1e-2),
                  lag_table(s$m, s$ob, init = s$init, tol = 1e-2,
                            window = 1))) {
    expect_identical(suppressWarnings(aci_range(tb))$status,
                     uncapped_range(tb)$status)
  }
})

test_that("cut tables reach the same status through aci() and saveRDS()", {
  s <- max_lag_setup()
  tb <- lag_table(s$m, s$ob, init = s$init, max_lag = 5)
  want <- suppressWarnings(aci_range(tb))
  res <- aci(s$m, s$ob, init = s$init, table = tb)
  expect_identical(suppressWarnings(aci_range(res)), want)
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(tb, f)
  expect_identical(suppressWarnings(aci_range(readRDS(f))), want)
})

test_that("a cut row with several hidden components is censored", {
  z <- .compiled_cir_matrix_setup()
  tb <- lag_table(z$model, z$obs, init = z$init, max_lag = 3)
  up <- expect_cut_status(tb)
  expect_identical(sum(up), 8L)
  expect_lt(last_kept(tb)[9], 1e-5)
  pp <- aci_predprey_model()
  ob <- as_obs(simulate(pp, seed = 2, t_end = 2, dt = 0.01))
  up <- expect_cut_status(lag_table(pp, ob, max_lag = 10))
  expect_true(any(up))
})
