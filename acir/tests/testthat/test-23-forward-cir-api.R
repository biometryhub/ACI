# Forward-CIR reporting window, objective functionals, threshold-grid
# compatibility mode, subjective read-out convention and status vocabulary.
# Ledger entry C2b.

.c2b_setup <- function(t_end = 0.5, dt = 0.005, seed = 271) {
  model <- aci_dyad_model()
  sim <- simulate(model, seed = seed, t_end = t_end, dt = dt, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  bundle <- .compile_cgns_run(model, obs)
  filter <- .cgns_filter_compiled(bundle, init)
  list(model = model, obs = obs, init = init, bundle = bundle,
       filter = filter,
       table = .lag_table_compiled(bundle, mode = "forward", filter = filter))
}


# ---------------------------------------------------------------- window ----

test_that("a reporting window returns exactly the rows the whole record would", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  want <- c(3L, 17L, 40L, 61L)
  for (method in c("exact", "l1_linf")) {
    full <- suppressWarnings(aci_range(
      ds$table, method = method, epsilon = c(1e-6, 1e-4), min_M = 0))
    win <- suppressWarnings(aci_range(
      ds$table, method = method, epsilon = c(1e-6, 1e-4), min_M = 0,
      anchors = want))
    expect_identical(win$t, full$t[want])
    expect_identical(win$tau, full$tau[want])
    expect_identical(win$M, full$M[want])
    expect_identical(win$subjective, full$subjective[want, , drop = FALSE])
    expect_identical(win$tail_bound, full$tail_bound[want])
    expect_identical(win$status, full$status[want])
    expect_identical(win$interval, full$interval[want, , drop = FALSE])
    expect_identical(win$meta$anchors, want)
    expect_identical(full$meta$anchors, seq_along(full$t))
  }
})


test_that("the streamed window agrees with the streamed whole record", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  want <- c(1L, 9L, 55L, 99L)
  full <- suppressWarnings(.forward_cir_compiled(
    ds$bundle, filter = ds$filter, method = "exact",
    epsilon = c(1e-6, 1e-4), min_M = 0))
  win <- suppressWarnings(.forward_cir_compiled(
    ds$bundle, filter = ds$filter, method = "exact",
    epsilon = c(1e-6, 1e-4), min_M = 0, anchors = want))
  expect_identical(win$tau, full$tau[want])
  expect_identical(win$M, full$M[want])
  expect_identical(win$subjective, full$subjective[want, , drop = FALSE])
  expect_identical(win$tail_bound, full$tail_bound[want])
  expect_identical(win$status, full$status[want])
  # the two routes still agree with each other under the window
  win_tab <- suppressWarnings(aci_range(
    ds$table, method = "exact", epsilon = c(1e-6, 1e-4), min_M = 0,
    anchors = want))
  expect_equal(win$tau, win_tab$tau, tolerance = 1e-12)
  expect_identical(win$status, win_tab$status)
})


test_that("a window skips the primitives before its earliest anchor", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  late <- .compiled_forward_primitives(ds$bundle, ds$filter, from = 60L)
  all_ <- .compiled_forward_primitives(ds$bundle, ds$filter)
  expect_null(late$E[[59L]])                     # never formed
  expect_identical(late$E[[60L]], all_$E[[60L]]) # and identical from there on
  expect_identical(late$one_R[[70L]], all_$one_R[[70L]])
  expect_identical(late$T2[70L], all_$T2[70L])
  expect_identical(late$Ub[70L], all_$Ub[70L])
})


test_that("windows keep their order and reject anything that is not an index", {
  ds <- .c2b_setup()
  N1 <- length(ds$table$t)
  rev_ <- suppressWarnings(aci_range(ds$table, min_M = 0,
                                     anchors = c(40L, 5L, 12L)))
  expect_identical(rev_$meta$anchors, c(40L, 5L, 12L))
  expect_identical(rev_$t, ds$table$t[c(40L, 5L, 12L)])
  for (bad in list(0L, -1L, N1 + 1L, c(2L, 2L), 2.5, NA_integer_, "a",
                   integer(0)))
    expect_error(aci_range(ds$table, anchors = bad),
                 class = "aci_error_dims")
  expect_error(.forward_cir_compiled(ds$bundle, anchors = N1 + 1L),
               class = "aci_error_dims")
})


test_that("aci_range on an aci result carries the window through", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  res <- aci(ds$model, ds$obs, init = ds$init, keep = "paths")
  want <- c(4L, 31L)
  got <- suppressWarnings(aci_range(res, method = "l1_linf", min_M = 0,
                                    anchors = want))
  expect_identical(got$meta$anchors, want)
  expect_equal(
    got$tau,
    suppressWarnings(aci_range(ds$table, method = "l1_linf",
                               min_M = 0))$tau[want],
    tolerance = 1e-12)
})


# ------------------------------------------------- definitional objective ----

test_that("method = 'exact' is the layer-cake sum, not a quadrature of it", {
  set.seed(404); dt <- 0.02
  for (i in 1:15) {
    p <- abs(rnorm(41)) * exp(-seq(0, 4, length.out = 41)); p[41] <- 0
    q <- rev(cummax(rev(p)))
    got <- unname(.fwd_lengths(p, dt, "exact")["tau"])
    expect_identical(got, dt * sum(q) / max(p))
    # the same value whichever time-axis quadrature is named, because there
    # is no time-axis quadrature left to choose
    expect_identical(
      got, unname(.fwd_lengths(p, dt, "exact", "sum")["tau"]))
    expect_identical(
      got,
      unname(.fwd_lengths(p, dt, "exact", "simpson", "trapezoid")["tau"]))
  }
})


test_that("the exact objective is the threshold average of the counting range", {
  set.seed(77); dt <- 0.01
  p <- abs(rnorm(80)) * exp(-seq(0, 3, length.out = 80)); p[80] <- 0
  M <- max(p)
  # a fine midpoint rule over the threshold axis: unbiased for a step
  # function, so it converges on the closed form from both sides
  e <- (seq_len(200000L) - 0.5) * M / 200000
  last <- vapply(e, function(ee) { k <- which(p > ee)
    if (!length(k)) 0L else max(k) }, integer(1))
  expect_equal(unname(.fwd_lengths(p, dt, "exact")["tau"]),
               mean(last) * dt, tolerance = 1e-6)
})


# --------------------------------------------------- unequal-spacing rule ----

test_that(".simpson_xy integrates a quadratic exactly on an uneven grid", {
  f <- function(x) 3 * x^2 - 2 * x + 5
  F_ <- function(x) x^3 - x^2 + 5 * x
  for (n in c(7L, 8L, 21L, 32L)) {          # odd and even point counts
    x <- sort(c(0, 4, stats::runif(n - 2L, 0, 4)))
    x <- x[!duplicated(x)]
    exact <- F_(max(x)) - F_(min(x))
    expect_equal(.simpson_xy(x, f(x)), exact, tolerance = 1e-10)
  }
  # a log grid, which is what the compatibility mode actually integrates over
  x <- 10^seq(-6, 0.5, length.out = 513L)
  expect_equal(.simpson_xy(x, f(x)), F_(max(x)) - F_(min(x)),
               tolerance = 1e-10)
  # degenerate lengths follow simps.m onto the trapezoid
  expect_identical(.simpson_xy(c(1, 3), c(2, 4)), 6)
  expect_identical(.simpson_xy(1, 2), 0)
})


test_that(".simpson_xy on a unit grid reproduces the uniform rule", {
  set.seed(5)
  for (n in c(9L, 10L, 51L, 52L)) {
    y <- abs(rnorm(n)) + 0.1
    expect_equal(.simpson_xy(seq_len(n), y), .simpson(y),
                 tolerance = 1e-11)
  }
})


# ------------------------------------------ MATLAB compatibility quadrature --

test_that("the compatibility mode reproduces the reference definitional route", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  want <- c(5L, 20L, 45L)
  got <- suppressWarnings(aci_range(
    ds$table, method = "exact", quadrature = "matlab_eps_grid", min_M = 0,
    anchors = want))
  grid <- 10^seq(-6, 0.5, length.out = 513L)
  expect_equal(got$meta$epsilon_grid, grid)
  expect_identical(got$bound, "eps_grid_objective")
  # independent route: the read-out at every node, integrated with the
  # golden helper's own transcription of simps(x, y)
  for (i in seq_along(want)) {
    p <- lt_row(ds$table, want[i], pad = "zero")
    counts <- vapply(grid, function(e) { k <- which(p > e)
      if (!length(k)) 0 else max(k) }, numeric(1))
    expect_equal(got$tau[i],
                 .simps_xy(grid, counts * ds$table$dt) / max(p),
                 tolerance = 1e-12)
  }
  # and it is a different number from the exact average it approximates
  ex <- suppressWarnings(aci_range(ds$table, method = "exact", min_M = 0,
                                   anchors = want))
  expect_gt(max(abs(got$tau - ex$tau)), 1e-6)
})


test_that("a caller's own epsilon_grid is honoured and validated", {
  ds <- .c2b_setup()
  grid <- 10^seq(-5, 0, length.out = 65L)
  got <- suppressWarnings(aci_range(
    ds$table, method = "exact", quadrature = "matlab_eps_grid",
    epsilon_grid = grid, min_M = 0, anchors = 9L))
  expect_equal(got$meta$epsilon_grid, grid)
  p <- lt_row(ds$table, 9L, pad = "zero")
  counts <- vapply(grid, function(e) { k <- which(p > e)
    if (!length(k)) 0 else max(k) }, numeric(1))
  expect_equal(got$tau[1L],
               .simps_xy(grid, counts * ds$table$dt) / max(p),
               tolerance = 1e-12)
  for (bad in list(1, c(2, 1), c(1, 1, 2), c(-1, 1), c(1, NA), "a",
                   numeric(0)))
    expect_error(
      aci_range(ds$table, method = "exact",
                quadrature = "matlab_eps_grid", epsilon_grid = bad),
      class = "aci_error_dims")
})


test_that("the two epsilon roles cannot be taken from one vector", {
  ds <- .c2b_setup()
  # a grid supplied to any other mode is refused rather than quietly ignored
  expect_error(aci_range(ds$table, epsilon_grid = 10^seq(-5, 0, len = 9)),
               class = "aci_error_dims")
  expect_error(aci_range(ds$table, method = "l1_linf",
                         epsilon_grid = 10^seq(-5, 0, len = 9)),
               class = "aci_error_dims")
  expect_error(.forward_cir_compiled(ds$bundle,
                                     epsilon_grid = 10^seq(-5, 0, len = 9)),
               class = "aci_error_dims")
  # and the compatibility mode is defined for the definitional objective only
  expect_error(aci_range(ds$table, method = "l1_linf",
                         quadrature = "matlab_eps_grid"),
               class = "aci_error_dims")
  expect_error(.forward_cir_compiled(ds$bundle, method = "l1_linf",
                                     quadrature = "matlab_eps_grid"),
               class = "aci_error_dims")
  # the reporting thresholds keep their own argument and their own meaning
  got <- suppressWarnings(aci_range(
    ds$table, method = "exact", quadrature = "matlab_eps_grid",
    epsilon = c(1e-3, 1e-5), min_M = 0, anchors = 11L))
  expect_identical(dim(got$subjective), c(1L, 2L))
  expect_equal(got$meta$epsilon_grid, 10^seq(-6, 0.5, length.out = 513L))
})


# --------------------------------------------------- read-out convention ----

test_that("the default subjective read-out counts cells, one step above G.7", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  epsilon <- c(1e-1, 1e-3, 1e-5, 1e-9)
  cnt <- suppressWarnings(aci_range(ds$table, epsilon = epsilon, min_M = 0))
  lag <- suppressWarnings(aci_range(ds$table, epsilon = epsilon, min_M = 0,
                                    convention = "lag_time"))
  expect_identical(cnt$meta$convention, "count")
  expect_identical(lag$meta$convention, "lag_time")
  zero <- cnt$subjective == 0
  expect_gt(sum(zero), 0L)
  expect_gt(sum(!zero), 100L)
  # zero cells stay zero under both, and every other cell differs by exactly
  # one grid step
  expect_true(all(lag$subjective[zero] == 0))
  # Exactly one grid step, verified on the grid index rather than on the
  # difference: dt * k and dt * (k - 1) are each correctly rounded, so their
  # difference is dt only to within the rounding of the two products.
  k <- round(cnt$subjective / ds$table$dt)
  expect_identical(cnt$subjective, ds$table$dt * k)
  expect_identical(lag$subjective, ds$table$dt * pmax(k - 1, 0))
  expect_lt(max(abs(cnt$subjective[!zero] - lag$subjective[!zero] -
                      ds$table$dt)), 1e-15)   # measured 5.12e-17
  # the objective is not touched by the convention
  expect_identical(cnt$tau, lag$tau)
  # and the streamed route reads out the same way
  st <- suppressWarnings(.forward_cir_compiled(
    ds$bundle, filter = ds$filter, epsilon = epsilon, min_M = 0,
    convention = "lag_time"))
  expect_equal(st$subjective, lag$subjective, tolerance = 1e-12)
  expect_error(aci_range(ds$table, convention = "nope"))
})


# ------------------------------------------------------------- statuses ----

test_that("the status vocabulary reports the four row outcomes", {
  expect_identical(.cir_status(c(1, 0.5, 0, 0, 0, 0, 0), 1, 1e-5), "resolved")
  expect_identical(.cir_status(c(1, 0.5, 0.5, 0.5, 0), 1, 1e-5), "censored")
  expect_identical(.cir_status(c(1e-9, 1e-9, 0), 1e-9, 1e-5),
                   "below_threshold")
  expect_identical(.cir_status(c(1, 0), 1, 1e-5), "insufficient")
  expect_identical(.cir_status(1, 1, 1e-5), "insufficient")
  # insufficient outranks everything: it is the one outcome with no row
  expect_identical(.cir_status(c(1e-9, 0), 1e-9, 1e-5), "insufficient")
  # the record-end test is against what the record has left, not the last
  # cell, which is identically zero on a complete row
  expect_identical(.cir_status(c(1, 1, 1, 0, 0, 0, 0), 1, 1e-5), "resolved")
  expect_identical(.cir_status(c(1, 1, 1, 1, 1, 0, 0), 1, 1e-5), "censored")
})


test_that("status agrees with the mask and with both routes", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .c2b_setup()
  got <- suppressWarnings(aci_range(ds$table))
  expect_s3_class(got$status, "factor")
  expect_identical(levels(got$status),
                   c("resolved", "censored", "below_threshold",
                     "insufficient"))
  expect_identical(length(got$status), length(got$t))
  N1 <- length(ds$table$t)
  # the last two anchors have fewer than three observations left
  expect_identical(as.character(got$status[c(N1 - 1L, N1)]),
                   c("insufficient", "insufficient"))
  # below_threshold is exactly the masking predicate under the default floor
  masked <- is.na(got$tau) & got$M >= 1e-14
  expect_identical(got$status[masked],
                   factor(rep("below_threshold", sum(masked)),
                          levels = levels(got$status)))
  # a floor nothing can clear puts every full row below it
  all_low <- suppressWarnings(aci_range(ds$table, min_M = 1e6))
  expect_true(all(as.character(all_low$status[seq_len(N1 - 2L)]) ==
                    "below_threshold"))
  # the streamed route says the same thing
  st <- suppressWarnings(.forward_cir_compiled(ds$bundle,
                                               filter = ds$filter))
  expect_identical(st$status, got$status)
  # min_M = NULL turns off the mask, not the status
  none <- suppressWarnings(aci_range(ds$table, min_M = NULL))
  expect_identical(none$status, got$status)
  expect_false(anyNA(none$tau[!is.na(got$M) & got$M >= 1e-14]))
})


test_that("printing a range reports its statuses", {
  ds <- .c2b_setup()
  out <- paste(utils::capture.output(
    print(suppressWarnings(aci_range(ds$table, anchors = c(1L, 5L))))),
    collapse = "\n")
  expect_match(out, "status: ")
  expect_match(out, "cir_result")
})
