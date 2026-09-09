# The regularization policy is audible, not only recorded ---------------------
#
# `regularize = "floor"` has always written a complete record onto the returned
# object. What it did not do was say so: a floor taken in the backward smoother
# raised no condition of any kind, because the stiffness diagnostic exists only
# in the explicit filter kernel. A run could therefore return a 6e21-nat
# information score with nothing to distinguish it from a converged one.
#
# The contract asserted here is that a call in which a floor fires raises
# exactly one `aci_warn_regularized`, that the printed object gains one line
# naming the count, and that a run which takes no floor is silent and prints
# exactly what it printed before. The deliberately unstable ENSO case checks
# reporting and consistency; its downstream values are not portable numerical
# references after an ill-conditioned covariance has been floored.

reg_capture <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(value = value,
       classes = vapply(warnings, function(w) class(w)[1L], character(1)),
       messages = vapply(warnings, conditionMessage, character(1)))
}

n_regularized <- function(cap) sum(cap$classes == "aci_warn_regularized")

## A site's only floor must record the value that the strict route refuses.
## The strict route stops before updating the regularization recorder.
expect_recorded_refusal <- function(row, error) {
  expect_identical(row$n, 1L)
  expect_s3_class(error, "aci_error_covariance_not_spd")
  expect_identical(row$site, error$site)
  expect_identical(row$first_index, error$index)
  expect_identical(row$first_time, error$time)
  expect_true(is.finite(error$value))
  expect_lt(error$value, 0)
  expect_equal(row$worst_value, error$value, tolerance = 1e-9)
}

## A legitimate but vanishingly tight scalar prior: the floor is taken in the
## backward smoother, which has no stiffness diagnostic, so before this change
## the call returned with an empty warning vector.
tight_scalar <- function() {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 11, t_end = 3, dt = 0.001)
  list(model = m, obs = s$obs, init = list(mean = 2, cov = matrix(1e-14, 1, 1)))
}

## The same dyad on a 200x coarser grid, where the explicit Riccati step
## overshoots: this one floors in the table and online routes too, and is cheap
## enough (16 points) to carry lag_table() and aci_range().
coarse_scalar <- function() {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 11, t_end = 3, dt = 0.001)
  idx <- seq.int(1L, length(s$obs$t), by = 200L)
  list(model = m,
       obs = observed_trajectory(s$obs$t[idx], s$obs$x[idx, , drop = FALSE]),
       init = list(mean = 2, cov = matrix(0.1, 1, 1)))
}

## The three-hidden-state ENSO partition under its automatic prior, which is
## about 380x too wide on two of its three components. This is the record the
## policy exists to make visible: it floors twice and scores about 6e21.
enso_record <- function() {
  m <- aci_enso_model(hidden = c("u", "hW", "tau"))
  ob <- as_obs(simulate(m, seed = 12, t_end = 4, dt = 0.005, burn_in = 0))
  list(model = m, obs = ob)
}


test_that("a floor taken only in the smoother is reported once", {
  d <- tight_scalar()
  a <- reg_capture(aci(d$model, d$obs, init = d$init, keep = "paths",
                       regularize = "floor"))
  expect_identical(n_regularized(a), 1L)
  ## the case's point: no other diagnostic fires here at all
  expect_false("aci_warn_riccati_stiff" %in% a$classes)
  expect_match(a$messages[a$classes == "aci_warn_regularized"],
               "smoother_backward", fixed = TRUE)
  expect_match(a$messages[a$classes == "aci_warn_regularized"],
               "meta$regularization", fixed = TRUE)
  r <- a$value$meta$regularization
  expect_true(r$fired)
  expect_identical(r$n_events, 1L)
  expect_identical(r$sites$site, "smoother_backward")
  expect_identical(r$sites$first_index, 1L)
  expect_equal(r$sites$first_time, 0)
  expect_equal(r$sites$worst_value, -4306.0395037011476, tolerance = 1e-9)
  expect_match(paste(capture.output(print(a$value)), collapse = "\n"),
               "regularized: 1 floor event(s)", fixed = TRUE)

  ## the forward pass alone takes no floor, so it stays silent
  f <- reg_capture(aci_filter(d$model, d$obs, init = d$init,
                              regularize = "floor"))
  expect_identical(n_regularized(f), 0L)
  expect_false(f$value$meta$regularization$fired)
  expect_identical(f$value$meta$regularization$n_events, 0L)
  expect_false(any(grepl("regularized",
                         capture.output(print(f$value)), fixed = TRUE)))

  ## and the smoother, which recomputes the filter, reports the same one event
  s <- reg_capture(aci_smoother(d$model, d$obs, init = d$init,
                                regularize = "floor"))
  expect_identical(n_regularized(s), 1L)
  expect_identical(s$value$meta$regularization$n_events, 1L)
})


test_that("the table and range verbs report and carry the same record", {
  d <- coarse_scalar()
  tb <- reg_capture(lag_table(d$model, d$obs, mode = "forward",
                              init = d$init, regularize = "floor"))
  expect_identical(n_regularized(tb), 1L)
  rt <- tb$value$meta$regularization
  expect_identical(rt$n_events, 5L)
  expect_identical(rt$sites$site, c("filter_explicit", "smoother_onelag"))
  expect_identical(rt$sites$n, c(2L, 3L))
  expect_match(paste(capture.output(print(tb$value)), collapse = "\n"),
               "regularized: 5 floor event(s)", fixed = TRUE)

  ## aci_range() reduces stored rows and runs no recursion, so it carries the
  ## table's record through to its own print without reporting a second time
  cr <- reg_capture(aci_range(tb$value, min_M = 0))
  expect_identical(n_regularized(cr), 0L)
  expect_identical(cr$value$meta$regularization$n_events, 5L)
  expect_match(paste(capture.output(print(cr$value)), collapse = "\n"),
               "regularized: 5 floor event(s)", fixed = TRUE)

  ## aci() on the same record floors six times, in the filter and the backward
  ## smoother, and still reports once
  a <- reg_capture(aci(d$model, d$obs, init = d$init, keep = "paths",
                       regularize = "floor"))
  expect_identical(n_regularized(a), 1L)
  expect_identical(a$value$meta$regularization$n_events, 6L)
  expect_identical(a$value$meta$regularization$sites$site,
                   c("filter_explicit", "smoother_backward"))

  ## a reused table runs no recursion, so it reports nothing new and shows the
  ## record its own rows were built under
  ra <- reg_capture(aci(d$model, d$obs, table = tb$value))
  expect_identical(n_regularized(ra), 0L)
  expect_true(ra$value$meta$regularization$fired)
  expect_identical(ra$value$meta$regularization$n_events, 5L)
  expect_match(paste(capture.output(print(ra$value)), collapse = "\n"),
               "regularized: 5 floor event(s)", fixed = TRUE)

  on <- reg_capture(aci_online(d$model, d$obs, lag = 5, init = d$init,
                               regularize = "floor"))
  expect_identical(n_regularized(on), 1L)
  expect_identical(on$value$meta$regularization$n_events, 5L)
  expect_identical(on$value$meta$regularization$sites$site,
                   c("filter_explicit", "smoother_onelag"))
  expect_match(paste(capture.output(print(on$value)), collapse = "\n"),
               "regularized: 5 floor event(s)", fixed = TRUE)
})


test_that("the floored ENSO record reports each failure and carries its metric", {
  d <- enso_record()
  a <- reg_capture(aci(d$model, d$obs, keep = "paths", regularize = "floor"))
  expect_identical(n_regularized(a), 1L)
  ## the two diagnostics this particular record already raised are retained
  expect_true("aci_warn_diffuse_init" %in% a$classes)
  expect_true("aci_warn_riccati_stiff" %in% a$classes)
  r <- a$value$meta$regularization
  expect_true(r$fired)
  expect_identical(r$n_events, 2L)
  expect_identical(r$sites$site, c("filter_explicit", "smoother_backward"))
  expect_identical(r$sites$n, c(1L, 1L))
  expect_identical(r$sites$first_index, c(2L, 2L))
  expect_equal(r$sites$first_time, c(0.005, 0.005))
  expect_equal(r$sites$worst_value[1L], -7.1493603591190933, tolerance = 1e-9)
  ## The first floor leaves a covariance with condition number about 1.7e11.
  ## Downstream values amplify rounding differences across BLAS/LAPACK builds.
  ## This bound identifies the inflated stress result; it is not an accuracy
  ## tolerance or evidence that flooring resolves the dynamics.
  expect_true(all(is.finite(a$value$aci)))
  expect_gt(max(a$value$aci), 1e20)
  expect_match(paste(capture.output(print(a$value)), collapse = "\n"),
               "regularized: 2 floor event(s)", fixed = TRUE)

  ## aci_metric() keeps no record of its own; the documented consequence is
  ## that the same 6e21 comes back from a bare data.frame
  km <- aci_metric(a$value$paths$smoother, a$value$paths$filter)
  expect_s3_class(km, "data.frame")
  expect_identical(names(km), c("t", "total", "signal", "dispersion"))
  expect_identical(km$t, a$value$t)
  expect_identical(km$total, a$value$aci)
  expect_identical(km$signal, a$value$signal)
  expect_identical(km$dispersion, a$value$dispersion)

  f <- reg_capture(aci_filter(d$model, d$obs, regularize = "floor"))
  expect_identical(n_regularized(f), 1L)
  expect_identical(f$value$meta$regularization$n_events, 1L)
  expect_identical(f$value$meta$regularization$sites$site, "filter_explicit")

  s <- reg_capture(aci_smoother(d$model, d$obs, regularize = "floor"))
  expect_identical(n_regularized(s), 1L)
  expect_identical(s$value$meta$regularization$n_events, 2L)
  expect_identical(s$value$meta$regularization$sites, r$sites)

  ## Reuse the same floored filter so the strict smoother reaches the second
  ## failure. Recomputing a strict filter would stop at the first one.
  es <- reg_capture(tryCatch(
    aci_smoother(d$model, d$obs, filter = f$value, regularize = "none"),
    error = identity))$value
  expect_recorded_refusal(r$sites[2L, , drop = FALSE], es)

  on <- reg_capture(aci_online(d$model, d$obs, lag = 20,
                               regularize = "floor"))
  expect_identical(n_regularized(on), 1L)
  ro <- on$value$meta$regularization
  expect_identical(ro$n_events, 2L)
  expect_identical(ro$sites$site, c("filter_explicit", "smoother_onelag"))
  expect_identical(ro$sites$n, c(1L, 1L))
  expect_identical(ro$sites$first_index, c(2L, 2L))
  expect_equal(ro$sites$first_time, c(0.005, 0.005))
  expect_equal(ro$sites$worst_value[1L], -7.1493603591190933, tolerance = 1e-9)
  eo <- reg_capture(tryCatch(
    aci_online(d$model, d$obs, lag = 20, filter = f$value,
               regularize = "none"), error = identity))$value
  expect_recorded_refusal(ro$sites[2L, , drop = FALSE], eo)
})


test_that("matrix floor records retain each site's minimum and first location", {
  rec <- .aci_reg_new("floor", c(0, 0.25, 0.5, 0.75))
  rec$j <- 2L
  ## The eigenvalues are 5 and -7; no measured package output is the reference.
  .cov_guard(matrix(c(-1, 6, 6, -1), 2, 2), rec, "smoother_backward")
  rec$j <- 3L
  .cov_guard(diag(c(2, -3)), rec, "smoother_backward")
  expect_equal(.aci_reg_freeze(rec)$sites$worst_value, -7, tolerance = 1e-12)
  rec$j <- 4L
  .cov_guard(diag(c(4, -5)), rec, "smoother_onelag")
  .cov_guard(diag(c(2, -9)), rec, "smoother_backward")
  r <- .aci_reg_freeze(rec)
  expect_true(r$fired)
  expect_identical(r$n_events, 4L)
  expect_identical(r$sites$site, c("smoother_backward", "smoother_onelag"))
  expect_identical(r$sites$n, c(3L, 1L))
  expect_identical(r$sites$first_index, c(2L, 4L))
  expect_identical(r$sites$first_time, c(0.25, 0.75))
  expect_equal(r$sites$worst_value, c(-9, -5), tolerance = 1e-12)
})


test_that("a run that takes no floor is silent and prints as before", {
  m <- aci_dyad_model()
  ob <- as_obs(simulate(m, seed = 11, t_end = 0.4, dt = 0.002))
  ini <- list(mean = 2, cov = matrix(0.1, 1, 1))
  for (policy in c("none", "floor")) {
    a <- reg_capture(aci(m, ob, init = ini, keep = "paths",
                         regularize = policy))
    expect_identical(n_regularized(a), 0L)
    expect_identical(a$value$meta$regularization$policy, policy)
    expect_false(a$value$meta$regularization$fired)
    expect_identical(a$value$meta$regularization$n_events, 0L)
    ## byte-identical to the pre-change capture, not a snapshot file
    expect_identical(
      capture.output(print(a$value)),
      "<aci_result> engine = cgns | peak ACI = 0.339 at t = 0.258")
    expect_identical(capture.output(print(a$value$paths$filter)),
                     "<da_path_gaussian> kind = filter, l = 1, N+1 = 201")
  }
  tb <- reg_capture(lag_table(m, ob, mode = "forward", init = ini,
                              regularize = "floor"))
  expect_identical(n_regularized(tb), 0L)
  expect_identical(
    capture.output(print(tb$value)),
    c("<lag_table> mode = forward, N+1 = 201, tol = 1e-08",
      paste("  mean retained lag: 100.0 steps;",
            "max heuristic tail estimate: 0.00e+00")))
  cr <- reg_capture(aci_range(tb$value, min_M = 0))
  expect_identical(n_regularized(cr), 0L)
  expect_identical(
    capture.output(print(cr$value)),
    c(paste("<cir_result> forward | method = exact (layer_cake_objective)",
            "| masked/NA: 1 of 201"),
      "  status: censored 199, insufficient 2"))
  on <- reg_capture(aci_online(m, ob, lag = 5, init = ini,
                               regularize = "floor"))
  expect_identical(n_regularized(on), 0L)
  expect_identical(capture.output(print(on$value)),
                   "<da_path_gaussian> kind = online, l = 1, N+1 = 201")
})
