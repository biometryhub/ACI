## Ledger C3d: the strict covariance policy (decision D3).
##
## `regularize = "none"` is the default: any covariance that leaves the
## positive-definite cone inside a state recursion, a metric input or the
## likelihood stops the run with a classed condition naming the site, the grid
## index and the time. `regularize = "floor"` is the previous behaviour, and
## every floor it takes is recorded in the result's `meta$regularization`.
##
## The three firing cases below are the only ones in this workspace where a
## floor is taken at all: O1 counted 29 events in 4 371 482 floor-capable
## invocations across 47 runs, and every one of the 29 was in one of these
## three deliberately broken probes.

`%||%` <- function(a, b) if (is.null(a)) b else a

quiet <- function(expr)
  withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))

## Local option scope, so no block leaks a policy into the rest of the suite.
with_reg <- function(policy, expr) {
  old <- options(aci.regularize = policy)
  on.exit(options(old), add = TRUE)
  force(expr)
}

pol_dyad <- function(t_end = 3, dt = 0.001, seed = 11) {
  m <- aci_dyad_model()
  s <- simulate(m, seed = seed, t_end = t_end, dt = dt)
  list(model = m, obs = s$obs, init = list(mean = 2, cov = matrix(0.1, 1, 1)))
}

## The S4 probe: the oracle dyad path resampled to a 200x coarser grid, where
## the explicit Riccati step overshoots.
probe_S4 <- function() {
  d <- pol_dyad(t_end = 3, dt = 0.001, seed = 11)
  idx <- seq.int(1L, length(d$obs$t), by = 200L)
  d$obs <- observed_trajectory(d$obs$t[idx],
                               d$obs$x[idx, , drop = FALSE])
  d
}

## The S5 probe: a legitimate but vanishingly tight scalar prior.
probe_S5 <- function() {
  d <- pol_dyad(t_end = 3, dt = 0.001, seed = 11)
  d$init <- list(mean = 2, cov = matrix(1e-14, 1, 1))
  d
}


test_that("the policy argument resolves, validates and is option-backed", {
  expect_identical(.aci_regularize(NULL), "none")
  expect_identical(.aci_regularize("none"), "none")
  expect_identical(.aci_regularize("floor"), "floor")
  with_reg("floor", {
    expect_identical(.aci_regularize(NULL), "floor")
    ## an explicit argument always beats the option
    expect_identical(.aci_regularize("none"), "none")
  })
  for (bad in list("ridge", "", NA_character_, 1, TRUE, c("none", "floor")))
    expect_error(.aci_regularize(bad), class = "aci_error_dims")
})


test_that("every site id used in the sources is in the site register", {
  ids <- names(.ACI_COV_SITES)
  expect_true(length(ids) >= 11L)
  expect_false(anyDuplicated(ids) > 0L)
  for (id in ids) {
    info <- .ACI_COV_SITES[[id]]
    expect_type(info$role, "character")
    expect_type(info$where, "character")
    expect_true(nzchar(info$role) && nzchar(info$where))
  }
  ## each id produces a message that names its own quantity and its own place
  for (id in ids) {
    e <- tryCatch(.aci_stop_cov(id, 7L, 0.35, -1.5), error = function(e) e)
    expect_s3_class(e, "aci_error_covariance_not_spd")
    expect_s3_class(e, "aci_error_spd")
    expect_s3_class(e, "aci_error")
    expect_identical(e$site, id)
    expect_identical(e$role, .ACI_COV_SITES[[id]]$role)
    expect_identical(e$index, 7L)
    expect_identical(e$time, 0.35)
    expect_identical(e$value, -1.5)
    expect_match(conditionMessage(e), .ACI_COV_SITES[[id]]$role, fixed = TRUE)
    expect_match(conditionMessage(e), .ACI_COV_SITES[[id]]$where, fixed = TRUE)
    expect_match(conditionMessage(e), "index 7 (time 0.35)", fixed = TRUE)
  }
})


test_that("a clean run always carries a zero-event record, never NULL", {
  d <- pol_dyad(t_end = 0.4, dt = 0.002)
  shape <- function(r, policy) {
    expect_type(r, "list")
    expect_identical(names(r),
                     c("policy", "fired", "n_events", "eps", "sites"))
    expect_identical(r$policy, policy)
    expect_false(r$fired)
    expect_identical(r$n_events, 0L)
    expect_identical(r$eps, 1e-12)
    expect_s3_class(r$sites, "data.frame")
    expect_identical(nrow(r$sites), 0L)
    expect_identical(names(r$sites),
                     c("site", "role", "n", "first_index", "first_time",
                       "worst_value"))
  }
  for (p in c("none", "floor")) {
    f <- aci_filter(d$model, d$obs, init = d$init, regularize = p)
    s <- aci_smoother(d$model, d$obs, filter = f, init = d$init, regularize = p)
    a <- aci(d$model, d$obs, init = d$init, keep = "paths", regularize = p)
    tb <- lag_table(d$model, d$obs, mode = "forward", init = d$init,
                    regularize = p)
    cr <- aci_range(tb, min_M = 0)
    shape(f$meta$regularization, p)
    shape(s$meta$regularization, p)
    shape(a$meta$regularization, p)
    shape(a$paths$filter$meta$regularization, p)
    shape(a$paths$smoother$meta$regularization, p)
    shape(tb$meta$regularization, p)
    shape(cr$meta$regularization, p)
  }
})


test_that("the multivariate and implicit routes carry the record too", {
  m <- aci_enso_model(hidden = c("u", "hW"), variant = "aci_code")
  s <- simulate(m, seed = 8, t_end = 1, dt = 0.005)
  ini <- list(mean = c(0, 0), cov = 0.05 * diag(2))
  for (st in c("explicit", "implicit")) {
    f <- aci_filter(m, s$obs, init = ini, stepper = st)
    sm <- aci_smoother(m, s$obs, filter = f, init = ini, stepper = st)
    expect_identical(f$meta$regularization$policy, "none")
    expect_identical(sm$meta$regularization$policy, "none")
    expect_false(f$meta$regularization$fired)
    expect_false(sm$meta$regularization$fired)
  }
  a <- aci(m, s$obs, init = ini, keep = "paths", regularize = "floor")
  expect_identical(a$meta$regularization$policy, "floor")
  expect_false(a$meta$regularization$fired)
})


test_that("S4: a coarse explicit grid stops, and names where and when", {
  d <- probe_S4()
  e <- tryCatch(quiet(aci(d$model, d$obs, init = d$init)),
                error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_s3_class(e, "aci_error_spd")
  expect_identical(e$site, "filter_explicit")
  expect_identical(e$role, "filter covariance")
  expect_type(e$index, "integer")
  expect_gt(e$index, 1L)
  expect_equal(e$time, d$obs$t[e$index])
  expect_lt(e$value, 0)
  ## the same case under the opt-in completes and records what it did
  a <- quiet(aci(d$model, d$obs, init = d$init, keep = "paths",
                 regularize = "floor"))
  r <- a$meta$regularization
  expect_true(r$fired)
  expect_gt(r$n_events, 0L)
  expect_true(all(r$sites$site %in% names(.ACI_COV_SITES)))
  expect_true("filter_explicit" %in% r$sites$site)
  expect_identical(sum(r$sites$n), r$n_events)
  expect_true(all(r$sites$worst_value < 0))
  expect_true(all(r$sites$first_index >= 1L))
  expect_equal(r$sites$first_time, d$obs$t[r$sites$first_index])
  ## every floored covariance is back inside the cone
  expect_true(all(a$paths$filter$cov > 0))
  expect_true(all(a$paths$smoother$cov > 0))
  expect_true(all(is.finite(a$aci)))
  ## and the filter path's own record is the prefix of the call's
  expect_identical(a$paths$filter$meta$regularization$policy, "floor")
  expect_lte(a$paths$filter$meta$regularization$n_events, r$n_events)
})


test_that("S5: a vanishingly tight scalar prior stops in the smoother", {
  d <- probe_S5()
  e <- tryCatch(quiet(aci(d$model, d$obs, init = d$init)),
                error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_identical(e$site, "smoother_backward")
  expect_identical(e$role, "smoother covariance")
  expect_lt(e$value, 0)
  a <- quiet(aci(d$model, d$obs, init = d$init, keep = "paths",
                 regularize = "floor"))
  r <- a$meta$regularization
  expect_true(r$fired)
  expect_identical(r$sites$site, "smoother_backward")
  expect_true(all(a$paths$smoother$cov > 0))
})


test_that("S5b: a vanishingly tight matrix prior stops in the smoother", {
  m <- aci_enso_model(hidden = c("u", "hW"), variant = "aci_code")
  s <- simulate(m, seed = 8, t_end = 1, dt = 0.005)
  ini <- list(mean = c(0, 0), cov = 1e-14 * diag(2))
  e <- tryCatch(quiet(aci(m, s$obs, init = ini)), error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_true(e$site %in% c("smoother_backward", "filter_explicit"))
  expect_lt(e$value, 0)
  a <- quiet(aci(m, s$obs, init = ini, keep = "paths", regularize = "floor"))
  expect_true(a$meta$regularization$fired)
  ## a floored matrix covariance is SPD again at every stored slice
  cv <- a$paths$smoother$cov
  for (j in seq_len(dim(cv)[3L]))
    expect_false(is.null(tryCatch(chol(cv[, , j]), error = function(e) NULL)))
})


test_that("aci_smoother carries one record across a recomputed filter", {
  d <- probe_S5()
  ## the smoother recomputes the filter here, so both halves share one record
  sm <- quiet(aci_smoother(d$model, d$obs, init = d$init, regularize = "floor"))
  expect_true(sm$meta$regularization$fired)
  expect_identical(sm$meta$regularization$policy, "floor")
  expect_error(quiet(aci_smoother(d$model, d$obs, init = d$init)),
               class = "aci_error_covariance_not_spd")
  ## a filter supplied from outside leaves only the backward pass to record
  f <- aci_filter(d$model, d$obs, init = d$init)
  expect_false(f$meta$regularization$fired)
  sm2 <- quiet(aci_smoother(d$model, d$obs, filter = f, init = d$init,
                            regularize = "floor"))
  expect_identical(sm2$meta$regularization$n_events,
                   sm$meta$regularization$n_events)
})


test_that("a non-finite covariance stops under both policies", {
  ## Flooring projects a negative covariance back into the cone, so it is a
  ## remedy for leaving the cone and not for diverging out of the arithmetic.
  ## A non-finite variance therefore aborts whatever the policy says, and the
  ## condition still names where and when.
  for (p in c("none", "floor")) for (bad in c(-Inf, Inf, NaN, NA_real_)) {
    rec <- .aci_reg_new(p, c(0, 0.5, 1))
    rec$j <- 2L
    e <- tryCatch(.cov_guard_scalar(bad, rec, "filter_explicit"),
                  error = function(e) e)
    expect_s3_class(e, "aci_error_covariance_not_spd")
    expect_s3_class(e, "aci_error_spd")
    expect_identical(e$site, "filter_explicit")
    expect_identical(e$index, 2L)
    expect_identical(e$time, 0.5)
    expect_match(conditionMessage(e), "not recoverable by regularisation",
                 fixed = TRUE)
    ## an unrecoverable failure is not a floor event
    expect_identical(rec$n, 0L)
  }
  ## the matrix guard refuses a non-finite covariance the same way: strictly,
  ## and under "floor" through spd_floor()'s own finiteness contract
  S <- matrix(c(1, NA, NA, 1), 2, 2)
  rec_n <- .aci_reg_new("none", c(0, 1)); rec_n$j <- 2L
  expect_s3_class(tryCatch(.cov_guard(S, rec_n, "smoother_backward"),
                           error = function(e) e),
                  "aci_error_covariance_not_spd")
  rec_f <- .aci_reg_new("floor", c(0, 1)); rec_f$j <- 2L
  expect_error(.cov_guard(S, rec_f, "smoother_backward"),
               class = "aci_error_spd")
  ## a large but finite excursion is a floor, not a divergence
  m <- aci_model(
    Lx = function(t, x) matrix(1e8, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(0, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(0, 1, 1), k = 1, l = 1)
  ob <- observed_trajectory(c(0, 1, 2, 3), matrix(c(0, 1, 0, 1), ncol = 1))
  ini <- list(mean = 0, cov = matrix(1, 1, 1))
  e <- tryCatch(quiet(aci_filter(m, ob, init = ini)), error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_equal(e$value, -1e16)
  f <- quiet(aci_filter(m, ob, init = ini, regularize = "floor"))
  expect_true(all(f$cov > 0))
  expect_true(f$meta$regularization$fired)
})


test_that("the guards are the previous arithmetic when nothing fires", {
  set.seed(707)
  rec_n <- .aci_reg_new("none", c(0, 1, 2))
  rec_f <- .aci_reg_new("floor", c(0, 1, 2))
  rec_n$j <- rec_f$j <- 2L
  for (i in 1:12) {
    A <- matrix(stats::rnorm(9), 3, 3)
    S <- crossprod(A) + diag(3)
    expect_identical(.cov_guard(S, rec_n, "filter_explicit"), sym(S))
    expect_identical(.cov_guard(S, rec_f, "filter_explicit"), .sym_floor(S))
    g <- .cov_guard_chol(S, rec_n, "metric_reference")
    expect_identical(g$R, S)
    expect_identical(g$ch, chol(S))
    expect_identical(.cov_guard_scalar(S[1, 1], rec_n, "filter_explicit"),
                     S[1, 1])
  }
  expect_identical(rec_n$n, 0L)
  expect_identical(rec_f$n, 0L)
})


test_that("a fired matrix guard floors exactly as spd_floor() does", {
  S <- diag(c(1, -2, 3))
  rec_f <- .aci_reg_new("floor", c(0, 0.25))
  rec_f$j <- 2L
  expect_identical(.cov_guard(S, rec_f, "smoother_backward"), spd_floor(sym(S)))
  expect_identical(rec_f$n, 1L)
  fr <- .aci_reg_freeze(rec_f)
  expect_identical(fr$sites$site, "smoother_backward")
  expect_identical(fr$sites$first_index, 2L)
  expect_identical(fr$sites$first_time, 0.25)
  expect_equal(fr$sites$worst_value, -2)
  ## and the same matrix stops under the strict recorder
  rec_n <- .aci_reg_new("none", c(0, 0.25)); rec_n$j <- 2L
  e <- tryCatch(.cov_guard(S, rec_n, "smoother_backward"),
                error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_equal(e$value, -2)
  ## the reference-covariance guard floors through the same primitive
  g <- .cov_guard_chol(S, .aci_reg_new("floor", 0), "metric_reference")
  expect_identical(g$R, spd_floor(S))
  expect_identical(g$ch, chol(spd_floor(S)))
})


test_that("the governed solves reach the ladder only under 'floor'", {
  ## chol_solve() keeps safe_chol()'s jitter ladder wherever it is called
  ## without a recorder; the four recursion sites pass one, and there the
  ## ladder is a regularisation the policy governs.
  S <- matrix(c(1, 1, 1, 1), 2, 2)          # singular, ladder-recoverable
  B <- diag(2)
  plain <- chol_solve(S, B, "gxx")
  rec_f <- .aci_reg_new("floor", c(0, 1)); rec_f$j <- 2L
  expect_identical(chol_solve(S, B, "gxx", rec_f, "onelag_filter_cov"), plain)
  expect_identical(rec_f$n, 1L)
  expect_identical(.aci_reg_freeze(rec_f)$sites$site, "onelag_filter_cov")
  rec_n <- .aci_reg_new("none", c(0, 1)); rec_n$j <- 2L
  e <- tryCatch(chol_solve(S, B, "gxx", rec_n, "onelag_filter_cov"),
                error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_identical(e$site, "onelag_filter_cov")
  ## a matrix the plain factorisation accepts never reaches either branch
  P <- matrix(c(4, 1, 1, 3), 2, 2)
  expect_identical(chol_solve(P, B, "gxx", rec_n, "onelag_filter_cov"),
                   chol_solve(P, B, "gxx"))
  expect_identical(rec_n$n, 0L)
})


test_that("the conditioning Gram sites are deliberately not governed", {
  ## O1 measured 0 ladder entries in 226 413 conditioning-Gram invocations.
  ## Those sites invert an observed-noise Gram, a model input whose positive
  ## definiteness validate_cgns() already contracts, so their remedy is "fix
  ## the model", not "reduce dt". They keep safe_chol()'s ladder under every
  ## policy, and this block is what fails if that boundary is ever moved.
  G <- matrix(c(1, 1, 1, 1), 2, 2)
  for (p in c("none", "floor"))
    with_reg(p, {
      expect_identical(masked_ginv(G, 1:2), {
        M <- matrix(0, 2, 2)
        M[1:2, 1:2] <- chol_solve(G, diag(2), "masked Gram")
        M
      })
      expect_false(is.null(safe_chol(G, "gxx")))
      expect_identical(chol_solve(G, diag(2), "gxx"),
                       backsolve(safe_chol(G, "gxx"),
                                 forwardsolve(t(safe_chol(G, "gxx")), diag(2))))
    })
})


test_that("the exported regularisers keep their documented behaviour", {
  A <- matrix(c(1, 2, 2, 1), 2, 2)
  for (p in c("none", "floor"))
    with_reg(p, {
      expect_identical(spd_floor(matrix(c(2, 0, 0, 3), 2, 2)),
                       matrix(c(2, 0, 0, 3), 2, 2))
      fl <- spd_floor(A)
      expect_false(is.null(tryCatch(chol(fl), error = function(e) NULL)))
      expect_error(safe_chol(A, "x"), class = "aci_error_spd")
      expect_identical(.scalar_spd_floor(2), 2)
      expect_equal(.scalar_spd_floor(-2), 2e-12, tolerance = 0)
    })
})


test_that("a saved result's own policy governs its recomputed forward CIR", {
  d <- pol_dyad(t_end = 0.3, dt = 0.002)
  a <- aci(d$model, d$obs, init = d$init, keep = "paths")
  expect_identical(a$meta$regularization$policy, "none")
  ## the option changes after the result was produced; the result does not
  with_reg("floor", {
    cr <- aci_range(a, min_M = 0, anchors = 1:5)
    expect_identical(cr$meta$regularization$policy, "none")
  })
  af <- aci(d$model, d$obs, init = d$init, keep = "paths",
            regularize = "floor")
  crf <- aci_range(af, min_M = 0, anchors = 1:5)
  expect_identical(crf$meta$regularization$policy, "floor")
  ## and the two agree numerically, because nothing fires here
  cr0 <- aci_range(a, min_M = 0, anchors = 1:5)
  expect_identical(cr0$tau, crf$tau)
  expect_identical(cr0$M, crf$M)
})


test_that("the option is a default only, and the argument always wins", {
  d <- probe_S5()
  with_reg("floor", {
    a <- quiet(aci(d$model, d$obs, init = d$init))
    expect_identical(a$meta$regularization$policy, "floor")
    expect_true(a$meta$regularization$fired)
    expect_error(quiet(aci(d$model, d$obs, init = d$init,
                           regularize = "none")),
                 class = "aci_error_covariance_not_spd")
  })
  with_reg("none",
    expect_error(quiet(aci(d$model, d$obs, init = d$init)),
                 class = "aci_error_covariance_not_spd"))
  for (f in list(aci, aci_filter, aci_smoother, lag_table))
    expect_error(f(d$model, d$obs, init = d$init, regularize = "ridge"),
                 class = "aci_error_dims")
})
