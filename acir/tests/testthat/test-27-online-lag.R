# Public fixed-lag online smoother -------------------------------------------
#
# Every number asserted below was measured in
# scratch/shootout-0.1.0/overdrive/o4-online/evidence/ and re-measured on the
# integrated build in scratch/shootout-0.1.0/sessionc/w3c/. The three claims
# worth a test are the two boundaries and the semantics of the lag itself:
#
#   * lag 0 returns the filter moments unchanged, value for value;
#   * lag L at index j is the complete Theorem 3 posterior on the prefix ending
#     at j + L, and is NOT the one ending at j + L - 1 or j + L + 1;
#   * lag Inf is the complete Theorem 3 smoother, which is NOT aci_smoother().
#
# The last of those is the sentence aciR's own roxygen gets wrong; assert the
# separation rather than the identity.

.online_dyad_run <- function(n = 2001L) {
  signal <- read.csv(
    testthat::test_path("fixtures", "oracles", "dyad_signal_x.csv"),
    header = FALSE
  )[seq_len(n), ]
  obs <- observed_trajectory(signal$V1, matrix(signal$V2, ncol = 1L))
  list(
    model = aci_dyad_model(), obs = obs,
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L))
  )
}


test_that("lag zero returns the filter moments unchanged", {
  r <- .online_dyad_run()
  filt <- aci_filter(r$model, r$obs, init = r$init)
  on0 <- aci_online(r$model, r$obs, lag = 0, filter = filt, init = r$init)
  expect_identical(on0$kind, "online")
  expect_identical(on0$meta$scheme, "theorem3_discrete")
  expect_identical(on0$meta$route, "thmD1_online_window")
  expect_identical(on0$mean, filt$mean)
  expect_identical(on0$cov, filt$cov)
  expect_identical(on0$meta$lag_effective, integer(length(filt$t)))
  expect_false(on0$meta$saturated)
  expect_identical(on0$meta$regularization$policy, "none")
  expect_false(on0$meta$regularization$fired)
})


test_that("a lag of L admits the record through index j + L and no further", {
  r <- .online_dyad_run()
  filt <- aci_filter(r$model, r$obs, init = r$init)
  prefix_mean <- function(j, n) {
    ob <- observed_trajectory(r$obs$t[seq_len(n)],
                              r$obs$x[seq_len(n), , drop = FALSE])
    b <- .compile_cgns_run(r$model, ob)
    f <- .cgns_filter_compiled(b, init = r$init, stepper = "explicit",
                               nsub = 1L, validate = FALSE)
    .smoother_thmD1_compiled(b, f, validate = FALSE, warn_cost = FALSE)$mean[j, 1L]
  }
  for (j in c(51L, 777L)) for (L in c(1L, 4L, 8L)) {
    got <- aci_online(r$model, r$obs, lag = L, filter = filt,
                      init = r$init)$mean[j, 1L]
    expect_equal(got, prefix_mean(j, j + L), tolerance = 1e-12)
    # the neighbouring prefixes are a thousand times further away than the
    # agreement above, so this pins the off-by-one, not just the value
    expect_gt(abs(got - prefix_mean(j, j + L - 1L)), 1e-3)
    expect_gt(abs(got - prefix_mean(j, j + L + 1L)), 1e-3)
  }
})


test_that("[authors source] fixed-lag moments match the pinned online oracle", {
  ref <- read.csv(
    testthat::test_path("fixtures", "oracles", "cir_online_reference.csv")
  )
  r <- .online_dyad_run()
  b <- .compile_cgns_run(r$model, r$obs)
  f <- .cgns_filter_compiled(b, init = r$init, stepper = "explicit", nsub = 1L,
                             validate = FALSE)
  aux <- .online_aux_compiled(b, f)
  got <- .online_at_compiled(b, f, aux, ref$j, ref$n - ref$j)
  expect_gt(max(ref$n - ref$j), 100L)
  expect_lt(max(abs(got$mean[, 1L] - ref$online_mean)), 1e-6)
  expect_lt(max(abs(got$cov[1L, 1L, ] - ref$online_cov)), 1e-6)
})


test_that("full lag is the complete Theorem 3 smoother, not the backward ODE", {
  r <- .online_dyad_run()
  filt <- aci_filter(r$model, r$obs, init = r$init)
  b <- .compile_cgns_run(r$model, r$obs)
  f <- .cgns_filter_compiled(b, init = r$init, stepper = "explicit", nsub = 1L,
                             validate = FALSE)
  thm <- .smoother_thmD1_compiled(b, f, validate = FALSE, warn_cost = FALSE)
  on_inf <- aci_online(r$model, r$obs, lag = Inf, filter = filt, init = r$init)
  expect_true(isTRUE(on_inf$meta$saturated))
  expect_identical(on_inf$meta$route, "thmD1_backward")
  expect_equal(on_inf$mean, thm$mean, tolerance = 1e-12)
  expect_equal(on_inf$cov, thm$cov, tolerance = 1e-12)

  # the scheme gap the documentation states, asserted as a separation
  ode <- aci_smoother(r$model, r$obs, filter = filt, init = r$init)
  expect_identical(ode$meta$scheme, "backward_ode_euler")
  expect_identical(thm$meta$scheme, "theorem3_discrete")
  expect_gt(max(abs(on_inf$mean - ode$mean)), 1e-3)

  # a saturating finite lag reaches the same fixed point by the other route
  forced <- .da_online_compiled(b, f, lag = b$N, method = "window")
  expect_identical(forced$meta$route, "thmD1_online_window")
  expect_equal(forced$mean, thm$mean, tolerance = 1e-10)
  expect_equal(forced$cov, thm$cov, tolerance = 1e-10)
})


test_that("the fixed-lag path is the object the lag table is built from", {
  r <- .online_dyad_run(801L)
  b <- .compile_cgns_run(r$model, r$obs)
  f <- .cgns_filter_compiled(b, init = r$init, stepper = "explicit", nsub = 1L,
                             validate = FALSE)
  aux <- .online_aux_compiled(b, f)
  thm <- .smoother_thmD1_compiled(b, f, validate = FALSE, warn_cost = FALSE)
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  tb <- lag_table(r$model, r$obs, mode = "forward", tol = 0, init = r$init)
  expect_identical(tb$meta$scheme, "theorem3_discrete")
  for (j in c(11L, 101L, 401L)) {
    cells <- lt_row(tb, j)
    K <- min(length(cells), 40L)
    got <- .online_at_compiled(b, f, aux, rep(j, K), 0:(K - 1L))
    kl <- vapply(seq_len(K), function(i) aci_metric_pair(
      thm$mean[j, ], matrix(thm$cov[, , j], b$l, b$l),
      got$mean[i, ], matrix(got$cov[, , i], b$l, b$l), decompose = FALSE
    ), numeric(1L))
    expect_lt(max(abs(kl - cells[seq_len(K)])), 1e-12)
  }
})


test_that("the online path is refused where a complete smoother is required", {
  r <- .online_dyad_run(401L)
  filt <- aci_filter(r$model, r$obs, init = r$init)
  on5 <- aci_online(r$model, r$obs, lag = 5, filter = filt, init = r$init)
  expect_error(
    lag_table(r$model, r$obs, mode = "forward", tol = 0, init = r$init,
              smoother = on5),
    class = "aci_error_dims"
  )
  expect_error(
    aci_smoother(r$model, r$obs, init = r$init, filter = on5),
    class = "aci_error_dims"
  )
  expect_error(
    aci_online(r$model, r$obs, lag = 2, init = r$init, filter = on5),
    class = "aci_error_dims"
  )
  # but it is a Gaussian path, so the divergence functionals accept it
  full <- aci_online(r$model, r$obs, lag = Inf, filter = filt, init = r$init)
  kl <- aci_metric(full, on5, decompose = FALSE)
  expect_true(all(kl$total >= 0))
  expect_gt(max(kl$total), 0)
  # and so do the read-outs that only need Gaussian moments
  expect_s3_class(as.data.frame(on5), "data.frame")
  expect_identical(nrow(as.data.frame(on5)), length(on5$t))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(on5), on5)
})


test_that("aci_online validates its lag and its filter", {
  r <- .online_dyad_run(201L)
  expect_error(aci_online(r$model, r$obs), class = "aci_error_dims")
  expect_error(aci_online(r$model, r$obs, lag = -1), class = "aci_error_dims")
  expect_error(aci_online(r$model, r$obs, lag = 2.5), class = "aci_error_dims")
  expect_error(aci_online(r$model, r$obs, lag = NA), class = "aci_error_dims")
  expect_error(aci_online(r$model, r$obs, lag = 1, bogus = 1),
               class = "aci_error_dims")
  expect_error(aci_online(r$model, r$obs, lag = 1, force_validate = NA),
               class = "aci_error_dims")
  expect_error(aci_online(r$model, r$obs, lag = 1, regularize = "clip"),
               class = "aci_error_dims")
  expect_error(
    aci_online(r$model, r$obs, lag = 1,
               filter = aci_smoother(r$model, r$obs, init = r$init)),
    class = "aci_error_dims"
  )
  # an implicit or sub-stepped filter is not the discretization Theorem 3 is
  # exact for, and is refused rather than silently accepted
  fi <- aci_filter(r$model, r$obs, init = r$init, stepper = "implicit")
  expect_error(aci_online(r$model, r$obs, lag = 1, filter = fi),
               class = "aci_error_stepper")
  expect_error(
    aci_online(r$model, r$obs, lag = 1,
               filter = aci_filter(r$model, r$obs, init = r$init),
               init = list(mean = 0, cov = matrix(1, 1L, 1L))),
    class = "aci_error_dims"
  )
  expect_error(aci_online(structure(list(k = 1L, l = 1L),
                                   class = "stochastic_model"),
                         r$obs, lag = 1),
               class = "aci_error_not_implemented")
})


test_that("[independent MATLAB] the matrix path matches the mv online oracle", {
  signal <- read.csv(testthat::test_path("fixtures", "oracles", "mv_signal.csv"))
  ref <- read.csv(
    testthat::test_path("fixtures", "oracles", "mv_online_reference.csv")
  )
  model <- .compiled_oracle_mv_model()
  obs <- observed_trajectory(signal$t, as.matrix(signal[, c("x1", "x2")]))
  b <- .compile_cgns_run(model, obs)
  f <- .cgns_filter_compiled(
    b, init = list(mean = c(0.8, 0.2), cov = 0.2 * diag(2L)),
    stepper = "explicit", nsub = 1L, validate = FALSE
  )
  got <- .online_at_compiled(b, f, .online_aux_compiled(b, f),
                             ref$j, ref$n - ref$j)
  got_cov <- t(vapply(seq_len(nrow(ref)), function(i)
    c(got$cov[1L, 1L, i], got$cov[1L, 2L, i], got$cov[2L, 2L, i]), numeric(3L)))
  expect_lt(max(abs(got$mean - as.matrix(ref[, c("om1", "om2")]))), 1e-6)
  expect_lt(
    max(abs(got_cov - as.matrix(ref[, c("oc11", "oc12", "oc22")]))), 1e-6
  )
})


test_that("the sliding window and the direct accumulator agree on a matrix model", {
  # l = 3 with a correlated, state-dependent Gram, so the window's associative
  # composition is exercised on genuinely non-commuting E products.
  signal <- read.csv(
    testthat::test_path("fixtures", "oracles", "enso_signal.csv")
  )[seq_len(801L), ]
  model <- aci_enso_model(hidden = c("u", "hW", "tau"), variant = "aci_code")
  obs <- observed_trajectory(signal$t,
                             as.matrix(signal[, c("T_C", "T_E", "I")]))
  init <- list(mean = c(6.9136e-04, -0.0028, -0.0256), cov = 0.01 * diag(3L))
  b <- .compile_cgns_run(model, obs)
  f <- .cgns_filter_compiled(b, init = init, stepper = "explicit", nsub = 1L,
                             validate = FALSE)
  aux <- .online_aux_compiled(b, f)
  anchors <- seq_len(b$N1)
  for (L in c(1L, 7L, 60L)) {
    win <- .online_window_compiled(b, f, aux, L)
    at <- .online_at_compiled(b, f, aux, anchors, L)
    expect_equal(win$mean, at$mean, tolerance = 1e-12)
    expect_equal(win$cov, at$cov, tolerance = 1e-12)
    expect_identical(win$lag_effective, at$lag_effective)
    # the window really did carry L steps of information away from the end
    expect_identical(win$lag_effective[1L], L)
    expect_identical(win$lag_effective[b$N1], 0L)
  }
  # and at saturation both agree with the backward sweep on the whole record
  thm <- .smoother_thmD1_compiled(b, f, validate = FALSE, warn_cost = FALSE)
  win <- .online_window_compiled(b, f, aux, b$N)
  expect_equal(win$mean, thm$mean, tolerance = 1e-10)
  expect_equal(win$cov, thm$cov, tolerance = 1e-10)
  # the accumulated covariances never left the positive-definite cone
  expect_gt(min(apply(win$cov, 3L, function(R)
    min(eigen(R, symmetric = TRUE, only.values = TRUE)$values))), 0)
})


test_that("aci_online carries the covariance policy of its call", {
  r <- .online_dyad_run(401L)
  on <- aci_online(r$model, r$obs, lag = 3, init = r$init)
  expect_identical(on$meta$regularization$policy, "none")
  expect_identical(on$meta$regularization$n_events, 0L)
  on_fl <- aci_online(r$model, r$obs, lag = 3, init = r$init,
                      regularize = "floor")
  expect_identical(on_fl$meta$regularization$policy, "floor")
  # nothing fires on this record, so the two policies are the same numbers
  expect_identical(on$mean, on_fl$mean)
  expect_identical(on$cov, on_fl$cov)
  # the option is the default, and the call argument overrides it
  old <- options(aci.regularize = "floor")
  on.exit(options(old), add = TRUE)
  expect_identical(
    aci_online(r$model, r$obs, lag = 3, init = r$init)$meta$regularization$policy,
    "floor"
  )
  expect_identical(
    aci_online(r$model, r$obs, lag = 3, init = r$init,
               regularize = "none")$meta$regularization$policy,
    "none"
  )
})


test_that("aci() reports the discretization scheme, not the route", {
  r <- .online_dyad_run(401L)
  a <- aci(r$model, r$obs, init = r$init, keep = "paths")
  expect_identical(a$meta$smoother_scheme, "backward_ode_euler")
  expect_identical(a$paths$smoother$meta$scheme, "backward_ode_euler")
  expect_identical(a$paths$smoother$meta$route, "backward_ode")
  # a reused table was built under the other scheme and says so
  at <- aci(r$model, r$obs, init = r$init, keep = "table")
  expect_identical(at$table$meta$scheme, "theorem3_discrete")
  expect_identical(at$table$meta$reference_smoother, "thmD1_online_complete")
  expect_identical(
    aci(r$model, r$obs, init = r$init, table = at$table)$meta$smoother_scheme,
    "theorem3_discrete"
  )
  # a table saved before meta$scheme existed still answers with its route tag
  legacy <- at$table
  legacy$meta$scheme <- NULL
  expect_identical(
    aci(r$model, r$obs, init = r$init, table = legacy)$meta$smoother_scheme,
    "thmD1_online_complete"
  )
})
