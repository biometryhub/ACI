# The realised observation-noise Gram contract --------------------------------
#
# The model contract on the observation-noise Gram, enforced on the realised
# path rather than at the constructor's five probe points: every public verb,
# the scalar and matrix routes, the specialised dyad realiser, both conditional
# routes, and cached and cache-disabled realisation.

quietly <- function(expr) withCallingHandlers(
  expr, warning = function(w) invokeRestart("muffleWarning"))

## Sx1 collapses at t = 0.01, which is not one of the constructor's probe
## times ((i - 1) * sqrt(2)), so the model is accepted at construction.
gram_scalar_model <- function(Lx_value) aci_model(
  Lx = function(t, x) matrix(Lx_value, 1, 1),
  fx = function(t, x) 0,
  Ly = function(t, x) matrix(-0.5, 1, 1),
  fy = function(t, x) 0,
  Sx1 = function(t, x) matrix(if (abs(t - 0.01) < 1e-12) 0 else 1, 1, 1),
  Sy2 = function(t, x) matrix(0.5, 1, 1),
  k = 1, l = 1)

gram_scalar_record <- function()
  observed_trajectory(c(0, 0.01, 0.02, 0.03),
                      matrix(c(0, 0.1, 0.2, 0.3), ncol = 1))

gram_scalar_init <- list(mean = 0, cov = matrix(1, 1, 1))


test_that(paste0("a Gram that collapses between probe points is refused, ",
                 "whatever the stepper or the covariance policy"), {
  ob <- gram_scalar_record()
  for (stepper in c("explicit", "implicit"))
    for (regularize in c("none", "floor")) {
      info <- paste(stepper, regularize)
      e <- expect_error(
        quietly(aci_filter(gram_scalar_model(1e-6), ob,
                           init = gram_scalar_init, stepper = stepper,
                           regularize = regularize)),
        class = "aci_error_gram_path", info = info)
      expect_s3_class(e, "aci_error_gram")
      expect_identical(e$index, 2L, info = info)
      expect_identical(e$time, 0.01, info = info)
      expect_identical(e$rcond, 0, info = info)
      ## The remedies for integration instability are named as remedies that
      ## do not apply, never offered.
      expect_false(grepl("raise nsub", conditionMessage(e), fixed = TRUE))
      expect_false(grepl("call with regularize", conditionMessage(e),
                         fixed = TRUE))
      expect_match(conditionMessage(e), "observation-model violation")
      expect_match(conditionMessage(e), "do not repair it")
    }
})


test_that("every public verb on the same record refuses identically", {
  m <- gram_scalar_model(1e-6)
  ob <- gram_scalar_record()
  ini <- gram_scalar_init
  calls <- list(
    aci = function() aci(m, ob, init = ini),
    aci_filter = function() aci_filter(m, ob, init = ini),
    aci_smoother = function() aci_smoother(m, ob, init = ini),
    online_lag2 = function() aci_online(m, ob, lag = 2, init = ini),
    online_full = function() aci_online(m, ob, lag = Inf, init = ini),
    lag_table = function() lag_table(m, ob, mode = "forward", init = ini)
  )
  for (nm in names(calls)) {
    e <- expect_error(quietly(calls[[nm]]()), class = "aci_error_gram_path",
                      info = nm)
    expect_identical(e$index, 2L, info = nm)
    expect_identical(e$time, 0.01, info = nm)
  }
})


test_that(paste0("strong coupling reports the Gram violation, not the ",
                 "covariance policy"), {
  ob <- gram_scalar_record()
  for (regularize in c("none", "floor"))
    for (nsub in list(1L, 4L)) {
      e <- expect_error(
        quietly(aci_filter(gram_scalar_model(1), ob, init = gram_scalar_init,
                           regularize = regularize, nsub = nsub)),
        class = "aci_error_gram_path")
      expect_false(inherits(e, "aci_error_covariance_not_spd"))
      expect_identical(e$index, 2L)
    }
})


test_that("matrix rank loss in the interior is refused at its own index", {
  m <- aci_model(
    Lx = function(t, x) diag(1e-6, 2),
    fx = function(t, x) c(0, 0),
    Ly = function(t, x) diag(-0.5, 2),
    fy = function(t, x) c(0, 0),
    Sx1 = function(t, x) if (abs(t - 0.01) < 1e-12)
      matrix(c(1, 1, 0, 0), 2, 2) else diag(2),
    Sy2 = function(t, x) diag(0.5, 2),
    k = 2, l = 2)
  ob <- observed_trajectory(
    c(0, 0.01, 0.02, 0.03),
    cbind(c(0, 0.1, 0.2, 0.3), c(0, -0.05, -0.1, -0.15)))
  ini <- list(mean = c(0, 0), cov = diag(2))
  for (stepper in c("explicit", "implicit")) {
    e <- expect_error(quietly(aci_filter(m, ob, init = ini, stepper = stepper)),
                      class = "aci_error_gram_path", info = stepper)
    expect_identical(e$index, 2L, info = stepper)
    expect_identical(e$time, 0.01, info = stepper)
  }
})


test_that("cached and cache-disabled realisation refuse identically", {
  m <- gram_scalar_model(1e-6)
  ob <- gram_scalar_record()
  ini <- gram_scalar_init
  for (cache in c(TRUE, FALSE)) {
    old <- options(aci.realiser_cache = cache)
    on.exit(options(old), add = TRUE)
    .realise_cache_clear()
    for (attempt in 1:2) {
      e <- expect_error(quietly(aci_filter(m, ob, init = ini)),
                        class = "aci_error_gram_path",
                        info = paste(cache, attempt))
      expect_identical(e$index, 2L)
    }
    options(old)
  }
  .realise_cache_clear()
})


## A two-channel model whose noise collapses on one channel exactly where that
## channel's observation is zero. The constructor's probe states are
## sin(seq_len(k) * (i + 0.25)), none of which is zero, so the model is
## accepted at construction and the collapse is a property of the record.
conditional_gram_model <- function(zero_channel) aci_model(
  Lx = function(t, x) matrix(c(0.7, -0.25), 2, 1),
  fx = function(t, x) c(-0.3 * x[1], -0.2 * x[2]),
  Ly = function(t, x) matrix(-0.8, 1, 1),
  fy = function(t, x) 0.15 * cos(t),
  Sx1 = function(t, x) {
    d <- c(0.7, 0.6)
    if (abs(x[zero_channel]) < 1e-12) d[zero_channel] <- 0
    diag(d, 2)
  },
  Sx2 = function(t, x) matrix(0, 2, 1),
  Sy1 = function(t, x) matrix(0, 1, 2),
  Sy2 = function(t, x) matrix(0.8, 1, 1),
  k = 2, l = 1)

conditional_gram_record <- function(zero_channel, zero_index) {
  tt <- seq(0, 0.02, by = 0.005)
  xx <- cbind(0.3 + 0.01 * seq_along(tt), -0.2 + 0.01 * seq_along(tt))
  xx[zero_index, zero_channel] <- 0
  observed_trajectory(tt, xx)
}

conditional_gram_init <- list(mean = 0.2, cov = matrix(0.4, 1, 1))


test_that("the check tests the whole realised slice, whatever the route", {
  ## Slice 1 is degenerate on channel 2 only; channel 1 is the target.
  k <- 2L
  N1 <- 4L
  gxx <- array(0, c(k, k, N1))
  for (j in seq_len(N1)) gxx[, , j] <- diag(c(0.49, 0.36))
  gxx[, , 1L] <- diag(c(0.49, 0))
  tt <- c(0, 0.005, 0.01, 0.015)
  N <- N1 - 1L
  ## Unmasked, masked, masked under the ACI_code first-slice convention, and
  ## masked on an empty target: the tested matrix is the whole slice in every
  ## case, so a degeneracy outside the target block is refused at the same
  ## index and the same time as it is unmasked.
  e <- expect_error(.compiled_precision_path(gxx, N, tgrid = tt),
                    class = "aci_error_gram_path")
  expect_identical(e$index, 1L)
  expect_identical(e$time, 0)
  e <- expect_error(
    .compiled_precision_path(gxx, N, target = 1L, tgrid = tt),
    class = "aci_error_gram_path")
  expect_identical(e$index, 1L)
  expect_identical(e$time, 0)
  e <- expect_error(
    .compiled_precision_path(gxx, N, target = 1L, first_step = "matlab",
                             tgrid = tt),
    class = "aci_error_gram_path")
  expect_identical(e$index, 1L)
  expect_identical(e$time, 0)
  e <- expect_error(
    .compiled_precision_path(gxx, N, target = integer(0), tgrid = tt),
    class = "aci_error_gram_path")
  expect_identical(e$index, 1L)
  expect_identical(e$time, 0)
  ## A degeneracy inside the target block is refused at its own index, and
  ## the message names no block.
  gxx2 <- gxx
  gxx2[, , 1L] <- diag(c(0.49, 0.36))
  gxx2[, , 3L] <- diag(c(0, 0.36))
  e <- expect_error(
    .compiled_precision_path(gxx2, N, target = 1L, tgrid = tt),
    class = "aci_error_gram_path")
  expect_identical(e$index, 3L)
  expect_identical(e$time, 0.01)
  expect_false(grepl("target block", conditionMessage(e), fixed = TRUE))
  ## The terminal slice is never inverted and is never tested.
  gxx3 <- gxx
  gxx3[, , 1L] <- diag(c(0.49, 0.36))
  gxx3[, , N1] <- diag(c(0, 0))
  expect_silent(.compiled_precision_path(gxx3, N, tgrid = tt))
  ## An empty target on an admissible record still inverts nothing.
  expect_identical(.compiled_precision_path(gxx3, N, target = integer(0),
                                            tgrid = tt),
                   array(0, c(k, k, N)))
})


test_that("the conditional routes refuse a degenerate target block", {
  ini <- conditional_gram_init
  m_target <- conditional_gram_model(1L)
  ob_target <- conditional_gram_record(1L, 3L)
  ## The whole realised Gram is tested at every interval start, so the refusal
  ## is at the same index and the same time whether the unconditioned
  ## precision path is realised in the cache or the conditional branch derives
  ## its own.
  for (cache in c(TRUE, FALSE)) {
    old <- options(aci.realiser_cache = cache)
    on.exit(options(old), add = TRUE)
    .realise_cache_clear()
    for (method in c("mask", "reduce")) {
      info <- paste(cache, method)
      e <- expect_error(
        quietly(aci_filter(m_target, ob_target, init = ini,
                           conditional = aci_conditional(2, method))),
        class = "aci_error_gram_path", info = info)
      expect_identical(e$index, 3L, info = info)
      expect_identical(e$time, 0.01, info = info)
      expect_false(grepl("target block", conditionMessage(e), fixed = TRUE),
                   info = info)
    }
    options(old)
  }
  .realise_cache_clear()
})


test_that("a degenerate non-target block is refused whatever the cache state", {
  ## The whole slice is what the model contract is on, so a degeneracy on a
  ## non-target channel is refused on the conditional routes exactly as it is
  ## unconditioned, and switching the realisation cache off does not make the
  ## same record admissible.
  ini <- conditional_gram_init
  m_nontarget <- conditional_gram_model(2L)
  ob_nontarget <- conditional_gram_record(2L, 1L)
  .realise_cache_clear()
  e <- expect_error(
    quietly(aci_filter(m_nontarget, ob_nontarget, init = ini,
                       conditional = aci_conditional(2, "mask"))),
    class = "aci_error_gram_path")
  expect_identical(e$index, 1L)
  e <- expect_error(quietly(aci_filter(m_nontarget, ob_nontarget, init = ini)),
                    class = "aci_error_gram_path")
  expect_identical(e$index, 1L)
  old <- options(aci.realiser_cache = FALSE)
  on.exit(options(old), add = TRUE)
  .realise_cache_clear()
  for (first_step in c("uniform", "matlab")) {
    e <- expect_error(
      quietly(aci_filter(m_nontarget, ob_nontarget, init = ini,
                         conditional = aci_conditional(
                           2, "mask", first_step = first_step))),
      class = "aci_error_gram_path", info = first_step)
    expect_identical(e$index, 1L, info = first_step)
    expect_identical(e$time, 0, info = first_step)
  }
  options(old)
  .realise_cache_clear()
})


test_that("every route and both cache states refuse the same record", {
  ## method = "reduce" drops the non-target channels before it inverts, so
  ## nothing downstream of the reduction can see this degeneracy; with the
  ## cache off nothing upstream would either. The six combinations are
  ## asserted together because the claim is that neither the route nor the
  ## cache state changes which records are admissible.
  ini <- conditional_gram_init
  m <- conditional_gram_model(2L)
  ob <- conditional_gram_record(2L, 1L)
  for (cache in c(TRUE, FALSE)) {
    old <- options(aci.realiser_cache = cache)
    on.exit(options(old), add = TRUE)
    .realise_cache_clear()
    for (method in list("mask", "reduce", NULL)) {
      info <- paste(cache, if (is.null(method)) "none" else method)
      cond <- if (is.null(method)) NULL else aci_conditional(2, method)
      e <- expect_error(
        quietly(aci_filter(m, ob, init = ini, conditional = cond)),
        class = "aci_error_gram_path", info = info)
      expect_identical(e$index, 1L, info = info)
      expect_identical(e$time, 0, info = info)
    }
    options(old)
  }
  .realise_cache_clear()
})


test_that(paste0("the specialised dyad realiser carries the same test, ",
                 "with the same k = 1 limitation"), {
  ## rcond of a non-zero 1 by 1 Gram is 1, so a vanishingly small but non-zero
  ## scalar observation noise is accepted here exactly as the constructor
  ## probe accepts it: the rule is relative conditioning, not a noise floor.
  ## The constructor already requires s_x > 0 and the realiser descriptor is
  ## locked, so the same test inside .compile_dyad_cgns() is defence in depth
  ## with no reachable failing input; what is testable is that it refuses
  ## nothing the constructor accepts.
  expect_error(aci_dyad_model(params = list(d_x = 0.5, gamma = 2, f_x = 0.5,
                                            s_x = 0, d_y = 0.5, f_y = 1,
                                            s_y = 1)),
               class = "aci_error_model_contract")
  small <- aci_dyad_model(params = list(d_x = 0.5, gamma = 2, f_x = 0.5,
                                        s_x = 1e-9, d_y = 0.5, f_y = 1,
                                        s_y = 1))
  ob <- as_obs(simulate(small, seed = 1, t_end = 0.2, dt = 0.01))
  bundle <- .compile_cgns_run(small, ob)
  expect_identical(bundle$realization, "dyad_directed")
  expect_equal(as.numeric(bundle$coefficients$gxx[1, 1, 1]), 1e-18,
               tolerance = 0)
  ## The ordinary dyad model is unaffected.
  m <- aci_dyad_model()
  expect_s3_class(.compile_cgns_run(m, as_obs(simulate(m, seed = 1,
                                                       t_end = 0.2,
                                                       dt = 0.01))),
                  "compiled_cgns")
})


test_that(paste0("an ill-conditioned but factorisable Gram is refused on ",
                 "the path as it is at the constructor"), {
  ## diag(1, 1e-14) has a Cholesky factor but rcond 1e-14, so a check that
  ## only reacted to a failed factorisation would accept it. The constructor
  ## already refuses it; the realised path must agree.
  expect_error(
    aci_model(Lx = function(t, x) diag(2), fx = function(t, x) c(0, 0),
              Ly = function(t, x) diag(-0.5, 2), fy = function(t, x) c(0, 0),
              Sx1 = function(t, x) diag(c(1, 1e-7)),
              Sy2 = function(t, x) diag(2), k = 2, l = 2),
    class = "aci_error_gram")
  m <- aci_model(
    Lx = function(t, x) diag(2), fx = function(t, x) c(0, 0),
    Ly = function(t, x) diag(-0.5, 2), fy = function(t, x) c(0, 0),
    Sx1 = function(t, x) if (abs(t - 0.01) < 1e-12)
      diag(c(1, 1e-7)) else diag(2),
    Sy2 = function(t, x) diag(2), k = 2, l = 2)
  ob <- observed_trajectory(
    c(0, 0.01, 0.02, 0.03),
    cbind(c(0, 0.1, 0.2, 0.3), c(0, -0.05, -0.1, -0.15)))
  e <- expect_error(
    quietly(aci_filter(m, ob, init = list(mean = c(0, 0), cov = diag(2)))),
    class = "aci_error_gram_path")
  expect_identical(e$index, 2L)
  expect_equal(e$rcond, 1e-14, tolerance = 1e-9)
  gxx <- array(0, c(2, 2, 2))
  gxx[, , 1] <- diag(2)
  gxx[, , 2] <- diag(c(1, 1e-14))
  expect_true(!is.null(tryCatch(chol(gxx[, , 2]), error = function(e) NULL)))
})


test_that("valid records are untouched by the check", {
  ## A positive-noise control of the weak-coupling case.
  pos <- aci_model(
    Lx = function(t, x) matrix(1e-6, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-0.5, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(0.5, 1, 1), k = 1, l = 1)
  f <- quietly(aci_filter(pos, gram_scalar_record(), init = gram_scalar_init))
  expect_equal(as.numeric(f$mean),
               c(0, 9.9999999999999995e-08, 1.9874999999999799e-07,
                 2.9626374999999405e-07),
               tolerance = 1e-12)
  expect_equal(as.numeric(f$cov),
               c(1, 0.99249999999998995, 0.98507499999998016,
                 0.97772424999997065),
               tolerance = 1e-12)
  expect_equal(f$meta$loglik, 2.6509396793681486, tolerance = 1e-12)
  expect_identical(f$meta$regularization$n_events, 0L)

  ## A dense correlated model with rectangular shared channels.
  dense <- aci_model(
    Lx = function(t, x) matrix(c(0.2, 0.1, -0.05, 0.3), 2, 2),
    fx = function(t, x) c(0.1, -0.2),
    Ly = function(t, x) matrix(c(-0.4, 0.05, 0.02, -0.6), 2, 2),
    fy = function(t, x) c(0, 0.1),
    Sx1 = function(t, x) matrix(c(0.5, 0.1, 0.2, 0.4, -0.1, 0.3), 2, 3),
    Sy1 = function(t, x) matrix(c(0.2, 0.05, -0.1, 0.3, 0.15, 0.1), 2, 3),
    Sx2 = function(t, x) matrix(c(0.3, 0.1), 2, 1),
    Sy2 = function(t, x) matrix(c(0.25, 0.35), 2, 1), k = 2, l = 2)
  obd <- observed_trajectory(c(0, 0.01, 0.02),
                             cbind(c(0, 0.1, 0.2), c(0, -0.05, -0.1)))
  fd <- quietly(aci_filter(dense, obd,
                           init = list(mean = c(0, 0), cov = diag(2))))
  expect_true(all(is.finite(fd$mean)))
  expect_equal(fd$meta$loglik, 2.5830017959972604, tolerance = 1e-12)
  expect_identical(fd$meta$regularization$n_events, 0L)

  ## Shipped library models stay orders of magnitude above the threshold.
  min_rcond <- function(model, obs) {
    full <- .realise_cgns_grid_once(model, obs)
    kk <- dim(full$gxx)[1L]
    min(vapply(seq_len(dim(full$gxx)[3L]), function(j) {
      g <- full$gxx[, , j]; dim(g) <- c(kk, kk); rcond(g)
    }, numeric(1)))
  }
  dyad <- aci_dyad_model()
  expect_gt(min_rcond(dyad, as_obs(simulate(dyad, seed = 1, t_end = 2,
                                            dt = 0.01))), 1e-4)
  for (hidden in list(c("hW", "tau"), c("u", "hW", "tau"))) {
    e <- aci_enso_model(hidden = hidden)
    ob <- as_obs(simulate(e, seed = 12, t_end = 4, dt = 0.005, burn_in = 0))
    expect_gt(min_rcond(e, ob), 1e-4)
  }
})
