# Theorem-3 auxiliaries for a scalar hidden state as vector arithmetic
# (.online_aux_scalar, .forward_primitives_scalar,
# .smoother_thmD1_scalar_moments) against the per-interval kernels they
# stand in for, on a k = 1 and a k = 5 record, and the covariance policy
# reached by the same route from both.

.scalar_aux_case <- function(which) {
  if (which == "dyad") {
    model <- aci_dyad_model()
    sim <- simulate(model, seed = 271, T = 3, dt = 0.005, burn_in = 0)
    init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  } else {
    model <- aci_enso_model(hidden = "u")
    sim <- stats::simulate(model, seed = 3, T = 4, dt = 0.005)
    init <- list(mean = 0, cov = matrix(0.1, 1, 1))
  }
  obs <- as_obs(sim)
  bundle <- .compile_cgns_run(model, obs)
  filt <- .cgns_filter_compiled(bundle, init = init, stepper = "explicit",
                                nsub = 1L, validate = FALSE)
  list(bundle = bundle, filt = filt)
}

## The per-interval kernels, one interval at a time.
.scalar_aux_by_interval <- function(bundle, filt, rec) {
  N <- bundle$N
  E <- b <- P <- numeric(N)
  F_ <- matrix(NA_real_, N, bundle$k)
  mu1 <- R1 <- numeric(N)
  for (n in seq_len(N)) {
    rec$j <- n
    co <- .compiled_co(bundle, n)
    aux <- .thmD1_aux_compiled(bundle, n, filt$cov[, , n], co = co, rec = rec)
    ol <- .onelag_stats(
      co, aux, filt$mean[n, ], filt$cov[, , n], filt$mean[n + 1L, ],
      filt$cov[, , n + 1L], bundle$x[n + 1L, ] - bundle$x[n, ], bundle$dt,
      1L, rec
    )
    E[n] <- aux$E[1L, 1L]
    F_[n, ] <- aux$F
    mu1[n] <- ol$mu
    R1[n] <- ol$R
  }
  list(E = E, F = F_, mu1 = mu1, R1 = R1)
}

.scalar_aux_gate <- 1e-12

for (which in c("dyad", "enso")) test_that(sprintf(
  "the scalar auxiliaries match the per-interval kernels (%s)", which), {
  cs <- .scalar_aux_case(which)
  bundle <- cs$bundle
  filt <- cs$filt
  expect_identical(bundle$l, 1L)
  ref <- .scalar_aux_by_interval(bundle, filt, .aci_reg_for(NULL, bundle$t))
  got <- .online_aux_scalar(bundle, filt)
  expect_identical(got$idx, seq_len(bundle$N))
  ## within a few roundings of the kernels' own arithmetic: the last bit of
  ## the inverse variance is the BLAS's triangular solve, and for k > 1 the
  ## gain row sums k products in the BLAS's order
  expect_lt(max(abs(got$E - ref$E) / abs(ref$E)), 1e-15)
  expect_lt(max(abs(got$F - ref$F) / pmax(1, abs(ref$F))), 1e-14)
  prim <- .forward_primitives_scalar(bundle, filt)
  expect_lt(max(abs(prim$one_mu[, 1L] - ref$mu1) / pmax(1, abs(ref$mu1))),
            1e-15)
  expect_lt(max(abs(prim$one_R_v - ref$R1) / abs(ref$R1)), 1e-15)
  ## the cells the matrix path reads
  expect_identical(prim$E_v, vapply(prim$E, function(e) e[1L, 1L], 1))
  expect_identical(prim$dR_v, vapply(prim$dR, function(e) e[1L, 1L], 1))
  expect_equal(prim$dmu[, 1L], ref$mu1 - filt$mean[seq_len(bundle$N), 1L])
  ## a window skips the intervals before it and agrees on the rest
  late <- .forward_primitives_scalar(bundle, filt, from = 40L)
  expect_null(late$E[[39L]])
  expect_identical(late$E[[40L]], prim$E[[40L]])
  expect_identical(late$s_n[41L:bundle$N], prim$s_n[41L:bundle$N])
  expect_identical(late$s_n[1L:39L], numeric(39L))
})


for (which in c("dyad", "enso")) test_that(sprintf(
  "the scalar smoother moments are the Theorem 3 recursion's (%s)", which), {
  cs <- .scalar_aux_case(which)
  bundle <- cs$bundle
  filt <- cs$filt
  ## the recursion, one interval at a time, on the kernels
  N1 <- bundle$N1
  rec <- .aci_reg_for(NULL, bundle$t)
  MU <- numeric(N1)
  CV <- numeric(N1)
  MU[N1] <- filt$mean[N1, 1L]
  CV[N1] <- filt$cov[1L, 1L, N1]
  for (j in (N1 - 1L):1L) {
    rec$j <- j
    co <- .compiled_co(bundle, j)
    aux <- .thmD1_aux_compiled(bundle, j, filt$cov[, , j], co = co, rec = rec)
    st <- .onelag_stats(
      co, aux, filt$mean[j, ], filt$cov[, , j], MU[j + 1L],
      matrix(CV[j + 1L], 1L, 1L), bundle$x[j + 1L, ] - bundle$x[j, ],
      bundle$dt, 1L, rec
    )
    MU[j] <- st$mu
    CV[j] <- st$R
  }
  got <- .smoother_thmD1_scalar_moments(bundle, filt)
  expect_lt(max(abs(got$mu - MU) / pmax(1, abs(MU))), 1e-15)
  expect_lt(max(abs(got$cv - CV) / abs(CV)), 1e-15)
  ## and the smoother verb returns them
  sm <- .smoother_thmD1_compiled(bundle, filt, validate = FALSE,
                                 warn_cost = FALSE)
  expect_identical(sm$mean[, 1L], got$mu)
  expect_identical(sm$cov[1L, 1L, ], got$cv)
  expect_identical(sm$meta$route, "thmD1")
})


test_that("a variance the policy must see sends both routes the same way", {
  cs <- .scalar_aux_case("dyad")
  bundle <- cs$bundle
  bad <- cs$filt
  bad$cov[1L, 1L, 10L] <- -1e-9
  expect_null(.online_aux_scalar(bundle, bad))
  expect_null(.forward_primitives_scalar(bundle, bad))
  expect_null(.smoother_thmD1_scalar_moments(bundle, bad))
  ## strict policy: the per-interval route aborts, and so does the verb
  expect_error(
    .smoother_thmD1_compiled(bundle, bad, validate = FALSE, warn_cost = FALSE,
                             regularize = "none"),
    class = "aci_error_covariance_not_spd"
  )
  ## floor policy: one recorded event at that interval, the floored value
  ## carried through the recursion
  rec <- .aci_reg_for("floor", bundle$t)
  sm <- .smoother_thmD1_compiled(bundle, bad, validate = FALSE,
                                 warn_cost = FALSE, regularize = rec)
  frozen <- .aci_reg_freeze(rec)
  expect_identical(frozen$n_events, 1L)
  expect_true(all(is.finite(sm$cov)))
  expect_true(all(sm$cov > 0))
})


test_that("the online window route reads the scalar auxiliaries", {
  cs <- .scalar_aux_case("dyad")
  bundle <- cs$bundle
  filt <- cs$filt
  rec <- .aci_reg_for(NULL, bundle$t)
  aux <- .online_aux_compiled(bundle, filt, rec)
  ref <- .scalar_aux_by_interval(bundle, filt, .aci_reg_for(NULL, bundle$t))
  expect_lt(max(abs(aux$E[1L, 1L, ] - ref$E) / abs(ref$E)), 1e-15)
  expect_lt(max(abs(aux$dmu[, 1L] - (ref$mu1 - filt$mean[-bundle$N1, 1L]))),
            1e-15)
  win <- .da_online_compiled(bundle, filt, lag = 25L, aux = aux,
                             regularize = rec)
  expect_identical(win$meta$route, "thmD1_online_window")
  expect_true(all(is.finite(win$mean)))
})
