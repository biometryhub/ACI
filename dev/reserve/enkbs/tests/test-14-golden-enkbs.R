## acir reserve file
## Origin: aci/tests/testthat/test-14-golden-enkbs.R
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: G3/G4 golden grade against the published EnKBS dyad experiment.
## Verbatim copy from the aci 0.0.30 sources; not modified.

# Golden alignment of the ensemble engine against a literal port of the
# published EnKBS dyad experiment (helper-golden-p3.R; EnKBS-main/dyad/utov.m).
# Both sides are driven by one set of increments and one initial ensemble, so
# every difference below is implementation, not randomness. The port adds its
# observation perturbation as the reference does; the package receives the
# negated draws and subtracts, which realises identical innovations.

.golden_p3_params <- function() {
  list(d_u = 0.5, d_v = 0.5, c = 2, F_u = 1, F_v = 0.8,
       sigma_u = 0.5, sigma_v = 1)
}

.golden_p3_setup <- function(n1 = 161L, m = 40L, dt = 1e-3, seed = 42L) {
  p <- .golden_p3_params()
  set.seed(seed)
  truth <- golden_p3_truth(n1, dt, p, stats::rnorm(n1 - 1), stats::rnorm(n1 - 1))
  z_u <- matrix(stats::rnorm(m * (n1 - 1)), m)
  z_v <- matrix(stats::rnorm(m * (n1 - 1)), m)
  u0 <- truth$u[1] + 0.1 * stats::rnorm(m)
  B <- array(NA_real_, c(1L, n1 - 1L, m))
  W <- array(NA_real_, c(1L, n1 - 1L, m))
  for (e in seq_len(m)) {
    B[1L, , e] <- z_u[e, ]
    W[1L, , e] <- -z_v[e, ]
  }
  list(p = p, dt = dt, n1 = n1, m = m, truth = truth, z_u = z_u, u0 = u0,
       model = model_dyad(variant = "p3", observe = "y"),
       obs = observed_trajectory(dt * (seq_len(n1) - 1L),
                                 matrix(truth$v, ncol = 1L)),
       noise = structure(list(B = B, W = W, seed = NULL),
                         class = "noise_store"),
       filt = golden_p3_enkbf(u0, truth$v, z_u, z_v, p, dt))
}

test_that("G3: EnKBF and EnKBS match the ported published dyad code", {
  gs <- .golden_p3_setup()
  run <- enkbf(gs$model, gs$obs, m = gs$m,
               ic_sampler = function(mm) matrix(gs$u0, 1L, mm),
               noise = gs$noise)
  got_f <- matrix(run$path$members[1L, , ], gs$n1, gs$m)
  expect_lt(max(abs(got_f - t(gs$filt))), 1e-10)

  port_s <- golden_p3_enkbs(gs$filt, gs$truth$v, gs$z_u, gs$p, gs$dt)
  sm <- enkbs(run$model, run$path, run$noise)
  got_s <- matrix(sm$members[1L, , ], gs$n1, gs$m)
  expect_lt(max(abs(got_s - t(port_s))), 1e-9)
})

test_that("G4: ensemble lag table and forward CIR match the port triangle", {
  gs <- .golden_p3_setup()
  run <- enkbf(gs$model, gs$obs, m = gs$m,
               ic_sampler = function(mm) matrix(gs$u0, 1L, mm),
               noise = gs$noise)
  sm <- enkbs(run$model, run$path, run$noise)
  tab <- suppressWarnings(
    .ensemble_lag_table_from_run(gs$model, gs$obs, run, sm)
  )
  gd <- golden_p3_delta(gs$filt, gs$truth$v, gs$z_u, gs$p, gs$dt)

  # The diagonal is the ACI metric of the ensemble route: the complete-record
  # posterior against the filter, which the reference stores at (j, j).
  expect_lt(max(abs(tab$diag - diag(gd$delta))), 1e-8)

  # Every retained lagged cell: row j holds horizons j..N1.
  for (j in seq_len(gs$n1)) {
    expect_lt(max(abs(tab$rows[[j]] - gd$delta[j, j:gs$n1])), 1e-8)
  }

  # The forward ensemble CIR: the reference's Simpson ratio per anchor. The
  # port zeroes anchors whose peak is at or below its threshold; the package
  # masks them, so the comparison is on the anchors both sides resolve.
  fc <- forward_cir(tab, method = "l1_linf", min_M = NULL)
  port_cir <- golden_p3_cir(gd$delta, gs$dt)
  peaks <- vapply(seq_len(gs$n1), function(r)
    max(gd$delta[r, r:gs$n1]), numeric(1))
  keep <- which(peaks > 1e-5)
  expect_gt(length(keep), 50L)
  expect_lt(max(abs(fc$M[keep] - peaks[keep])), 1e-10)
  expect_lt(max(abs(fc$tau[keep] - port_cir[keep])), 1e-9)
})
