## acir reserve file
## Origin: aci/tests/testthat/test-03-engine.R:140-159
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: T5 one-lag column check; requires mode='one_lag'.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: T5: the one-lag column of lag_table(mode='one_lag') against brute-force truncated smoothing.

test_that("T5: one-lag column matches brute-force truncated smoothing", {
  m <- model_dyad(); dt <- 5e-3
  s <- simulate(m, seed = 9, T = 2, dt = dt)
  filt <- suppressWarnings(da_filter(m, s$obs))
  smoo <- online_smoother(m, s$obs, filt)
  lt <- lag_table(m, s$obs, mode = "one_lag", filter = filt, smoother = smoo)
  N1 <- length(s$obs$t)
  obs_tr <- observed_trajectory(s$obs$t[-N1], s$obs$x[-N1, , drop = FALSE])
  filt_tr <- suppressWarnings(da_filter(m, obs_tr))
  smoo_tr <- online_smoother(m, obs_tr, filt_tr)
  brute <- sapply(1:(N1 - 1), function(j)
    unname(gaussian_kl(smoo$mean[j, ], smoo$cov[, , j],
                       smoo_tr$mean[j, ], smoo_tr$cov[, , j], decompose = FALSE)))
  ol <- lt_onelag(lt)[1:(N1 - 1)]
  keep <- brute > 1e-7          # compare where signal exists
  expect_gt(sum(keep), 20)
  relerr <- abs(ol[keep] - brute[keep]) / brute[keep]
  expect_lt(stats::median(relerr), 0.05)
  expect_gt(stats::cor(ol[keep], brute[keep]), 0.999)
})


## Reserve staging block
## Origin: aci 0.0.30 (git tree 97f6b124) tests/testthat/test-03-engine.R:235-248
## Category: fbcir
## Note: One-lag propagation / operator-norm truncation block.

test_that("one-lag propagation does not truncate on an unscaled operator norm", {
  dt <- 0.01; ax <- (1 - 5e-14) / dt
  m <- cgns_model(
    Lx = function(t, x) matrix(ax, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(0, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
  ob <- observed_trajectory(dt * 0:3, matrix(c(0, 0, 0, 1e6), ncol = 1))
  tab <- lag_table(m, ob, mode = "one_lag",
                   init = list(mean = 0, cov = matrix(1 / ax, 1, 1)))
  expect_true(all(is.finite(lt_onelag(tab))))
  expect_gt(lt_onelag(tab)[2], 1)
  expect_true(is.na(tab$meta$stop_index))
})
