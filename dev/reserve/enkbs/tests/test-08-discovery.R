## acir reserve file
## Origin: aci/tests/testthat/test-08-discovery.R
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Whole test file for the discovery family.
## Verbatim copy from the aci 0.0.30 sources; not modified.

test_that("D1: FFBS sample moments match the smoother marginals (linear)", {
  lin <- cgns_model(Lx = function(t,x) matrix(1,1,1), fx = function(t,x) 0.2,
                    Ly = function(t,x) matrix(-1,1,1), fy = function(t,x) 0.3,
                    Sx1 = function(t,x) matrix(0.5,1,1),
                    Sy2 = function(t,x) matrix(0.6,1,1), k = 1, l = 1)
  s <- simulate(lin, seed = 1, T = 2, dt = 5e-3, ic = list(x0 = 0, y0 = 0.5))
  filt <- da_filter(lin, s$obs, init = list(mean = 0, cov = 1))
  smoo <- da_smooth(lin, s$obs, filter = filt)
  Y <- sample_paths(lin, s$obs, n_samples = 600, seed = 2, filter = filt)
  idx <- c(50, 150, 250, 350)
  for (j in idx) {
    expect_lt(abs(mean(Y[1, j, ]) - smoo$mean[j, 1]),
              4 * sqrt(smoo$cov[1, 1, j] / 600) + 0.02)
    expect_lt(abs(var(Y[1, j, ]) / smoo$cov[1, 1, j] - 1), 0.25)
  }
  # path continuity: consecutive-time correlation should be high
  expect_gt(cor(Y[1, 200, ], Y[1, 201, ]), 0.9)
})

test_that("polynomial library honours degree and includes P3 cubic terms", {
  l0 <- cgns_library(k = 1, l = 1, degree = 0)
  expect_equal(l0$names, "1")
  l3 <- cgns_library(k = 1, l = 2, degree = 3)
  expect_true(all(c("x1*y1^2", "x1*y2^2", "x1*y1*y2") %in% l3$names))
})

test_that("candidate-only causation entropy is invariant to unrequested nuisance columns", {
  set.seed(17)
  a <- rnorm(100); b <- rnorm(100); nuisance <- a + b
  target <- a + 0.2 * rnorm(100)
  base <- causation_entropy(cbind(a = a, b = b), target, c("a", "b"))
  wide <- causation_entropy(cbind(a = a, b = b, nuisance = nuisance),
                            target, c("a", "b"))
  expect_equal(wide, base, tolerance = 0)
})

test_that("learned model preserves shared diffusion channels", {
  lib <- cgns_library(1, 1, degree = 1)
  template <- list(k = 1, l = 1,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sx2 = function(t, x) matrix(0, 1, 1),
    Sy1 = function(t, x) matrix(0.5, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1))
  rebuilt <- model_from_learned(list(c("1" = 0, "y1" = 1)),
                                list(c("y1" = -1)), lib, template)
  expect_equal(cgns_grams(rebuilt, 0, 0)$gyx[1, 1], 0.5)
})

test_that("D2/D3: dyad structure + coefficients recovered; nil term pruned", {
  m <- model_dyad()             # dx = (0.5 - 0.5 x + 2 x y) dt; dy = (1 - 0.5 y - 2 x^2) dt
  s <- simulate(m, seed = 6, T = 30, dt = 2e-3, burn_in = 2)
  lib <- cgns_library(k = 1, l = 1, degree = 2)
  lm_ <- suppressWarnings(learn_model(m, s$obs, lib, n_samples = 20,
                                      ce_threshold = "auto", seed = 3))
  ex <- lm_$equations$dx1; ey <- lm_$equations$dy1
  # structure: true terms recovered in both equations
  expect_true(all(c("1", "x1", "x1*y1") %in% ex$kept))
  expect_true(all(c("1", "x1^2") %in% ey$kept))
  # nil term stays below threshold (D3)
  expect_lt(ex$ce[["y1^2"]], attr(threshold_structure(ex$ce, "auto"), "threshold"))
  # observed equation coefficients: tight (regression on true observations)
  expect_lt(abs(ex$coef[["x1*y1"]] - 2) / 2, 0.25)
  # A single posterior-sampling M-step has the documented h-transform bias;
  # test signs/order rather than claiming oracle coefficient accuracy.
  expect_lt(ex$coef[["x1"]], 0)
  expect_gt(ex$coef[["1"]], 0)
  # hidden equation coefficients: single-pass h-transform bias documented in
  # learn_model(); assert sign + order of magnitude, not tight recovery
  expect_lt(ey$coef[["x1^2"]], -0.5)
  expect_gt(ey$coef[["1"]], 0.4)
  # Diffusion amplitudes use increment variance times dt (not derivative SD).
  expect_equal(ex$sigma_hat, 0.5, tolerance = 0.08)
  expect_equal(ey$sigma_hat, 1, tolerance = 0.15)
})

test_that("D4: constrained_mle satisfies KKT conditions exactly", {
  set.seed(4)
  Th <- matrix(rnorm(200 * 5), 200); colnames(Th) <- paste0("t", 1:5)
  xi_true <- c(1, -2, 0.5, 0, 3)
  y <- drop(Th %*% xi_true) + rnorm(200, 0, 0.1)
  A <- matrix(c(1, 1, 0, 0, 0), 1); b <- -1        # xi1 + xi2 = -1
  f <- constrained_mle(Th, y, constraints = list(A = A, b = b))
  expect_equal(drop(A %*% f$coef), -1, tolerance = 1e-10)
  grad <- drop(crossprod(Th) %*% f$coef - crossprod(Th, y))
  # stationarity: gradient lies in row space of A
  proj <- drop(t(A) %*% solve(A %*% t(A)) %*% (A %*% grad))
  expect_lt(max(abs(grad - proj)), 1e-6 * max(abs(grad), 1))
  f0 <- constrained_mle(Th, y)
  expect_lt(max(abs(f0$coef - xi_true)), 0.05)
})

test_that("D5: joint energy constraint holds exactly and preserves recovery", {
  m <- model_dyad(); s <- simulate(m, seed = 6, T = 30, dt = 2e-3, burn_in = 2)
  lib <- cgns_library(k = 1, l = 1, degree = 2)
  fit <- suppressWarnings(learn_model(m, s$obs, lib, n_samples = 15, seed = 3,
                                      energy_pairs = list(c(obs = 1, hid = 1))))
  cxy <- fit$equations$dx1$coef[["x1*y1"]]
  cx2 <- fit$equations$dy1$coef[["x1^2"]]
  expect_lt(abs(cxy + cx2), 1e-8)                 # exact conservation
  expect_lt(abs(cxy - 2) / 2, 0.3)                # recovery intact
  expect_true(fit$equations$dx1$energy_constrained)
})

test_that("discovery rejects unnamed equation constraints and unused FFBS arguments", {
  m <- model_dyad(); s <- simulate(m, seed = 1, T = 0.1, dt = 0.01)
  lib <- cgns_library(1, 1, degree = 1)
  expect_error(learn_model(m, s$obs, lib, n_samples = 2,
                           constraints = list(list(A = matrix(1), b = 0))),
               class = "aci_error_dims")
  expect_error(sample_paths(m, s$obs, n_samples = 2, typo = TRUE),
               class = "aci_error_dims")
})
