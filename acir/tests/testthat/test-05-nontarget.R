chain_model <- function(gA = 0) aci_model(
  Lx = function(t, x) matrix(c(gA, 2), 2, 1),
  fx = function(t, x) c(-x[1] + 1.5 * x[2], -x[2]),
  Ly = function(t, x) matrix(-1, 1, 1),
  fy = function(t, x) 0,
  Sx1 = function(t, x) diag(0.3, 2),
  Sy2 = function(t, x) matrix(0.5, 1, 1), k = 2, l = 1)

test_that("T7a: chain y -> xB -> xA gives near-zero conditional ACI", {
  m <- chain_model(gA = 0)
  s <- simulate(m, seed = 33, T = 4, dt = 2e-3, ic = list(x0 = c(0, 0), y0 = 0.5))
  a_un <- suppressWarnings(aci(m, s$obs))
  a_c  <- suppressWarnings(
    aci(m, s$obs, conditional = aci_conditional(given = 2, method = "reduce")))
  burn <- seq_len(200)
  expect_gt(max(a_un$aci[-burn]), 20 * max(a_c$aci[-burn]))
})
test_that("T7b: 'reduce' and 'mask' agree when both defined", {
  m <- chain_model(gA = 1)
  s <- simulate(m, seed = 34, T = 3, dt = 2e-3, ic = list(x0 = c(0, 0), y0 = 0.5))
  a1 <- suppressWarnings(aci(m, s$obs,
          conditional = aci_conditional(given = 2, method = "reduce")))
  a2 <- suppressWarnings(aci(m, s$obs,
          conditional = aci_conditional(given = 2, method = "mask")))
  expect_equal(a1$aci, a2$aci, tolerance = 1e-8)
  expect_gt(max(a1$aci[-(1:100)]), 1e-3)     # genuinely nonzero conditional link
})

test_that("prescribed-forcing paths retain stable source-model provenance", {
  m <- chain_model(gA = 1)
  s <- simulate(m, seed = 35, T = 0.2, dt = 0.01,
                ic = list(x0 = c(0, 0), y0 = 0.5))
  nt <- aci_conditional(given = 2, method = "reduce")
  f <- suppressWarnings(aci_filter(m, s$obs, conditional = nt))
  sm <- expect_no_error(aci_smoother(m, s$obs, filter = f, conditional = nt))
  expect_s3_class(sm, "da_path_gaussian")
  expect_no_error(lag_table(m, s$obs, mode = "forward", filter = f,
                            conditional = nt))
  expect_s3_class(suppressWarnings(aci(m, s$obs, conditional = nt)),
                  "aci_result")
})

test_that("non-target block specifications reject empty, fractional, and duplicate indices", {
  m <- chain_model(); s <- simulate(m, seed = 1, T = 0.1, dt = 0.01)
  for (bad in list(numeric(0), 1.5, c(1, 1), 3))
    expect_error(aci_filter(m, s$obs,
                            conditional = aci_conditional(bad, "mask")),
                 class = "aci_error_nontarget")
})

test_that("masked-innovation likelihood scores only target components", {
  m <- aci_model(
    Lx = function(t, x) matrix(c(1, 2), 2, 1),
    fx = function(t, x) c(0, 0),
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) diag(1, 2),
    Sy2 = function(t, x) matrix(0.5, 1, 1), k = 2, l = 1)
  tt <- 0.01 * 0:20
  x1 <- sin(tt)
  ob1 <- observed_trajectory(tt, cbind(x1, cos(tt)))
  ob2 <- observed_trajectory(tt, cbind(x1, 100 * cos(7 * tt)))
  nt <- aci_conditional(2, "mask")
  ini <- list(mean = 0, cov = matrix(1, 1, 1))
  f1 <- aci_filter(m, ob1, conditional = nt, init = ini)
  f2 <- aci_filter(m, ob2, conditional = nt, init = ini)
  expect_equal(f1$mean, f2$mean, tolerance = 0)
  expect_equal(f1$cov, f2$cov, tolerance = 0)
  expect_equal(f1$meta$loglik, f2$meta$loglik, tolerance = 0)
  expect_equal(f1$meta$likelihood_idx, 1L)
})
