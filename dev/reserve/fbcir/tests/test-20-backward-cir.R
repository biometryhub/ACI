## acir reserve file
## Origin: aci/tests/testthat/test-20-compiled-cir.R:144-183
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Backward CIR prefix-sharing assertions excised from test-20-compiled-cir.R.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Backward CIR prefix sharing across reference times (.slice_compiled_cgns / .slice_compiled_filter).

test_that("backward CIR prefixes reuse one compiled coefficient path", {
  calls <- 0L
  model <- cgns_model(
    Lx = function(t, x) {
      calls <<- calls + 1L
      matrix(0.7 + 0.05 * x[1], 1, 1)
    },
    fx = function(t, x) -0.3 * x,
    Ly = function(t, x) matrix(-0.8, 1, 1),
    fy = function(t, x) 0.1 * sin(t),
    Sx1 = function(t, x) matrix(0.5, 1, 1),
    Sy2 = function(t, x) matrix(0.7, 1, 1),
    k = 1L, l = 1L,
    meta = list(ic_default = list(x0 = 0, y0 = 0.2))
  )
  tt <- seq(0, 0.3, by = 0.005)
  obs <- observed_trajectory(tt, matrix(0.1 * sin(2 * tt), ncol = 1L))
  init <- list(mean = 0.2, cov = matrix(0.4, 1, 1))
  result <- aci(model, obs, init = init, keep = "paths")
  calls <- 0L
  refs <- c(0.15, 0.225, 0.3)
  reused <- backward_cir(result, T = refs, method = "exact", min_M = 0)
  expect_identical(calls, as.integer(length(tt)))

  separate <- lapply(refs, function(Tv) {
    keep <- tt <= Tv + 1e-12
    table <- lag_table(
      model,
      observed_trajectory(tt[keep], obs$x[keep, , drop = FALSE]),
      mode = "one_lag", init = init
    )
    backward_cir(table, method = "exact", min_M = 0)
  })
  expect_equal(reused$t, vapply(separate, `[[`, numeric(1L), "t"),
               tolerance = 0)
  expect_equal(reused$tau, vapply(separate, `[[`, numeric(1L), "tau"),
               tolerance = 1e-12)
  expect_equal(reused$M, vapply(separate, `[[`, numeric(1L), "M"),
               tolerance = 1e-12)
})


## Reserve staging block
## Origin: aci 0.0.30 (git tree 97f6b124) tests/testthat/test-20-compiled-cir.R:186-216
## Category: fbcir
## Note: Backward CIR coefficient-realisation bound.

test_that("backward CIR does not realise coefficients beyond requested prefixes", {
  calls <- 0L
  model <- cgns_model(
    Lx = function(t, x) {
      calls <<- calls + 1L
      if (t > 0.151 && t < 0.4) return(matrix(NaN, 1L, 1L))
      matrix(0.6, 1L, 1L)
    },
    fx = function(t, x) -0.2 * x,
    Ly = function(t, x) matrix(-0.7, 1L, 1L),
    fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.5, 1L, 1L),
    Sy2 = function(t, x) matrix(0.8, 1L, 1L),
    k = 1L, l = 1L,
    meta = list(ic_default = list(x0 = 0, y0 = 0))
  )
  tt <- seq(0, 0.3, by = 0.005)
  obs <- observed_trajectory(tt, matrix(0.1 * sin(tt), ncol = 1L))
  init <- list(mean = 0, cov = matrix(0.3, 1L, 1L))
  calls <- 0L

  got <- expect_no_error(backward_cir(
    model, obs, T = c(0.1, 0.15), init = init, min_M = 0
  ))
  expect_equal(got$t, c(0.1, 0.15), tolerance = 0)
  expect_identical(calls, 31L)
  expect_error(
    backward_cir(model, obs, T = 0.2, init = init, min_M = 0),
    class = "aci_error_model_contract"
  )
})
