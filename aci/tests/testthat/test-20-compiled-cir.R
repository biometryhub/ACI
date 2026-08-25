.expect_compiled_forward_equal <- function(actual, expected,
                                           tolerance = 1e-12) {
  expect_s3_class(actual, "cir_result")
  expect_identical(names(actual), names(expected))
  expect_equal(actual$t, expected$t, tolerance = 0)
  expect_equal(actual$tau, expected$tau, tolerance = tolerance)
  expect_equal(actual$M, expected$M, tolerance = tolerance)
  expect_equal(actual$interval, expected$interval, tolerance = tolerance)
  expect_identical(actual$direction, expected$direction)
  expect_identical(actual$method, expected$method)
  expect_identical(actual$bound, expected$bound)
  expect_equal(actual$tail_bound, expected$tail_bound,
               tolerance = tolerance)
  expect_equal(actual$subjective, expected$subjective,
               tolerance = tolerance)
  expect_equal(actual$dt, expected$dt, tolerance = 0)
  expect_identical(actual$meta, expected$meta)
}


.compiled_cir_matrix_setup <- function() {
  model <- cgns_model(
    Lx = function(t, x) matrix(c(
      0.65 + 0.04 * x[1], -0.12,
      0.08, 0.45 + 0.03 * x[2]
    ), 2, 2, byrow = TRUE),
    fx = function(t, x) c(-0.25 * x[1] + 0.08,
                           -0.18 * x[2] - 0.03),
    Ly = function(t, x) matrix(c(-0.7, -0.04, 0.09, -0.55),
                                2, 2, byrow = TRUE),
    fy = function(t, x) c(0.1 * sin(t) - 0.06 * x[1],
                           0.12 * cos(t) + 0.04 * x[2]),
    Sx1 = function(t, x) matrix(c(0.65, 0.05, 0, 0.55), 2, 2),
    Sy1 = function(t, x) matrix(c(0.05, 0.01, 0, 0.04), 2, 2),
    Sy2 = function(t, x) matrix(c(0.72, 0.03, 0, 0.64), 2, 2),
    k = 2, l = 2, name = "compiled-cir-matrix"
  )
  t <- seq(0, 0.12, by = 0.01)
  obs <- observed_trajectory(
    t, cbind(0.25 + 0.08 * sin(4 * t), -0.15 + 0.06 * cos(3 * t))
  )
  list(
    model = model, obs = obs,
    init = list(
      mean = c(0.15, -0.08),
      cov = matrix(c(0.42, 0.025, 0.025, 0.31), 2, 2)
    )
  )
}


test_that("streamed scalar forward CIR matches reduction of retained rows", {
  old <- options(aci.default_tol = 1e-5)
  on.exit(options(old), add = TRUE)
  model <- model_dyad()
  sim <- simulate(model, seed = 271, T = 0.5, dt = 0.005, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  bundle <- .compile_cgns_run(model, obs)
  filter <- .cgns_filter_compiled(bundle, init)
  table <- .lag_table_compiled(bundle, mode = "forward", filter = filter)

  for (method in c("exact", "l1_linf")) {
    for (quadrature in c("simpson", "sum")) {
      expected <- suppressWarnings(forward_cir(
        table, method = method, eps = c(1e-6, 1e-4), min_M = 0,
        quadrature = quadrature, simpson_close = "quadratic"
      ))
      actual <- suppressWarnings(.forward_cir_compiled(
        bundle, filter = filter, method = method,
        eps = c(1e-6, 1e-4), min_M = 0,
        quadrature = quadrature, simpson_close = "quadratic"
      ))
      .expect_compiled_forward_equal(actual, expected)
    }
  }
})


test_that("streamed correlated matrix forward CIR matches retained rows", {
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ds <- .compiled_cir_matrix_setup()
  bundle <- .compile_cgns_run(ds$model, ds$obs)
  filter <- .cgns_filter_compiled(bundle, ds$init)
  table <- .lag_table_compiled(bundle, mode = "forward", filter = filter)
  expected <- forward_cir(
    table, method = "exact", eps = 1e-8, min_M = 0,
    simpson_close = "trapezoid"
  )
  actual <- .forward_cir_compiled(
    bundle, filter = filter, method = "exact", eps = 1e-8, min_M = 0,
    simpson_close = "trapezoid"
  )
  .expect_compiled_forward_equal(actual, expected, tolerance = 1e-11)
})


test_that("forward_cir on an unstored ACI result uses the streamed route", {
  old <- options(aci.default_tol = 1e-6)
  on.exit(options(old), add = TRUE)
  model <- model_dyad()
  sim <- simulate(model, seed = 81, T = 0.25, dt = 0.005, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  result <- aci(model, obs, init = init, keep = "paths")
  expect_null(result$table)
  direct <- suppressWarnings(forward_cir(
    result, method = "l1_linf", min_M = 0
  ))
  table <- lag_table(
    model, obs, mode = "forward", filter = result$paths$filter,
    init = result$handles$init
  )
  expected <- suppressWarnings(forward_cir(
    table, method = "l1_linf", min_M = 0
  ))
  .expect_compiled_forward_equal(direct, expected)
})


test_that("streamed forward CIR preserves warnings and argument errors", {
  model <- model_dyad()
  sim <- simulate(model, seed = 9, T = 0.08, dt = 0.01, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  bundle <- .compile_cgns_run(model, obs)
  implicit <- .cgns_filter_compiled(
    bundle, init, stepper = "implicit", nsub = 1L
  )
  expect_warning(
    .forward_cir_compiled(bundle, filter = implicit, min_M = 0),
    class = "aci_warn_stepper"
  )
  expect_error(
    .forward_cir_compiled(bundle, eps = -1), class = "aci_error_dims"
  )
  expect_error(
    .forward_cir_compiled(bundle, nonsense = TRUE), class = "aci_error_dims"
  )
})


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
