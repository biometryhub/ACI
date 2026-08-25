compiled_conditioning_fixture <- function(counter = NULL) {
  counted <- function(name, fun) {
    force(name)
    force(fun)
    function(t, x) {
      if (!is.null(counter)) counter[[name]] <- counter[[name]] + 1L
      fun(t, x)
    }
  }
  model <- cgns_model(
    Lx = counted("Lx", function(t, x)
      matrix(c(0.7 + 0.1 * x[1], -0.25 + 0.05 * x[2]), 2, 1)),
    fx = counted("fx", function(t, x)
      c(-0.3 * x[1] + 0.1 * sin(t), -0.2 * x[2] - 0.05)),
    Ly = counted("Ly", function(t, x) matrix(-0.8 + 0.02 * x[1], 1, 1)),
    fy = counted("fy", function(t, x) 0.15 * cos(t) - 0.1 * x[1]),
    Sx1 = counted("Sx1", function(t, x) diag(c(0.7, 0.6), 2)),
    Sx2 = counted("Sx2", function(t, x) matrix(0, 2, 1)),
    Sy1 = counted("Sy1", function(t, x) matrix(c(0.04, 0), 1, 2)),
    Sy2 = counted("Sy2", function(t, x) matrix(0.8, 1, 1)),
    k = 2, l = 1, name = "compiled-conditioning-test",
    meta = list(ic_default = list(x0 = c(0.3, -0.2), y0 = 0.1))
  )
  tt <- seq(0, 0.2, by = 0.005)
  xx <- cbind(
    target = 0.3 + 0.1 * sin(3 * tt),
    prescribed = -0.2 + 0.08 * cos(2 * tt)
  )
  list(
    model = model,
    obs = observed_trajectory(tt, xx),
    init = list(mean = 0.2, cov = matrix(0.4, 1, 1))
  )
}


expect_compiled_conditioning_fields <- function(strategy) {
  ds <- compiled_conditioning_fixture()
  nt <- nontarget(2, strategy)
  compiled <- .compile_cgns_complete(ds$model, ds$obs, nt)
  ix <- .nt_indices(nt, ds$obs)
  prescribed <- identical(strategy, "prescribed_forcing")
  target <- ix$A
  expected_k <- if (prescribed) length(target) else ds$model$k

  expect_identical(compiled$source_model, ds$model)
  expect_identical(compiled$source_obs, ds$obs)
  expect_identical(compiled$provenance$source_model, ds$model)
  expect_identical(compiled$provenance$source_obs, ds$obs)
  expect_identical(compiled$nontarget, nt)
  expect_identical(compiled$likelihood_idx,
                   if (prescribed) 1L else as.integer(target))
  expect_identical(compiled$k, as.integer(expected_k))
  expect_equal(
    compiled$x,
    if (prescribed) ds$obs$x[, target, drop = FALSE] else ds$obs$x,
    tolerance = 0
  )

  path_cross <- FALSE
  for (j in seq_along(ds$obs$t)) {
    co <- eval_coefs(ds$model, ds$obs$t[j], ds$obs$x[j, ])
    expected <- if (prescribed) list(
      Lx = co$Lx[target, , drop = FALSE], fx = co$fx[target],
      Ly = co$Ly, fy = co$fy,
      gxx = co$gxx[target, target, drop = FALSE],
      gyy = co$gyy, gyx = co$gyx[, target, drop = FALSE]
    ) else co[c("Lx", "fx", "Ly", "fy", "gxx", "gyy", "gyx")]
    expect_equal(matrix(compiled$coefficients$Lx[, , j], expected_k, 1),
                 expected$Lx, tolerance = 0)
    expect_equal(compiled$coefficients$fx[j, ], expected$fx, tolerance = 0)
    expect_equal(matrix(compiled$coefficients$Ly[, , j], 1, 1),
                 expected$Ly, tolerance = 0)
    expect_equal(compiled$coefficients$fy[j, ], expected$fy, tolerance = 0)
    expect_equal(matrix(compiled$coefficients$gxx[, , j],
                        expected_k, expected_k), expected$gxx, tolerance = 0)
    expect_equal(matrix(compiled$coefficients$gyy[, , j], 1, 1),
                 expected$gyy, tolerance = 0)
    expect_equal(matrix(compiled$coefficients$gyx[, , j], 1, expected_k),
                 expected$gyx, tolerance = 0)
    path_cross <- path_cross || .has_cross_noise(expected)
    if (j < length(ds$obs$t)) {
      expected_weight <- if (prescribed) {
        chol_solve(expected$gxx, diag(expected_k), "gxx")
      } else {
        masked_ginv(co$gxx, target)
      }
      expect_equal(
        matrix(compiled$coefficients$gxx_weight[, , j],
               expected_k, expected_k),
        expected_weight, tolerance = 0
      )
    }
  }
  expect_identical(compiled$correlated_noise,
                   isTRUE(compiled$model$meta$correlated_noise) || path_cross)

  compiled_filter <- .cgns_filter_compiled(
    compiled, ds$init, stepper = "explicit", nsub = 2L
  )
  public_filter <- da_filter(
    ds$model, ds$obs, init = ds$init, nontarget = nt,
    stepper = "explicit", nsub = 2L
  )
  expect_equal(compiled_filter$mean, public_filter$mean, tolerance = 0)
  expect_equal(compiled_filter$cov, public_filter$cov, tolerance = 0)
  expect_equal(compiled_filter$meta$loglik, public_filter$meta$loglik,
               tolerance = 1e-13)
  expect_identical(compiled_filter$meta$likelihood_idx,
                   compiled$likelihood_idx)
  expect_identical(compiled_filter$meta$source_model, ds$model)

  compiled_smoother <- .cgns_smoother_compiled(
    compiled, compiled_filter
  )
  public_smoother <- da_smooth(
    ds$model, ds$obs, filter = public_filter, nontarget = nt
  )
  expect_equal(compiled_smoother$mean, public_smoother$mean, tolerance = 0)
  expect_equal(compiled_smoother$cov, public_smoother$cov, tolerance = 0)
  expect_identical(compiled_smoother$meta$route, public_smoother$meta$route)
  expect_identical(compiled_smoother$meta$source_model, ds$model)
}


test_that("one-pass inflate conditioning matches direct coefficient equations", {
  expect_compiled_conditioning_fields("inflate")
})


test_that("one-pass prescribed forcing matches direct coefficient equations", {
  expect_compiled_conditioning_fields("prescribed_forcing")
})


.direct_compiled_loglik <- function(bundle, filter) {
  ia <- bundle$likelihood_idx
  ll <- 0
  for (j in seq_len(bundle$N)) {
    Lx <- matrix(
      bundle$coefficients$Lx[ia, , j, drop = FALSE],
      length(ia), bundle$l
    )
    innovation <- (
      bundle$rate[j, ia] -
        bundle$coefficients$fx[j, ia] -
        drop(Lx %*% filter$mean[j, ])
    ) * bundle$dt
    R <- matrix(filter$cov[, , j], bundle$l, bundle$l)
    gxx <- matrix(
      bundle$coefficients$gxx[ia, ia, j, drop = FALSE],
      length(ia), length(ia)
    )
    predictive <- Lx %*% R %*% t(Lx) * bundle$dt + gxx
    predictive <- (predictive + t(predictive)) * bundle$dt / 2
    logdet <- as.numeric(determinant(
      predictive, logarithm = TRUE
    )$modulus)
    ll <- ll - 0.5 * drop(crossprod(
      innovation, solve(predictive, innovation)
    )) - 0.5 * logdet - 0.5 * length(ia) * log(2 * pi)
  }
  ll
}


test_that("predictive likelihood matches its direct Gaussian equation", {
  ds <- compiled_conditioning_fixture()
  specs <- list(
    unconditioned = NULL,
    inflate = nontarget(2, "inflate"),
    prescribed = nontarget(2, "prescribed_forcing")
  )
  for (label in names(specs)) {
    bundle <- .compile_cgns_run(ds$model, ds$obs, specs[[label]])
    filter <- .cgns_filter_compiled(
      bundle, ds$init, stepper = "explicit", nsub = 1L
    )
    expect_equal(
      filter$meta$loglik,
      .direct_compiled_loglik(bundle, filter),
      tolerance = 1e-12,
      info = label
    )
  }
})


test_that("complete generic compilation covers the unconditioned contract", {
  ds <- compiled_conditioning_fixture()
  compiled <- .compile_cgns_complete(ds$model, ds$obs)

  expect_identical(compiled$likelihood_idx, 1:2)
  expect_null(compiled$nontarget)
  expect_identical(compiled$source_model, ds$model)
  expect_identical(compiled$source_obs, ds$obs)
  for (j in seq_along(ds$obs$t)) {
    co <- eval_coefs(ds$model, ds$obs$t[j], ds$obs$x[j, ])
    expect_equal(matrix(compiled$coefficients$Lx[, , j], 2, 1),
                 co$Lx, tolerance = 0)
    expect_equal(compiled$coefficients$fx[j, ], co$fx, tolerance = 0)
    expect_equal(matrix(compiled$coefficients$gxx[, , j], 2, 2),
                 co$gxx, tolerance = 0)
  }
})


test_that("the source coefficient closures are called once per grid point", {
  fields <- c("Lx", "fx", "Ly", "fy", "Sx1", "Sx2", "Sy1", "Sy2")
  for (strategy in list(NULL, "inflate", "prescribed_forcing")) {
    counter <- new.env(parent = emptyenv())
    for (field in fields) counter[[field]] <- 0L
    ds <- compiled_conditioning_fixture(counter)
    for (field in fields) counter[[field]] <- 0L
    nt <- if (is.null(strategy)) NULL else nontarget(2, strategy)

    bundle <- .compile_cgns_complete(ds$model, ds$obs, nt)

    expect_s3_class(bundle, "compiled_cgns")
    expect_identical(
      unname(vapply(fields, function(field) counter[[field]], integer(1))),
      rep.int(as.integer(length(ds$obs$t)), length(fields)),
      info = strategy %||% "unconditioned"
    )
  }
})


test_that("prescribed forcing keeps the established A-B cross-block tolerance", {
  cross_model <- function(epsilon) cgns_model(
    Lx = function(t, x) matrix(c(1, 0.5), 2, 1),
    fx = function(t, x) c(0, 0),
    Ly = function(t, x) matrix(-1, 1, 1),
    fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(c(1, epsilon, 0, 1), 2, 2,
                                 byrow = TRUE),
    Sy2 = function(t, x) matrix(0.5, 1, 1),
    k = 2, l = 1
  )
  tt <- 0.01 * 0:10
  obs <- observed_trajectory(tt, cbind(sin(tt), cos(tt)))
  nt <- nontarget(2, "prescribed_forcing")

  expect_no_error(.compile_cgns_complete(cross_model(0.5e-12), obs, nt))
  expect_error(
    .compile_cgns_complete(cross_model(2e-12), obs, nt),
    class = "aci_error_nontarget_crossnoise"
  )
})


test_that("authenticated affine batches preserve exact coefficient fields", {
  calls <- new.env(parent = emptyenv())
  calls$f <- 0L
  calls$g <- 0L
  f_full <- function(t, x, y) {
    calls$f <- calls$f + 1L
    c(0.2 - 0.3 * x[1] + (0.7 + 0.1 * x[1]) * y[1] - 0.2 * y[2],
      -0.1 - 0.2 * x[2] + 0.05 * y[1] + 0.4 * y[2])
  }
  g_full <- function(t, x, y) {
    calls$g <- calls$g + 1L
    c(0.1 * sin(t) - 0.6 * y[1] + 0.08 * y[2],
      -0.03 * x[1] - 0.04 * y[1] - 0.5 * y[2])
  }
  model <- cgns_from_affine(
    f_full, g_full,
    Sx = function(t, x) matrix(c(0.7, 0.05, 0, 0.6), 2, 2),
    Sy_hidden = function(t, x) matrix(c(0.8, 0.02, 0, 0.7), 2, 2),
    Sy_shared = function(t, x) matrix(c(0.03, 0, 0.01, 0.02), 2, 2),
    k = 2L, l = 2L
  )
  tt <- seq(0, 0.1, by = 0.01)
  obs <- observed_trajectory(
    tt, cbind(0.2 + 0.05 * sin(tt), -0.1 + 0.04 * cos(tt))
  )

  calls$f <- calls$g <- 0L
  generic <- .compile_cgns_complete(model, obs)
  generic_calls <- c(f = calls$f, g = calls$g)
  calls$f <- calls$g <- 0L
  batch <- .compile_cgns_run(model, obs)
  batch_calls <- c(f = calls$f, g = calls$g)

  expect_identical(batch$realization, "affine_batch")
  for (field in names(generic$coefficients))
    expect_identical(batch$coefficients[[field]],
                     generic$coefficients[[field]], info = field)
  expect_identical(batch$correlated_noise, generic$correlated_noise)
  expect_identical(batch_calls,
                   c(f = length(tt) * 3L, g = length(tt) * 3L))
  expect_identical(generic_calls,
                   c(f = length(tt) * 4L, g = length(tt) * 4L))

  changed <- model
  changed$fx <- function(t, x) c(0, 0)
  expect_null(.cgns_realizer_descriptor(changed))
  expect_identical(
    .compile_cgns_run(changed, obs)$realization,
    "generic_closure_one_pass"
  )
})


test_that("affine batches condition realised arrays without changing estimands", {
  model <- cgns_from_affine(
    f_full = function(t, x, y)
      c(-0.2 * x[1] + 0.7 * y, -0.1 * x[2] - 0.25 * y),
    g_full = function(t, x, y) 0.05 * x[1] - 0.6 * y,
    Sx = function(t, x) diag(c(0.7, 0.6)),
    Sy_hidden = function(t, x) matrix(0.8, 1, 1),
    k = 2L, l = 1L
  )
  tt <- seq(0, 0.1, by = 0.01)
  obs <- observed_trajectory(tt, cbind(sin(tt), cos(tt)))
  for (strategy in c("inflate", "prescribed_forcing")) {
    nt <- nontarget(2L, strategy)
    generic <- .compile_cgns_complete(model, obs, nt)
    batch <- .compile_cgns_run(model, obs, nt)
    expect_identical(batch$realization, "affine_batch")
    expect_identical(batch$likelihood_idx, generic$likelihood_idx)
    for (field in names(generic$coefficients))
      expect_identical(batch$coefficients[[field]],
                       generic$coefficients[[field]], info = strategy)
  }
})
