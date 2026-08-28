compiled_conditioning_fixture <- function(counter = NULL) {
  counted <- function(name, fun) {
    force(name)
    force(fun)
    function(t, x) {
      if (!is.null(counter)) counter[[name]] <- counter[[name]] + 1L
      fun(t, x)
    }
  }
  model <- aci_model(
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


expect_compiled_conditioning_fields <- function(method) {
  ds <- compiled_conditioning_fixture()
  nt <- aci_conditional(2, method)
  compiled <- .compile_cgns_complete(ds$model, ds$obs, nt)
  ix <- .nt_indices(nt, ds$obs)
  prescribed <- identical(method, "reduce")
  target <- ix$A
  expected_k <- if (prescribed) length(target) else ds$model$k

  expect_identical(compiled$source_model, ds$model)
  expect_identical(compiled$source_obs, ds$obs)
  expect_identical(compiled$provenance$source_model, ds$model)
  expect_identical(compiled$provenance$source_obs, ds$obs)
  expect_identical(compiled$conditional, nt)
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
  public_filter <- aci_filter(
    ds$model, ds$obs, init = ds$init, conditional = nt,
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
  public_smoother <- aci_smoother(
    ds$model, ds$obs, filter = public_filter, conditional = nt
  )
  expect_equal(compiled_smoother$mean, public_smoother$mean, tolerance = 0)
  expect_equal(compiled_smoother$cov, public_smoother$cov, tolerance = 0)
  expect_identical(compiled_smoother$meta$route, public_smoother$meta$route)
  expect_identical(compiled_smoother$meta$source_model, ds$model)
}


test_that("one-pass mask conditioning matches direct coefficient equations", {
  expect_compiled_conditioning_fields("mask")
})


test_that("one-pass prescribed forcing matches direct coefficient equations", {
  expect_compiled_conditioning_fields("reduce")
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
    mask = aci_conditional(2, "mask"),
    reduce = aci_conditional(2, "reduce")
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
  expect_null(compiled$conditional)
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
  for (method in list(NULL, "mask", "reduce")) {
    counter <- new.env(parent = emptyenv())
    for (field in fields) counter[[field]] <- 0L
    ds <- compiled_conditioning_fixture(counter)
    for (field in fields) counter[[field]] <- 0L
    nt <- if (is.null(method)) NULL else aci_conditional(2, method)

    bundle <- .compile_cgns_complete(ds$model, ds$obs, nt)

    expect_s3_class(bundle, "compiled_cgns")
    expect_identical(
      unname(vapply(fields, function(field) counter[[field]], integer(1))),
      rep.int(as.integer(length(ds$obs$t)), length(fields)),
      info = method %||% "unconditioned"
    )
  }
})


test_that("prescribed forcing keeps the established A-B cross-block tolerance", {
  cross_model <- function(epsilon) aci_model(
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
  nt <- aci_conditional(2, "reduce")

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
  model <- aci_model_from_affine(
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
  model <- aci_model_from_affine(
    f_full = function(t, x, y)
      c(-0.2 * x[1] + 0.7 * y, -0.1 * x[2] - 0.25 * y),
    g_full = function(t, x, y) 0.05 * x[1] - 0.6 * y,
    Sx = function(t, x) diag(c(0.7, 0.6)),
    Sy_hidden = function(t, x) matrix(0.8, 1, 1),
    k = 2L, l = 1L
  )
  tt <- seq(0, 0.1, by = 0.01)
  obs <- observed_trajectory(tt, cbind(sin(tt), cos(tt)))
  for (method in c("mask", "reduce")) {
    nt <- aci_conditional(2L, method)
    generic <- .compile_cgns_complete(model, obs, nt)
    batch <- .compile_cgns_run(model, obs, nt)
    expect_identical(batch$realization, "affine_batch")
    expect_identical(batch$likelihood_idx, generic$likelihood_idx)
    for (field in names(generic$coefficients))
      expect_identical(batch$coefficients[[field]],
                       generic$coefficients[[field]], info = method)
  }
})


test_that("the affine batch realiser reproduces the generic arrays on every shape", {
  ## The batch realiser assembles its arrays from flat contiguous writes and
  ## differences a whole coefficient block at a time.  The generic one-pass
  ## realiser, which evaluates the model's own coefficient closures, is the
  ## reference for all of it.
  shapes <- list(
    list(k = 1L, l = 1L), list(k = 3L, l = 1L),
    list(k = 1L, l = 3L), list(k = 3L, l = 3L), list(k = 2L, l = 3L)
  )
  for (sh in shapes) {
    k <- sh$k
    l <- sh$l
    model <- aci_model_from_affine(
      f_full = function(t, x, y)
        as.numeric(-0.4 * x + 0.1 * sin(t + seq_len(k)) +
                   outer(seq_len(k) / k, seq_len(l) / l) %*% y),
      g_full = function(t, x, y)
        as.numeric(-0.6 * y + 0.05 * cos(t) * sum(x) +
                   outer(seq_len(l), seq_len(l), function(a, b)
                     0.02 * (a - b)) %*% y),
      Sx = function(t, x) diag(0.5 + 0.02 * abs(x), k, k),
      Sy_hidden = function(t, x) diag(0.7 + 0.01 * t, l, l),
      Sy_shared = function(t, x) matrix(0.01 * (1 + t), l, k),
      k = k, l = l
    )
    tt <- seq(0, 0.3, by = 0.005)
    xx <- matrix(0, length(tt), k)
    for (i in seq_len(k)) xx[, i] <- 0.3 * sin(3 * tt + i) - 0.1 * i
    obs <- observed_trajectory(tt, xx)
    lab <- sprintf("k=%d l=%d", k, l)

    generic <- .realise_cgns_grid_once(model, as_obs(obs))
    batch <- .realise_affine_cgns_grid_once(model, as_obs(obs))
    for (field in names(generic))
      expect_identical(batch[[field]], generic[[field]],
                       info = paste(lab, field))
    expect_identical(attr(batch, "path_cross_noise"),
                     attr(generic, "path_cross_noise"), info = lab)
    expect_identical(dim(batch$Lx), c(k, l, length(tt)), info = lab)
    expect_identical(dim(batch$fx), c(length(tt), k), info = lab)
    expect_identical(dim(batch$fy), c(length(tt), l), info = lab)
    expect_identical(dim(batch$gyx), c(l, k, length(tt)), info = lab)

    ## and slice by slice against a direct evaluation of the affine drifts
    d <- .cgns_realizer_descriptor(model)
    for (j in c(1L, 7L, length(tt))) {
      base_f <- as.numeric(d$spec$f_full(tt[j], xx[j, ], numeric(l)))
      cols <- vapply(seq_len(l), function(i)
        as.numeric(d$spec$f_full(tt[j], xx[j, ],
                                 replace(numeric(l), i, 1))) - base_f,
        numeric(k))
      expect_identical(as.numeric(batch$Lx[, , j]), as.numeric(cols),
                       info = paste(lab, j))
      expect_identical(batch$fx[j, ], base_f, info = paste(lab, j))
    }
  }
})


test_that("the whole-path cross-noise reduction matches the per-slice test", {
  make <- function(gyx_scale) {
    N1 <- 6L
    gxx <- array(0, c(2L, 2L, N1))
    gyy <- array(0, c(2L, 2L, N1))
    gyx <- array(0, c(2L, 2L, N1))
    for (j in seq_len(N1)) {
      gxx[, , j] <- diag(c(1, 2)) * j
      gyy[, , j] <- diag(c(3, 0.5)) * j
      gyx[, , j] <- gyx_scale[j] * matrix(c(1, 0, 0, 1), 2, 2)
    }
    list(gxx = gxx, gyy = gyy, gyx = gyx)
  }
  per_slice <- function(co) {
    N1 <- dim(co$gxx)[3L]
    any(vapply(seq_len(N1), function(j) .has_cross_noise(list(
      gxx = matrix(co$gxx[, , j], 2, 2),
      gyy = matrix(co$gyy[, , j], 2, 2),
      gyx = matrix(co$gyx[, , j], 2, 2))), logical(1)))
  }
  tol <- 100 * .Machine$double.eps
  ## exactly zero cross-noise; just under the tolerance on every slice; over
  ## the tolerance on one interior slice only; over on the last slice only.
  scales <- list(
    rep(0, 6L),
    tol * sqrt(seq_len(6L) * 2 * seq_len(6L) * 3) * 0.5,
    replace(rep(0, 6L), 4L, 1),
    replace(rep(0, 6L), 6L, 1e-12)
  )
  for (i in seq_along(scales)) {
    co <- make(scales[[i]])
    expect_identical(.compiled_path_has_cross_noise(co), per_slice(co),
                     info = paste("case", i))
  }
  ## non-finite scale is rejected by both forms
  co <- make(rep(1, 6L))
  co$gyy[1L, 1L, 3L] <- Inf
  expect_identical(.compiled_path_has_cross_noise(co), per_slice(co))
})


test_that("the precision path equals a per-slice solve on both branches", {
  for (k in c(1L, 2L, 4L)) {
    N <- 9L
    gxx <- array(0, c(k, k, N + 1L))
    for (j in seq_len(N + 1L)) {
      A <- matrix(0.05 * (seq_len(k * k) + j), k, k)
      gxx[, , j] <- diag(k) * (1 + 0.1 * j) + A %*% t(A)
    }
    got <- .compiled_precision_path(gxx, N)
    expect_identical(dim(got), c(k, k, N), info = paste("k", k))
    for (j in seq_len(N)) {
      gram <- matrix(gxx[, , j], k, k)
      ## These Grams are dense, so the unmasked route's chol2inv() and the
      ## definition's two triangular solves differ in the last bit or two.
      ## The equality carries the budget; the identity below keeps this a pin
      ## on the route the function actually takes, so it still fails if the
      ## guarded pair is ever entered on a Gram the fast pair should handle.
      expect_equal(matrix(got[, , j], k, k),
                   chol_solve(gram, diag(k), "gxx"),
                   tolerance = 1e-12, info = paste("k", k, "j", j))
      expect_identical(matrix(got[, , j], k, k),
                       chol2inv(safe_chol(gram, "gxx")) + 0,
                       info = paste("k", k, "j", j))
    }
    if (k > 1L) {
      target <- seq_len(k - 1L)
      gotm <- .compiled_precision_path(gxx, N, target = target)
      for (j in seq_len(N))
        expect_identical(matrix(gotm[, , j], k, k),
                         masked_ginv(matrix(gxx[, , j], k, k), target),
                         info = paste("masked k", k, "j", j))
    }
  }
})


test_that("a diagonal Gram inverts to the same bits on either route", {
  ## This is the structural reason the chol2inv route costs the ACI_code scope
  ## nothing: every observation Gram reachable there is diagonal, and on a
  ## diagonal Gram the Cholesky factor is exact, so chol2inv() and two
  ## triangular solves against an identity are the same divisions.
  ##
  ## The comparison is made on the serialized bytes as well as with
  ## `identical()`, because `identical()` and a max-absolute-difference both
  ## report `0 == -0` and `chol2inv()` writes `-0.0` into the off-diagonal
  ## zeros that the triangular route writes as `+0.0`. The `+ 0` in
  ## `.compiled_precision_path()`'s fast pair is what makes this pass.
  bitwise <- function(a, b) identical(serialize(a, NULL, version = 3L),
                                      serialize(b, NULL, version = 3L))
  for (k in c(1L, 2L, 3L, 5L)) {
    N <- 12L
    gxx <- array(0, c(k, k, N + 1L))
    for (j in seq_len(N + 1L))
      gxx[, , j] <- diag(seq_len(k) * (0.3 + 0.05 * j), k)
    got <- .compiled_precision_path(gxx, N)
    for (j in seq_len(N)) {
      gram <- matrix(gxx[, , j], k, k)
      ref <- chol_solve(gram, diag(k), "gxx")
      expect_identical(matrix(got[, , j], k, k), ref,
                       info = paste("diagonal k", k, "j", j))
      expect_true(bitwise(matrix(got[, , j], k, k), ref),
                  info = paste("diagonal bytes k", k, "j", j))
    }
  }
  ## and the realised ENSO Gram is in fact diagonal, on every slice
  skip_if_not(exists("aci_enso_model"))
  m <- aci_enso_model(hidden = c("u", "hW", "tau"), variant = "aci_code")
  on <- m$meta$vars$observed
  set.seed(4L)
  n <- 121L
  x <- matrix(0.2 * sin(seq_len(n * length(on)) / 7), n, length(on))
  ob <- observed_trajectory(seq(0, by = 1 / 12, length.out = n), x)
  b <- .compile_cgns_run(m, ob)
  g <- b$coefficients$gxx
  off <- max(vapply(seq_len(dim(g)[3L]), function(j) {
    G <- matrix(g[, , j], b$k, b$k)
    max(abs(G - diag(diag(G), b$k)))
  }, 1))
  expect_identical(off, 0)
})


test_that("the unmasked precision route is chol2inv, within budget of a solve", {
  ## Dense, deliberately ill-conditioned Grams: the only surface on which the
  ## two routes can disagree at all. The budget is 1e-13 RELATIVE to the
  ## magnitude of the inverse, measured against the guarded per-slice
  ## definition. O3 (scratch/shootout-0.1.0/overdrive/o3-chol2inv) measured the
  ## worst deviation at 9.6e-14 relative on a kappa 7.5e5 fixture.
  worst <- 0
  for (k in c(2L, 3L, 5L)) {
    for (lk in c(2, 4, 6)) {
      set.seed(700L + 10L * k + lk)
      N <- 12L
      ev <- 10^seq(0, -lk, length.out = k)
      gxx <- array(0, c(k, k, N + 1L))
      for (j in seq_len(N + 1L)) {
        Q <- qr.Q(qr(matrix(stats::rnorm(k * k), k, k)))
        S <- Q %*% diag(ev * (1 + 0.01 * j), k) %*% t(Q)
        gxx[, , j] <- (S + t(S)) / 2
      }
      got <- .compiled_precision_path(gxx, N)
      for (j in seq_len(N)) {
        gram <- matrix(gxx[, , j], k, k)
        ref <- chol_solve(gram, diag(k), "gxx")
        new <- matrix(got[, , j], k, k)
        ## the route pin: exactly the inverse chol2inv() gives on that factor
        expect_identical(new, chol2inv(safe_chol(gram, "gxx")) + 0,
                         info = paste("k", k, "log10kappa", lk, "j", j))
        ## the budget: 1e-13 relative to the scale of the inverse
        rel <- max(abs(new - ref)) / max(abs(ref))
        worst <- max(worst, rel)
        expect_lt(rel, 1e-13)
        ## chol2inv() is exactly symmetric where a triangular solve is not
        expect_identical(new, t.default(new),
                         info = paste("symmetry k", k, "j", j))
      }
    }
  }
  expect_lt(worst, 1e-13)
})


test_that("auto zero cross-channels keep the partner width they were built with", {
  probes <- 0L
  model <- aci_model(
    Lx = function(t, x) matrix(0.5, 1, 1),
    fx = function(t, x) -0.3 * x,
    Ly = function(t, x) matrix(-0.7, 1, 1),
    fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(c(0.4, 0.1), 1, 2),
    Sy2 = function(t, x) { probes <<- probes + 1L; matrix(c(0.8, 0.2, 0.05), 1, 3) },
    k = 1, l = 1
  )
  ## the auto Sx2 must be the k by ncol(Sy2) zero block at every point
  for (pt in list(c(0, 0.2), c(1.7, -3), c(-2.5, 11))) {
    expect_identical(model$Sx2(pt[1L], pt[2L]), matrix(0, 1L, 3L))
    expect_identical(model$Sy1(pt[1L], pt[2L]), matrix(0, 1L, 2L))
  }
  ## and it must not re-evaluate the partner to learn that width again
  before <- probes
  for (i in 1:20) model$Sx2(0.1 * i, 0.5)
  expect_identical(probes, before)

  ## a mid-path change in the partner's channel count is still rejected: the
  ## realiser compares the zero block's width against the partner every step.
  wobbly <- aci_model_from_affine(
    f_full = function(t, x, y) -0.3 * x + 0.5 * y,
    g_full = function(t, x, y) -0.7 * y,
    Sx = function(t, x) matrix(0.4, 1, 1),
    Sy_hidden = function(t, x)
      if (abs(t - 0.02) < 1e-12) matrix(c(0.8, 0.1), 1, 2) else matrix(0.8, 1, 1),
    k = 1, l = 1
  )
  obs <- observed_trajectory(seq(0, 0.05, by = 0.01), sin(seq(0, 0.05, by = 0.01)))
  expect_error(.compile_cgns_run(wobbly, obs),
               class = "aci_error_model_contract")
})


test_that("the affine batch realiser rejects a drift whose width changes on the grid", {
  bad <- aci_model_from_affine(
    f_full = function(t, x, y)
      if (abs(t - 0.02) < 1e-12 && y[1L] != 0) c(-0.3 * x + 0.5 * y, 0)
      else -0.3 * x + 0.5 * y,
    g_full = function(t, x, y) -0.7 * y,
    Sx = function(t, x) matrix(0.4, 1, 1),
    Sy_hidden = function(t, x) matrix(0.8, 1, 1),
    k = 1, l = 1
  )
  obs <- observed_trajectory(seq(0, 0.05, by = 0.01), sin(seq(0, 0.05, by = 0.01)))
  expect_error(.compile_cgns_run(bad, obs), class = "aci_error_model_contract")
  expect_error(.compile_cgns_run(bad, obs), "incompatible dimensions at index 3")
})
