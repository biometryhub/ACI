compiled_matrix_setup <- function() {
  model <- aci_model(
    Lx = function(t, x) matrix(c(
      0.7 + 0.05 * x[1], -0.2,
      0.1, 0.5 + 0.03 * x[2]
    ), 2, 2, byrow = TRUE),
    fx = function(t, x) c(-0.3 * x[1] + 0.1, -0.2 * x[2] - 0.05),
    Ly = function(t, x) matrix(c(-0.8, 0.1, -0.05, -0.6), 2, 2),
    fy = function(t, x) c(0.2 * sin(t) - 0.1 * x[1],
                           0.15 * cos(t) + 0.05 * x[2]),
    Sx1 = function(t, x) matrix(c(0.7, 0.1, 0, 0.6), 2, 2),
    Sy1 = function(t, x) matrix(c(0.04, 0, 0.01, 0.03), 2, 2),
    Sy2 = function(t, x) matrix(c(0.8, 0.05, 0, 0.7), 2, 2),
    k = 2, l = 2, name = "compiled-matrix-test"
  )
  t <- seq(0, 0.2, by = 0.005)
  x <- cbind(0.3 + 0.1 * sin(3 * t), -0.2 + 0.08 * cos(2 * t))
  obs <- observed_trajectory(t, x)
  init <- list(
    mean = c(0.2, -0.1),
    cov = matrix(c(0.4, 0.03, 0.03, 0.3), 2, 2)
  )
  list(
    model = model, obs = obs, init = init,
    bundle = .compile_cgns_complete(model, obs)
  )
}


# Frozen from the untouched 0.0.21 matrix equations on indices 1, 11, 21 and
# 41. The implicit stepper is an aci policy with no supplied MATLAB
# implementation, so this is a version regression, not a MATLAB-parity claim.
# Explicit matrix execution has an independent MATLAB grade in test-19.
.implicit_0021_regression <- list(
  `1` = list(
    filter_mean = matrix(c(
      0.2, 0.19518629327433684, 0.19089427639676157, 0.18302454533278029,
      -0.1, -0.08887487590904845, -0.07863256806157197,
      -0.06054855899486391
    ), 4, 2),
    filter_cov = matrix(c(
      0.4, 0.03, 0.3,
      0.39068800740510845, 0.031740314040495154, 0.30238617870456463,
      0.38245029120445495, 0.033247455348803273, 0.30459787226032314,
      0.36869065836718162, 0.035681259828022199, 0.30854253733223302
    ), 3, 4),
    loglik = 173.05021205003564,
    smoother_mean = matrix(c(
      0.21311499262957756, 0.20495311095082192, 0.19724694729561756,
      0.18302454533278029, -0.10241206734558872, -0.091341203043578137,
      -0.080696356362158689, -0.060548558994863909
    ), 4, 2),
    smoother_cov = matrix(c(
      0.37468439794488873, 0.032212163359455454, 0.28748580260762208,
      0.37146693581502577, 0.033364233703737728, 0.29247838920445318,
      0.3694069766583708, 0.034318447707379608, 0.29762726678567764,
      0.36869065836718162, 0.035681259828022199, 0.30854253733223302
    ), 3, 4)
  ),
  `2` = list(
    filter_mean = matrix(c(
      0.2, 0.19517649710550888, 0.19087676153397423, 0.18299593957408755,
      -0.1, -0.088861812155334882, -0.07860831737757859,
      -0.060506963144165246
    ), 4, 2),
    filter_cov = matrix(c(
      0.4, 0.03, 0.3,
      0.39073218096246309, 0.031740001777592647, 0.30241808861618813,
      0.38253430050577525, 0.033247046497050584, 0.30465896453580232,
      0.36884310337743698, 0.035681131787158704, 0.30865469249645616
    ), 3, 4),
    loglik = 173.05019787755526,
    smoother_mean = matrix(c(
      0.2130534208734006, 0.20489835941063678, 0.19719996179049726,
      0.18299593957408755, -0.10232570812645503, -0.091264706752629404,
      -0.080630643274903963, -0.060506963144165246
    ), 4, 2),
    smoother_cov = matrix(c(
      0.37482432193594412, 0.032215780491834843, 0.28763319685040123,
      0.37162178961148323, 0.033366319778633052, 0.29262557540623402,
      0.36956935618971398, 0.034319324175923108, 0.29776870095394126,
      0.36884310337743698, 0.035681131787158704, 0.30865469249645616
    ), 3, 4)
  )
)


test_that("compiled matrix dispatch covers explicit substepping", {
  ds <- compiled_matrix_setup()
  for (nsub in c(1L, 2L)) {
    direct <- .cgns_filter_matrix_compiled(
      ds$bundle, ds$init, stepper = "explicit", nsub = nsub
    )
    dispatched <- .cgns_filter_compiled(
      ds$bundle, ds$init, stepper = "explicit", nsub = nsub
    )
    expect_identical(dispatched, direct)
    expect_true(all(vapply(seq_len(ds$bundle$N1), function(j)
      !is.null(tryCatch(chol(direct$cov[, , j]), error = function(e) NULL)),
      logical(1))))
    expect_identical(direct$meta$stepper, "explicit")
    expect_identical(direct$meta$nsub, nsub)
  }
})


test_that("implicit matrix filter and smoother retain 0.0.21 policy", {
  ds <- compiled_matrix_setup()
  idx <- c(1L, 11L, 21L, 41L)
  for (nsub in c(1L, 2L)) {
    reference <- .implicit_0021_regression[[as.character(nsub)]]
    filter <- .cgns_filter_compiled(
      ds$bundle, ds$init, stepper = "implicit", nsub = nsub
    )
    smoother <- .cgns_smoother_compiled(ds$bundle, filter)
    filter_cov <- vapply(idx, function(i) c(
      filter$cov[1, 1, i], filter$cov[1, 2, i], filter$cov[2, 2, i]
    ), numeric(3))
    smoother_cov <- vapply(idx, function(i) c(
      smoother$cov[1, 1, i], smoother$cov[1, 2, i], smoother$cov[2, 2, i]
    ), numeric(3))
    expect_equal(unname(filter$mean[idx, ]), reference$filter_mean,
                 tolerance = 1e-14)
    expect_equal(unname(filter_cov), reference$filter_cov, tolerance = 1e-14)
    expect_equal(filter$meta$loglik, reference$loglik, tolerance = 1e-14)
    expect_equal(unname(smoother$mean[idx, ]), reference$smoother_mean,
                 tolerance = 1e-14)
    expect_equal(unname(smoother_cov), reference$smoother_cov,
                 tolerance = 1e-14)
    expect_identical(smoother$meta$route, "backward_ode_correlated")
  }
})


test_that("compiled dispatcher covers scalar implicit execution", {
  model <- aci_dyad_model()
  sim <- simulate(model, seed = 71, T = 0.2, dt = 0.002, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  bundle <- .compile_cgns_run(model, obs)
  matrix_route <- .cgns_filter_matrix_compiled(
    bundle, init, stepper = "implicit", nsub = 2L
  )
  dispatched <- .cgns_filter_compiled(
    bundle, init, stepper = "implicit", nsub = 2L
  )
  expect_identical(dispatched, matrix_route)
  expect_true(all(dispatched$cov > 0))
})


test_that("complete compiled matrix execution retains the aci result contract", {
  ds <- compiled_matrix_setup()
  filter <- .cgns_filter_compiled(ds$bundle, ds$init)
  smoother <- .cgns_smoother_compiled(ds$bundle, filter)
  metric <- aci_metric(smoother, filter, decompose = TRUE)
  compiled <- .aci_cgns_compiled(ds$bundle, ds$init)
  public <- aci(ds$model, ds$obs, init = ds$init)

  expect_s3_class(compiled, "aci_result")
  expect_equal(compiled$aci, metric$total, tolerance = 0)
  expect_equal(compiled$signal, metric$signal, tolerance = 0)
  expect_equal(compiled$dispersion, metric$dispersion, tolerance = 0)
  expect_identical(compiled$paths$filter, filter)
  expect_identical(compiled$paths$smoother, smoother)
  expect_identical(compiled$handles$model, ds$model)
  expect_identical(compiled$handles$obs, ds$obs)
  expect_equal(public$aci, compiled$aci, tolerance = 0)
  expect_equal(public$paths$filter$meta$loglik,
               compiled$paths$filter$meta$loglik, tolerance = 0)
})


test_that("matrix filter honours loglik = FALSE without moving the moments", {
  ds <- compiled_matrix_setup()
  carried <- setdiff(names(.cgns_filter_compiled(ds$bundle, ds$init)$meta),
                     "loglik")
  for (stepper in c("explicit", "implicit")) {
    for (nsub in c(1L, 2L)) {
      scored <- .cgns_filter_compiled(ds$bundle, ds$init, stepper = stepper,
                                      nsub = nsub)
      skipped <- .cgns_filter_compiled(ds$bundle, ds$init, stepper = stepper,
                                       nsub = nsub, loglik = FALSE)
      expect_identical(skipped$t, scored$t)
      expect_identical(skipped$mean, scored$mean)
      expect_identical(skipped$cov, scored$cov)
      expect_false(is.null(scored$meta$loglik))
      expect_null(skipped$meta$loglik)
      expect_identical(skipped$meta[carried], scored$meta[carried])
    }
  }
  ## The frozen 0.0.21 values still pin the scored arm on the implicit path.
  for (nsub in c(1L, 2L))
    expect_equal(
      .cgns_filter_compiled(ds$bundle, ds$init, stepper = "implicit",
                            nsub = nsub)$meta$loglik,
      .implicit_0021_regression[[as.character(nsub)]]$loglik,
      tolerance = 1e-14
    )
  expect_error(.cgns_filter_compiled(ds$bundle, ds$init, loglik = "yes"),
               class = "aci_error_dims")
})


test_that("aci(loglik = FALSE) leaves every matrix metric untouched", {
  ds <- compiled_matrix_setup()
  scored <- aci(ds$model, ds$obs, init = ds$init)
  skipped <- aci(ds$model, ds$obs, init = ds$init, loglik = FALSE)

  expect_identical(skipped$aci, scored$aci)
  expect_identical(skipped$signal, scored$signal)
  expect_identical(skipped$dispersion, scored$dispersion)
  expect_identical(skipped$paths$filter$mean, scored$paths$filter$mean)
  expect_identical(skipped$paths$filter$cov, scored$paths$filter$cov)
  expect_identical(skipped$paths$smoother$mean, scored$paths$smoother$mean)
  expect_identical(skipped$paths$smoother$cov, scored$paths$smoother$cov)
  expect_false(is.null(scored$paths$filter$meta$loglik))
  expect_null(skipped$paths$filter$meta$loglik)
  ## A loglik-free filter is still a valid input to the public smoother.
  expect_identical(
    aci_smoother(ds$model, ds$obs, filter = skipped$paths$filter)$cov,
    scored$paths$smoother$cov
  )
})


test_that(".sym_floor() reproduces spd_floor(sym()) on the recursion inputs", {
  set.seed(31L)
  cases <- list(
    spd = crossprod(matrix(stats::rnorm(9), 3, 3)) + diag(3),
    near_asym = matrix(c(2, 0.3 + 1e-13, 0.1, 0.3, 1.5, 0.2,
                         0.1, 0.2, 1.1), 3, 3),
    singular = diag(c(1, 1, 0)),
    indefinite = matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), 3, 3),
    scalar = matrix(0.25, 1, 1),
    negative_scalar = matrix(-0.25, 1, 1)
  )
  for (nm in names(cases)) {
    m <- cases[[nm]]
    expect_identical(.sym_floor(m), spd_floor(sym(m)), info = nm)
  }
  expect_error(.sym_floor(matrix(c(1, NA, NA, 1), 2, 2)),
               class = "aci_error_spd")
})


test_that("the matrix kernels read the same coefficients as .compiled_co()", {
  ds <- compiled_matrix_setup()
  b <- .compile_cgns_run(ds$model, ds$obs)
  ## The filter and smoother slice the coefficient arrays inline rather than
  ## through .compiled_co(); the two readings must be the same matrices.
  for (j in c(1L, 2L, b$N)) {
    co <- .compiled_co(b, j)
    Lxj <- b$coefficients$Lx[, , j]; dim(Lxj) <- c(b$k, b$l)
    Lyj <- b$coefficients$Ly[, , j]; dim(Lyj) <- c(b$l, b$l)
    gyyj <- b$coefficients$gyy[, , j]; dim(gyyj) <- c(b$l, b$l)
    gyxj <- b$coefficients$gyx[, , j]; dim(gyxj) <- c(b$l, b$k)
    gxxj <- b$coefficients$gxx[, , j]; dim(gxxj) <- c(b$k, b$k)
    Gi <- b$coefficients$gxx_weight[, , j]; dim(Gi) <- c(b$k, b$k)
    expect_identical(Lxj, co$Lx)
    expect_identical(Lyj, co$Ly)
    expect_identical(gyyj, co$gyy)
    expect_identical(gyxj, co$gyx)
    expect_identical(gxxj, co$gxx)
    expect_identical(Gi, .compiled_ginv(b, j))
    expect_identical(as.numeric(b$coefficients$fx[j, ]), co$fx)
    expect_identical(as.numeric(b$coefficients$fy[j, ]), co$fy)
  }
})


test_that("the smoother's filtered-covariance inverse is chol2inv on the factor", {
  ## The backward recursion inverts the filtered covariance with chol2inv()
  ## rather than two triangular solves against an identity. For l == 1 the two
  ## are one division and agree bit for bit; above that they differ in the last
  ## bit or two and the recursion is budgeted for it.
  set.seed(11L)
  for (l in c(1L, 2L, 3L, 5L)) {
    A <- matrix(stats::rnorm(l * l), l, l)
    R <- crossprod(A) + l * diag(l)
    ch <- chol(R)
    tri <- backsolve(ch, forwardsolve(t(ch), diag(l)))
    if (l == 1L)
      expect_identical(chol2inv(ch), tri, info = paste("l", l))
    else
      expect_equal(chol2inv(ch), tri, tolerance = 1e-12, info = paste("l", l))
    ## chol2inv() is exactly symmetric; the triangular route need not be
    expect_identical(chol2inv(ch), t.default(chol2inv(ch)), info = paste("l", l))
  }
})
