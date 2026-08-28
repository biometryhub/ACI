compiled_dyad_setup <- function(T = 0.4, dt = 0.001, seed = 41) {
  model <- aci_dyad_model()
  sim <- simulate(model, seed = seed, T = T, dt = dt, burn_in = 0)
  obs <- as_obs(sim)
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  list(model = model, obs = obs, init = init,
       bundle = .compile_cgns_complete(model, obs))
}


expect_scalar_paths_equal <- function(actual, expected, tolerance = 1e-13) {
  expect_equal(actual$t, expected$t, tolerance = 0)
  expect_equal(actual$mean, expected$mean, tolerance = tolerance)
  expect_equal(actual$cov, expected$cov, tolerance = tolerance)
  expect_identical(dim(actual$mean), dim(expected$mean))
  expect_identical(dim(actual$cov), dim(expected$cov))
}


capture_aci_warnings <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(
    value = value,
    classes = vapply(warnings, function(w) class(w)[1L], character(1)),
    messages = vapply(warnings, conditionMessage, character(1)))
}


test_that("compiled CGNS generic fallback realizes each grid point once", {
  calls <- 0L
  model <- aci_model(
    Lx = function(t, x) { calls <<- calls + 1L; matrix(1 + t + x[1], 1, 1) },
    fx = function(t, x) 0.2 - x,
    Ly = function(t, x) matrix(-0.7 + 0.1 * t, 1, 1),
    fy = function(t, x) sin(t) - 0.2 * x^2,
    Sx1 = function(t, x) matrix(0.4 + 0.01 * t, 1, 1),
    Sy1 = function(t, x) matrix(0.05 * x[1], 1, 1),
    Sy2 = function(t, x) matrix(0.8, 1, 1), k = 1, l = 1)
  obs <- observed_trajectory(seq(0, 0.1, by = 0.01),
                             matrix(seq(0.2, 0.3, length.out = 11), ncol = 1))
  calls <- 0L
  bundle <- .compile_cgns_complete(model, obs)

  expect_s3_class(bundle, "compiled_cgns")
  expect_identical(bundle$realization, "generic_closure_one_pass")
  expect_equal(calls, length(obs$t))
  expect_true(bundle$correlated_noise)
  expect_null(bundle$init)
  expect_identical(dim(bundle$rate), c(10L, 1L))
  expect_equal(bundle$rate[, 1], diff(obs$x[, 1]) / obs$dt, tolerance = 0)

  for (j in c(1L, 6L, 11L)) {
    co <- eval_coefs(model, obs$t[j], obs$x[j, ])
    expect_equal(matrix(bundle$coefficients$Lx[, , j], 1L, 1L),
                 co$Lx, tolerance = 0)
    expect_equal(bundle$coefficients$fx[j, ], co$fx, tolerance = 0)
    expect_equal(matrix(bundle$coefficients$Ly[, , j], 1L, 1L),
                 co$Ly, tolerance = 0)
    expect_equal(bundle$coefficients$fy[j, ], co$fy, tolerance = 0)
    expect_equal(matrix(bundle$coefficients$gxx[, , j], 1L, 1L),
                 co$gxx, tolerance = 0)
    expect_equal(matrix(bundle$coefficients$gyy[, , j], 1L, 1L),
                 co$gyy, tolerance = 0)
    expect_equal(matrix(bundle$coefficients$gyx[, , j], 1L, 1L),
                 co$gyx, tolerance = 0)
  }
})


test_that("directed dyad realiser fills the same neutral compiled contract", {
  ds <- compiled_dyad_setup(T = 0.2, dt = 0.001)
  directed <- .compile_dyad_cgns(ds$model, ds$obs)
  expect_identical(ds$bundle$realization, "generic_closure_one_pass")
  expect_identical(directed$realization, "dyad_directed")
  for (field in names(ds$bundle$coefficients))
    expect_identical(directed$coefficients[[field]],
                     ds$bundle$coefficients[[field]])
  expect_identical(directed$scalar, ds$bundle$scalar)
  expect_no_error(.validate_compiled_cgns(
    directed, ds$model, ds$obs, scalar = TRUE))

  generic_result <- .aci_scalar_compiled(ds$bundle, ds$init)
  directed_result <- .aci_scalar_compiled(directed, ds$init)
  expect_equal(directed_result$aci, generic_result$aci, tolerance = 0)
  expect_equal(directed_result$paths$filter$mean,
               generic_result$paths$filter$mean, tolerance = 0)
  expect_equal(directed_result$paths$filter$cov,
               generic_result$paths$filter$cov, tolerance = 0)
  expect_equal(directed_result$paths$smoother$mean,
               generic_result$paths$smoother$mean, tolerance = 0)
  expect_equal(directed_result$paths$smoother$cov,
               generic_result$paths$smoother$cov, tolerance = 0)

  changed <- ds$model
  changed$fx <- function(t, x) -0.6 * x + 0.5
  expect_error(.compile_dyad_cgns(changed, ds$obs),
               class = "aci_error_compiled_contract")
  expect_error(.compile_dyad_cgns(
    ds$model, ds$obs, conditional = aci_conditional(1, "mask")),
    class = "aci_error_compiled_contract")
})


test_that("automatic realiser selection uses sealed constructor identity", {
  ds <- compiled_dyad_setup(T = 0.05, dt = 0.005)
  model <- ds$model
  expect_true(is.environment(.cgns_realizer_descriptor(model)))
  expect_identical(
    .compile_cgns_run(model, ds$obs)$realization,
    "dyad_directed"
  )

  metadata_only <- model
  metadata_only$name <- "a harmless presentation label"
  metadata_only$meta$params$d_x <- 999
  expect_identical(
    .compile_cgns_run(metadata_only, ds$obs)$realization,
    "dyad_directed"
  )

  changed <- model
  changed$fx <- function(t, x) -0.9 * x + 0.25
  expect_null(.cgns_realizer_descriptor(changed))
  expect_identical(
    .compile_cgns_run(changed, ds$obs)$realization,
    "generic_closure_one_pass"
  )

  expect_error(environment(model$fx)$p$d_x <- 999)
  roundtrip <- unserialize(serialize(model, NULL))
  expect_null(.cgns_realizer_descriptor(roundtrip))
  expect_identical(
    .compile_cgns_run(roundtrip, ds$obs)$realization,
    "generic_closure_one_pass"
  )
})


test_that("compiled identity excludes init but binds model, observations, and policy", {
  ds <- compiled_dyad_setup(T = 0.05, dt = 0.005)
  expect_no_error(.validate_compiled_cgns(
    ds$bundle, ds$model, ds$obs, scalar = TRUE))

  changed_model <- ds$model
  changed_model$fx <- function(t, x) -0.6 * x + 0.5
  expect_error(.validate_compiled_cgns(ds$bundle, changed_model, ds$obs),
               class = "aci_error_model_contract")
  changed_obs <- ds$obs; changed_obs$x[2, 1] <- changed_obs$x[2, 1] + 0.1
  expect_error(.validate_compiled_cgns(ds$bundle, ds$model, changed_obs),
               class = "aci_error_obs_contract")
  expect_error(.validate_compiled_cgns(
    ds$bundle, conditional = aci_conditional(1, "mask")),
    class = "aci_error_nontarget")

  f1 <- .cgns_filter_scalar(ds$bundle, ds$init)
  f2 <- .cgns_filter_scalar(
    ds$bundle, list(mean = -1, cov = matrix(0.5, 1, 1)))
  expect_false(isTRUE(all.equal(f1$mean, f2$mean)))
  expect_identical(f1$meta$source_model, ds$model)
  expect_identical(f1$meta$obs_x, ds$obs$x)
})


test_that("compiled bundles reject internal mutation and mismatched metric grids", {
  ds <- compiled_dyad_setup(T = 0.05, dt = 0.005)
  mutations <- list(
    function(x) { x$t[2] <- x$t[2] + 0.001; x },
    function(x) { x$x[2, 1] <- x$x[2, 1] + 0.1; x },
    function(x) { x$dt <- x$dt * 2; x },
    function(x) { x$model$k <- 2L; x },
    function(x) { x$source_obs$x[2, 1] <- x$source_obs$x[2, 1] + 0.1; x },
    function(x) { x$coefficients$Lx <- x$coefficients$Lx[, , -1]; x },
    function(x) { x$scalar$Lx[1] <- x$scalar$Lx[1] + 1; x })
  for (mutate in mutations)
    expect_error(.validate_compiled_cgns(
      mutate(ds$bundle), conditional = NULL, scalar = TRUE),
      class = "aci_error_compiled_contract")

  filter <- .cgns_filter_scalar(ds$bundle, ds$init)
  smoother <- .cgns_smoother_scalar(ds$bundle, filter)
  shifted_filter <- filter; shifted_filter$t <- shifted_filter$t + 1
  shifted_smoother <- smoother; shifted_smoother$t <- shifted_smoother$t + 1
  expect_error(.gaussian_kl_path_scalar(
    ds$bundle, shifted_smoother, shifted_filter),
    class = "aci_error_dims")
})


test_that("prescribed forcing compiles after resolution to a scalar contract", {
  model <- aci_model(
    Lx = function(t, x) matrix(c(1, 0.5), 2, 1),
    fx = function(t, x) c(-x[1] + 0.2 * x[2], -0.4 * x[2]),
    Ly = function(t, x) matrix(-0.8, 1, 1),
    fy = function(t, x) 0.1 * x[2],
    Sx1 = function(t, x) diag(c(0.4, 0.3)),
    Sy2 = function(t, x) matrix(0.7, 1, 1), k = 2, l = 1)
  sim <- simulate(model, seed = 9, T = 0.1, dt = 0.005,
                  ic = list(x0 = c(0, 0), y0 = 0.2))
  nt <- aci_conditional(2, "reduce")
  init <- list(mean = 0.2, cov = matrix(0.3, 1, 1))
  bundle <- .compile_cgns_complete(model, sim$obs, conditional = nt)
  public <- aci_filter(model, sim$obs, init = init, conditional = nt)
  scalar <- .cgns_filter_scalar(bundle, init)

  expect_identical(c(bundle$k, bundle$l), c(1L, 1L))
  expect_identical(bundle$source_model, model)
  expect_identical(bundle$conditional, nt)
  expect_scalar_paths_equal(scalar, public)
  expect_equal(scalar$meta$loglik, public$meta$loglik, tolerance = 1e-13)
})


test_that("private scalar kernels and production public routing agree", {
  ds <- compiled_dyad_setup()
  for (nsub in c(1L, 3L)) {
    public_filter <- aci_filter(ds$model, ds$obs, init = ds$init,
                                stepper = "explicit", nsub = nsub)
    scalar_filter <- .cgns_filter_scalar(ds$bundle, ds$init, nsub = nsub)
    expect_scalar_paths_equal(scalar_filter, public_filter)
    expect_equal(scalar_filter$meta$loglik, public_filter$meta$loglik,
                 tolerance = 1e-13)
    expect_identical(scalar_filter$meta$stepper, "explicit")
    expect_identical(scalar_filter$meta$nsub, nsub)

    public_smoother <- aci_smoother(ds$model, ds$obs, filter = public_filter)
    scalar_smoother <- .cgns_smoother_scalar(ds$bundle, scalar_filter)
    expect_scalar_paths_equal(scalar_smoother, public_smoother)
    expect_identical(scalar_smoother$meta$route, public_smoother$meta$route)
    expect_identical(scalar_smoother$mean[nrow(scalar_smoother$mean), ],
                     scalar_filter$mean[nrow(scalar_filter$mean), ])

    for (decompose in c(FALSE, TRUE)) {
      public_metric <- aci_metric(public_smoother, public_filter,
                                  decompose = decompose)
      scalar_metric <- .gaussian_kl_path_scalar(
        ds$bundle, scalar_smoother, scalar_filter, decompose = decompose)
      expect_equal(scalar_metric, public_metric, tolerance = 1e-13)
    }
  }

  public <- aci(ds$model, ds$obs, init = ds$init)
  scalar <- .aci_scalar_compiled(ds$bundle, init = ds$init)
  expect_equal(scalar$aci, public$aci, tolerance = 1e-13)
  expect_equal(scalar$signal, public$signal, tolerance = 1e-13)
  expect_equal(scalar$dispersion, public$dispersion, tolerance = 1e-13)
  expect_identical(scalar$handles$model, ds$model)
  expect_identical(scalar$handles$obs, ds$obs)
  expect_identical(scalar$meta$smoother_scheme,
                   public$meta$smoother_scheme)
  expect_null(.aci_scalar_compiled(
    ds$bundle, init = ds$init, keep = "none")$paths)
})


test_that("scalar proof preserves warnings and covariance policy", {
  ds <- compiled_dyad_setup(T = 0.1, dt = 0.005)
  public_default <- capture_aci_warnings(aci_filter(ds$model, ds$obs))
  scalar_default <- capture_aci_warnings(.cgns_filter_scalar(ds$bundle))
  expect_identical(scalar_default$classes, public_default$classes)
  expect_identical(scalar_default$messages, public_default$messages)
  expect_identical(scalar_default$classes, "aci_warn_diffuse_init")
  expect_scalar_paths_equal(scalar_default$value, public_default$value)
  for (bad in c(0, -1))
    expect_error(.cgns_filter_scalar(
      ds$bundle, list(mean = 0, cov = matrix(bad, 1, 1))),
      class = "aci_error_spd")
  expect_error(.cgns_filter_scalar(ds$bundle, ds$init, nsub = 0),
               class = "aci_error_dims")

  expect_identical(.scalar_spd_floor(2), 2)
  expect_equal(.scalar_spd_floor(-2), 2e-12, tolerance = 0)
  expect_equal(.scalar_spd_floor(0), 1e-300, tolerance = 0)
  expect_error(.scalar_spd_floor(Inf), class = "aci_error_spd")

  # Updated for the strict covariance policy (decision D3): the stiff case
  # below floors the variance, so both routes are asked for the opt-in, and
  # both are also pinned to stop identically under the default.
  stiff <- compiled_dyad_setup(T = 0.2, dt = 0.005, seed = 55)
  big <- list(mean = 0, cov = matrix(100, 1, 1))
  public_strict <- tryCatch(aci_filter(stiff$model, stiff$obs, init = big),
                            error = function(e) e)
  scalar_strict <- tryCatch(.cgns_filter_scalar(stiff$bundle, big),
                            error = function(e) e)
  expect_s3_class(public_strict, "aci_error_covariance_not_spd")
  expect_s3_class(scalar_strict, "aci_error_covariance_not_spd")
  expect_identical(conditionMessage(scalar_strict),
                   conditionMessage(public_strict))
  expect_identical(scalar_strict$index, public_strict$index)
  expect_identical(scalar_strict$site, "filter_explicit")
  public_stiff <- capture_aci_warnings(
    aci_filter(stiff$model, stiff$obs, init = big, regularize = "floor"))
  scalar_stiff <- capture_aci_warnings(
    .cgns_filter_scalar(stiff$bundle, big, regularize = "floor"))
  expect_identical(scalar_stiff$classes, public_stiff$classes)
  expect_identical(scalar_stiff$messages, public_stiff$messages)
  expect_true("aci_warn_riccati_stiff" %in% scalar_stiff$classes)
  expect_scalar_paths_equal(scalar_stiff$value, public_stiff$value)
  expect_identical(scalar_stiff$value$meta$regularization,
                   public_stiff$value$meta$regularization)
  expect_true(public_stiff$value$meta$regularization$fired)

  dt_model <- aci_dyad_model(params = within(ds$model$meta$params, d_y <- 25))
  dt_obs <- observed_trajectory(
    c(0, 0.005), matrix(c(0, 0), ncol = 1))
  dt_init <- list(mean = 0, cov = matrix(0.01, 1, 1))
  dt_bundle <- .compile_cgns_complete(dt_model, dt_obs)
  public_dt <- capture_aci_warnings(
    aci_filter(dt_model, dt_obs, init = dt_init))
  scalar_dt <- capture_aci_warnings(
    .cgns_filter_scalar(dt_bundle, dt_init))
  expect_identical(scalar_dt$classes, public_dt$classes)
  expect_identical(scalar_dt$messages, public_dt$messages)
  expect_identical(scalar_dt$classes, "aci_warn_dt_stability")
  expect_scalar_paths_equal(scalar_dt$value, public_dt$value)

  floor_model <- aci_model(
    Lx = function(t, x) matrix(10, 1, 1),
    fx = function(t, x) 0,
    Ly = function(t, x) matrix(0, 1, 1),
    fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(0, 1, 1), k = 1, l = 1)
  floor_obs <- observed_trajectory(
    c(0, 0.02), matrix(c(0, 0), ncol = 1))
  floor_init <- list(mean = 0, cov = matrix(1, 1, 1))
  floor_bundle <- .compile_cgns_complete(floor_model, floor_obs)
  # Updated for D3: floor_model is built to drive the variance negative, so it
  # is the package's direct test of the floor value and now names the policy.
  expect_error(aci_filter(floor_model, floor_obs, init = floor_init),
               class = "aci_error_covariance_not_spd")
  public_floor <- capture_aci_warnings(
    aci_filter(floor_model, floor_obs, init = floor_init,
               regularize = "floor"))
  scalar_floor <- capture_aci_warnings(
    .cgns_filter_scalar(floor_bundle, floor_init, regularize = "floor"))
  expect_identical(scalar_floor$classes, public_floor$classes)
  expect_identical(scalar_floor$messages, public_floor$messages)
  expect_equal(scalar_floor$value$cov, public_floor$value$cov, tolerance = 0)
  expect_gt(scalar_floor$value$cov[1, 1, 2], 0)
  expect_equal(scalar_floor$value$cov[1, 1, 2], 1e-12, tolerance = 0)
  fr <- public_floor$value$meta$regularization
  expect_identical(fr$policy, "floor")
  expect_identical(fr$n_events, 1L)
  expect_identical(fr$eps, 1e-12)
  expect_identical(fr$sites$site, "filter_explicit")
  expect_identical(fr$sites$first_index, 2L)
  expect_equal(fr$sites$first_time, 0.02)
  expect_lt(fr$sites$worst_value, 0)

  warned <- FALSE
  expect_error(withCallingHandlers(
    .aci_scalar_compiled(ds$bundle, decompose = NA),
    warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
    class = "aci_error_dims")
  expect_false(warned)
})


test_that("compiled scalar path matches the bundled independent MATLAB P1 port", {
  ds <- compiled_dyad_setup(T = 0.3, dt = 0.001, seed = 333)
  p <- ds$model$meta$params
  golden <- golden_p1_moments(
    ds$obs$x[, 1], ds$obs$dt,
    list(gamma = p$gamma, F_x = p$f_x, d_x = p$d_x,
         F_y = p$f_y, d_y = p$d_y, s_x = p$s_x, s_y = p$s_y),
    mu0 = ds$init$mean, R0 = ds$init$cov[1, 1])
  filter <- .cgns_filter_scalar(ds$bundle, ds$init)
  smoother <- .cgns_smoother_scalar(ds$bundle, filter, validate = FALSE)
  metric <- .gaussian_kl_path_scalar(ds$bundle, smoother, filter)

  expect_equal(filter$mean[, 1], golden$fm, tolerance = 1e-13)
  expect_equal(filter$cov[1, 1, ], golden$fc, tolerance = 1e-13)
  expect_equal(smoother$mean[, 1], golden$sm, tolerance = 1e-13)
  expect_equal(smoother$cov[1, 1, ], golden$sc, tolerance = 1e-13)
  expect_equal(metric$total, golden$aci, tolerance = 1e-12)
  expect_equal(metric$signal, golden$signal, tolerance = 1e-12)
  expect_equal(metric$dispersion, golden$dispersion, tolerance = 1e-12)
})


test_that("time-varying correlated scalar coefficients retain route and parity", {
  model <- aci_model(
    Lx = function(t, x) matrix(0.7 + 0.1 * t, 1, 1),
    fx = function(t, x) 0.1 - 0.05 * x,
    Ly = function(t, x) matrix(-1.1 + 0.05 * t, 1, 1),
    fy = function(t, x) -0.05 + 0.03 * x,
    Sx1 = function(t, x) matrix(0.4, 1, 1),
    Sx2 = function(t, x) matrix(0.1 + 0.02 * t, 1, 1),
    Sy1 = function(t, x) matrix(0.2 * x[1], 1, 1),
    Sy2 = function(t, x) matrix(0.5, 1, 1), k = 1, l = 1)
  obs <- observed_trajectory(0.01 * 0:50,
    matrix(cumsum(c(0, rep(c(0.01, -0.006), 25))), ncol = 1))
  init <- list(mean = 0.2, cov = matrix(0.8, 1, 1))
  bundle <- .compile_cgns_complete(model, obs)
  public_filter <- aci_filter(model, obs, init = init, nsub = 2)
  scalar_filter <- .cgns_filter_scalar(bundle, init, nsub = 2)
  public_smoother <- aci_smoother(model, obs, filter = public_filter)
  scalar_smoother <- .cgns_smoother_scalar(bundle, scalar_filter)

  expect_true(bundle$correlated_noise)
  expect_identical(scalar_smoother$meta$route, "backward_ode_correlated")
  expect_scalar_paths_equal(scalar_filter, public_filter)
  expect_scalar_paths_equal(scalar_smoother, public_smoother)
})


test_that("terminal-only cross-noise is retained in scalar route metadata", {
  model <- aci_model(
    Lx = function(t, x) matrix(0.5, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-0.8, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.4, 1, 1),
    Sy1 = function(t, x) matrix(if (abs(t - 0.1) < 1e-14) 0.2 else 0,
                                1, 1),
    Sy2 = function(t, x) matrix(0.7, 1, 1), k = 1, l = 1)
  expect_false(model$meta$correlated_noise)
  obs <- observed_trajectory(seq(0, 0.1, by = 0.01), matrix(0, 11, 1))
  init <- list(mean = 0, cov = matrix(0.4, 1, 1))
  bundle <- .compile_cgns_complete(model, obs)
  filter <- .cgns_filter_scalar(bundle, init)
  smoother <- .cgns_smoother_scalar(bundle, filter)
  public <- aci_smoother(model, obs, filter = aci_filter(model, obs, init = init))

  expect_true(bundle$correlated_noise)
  expect_identical(smoother$meta$route, "backward_ode_correlated")
  expect_identical(smoother$meta$route, public$meta$route)
})


test_that("valid named supplied filters retain the terminal invariant", {
  ds <- compiled_dyad_setup(T = 0.05, dt = 0.005)
  filter <- .cgns_filter_scalar(ds$bundle, ds$init)
  colnames(filter$mean) <- "hidden"
  smoother <- expect_no_error(.cgns_smoother_scalar(ds$bundle, filter))
  expect_lt(abs(smoother$mean[nrow(smoother$mean), 1] -
                  filter$mean[nrow(filter$mean), 1]), 1e-12)
})


external_scalar_model <- function(name) {
  if (identical(name, "dyad_reference_head")) return(aci_dyad_model())
  aci_model(
    Lx = function(t, x) matrix(0.8 + 0.35 * sin(3 * t), 1, 1),
    fx = function(t, x) 0.2 * cos(2 * t) - 0.1 * x,
    Ly = function(t, x) matrix(-1.3, 1, 1),
    fy = function(t, x) 0.45 - 0.15 * x^2,
    Sx1 = function(t, x) matrix(0.5, 1, 1),
    Sx2 = function(t, x) matrix(0.3, 1, 1),
    Sy1 = function(t, x) matrix(0.4, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
}


test_that("compiled scalar paths retain full external MATLAB-fixture parity", {
  oracle_root <- Sys.getenv("ACI_ORACLE_PARITY_ROOT", unset = "")
  if (!nzchar(oracle_root) || !dir.exists(oracle_root))
    skip("Set ACI_ORACLE_PARITY_ROOT to run the full external MATLAB fixtures.")

  for (name in c("dyad_reference_head", "arbitrary_cross_noise")) {
    dataset_dir <- file.path(oracle_root, "datasets", name)
    report_path <- file.path(oracle_root, "reports", paste0("matlab_", name, ".csv"))
    expect_true(file.exists(file.path(dataset_dir, "meta.dcf")))
    expect_true(file.exists(file.path(dataset_dir, "arrays.csv")))
    expect_true(file.exists(report_path))
    meta <- read.dcf(file.path(dataset_dir, "meta.dcf"))
    arrays <- utils::read.csv(file.path(dataset_dir, "arrays.csv"))
    reference <- utils::read.csv(report_path, check.names = FALSE)
    num <- function(field) as.numeric(meta[1L, field])
    dt <- num("dt")
    obs <- observed_trajectory(seq(0, by = dt, length.out = nrow(arrays)),
                               matrix(arrays$x, ncol = 1))
    model <- external_scalar_model(name)
    init <- list(mean = num("mu0"), cov = matrix(num("R0"), 1, 1))
    bundle <- .compile_cgns_complete(model, obs)
    filter <- .cgns_filter_scalar(bundle, init)
    smoother <- .cgns_smoother_scalar(bundle, filter, validate = FALSE)
    metric <- .gaussian_kl_path_scalar(bundle, smoother, filter)

    expect_lt(max(abs(filter$mean[, 1] - reference$filter_mean)), 1e-12)
    expect_lt(max(abs(filter$cov[1, 1, ] - reference$filter_cov)), 1e-12)
    expect_lt(max(abs(smoother$mean[, 1] - reference$smoother_mean)), 1e-12)
    expect_lt(max(abs(smoother$cov[1, 1, ] - reference$smoother_cov)), 1e-12)
    expect_lt(max(abs(metric$total - reference$ACI_metric)), 1e-11)
    expect_lt(max(abs(metric$signal - reference$signal_part)), 1e-11)
    expect_lt(max(abs(metric$dispersion - reference$dispersion_part)), 1e-11)
  }
})


test_that("scalar filter honours loglik = FALSE without moving the moments", {
  ds <- compiled_dyad_setup()
  carried <- setdiff(names(.cgns_filter_scalar(ds$bundle, ds$init)$meta),
                     "loglik")
  for (nsub in c(1L, 3L)) {
    scored <- .cgns_filter_scalar(ds$bundle, ds$init, nsub = nsub)
    skipped <- .cgns_filter_scalar(ds$bundle, ds$init, nsub = nsub,
                                   loglik = FALSE)
    expect_identical(skipped$t, scored$t)
    expect_identical(skipped$mean, scored$mean)
    expect_identical(skipped$cov, scored$cov)
    expect_false(is.null(scored$meta$loglik))
    expect_null(skipped$meta$loglik)
    expect_identical(skipped$meta[carried], scored$meta[carried])
  }
})


test_that("aci_filter defaults to loglik = TRUE and skips it on request", {
  ds <- compiled_dyad_setup()
  default <- aci_filter(ds$model, ds$obs, init = ds$init)
  scored <- aci_filter(ds$model, ds$obs, init = ds$init, loglik = TRUE)
  skipped <- aci_filter(ds$model, ds$obs, init = ds$init, loglik = FALSE)

  expect_identical(default, scored)
  expect_equal(default$meta$loglik,
               .cgns_likelihood_scalar_kernel(ds$bundle, default$mean[, 1],
                                              default$cov[1, 1, ]),
               tolerance = 0)
  expect_identical(skipped$mean, default$mean)
  expect_identical(skipped$cov, default$cov)
  expect_null(skipped$meta$loglik)
  expect_identical(skipped$meta$likelihood_idx, default$meta$likelihood_idx)
  expect_error(aci_filter(ds$model, ds$obs, init = ds$init, loglik = NA),
               class = "aci_error_dims")
  expect_error(aci_filter(ds$model, ds$obs, init = ds$init, loglik = c(TRUE, TRUE)),
               class = "aci_error_dims")
})


test_that("aci(loglik = FALSE) leaves every scalar metric untouched", {
  ds <- compiled_dyad_setup()
  scored <- aci(ds$model, ds$obs, init = ds$init, keep = "paths")
  skipped <- aci(ds$model, ds$obs, init = ds$init, keep = "paths",
                 loglik = FALSE)

  expect_identical(skipped$aci, scored$aci)
  expect_identical(skipped$signal, scored$signal)
  expect_identical(skipped$dispersion, scored$dispersion)
  expect_identical(skipped$paths$filter$mean, scored$paths$filter$mean)
  expect_identical(skipped$paths$filter$cov, scored$paths$filter$cov)
  expect_identical(skipped$paths$smoother$mean, scored$paths$smoother$mean)
  expect_identical(skipped$paths$smoother$cov, scored$paths$smoother$cov)
  expect_false(is.null(scored$paths$filter$meta$loglik))
  expect_null(skipped$paths$filter$meta$loglik)
})
