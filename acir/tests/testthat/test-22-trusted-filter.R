## Ledger C2a: the trusted supplied-filter route in aci_smoother().
## Two claims are under test. The authenticated route and the full-validation
## route return the same smoother to the bit; and every way a supplied filter
## can stop being the object this namespace produced sends it back to the
## unchanged full-validation route.

trusted_scalar_setup <- function(T = 0.3, dt = 0.001, seed = 19) {
  model <- aci_dyad_model()
  obs <- as_obs(simulate(model, seed = seed, T = T, dt = dt, burn_in = 0))
  init <- list(mean = 2, cov = matrix(0.1, 1, 1))
  list(model = model, obs = obs, init = init,
       filter = aci_filter(model, obs, init = init))
}


trusted_matrix_setup <- function(T = 0.3, dt = 0.002, seed = 23) {
  model <- aci_model(
    Lx = function(t, x) matrix(c(0.8, -0.3, 0.15, 0.6), 2, 2),
    fx = function(t, x) c(-0.4 * x[1] + 0.1, -0.25 * x[2]),
    Ly = function(t, x) matrix(c(-0.9, 0.2, 0.1, -0.7), 2, 2),
    fy = function(t, x) c(0.2, -0.1),
    Sx1 = function(t, x) diag(c(0.6, 0.5), 2),
    Sy2 = function(t, x) diag(c(0.7, 0.4), 2),
    k = 2, l = 2, name = "trusted-matrix")
  obs <- as_obs(simulate(model, seed = seed, T = T, dt = dt, burn_in = 0))
  init <- list(mean = c(0, 0), cov = diag(0.1, 2))
  list(model = model, obs = obs, init = init,
       filter = aci_filter(model, obs, init = init))
}


expect_same_smoother <- function(a, b) {
  expect_identical(a$t, b$t)
  expect_identical(a$mean, b$mean)
  expect_identical(a$cov, b$cov)
  expect_identical(a$meta, b$meta)
}


test_that("the authenticated route returns the validated smoother to the bit", {
  for (ds in list(trusted_scalar_setup(), trusted_matrix_setup())) {
    trusted <- aci_smoother(ds$model, ds$obs, filter = ds$filter, init = ds$init)
    forced <- aci_smoother(ds$model, ds$obs, filter = ds$filter, init = ds$init,
                           force_validate = TRUE)
    recomputed <- aci_smoother(ds$model, ds$obs, init = ds$init)

    expect_true(.da_filter_authenticated(
      ds$filter, .compile_cgns_run(ds$model, ds$obs), ds$model))
    expect_same_smoother(trusted, forced)
    expect_same_smoother(trusted, recomputed)
  }
})


test_that("the authenticated route survives the substepped and implicit filters", {
  ds <- trusted_scalar_setup()
  for (args in list(list(stepper = "explicit", nsub = 3L),
                    list(stepper = "implicit", nsub = 1L),
                    list(stepper = "implicit", nsub = 2L))) {
    f <- do.call(aci_filter, c(list(ds$model, ds$obs, init = ds$init), args))
    trusted <- do.call(aci_smoother,
                       c(list(ds$model, ds$obs, filter = f, init = ds$init), args))
    forced <- do.call(aci_smoother,
                      c(list(ds$model, ds$obs, filter = f, init = ds$init,
                             force_validate = TRUE), args))
    expect_same_smoother(trusted, forced)
    expect_identical(trusted$meta$nsub, args$nsub)
  }
})


test_that("the authenticated route carries a non-target reduction", {
  model <- aci_enso_model(hidden = "hW", variant = "aci_code")
  observed <- model$meta$vars$observed
  sim <- simulate(model, seed = 31, T = 0.4, dt = 0.002, burn_in = 0)
  x <- as_obs(sim)$x
  colnames(x) <- observed
  obs <- observed_trajectory(as_obs(sim)$t, x)
  init <- list(mean = 0, cov = matrix(0.1, 1, 1))
  nt <- aci_conditional(given = setdiff(observed, "TC"), method = "mask")

  f <- aci_filter(model, obs, conditional = nt, init = init)
  trusted <- aci_smoother(model, obs, filter = f, conditional = nt, init = init)
  forced <- aci_smoother(model, obs, filter = f, conditional = nt, init = init,
                         force_validate = TRUE)
  expect_same_smoother(trusted, forced)

  ## the same filter offered for the unconditioned run must not authenticate,
  ## and must still raise the unchanged non-target condition
  expect_false(.da_filter_authenticated(
    f, .compile_cgns_run(model, obs), model))
  expect_error(aci_smoother(model, obs, filter = f, init = init),
               class = "aci_error_nontarget")
})


test_that("a hand-mutated covariance loses authentication and is validated", {
  ds <- trusted_scalar_setup()
  bundle <- .compile_cgns_run(ds$model, ds$obs)

  ## a value change that the SPD check catches
  broken <- ds$filter
  broken$cov[1, 1, 5] <- -1
  expect_false(.da_filter_authenticated(broken, bundle, ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = broken, init = ds$init),
               class = "aci_error_spd")

  ## a value change that leaves the path valid: the validated route accepts it,
  ## and the authenticated route must not be reachable to short-circuit that
  nudged <- ds$filter
  nudged$cov[1, 1, 5] <- nudged$cov[1, 1, 5] * 1.5
  expect_false(.da_filter_authenticated(nudged, bundle, ds$model))
  expect_s3_class(
    aci_smoother(ds$model, ds$obs, filter = nudged, init = ds$init),
    "da_path_gaussian")

  ## a non-finite mean, which only the validation reports
  nan_mean <- ds$filter
  nan_mean$mean[3, 1] <- NaN
  expect_false(.da_filter_authenticated(nan_mean, bundle, ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = nan_mean, init = ds$init),
               class = "aci_error_dims")
})


test_that("the invariants the token stands in for hold on every sealed path", {
  ## Skipping .strict_chol() is licensed by the filter kernels, not by the
  ## token: every covariance they emit is finite, exactly symmetric and
  ## Cholesky-decomposable, and a path whose moments are not finite is left
  ## unsealed so that it goes back through the full validation.
  for (ds in list(trusted_scalar_setup(), trusted_matrix_setup())) {
    f <- ds$filter
    expect_false(is.null(attr(f, ".aci_trusted", exact = TRUE)))
    expect_true(all(is.finite(f$mean)))
    expect_true(all(is.finite(f$cov)))
    l <- dim(f$cov)[1L]
    slice <- function(j) matrix(f$cov[, , j], l, l)
    expect_true(all(vapply(seq_len(dim(f$cov)[3L]), function(j) {
      R <- slice(j)
      identical(R, t(R)) &&
        !is.null(tryCatch(chol(R), error = function(e) NULL))
    }, logical(1))))
  }

  ## a path carrying a non-finite moment is not sealed at all
  ds <- trusted_scalar_setup()
  bare <- .cgns_filter_compiled(.compile_cgns_run(ds$model, ds$obs), ds$init)
  expect_false(is.null(attr(.attach_da_trusted(bare), ".aci_trusted",
                            exact = TRUE)))
  bare$mean[2, 1] <- Inf
  expect_null(attr(.attach_da_trusted(bare), ".aci_trusted", exact = TRUE))
  bare$mean[2, 1] <- 0
  bare$cov[1, 1, 2] <- NaN
  expect_null(attr(.attach_da_trusted(bare), ".aci_trusted", exact = TRUE))
})


test_that("a hand-mutated model or provenance loses authentication", {
  ds <- trusted_scalar_setup()
  bundle <- .compile_cgns_run(ds$model, ds$obs)
  other <- aci_dyad_model(params = list(d_x = 0.5, gamma = 0.5, f_x = 0.5,
                                    s_x = 0.5, d_y = 0.5, f_y = 1, s_y = 1))

  swapped <- ds$filter
  swapped$meta$source_model <- other
  expect_false(.da_filter_authenticated(swapped, bundle, ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = swapped, init = ds$init),
               class = "aci_error_model_contract")

  ## the caller's model, not the token's, decides: an untouched filter offered
  ## against a different model must not authenticate
  expect_false(.da_filter_authenticated(
    ds$filter, .compile_cgns_run(other, ds$obs), other))

  moved <- ds$filter
  moved$meta$obs_x[2, 1] <- moved$meta$obs_x[2, 1] + 1
  expect_false(.da_filter_authenticated(moved, bundle, ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = moved, init = ds$init),
               class = "aci_error_dims")

  stripped <- ds$filter
  stripped$meta$obs_x <- NULL
  expect_false(.da_filter_authenticated(stripped, bundle, ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = stripped, init = ds$init),
               class = "aci_error_dims")

  relabelled <- ds$filter
  relabelled$kind <- "smoother"
  expect_false(.da_filter_authenticated(relabelled, bundle, ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = relabelled, init = ds$init),
               class = "aci_error_dims")
})


test_that("a foreign or forged token is not authenticated", {
  ds <- trusted_scalar_setup()
  bundle <- .compile_cgns_run(ds$model, ds$obs)
  token <- attr(ds$filter, ".aci_trusted", exact = TRUE)

  ## no token at all: a filter built by the internal kernels never carries one
  bare <- .cgns_filter_compiled(bundle, ds$init)
  expect_null(attr(bare, ".aci_trusted", exact = TRUE))
  expect_false(.da_filter_authenticated(bare, bundle, ds$model))

  ## the sentinel does not survive serialization, so a token that has been
  ## through saveRDS() is dead even when every value in it is right
  dead <- ds$filter
  attr(dead, ".aci_trusted") <- unserialize(serialize(token, NULL))
  expect_false(.da_filter_authenticated(dead, bundle, ds$model))
  expect_same_smoother(
    aci_smoother(ds$model, ds$obs, filter = dead, init = ds$init),
    aci_smoother(ds$model, ds$obs, filter = ds$filter, init = ds$init))

  ## and it cannot be forged from a lookalike
  forged <- ds$filter
  attr(forged, ".aci_trusted") <- utils::modifyList(
    token, list(sentinel = new.env(parent = emptyenv())))
  expect_false(.da_filter_authenticated(forged, bundle, ds$model))

  wrong_abi <- ds$filter
  attr(wrong_abi, ".aci_trusted") <- utils::modifyList(token, list(abi = 2L))
  expect_false(.da_filter_authenticated(wrong_abi, bundle, ds$model))

  not_a_token <- ds$filter
  attr(not_a_token, ".aci_trusted") <- "trust me"
  expect_false(.da_filter_authenticated(not_a_token, bundle, ds$model))

  ## a token lifted onto a different path: the token's own components no longer
  ## match the path carrying it
  other_ds <- trusted_scalar_setup(seed = 77)
  transplanted <- other_ds$filter
  attr(transplanted, ".aci_trusted") <- token
  expect_false(.da_filter_authenticated(
    transplanted, .compile_cgns_run(other_ds$model, other_ds$obs),
    other_ds$model))
})


test_that("a filter built for a different grid is rejected as before", {
  ds <- trusted_scalar_setup()
  longer <- trusted_scalar_setup(T = 0.5)
  expect_false(.da_filter_authenticated(
    longer$filter, .compile_cgns_run(ds$model, ds$obs), ds$model))
  expect_error(aci_smoother(ds$model, ds$obs, filter = longer$filter,
                            init = ds$init),
               class = "aci_error_dims")
  expect_error(aci_smoother(ds$model, ds$obs, filter = list(a = 1),
                            init = ds$init),
               class = "aci_error_dims")
})


test_that("force_validate is checked and does not move the result", {
  ds <- trusted_scalar_setup()
  expect_error(aci_smoother(ds$model, ds$obs, filter = ds$filter, init = ds$init,
                            force_validate = NA),
               class = "aci_error_dims")
  expect_error(aci_smoother(ds$model, ds$obs, filter = ds$filter, init = ds$init,
                            force_validate = c(TRUE, TRUE)),
               class = "aci_error_dims")
  expect_error(aci_smoother(ds$model, ds$obs, filter = ds$filter, init = ds$init,
                            force_validate = "yes"),
               class = "aci_error_dims")
  ## force_validate = TRUE must still reject a broken filter
  broken <- ds$filter
  broken$cov[1, 1, 5] <- -1
  expect_error(aci_smoother(ds$model, ds$obs, filter = broken, init = ds$init,
                            force_validate = TRUE),
               class = "aci_error_spd")
})


test_that("the token leaves the filter path's public contract unchanged", {
  ds <- trusted_scalar_setup()
  f <- ds$filter

  ## the token is pure data apart from the shared sentinel, so two filters
  ## computed from the same inputs are still identical()
  expect_identical(f, aci_filter(ds$model, ds$obs, init = ds$init))
  expect_identical(names(f), c("t", "mean", "cov", "kind", "meta"))
  expect_identical(class(f), c("da_path_gaussian", "da_path"))
  expect_identical(sort(names(attr(f, ".aci_trusted", exact = TRUE))),
                   sort(c("sentinel", "abi", "mean", "cov")))

  ## it carries only what the run cannot re-establish live, and holds the
  ## path's own arrays rather than copies of them; object.size() double-counts
  ## those references, but nothing is allocated twice
  token <- attr(f, ".aci_trusted", exact = TRUE)
  expect_identical(token$mean, f$mean)
  expect_identical(token$cov, f$cov)

  ## smoother paths are never fed back in as a filter and carry no token
  expect_null(attr(aci_smoother(ds$model, ds$obs, init = ds$init),
                   ".aci_trusted", exact = TRUE))

  ## nor do the filter paths aci() keeps, which never re-enter aci_smoother()
  expect_null(attr(aci(ds$model, ds$obs, init = ds$init,
                       keep = "paths")$paths$filter,
                   ".aci_trusted", exact = TRUE))
})


## Ledger C3c: .validate_gaussian_path()'s l == 1 positivity test.
## The claim is an equivalence, not an approximation: for a 1x1 slice the
## vector test accepts exactly what .strict_chol() accepts, and rejects exactly
## what it rejects. Both directions are checked against .strict_chol() itself.

test_that("for a 1x1 matrix chol() succeeds exactly when the value is finite and positive", {
  accepts <- function(v) {
    ch <- tryCatch(.strict_chol(matrix(v, 1L, 1L), "probe"),
                   error = function(e) NULL)
    !is.null(ch)
  }
  predicate <- function(v) isTRUE(is.finite(v) & v > 0)
  values <- c(1, 0.5, 1e-300, 1e300, .Machine$double.xmin, 2, 1e-12,
              0, -0, -1e-300, -1, -1e300, NA_real_, NaN, Inf, -Inf)
  for (v in values)
    expect_identical(accepts(v), predicate(v),
                     info = paste("value", format(v)))
  ## and the two directions stated separately, so a change to either shows up
  ## as the direction it broke
  for (v in c(1, 0.5, 1e-300, 1e300, 2, 1e-12))
    expect_true(accepts(v), info = paste("accepts", format(v)))
  for (v in c(0, -0, -1e-300, -1, -1e300, NA_real_, NaN, Inf, -Inf))
    expect_false(accepts(v), info = paste("rejects", format(v)))
  ## a random sweep, so the equivalence is not a property of the listed values
  set.seed(31L)
  for (v in c(stats::rnorm(200L), stats::rnorm(200L, sd = 1e-8),
              stats::rnorm(50L, sd = 1e8)))
    expect_identical(accepts(v), predicate(v), info = paste("sweep", format(v)))
})


test_that("the l == 1 validator rejects the same paths, at the same index, with the same error", {
  ds <- trusted_scalar_setup()
  n <- length(ds$obs$t)
  ## a path this namespace did not produce, so force_validate is not needed to
  ## reach the loop: strip the token and let the full route run
  strip <- function(f) { attr(f, ".aci_trusted") <- NULL; f }

  ## the accepting case
  expect_silent(.validate_gaussian_path(strip(ds$filter), ds$obs, 1L, "filter"))

  ## every rejecting case the loop is responsible for
  for (bad in c(0, -1, -1e-300)) {
    for (at in c(1L, 7L, n)) {
      f <- strip(ds$filter)
      f$cov[1L, 1L, at] <- bad
      err <- tryCatch(.validate_gaussian_path(f, ds$obs, 1L, "filter"),
                      error = function(e) e)
      expect_s3_class(err, "aci_error_spd")
      expect_identical(conditionMessage(err), sprintf(
        "Matrix (filter covariance at index %d) must be positive definite.", at))
    }
  }
  ## the first bad index is the one reported, not the last
  f <- strip(ds$filter)
  f$cov[1L, 1L, 5L] <- -1
  f$cov[1L, 1L, 9L] <- -1
  expect_error(.validate_gaussian_path(f, ds$obs, 1L, "filter"),
               "index 5", class = "aci_error_spd")

  ## a non-finite covariance is still caught earlier, by the dimension check,
  ## with its own class; the positivity branch is not what reports it
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    f <- strip(ds$filter)
    f$cov[1L, 1L, 4L] <- bad
    expect_error(.validate_gaussian_path(f, ds$obs, 1L, "filter"),
                 class = "aci_error_dims")
  }

  ## and l > 1 still walks .strict_chol() per slice
  dm <- trusted_matrix_setup()
  fm <- strip(dm$filter)
  expect_silent(.validate_gaussian_path(fm, dm$obs, 2L, "filter"))
  fm$cov[, , 6L] <- matrix(c(1, 2, 2, 1), 2L, 2L)
  expect_error(.validate_gaussian_path(fm, dm$obs, 2L, "filter"),
               class = "aci_error_spd")
})
