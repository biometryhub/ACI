## Ledger entry C2d.  aci_enso_model() declares a whole-path realisation of its own
## coefficient expressions; every assertion here is that taking it changes
## nothing but the time it takes.  `identical()` throughout, never a tolerance.

enso6_partitions <- function() list(
  joint = c("u", "hW", "tau"), u = "u", hW = "hW",
  tau_full = "tau", te = "TE", u_hW = c("u", "hW"),
  all4 = c("u", "hW", "TE", "tau")
)

enso6_model_for <- function(hidden) {
  ## hidden = "tau" defaults to the reduced observation set, which declares an
  ## estimand and so refuses composition with a caller's conditional spec (C2c).
  ## The conditioning cases below need the five-channel build.
  if (identical(hidden, "tau"))
    aci_enso_model(hidden = "tau", variant = "aci_code", observations = "full")
  else aci_enso_model(hidden = hidden, variant = "aci_code")
}

enso6_path <- function(model, n = 61L, from = 0) {
  tt <- from + seq(0, by = 1 / 24, length.out = n)
  nm <- model$meta$vars$observed
  x <- matrix(0, n, length(nm))
  base <- c(u = 0.01, hW = -0.02, TC = 0.03, TE = 0.05, tau = -0.03, I = 1.6)
  amp <- c(u = 0.02, hW = 0.03, TC = 0.04, TE = 0.06, tau = 0.05, I = 0.7)
  for (i in seq_along(nm))
    x[, i] <- base[[nm[i]]] + amp[[nm[i]]] * sin(1.3 * tt + i)
  observed_trajectory(tt, x, names = nm)
}


test_that("the ENSO6 batch realiser reproduces the grid realiser on every partition", {
  for (nm in names(enso6_partitions())) {
    model <- enso6_model_for(enso6_partitions()[[nm]])
    for (n in c(2L, 3L, 61L)) {
      obs <- as_obs(enso6_path(model, n = n))
      batch <- .realise_affine_cgns_grid_once(model, obs, batch = TRUE)
      grid <- .realise_affine_cgns_grid_once(model, obs, batch = FALSE)
      generic <- .realise_cgns_grid_once(model, obs)
      lab <- sprintf("%s n=%d", nm, n)
      for (field in names(grid)) {
        expect_identical(batch[[field]], grid[[field]],
                         info = paste(lab, field, "vs grid"))
        expect_identical(batch[[field]], generic[[field]],
                         info = paste(lab, field, "vs generic closures"))
      }
      expect_identical(attr(batch, "path_cross_noise"),
                       attr(grid, "path_cross_noise"), info = lab)
      expect_identical(attr(batch, "path_cross_noise"),
                       attr(generic, "path_cross_noise"), info = lab)
      expect_identical(dim(batch$Lx), c(model$k, model$l, n), info = lab)
      expect_identical(dim(batch$fx), c(n, model$k), info = lab)
      expect_identical(dim(batch$fy), c(n, model$l), info = lab)
      expect_identical(dim(batch$gyx), c(model$l, model$k, n), info = lab)
    }
  }
})


test_that("the ENSO6 batch arrays match the model's own coefficient closures slice by slice", {
  for (nm in c("joint", "hW", "u_hW")) {
    model <- enso6_model_for(enso6_partitions()[[nm]])
    obs <- as_obs(enso6_path(model, n = 61L))
    batch <- .realise_affine_cgns_grid_once(model, obs)
    for (j in c(1L, 2L, 30L, 61L)) {
      co <- eval_coefs(model, obs$t[j], obs$x[j, ])
      lab <- paste(nm, j)
      ## `[ , , j]` drops the unit margins a scalar partition has, so the
      ## slices are compared as the doubles they are; the array dimensions are
      ## asserted on the whole realisation above.
      expect_identical(as.numeric(batch$Lx[, , j]), as.numeric(co$Lx), info = lab)
      expect_identical(batch$fx[j, ], as.numeric(co$fx), info = lab)
      expect_identical(as.numeric(batch$Ly[, , j]), as.numeric(co$Ly), info = lab)
      expect_identical(batch$fy[j, ], as.numeric(co$fy), info = lab)
      expect_identical(as.numeric(batch$gxx[, , j]), as.numeric(co$gxx), info = lab)
      expect_identical(as.numeric(batch$gyy[, , j]), as.numeric(co$gyy), info = lab)
      expect_identical(as.numeric(batch$gyx[, , j]), as.numeric(co$gyx), info = lab)
    }
  }
})


test_that("ENSO6 bundles are unchanged under both conditioning strategies", {
  grid_compile <- function(model, obs, nt = NULL) {
    nt <- .model_estimand_spec(model, nt)
    src <- as_obs(obs)
    .compile_cgns_complete(
      model, src, conditional = nt,
      full = .realise_affine_cgns_grid_once(model, src, batch = FALSE),
      realization = "affine_batch")
  }
  for (nm in names(enso6_partitions())) {
    model <- enso6_model_for(enso6_partitions()[[nm]])
    obs <- enso6_path(model, n = 61L)
    observed <- model$meta$vars$observed
    given <- setdiff(observed, "TC")
    specs <- list(
      unconditioned = NULL,
      mask = aci_conditional(given = given, method = "mask"),
      reduce = aci_conditional(given = given, method = "reduce"),
      matlab = aci_conditional(target = "TC", method = "mask",
                               first_step = "matlab"))
    for (sn in names(specs)) {
      got <- .compile_cgns_run(model, obs, specs[[sn]])
      want <- grid_compile(model, obs, specs[[sn]])
      lab <- paste(nm, sn)
      expect_identical(got$realization, "affine_batch", info = lab)
      expect_identical(names(got$coefficients), names(want$coefficients), info = lab)
      for (field in names(want$coefficients))
        expect_identical(got$coefficients[[field]], want$coefficients[[field]],
                         info = paste(lab, field))
      expect_identical(got$correlated_noise, want$correlated_noise, info = lab)
      expect_identical(got$likelihood_idx, want$likelihood_idx, info = lab)
      expect_identical(got$scalar, want$scalar, info = lab)
    }
  }
})


test_that("the declared tau reduction is realised through the full model and then sliced", {
  full <- aci_enso_model(hidden = "tau", variant = "aci_code", observations = "full")
  reduced <- aci_enso_model(hidden = "tau", variant = "aci_code")
  obs <- enso6_path(full, n = 61L)
  declared <- .compile_cgns_run(reduced, obs)
  explicit <- .compile_cgns_run(
    full, obs, aci_conditional(given = c("u", "hW"), method = "reduce"))
  expect_identical(declared$coefficients, explicit$coefficients)
  expect_identical(declared$likelihood_idx, explicit$likelihood_idx)
  expect_identical(declared$k, 3L)
  ## and the reduced arrays really are slices of the realised five-channel ones
  wide <- .realise_affine_cgns_grid_once(full, as_obs(obs))
  expect_identical(declared$coefficients$Lx,
                   array(wide$Lx[3:5, , , drop = FALSE], c(3L, 1L, 61L)))
  expect_identical(declared$coefficients$gxx,
                   array(wide$gxx[3:5, 3:5, , drop = FALSE], c(3L, 3L, 61L)))
})


test_that("ENSO6 realiser selection uses sealed constructor identity", {
  model <- aci_enso_model(hidden = c("u", "hW", "tau"), variant = "aci_code")
  obs <- enso6_path(model, n = 41L)
  descriptor <- .cgns_realizer_descriptor(model)
  expect_true(is.environment(descriptor))
  expect_identical(descriptor$id, "enso6_aci_code_v1")
  expect_true(is.function(descriptor$spec$f_full))
  expect_true(is.function(descriptor$spec$g_full))
  expect_identical(.compile_cgns_run(model, obs)$realization, "affine_batch")

  ## metadata is not a routing key, and does not reach the arrays
  metadata_only <- model
  metadata_only$name <- "a harmless presentation label"
  metadata_only$meta$params$r <- 999
  expect_identical(.compile_cgns_run(metadata_only, obs)$coefficients,
                   .compile_cgns_run(model, obs)$coefficients)

  ## a replaced coefficient function authenticates as nothing
  changed <- model
  changed$fx <- function(t, x) rep(0, 3)
  expect_null(.cgns_realizer_descriptor(changed))
  expect_identical(.compile_cgns_run(changed, obs)$realization,
                   "generic_closure_one_pass")
  expect_error(.realise_affine_cgns_grid_once(changed, as_obs(obs)),
               class = "aci_error_compiled_contract")

  roundtrip <- unserialize(serialize(model, NULL))
  expect_null(.cgns_realizer_descriptor(roundtrip))
  expect_identical(.compile_cgns_run(roundtrip, obs)$realization,
                   "generic_closure_one_pass")

  ## the captured constructor parameters cannot be mutated behind the
  ## descriptor's back: both routes read one locked environment
  ce <- environment(descriptor$spec$f_full)
  expect_true(environmentIsLocked(ce))
  expect_identical(ce, environment(descriptor$spec$enso6$drift))
  expect_identical(ce, environment(descriptor$spec$enso6$sd))
  expect_identical(ce, environment(model$Sx1))
  expect_error(assign("p", list(), envir = ce))
  expect_error(eval(quote(p$r <- 999), ce))
  expect_error(eval(quote(sigma_E <- 1), ce))
})


test_that("an aci_model_from_affine model with no batch declaration still realises", {
  model <- aci_model_from_affine(
    f_full = function(t, x, y) c(-0.2 * x[1] + 0.7 * y, -0.1 * x[2] - 0.25 * y),
    g_full = function(t, x, y) 0.05 * x[1] - 0.6 * y,
    Sx = function(t, x) diag(c(0.7, 0.6)),
    Sy_hidden = function(t, x) matrix(0.8, 1, 1), k = 2L, l = 1L)
  tt <- seq(0, 0.2, by = 0.01)
  obs <- as_obs(observed_trajectory(tt, cbind(sin(tt), cos(tt))))
  descriptor <- .cgns_realizer_descriptor(model)
  expect_identical(descriptor$id, "affine_model_v1")
  expect_null(descriptor$spec$enso6)
  expect_identical(.realise_affine_cgns_grid_once(model, obs, batch = TRUE),
                   .realise_affine_cgns_grid_once(model, obs, batch = FALSE))
})


test_that("the ENSO6 batch realiser keeps the per-point contract errors", {
  model <- aci_enso_model(hidden = "hW", variant = "aci_code")
  spec <- .cgns_realizer_descriptor(model)$spec$enso6
  obs <- as_obs(enso6_path(model, n = 9L))

  ## a non-finite realised coefficient stops at the same index either way
  nm <- model$meta$vars$observed
  bad_x <- obs$x
  bad_x[4L, match("TC", nm)] <- 1e200
  bad <- as_obs(observed_trajectory(obs$t, bad_x))
  expect_error(.realise_affine_cgns_grid_once(model, bad, batch = TRUE),
               "non-finite values at index 4", class = "aci_error_model_contract")
  expect_error(.realise_affine_cgns_grid_once(model, bad, batch = FALSE),
               "non-finite values at index 4", class = "aci_error_model_contract")

  ## and a declaration whose drift does not return one value per state and
  ## grid point is a model-contract error, not a silent reshape
  short <- spec
  short$drift <- function(t, s) numeric(3)
  expect_error(.realise_enso6_cgns_grid_once(short, obs, model$k, model$l),
               class = "aci_error_model_contract")
  short2 <- spec
  short2$sd <- function(t, s) numeric(3)
  expect_error(.realise_enso6_cgns_grid_once(short2, obs, model$k, model$l),
               class = "aci_error_model_contract")
})


test_that("the ENSO6 drift and diffusion hold elementwise over a whole path", {
  ## This is the property the batch realiser rests on: the constructor's
  ## expressions give the same doubles whether they are handed one grid point
  ## or the whole grid.
  model <- aci_enso_model(hidden = c("u", "hW", "tau"), variant = "aci_code")
  spec <- .cgns_realizer_descriptor(model)$spec$enso6
  obs <- as_obs(enso6_path(model, n = 17L))
  n <- length(obs$t)
  state <- rep(list(numeric(n)), 6L)
  for (i in seq_len(model$k)) state[[spec$obs_i[i]]] <- obs$x[, i]
  whole <- spec$drift(obs$t, state)
  dim(whole) <- c(n, 6L)
  amps <- spec$sd(obs$t, state)
  dim(amps) <- c(n, 6L)
  for (j in seq_len(n)) {
    point <- vapply(state, function(z) z[j], numeric(1))
    expect_identical(whole[j, ], spec$drift(obs$t[j], point), info = j)
    expect_identical(amps[j, ], spec$sd(obs$t[j], point), info = j)
  }
})


test_that("the precision path keeps the per-slice arithmetic on both branches", {
  ## The guarded per-slice form IS the definition.  The shipped path lifts
  ## safe_chol()'s re-entry checks out of the loop, and takes the unmasked
  ## inverse with chol2inv() on the same factor.  These Grams are dense, so the
  ## unmasked route agrees with the definition to the last bit or two and is
  ## pinned exactly against its own route; the masked branch is unchanged and
  ## must not move a bit.
  set.seed(19)
  k <- 4L
  N <- 25L
  gxx <- array(0, c(k, k, N + 1L))
  for (j in seq_len(N + 1L)) {
    A <- matrix(0.05 * (seq_len(k * k) + j) + 0.02 * rnorm(k * k), k, k)
    gxx[, , j] <- diag(k) * (1 + 0.1 * j) + A %*% t(A)
  }
  guarded <- function(gm) chol_solve(gm, diag(k), "gxx")
  route <- function(gm) chol2inv(safe_chol(gm, "gxx")) + 0
  definition <- function(g, N, target = NULL, matlab = FALSE,
                         full_inv = guarded) {
    kk <- k * k
    pr <- numeric(kk * N)
    from <- 1L
    if (!is.null(target) && matlab) {
      gm <- g[, , 1L]; dim(gm) <- c(k, k)
      pr[seq_len(kk)] <- full_inv(gm)
      from <- 2L
    }
    if (from <= N) for (j in from:N) {
      gm <- g[, , j]; dim(gm) <- c(k, k)
      at <- (j - 1L) * kk
      pr[(at + 1L):(at + kk)] <- if (is.null(target))
        full_inv(gm) else masked_ginv(gm, target)
    }
    dim(pr) <- c(k, k, N)
    pr
  }
  expect_equal(.compiled_precision_path(gxx, N), definition(gxx, N),
               tolerance = 1e-12)
  expect_identical(.compiled_precision_path(gxx, N),
                   definition(gxx, N, full_inv = route))
  for (tg in list(1L, c(1L, 3L), c(2L, 4L), seq_len(k))) {
    ## masked slices are untouched by the route change
    expect_identical(.compiled_precision_path(gxx, N, target = tg),
                     definition(gxx, N, tg), info = paste(tg, collapse = ","))
    expect_equal(
      .compiled_precision_path(gxx, N, target = tg, first_step = "matlab"),
      definition(gxx, N, tg, matlab = TRUE), tolerance = 1e-12,
      info = paste(tg, collapse = ","))
    expect_identical(
      .compiled_precision_path(gxx, N, target = tg, first_step = "matlab"),
      definition(gxx, N, tg, matlab = TRUE, full_inv = route),
      info = paste(tg, collapse = ","))
  }

  ## a slice that does not factor sends the whole path through the guarded
  ## pair, jitter ladder and all
  ladder <- gxx
  ladder[, , 7L] <- matrix(1, k, k)
  expect_identical(.compiled_precision_path(ladder, N), definition(ladder, N))
  expect_identical(.compiled_precision_path(ladder, N, target = c(1L, 3L)),
                   definition(ladder, N, c(1L, 3L)))

  ## a non-finite Gram raises the same classed error from the same branch
  nonfinite <- gxx
  nonfinite[1L, 1L, 5L] <- NA_real_
  expect_error(.compiled_precision_path(nonfinite, N), class = "aci_error_spd")
  expect_error(.compiled_precision_path(nonfinite, N, target = c(1L, 3L)),
               class = "aci_error_dims")

  ## an invalid target is rejected once, up front, with masked_ginv()'s message
  for (bad in list(0L, k + 1L, c(1L, 1L), 1.5, -1L, NA_integer_))
    expect_error(.compiled_precision_path(gxx, N, target = bad),
                 class = "aci_error_dims")
  ## an empty target is the zero weighting, as masked_ginv() defines it
  expect_identical(.compiled_precision_path(gxx, N, target = integer(0)),
                   array(0, c(k, k, N)))
})
