## Ledger C2c: the conditional and tau-partition API surface.
##
## Three things are graded here.  Target-by-name must be exactly the
## complement of block-by-name and change no number.  The MATLAB first-slice
## convention must be inert where the mask is inert and must otherwise
## reproduce a measured, non-round-off difference.  And the tau partition's
## two observation sets must be two named estimands, each reproducible on its
## own terms.

## --------------------------------------------------------------------------
## The ENSO path used by the measured constants below.  It is the same
## source-derived realisation the C-series comparison reports were computed on
## (`aci_enso_model(hidden = c("u","hW","tau"))`, seed 42, t_end = 20, dt = 0.005,
## 4001 points), regenerated here rather than shipped, and it agrees with the
## stored copy of that path to 5e-15.  Built once for the whole file.
enso_c2c_cache <- new.env(parent = emptyenv())

enso_c2c_path <- function() {
  if (!is.null(enso_c2c_cache$path)) return(enso_c2c_cache$path)
  sim <- simulate(aci_enso_model(hidden = c("u", "hW", "tau"), variant = "aci_code"),
                  seed = 42, t_end = 20, dt = 0.005)
  x <- as.matrix(sim$obs$x)
  y <- as.matrix(sim$hidden)
  colnames(x) <- c("TC", "TE", "I")
  colnames(y) <- c("u", "hW", "tau")
  enso_c2c_cache$path <- cbind(t = sim$obs$t, x, y)
  enso_c2c_cache$path
}

enso_c2c_case <- function(hidden, ...) {
  p <- enso_c2c_path()
  m <- aci_enso_model(hidden = hidden, variant = "aci_code", ...)
  list(model = m,
       obs = observed_trajectory(p[, "t"],
                                 p[, m$meta$vars$observed, drop = FALSE]),
       init = list(mean = p[1L, hidden], cov = matrix(0.1, 1, 1)))
}

## Every quantity an ACI run reports, as plain numerics.
c2c_all <- function(a) list(filter_mean = as.numeric(a$paths$filter$mean),
                            filter_cov  = as.numeric(a$paths$filter$cov),
                            smooth_mean = as.numeric(a$paths$smoother$mean),
                            smooth_cov  = as.numeric(a$paths$smoother$cov),
                            aci = as.numeric(a$aci),
                            signal = as.numeric(a$signal),
                            dispersion = as.numeric(a$dispersion))

expect_c2c_identical <- function(a, b, info = NULL) {
  for (nm in names(a))
    expect_identical(a[[nm]], b[[nm]], info = paste(info, nm))
}


################################################################################
# 1. target = names the target block; given = names its complement
################################################################################

test_that("target and given name the two sides of the same split", {
  obs <- observed_trajectory(0.01 * 0:5,
                             cbind(a = sin(0.01 * 0:5), b = cos(0.01 * 0:5),
                                   c = 0.01 * 0:5))
  for (pair in list(list(t = "a",             b = c("b", "c")),
                    list(t = c("a", "c"),     b = "b"),
                    list(t = 2L,              b = c(1L, 3L)),
                    list(t = c(3L, 1L),       b = 2L))) {
    by_target <- .nt_indices(aci_conditional(target = pair$t, method = "mask"), obs)
    by_given <- .nt_indices(aci_conditional(given = pair$b, method = "mask"), obs)
    expect_identical(by_target, by_given)
  }
  ## indices need no column names; names do
  bare <- observed_trajectory(0.01 * 0:5, unname(obs$x))
  expect_identical(.nt_indices(aci_conditional(target = 2L, method = "mask"), bare),
                   list(A = 2L, B = c(1L, 3L)))
  expect_error(.nt_indices(aci_conditional(target = "a", method = "mask"), bare),
               class = "aci_error_nontarget")
})


test_that("target = changes no assimilated number", {
  ds <- enso_c2c_case("hW")
  observed <- ds$model$meta$vars$observed
  for (method in c("mask", "reduce")) {
    by_t <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                conditional = aci_conditional(target = "TC", method = method))
    by_b <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                conditional = aci_conditional(given = setdiff(observed, "TC"),
                                              method = method))
    expect_c2c_identical(c2c_all(by_t), c2c_all(by_b), method)
    expect_identical(by_t$paths$filter$meta$likelihood_idx,
                     by_b$paths$filter$meta$likelihood_idx)
  }
})


test_that("a split needs exactly one side, and both sides must be non-empty", {
  expect_error(aci_conditional(), class = "aci_error_nontarget")
  expect_error(aci_conditional(method = "mask"), class = "aci_error_nontarget")
  expect_error(aci_conditional(given = 1, target = 2), class = "aci_error_nontarget")
  obs <- observed_trajectory(0.01 * 0:5, cbind(sin(0.01 * 0:5), cos(0.01 * 0:5)))
  ## target covering every channel leaves nothing to condition on
  expect_error(.nt_indices(aci_conditional(target = c(1L, 2L), method = "mask"), obs),
               class = "aci_error_nontarget")
  ## and given covering every channel leaves no target, as before
  expect_error(.nt_indices(aci_conditional(given = c(1L, 2L), method = "mask"), obs),
               class = "aci_error_nontarget")
  for (bad in list(numeric(0), 1.5, c(1, 1), 3, NA_real_, Inf))
    expect_error(.nt_indices(aci_conditional(target = bad, method = "mask"), obs),
                 class = "aci_error_nontarget")
})


test_that("printing names whichever side was given", {
  expect_output(print(aci_conditional(2, "mask")), "x_B = \\{2\\}, method = mask")
  expect_output(print(aci_conditional(target = 1, method = "mask")),
                "x_A = \\{1\\}, method = mask")
  expect_output(print(aci_conditional(target = 1, method = "mask",
                                      first_step = "matlab")),
                "first_step = matlab")
  ## the default convention stays out of the line it did not appear on before
  expect_false(grepl("first_step",
                     paste(capture.output(print(aci_conditional(2, "mask"))),
                           collapse = "")))
})


################################################################################
# 2. the first-slice convention
################################################################################

test_that("the masked precision path differs from uniform only in its first slice", {
  k <- 4L
  N <- 9L
  gxx <- array(0, c(k, k, N + 1L))
  for (j in seq_len(N + 1L)) {
    A <- matrix(0.05 * (seq_len(k * k) + j), k, k)
    gxx[, , j] <- diag(k) * (1 + 0.1 * j) + A %*% t(A)
  }
  target <- c(1L, 3L)
  uniform <- .compiled_precision_path(gxx, N, target = target)
  matlab <- .compiled_precision_path(gxx, N, target = target,
                                     first_step = "matlab")
  ## slice 1 is the full inverse; the target block coincides with the mask's
  ## only where the Gram is block-diagonal, which this one is not.  The full
  ## inverse is taken by chol2inv() on a dense Gram here, so it agrees with a
  ## per-slice solve to the last bit or two and is pinned exactly against its
  ## own route.
  expect_equal(matrix(matlab[, , 1L], k, k),
               chol_solve(matrix(gxx[, , 1L], k, k), diag(k), "gxx"),
               tolerance = 1e-12)
  expect_identical(matrix(matlab[, , 1L], k, k),
                   chol2inv(safe_chol(matrix(gxx[, , 1L], k, k), "gxx")) + 0)
  expect_false(identical(matrix(matlab[, , 1L], k, k),
                         matrix(uniform[, , 1L], k, k)))
  for (j in 2:N)
    expect_identical(matrix(matlab[, , j], k, k), matrix(uniform[, , j], k, k))
  ## the default is uniform, and the unconditioned branch has no mask to move
  expect_identical(.compiled_precision_path(gxx, N, target = target,
                                            first_step = "uniform"), uniform)
  expect_identical(.compiled_precision_path(gxx, N), .compiled_precision_path(gxx, N))
  ## a single-interval record still gets its one slice
  expect_equal(
    matrix(.compiled_precision_path(gxx, 1L, target = target,
                                    first_step = "matlab")[, , 1L], k, k),
    chol_solve(matrix(gxx[, , 1L], k, k), diag(k), "gxx"),
    tolerance = 1e-12)
  expect_identical(
    matrix(.compiled_precision_path(gxx, 1L, target = target,
                                    first_step = "matlab")[, , 1L], k, k),
    chol2inv(safe_chol(matrix(gxx[, , 1L], k, k), "gxx")) + 0)
})


test_that("first_step = 'matlab' belongs to the masked method only", {
  expect_error(aci_conditional(2, "reduce", first_step = "matlab"),
               class = "aci_error_nontarget")
  expect_error(aci_conditional(target = 1, method = "reduce",
                               first_step = "matlab"),
               class = "aci_error_nontarget")
  expect_no_error(aci_conditional(2, "reduce", first_step = "uniform"))
  expect_error(aci_conditional(2, "mask", first_step = "nope"))
})


test_that("first_step = 'uniform' is the shipped behaviour to the bit", {
  ds <- enso_c2c_case("hW")
  spec <- aci_conditional(target = "TC", method = "mask")
  explicit <- aci_conditional(target = "TC", method = "mask",
                              first_step = "uniform")
  expect_identical(spec$first_step, "uniform")
  a1 <- aci(ds$model, ds$obs, init = ds$init, conditional = spec,
            keep = "paths", decompose = TRUE)
  a2 <- aci(ds$model, ds$obs, init = ds$init, conditional = explicit,
            keep = "paths", decompose = TRUE)
  expect_c2c_identical(c2c_all(a1), c2c_all(a2))
})


test_that("first_step = 'matlab' reproduces the measured hW transient", {
  ## C4 section 2.3, measured independently through aciR with a supplied
  ## S_xoS_x_inv, on this same 4001-point path.  These are that measurement to
  ## the digits it reported: signed -1.0798e-01 at step 2, peak 1.5619e-01 at
  ## step 4, filter covariance 3.6563e-02, ACI 5.7434e-01.
  ds <- enso_c2c_case("hW")
  uniform <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                 conditional = aci_conditional(target = "TC", method = "mask"))
  matlab <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                conditional = aci_conditional(target = "TC", method = "mask",
                                              first_step = "matlab"))
  d <- as.numeric(matlab$paths$filter$mean) - as.numeric(uniform$paths$filter$mean)
  expect_equal(d[1L], 0)                       # the prior is untouched
  expect_equal(d[2L], -1.079777623214551e-01, tolerance = 1e-9)
  expect_identical(which.max(abs(d)), 4L)
  expect_equal(max(abs(d)), 1.561914587488461e-01, tolerance = 1e-9)
  expect_equal(max(abs(matlab$paths$filter$cov - uniform$paths$filter$cov)),
               3.656250000000001e-02, tolerance = 1e-9)
  expect_equal(max(abs(matlab$aci - uniform$aci)),
               5.743363354542981e-01, tolerance = 1e-9)
  ## and it is a transient on a live scale, not a round-off: the step-2 filter
  ## mean itself is 0.1438, so this is a 75% move
  expect_gt(abs(d[2L]) / abs(as.numeric(uniform$paths$filter$mean)[2L]), 0.7)
})


test_that("first_step is exactly inert where the mask itself is inert", {
  ## The u-hidden conditional case has one non-zero Lx row, the target's, so
  ## zeroing the rest of the Gram inverse cannot move the gain - with or
  ## without the first slice (C4 section 2.3, "Same experiment, u hidden").
  ds <- enso_c2c_case("u")
  uniform <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                 conditional = aci_conditional(target = "TC", method = "mask"))
  matlab <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                conditional = aci_conditional(target = "TC", method = "mask",
                                              first_step = "matlab"))
  expect_c2c_identical(c2c_all(uniform), c2c_all(matlab))
})


################################################################################
# 3. the tau partition's two observation sets
################################################################################

test_that("hidden = 'tau' defaults to the reduced observation set", {
  reduced <- aci_enso_model(hidden = "tau", variant = "aci_code")
  full <- aci_enso_model(hidden = "tau", variant = "aci_code", observations = "full")
  expect_identical(reduced$meta$observations, "reduced")
  expect_identical(full$meta$observations, "full")
  expect_s3_class(reduced$meta$estimand_nontarget, "aci_conditional_spec")
  expect_identical(reduced$meta$estimand_nontarget$given, c(1L, 2L))
  expect_identical(reduced$meta$estimand_nontarget$method, "reduce")
  expect_null(full$meta$estimand_nontarget)
  ## the model itself still observes all five channels
  expect_identical(reduced$k, 5L)
  expect_identical(reduced$meta$vars$observed, c("u", "hW", "TC", "TE", "I"))
  ## every other partition is full, and cannot ask for the reduction
  for (h in list("u", "hW", c("u", "hW"), c("u", "hW", "tau"))) {
    expect_identical(aci_enso_model(hidden = h, variant = "aci_code")$meta$observations,
                     "full")
    expect_error(aci_enso_model(hidden = h, variant = "aci_code",
                                observations = "reduced"),
                 class = "aci_error_model_contract")
  }
})


test_that("the reduced tau estimand is the prescribed-forcing reduction", {
  ds <- enso_c2c_case("tau")                    # observations = "reduced"
  full <- aci_enso_model(hidden = "tau", variant = "aci_code", observations = "full")
  declared <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE)
  explicit <- aci(full, ds$obs, init = ds$init, keep = "paths", decompose = TRUE,
                  conditional = aci_conditional(given = c("u", "hW"),
                                                method = "reduce"))
  expect_c2c_identical(c2c_all(declared), c2c_all(explicit))
  ## and it reaches the reduced three-channel system
  f <- aci_filter(ds$model, ds$obs, init = ds$init)
  expect_identical(f$meta$likelihood_idx, 1:3)
  bundle <- .compile_cgns_run(ds$model, ds$obs)
  expect_identical(bundle$k, 3L)
  expect_identical(bundle$model$meta$nontarget_reduction,
                   list(target = 3:5, prescribed = 1:2))
  ## the reduced model must not carry the declaration onward
  expect_null(bundle$model$meta$estimand_nontarget)
  ## a declared estimand needs no column names
  bare <- observed_trajectory(ds$obs$t, unname(ds$obs$x))
  expect_identical(aci_filter(ds$model, bare, init = ds$init)$mean, f$mean)
})


test_that("the two tau estimands differ by the amount C3 measured", {
  ## C3 section 5: the reference script asserts the reduced and full
  ## observation sets agree.  On this path they do not.  Filter mean 2.470e-01,
  ## filter covariance 8.492e-03, smoother mean 1.906e-01, ACI 7.765e-01,
  ## Pearson 0.905, with the time-averaged ACI within 0.5%.
  ds <- enso_c2c_case("tau")
  full <- aci_enso_model(hidden = "tau", variant = "aci_code", observations = "full")
  r <- aci(ds$model, ds$obs, init = ds$init, keep = "paths", decompose = TRUE)
  f <- aci(full, ds$obs, init = ds$init, keep = "paths", decompose = TRUE)
  expect_equal(max(abs(r$paths$filter$mean - f$paths$filter$mean)),
               2.470215324439298e-01, tolerance = 1e-9)
  expect_equal(max(abs(r$paths$filter$cov - f$paths$filter$cov)),
               8.491541936276825e-03, tolerance = 1e-9)
  expect_equal(max(abs(r$paths$smoother$mean - f$paths$smoother$mean)),
               1.905819377600842e-01, tolerance = 1e-9)
  expect_equal(max(abs(r$aci - f$aci)), 7.764861042153446e-01, tolerance = 1e-9)
  expect_equal(cor(as.numeric(r$aci), as.numeric(f$aci)), 0.904551,
               tolerance = 1e-5)
  ## the average level is nearly preserved while the curve is not - the mean
  ## ACI moves 0.51%, the pointwise curve by three times its own mean level.
  ## That is why a fidelity claim has to name its observation set.
  expect_equal((mean(r$aci) - mean(f$aci)) / mean(f$aci), -0.0050747380,
               tolerance = 1e-6)
  expect_gt(max(abs(r$aci - f$aci)) / mean(f$aci), 3)
})


test_that("a declared estimand refuses to compose with a second one", {
  ds <- enso_c2c_case("tau")
  for (spec in list(aci_conditional(target = "TC", method = "mask"),
                    aci_conditional(given = c("u", "hW"),
                                    method = "reduce"))) {
    expect_error(aci_filter(ds$model, ds$obs, init = ds$init,
                            conditional = spec),
                 class = "aci_error_nontarget")
    expect_error(aci(ds$model, ds$obs, init = ds$init, conditional = spec),
                 class = "aci_error_nontarget")
  }
  ## a malformed declaration is caught rather than silently ignored
  broken <- ds$model
  broken$meta$estimand_nontarget <- list(given = 1L)
  expect_error(aci_filter(broken, ds$obs, init = ds$init),
               class = "aci_error_nontarget")
})


################################################################################
# 5. The recorded conditional estimand is the reference script's (W3e)
################################################################################

test_that("the recorded conditioning set is the script's, not the inert one", {
  ## Found while drafting the vignettes (O6 finding D1).  aci_enso_model() used to
  ## record the set of observed channels the hidden variable REACHES as the
  ## target - (TC, TE, I) - which made the complement {u, tau}.  Masking that
  ## complement is exactly inert, because the u and tau rows of Lx are
  ## identically zero, so a caller who built a mask from the field got the
  ## UNCONDITIONAL run under a conditional name.  The reference scripts mask
  ## down to T_C alone (ENSO_model_cond_ACI_h_W_unobs.m:1202,
  ## ENSO_model_cond_ACI_u_unobs.m:1205; the T_E and I lines are commented
  ## out), and that complement is emphatically not inert.
  ds <- enso_c2c_case("hW")
  m <- ds$model
  observed <- m$meta$vars$observed
  expect_identical(observed, c("u", "TC", "TE", "tau", "I"))
  expect_identical(observed[m$meta$target_obs_idx], "TC")
  expect_identical(observed[m$meta$conditioning_obs_idx],
                   c("u", "TE", "tau", "I"))
  expect_identical(m$meta$causal_link, "(hW) -> (TC) | (u,TE,tau,I)")
  expect_match(m$meta$estimand_provenance,
               "ENSO_model_cond_ACI_h_W_unobs.m:1202", fixed = TRUE)
  ## The field is a record, not a default: nothing is conditioned unless the
  ## caller says so, and this partition declares no estimand of its own.
  expect_null(m$meta$estimand_nontarget)
  expect_match(m$meta$conditioning_note, "a record, not a default",
               fixed = TRUE)

  ## The two complements, measured, so the fix is a fact and not an assertion
  ## about intent.  Only the TC and TE rows of Lx are non-zero.
  lx <- vapply(seq_along(ds$obs$t),
               function(j) abs(drop(eval_coefs(m, ds$obs$t[j],
                                               ds$obs$x[j, ])$Lx)),
               numeric(m$k))
  rownames(lx) <- observed
  reach <- apply(lx, 1L, max)
  expect_identical(unname(reach[c("u", "tau")]), c(0, 0))
  expect_gt(min(reach[c("TC", "TE")]), 0.4)

  base <- aci(m, ds$obs, init = ds$init, decompose = TRUE, keep = "paths")
  inert <- aci(m, ds$obs, init = ds$init, decompose = TRUE, keep = "paths",
               conditional = aci_conditional(given = c("u", "tau"),
                                             method = "mask"))
  script <- aci(m, ds$obs, init = ds$init, decompose = TRUE, keep = "paths",
                conditional = aci_conditional(target = "TC", method = "mask"))
  ## The old recorded complement: identical, not merely close.
  expect_c2c_identical(c2c_all(base), c2c_all(inert), "old recorded complement")
  expect_identical(max(abs(base$aci - inert$aci)), 0)
  ## The script's complement: a difference of the size of the answer itself.
  expect_false(identical(as.numeric(base$aci), as.numeric(script$aci)))
  expect_gt(max(abs(base$aci - script$aci)) / max(base$aci), 0.2)
  ## And the complement read straight off the field reproduces the script arm.
  from_meta <- aci(m, ds$obs, init = ds$init, decompose = TRUE, keep = "paths",
                   conditional = aci_conditional(given = observed[
                       m$meta$conditioning_obs_idx], method = "mask"))
  expect_c2c_identical(c2c_all(script), c2c_all(from_meta), "field-built mask")
})

test_that("the unconditional partitions record an empty conditioning set", {
  ## The tau and joint scripts leave every masking line commented
  ## (tau_unobs.m:1205-1211, u_h_W_tau_unobs.m:1218-1226), so their runs are
  ## unconditional and their masked complement is empty.  For the reduced tau
  ## estimand the two prescribed channels are still conditioned upon, by
  ## prescription rather than by masking, and the field says so.
  joint <- aci_enso_model(hidden = c("u", "hW", "tau"), variant = "aci_code")
  expect_identical(joint$meta$vars$observed[joint$meta$target_obs_idx],
                   c("TC", "TE", "I"))
  expect_length(joint$meta$conditioning_obs_idx, 0L)
  expect_identical(joint$meta$causal_link, "(u,hW,tau) -> (TC,TE,I)")

  full <- aci_enso_model(hidden = "tau", variant = "aci_code",
                         observations = "full")
  expect_identical(full$meta$vars$observed[full$meta$target_obs_idx],
                   c("u", "hW", "TC", "TE", "I"))
  expect_length(full$meta$conditioning_obs_idx, 0L)

  red <- aci_enso_model(hidden = "tau", variant = "aci_code",
                        observations = "reduced")
  expect_identical(red$meta$vars$observed[red$meta$target_obs_idx],
                   c("TC", "TE", "I"))
  expect_identical(red$meta$vars$observed[red$meta$conditioning_obs_idx],
                   c("u", "hW"))
  expect_identical(red$meta$causal_link, "(tau) -> (TC,TE,I) | (u,hW)")
  ## and those two are exactly the channels the declared estimand prescribes
  expect_identical(red$meta$estimand_nontarget$given,
                   red$meta$conditioning_obs_idx)
})
