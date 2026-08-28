# Always-on source-derived grades for the scalar ENSO partitions ---------------
#
# ACI_code carries five ENSO scripts, not one joint ENSO case.  The joint
# three-hidden case is graded by the `enso` fixture in test-19; these fixtures
# grade the three SCALAR-hidden partitions (u, h_W, tau) and the reduced
# three-dimensional observation the tau script uses as its operational
# shortcut.  `oracle-manifest-partitions.yml` is the authority for their
# provenance and scope.  In particular:
#
# * these are SOURCE-DERIVED, not authors-source.  The driving path is an R
#   simulation of the ACI_code ENSO system, not a MATLAB realisation, and
#   aci_enso_model()'s own metadata records matlab_simulator_parity = FALSE;
# * the pinned values were produced by aci 0.0.30 and independently reproduced
#   by aciR 0.2.3's multivariate recursions on the identical realised
#   coefficient arrays, to 2.11e-15.  Two implementations agreeing is useful
#   evidence and is NOT an authors'-reference grounding;
# * no aciR constructor for these partitions exists, so nothing here grades a
#   second independent coefficient realisation;
# * the upstream u, h_W and tau scripts run conditional ACI with the first
#   filter step's inverse left unmasked.  These fixtures grade the
#   UNCONDITIONAL observation.  A conditional grade needs its own fixture built
#   to match that convention deliberately.
#
# Keep those distinctions in test descriptions and failure messages.
#
# The observation set is named on EVERY construction below.  Since C2c,
# aci_enso_model(hidden = "tau") defaults to the reduced three-channel estimand,
# and the two estimands differ by up to 0.247 in the filter mean: a test that
# leant on the default would compare one arm against the other's pinned bytes
# the next time that default moved.  Two different fixtures grade the two arms
# here, and each names what it is.

.partition_oracle_dir <- testthat::test_path("fixtures", "oracles")
.partition_oracle_manifest <- "oracle-manifest-partitions.yml"

# Gate rationale, recorded in the manifest's `validation:` block and restated
# here because a reader of the test should not have to go and find it.  This is
# a regression gate on a pinned 15-significant-digit CSV, not a claim about
# accuracy.  The file's own quantisation floor is 1.7e-15 (the largest per-step
# quantity pinned is the wind-burst Gram diagonal at 3.383, quantised at
# 5 * 3.383 * 1e-16).  Measured agreement sits at that floor: 5.33e-15,
# dominated by the coefficient columns, while the producer and this package
# agree BIT-FOR-BIT in memory on every coefficient array and every filter and
# smoother moment.  So the number gated is a file round trip.  1e-8 is six
# orders above it and seven below the smallest scientifically meaningful
# divergence, the reduced tau arm's 0.247 departure from the full one.  The
# rest of this directory uses 1e-6; this is two decades tighter because these
# fixtures grade a same-project producer, with no cross-language printing slack
# to absorb.  It is deliberately not tighter than 1e-8: 4001 sequential steps
# under a different BLAS could move the accumulated sums past the floor with no
# code change, and a genuine realiser or recursion regression would move a
# coefficient far further than 1e-8 in any case.
.partition_oracle_tolerance <- 1e-8

.partition_oracle_file <- function(name) {
  path <- file.path(.partition_oracle_dir, name)
  expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("always-on partition fixture %s is missing", name)
  )
  path
}

# Shared with test-19 through helper-oracle-manifest.R: one regex reads both
# manifests, and test-19 asserts that neither manifest has an entry the reader
# cannot see.
.partition_manifest_hashes <- function(algorithm = "md5") {
  .oracle_manifest_hashes(
    .partition_oracle_file(.partition_oracle_manifest), algorithm
  )
}

.partition_signal <- function() {
  read.csv(.partition_oracle_file("enso6_partition_signal.csv"))
}

# One partition rebuilt from the pinned signal through the PUBLIC entry points.
# The compiled bundle is taken as well, since the realised coefficient arrays
# have no public accessor and are half of what these fixtures grade.
#
# `observations` is mandatory, not defaulted: see the file header.  A model
# that declares its own observation set refuses a second specification, so
# `spec` is only ever supplied alongside observations = "full".
.partition_run <- function(hidden, observations, signal, spec = NULL) {
  model <- aci_enso_model(variant = "aci_code", hidden = hidden,
                          observations = observations)
  ov <- model$meta$vars$observed
  obs <- observed_trajectory(signal$t, as.matrix(signal[, ov, drop = FALSE]))
  init <- list(mean = as.numeric(signal[[hidden]][1L]),
               cov = matrix(0.1, 1L, 1L))
  bundle <- .compile_cgns_run(model, obs, spec)
  filter <- aci_filter(model, obs, init = init, conditional = spec)
  smoother <- aci_smoother(model, obs, filter = filter, init = init,
                           conditional = spec)
  metric <- aci(model, obs, init = init, conditional = spec, decompose = TRUE)
  list(model = model, obs = obs, init = init, bundle = bundle,
       filter = filter, smoother = smoother, metric = metric,
       observed_names = colnames(bundle$x))
}

# The realised coefficient arrays laid out as the fixture stores them: gxx by
# its diagonal only, which is lossless because the off-diagonal is exactly
# zero on every step and the test asserts it separately.
.partition_coefficients <- function(run) {
  co <- run$bundle$coefficients
  ov <- run$observed_names
  k <- run$bundle$k
  out <- list()
  for (i in seq_len(k)) out[[paste0("coef_Lx_", ov[i])]] <- co$Lx[i, 1L, ]
  for (i in seq_len(k)) out[[paste0("coef_fx_", ov[i])]] <- co$fx[, i]
  out[["coef_Ly"]] <- co$Ly[1L, 1L, ]
  out[["coef_fy"]] <- co$fy[, 1L]
  for (i in seq_len(k)) {
    out[[paste0("coef_gxx_diag_", ov[i])]] <- co$gxx[i, i, ]
  }
  out[["coef_gyy"]] <- co$gyy[1L, 1L, ]
  for (i in seq_len(k)) out[[paste0("coef_gyx_", ov[i])]] <- co$gyx[1L, i, ]
  out
}

.partition_series <- function(run) {
  list(
    ref_filter_mean = as.numeric(run$filter$mean),
    ref_filter_cov = as.numeric(run$filter$cov),
    ref_smooth_mean = as.numeric(run$smoother$mean),
    ref_smooth_cov = as.numeric(run$smoother$cov),
    ref_aci = as.numeric(run$metric$aci)
  )
}

.partition_max_error <- function(computed, ref, idx) {
  vapply(
    names(computed),
    function(k) max(abs(computed[[k]][idx] - ref[[k]])),
    numeric(1L)
  )
}

# The fixtures sample every twentieth step.  A defect at one of the 3800
# unsampled steps propagates forward through a sequential recursion and would
# be caught anyway, but a per-step output defect need not.  The full-record
# summary closes that gap: every series is reduced over all 4001 steps and
# graded against pinned scalars.
#
# These reductions are scored against max(1, |pinned|) rather than absolutely.
# A sum of absolute values over 4001 steps reaches 8002 here, and fifteen
# significant digits quantise a number that size at 4e-11, coarser than the
# gate.  Dividing by the pinned magnitude puts these on the same footing as
# the per-step comparisons instead of forcing the gate open to accommodate
# the file format.  Below unit magnitude the divisor is one, so a small mean
# is still graded absolutely.
.partition_summary_rows <- function(partition, computed, summary) {
  rows <- summary[summary$partition == partition, , drop = FALSE]
  expect_gt(nrow(rows), 0L)
  scaled <- function(a, b) abs(a - b) / max(1, abs(b))
  err <- numeric(0L)
  for (i in seq_len(nrow(rows))) {
    key <- rows$series[i]
    v <- computed[[key]]
    if (is.null(v)) v <- computed[[paste0("coef_", key)]]
    if (is.null(v)) v <- computed[[paste0("ref_", key)]]
    expect_false(is.null(v), info = sprintf("series %s absent", key))
    expect_identical(length(v), as.integer(rows$n[i]))
    err <- c(err, scaled(min(v), rows$min[i]), scaled(max(v), rows$max[i]),
             scaled(mean(v), rows$mean[i]),
             scaled(sum(abs(v)), rows$sum_abs[i]))
  }
  max(err)
}


test_that("partition oracle assets retain every manifest-pinned byte", {
  # The T_C-hidden zeroth-order entries were added to this manifest in W3e and
  # are graded by test-29-tc-zeroth-order.R.  They are pinned here, with the
  # rest of the file, because the hazard this block guards against is a fixture
  # nobody pins, and that hazard does not respect a test-file boundary.
  expected <- c(
    "enso6_partition_signal.csv",
    "enso6_partition_u_reference.csv",
    "enso6_partition_hW_reference.csv",
    "enso6_partition_tau_reference.csv",
    "enso6_partition_tau_reduced3_reference.csv",
    "enso6_partition_fullpath_summary.csv",
    "tc_coefficients_caseA.csv",
    "tc_outputs_caseA_intended.csv",
    "tc_outputs_caseA_literal.csv",
    "tc_outputs_caseB_window_monthly_intended.csv",
    "tc_outputs_caseB_window_monthly_literal.csv",
    "tc_outputs_caseB_window_head_intended.csv",
    "tc_outputs_caseB_window_head_literal.csv",
    "tc_d1_divergence.csv"
  )
  md5 <- .partition_manifest_hashes("md5")
  sha256 <- .partition_manifest_hashes("sha256")
  expect_setequal(names(md5), expected)
  expect_setequal(names(sha256), expected)
  expect_true(all(nchar(md5) == 32L))
  expect_true(all(nchar(sha256) == 64L))

  # This manifest is deliberately NOT merged into oracle-manifest.yml.  The two
  # describe different evidence classes, and merging them would put one
  # expected-name list in front of two readers with different expectations.
  # test-19 asserts that both files are fully pinned; this asserts they stay
  # disjoint, so a fixture cannot be pinned twice with two different hashes.
  expect_length(
    intersect(names(md5),
              .oracle_manifest_declared_csvs(
                .partition_oracle_file("oracle-manifest.yml"))),
    0L
  )

  for (name in names(md5)) {
    expect_identical(
      unname(tools::md5sum(.partition_oracle_file(name))),
      unname(md5[[name]]),
      info = sprintf(
        "%s changed without a partition-manifest refresh", name
      )
    )
  }
})


test_that("the pinned partition signal is the documented seeded simulation", {
  signal <- .partition_signal()
  expect_identical(nrow(signal), 4001L)
  expect_setequal(names(signal), c("t", "TC", "TE", "I", "u", "hW", "tau"))
  expect_equal(signal$t[1L], 0)
  expect_equal(signal$t[4001L], 20)
  expect_equal(diff(signal$t[1:2]), 0.005)
  # The generator is recorded, not re-run: aci_enso_model()'s own metadata says
  # the R simulator does not reproduce the MATLAB mixed scheme, so the path is
  # source-derived and the fixture pins its bytes rather than its seed.
  expect_identical(
    aci_enso_model(variant = "aci_code",
                   hidden = c("u", "hW", "tau"))$meta$matlab_simulator_parity,
    FALSE
  )
})


test_that("the graded observation sets are named, not inherited from a default", {
  # Default-proofing, asserted rather than assumed.  The tau default moved in
  # C2c; if it moves again, or if `observations` were to acquire a default for
  # the other partitions, these fixtures would silently start grading a
  # different estimand.  Nothing below reads a default -- this is the only
  # place the defaults are looked at, and it looks at them to record them.
  expect_identical(
    aci_enso_model(variant = "aci_code", hidden = "tau")$meta$observations,
    "reduced"
  )
  expect_identical(
    aci_enso_model(variant = "aci_code", hidden = "u")$meta$observations,
    "full"
  )
  # A model that declares its own observation set refuses a second
  # specification, which is why the reduced arm below is built two ways and
  # never one way with an override.
  signal <- .partition_signal()
  model <- aci_enso_model(variant = "aci_code", hidden = "tau",
                          observations = "reduced")
  obs <- observed_trajectory(
    signal$t,
    as.matrix(signal[, model$meta$vars$observed, drop = FALSE])
  )
  expect_error(
    .compile_cgns_run(model, obs,
                      aci_conditional(given = c("TC"), method = "mask")),
    class = "aci_error_nontarget"
  )
})


for (part in c("u", "hW", "tau")) {
  local({
    hidden <- part
    test_that(sprintf(
      "[source-derived harness] %s-hidden ENSO partition matches its pinned coefficients and recursions",
      hidden
    ), {
      signal <- .partition_signal()
      ref <- read.csv(.partition_oracle_file(
        sprintf("enso6_partition_%s_reference.csv", hidden)
      ))
      run <- .partition_run(hidden, "full", signal)

      # dispatch: five observed, one hidden, so the matrix branch, and the
      # affine batch realiser.  A scalar-branch regression would be a
      # different code path scoring the same numbers.
      expect_identical(run$bundle$realization, "affine_batch")
      expect_identical(run$bundle$k, 5L)
      expect_identical(run$bundle$l, 1L)
      expect_identical(run$model$meta$observations, "full")
      expect_setequal(run$observed_names,
                      setdiff(c("u", "hW", "TC", "TE", "tau", "I"), hidden))

      idx <- ref$index
      expect_identical(length(idx), 201L)
      expect_identical(idx[1L], 1L)
      expect_identical(idx[201L], 4001L)

      coefficients <- .partition_coefficients(run)
      series <- .partition_series(run)
      # `truth` is the signal's own hidden column, not a computed quantity.
      # Grading it ties the reference file and the summary file back to the
      # signal file, so a mismatched trio cannot pass.
      truth <- list(truth = signal[[hidden]])
      errors <- c(
        .partition_max_error(coefficients, ref, idx),
        .partition_max_error(series, ref, idx),
        .partition_max_error(truth, ref, idx)
      )
      expect_lt(max(errors), .partition_oracle_tolerance)

      # Structural invariants the fixture layout depends on.  gxx is stored by
      # its diagonal, so a non-zero off-diagonal would make the fixture lossy
      # rather than merely wrong, and gyx is exactly zero on every ACI_code
      # model, which is why the cross-noise terms are declared ungraded.
      co <- run$bundle$coefficients
      expect_identical(max(abs(co$gyx)), 0)
      expect_identical(
        max(abs(co$gxx - array(apply(co$gxx, 3L, function(z) diag(diag(z))),
                               dim(co$gxx)))),
        0
      )
      # The Gram inverse is realised on the N intervals, one slice short of
      # the N1 coefficient slices.  Any consumer that assumes N1 must pad.
      expect_identical(dim(co$gxx_weight)[3L], run$bundle$N1 - 1L)

      # Analytic boundary: the smoother is the filter at the final index.
      expect_lt(
        max(abs(as.numeric(run$smoother$mean)[4001L] -
                  as.numeric(run$filter$mean)[4001L]),
            abs(as.numeric(run$smoother$cov)[4001L] -
                  as.numeric(run$filter$cov)[4001L])),
        .partition_oracle_tolerance
      )

      summary <- read.csv(
        .partition_oracle_file("enso6_partition_fullpath_summary.csv")
      )
      expect_lt(
        .partition_summary_rows(hidden, c(coefficients, series, truth),
                                summary),
        .partition_oracle_tolerance
      )
    })
  })
}


test_that("[source-derived harness] the reduced tau partition matches its pinned run", {
  signal <- .partition_signal()
  ref <- read.csv(
    .partition_oracle_file("enso6_partition_tau_reduced3_reference.csv")
  )
  # The declared route: aci_enso_model(observations = "reduced") carries its own
  # prescribed-forcing specification (C2c), and every public entry point
  # applies it.  The fixture was produced before that route existed, by the
  # explicit specification on the five-channel model; the next block asserts
  # the two are the same run rather than two runs that happen to pass one gate.
  run <- .partition_run("tau", "reduced", signal)

  expect_identical(run$bundle$k, 3L)
  expect_identical(run$bundle$l, 1L)
  expect_identical(run$observed_names, c("TC", "TE", "I"))
  expect_identical(run$model$meta$observations, "reduced")

  idx <- ref$index
  coefficients <- .partition_coefficients(run)
  series <- .partition_series(run)
  truth <- list(truth = signal$tau)
  errors <- c(
    .partition_max_error(coefficients, ref, idx),
    .partition_max_error(series, ref, idx),
    .partition_max_error(truth, ref, idx)
  )
  expect_lt(max(errors), .partition_oracle_tolerance)

  summary <- read.csv(
    .partition_oracle_file("enso6_partition_fullpath_summary.csv")
  )
  expect_lt(
    .partition_summary_rows("tau_reduced3", c(coefficients, series, truth),
                            summary),
    .partition_oracle_tolerance
  )
})


test_that("the declared reduced estimand is the explicit prescribed-forcing run", {
  # The fixture's own construction: the five-channel model with the reduction
  # supplied by the caller.  Identity, not a tolerance -- the declaration is
  # meant to BE this specification, so anything but bit-equality is a defect,
  # and a tolerance here would hide the case where the declaration drifts to a
  # nearby-but-different reduction.
  signal <- .partition_signal()
  declared <- .partition_run("tau", "reduced", signal)
  explicit <- .partition_run(
    "tau", "full", signal,
    aci_conditional(given = c("u", "hW"), method = "reduce")
  )
  expect_identical(explicit$observed_names, c("TC", "TE", "I"))
  expect_identical(.partition_coefficients(declared),
                   .partition_coefficients(explicit))
  expect_identical(.partition_series(declared), .partition_series(explicit))
})


test_that("the two conditional strategies are the same operator on this model", {
  # gxx is diagonal and gyx is zero here, so masking the non-target rows of the
  # Gram inverse gives exactly the gain deleting those channels gives.  Pinned
  # as an identity rather than a tolerance: the C3 measurement was identical(),
  # max difference 0, on every one of the five series.
  #
  # Both arms are built on the five-channel model, because a model that
  # declares its own reduction will not accept either specification.
  signal <- .partition_signal()
  reduced <- .partition_run(
    "tau", "full", signal,
    aci_conditional(given = c("u", "hW"), method = "reduce")
  )
  inflated <- .partition_run(
    "tau", "full", signal,
    aci_conditional(given = c("u", "hW"), method = "mask")
  )
  expect_identical(.partition_series(reduced), .partition_series(inflated))
})


test_that("the reduced tau observation is a DIFFERENT estimand from the full one", {
  # ACI_code's tau script asserts that assimilating (T_C, T_E, I) with u and
  # h_W prescribed equals assimilating all five.  It does not, on this path.
  # The gap is pinned so that a future change cannot quietly make the two
  # agree and be read as a fix.  Recorded C3 measurements: filter mean 2.470e-1
  # max, ACI 7.765e-1 max, ACI Pearson 0.90455.
  #
  # This is also the assertion that makes the two tau fixtures worth shipping
  # separately: they are not two spellings of one run.
  signal <- .partition_signal()
  full <- .partition_run("tau", "full", signal)
  reduced <- .partition_run("tau", "reduced", signal)
  expect_identical(full$bundle$k, 5L)
  expect_identical(reduced$bundle$k, 3L)
  gap <- max(abs(as.numeric(full$metric$aci) - as.numeric(reduced$metric$aci)))
  expect_gt(gap, 0.5)
  expect_lt(
    abs(max(abs(as.numeric(full$filter$mean) -
                  as.numeric(reduced$filter$mean))) - 0.2470),
    1e-3
  )
})
