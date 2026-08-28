# Always-on external oracle grades for the compiled CGNS kernels ----------------
#
# The numerical fixtures in fixtures/oracles are copied byte-for-byte from the
# hash-pinned aciR oracle set.  `oracle-manifest.yml` is the authority for their
# provenance and scope.  In particular:
#
# * dyad and predator-prey are derived from the authors' supplied MATLAB;
# * ENSO grades source-derived coefficients and recursions on the harness path,
#   not the authors' particular stochastic realisation;
# * cross, mv, cir_cross, and mv_online are independent MATLAB transcriptions,
#   not authors-source numerical groundings.
#
# Keep those distinctions in test descriptions and failure messages.  A shared
# equation transcription is useful evidence, but it is not stronger evidence
# merely because it ran in MATLAB.
#
# The source-derived scalar ENSO partitions live in a second manifest,
# oracle-manifest-partitions.yml, and are graded by
# test-28-partition-oracles.R.  The byte-pinning test below covers BOTH
# manifests, because the hazard it guards against is a fixture nobody pins,
# and that hazard does not respect a file boundary.

.compiled_oracle_dir <- testthat::test_path("fixtures", "oracles")
.compiled_oracle_tolerance <- 1e-6

.compiled_oracle_file <- function(name) {
  path <- file.path(.compiled_oracle_dir, name)
  expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("always-on oracle fixture %s is missing", name)
  )
  path
}

# The reader lives in helper-oracle-manifest.R so that this file and
# test-28-partition-oracles.R share one regex.  Its character class admits
# digits; the test below proves it does.
.compiled_oracle_manifest_hashes <- function(algorithm = "md5") {
  .oracle_manifest_hashes(
    .compiled_oracle_file("oracle-manifest.yml"), algorithm
  )
}

.compiled_oracle_run <- function(model, obs, init) {
  # Grade production automatic realisation and scalar/matrix dispatch, not a
  # development-only compiler entry point.
  bundle <- .compile_cgns_run(model, obs)
  filter <- .cgns_filter_compiled(
    bundle, init = init, stepper = "explicit", nsub = 1L
  )
  smoother <- .cgns_smoother_compiled(bundle, filter)
  metric <- .gaussian_kl_path_compiled(
    bundle, smoother, filter, decompose = TRUE
  )
  list(bundle = bundle, filter = filter, smoother = smoother, metric = metric)
}

.compiled_oracle_scalar_error <- function(run, ref, idx) {
  max(
    abs(run$filter$mean[idx, 1L] - ref$filter_mean),
    abs(run$filter$cov[1L, 1L, idx] - ref$filter_cov),
    abs(run$smoother$mean[idx, 1L] - ref$smoother_mean),
    abs(run$smoother$cov[1L, 1L, idx] - ref$smoother_cov),
    abs(run$metric$total[idx] - ref$ACI_metric)
  )
}

.compiled_oracle_cross_model <- function() {
  aci_model(
    Lx = function(t, x) matrix(2 * x[1L], 1L, 1L),
    fx = function(t, x) 0.5 - 0.5 * x[1L],
    Ly = function(t, x) matrix(-0.5, 1L, 1L),
    fy = function(t, x) 1 - 2 * x[1L]^2,
    Sx1 = function(t, x) matrix(0.6, 1L, 1L),
    Sx2 = function(t, x) matrix(0.3, 1L, 1L),
    Sy1 = function(t, x) matrix(0.5, 1L, 1L),
    Sy2 = function(t, x) matrix(0.8, 1L, 1L),
    k = 1L, l = 1L, name = "independent-cross-oracle"
  )
}

.compiled_oracle_online_moments <- function(model, t, x, init, ref) {
  # The fixed-lag distribution p(y_j | x_0:n) an oracle row (j, n) requests is
  # what the online accumulator returns at anchor j and lag n - j, so the whole
  # fixture is graded from one forward auxiliary pass over the full record
  # instead of one compile + filter + backward sweep per distinct prefix.
  # The equivalence of the two routes is not assumed here: test-27's
  # "a lag of L admits the record through index j + L and no further" block
  # asserts it directly against .smoother_thmD1_compiled() on the prefixes,
  # and is the independent check that keeps this shortcut honest.
  obs <- observed_trajectory(t, x)
  bundle <- .compile_cgns_run(model, obs)
  filter <- .cgns_filter_compiled(
    bundle, init = init, stepper = "explicit", nsub = 1L, validate = FALSE
  )
  got <- .online_at_compiled(
    bundle, filter, .online_aux_compiled(bundle, filter), ref$j, ref$n - ref$j
  )
  list(
    mean = lapply(seq_len(nrow(ref)), function(i) as.numeric(got$mean[i, ])),
    cov = lapply(seq_len(nrow(ref)), function(i)
      matrix(got$cov[, , i], bundle$l, bundle$l))
  )
}


test_that("compiled oracle assets retain every manifest-pinned byte", {
  expected <- c(
    "dyad_signal_x.csv", "dyad_reference.csv",
    "cross_signal_x.csv", "cross_reference.csv",
    "cir_online_reference.csv", "cir_range_reference.csv",
    "predprey_signal.csv",
    "predprey_reference_predator_to_prey.csv",
    "predprey_reference_prey_to_predator.csv",
    "mv_signal.csv", "mv_reference.csv",
    "cir_cross_online_reference.csv", "cir_cross_range_reference.csv",
    "enso_signal.csv", "enso_reference.csv",
    "mv_online_reference.csv"
  )
  md5 <- .compiled_oracle_manifest_hashes("md5")
  sha256 <- .compiled_oracle_manifest_hashes("sha256")
  expect_setequal(names(md5), expected)
  expect_setequal(names(sha256), expected)
  expect_true(all(nchar(md5) == 32L))
  expect_true(all(nchar(sha256) == 64L))

  for (name in names(md5)) {
    expect_identical(
      unname(tools::md5sum(.compiled_oracle_file(name))),
      unname(md5[[name]]),
      info = sprintf(
        "%s changed without an oracle-manifest refresh", name
      )
    )
  }

  # Every entry of every manifest is actually pinned by a reader --------------
  #
  # The list above is hand-maintained, so on its own it cannot detect a
  # manifest entry no reader ever sees: the strict pattern skips the entry,
  # `names(md5)` never mentions it, and expect_setequal() compares two sets
  # that agree about a fixture neither of them pins.  That is exactly what the
  # digit-free character class did to `enso6_partition_*`.  The permissive
  # scan below finds block headers by shape rather than by name class, so the
  # two disagree the moment the strict reader goes blind to a name.
  manifests <- c("oracle-manifest.yml", "oracle-manifest-partitions.yml")
  for (mf in manifests) {
    path <- .compiled_oracle_file(mf)
    declared <- .oracle_manifest_declared_csvs(path)
    for (algorithm in c("md5", "sha256")) {
      pinned <- names(.oracle_manifest_hashes(path, algorithm))
      expect_setequal(pinned, declared)
      expect_identical(
        length(pinned), length(declared),
        info = sprintf("%s: %s reader saw %d of %d declared entries",
                       mf, algorithm, length(pinned), length(declared))
      )
    }
  }

  # The regex fix itself, stated as a fact rather than assumed.  The shared
  # reader returns all fourteen entries of the partitions manifest: the six
  # `enso6_partition_*` names, every one of which bears a digit, and the eight
  # T_C zeroth-order names added in W3e, none of which does.  Under the
  # previous `[A-Za-z_]` character class the first six returned character(0)
  # silently, which is what the fix and this assertion are about.
  partition_pinned <- .oracle_manifest_hashes(
    .compiled_oracle_file("oracle-manifest-partitions.yml"), "sha256"
  )
  expect_length(partition_pinned, 14L)
  scalar_partitions <- grep("^enso6_partition_", names(partition_pinned),
                            value = TRUE)
  expect_length(scalar_partitions, 6L)
  expect_true(all(grepl("[0-9]", scalar_partitions)))
  expect_length(grep("^tc_", names(partition_pinned), value = TRUE), 8L)

  # No fixture in the directory is unpinned by both manifests.  The reciprocal
  # of the check above: an entry with no file, and a file with no entry, are
  # different failures and both are worth catching.
  on_disk <- list.files(.compiled_oracle_dir, pattern = "\\.csv$")
  in_manifests <- unlist(lapply(manifests, function(mf)
    .oracle_manifest_declared_csvs(.compiled_oracle_file(mf))))
  expect_setequal(on_disk, in_manifests)
  expect_identical(anyDuplicated(in_manifests), 0L)
})


test_that("[authors source] compiled scalar dyad matches the pinned MATLAB output", {
  signal <- read.csv(
    .compiled_oracle_file("dyad_signal_x.csv"), header = FALSE
  )
  ref <- read.csv(.compiled_oracle_file("dyad_reference.csv"))
  obs <- observed_trajectory(signal$V1, matrix(signal$V2, ncol = 1L))
  run <- .compiled_oracle_run(
    aci_dyad_model(), obs,
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L))
  )
  expect_identical(run$bundle$realization, "dyad_directed")
  idx <- seq.int(1L, nrow(signal), by = 100L)
  expect_equal(length(idx), nrow(ref))
  expect_lt(
    .compiled_oracle_scalar_error(run, ref, idx),
    .compiled_oracle_tolerance
  )
})


test_that("[independent MATLAB] compiled correlated scalar path matches its oracle", {
  signal <- read.csv(
    .compiled_oracle_file("cross_signal_x.csv"), header = FALSE
  )
  ref <- read.csv(.compiled_oracle_file("cross_reference.csv"))
  obs <- observed_trajectory(signal$V1, matrix(signal$V2, ncol = 1L))
  run <- .compiled_oracle_run(
    .compiled_oracle_cross_model(), obs,
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L))
  )
  expect_identical(run$bundle$realization, "generic_closure_one_pass")
  idx <- seq.int(1L, nrow(signal), by = 100L)
  expect_true(run$bundle$correlated_noise)
  expect_lt(
    .compiled_oracle_scalar_error(run, ref, idx),
    .compiled_oracle_tolerance
  )
})


test_that("[authors source] compiled predator-prey paths match both directions", {
  signal <- read.csv(.compiled_oracle_file("predprey_signal.csv"))
  cases <- list(
    predator_to_prey = list(observed = signal$prey, hidden = "predator"),
    prey_to_predator = list(observed = signal$predator, hidden = "prey")
  )

  for (direction in names(cases)) {
    case <- cases[[direction]]
    ref <- read.csv(.compiled_oracle_file(sprintf(
      "predprey_reference_%s.csv", direction
    )))
    obs <- observed_trajectory(
      signal$t, matrix(case$observed, ncol = 1L)
    )
    run <- .compiled_oracle_run(
      aci_predprey_model(hidden = case$hidden), obs,
      init = list(mean = 4, cov = matrix(0.1, 1L, 1L))
    )
    expect_identical(run$bundle$realization, "generic_closure_one_pass")
    expect_lt(
      .compiled_oracle_scalar_error(run, ref, ref$index),
      .compiled_oracle_tolerance
    )
  }
})


test_that("[independent MATLAB] compiled correlated matrix path matches its oracle", {
  signal <- read.csv(.compiled_oracle_file("mv_signal.csv"))
  ref <- read.csv(.compiled_oracle_file("mv_reference.csv"))
  obs <- observed_trajectory(
    signal$t, as.matrix(signal[, c("x1", "x2")])
  )
  run <- .compiled_oracle_run(
    .compiled_oracle_mv_model(), obs,
    init = list(mean = c(0.8, 0.2), cov = 0.2 * diag(2L))
  )
  expect_identical(run$bundle$realization, "generic_closure_one_pass")
  idx <- ref$index
  errors <- c(
    filter_mean = max(
      abs(run$filter$mean[idx, 1L] - ref$fm1),
      abs(run$filter$mean[idx, 2L] - ref$fm2)
    ),
    filter_cov = max(
      abs(run$filter$cov[1L, 1L, idx] - ref$fc11),
      abs(run$filter$cov[1L, 2L, idx] - ref$fc12),
      abs(run$filter$cov[2L, 2L, idx] - ref$fc22)
    ),
    smoother_mean = max(
      abs(run$smoother$mean[idx, 1L] - ref$sm1),
      abs(run$smoother$mean[idx, 2L] - ref$sm2)
    ),
    smoother_cov = max(
      abs(run$smoother$cov[1L, 1L, idx] - ref$sc11),
      abs(run$smoother$cov[1L, 2L, idx] - ref$sc12),
      abs(run$smoother$cov[2L, 2L, idx] - ref$sc22)
    ),
    metric = max(abs(run$metric$total[idx] - ref$ACI_metric))
  )
  expect_true(run$bundle$correlated_noise)
  expect_gt(max(abs(ref$fc12)), 1e-3)
  expect_lt(max(errors), .compiled_oracle_tolerance)
})


test_that("[source-derived harness] compiled ENSO recursions match the pinned path", {
  signal <- read.csv(.compiled_oracle_file("enso_signal.csv"))
  ref <- read.csv(.compiled_oracle_file("enso_reference.csv"))
  model <- aci_enso_model(
    hidden = c("u", "hW", "tau"), variant = "aci_code"
  )
  obs <- observed_trajectory(
    signal$t, as.matrix(signal[, c("T_C", "T_E", "I")])
  )
  run <- .compiled_oracle_run(
    model, obs,
    init = list(
      mean = c(6.9136e-04, -0.0028, -0.0256),
      cov = 0.01 * diag(3L)
    )
  )
  expect_identical(run$bundle$realization, "affine_batch")
  idx <- ref$index
  errors <- c(
    filter_mean = max(
      abs(run$filter$mean[idx, 1L] - ref$fm1),
      abs(run$filter$mean[idx, 2L] - ref$fm2),
      abs(run$filter$mean[idx, 3L] - ref$fm3)
    ),
    filter_cov = max(
      abs(run$filter$cov[1L, 1L, idx] - ref$fc11),
      abs(run$filter$cov[3L, 3L, idx] - ref$fc33),
      abs(run$filter$cov[1L, 3L, idx] - ref$fc13)
    ),
    smoother_mean = max(
      abs(run$smoother$mean[idx, 1L] - ref$sm1),
      abs(run$smoother$mean[idx, 3L] - ref$sm3)
    ),
    smoother_cov = max(
      abs(run$smoother$cov[1L, 1L, idx] - ref$sc11),
      abs(run$smoother$cov[3L, 3L, idx] - ref$sc33)
    ),
    metric = max(abs(run$metric$total[idx] - ref$ACI_metric))
  )
  expect_identical(model$meta$matlab_simulator_parity, FALSE)
  expect_gt(max(run$bundle$coefficients$gxx[3L, 3L, ]) /
              min(run$bundle$coefficients$gxx[3L, 3L, ]), 2)
  expect_lt(max(errors), .compiled_oracle_tolerance)
})


test_that("[authors source] compiled scalar online smoother matches its oracle", {
  signal <- read.csv(
    .compiled_oracle_file("dyad_signal_x.csv"), header = FALSE
  )[seq_len(2001L), ]
  ref <- read.csv(.compiled_oracle_file("cir_online_reference.csv"))
  got <- .compiled_oracle_online_moments(
    aci_dyad_model(), signal$V1, matrix(signal$V2, ncol = 1L),
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L)), ref = ref
  )
  got_mean <- vapply(got$mean, `[[`, numeric(1L), 1L)
  got_cov <- vapply(got$cov, `[[`, numeric(1L), 1L)
  expect_gt(max(ref$n - ref$j), 100L)
  expect_lt(max(abs(got_mean - ref$online_mean)),
            .compiled_oracle_tolerance)
  expect_lt(max(abs(got_cov - ref$online_cov)),
            .compiled_oracle_tolerance)
})


test_that("[independent MATLAB] compiled correlated scalar online path matches", {
  signal <- read.csv(
    .compiled_oracle_file("cross_signal_x.csv"), header = FALSE
  )[seq_len(2001L), ]
  ref <- read.csv(
    .compiled_oracle_file("cir_cross_online_reference.csv")
  )
  got <- .compiled_oracle_online_moments(
    .compiled_oracle_cross_model(), signal$V1,
    matrix(signal$V2, ncol = 1L),
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L)), ref = ref
  )
  got_mean <- vapply(got$mean, `[[`, numeric(1L), 1L)
  got_cov <- vapply(got$cov, `[[`, numeric(1L), 1L)
  expect_gt(max(ref$n - ref$j), 100L)
  expect_lt(max(abs(got_mean - ref$online_mean)),
            .compiled_oracle_tolerance)
  expect_lt(max(abs(got_cov - ref$online_cov)),
            .compiled_oracle_tolerance)
})


test_that("[independent MATLAB] compiled correlated matrix online path matches", {
  signal <- read.csv(.compiled_oracle_file("mv_signal.csv"))
  ref <- read.csv(.compiled_oracle_file("mv_online_reference.csv"))
  got <- .compiled_oracle_online_moments(
    .compiled_oracle_mv_model(), signal$t,
    as.matrix(signal[, c("x1", "x2")]),
    init = list(mean = c(0.8, 0.2), cov = 0.2 * diag(2L)), ref = ref
  )
  got_mean <- do.call(rbind, got$mean)
  got_cov <- vapply(got$cov, function(R) c(R[1L, 1L], R[1L, 2L],
                                           R[2L, 2L]), numeric(3L))
  ref_mean <- as.matrix(ref[, c("om1", "om2")])
  ref_cov <- t(as.matrix(ref[, c("oc11", "oc12", "oc22")]))
  expect_gt(max(ref$n - ref$j), 1000L)
  expect_gt(max(abs(ref$oc12)), 1e-3)
  expect_lt(max(abs(got_mean - ref_mean)), .compiled_oracle_tolerance)
  expect_lt(max(abs(got_cov - ref_cov)), .compiled_oracle_tolerance)
})


.compiled_oracle_cir <- function(model, signal_file, range_file) {
  signal <- read.csv(
    .compiled_oracle_file(signal_file), header = FALSE
  )[seq_len(2001L), ]
  reference <- read.csv(.compiled_oracle_file(range_file))
  obs <- observed_trajectory(signal$V1, matrix(signal$V2, ncol = 1L))
  bundle <- .compile_cgns_run(model, obs)
  old <- options(aci.default_tol = 0)
  on.exit(options(old), add = TRUE)
  ## Ledger entry C2b: the fixture pins exactly these anchors, so the run is
  ## windowed to them rather than reducing all 2001 rows and discarding all
  ## but 15.  The graded comparison is unchanged; what it costs is not.
  got <- .forward_cir_compiled(
    bundle,
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L)),
    method = "l1_linf", epsilon = c(1e-1, 1e-2, 1e-3, 1e-4),
    min_M = NULL, quadrature = "simpson", simpson_close = "quadratic",
    anchors = reference$j
  )
  ## Ledger entry C2b: the default subjective read-out is now the reference
  ## script's counting convention (subj_CIR_idx * dt), so the pinned MATLAB
  ## columns are compared directly instead of being shifted down one grid
  ## step first.
  expected_subjective <- as.matrix(reference[, 4:7, drop = FALSE])
  list(
    anchors_ok = identical(got$meta$anchors, as.integer(reference$j)),
    peak = max(abs(got$M - reference$peak)),
    objective = max(abs(got$tau - reference$objective)),
    subjective = max(abs(got$subjective - expected_subjective)),
    lag_time_shift = max(abs(
      .forward_cir_compiled(
        bundle, init = list(mean = 2, cov = matrix(0.1, 1L, 1L)),
        method = "l1_linf", epsilon = c(1e-1, 1e-2, 1e-3, 1e-4),
        min_M = NULL, anchors = reference$j, convention = "lag_time"
      )$subjective -
      pmax(as.matrix(reference[, 4:7, drop = FALSE]) - bundle$dt, 0)
    ))
  )
}


test_that("[authors source] streamed dyad CIR matches the pinned MATLAB ranges", {
  errors <- .compiled_oracle_cir(
    aci_dyad_model(), "dyad_signal_x.csv", "cir_range_reference.csv"
  )
  expect_true(errors$anchors_ok)
  expect_lt(errors$peak, 1e-10)
  expect_lt(errors$objective, 1e-10)
  expect_lt(errors$subjective, 1e-12)
  expect_lt(errors$lag_time_shift, 1e-12)
})


test_that("[independent MATLAB] streamed correlated CIR matches its ranges", {
  errors <- .compiled_oracle_cir(
    .compiled_oracle_cross_model(), "cross_signal_x.csv",
    "cir_cross_range_reference.csv"
  )
  expect_true(errors$anchors_ok)
  expect_lt(errors$peak, 1e-10)
  expect_lt(errors$objective, 1e-10)
  expect_lt(errors$subjective, 1e-12)
  expect_lt(errors$lag_time_shift, 1e-12)
})
