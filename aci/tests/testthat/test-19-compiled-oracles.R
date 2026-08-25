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

.compiled_oracle_manifest_hashes <- function(algorithm = "md5") {
  lines <- trimws(readLines(
    .compiled_oracle_file("oracle-manifest.yml"), warn = FALSE
  ))
  entries <- grep("^[A-Za-z_]+\\.csv:$", lines)
  hashes <- character(0L)
  for (i in entries) {
    block <- lines[seq.int(i + 1L, min(i + 3L, length(lines)))]
    hit <- grep(sprintf("^%s: [0-9a-f]+$", algorithm), block, value = TRUE)
    expect_length(hit, 1L)
    hashes[sub(":$", "", lines[i])] <- sub(
      sprintf("^%s: ", algorithm), "", hit[1L]
    )
  }
  hashes
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
  cgns_model(
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

.compiled_oracle_mv_model <- function() {
  noise <- matrix(
    c(
      0.60, 0.10, 0.25, 0.05,
      0.20, 0.50, 0.10, 0.30,
      0.15, 0.05, 0.70, 0.10,
      0.05, 0.10, 0.15, 0.55
    ),
    4L, 4L
  )
  sx <- noise[1:2, , drop = FALSE]
  sy <- noise[3:4, , drop = FALSE]
  cgns_model(
    Lx = function(t, x) matrix(
      c(0.8 + 0.3 * x[1L], 0.1 * x[2L],
        0.2 * sin(t), 0.6 - 0.2 * x[1L]),
      2L, 2L
    ),
    fx = function(t, x) c(0.4 - 0.5 * x[1L], -0.3 * x[2L] + 0.2),
    Ly = function(t, x) matrix(
      c(-1.2 + 0.1 * x[1L], 0.2,
        0.3, -0.9 - 0.1 * x[2L]),
      2L, 2L
    ),
    fy = function(t, x) c(0.5 - 0.2 * x[1L]^2,
                           0.1 - 0.15 * x[2L]),
    Sx1 = function(t, x) sx[, 1:2, drop = FALSE],
    Sx2 = function(t, x) sx[, 3:4, drop = FALSE],
    Sy1 = function(t, x) sy[, 1:2, drop = FALSE],
    Sy2 = function(t, x) sy[, 3:4, drop = FALSE],
    k = 2L, l = 2L, name = "independent-matrix-oracle"
  )
}

.compiled_oracle_online_moments <- function(model, t, x, init, ref) {
  # A complete Theorem-3 smoother on the prefix ending at n is exactly the
  # fixed-lag distribution p(y_j | x_0:n) requested by an oracle row (j, n).
  # Use the compiled online entry point itself for each distinct prefix; do not
  # restate its E/F or ordered-product recurrence in this test.
  means <- vector("list", nrow(ref))
  covs <- vector("list", nrow(ref))
  for (n in sort(unique(ref$n))) {
    rows <- which(ref$n == n)
    obs <- observed_trajectory(t[seq_len(n)], x[seq_len(n), , drop = FALSE])
    bundle <- .compile_cgns_run(model, obs)
    filter <- .cgns_filter_compiled(
      bundle, init = init, stepper = "explicit", nsub = 1L,
      validate = FALSE
    )
    online <- .smoother_thmD1_compiled(
      bundle, filter, validate = FALSE, warn_cost = FALSE
    )
    for (i in rows) {
      means[[i]] <- as.numeric(online$mean[ref$j[i], ])
      covs[[i]] <- matrix(
        online$cov[, , ref$j[i]], bundle$l, bundle$l
      )
    }
  }
  list(mean = means, cov = covs)
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
})


test_that("[authors source] compiled scalar dyad matches the pinned MATLAB output", {
  signal <- read.csv(
    .compiled_oracle_file("dyad_signal_x.csv"), header = FALSE
  )
  ref <- read.csv(.compiled_oracle_file("dyad_reference.csv"))
  obs <- observed_trajectory(signal$V1, matrix(signal$V2, ncol = 1L))
  run <- .compiled_oracle_run(
    model_dyad(), obs,
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
      model_predator_prey(hidden = case$hidden), obs,
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
  model <- model_enso6(
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
    model_dyad(), signal$V1, matrix(signal$V2, ncol = 1L),
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
  got <- .forward_cir_compiled(
    bundle,
    init = list(mean = 2, cov = matrix(0.1, 1L, 1L)),
    method = "l1_linf", eps = c(1e-1, 1e-2, 1e-3, 1e-4),
    min_M = NULL, quadrature = "simpson", simpson_close = "quadratic"
  )
  expected_subjective <- pmax(
    as.matrix(reference[, 4:7, drop = FALSE]) - bundle$dt,
    0
  )
  list(
    peak = max(abs(got$M[reference$j] - reference$peak)),
    objective = max(abs(got$tau[reference$j] - reference$objective)),
    subjective = max(abs(
      got$subjective[reference$j, , drop = FALSE] - expected_subjective
    ))
  )
}


test_that("[authors source] streamed dyad CIR matches the pinned MATLAB ranges", {
  errors <- .compiled_oracle_cir(
    model_dyad(), "dyad_signal_x.csv", "cir_range_reference.csv"
  )
  expect_lt(errors$peak, 1e-10)
  expect_lt(errors$objective, 1e-10)
  expect_lt(errors$subjective, 1e-12)
})


test_that("[independent MATLAB] streamed correlated CIR matches its ranges", {
  errors <- .compiled_oracle_cir(
    .compiled_oracle_cross_model(), "cross_signal_x.csv",
    "cir_cross_range_reference.csv"
  )
  expect_lt(errors$peak, 1e-10)
  expect_lt(errors$objective, 1e-10)
  expect_lt(errors$subjective, 1e-12)
})
