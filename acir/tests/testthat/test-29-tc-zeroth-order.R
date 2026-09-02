# Always-on source-derived grades for the T_C-hidden zeroth-order partition ----
#
# ACI_code-main/ENSO_model_cond_ACI_T_C_unobs.m is the fifth and last ENSO
# script, and the only one whose partition is not an exact split of the
# six-state system: with T_C hidden the damping c_1(t, T_C) T_C is cubic in the
# hidden state, and the script restores conditional linearity by a zeroth-order
# Taylor expansion of c_1 about the climatology T_C = 0.  What it assimilates
# is a two-channel reduced model, observed (T_E, I), hidden T_C, with u, h_W
# and tau prescribed from their observed series (:999-1021).
#
# `oracle-manifest-partitions.yml` is the authority for the provenance of the
# fixtures graded here.  In particular:
#
# * these are SOURCE-DERIVED from an INDEPENDENT TRANSCRIPTION.  The pinned
#   values come from tools/fixtures/make-tc-fixtures.R, a statement-by-
#   statement R rendering of the script's filter (:1198-1257), smoother
#   (:1263-1324) and ACI (:1394-1397) sections, written without reference to
#   this package's kernels.  That is one step stronger than package-to-package
#   agreement and one step weaker than an authors-source fixture.  No MATLAB
#   was executed anywhere: nothing here grades the MATLAB realisation;
# * the driving paths are R simulations of the ACI_code ENSO system, not
#   MATLAB realisations, and aci_enso_model()'s own metadata records
#   matlab_simulator_parity = FALSE.  A path is an input to a filter, not an
#   oracle for one;
# * the seasonal coefficients are graded on the STATE-TIME phase, the
#   convention every acir constructor uses.  The MATLAB stored arrays hold
#   element 1 at state time and elements 2..N+1 one step ahead
#   (:1046/:1052 against :1140/:1150); the coefficient fixture carries both
#   columns, and the difference on a 4001-point path is up to 9.3e-04 in the
#   filter mean and 1.4e-02 in ACI;
# * the T_E-target conditional mask is EXACTLY INERT on this partition, so
#   nothing here grades the conditioning machinery.  It is asserted as an
#   identity below precisely because it is a negative result: the I row of Lx
#   is zero and the noise cross-Gram is zero, so the I precision multiplies
#   zero in the filter gain, and both the mask and the MATLAB first-slice
#   convention move nothing at all.  h_W remains the only live conditional
#   test in this package.
#
# Keep those distinctions in test descriptions and failure messages.

.tc_oracle_dir <- testthat::test_path("fixtures", "oracles")

# Gate rationale, recorded in the manifest's `tc_zeroth_order: validation:`
# block and restated here so a reader of the test need not go and find it.
# This is a regression gate on a pinned 15-significant-digit CSV, not a claim
# about accuracy.  The largest per-step quantity pinned is the ACI total, which
# reaches 4.2 on the case B window and is quantised at 5 * 4.2 * 1e-16 = 2.1e-15
# by the file format alone.  The measured agreement between this package's
# public chain and the transcription is 4.7e-15 on case A and 1.9e-14 on the
# 22001-step case B record, so the number gated is a file round trip plus the
# accumulation of one long sequential recursion.  1e-12 sits about fifty times
# above the worst measured value and far below any difference that would mean
# something: the smallest scientifically interesting divergence in this file is
# the D1 defect's 0.105 on the filter mean.  It is deliberately not tighter:
# 22000 sequential steps under a different BLAS could move the accumulated sums
# past the floor with no code change.
.tc_oracle_tolerance <- 1e-12

.tc_oracle_file <- function(name) {
  path <- file.path(.tc_oracle_dir, name)
  expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("always-on T_C fixture %s is missing", name)
  )
  path
}

# ACI_code-main/ENSO_model_cond_ACI_T_C_unobs.m parameters, transcribed here so
# the analytic identities below are checked against the source and not against
# the constructor's own arithmetic.
.tc_p <- local({
  fct <- 0.65                                  # :913
  list(factor = fct, b_0 = 2.5, mu = 0.5,      # :920, :922
       gamma_C = 0.75 * fct, gamma_E = 0.75 * fct,          # :938, :940
       r_C = 0.75 * fct * 2.5 * 0.5 / 2,                    # :942
       r_E = 3 * 0.75 * fct * 2.5 * 0.5 / 2,                # :944
       zeta_C = 0.75 * fct * 2.5 * 0.5 / 2,                 # :946
       zeta_E = 0.75 * fct * 2.5 * 0.5 / 2,                 # :948
       C_u = 0.03 * fct, lambda = 2 / 60, m = 2,            # :951, :955, :957
       sigma_C = 0.04 * sqrt(fct),                          # :962
       sigma_E = sqrt(5) * 1e-2 * sqrt(fct),                # :973
       dt = 0.005, k_dt = 100)                              # :894, :903
})

.tc_c1_zero <- function(t) {
  # c_1(t, 0), the zeroth-order substitution: 25*(0 + 0.75/7.5)^2 + 0.9.
  (25 * (0 + 0.75 / 7.5)^2 + 0.9) *
    (1 + 0.3 * sin(t * 2 * pi / 6 - pi / 6)) * .tc_p$factor
}

.tc_signal <- function() {
  d <- read.csv(.tc_oracle_file("enso6_partition_signal.csv"))
  data.frame(t = d$t, u = d$u, hW = d$hW, TC = d$TC, TE = d$TE,
             tau = d$tau, I = d$I)
}

.tc_model <- function(path, defect = FALSE) {
  aci_enso_model(hidden = "TC", variant = "aci_code",
                 approximation = "zeroth_order_c1", prescribed = path,
                 matlab_defect_compat = defect)
}

# One run through the PUBLIC entry point.  aci(keep = "paths") compiles once
# and returns the filter and the smoother it used, which is the whole graded
# set; three separate calls would compile the same 22001-step bundle three
# times.
.tc_run <- function(path, defect = FALSE) {
  m <- .tc_model(path, defect)
  ob <- observed_trajectory(path$t, cbind(TE = path$TE, I = path$I))
  init <- list(mean = path$TC[1L], cov = matrix(0.1, 1L, 1L))
  a <- aci(m, ob, init = init, decompose = TRUE, keep = "paths")
  list(model = m, obs = ob, init = init,
       filter_mean    = as.numeric(a$paths$filter$mean),
       filter_cov     = as.numeric(a$paths$filter$cov),
       smoother_mean  = as.numeric(a$paths$smoother$mean),
       smoother_cov   = as.numeric(a$paths$smoother$cov),
       aci_signal     = as.numeric(a$signal),
       aci_dispersion = as.numeric(a$dispersion),
       aci_total      = as.numeric(a$aci))
}

.tc_graded <- c("filter_mean", "filter_cov", "smoother_mean", "smoother_cov",
                "aci_signal", "aci_dispersion", "aci_total")

# Case B: a 110-model-year path, so that the script's own analysis window - 14
# model years with a two-year buffer at each end (:1351-1354) - actually
# exists.  Case A's 4000 steps are 40 model months and hold no such window.
#
# The window path is NOT shipped.  It is 2.0 MB, it regenerates byte-for-byte
# from its seed, and its SHA-256 is pinned in the manifest; the assertion below
# is what makes the seed as strong as the bytes would be, and it fails loudly
# and correctly if simulate() ever moves.
.tc_case_b <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    m <- aci_enso_model(variant = "aci_code", hidden = c("u", "hW", "tau"))
    s <- stats::simulate(m, seed = 4242L, t_end = 110, dt = .tc_p$dt)
    x <- as.matrix(s$obs$x); y <- as.matrix(s$hidden)
    path <- data.frame(t = s$obs$t, u = y[, 1L], hW = y[, 2L], TC = x[, 1L],
                       TE = x[, 2L], tau = y[, 3L], I = x[, 3L])
    a <- 3 * 12 * .tc_p$k_dt + 1
    cached <<- list(path = path, idx = a:(a + 14 * 12 * .tc_p$k_dt))
    cached
  }
})

.tc_case_b_runs <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    b <- .tc_case_b()
    cached <<- list(intended = .tc_run(b$path, FALSE),
                    literal  = .tc_run(b$path, TRUE))
    cached
  }
})


# 1. Constructor contract ------------------------------------------------------

test_that("the TC partition is admissible only as the zeroth-order approximation", {
  path <- .tc_signal()

  # No approximation named: the cubic damping is not conditionally linear, and
  # the message has to say which substitution makes it so.
  expect_error(aci_enso_model(hidden = "TC"), class = "aci_error_model_contract")
  expect_error(aci_enso_model(hidden = "TC"), "zeroth_order_c1")

  # A TC-containing hidden set that is not exactly {TC} keeps the pre-existing
  # exact-split error: no MATLAB script and no closed form covers it.
  for (h in list(c("TC", "hW"), c("TC", "u", "tau"), c("TC", "I"))) {
    expect_error(aci_enso_model(hidden = h, approximation = "zeroth_order_c1"),
                 "must remain observed", class = "aci_error_model_contract")
  }
  expect_error(aci_enso_model(hidden = "I"), "must remain observed",
               class = "aci_error_model_contract")
  expect_error(aci_enso_model(hidden = c("hW", "I")), "must remain observed",
               class = "aci_error_model_contract")

  # And the approximation is refused where there is nothing to approximate.
  expect_error(aci_enso_model(hidden = "hW", approximation = "zeroth_order_c1"),
               class = "aci_error_model_contract")
  expect_error(aci_enso_model(hidden = c("u", "hW", "tau"),
                              approximation = "zeroth_order_c1"),
               class = "aci_error_model_contract")

  # Unlike every other partition this model is not self-contained.
  expect_error(aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1"),
               "prescribed", class = "aci_error_model_contract")
  expect_error(
    aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1",
                   prescribed = path[c("t", "u", "hW")]),
    class = "aci_error_model_contract"
  )
  # A non-uniform grid is refused rather than interpolated over.
  bent <- path; bent$t[10L] <- bent$t[10L] + 1e-4
  expect_error(
    aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1",
                   prescribed = bent),
    "uniform", class = "aci_error_model_contract"
  )

  # The defect flag and the prescribed forcings belong to this branch alone.
  expect_error(aci_enso_model(hidden = "hW", matlab_defect_compat = TRUE),
               class = "aci_error_model_contract")
  expect_error(aci_enso_model(hidden = "hW", prescribed = path),
               class = "aci_error_model_contract")
  expect_error(aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1",
                              prescribed = path, matlab_defect_compat = NA),
               class = "aci_error_model_contract")

  # There is no five-channel TC estimand to ask for.
  expect_error(
    aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1",
                   prescribed = path, observations = "full"),
    class = "aci_error_model_contract"
  )

  m <- .tc_model(path)
  expect_s3_class(m, "cgns_model")
  expect_identical(m$k, 2L)
  expect_identical(m$l, 1L)
  expect_identical(m$meta$vars$observed, c("TE", "I"))
  expect_identical(m$meta$vars$hidden, "TC")
  expect_identical(m$meta$vars$prescribed, c("u", "hW", "tau"))
  expect_identical(m$meta$source_partition, "TC")
  expect_identical(m$meta$approximation, "zeroth_order_c1")
  expect_identical(m$meta$observations, "reduced")
  expect_identical(m$meta$coefficient_phase, "state_time")
  expect_false(m$meta$matlab_defect_compat$active)
  expect_true(.tc_model(path, TRUE)$meta$matlab_defect_compat$active)
  expect_identical(m$meta$target_obs_idx, 1L)
  expect_identical(m$meta$conditioning_obs_idx, 2L)
  expect_identical(m$meta$causal_link, "(TC) -> (TE) | (u,hW,tau,I)")
  expect_identical(m$meta$prescribed_grid$n, nrow(path))
  expect_equal(m$meta$prescribed_grid$dt, .tc_p$dt)

  # No batch realiser is attached: this is a distinct two-channel construction
  # and it compiles through the generic route.  Asserting that keeps the other
  # partitions' descriptor authentication meaningful.
  expect_null(.cgns_realizer_descriptor(m))
  ob <- observed_trajectory(path$t, cbind(TE = path$TE, I = path$I))
  expect_identical(.compile_cgns_run(m, ob, NULL)$realization,
                   "generic_closure_one_pass")
})

test_that("the TC model refuses to simulate and refuses a foreign grid", {
  path <- .tc_signal()
  m <- .tc_model(path)

  # Falling through to the generic simulator would integrate the linearised
  # hidden drift and return a plausible-looking path of a system nobody wrote
  # down.  The message has to name the constructor that can generate one.
  expect_false(m$meta$simulate_supported)
  expect_error(stats::simulate(m, seed = 1L, t_end = 1, dt = .tc_p$dt),
               class = "aci_error_model_contract")
  expect_error(stats::simulate(m, seed = 1L, t_end = 1, dt = .tc_p$dt),
               "aci_enso_model(hidden = c(\"u\", \"hW\", \"tau\")", fixed = TRUE)

  # The prescribed forcings are looked up by index, so an observation grid that
  # is not the grid they were supplied on is refused rather than resolved to a
  # neighbouring step.
  init <- list(mean = path$TC[1L], cov = matrix(0.1, 1L, 1L))
  n <- 200L
  shifted <- observed_trajectory(path$t[seq_len(n)] + 1e-4,
                                 cbind(TE = path$TE[seq_len(n)],
                                       I = path$I[seq_len(n)]))
  expect_error(aci_filter(m, shifted, init = init), class = "aci_error_dims")
  expect_error(aci(m, shifted, init = init), class = "aci_error_dims")
  expect_error(aci_smoother(m, shifted, init = init), class = "aci_error_dims")
  # Times past the end of the stored forcing would clamp to its last value.
  beyond <- observed_trajectory(path$t[seq_len(n)] + 20,
                                cbind(TE = path$TE[seq_len(n)],
                                      I = path$I[seq_len(n)]))
  expect_error(aci_filter(m, beyond, init = init), class = "aci_error_dims")

  # A leading sub-record is on the grid and is accepted, and so is a coarser
  # observation grid whose every time is still a forcing grid point: each
  # observation is then paired with the forcing actually stored at its time.
  head_obs <- observed_trajectory(path$t[seq_len(n)],
                                  cbind(TE = path$TE[seq_len(n)],
                                        I = path$I[seq_len(n)]))
  expect_s3_class(aci_filter(m, head_obs, init = init), "da_path")
  every2 <- seq(1L, 401L, by = 2L)
  coarse <- observed_trajectory(path$t[every2],
                                cbind(TE = path$TE[every2], I = path$I[every2]))
  expect_s3_class(aci_filter(m, coarse, init = init), "da_path")

  # Every other model is untouched by the check.
  d <- aci_dyad_model()
  sim <- stats::simulate(d, seed = 1L, t_end = 1, dt = 0.01)
  expect_null(d$meta$prescribed_grid)
  expect_s3_class(aci_filter(d, as_obs(sim), init = list(mean = 2, cov = 0.1)),
                  "da_path")
})


# 2. Coefficient identities ----------------------------------------------------

test_that("the TC coefficients are the six-state drift split about TC = 0", {
  path <- .tc_signal()
  m <- .tc_model(path)
  md <- .tc_model(path, TRUE)
  ts <- path$t[c(1L, 2L, 137L, 1500L, 4001L)]
  xs <- list(c(0, 1.6), c(-0.4, 0.2), c(0.9, 3.7))

  for (t in ts) for (x in xs) {
    # L_x = [-zeta_E; 0], constant in both arguments (:1024).
    expect_equal(m$Lx(t, x), matrix(c(-.tc_p$zeta_E, 0), 2L, 1L))
    # L_y is a time-only series: the whole point of the substitution.
    expect_equal(drop(m$Ly(t, x)), .tc_p$r_C - .tc_c1_zero(t))
    expect_identical(m$Ly(t, x), m$Ly(t, xs[[1L]]))
    # S_y = sigma_C (:1037); S_x = diag(sigma_E, sigma_I(I)) (:1046, :1144),
    # with the documented variance floor on the I channel.
    expect_equal(drop(m$Sy2(t, x)), .tc_p$sigma_C)
    Sx1 <- m$Sx1(t, x)
    expect_equal(Sx1[1L, 1L], .tc_p$sigma_E)
    expect_equal(Sx1[2L, 2L],
                 sqrt(.tc_p$lambda * x[2L] * (4 - x[2L]) + 1e-3 * .tc_p$lambda))
    expect_equal(Sx1[1L, 2L], 0)
    # The I row of f_x is -lambda (I - m), free of every prescribed channel.
    expect_equal(m$fx(t, x)[2L], -.tc_p$lambda * (x[2L] - .tc_p$m))
    # The defect flips exactly one term of f_y and nothing else.
    j <- round(t / .tc_p$dt) + 1L
    expect_equal(m$fy(t, x) - md$fy(t, x), .tc_p$gamma_C * path$hW[j])
    expect_identical(m$fx(t, x), md$fx(t, x))
    expect_identical(m$Ly(t, x), md$Ly(t, x))
    expect_identical(m$Lx(t, x), md$Lx(t, x))
  }

  # The construction is the exact partition's own drift with c1 frozen at
  # TC = 0, so it must agree with aci_enso_model(hidden = c("u","hW","tau"))
  # evaluated there: same algebra, one substitution, not two transcriptions.
  joint <- aci_enso_model(variant = "aci_code", hidden = c("u", "hW", "tau"))
  for (t in ts) for (x in xs) {
    j <- round(t / .tc_p$dt) + 1L
    y <- c(path$u[j], path$hW[j], path$tau[j])
    full <- joint$fx(t, c(0, x[1L], x[2L])) +
      drop(joint$Lx(t, c(0, x[1L], x[2L])) %*% y)
    # joint observes (TC, TE, I): rows 2 and 3 are the TE and I drifts, and
    # row 1 at TC = 0 is the TC drift with the zeroth-order damping already
    # applied, which is exactly f_y.
    expect_equal(m$fx(t, x), full[2:3])
    expect_equal(m$fy(t, x), full[1L])
  }
})


# 3. Fixture regression --------------------------------------------------------

test_that("[source-derived] the TC coefficients match the pinned transcription", {
  path <- .tc_signal()
  ref <- read.csv(.tc_oracle_file("tc_coefficients_caseA.csv"))
  expect_identical(nrow(ref), nrow(path))
  expect_equal(ref$t, path$t, tolerance = .tc_oracle_tolerance)

  m <- .tc_model(path)
  md <- .tc_model(path, TRUE)
  ob <- observed_trajectory(path$t, cbind(TE = path$TE, I = path$I))
  b <- .compile_cgns_run(m, ob, NULL)
  bd <- .compile_cgns_run(md, ob, NULL)
  co <- b$coefficients

  err <- function(a, e) max(abs(a - e))
  # The shipped state-time phase.  L_y and the T_E row of f_x are the only
  # phase-dependent quantities; they are compared against the state-time
  # columns, and the MATLAB-phase columns in the same file record what the
  # stored-array convention would have given.
  expect_lt(err(as.numeric(co$Ly), ref$L_y_state_time), .tc_oracle_tolerance)
  expect_lt(err(co$fx[, 1L], ref$f_x_TE_state_time), .tc_oracle_tolerance)
  expect_lt(err(co$fx[, 2L], ref$f_x_I), .tc_oracle_tolerance)
  expect_lt(err(as.numeric(co$fy), ref$f_y_intended), .tc_oracle_tolerance)
  expect_lt(err(as.numeric(bd$coefficients$fy), ref$f_y_literal),
            .tc_oracle_tolerance)
  # sigma_I carries the documented variance floor, so it grades against the
  # floored column.  The unfloored MATLAB column is in the same file, and the
  # gap between them is bounded here rather than left implicit.
  sigma_I <- sqrt(co$gxx[2L, 2L, ])
  expect_lt(err(sigma_I, ref$sigma_I_floored), .tc_oracle_tolerance)
  # The whole of the deviation from the MATLAB expression, stated as the exact
  # identity rather than as a bound: the floor adds 1e-3 * lambda to the
  # variance and nothing else.  Its size on this record is 3.7e-4 at most, and
  # it is inert - the I row of Lx is zero, so this entry never reaches the
  # posterior, and the T_E-target mask zeroes its precision in any case.
  expect_lt(err(sigma_I, sqrt(ref$sigma_I_matlab^2 + 1e-3 * .tc_p$lambda)),
            .tc_oracle_tolerance)
  expect_true(all(sigma_I > ref$sigma_I_matlab))
  expect_gt(max(abs(sigma_I - ref$sigma_I_matlab)), .tc_oracle_tolerance)
  expect_lt(max(abs(sigma_I - ref$sigma_I_matlab)), 1e-3)
  # The Gram is diagonal, which is what makes the T_E-target mask inert.
  expect_identical(max(abs(co$gxx[1L, 2L, ])), 0)
  expect_identical(max(abs(co$gyx)), 0)
  expect_identical(max(abs(co$Lx[2L, 1L, ])), 0)

  # The two phase conventions differ by more than round-off, which is why the
  # constructor documents the choice instead of leaving it to be discovered.
  expect_gt(max(abs(ref$L_y_state_time - ref$L_y_matlab_phase)), 1e-6)
  expect_gt(max(abs(ref$f_x_TE_state_time - ref$f_x_TE_matlab_phase)), 1e-6)
})

test_that("[source-derived] the TC case A moments match the pinned transcription", {
  path <- .tc_signal()
  for (defect in c(FALSE, TRUE)) {
    ref <- read.csv(.tc_oracle_file(sprintf(
      "tc_outputs_caseA_%s.csv", if (defect) "literal" else "intended")))
    run <- .tc_run(path, defect)
    expect_identical(nrow(ref), length(run$filter_mean))
    for (q in .tc_graded) {
      expect_lt(
        max(abs(run[[q]] - ref[[if (q == "aci_total") "aci_total" else q]])),
        .tc_oracle_tolerance,
        label = sprintf("case A %s, %s", if (defect) "literal" else "intended", q)
      )
    }
  }
})

test_that("[source-derived] the TC case B window matches the pinned transcription", {
  b <- .tc_case_b()
  # The 22001-point record regenerates from its seed; the manifest pins the
  # SHA-256 of the window it produces, and the window arithmetic is the
  # script's own (:1351-1354) with sim_year_start = 3, ACI_period_years = 14.
  expect_identical(nrow(b$path), 22001L)
  expect_identical(length(b$idx), 16801L)
  expect_equal(b$path$t[min(b$idx)], 18)
  expect_equal(b$path$t[max(b$idx)], 102)

  runs <- .tc_case_b_runs()
  monthly <- b$idx[seq(1L, length(b$idx), by = .tc_p$k_dt)]
  headwin <- b$idx[seq_len(200L)]
  for (defect in c(FALSE, TRUE)) {
    tag <- if (defect) "literal" else "intended"
    run <- runs[[tag]]
    for (part in c("monthly", "head")) {
      idx <- if (identical(part, "monthly")) monthly else headwin
      ref <- read.csv(.tc_oracle_file(sprintf(
        "tc_outputs_caseB_window_%s_%s.csv", part, tag)))
      expect_identical(nrow(ref), length(idx))
      expect_lt(max(abs(ref$t - b$path$t[idx])), .tc_oracle_tolerance)
      for (q in .tc_graded) {
        expect_lt(max(abs(run[[q]][idx] - ref[[q]])), .tc_oracle_tolerance,
                  label = sprintf("case B %s %s, %s", part, tag, q))
      }
    }
  }
})


# 4. The D1 divergence, pinned as a measurement --------------------------------

test_that("[source-derived] the omitted gamma_C h_W term moves only the mean channel", {
  # The reference script's f_y omits gamma_C * h_W (:1053, :1151).  The same
  # script's simulator drift includes it (:1124), h_W is prescribed and
  # observed here, and the sibling scripts carry the term in the corresponding
  # T_C coefficient rows (tau_unobs.m:1138; u_h_W_tau_unobs.m:1016,:1131), so
  # this is a source defect and not a convention.  acir includes the term by
  # default and reproduces the omission under matlab_defect_compat = TRUE.
  #
  # The sharpest available guard that the flag touches only the mean channel:
  # f_y enters neither the Riccati equation nor the backward covariance ODE, so
  # both covariances must be BIT-identical, not merely close.
  ref <- read.csv(.tc_oracle_file("tc_d1_divergence.csv"))
  path <- .tc_signal()
  b <- .tc_case_b()
  runs <- .tc_case_b_runs()

  cases <- list(
    list(name = "A whole record", a = .tc_run(path, FALSE),
         b = .tc_run(path, TRUE), idx = seq_len(nrow(path)), t = path$t),
    list(name = "B 14-year window", a = runs$intended, b = runs$literal,
         idx = b$idx, t = b$path$t))

  trapz <- function(t, v) sum(diff(t) * (utils::head(v, -1L) + v[-1L]) / 2)

  for (case in cases) {
    rows <- ref[ref$case == case$name, ]
    expect_identical(nrow(rows), 8L)
    for (q in .tc_graded) {
      row <- rows[rows$quantity == q, ]
      expect_identical(nrow(row), 1L)
      va <- case$a[[q]][case$idx]; vb <- case$b[[q]][case$idx]
      if (isTRUE(row$identical)) {
        # Asserted as identity, not as a tolerance: anything non-zero here
        # means f_y has reached a covariance recursion.
        expect_identical(va, vb, label = sprintf("%s %s", case$name, q))
        expect_identical(max(abs(va - vb)), 0)
      } else {
        expect_false(identical(va, vb))
        expect_equal(max(abs(va - vb)), row$max_abs, tolerance = 1e-6,
                     label = sprintf("%s %s max abs", case$name, q))
        expect_equal(sqrt(mean((va - vb)^2)), row$rmsd, tolerance = 1e-6,
                     label = sprintf("%s %s rmsd", case$name, q))
      }
    }
    # The headline figure: the time-integrated ACI over the analysis window.
    row <- rows[rows$quantity == "aci_total_integrated", ]
    ia <- trapz(case$t[case$idx], case$a$aci_total[case$idx])
    ib <- trapz(case$t[case$idx], case$b$aci_total[case$idx])
    expect_equal(ia, row$max_abs, tolerance = 1e-6)
    expect_equal(ib, row$rmsd, tolerance = 1e-6)
    expect_gt(ib, ia)
  }

  # On the script's own 14-model-year window the omission roughly doubles the
  # time-integrated ACI.  That is the number that makes the defect a decision
  # rather than a footnote.
  row <- ref[ref$case == "B 14-year window" &
               ref$quantity == "aci_total_integrated", ]
  expect_gt(row$rmsd / row$max_abs, 1.9)
  expect_lt(row$rmsd / row$max_abs, 2.1)
})


# 5. Conditioning is exactly inert on this partition ---------------------------

test_that("[negative result] the TE-target mask moves nothing on the TC partition", {
  # A negative regression, asserted as exact equality rather than as a bound.
  # The I row of L_x is zero and the noise cross-Gram is zero, so the I-channel
  # precision multiplies zero in the filter gain: masking it, and leaving the
  # first slice unmasked as the reference scripts do, are both structurally
  # incapable of moving any quantity.  If any of these stops being an identity,
  # the construction has changed and the fixtures above no longer grade what
  # they claim to.  This partition is therefore NOT a discriminating test of
  # the conditioning API; h_W remains the only live one.
  path <- .tc_signal()
  m <- .tc_model(path)
  ob <- observed_trajectory(path$t, cbind(TE = path$TE, I = path$I))
  init <- list(mean = path$TC[1L], cov = matrix(0.1, 1L, 1L))

  declared <- m$meta$estimand_nontarget
  expect_s3_class(declared, "aci_conditional_spec")
  expect_identical(declared$target, 1L)
  expect_identical(declared$method, "mask")
  expect_identical(declared$first_step, "uniform")
  # The declaration is not composable, as C2c settled for the tau partition.
  expect_error(aci_filter(m, ob, init = init,
                          conditional = aci_conditional(given = 2L,
                                                        method = "mask")),
               class = "aci_error_nontarget")

  variants <- list(
    unmasked = NULL,
    matlab_first_step = aci_conditional(target = 1L, method = "mask",
                                        first_step = "matlab"))
  base <- .tc_run(path, FALSE)
  for (nm in names(variants)) {
    alt <- m
    alt$meta$estimand_nontarget <- variants[[nm]]
    a <- aci(alt, ob, init = init, decompose = TRUE, keep = "paths")
    got <- list(filter_mean = as.numeric(a$paths$filter$mean),
                filter_cov = as.numeric(a$paths$filter$cov),
                smoother_mean = as.numeric(a$paths$smoother$mean),
                smoother_cov = as.numeric(a$paths$smoother$cov),
                aci_signal = as.numeric(a$signal),
                aci_dispersion = as.numeric(a$dispersion),
                aci_total = as.numeric(a$aci))
    for (q in .tc_graded) {
      expect_identical(got[[q]], base[[q]], label = sprintf("%s, %s", nm, q))
      expect_identical(max(abs(got[[q]] - base[[q]])), 0)
    }
  }
})


# 6. The retired claim cannot come back ----------------------------------------

test_that("the TC partition is no longer recorded as unsupported", {
  # Before this partition existed, aci_enso_model() carried
  # meta$unsupported_partitions$TC saying the zeroth-order case was "not the
  # exact six-state CGNS split constructed here".  It is constructed here now,
  # so the sentence has to go and has to stay gone.
  path <- .tc_signal()
  models <- list(aci_enso_model(), aci_enso_model(hidden = "u"),
                 aci_enso_model(hidden = "hW"),
                 aci_enso_model(hidden = "tau", observations = "reduced"),
                 aci_enso_model(hidden = "tau", observations = "full"),
                 aci_enso_model(hidden = c("u", "hW", "tau")),
                 .tc_model(path))
  for (m in models) {
    expect_null(m$meta$unsupported_partitions)
    expect_null(m$meta$unsupported_partitions$TC)
  }
  # The claim itself, wherever a constructor might still deparse it.
  src <- paste(unlist(lapply(
    c("aci_enso_model", ".enso6_tc_zeroth_order", ".enso6_prescribed_grid"),
    function(f) deparse(get(f, envir = asNamespace("acir"))))), collapse = " ")
  expect_false(grepl("not the exact six-state", src, fixed = TRUE))
  expect_false(grepl("unsupported_partitions", src, fixed = TRUE))
})


# 7. Strategy equivalence ------------------------------------------------------

test_that("inflate and prescribed forcing agree on the TC partition", {
  # The same family result C3 and C4 found for the other ENSO partitions.
  # Here it is a weak check by construction, because the mask is inert, but it
  # does exercise the reduction machinery on a k = 2 model.
  path <- .tc_signal()
  m <- .tc_model(path)
  ob <- observed_trajectory(path$t, cbind(TE = path$TE, I = path$I))
  init <- list(mean = path$TC[1L], cov = matrix(0.1, 1L, 1L))

  pf <- m
  pf$meta$estimand_nontarget <- aci_conditional(given = 2L,
                                                method = "reduce")
  a_inflate <- aci(m, ob, init = init, decompose = TRUE)
  a_prescribed <- aci(pf, ob, init = init, decompose = TRUE)
  expect_lt(max(abs(a_inflate$aci - a_prescribed$aci)), 1e-14)
  expect_lt(max(abs(a_inflate$signal - a_prescribed$signal)), 1e-14)
  expect_lt(max(abs(a_inflate$dispersion - a_prescribed$dispersion)), 1e-14)

  red <- aci_conditional_reduce(m, ob, aci_conditional(given = 2L,
                                                       method = "reduce"))
  expect_identical(red$model$k, 1L)
  expect_identical(red$obs$k, 1L)
  expect_null(red$model$meta$estimand_nontarget)
})
