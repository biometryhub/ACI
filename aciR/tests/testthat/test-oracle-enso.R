# Independent-oracle grade of the stochastic ENSO model ------------------------
#
# The largest system the package expresses, and the only one whose observation
# noise varies in time. That last property is why this grade matters beyond the
# model itself: the Walker circulation's observation-noise variance is
# multiplicative in its own state and swings by more than a factor of three
# across the record, so a filter that read it at its first step alone would be
# integrating a different system and would still produce a plausible answer.
#
# Scope. The harness transcribes the reference's coefficient construction and
# the vector recursions, but generates its observed path by plain
# Euler-Maruyama rather than the reference's mixed scheme, because the
# reference's wind-burst update carries a correction this package does not
# reproduce, for reasons recorded in design/2026-08-13_milstein_anomaly.md.
# The filter is graded on whatever signal it is driven with, so this grades the
# coefficients and the recursions, NOT the reference's particular realisation.

.aci_enso_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

.aci_enso_pieces <- function() {
  signal <- read.csv(.aci_enso_fixture("enso_signal.csv"))
  list(
    signal = signal,
    comp = aci_enso_components(
      signal$T_C, signal$T_E, signal$I, time = signal$t
    ),
    x = rbind(signal$T_C, signal$T_E, signal$I),
    dt = 0.005,
    mu0 = c(6.9136e-04, -0.0028, -0.0256),
    R0 = 0.01 * diag(3L)
  )
}

test_that("the ENSO fixture exercises a genuinely time-varying noise", {
  # Without this the grade would say nothing about the property that makes
  # this model different from every other the package carries.
  p <- .aci_enso_pieces()
  observation_variance <- p$comp$S_xoS_x[3L, 3L, ]
  expect_gt(
    max(observation_variance) / min(observation_variance), 2
  )
  latent_variance <- p$comp$S_yoS_y[3L, 3L, ]
  expect_gt(max(latent_variance) / min(latent_variance), 2)

  # The coupling and self-drift move too, and the seasonal terms are alive.
  expect_gt(diff(range(p$comp$L_x[1L, 1L, ])), 1e-3)
  expect_gt(diff(range(p$comp$L_y[1L, 3L, ])), 1e-4)
  expect_gt(diff(range(p$comp$f_x[2L, ])), 1e-4)
})

test_that("the ENSO model reproduces the MATLAB oracle to 1e-6", {
  p <- .aci_enso_pieces()
  ref <- read.csv(.aci_enso_fixture("enso_reference.csv"))

  filt <- aci_filter(p$x, p$comp, p$dt, mu0 = p$mu0, R0 = p$R0)
  smooth <- aci_smoother(p$x, p$comp, p$dt, filt)
  metric <- aci_metric(filt, smooth)
  idx <- ref$index

  expect_lt(max(abs(filt$mean[1L, idx] - ref$fm1)), 1e-6)
  expect_lt(max(abs(filt$mean[2L, idx] - ref$fm2)), 1e-6)
  expect_lt(max(abs(filt$mean[3L, idx] - ref$fm3)), 1e-6)
  expect_lt(max(abs(filt$cov[1L, 1L, idx] - ref$fc11)), 1e-6)
  expect_lt(max(abs(filt$cov[3L, 3L, idx] - ref$fc33)), 1e-6)
  expect_lt(max(abs(filt$cov[1L, 3L, idx] - ref$fc13)), 1e-6)
  expect_lt(max(abs(smooth$mean[1L, idx] - ref$sm1)), 1e-6)
  expect_lt(max(abs(smooth$mean[3L, idx] - ref$sm3)), 1e-6)
  expect_lt(max(abs(smooth$cov[1L, 1L, idx] - ref$sc11)), 1e-6)
  expect_lt(max(abs(smooth$cov[3L, 3L, idx] - ref$sc33)), 1e-6)
  expect_lt(max(abs(metric[idx] - ref$ACI_metric)), 1e-6)

  expect_gt(max(ref$ACI_metric), 0.1)
})

test_that("the conditional question runs on the ENSO system", {
  # The configuration the case study is built around: what the central-Pacific
  # temperature says about the latent state, given that the eastern-Pacific
  # temperature and the Walker circulation are also observed.
  p <- .aci_enso_pieces()
  conditioned <- aci_conditional(p$comp, target = 1L)

  filt <- aci_filter(p$x, conditioned, p$dt, mu0 = p$mu0, R0 = p$R0)
  smooth <- aci_smoother(p$x, conditioned, p$dt, filt)
  metric <- aci_metric(filt, smooth)

  expect_true(all(is.finite(metric)))
  expect_true(all(metric >= 0))

  # Conditioning changes the estimand, so the answer must differ from the
  # unconditional one; equality would mean the masking did nothing.
  base_filt <- aci_filter(p$x, p$comp, p$dt, mu0 = p$mu0, R0 = p$R0)
  base <- aci_metric(
    base_filt, aci_smoother(p$x, p$comp, p$dt, base_filt)
  )
  expect_gt(max(abs(metric - base)), 1e-3)

  # The masked weight tracks the time-variation of the covariance it is
  # derived from. Conditioning on the central-Pacific temperature, whose own
  # observation noise is constant, correctly gives a constant weight.
  # Conditioning on the Walker circulation, whose noise is multiplicative in
  # its own state, must give a weight that moves. Building either from the
  # first step alone would be the defect this package already found once.
  on_temperature <- aci_conditional(p$comp, target = 1L)
  on_circulation <- aci_conditional(p$comp, target = 3L)
  expect_length(dim(on_circulation$S_xoS_x_inv), 3L)
  weight <- on_circulation$S_xoS_x_inv[3L, 3L, ]
  expect_gt(max(weight) / min(weight), 2)
  expect_equal(
    on_temperature$S_xoS_x_inv[1L, 1L, 1L],
    1 / p$comp$S_xoS_x[1L, 1L, 1L], tolerance = 1e-12
  )
})

test_that("an out-of-domain Walker circulation is refused", {
  # Its observation noise vanishes at both ends of (0, 4), leaving the
  # covariance the filter inverts singular. Refused rather than regularised.
  p <- aci_enso_parameters()
  tt <- seq(0, by = 0.005, length.out = 50L)
  for (bad in c(0, 4, -0.1, 4.2)) {
    expect_error(
      aci_enso_components(
        T_C = rep(0.01, 50L), T_E = rep(0.01, 50L),
        I = c(rep(1.5, 49L), bad), p = p, time = tt
      ),
      "strictly inside"
    )
  }
})

test_that("the ENSO parameters are derived, not restated", {
  p <- aci_enso_parameters()
  # The reference derives most parameters from a few physical ones; a caller
  # changing the diversity factor must see the derived ones move with it.
  q <- aci_enso_parameters(factor = 0.5)
  expect_false(isTRUE(all.equal(p$r, q$r)))
  expect_false(isTRUE(all.equal(p$gamma_C, q$gamma_C)))
  expect_false(isTRUE(all.equal(p$sigma_u, q$sigma_u)))
  expect_equal(p$r_C, p$gamma_C * p$b_0 * p$mu / 2)
  expect_equal(p$delta_h, p$alpha_2 * p$b_0 * p$mu)
})
