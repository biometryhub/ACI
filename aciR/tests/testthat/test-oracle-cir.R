# -- binding independent-oracle grade: the causal influence range --------------
#
# The online smoother and the causal influence range are graded against the
# reference MATLAB implementation, transcribed in oracle/aci_oracle_cir.m.
#
# The harness reuses the dyad signal fixture rather than generating its own, so
# this grade and the dyad grade are driven by the identical hash-pinned path.
# It uses the leading 2001 steps: the reference algorithm stores a staggered
# triangle quadratic in the record, which at the full length would need several
# gigabytes to hold quantities immediately reduced to scalars. aciR does not
# form that triangle at all, so the truncation is a limit of the harness and
# not of the package.
#
# Like the other oracle tests this one runs from the installed fixtures and
# never skips.

.aci_cir_oracle_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

.aci_cir_oracle_setup <- function() {
  signal <- read.csv(
    .aci_cir_oracle_fixture("dyad_signal_x.csv"), header = FALSE
  )
  x <- signal$V2[seq_len(2001L)]
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  comp <- aci_dyad_components(x, p)
  list(
    x = x, comp = comp, dt = 0.001,
    filt = aci_filter(x, comp, 0.001, mu0 = p$F_y / p$d_y, R0 = 0.1)
  )
}

test_that("the online smoother reproduces the MATLAB oracle to 1e-6", {
  s <- .aci_cir_oracle_setup()
  ref <- read.csv(.aci_cir_oracle_fixture("cir_online_reference.csv"))
  aux <- aciR:::.aci_online_aux(s$x, s$comp, s$dt, s$filt)

  got <- vapply(seq_len(nrow(ref)), function(i) {
    j <- ref$j[i]
    n <- ref$n[i]
    if (n == j) {
      return(c(s$filt$mean[j], s$filt$cov[j]))
    }
    k <- seq.int(j, n - 1L)
    d <- aux$cum_sign[k] * aux$cum_sign[j] *
      exp(aux$cum_log[k] - aux$cum_log[j])
    c(
      s$filt$mean[j] + sum(d * aux$innov_mean[k]),
      s$filt$cov[j] + sum(d * d * aux$innov_cov[k])
    )
  }, numeric(2))

  # The fixture must exercise a genuine range of lags, or this grades only the
  # trivial diagonal where the online smoother is the filter.
  expect_gt(length(unique(ref$n - ref$j)), 3L)
  expect_gt(max(ref$n - ref$j), 100L)

  expect_lt(max(abs(got[1L, ] - ref$online_mean)), 1e-6)
  expect_lt(max(abs(got[2L, ] - ref$online_cov)), 1e-6)
})

test_that("the causal influence range reproduces the MATLAB oracle to 1e-6", {
  s <- .aci_cir_oracle_setup()
  ref <- read.csv(.aci_cir_oracle_fixture("cir_range_reference.csv"))
  epsilon <- c(1e-1, 1e-2, 1e-3, 1e-4)

  # The saturation margin is stood down for the grade: the reference reports
  # every value it computes, so comparing like with like means comparing the
  # raw quantity. The margin is graded by its own tests.
  rng <- aci_cir(
    s$x, s$comp, s$dt, s$filt, window = ref$j, epsilon = epsilon,
    margin = 1e-9
  )

  expect_lt(max(abs(rng$peak - ref$peak)), 1e-6)
  expect_lt(max(abs(rng$objective - ref$objective)), 1e-6)

  subjective_ref <- as.matrix(
    ref[, c("subj_0.1", "subj_0.01", "subj_0.001", "subj_0.0001")]
  )
  expect_lt(max(abs(t(rng$subjective) - subjective_ref)), 1e-6)

  # Non-degeneracy: a fixture of zeros would pass every comparison above while
  # grading nothing at all.
  expect_gt(max(subjective_ref), 0)
  expect_gt(min(ref$peak), 0)
  expect_gt(max(ref$objective), 0)
})
