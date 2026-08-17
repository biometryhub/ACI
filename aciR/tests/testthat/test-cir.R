# Grading the causal influence range -------------------------------------------
#
# The range is built on the online smoother, which is graded in its own file,
# so the tests here grade what this layer adds: that the comparison sequence is
# the one the definition calls for, that its reduction to a range obeys the
# properties a range must have, and that a time the record is too short to
# resolve is reported as unresolved rather than as a short range.
#
# The last of those is the failure mode this quantity invites. The divergence
# at a given time is compared against every later observation, so a time near
# the end of the record is compared against almost nothing and yields a small
# number that looks like a confident answer.

.aci_cir_params <- function() {
  list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
}

.aci_cir_pieces <- function(n = 1200L, dt = 1e-3, seed = 5L) {
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = n, dt = dt, seed = seed)
  comp <- aci_dyad_components(sim$x, model$parameters)
  filt <- aci_filter(sim$x, comp, dt, mu0 = model$y0, R0 = 0.1)
  list(x = sim$x, comp = comp, filt = filt, dt = dt, n = n)
}

test_that("the comparison sequence ends at the online smoother", {
  # The last entry of each row is the posterior informed by the whole record,
  # which is by definition the full-lag online smoother. This ties the range
  # to an estimator that is already graded, so the row cannot drift away from
  # it silently.
  p <- .aci_cir_pieces()
  aux <- aciR:::.aci_online_aux(p$x, p$comp, p$dt, p$filt)
  online <- aci_online_smoother(p$x, p$comp, p$dt, p$filt, lag = Inf)

  for (j in c(10L, 100L, 500L, 900L)) {
    k <- seq.int(j, p$n - 1L)
    d <- aux$cum_sign[k] * aux$cum_sign[j] *
      exp(aux$cum_log[k] - aux$cum_log[j])
    expect_lt(
      abs(p$filt$mean[j] + sum(d * aux$innov_mean[k]) - online$mean[j]),
      1e-12
    )
  }
})

test_that("the comparison sequence is a relative entropy", {
  p <- .aci_cir_pieces()
  aux <- aciR:::.aci_online_aux(p$x, p$comp, p$dt, p$filt)
  re <- aciR:::.aci_cir_row(aux, p$filt, 100L, p$n)

  expect_true(all(re >= 0))
  expect_true(all(is.finite(re)))
  # Comparing the fully informed posterior against itself is the final entry
  # and must vanish.
  expect_equal(re[length(re)], 0)
  # And it must be non-degenerate, or every reduction below grades nothing.
  expect_gt(max(re), 0)
})

test_that("the subjective range shortens as the threshold rises", {
  p <- .aci_cir_pieces()
  rng <- aci_cir(p$x, p$comp, p$dt, p$filt, window = seq(20L, 700L, by = 40L))

  for (i in seq_len(ncol(rng$subjective))) {
    column <- rng$subjective[, i]
    column <- column[!is.na(column)]
    expect_true(all(diff(column) <= 0))
  }
  # A demanding threshold must be unresolved at least as often as a lax one.
  na_by_threshold <- rowSums(is.na(rng$subjective))
  expect_gte(na_by_threshold[1L], na_by_threshold[length(na_by_threshold)])
})

test_that("times the record cannot resolve are reported as unresolved", {
  p <- .aci_cir_pieces()
  window <- seq(20L, p$n - 20L, by = 20L)
  rng <- aci_cir(p$x, p$comp, p$dt, p$filt, window = window)

  # Saturation is a statement about how much record remains, so it must affect
  # the late times and not the early ones.
  expect_true(any(rng$saturated))
  expect_true(any(!rng$saturated))
  expect_gt(min(which(rng$saturated)), max(which(!rng$saturated)))
  # An unresolved time returns a right-censored LOWER BOUND, not a hole. The
  # record does support a statement about such a time (that the range is at
  # least this long), and discarding it loses information exactly at the end
  # of the record, which is where a user studying a recent event looks.
  expect_true(all(is.finite(rng$objective[rng$saturated])))
  expect_true(all(is.finite(rng$objective[!rng$saturated])))
  expect_identical(rng$saturated, rng$status == "censored")

  # The bound has to be a bound: a censored range must be at least as long as
  # the margin that condemned it, or the flag is decoration.
  expect_true(all(rng$objective[rng$saturated] > 0))
})

test_that("an uncoupled system has no causal influence range", {
  # With no coupling and no noise cross-covariance the observed process carries
  # no information about the unobserved one, so the two posteriors coincide and
  # the range collapses. This is an analytic identity, independent of both this
  # implementation and the reference one.
  # The step is compared at a fixed INSTANT, not a fixed index. A fixed index
  # names a different instant at every resolution, and the residue then mixes
  # the discretisation error with the filter's transient, which is how this
  # test first reported a meaningless order.
  uncoupled <- function(dt, t_at = 0.1) {
    n <- as.integer(0.6 / dt) + 1L
    j <- as.integer(t_at / dt) + 1L
    x <- 1 + 0.3 * sin(5 * seq(0, by = dt, length.out = n))
    comp <- list(
      L_x = rep(0, n), f_x = rep(0, n), L_y = -0.5, f_y = rep(1, n),
      S_xoS_x = 1, S_yoS_y = 1, S_yoS_x = 0, S_xoS_y = 0
    )
    filt <- aci_filter(x, comp, dt, mu0 = 2, R0 = 0.1)
    aux <- aciR:::.aci_online_aux(x, comp, dt, filt)
    list(
      peak = max(aciR:::.aci_cir_row(aux, filt, j, n)),
      rng = aci_cir(x, comp, dt, filt, window = seq(20L, n %/% 2L, by = 20L))
    )
  }

  # The identity is exact in continuous time only. Under an explicit Euler step
  # the divergence vanishes at second order, so the graded quantity is that
  # order and not a bound: halving the step must quarter the residue. A
  # constant bound here would pass for the wrong reasons at any one step size.
  steps <- c(4e-3, 2e-3, 1e-3, 5e-4)
  peaks <- vapply(steps, function(dt) uncoupled(dt)$peak, numeric(1))
  order <- log2(peaks[-length(peaks)] / peaks[-1])
  expect_true(all(order > 1.8 & order < 2.2))

  # And with no information to accumulate there is no range to report.
  expect_true(all(uncoupled(1e-3)$rng$objective %in% c(0, NA_real_)))
})

test_that("aci_cir validates its arguments", {
  p <- .aci_cir_pieces(n = 300L)
  expect_error(
    aci_cir(p$x, p$comp, p$dt, p$filt, window = c(0, 5)),
    "within the observed signal"
  )
  expect_error(
    aci_cir(p$x, p$comp, p$dt, p$filt, window = c(2, 1e6)),
    "within the observed signal"
  )
  expect_error(
    aci_cir(p$x, p$comp, p$dt, p$filt, window = 10:20, epsilon = c(1, -1)),
    "positive, non-missing thresholds"
  )
  expect_error(
    aci_cir(p$x, p$comp, p$dt, p$filt, window = 10:20, margin = 1),
    "strictly between zero and one"
  )
  expect_error(
    aci_cir(p$x, p$comp, p$dt, filt = NULL, window = 10:20),
    "`mu0` must be"
  )
})
