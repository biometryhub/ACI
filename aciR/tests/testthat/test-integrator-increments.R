# -- grading the integrator against the reference, by inverting the scheme -----
#
# aci_simulate() used to be ungraded, and the stated reason was that R and
# MATLAB draw normal variates by different algorithms, so a simulated path
# cannot coincide with the reference path even at a matching seed. That is
# true. It is also a statement about two random number generators, not about
# the integrator, and grading path equality under a shared seed would have
# graded the wrong thing.
#
# Euler-Maruyama is invertible, which makes the integrator itself gradeable.
# For the dyad model
#
#   dx = (- d_x x + gamma x y + F_x) dt + sigma_x dW_x
#   dy = (- d_y y - gamma x^2 + F_y) dt + sigma_y dW_y
#
# the variates the reference actually drew are recoverable from its own
# captured path by subtracting the drift and dividing by the noise
# coefficient. Feeding them back through aci_simulate(increments = ) must
# return the reference path itself -- not a path with the same statistics, the
# same path. That is an exact grade of the scheme, its coefficient assembly and
# its step ordering, and it is what these tests run.
#
# The inverse below is written from the two equations above, with the
# parameters read off the model object, rather than by calling anything the
# simulator uses. An oracle assembled from the code it grades tests
# self-consistency; this one does not share the integrator's source.
#
# -- the fixture ---------------------------------------------------------------
#
# The inversion needs BOTH components at consecutive steps of the integration
# grid, and neither packaged oracle fixture supplies that. inst/extdata's
# dyad_signal_x.csv carries the observed path only, and dyad_reference.csv
# carries x and y at every hundredth step; from x alone the recursion is one
# equation short at every step, and the unobserved path is not recoverable.
#
# fixtures/dyad_true_path_head.csv is therefore the first 3001 steps of both
# components of the published run, taken at full double precision from the
# parity harness's captured MATLAB workspace for the reduced (N = 3000) dyad
# capture, which runs the reference's own rng(333). The first test below
# checks its observed column against the packaged signal, so the unobserved
# column is tied to the published run rather than merely asserted to belong
# to it.
#
# -- the declared tolerance ----------------------------------------------------
#
# Declared from the arithmetic before any comparison was run, not read off the
# disagreement. Two quantities set it.
#
#   Local injection. The recovery re-associates the forward expression, so the
#   variate handed back is not the reference's original bit pattern: the
#   subtraction, the division and the re-multiplication each round once, and
#   the re-driven addition rounds once more. That last term dominates, at half
#   an ulp of the state; five ulps of the largest state the window visits
#   bounds the whole chain with a factor of seven to spare.
#
#   Amplification. A deviation obeys d[j + 1] = (I + dt J[j]) d[j] + e, with J
#   the Jacobian of the drift at the reference state -- the noise is additive,
#   so it does not enter J. An injection at one step is amplified by the
#   product of the largest singular values of (I + dt J) over every later step.
#
# Summing the injection over the steps, each weighted by its own amplification,
# gives 1.15e-10 over this window. The gate is that derived bound, recomputed
# by .aci_roundoff_bound() below so it tracks whichever path it is applied to,
# under a declared ceiling of 1e-9 that stops it inflating unnoticed.
#
# Measured afterwards, for the record and not as a gate: 5.1e-15 on the
# observed component and 1.8e-15 on the unobserved, some four orders of
# magnitude inside the bound. The bound is conservative because the
# singular-value product assumes every injection aligns with the locally most
# expansive direction at every later step, which no round-off sequence does.
# It stands as the gate because it was derived in advance.
#
# One property here is exact rather than bounded, and one is not. Supplying the
# variates a seeded call would have drawn reproduces that call BIT FOR BIT,
# because it is the same code integrating the same doubles; that is asserted
# with expect_identical(), and it is what shows the new argument leaves the
# existing path untouched. Recovering variates from a path and re-driving it is
# NOT bit-exact, and cannot be: the recovery divides by the noise coefficient
# and the integrator multiplies by it again, and fl(s * fl(v / s)) is not v for
# every v. It returns to within one or two ulps, which is what that test
# asserts.

.aci_dyad_true_path <- function() {
  path <- testthat::test_path("fixtures", "dyad_true_path_head.csv")
  testthat::expect_true(file.exists(path))
  read.csv(path)
}

# The inverse of one Euler-Maruyama step, written from the model equations.
# Returns the standard normal variates, in the shape aci_simulate() consumes.
.aci_recover_increments <- function(x, y, dt, p) {
  n <- length(x)
  from <- seq_len(n - 1L)
  drift_x <- -p$d_x * x[from] + p$gamma * x[from] * y[from] + p$F_x
  drift_y <- -p$d_y * y[from] - p$gamma * x[from]^2 + p$F_y
  list(
    dW_x = (x[-1L] - x[from] - drift_x * dt) / (p$sigma_x * sqrt(dt)),
    dW_y = (y[-1L] - y[from] - drift_y * dt) / (p$sigma_y * sqrt(dt))
  )
}

# The declared round-off gate: the per-step injection, summed over the steps
# and weighted by the amplification each one still has to pass through. See the
# derivation in the header.
.aci_roundoff_bound <- function(x, y, dt, p) {
  n <- length(x)
  from <- seq_len(n - 1L)
  injection <- 5 * .Machine$double.eps * max(abs(x), abs(y), 1)
  m11 <- 1 + dt * (-p$d_x + p$gamma * y[from])
  m12 <- dt * p$gamma * x[from]
  m21 <- dt * (-2 * p$gamma * x[from])
  m22 <- rep.int(1 - dt * p$d_y, n - 1L)
  # Largest singular value of a 2 x 2 matrix, in closed form.
  frobenius <- m11^2 + m12^2 + m21^2 + m22^2
  squared_det <- (m11 * m22 - m12 * m21)^2
  growth <- sqrt(
    (frobenius + sqrt(pmax(frobenius^2 - 4 * squared_det, 0))) / 2
  )
  # Tail products: an injection at step k is amplified by every step after it.
  tail_log <- rev(cumsum(rev(log(growth))))
  injection * sum(exp(c(tail_log[-1L], 0)))
}

test_that("the captured true path is the published dyad run", {
  ref <- .aci_dyad_true_path()
  expect_identical(nrow(ref), 3001L)
  expect_identical(names(ref), c("x", "y"))

  signal_csv <- system.file("extdata", "dyad_signal_x.csv", package = "aciR")
  expect_true(nzchar(signal_csv))
  sig <- read.csv(signal_csv, header = FALSE)
  expect_gte(nrow(sig), nrow(ref))

  # The packaged signal is decimal text at fifteen significant digits, so over
  # a window reaching 2.01 it can agree with the captured doubles to about
  # 1e-14 and no closer. 1e-13 is the gate that file format justifies.
  expect_lt(max(abs(ref$x - sig$V2[seq_len(nrow(ref))])), 1e-13)

  # The capture starts at the model's own initial state, exactly.
  model <- aci_dyad_model()
  expect_identical(ref$x[1L], model$x0)
  expect_identical(ref$y[1L], model$y0)
})

test_that("the recovered increments are the reference's own normal draws", {
  ref <- .aci_dyad_true_path()
  p <- aci_dyad_model()$parameters
  z <- .aci_recover_increments(ref$x, ref$y, 0.001, p)

  expect_length(z$dW_x, nrow(ref) - 1L)
  expect_length(z$dW_y, nrow(ref) - 1L)
  expect_true(all(is.finite(c(z$dW_x, z$dW_y))))

  # These gates check that the recovered quantity is a plausible normal sample,
  # which catches a wrong `sigma` or a mis-scaled `dt`. They do NOT check the
  # drift, and an earlier version of this comment claimed they did.
  #
  # Measured, not argued: recovering with `F_x`, `F_y`, `d_x` or `gamma`
  # deliberately wrong passes every gate below. A constant drift error `delta`
  # displaces each variate by only `delta * sqrt(dt) / sigma`, which at
  # `dt = 0.001` is 0.063 per unit for x and 0.032 for y, so the mean gate
  # tolerates forcing errors up to roughly 1.6 and 3.2 respectively.
  #
  # Nor would a wrong drift here still reproduce the path: only an error shared
  # by this recovery and the simulator cancels, and that shared error is what
  # `test-oracle-dyad.R` catches by grading the filter against the reference's
  # own output.
  expect_lt(abs(mean(z$dW_x)), 0.1)
  expect_lt(abs(mean(z$dW_y)), 0.1)
  expect_lt(abs(stats::sd(z$dW_x) - 1), 0.07)
  expect_lt(abs(stats::sd(z$dW_y) - 1), 0.07)
  expect_lt(max(abs(c(z$dW_x, z$dW_y))), 6)
})

test_that("the integrator reproduces the reference path from its increments", {
  ref <- .aci_dyad_true_path()
  model <- aci_dyad_model()
  p <- model$parameters
  dt <- 0.001

  z <- .aci_recover_increments(ref$x, ref$y, dt, p)
  driven <- aci_simulate(model, n = nrow(ref), dt = dt, increments = z)

  gate <- .aci_roundoff_bound(ref$x, ref$y, dt, p)
  # The derived bound must stay under the ceiling declared in the header, so a
  # future change to the fixture or the model cannot widen the grade in
  # silence.
  expect_lt(gate, 1e-9)
  expect_lt(max(abs(driven$x - ref$x)), gate)
  expect_lt(max(abs(driven$y - ref$y)), gate)
})

test_that("supplied increments reproduce a seeded path bit for bit", {
  model <- aci_dyad_model()
  n <- 2000L
  # aci_simulate() seeds, then draws the observed increments and the
  # unobserved ones in that order, so the same two calls outside reproduce the
  # exact doubles it integrated.
  set.seed(42)
  z <- list(dW_x = stats::rnorm(n - 1L), dW_y = stats::rnorm(n - 1L))

  seeded <- aci_simulate(model, n = n, seed = 42)
  driven <- aci_simulate(model, n = n, increments = z)

  expect_identical(driven$x, seeded$x)
  expect_identical(driven$y, seeded$y)
  expect_identical(driven$t, seeded$t)
})

test_that("increments recovered from a simulated path re-drive it", {
  model <- aci_dyad_model()
  p <- model$parameters
  dt <- 0.001
  sim <- aci_simulate(model, n = 2000L, dt = dt, seed = 333)

  z <- .aci_recover_increments(sim$x, sim$y, dt, p)
  again <- aci_simulate(model, n = 2000L, dt = dt, increments = z)

  gate <- .aci_roundoff_bound(sim$x, sim$y, dt, p)
  expect_lt(gate, 1e-9)
  expect_lt(max(abs(again$x - sim$x)), gate)
  expect_lt(max(abs(again$y - sim$y)), gate)
})

test_that("a driven path consumes none of the caller's random stream", {
  model <- aci_dyad_model()
  z <- list(dW_x = rep(0.5, 9L), dW_y = rep(-0.25, 9L))

  set.seed(11)
  before <- stats::runif(3L)
  set.seed(11)
  invisible(aci_simulate(model, n = 10L, increments = z))
  expect_identical(before, stats::runif(3L))
})

test_that("the increments contract is enforced rather than coerced", {
  model <- aci_dyad_model()
  ok <- list(dW_x = rep(0, 9L), dW_y = rep(0, 9L))

  expect_error(
    aci_simulate(model, n = 10L, seed = 1, increments = ok),
    "not both"
  )
  expect_error(
    aci_simulate(model, n = 10L, increments = c(1, 2)),
    "named list"
  )
  expect_error(
    aci_simulate(model, n = 10L, increments = list(dW_x = rep(0, 9L))),
    "`dW_y`"
  )
  expect_error(
    aci_simulate(model, n = 10L, increments = list()),
    "`dW_x`, `dW_y`"
  )
  # A scalar is not recycled to the run length, and a short vector is not
  # padded: the increments are data, not a default.
  expect_error(
    aci_simulate(model, n = 10L, increments = list(dW_x = 0, dW_y = 0)),
    "not recycled"
  )
  expect_error(
    aci_simulate(
      model, n = 10L,
      increments = list(dW_x = rep(0, 9L), dW_y = rep(0, 8L))
    ),
    "not recycled"
  )
  expect_error(
    aci_simulate(
      model, n = 10L,
      increments = list(dW_x = as.character(rep(0, 9L)), dW_y = rep(0, 9L))
    ),
    "plain numeric vector"
  )
  expect_error(
    aci_simulate(
      model, n = 10L,
      increments = list(dW_x = matrix(0, nrow = 3L, ncol = 3L),
                        dW_y = rep(0, 9L))
    ),
    "plain numeric vector"
  )
  expect_error(
    aci_simulate(
      model, n = 10L,
      increments = list(dW_x = c(rep(0, 8L), NA_real_), dW_y = rep(0, 9L))
    ),
    "complete and finite"
  )
  expect_error(
    aci_simulate(
      model, n = 10L,
      increments = list(dW_x = rep(0, 9L), dW_y = c(rep(0, 8L), Inf))
    ),
    "complete and finite"
  )

  # A rejected call must leave the caller's generator exactly as it found it.
  set.seed(5)
  before <- stats::runif(2L)
  set.seed(5)
  expect_error(aci_simulate(model, n = 10L, increments = list()))
  expect_identical(before, stats::runif(2L))
})

test_that("supplied increments compose with the Milstein scheme", {
  # Geometric Brownian motion, whose exact solution depends on the driving
  # path only through the sum of its increments, so the reproduced path can be
  # checked against the closed form of the very realisation supplied.
  a <- 0.6
  b <- 0.4
  model <- aci_cgns_model(
    L_x = 0, f_x = function(x) a * x, L_y = -1, f_y = 0,
    S_xoS_x = b^2, S_yoS_y = 1, x0 = 1, y0 = 0
  )
  n <- 4001L
  dt <- 1 / 4000
  set.seed(21)
  z <- list(dW_x = stats::rnorm(n - 1L), dW_y = stats::rnorm(n - 1L))

  driven <- aci_simulate(
    model, n = n, dt = dt, scheme = "milstein", increments = z,
    sigma_x = function(x) b * x, d_sigma_x = function(x) rep(b, length(x))
  )
  exact <- exp((a - b^2 / 2) * (n - 1L) * dt + b * sqrt(dt) * sum(z$dW_x))

  # Milstein is strongly order one, so at dt = 1/4000 the terminal value
  # tracks the exact solution of this realisation to a few parts in a
  # thousand. The point of the assertion is that the supplied path is the one
  # integrated, not that the scheme is exact.
  expect_lt(abs(driven$x[n] - exact) / exact, 0.005)
})
