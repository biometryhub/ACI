# Grading the Milstein scheme --------------------------------------------------
#
# The oracle here is a theorem rather than a reference implementation. For a
# stochastic differential equation with a state-dependent diffusion, the
# Euler-Maruyama scheme converges strongly at order one half and the Milstein
# scheme at order one. Those rates are properties of the schemes, established
# independently of any code, so measuring them grades the implementation
# against something this package did not author.
#
# The test system is geometric Brownian motion, whose solution is known in
# closed form, driven by the SAME Wiener path the simulator uses, so the
# comparison is against the exact solution of the very realisation being
# integrated, not against a different draw from the same law.
#
# A rate is the right assertion here for the same reason it was for the online
# smoother: a scheme that is merely close at one step size proves nothing, and
# a wrongly stated correction term does not converge at a clean order.

.aci_gbm_pieces <- function(a = 0.6, b = 0.4, x0 = 1) {
  list(
    a = a, b = b, x0 = x0,
    model = aci_cgns_model(
      L_x = 0, f_x = function(x) a * x, L_y = -1, f_y = 0,
      S_xoS_x = b^2, S_yoS_y = 1, x0 = x0, y0 = 0
    ),
    sigma_x = function(x) b * x,
    d_sigma_x = function(x) rep(b, length(x))
  )
}

# Mean absolute deviation from the exact solution over many realisations: the
# strong error, which is what the convergence orders above are stated for.
.aci_strong_error <- function(p, dt, scheme, paths = 200L, t_end = 1) {
  n <- as.integer(t_end / dt) + 1L
  errs <- vapply(seq_len(paths), function(s) {
    sim <- aci_simulate(
      p$model, n = n, dt = dt, seed = s, scheme = scheme,
      sigma_x = p$sigma_x, d_sigma_x = p$d_sigma_x
    )
    # The simulator draws the observed increments first, so re-drawing them
    # under the same seed recovers the very path it integrated.
    set.seed(s)
    dw <- stats::rnorm(n - 1L)
    exact <- p$x0 * exp(
      (p$a - p$b^2 / 2) * (n - 1L) * dt + p$b * sqrt(dt) * sum(dw)
    )
    abs(sim$x[n] - exact)
  }, numeric(1))
  mean(errs)
}

test_that("Euler-Maruyama converges strongly at order one half", {
  p <- .aci_gbm_pieces()
  steps <- c(1 / 32, 1 / 64, 1 / 128, 1 / 256)
  err <- vapply(
    steps, .aci_strong_error, numeric(1), p = p, scheme = "euler_maruyama"
  )
  order <- mean(log2(err[-length(err)] / err[-1]))
  expect_gt(order, 0.3)
  expect_lt(order, 0.75)
})

test_that("Milstein converges strongly at order one", {
  p <- .aci_gbm_pieces()
  steps <- c(1 / 32, 1 / 64, 1 / 128, 1 / 256)
  err <- vapply(
    steps, .aci_strong_error, numeric(1), p = p, scheme = "milstein"
  )
  order <- mean(log2(err[-length(err)] / err[-1]))
  expect_gt(order, 0.8)
  expect_lt(order, 1.3)
})

test_that("Milstein beats Euler-Maruyama on the same paths", {
  # The orders above are the substantive claim; this is the consequence a user
  # would notice, and it must hold at every step size rather than on average.
  p <- .aci_gbm_pieces()
  for (dt in c(1 / 32, 1 / 128)) {
    euler <- .aci_strong_error(p, dt, "euler_maruyama")
    milstein <- .aci_strong_error(p, dt, "milstein")
    expect_lt(milstein, euler)
  }
})

test_that("the schemes coincide when the diffusion is constant", {
  # With a diffusion that does not vary with the state, the correction term is
  # identically zero and the two schemes must agree exactly, not approximately.
  p <- .aci_gbm_pieces()
  constant <- function(x) rep(0.4, length(x))
  zero <- function(x) rep(0, length(x))
  a <- aci_simulate(
    p$model, n = 500L, dt = 1e-3, seed = 7L, scheme = "euler_maruyama",
    sigma_x = constant, d_sigma_x = zero
  )
  b <- aci_simulate(
    p$model, n = 500L, dt = 1e-3, seed = 7L, scheme = "milstein",
    sigma_x = constant, d_sigma_x = zero
  )
  expect_identical(a$x, b$x)
})

test_that("the Milstein scheme refuses to guess its correction term", {
  p <- .aci_gbm_pieces()
  expect_error(
    aci_simulate(p$model, n = 100L, scheme = "milstein"),
    "Supply `sigma_x`"
  )
  expect_error(
    aci_simulate(
      p$model, n = 100L, scheme = "milstein", sigma_x = p$sigma_x
    ),
    "`d_sigma_x` must be supplied"
  )
  expect_error(
    aci_simulate(p$model, n = 100L, sigma_x = 0.4),
    "must be a function"
  )
})
