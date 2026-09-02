## Conditioning of the scalar dispersion term.
##
## dispersion = 0.5 * (R_p/R_q - 1 - log(R_p/R_q)) is an O(delta^2) quantity
## built from O(delta) pieces, so how it is written decides how much of it
## survives near R_p == R_q -- the operating point of every converged ACI run.
##
## Reference: the Maclaurin series of 0.5 * (d - log1p(d)),
##   0.5 * (d^2/2 - d^3/3 + d^4/4 - d^5/5 + ...).
## At |d| <= 1e-4 the d^6 and d^7 terms are below the double-precision
## resolution of the sum: adding them changes no bit of the four-term value at
## any delta used here (verified, max relative change exactly 0). The four-term
## series is therefore the exact answer for these inputs, not an approximation.
##
## Measured relative error against that reference (scratch, R 4.5.2,
## aarch64-apple-darwin20), R_q = 1, R_p = 1 + delta:
##
##   delta   direct form   log1p form
##   1e-04   6.83e-09      3.48e-13
##   1e-06   4.47e-05      3.16e-11
##   1e-08   5.00e-01      2.60e-09
##   1e-09   5.00e-01      1.56e-08
##   1e-10   5.00e-01      8.43e-07
##   1e-11   5.00e-01      1.60e-05
##   1e-12   5.00e-01      1.29e-04
##   1e-13   4.99e-01      1.24e-03
##
## The direct form loses the whole term from delta = 1e-8 down (relative error
## 0.5 = the value is off by a factor of two, and it is only that well behaved
## because R_q = 1 makes log(Lq) vanish exactly). At R_q = 3.7 the same deltas
## give relative errors of 1.0, 4.44e+04 and 4.44e+08 -- wrong by orders of
## magnitude, not merely clamped.
##
## Tolerance. The log1p form inherits one rounding of log1p at magnitude delta
## into a result of magnitude delta^2/2, so its relative error is bounded by
## ~eps/delta. Measured, the ratio (relative error) / (eps/delta) stays in
## [0.07, 0.73] across the whole grid and both signs. The test uses 4 * eps /
## delta: at least a 5x margin on every measured point, which leaves room for a
## platform log1p that is a couple of ulps short of correctly rounded, while
## still being ~1/delta tighter than anything the direct form could pass.

kl_disp <- function(R_p, R_q) {
  acir:::.gaussian_kl_scalar_kernel(
    numeric(length(R_p)), R_p, numeric(length(R_p)), R_q
  )$dispersion
}

## Four-term series; see the header for why it is exact at these deltas.
disp_series <- function(d) {
  0.5 * (d^2 / 2 - d^3 / 3 + d^4 / 4 - d^5 / 5)
}

## Relative-error budget of the log1p form at a given delta.
disp_tol <- function(d) {
  pmax(4 * .Machine$double.eps / abs(d), 8 * .Machine$double.eps)
}

test_that("scalar dispersion survives the near-tie, R_q = 1", {
  R_q <- rep(1, 10L)
  R_p <- 1 + 10^-(4:13)
  ## Exact by Sterbenz: the stored delta, free of input representation error.
  d <- R_p - R_q
  expect_true(all(d > 0))

  got <- kl_disp(R_p, R_q)
  ref <- disp_series(d)
  rel <- abs(got - ref) / ref

  expect_true(all(is.finite(rel)))
  expect_true(all(rel <= disp_tol(d)))
  ## Nothing was clamped away: every value is still strictly positive.
  expect_true(all(got > 0))
})

test_that("scalar dispersion survives the near-tie, R_p < R_q", {
  ## The operative branch: the smoother variance sits below the filter
  ## variance, so a converged ACI run spends its time at delta slightly < 0.
  R_q <- rep(1, 10L)
  R_p <- 1 - 10^-(4:13)
  d <- R_p - R_q
  expect_true(all(d < 0))

  got <- kl_disp(R_p, R_q)
  ref <- disp_series(d)
  rel <- abs(got - ref) / ref

  expect_true(all(is.finite(rel)))
  expect_true(all(rel <= disp_tol(d)))
  expect_true(all(got > 0))
})

test_that("scalar dispersion near-tie accuracy is scale free", {
  ## R_q away from 1 is where the direct form failed hardest, because
  ## log(R_q) no longer vanishes exactly.
  for (R_qv in c(3.7, 1e-3, 5e2)) {
    R_q <- rep(R_qv, 10L)
    R_p <- R_qv * (1 + 10^-(4:13))
    d <- R_p / R_q - 1
    got <- kl_disp(R_p, R_q)
    ref <- disp_series(d)
    rel <- abs(got - ref) / ref
    expect_true(all(rel <= disp_tol(d)),
                info = sprintf("R_q = %g, max ratio %.3g", R_qv,
                               max(rel / disp_tol(d))))
  }
})

test_that("scalar dispersion still matches the direct form where both are sound", {
  ## Away from the tie the two are numerically interchangeable, so the change
  ## cannot have moved any well-conditioned value.
  set.seed(21)
  R_q <- exp(rnorm(500L))
  R_p <- R_q * exp(rnorm(500L, 0, 0.5))
  ## |delta| > 0.05 puts the direct form's own eps/delta^2 error near 1e-13,
  ## so a disagreement above that would be the new form's, not the old one's.
  keep <- abs(R_p / R_q - 1) > 0.05
  expect_gt(sum(keep), 400L)

  Lq <- sqrt(R_q); Lp <- sqrt(R_p); A <- Lp / Lq
  direct <- 0.5 * (A * A - 1 + 2 * log(Lq) - 2 * log(Lp))

  expect_equal(kl_disp(R_p, R_q)[keep], pmax(direct, 0)[keep],
               tolerance = 1e-11)
})

test_that("kernel contract is unchanged by the dispersion rewrite", {
  set.seed(22)
  n <- 200L
  mu_p <- rnorm(n); mu_q <- rnorm(n)
  R_q <- exp(rnorm(n)); R_p <- R_q * exp(rnorm(n, 0, 0.4))

  out <- acir:::.gaussian_kl_scalar_kernel(mu_p, R_p, mu_q, R_q)
  expect_named(out, c("total", "signal", "dispersion"))
  expect_identical(out$total, out$signal + out$dispersion)
  expect_true(all(out$signal >= 0), label = "signal clamped at zero")
  expect_true(all(out$dispersion >= 0), label = "dispersion clamped at zero")

  ## Signal is untouched: it is still 0.5 * (dmu)^2 / R_q.
  expect_identical(out$signal, pmax(0.5 * ((mu_q - mu_p) / sqrt(R_q))^2, 0))

  ## Identity of distributions gives an exact zero, both parts.
  id <- acir:::.gaussian_kl_scalar_kernel(mu_p, R_q, mu_p, R_q)
  expect_identical(id$total, rep(0, n))
  expect_identical(id$dispersion, rep(0, n))

  ## decompose = FALSE contract.
  bare <- acir:::.gaussian_kl_scalar_kernel(mu_p, R_p, mu_q, R_q,
                                            decompose = FALSE)
  expect_named(bare, "total")
  expect_identical(bare$total, out$total)

  ## Positive-definiteness aborts, reference covariance checked first.
  expect_error(acir:::.gaussian_kl_scalar_kernel(0, 1, 0, 0),
               class = "aci_error_spd")
  expect_error(acir:::.gaussian_kl_scalar_kernel(0, -1, 0, 1),
               class = "aci_error_spd")
  expect_error(acir:::.gaussian_kl_scalar_kernel(0, -1, 0, -1),
               class = "aci_error_spd")
})
