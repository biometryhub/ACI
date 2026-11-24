test_that(paste0("scalar dispersion remains finite when only the full ",
                 "ratio overflows"), {
  rp <- c(1e308, .Machine$double.xmax, 1e-15)
  rq <- c(0.5, 0.5, .Machine$double.xmin * .Machine$double.eps)
  n <- length(rp)
  p <- new_da_path(seq_len(n), matrix(0, n, 1),
                   array(rp, c(1, 1, n)), "smoother")
  q <- new_da_path(seq_len(n), matrix(0, n, 1),
                   array(rq, c(1, 1, n)), "filter")
  # Rounded from 800-digit evaluations on the exact represented variances.
  expected <- c(1e308, .Machine$double.xmax, 1.0120112665365532e308)
  got <- aci_metric(p, q)
  expect_true(all(is.infinite(rp / rq)))
  expect_true(all(is.finite(got$dispersion)))
  expect_true(all(abs(got$dispersion / expected - 1) < 1e-15))
  expect_identical(got$signal, numeric(n))
  expect_identical(got$total, got$dispersion)
  expect_identical(aci_metric(p, q, decompose = FALSE)$total, got$total)
})

test_that("unrepresentable scalar dispersion is infinite rather than NaN", {
  p <- new_da_path(1, matrix(0), array(1e308, c(1, 1, 1)), "smoother")
  q <- new_da_path(1, matrix(0), array(0.25, c(1, 1, 1)), "filter")
  got <- aci_metric(p, q)
  expect_identical(got$dispersion, Inf)
  expect_identical(got$total, Inf)
  expect_identical(got$signal, 0)
})

test_that("two infinite variances give a NaN dispersion, not an error", {
  # Exported functions refuse infinite moments, so the kernel is called
  # directly.
  got <- .gaussian_kl_scalar_kernel(c(0, 0), c(1, Inf), c(1, 1), c(2, Inf))
  one <- .gaussian_kl_scalar_kernel(0, 1, 1, 2)
  expect_identical(got$signal, c(one$signal, 0))
  expect_identical(got$dispersion[1], one$dispersion)
  expect_true(is.nan(got$dispersion[2]))
  expect_true(is.nan(got$total[2]))
  # The same pair beside subnormal variance ratios.
  got <- .gaussian_kl_scalar_kernel(c(0, 0, 0), c(Inf, 1e-300, 1e-300),
                                    c(1, 1, 1), c(Inf, 1e10, 1e10))
  rest <- .gaussian_kl_scalar_kernel(c(0, 0), c(1e-300, 1e-300), c(1, 1),
                                     c(1e10, 1e10))
  expect_true(is.nan(got$total[1]))
  expect_identical(got$total[2:3], rest$total)
})

test_that("the scalar kernel recycles arguments of different lengths", {
  recycled <- function(mu_p, R_p, mu_q, R_q) {
    n <- max(lengths(list(mu_p, R_p, mu_q, R_q)))
    .gaussian_kl_scalar_kernel(rep_len(mu_p, n), rep_len(R_p, n),
                               rep_len(mu_q, n), rep_len(R_q, n))
  }
  # One case for each branch that indexes single elements: the exact
  # overflow decision, signal recovery, a subnormal variance ratio and an
  # overflowing one. The first case's mean and variance put the signal of
  # mu_p = 0 exactly at the overflow midpoint, as in test-45.
  tie_variance <- as.numeric("0x1.948b0fcd6e9e0p+0")
  tie_mean <- as.numeric("0x1.c71c71c71c71cp+512")
  cases <- list(list(c(0, 2^300), tie_variance, tie_mean, tie_variance),
                list(c(0, 0), c(1, 1), c(1, 1.5e154), 1),
                list(c(0, 0), c(1e-300, 1e-300), c(0, 0), 1e10),
                list(c(0, 0), 1e308, c(0, 0), c(0.5, 0.5)))
  for (args in cases)
    expect_identical(do.call(.gaussian_kl_scalar_kernel, args),
                     do.call(recycled, args))
})
