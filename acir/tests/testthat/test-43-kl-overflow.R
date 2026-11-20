test_that("scalar dispersion remains finite when only the full ratio overflows", {
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
