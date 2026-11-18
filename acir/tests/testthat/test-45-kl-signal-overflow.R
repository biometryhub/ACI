test_that("scalar signal avoids intermediate overflow and respects its limit", {
  mu_p <- c(0, 0, -1e308, 0, 0, 0, 0)
  mu_q <- c(1.5e154, -1.5e154, 1e308, 2^-25, 2^512,
            as.numeric("0x1.255c5eb8e2014p+467"),
            as.numeric("0x1.ebe9a35113493p+730"))
  variance <- c(1, 1, 1.5e308, 2^-1074, 0.5,
                as.numeric("0x1.502c92276d45ap-91"),
                as.numeric("0x1.d89d06d9d8e02p+436"))
  # Rounded from high-precision evaluation on these exact binary inputs.
  # The last two cases straddle the overflow rounding midpoint. Merely
  # changing the order of squaring and division misclassifies both.
  expected <- c(1.1250000000000002e308, 1.1250000000000002e308,
                1.3333333333333333e308, 8.98846567431158e307,
                Inf, .Machine$double.xmax, Inf)
  metric <- function(indices, decompose = TRUE) {
    n <- length(indices)
    p <- new_da_path(seq_len(n), matrix(mu_p[indices], n, 1),
                     array(variance[indices], c(1, 1, n)), "smoother")
    q <- new_da_path(seq_len(n), matrix(mu_q[indices], n, 1),
                     array(variance[indices], c(1, 1, n)), "filter")
    aci_metric(p, q, decompose = decompose)
  }
  indices <- seq_along(mu_p)
  got <- metric(indices)
  finite <- is.finite(expected)
  expect_true(all(is.infinite(0.5 * ((mu_q - mu_p) / sqrt(variance))^2)))
  expect_identical(is.finite(got$signal), finite)
  expect_true(all(abs(got$signal[finite] / expected[finite] - 1) < 1e-15))
  expect_identical(got$signal[!finite], expected[!finite])
  expect_identical(got$dispersion, numeric(length(indices)))
  expect_identical(got$total, got$signal)
  expect_identical(metric(indices, FALSE)$total, got$total)
  expect_identical(metric(rev(indices))$signal, rev(got$signal))
  expect_identical(vapply(indices, function(i) metric(i)$signal, numeric(1)),
                   got$signal)
})

test_that("finite scalar signal arithmetic is preserved in mixed calls", {
  mu_p <- c(0, 0, 0, 0, 0, 0)
  mu_q <- c(0, 1, -1, 1e150, 1e-160, 1.5e154)
  variance <- c(1, 2, 3, 1, 1e-320, 1)
  n <- length(mu_q)
  p <- new_da_path(seq_len(n), matrix(mu_p, n, 1),
                   array(variance, c(1, 1, n)), "smoother")
  q <- new_da_path(seq_len(n), matrix(mu_q, n, 1),
                   array(variance, c(1, 1, n)), "filter")
  previous <- 0.5 * ((mu_q - mu_p) / sqrt(variance))^2
  got <- aci_metric(p, q)
  finite <- is.finite(previous)
  expect_identical(got$signal[finite], previous[finite])
  expect_identical(got$dispersion, numeric(n))
  expect_true(is.finite(got$signal[n]))
})
