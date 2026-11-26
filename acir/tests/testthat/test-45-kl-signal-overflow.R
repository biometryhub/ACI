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

test_that(paste0("scalar signal overflow is decided exactly at the ",
                 "rounding midpoint"), {
  # 2001599834386887^2 = 222399981598543 * (2^54 - 1), so for mu_p = 0 the
  # first mean and variance give a signal equal to the overflow midpoint
  # 2^1024 - 2^970, which rounds to Inf. A small positive mu_p takes the exact
  # signal just below the midpoint, where it rounds to the largest double,
  # and a small negative one takes it above. The sixth and seventh cases move
  # mu_q instead, by an amount whose half underflows. The last two cases
  # reach the midpoint only through the square of the rounding error of the
  # scaled mean difference, and lie above it by about 5e-34 and 3e-34 of its
  # value. Expected values come from exact rational arithmetic.
  tie_variance <- as.numeric("0x1.948b0fcd6e9e0p+0")
  tie_mean <- as.numeric("0x1.c71c71c71c71cp+512")
  mu_p <- c(0, 2^300, -2^300, 2^-1074, -2^-1074, -tie_mean, -tie_mean,
            as.numeric("-0x1.d1fabb5cf87c4p+457"),
            as.numeric("-0x1.743c1946018ffp+458"))
  mu_q <- c(rep(tie_mean, 5), 2^-1074, -2^-1074,
            as.numeric("0x1.a0f0f006a73fep+512"),
            as.numeric("0x1.9bdbab49bc313p+512"))
  variance <- c(rep(tie_variance, 7), as.numeric("0x1.5387f76c468aep+0"),
                as.numeric("0x1.4b4d8a3ea284dp+0"))
  xmax <- .Machine$double.xmax
  expected <- c(Inf, xmax, Inf, xmax, Inf, Inf, xmax, Inf, Inf)
  n <- length(mu_p)
  p <- new_da_path(seq_len(n), matrix(mu_p, n, 1),
                   array(variance, c(1, 1, n)), "smoother")
  q <- new_da_path(seq_len(n), matrix(mu_q, n, 1),
                   array(variance, c(1, 1, n)), "filter")
  expect_identical(aci_metric(p, q)$signal, expected)
})
