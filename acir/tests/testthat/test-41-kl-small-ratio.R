test_that("scalar KL stays finite for small positive variance ratios", {
  rp <- c(0.49, 1e-8, 1e-12, 1e-16, 1e-20, 1e-100, 1e-300, 1e-300)
  rq <- c(rep(1, 7), 1e100)
  n <- length(rp)
  p <- new_da_path(seq_len(n), matrix(0, n, 1),
                   array(rp, c(1, 1, n)), "smoother")
  q <- new_da_path(seq_len(n), matrix(0, n, 1),
                   array(rq, c(1, 1, n)), "filter")
  expected <- vapply(seq_len(n), function(j)
    aci_metric_pair(0, matrix(rp[j]), 0, matrix(rq[j]))[["total"]],
    numeric(1))
  got <- aci_metric(p, q)
  expect_true(all(is.finite(got$total)))
  expect_equal(got$total, expected, tolerance = 1e-12)
  expect_identical(got$signal, numeric(n))
  expect_identical(got$total, got$dispersion)
  expect_identical(aci_metric(p, q, decompose = FALSE)$total, got$total)
  expect_equal(got$total[5], 0.5 * (1e-20 - 1 - log(1e-20)),
               tolerance = 1e-14)
  # The final positive input ratio underflows, but its logarithm is finite.
  expect_identical(rp[n] / rq[n], 0)
  expect_equal(got$total[n], 0.5 * (-1 + 400 * log(10)),
               tolerance = 1e-12)
})

test_that("small-ratio dispersion retains accuracy across variance scales", {
  # Binary fractions make the input ratio exact at every scale.
  rq <- 2^c(-900, 0, 900)
  rp <- rq / 4
  n <- length(rp)
  p <- new_da_path(seq_len(n), matrix(0, n, 1),
                   array(rp, c(1, 1, n)), "smoother")
  q <- new_da_path(seq_len(n), matrix(0, n, 1),
                   array(rq, c(1, 1, n)), "filter")
  expected <- 0.5 * (-0.75 + log(4))
  expect_true(all(abs(aci_metric(p, q)$dispersion - expected) < 1e-15))
})

test_that("subnormal positive ratios retain their logarithmic contribution", {
  p <- new_da_path(1, matrix(0), array(1e-300, c(1, 1, 1)), "smoother")
  q <- new_da_path(1, matrix(0), array(1e23, c(1, 1, 1)), "filter")
  # 800-digit reference evaluated on the exact binary input values.
  expected <- 371.36749251853837
  expect_equal(aci_metric(p, q)$dispersion, expected, tolerance = 1e-12)
})
