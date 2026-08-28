## acir reserve file
## Origin: aci/tests/testthat/test-01-kl.R:23-31
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: projected_kl identity and validation assertions.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: I3 projected_kl identity; projected_kl's only in-package consumer is the excluded extremes family.

test_that("I3: projected_kl equals 1-D KL of projections", {
  set.seed(3); l <- 3
  R0 <- crossprod(matrix(rnorm(9), 3)) + diag(3); RE <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  m0 <- rnorm(3); mE <- rnorm(3); v <- rnorm(3); v <- v / sqrt(sum(v^2))
  expect_equal(projected_kl(v, m0, R0, mE, RE),
               unname(gaussian_kl(drop(v %*% mE), drop(t(v) %*% RE %*% v),
                                  drop(v %*% m0), drop(t(v) %*% R0 %*% v),
                                  decompose = FALSE)), tolerance = 1e-12)
})


## Reserve staging block
## Origin: aci 0.0.30 (git tree 97f6b124) tests/testthat/test-01-kl.R:51-55
## Category: extensions
## Note: projected_kl validation assertions from the mixed block at :44.

  expect_error(projected_kl(0, 0, matrix(1), 1, matrix(1)),
               class = "aci_error_dims")
  expect_error(projected_kl(c(1, 0), c(0, 0), diag(c(1, -1)),
                            c(1, 1), diag(2)),
               class = "aci_error_spd")
