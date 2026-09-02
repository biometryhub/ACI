## acir reserve file
## Origin: aci/tests/testthat/test-01-kl.R:66-70
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: empirical_kl rank assertions.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: empirical_kl rank assertions from the mixed block at :58.

  A <- cbind(1:4, 2 * (1:4))
  B <- cbind(1:5, (1:5)^2)
  expect_error(empirical_kl(A, B), class = "aci_error_ensemble_rank")
  expect_error(empirical_kl(matrix(1:6, 2, 3), matrix(1:9, 3, 3)),
               class = "aci_error_ensemble_rank")
