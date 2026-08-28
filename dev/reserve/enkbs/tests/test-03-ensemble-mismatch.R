## acir reserve file
## Origin: aci/tests/testthat/test-03-engine.R:106-111
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Mismatched-observation guard exercised through enkbf on model_l84.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Mismatched-observation guard exercised on model_l84 (E) and enkbf; the mainline block is retargeted to model_dyad and the closed-form route only.

test_that("closed-form and ensemble engines reject mismatched observations", {
  m <- model_l84()
  bad <- observed_trajectory(seq(0, 0.1, by = 0.01), matrix(0, 11, 1))
  expect_error(da_filter(m, bad), class = "aci_error_dims")
  expect_error(enkbf(m, bad, m = 5), class = "aci_error_dims")
})
