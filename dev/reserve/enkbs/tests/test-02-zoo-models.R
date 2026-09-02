## acir reserve file
## Origin: aci/tests/testthat/test-02-models.R:50-58
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Cross-family Z6-lite block (F/E/P constructors); filed with the EnKBS majority, see DISPOSITIONS.md.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Z6-lite asserted model_tipping_triad (F), model_pathways (P), model_l84 (E), model_l96 (E) and model_topographic (P). The acir mainline block is rewritten for the three ACI_code constructors.

test_that("Z6-lite: zoo models construct and validate", {
  expect_s3_class(model_tipping_triad(0.1), "cgns_model")
  expect_s3_class(model_pathways(), "cgns_model")
  expect_s3_class(model_l84(), "cgns_model")
  expect_s3_class(model_l96(n = 12), "stochastic_model")
  expect_s3_class(model_enso6(), "cgns_model")            # (hW, tau) hidden split IS CGNS
  expect_error(model_enso6(hidden = c("I")), class = "aci_error_model_contract")
  expect_s3_class(model_topographic(), "stochastic_model")
})
