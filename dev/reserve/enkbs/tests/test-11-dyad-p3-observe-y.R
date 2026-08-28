## acir reserve file
## Origin: aci/tests/testthat/test-11-model-parity.R:8-16
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: p3 preset and reverse observe='y' partition assertions.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: Reverse observe='y' partition and the p3 (EnKBS) preset assertions.

  reverse_ic <- model_dyad(observe = "y")$meta$ic_default
  expect_equal(reverse_ic$x0, 2)
  expect_equal(reverse_ic$y0, 1)
  expect_match(model_dyad(observe = "y")$meta$source_status,
               "package extension")
  expect_equal(model_dyad(variant = "p3", observe = "y")$meta$source_status,
               "paper + MATLAB checked (published EnKBS dyad experiment)")
  expect_match(model_dyad(variant = "p3")$meta$provenance,
               "EnKBS causal-inference")
