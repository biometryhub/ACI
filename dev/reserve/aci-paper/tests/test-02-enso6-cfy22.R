## acir reserve file
## Origin: aci/tests/testthat/test-02-models.R:79-87
## Source package: aci 0.0.30, git tree 97f6b124
## Category: aci-paper
## Intended release: 0.1.x, TBD with the supervisor/collaborators
## Reason: cfy22 arm of the enso6 variant loop.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The enso6 block looped over variants c('cfy22','aci_code'); the mainline keeps the aci_code arm only.

test_that("enso6 golden variant and u-hidden partitions stay CGNS", {
  for (v in c("cfy22", "aci_code")) {
    m <- model_enso6(hidden = c("u", "hW", "tau"), variant = v)
    expect_s3_class(m, "cgns_model")
    expect_true(all(is.finite(simulate(m, seed = 12, T = 1, dt = 5e-3)$obs$x)))
  }
  expect_equal(model_enso6(variant = "aci_code")$meta$params$factor, 0.65)
  expect_error(model_enso6(hidden = "TC"), class = "aci_error_model_contract")
})
