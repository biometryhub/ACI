## acir reserve file
## Origin: aci/tests/testthat/test-11-model-parity.R:223-226
## Source package: aci 0.0.30, git tree 97f6b124
## Category: aci-paper
## Intended release: 0.1.x, TBD with the supervisor/collaborators
## Reason: cfy22 initial-condition assertions.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: cfy22 initial-condition assertions from the ENSO aci_code block.

  cfy <- model_enso6(hidden = "u", variant = "cfy22")
  cfy_ic <- c(u = 0, hW = 0, TC = 0.1, TE = 0.1, tau = 0, I = 2)
  expect_equal(cfy$meta$ic_default$x0, cfy_ic[cfy$meta$vars$observed])
  expect_equal(cfy$meta$ic_default$y0, cfy_ic[cfy$meta$vars$hidden])
