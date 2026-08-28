## acir reserve file
## Origin: aci/tests/testthat/test-18-compiled-lag.R:136-136
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: one_lag storage-mode rows excised from test-18-compiled-lag.R.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The one_lag storage-mode row of the settings list.

    one_lag = list(tol = 1e-7, window = 2L, max_lag = Inf),


## Reserve staging block
## Origin: aci 0.0.30 (git tree 97f6b124) tests/testthat/test-18-compiled-lag.R:178-195
## Category: dead
## Note: The smoother_only mode was unreachable in aci 0.0.30 (.lag_table_compiled never passes it); the mode and this block are removed from the acir mainline.

test_that("compiled lag core retains internal smoother-only state", {
  ds <- .compiled_lag_matrix_setup()
  compiled_filter <- .cgns_filter_compiled(ds$bundle, ds$init)
  compiled_smoother <- .smoother_thmD1_compiled(
    ds$bundle, compiled_filter
  )
  compiled <- .lagtable_core_compiled(
    ds$bundle, compiled_filter, compiled_smoother,
    mode = "smoother_only", tol = 0, window = Inf, max_lag = Inf
  )
  expect_null(compiled$rows)
  expect_equal(compiled$diag,
               gaussian_kl_path(compiled_smoother, compiled_filter)$total,
               tolerance = 1e-11)
  expect_equal(compiled$smu, compiled_smoother$mean, tolerance = 1e-11)
  expect_equal(compiled$scov, compiled_smoother$cov, tolerance = 1e-11)
  expect_true(all(compiled$tailbnd == 0))
})
