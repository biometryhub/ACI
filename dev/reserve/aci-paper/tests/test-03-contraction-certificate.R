## acir reserve file
## Origin: aci/tests/testthat/test-03-engine.R:353-378
## Source package: aci 0.0.30, git tree 97f6b124
## Category: aci-paper
## Intended release: 0.1.x, TBD with the supervisor/collaborators
## Reason: Tests lt_contraction_certificate, which leaves the 0.1.0 surface.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: lt_contraction_certificate diagnostics (andreou2026smoother eqs. 3.18-3.19); the function is dropped from the 0.1.0 surface.

test_that("contraction certificate reports per-step E diagnostics", {
  m <- model_dyad()
  s <- simulate(m, seed = 1, T = 2, dt = 0.01)
  init <- list(mean = 0, cov = diag(1, 1))
  tb <- lag_table(m, as_obs(s), mode = "forward", init = init)
  cert <- lt_contraction_certificate(tb)
  expect_identical(nrow(cert), length(tb$t) - 1L)
  expect_true(all(is.finite(cert$lambda_min)))
  expect_true(all(is.finite(cert$enorm)))
  # scalar hidden state: operator norm and spectral radius coincide with |E|
  expect_identical(cert$enorm, cert$rho_E)
  expect_equal(attr(cert, "gamma"), max(cert$enorm), tolerance = 1e-15)
  expect_identical(attr(cert, "condition_318"), all(cert$lambda_min > 0))
  # independent plumbing check at one step: rebuild E from a fresh filter
  filt <- da_filter(m, as_obs(s), init = init)
  co <- eval_coefs(m, s$obs$t[5], s$obs$x[5, ])
  Rf <- filt$cov[, , 5]
  Rfi <- chol_solve(Rf, diag(1), "Rf")
  Gi <- chol_solve(co$gxx, diag(1), "gxx")
  Gx <- co$Lx + t(co$gyx) %*% Rfi
  Gy <- co$Ly + co$gyy %*% Rfi
  E <- diag(1) + (co$gyx %*% Gi %*% Gx - Gy) * 0.01
  expect_equal(cert$enorm[5], abs(E[1, 1]), tolerance = 1e-12)
  expect_equal(cert$lambda_min[5], (1 - E[1, 1]) / 0.01,
               tolerance = 1e-8)
})
