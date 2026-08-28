## acir reserve file
## Origin: aci/tests/testthat/test-13-fbcir-models.R
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Whole test file; drives model_multiscale_fbcir and model_topographic_layered_fbcir only.
## Verbatim copy from the aci 0.0.30 sources; not modified.

test_that("FBCIR multiscale joint coefficients reproduce Eq. (25)", {
  m <- model_multiscale_fbcir(partition = "joint")
  x <- c(x1 = 0.4, x2 = -0.7)
  co <- eval_coefs(m, 9, x)

  expect_equal(co$Lx, matrix(c(1.24, 0, 0, 0.1), 2, 2, byrow = TRUE),
               tolerance = 1e-14)
  expect_equal(co$fx, c(1.603666666666667, 5.48), tolerance = 1e-13)
  expect_equal(co$Ly, matrix(c(-5, 4, -4, -12), 2, 2, byrow = TRUE))
  expect_equal(co$fy, c(0.504, -0.93), tolerance = 1e-14)
  expect_equal(m$Sx1(9, x), diag(c(0.15, 0.3)))
  expect_equal(m$Sx2(9, x), matrix(c(1.52, 0, 0, 29 / 6),
                                          2, 2, byrow = TRUE),
               tolerance = 1e-14)
  expect_equal(m$Sy2(9, x), diag(c(1, 2) / sqrt(0.1)))
  # W_y is shared by x and y, so the cross Gram must not disappear.
  expect_equal(co$gyx, m$Sy2(9, x) %*% t(m$Sx2(9, x)), tolerance = 1e-14)
  expect_true(isTRUE(m$meta$correlated_noise))
  expect_true(isTRUE(m$meta$shared_noise_preserved))
})

test_that("FBCIR multiscale conditional y1 split preserves ordering and noise", {
  m <- model_multiscale_fbcir(partition = "y1")
  x <- c(x1 = 0.4, x2 = -0.7, y2 = -0.2)
  co <- eval_coefs(m, 9, x)

  expect_equal(co$Lx, matrix(c(1.24, 0, -4), 3, 1), tolerance = 1e-14)
  expect_equal(co$fx, c(1.603666666666667, 5.46, 1.47), tolerance = 1e-13)
  expect_equal(co$Ly, matrix(-5, 1, 1))
  expect_equal(co$fy, -0.296, tolerance = 1e-14)
  expect_equal(m$Sx1(9, x), matrix(c(
    0.15, 0, 0,
    0, 0.3, 29 / 6,
    0, 0, 2 / sqrt(0.1)
  ), 3, 3, byrow = TRUE), tolerance = 1e-14)
  expect_equal(m$Sx2(9, x), matrix(c(1.52, 0, 0), 3, 1), tolerance = 1e-14)
  expect_equal(m$Sy2(9, x), matrix(1 / sqrt(0.1), 1, 1))
  expect_equal(m$meta$target_obs_idx, 1L)
  expect_equal(m$meta$conditioning_obs_idx, c(2L, 3L))
  expect_true(isTRUE(m$meta$correlated_noise))
})

test_that("FBCIR multiscale conditional y2 split preserves ordering and noise", {
  m <- model_multiscale_fbcir(partition = "y2")
  x <- c(x1 = 0.4, x2 = -0.7, y1 = 0.3)
  co <- eval_coefs(m, 9, x)

  expect_equal(co$Lx, matrix(c(0, 0.1, 4), 3, 1), tolerance = 1e-14)
  expect_equal(co$fx, c(1.975666666666667, 5.48, -0.996), tolerance = 1e-13)
  expect_equal(co$Ly, matrix(-12, 1, 1))
  expect_equal(co$fy, -2.13, tolerance = 1e-14)
  expect_equal(m$Sx1(9, x), matrix(c(
    0.15, 0, 1.52,
    0, 0.3, 0,
    0, 0, 1 / sqrt(0.1)
  ), 3, 3, byrow = TRUE), tolerance = 1e-14)
  expect_equal(m$Sx2(9, x), matrix(c(0, 29 / 6, 0), 3, 1), tolerance = 1e-14)
  expect_equal(m$Sy2(9, x), matrix(2 / sqrt(0.1), 1, 1))
  expect_equal(m$meta$target_obs_idx, 2L)
  expect_equal(m$meta$conditioning_obs_idx, c(1L, 3L))
  expect_true(isTRUE(m$meta$correlated_noise))
})

test_that("all multiscale partitions encode the same four-state SDE", {
  t <- 2.75
  z <- c(x1 = 0.4, x2 = -0.7, y1 = 0.3, y2 = -0.2)
  joint <- model_multiscale_fbcir(partition = "joint")
  y1 <- model_multiscale_fbcir(partition = "y1")
  y2 <- model_multiscale_fbcir(partition = "y2")

  drift_joint <- c(joint$f(t, z[1:2], z[3:4]),
                   joint$g(t, z[1:2], z[3:4]))
  D_joint <- rbind(cbind(joint$Sx1(t, z[1:2]), joint$Sx2(t, z[1:2])),
                   cbind(joint$Sy1(t, z[1:2]), joint$Sy2(t, z[1:2])))

  x_y1 <- z[c("x1", "x2", "y2")]
  local_drift_y1 <- c(y1$f(t, x_y1, z["y1"]),
                      y1$g(t, x_y1, z["y1"]))
  local_D_y1 <- rbind(cbind(y1$Sx1(t, x_y1), y1$Sx2(t, x_y1)),
                      cbind(y1$Sy1(t, x_y1), y1$Sy2(t, x_y1)))
  full_order_y1 <- c(1, 2, 4, 3)

  x_y2 <- z[c("x1", "x2", "y1")]
  drift_y2 <- c(y2$f(t, x_y2, z["y2"]),
                y2$g(t, x_y2, z["y2"]))
  D_y2 <- rbind(cbind(y2$Sx1(t, x_y2), y2$Sx2(t, x_y2)),
                cbind(y2$Sy1(t, x_y2), y2$Sy2(t, x_y2)))

  expect_equal(local_drift_y1[full_order_y1], drift_joint, tolerance = 1e-14)
  expect_equal(drift_y2, drift_joint, tolerance = 1e-14)
  expect_equal((local_D_y1 %*% t(local_D_y1))[full_order_y1, full_order_y1],
               D_joint %*% t(D_joint), tolerance = 1e-14)
  expect_equal(D_y2 %*% t(D_y2), D_joint %*% t(D_joint), tolerance = 1e-14)
})

test_that("FBCIR constructors validate without advancing the caller RNG", {
  set.seed(2468)
  invisible(stats::runif(1))
  invisible(model_multiscale_fbcir(partition = "joint"))
  invisible(model_multiscale_fbcir(partition = "y1"))
  invisible(model_multiscale_fbcir(partition = "y2"))
  invisible(model_topographic_layered_fbcir())
  actual <- stats::runif(1)

  set.seed(2468)
  invisible(stats::runif(1))
  expected <- stats::runif(1)
  expect_equal(actual, expected)
})

test_that("layered topographic repository model is distinct and transcribed", {
  m <- model_topographic_layered_fbcir(target_mode = 2)
  x <- c(0.2, -0.4, 0.6, -0.8, 1.0, -1.2)
  co <- eval_coefs(m, 7, x)
  omega <- c(1, 0.5, 0.25)

  expect_equal(co$Lx, matrix(c(
    -0.4 - omega[1], -0.2,
    2 * -0.8 - omega[2] / 2, -2 * 0.6,
    3 * -1.2 - omega[3] / (sqrt(2) * 3^2), -3 * 1.0
  ), 6, 1), tolerance = 1e-14)
  expect_equal(co$fx, c(
    -(1 / 40) * 0.2 - 2 * -0.4,
    -(1 / 40) * -0.4 + 2 * 0.2,
    -(1 / 40) * 0.6 - 2 * -0.8 / 2,
    -(1 / 40) * -0.8 + 2 * 0.6 / 2,
    -(1 / 40) * 1.0 - 2 * -1.2 / 3,
    -(1 / 40) * -1.2 + 2 * 1.0 / 3
  ), tolerance = 1e-14)
  expect_equal(co$Ly, matrix(-1 / 40, 1, 1))
  expect_equal(co$fy,
               omega[1] * 0.2 + 2 * omega[2] * 0.6 +
                 omega[3] * 1.0 / sqrt(2), tolerance = 1e-14)
  expect_equal(diag(m$Sx1(7, x)), c(rep(1 / 40, 4), rep(1 / 120, 2)),
               tolerance = 1e-14)
  expect_equal(m$Sy2(7, x), matrix(1 / (20 * sqrt(2)), 1, 1))
  expect_equal(m$meta$target_obs_idx, 3:4)
  expect_equal(m$meta$provenance,
               "repository_case_study_named_in_andreou2026cir_code_availability")
  expect_false(m$meta$paper_scope)
})

test_that("FBCIR constructors reject invalid scientific parameters", {
  expect_error(model_multiscale_fbcir(epsilon = 0),
               class = "aci_error_model_contract")
  expect_error(model_multiscale_fbcir(params = list(I = diag(3))),
               class = "aci_error_model_contract")
  expect_error(model_multiscale_fbcir(params = list(gamma_1 = 0.5)),
               class = "aci_error_model_contract")
  expect_error(model_topographic_layered_fbcir(target_mode = 0),
               class = "aci_error_model_contract")
})
