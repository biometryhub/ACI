test_that("I1/I5: KL nonnegative, zero at identity, decomposition adds up", {
  set.seed(1)
  for (i in 1:50) {
    l <- sample(1:4, 1)
    A <- matrix(rnorm(l * l), l); Rp <- crossprod(A) + diag(l) * 0.1
    B <- matrix(rnorm(l * l), l); Rq <- crossprod(B) + diag(l) * 0.1
    mp <- rnorm(l); mq <- rnorm(l)
    v <- gaussian_kl(mp, Rp, mq, Rq)
    expect_gte(v["total"], -1e-12)
    expect_equal(unname(v["total"]), unname(v["signal"] + v["dispersion"]), tolerance = 1e-12)
    expect_lt(abs(gaussian_kl(mp, Rp, mp, Rp, decompose = FALSE)), 1e-10)
  }
})
test_that("I2: 1-D closed form matches quadrature", {
  set.seed(2)
  for (i in 1:20) {
    mp <- rnorm(1); sp <- runif(1, .3, 2); mq <- rnorm(1); sq_ <- runif(1, .3, 2)
    f <- function(z) dnorm(z, mp, sp) * (dnorm(z, mp, sp, log = TRUE) - dnorm(z, mq, sq_, log = TRUE))
    num <- integrate(f, mp - 12 * sp, mp + 12 * sp, rel.tol = 1e-10)$value
    expect_equal(unname(gaussian_kl(mp, sp^2, mq, sq_^2, decompose = FALSE)), num, tolerance = 1e-6)
  }
})
test_that("I3: projected_kl equals 1-D KL of projections", {
  set.seed(3); l <- 3
  R0 <- crossprod(matrix(rnorm(9), 3)) + diag(3); RE <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  m0 <- rnorm(3); mE <- rnorm(3); v <- rnorm(3); v <- v / sqrt(sum(v^2))
  expect_equal(projected_kl(v, m0, R0, mE, RE),
               unname(gaussian_kl(drop(v %*% mE), drop(t(v) %*% RE %*% v),
                                  drop(v %*% m0), drop(t(v) %*% R0 %*% v),
                                  decompose = FALSE)), tolerance = 1e-12)
})
test_that("I4: coordinate-freeness under affine maps", {
  set.seed(4); l <- 3
  Rp <- crossprod(matrix(rnorm(9), 3)) + diag(3); Rq <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  mp <- rnorm(3); mq <- rnorm(3)
  A <- matrix(rnorm(9), 3); while (abs(det(A)) < .3) A <- matrix(rnorm(9), 3)
  b <- rnorm(3)
  v1 <- gaussian_kl(mp, Rp, mq, Rq, decompose = FALSE)
  v2 <- gaussian_kl(A %*% mp + b, A %*% Rp %*% t(A), A %*% mq + b, A %*% Rq %*% t(A),
                    decompose = FALSE)
  expect_equal(unname(v1), unname(v2), tolerance = 1e-9)
})

test_that("KL path honours decomposition flag and projection validates v", {
  p <- new_da_path(c(0, 1), matrix(c(0, 1), 2, 1),
                   array(c(1, 1), c(1, 1, 2)), "p")
  q <- new_da_path(c(0, 1), matrix(c(0, 0), 2, 1),
                   array(c(1, 1), c(1, 1, 2)), "q")
  expect_equal(names(gaussian_kl_path(p, q, decompose = FALSE)),
               c("t", "total"))
  expect_error(projected_kl(0, 0, matrix(1), 1, matrix(1)),
               class = "aci_error_dims")
  expect_error(projected_kl(c(1, 0), c(0, 0), diag(c(1, -1)),
                            c(1, 1), diag(2)),
               class = "aci_error_spd")
})

test_that("public KL estimators reject singular or asymmetric covariance inputs", {
  singular <- diag(c(1, 0))
  asymmetric <- matrix(c(1, 0, 0.25, 1), 2, 2)
  expect_error(gaussian_kl(c(0, 0), singular, c(0, 0), diag(2)),
               class = "aci_error_spd")
  expect_error(gaussian_kl(c(0, 0), diag(2), c(0, 0), asymmetric),
               class = "aci_error_spd")

  A <- cbind(1:4, 2 * (1:4))
  B <- cbind(1:5, (1:5)^2)
  expect_error(empirical_kl(A, B), class = "aci_error_ensemble_rank")
  expect_error(empirical_kl(matrix(1:6, 2, 3), matrix(1:9, 3, 3)),
               class = "aci_error_ensemble_rank")
})
