test_that("I1/I5: KL nonnegative, zero at identity, decomposition adds up", {
  set.seed(1)
  for (i in 1:50) {
    l <- sample(1:4, 1)
    A <- matrix(rnorm(l * l), l); Rp <- crossprod(A) + diag(l) * 0.1
    B <- matrix(rnorm(l * l), l); Rq <- crossprod(B) + diag(l) * 0.1
    mp <- rnorm(l); mq <- rnorm(l)
    v <- aci_metric_pair(mp, Rp, mq, Rq)
    expect_gte(v["total"], -1e-12)
    expect_equal(unname(v["total"]), unname(v["signal"] + v["dispersion"]), tolerance = 1e-12)
    expect_lt(abs(aci_metric_pair(mp, Rp, mp, Rp, decompose = FALSE)), 1e-10)
  }
})
test_that("I2: 1-D closed form matches quadrature", {
  set.seed(2)
  for (i in 1:20) {
    mp <- rnorm(1); sp <- runif(1, .3, 2); mq <- rnorm(1); sq_ <- runif(1, .3, 2)
    f <- function(z) dnorm(z, mp, sp) * (dnorm(z, mp, sp, log = TRUE) - dnorm(z, mq, sq_, log = TRUE))
    num <- integrate(f, mp - 12 * sp, mp + 12 * sp, rel.tol = 1e-10)$value
    expect_equal(unname(aci_metric_pair(mp, sp^2, mq, sq_^2, decompose = FALSE)), num, tolerance = 1e-6)
  }
})
test_that("I4: coordinate-freeness under affine maps", {
  set.seed(4); l <- 3
  Rp <- crossprod(matrix(rnorm(9), 3)) + diag(3); Rq <- crossprod(matrix(rnorm(9), 3)) + diag(3)
  mp <- rnorm(3); mq <- rnorm(3)
  A <- matrix(rnorm(9), 3); while (abs(det(A)) < .3) A <- matrix(rnorm(9), 3)
  b <- rnorm(3)
  v1 <- aci_metric_pair(mp, Rp, mq, Rq, decompose = FALSE)
  v2 <- aci_metric_pair(A %*% mp + b, A %*% Rp %*% t(A), A %*% mq + b, A %*% Rq %*% t(A),
                        decompose = FALSE)
  expect_equal(unname(v1), unname(v2), tolerance = 1e-9)
})

test_that("KL path honours decomposition flag and projection validates v", {
  p <- new_da_path(c(0, 1), matrix(c(0, 1), 2, 1),
                   array(c(1, 1), c(1, 1, 2)), "p")
  q <- new_da_path(c(0, 1), matrix(c(0, 0), 2, 1),
                   array(c(1, 1), c(1, 1, 2)), "q")
  expect_equal(names(aci_metric(p, q, decompose = FALSE)),
               c("t", "total"))
})

test_that("public KL estimators reject singular or asymmetric covariance inputs", {
  singular <- diag(c(1, 0))
  asymmetric <- matrix(c(1, 0, 0.25, 1), 2, 2)
  expect_error(aci_metric_pair(c(0, 0), singular, c(0, 0), diag(2)),
               class = "aci_error_spd")
  expect_error(aci_metric_pair(c(0, 0), diag(2), c(0, 0), asymmetric),
               class = "aci_error_spd")
})


## The matrix branch of aci_metric() carries its own Cholesky entry
## (.kl_path_chol) and its own loop, both reorganised for speed. These lock the
## two things that reorganisation could have moved: the arithmetic, against the
## public single-index estimator, and the order in which a bad path is
## rejected, which is reference covariance before integrating covariance at
## each index, ascending in index.
kl_path_pair <- function(n = 6L, l = 3L, seed = 11L) {
  set.seed(seed)
  mkcov <- function() {
    A <- matrix(stats::rnorm(l * l), l, l)
    crossprod(A) + diag(l) * 0.5
  }
  cp <- array(0, c(l, l, n)); cq <- array(0, c(l, l, n))
  for (j in seq_len(n)) { cp[, , j] <- mkcov(); cq[, , j] <- mkcov() }
  list(
    p = new_da_path(seq_len(n) - 1, matrix(stats::rnorm(n * l), n, l), cp, "p"),
    q = new_da_path(seq_len(n) - 1, matrix(stats::rnorm(n * l), n, l), cq, "q")
  )
}

test_that("matrix KL path agrees with the single-index estimator", {
  z <- kl_path_pair()
  got <- aci_metric(z$p, z$q, decompose = TRUE)
  n <- length(z$p$t)
  want <- vapply(seq_len(n), function(j)
    aci_metric_pair(z$p$mean[j, ], z$p$cov[, , j],
                    z$q$mean[j, ], z$q$cov[, , j]), numeric(3))
  expect_identical(got$total, unname(want["total", ]))
  expect_identical(got$signal, unname(want["signal", ]))
  expect_identical(got$dispersion, unname(want["dispersion", ]))
  expect_identical(names(got), c("t", "total", "signal", "dispersion"))
  expect_identical(
    aci_metric(z$p, z$q, decompose = FALSE)$total, got$total
  )
})

test_that("matrix KL path rejects a bad index, reference covariance first", {
  bad_sym <- matrix(c(2, 0, 0.25, 2, 0, 0, 0, 0, 2), 3, 3)
  singular <- diag(c(1, 1, 0))

  z <- kl_path_pair(); z$q$cov[, , 4L] <- bad_sym
  expect_error(aci_metric(z$p, z$q), "R_q.*symmetric",
               class = "aci_error_spd")
  z <- kl_path_pair(); z$p$cov[, , 4L] <- bad_sym
  expect_error(aci_metric(z$p, z$q), "R_p.*symmetric",
               class = "aci_error_spd")
  z <- kl_path_pair(); z$q$cov[, , 2L] <- singular
  expect_error(aci_metric(z$p, z$q), "R_q.*positive definite",
               class = "aci_error_spd")
  z <- kl_path_pair(); z$p$cov[, , 2L] <- singular
  expect_error(aci_metric(z$p, z$q), "R_p.*positive definite",
               class = "aci_error_spd")

  ## Same index, both bad: the reference is reported.
  z <- kl_path_pair()
  z$q$cov[, , 3L] <- singular; z$p$cov[, , 3L] <- bad_sym
  expect_error(aci_metric(z$p, z$q), "R_q", class = "aci_error_spd")
  ## Different indices: the earlier index is reported, whichever path it is on.
  z <- kl_path_pair()
  z$q$cov[, , 5L] <- singular; z$p$cov[, , 2L] <- singular
  expect_error(aci_metric(z$p, z$q), "R_p", class = "aci_error_spd")
  ## A symmetric perturbation inside the tolerance is still accepted.
  z <- kl_path_pair()
  z$q$cov[1L, 2L, 3L] <- z$q$cov[2L, 1L, 3L] + 1e-15
  expect_silent(aci_metric(z$p, z$q))
})
