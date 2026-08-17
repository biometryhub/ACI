# Analytic identities for the vector-valued core -------------------------------
#
# These are the PRIMARY grounding for the matrix noise-cross-covariance terms,
# because no other grounding for them exists. Every scalar model in the
# reference implementation sets the cross-covariance to zero, and all five ENSO
# scripts say in as many words that the cross-interaction terms are absent, so
# the matrix-valued terms aciR exposes have no upstream counterpart at all.
#
# oracle/aci_oracle_mv.m grades them against a second independent
# implementation, which refutes a transcription error but is not an authors'
# reference. What follows depends on neither implementation: it is algebra on
# the governing equations, solved in closed form.

# The stationary filtered covariance of a scalar constant-coefficient CGNS is
# the positive root of the algebraic Riccati equation. Identical to the helper
# in test-identities.R; repeated rather than shared so that each file states
# the oracle it grades against.
.riccati_cov_mv <- function(L_x, L_y, S_xoS_x, S_yoS_y, S_yoS_x) {
  k <- L_y * S_xoS_x - S_yoS_x * L_x
  (k + sqrt(k^2 + L_x^2 * (S_yoS_y * S_xoS_x - S_yoS_x^2))) / L_x^2
}

test_that("the vector filter sits at the block-diagonal Riccati fixed point", {
  # Two scalar systems, each with a NON-ZERO noise cross-covariance, assembled
  # block-diagonally into one matrix system. Block-diagonal means the joint
  # Riccati equation decouples into the two scalar ones, whose roots are known
  # in closed form, so the matrix recursion has an exact target to hit, and
  # the cross-noise terms are exercised while it hits it.
  a <- list(L_x = 1.3, L_y = -0.7, S_xoS_x = 0.8, S_yoS_y = 1.1, S_yoS_x = 0.35)
  b <- list(L_x = 0.9, L_y = -1.1, S_xoS_x = 1.2, S_yoS_y = 0.6, S_yoS_x = -0.4)

  for (q in list(a, b)) {
    expect_gt(q$S_xoS_x * q$S_yoS_y - q$S_yoS_x^2, 0)
    expect_true(q$S_yoS_x != 0)
  }

  r_a <- do.call(.riccati_cov_mv, a)
  r_b <- do.call(.riccati_cov_mv, b)
  # Each root must annihilate its own Riccati residual: the oracle checks
  # itself before anything is graded against it.
  for (q in list(list(q = a, r = r_a), list(q = b, r = r_b))) {
    residual <- 2 * q$q$L_y * q$r + q$q$S_yoS_y -
      (q$q$S_yoS_x + q$r * q$q$L_x)^2 / q$q$S_xoS_x
    expect_lt(abs(residual), 1e-12)
    expect_gt(q$r, 0)
  }

  n <- 4000L
  dt <- 1e-3
  # A constant observed signal keeps the drift terms out of the covariance
  # recursion, which is what the Riccati equation describes.
  x <- matrix(0, nrow = 2L, ncol = n)
  comp <- list(
    L_x = diag(c(a$L_x, b$L_x)),
    f_x = c(0, 0),
    L_y = diag(c(a$L_y, b$L_y)),
    f_y = c(0, 0),
    S_xoS_x = diag(c(a$S_xoS_x, b$S_xoS_x)),
    S_yoS_y = diag(c(a$S_yoS_y, b$S_yoS_y)),
    S_yoS_x = diag(c(a$S_yoS_x, b$S_yoS_x))
  )
  r_star <- diag(c(r_a, r_b))
  filt <- aci_filter(x, comp, dt, mu0 = c(0, 0), R0 = r_star)

  # Started at the fixed point, the recursion must stay there: the fixed point
  # of the differential equation is exactly the fixed point of its Euler map.
  expect_lt(max(abs(filt$cov[, , n] - r_star)), 1e-10)
  expect_lt(max(abs(filt$cov[1L, 2L, ])), 1e-12)
})

test_that("the vector metric agrees with the scalar metric at one dimension", {
  # The multivariate relative entropy must reduce to the scalar one. This is an
  # identity between two independently written expressions, so it grades the
  # matrix formula rather than merely exercising it.
  set.seed(4)
  m_f <- 0.4
  m_s <- 0.9
  r_f <- 1.7
  r_s <- 0.6
  scalar <- 0.5 * (m_s - m_f)^2 / r_f +
    0.5 * (r_s / r_f - 1 - log(r_s / r_f))
  mv <- aciR:::.aci_metric_mv(
    list(mean = matrix(m_f, 1L, 1L), cov = array(r_f, c(1L, 1L, 1L))),
    list(mean = matrix(m_s, 1L, 1L), cov = array(r_s, c(1L, 1L, 1L)))
  )
  expect_equal(mv$value, scalar, tolerance = 1e-12)
})

test_that("the vector metric is right on a correlated pair", {
  # A regression test for a real defect. The trace term of the multivariate
  # relative entropy is tr(R_f^-1 R_s). Written as tr(c_f^-1 R_s c_f^-1) it
  # agrees for DIAGONAL covariances and disagrees otherwise, so no
  # block-diagonal or scalar test can distinguish the two. This one can: the
  # covariances below are genuinely correlated, and the expected value is
  # computed from the definition using solve() and determinant(), which share
  # no code with the Cholesky path under test.
  r_f <- matrix(c(1.4, 0.6, 0.6, 0.9), 2L, 2L)
  r_s <- matrix(c(0.7, -0.3, -0.3, 1.1), 2L, 2L)
  d <- c(0.5, -0.8)

  expected <- 0.5 * (
    as.numeric(t(d) %*% solve(r_f) %*% d) +
      sum(diag(solve(r_f) %*% r_s)) - 2 -
      log(det(r_s) / det(r_f))
  )
  got <- aciR:::.aci_metric_mv(
    list(mean = matrix(0, 2L, 1L), cov = array(r_f, c(2L, 2L, 1L))),
    list(mean = matrix(d, 2L, 1L), cov = array(r_s, c(2L, 2L, 1L)))
  )
  expect_equal(got$value, expected, tolerance = 1e-12)
  expect_gt(got$value, 0)

  # And the defect it guards against would have produced a different number:
  # if the two expressions agreed here the test would be grading nothing.
  c_f <- chol(r_f)
  wrong_trace <- sum(diag(
    backsolve(c_f, t(backsolve(c_f, t(r_s), transpose = TRUE)))
  ))
  expect_gt(abs(wrong_trace - sum(diag(solve(r_f) %*% r_s))), 1e-6)
})

test_that("the vector core rejects an inadmissible system", {
  x <- matrix(0, 2L, 10L)
  base <- list(
    L_x = diag(2), f_x = c(0, 0), L_y = -diag(2), f_y = c(0, 0),
    S_xoS_x = diag(2), S_yoS_y = diag(2), S_yoS_x = matrix(0, 2L, 2L)
  )

  bad_dim <- base
  bad_dim$L_x <- matrix(1, 3L, 2L)
  expect_error(aci_filter(x, bad_dim, 1e-3, c(0, 0), diag(2)), "3 by 2|2 by 2")

  singular <- base
  singular$S_xoS_x <- matrix(c(1, 1, 1, 1), 2L, 2L)
  expect_error(
    aci_filter(x, singular, 1e-3, c(0, 0), diag(2)),
    "symmetric positive definite"
  )

  asymmetric <- base
  asymmetric$S_yoS_y <- matrix(c(1, 0.3, -0.3, 1), 2L, 2L)
  expect_error(
    aci_filter(x, asymmetric, 1e-3, c(0, 0), diag(2)), "must be symmetric"
  )
})
