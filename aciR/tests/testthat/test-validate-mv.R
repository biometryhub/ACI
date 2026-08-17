# The vector-valued validator's guards, driven through the public entry point
# rather than by calling the internal checkers directly. What matters is not
# that a helper rejects a malformed argument, but that a user handing a
# malformed components list to `aci_filter()` is stopped with a message naming
# the component and saying what shape was wanted. A guard that fires only when
# called directly is not protecting anybody.
#
# Each case perturbs exactly one thing about a known-good list, so a failure
# here names the specific guard that stopped working.

.mv_valid <- function(n = 40L, n_x = 2L, n_y = 2L) {
  list(
    L_x = array(rep(diag(-1, n_x, n_y), n), c(n_x, n_y, n)),
    f_x = matrix(0, n_x, n),
    L_y = array(rep(diag(-1, n_y, n_y), n), c(n_y, n_y, n)),
    f_y = matrix(0, n_y, n),
    S_xoS_x = diag(0.1, n_x),
    S_yoS_y = diag(0.1, n_y),
    S_yoS_x = matrix(0, n_y, n_x),
    S_xoS_y = matrix(0, n_x, n_y)
  )
}

.mv_signal <- function(n = 40L, n_x = 2L) {
  set.seed(4L)
  matrix(stats::rnorm(n_x * n), n_x, n)
}

.mv_filter <- function(comp, x = .mv_signal(), n_y = 2L) {
  aci_filter(x, comp, dt = 0.01, mu0 = rep(0, n_y), R0 = diag(0.1, n_y))
}

test_that("a well-formed vector-valued components list runs", {
  f <- .mv_filter(.mv_valid())
  expect_identical(dim(f$mean), c(2L, 40L))
  expect_identical(dim(f$cov), c(2L, 2L, 40L))
  expect_true(all(is.finite(f$mean)))
})

test_that("a missing component is named up front", {
  comp <- .mv_valid()
  comp$S_yoS_y <- NULL
  expect_error(.mv_filter(comp), "missing the component")
  expect_error(.mv_filter(comp), "S_yoS_y")

  # The required set is L_x, f_x, L_y, f_y, S_xoS_x, S_yoS_y and S_yoS_x.
  # S_xoS_y is not among them, so removing it is not an absence.
  optional <- .mv_valid()
  optional$S_xoS_y <- NULL
  expect_no_error(.mv_filter(optional))

  two <- .mv_valid()
  two$S_yoS_y <- NULL
  two$S_yoS_x <- NULL
  # The message pluralises, so that two absences do not read as one.
  expect_error(.mv_filter(two), "components")
})

test_that("the Grammian must be square in the observed dimension", {
  comp <- .mv_valid()
  comp$S_xoS_x <- diag(0.1, 3L)
  expect_error(.mv_filter(comp), "S_xoS_x")
  expect_error(.mv_filter(comp), "square in the observed dimension")
})

test_that("a coefficient array is rejected for shape and for non-finiteness", {
  wrong_shape <- .mv_valid()
  wrong_shape$L_x <- array(0, c(5L, 5L, 40L))
  expect_error(.mv_filter(wrong_shape), "L_x")
  expect_error(.mv_filter(wrong_shape), "time in columns|numeric matrix")

  # Shape is right, contents are not. This is a separate guard, because a
  # correctly shaped array of NA passes every dimension test.
  not_finite <- .mv_valid()
  not_finite$L_x[1L, 1L, 1L] <- NA_real_
  expect_error(.mv_filter(not_finite), "L_x.*finite throughout")
})

test_that("a drift term is rejected for shape and for non-finiteness", {
  wrong_shape <- .mv_valid()
  wrong_shape$f_x <- matrix(0, 7L, 40L)
  expect_error(.mv_filter(wrong_shape), "f_x")

  not_finite <- .mv_valid()
  not_finite$f_x[1L, 1L] <- Inf
  expect_error(.mv_filter(not_finite), "f_x.*finite throughout")
})

test_that("a supplied inverse Grammian is checked for shape", {
  comp <- .mv_valid()
  # The conditional construction supplies this deliberately, and it need not be
  # the inverse of anything, so the only thing that can be checked is its shape.
  comp$S_xoS_x_inv <- diag(1, 4L)
  expect_error(.mv_filter(comp), "S_xoS_x_inv")
})

test_that("the signal must be a complete finite matrix, two steps or more", {
  comp <- .mv_valid()

  # One column is not a record. A single observation cannot be smoothed.
  expect_error(
    aci_filter(matrix(1, 2L, 1L), comp, dt = 0.01,
               mu0 = rep(0, 2L), R0 = diag(0.1, 2L)),
    "at least two columns"
  )

  bad <- .mv_signal()
  bad[1L, 5L] <- NA_real_
  expect_error(.mv_filter(comp, x = bad), "complete and finite")
  # The message says where, because a long record makes a bare rejection
  # useless for finding the offending value.
  expect_error(.mv_filter(comp, x = bad), "position")

  inf <- .mv_signal()
  inf[2L, 9L] <- Inf
  expect_error(.mv_filter(comp, x = inf), "complete and finite")
})

test_that("a one-row signal is promoted rather than rejected", {
  # A single observed component supplied as a bare vector is a matrix with one
  # row, and treating it as an error would make the scalar case awkward to
  # express in the vector schema.
  n <- 30L
  comp <- .mv_valid(n = n, n_x = 1L, n_y = 1L)
  x <- as.numeric(stats::rnorm(n))
  f <- aci_filter(x, comp, dt = 0.01, mu0 = 0, R0 = 0.1)
  expect_identical(dim(f$mean), c(1L, n))
})
