# -- the truncation path, on a system fast enough to reach it ------------------
#
# `tol` stops the accumulated update product once it can no longer move the
# estimate. On every system this package ships the product never gets there:
# a per-step contraction of ~0.998 needs some 26,000 steps to reach 1e-18 and
# the records are shorter, so the branch is dead on all the real fixtures. A
# path that is never exercised is a path nobody has checked.
#
# The fixture below is synthetic on purpose. The identity under test is
# algebraic -- once the running product is negligible, filling the rest of the
# sequence forward must leave the fully informed posterior where a `tol = 0`
# walk leaves it -- and that identity does not care where the contraction came
# from. Waiting for a "realistic" fast contraction is how this path stays
# untested.
#
# It is an admissible CGNS driven through `.aci_online_aux_mv`, not a
# hand-built array of update matrices; testing a hand-built array would be
# testing a different function. It is deliberately not in the oracle manifest:
# this is a path test, not a grade against the authors.

.fast_contraction_system <- function(n = 400L, dt = 0.01) {
  set.seed(7L)
  x <- rbind(cumsum(stats::rnorm(n + 1L, sd = 0.05)),
             cumsum(stats::rnorm(n + 1L, sd = 0.05)))
  # A strong latent self-drift is what makes E contract fast: with a diagonal
  # L_y = -a and the filter covariance at its fixed point, E is about
  # 1 - a * dt per step.
  comp <- list(
    L_x = diag(c(0.5, 0.5)),
    L_y = diag(c(-30, -30)),
    f_x = matrix(0.02, 2L, n + 1L),
    f_y = matrix(0.05, 2L, n + 1L),
    S_xoS_x = diag(2L),
    S_yoS_y = diag(2L),
    S_yoS_x = matrix(0, 2L, 2L)
  )
  filt <- aci_filter(x, comp, dt, mu0 = c(0, 0), R0 = diag(2L) * 0.1)
  # The auxiliaries want the validated components: the observation-noise
  # inverse is derived there, not supplied. Calling the internal with a raw
  # list would test a state the public path never produces.
  list(x = x, comp = aciR:::.aci_check_components_mv(comp, x), dt = dt,
       filt = filt, n = n + 1L)
}

test_that("the fixture contracts once settled, and grows before it does", {
  # A precondition, asserted rather than assumed. If the system stops
  # contracting, every assertion below would pass vacuously against a branch
  # that never ran -- which is the failure this whole file exists to prevent.
  f <- .fast_contraction_system()
  aux <- aciR:::.aci_online_aux_mv(f$x, f$comp, f$dt, f$filt)
  radius <- vapply(seq_len(dim(aux$E_j)[3L]), function(j) {
    max(Mod(eigen(aux$E_j[, , j], only.values = TRUE)$values))
  }, numeric(1L))

  # Settled behaviour: a strong contraction, which is what lets `tol` fire.
  expect_lt(stats::median(radius), 0.85)
  expect_lt(max(radius[-seq_len(20L)]), 1)

  # And the part worth keeping. The paper bounds the spectral radius of each
  # factor below one, and that bound is a statement about the SETTLED filter,
  # not about every step. Here the first few factors exceed one, because the
  # filter covariance starts at R0 = 0.1 rather than at its fixed point, so
  # S_yoS_y R_f^{-1} is momentarily large enough to flip the sign of the drift
  # term in E.
  #
  # This is asserted, not tolerated. Any lag bound estimated from a contraction
  # rate measured over all steps would be misled by this transient, which is
  # the failure mode a max-entry test on the accumulated product cannot have:
  # the product is what contracts, and it is what the code tests.
  expect_gt(max(radius[seq_len(20L)]), 1)
})

test_that("`tol` truncates the online smoother on a fast-contracting system", {
  f <- .fast_contraction_system()
  cut <- aci_online_smoother(f$x, f$comp, f$dt, f$filt, lag = Inf,
                             tol = 1e-18)
  # The whole point: the loop stopped early. On every shipped fixture this is
  # equal to the record length and the branch is dead.
  expect_lt(cut$lag_effective, f$n - 1L)

  # `tol` must be positive by contract, so the reference walk uses the
  # smallest value that cannot fire inside this record: at a contraction of
  # ~0.7 the product needs some 2,000 steps to reach 1e-300 and the record is
  # 401 long.
  full <- aci_online_smoother(f$x, f$comp, f$dt, f$filt, lag = Inf,
                              tol = 1e-300)
  expect_identical(full$lag_effective, f$n - 1L)

  # Truncating must not move the answer.
  expect_equal(cut$mean, full$mean, tolerance = 1e-14)
  expect_equal(cut$cov, full$cov, tolerance = 1e-14)
})

test_that("the CIR row's fill-forward leaves the informed posterior alone", {
  f <- .fast_contraction_system()
  aux <- aciR:::.aci_online_aux_mv(f$x, f$comp, f$dt, f$filt)

  for (j in c(10L, 80L, 200L)) {
    cut <- aciR:::.aci_cir_row_mv(aux, f$filt, j, f$n, tol = 1e-18,
                                  horizon = f$n)
    walk <- aciR:::.aci_cir_row_mv(aux, f$filt, j, f$n, tol = 1e-300,
                                   horizon = f$n)
    expect_length(cut, length(walk))
    # Two things land together if the fill-forward is wrong: the target the row
    # is measured against moves, and the tail is written as exact zero. The
    # whole row is compared, not just its head, so either would show.
    expect_equal(cut, walk, tolerance = 1e-13)
    expect_equal(cut[length(cut)], walk[length(walk)], tolerance = 1e-13)
  }
})
