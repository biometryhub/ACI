# -- the online smoother and causal influence range, vector case ---------------
#
# One thing does NOT carry over from the scalar implementation, and it is the
# thing that made the scalar one fast. In one dimension the ordered product of
# the per-step auxiliary matrices reduces to a difference of cumulative
# logarithms, so any range is recoverable in constant time. Matrices do not
# commute and there is no such reduction. What survives is the reason the
# reduction was worth having: the products decay geometrically, so the
# accumulation is truncated rather than reconstructed.
#
# Three gradings, in increasing strength:
#
#   1. Lag zero is the filter, exactly, as in one dimension.
#   2. At one dimension the vector path must reproduce the SCALAR online
#      smoother and range at every lag -- and that scalar path is graded
#      against the authors' MATLAB, so this ties the new code to the package's
#      oldest oracle through an intermediary.
#   3. With a matrix cross-covariance it is graded against an independent
#      MATLAB transcription of equations (3.5) to (3.7) in FULL generality.
#      The reference implementation carries only the zero-cross-noise
#      specialisation, equation (3.8), so there is nothing upstream to grade
#      the general form against and this is not an authors' grounding.

.aci_mvos_fixture <- function(name) {
  path <- system.file("extdata", name, package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = sprintf("oracle fixture %s must ship in inst/extdata", name)
  )
  path
}

.aci_mvos_collapse <- function(n_keep = 1200L) {
  x <- read.csv(
    .aci_mvos_fixture("dyad_signal_x.csv"), header = FALSE
  )$V2[seq_len(n_keep)]
  p <- list(
    d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
    sigma_x = 0.5, sigma_y = 1
  )
  n <- length(x)
  list(
    x = x, dt = 0.001,
    scalar = aci_dyad_components(x, p),
    vector = list(
      L_x = array(2 * x, c(1L, 1L, n)),
      f_x = matrix(0.5 - 0.5 * x, 1L, n),
      L_y = matrix(-0.5, 1L, 1L),
      f_y = matrix(1 - 2 * x^2, 1L, n),
      S_xoS_x = matrix(0.25, 1L, 1L),
      S_yoS_y = matrix(1, 1L, 1L),
      S_yoS_x = matrix(0, 1L, 1L)
    ),
    X = matrix(x, 1L, n)
  )
}

test_that("a lag of zero is exactly the filter in the vector case", {
  p <- .aci_mvos_collapse(400L)
  filt <- aci_filter(p$X, p$vector, p$dt, mu0 = 2, R0 = matrix(0.1))
  online <- aci_online_smoother(p$X, p$vector, p$dt, filt, lag = 0)
  expect_identical(online$mean, filt$mean)
  expect_identical(online$cov, filt$cov)
})

test_that("the vector online smoother collapses onto the scalar one", {
  # At every lag, not merely at full lag: a lag off by one agrees at the
  # boundary, where both saturate at the end of the record, and nowhere else.
  # That is exactly the defect this test caught while it was being written.
  p <- .aci_mvos_collapse()
  fs <- aci_filter(p$x, p$scalar, p$dt, mu0 = 2, R0 = 0.1)
  fm <- aci_filter(p$X, p$vector, p$dt, mu0 = 2, R0 = matrix(0.1))

  for (lag in c(0, 1, 7, 50, 400, Inf)) {
    os <- aci_online_smoother(p$x, p$scalar, p$dt, fs, lag = lag)
    om <- aci_online_smoother(p$X, p$vector, p$dt, fm, lag = lag)
    expect_lt(max(abs(om$mean[1L, ] - os$mean)), 1e-12)
    expect_lt(max(abs(om$cov[1L, 1L, ] - os$cov)), 1e-12)
  }
})

test_that("the vector range collapses onto the scalar range", {
  p <- .aci_mvos_collapse()
  fs <- aci_filter(p$x, p$scalar, p$dt, mu0 = 2, R0 = 0.1)
  fm <- aci_filter(p$X, p$vector, p$dt, mu0 = 2, R0 = matrix(0.1))
  window <- seq(20L, 700L, by = 40L)
  epsilon <- c(1e-1, 1e-2, 1e-3)

  a <- aci_cir(
    p$x, p$scalar, p$dt, fs, window = window, margin = 1e-9,
    epsilon = epsilon
  )
  b <- aci_cir(
    p$X, p$vector, p$dt, fm, window = window, margin = 1e-9,
    epsilon = epsilon
  )
  expect_lt(max(abs(a$peak - b$peak)), 1e-12)
  expect_lt(max(abs(a$objective - b$objective)), 1e-12)
  expect_lt(max(abs(a$subjective - b$subjective)), 1e-12)
  # Non-degenerate, or the agreement above is between two sets of zeros.
  expect_gt(min(b$peak), 1e-3)
})

test_that("the vector online smoother reproduces the oracle with cross-noise", {
  signal <- read.csv(.aci_mvos_fixture("mv_signal.csv"))
  ref <- read.csv(.aci_mvos_fixture("mv_online_reference.csv"))
  # Built locally rather than borrowed from test-oracle-mv.R: helpers are not
  # reliably shared between test files, and a grade that depends on another
  # file's private function is a grade that can vanish when that file moves.
  pieces <- local({
    n <- nrow(signal)
    x <- rbind(signal$x1, signal$x2)
    s_mat <- matrix(
      c(0.60, 0.10, 0.25, 0.05, 0.20, 0.50, 0.10, 0.30,
        0.15, 0.05, 0.70, 0.10, 0.05, 0.10, 0.15, 0.55), 4L, 4L
    )
    s_x <- s_mat[1:2, ]
    s_y <- s_mat[3:4, ]
    L_x <- array(0, c(2L, 2L, n))
    L_y <- array(0, c(2L, 2L, n))
    f_x <- matrix(0, 2L, n)
    f_y <- matrix(0, 2L, n)
    for (j in seq_len(n)) {
      x1 <- x[1L, j]
      x2 <- x[2L, j]
      tt <- signal$t[j]
      L_x[, , j] <- matrix(
        c(0.8 + 0.3 * x1, 0.1 * x2, 0.2 * sin(tt), 0.6 - 0.2 * x1), 2L, 2L
      )
      L_y[, , j] <- matrix(
        c(-1.2 + 0.1 * x1, 0.2, 0.3, -0.9 - 0.1 * x2), 2L, 2L
      )
      f_x[, j] <- c(0.4 - 0.5 * x1, -0.3 * x2 + 0.2)
      f_y[, j] <- c(0.5 - 0.2 * x1^2, 0.1 - 0.15 * x2)
    }
    list(x = x, comp = list(
      L_x = L_x, L_y = L_y, f_x = f_x, f_y = f_y,
      S_xoS_x = s_x %*% t(s_x), S_yoS_y = s_y %*% t(s_y),
      S_yoS_x = s_y %*% t(s_x)
    ))
  })
  dt <- 0.002

  comp <- aciR:::.aci_check_components_mv(pieces$comp, pieces$x)
  filt <- aci_filter(
    pieces$x, pieces$comp, dt, mu0 = c(0.8, 0.2), R0 = 0.2 * diag(2L)
  )
  aux <- aciR:::.aci_online_aux_mv(pieces$x, comp, dt, filt)

  err_mean <- numeric(nrow(ref))
  err_cov <- numeric(nrow(ref))
  for (i in seq_len(nrow(ref))) {
    j <- ref$j[i]
    n <- ref$n[i]
    mu <- filt$mean[, j]
    r <- filt$cov[, , j]
    d <- diag(2L)
    if (n > j) {
      for (k in seq.int(j, n - 1L)) {
        mu <- mu + d %*% aux$innov_mean[, k]
        r <- r + d %*% aux$innov_cov[, , k] %*% t(d)
        r <- (r + t(r)) / 2
        d <- d %*% aux$E_j[, , k]
      }
    }
    err_mean[i] <- max(abs(as.numeric(mu) - c(ref$om1[i], ref$om2[i])))
    err_cov[i] <- max(
      abs(r[1L, 1L] - ref$oc11[i]), abs(r[1L, 2L] - ref$oc12[i]),
      abs(r[2L, 2L] - ref$oc22[i])
    )
  }

  expect_lt(max(err_mean), 1e-6)
  expect_lt(max(err_cov), 1e-6)
  # The fixture must span real lags, or it grades only the trivial diagonal.
  expect_gt(max(ref$n - ref$j), 1000L)
  expect_gt(length(unique(ref$n - ref$j)), 3L)
})
