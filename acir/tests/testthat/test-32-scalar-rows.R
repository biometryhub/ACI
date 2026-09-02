# Vectorised scalar forward rows (.cir_scalar_row) against the cell-by-cell
# recursion they replace, on the worst case the rows face: a long record of
# contracting update factors, near-equal variances, a one-cell freeze window.
# The last block pins the cancellation the blocked cumulative logarithms exist
# to avoid, measured on the same record.

.scalar_rows_record <- function(N, seed = 20260901L) {
  set.seed(seed)
  E <- runif(N, 0.90, 0.9995)
  dmu <- stats::rnorm(N, sd = 0.05)
  dR <- stats::runif(N, -1e-4, 1e-4)
  list(
    E = E, dmu = dmu, dR = dR,
    sp = .cir_scalar_primitives(
      lapply(E, function(e) matrix(e, 1L, 1L)), matrix(dmu, ncol = 1L),
      lapply(dR, function(r) matrix(r, 1L, 1L))
    )
  )
}

## The suffix estimates as .compiled_forward_primitives() forms them for l = 1.
.scalar_rows_suffix <- function(rec) {
  N <- length(rec$E)
  T2 <- Ub <- numeric(N + 1L)
  for (n in N:1) {
    eh2 <- max(1, abs(rec$E[n]))^2
    T2[n] <- if (n < N) rec$dmu[n + 1L]^2 + eh2 * T2[n + 1L] else 0
    Ub[n] <- if (n < N) abs(rec$dR[n + 1L]) + eh2 * Ub[n + 1L] else 0
    if (!is.finite(T2[n]) || T2[n] > 1e12) T2[n] <- Inf
    if (!is.finite(Ub[n]) || Ub[n] > 1e12) Ub[n] <- Inf
  }
  list(T2 = T2, Ub = Ub)
}

## The recursion .compiled_forward_reduce() runs for l = 1: sequential
## products, the relative entropy of .kl_fast() and the freeze rule, cell by
## cell. The reference posterior is the row's own end point, which is what
## the Theorem 3 smoother is at the anchor.
.scalar_rows_recursion <- function(j, rec, sfx, mu_f, R_f, tol, window,
                                   max_lag) {
  E <- rec$E; dmu <- rec$dmu; dR <- rec$dR
  N <- length(E)
  d <- c(1, cumprod(E[j:(N - 1L)]))
  mu_s <- mu_f + sum(d * dmu[j:N])
  R_s <- R_f + sum(d * d * dR[j:N])
  lam <- R_s
  ## .kl_fast() for a scalar state: Cholesky factors are square roots, the
  ## forward solves are divisions
  chp <- sqrt(R_s)
  kl <- function(mu, R) {
    chq <- sqrt(R)
    w <- (mu - mu_s) / chq
    A <- chp / chq
    max(0.5 * (w * w + A * A - 1) + log(chq) - log(chp), 0)
  }
  diag0 <- kl(mu_f, R_f)
  mu <- mu_f + dmu[j]
  R <- R_f + dR[j]
  D <- 1
  row <- c(diag0, kl(mu, R))
  lag_kept <- 1L
  below <- 0L
  tail <- 0
  if (j < N) for (n in seq.int(j + 1L, N)) {
    lag_length <- n + 1L - j
    if (is.finite(max_lag) && lag_length > max_lag) {
      lag_kept <- as.integer(max_lag)
      tail <- row[length(row)]
      break
    }
    D <- D * E[n - 1L]
    mu <- mu + D * dmu[n]
    R <- R + D * dR[n] * D
    Pval <- kl(mu, R)
    row <- c(row, Pval)
    lag_kept <- lag_length
    if (tol > 0) {
      fut <- 1.5 * (D^2 * sfx$T2[n] / (2 * lam) + D^2 * sfx$Ub[n] / (2 * lam))
      ok <- is.finite(fut) && (Pval + fut) < tol
      below <- if (ok) below + 1L else 0L
      if (below >= window || lag_length >= max_lag) {
        tail <- if (is.finite(fut)) Pval + fut else Pval
        break
      }
    } else if (lag_length >= max_lag) {
      tail <- Pval
      break
    }
  }
  list(row = row, lag_kept = lag_kept, tail = tail, mu_s = mu_s, R_s = R_s,
       lam = lam, diag0 = diag0)
}

.scalar_rows_builder <- function(j, rec, sfx, ref, mu_f, R_f, tol, window,
                                 max_lag) {
  .cir_scalar_row(
    j, rec$sp, sfx$T2, sfx$Ub, length(rec$E), mu_f, R_f, ref$mu_s, ref$R_s,
    ref$diag0, ref$lam, tol, window, max_lag, .aci_reg_for(NULL, NA_real_)
  )
}

.scalar_rows_gate <- 1e-12


test_that("blocked rows match the recursion on a 100,000-step record", {
  rec <- .scalar_rows_record(100000L)
  sfx <- .scalar_rows_suffix(rec)
  for (j in c(1L, 50000L, 95000L, 99999L, 100000L)) {
    ref <- .scalar_rows_recursion(j, rec, sfx, 0.3, 1, tol = 0,
                                  window = Inf, max_lag = Inf)
    got <- .scalar_rows_builder(j, rec, sfx, ref, 0.3, 1, tol = 0,
                                window = Inf, max_lag = Inf)
    expect_identical(got$lag_kept, ref$lag_kept)
    expect_identical(got$tail, 0)
    expect_false(got$frozen)
    expect_length(got$row, length(ref$row))
    ## cell by cell, within the gate of the row's peak. The objectives are
    ## not compared on a record this long: once the increments fall below
    ## half an ulp each route's posterior freezes at its own last value, and
    ## the constant ulp-sized cell offset that leaves against the smoother
    ## sums over tens of thousands of flat cells, for any two arithmetic
    ## routes alike. The 20,001-step test below compares them.
    expect_lt(max(abs(got$row - ref$row)),
              .scalar_rows_gate * max(abs(ref$row)))
  }
})


test_that("the freeze cut lands on the recursion's cell, tail included", {
  rec <- .scalar_rows_record(20001L)
  sfx <- .scalar_rows_suffix(rec)
  cases <- list(
    list(tol = 1e-8, window = 1L, max_lag = Inf),
    list(tol = 1e-8, window = 3L, max_lag = Inf),
    list(tol = 1e-5, window = 1L, max_lag = 500),
    list(tol = 1e-5, window = 3L, max_lag = 2),
    list(tol = 1e-5, window = 3L, max_lag = 1),
    list(tol = 0, window = Inf, max_lag = 700),
    list(tol = 0, window = Inf, max_lag = 1)
  )
  frozen_by_window <- 0L
  for (j in c(1L, 7000L, 19000L, 20000L)) for (cs in cases) {
    ref <- .scalar_rows_recursion(j, rec, sfx, 0.3, 1, cs$tol, cs$window,
                                  cs$max_lag)
    got <- .scalar_rows_builder(j, rec, sfx, ref, 0.3, 1, cs$tol, cs$window,
                                cs$max_lag)
    expect_identical(got$lag_kept, ref$lag_kept)
    expect_length(got$row, length(ref$row))
    expect_lt(max(abs(got$row - ref$row)),
              .scalar_rows_gate * max(abs(ref$row)))
    ## the tail estimate is the divergence at the cut plus the bound; the
    ## divergence there is a difference of nearly equal posteriors, so the
    ## tail is held on the row's scale, not its own
    expect_lt(abs(got$tail - ref$tail),
              .scalar_rows_gate * max(abs(ref$row)))
    ## both objectives of the reported range, on the zero-padded row
    if (length(ref$row) > 3L) for (method in c("exact", "l1_linf")) {
      tau <- vapply(list(got$row, ref$row), function(p)
        .fwd_lengths(p, 0.005, method)[["tau"]], numeric(1L))
      expect_equal(tau[1L], tau[2L], tolerance = .scalar_rows_gate)
    }
    if (got$frozen && cs$tol > 0 && got$lag_kept < cs$max_lag)
      frozen_by_window <- frozen_by_window + 1L
  }
  ## the window rule, not only max_lag, was exercised
  expect_gt(frozen_by_window, 0L)
})


test_that("a non-positive variance reaches the covariance policy by cell", {
  rec <- .scalar_rows_record(2000L)
  ## one negative increment large enough to take the variance below zero at
  ## a single cell, restored by the next
  rec$dR[1005L] <- -3
  rec$dR[1006L] <- 3.2
  rec$sp <- .cir_scalar_primitives(
    lapply(rec$E, function(e) matrix(e, 1L, 1L)), matrix(rec$dmu, ncol = 1L),
    lapply(rec$dR, function(r) matrix(r, 1L, 1L))
  )
  sfx <- .scalar_rows_suffix(rec)
  j <- 1000L
  ref <- .scalar_rows_recursion(j, rec, sfx, 0.3, 1, 0, Inf, Inf)
  strict <- .aci_reg_for("none", NA_real_)
  expect_error(
    .cir_scalar_row(j, rec$sp, sfx$T2, sfx$Ub, 2000L, 0.3, 1, ref$mu_s,
                    ref$R_s, ref$diag0, ref$lam, 0, Inf, Inf, strict),
    class = "aci_error_spd"
  )
  lenient <- .aci_reg_for("floor", NA_real_)
  got <- .cir_scalar_row(j, rec$sp, sfx$T2, sfx$Ub, 2000L, 0.3, 1, ref$mu_s,
                         ref$R_s, ref$diag0, ref$lam, 0, Inf, Inf, lenient)
  expect_length(got$row, length(ref$row))
  ## the policy fired once, at that cell, and every other cell agrees; the
  ## reference recursion here does not floor, so its own value there is NaN
  bad <- 1005L - j + 2L
  expect_true(is.nan(ref$row[bad]))
  expect_true(is.finite(got$row[bad]))
  expect_identical(.aci_reg_freeze(lenient)$n_events, 1L)
  expect_lt(max(abs(got$row[-bad] - ref$row[-bad])),
            .scalar_rows_gate * max(abs(ref$row[-bad])))
})


test_that("differenced cumulative logarithms cancel, blocked ones do not", {
  ## The ordered products d_k = prod_{i=j}^{k-1} E_i for cells close to a
  ## late anchor, three ways: the sequential product (the recursion's own
  ## arithmetic), the difference of a record-length cumulative logarithm
  ## (aciR's form), and the blocked form as .cir_scalar_row() reads it back
  ## through a unit increment at k with the reference posterior at zero
  ## (cell k then holds 0.5 * d_k^2).
  rec <- .scalar_rows_record(100000L)
  j <- 95000L
  N <- length(rec$E)
  sfx <- .scalar_rows_suffix(rec)
  cum <- c(0, cumsum(log(rec$E)))
  d_seq <- c(1, cumprod(rec$E[j:(N - 1L)]))
  ks <- j + 1:40
  d_naive <- exp(cum[ks] - cum[j])
  d_blocked <- vapply(ks, function(k) {
    one <- rec
    one$dmu[] <- 0
    one$dmu[k] <- 1
    one$dR[] <- 0
    one$sp <- .cir_scalar_primitives(
      lapply(one$E, function(e) matrix(e, 1L, 1L)),
      matrix(one$dmu, ncol = 1L),
      lapply(one$dR, function(r) matrix(r, 1L, 1L))
    )
    got <- .cir_scalar_row(j, one$sp, sfx$T2, sfx$Ub, N, 0, 1, 0, 1, 0, 1,
                           0, Inf, Inf, .aci_reg_for(NULL, NA_real_))
    sqrt(2 * got$row[k - j + 2L])
  }, numeric(1L))
  truth <- d_seq[ks - j + 1L]
  err_naive <- max(abs(d_naive - truth) / truth)
  err_blocked <- max(abs(d_blocked - truth) / truth)
  ## The differenced form errs by at least the rounding of the
  ## record-length sum itself, eps * |sum log E| (1.1e-12 on this record),
  ## on every platform; where cumsum accumulates in double rather than in
  ## extended precision (arm64) the accumulation adds to it (2.9e-12 here,
  ## 6e-13 on x86). The blocked form measures 2e-15.
  expect_gt(err_naive, 1e-13)
  expect_lt(err_blocked, 1e-14)
  expect_lt(err_blocked, err_naive / 10)
})
