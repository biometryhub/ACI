################################################################################
## aci-cir-rows.R - vectorised forward divergence rows for a scalar hidden state
################################################################################


#' Scalar views of the per-interval forward primitives (internal)
#'
#' The compiled primitives keep every per-interval quantity as a matrix so that
#' one code path serves any hidden dimension. For a scalar hidden state the
#' row builder below wants plain vectors, read once per call rather than once
#' per cell.
#'
#' @param E List of `1 x 1` update matrices, one per interval; entries before
#'   the first computed interval may be `NULL`.
#' @param dmu `N x 1` matrix of one-lag mean increments.
#' @param dR List of `1 x 1` one-lag covariance increments, `NULL`-padded like
#'   `E`.
#' @returns A list with the vectors `E`, `dmu`, `dR`, `logE` (logarithm of
#'   `abs(E)`) and `sgnE` (sign of `E`).
#' @noRd
.cir_scalar_primitives <- function(E, dmu, dR) {
  take <- function(x) if (is.null(x)) NA_real_ else x[1L]
  E_v <- vapply(E, take, numeric(1L))
  list(
    E = E_v, dmu = as.numeric(dmu[, 1L]), dR = vapply(dR, take, numeric(1L)),
    logE = log(abs(E_v)), sgnE = sign(E_v)
  )
}


#' Cumulative sum restarted every `block` terms (internal)
#'
#' Each term's cumulative value is the sum of the complete blocks before its
#' own block, a sum of at most `length(x) / block` block totals, plus the
#' partial sum inside its block. No value is the difference of two long
#' running sums, which is what keeps the ordered products below accurate on a
#' record of tens of thousands of steps.
#'
#' @param x Numeric vector.
#' @param block Positive block length.
#' @returns Numeric vector of the same length as `x`.
#' @noRd
.cir_blocked_cumsum <- function(x, block) {
  n <- length(x)
  if (n <= block) return(cumsum(x))
  out <- numeric(n)
  base <- 0
  from <- 1L
  while (from <= n) {
    to <- min(from + block - 1L, n)
    part <- cumsum(x[from:to])
    out[from:to] <- base + part
    base <- base + part[to - from + 1L]
    from <- to + 1L
  }
  out
}


#' One anchor's forward divergence row as vector arithmetic (internal)
#'
#' For a scalar hidden state every fixed-lag posterior at anchor `j` is an
#' affine function of the per-interval primitives, so a whole row is a
#' cumulative sum rather than a walk over cells. Ported from aciR 0.2.3
#' (`.aci_cir_row()` in `R/aci-cir.R` and the cumulative logarithms of
#' `.aci_online_aux()` in `R/aci-online-smoother.R`), tagged `parents-final`.
#'
#' The ordered products of the update factors are formed as exponentials of
#' summed logarithms, never as a product of thousands of factors below one.
#' Unlike the parent, the logarithms are summed from the anchor outward in
#' blocks of `block` cells (`.cir_blocked_cumsum()`), so no cell is ever the
#' difference of two record-length cumulative sums. Differencing them errs by
#' at least the rounding of that sum, `eps * |sum(log E)|`, on any platform
#' (6e-13 at 100,000 steps of contracting factors on x86, whose `cumsum`
#' accumulates in extended precision), and by 1.7e-12 at 20,000 steps and
#' 2.8e-11 at 100,000 where it accumulates in double (arm64), against a 1e-12
#' gate; the blocked form stays below 3e-13 on the far cells and 2e-15 on the
#' cells near the anchor. The mean
#' and variance sums run in the same order as the cell-by-cell recursion they
#' replace, and the relative entropy of each cell is evaluated in the
#' operations of `.kl_fast()` for a scalar state (Cholesky factor `sqrt`,
#' the two forward solves as divisions, the same summation order), so a cell
#' is bit-identical to the recursion's given the same posterior.
#'
#' The adaptive freeze of `.compiled_forward_reduce()` and
#' `.lagtable_core_compiled()` is reproduced as a cut of the row by index. The
#' cells are formed in chunks that grow fourfold from `chunk`, so a row that
#' freezes early costs about its own length and a long row costs a handful of
#' vector passes; the freeze predicate is evaluated on each chunk and the row
#' stops at the first cell where it has held for `window` consecutive cells,
#' or at `max_lag`. The tail estimate recorded at the cut is the one the
#' recursion records.
#'
#' Cell `0` is the caller's lag-zero divergence. Cell `c >= 1` scores the
#' posterior informed through interval `k = j + c - 1` against the reference
#' smoother at the anchor.
#'
#' @param j Anchor index, `1 <= j <= N`.
#' @param sp Scalar primitives from `.cir_scalar_primitives()`.
#' @param T2,Ub Suffix estimates from the primitives, length `N + 1`.
#' @param N Number of intervals.
#' @param mu_f,R_f Filtered mean and variance at the anchor.
#' @param mu_s,R_s Reference smoother mean and (guarded) variance at the
#'   anchor.
#' @param diag0 The lag-zero divergence, cell `0` of the row.
#' @param lam Smallest eigenvalue of the guarded reference variance.
#' @param tol,window,max_lag Adaptive freeze controls, as in the recursions.
#' @param rec Covariance-policy recorder, its `j` already set to the anchor.
#' @param block Cells per numerical block of the cumulative logarithms.
#' @param chunk Cells in the first evaluation chunk; later chunks grow
#'   fourfold.
#' @returns A list with `row` (cells `0` to the cut, no padding), `lag_kept`
#'   (the index of the last cell kept), `tail` (the tail estimate at the cut,
#'   `0` for an uncut row) and `frozen`.
#' @noRd
.cir_scalar_row <- function(j, sp, T2, Ub, N, mu_f, R_f, mu_s, R_s, diag0,
                            lam, tol, window, max_lag, rec, block = 512L,
                            chunk = 2048L) {
  C <- N + 1L - j
  k_last <- if (is.finite(max_lag)) as.integer(min(N, j + max_lag - 1)) else N
  row <- numeric(C + 1L)
  row[1L] <- diag0
  cut <- C
  tail <- 0
  frozen <- FALSE
  carry_log <- 0
  carry_sgn <- 1
  mu0 <- mu_f
  rr0 <- R_f
  below <- 0L
  chp <- sqrt(R_s)
  log_chp <- log(chp)
  ## The recursion's tail bound, in its arithmetic order, for a scalar state
  ## (`sqrt(l)` is exactly one).
  bound <- function(Dn2, k) {
    1.5 * (Dn2 * T2[k] / (2 * lam) + Dn2 * Ub[k] / (2 * lam))
  }

  ## The cells k .. k1 of the row: ordered products from the carried total,
  ## then the running mean and variance sums continued from the carried
  ## values.
  form <- function(k, k1) {
    if (k1 == k) {
      lg <- carry_log
      sg <- carry_sgn
    } else {
      lg <- carry_log + c(0, .cir_blocked_cumsum(sp$logE[k:(k1 - 1L)], block))
      sg <- carry_sgn * c(1, cumprod(sp$sgnE[k:(k1 - 1L)]))
    }
    d <- sg * exp(lg)
    kk <- k:k1
    list(
      lg = lg, sg = sg, d = d,
      mu = mu0 + cumsum(d * sp$dmu[kk]),
      rr = rr0 + cumsum(d * d * sp$dR[kk])
    )
  }

  k <- j
  size <- chunk
  while (k <= k_last) {
    k1 <- min(k + size - 1L, k_last)
    ch <- form(k, k1)
    ## A variance that is not positive goes through the covariance policy one
    ## cell at a time, as the recursion does; the chunk is shortened so that
    ## every cell before it is formed in one piece and the offending cell is
    ## formed alone.
    rq <- ch$rr
    lim <- range(rq)
    if (!(lim[1L] > 0) || !is.finite(lim[2L])) {
      bad <- which(!is.finite(rq) | rq <= 0)
      k1 <- if (bad[1L] > 1L) k + bad[1L] - 2L else k
      ch <- form(k, k1)
      rq <- ch$rr
      if (bad[1L] == 1L)
        rq[1L] <- .cov_guard_chol(
          matrix(rq[1L], 1L, 1L), rec, "metric_reference"
        )$R[1L, 1L]
    }
    m <- k1 - k + 1L
    chq <- sqrt(rq)
    w <- (ch$mu - mu_s) / chq
    A <- chp / chq
    val <- pmax(0.5 * (w * w + A * A - 1) + log(chq) - log_chp, 0)
    cc <- seq.int(k - j + 1L, k1 - j + 1L)
    row[cc + 1L] <- val

    ## The freeze predicate on this chunk's cells of positive lag two and
    ## beyond (the first chunk starts at lag one, which is never tested),
    ## continuing the count of consecutive holds carried in.
    if (tol > 0 && k1 > j) {
      first <- if (k == j) 2L else 1L
      kt <- seq.int(k + first - 1L, k1)
      vt <- val[first:m]
      fut <- bound(ch$d[first:m]^2, kt)
      ok <- is.finite(fut) & (vt + fut) < tol
      n_ok <- length(ok)
      idx <- seq_len(n_ok)
      last_fail <- cummax(idx * !ok)
      run <- idx - last_fail
      run[last_fail == 0L] <- run[last_fail == 0L] + below
      h <- match(TRUE, run >= window)
      if (!is.na(h)) {
        cut <- kt[h] - j + 1L
        tail <- vt[h] + fut[h]
        frozen <- TRUE
        break
      }
      below <- run[n_ok]
    }

    mu0 <- ch$mu[m]
    rr0 <- ch$rr[m]
    carry_log <- ch$lg[m] + sp$logE[k1]
    carry_sgn <- ch$sg[m] * sp$sgnE[k1]
    k <- k1 + 1L
    size <- 4L * size
  }

  if (!frozen && is.finite(max_lag) && max_lag <= C) {
    ## Reached max_lag: the recursion freezes at that cell and records the
    ## tail estimate when it tested the cell, and the cell's own value when it
    ## did not (a lag of one is never tested).
    cut <- k_last - j + 1L
    frozen <- TRUE
    last <- row[cut + 1L]
    if (tol > 0 && cut >= 2L) {
      fut <- bound(ch$d[length(ch$d)]^2, k_last)
      tail <- if (is.finite(fut)) last + fut else last
    } else {
      tail <- last
    }
  }

  list(
    row = row[seq_len(cut + 1L)], lag_kept = cut, tail = tail, frozen = frozen
  )
}
