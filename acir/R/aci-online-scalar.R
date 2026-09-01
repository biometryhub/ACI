################################################################################
## aci-online-scalar.R - Theorem-3 auxiliaries for a scalar hidden state, as
## vector arithmetic over the record
################################################################################


#' Products of per-interval matrices with per-interval vectors (internal)
#'
#' `A` holds one `k x k` matrix per interval, `v` one `k`-vector per interval
#' (as a row). `.scalar_aux_matvec()` returns `A[, , n] %*% v[n, ]` as row
#' `n`; `.scalar_aux_vecmat()` returns `v[n, ] %*% A[, , n]`. For `k = 1` both
#' are the plain product, so the value is the one the matrix kernels form.
#'
#' @param A Array `k x k x n`.
#' @param v Matrix `n x k`.
#' @returns Matrix `n x k`.
#' @noRd
.scalar_aux_matvec <- function(A, v) {
  k <- dim(A)[1L]
  n <- dim(A)[3L]
  if (k == 1L) return(matrix(as.numeric(A) * v[, 1L], n, 1L))
  W <- array(rep(t(v), each = k), c(k, k, n))
  t(colSums(aperm(A * W, c(2L, 1L, 3L))))
}


#' @rdname dot-scalar_aux_matvec
#' @noRd
.scalar_aux_vecmat <- function(v, A) {
  k <- dim(A)[1L]
  n <- dim(A)[3L]
  if (k == 1L) return(matrix(v[, 1L] * as.numeric(A), n, 1L))
  W <- aperm(array(t(v), c(k, n, k)), c(1L, 3L, 2L))
  t(colSums(A * W))
}


#' Theorem-3 auxiliaries over the whole record for a scalar hidden state
#' (internal)
#'
#' For `l = 1` and any `k`, the per-interval quantities that
#' `.thmD1_aux_compiled()` and `.onelag_stats()` form one interval at a time
#' (the update factor `E`, the gain row `F`, the offset `b` and the
#' covariance offset `P` of the one-lag posterior) are vector expressions
#' over the record. Ported from aciR 0.2.3 (`.aci_online_aux()` in
#' `R/aci-online-smoother.R`, tag `parents-final`), written in the operations
#' of the matrix kernels in their order (the inverse of the filtered variance
#' as `chol_solve()` forms it, two divisions by its square root), so that for
#' `k = 1` every value agrees with the kernels' to the last bit up to the one
#' rounding a BLAS's triangular solve may take differently from a division
#' (one ulp on Accelerate, measured), and for `k > 1` it differs only by the
#' summation order inside the products the BLAS took.
#'
#' The covariance policy is not applied here. The kernels guard the filtered
#' variance at every interval; when every one of them is positive and finite
#' the guard never fires and the two routes agree, and when any is not this
#' returns `NULL` so that the caller takes the per-interval route, which
#' records the event and floors the value as the policy says.
#'
#' @param bundle A compiled CGNS bundle with `l = 1`.
#' @param filt A compatible explicit single-step filter path.
#' @param from First interval to form; intervals before it are not computed.
#' @returns A list with `idx` (the intervals formed), `E`, `b`, `P` (numeric,
#'   one per formed interval) and `F` (a matrix, one row per formed
#'   interval), or `NULL`.
#' @noRd
.online_aux_scalar <- function(bundle, filt, from = 1L) {
  N <- bundle$N
  k <- bundle$k
  dt <- bundle$dt
  idx <- seq.int(from, N)
  n <- length(idx)
  Rf <- as.numeric(filt$cov[1L, 1L, idx])
  if (!all(is.finite(Rf)) || any(Rf <= 0)) return(NULL)
  co <- bundle$coefficients
  Lx <- t(matrix(co$Lx[, 1L, idx], k, n))
  fx <- matrix(co$fx[idx, ], n, k)
  Ly <- as.numeric(co$Ly[1L, 1L, idx])
  fy <- as.numeric(co$fy[idx, 1L])
  gyy <- as.numeric(co$gyy[1L, 1L, idx])
  gyx <- t(matrix(co$gyx[1L, , idx], k, n))
  Gi <- array(co$gxx_weight[, , idx], c(k, k, n))
  muf <- as.numeric(filt$mean[idx, 1L])

  ## .thmD1_aux_compiled(), interval by interval, as vectors
  ch <- sqrt(Rf)
  Rfi <- (1 / ch) / ch
  Gx <- Lx + gyx * Rfi
  Gy <- Ly + gyy * Rfi
  K <- .scalar_aux_matvec(Gi, Gx)
  H <- Rfi * (Ly * Rf + Rf * Ly + gyy)
  E <- 1 + (rowSums(.scalar_aux_vecmat(gyx, Gi) * Gx) - Gy) * dt
  KR <- K * Rf
  GxKR <- rowSums(Gx * KR)
  KRK <- aperm(array(t(KR), c(k, n, k)), c(1L, 3L, 2L)) *
    array(rep(t(K), each = k), c(k, k, n))
  LxG <- .scalar_aux_vecmat(Lx, Gi + KRK * dt)
  F_ <- -Rf * (
    K + (GxKR * K - ((Rfi * H) * Rf) * K + Ly * K) * dt - LxG
  )

  ## .onelag_stats(), the offsets that do not depend on the later moments
  ILy <- 1 + Ly * dt
  dx <- bundle$x[idx + 1L, , drop = FALSE] - bundle$x[idx, , drop = FALSE]
  b <- muf - E * (ILy * muf + fy * dt) +
    rowSums(F_ * (dx - (Lx * muf + fx) * dt))
  P <- Rf - (E * ILy) * Rf - (rowSums(F_ * Lx) * Rf) * dt

  list(idx = idx, E = E, F = F_, b = b, P = P)
}


#' One-lag forward primitives for a scalar hidden state (internal)
#'
#' The per-interval fields that `.compiled_forward_primitives()` and the
#' primitives loop of `.lagtable_core_compiled()` collect, formed from
#' `.online_aux_scalar()` against the filter at the next time. The one-lag
#' variance goes through the same policy as `.onelag_stats()` applies: when
#' any is not positive and finite the result is `NULL` and the caller takes
#' the per-interval route.
#'
#' @param bundle A compiled CGNS bundle with `l = 1`.
#' @param filt A compatible explicit single-step filter path.
#' @param from First interval needed.
#' @returns A list with `E` (`1 x 1` matrices, `NULL` before `from`), `dmu`
#'   (`N x 1`), `dR` (`1 x 1` matrices), `one_mu` (`N x 1`), `one_R` (`1 x 1`
#'   matrices) and the norms `s_n`, `r_n`, `e_n` (length `N`, zero before
#'   `from`), or `NULL`.
#' @noRd
.forward_primitives_scalar <- function(bundle, filt, from = 1L) {
  aux <- .online_aux_scalar(bundle, filt, from = from)
  if (is.null(aux)) return(NULL)
  N <- bundle$N
  idx <- aux$idx
  muf <- as.numeric(filt$mean[idx, 1L])
  Rf <- as.numeric(filt$cov[1L, 1L, idx])
  mu1 <- aux$E * as.numeric(filt$mean[idx + 1L, 1L]) + aux$b
  R1 <- (aux$E * as.numeric(filt$cov[1L, 1L, idx + 1L])) * aux$E + aux$P
  if (!all(is.finite(R1)) || any(R1 <= 0)) return(NULL)
  dmu_v <- mu1 - muf
  dR_v <- R1 - Rf
  as_cells <- function(v) {
    out <- vector("list", N)
    out[idx] <- lapply(v, matrix, 1L, 1L)
    out
  }
  dmu <- matrix(NA_real_, N, 1L)
  dmu[idx, 1L] <- dmu_v
  one_mu <- matrix(NA_real_, N, 1L)
  one_mu[idx, 1L] <- mu1
  s_n <- r_n <- e_n <- numeric(N)
  s_n[idx] <- abs(dmu_v)
  r_n[idx] <- abs(dR_v)
  e_n[idx] <- abs(aux$E)
  list(
    E = as_cells(aux$E), dmu = dmu, dR = as_cells(dR_v),
    one_mu = one_mu, one_R = as_cells(R1), s_n = s_n, r_n = r_n, e_n = e_n,
    E_v = aux$E, dR_v = dR_v, one_R_v = R1
  )
}


#' Theorem-3 smoother moments for a scalar hidden state (internal)
#'
#' The backward recursion of `.smoother_thmD1_compiled()` on the auxiliaries
#' of `.online_aux_scalar()`: one scalar step per interval, in the operations
#' of `.onelag_stats()`. The smoothed variances go through the same policy
#' check as there; when any is not positive and finite the result is `NULL`
#' and the caller runs the per-interval route.
#'
#' @param bundle A compiled CGNS bundle with `l = 1`.
#' @param filt A compatible explicit single-step filter path.
#' @returns A list with `mu` and `cv` (length `N + 1`), or `NULL`.
#' @noRd
.smoother_thmD1_scalar_moments <- function(bundle, filt) {
  aux <- .online_aux_scalar(bundle, filt)
  if (is.null(aux)) return(NULL)
  N1 <- bundle$N1
  mu <- cv <- numeric(N1)
  mu[N1] <- filt$mean[N1, 1L]
  cv[N1] <- filt$cov[1L, 1L, N1]
  E <- aux$E
  b <- aux$b
  P <- aux$P
  for (j in (N1 - 1L):1L) {
    mu[j] <- E[j] * mu[j + 1L] + b[j]
    cv[j] <- (E[j] * cv[j + 1L]) * E[j] + P[j]
  }
  if (!all(is.finite(cv)) || any(cv <= 0)) return(NULL)
  list(mu = mu, cv = cv)
}
