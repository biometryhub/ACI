# Fixed-lag online smoother, vector case ---------------------------------------
#
# The matrix counterpart of aci-online-smoother.R. One thing does not carry
# over, and it is the thing that made the scalar version fast.
#
# In one dimension the ordered product of the per-step auxiliary matrices is a
# product of numbers, so it reduces to a difference of cumulative logarithms
# and any range is recoverable in constant time. Matrices do not commute, so
# there is no such reduction. The product over a range depends on the order of
# its factors and cannot be recovered from endpoint summaries.
#
# What does carry over is the reason the reduction was worth having. The source
# paper bounds the spectral radius of each factor below one, so the products
# decay geometrically whatever their order. The accumulation is therefore
# truncated once its norm falls below tolerance, which bounds the work by the
# effective lag rather than the record length. That was the substance of the
# scalar design; the logarithms were a constant factor on top of it.
#
# Symbols follow the source paper, with E_j, F_j, G, H and K as equations (3.5)
# to (3.7), b and P as (3.15) and (3.16), and the update matrix D accumulating
# on the RIGHT as (3.12).

#' Per-step auxiliary matrices of the vector online smoother
#'
#' Computes `E_j`, `F_j`, `b` and `P` at every step, and from them the
#' innovation each new observation contributes to earlier retained steps.
#'
#' The expressions are equations (3.5) to (3.7), (3.15) and (3.16) of the
#' source paper, in full generality. The reference implementation carries only
#' the specialisation for a zero noise cross-covariance, equation (3.8), which
#' is what its ENSO scripts use; the general form is implemented here because
#' the components schema exposes the cross-covariance and would otherwise offer
#' a path the recursion could not take.
#'
#' @param x Numeric matrix. The observed signal, one row per component.
#' @param comp A validated vector components list.
#' @param dt Numeric scalar. The integration time step.
#' @param filt A list with `mean` and `cov`, the filtered posterior.
#'
#' @returns A list with the arrays `E_j` and `F_j`, the matrix `innov_mean`
#'   and the array `innov_cov`.
#'
#' @noRd
#' @keywords internal
.aci_online_aux_mv <- function(x, comp, dt, filt) {
  n <- ncol(x)
  n_y <- nrow(filt$mean)
  identity <- diag(n_y)

  E_j <- array(0, c(n_y, n_y, n))
  F_j <- array(0, c(n_y, nrow(x), n))
  innov_mean <- matrix(0, n_y, n)
  innov_cov <- array(0, c(n_y, n_y, n))

  for (j in seq_len(n)) {
    L_x <- .aci_slice(comp$L_x, j)
    L_y <- .aci_slice(comp$L_y, j)
    S_yoS_y <- .aci_slice(comp$S_yoS_y, j)
    S_yoS_x <- .aci_slice(comp$S_yoS_x, j)
    S_xoS_y <- t(S_yoS_x)
    inv <- .aci_slice(comp$S_xoS_x_inv, j)
    R_f <- .aci_slice(filt$cov, j)
    R_f_inv <- chol2inv(.aci_chol(R_f))

    # (3.7)
    G_x <- L_x + S_xoS_y %*% R_f_inv
    G_y <- L_y + S_yoS_y %*% R_f_inv
    H_j <- R_f_inv %*% (L_y %*% R_f + R_f %*% t(L_y) + S_yoS_y)
    K_j <- inv %*% G_x

    # (3.5)
    E_j[, , j] <- identity + (S_yoS_x %*% inv %*% G_x - G_y) * dt

    # (3.6). Every term is l-by-k; the transposes are what make it so, and
    # they are invisible in one dimension, which is why the scalar form could
    # not simply be widened.
    F_j[, , j] <- -R_f %*% (
      t(K_j) +
        (t(G_x) %*% K_j %*% R_f %*% t(K_j) -
           R_f_inv %*% t(H_j) %*% R_f %*% t(K_j) +
           t(L_y) %*% t(K_j)) * dt -
        t(L_x) %*% (inv + K_j %*% R_f %*% t(K_j) * dt)
    )
  }

  # (3.13) to (3.16). The displacement each new observation imposes on the
  # step before it, which every earlier step then inherits, damped by D.
  for (k in seq_len(n - 1L)) {
    L_x <- .aci_slice(comp$L_x, k)
    L_y <- .aci_slice(comp$L_y, k)
    R_f <- .aci_slice(filt$cov, k)
    e_k <- E_j[, , k]
    f_k <- F_j[, , k]

    b_k <- filt$mean[, k] -
      e_k %*% ((diag(nrow(L_y)) + L_y * dt) %*% filt$mean[, k] +
                 .aci_column(comp$f_y, k) * dt) +
      f_k %*% (x[, k + 1L] - x[, k] -
                 (L_x %*% filt$mean[, k] + .aci_column(comp$f_x, k)) * dt)
    p_k <- R_f - e_k %*% (diag(nrow(L_y)) + L_y * dt) %*% R_f -
      f_k %*% L_x %*% R_f * dt

    innov_mean[, k] <- e_k %*% filt$mean[, k + 1L] + b_k - filt$mean[, k]
    innov_cov[, , k] <- e_k %*% .aci_slice(filt$cov, k + 1L) %*% t(e_k) +
      p_k - R_f
  }

  list(
    E_j = E_j, F_j = F_j,
    innov_mean = innov_mean, innov_cov = innov_cov
  )
}

#' Vector online smoother at a fixed lag
#'
#' @param x Numeric matrix. The observed signal.
#' @param comp A validated vector components list.
#' @param dt Numeric scalar. The integration time step.
#' @param filt A list with `mean` and `cov`.
#' @param lag Numeric scalar. Future steps admitted, or `Inf`.
#' @param tol Numeric scalar. Truncation tolerance on the update norm.
#'
#' @returns A list with `mean`, `cov` and `lag_effective`.
#'
#' @noRd
#' @keywords internal
.aci_online_smoother_mv <- function(x, comp, dt, filt, lag, tol) {
  n <- ncol(x)
  n_y <- nrow(filt$mean)
  aux <- .aci_online_aux_mv(x, comp, dt, filt)

  mean <- filt$mean
  cov <- filt$cov
  lag_effective <- 0L

  for (j in seq_len(n - 1L)) {
    # A lag of L admits the observations at k = j to j + L - 1, so it
    # contributes L terms and a lag of zero contributes none, leaving the
    # filter. Running to j + L instead would silently make every lag one step
    # longer than asked for, and would agree with the scalar path only at full
    # lag, where both saturate at the end of the record.
    last <- min(j + lag - 1L, n - 1L)
    if (last < j) {
      next
    }
    d <- diag(n_y)
    for (k in seq.int(j, last)) {
      mean[, j] <- mean[, j] + d %*% aux$innov_mean[, k]
      cov[, , j] <- cov[, , j] + d %*% aux$innov_cov[, , k] %*% t(d)
      lag_effective <- max(lag_effective, k - j + 1L)
      # D accumulates on the right, per equation (3.12).
      d <- d %*% aux$E_j[, , k]
      if (max(abs(d)) < tol) {
        break
      }
    }
    cov[, , j] <- (cov[, , j] + t(cov[, , j])) / 2
  }

  list(mean = mean, cov = cov, lag_effective = as.integer(lag_effective))
}

#' One row of the vector causal-influence-range divergence sequence
#'
#' @param aux Auxiliary matrices from [.aci_online_aux_mv()].
#' @param filt A list with `mean` and `cov`.
#' @param j Integer scalar. The reported step.
#' @param n Integer scalar. Length of the observed signal.
#' @param tol Numeric scalar. Truncation tolerance.
#' @param horizon Integer scalar. How far forward the comparison may look. The
#'   fully informed posterior it is compared against is unaffected and is
#'   always taken over the whole record.
#'
#' @returns A numeric vector of divergences, or `NULL` at the final step.
#'
#' @noRd
#' @keywords internal
.aci_cir_row_mv <- function(aux, filt, j, n, tol, horizon = n) {
  if (j >= n) {
    return(NULL)
  }
  n_y <- nrow(filt$mean)
  count <- n - j + 1L
  mu <- matrix(0, n_y, count)
  rr <- array(0, c(n_y, n_y, count))
  mu[, 1L] <- filt$mean[, j]
  rr[, , 1L] <- .aci_slice(filt$cov, j)

  d <- diag(n_y)
  for (k in seq.int(j, n - 1L)) {
    i <- k - j + 2L
    mu[, i] <- mu[, i - 1L] + d %*% aux$innov_mean[, k]
    rr[, , i] <- rr[, , i - 1L] + d %*% aux$innov_cov[, , k] %*% t(d)
    rr[, , i] <- (rr[, , i] + t(rr[, , i])) / 2
    d <- d %*% aux$E_j[, , k]
    if (max(abs(d)) < tol && i < count) {
      # The remaining updates cannot move the estimate, so the rest of the
      # sequence is the fully informed value it has already reached.
      for (m in seq.int(i + 1L, count)) {
        mu[, m] <- mu[, i]
        rr[, , m] <- rr[, , i]
      }
      break
    }
  }

  mu_end <- mu[, count]
  r_end <- rr[, , count]
  value <- .aci_metric_mv(
    list(mean = mu, cov = rr),
    list(
      mean = matrix(rep(mu_end, count), n_y, count),
      cov = array(rep(as.numeric(r_end), count), c(n_y, n_y, count))
    )
  )$value
  # Truncated only after mu_end and r_end have been taken over the whole
  # record; see the note in the scalar .aci_cir_row().
  value[seq_len(min(count, horizon - j + 1L))]
}
