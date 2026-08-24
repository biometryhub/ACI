################################################################################
## utils.R - foundations: conditions, linear algebra, observation contract
################################################################################


#' Signal a classed aci error (internal)
#'
#' Raises an error whose condition object carries the supplied subclasses ahead
#' of `aci_error`, so callers can handle failures programmatically.
#'
#' @param class Character vector of condition subclasses to prepend.
#' @param msg 1-length character message.
#' @param ... Further named fields stored on the condition object.
#' @returns Never returns; called for the condition it raises.
#' @noRd
aci_abort <- function(class, msg, ...) {
  stop(structure(class = c(class, "aci_error", "error", "condition"),
                 list(message = msg, call = sys.call(-1), ...)))
}


#' Signal a classed aci warning (internal)
#'
#' Raises a warning whose condition object carries the supplied subclasses ahead
#' of `aci_warning`.
#'
#' @param class Character vector of condition subclasses to prepend.
#' @param msg 1-length character message.
#' @param ... Further named fields stored on the condition object.
#' @returns Invisibly `NULL`; called for the condition it raises.
#' @noRd
aci_warn <- function(class, msg, ...) {
  warning(structure(class = c(class, "aci_warning", "warning", "condition"),
                    list(message = msg, call = sys.call(-1), ...)))
}


#' Default value for `NULL` (internal)
#'
#' @param a Value to use when it is not `NULL`.
#' @param b Fallback value used when `a` is `NULL`.
#' @returns `a` unless it is `NULL`, otherwise `b`.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Symmetric part of a square matrix (internal)
#'
#' @param R Square numeric matrix.
#' @returns The symmetric part `(R + t(R)) / 2`.
#' @noRd
sym <- function(R) (R + t(R)) / 2


#' Cholesky factor without regularization (internal)
#'
#' Validates that a matrix is finite, square, symmetric and positive definite,
#' and returns its Cholesky factor.
#'
#' @param R Square numeric matrix.
#' @param where 1-length character naming the matrix in error messages.
#' @returns Upper-triangular Cholesky factor of `R`.
#' @noRd
.strict_chol <- function(R, where = "covariance") {
  R <- as.matrix(R)
  if (!is.numeric(R) || nrow(R) < 1L || nrow(R) != ncol(R) ||
      any(!is.finite(R)))
    aci_abort("aci_error_spd",
              sprintf("Matrix (%s) must be finite, numeric, square, and non-empty.", where))
  if (max(abs(R - t(R))) > 1e-12 * max(1, max(abs(R))))
    aci_abort("aci_error_spd", sprintf("Matrix (%s) must be symmetric.", where))
  R <- sym(R)
  ch <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(ch))
    aci_abort("aci_error_spd",
              sprintf("Matrix (%s) must be positive definite.", where))
  ch
}


#' Cholesky with escalating jitter ladder
#'
#' Returns the Cholesky factor of a square numeric matrix, adding an escalating
#' ridge to the diagonal when the unmodified matrix is not numerically positive
#' definite. The ladder runs from `1e-12` to `1e-6` relative to the mean
#' diagonal entry.
#'
#' @param R Symmetric covariance matrix.
#' @param where Context label used in numerical error messages.
#' @returns Upper-triangular Cholesky factor of `R`, possibly of `R` plus a
#'   diagonal ridge.
#'
#' @seealso [spd_floor()]
#'
#' @examples
#' safe_chol(matrix(c(2, 0.5, 0.5, 1), 2, 2))
#'
#' @export
safe_chol <- function(R, where = "covariance") {
  R <- as.matrix(R)
  if (nrow(R) < 1L || nrow(R) != ncol(R) || any(!is.finite(R)))
    aci_abort("aci_error_spd", sprintf("Matrix (%s) must be finite, square, and non-empty.", where))
  R <- sym(R)
  ch <- tryCatch(chol(R), error = function(e) NULL)
  if (!is.null(ch)) return(ch)
  l <- nrow(R); base <- sum(diag(R)) / max(l, 1)
  if (!is.finite(base) || base <= 0) base <- 1
  for (eta in 10^seq(-12, -6)) {
    ch <- tryCatch(chol(R + eta * base * diag(l)), error = function(e) NULL)
    if (!is.null(ch)) return(ch)
  }
  aci_abort("aci_error_spd", sprintf(
    "Matrix (%s) is not positive definite even after jitter ladder; min eig = %.3e.",
    where, tryCatch(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values),
                    error = function(e) NA_real_)))
}


#' Project to SPD by eigenvalue flooring (used after Euler covariance updates)
#'
#' Returns the matrix unchanged when it already admits a Cholesky factor, and
#' otherwise floors its eigenvalues at a small positive multiple of the largest
#' absolute eigenvalue before reassembling it.
#'
#' @param R Symmetric matrix.
#' @param eps Minimum eigenvalue.
#' @returns A symmetric positive-definite matrix of the same dimension as `R`.
#'
#' @seealso [safe_chol()]
#'
#' @examples
#' spd_floor(matrix(c(1, 2, 2, 1), 2, 2))
#'
#' @export
spd_floor <- function(R, eps = 1e-12) {
  R <- as.matrix(R)
  if (nrow(R) < 1L || nrow(R) != ncol(R) || any(!is.finite(R)) ||
      length(eps) != 1L || !is.finite(eps) || eps <= 0)
    aci_abort("aci_error_spd", "spd_floor() needs a finite non-empty square matrix and eps > 0.")
  R <- sym(R)
  ok <- tryCatch({ chol(R); TRUE }, error = function(e) FALSE)
  if (ok) return(R)
  e <- eigen(R, symmetric = TRUE)
  fl <- max(eps * max(abs(e$values), 1e-300), 1e-300)
  sym(e$vectors %*% (pmax(e$values, fl) * t(e$vectors)))
}


#' Solve a symmetric positive-definite system by Cholesky (internal)
#'
#' @param R Symmetric positive-definite matrix, the left-hand side.
#' @param B Numeric matrix or vector, the right-hand side.
#' @param where 1-length character naming the system in error messages.
#' @returns The solution of `R x = B`.
#' @noRd
chol_solve <- function(R, B, where = "solve") {
  ch <- safe_chol(R, where)
  backsolve(ch, forwardsolve(t(ch), B))
}


#' Log-determinant via the Cholesky factor (internal)
#'
#' @param R Symmetric positive-definite matrix.
#' @param where 1-length character naming the matrix in error messages.
#' @returns 1-length numeric, the log-determinant of `R`.
#' @noRd
logdet_chol <- function(R, where = "logdet") 2 * sum(log(diag(safe_chol(R, where))))


#' Masked Gram inverse for non-target inflation (internal)
#'
#' Builds an inverse observation Gram matrix supported only on the target block,
#' giving zero filter weight to the innovations of the remaining channels.
#'
#' @param gxx Square numeric observation Gram matrix.
#' @param idxA Integer vector of target row/column indices within `gxx`.
#' @returns A matrix of the same dimension as `gxx`, zero outside the target
#'   block and equal to the inverse of the target block inside it.
#' @noRd
masked_ginv <- function(gxx, idxA) {
  gxx <- as.matrix(gxx)
  if (!is.numeric(gxx) || nrow(gxx) < 1L || nrow(gxx) != ncol(gxx) ||
      any(!is.finite(gxx)))
    aci_abort("aci_error_dims", "gxx must be a finite non-empty square matrix.")
  k <- nrow(gxx); M <- matrix(0, k, k)
  if (!is.numeric(idxA) || any(!is.finite(idxA)) ||
      any(idxA != floor(idxA)) || any(idxA < 1L) || any(idxA > k) ||
      anyDuplicated(idxA))
    aci_abort("aci_error_dims", "idxA contains invalid Gram-matrix indices.")
  if (length(idxA) == 0) return(M)
  idxA <- as.integer(idxA)
  M[idxA, idxA] <- chol_solve(gxx[idxA, idxA, drop = FALSE],
                              diag(length(idxA)), "masked Gram")
  M
}


#' Observed trajectory on a uniform time grid
#'
#' Constructs the observation object consumed throughout the package. The grid
#' must be strictly increasing and uniformly spaced, and the observations must
#' be finite; version 0 additionally assumes the observations are effectively
#' noise-free. The noise-free restriction follows the method's current
#' published scope: andreou2026cir Section 2.1 leaves noise-contaminated
#' observations to future work, and its Discussion lists them as an open
#' direction, so this is a limitation of the framework as published, not a
#' modelling assumption added by the package.
#'
#' @param t Numeric vector of observation times, strictly increasing and
#'   uniformly spaced.
#' @param x Numeric matrix of observations with one row per time and at least
#'   one column, or a vector coerced to a single column.
#' @param noise_free Logical; must be `TRUE`, since noisy observations are not
#'   supported in this version.
#' @param names Optional character vector of unique, non-empty column names,
#'   one per observed channel.
#' @returns An object of class `obs_traj`: a list with the time vector `t`, the
#'   observation matrix `x`, the step `dt`, the observed dimension `k` and the
#'   flag `noise_free`.
#'
#' @seealso [as_obs()], [da_filter()]
#'
#' @examples
#' observed_trajectory(t = seq(0, 1, by = 0.1),
#'                     x = matrix(rnorm(11), ncol = 1))
#'
#' @export
observed_trajectory <- function(t, x, noise_free = TRUE, names = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(t) || !is.numeric(x))
    aci_abort("aci_error_obs_contract", "t and x must be numeric.")
  if (length(t) != nrow(x))
    aci_abort("aci_error_dims", "length(t) must equal nrow(x).")
  if (ncol(x) < 1L)
    aci_abort("aci_error_dims", "x must contain at least one observed column.")
  if (any(!is.finite(t)) || any(!is.finite(x)))
    aci_abort("aci_error_obs_contract", "Non-finite values in observations.")
  dtv <- diff(t)
  if (length(dtv) < 1 || any(dtv <= 0))
    aci_abort("aci_error_obs_contract", "t must be strictly increasing.")
  if (max(abs(dtv - dtv[1])) > 1e-8 * max(abs(dtv[1]), 1e-12))
    aci_abort("aci_error_obs_contract",
      "v0 requires a uniform time grid (invariant #11); resample first.")
  if (!isTRUE(noise_free))
    aci_abort("aci_error_obs_contract",
      "v0 assumes effectively noise-free observations (invariant #11); noisy obs are a roadmap item.")
  if (!is.null(names)) {
    if (!is.character(names) || length(names) != ncol(x) || anyNA(names) ||
        any(!nzchar(names)) || anyDuplicated(names))
      aci_abort("aci_error_dims", "names must contain one unique non-empty name per observed column.")
    colnames(x) <- names
  }
  structure(list(t = as.numeric(t), x = x, dt = dtv[1], k = ncol(x),
                 noise_free = TRUE),
            class = "obs_traj")
}


#' Coerce to an observed trajectory
#'
#' Generic converting the supported observation representations to an
#' `obs_traj` object.
#'
#' @param x Object to coerce: an `obs_traj`, a numeric matrix, or a simulation
#'   of class `aci_sim`.
#' @param ... Arguments passed to methods.
#' @returns An object of class `obs_traj`; see [observed_trajectory()].
#'
#' @seealso [observed_trajectory()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' as_obs(sim)
#'
#' @export
as_obs <- function(x, ...) UseMethod("as_obs")


#' @describeIn as_obs Returns the trajectory unchanged.
#' @export
as_obs.obs_traj <- function(x, ...) x


#' @describeIn as_obs Builds a uniform grid from `dt` and `t0` for a matrix of
#'   observations.
#' @param dt Positive 1-length numeric step used to build the time grid.
#' @param t0 1-length numeric time of the first observation.
#' @export
as_obs.matrix <- function(x, dt = 1, t0 = 0, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied while coercing a matrix observation.")
  if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0 ||
      !is.numeric(t0) || length(t0) != 1L || !is.finite(t0))
    aci_abort("aci_error_obs_contract",
              "dt must be one finite positive value and t0 one finite value.")
  observed_trajectory(t0 + dt * (seq_len(nrow(x)) - 1), x)
}


#' @describeIn as_obs Extracts the observation component of a simulation.
#' @export
as_obs.aci_sim <- function(x, ...) x$obs


#' Print an observed trajectory
#'
#' @param x An `obs_traj` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.obs_traj <- function(x, ...) {
  cat(sprintf("<obs_traj> k = %d, N+1 = %d, dt = %g, span [%g, %g]\n",
              x$k, length(x$t), x$dt, x$t[1], x$t[length(x$t)]))
  invisible(x)
}


#' Coerce an observed trajectory to a data frame
#'
#' @param x An `obs_traj` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns A data frame in long form with columns `t`, `var` and `value`.
#' @export
as.data.frame.obs_traj <- function(x, ...) {
  nm <- colnames(x$x) %||% paste0("x", seq_len(x$k))
  data.frame(t = rep(x$t, x$k), var = rep(nm, each = length(x$t)),
             value = as.vector(x$x))
}
