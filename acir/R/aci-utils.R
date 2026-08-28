################################################################################
## aci-utils.R - foundations: conditions, linear algebra, observation contract
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
#' Called once or twice per step inside the recursions, so the transpose goes
#' straight to `t.default()`: every caller passes a plain numeric matrix, which
#' is the method the generic dispatches to anyway.
#'
#' @param R Square numeric matrix.
#' @returns The symmetric part `(R + t(R)) / 2`.
#' @noRd
sym <- function(R) (R + t.default(R)) / 2


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


#' `spd_floor()` of the symmetric part, without the guarded entry (internal)
#'
#' `spd_floor(sym(R))` is the covariance policy of the Euler recursions, and
#' those recursions apply it at every substep. This is the same policy written
#' for that position: `spd_floor()` re-symmetrises its argument, and `sym()` is
#' idempotent to the bit, so the entry checks and the second symmetrisation are
#' redundant when the caller has just formed `R` arithmetically. A matrix that
#' still fails is handed to [spd_floor()] itself, which reproduces the
#' eigenvalue floor and the non-finite error unchanged.
#'
#' @param R Square numeric matrix.
#' @returns `spd_floor(sym(R))`.
#' @noRd
.sym_floor <- function(R) {
  S <- (R + t.default(R)) / 2
  if (is.null(tryCatch(chol.default(S), error = function(e) NULL)))
    spd_floor(S) else S
}


# =========================================================================
# covariance policy: resolution, recording, guards
# =========================================================================

#' Resolve the covariance policy for one call (internal)
#'
#' `NULL` reads `getOption("aci.regularize")`, the pattern [lag_table()]
#' already uses for `aci.default_tol`. The resolved value, never the option,
#' is what a result records, so a saved result stays self-describing when the
#' option later changes.
#'
#' @param regularize `NULL`, `"none"` or `"floor"`.
#' @returns `"none"` or `"floor"`.
#' @noRd
.aci_regularize <- function(regularize = NULL) {
  x <- regularize %||% getOption("aci.regularize", "none")
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !x %in% c("none", "floor"))
    aci_abort("aci_error_dims", "regularize must be \"none\" or \"floor\".")
  x
}


#' Per-run regularisation recorder (internal)
#'
#' One environment per public call. `j` is set once per outer grid step by the
#' kernel, so the guards need no index argument in the hot loop; everything
#' else is touched only when a floor actually fires.
#'
#' @param policy `"none"` or `"floor"`.
#' @param tgrid Numeric time grid of the run.
#' @returns An environment.
#' @noRd
.aci_reg_new <- function(policy, tgrid) {
  e <- new.env(parent = emptyenv())
  e$policy <- policy
  e$strict <- identical(policy, "none")
  e$tgrid <- tgrid
  e$j <- NA_integer_
  e$n <- 0L
  e$sites <- list()
  e
}


#' Recorder for one kernel, created or inherited (internal)
#'
#' Internal kernels take their policy in one argument that is either the
#' user-facing `"none"`/`"floor"`/`NULL`, or a recorder a caller has already
#' created. Passing the recorder is how one public call shares a single record
#' across the filter, the smoother and the metric it feeds.
#'
#' @param regularize `NULL`, a policy string, or a recorder.
#' @param tgrid Numeric time grid used when a recorder must be created.
#' @returns An environment from [.aci_reg_new()].
#' @noRd
.aci_reg_for <- function(regularize, tgrid) {
  if (is.environment(regularize)) return(regularize)
  .aci_reg_new(.aci_regularize(regularize), tgrid)
}


#' The zero-row site frame every clean run reports (internal)
#'
#' Built once at load. A kernel that took no floor is the overwhelmingly
#' common case and is on the timed path, and constructing this frame per call
#' is about 30 microseconds - a measurable share of a millisecond-scale dyad
#' filter. One shared immutable value costs nothing and is `identical()` to
#' the constructed frame.
#'
#' @noRd
.ACI_REG_NO_SITES <- data.frame(
  site = character(0), role = character(0), n = integer(0),
  first_index = integer(0), first_time = numeric(0),
  worst_value = numeric(0), stringsAsFactors = FALSE)


#' Freeze a recorder into the `meta$regularization` contract (internal)
#'
#' Always returns a list, never `NULL`, so a downstream consumer never has to
#' test for absence.
#'
#' @param rec A recorder from [.aci_reg_new()].
#' @returns A list with `policy`, `fired`, `n_events`, `eps` and `sites`.
#' @noRd
.aci_reg_freeze <- function(rec) {
  if (!length(rec$sites))
    return(list(policy = rec$policy, fired = FALSE, n_events = 0L,
                eps = 1e-12, sites = .ACI_REG_NO_SITES))
  sites <- do.call(rbind, lapply(rec$sites, function(s)
    data.frame(site = s$site, role = s$role, n = s$n,
               first_index = s$first_index, first_time = s$first_time,
               worst_value = s$worst_value, stringsAsFactors = FALSE)))
  rownames(sites) <- NULL
  list(policy = rec$policy, fired = rec$n > 0L, n_events = rec$n,
       eps = 1e-12, sites = sites)
}


#' Record one floor event (internal)
#'
#' @param rec Recorder.
#' @param site Short site id.
#' @param index Grid index.
#' @param time Grid time.
#' @param value Offending value.
#' @returns Invisibly `NULL`.
#' @noRd
.aci_reg_note <- function(rec, site, index, time, value) {
  rec$n <- rec$n + 1L
  s <- rec$sites[[site]]
  if (is.null(s))
    s <- list(site = site, role = .ACI_COV_SITES[[site]]$role, n = 0L,
              first_index = index, first_time = time, worst_value = value)
  s$n <- s$n + 1L
  if (is.finite(value) &&
      (!is.finite(s$worst_value) || value < s$worst_value))
    s$worst_value <- value
  rec$sites[[site]] <- s
  invisible(NULL)
}


#' Site register: role and prose used in errors and records (internal)
#'
#' Keyed by the short site id passed at each guarded call site, so the wording
#' of a covariance failure lives in one place and a later release can refine
#' the policy per role without touching the kernels.
#'
#' @noRd
.ACI_COV_SITES <- list(
  filter_explicit      = list(role = "filter covariance",
                              where = "explicit Riccati step"),
  filter_implicit_q    = list(role = "filter covariance",
                              where = "implicit predicted noise term"),
  filter_implicit_p    = list(role = "filter covariance",
                              where = "implicit predicted covariance"),
  filter_implicit_info = list(role = "implicit information matrix",
                              where = "implicit information-matrix inverse"),
  filter_implicit_r    = list(role = "filter covariance",
                              where = "implicit posterior covariance"),
  smoother_backward    = list(role = "smoother covariance",
                              where = "backward Riccati step"),
  smoother_filter_cov  = list(role = "filter covariance",
                              where = "backward smoother's filtered-covariance inverse"),
  smoother_onelag      = list(role = "one-lag smoothed covariance",
                              where = "Theorem 3 one-lag statistics"),
  onelag_filter_cov    = list(role = "filter covariance",
                              where = "Theorem 3 auxiliary matrices"),
  likelihood_innov     = list(role = "innovation covariance",
                              where = "predictive log-likelihood"),
  metric_reference     = list(role = "reference covariance",
                              where = "relative-entropy denominator")
)


#' Abort on a covariance that left the positive-definite cone (internal)
#'
#' Name the quantity, say what it became, say where, then give the remedies -
#' the shape `aciR`'s `.aci_stop_covariance()` uses - plus the sentence this
#' package needs and `aciR` does not: how to get the previous behaviour back.
#'
#' @param site Short site id, a name of `.ACI_COV_SITES`.
#' @param index Grid index at which the covariance failed.
#' @param time Grid time at `index`.
#' @param value Offending value: the smallest eigenvalue on the matrix path,
#'   the variance itself on the scalar path.
#' @returns Never returns; called for the condition it raises.
#' @noRd
.aci_stop_cov <- function(site, index, time, value) {
  info <- .ACI_COV_SITES[[site]]
  aci_abort(
    c("aci_error_covariance_not_spd", "aci_error_spd"),
    sprintf(paste0(
      "The %s must stay finite and positive definite; it reached %s at ",
      "index %d (time %g), in the %s. Reduce dt, raise nsub, or use ",
      "stepper = \"implicit\", which preserves positivity. To keep the ",
      "previous behaviour, call with regularize = \"floor\"; every floored ",
      "step is then recorded in the result's meta$regularization."),
      info$role, format(value), index, time, info$where),
    site = site, role = info$role, index = index, time = time, value = value)
}


#' Cold half of the matrix guards: stop, or record and report (internal)
#'
#' Reached only after a factorisation has already failed, which the fire-count
#' evidence puts at 29 events in 4.37e6 invocations across every standard case
#' and every stress probe in the workspace.
#'
#' @param S Symmetric matrix that failed its Cholesky factorisation.
#' @param rec Recorder.
#' @param site Short site id.
#' @returns Invisibly `NULL` under `"floor"`; never returns under `"none"`.
#' @noRd
.cov_guard_event <- function(S, rec, site) {
  lam <- tryCatch(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values),
                  error = function(e) NA_real_)
  tt <- rec$tgrid[rec$j]
  if (isTRUE(rec$strict)) .aci_stop_cov(site, rec$j, tt, lam)
  .aci_reg_note(rec, site, rec$j, tt, lam)
  invisible(NULL)
}


#' Guarded symmetric covariance update, matrix path (internal)
#'
#' Replaces `.sym_floor()` at the recursion call sites. The non-firing path is
#' `.sym_floor()`'s own arithmetic plus two argument promises; everything
#' policy-related happens only after `chol.default()` has already failed.
#'
#' @param R Square numeric matrix, freshly formed by the caller.
#' @param rec Recorder from [.aci_reg_new()].
#' @param site Short site id.
#' @returns The symmetric part of `R`, or its eigenvalue-floored projection
#'   under `regularize = "floor"`.
#' @noRd
.cov_guard <- function(R, rec, site) {
  S <- (R + t.default(R)) / 2
  if (!is.null(tryCatch(chol.default(S), error = function(e) NULL))) return(S)
  .cov_guard_event(S, rec, site)
  spd_floor(S)
}


#' Guarded variance update, scalar path (internal)
#'
#' Replaces the four inline `if (R <= 0) R <- max(...)` floors in
#' `aci-kernels-scalar.R` and the likelihood's `.scalar_spd_floor()` call. The
#' callers keep today's `!is.finite(v) || v <= 0` entry test, so the
#' non-firing path is the two comparisons it always was and this body is cold.
#'
#' @param value One variance.
#' @param rec Recorder.
#' @param site Short site id.
#' @returns `value` when positive; the scalar floor under
#'   `regularize = "floor"`.
#' @noRd
.cov_guard_scalar <- function(value, rec, site) {
  tt <- rec$tgrid[rec$j]
  ## A non-finite variance is not recoverable by regularisation, so it aborts
  ## under both policies.
  if (!is.finite(value))
    aci_abort(c("aci_error_covariance_not_spd", "aci_error_spd"),
              sprintf(paste0("The %s became %s at index %d (time %g). ",
                             "This is not recoverable by regularisation; ",
                             "reduce dt or raise nsub."),
                      .ACI_COV_SITES[[site]]$role, format(value), rec$j, tt),
              site = site, role = .ACI_COV_SITES[[site]]$role,
              index = rec$j, time = tt, value = value)
  if (value > 0) return(value)
  if (isTRUE(rec$strict)) .aci_stop_cov(site, rec$j, tt, value)
  .aci_reg_note(rec, site, rec$j, tt, value)
  max(1e-12 * max(abs(value), 1e-300), 1e-300)
}


#' Guarded Cholesky of a reference covariance (internal)
#'
#' Replaces `.chol_or_floor()` at the metric and likelihood call sites, and
#' returns the same `list(R = , ch = )` shape.
#'
#' @param R Square numeric matrix.
#' @param rec Recorder.
#' @param site Short site id.
#' @returns A list with the possibly floored `R` and its Cholesky factor `ch`.
#' @noRd
.cov_guard_chol <- function(R, rec, site) {
  ch <- tryCatch(chol(R), error = function(e) NULL)
  if (!is.null(ch)) return(list(R = R, ch = ch))
  .cov_guard_event(R, rec, site)
  Rf <- spd_floor(R)
  list(R = Rf, ch = chol(Rf))
}


#' Guarded Cholesky factor for a governed linear solve (internal)
#'
#' The jitter ladder is a second silent regularisation inside the recursions,
#' so the four recursion call sites route their factorisation through the
#' policy. Under `"floor"` the failing branch hands straight back to
#' [safe_chol()], so the ladder and its numbers are exactly what they were.
#'
#' @param R Square numeric matrix.
#' @param rec Recorder.
#' @param site Short site id.
#' @param where 1-length character naming the system in error messages.
#' @returns The upper-triangular Cholesky factor used by [chol_solve()].
#' @noRd
.cov_guard_factor <- function(R, rec, site, where) {
  S <- sym(as.matrix(R))
  ch <- tryCatch(chol.default(S), error = function(e) NULL)
  if (!is.null(ch)) return(ch)
  .cov_guard_event(S, rec, site)
  safe_chol(R, where)
}


#' Solve a symmetric positive-definite system by Cholesky (internal)
#'
#' @param R Symmetric positive-definite matrix, the left-hand side.
#' @param B Numeric matrix or vector, the right-hand side.
#' @param where 1-length character naming the system in error messages.
#' @param rec Optional recorder. Supplied at the four recursion call sites, so
#'   that the jitter ladder is reached only under `regularize = "floor"`;
#'   `NULL` everywhere else keeps [safe_chol()]'s behaviour unchanged.
#' @param site Short site id, required when `rec` is supplied.
#' @returns The solution of `R x = B`.
#' @noRd
chol_solve <- function(R, B, where = "solve", rec = NULL, site = NULL) {
  ch <- if (is.null(rec)) safe_chol(R, where) else
    .cov_guard_factor(R, rec, site, where)
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
#' @seealso [as_obs()], [aci_filter()]
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
#' m <- aci_dyad_model()
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
