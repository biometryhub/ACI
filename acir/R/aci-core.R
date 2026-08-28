################################################################################
## aci-core.R - causal layer: ACI, CIR, KL metric
## ########################################################################## ##
##
## Contents:
##   * Gaussian relative entropy, signal/dispersion decomposition, paths:
##       - aci_metric_pair, aci_metric
##
##   * the ACI metric front-end aci():
##       - aci, print.aci_result, as.data.frame.aci_result
##
##   * forward causal influence ranges (exact + l1_linf):
##       - .fwd_lengths, new_cir_result, print.cir_result, aci_range,
##         aci_range.lag_table, aci_range.aci_result
##
################################################################################


################################################################################
# Gaussian relative entropy, signal/dispersion decomposition, paths
################################################################################

#' Gaussian relative entropy
#'
#' Relative entropy of one multivariate Gaussian from another, optionally split
#' into its signal and dispersion parts.
#'
#' Public KL values never apply a covariance ridge; callers wanting
#' regularisation opt in explicitly with [spd_floor()]. andreou2026cir
#' (Section 2.2, closing paragraph) states that regularization is expected
#' only in the degenerate limit, which is the published basis for this
#' function's strictness. The state recursions are strict by default too, and
#' regularise only when a call asks for it with `regularize = "floor"`; see
#' [aci_filter()].
#'
#' @param mu_p Numeric vector, mean of the first distribution.
#' @param R_p Covariance matrix of the first distribution.
#' @param mu_q Numeric vector, mean of the second distribution.
#' @param R_q Covariance matrix of the second distribution.
#' @param decompose `TRUE` to return the signal and dispersion parts alongside
#'   the total.
#' @returns A named numeric vector with the `total` and, when `decompose` is
#'   `TRUE`, the `signal` and `dispersion` parts.
#'
#' @seealso [aci_metric()], [aci()]
#'
#' @examples
#' aci_metric_pair(mu_p = 0, R_p = matrix(1), mu_q = 1, R_q = matrix(2))
#'
#' @export
aci_metric_pair <- function(mu_p, R_p, mu_q, R_q, decompose = TRUE) {
  mu_p <- as.numeric(mu_p); mu_q <- as.numeric(mu_q)
  R_p <- as.matrix(R_p); R_q <- as.matrix(R_q)
  l <- length(mu_p)
  if (!is.numeric(R_p) || !is.numeric(R_q) || l < 1L || length(mu_q) != l ||
      !identical(dim(R_p), c(l, l)) || !identical(dim(R_q), c(l, l)) ||
      any(!is.finite(c(mu_p, mu_q, R_p, R_q))))
    aci_abort("aci_error_dims", "aci_metric_pair: dimension mismatch.")
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  # Public KL values must not depend on an undocumented covariance ridge.
  # Callers wanting regularisation can opt in explicitly with spd_floor().
  Lq <- .strict_chol(R_q, "R_q"); Lp <- .strict_chol(R_p, "R_p")
  d <- mu_q - mu_p
  w <- forwardsolve(t(Lq), d)                     # Lq' \ d
  signal <- 0.5 * sum(w * w)
  A <- forwardsolve(t(Lq), t(Lp))                 # tr(Rq^{-1}Rp) = ||Lq'^{-1} Lp'||_F^2
  tr <- sum(A * A)
  dispersion <- 0.5 * (tr - l + 2 * sum(log(diag(Lq))) - 2 * sum(log(diag(Lp))))
  signal <- max(signal, 0); dispersion <- max(dispersion, 0)
  tot <- signal + dispersion
  if (isTRUE(getOption("aci.debug_assert", FALSE)) && tot < -1e-10)
    aci_abort("aci_error_internal", "Negative KL encountered.")
  if (decompose) c(total = tot, signal = signal, dispersion = dispersion) else tot
}


#' Cholesky of one path covariance slice on the KL-path contract (internal)
#'
#' [aci_metric()] establishes the shape, dimension and finiteness of the
#' whole covariance array once, before the loop, so per index only symmetry and
#' positive definiteness remain to check. This is `.strict_chol()` with the
#' checks the caller has already discharged removed and the transpose shared
#' between the symmetry test and the symmetrisation. The arithmetic, the
#' tolerance and the error classes and messages are unchanged, so an index that
#' `.strict_chol()` rejects is rejected here in the same way.
#'
#' @param R Square finite numeric matrix.
#' @param where 1-length character naming the matrix in error messages.
#' @returns Upper-triangular Cholesky factor of the symmetric part of `R`.
#' @noRd
.kl_path_chol <- function(R, where) {
  tR <- t.default(R)
  if (max(abs(R - tR)) > 1e-12 * max(1, max(abs(R))))
    aci_abort("aci_error_spd", sprintf("Matrix (%s) must be symmetric.", where))
  ch <- tryCatch(chol.default((R + tR) / 2), error = function(e) NULL)
  if (is.null(ch))
    aci_abort("aci_error_spd",
              sprintf("Matrix (%s) must be positive definite.", where))
  ch
}


#' Gaussian relative entropy along a pair of paths
#'
#' Evaluates the relative entropy at every time of two assimilation paths on a
#' common grid. Gaussian relative entropy is oriented as smoother relative to
#' filter.
#'
#' @param p A `da_path_gaussian` object, the integrating distribution.
#' @param q A `da_path_gaussian` object on the same time grid.
#' @param decompose `TRUE` to return the signal and dispersion parts alongside
#'   the total.
#' @returns A data frame with the time column `t` and either `total` alone or
#'   `total`, `signal` and `dispersion`.
#'
#' @seealso [aci_metric_pair()], [aci()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- aci_filter(m, ob)
#' s <- aci_smoother(m, ob, filter = f)
#' head(aci_metric(s, f))
#'
#' @export
aci_metric <- function(p, q, decompose = TRUE) {
  if (!inherits(p, "da_path_gaussian") || !inherits(q, "da_path_gaussian"))
    aci_abort("aci_error_dims", "p and q must be Gaussian assimilation paths.")
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  if (length(p$t) != length(q$t) || max(abs(p$t - q$t)) > 1e-10)
    aci_abort("aci_error_dims", "aci_metric: mismatched time grids.")
  n <- length(p$t)
  if (!is.matrix(p$mean) || !is.matrix(q$mean))
    aci_abort("aci_error_dims", "aci_metric: path means must be matrices.")
  lp <- ncol(p$mean); lq <- ncol(q$mean)
  if (n < 1L || lp < 1L || lq != lp ||
      !identical(dim(p$mean), c(n, lp)) ||
      !identical(dim(q$mean), c(n, lp)) ||
      !identical(dim(p$cov), c(lp, lp, n)) ||
      !identical(dim(q$cov), c(lp, lp, n)) ||
      any(!is.finite(c(p$t, q$t, p$mean, q$mean, p$cov, q$cov))))
    aci_abort("aci_error_dims",
              "aci_metric: incompatible or non-finite moments.")
  if (lp == 1L) {
    R_p <- as.numeric(p$cov)
    R_q <- as.numeric(q$cov)
    # Preserve aci_metric_pair()'s per-index validation order: reference
    # covariance first, integrating covariance second, before the vectorised
    # arithmetic.
    for (j in seq_len(n)) {
      if (R_q[j] <= 0)
        aci_abort("aci_error_spd", "Matrix (R_q) must be positive definite.")
      if (R_p[j] <= 0)
        aci_abort("aci_error_spd", "Matrix (R_p) must be positive definite.")
    }
    scalar <- .gaussian_kl_scalar_kernel(
      as.numeric(p$mean), R_p, as.numeric(q$mean), R_q,
      decompose = decompose
    )
    return(data.frame(t = p$t, scalar, check.names = FALSE))
  }
  ## Matrix path. The work that survives per index is the two Cholesky factors
  ## and the two triangular solves; the array components, the diagonal indices,
  ## the debug-assert option and the decomposition branch are read once here
  ## rather than at every index. The reference covariance is still factorised
  ## before the integrating one at each index, so the order in which a bad path
  ## is rejected is unchanged. `t.default()` stands in for the generic on the
  ## plain numeric matrices this loop handles, as in the compiled kernels.
  pm <- p$mean; qm <- q$mean; pc <- p$cov; qc <- q$cov
  dm <- c(lp, lp)
  dg <- seq.int(1L, lp * lp, by = lp + 1L)
  dbg <- isTRUE(getOption("aci.debug_assert", FALSE))
  dec <- isTRUE(decompose)
  total <- numeric(n)
  signal <- if (dec) numeric(n) else NULL
  dispersion <- if (dec) numeric(n) else NULL
  for (j in seq_len(n)) {
    Rq <- qc[, , j]; dim(Rq) <- dm
    Rp <- pc[, , j]; dim(Rp) <- dm
    Lq <- .kl_path_chol(Rq, "R_q")
    Lp <- .kl_path_chol(Rp, "R_p")
    tLq <- t.default(Lq)
    w <- forwardsolve(tLq, qm[j, ] - pm[j, ])
    sg <- 0.5 * sum(w * w)
    A <- forwardsolve(tLq, t.default(Lp))
    ds <- 0.5 * (
      sum(A * A) - lp + 2 * sum(log(Lq[dg])) - 2 * sum(log(Lp[dg]))
    )
    sg <- max(sg, 0)
    ds <- max(ds, 0)
    tt <- sg + ds
    if (dbg && tt < -1e-10)
      aci_abort("aci_error_internal", "Negative KL encountered.")
    total[j] <- tt
    if (dec) { signal[j] <- sg; dispersion[j] <- ds }
  }
  if (!dec) return(data.frame(t = p$t, total = total))
  data.frame(t = p$t, total = total, signal = signal,
             dispersion = dispersion)
}


################################################################################
# ACI / aci()
################################################################################

#' Assimilative causal inference
#'
#' Runs the filter and smoother for a model and observed record and scores each
#' time by the relative entropy of the smoother from the filter. A normal
#' `aci()` call uses the supplied-code backward-ODE headline smoother, including
#' its correlated-noise correction, independently of `keep`. [lag_table()] and
#' `aci(table = ...)` instead use the complete online Theorem 3 smoother;
#' their finite-grid diagonal can therefore differ from headline ACI.
#'
#' @param model A `cgns_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param engine One of `"auto"` or `"cgns"`; `"auto"` selects the closed-form
#'   engine for a conditional-Gaussian model.
#' @param conditional Optional `aci_conditional_spec`; see [aci_conditional()].
#' @param table Optional precomputed `lag_table`, whose online diagonal is used
#'   in place of the headline smoother.
#' @param keep One of `"paths"`, `"table"` or `"none"`, selecting which objects
#'   are retained on the result. It does not select a smoother.
#' @param decompose `TRUE` to retain the signal and dispersion parts.
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param regularize Covariance policy for this call; see [aci_filter()]. One
#'   record covers the filter, the smoother and any table this call builds, and
#'   is returned in `meta$regularization`.
#' @param loglik `TRUE` (the default) accumulates the predictive
#'   log-likelihood on the internal filter, where `keep = "paths"` exposes it as
#'   `paths$filter$meta$loglik`. ACI itself never uses it, so `FALSE` skips that
#'   work and leaves `paths$filter$meta$loglik` `NULL`; every ACI quantity is
#'   unchanged.
#' @param ... Must be empty; unused arguments are an error.
#' @returns An object of class `aci_result`.
#'
#' @seealso [aci_range()], [lag_table()], [aci_metric()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' a <- aci(m, ob)
#' a
#'
#' @export
aci <- function(model, obs, engine = c("auto", "cgns"),
                conditional = NULL, table = NULL,
                keep = c("paths", "table", "none"), decompose = TRUE,
                init = NULL,
                stepper = c("explicit", "implicit"), nsub = 1L,
                regularize = NULL, loglik = TRUE, ...) {
  ## R partial-matches `m =` to `model` before positional filling, so
  ## aci(mod, ob, m = 50) silently becomes model = 50, obs = mod. `m` is
  ## the ensemble size of the deferred ensemble engine, not an argument
  ## here; sys.call() still sees the literal name match.call() rewrites.
  if ("m" %in% names(sys.call()))
    aci_abort("aci_error_not_implemented", paste(
      "aci() has no argument m in this release; the ensemble engine and",
      "its ensemble size m are out of scope. R would partial-match m to",
      "model, silently shifting the other arguments."))
  dots <- list(...)
  stepper <- match.arg(stepper); keep <- match.arg(keep)
  policy <- .aci_regularize(regularize)
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose))
    aci_abort("aci_error_dims", "decompose must be TRUE or FALSE.")
  engine <- match.arg(engine)
  obs <- as_obs(obs)
  if (!inherits(model, "stochastic_model") || obs$k != model$k)
    aci_abort("aci_error_model_contract",
              "model must be a stochastic_model matching the observation dimension.")
  if (engine == "auto") engine <- if (inherits(model, "cgns_model")) "cgns" else
    aci_abort("aci_error_model_contract",
              "Only cgns_model is supported in this release; the ensemble engine is out of scope.")
  if (!is.null(table)) {
    if (length(dots))
      aci_abort("aci_error_dims", "Unused arguments were supplied with a reused lag table.")
    .validate_lag_table_source(table, model, obs, conditional, init)
    table_init <- table$meta$init %||% init
    return(structure(list(t = table$t, aci = table$diag,
                          signal = if (decompose) table$diag_signal else NULL,
                          dispersion = if (decompose) table$diag_dispersion else NULL,
                          paths = NULL, table = table,
                          handles = list(model = model, obs = obs,
                                         conditional = table$meta$conditional,
                                         init = table_init),
                          meta = list(engine = "reused_table",
                                      conditional = table$meta$conditional,
                                      ## the discretization scheme, with the
                                      ## finer route tag as the fallback a
                                      ## table saved before `scheme` existed
                                      ## still answers with.
                                      smoother_scheme =
                                        table$meta$scheme %||%
                                        table$meta$reference_smoother %||%
                                        "unspecified_lag_table_reference",
                                      ## no recursion runs here: the record is
                                      ## the one the table's own rows were
                                      ## built under.
                                      regularization =
                                        table$meta$regularization %||%
                                        .aci_reg_freeze(
                                          .aci_reg_new(policy, obs$t)))),
                     class = "aci_result"))
  }
  bundle <- NULL
  rec <- .aci_reg_new(policy, obs$t)
  if (engine == "cgns") {
    if (length(dots))
      aci_abort("aci_error_dims", "Unused arguments were supplied to the closed-form ACI engine.")
    if (!inherits(model, "cgns_model"))
      aci_abort("aci_error_model_contract", "engine='cgns' requires a cgns_model.")
    bundle <- .compile_cgns_run(model, obs, conditional)
    filt <- .cgns_filter_compiled(
      bundle, init = init, stepper = stepper, nsub = nsub, validate = FALSE,
      loglik = loglik, regularize = rec
    )
    smoo <- .cgns_smoother_compiled(bundle, filt, validate = FALSE,
                                    regularize = rec)
  }
  klp <- .gaussian_kl_path_compiled(
    bundle, smoo, filt, decompose = decompose, validate = FALSE
  )
  tab <- NULL
  if (keep == "table") {
    tab <- .lag_table_compiled(
      bundle, mode = "forward", filter = filt, init = filt$meta$init,
      validate = FALSE, regularize = rec
    )
  }
  actual_init <- filt$meta$init
  structure(list(t = klp$t, aci = klp$total,
                 signal = if (decompose) klp$signal else NULL,
                 dispersion = if (decompose) klp$dispersion else NULL,
                 paths = if (keep != "none") list(filter = filt, smoother = smoo) else NULL,
                 table = tab,
                 handles = list(model = model, obs = obs,
                                conditional = conditional,
                                init = actual_init),
                 meta = list(engine = engine, conditional = conditional,
                             smoother_scheme = smoo$meta$scheme %||% "unspecified",
                             table_reference = if (!is.null(tab))
                               tab$meta$reference_smoother else NULL,
                             regularization = .aci_reg_freeze(rec))),
            class = "aci_result")
}


#' Print aci() result
#'
#' @param x An `aci_result` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.aci_result <- function(x, ...) {
  pk <- which.max(x$aci)
  cat(sprintf("<aci_result> engine = %s%s | peak ACI = %.4g at t = %.4g\n",
              x$meta$engine,
              if (!is.null(x$meta$conditional)) " (conditional)" else "",
              x$aci[pk], x$t[pk]))
  invisible(x)
}


#' aci() result to data.frame
#'
#' @param x An `aci_result` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns A data frame with one row per time, carrying the metric and, when
#'   retained, its signal and dispersion parts.
#' @export
as.data.frame.aci_result <- function(x, ...) {
  d <- data.frame(t = x$t, aci = x$aci)
  if (!is.null(x$signal)) { d$signal <- x$signal; d$dispersion <- x$dispersion }
  d
}


################################################################################
# forward causal influence ranges (exact + l1_linf)
################################################################################

#' Composite Simpson's rule on a unit-spaced grid (internal)
#'
#' @param y Numeric vector of samples on a uniformly spaced grid.
#' @returns 1-length numeric, the integral in units of the sample spacing.
#' @noRd
.simpson <- function(y, close = c("quadratic", "trapezoid")) {
  close <- match.arg(close)
  n <- length(y)
  if (n < 3L) {
    return(if (n == 2L) 0.5 * sum(y) else 0)
  }
  m <- if (n %% 2L == 1L) n else n - 1L
  i <- seq(1L, m - 2L, by = 2L)
  s <- sum(y[i] + 4 * y[i + 1L] + y[i + 2L]) / 3
  ## Composite Simpson needs an odd point count; an even-length grid leaves
  ## one interval over.  ACI_code-main/simps.m:106-115 closes it by fitting a
  ## quadratic through the last three points and integrating that over the
  ## final interval, which reduces on a uniform grid to (5 y_n + 8 y_{n-1} -
  ## y_{n-2}) / 12.  The "trapezoid" close is the package's pre-0.0.21 rule,
  ## retained so earlier results stay reproducible; the two differ by
  ## -(y_n - 2 y_{n-1} + y_{n-2}) / 12, zero on a flat or linear tail.
  if (n %% 2L == 0L)
    s <- s + if (identical(close, "trapezoid")) 0.5 * (y[n - 1L] + y[n]) else
      (5 * y[n] + 8 * y[n - 1L] - y[n - 2L]) / 12
  return(s)
}


#' Composite Simpson's rule on an arbitrary grid (internal)
#'
#' The unequal-spacing rule of `ACI_code-main/simps.m:93-115`: a parabola
#' through each consecutive triple, integrated over the pair of intervals it
#' spans, and, when the point count is even, the quadratic through the last
#' three points integrated over the final interval alone. Unlike `.simpson()`,
#' which is the uniform rule and returns the integral in units of the sample
#' spacing, this one carries its abscissae and returns the integral itself.
#'
#' @param x Numeric vector of strictly increasing abscissae.
#' @param y Numeric vector of samples at `x`.
#' @returns 1-length numeric, the integral of `y` over `x`.
#' @noRd
.simpson_xy <- function(x, y) {
  m <- length(x)
  ## simps.m falls back on trapz below three points.
  if (m < 3L)
    return(if (m == 2L) 0.5 * (x[2L] - x[1L]) * (y[1L] + y[2L]) else 0)
  i <- seq.int(1L, m - 2L, by = 2L)
  h1 <- x[i + 1L] - x[i]
  h2 <- x[i + 2L] - x[i + 1L]
  al <- (h1 + h2) / h1 / 6
  s <- sum(al * (2 * h1 - h2) * y[i] +
           al * (h1 + h2)^2 / h2 * y[i + 1L] +
           al * h1 / h2 * (2 * h2 - h1) * y[i + 2L])
  if (m %% 2L == 0L) {
    ## An even point count leaves one interval over.  simps.m closes it with
    ## the quadratic through the last three points, fitted by a Vandermonde
    ## solve and integrated over that interval alone.
    k <- c(m - 2L, m - 1L, m)
    cf <- solve(cbind(x[k]^2, x[k], 1), y[k])
    s <- s + cf[1L] * (x[m]^3 - x[m - 1L]^3) / 3 +
             cf[2L] * (x[m]^2 - x[m - 1L]^2) / 2 +
             cf[3L] * (x[m] - x[m - 1L])
  }
  return(s)
}


#' Resolve the MATLAB-compatibility threshold grid (internal)
#'
#' @param epsilon_grid `NULL` for the reference grid, or strictly
#'   increasing finite
#'   non-negative quadrature nodes.
#' @returns Numeric vector of quadrature nodes, ascending.
#' @noRd
.cir_eps_grid <- function(epsilon_grid) {
  ## ACI_code-main/dyad_interaction_model.m:441-444 takes 513 points linear in
  ## the order of magnitude from -6 to 0.5, and integrates over 10^that in
  ## ascending order.
  if (is.null(epsilon_grid)) return(10^seq(-6, 0.5, length.out = 513L))
  if (!is.numeric(epsilon_grid) || length(epsilon_grid) < 2L ||
      any(!is.finite(epsilon_grid)) || any(epsilon_grid < 0) ||
      is.unsorted(epsilon_grid, strictly = TRUE))
    aci_abort(
      "aci_error_dims",
      paste("epsilon_grid must be at least two finite non-negative",
            "thresholds in strictly increasing order.")
    )
  return(as.numeric(epsilon_grid))
}


#' Settle the MATLAB-compatibility mode's arguments (internal)
#'
#' Keeps the two epsilon roles apart: reporting thresholds reach the result
#' through `epsilon`, quadrature nodes only through `epsilon_grid`, and
#' `epsilon_grid` is
#' accepted only by the mode that has nodes to place.
#'
#' @param quadrature Resolved quadrature name.
#' @param method Resolved objective functional.
#' @param epsilon_grid The caller's `epsilon_grid`, possibly `NULL`.
#' @returns The quadrature nodes, or `NULL` when the mode is not in use.
#' @noRd
.cir_compat_grid <- function(quadrature, method, epsilon_grid) {
  if (identical(quadrature, "matlab_eps_grid")) {
    if (!identical(method, "exact"))
      aci_abort(
        "aci_error_dims",
        paste("quadrature = 'matlab_eps_grid' reproduces the reference",
              "definitional objective and needs method = 'exact'.")
      )
    return(.cir_eps_grid(epsilon_grid))
  }
  if (!is.null(epsilon_grid))
    aci_abort(
      "aci_error_dims",
      paste("epsilon_grid supplies quadrature nodes for",
            "quadrature = 'matlab_eps_grid' only; use epsilon for reporting",
            "thresholds.")
    )
  return(NULL)
}


#' Calculate tau from one row of divergences (internal)
#'
#' Tau is the causal influence range at one anchor time.
#'
#' @param dt Positive 1-length numeric step.
#' @param p Numeric vector of finite, non-negative one-lag divergences.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param M Max p or a bias-corrected metric.
#' @param direction Character, "fwd" or "bwd".
#' @param quadrature One of `"simpson"`, `"sum"` or `"matlab_eps_grid"`.
#' @param simpson_close Even-grid Simpson closing convention.
#' @param epsilon_grid Threshold nodes for `"matlab_eps_grid"`.
#' @returns 1-length numeric `tau`
#' @noRd
.calc_tau <- function(dt, p, method, M, direction, quadrature,
                      simpson_close = "quadratic", epsilon_grid = NULL){
  if (method != "l1_linf") {
    ## The definitional objective: the subjective range averaged over every
    ## threshold.  Reading the range off the running extremum makes that
    ## average a finite sum with no quadrature error at all, because the
    ## count of cells above a threshold, integrated over the threshold, is
    ## identically the sum of the running extremum.  So the definitional
    ## objective on the time grid is dt * sum(suffix max) / M, exactly.
    if(direction == "fwd"){
      p <- rev(cummax(rev(p)))
    } else{ # direction = "bwd"
      p <- rev(cummin(rev(p)))
    }
    if (quadrature == "matlab_eps_grid") {
      ## MATLAB compatibility: the reference defn_objective_CIR, a Simpson
      ## quadrature of the subjective read-out over a threshold grid rather
      ## than the exact average above.  The count of cells strictly above a
      ## node is the position of the last exceedance, so one sorted search
      ## over the non-decreasing reversed running extremum answers the whole
      ## grid.  This read-out counts cells (the reference convention) whatever
      ## `convention` the subjective matrix is reported in, so that the
      ## compatibility objective stays comparable with the exact one.
      counts <- length(p) - findInterval(epsilon_grid, rev(p))
      return(dt * .simpson_xy(epsilon_grid, counts) / M)
    }
    return(dt * sum(p) / M)
  }
  ## The two published conventions for the ratio: the plain andreou2026cir
  ## Appendix G L1 grid-function sum (eqs. G.8/G.14, the FBCIR scripts'
  ## active lines) and composite Simpson (the ACI and EnKBS codebases'
  ## active lines).
  if (quadrature == "sum") {
    tau <- dt * sum(p) / M
  } else {
    tau <- dt * .simpson(p, close = simpson_close) / M
  }
  return(tau)
}


#' Forward influence duration from one table row (internal)
#'
#' @param p Numeric vector of finite, non-negative forward divergences.
#' @param dt Positive 1-length numeric step.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param quadrature One of `"simpson"`, `"sum"` or `"matlab_eps_grid"`.
#' @param simpson_close Even-grid Simpson closing convention.
#' @param epsilon_grid Threshold nodes for `"matlab_eps_grid"`.
#' @returns Named numeric vector with the range `tau` and the strength `M`.
#' @noRd
.fwd_lengths <- function(p, dt, method, quadrature = "simpson",
                        simpson_close = "quadratic", epsilon_grid = NULL) {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p)) || any(p < 0) || length(dt) != 1L || !is.finite(dt) || dt <= 0){
    aci_abort("aci_error_dims", "Forward discrepancy values must be finite and non-negative, with dt > 0.")
  }
  M <- max(p)
  if (!is.finite(M) || M < 1e-14){ return(c(tau = NA_real_, M = M)) }
  # the l1_linf ratio is integrated with composite Simpson, following the ACI
  # reference code; the FBCIR scripts' active ratio lines use the plain
  # andreou2026cir Appendix G L1 sum, with Simpson retained there as a
  # commented alternative.  The definitional objective is a finite sum.
  tau <- .calc_tau(p = p, dt = dt, method = method, M = M, direction = "fwd",
                   quadrature = quadrature, simpson_close = simpson_close,
                   epsilon_grid = epsilon_grid)
  return(c(tau = tau, M = M))
}


#' Resolve the minimum influence strength (internal)
#'
#' @param min_M Either `"auto"` (reads `aci.cir_min_strength` option) or a finite non-negative number.
#' @returns 1-length numeric strength threshold.
#' @noRd
.cir_min_strength <- function(min_M) {
  if (is.character(min_M) && length(min_M) == 1L &&
      !is.na(min_M) && identical(min_M, "auto")) {
    min_M <- getOption("aci.cir_min_strength", 1e-5)
  }
  if (!is.numeric(min_M) || length(min_M) != 1L ||
      !is.finite(min_M) || min_M < 0)
    aci_abort("aci_error_dims", "min_M must be 'auto' or one finite non-negative number.")
  as.numeric(min_M)
}


#' Resolve the requested anchor times (internal)
#'
#' @param anchors `NULL` for every anchor time, or whole numbers indexing the
#'   anchor grid.
#' @param N1 Number of anchor times available.
#' @returns Integer vector of anchor indices, in the order requested.
#' @noRd
.cir_anchors <- function(anchors, N1) {
  if (is.null(anchors)) return(seq_len(N1))
  if (!is.numeric(anchors) || !length(anchors) || any(!is.finite(anchors)) ||
      any(anchors != floor(anchors)) || any(anchors < 1) ||
      any(anchors > N1) || anyDuplicated(anchors))
    aci_abort("aci_error_dims", sprintf(
      "anchors must be distinct whole numbers in 1:%d, or NULL for all of them.",
      N1))
  as.integer(anchors)
}


#' Resolution status of one forward-CIR row (internal)
#'
#' Four outcomes, checked in this order.
#'
#' `"insufficient"`: fewer than three observations from the anchor onwards, so
#' no range is supported at all.
#'
#' `"below_threshold"`: the peak divergence is under the strength floor - the
#' same predicate that masks the range.
#'
#' `"censored"`: the exceedance runs into the end of the record. The last cell
#' of a complete row is identically zero, because the online estimate given the
#' whole record is the smoother it is scored against, so a range can never
#' formally reach the end of the row and testing the last cell alone would
#' never fire (measured: on the pinned 2001-point dyad record the last
#' exceedance of the 1e-5 floor is never later than the second-to-last cell,
#' at any anchor). The test is instead between two parts of the same row: the
#' record left after the last exceedance, against the exceedance itself. A
#' record that does not outlast the influence it measured did not resolve it,
#' and the reported range is a lower bound. No margin or horizon is taken from
#' the caller; the comparison is scale-free and entirely row-derived.
#'
#' `"resolved"`: everything else.
#'
#' @param p Numeric vector, one zero-padded forward row.
#' @param M Peak divergence of that row.
#' @param floor_M Strength floor.
#' @returns 1-length character.
#' @noRd
.cir_status <- function(p, M, floor_M) {
  n <- length(p)
  if (n < 3L) return("insufficient")
  if (M < floor_M) return("below_threshold")
  above <- which(p > floor_M)
  settled <- if (length(above)) above[length(above)] else 0L
  if (2L * settled > n) return("censored")
  return("resolved")
}


#' Construct a causal influence range (internal)
#'
#' @param t Numeric vector of anchor times.
#' @param tau Numeric vector of ranges, in the model's time units.
#' @param M Numeric vector of influence strengths.
#' @param direction Either `"forward"` or `"backward"`.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param dt Positive 1-length numeric step.
#' @param interval Two-column matrix of range endpoints.
#' @param tail_bound Optional heuristic tail estimates carried from the table.
#' @param subjective Optional subjective ranges at the requested thresholds.
#' @param status Optional per-anchor resolution status factor.
#' @param meta Optional named list of metadata carried on the object.
#' @param bound Optional label recording the estimator's discrete semantics.
#' @returns An object of class `cir_result`.
#' @noRd
new_cir_result <- function(t, tau, M, direction, method, dt, interval,
                           tail_bound = NULL, subjective = NULL, status = NULL,
                           meta = list(), bound = NULL) {
  structure(list(t = t, tau = tau, M = M, interval = interval,
                 direction = direction, method = method,
                 bound = bound %||% switch(method, l1_linf =
                   if (direction == "forward") "underestimate" else "overestimate",
                   exact = "layer_cake_objective"),
                 tail_bound = tail_bound, subjective = subjective,
                 status = status, dt = dt, meta = meta),
            class = "cir_result")
}


#' Levels of the per-anchor resolution status, in reporting order (internal)
#' @noRd
.cir_status_levels <- c("resolved", "censored", "below_threshold",
                        "insufficient")


#' Print the causal influence range
#'
#' @param x A `cir_result` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.cir_result <- function(x, ...) {
  cat(sprintf("<cir_result> %s | method = %s (%s) | masked/NA: %d of %d\n",
              x$direction, x$method, x$bound, sum(is.na(x$tau)), length(x$tau)))
  ab <- x$meta$above_baseline
  if (is.null(ab) && !is.null(x$meta$per_reference))
    ab <- vapply(x$meta$per_reference,
                 function(mm) isTRUE(mm$above_baseline), logical(1))
  if (!is.null(ab) && length(ab))
    cat(sprintf("  above initial-time baseline: %d of %d anchors\n",
                sum(ab), length(ab)))
  if (!is.null(x$status) && length(x$status)) {
    tb <- table(x$status)
    tb <- tb[tb > 0L]
    cat(sprintf("  status: %s\n",
                paste(sprintf("%s %d", names(tb), as.integer(tb)),
                      collapse = ", ")))
  }
  invisible(x)
}


#' Causal influence range
#'
#' Summarizes the duration of influence on the discrete time grid, forward from
#' each anchor time. A finite adaptive table is labelled
#' `objective_on_truncated_table`; its `tail_bound` field is a heuristic tail
#' estimate and must not be interpreted as a certified error bound. The
#' `l1_linf` estimators are ratios, integrated with composite Simpson over the
#' whole span, following the ACI reference code.
#'
#' Only `direction = "forward"` is in this release. The backward range is a
#' `FBCIR_code-main` feature: it appears in no ACI_code script, and it is held
#' in the development reserve (`reserve/fbcir/`) for the release that brings
#' that family in. `direction = "backward"` raises `aci_error_not_implemented`
#' rather than returning a forward answer under a backward label.
#'
#' @param x A `lag_table` or `aci_result` object.
#' @param direction `"forward"`, the only direction in this release. See above.
#' @param ... Arguments passed to methods.
#' @returns An object of class `cir_result` with `direction` `"forward"`.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' @seealso [lag_table()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' aci_range(tb)
#'
#' @export
aci_range <- function(x, direction = c("forward", "backward"), ...)
  UseMethod("aci_range")


#' Resolve and gate the range direction (internal)
#'
#' @param direction The `direction` argument as the caller supplied it.
#' @returns `"forward"`; called for its error condition.
#' @noRd
.aci_range_direction <- function(direction) {
  direction <- match.arg(direction, c("forward", "backward"))
  if (identical(direction, "backward"))
    aci_abort("aci_error_not_implemented", paste(
      "The backward causal influence range is not in this release: it is a",
      "FBCIR_code-main feature with no ACI_code counterpart, and is held in",
      "the development reserve (reserve/fbcir/) for the release that brings",
      "that family in."))
  return(direction)
}


#' @describeIn aci_range Forward range from a precomputed lag table.
#' @param method The objective functional. `"exact"` is the definitional
#'   objective range, the subjective range averaged over every threshold;
#'   reading the range off the running maximum makes that average a finite
#'   sum, `dt * sum(suffix max) / M`, with no quadrature error. `"l1_linf"` is
#'   the efficient ratio the ACI reference script computes, `dt * integral(row)
#'   / M`. They are different functionals, not two quadratures of one: they
#'   coincide only where the divergence decreases with lag.
#' @param epsilon Optional numeric vector of finite non-negative
#'   thresholds at which
#'   the subjective range is also reported. Its read-out convention is set by
#'   `convention`. This vector reports thresholds and nothing else; the
#'   MATLAB-compatibility quadrature takes its nodes from `epsilon_grid`, a
#'   separate argument, so that one vector is never asked to do both jobs.
#' @param min_M Either `"auto"` or one finite non-negative number; ranges whose
#'   strength falls below it are masked. It also sets the strength floor the
#'   `status` vocabulary is judged against, and `min_M = NULL` turns off the
#'   masking but not the status, which then falls back to the
#'   `aci.cir_min_strength` option.
#' @param masked_value What a masked range reports: `"na"` (the default) keeps
#'   the masking visible as `NA`; `"zero"` follows andreou2026cir Remarks B.4
#'   and C.4 and the FBCIR scripts, which set the length to 0. The paper's
#'   forward figures are computed on an untruncated table; the adaptive
#'   truncation of the default `lag_table()` is a package storage device whose
#'   use is recorded in the result's `bound` label, and `tol = 0` reproduces
#'   the untruncated convention.
#' @param quadrature How a row is reduced to its objective range, recorded in
#'   `meta$quadrature`. For `method = "l1_linf"`: `"simpson"` (the default)
#'   follows the ACI and EnKBS codebases' active lines, and `"sum"` is the
#'   literal andreou2026cir eq. G.8 L1 grid-function sum, which the FBCIR
#'   scripts' active lines use. For `method = "exact"` the threshold average is
#'   a finite sum, so there is no time-axis quadrature to choose and
#'   `"simpson"` and `"sum"` return the same exact value; `"matlab_eps_grid"`
#'   instead reproduces the reference script's `defn_objective_CIR`, a Simpson
#'   quadrature of the subjective read-out over the threshold grid
#'   `epsilon_grid`.
#'   `"matlab_eps_grid"` is defined for `method = "exact"` only.
#' @param simpson_close Closing rule for the leftover interval when a grid has
#'   an even number of points: `"quadratic"` (the default) fits a quadratic
#'   through the last three points, following `simps.m` in the ACI reference
#'   code; `"trapezoid"` is the package's pre-0.0.21 rule and reproduces
#'   results reported by earlier versions. Odd-length grids are unaffected.
#' @param anchors Optional whole numbers indexing the anchor times to report;
#'   `NULL` (the default) reports every anchor time. Only the requested rows
#'   are formed and reduced, so a handful of anchors on a long record costs a
#'   fraction of the whole-record computation. The result is subset to the
#'   anchors asked for, in the order asked for, and records them in
#'   `meta$anchors`.
#' @param convention The subjective read-out. `"count"` (the default) is the
#'   ACI reference script's `subj_CIR_idx * dt`, the 1-based index of the last
#'   exceedance times the step. `"lag_time"` is andreou2026cir eq. G.7, the lag
#'   time of that exceedance with the first cell at lag 0, one grid step below.
#'   Cells with no exceedance report 0 under both. The objective range is not
#'   affected: `"exact"` and `"matlab_eps_grid"` both average the counting
#'   read-out, which is the convention the definition is written in.
#' @param epsilon_grid Threshold quadrature nodes for
#'   `quadrature = "matlab_eps_grid"`, strictly increasing. `NULL` (the
#'   default) uses the reference script's 513-point grid, `10^seq(-6, 0.5,
#'   length.out = 513)`. Supplying it with any other `quadrature` is an error,
#'   so the reporting thresholds in `epsilon` and these quadrature nodes
#'   can never
#'   be taken from one vector.
#' @export
aci_range.lag_table <- function(x, direction = c("forward", "backward"),
                                method = c("exact", "l1_linf"),
                                epsilon = NULL, min_M = "auto",
                                masked_value = c("na", "zero"),
                                quadrature = c("simpson", "sum",
                                "matlab_eps_grid"),
                                simpson_close = c("quadratic", "trapezoid"),
                                anchors = NULL,
                                convention = c("count", "lag_time"),
                                epsilon_grid = NULL,
                                ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims",
              "Unused arguments were supplied to aci_range().")
  .aci_range_direction(direction)
  masked_value <- match.arg(masked_value)
  quadrature <- match.arg(quadrature)
  simpson_close <- match.arg(simpson_close)
  convention <- match.arg(convention)
  method <- match.arg(method)
  grid <- .cir_compat_grid(quadrature, method, epsilon_grid)
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  if (!x$mode %in% c("forward", "full"))
    aci_abort("aci_error_dims",
              "aci_range needs a 'forward' or 'full' mode table.")
  if (!is.null(epsilon) &&
      (!is.numeric(epsilon) || !length(epsilon) ||
       any(!is.finite(epsilon)) || any(epsilon < 0)))
    aci_abort("aci_error_dims",
              "epsilon must contain finite non-negative thresholds.")
  N1 <- length(x$t)
  idx <- .cir_anchors(anchors, N1)
  n <- length(idx)
  ## The strength floor drives both the mask and the status vocabulary. It is
  ## resolved even when masking is off, so that `status` keeps a scale.
  floor_M <- .cir_min_strength(if (is.null(min_M)) "auto" else min_M)
  structurally_truncated <- any(vapply(idx, function(j)
    length(x$rows[[j]]) < (N1 - j + 1L), logical(1)))
  tau <- M <- rep(NA_real_, n)
  status <- character(n)
  subj <- if (is.null(epsilon)) NULL else
    matrix(0, nrow = n, ncol = length(epsilon))
  ## idx * dt counts cells, the ACI reference script's subj_CIR_idx * dt;
  ## (idx - 1) * dt is the lag time of the last exceedance, andreou2026cir
  ## eq. G.7 with the first cell at lag 0, one grid step below.
  off <- if (identical(convention, "count")) 0 else 1
  for (i in seq_len(n)) {
    ## One pass per requested anchor: the row is formed once and every
    ## quantity read off it, rather than once per threshold as well.
    p <- lt_row(x, idx[i], pad = "zero")
    r <- .fwd_lengths(p, x$dt, method, quadrature, simpson_close, grid)
    tau[i] <- r[["tau"]]; M[i] <- r[["M"]]
    status[i] <- .cir_status(p, M[i], floor_M)
    if (!is.null(subj)) for (q in seq_along(epsilon)) {
      hit <- which(p > epsilon[q])
      subj[i, q] <- if (!length(hit)) 0 else x$dt * (max(hit) - off)
    }
  }
  if (!is.null(min_M)) {
    mm <- .cir_min_strength(min_M)
    bad <- which(M < mm)
    ## The published convention (andreou2026cir Rmks B.4 and C.4, and the
    ## FBCIR scripts) sets a low-strength length to 0; the package default
    ## reports NA so that masking stays visible in downstream summaries.
    if (length(bad)) { tau[bad] <- if (masked_value == "zero") 0 else NA_real_
      aci_warn("aci_warn_low_signal", sprintf(
        "%d forward CIR values masked (M < %.3g); interpret CIRs jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).",
        length(bad), mm)) }
  }
  tt <- x$t[idx]
  new_cir_result(tt, tau, M, "forward", method, x$dt,
                 interval = cbind(tt, tt + tau, deparse.level = 0L),
                 tail_bound = x$tailbnd[idx], subjective = subj,
                 status = factor(status, levels = .cir_status_levels),
                 meta = list(quadrature = quadrature,
                             simpson_close = simpson_close,
                             convention = convention,
                             epsilon_grid = grid,
                             table_conditional = x$meta$conditional,
                             anchors = idx,
                             structurally_truncated = structurally_truncated,
                             ## Reducing stored rows forms no covariance, so
                             ## this method has no policy of its own; it
                             ## reports the one the table was built under.
                             regularization = x$meta$regularization),
                 bound = .cir_bound(method, quadrature,
                                    structurally_truncated))
}


#' Label the estimator's discrete semantics (internal)
#'
#' @param method Resolved objective functional.
#' @param quadrature Resolved quadrature name.
#' @param truncated `TRUE` when some reported row was not stored in full.
#' @returns 1-length character, or `NULL` for the default label.
#' @noRd
.cir_bound <- function(method, quadrature, truncated) {
  if (isTRUE(truncated))
    return(if (method == "exact") "objective_on_truncated_table"
           else "lower_ratio_on_truncated_table_only")
  if (method == "exact" && quadrature == "matlab_eps_grid")
    return("eps_grid_objective")
  return(NULL)
}


#' @describeIn aci_range Forward range from an ACI result, building the table
#'   when it was not retained. The covariance policy is not re-resolved here:
#'   it is read off the result's own `meta$regularization`, so re-analysing a
#'   saved result in a new session cannot silently change it.
#' @export
aci_range.aci_result <- function(x, direction = c("forward", "backward"),
                                 ...) {
  .aci_range_direction(direction)
  if (!is.null(x$table)) return(aci_range(x$table, ...))
  bundle <- .compile_cgns_run(
    x$handles$model, x$handles$obs, x$handles$conditional
  )
  .forward_cir_compiled(
    bundle, filter = x$paths$filter, init = x$handles$init,
    regularize = x$meta$regularization$policy, ...
  )
}
