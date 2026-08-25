################################################################################
## assimilation.R - closed-form engine: filter, smoother, lag table, conditioning, regularization
## ########################################################################## ##
##
## Contents:
##   * closed-form CGNS filter/smoother and the Theorem 3 lag table:
##       - da_filter, da_smooth, new_da_path, print.da_path_gaussian,
##         as.data.frame.da_path_gaussian, .default_init,
##         da_filter.cgns_model, da_smooth.cgns_model, .chol_or_floor,
##         .kl_fast, .onelag_stats, lag_table, lt_diag, lt_onelag, lt_row,
##         lt_tail_bound, truncation_profile, lt_contraction_certificate,
##         print.lag_table, as.data.frame.lag_table
##
##   * conditional ACI strategies (prescribed forcing / masked innovations):
##       - nontarget, print.nontarget_spec, .nt_indices, reduce_nontarget, .resolve_nontarget
##
##   * Gaspari-Cohn localization and inflation:
##       - gaspari_cohn, localization_spec, apply_inflation
##
################################################################################


################################################################################
# closed-form CGNS filter/smoother and the Theorem 3 lag table
################################################################################

#' Data assimilation filter
#'
#' Generic reconstructing the hidden state from the observed record up to each
#' time. The closed-form method is used for a `cgns_model`; a general
#' `stochastic_model` is filtered by the ensemble engine.
#'
#' @param model A `cgns_model` or `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param ... Arguments passed to methods.
#' @returns An assimilation path: `da_path_gaussian` for the closed-form engine,
#'   `da_path_ensemble` for the ensemble engine.
#'
#' @seealso [da_smooth()], [aci()], [lag_table()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- da_filter(m, ob)
#' f
#'
#' @export
da_filter <- function(model, obs, ...) UseMethod("da_filter")


#' Data assimilation smoother
#'
#' Generic reconstructing the hidden state from the whole observed record. The
#' closed-form method is used for a `cgns_model`; a general `stochastic_model`
#' is smoothed by the ensemble engine.
#'
#' @param model A `cgns_model` or `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param ... Arguments passed to methods.
#' @returns An assimilation path: `da_path_gaussian` for the closed-form engine,
#'   `da_path_ensemble` for the ensemble engine.
#'
#' @seealso [da_filter()], [aci()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- da_filter(m, ob)
#' da_smooth(m, ob, filter = f)
#'
#' @export
da_smooth <- function(model, obs, ...) UseMethod("da_smooth")


#' Construct a Gaussian assimilation path (internal)
#'
#' @param t Numeric vector of times.
#' @param mean Numeric matrix of hidden means, one row per time.
#' @param cov Numeric array of hidden covariances, `l` by `l` by time.
#' @param kind 1-length character, `"filter"` or `"smoother"`.
#' @param meta Optional named list of metadata carried on the object.
#' @returns An object of class `da_path_gaussian`.
#' @noRd
new_da_path <- function(t, mean, cov, kind, meta = list()) {
  structure(list(t = t, mean = mean, cov = cov, kind = kind, meta = meta),
            class = c("da_path_gaussian", "da_path"))
}


#' Print a Gaussian assimilation path
#'
#' @param x A `da_path_gaussian` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.da_path_gaussian <- function(x, ...) {
  cat(sprintf("<da_path_gaussian> kind = %s, l = %d, N+1 = %d\n",
              x$kind, ncol(x$mean), length(x$t))); invisible(x)
}


#' Coerce a Gaussian assimilation path to a data frame
#'
#' @param x A `da_path_gaussian` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns A data frame in long form with one row per time and hidden
#'   component, carrying the mean and its marginal variance.
#' @export
as.data.frame.da_path_gaussian <- function(x, ...) {
  l <- ncol(x$mean)
  data.frame(t = rep(x$t, l), var = rep(paste0("y", 1:l), each = length(x$t)),
             mean = as.vector(x$mean),
             sd = as.vector(t(apply(x$cov, 3, function(R) sqrt(pmax(diag(R), 0))))))
}


#' Default diffuse prior for the hidden state (internal)
#'
#' @param m A `cgns_model` object.
#' @param co0 Coefficient list at the first observation, from `eval_coefs()`.
#' @returns A list with the initial `mean` and `cov`.
#' @noRd
.default_init <- function(m, co0) {
  # Weakly-informative default: ~5x the OU stationary scale gyy / (2|Ly|).
  # A 100x-diffuse prior destabilizes the explicit-Euler Riccati when Lx is
  # strongly state-dependent (R^2 Lx' gxx^{-1} Lx dt >> R); see SPEC-01 s8.
  c0 <- 5 * max(mean(diag(co0$gyy)), 1e-8) /
        max(2 * mean(abs(diag(as.matrix(co0$Ly)))), 1e-4)
  list(mean = rep(0, m$l), cov = diag(c0, m$l), diffuse = TRUE)
}

#' Check a supplied Gaussian path against its model and observations (internal)
#'
#' @param path A `da_path_gaussian` object.
#' @param ob An `obs_traj` object.
#' @param l Hidden dimension the path must carry.
#' @param label 1-length character naming the path in error messages.
#' @param nontarget Optional non-target tag the path must have been built with.
#' @param model Optional reduced model the path must match.
#' @param source_model Optional original model the path must match.
#' @returns Invisibly `TRUE`; called for its error conditions.
#' @noRd
.validate_gaussian_path <- function(path, ob, l, label, nontarget = NULL,
                                    model = NULL, source_model = NULL) {
  if (!inherits(path, "da_path_gaussian"))
    aci_abort("aci_error_dims", sprintf("%s must be a da_path_gaussian object.", label))
  if (!identical(path$kind, label))
    aci_abort("aci_error_dims",
              sprintf("%s has kind '%s', not '%s'.", label,
                      path$kind %||% "<missing>", label))
  N1 <- length(ob$t)
  if (length(path$t) != N1 || any(!is.finite(path$t)) ||
      max(abs(path$t - ob$t)) > 1e-10 * max(1, max(abs(ob$t))))
    aci_abort("aci_error_dims", sprintf("%s uses a different time grid.", label))
  if (!identical(dim(path$mean), c(N1, l)) ||
      !identical(dim(path$cov), c(l, l, N1)) ||
      any(!is.finite(c(path$mean, path$cov))))
    aci_abort("aci_error_dims", sprintf("%s has incompatible or non-finite moments.", label))
  for (j in seq_len(N1))
    .strict_chol(path$cov[, , j], sprintf("%s covariance at index %d", label, j))
  if (!identical(path$meta$nontarget %||% NULL, nontarget %||% NULL))
    aci_abort("aci_error_nontarget",
              sprintf("%s was computed under a different nontarget specification.", label))
  stored_obs <- path$meta$obs_x
  if (is.null(stored_obs))
    aci_abort("aci_error_dims",
              sprintf("%s lacks observation provenance; recompute it with the package.", label))
  if (!identical(dim(stored_obs), dim(ob$x)) ||
      any(abs(stored_obs - ob$x) > 1e-12 * pmax(1, abs(ob$x))))
    aci_abort("aci_error_dims",
              sprintf("%s was computed from different observation values.", label))
  # A prescribed-forcing reduction contains newly-created lookup closures on
  # every call, so two valid reductions of the same source model are not
  # `identical()`. Public paths therefore carry the stable original model;
  # fall back to the resolved model only for reduced/internal paths.
  expected_model <- source_model %||% model
  stored_model <- path$meta$source_model %||% path$meta$model
  if (!is.null(expected_model) && is.null(stored_model))
    aci_abort("aci_error_model_contract",
              sprintf("%s lacks model provenance; recompute it with the package.", label))
  if (!is.null(expected_model) && !identical(stored_model, expected_model))
    aci_abort("aci_error_model_contract",
              sprintf("%s was computed with a different model object.", label))
  invisible(TRUE)
}


#' Compare two hidden-state priors for equality (internal)
#'
#' @param a First prior, a list with `mean` and `cov`.
#' @param b Second prior, a list with `mean` and `cov`.
#' @param l Hidden dimension.
#' @returns 1-length logical, `TRUE` when the two priors agree.
#' @noRd
.same_gaussian_init <- function(a, b, l) {
  if (is.null(a) || is.null(b) || is.null(a$cov) || is.null(b$cov)) return(FALSE)
  am <- as.numeric(a$mean %||% rep(0, l)); bm <- as.numeric(b$mean %||% rep(0, l))
  ac <- as.matrix(a$cov); bc <- as.matrix(b$cov)
  length(am) == l && length(bm) == l &&
    identical(dim(ac), c(l, l)) && identical(dim(bc), c(l, l)) &&
    max(abs(am - bm)) <= 1e-12 * max(1, max(abs(c(am, bm)))) &&
    max(abs(ac - bc)) <= 1e-12 * max(1, max(abs(c(ac, bc))))
}


#' @describeIn da_filter Closed-form filter for a conditional-Gaussian model.
#' @param init Optional list with the initial hidden `mean` and `cov`; `NULL`
#'   uses a diffuse prior and warns.
#' @param nontarget Optional `nontarget_spec` selecting a conditional ACI
#'   reduction; see [nontarget()].
#' @param stepper Either `"explicit"` or `"implicit"`; the implicit Riccati step
#'   preserves positivity.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @export
da_filter.cgns_model <- function(model, obs, init = NULL, nontarget = NULL,
                                 stepper = c("explicit", "implicit"),
                                 nsub = 1L, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to da_filter().")
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  bundle <- .compile_cgns_run(model, obs, nontarget)
  .cgns_filter_compiled(
    bundle, init, stepper = match.arg(stepper), nsub = nsub,
    validate = FALSE
  )
}


#' @describeIn da_smooth Closed-form backward-ODE smoother for a
#'   conditional-Gaussian model.
#' @param filter Optional precomputed filter path; recomputed when `NULL`.
#' @param nontarget Optional `nontarget_spec` selecting a conditional ACI
#'   reduction; see [nontarget()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @export
da_smooth.cgns_model <- function(model, obs, filter = NULL, nontarget = NULL,
                                 init = NULL,
                                 stepper = c("explicit", "implicit"),
                                 nsub = 1L, ...) {
  stepper_supplied <- !missing(stepper); nsub_supplied <- !missing(nsub)
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to da_smooth().")
  stepper <- match.arg(stepper)
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  bundle <- .compile_cgns_run(model, obs, nontarget)
  if (!is.null(filter))
    .validate_gaussian_path(
      filter, bundle$obs, bundle$l, "filter", bundle$nontarget,
      model = bundle$model, source_model = model
    )
  if (!is.null(filter) && !is.null(init) &&
      !.same_gaussian_init(init, filter$meta$init, bundle$l))
    aci_abort("aci_error_dims",
              "init conflicts with the prior stored on the supplied filter.")
  if (!is.null(filter) && stepper_supplied &&
      !identical(filter$meta$stepper, stepper))
    aci_abort("aci_error_stepper", "The supplied filter uses a different stepper.")
  if (!is.null(filter) && nsub_supplied &&
      !identical(as.integer(filter$meta$nsub), as.integer(nsub)))
    aci_abort("aci_error_stepper", "The supplied filter uses a different nsub value.")
  filt <- filter %||% .cgns_filter_compiled(
    bundle, init, stepper = stepper, nsub = nsub, validate = FALSE
  )
  p <- .cgns_smoother_compiled(bundle, filt, validate = FALSE)
  # invariant #3: terminal condition equals filter terminal (asserted)
  stopifnot(max(abs(p$mean[length(p$t), ] - filt$mean[length(p$t), ])) < 1e-12)
  p
}


#' Cholesky factor with a floor-on-failure fallback (internal)
#'
#' A lean alternative to [safe_chol()] for the innermost loops, where the
#' guarded stack dominates the cost.
#'
#' @param R Square numeric matrix.
#' @returns A list with the possibly floored matrix `R` and its Cholesky factor
#'   `ch`.
#' @noRd
.chol_or_floor <- function(R) {
  ch <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(ch)) { R <- spd_floor(R); ch <- chol(R) }
  list(R = R, ch = ch)
}


#' Gaussian relative entropy from a precomputed Cholesky factor (internal)
#'
#' @param mu_p Numeric vector, mean of the first distribution.
#' @param chp Cholesky factor of the covariance of the first distribution.
#' @param mu_q Numeric vector, mean of the second distribution.
#' @param Rq Covariance matrix of the second distribution.
#' @returns 1-length numeric, the relative entropy of `p` from `q`, floored at
#'   zero.
#' @noRd
.kl_fast <- function(mu_p, chp, mu_q, Rq) {
  cq <- .chol_or_floor(Rq)
  d <- mu_q - mu_p
  w <- forwardsolve(t(cq$ch), d)
  A <- forwardsolve(t(cq$ch), t(chp))
  max(0.5 * (sum(w * w) + sum(A * A) - length(mu_p)) +
      sum(log(diag(cq$ch))) - sum(log(diag(chp))), 0)
}


#' One-lag smoothed statistics (internal)
#'
#' @param co Coefficient list from `eval_coefs()`.
#' @param aux Theorem-3 auxiliary matrices `E` and `F`.
#' @param muf_nm1 Filtered mean at the earlier index.
#' @param Rf_nm1 Filtered covariance at the earlier index.
#' @param muf_n Smoothed or filtered mean at the later index.
#' @param Rf_n Smoothed or filtered covariance at the later index.
#' @param dx Observed increment across the step.
#' @param dt Positive 1-length numeric step.
#' @param l Hidden dimension.
#' @returns A list with the updated `mu` and `R`.
#' @noRd
.onelag_stats <- function(co, aux, muf_nm1, Rf_nm1, muf_n, Rf_n, dx, dt, l) {
  ILy <- diag(l) + co$Ly * dt
  b   <- muf_nm1 - drop(aux$E %*% (drop(ILy %*% muf_nm1) + co$fy * dt)) +
         drop(aux$F %*% (dx - (drop(co$Lx %*% muf_nm1) + co$fx) * dt))
  Pt  <- Rf_nm1 - aux$E %*% ILy %*% Rf_nm1 - aux$F %*% co$Lx %*% Rf_nm1 * dt
  mus <- drop(aux$E %*% muf_n) + b
  Rs  <- spd_floor(sym(aux$E %*% Rf_n %*% t(aux$E) + Pt))
  list(mu = mus, R = Rs)
}


#' Finite-lag divergence table
#'
#' Builds the table of finite-lag divergences consumed by the causal
#' influence range estimators. A lag table uses the complete online Theorem 3
#' smoother as its reference. That reference costs O(N) time-point work; table
#' construction then costs work proportional to the retained lag cells, with
#' O(N^2) cells for a full table in the worst case.
#'
#' @param model A `cgns_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param mode One of `"forward"`, `"one_lag"` or `"full"`, selecting which
#'   cells are retained.
#' @param tol Positive tolerance below which a row is frozen by the adaptive
#'   storage rule.
#' @param window Number of consecutive steps a row must stay below `tol` before
#'   it is frozen.
#' @param max_lag Maximum positive lag retained, or `Inf` for no cap.
#' @param filter Optional precomputed filter path.
#' @param smoother Optional precomputed smoother path.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param ... Must be empty; unused arguments are an error.
#' @returns An object of class `lag_table`.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
#' **36**(4), 71. \doi{10.1007/s00332-026-10271-x}
#'
#' @seealso [forward_cir()], [backward_cir()], [lt_row()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' lag_table(m, ob, mode = "forward")
#'
#' @export
lag_table <- function(model, obs, mode = c("forward", "one_lag", "full"),
                      tol = getOption("aci.default_tol", 1e-8), window = 3L,
                      max_lag = Inf, filter = NULL, smoother = NULL,
                      nontarget = NULL, init = NULL,
                      stepper = "explicit", nsub = 1L, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to lag_table().")
  if (!is.character(stepper) || length(stepper) != 1L ||
      !identical(stepper, "explicit") || !is.numeric(nsub) ||
      length(nsub) != 1L || !is.finite(nsub) || nsub != 1)
    aci_abort("aci_error_stepper",
              "lag_table requires stepper = 'explicit' and nsub = 1.")
  mode <- match.arg(mode)
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_not_implemented",
              paste(
                "lag_table() is the closed-form CGNS route; use",
                "ensemble_lag_table() or aci(..., engine = 'ensemble',",
                "keep = 'table') for the andreou2026cir/jiang2026enkbs forward ensemble route."))
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    aci_abort("aci_error_dims", "tol must be one finite non-negative number.")
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window < 1L || window != floor(window))
    aci_abort("aci_error_dims", "window must be a positive integer.")
  if (!is.numeric(max_lag) || length(max_lag) != 1L || is.na(max_lag) ||
      max_lag < 1 || (!is.infinite(max_lag) && max_lag != floor(max_lag)))
    aci_abort("aci_error_dims", "max_lag must be a positive integer or Inf.")
  if (mode == "full" && is.finite(max_lag))
    aci_abort("aci_error_dims", "mode = 'full' requires max_lag = Inf.")
  bundle <- .compile_cgns_run(model, obs, nontarget)
  .lag_table_compiled(
    bundle, mode = mode, tol = tol, window = window, max_lag = max_lag,
    filter = filter, smoother = smoother, init = init, validate = FALSE
  )
}


#' Check a lag table against the model and observations it is reused with (internal)
#'
#' @param table A `lag_table` object.
#' @param model A `cgns_model` object.
#' @param obs An `obs_traj` object.
#' @param nontarget Optional non-target tag the table must have been built with.
#' @param init Optional prior the table must have been built with.
#' @returns Invisibly `TRUE`; called for its error conditions.
#' @noRd
.validate_lag_table_source <- function(table, model, obs,
                                       nontarget = NULL, init = NULL) {
  if (!inherits(table, "lag_table"))
    aci_abort("aci_error_dims", "table must be a lag_table.")
  obs <- as_obs(obs)
  if (length(table$t) != length(obs$t) ||
      max(abs(table$t - obs$t)) > 1e-10 * max(1, max(abs(obs$t))))
    aci_abort("aci_error_dims", "The supplied lag table uses a different observation grid.")
  sx <- table$meta$source_obs_x
  if (is.null(sx) || !identical(dim(sx), dim(obs$x)) ||
      any(abs(sx - obs$x) > 1e-12 * pmax(1, abs(obs$x))))
    aci_abort("aci_error_dims",
              "The supplied lag table was computed from different observation values.")
  if (is.null(table$meta$source_model) ||
      !identical(table$meta$source_model, model))
    aci_abort("aci_error_model_contract",
              "The supplied lag table was computed with a different model object.")
  if (!identical(nontarget %||% NULL, table$meta$nontarget %||% NULL))
    aci_abort("aci_error_nontarget",
              "The supplied lag table used a different nontarget specification.")
  if (!is.null(init) &&
      !.same_gaussian_init(init, table$meta$init, model$l))
    aci_abort("aci_error_dims",
              "The supplied lag table was computed with a different initialization.")
  invisible(TRUE)
}


#' Diagonal of a lag table
#'
#' Accesses the zero-lag entries of a table without depending on its storage
#' representation.
#'
#' @param x A `lag_table` object.
#' @returns Numeric vector of the zero-lag divergences, one per time.
#'
#' @seealso [lt_row()], [lt_onelag()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' head(lt_diag(tb))
#'
#' @export
lt_diag <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  x$diag
}


#' One-lag sequence of a lag table
#'
#' Accesses the single-observation-lag entries used by the backward causal
#' influence range.
#'
#' @param x A `lag_table` object.
#' @returns Numeric vector of the one-lag divergences.
#'
#' @seealso [backward_cir()], [lt_row()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "one_lag")
#' head(lt_onelag(tb))
#'
#' @export
lt_onelag <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  out <- x$onelag %||% aci_abort("aci_error_dims", "Table was not built in 'one_lag' mode.")
  # A one-lag discrepancy exists for the N intervals only. The internal N+1
  # allocation has a terminal zero sentinel used by the sweep; exposing that
  # sentinel makes every suffix minimum zero and corrupts exact backward CIR.
  if (length(out) == length(x$t)) out <- out[-length(out)]
  out
}


#' One row of a lag table
#'
#' Accesses the divergences at one anchor time across increasing positive lag,
#' padding the cells that adaptive storage did not retain.
#'
#' @param x A `lag_table` object.
#' @param j Integer index of the anchor time.
#' @param pad Either `"zero"` or `"na"`, the value used for cells the table did
#'   not retain.
#' @returns Numeric vector of divergences at increasing positive lag.
#'
#' @seealso [lt_diag()], [forward_cir()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' head(lt_row(tb, 1))
#'
#' @export
lt_row <- function(x, j, pad = c("zero", "na")) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  pad <- match.arg(pad)
  if (is.null(x$rows)) aci_abort("aci_error_dims", "Table has no stored rows (mode = 'one_lag').")
  if (length(j) != 1L || !is.numeric(j) || !is.finite(j) ||
      j != as.integer(j) || j < 1L || j > length(x$t))
    aci_abort("aci_error_dims", "j must be one valid integer row index.")
  j <- as.integer(j)
  r <- x$rows[[j]]; want <- length(x$t) - j + 1
  if (length(r) < want) r <- c(r, rep(if (pad == "zero") 0 else NA_real_, want - length(r)))
  r
}


#' Heuristic tail estimate of a lag table
#'
#' The historical `lt_tail_bound()` name is retained for compatibility, but its
#' value is a heuristic tail estimate, not a certified mathematical error bound.
#'
#' @param x A `lag_table` object.
#' @param j Optional integer index of a single anchor time; `NULL` returns the
#'   estimate at every time.
#' @returns Numeric vector of heuristic tail estimates.
#'
#' @seealso [truncation_profile()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' head(lt_tail_bound(tb))
#'
#' @export
lt_tail_bound <- function(x, j = NULL) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  if (is.null(x$tailbnd))
    aci_abort("aci_error_dims", "This lag table mode has no truncation estimates.")
  if (is.null(j)) return(x$tailbnd)
  if (length(j) != 1L || !is.numeric(j) || !is.finite(j) ||
      j != as.integer(j) || j < 1L || j > length(x$t))
    aci_abort("aci_error_dims", "j must be one valid integer row index.")
  x$tailbnd[as.integer(j)]
}


#' Retained-cell profile of a lag table
#'
#' Reports how many lag cells adaptive storage retained at each anchor time.
#'
#' @param x A `lag_table` object.
#' @returns A data frame with one row per anchor time, giving the retained cell
#'   count and the corresponding lag reach.
#'
#' @seealso [lt_tail_bound()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' head(truncation_profile(tb))
#'
#' @export
truncation_profile <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  if (is.null(x$L) || is.null(x$tailbnd))
    aci_abort("aci_error_dims", "This lag table mode has no truncation profile.")
  data.frame(j = seq_along(x$t), t = x$t, L = x$L,
             tail_estimate = x$tailbnd)
}


#' Per-step contraction certificate for the Theorem 3 update
#'
#' Evaluates the checkable per-step contraction condition of the online
#' smoother update matrices (andreou2026smoother eqs. 3.18-3.19): when the
#' Hermitian part of `(I - E^j) / dt` is positive definite at every step,
#' the spectral radius of every `E^j` falls below one for a small enough
#' step. Reported per step: `lambda_min`, the smallest eigenvalue of that
#' Hermitian part (positive at the step means condition 3.18 holds there);
#' `enorm`, the operator 2-norm of `E^j`; and `rho_E`, its spectral radius.
#'
#' The condition is stated for the continuous-time generator and holds for
#' `dt` sufficiently small, so a positive `lambda_min` alongside an `enorm`
#' at or above one flags a step-size margin rather than a broken model. The
#' spectral radius is not submultiplicative, so no bound on the accumulated
#' update products follows from `rho_E < 1` alone; [lt_tail_bound()] remains
#' a heuristic estimate.
#'
#' @param x A `lag_table` object.
#' @returns A data frame with one row per grid step, carrying `j`, `t`,
#'   `lambda_min`, `enorm` and `rho_E`, with attributes `gamma` (the largest
#'   `enorm`) and `condition_318` (`TRUE` when every `lambda_min` is
#'   positive).
#'
#' @references
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
#' **36**(4), 71. \doi{10.1007/s00332-026-10271-x}
#'
#' @seealso [lag_table()], [lt_tail_bound()], [truncation_profile()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' tb <- lag_table(m, as_obs(sim), mode = "forward",
#'                 init = list(mean = 0, cov = diag(1, 1)))
#' cert <- lt_contraction_certificate(tb)
#' attr(cert, "gamma"); attr(cert, "condition_318")
#'
#' @export
lt_contraction_certificate <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  model <- x$meta$source_model
  if (is.null(model) || is.null(x$meta$source_obs_x))
    aci_abort("aci_error_dims",
              "This lag table carries no model/observation handles.")
  obs <- observed_trajectory(x$t, x$meta$source_obs_x)
  bundle <- .compile_cgns_run(model, obs, x$meta$nontarget)
  filt <- .cgns_filter_compiled(
    bundle, x$meta$init, stepper = "explicit", nsub = 1L,
    validate = FALSE
  )
  l <- bundle$l; dt <- x$dt; N1 <- length(x$t)
  lam <- en <- rh <- numeric(N1 - 1L)
  for (j in seq_len(N1 - 1L)) {
    co <- .compiled_co(bundle, j)
    aux <- .thmD1_aux_compiled(bundle, j, filt$cov[, , j], co = co)
    Hm  <- (diag(l) - aux$E) / dt
    Hm  <- (Hm + t(Hm)) / 2
    lam[j] <- min(eigen(Hm, symmetric = TRUE, only.values = TRUE)$values)
    en[j]  <- if (l == 1L) abs(aux$E[1, 1]) else norm(aux$E, "2")
    rh[j]  <- if (l == 1L) abs(aux$E[1, 1]) else
      max(abs(eigen(aux$E, only.values = TRUE)$values))
  }
  out <- data.frame(j = seq_len(N1 - 1L), t = x$t[seq_len(N1 - 1L)],
                    lambda_min = lam, enorm = en, rho_E = rh)
  attr(out, "gamma") <- max(en)
  attr(out, "condition_318") <- all(lam > 0)
  out
}


#' Print a lag table
#'
#' @param x A `lag_table` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.lag_table <- function(x, ...) {
  cat(sprintf("<lag_table> mode = %s, N+1 = %d, tol = %g\n", x$mode,
              length(x$t), x$meta$tol))
  if (!is.null(x$L)) cat(sprintf("  mean retained lag: %.1f steps; max heuristic tail estimate: %.2e\n",
                                 mean(x$L, na.rm = TRUE), max(x$tailbnd)))
  invisible(x)
}


#' Coerce a lag table to a data frame
#'
#' @param x A `lag_table` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns A data frame in long form with one row per retained cell, carrying
#'   the anchor time, the lag and the divergence.
#' @export
as.data.frame.lag_table <- function(x, ...) {
  if (is.null(x$rows)) {
    ol <- lt_onelag(x)
    return(data.frame(j = seq_along(ol), t = x$t[seq_along(ol)], onelag = ol))
  }
  do.call(rbind, lapply(seq_along(x$rows), function(j)
    data.frame(j = j, n = j + seq_along(x$rows[[j]]) - 1, P = x$rows[[j]])))
}

################################################################################
# Conditional ACI strategies (prescribed forcing & masked innovations)
################################################################################

#' Conditional ACI non-target specification
#'
#' Describes which observed channels are excluded from the causal question, and
#' by which strategy. `nontarget()` describes conditional ACI masking;
#' [reduce_nontarget()] is the prescribed-forcing reduction used when requested.
#'
#' @param blocks Integer or character vector naming the non-target observed
#'   channels.
#' @param strategy Either `"prescribed_forcing"`, which substitutes the
#'   non-target channels as known forcing, or `"inflate"`, which gives their
#'   innovations zero weight in the filter.
#' @returns An object of class `nontarget_spec`.
#'
#' @seealso [reduce_nontarget()], [da_filter()]
#'
#' @examples
#' nontarget(blocks = 2, strategy = "inflate")
#'
#' @export
nontarget <- function(blocks, strategy = c("prescribed_forcing", "inflate")) {
  strategy <- match.arg(strategy)
  structure(list(blocks = blocks, strategy = strategy), class = "nontarget_spec")
}


#' Print a non-target specification
#'
#' @param x A `nontarget_spec` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.nontarget_spec <- function(x, ...) {
  cat(sprintf("<nontarget_spec> x_B = {%s}, strategy = %s\n",
              paste(x$blocks, collapse = ", "), x$strategy)); invisible(x)
}


#' Resolve non-target channel names to indices (internal)
#'
#' @param spec A `nontarget_spec` object.
#' @param obs An `obs_traj` object.
#' @returns A list with the non-target indices and their target complement.
#' @noRd
.nt_indices <- function(spec, obs) {
  b <- spec$blocks
  if (is.character(b)) {
    if (!length(b) || anyNA(b) || any(!nzchar(b)) || anyDuplicated(b))
      aci_abort("aci_error_nontarget",
                "Named non-target blocks must be unique, non-empty names.")
    nm <- colnames(obs$x)
    if (is.null(nm)) aci_abort("aci_error_nontarget", "Named blocks need named obs columns.")
    b <- match(b, nm)
  } else {
    if (!is.numeric(b) || !length(b) || any(!is.finite(b)) ||
        any(b != floor(b)) || anyDuplicated(b))
      aci_abort("aci_error_nontarget",
                "Numeric non-target blocks must be unique finite integer indices.")
  }
  if (any(is.na(b)) || any(b < 1) || any(b > obs$k))
    aci_abort("aci_error_nontarget", "non-target blocks out of range.")
  idxB <- sort(as.integer(b)); idxA <- setdiff(seq_len(obs$k), idxB)
  if (length(idxA) == 0)
    aci_abort("aci_error_nontarget", "At least one observed target component (x_A) is required.")
  list(A = idxA, B = idxB)
}


#' Reduce a model by prescribing the non-target channels
#'
#' Rewrites a model so that the non-target observed channels enter as known
#' forcing rather than as observations to be assimilated, leaving the target
#' channels as the only observed process.
#' x_B becomes prescribed forcing in the (x_A, y) system.
#'
#' @param model A `cgns_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param spec A `nontarget_spec` object; see [nontarget()].
#' @returns A list with the reduced `model` and the reduced `obs`.
#'
#' @seealso [nontarget()]
#'
#' @examples
#' m2 <- cgns_model(
#'   Lx = function(t, x) matrix(c(1, 0), 2, 1),
#'   fx = function(t, x) -0.5 * x,
#'   Ly = function(t, x) matrix(-0.5, 1, 1),
#'   fy = function(t, x) 0,
#'   Sx1 = function(t, x) diag(0.5, 2),
#'   Sy2 = function(t, x) matrix(1, 1, 1),
#'   k = 2, l = 1)
#' sim2 <- simulate(m2, seed = 1, T = 1, dt = 0.01)
#' reduce_nontarget(m2, as_obs(sim2), nontarget(2))
#'
#' @export
reduce_nontarget <- function(model, obs, spec) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_nontarget", "prescribed_forcing requires a cgns_model (use 'inflate' otherwise).")
  obs <- as_obs(obs)
  if (!inherits(spec, "nontarget_spec") ||
      !identical(spec$strategy, "prescribed_forcing"))
    aci_abort("aci_error_nontarget",
              "reduce_nontarget() requires nontarget(..., strategy = 'prescribed_forcing').")
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  ix <- .nt_indices(spec, obs); A <- ix$A; B <- ix$B
  # Admissibility (SPEC-01 s6 S1): the A-B cross-block must vanish along
  # the entire prescribed path, not only at its initial point.
  bad_cross <- any(vapply(seq_along(obs$t), function(j) {
    co <- eval_coefs(model, obs$t[j], obs$x[j, ])
    max(abs(co$gxx[A, B, drop = FALSE])) >
      1e-12 * max(abs(co$gxx), 1e-300)
  }, logical(1)))
  if (bad_cross)
    aci_abort("aci_error_nontarget_crossnoise",
      "gxx has a nonzero A-B cross-block; use nontarget(strategy = 'inflate') (SPEC-01 s6, pending SI equivalence transcription).")
  t0 <- obs$t[1]; dt <- obs$dt; XB <- obs$x[, B, drop = FALSE]; Nrow <- nrow(XB)
  lookupB <- function(t) {
    j <- as.integer(round((t - t0) / dt)) + 1L
    XB[min(max(j, 1L), Nrow), ]
  }
  assemble <- function(t, xA) { full <- numeric(length(A) + length(B))
    full[A] <- xA; full[B] <- lookupB(t); full }
  mk <- function(fun, rows = NULL) function(t, xA) {
    v <- fun(t, assemble(t, xA)); if (is.null(rows)) v else
      if (is.matrix(v)) v[rows, , drop = FALSE] else v[rows] }
  red <- cgns_model(
    Lx = mk(model$Lx, A), fx = mk(model$fx, A),
    Ly = mk(model$Ly),    fy = mk(model$fy),
    Sx1 = mk(model$Sx1, A), Sx2 = mk(model$Sx2, A),
    Sy1 = mk(model$Sy1),    Sy2 = mk(model$Sy2),
    k = length(A), l = model$l,
    name = paste0(model$name, "|reduced"),
    meta = utils::modifyList(model$meta,
                             list(nontarget_reduction = list(target = A,
                                                             prescribed = B))))
  if (!is.null(model$meta$ic_default)) {
    ic <- model$meta$ic_default
    if (is.list(ic) && length(ic$x0) == model$k && length(ic$y0) == model$l &&
        all(is.finite(c(ic$x0, ic$y0)))) {
      red$meta$ic_default <- list(x0 = as.numeric(ic$x0)[A],
                                  y0 = as.numeric(ic$y0))
    } else {
      red$meta$ic_default <- NULL
    }
  }
  list(model = red,
       obs = observed_trajectory(obs$t, obs$x[, A, drop = FALSE]),
       map = ix)
}


# Resolve a nontarget spec into what the CGNS kernels need:
#   list(model, obs, ginv) where ginv(gxx) is the (possibly masked) Gram inverse.
#' Apply a non-target strategy before assimilation (internal)
#'
#' @param model A `cgns_model` object.
#' @param obs An `obs_traj` object.
#' @param spec Optional `nontarget_spec` object; `NULL` leaves the model alone.
#' @returns A list with the resolved `model`, `obs`, Gram inverse `ginv`,
#'   `likelihood_idx` and the metadata `tag`.
#' @noRd
.resolve_nontarget <- function(model, obs, spec) {
  if (is.null(spec))
    return(list(model = model, obs = obs,
                ginv = function(gxx) chol_solve(gxx, diag(nrow(gxx)), "gxx"),
                likelihood_idx = seq_len(obs$k), tag = NULL))
  if (!inherits(spec, "nontarget_spec"))
    aci_abort("aci_error_nontarget", "nontarget must be created by nontarget().")
  if (spec$strategy == "prescribed_forcing") {
    red <- reduce_nontarget(model, obs, spec)
    return(list(model = red$model, obs = red$obs,
                ginv = function(gxx) chol_solve(gxx, diag(nrow(gxx)), "gxx"),
                likelihood_idx = seq_len(red$obs$k), tag = spec))
  }
  ix <- .nt_indices(spec, obs)
  list(model = model, obs = obs,
       ginv = function(gxx) masked_ginv(gxx, ix$A),
       likelihood_idx = ix$A, tag = spec)
}


################################################################################
# Gaspari-Cohn localization and inflation
################################################################################

#' Gaspari-Cohn localization weight
#'
#' Evaluates the compactly supported fifth-order piecewise-rational correlation
#' function used to taper ensemble covariances with distance.
#'
#' @param r Numeric vector of non-negative distances, scaled by the
#'   localization radius.
#' @returns Numeric vector of weights in the unit interval, zero at scaled
#'   distances beyond two.
#'
#' @references
#' Gaspari, G. and Cohn, S. E. (1999). Construction of correlation functions
#' in two and three dimensions. *Quarterly Journal of the Royal Meteorological
#' Society* **125**, 723-757. \doi{10.1002/qj.49712555417}
#'
#' @seealso [localization_spec()]
#'
#' @examples
#' gaspari_cohn(seq(0, 2.5, by = 0.5))
#'
#' @export
gaspari_cohn <- function(r) {
  if (!is.numeric(r) || any(!is.finite(r)))
    aci_abort("aci_error_dims", "r must contain finite numeric distances.")
  attrs <- attributes(r)
  r <- abs(r); out <- numeric(length(r))
  i1 <- r < 1; i2 <- r >= 1 & r < 2
  z <- r[i1]
  out[i1] <- 1 - (5/3) * z^2 + (5/8) * z^3 + (1/2) * z^4 - (1/4) * z^5
  z <- r[i2]
  out[i2] <- 4 - 5 * z + (5/3) * z^2 + (5/8) * z^3 - (1/2) * z^4 +
             (1/12) * z^5 - 2 / (3 * z)
  out <- pmax(out, 0)
  attributes(out) <- attrs
  out
}


#' Localization specification for the ensemble engine
#'
#' Builds the taper matrix applied to ensemble cross-covariances, from state
#' and observation coordinates and a localization radius.
#'
#' @param coords Numeric vector or matrix of hidden-state coordinates.
#' @param radius Positive 1-length numeric localization radius.
#' @param coords_obs Optional coordinates of the observed channels; `NULL`
#'   reuses `coords`.
#' @param distance Either `"euclidean"` or `"cyclic"`.
#' @param period Positive 1-length numeric domain period, required when
#'   `distance` is `"cyclic"`.
#' @returns A list carrying the taper matrix and the settings used to build it.
#'
#' @seealso [gaspari_cohn()], [enkbf()]
#'
#' @examples
#' localization_spec(coords = 1:5, radius = 2)
#'
#' @export
localization_spec <- function(coords, radius, coords_obs = NULL,
                              distance = c("euclidean", "cyclic"), period = NULL) {
  distance <- match.arg(distance)
  if (length(radius) != 1L || !is.finite(radius) || radius <= 0)
    aci_abort("aci_error_dims", "radius must be a finite positive scalar.")
  if (distance == "cyclic" && is.null(period))
    aci_abort("aci_error_dims", "cyclic distance needs `period`.")
  co <- as.matrix(coords); cobs <- if (is.null(coords_obs)) co else as.matrix(coords_obs)
  if (!nrow(co) || !ncol(co) || ncol(cobs) != ncol(co) ||
      any(!is.finite(c(co, cobs))))
    aci_abort("aci_error_dims", "coords and coords_obs must be finite matrices with matching columns.")
  if (distance == "cyclic" &&
      (length(period) != 1L || !is.finite(period) || period <= 0))
    aci_abort("aci_error_dims", "period must be a finite positive scalar.")
  if (distance == "cyclic" && ncol(co) != 1L)
    aci_abort("aci_error_dims", "cyclic distance currently requires one-dimensional coordinates.")
  dfun <- function(a, b) {
    d <- outer(seq_len(nrow(a)), seq_len(nrow(b)),
               Vectorize(function(i, j) sqrt(sum((a[i, ] - b[j, ])^2))))
    if (distance == "cyclic") {
      d <- outer(a[, 1], b[, 1], function(x, y) {
        dd <- abs((x - y) %% period)
        pmin(dd, period - dd)
      })
    }
    d
  }
  structure(list(radius = radius,
                 C2 = gaspari_cohn(dfun(co, co) / radius),        # l x l (hidden-hidden)
                 C1 = gaspari_cohn(dfun(co, cobs) / radius)),     # l x k (hidden-obs)
            class = "localization_spec")
}


#' Apply multiplicative inflation to ensemble members
#'
#' Scales the anomalies of an ensemble about its mean. The inflation parameter
#' is named for the variance factor, while the anomalies are multiplied by its
#' square root.
#'
#' @param members Numeric matrix of ensemble members.
#' @param delta2 Positive 1-length numeric variance inflation factor.
#' @returns A matrix of the same dimension as `members`, with inflated
#'   anomalies about the unchanged ensemble mean.
#'
#' @seealso [enkbf()]
#'
#' @examples
#' apply_inflation(matrix(rnorm(20), nrow = 2), delta2 = 1.1)
#'
#' @export
apply_inflation <- function(members, delta2) {
  members <- as.matrix(members)
  if (ncol(members) < 2L || any(!is.finite(members)))
    aci_abort("aci_error_dims", "members must be a finite matrix with at least two columns.")
  if (length(delta2) != 1L || !is.finite(delta2) || delta2 < 1)
    aci_abort("aci_error_model_contract", "inflation delta^2 must be a finite scalar >= 1.")
  if (delta2 == 1) return(members)
  yb <- rowMeans(members)
  # jiang2026enkbs eq. 15 anomaly rescaling
  yb + sqrt(delta2) * (members - yb)
}
