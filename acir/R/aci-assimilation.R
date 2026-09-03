################################################################################
## aci-assimilation.R - closed-form filter, smoother, lag table, conditioning
## ########################################################################## ##
##
## Contents:
##   * closed-form CGNS filter/smoother and the Theorem 3 lag table:
##       - aci_filter, aci_smoother, aci_online, new_da_path,
##         print.da_path_gaussian, as.data.frame.da_path_gaussian,
##         .default_init, .attach_da_trusted, .da_filter_authenticated,
##         aci_filter.cgns_model, aci_smoother.cgns_model,
##         aci_online.cgns_model,
##         .kl_fast, .onelag_stats, lag_table, lt_diag, lt_row,
##         lt_tail_bound, print.lag_table, as.data.frame.lag_table
##
##   * conditional ACI strategies (prescribed forcing / masked innovations):
##       - aci_conditional, print.aci_conditional_spec, .nt_indices,
##         .model_estimand_spec,
##         .check_prescribed_grid,
##         aci_conditional_reduce, .resolve_nontarget
##
################################################################################


################################################################################
# closed-form CGNS filter/smoother and the Theorem 3 lag table
################################################################################

#' Data assimilation filter
#'
#' Generic reconstructing the hidden state from the observed record up to each
#' time. The closed-form method is used for a `cgns_model`; a general
#' `stochastic_model` is out of scope in this release.
#'
#' @param model A `cgns_model` or `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param ... Arguments passed to methods.
#' @returns An assimilation path: `da_path_gaussian` for the closed-form engine.
#'   Its `meta$loglik` holds the predictive log-likelihood of the observed
#'   record, or `NULL` when the method was called with `loglik = FALSE`.
#'
#' @seealso [aci_smoother()], [aci()], [lag_table()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- aci_filter(m, ob)
#' f
#'
#' @export
aci_filter <- function(model, obs, ...) UseMethod("aci_filter")


#' Data assimilation smoother
#'
#' Generic reconstructing the hidden state from the whole observed record. The
#' closed-form method is used for a `cgns_model`; a general `stochastic_model`
#' is out of scope in this release.
#'
#' @param model A `cgns_model` or `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param ... Arguments passed to methods.
#' @returns An assimilation path: `da_path_gaussian` for the closed-form engine.
#'
#' @seealso [aci_filter()], [aci()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- aci_filter(m, ob)
#' aci_smoother(m, ob, filter = f)
#'
#' @export
aci_smoother <- function(model, obs, ...) UseMethod("aci_smoother")


#' Fixed-lag online data assimilation
#'
#' Generic reconstructing the hidden state at each time from the observed
#' record up to a fixed number of steps ahead of it. The closed-form method is
#' used for a `cgns_model`; a general `stochastic_model` is out of scope in
#' this release.
#'
#' `lag` is the number of future observations each estimate may condition on:
#' the estimate at index `j` uses the observed record through index `j + lag`,
#' and saturates at the end of the record. `lag = 0` returns the filter
#' moments unchanged. `lag = Inf` returns the complete Theorem 3 posterior
#' given the whole record.
#'
#' @section Scheme:
#' `aci_online()` computes the **discrete** Theorem 3 posterior: the exact
#' conditional law of the hidden state given the observed increments on the
#' sampling grid under the explicit single-step discretization. [aci_smoother()]
#' integrates the **continuous** backward smoothing equations with an Euler
#' step of the same size. These are two discretizations of the same
#' continuous-time object and they agree only to first order in the step, so at
#' full lag `aci_online()` does **not** reproduce [aci_smoother()]. The gap
#' grows
#' with the length of the record; it is not a constant offset. On the packaged
#' ENSO partition (`l = 3`, `dt = 0.005`) the smoothed means differ by up to
#' 1.89e-02 against a mean scale of 0.388 over 401 steps, and by up to 9.58e-02
#' against a scale of 2.22 over 4001 steps. The resulting ACI values differ by
#' up to 0.104 against a scale of 1.093 at 401 steps and 0.482 against 2.347 at
#' 4001. [lag_table()] and the causal influence range estimators use the
#' discrete scheme throughout, which is why `lt_diag()` and `aci()$aci` can
#' differ by that same amount. `meta$scheme` records which scheme produced a
#' path: `"theorem3_discrete"` here and on the lag table's reference smoother,
#' `"backward_ode_euler"` on [aci_smoother()].
#'
#' The one exact boundary is the other end: at `lag = 0` the returned moments
#' are the filter moments, unchanged value for value.
#'
#' @param model A `cgns_model` or `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param lag Number of future steps each estimate may condition on. A
#'   non-negative whole number, or `Inf` for the whole record. No default: the
#'   lag is the argument the function exists for, and defaulting it invites the
#'   full-lag result to be mistaken for [aci_smoother()].
#' @param ... Arguments passed to methods.
#' @returns An assimilation path of kind `"online"`, carrying `meta$lag`, the
#'   per-anchor `meta$lag_effective`, `meta$saturated` and `meta$scheme`. Its
#'   kind is what keeps it out of the places a complete smoother is required:
#'   `lag_table(smoother = )` rejects it with `aci_error_dims`.
#'
#' @references
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
#' **36**(4), 71. \doi{10.1007/s00332-026-10271-x}
#'
#' @seealso [aci_filter()], [aci_smoother()], [lag_table()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' aci_online(m, ob, lag = 5)
#'
#' @export
aci_online <- function(model, obs, lag, ...) UseMethod("aci_online")


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
#' @param conditional Optional conditional tag the path must have been
#'   built with.
#' @param model Optional reduced model the path must match.
#' @param source_model Optional original model the path must match.
#' @returns Invisibly `TRUE`; called for its error conditions.
#' @noRd
.validate_gaussian_path <- function(path, ob, l, label, conditional = NULL,
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
  if (l == 1L) {
    ## For a 1x1 slice `.strict_chol()` reduces exactly to a positivity test.
    ## The dimension check above has already established square, non-empty and
    ## finite (and rejects non-finite moments with `aci_error_dims` before this
    ## point, which is why the finiteness term below is belt and braces rather
    ## than a second gate); `max(abs(R - t(R)))` is 0 for a 1x1 matrix, so the
    ## symmetry check cannot fire; `sym()` returns it unchanged; and `chol()` on
    ## a 1x1 matrix succeeds if and only if its single value is finite and
    ## strictly positive. The accept set is therefore identical, and `test-22`
    ## pins that in both directions on finite, non-finite, zero and negative
    ## values, against `.strict_chol()` itself.
    cvv <- as.numeric(path$cov)
    ok <- is.finite(cvv) & cvv > 0
    if (!all(ok)) {
      j <- which.max(!ok)
      aci_abort("aci_error_spd", sprintf(
        "Matrix (%s covariance at index %d) must be positive definite.",
        label, j))
    }
  } else {
    for (j in seq_len(N1))
      .strict_chol(path$cov[, , j],
                   sprintf("%s covariance at index %d", label, j))
  }
  if (!identical(path$meta$conditional %||% NULL, conditional %||% NULL))
    aci_abort("aci_error_nontarget", sprintf(
      "%s was computed under a different conditional specification.", label))
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


# Private identity used to recognize a filter path this namespace produced.
# Like the CGNS realiser sentinel it is deliberately not reconstructible from
# outside: a path that arrives from another namespace, from `unserialize()`, or
# from a hand-built list carries no live reference to it and is simply
# revalidated in full.
.da_trusted_sentinel <- new.env(parent = emptyenv())


#' Seal a package-produced filter path against its own moments (internal)
#'
#' The token holds the two things the run cannot re-establish from anywhere
#' else: the mean and covariance paths this namespace produced. Everything the
#' smoother also needs to agree on - the time grid, the observations, the
#' non-target tag, the model - already has an authoritative live copy on the
#' compiled bundle or in the caller's own arguments, and is checked against
#' that instead of being duplicated here.
#'
#' It holds references, not copies, so the sealed path allocates nothing new,
#' and R's copy-on-modify makes any later edit to the path visible as a
#' mismatch. The token is pure data apart from the namespace-local sentinel, so
#' two filter paths computed from the same inputs stay `identical()`, which the
#' package's own dispatch contract relies on.
#'
#' A path whose moments are not finite is left unsealed. That is the one
#' precondition of `.strict_chol()` the filter kernels do not already
#' guarantee, so refusing to seal it keeps the two routes exactly equivalent
#' rather than nearly so. The test is on the sums: any non-finite element makes
#' its sum non-finite, and the one false negative, a finite path whose sum
#' overflows, merely declines to seal and is validated in full as before. That
#' is one pass and no allocation, where `all(is.finite())` is a pass plus a
#' logical vector as long as the path.
#'
#' @param path A `da_path_gaussian` filter built by this namespace.
#' @returns `path`, sealed when its moments are finite.
#' @noRd
.attach_da_trusted <- function(path) {
  if (!is.finite(sum(path$mean)) || !is.finite(sum(path$cov))) return(path)
  attr(path, ".aci_trusted") <- list(
    sentinel = .da_trusted_sentinel,
    abi      = 1L,
    mean     = path$mean,
    cov      = path$cov
  )
  path
}


#' Authenticate a supplied filter path against the current run (internal)
#'
#' Two questions, in order. Are these the moments this namespace returned,
#' unedited? And was the path built for the run now being smoothed? Only when
#' both hold may `.validate_gaussian_path()` be skipped; every other outcome,
#' including a token that cannot be read at all, returns `FALSE` and leaves the
#' caller on the unchanged full-validation path.
#'
#' Each check below stands in for one the full validation performs. The
#' per-step Cholesky loop is the one with no cheap live equivalent, and it is
#' what the token exists to retire: every covariance the filter kernels store
#' has passed the covariance policy, so a covariance that is still ours is
#' positive definite because the run would have stopped, or been floored back
#' into the cone on request, otherwise.
#'
#' @param filter A supplied path, of any shape.
#' @param bundle The compiled CGNS bundle for the current run.
#' @param source_model The model the caller passed to [aci_smoother()].
#' @returns 1-length logical.
#' @noRd
.da_filter_authenticated <- function(filter, bundle, source_model) {
  tok <- attr(filter, ".aci_trusted", exact = TRUE)
  if (!is.list(tok) ||
      !identical(tok$sentinel, .da_trusted_sentinel) ||
      !identical(tok$abi, 1L))
    return(FALSE)
  ## the moments are untouched since we built them
  if (!inherits(filter, "da_path_gaussian") ||
      !identical(tok$mean, filter$mean) ||
      !identical(tok$cov, filter$cov))
    return(FALSE)
  ## and the path belongs to this run
  N1 <- as.integer(bundle$N1); l <- as.integer(bundle$l)
  identical(filter$kind, "filter") &&
    identical(dim(filter$mean), c(N1, l)) &&
    identical(dim(filter$cov), c(l, l, N1)) &&
    identical(filter$t, bundle$t) &&
    identical(filter$meta$obs_x, bundle$x) &&
    identical(filter$meta$conditional %||% NULL,
              bundle$conditional %||% NULL) &&
    identical(filter$meta$source_model, source_model)
}


#' @describeIn aci_filter Closed-form filter for a conditional-Gaussian model.
#' @param init Optional list with the initial hidden `mean` and `cov`; `NULL`
#'   uses a diffuse prior and warns.
#' @param conditional Optional `aci_conditional_spec` selecting a
#'   conditional ACI
#'   reduction; see [aci_conditional()].
#' @param stepper Either `"explicit"` or `"implicit"`; the implicit Riccati step
#'   preserves positivity.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param regularize Covariance policy for this call. `"none"` (the default,
#'   and the value of `getOption("aci.regularize")` when it is unset) stops
#'   with a classed `aci_error_covariance_not_spd` naming the site, grid index
#'   and time as soon as a covariance leaves the positive-definite cone.
#'   `"floor"` is the previous behaviour: the covariance is projected back by
#'   [spd_floor()] and every such event is recorded in the result's
#'   `meta$regularization`.
#' @param loglik `TRUE` (the default) accumulates the predictive
#'   log-likelihood into `meta$loglik`. `FALSE` skips that work; the filter
#'   moments are unchanged and `meta$loglik` is `NULL`. The likelihood is not
#'   used by ACI, so `FALSE` is the cheaper choice when only the state estimate
#'   is wanted.
#' @export
aci_filter.cgns_model <- function(model, obs, init = NULL, conditional = NULL,
                                  stepper = c("explicit", "implicit"),
                                  nsub = 1L, regularize = NULL, loglik = TRUE,
                                  ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims",
              "Unused arguments were supplied to aci_filter().")
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  bundle <- .compile_cgns_run(model, obs, conditional)
  .attach_da_trusted(
    .cgns_filter_compiled(
      bundle, init, stepper = match.arg(stepper), nsub = nsub,
      validate = FALSE, loglik = loglik,
      regularize = .aci_regularize(regularize)
    )
  )
}


#' @describeIn aci_smoother Closed-form backward-ODE smoother for a
#'   conditional-Gaussian model.
#' @param filter Optional precomputed filter path; recomputed when `NULL`.
#' @param conditional Optional `aci_conditional_spec` selecting a
#'   conditional ACI
#'   reduction; see [aci_conditional()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param force_validate `FALSE` (the default) lets a `filter` that
#'   [aci_filter()] produced for this same model, observations and conditional
#'   specification, and that has not been altered since, skip the per-step
#'   re-validation of its covariances. Any other supplied path, including one
#'   that has been through `saveRDS()`, is validated in full as before. `TRUE`
#'   validates unconditionally. The smoother result is the same either way.
#' @param regularize Covariance policy for this call; see [aci_filter()]. One
#'   record covers the whole call, so a filter recomputed here and the backward
#'   recursion that consumes it share the `meta$regularization` on the returned
#'   smoother.
#' @export
aci_smoother.cgns_model <- function(model, obs, filter = NULL,
                                    conditional = NULL,
                                    init = NULL,
                                    stepper = c("explicit", "implicit"),
                                    nsub = 1L, regularize = NULL,
                                    force_validate = FALSE, ...) {
  stepper_supplied <- !missing(stepper); nsub_supplied <- !missing(nsub)
  if (length(list(...)))
    aci_abort("aci_error_dims",
              "Unused arguments were supplied to aci_smoother().")
  if (length(force_validate) != 1L || !is.logical(force_validate) ||
      is.na(force_validate))
    aci_abort("aci_error_dims", "force_validate must be TRUE or FALSE.")
  stepper <- match.arg(stepper)
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  bundle <- .compile_cgns_run(model, obs, conditional)
  rec <- .aci_reg_new(.aci_regularize(regularize), bundle$t)
  ## A path this namespace built for this exact run, still untouched, has
  ## already been through every check below; re-deriving N+1 Cholesky factors
  ## of covariances the filter kernel itself floored proves nothing new.
  if (!is.null(filter) &&
      !(isFALSE(force_validate) &&
        .da_filter_authenticated(filter, bundle, model)))
    .validate_gaussian_path(
      filter, bundle$obs, bundle$l, "filter", bundle$conditional,
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
    bundle, init, stepper = stepper, nsub = nsub, validate = FALSE,
    regularize = rec
  )
  p <- .cgns_smoother_compiled(bundle, filt, validate = FALSE, regularize = rec)
  # invariant #3: terminal condition equals filter terminal (asserted)
  stopifnot(max(abs(p$mean[length(p$t), ] - filt$mean[length(p$t), ])) < 1e-12)
  p
}


#' @describeIn aci_online Closed-form fixed-lag online Theorem 3 smoother for a
#'   conditional-Gaussian model.
#' @param filter Optional precomputed filter path; recomputed when `NULL`. It
#'   must be the explicit single-step filter, which is the discretization the
#'   Theorem 3 recursions are exact for.
#' @param conditional Optional `aci_conditional_spec` selecting a
#'   conditional ACI
#'   reduction; see [aci_conditional()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param regularize Covariance policy for this call; see [aci_filter()]. One
#'   record covers the whole call, and is returned in `meta$regularization`.
#' @param force_validate `FALSE` (the default) lets a `filter` that
#'   [aci_filter()] produced for this same run, unaltered since, skip per-step
#'   re-validation, as in [aci_smoother()].
#' @export
aci_online.cgns_model <- function(model, obs, lag, filter = NULL,
                                  conditional = NULL, init = NULL,
                                  regularize = NULL, force_validate = FALSE,
                                  ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims",
              "Unused arguments were supplied to aci_online().")
  if (missing(lag))
    aci_abort("aci_error_dims", "aci_online() requires a lag.")
  if (!is.numeric(lag) || length(lag) != 1L || is.na(lag) || lag < 0 ||
      (is.finite(lag) && lag != floor(lag)))
    aci_abort("aci_error_dims",
              "lag must be one non-negative whole number, or Inf.")
  if (length(force_validate) != 1L || !is.logical(force_validate) ||
      is.na(force_validate))
    aci_abort("aci_error_dims", "force_validate must be TRUE or FALSE.")
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims",
              "Observation dimension does not match the CGNS model.")
  bundle <- .compile_cgns_run(model, obs, conditional)
  rec <- .aci_reg_new(.aci_regularize(regularize), bundle$t)
  if (!is.null(filter) &&
      !(isFALSE(force_validate) &&
        .da_filter_authenticated(filter, bundle, model)))
    .validate_gaussian_path(
      filter, bundle$obs, bundle$l, "filter", bundle$conditional,
      model = bundle$model, source_model = model
    )
  if (!is.null(filter) && !is.null(init) &&
      !.same_gaussian_init(init, filter$meta$init, bundle$l))
    aci_abort("aci_error_dims",
              "init conflicts with the prior stored on the supplied filter.")
  filt <- filter %||% .cgns_filter_compiled(
    bundle, init, stepper = "explicit", nsub = 1L, validate = FALSE,
    loglik = FALSE, regularize = rec
  )
  .da_online_compiled(bundle, filt, lag, regularize = rec)
}


#' @describeIn aci_filter Classed not-implemented condition for a general
#'   (non-CGNS) stochastic model.
#' @export
aci_filter.stochastic_model <- function(model, obs, ...)
  aci_abort("aci_error_not_implemented",
            "Only cgns_model is supported in this release; the ensemble engine is out of scope.")


#' @describeIn aci_smoother Classed not-implemented condition for a general
#'   (non-CGNS) stochastic model.
#' @export
aci_smoother.stochastic_model <- function(model, obs, ...)
  aci_abort("aci_error_not_implemented",
            "Only cgns_model is supported in this release; the ensemble engine is out of scope.")


#' @describeIn aci_online Classed not-implemented condition for a general
#'   (non-CGNS) stochastic model.
#' @export
aci_online.stochastic_model <- function(model, obs, lag, ...)
  aci_abort("aci_error_not_implemented",
            "Only cgns_model is supported in this release; the ensemble engine is out of scope.")


#' Gaussian relative entropy from a precomputed Cholesky factor (internal)
#'
#' @param mu_p Numeric vector, mean of the first distribution.
#' @param chp Cholesky factor of the covariance of the first distribution.
#' @param mu_q Numeric vector, mean of the second distribution.
#' @param Rq Covariance matrix of the second distribution.
#' @param rec Covariance-policy recorder; its `j` names the anchor row the
#'   caller is scoring.
#' @returns 1-length numeric, the relative entropy of `p` from `q`, floored at
#'   zero.
#' @noRd
.kl_fast <- function(mu_p, chp, mu_q, Rq, rec) {
  cq <- .cov_guard_chol(Rq, rec, "metric_reference")
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
#' @param rec Covariance-policy recorder; its `j` names the interval index.
#' @returns A list with the updated `mu` and `R`.
#' @noRd
.onelag_stats <- function(co, aux, muf_nm1, Rf_nm1, muf_n, Rf_n, dx, dt, l,
                          rec = .aci_reg_for(NULL, NA_real_)) {
  ILy <- diag(l) + co$Ly * dt
  b   <- muf_nm1 - drop(aux$E %*% (drop(ILy %*% muf_nm1) + co$fy * dt)) +
         drop(aux$F %*% (dx - (drop(co$Lx %*% muf_nm1) + co$fx) * dt))
  Pt  <- Rf_nm1 - aux$E %*% ILy %*% Rf_nm1 - aux$F %*% co$Lx %*% Rf_nm1 * dt
  mus <- drop(aux$E %*% muf_n) + b
  ## `.cov_guard()` forms the symmetric part itself, exactly as `spd_floor()`
  ## did after `sym()`; the arithmetic is unchanged.
  Rs  <- .cov_guard(aux$E %*% Rf_n %*% t(aux$E) + Pt, rec, "smoother_onelag")
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
#' @param mode Either `"forward"` or `"full"`, selecting which cells are
#'   retained.
#' @param tol Positive tolerance below which a row is frozen by the adaptive
#'   storage rule.
#' @param window Number of consecutive steps a row must stay below `tol` before
#'   it is frozen.
#' @param max_lag Maximum positive lag retained, or `Inf` for no cap.
#' @param filter Optional precomputed filter path.
#' @param smoother Optional precomputed smoother path.
#' @param conditional Optional `aci_conditional_spec`; see [aci_conditional()].
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param regularize Covariance policy for this call; see [aci_filter()]. One
#'   record covers the filter, the Theorem 3 reference smoother and every
#'   relative-entropy denominator the table forms, and is returned in
#'   `meta$regularization`.
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
#' @seealso [aci_range()], [lt_row()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' lag_table(m, ob, mode = "forward")
#'
#' @export
lag_table <- function(model, obs, mode = c("forward", "full"),
                      tol = getOption("aci.default_tol", 1e-8), window = 3L,
                      max_lag = Inf, filter = NULL, smoother = NULL,
                      conditional = NULL, init = NULL,
                      stepper = "explicit", nsub = 1L, regularize = NULL, ...) {
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
              "lag_table() requires a cgns_model.")
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
  bundle <- .compile_cgns_run(model, obs, conditional)
  .lag_table_compiled(
    bundle, mode = mode, tol = tol, window = window, max_lag = max_lag,
    filter = filter, smoother = smoother, init = init, validate = FALSE,
    regularize = .aci_regularize(regularize)
  )
}


#' Check a lag table against the model and observations it is reused with (internal)
#'
#' @param table A `lag_table` object.
#' @param model A `cgns_model` object.
#' @param obs An `obs_traj` object.
#' @param conditional Optional conditional tag the table must have been
#'   built with.
#' @param init Optional prior the table must have been built with.
#' @returns Invisibly `TRUE`; called for its error conditions.
#' @noRd
.validate_lag_table_source <- function(table, model, obs,
                                       conditional = NULL, init = NULL) {
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
  if (!identical(conditional %||% NULL, table$meta$conditional %||% NULL))
    aci_abort("aci_error_nontarget",
              paste("The supplied lag table used a different",
                    "conditional specification."))
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
#' @seealso [lt_row()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
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
#' @seealso [lt_diag()], [aci_range()]
#'
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "forward")
#' head(lt_row(tb, 1))
#'
#' @export
lt_row <- function(x, j, pad = c("zero", "na")) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  pad <- match.arg(pad)
  if (is.null(x$rows)) aci_abort("aci_error_dims", "Table has no stored rows.")
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
#' @examples
#' m <- aci_dyad_model()
#' sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
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
  do.call(rbind, lapply(seq_along(x$rows), function(j)
    data.frame(j = j, n = j + seq_along(x$rows[[j]]) - 1, P = x$rows[[j]])))
}

################################################################################
# Conditional ACI strategies (prescribed forcing & masked innovations)
################################################################################

#' Conditional ACI specification
#'
#' Describes which observed channels carry the causal question and which are
#' conditioned out of it, and by which method. `aci_conditional()` describes
#' conditional ACI masking; [aci_conditional_reduce()] is the model reduction
#' `method = "reduce"` asks for.
#'
#' # Estimand
#'
#' Split the observed process into target and conditioning channels,
#' `x = (x_A, x_B)`, with hidden process `y`. Conditional ACI is the estimand
#' `y(t) -> x_A | x_B`: only the target channels `x_A` transfer information
#' into the hidden posterior, and the conditioning channels `x_B` are
#' conditioned upon rather than assimilated. `"mask"` realises it by giving
#' the `x_B` innovations zero weight in the filter gain, the Riccati term and
#' the online-smoother gain, through an observation-precision matrix supported
#' only on the `A` block; `"reduce"` realises it by rewriting the
#' model so `x_B` enters as a known time series (prescribed forcing). The two
#' coincide when the `A`-`B` noise cross-block vanishes, which the reduction
#' checks along the whole path.
#'
#' `target` names `x_A` directly, which is how the reference scripts write the
#' question (`h_W(t) -> T_C | (u, T_E, tau, I)`,
#' `ENSO_model_cond_ACI_h_W_unobs.m:1199-1202`). `given` names the complement
#' `x_B`. Supply exactly one; the other side is derived.
#'
#' # First-slice convention
#'
#' The reference scripts fill the first slice of the observation-precision
#' array with the full Gram inverse before the target-only overwrite, and mask
#' only the later slices
#' (`ENSO_model_cond_ACI_h_W_unobs.m:1197` against `:1250`). `first_step`
#' selects between masking every slice, `"uniform"`, and reproducing that
#' asymmetry, `"matlab"`. It is not a round-off-level choice: on a 4001-point
#' ENSO path with `h_W` hidden and `T_C` the target it moves the step-2 filter
#' mean by 0.108 and the peak ACI by 0.574, and the difference decays through
#' the record rather than vanishing. It is inert wherever the mask itself is
#' inert. `"matlab"` applies only to `"mask"`, which is where a masked
#' precision path exists.
#'
#' @param given Integer or character vector naming the conditioning observed
#'   channels `x_B`. Supply this or `target`, not both.
#' @param method Either `"mask"`, which gives the conditioning channels'
#'   innovations zero weight in the filter, or `"reduce"`, which substitutes
#'   them as known forcing.
#' @param target Integer or character vector naming the target observed
#'   channels `x_A`. Supply this or `given`, not both.
#' @param first_step Either `"uniform"`, which masks the observation precision
#'   at every step, or `"matlab"`, which leaves the first slice unmasked as the
#'   reference scripts do. `"matlab"` requires `method = "mask"`.
#' @returns An object of class `aci_conditional_spec`.
#'
#' @seealso [aci_conditional_reduce()], [aci_filter()]
#'
#' @examples
#' aci_conditional(given = 2, method = "mask")
#' aci_conditional(target = 1, method = "mask")
#'
#' @export
aci_conditional <- function(given = NULL, method = c("mask", "reduce"),
                            target = NULL,
                            first_step = c("uniform", "matlab")) {
  method     <- match.arg(method)
  first_step <- match.arg(first_step)
  if (is.null(given) && is.null(target))
    aci_abort("aci_error_nontarget",
              paste("Supply the conditioning channels as given, or the",
                    "target channels as target."))
  if (!is.null(given) && !is.null(target))
    aci_abort("aci_error_nontarget",
              paste("given and target name the two sides of one split;",
                    "supply exactly one."))
  if (identical(first_step, "matlab") && !identical(method, "mask"))
    aci_abort("aci_error_nontarget", paste(
      "first_step = 'matlab' is a convention of the masked observation",
      "precision; it requires method = 'mask'."))
  structure(list(given = given, method = method, target = target,
                 first_step = first_step), class = "aci_conditional_spec")
}


#' Print a conditional ACI specification
#'
#' @param x An `aci_conditional_spec` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.aci_conditional_spec <- function(x, ...) {
  side <- if (is.null(x$target))
    sprintf("x_B = {%s}", paste(x$given, collapse = ", ")) else
    sprintf("x_A = {%s}", paste(x$target, collapse = ", "))
  cat(sprintf("<aci_conditional_spec> %s, method = %s%s\n", side, x$method,
              if (identical(x$first_step %||% "uniform", "uniform")) "" else
                sprintf(", first_step = %s", x$first_step)))
  invisible(x)
}


#' Resolve conditioning channel names to indices (internal)
#'
#' @param spec A `aci_conditional_spec` object.
#' @param obs An `obs_traj` object.
#' @returns A list with the conditioning indices and their target complement.
#' @noRd
.nt_indices <- function(spec, obs) {
  resolve <- function(v, label) {
    if (is.character(v)) {
      if (!length(v) || anyNA(v) || any(!nzchar(v)) || anyDuplicated(v))
        aci_abort("aci_error_nontarget",
                  sprintf("Named %s channels must be unique, non-empty names.",
                          label))
      nm <- colnames(obs$x)
      if (is.null(nm)) aci_abort("aci_error_nontarget",
                                 "Named channels need named obs columns.")
      v <- match(v, nm)
    } else {
      if (!is.numeric(v) || !length(v) || any(!is.finite(v)) ||
          any(v != floor(v)) || anyDuplicated(v))
        aci_abort("aci_error_nontarget",
                  sprintf(paste("Numeric %s channels must be unique finite",
                                "integer indices."), label))
    }
    if (any(is.na(v)) || any(v < 1) || any(v > obs$k))
      aci_abort("aci_error_nontarget",
                sprintf("%s channels out of range.", label))
    sort(as.integer(v))
  }
  if (is.null(spec$target)) {
    idxB <- resolve(spec$given, "conditioning")
    idxA <- setdiff(seq_len(obs$k), idxB)
    if (length(idxA) == 0)
      aci_abort("aci_error_nontarget",
                "At least one observed target component (x_A) is required.")
  } else {
    idxA <- resolve(spec$target, "target")
    idxB <- setdiff(seq_len(obs$k), idxA)
    if (length(idxB) == 0)
      aci_abort("aci_error_nontarget",
                paste("At least one observed conditioning component (x_B)",
                      "is required."))
  }
  list(A = idxA, B = idxB)
}


#' Apply a model's declared assimilation estimand (internal)
#'
#' A constructor may record the observation set its estimand is defined on, as
#' an `aci_conditional_spec` in `meta$estimand_nontarget`. Assimilation then
#' runs that
#' specification whenever the caller supplies none. Composing a second
#' specification on top of a declared one is not supported.
#'
#' @param model The model handed to a public entry point.
#' @param spec The caller's `conditional` argument, possibly `NULL`.
#' @returns The specification to compile with.
#' @noRd
.model_estimand_spec <- function(model, spec) {
  declared <- if (inherits(model, "cgns_model"))
    model$meta$estimand_nontarget else NULL
  if (is.null(declared)) return(spec)
  if (!inherits(declared, "aci_conditional_spec"))
    aci_abort("aci_error_nontarget",
              "meta$estimand_nontarget must be built by aci_conditional().")
  if (!is.null(spec))
    aci_abort("aci_error_nontarget", paste(
      "This model already declares its own observation set through",
      "meta$estimand_nontarget, so a second conditional specification",
      "cannot be composed with it. Build the model on the full observation",
      "set instead."))
  declared
}


#' Refuse an observation grid a prescribed forcing was not stored on (internal)
#'
#' A model whose coefficients read constructor-supplied forcings by index
#' records that grid in `meta$prescribed_grid`. Running it on any other grid
#' would silently pair each observation with a neighbouring forcing value, so
#' the grids are compared once, where the model meets the observations, rather
#' than interpolated. Models without such a declaration return immediately.
#'
#' @param model The model handed to a public entry point.
#' @param obs An `obs_traj` object.
#' @returns `NULL`, invisibly; called for its error.
#' @noRd
.check_prescribed_grid <- function(model, obs) {
  g <- if (inherits(model, "cgns_model")) model$meta$prescribed_grid else NULL
  if (is.null(g)) return(invisible(NULL))
  if (!is.list(g) || !all(c("t0", "dt", "n") %in% names(g)))
    aci_abort("aci_error_model_contract",
              "meta$prescribed_grid must record t0, dt and n.")
  tol <- 1e-9 * g$dt
  j <- round((obs$t - g$t0) / g$dt) + 1
  if (length(obs$t) > g$n || any(j < 1) || any(j > g$n) ||
      max(abs(obs$t - (g$t0 + (j - 1) * g$dt))) > tol)
    aci_abort("aci_error_dims", sprintf(paste(
      "The observation times are not the grid this model's prescribed",
      "forcing (%s) was supplied on: t = %.10g to %.10g by %.10g, %d points.",
      "The forcings are looked up by index and never interpolated, so",
      "assimilate on that grid or rebuild the model on this one."),
      paste(g$channels %||% "prescribed", collapse = ", "),
      g$t0, g$t0 + (g$n - 1) * g$dt, g$dt, g$n))
  invisible(NULL)
}


#' Reduce a model by prescribing the conditioning channels
#'
#' Rewrites a model so that the conditioning observed channels enter as known
#' forcing rather than as observations to be assimilated, leaving the target
#' channels as the only observed process.
#' x_B becomes prescribed forcing in the (x_A, y) system.
#'
#' @param model A `cgns_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param spec A `aci_conditional_spec` object; see [aci_conditional()].
#' @returns A list with the reduced `model` and the reduced `obs`.
#'
#' @seealso [aci_conditional()]
#'
#' @examples
#' m2 <- aci_model(
#'   Lx = function(t, x) matrix(c(1, 0), 2, 1),
#'   fx = function(t, x) -0.5 * x,
#'   Ly = function(t, x) matrix(-0.5, 1, 1),
#'   fy = function(t, x) 0,
#'   Sx1 = function(t, x) diag(0.5, 2),
#'   Sy2 = function(t, x) matrix(1, 1, 1),
#'   k = 2, l = 1)
#' sim2 <- simulate(m2, seed = 1, t_end = 1, dt = 0.01)
#' aci_conditional_reduce(m2, as_obs(sim2), aci_conditional(2, "reduce"))
#'
#' @export
aci_conditional_reduce <- function(model, obs, spec) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_nontarget", paste(
      "method = 'reduce' requires an aci_model (use 'mask' otherwise)."))
  obs <- as_obs(obs)
  if (!inherits(spec, "aci_conditional_spec") ||
      !identical(spec$method, "reduce"))
    aci_abort("aci_error_nontarget",
              paste("aci_conditional_reduce() requires",
                    "aci_conditional(..., method = 'reduce')."))
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
      paste("The observation-noise Gram gxx couples the target and",
            "non-target channels; method = 'reduce' needs a vanishing",
            "cross-block, so use aci_conditional(method = 'mask')."))
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
  red <- aci_model(
    Lx = mk(model$Lx, A), fx = mk(model$fx, A),
    Ly = mk(model$Ly),    fy = mk(model$fy),
    Sx1 = mk(model$Sx1, A), Sx2 = mk(model$Sx2, A),
    Sy1 = mk(model$Sy1),    Sy2 = mk(model$Sy2),
    k = length(A), l = model$l,
    name = paste0(model$name, "|reduced"),
    meta = utils::modifyList(model$meta,
                             list(estimand_nontarget = NULL,
                                  nontarget_reduction = list(target = A,
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


# Resolve a conditional spec into what the CGNS kernels need:
#   list(model, obs, ginv) where ginv(gxx) is the (possibly masked) Gram inverse.
#' Apply a conditional ACI method before assimilation (internal)
#'
#' @param model A `cgns_model` object.
#' @param obs An `obs_traj` object.
#' @param spec Optional `aci_conditional_spec` object; `NULL` leaves the model
#'   alone.
#' @returns A list with the resolved `model`, `obs`, Gram inverse `ginv`,
#'   `likelihood_idx` and the metadata `tag`.
#' @noRd
.resolve_nontarget <- function(model, obs, spec) {
  if (is.null(spec))
    return(list(model = model, obs = obs,
                ginv = function(gxx) chol_solve(gxx, diag(nrow(gxx)), "gxx"),
                likelihood_idx = seq_len(obs$k), tag = NULL))
  if (!inherits(spec, "aci_conditional_spec"))
    aci_abort("aci_error_nontarget",
              "conditional must be created by aci_conditional().")
  if (spec$method == "reduce") {
    red <- aci_conditional_reduce(model, obs, spec)
    return(list(model = red$model, obs = red$obs,
                ginv = function(gxx) chol_solve(gxx, diag(nrow(gxx)), "gxx"),
                likelihood_idx = seq_len(red$obs$k), tag = spec))
  }
  ix <- .nt_indices(spec, obs)
  list(model = model, obs = obs,
       ginv = function(gxx) masked_ginv(gxx, ix$A),
       likelihood_idx = ix$A, tag = spec)
}
