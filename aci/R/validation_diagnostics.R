################################################################################
## validation_diagnostics.R - statistical validation and significance
################################################################################

#' Structural zero-coupling check
#'
#' `nil_causality_check()` is a structural zero-coupling check. It is a
#' diagnostic, not a substitute for checking model assumptions or source
#' provenance.
#'
#' @param model A `cgns_model` object.
#' @param obs Optional observed trajectory; `NULL` simulates one.
#' @param direction Character naming the causal direction to check.
#' @param floor Positive numeric level below which the metric counts as nil.
#' @param T Positive total simulated time used when `obs` is `NULL`.
#' @param dt Positive integration step used when `obs` is `NULL`.
#' @param seed Non-negative whole number seeding the simulation.
#' @returns A list carrying the realised metric and the verdict of the check.
#'
#' @seealso [nil_surrogate_test()], [cross_validate()]
#'
#' @examples
#' nil_causality_check(model_dyad(), direction = list(cause = 1, effect = 1),
#'                     T = 5, dt = 0.01)
#'
#' @export
nil_causality_check <- function(model, obs = NULL, direction,
                                floor = 1e-4, T = 50, dt = 1e-3, seed = 1) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract",
              "nil_causality_check() needs a cgns_model for direction-specific marginals.")
  if (!is.list(direction) || is.null(direction$cause) || is.null(direction$effect))
    aci_abort("aci_error_direction",
              "direction must be a list with `cause` and `effect` indices.")
  raw_cause <- direction$cause; raw_effect <- direction$effect
  valid_idx <- function(z, upper)
    is.numeric(z) && length(z) && all(is.finite(z)) &&
    all(z == base::floor(z)) && !anyDuplicated(z) && all(z >= 1L & z <= upper)
  if (!valid_idx(raw_cause, model$l) || !valid_idx(raw_effect, model$k))
    aci_abort("aci_error_direction",
              "direction cause/effect indices must be unique finite integers inside the model dimensions.")
  cause <- sort(as.integer(raw_cause)); effect <- sort(as.integer(raw_effect))
  if (!is.numeric(floor) || length(floor) != 1L || !is.finite(floor) || floor < 0)
    aci_abort("aci_error_dims", "floor must be one finite non-negative value.")
  if (is.null(obs)) obs <- simulate(model, seed = seed, T = T, dt = dt)$obs
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "The observation width does not match the model.")
  if (length(obs$t) < 2L)
    aci_abort("aci_error_obs_contract",
              "nil_causality_check() needs at least two observation times.")

  # For a CGNS the exact direct sensitivity is the selected Lx block. Evaluate
  # it on the same realised trajectory and absolute clock as the empirical
  # check; random probes on [0, 10] can miss seasonal/path-dependent coupling
  # and used to mutate the caller's RNG even when `obs` was supplied.
  direct <- max(vapply(seq_along(obs$t), function(j) {
    L <- as.matrix(model$Lx(obs$t[j], obs$x[j, ]))
    max(abs(L[effect, cause, drop = FALSE]))
  }, numeric(1)))

  # Use only the selected effect innovations while retaining every other
  # observed component as prescribed state information.  Then take the same
  # selected hidden marginal in both Gaussian paths.  The previous
  # implementation scored the full model, so `direction` affected only the
  # structural probe and could report a false failure due to another channel.
  other_effects <- setdiff(seq_len(model$k), effect)
  nt <- if (length(other_effects))
    nontarget(other_effects, strategy = "inflate") else NULL
  filt <- suppressWarnings(da_filter(model, obs, nontarget = nt))
  smoo <- da_smooth(model, obs, filter = filt, nontarget = nt)
  marginal <- function(path) {
    cv <- path$cov[cause, cause, , drop = FALSE]
    dim(cv) <- c(length(cause), length(cause), length(path$t))
    new_da_path(path$t, path$mean[, cause, drop = FALSE], cv, path$kind,
                meta = path$meta)
  }
  score <- gaussian_kl_path(marginal(smoo), marginal(filt))
  burn <- min(ceiling(0.1 * nrow(score)), nrow(score) - 1L)
  peak <- max(score$total[-seq_len(burn)])
  a <- structure(list(t = score$t, aci = score$total,
                      signal = score$signal, dispersion = score$dispersion,
                      paths = list(filter = filt, smoother = smoo),
                      meta = list(engine = "cgns", direction = direction,
                                  source_status = "validation diagnostic")),
                 class = "aci_result")
  structure(list(direct_sensitivity = direct, empirical_peak = peak,
                 floor = floor,
                 structural_pass = direct < 1e-10,
                 empirical_pass = peak < floor,
                 direction = list(cause = cause, effect = effect),
                 series = a,
                 meta = list(source_status = paste(
                   "Direction-specific validation diagnostic; this wrapper is",
                   "not an estimand defined in the supplied ACI papers."))),
            class = "nilcheck")
}


#' Print a nil-causality check
#'
#' @param x A `nilcheck` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.nilcheck <- function(x, ...) {
  cat(sprintf(
    "<nilcheck> cause %s -> effect %s | direct sensitivity %.3g | empirical peak %.3g (floor %.3g)\n",
    paste(x$direction$cause, collapse = ","),
    paste(x$direction$effect, collapse = ","),
    x$direct_sensitivity, x$empirical_peak, x$floor))
  cat(sprintf(
    "  structural_pass: %s | empirical_pass: %s (TRUE flags the direction as nil)\n",
    x$structural_pass, x$empirical_pass))
  invisible(x)
}


#' Cross-engine validation of ensemble against closed-form moments
#'
#' `cross_validate()` compares ensemble moments with a compatible closed-form
#' CGNS result. This is a diagnostic, not a substitute for checking model
#' assumptions or source provenance.
#'
#' @param model A `cgns_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param m_grid Numeric vector of ensemble sizes to compare.
#' @param n_rep Number of repetitions at each ensemble size.
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param inflation Positive multiplicative variance inflation factor.
#' @param seed Non-negative whole number seeding the ensembles.
#' @param init Optional list with the initial hidden `mean` and `cov`.
#' @param ... Must be empty; unused arguments are an error.
#' @returns A data frame with one row per ensemble size and repetition, carrying
#'   the discrepancy against the closed-form moments.
#'
#' @seealso [nil_causality_check()], [as_gaussian()]
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' cross_validate(m, ob, m_grid = c(10, 20), n_rep = 1)
#' }
#'
#' @export
cross_validate <- function(model, obs, m_grid = c(10, 20, 50, 100), n_rep = 3,
                           localization = NULL, inflation = 1, seed = 1,
                           init = NULL, ...) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "cross_validate needs a cgns_model oracle.")
  dots <- list(...)
  if ("ic_sampler" %in% names(dots))
    aci_abort("aci_error_model_contract", paste(
      "cross_validate() controls the ensemble initial sampler so it matches",
      "the oracle prior; pass that Gaussian prior as `init`, not `ic_sampler`."))
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "The observation width does not match the model.")
  if (!is.numeric(m_grid) || !length(m_grid) || any(!is.finite(m_grid)) ||
      any(m_grid != floor(m_grid)) || any(m_grid <= model$l) ||
      anyDuplicated(m_grid))
    aci_abort("aci_error_ensemble_rank",
              "m_grid must contain unique integer ensemble sizes greater than l.")
  m_grid <- as.integer(m_grid)
  if (!is.numeric(n_rep) || length(n_rep) != 1L || !is.finite(n_rep) ||
      n_rep < 1L || n_rep != floor(n_rep))
    aci_abort("aci_error_dims", "n_rep must be a positive integer.")
  n_rep <- as.integer(n_rep)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) ||
      seed + 1000 * max(m_grid) + n_rep > .Machine$integer.max)
    aci_abort("aci_error_dims", "seed must keep every derived replicate seed inside R's integer range.")
  filt <- suppressWarnings(da_filter(model, obs, init = init))
  smoo <- da_smooth(model, obs, filter = filt)
  aci_o <- gaussian_kl_path(smoo, filt)$total
  common_init <- list(mean = as.numeric(filt$mean[1, ]),
                      cov = as.matrix(filt$cov[, , 1]))
  out <- do.call(rbind, lapply(m_grid, function(mm) do.call(rbind, lapply(
    seq_len(n_rep), function(r) {
      fr <- do.call(enkbf, c(list(model = model, obs = obs, m = mm,
                    seed = seed + 1000 * mm + r,
                    ic_sampler = common_init,
                    localization = localization, inflation = inflation), dots))
      # ensemble smoother needs obs_x for coefficient evaluation
      fr$path$meta$obs_x <- obs$x
      sm <- as_gaussian(enkbs(model, fr$path, fr$noise, localization = localization))
      fg <- as_gaussian(fr$path)
      aci_e <- gaussian_kl_path(sm, fg)$total
      data.frame(m = mm, rep = r,
                 smoother_rmse = sqrt(mean((sm$mean - smoo$mean)^2)),
                 aci_rmse = sqrt(mean((aci_e - aci_o)^2)))
    }))))
  attr(out, "oracle_peak_aci") <- max(aci_o)
  attr(out, "initial_prior") <- common_init
  attr(out, "source_status") <- paste(
    "Paper-inspired ensemble/oracle diagnostic; not a prescribed ACI",
    "calibration procedure. Both engines use the same initial Gaussian prior.")
  class(out) <- c("cv_result", class(out))
  out
}


#' Surrogate-null significance for the hidden-to-observed verdict
#'
#' Diagnostic layer adopted from the Opus-4.8 demo, with the statistic
#' replaced for the specified-model setting: the demo scores surrogates under
#' the null model (both sides ~ 0, a weak null tied to its fit-from-data
#' pipeline), while the raw ACI metric is power-inverted here (model error on
#' decoupled data inflates it). This test uses the predictive
#' likelihood-ratio LR = loglik(x | coupled candidate) - loglik(x | decoupled
#' candidate), whose parametric-bootstrap null comes from B simulations of
#' the decoupled model (Lx = 0, all other dynamics and noise kept). The
#' one-sided Monte Carlo p-value is one plus the number of null likelihood
#' ratios at least as large as the observed ratio, divided by `B + 1`.
#'
#' When the candidate model was learned from the same series (the
#' [aci_fit()] route), the test is not out of sample: the coupled model was
#' tuned to the observations, and neither model is re-estimated on the
#' surrogate replicates, so the null ratios carry no matching fitting
#' advantage and small p-values are optimistic there.
#'
#' @param model Coupled candidate model.
#' @param obs Observed trajectory.
#' @param B Number of parametric-bootstrap replicates.
#' @param burn_frac Fraction discarded before likelihood scoring.
#' @param seed Reproducibility seed.
#' @param init Optional Gaussian initialization.
#' @param ... Additional assimilation arguments.
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' nil_surrogate_test(m, ob, B = 4)
#' }
#'
#' @export
nil_surrogate_test <- function(model, obs, B = 49, burn_frac = 0.1,
                               seed = 1, init = NULL, ...) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "nil_surrogate_test needs a cgns_model in v0.")
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1L || B != as.integer(B))
    aci_abort("aci_error_dims", "B must be a positive integer.")
  B <- as.integer(B)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) ||
      seed + 7919 * B + 262 > .Machine$integer.max)
    aci_abort("aci_error_dims",
              "seed must keep every derived bootstrap seed inside R's integer range.")
  if (!is.numeric(burn_frac) || length(burn_frac) != 1L || !is.finite(burn_frac) ||
      burn_frac < 0 || burn_frac >= 1)
    aci_abort("aci_error_dims", "burn_frac must lie in [0, 1).")
  obs <- as_obs(obs)
  k <- model$k; l <- model$l
  null_model <- cgns_model(
    Lx = function(t, x) matrix(0, k, l),
    fx = model$fx, Ly = model$Ly, fy = model$fy,
    Sx1 = model$Sx1, Sx2 = model$Sx2, Sy1 = model$Sy1, Sy2 = model$Sy2,
    k = k, l = l, name = paste0(model$name, "|decoupled"), meta = model$meta)
  drop_burn <- function(ob) {
    j0 <- floor(burn_frac * (length(ob$t) - 1L))
    keep <- seq.int(j0 + 1L, length(ob$t))
    if (length(keep) < 2L)
      aci_abort("aci_error_obs_contract",
                "burn_frac leaves fewer than two observation times.")
    # Keep the original clock. Resetting the retained times to zero changes a
    # seasonal/time-inhomogeneous model and therefore changes the null being
    # tested.
    observed_trajectory(ob$t[keep], ob$x[keep, , drop = FALSE],
                        names = colnames(ob$x))
  }
  lr_of <- function(ob) {
    ob <- drop_burn(ob)
    llc <- suppressWarnings(da_filter(model, ob, init = init, ...))$meta$loglik
    ll0 <- suppressWarnings(da_filter(null_model, ob, init = init, ...))$meta$loglik
    llc - ll0
  }
  lr_obs <- lr_of(obs)
  Tspan <- obs$t[length(obs$t)] - obs$t[1]
  time_origin <- obs$t[1]
  y0 <- model$meta$ic_default$y0 %||% rep(0, l)
  shifted_for_simulation <- function(m) {
    if (abs(time_origin) < .Machine$double.eps) return(m)
    shift <- function(fun) function(t, x) fun(t + time_origin, x)
    cgns_model(Lx = shift(m$Lx), fx = shift(m$fx),
               Ly = shift(m$Ly), fy = shift(m$fy),
               Sx1 = shift(m$Sx1), Sx2 = shift(m$Sx2),
               Sy1 = shift(m$Sy1), Sy2 = shift(m$Sy2),
               k = m$k, l = m$l, name = paste0(m$name, "|time-shifted"),
               meta = m$meta)
  }
  absolute_obs <- function(ob) observed_trajectory(
    ob$t + time_origin, ob$x, names = colnames(obs$x))
  sim_null <- shifted_for_simulation(null_model)
  lr_null <- vapply(seq_len(B), function(b) {
    for (att in 0:2) {                # learned/decoupled nulls can diverge
      sb <- tryCatch(
        simulate(sim_null, seed = seed + 7919 * b + 131 * att, T = Tspan,
                 dt = obs$dt, ic = list(x0 = obs$x[1, ], y0 = y0)),
        error = function(e) NULL)
      if (!is.null(sb)) return(lr_of(absolute_obs(sb$obs)))
    }
    NA_real_
  }, numeric(1))
  lr_null <- lr_null[is.finite(lr_null)]
  null_type <- "decoupled"
  if (length(lr_null) < ceiling(B / 2)) {
    # The decoupled candidate is unstable on its own (learned polynomial
    # drifts often are once the restraining coupling is removed). Fall back
    # to a *linearized* null in the same family: the observed drift is
    # replaced by its best linear fit to the data (per channel, slope clamped
    # negative), noise and hidden dynamics kept. Flagged in the result.
    aci_warn("aci_warn_surrogate_linearized", paste(
      "Decoupled-null simulations diverged; using the linearized",
      "stationary null (see ?nil_surrogate_test)."))
    ab <- vapply(seq_len(k), function(i) {
      xi <- obs$x[-nrow(obs$x), i]; dxi <- diff(obs$x[, i]) / obs$dt
      cf <- stats::coef(stats::lm(dxi ~ xi))
      c(cf[1], -abs(cf[2]))
    }, numeric(2))
    lin_null <- cgns_model(
      Lx = function(t, x) matrix(0, k, l),
      fx = function(t, x) ab[1, ] + ab[2, ] * x,
      Ly = model$Ly, fy = model$fy,
      Sx1 = model$Sx1, Sx2 = model$Sx2, Sy1 = model$Sy1, Sy2 = model$Sy2,
      k = k, l = l, name = paste0(model$name, "|linearized-null"),
      meta = model$meta)
    lr_of2 <- function(ob) {
      ob2 <- drop_burn(ob)
      suppressWarnings(da_filter(model, ob2, init = init, ...))$meta$loglik -
        suppressWarnings(da_filter(lin_null, ob2, init = init, ...))$meta$loglik
    }
    lr_obs <- lr_of2(obs)
    sim_lin_null <- shifted_for_simulation(lin_null)
    lr_null <- vapply(seq_len(B), function(b) {
      sb <- tryCatch(
        simulate(sim_lin_null, seed = seed + 7919 * b, T = Tspan, dt = obs$dt,
                 ic = list(x0 = obs$x[1, ], y0 = y0)),
        error = function(e) NULL)
      if (is.null(sb)) NA_real_ else lr_of2(absolute_obs(sb$obs))
    }, numeric(1))
    lr_null <- lr_null[is.finite(lr_null)]
    null_type <- "linearized"
    if (length(lr_null) < ceiling(B / 2))
      aci_abort("aci_error_surrogate_unstable",
                "Even the linearized null diverged; test not applicable.")
  }
  complete <- length(lr_null) == B
  if (!complete)
    aci_warn("aci_warn_surrogate_incomplete", sprintf(paste(
      "Only %d of %d requested null simulations were usable; the reported",
      "Monte Carlo p-value is conditional on successful simulations and must",
      "be treated as exploratory."), length(lr_null), B))
  structure(list(observed_lr = lr_obs, null_lr = lr_null, B = length(lr_null),
                 requested_B = B,
                 null_type = null_type,
                 p_value = (1 + sum(lr_null >= lr_obs)) / (length(lr_null) + 1),
                 calibration = list(
                   source_status = paste(
                     "Experimental parametric-bootstrap validation extension;",
                     "not defined in the supplied ACI papers or MATLAB code."),
                   reference_null = "Lx = 0 with all other coefficients retained",
                   effective_null = null_type,
                   exact_reference_null = identical(null_type, "decoupled"),
                   all_requested_replicates_used = complete,
                   failed_replicates = B - length(lr_null),
                   time_origin_preserved = TRUE,
                   time_origin = time_origin)),
            class = "aci_surrogate_test")
}


#' Print a surrogate-null test result
#'
#' @param x An `aci_surrogate_test` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.aci_surrogate_test <- function(x, ...) {
  cat(sprintf("<aci_surrogate_test> LR = %.4g | %s null (B = %d/%d): median %.4g, max %.4g | p = %.3f%s\n",
              x$observed_lr, x$null_type, x$B, x$requested_B %||% x$B,
              stats::median(x$null_lr), max(x$null_lr), x$p_value,
              if (isFALSE(x$calibration$all_requested_replicates_used))
                " (incomplete calibration)" else ""))
  invisible(x)
}


#' Check a precomputed range against its model and observations (internal)
#'
#' A precomputed CIR is safe to reuse only when it came from the same exact
#' model, observations, conditioning choice, and Gaussian prior.  A bare data
#' frame has the right shape but cannot establish any of those facts.
#'
#' @param x A `cir_result` object.
#' @param model A `cgns_model` object.
#' @param obs An `obs_traj` object.
#' @param actual_init The prior the accompanying filter was run from.
#' @returns Invisibly `TRUE`; called for its error conditions.
#' @noRd
.validate_precomputed_cir <- function(x, model, obs, actual_init) {
  need <- c("t", "tau", "strength", "direction")
  if (!is.data.frame(x) || !nrow(x) || !all(need %in% names(x)))
    aci_abort("aci_error_dims",
              "Precomputed cir data must be a non-empty cir_table() result with t, tau, strength, and direction columns.")
  if (!is.numeric(x$t) || !is.numeric(x$tau) || !is.numeric(x$strength) ||
      any(!is.finite(x$t)) || any(!is.finite(x$strength)) ||
      any(x$strength < 0) || any(!is.finite(x$tau[!is.na(x$tau)])) ||
      !is.character(x$direction) || anyNA(x$direction) ||
      any(!x$direction %in% c("forward", "backward")) ||
      !any(x$direction == "forward"))
    aci_abort("aci_error_dims", paste(
      "Precomputed cir data must contain finite times and non-negative strengths,",
      "finite-or-NA ranges, and at least one forward row."))
  p <- attr(x, "aci_provenance", exact = TRUE)
  if (!is.list(p) || is.null(p$source_model) || is.null(p$source_obs_t) ||
      is.null(p$source_obs_x) || is.null(p$actual_init))
    aci_abort("aci_error_dims", paste(
      "Precomputed cir data lacks intact provenance; pass the unmodified result",
      "returned by cir_table()."))
  if (!identical(p$source_model, model))
    aci_abort("aci_error_model_contract",
              "Precomputed cir data was computed with a different model object.")
  if (length(p$source_obs_t) != length(obs$t) ||
      max(abs(p$source_obs_t - obs$t)) > 1e-10 * max(1, max(abs(obs$t))) ||
      !identical(dim(p$source_obs_x), dim(obs$x)) ||
      any(abs(p$source_obs_x - obs$x) > 1e-12 * pmax(1, abs(obs$x))))
    aci_abort("aci_error_dims",
              "Precomputed cir data was computed from different observations.")
  if (!is.null(p$nontarget))
    aci_abort("aci_error_nontarget", paste(
      "aci_check() is unconditional, but the precomputed cir data used a",
      "non-target conditioning specification."))
  if (!.same_gaussian_init(p$actual_init, actual_init, model$l))
    aci_abort("aci_error_dims",
              "Precomputed cir data was computed with a different Gaussian prior.")
  invisible(TRUE)
}


#' One-call diagnostic for a fitted or specified ACI model
#'
#' The counterpart to `gam.check()`: draws the reconstructed latent with its
#' uncertainty band and the influence metric, then prints the verdict, the
#' forward CIR, and a reconstruction-consistency check. The latter tests the
#' full matrix ordering `filter covariance - smoother covariance >= 0` at every
#' step, not only the first latent marginal. The compact plot still shows latent
#' component 1; the returned consistency diagnostics cover every component.
#' `test = TRUE` runs the surrogate-null likelihood-ratio bootstrap; `FALSE`
#' reports the effect-size screen only.
#'
#' @param model A fitted or specified CGNS model.
#' @param obs An observed trajectory coercible with `as_obs()`.
#' @param init Optional Gaussian initialization list with `mean` and `cov`.
#' @param test Whether to run the surrogate-null diagnostic.
#' @param B Number of bootstrap replicates.
#' @param floor Numerical floor used by the effect-size screen.
#' @param plot Whether to draw diagnostic plots.
#' @param scale Multiplicative scale for the displayed hidden state.
#' @param unit Label for the hidden-state scale.
#' @param seed Reproducibility seed.
#' @param cir Whether to calculate CIR summaries.
#' @param ... Additional arguments forwarded to lower-level routines.
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' aci_check(m, ob, test = FALSE, plot = FALSE)
#' }
#'
#' @export
aci_check <- function(model, obs, init = NULL, test = TRUE, B = 99L,
                      floor = 1e-3, plot = TRUE, scale = 1, unit = "",
                      seed = NULL, cir = TRUE, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to aci_check().")
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "aci_check() needs a cgns_model.")
  if (!is.logical(test) || length(test) != 1L || is.na(test) ||
      !is.logical(plot) || length(plot) != 1L || is.na(plot))
    aci_abort("aci_error_dims", "test and plot must each be TRUE or FALSE.")
  if (isTRUE(test) && (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1L || B != floor(B)))
    aci_abort("aci_error_dims", "B must be a positive integer when test = TRUE.")
  if (!is.numeric(floor) || length(floor) != 1L || !is.finite(floor) || floor < 0)
    aci_abort("aci_error_dims", "floor must be one finite non-negative number.")
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0)
    aci_abort("aci_error_dims", "scale must be one finite positive number.")
  if (!is.character(unit) || length(unit) != 1L || is.na(unit))
    aci_abort("aci_error_dims", "unit must be one character string.")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)))
    aci_abort("aci_error_dims", "seed must be NULL or one non-negative integer.")
  if (!(is.logical(cir) && length(cir) == 1L && !is.na(cir)) &&
      !is.data.frame(cir))
    aci_abort("aci_error_dims", "cir must be TRUE, FALSE, or a precomputed cir_table() data frame.")
  obs <- as_obs(obs)
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the CGNS model.")
  fl  <- suppressWarnings(da_filter(model, obs, init = init))
  sm  <- da_smooth(model, obs, filter = fl)
  kl  <- gaussian_kl_path(sm, fl)
  burn <- seq_len(ceiling(0.10 * nrow(kl)))
  post <- kl[-burn, ]
  sur <- if (isTRUE(test))
    suppressWarnings(nil_surrogate_test(model, obs, B = B, init = init,
                                        seed = seed %||% 1L)) else NULL
  ## Building a forward lag table is the expensive part of this function. When
  ## the caller already has one (from cir_table()), pass it in as `cir`, in a
  ## multi-site workflow that is the difference between one table per site and
  ## several.
  ct <- if (is.data.frame(cir)) cir else if (isTRUE(cir))
    cir_table(model, obs, init = init, direction = "forward") else NULL
  if (!is.null(ct)) .validate_precomputed_cir(ct, model, obs, fl$meta$init)

  # Conditioning on future observations cannot increase Gaussian covariance:
  # R_filter(t) - R_smoother(t) must be positive semidefinite. Test the entire
  # l by l ordering with a scale-aware numerical tolerance at every time.
  consistency_min_eigen <- vapply(seq_along(obs$t), function(j)
    min(eigen(sym(fl$cov[, , j] - sm$cov[, , j]), symmetric = TRUE,
              only.values = TRUE)$values), numeric(1))
  consistency_tolerance <- vapply(seq_along(obs$t), function(j)
    1e-10 * max(1, abs(fl$cov[, , j]), abs(sm$cov[, , j])), numeric(1))
  consistency_by_step <- consistency_min_eigen >= -consistency_tolerance
  ok <- mean(consistency_by_step)
  verdict <- if (!is.null(sur)) isTRUE(sur$p_value <= 0.05) && mean(post$total) > floor
             else mean(post$total) > floor

  if (isTRUE(plot)) {
    op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4.2, 2, 1)); on.exit(graphics::par(op))
    mu <- sm$mean[, 1] * scale; sd <- sqrt(pmax(sm$cov[1, 1, ], 0)) * scale
    plot(obs$t, mu, type = "n", xlab = "time", ylab = paste0("latent ", unit),
         ylim = range(c(mu - 2 * sd, mu + 2 * sd)), main = "Reconstructed latent cause")
    graphics::polygon(c(obs$t, rev(obs$t)), c(mu - 2 * sd, rev(mu + 2 * sd)),
                      border = NA, col = grDevices::adjustcolor("steelblue", 0.25))
    graphics::lines(obs$t, mu, col = "steelblue", lwd = 2)
    graphics::lines(obs$t, fl$mean[, 1] * scale, col = "grey40", lty = 2)
    graphics::legend("topright", c("smoother +/- 2sd", "filter"),
                     col = c("steelblue", "grey40"), lty = c(1, 2), bty = "n", cex = .8)
    plot(kl$t, kl$total, type = "l", lwd = 2, xlab = "time", ylab = "ACI (nats)",
         main = "Influence over time: KL(smoother || filter)")
    graphics::lines(kl$t, kl$signal, col = 4, lty = 2)
    graphics::lines(kl$t, kl$dispersion, col = 2, lty = 3)
    if (!is.null(ct)) {
      ft <- ct$t[ct$direction == "forward"]
      graphics::rug(ft[unique(round(seq(1, length(ft), length.out = 12)))],
                    col = "firebrick")
    }
  }

  cat("Assimilative causal inference: model check\n")
  cat(sprintf("  model  : %s (k = %d observed, l = %d hidden, %d steps)\n",
              model$name, model$k, model$l, length(obs$t)))
  cat(sprintf("  ACI    : peak %.4g, post-burn mean %.4g (signal %.3g / dispersion %.3g)\n",
              max(post$total), mean(post$total), mean(post$signal), mean(post$dispersion)))
  if (!is.null(sur))
    cat(sprintf("  verdict: %s   (surrogate LR p = %.3f, B = %d)\n",
                if (verdict) "cause drives the effects" else "not distinguishable from the null",
                sur$p_value, sur$B))
  else
    cat(sprintf("  verdict: %s   (effect-size screen, floor %.3g)\n",
                if (verdict) "cause drives the effects" else "below floor", floor))
  if (!is.null(ct)) {
    ftau <- ct$tau[ct$direction == "forward" & is.finite(ct$tau)]
    if (length(ftau))
      cat(sprintf("  forward CIR: median %.4g, range %.4g-%.4g time units\n",
                  stats::median(ftau), min(ftau), max(ftau)))
    else
      cat("  forward CIR: all ranges masked by the minimum-strength rule\n")
  }
  cat(sprintf("  reconstruction: full covariance PSD order at %.0f%% of steps%s\n",
              100 * ok, if (isTRUE(all.equal(ok, 1))) "  (consistent)" else "  (CHECK)"))
  invisible(list(verdict = verdict, surrogate = sur, cir = ct, kl = kl,
                 filter = fl, smoother = sm, consistency = ok,
                 consistency_by_step = consistency_by_step,
                 consistency_min_eigen = consistency_min_eigen,
                 consistency_tolerance = consistency_tolerance,
                 consistency_definition =
                   "filter covariance - smoother covariance is positive semidefinite"))
}


#' Observing-system simulation experiment for a CGNS model
#'
#' Simulate from the model, WITHHOLD the hidden path, re-assimilate the
#' simulated observations, and score the reconstruction against the output truth.
#' The perfect-model twin is the ACI papers' own validation mode
#' (andreou2026aci Results: the case studies use "a perfect-model setup ...
#' meaning the system generates the synthetic observational data and the
#' forecast model in ACI at the same"; jiang2026enkbs sec 3.2 benchmarks the
#' ensemble smoother against the closed form by RMSE), and the standard
#' observing-system simulation experiment of data assimilation
#' (Arnold & Dey 1986).
#'
#' It complements the two checks already available. `nil_causality_check()` asks
#' whether the metric vanishes when the coupling is removed; the smoother-vs-
#' filter variance comparison in `aci_check()` is a theoretical consistency test
#' (the smoother conditions on strictly more information, so it cannot be less
#' certain). This asks the empirical question neither answers: given this model
#' and this observation design, how well is the latent actually recovered?
#'
#' Reports the filter and smoother RMSE, the smoother's improvement over the
#' filter, and the share of latent-component/time cells where the truth lies
#' inside its marginal posterior 2-sd band. All scalar scores aggregate every
#' hidden component; component-level summaries are returned in `per_component`.
#' That last score is a calibration check from general uncertainty
#' quantification, not something the ACI papers prescribe: read it as "is the
#' posterior spread honest", and treat a value far from ~0.95 as a warning about
#' the noise amplitudes rather than about the causal verdict.
#'
#' @param model CGNS model.
#' @param init Gaussian assimilation initialization.
#' @param T Simulation duration.
#' @param dt Time step.
#' @param x0 Optional observed initial state.
#' @param y0 Optional hidden initial state.
#' @param nrep Number of twin replicates.
#' @param seed Reproducibility seed.
#' @param stepper Filter covariance stepper.
#' @param nsub Number of filter substeps.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications* **17**, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#'
#' Arnold, C. P., Jr. and Dey, C. H. (1986). Observing-systems simulation
#' experiments: past, present, and future. *Bulletin of the American
#' Meteorological Society* **67**(6), 687-695.
#' \doi{10.1175/1520-0477(1986)067<0687:OSSEPP>2.0.CO;2}
#' @examples
#' \donttest{
#' osse_twin(model_dyad(), init = list(mean = 2, cov = matrix(0.1, 1, 1)),
#'           T = 2, dt = 0.01, nrep = 2)
#' }
#'
#' @export
osse_twin <- function(model, init, T, dt, x0 = NULL, y0 = NULL, nrep = 12L,
                      seed = 1L, stepper = c("explicit", "implicit"), nsub = 1L) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "osse_twin() needs a cgns_model.")
  stepper <- match.arg(stepper)
  if (!is.numeric(T) || length(T) != 1L || !is.finite(T) || T <= 0 ||
      !is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0)
    aci_abort("aci_error_dims", "T and dt must be finite positive scalars.")
  N <- round(T / dt)
  if (N < 1L || abs(N * dt - T) > 1e-8 * max(T, dt))
    aci_abort("aci_error_dims", "T must be an integer multiple of dt.")
  if (!is.numeric(nrep) || length(nrep) != 1L || !is.finite(nrep) ||
      nrep < 1L || nrep != floor(nrep))
    aci_abort("aci_error_dims", "nrep must be a positive integer.")
  nrep <- as.integer(nrep)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed + nrep > .Machine$integer.max)
    aci_abort("aci_error_dims", paste(
      "seed must be a non-negative integer whose per-twin seed sequence",
      "does not exceed the R integer range."))
  if (!is.numeric(nsub) || length(nsub) != 1L || !is.finite(nsub) ||
      nsub < 1L || nsub != floor(nsub))
    aci_abort("aci_error_dims", "nsub must be a positive integer.")
  nsub <- as.integer(nsub)
  if (!is.null(x0) && (!is.numeric(x0) || length(x0) != model$k ||
      any(!is.finite(x0))))
    aci_abort("aci_error_dims", "x0 must contain one finite value per observed state.")
  if (!is.null(y0) && (!is.numeric(y0) || length(y0) != model$l ||
      any(!is.finite(y0))))
    aci_abort("aci_error_dims", "y0 must contain one finite value per hidden state.")
  if (!is.null(init)) {
    if (!is.list(init))
      aci_abort("aci_error_dims", "init must be NULL or a Gaussian-prior list.")
    if (!is.null(init$mean) && (!is.numeric(init$mean) ||
        length(init$mean) != model$l || any(!is.finite(init$mean))))
      aci_abort("aci_error_dims", "init$mean must contain one finite value per hidden state.")
    if (!is.null(init$cov)) {
      R0 <- as.matrix(init$cov)
      if (!is.numeric(R0) || !identical(dim(R0), c(model$l, model$l)) ||
          any(!is.finite(R0)))
        aci_abort("aci_error_dims", "init$cov must be a finite l by l matrix.")
      .strict_chol(R0, "osse_twin init$cov")
    }
  }
  m <- model
  if (!is.null(x0) || !is.null(y0)) {
    base_ic <- m$meta$ic_default %||% list()
    if (!is.list(base_ic))
      aci_abort("aci_error_model_contract", "model$meta$ic_default must be a list.")
    m$meta$ic_default <- list(x0 = x0 %||% base_ic$x0,
                              y0 = y0 %||% base_ic$y0)
  }
  twins <- lapply(seq_len(nrep), function(i) {
    sim <- simulate(m, seed = seed + i, T = T, dt = dt, burn_in = 0)
    fl  <- suppressWarnings(da_filter(m, sim$obs, init = init, stepper = stepper,
                                      nsub = nsub))
    sm  <- da_smooth(m, sim$obs, filter = fl)
    ef <- fl$mean - sim$hidden
    es <- sm$mean - sim$hidden
    # vapply() returns a bare vector when l == 1, whose transpose is 1 x n
    # rather than the n x l `es` is; the reshape fixes it for every l.
    marginal_var <- t(matrix(vapply(seq_along(sm$t), function(j)
      diag(as.matrix(sm$cov[, , j])), numeric(m$l)), nrow = m$l))
    covered <- abs(es) <= 2 * sqrt(pmax(marginal_var, 0))
    list(
      scalar = c(filter = sqrt(mean(ef^2)),
                 smoother = sqrt(mean(es^2)),
                 coverage = mean(covered)),
      filter_by_component = sqrt(colMeans(ef^2)),
      smoother_by_component = sqrt(colMeans(es^2)),
      coverage_by_component = colMeans(covered))
  })
  res <- do.call(rbind, lapply(twins, `[[`, "scalar"))
  fcomp <- do.call(rbind, lapply(twins, `[[`, "filter_by_component"))
  scomp <- do.call(rbind, lapply(twins, `[[`, "smoother_by_component"))
  ccomp <- do.call(rbind, lapply(twins, `[[`, "coverage_by_component"))
  component_names <- paste0("y", seq_len(m$l))
  colnames(fcomp) <- colnames(scomp) <- colnames(ccomp) <- component_names
  out <- list(nrep = nrep, per_twin = res, filter = mean(res[, "filter"]),
              smoother = mean(res[, "smoother"]), coverage = mean(res[, "coverage"]))
  out$gain <- if (out$filter > 0) 1 - out$smoother / out$filter else NA_real_
  out$l <- m$l
  out$per_component <- data.frame(
    component = component_names,
    filter = colMeans(fcomp),
    smoother = colMeans(scomp),
    coverage = colMeans(ccomp),
    row.names = NULL)
  out$per_twin_by_component <- list(
    filter = fcomp, smoother = scomp, coverage = ccomp)
  out$meta <- list(
    aggregation = paste(
      "RMSE and marginal 2-sd coverage aggregate all hidden components and",
      "all retained time points within each twin, then average over twins."),
    stepper = stepper, nsub = nsub,
    seeds = seed + seq_len(nrep),
    source_status = paste(
      "Perfect-model validation diagnostic; scalar coverage is a general",
      "uncertainty-calibration extension, not an ACI-paper estimand."))
  structure(out, class = "osse_twin")
}


#' Print a perfect-model OSSE result
#'
#' @param x An `osse_twin` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.osse_twin <- function(x, ...) {
  cat(sprintf("<osse_twin> %d twins, %d latent component%s | RMSE: filter %.4g -> smoother %.4g (%.0f%% better) | marginal 2sd coverage %.0f%%\n",
              x$nrep, x$l %||% 1L, if ((x$l %||% 1L) == 1L) "" else "s",
              x$filter, x$smoother, 100 * x$gain, 100 * x$coverage))
  invisible(x)
}
