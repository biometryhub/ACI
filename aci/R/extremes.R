################################################################################
## extremes.R - extreme-event workflow
## ########################################################################## ##
##
## Contents:
##   * extreme-event detection, onsets, composites, sensitive directions:
##       - detect_events, print.event_set, event_onset, event_influence, event_stats,
##         sensitive_directions
##
##   * moser2026extremes feature sets and mechanism clustering:
##       - .win_idx, .slope, features_pathways, features_topographic, feature_set_default,
##         classify_events, print.event_classes,
##
################################################################################


################################################################################
# extreme-event detection, onsets, composites, sensitive directions
################################################################################

#' Detect extreme events in an observed component (moser2026extremes s4.1-4.3
#' conventions)
#'
#' `threshold = "event_peak_quantile"` applies the quantile to amplitudes of
#' sign-consistent excursion peaks, matching the moser2026extremes dyad
#' convention. The raw trajectory quantile route (`threshold = "quantile"`) is
#' retained for the pathways-style workflow; the two definitions are not
#' interchangeable.
#'
#' @param obs Observed trajectory.
#' @param component Observed component index.
#' @param threshold Threshold mode: a raw-sample empirical quantile,
#'   a quantile of sign-excursion peak amplitudes, or an absolute magnitude.
#' @param q Quantile probability when `threshold = "quantile"` or
#'   `"event_peak_quantile"`.
#' @param value Absolute threshold when `threshold = "absolute"`.
#' @param two_sided Whether to detect both upper and lower extremes.
#' @param min_separation Minimum separation between retained event peaks.
#' @param burn_in Ignore observations before this time.
#' @param event_unit Whether events are sign-consistent excursions or complete threshold runs.
#' @returns An `event_set` data frame with one row per event: `peak_index`,
#'   `t_star` and `peak_value` (the peak), `sign`, `run_start` and `run_end`
#'   (the containing excursion or threshold run), and `duration`, the
#'   half-peak width (moser2026extremes Appendix A feature 3). `duration`
#'   measures the peak, not the run, so it is not `run_end - run_start`.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' detect_events(ob, q = 0.9)
#'
#' @export
detect_events <- function(obs, component = 1,
                          threshold = c("quantile", "event_peak_quantile", "absolute"), q = 0.9,
                          value = NULL, two_sided = FALSE,
                          min_separation = 0, burn_in = 0,
                          event_unit = c("sign_excursion", "threshold_run")) {
  obs <- as_obs(obs); threshold <- match.arg(threshold)
  event_unit <- match.arg(event_unit)
  if (!is.numeric(component) || length(component) != 1L ||
      !is.finite(component) || component != as.integer(component) ||
      component < 1L || component > obs$k)
    aci_abort("aci_error_dims", "component must index one observed variable.")
  if (!is.logical(two_sided) || length(two_sided) != 1L || is.na(two_sided))
    aci_abort("aci_error_dims", "two_sided must be TRUE or FALSE.")
  if (threshold %in% c("quantile", "event_peak_quantile") &&
      (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1))
    aci_abort("aci_error_dims", "q must be a finite value strictly between 0 and 1.")
  if (threshold == "quantile" && two_sided && q <= 0.5)
    aci_abort("aci_error_dims", "Two-sided sample quantiles require q > 0.5.")
  if (threshold == "absolute" &&
      (length(value) != 1L || !is.finite(value) || value <= 0))
    aci_abort("aci_error_dims", "value must be supplied as a finite positive magnitude for an absolute threshold.")
  if (threshold == "event_peak_quantile" && event_unit != "sign_excursion")
    aci_abort("aci_error_dims",
              "event_peak_quantile requires event_unit = 'sign_excursion'.")
  if (length(min_separation) != 1L || !is.finite(min_separation) || min_separation < 0 ||
      length(burn_in) != 1L || !is.finite(burn_in) || burn_in < 0)
    aci_abort("aci_error_dims", "min_separation and burn_in must be finite and non-negative.")
  x <- obs$x[, component]; t <- obs$t
  keep <- t >= burn_in
  if (!any(keep)) aci_abort("aci_error_dims", "burn_in excludes the complete trajectory.")
  peak_mode <- threshold == "event_peak_quantile"
  thr_hi <- if (peak_mode) Inf else if (threshold == "quantile")
    stats::quantile(x[keep], q) else abs(value)
  thr_lo <- if (peak_mode) -Inf else if (two_sided) {
    if (threshold == "quantile") stats::quantile(x[keep], 1 - q) else -abs(value)
  } else -Inf
  extreme <- if (peak_mode) keep else keep & (x > thr_hi | x < thr_lo)
  if (!peak_mode && !any(extreme)) {
    out <- structure(data.frame(), class = c("event_set", "data.frame"))
    attr(out, "threshold") <- c(hi = unname(thr_hi), lo = unname(thr_lo))
    attr(out, "component") <- component; attr(out, "event_unit") <- event_unit
    return(out)
  }
  membership <- if (event_unit == "threshold_run") extreme else {
    s <- sign(x)
    # Assign exact zeros to the nearest previous nonzero sign (and then the
    # following sign at the left boundary), so zero crossings delimit physical
    # excursions without fragmenting a flat crossing.
    for (i in seq_along(s)) if (s[i] == 0 && i > 1L) s[i] <- s[i - 1L]
    for (i in rev(seq_along(s))) if (s[i] == 0 && i < length(s)) s[i] <- s[i + 1L]
    s
  }
  runs <- rle(membership); ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  run_ids <- if (event_unit == "threshold_run") which(runs$values) else seq_along(runs$values)
  candidates <- lapply(run_ids, function(i) {
    idx <- starts[i]:ends[i]
    if (event_unit == "sign_excursion") idx <- idx[keep[idx]]
    if (!length(idx)) return(NULL)
    pk <- idx[which.max(abs(x[idx]))]
    if (!peak_mode && !extreme[pk]) return(NULL)
    if (peak_mode && !two_sided && x[pk] <= 0) return(NULL)
    half <- abs(x[pk]) / 2; sgn <- sign(x[pk])
    lo <- pk; while (lo > 1 && sign(x[lo - 1]) == sgn && abs(x[lo - 1]) >= half) lo <- lo - 1
    hi <- pk; while (hi < length(x) && sign(x[hi + 1]) == sgn && abs(x[hi + 1]) >= half) hi <- hi + 1
    data.frame(peak_index = pk, t_star = t[pk], peak_value = x[pk],
               sign = sgn, run_start = t[starts[i]], run_end = t[ends[i]],
               duration = t[hi] - t[lo])          # App. A feature 3 (half-peak width)
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (peak_mode && length(candidates)) {
    amp <- vapply(candidates, function(z)
      if (two_sided) abs(z$peak_value[1]) else z$peak_value[1], numeric(1))
    peak_threshold <- unname(stats::quantile(amp, q))
    # Inclusive selection retains the sole event in a one-excursion record and
    # treats tied peaks at the requested empirical quantile consistently.
    candidates <- candidates[amp >= peak_threshold]
    thr_hi <- peak_threshold
    thr_lo <- if (two_sided) -peak_threshold else -Inf
  }
  if (!length(candidates)) {
    out <- structure(data.frame(), class = c("event_set", "data.frame"))
    attr(out, "threshold") <- c(hi = unname(thr_hi), lo = unname(thr_lo))
    attr(out, "component") <- component; attr(out, "event_unit") <- event_unit
    return(out)
  }
  ev <- do.call(rbind, candidates)
  if (min_separation > 0 && nrow(ev) > 1) {
    ord <- order(-abs(ev$peak_value)); chosen <- integer()
    for (i in ord)
      if (!length(chosen) || all(abs(ev$t_star[i] - ev$t_star[chosen]) >= min_separation))
        chosen <- c(chosen, i)
    ev <- ev[chosen, , drop = FALSE]; ev <- ev[order(ev$t_star), , drop = FALSE]
  }
  rownames(ev) <- NULL
  attr(ev, "threshold") <- c(hi = unname(thr_hi), lo = unname(thr_lo))
  attr(ev, "component") <- component
  attr(ev, "event_unit") <- event_unit
  structure(ev, class = c("event_set", class(ev)))
}


#' Print an event set
#'
#' @param x An `event_set` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.event_set <- function(x, ...) {
  cat(sprintf("<event_set> %d events; thresholds hi=%.4g lo=%.4g\n",
              nrow(x), attr(x, "threshold")["hi"], attr(x, "threshold")["lo"]))
  if (nrow(x)) print.data.frame(utils::head(x, 5))
  invisible(x)
}


#' Event onset via the ACI (KL) threshold rule (moser2026extremes eqs. 3.1-3.2)
#'
#' Return the first time KL(smoother || filter) exceeds `kappa` between
#' `t_star - T_pre` and `t_star`.
#'
#' @param aci_result An object returned by `aci()`.
#' @param events Event set from `detect_events()`.
#' @param kappa ACI threshold.
#' @param T_pre Length of the pre-event search window.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' ev <- detect_events(ob, q = 0.9)
#' a <- aci(m, ob)
#' event_onset(a, ev)
#'
#' @export
event_onset <- function(aci_result, events, kappa = 0.8, T_pre = 1.5) {
  if (!inherits(aci_result, "aci_result"))
    aci_abort("aci_error_dims", "aci_result must be an aci_result object.")
  if (!is.data.frame(events) || !"t_star" %in% names(events) ||
      any(!is.finite(events$t_star)) ||
      any(events$t_star < min(aci_result$t) | events$t_star > max(aci_result$t)))
    aci_abort("aci_error_dims",
              "events must contain finite t_star values within the ACI time span.")
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) || kappa < 0)
    aci_abort("aci_error_dims", "kappa must be one finite non-negative threshold.")
  if (!is.numeric(T_pre) || length(T_pre) != 1L || !is.finite(T_pre) || T_pre <= 0)
    aci_abort("aci_error_dims", "T_pre must be one finite positive duration.")
  t <- aci_result$t; a <- aci_result$aci
  onset <- vapply(seq_len(nrow(events)), function(i) {
    win <- which(t >= events$t_star[i] - T_pre & t <= events$t_star[i])
    hit <- win[a[win] >= kappa]
    if (length(hit)) t[hit[1]] else NA_real_
  }, numeric(1))
  cbind(events, onset = onset, lead_time = events$t_star - onset)
}


#' Per-event forward influence and backward attribution windows
#'
#' The forward window applies the package's andreou2026cir
#' last-exceedance/layer-cake CIR to moser2026extremes event anchors. It agrees
#' with moser2026extremes eqs. 3.3-3.6 when the finite-lag KL decreases
#' monotonically, as moser2026extremes assumes, but moser2026extremes's
#' first-crossing definition differs on a non-monotone finite-grid row. The
#' paired backward window applies andreou2026cir backward CIR. Both are package
#' syntheses rather than a paired formula in moser2026extremes. The backward
#' lagged grid ends one step before the event reference time.
#'
#' @param model Model used for assimilation.
#' @param obs Observed trajectory.
#' @param events Event set from `detect_events()`.
#' @param forward_table Optional precomputed forward lag table.
#' @param nontarget Optional conditional ACI masking specification.
#' @param method CIR method. The forward value is the andreou2026cir
#'   last-exceedance objective applied to moser2026extremes event anchors
#'   (equivalent to moser2026extremes's first-crossing construction only for
#'   monotone lag-discrepancy rows); the backward andreou2026cir attribution is
#'   also a package synthesis rather than a moser2026extremes equation.
#' @param init Optional Gaussian initialization.
#' @param tol,window,max_lag Forward lag table controls, used only when
#'   `forward_table` is not supplied.
#' @param eps Optional subjective CIR discrepancy thresholds.
#' @param min_M Minimum discrepancy strength; `"auto"` uses the package option
#'   `aci.cir_min_strength` (default `1e-5`) and masks low-signal ranges.
#' @param ... Unused arguments are rejected.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#'
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' ev <- detect_events(ob, q = 0.9)
#' event_influence(m, ob, ev)
#' }
#'
#' @export
event_influence <- function(model, obs, events, forward_table = NULL,
                            nontarget = NULL,
                            method = c("exact", "l1_linf"), init = NULL,
                            tol = getOption("aci.default_tol", 1e-8),
                            window = 3L, max_lag = Inf,
                            eps = NULL, min_M = "auto", ...) {
  tol_supplied <- !missing(tol); window_supplied <- !missing(window)
  max_lag_supplied <- !missing(max_lag)
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to event_influence().")
  method <- match.arg(method)
  obs <- as_obs(obs)
  if (!is.data.frame(events) || !nrow(events) || !"t_star" %in% names(events) ||
      any(!is.finite(events$t_star)) ||
      any(events$t_star < min(obs$t) | events$t_star > max(obs$t)))
    aci_abort("aci_error_dims",
              "event_influence() needs finite event times within the observation span.")
  if (is.null(forward_table)) {
    nt_eff <- nontarget; init_eff <- init
    ft <- lag_table(model, obs, mode = "forward",
                    nontarget = nt_eff, init = init_eff, tol = tol,
                    window = window, max_lag = max_lag)
  } else {
    if (tol_supplied || window_supplied || max_lag_supplied)
      aci_abort("aci_error_dims", paste(
        "tol, window, and max_lag cannot be supplied when forward_table",
        "is reused; those settings are already fixed by the table."))
    nt_eff <- nontarget %||% forward_table$meta$nontarget
    init_eff <- init %||% forward_table$meta$init
    .validate_lag_table_source(forward_table, model, obs, nt_eff, init_eff)
    ft <- forward_table
  }
  fc <- forward_cir(ft, method = method, eps = eps, min_M = min_M)
  jmap <- vapply(events$t_star, function(ts) which.min(abs(obs$t - ts)), integer(1))
  bwd <- backward_cir(structure(list(handles = list(model = model, obs = obs,
                                                    nontarget = nt_eff,
                                                    init = init_eff)),
                                class = "aci_result"),
                      T = events$t_star, method = method,
                      eps = eps, min_M = min_M)
  tau_b <- if (length(jmap) == 1) bwd$tau else bwd$tau
  bint <- if (length(jmap) == 1L) matrix(bwd$interval, nrow = 1L) else bwd$interval
  data.frame(t_star = events$t_star,
             tau_forward = fc$tau[jmap], tau_backward = tau_b,
             window_lo = bint[, 1], backward_grid_end = bint[, 2],
             window_hi = events$t_star + fc$tau[jmap],
             forward_bound = rep(fc$bound, length(jmap)),
             forward_tail_estimate = fc$tail_bound[jmap],
             backward_bound = rep(bwd$bound, length(jmap)))
}


#' Composite conditional (event-locked) hidden statistics vs climatology
#'
#' @param paths A smoother path or list containing one.
#' @param events Event set from `detect_events()`.
#' @param w_pre Pre-event composite width.
#' @param w_post Post-event composite width.
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 30, dt = 0.01)
#' ob <- as_obs(sim)
#' ev <- detect_events(ob, q = 0.95, min_separation = 4, burn_in = 2)
#' f <- da_filter(m, ob)
#' s <- da_smooth(m, ob, filter = f)
#' event_stats(list(filter = f, smoother = s), ev)
#' }
#'
#' @export
event_stats <- function(paths, events, w_pre = 1.5, w_post = 1.5) {
  sm <- paths$smoother %||% paths
  if (!inherits(sm, "da_path_gaussian") || !identical(sm$kind, "smoother"))
    aci_abort("aci_error_dims", "paths must contain a Gaussian smoother path.")
  if (!is.data.frame(events) || !nrow(events) || !"t_star" %in% names(events))
    aci_abort("aci_error_dims", "event_stats() needs at least one event.")
  if (!is.numeric(w_pre) || length(w_pre) != 1L || !is.finite(w_pre) || w_pre < 0 ||
      !is.numeric(w_post) || length(w_post) != 1L || !is.finite(w_post) || w_post < 0)
    aci_abort("aci_error_dims", "w_pre and w_post must be finite non-negative durations.")
  if (any(!is.finite(events$t_star)))
    aci_abort("aci_error_dims", "Event times must be finite.")
  if (any(events$t_star < min(sm$t) | events$t_star > max(sm$t)))
    aci_abort("aci_error_dims", "Event times must lie within the smoother time span.")
  dt <- sm$t[2] - sm$t[1]; l <- ncol(sm$mean)
  off <- seq(-round(w_pre / dt), round(w_post / dt))
  rel_t <- off * dt
  centers <- vapply(events$t_star, function(ts) which.min(abs(sm$t - ts)), integer(1))
  complete <- centers + min(off) >= 1L & centers + max(off) <= length(sm$t)
  dropped <- sum(!complete)
  centers <- centers[complete]
  if (!length(centers))
    aci_abort("aci_error_dims",
              "No event has a complete requested pre/post window on the smoother grid.")
  cond_mean <- matrix(NA_real_, length(off), l)
  cond_cov <- array(NA_real_, c(l, l, length(off)))
  for (r in seq_along(off)) {
    idx <- centers + off[r]
    mus <- sm$mean[idx, , drop = FALSE]
    mu <- colMeans(mus)
    R <- matrix(0, l, l)
    for (ii in seq_along(idx)) {
      d <- mus[ii, ] - mu
      R <- R + sm$cov[, , idx[ii]] + tcrossprod(d)
    }
    cond_mean[r, ] <- mu
    cond_cov[, , r] <- sym(R / length(idx))
  }
  mu0 <- colMeans(sm$mean)
  R0 <- matrix(0, l, l)
  for (j in seq_along(sm$t)) {
    d <- sm$mean[j, ] - mu0
    R0 <- R0 + sm$cov[, , j] + tcrossprod(d)
  }
  R0 <- sym(R0 / length(sm$t))
  # vapply() returns a bare vector when l == 1, whose transpose is 1 x n rather
  # than the n x l this must be; the reshape fixes the orientation for every l.
  comp_sd <- t(matrix(vapply(seq_along(off), function(r)
    sqrt(pmax(diag(as.matrix(cond_cov[, , r])), 0)), numeric(l)), nrow = l))
  structure(list(rel_t = rel_t,
       conditional_mean = cond_mean, conditional_cov = cond_cov,
       composite_mean = cond_mean, composite_sd = comp_sd,
       climatology_mean = mu0, climatology_cov = R0,
       climatology_sd = sqrt(pmax(diag(R0), 0)),
       n_events = length(centers), n_events_dropped = dropped,
       event_indices_used = which(complete)), class = "aci_event_stats")
}


#' Event-sensitive directions
#'
#' Compare event-conditioned and unconditional hidden-state mixtures. The
#' primary direction uses a deterministic multistart numerical search for
#' moser2026extremes's full one-dimensional projected Gaussian KL on the
#' climatology-whitened unit sphere; it is not a certificate of the global
#' optimum. Separate mean-shift and covariance directions are returned as
#' diagnostics. This function does not construct moser2026extremes's
#' representative/MAP event paths or the weighted uncertainty quantities in
#' Eqs. 3.21-3.29.
#'
#' @param stats Event statistics returned by `event_stats()`.
#' @param at Relative event time or times.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 30, dt = 0.01)
#' ob <- as_obs(sim)
#' ev <- detect_events(ob, q = 0.95, min_separation = 4, burn_in = 2)
#' f <- da_filter(m, ob)
#' s <- da_smooth(m, ob, filter = f)
#' sensitive_directions(event_stats(list(filter = f, smoother = s), ev))
#' }
#'
#' @export
sensitive_directions <- function(stats, at = 0) {
  if (!inherits(stats, "aci_event_stats")) {
    aci_warn("aci_warn_deprecated", paste(
      "Passing filter/smoother paths to sensitive_directions() computes a",
      "posterior-contraction diagnostic, not the moser2026extremes event statistic; use",
      "posterior_contraction_directions() explicitly."))
    return(posterior_contraction_directions(stats, at))
  }
  if (!is.numeric(at) || !length(at) || any(!is.finite(at)))
    aci_abort("aci_error_dims", "at must contain finite relative times.")
  needed <- c("climatology_mean", "climatology_cov", "conditional_mean",
              "conditional_cov", "rel_t")
  if (!all(needed %in% names(stats)))
    aci_abort("aci_error_dims", "stats is missing event-mixture components.")
  mu0 <- stats$climatology_mean; R0 <- spd_floor(stats$climatology_cov)
  l <- length(mu0); C0 <- safe_chol(R0)
  canonical_sign <- function(V) {
    for (j in seq_len(ncol(V))) {
      ii <- which.max(abs(V[, j]))
      if (V[ii, j] < 0) V[, j] <- -V[, j]
    }
    V
  }
  full_direction <- function(RE, dm) {
    # Whiten v'R0v to one. With R0 = U' U and z = U v, the full projected
    # event KL is 1/2 { z'Az + (z'c)^2 - 1 - log(z'Az) } on ||z|| = 1,
    # where A = U'^-1 RE U^-1 and c = U'^-1 dm
    # (moser2026extremes eqs. 3.12, 3.17).
    L0 <- t(C0)
    left <- forwardsolve(L0, RE)
    A <- sym(t(forwardsolve(L0, t(left))))
    cvec <- drop(forwardsolve(L0, dm))
    score <- function(z) {
      nz <- sqrt(sum(z^2)); if (!is.finite(nz) || nz < 1e-14) return(-Inf)
      z <- z / nz; aa <- drop(crossprod(z, A %*% z))
      if (!is.finite(aa) || aa <= 0) return(-Inf)
      0.5 * (aa + sum(z * cvec)^2 - 1 - log(aa))
    }
    eeA <- eigen(A, symmetric = TRUE)
    starts <- eeA$vectors
    nc <- sqrt(sum(cvec^2))
    if (nc > 0) starts <- cbind(starts, cvec / nc)
    ns <- min(40L, max(8L, 4L * l))
    # Reshaped because vapply() returns a bare vector at l == 1, which cbind()
    # would then recycle down the column rather than append as starts.
    deterministic <- matrix(vapply(seq_len(ns), function(r)
      sin(seq_len(l) * (r + sqrt(2))), numeric(l)), nrow = l)
    starts <- cbind(starts, deterministic)
    best_z <- starts[, 1]; best <- score(best_z)
    fn <- function(w) -score(w)
    for (ii in seq_len(ncol(starts))) {
      fit <- tryCatch(stats::optim(starts[, ii], fn, method = "BFGS",
                                  control = list(maxit = 1000, reltol = 1e-12)),
                      error = function(e) NULL)
      z <- if (is.null(fit) || any(!is.finite(fit$par))) starts[, ii] else fit$par
      sc <- score(z)
      if (is.finite(sc) && sc > best) { best <- sc; best_z <- z }
    }
    best_z <- best_z / sqrt(sum(best_z^2))
    v <- drop(backsolve(C0, best_z))
    v <- v / sqrt(sum(v^2))
    if (v[which.max(abs(v))] < 0) v <- -v
    list(direction = v, score = best)
  }
  out <- lapply(at, function(rt) {
    j <- which.min(abs(stats$rel_t - rt))
    muE <- stats$conditional_mean[j, ]
    RE <- spd_floor(stats$conditional_cov[, , j])
    if (any(!is.finite(c(muE, RE))))
      aci_abort("aci_error_dims", "No complete event windows are available at the requested relative time.")
    dm <- muE - mu0
    fs <- full_direction(RE, dm)
    full_kl <- projected_kl(fs$direction, mu0, R0, muE, RE)
    vm <- drop(chol_solve(R0, dm, "event climatology covariance"))
    if (sqrt(sum(vm^2)) > 0) {
      vm <- vm / sqrt(sum(vm^2))
      if (vm[which.max(abs(vm))] < 0) vm <- -vm
      mean_kl <- projected_kl(vm, mu0, R0, muE, RE)
    } else {
      vm[] <- NA_real_; mean_kl <- NA_real_
    }
    L0 <- t(C0)
    left <- forwardsolve(L0, RE)
    whitened <- sym(t(forwardsolve(L0, t(left))))
    ee <- eigen(whitened, symmetric = TRUE)
    V <- backsolve(C0, ee$vectors)
    V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
    V <- canonical_sign(V)
    cov_score <- 0.5 * pmax(ee$values - 1 - log(ee$values), 0)
    ord <- order(cov_score, decreasing = TRUE)
    V <- V[, ord, drop = FALSE]; lam <- ee$values[ord]
    pk <- vapply(seq_len(l), function(ii)
      projected_kl(V[, ii], mu0, R0, muE, RE), numeric(1))
    list(relative_time = stats$rel_t[j],
         full_direction = fs$direction,
         full_projected_kl = full_kl,
         mean_direction = vm,
         mean_projected_kl = mean_kl,
         covariance_directions = V,
         directions = matrix(fs$direction, ncol = 1L,
                             dimnames = list(NULL, "full")),
         generalized_eigenvalues = lam,
         covariance_score = cov_score[ord], projected_kl = pk)
  })
  if (length(at) == 1L) out[[1]] else out
}


#' Filter-to-smoother uncertainty-contraction directions (package diagnostic,
#' not the moser2026extremes event-conditioned sensitive-direction estimand).
#'
#' @param paths List containing compatible Gaussian `filter` and `smoother` paths.
#' @param at Absolute time or times at which to calculate directions.
#' @returns A list of generalized eigen-decompositions, one per requested time,
#'   or that list's single element when one time is requested.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @seealso [sensitive_directions()], [event_stats()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' f <- da_filter(m, ob)
#' s <- da_smooth(m, ob, filter = f)
#' posterior_contraction_directions(list(filter = f, smoother = s), at = 1)
#'
#' @export
posterior_contraction_directions <- function(paths, at) {
  if (!is.list(paths) ||
      !inherits(paths$filter, "da_path_gaussian") ||
      !inherits(paths$smoother, "da_path_gaussian") ||
      !identical(paths$filter$kind, "filter") ||
      !identical(paths$smoother$kind, "smoother"))
    aci_abort("aci_error_dims",
              "paths must contain Gaussian filter and smoother paths.")
  filt <- paths$filter; smoo <- paths$smoother
  if (!is.numeric(at) || !length(at) || any(!is.finite(at)) ||
      any(at < min(filt$t) | at > max(filt$t)))
    aci_abort("aci_error_dims", "at must contain finite times inside the path span.")
  # Reuse the public path validator for grid, dimension, finiteness, and SPD
  # checks before solving a generalized eigenproblem.
  invisible(gaussian_kl_path(smoo, filt, decompose = FALSE))
  l <- ncol(filt$mean)
  out <- lapply(at, function(ts) {
    j <- which.min(abs(filt$t - ts))
    Rf <- filt$cov[, , j]; Rs <- smoo$cov[, , j]
    Wc <- safe_chol(Rs); A <- forwardsolve(t(Wc), t(forwardsolve(t(Wc), Rf)))
    e <- eigen(sym(A), symmetric = TRUE)
    V <- backsolve(Wc, e$vectors)                 # generalized eigvecs Rf v = lam Rs v
    V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
    pk <- vapply(seq_len(l), function(c1)
      projected_kl(V[, c1], filt$mean[j, ], Rf, smoo$mean[j, ], Rs), numeric(1))
    ord <- order(pk, decreasing = TRUE)
    list(t = filt$t[j], directions = V[, ord, drop = FALSE],
         projected_kl = pk[ord], gen_eigenvalues = e$values[ord])
  })
  if (length(at) == 1) out[[1]] else out
}


################################################################################
# moser2026extremes App A/B feature sets and mechanism clustering
################################################################################

#' Indices of a closed time window (internal)
#'
#' @param t Numeric vector of times.
#' @param lo Lower bound of the window.
#' @param hi Upper bound of the window.
#' @returns Integer vector of indices of `t` inside the window.
#' @noRd
.win_idx <- function(t, lo, hi) which(t >= lo - 1e-12 & t <= hi + 1e-12)


#' Ordinary least-squares slope over a window (internal)
#'
#' @param t Numeric vector of times.
#' @param y Numeric vector of values, one per time.
#' @returns 1-length numeric slope, or `NA` when fewer than two points.
#' @noRd
.slope <- function(t, y) if (length(t) > 1) unname(stats::coef(stats::lm(y ~ t))[2]) else NA_real_


#' Check that feature windows fit inside the record (internal)
#'
#' @param sim An `aci_sim` object.
#' @param events An `event_set` object.
#' @param pre Length of the pre-event window.
#' @param post Length of the post-event window.
#' @returns The retained events, as an integer index vector.
#' @noRd
.validate_feature_events <- function(sim, events, pre, post) {
  if (is.null(sim$obs) || !inherits(sim$obs, "obs_traj") ||
      !is.data.frame(events) || !nrow(events) ||
      !all(c("peak_index", "t_star") %in% names(events)))
    aci_abort("aci_error_dims",
              "Feature extraction needs a simulation and non-empty events with peak_index/t_star.")
  pk <- events$peak_index
  if (!is.numeric(pk) || any(!is.finite(pk)) || any(pk != floor(pk)) ||
      any(pk < 1L) || any(pk > length(sim$obs$t)) ||
      any(!is.finite(events$t_star)))
    aci_abort("aci_error_dims", "Event peak indices/times are invalid.")
  pk <- as.integer(pk)
  if (any(abs(sim$obs$t[pk] - events$t_star) > sim$obs$dt / 2 + 1e-12))
    aci_abort("aci_error_dims", "Event peak indices do not match this simulation grid.")
  if (any(events$t_star - pre < min(sim$obs$t) - 1e-12) ||
      any(events$t_star + post > max(sim$obs$t) + 1e-12))
    aci_abort("aci_error_dims",
              "Every event must have the complete feature pre/post window on the simulation grid.")
  invisible(pk)
}


#' Hidden-state matrix used for feature extraction (internal)
#'
#' @param sim An `aci_sim` object.
#' @param source Either `"hidden"` for the true path or `"smoother"` for the
#'   reconstructed one.
#' @param paths Optional assimilation paths, required when `source` is
#'   `"smoother"`.
#' @returns Numeric matrix of hidden values, one row per time.
#' @noRd
.feature_hidden_matrix <- function(sim, source, paths) {
  if (is.null(sim$model) || !inherits(sim$model, "stochastic_model"))
    aci_abort("aci_error_dims", "Feature extraction needs the simulation model metadata.")
  if (source == "hidden") {
    H <- sim$hidden
    if (!is.matrix(H) || !identical(dim(H), c(length(sim$obs$t), sim$model$l)) ||
        any(!is.finite(H)))
      aci_abort("aci_error_dims",
                "The simulation does not contain a compatible finite hidden truth path.")
    return(H)
  }
  sm <- paths$smoother %||% paths
  if (!inherits(sm, "da_path_gaussian"))
    aci_abort("aci_error_dims", "source = 'smoother' needs a Gaussian smoother path.")
  .validate_gaussian_path(sm, sim$obs, sim$model$l, "smoother",
                          sm$meta$nontarget %||% NULL,
                          source_model = sim$model)
  sm$mean
}


#' Eighteen-feature vector for the pathways model (moser2026extremes Appendix A)
#'
#' The pre-event window runs from `t_star - 1.5` to `t_star`.
#'
#' @param sim Simulation object containing observations and, when requested, hidden truth.
#' @param events Event set.
#' @param source Hidden truth or smoother reconstruction.
#' @param paths Assimilation paths required for smoother-derived features.
#' @param params Optional model parameters overriding simulation metadata.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' sp <- simulate(model_pathways(), seed = 1, T = 20, dt = 0.01)
#' ev <- detect_events(as_obs(sp), q = 0.95, min_separation = 3)
#' head(features_pathways(sp, ev))
#'
#' @export
features_pathways <- function(sim, events, source = c("hidden", "smoother"),
                              paths = NULL, params = NULL) {
  source <- match.arg(source)
  .validate_feature_events(sim, events, pre = 1.5, post = 0)
  p <- params %||% sim$model$meta$params
  t <- sim$obs$t; u <- sim$obs$x[, 1]; dt <- sim$obs$dt
  H <- .feature_hidden_matrix(sim, source, paths)
  gam <- H[, 1]; b <- H[, 2]
  ft <- vapply(seq_len(nrow(events)), function(i) {
    pk <- events$peak_index[i]; ts <- events$t_star[i]; sgn <- events$sign[i]
    w <- .win_idx(t, ts - 1.5, ts)
    G <- (-p$d_u + p$c * gam[w]) * u[w]; B <- b[w]
    iG <- sum(G) * dt; iB <- sum(B) * dt
    aG <- sum(abs(G)) * dt; aB <- sum(abs(B)) * dt
    j1 <- which.min(abs(t - (ts - 1)))
    den <- max(aG + aB, .Machine$double.eps)
    c(peak = u[pk], abs_peak = abs(u[pk]), duration = events$duration[i],
      pre_slope_u = .slope(t[w], u[w]), pre_mean_u = mean(u[w]),
      pre_mean_gamma = mean(gam[w]), pre_mean_b = mean(b[w]),
      pre_mean_G = mean(G), pre_mean_B = mean(B),
      int_G = iG, int_B = iB, signed_int_G = sgn * iG, signed_int_B = sgn * iB,
      rel_damping = aG / den, rel_forcing = aB / den,
      u_lag1 = u[j1], gamma_lag1 = gam[j1], b_lag1 = b[j1])
  }, numeric(18))
  out <- t(ft)
  if (any(!is.finite(out)))
    aci_abort("aci_error_dims", "Pathway features are non-finite; check event windows and model metadata.")
  out
}


#' Twenty-nine-feature vector for the topographic model (moser2026extremes
#' Appendix B)
#'
#' Relative to `t_star`, the pre window is -1 to -0.3, the peak window is
#' -0.1 to 0.1, and the post window is 0.2 to 1.
#' S_pre^(s)/S_pre^(l) growth-rate windows are chosen as (t*-0.3, t*) /
#' (t*-1, t*): an implementation choice pending exact source definitions
#' (transcription log); family energies use the listed positive-k members.
#'
#' @param sim Simulation object containing observations and modal states.
#' @param events Event set.
#' @param source Hidden truth or smoother reconstruction.
#' @param paths Assimilation paths required for smoother-derived features.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' \donttest{
#' st <- simulate(model_topographic(), seed = 1, T = 20, dt = 0.01)
#' ev <- detect_events(as_obs(st), q = 0.9, min_separation = 3, burn_in = 2)
#' head(features_topographic(st, ev))
#' }
#'
#' @export
features_topographic <- function(sim, events, source = c("hidden", "smoother"),
                                 paths = NULL) {
  source <- match.arg(source)
  .validate_feature_events(sim, events, pre = 1, post = 1)
  t <- sim$obs$t; V <- sim$obs$x[, 1]
  H <- .feature_hidden_matrix(sim, source, paths)
  Kp <- sim$model$meta$modes; nK <- nrow(Kp)
  key <- apply(Kp, 1, function(k) paste(k[1], k[2], sep = "_"))
  Z <- H[, 1:nK, drop = FALSE]; Zi <- H[, nK + 1:nK, drop = FALSE]
  # The state stores one member of each conjugate Fourier pair. Full modal
  # energy therefore has a factor of two (there is no zero mode in this set).
  E_k <- 2 * (Z^2 + Zi^2)
  fam <- sim$model$meta$mode_families
  fidx <- lapply(fam, function(nm) match(nm, key))
  Etot <- rowSums(E_k)
  ft <- vapply(seq_len(nrow(events)), function(i) {
    ts <- events$t_star[i]; pk <- events$peak_index[i]
    wpre  <- .win_idx(t, ts - 1,   ts - 0.3)
    wpost <- .win_idx(t, ts + 0.2, ts + 1)
    wsg   <- .win_idx(t, ts - 0.3, ts)      # short growth window (impl. choice)
    wlg   <- .win_idx(t, ts - 1,   ts)      # long growth window (impl. choice)
    Efam <- vapply(fidx, function(ii) mean(rowSums(E_k[wpre, ii, drop = FALSE])), numeric(1))
    Emax <- vapply(fidx[c("upper", "zonal", "lower")],
                   function(ii) max(rowSums(E_k[wpre, ii, drop = FALSE])), numeric(1))
    re_pre  <- colMeans(Z[wpre, match(c("1_1", "2_1", "1_0", "1_-1", "2_-1"), key), drop = FALSE])
    c(V_peak = V[pk], abs_V_peak = abs(V[pk]),
      S_pre_s = .slope(t[wsg], V[wsg]), S_pre_l = .slope(t[wlg], V[wlg]),
      S_post = .slope(t[wpost], V[wpost]),
      Std_pre_V = stats::sd(V[wpre]), Mean_pre_V = mean(V[wpre]),
      Mean_post_V = mean(V[wpost]),
      E_tot_peak = Etot[pk], E_tot_pre_mean = mean(Etot[wpre]),
      E_tot_pre_max = max(Etot[wpre]),
      E_upper = Efam[["upper"]], E_zonal = Efam[["zonal"]],
      E_lower = Efam[["lower"]], E_merid = Efam[["merid"]],
      E_upper_max = Emax[["upper"]], E_zonal_max = Emax[["zonal"]],
      E_lower_max = Emax[["lower"]],
      R1 = Efam[["upper"]] / Efam[["lower"]],
      R2 = Efam[["zonal"]] / (Efam[["upper"]] + Efam[["lower"]]),
      R3 = (Efam[["upper"]] - Efam[["lower"]]) / (Efam[["upper"]] + Efam[["lower"]]),
      re_pre_11 = re_pre[1], re_pre_21 = re_pre[2], re_pre_10 = re_pre[3],
      re_pre_1m1 = re_pre[4], re_pre_2m1 = re_pre[5],
      re_peak_11 = Z[pk, match("1_1", key)],
      re_peak_10 = Z[pk, match("1_0", key)],
      re_peak_1m1 = Z[pk, match("1_-1", key)])
  }, numeric(29))
  out <- t(ft)
  if (any(!is.finite(out)))
    aci_abort("aci_error_dims",
              "Topographic features are non-finite; check complete windows and energy denominators.")
  out
}


#' Default feature names for the event classifiers
#'
#' Request the default feature names used by the event classifiers.
#'
#' @param kind Either `"pathways"` or `"topographic"`.
#' @returns Character vector of feature names.
#'
#' @seealso [features_pathways()], [features_topographic()], [classify_events()]
#'
#' @examples
#' feature_set_default("pathways")
#'
#' @export
feature_set_default <- function(kind = c("pathways", "topographic"))
  switch(match.arg(kind), pathways = features_pathways,
         topographic = features_topographic)


#' Standardize + k-means + PCA projection (moser2026extremes s3.5 / Figs. 4.3d)
#'
#' @param features Numeric event-by-feature matrix.
#' @param k Number of clusters.
#' @param standardize Whether to standardize feature columns.
#' @param seed Reproducibility seed.
#' @param nstart Number of k-means random starts.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' set.seed(1)
#' classify_events(matrix(rnorm(60), ncol = 3), k = 2)
#'
#' @export
classify_events <- function(features, k = 3, standardize = TRUE, seed = 1,
                            nstart = 25) {
  X <- as.matrix(features)
  if (nrow(X) < 2L || ncol(X) < 1L || any(!is.finite(X)))
    aci_abort("aci_error_dims", "features must be a finite matrix with at least two events.")
  X <- X[, apply(X, 2, function(c1) all(is.finite(c1)) && stats::sd(c1) > 0), drop = FALSE]
  if (!ncol(X)) aci_abort("aci_error_dims", "No non-constant finite features remain.")
  if (length(k) != 1L || !is.finite(k) || k != as.integer(k) || k < 1L || k >= nrow(X))
    aci_abort("aci_error_dims", "k must be a positive integer smaller than the number of events.")
  if (length(nstart) != 1L || !is.finite(nstart) || nstart < 1L || nstart != as.integer(nstart))
    aci_abort("aci_error_dims", "nstart must be a positive integer.")
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize))
    aci_abort("aci_error_dims", "standardize must be TRUE or FALSE.")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    aci_abort("aci_error_dims", "seed must be one non-negative integer.")
  Xs <- if (standardize) scale(X) else X
  set.seed(seed)
  km <- stats::kmeans(Xs, centers = k, nstart = nstart)
  pc <- stats::prcomp(Xs)
  structure(list(cluster = km$cluster, centers = km$centers, kmeans = km,
                 pca = pc, projection = pc$x[, 1:min(2, ncol(pc$x)), drop = FALSE],
                 features_used = colnames(X)),
            class = "event_classes")
}


#' Print an event classification
#'
#' @param x An `event_classes` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.event_classes <- function(x, ...) {
  cat(sprintf("<event_classes> k = %d clusters over %d events, %d features\n",
              nrow(x$centers), length(x$cluster), length(x$features_used)))
  print(table(cluster = x$cluster)); invisible(x)
}
