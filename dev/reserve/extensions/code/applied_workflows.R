## acir reserve file
## Origin: aci/R/applied_workflows.R
## Source package: aci 0.0.30, git tree 97f6b124
## Category: extensions
## Intended release: unscheduled (package-only, no paper or MATLAB backing)
## Reason: Applied on-ramps (screen_outliers, calibrate_cgns, ...); package-only infrastructure.
## Verbatim copy from the aci 0.0.30 sources; not modified.

################################################################################
## applied_workflows.R - reusable, non-estimand-specific application workflows
## ########################################################################## ##
##
## Contents:
##   * observation preprocessing for the engine's uniform trajectory contract:
##       - screen_outliers, as_uniform_trajectory
##
##   * deterministic CGNS simulation and bounded calibration:
##       - simulate_deterministic, calibrate_cgns, print.cgns_calibration
##
################################################################################


#' Screen one-sided dropouts from an irregular observation series
#'
#' Rejects points that fall implausibly far BELOW a deliberately stiff fit, plus
#' anything outside the physical range `[lo, hi]`. The stiffness matters: left to
#' generalised cross-validation a spline with more than a dozen dates interpolates
#' straight through the dropouts, residuals stay small and nothing is rejected.
#' `df_max` therefore caps the flexibility of the screening fit only.
#'
#' Returns a logical vector over the input, with the fitted values attached.
#' @param t Observation times.
#' @param value Numeric observation series.
#' @param lo,hi Physical lower and upper bounds.
#' @param k Robust residual cutoff multiplier.
#' @param floor_abs Minimum absolute residual cutoff.
#' @param iter Number of screen/refit iterations.
#' @param df_max Maximum screening-spline degrees of freedom.
#'
#' @examples
#' set.seed(1)
#' tt <- seq(0, 10, length.out = 30)
#' v <- sin(tt) + rnorm(30, sd = 0.05)
#' v[c(5, 20)] <- v[c(5, 20)] - 1
#' screen_outliers(tt, v)
#'
#' @export
screen_outliers <- function(t, value, lo = -Inf, hi = Inf, k = 3, floor_abs = 0.10,
                            iter = 3L, df_max = 5) {
  t <- as.numeric(t); value <- as.numeric(value)
  if (length(t) != length(value))
    aci_abort("aci_error_dims", "screen_outliers(): t and value must be the same length.")
  if (!is.numeric(iter) || length(iter) != 1L || !is.finite(iter) ||
      iter < 0L || iter != floor(iter))
    aci_abort("aci_error_dims", "iter must be a non-negative integer.")
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0 ||
      !is.numeric(floor_abs) || length(floor_abs) != 1L ||
      !is.finite(floor_abs) || floor_abs < 0)
    aci_abort("aci_error_dims", "k and floor_abs must be finite non-negative scalars.")
  if (!is.numeric(lo) || length(lo) != 1L || is.na(lo) ||
      !is.numeric(hi) || length(hi) != 1L || is.na(hi) || lo > hi)
    aci_abort("aci_error_dims", "lo and hi must be ordered numeric scalar bounds.")
  if (!is.numeric(df_max) || length(df_max) != 1L || is.na(df_max) ||
      df_max < 3)
    aci_abort("aci_error_dims", "df_max must be at least 3 or Inf.")
  iter <- as.integer(iter)
  keep <- is.finite(t) & is.finite(value) & value >= lo & value <= hi
  if (sum(keep) < 4L) return(keep)
  for (i in seq_len(iter)) {
    fit <- .aci_spline(t[keep], value[keep], df_max = df_max)
    r   <- value[keep] - stats::predict(fit, t[keep])$y
    thr <- max(k * stats::mad(r), floor_abs)
    drop <- r < -thr
    if (!any(drop)) break
    keep[which(keep)[drop]] <- FALSE
    if (sum(keep) < 4L) break
  }
  attr(keep, "n_dropped") <- sum(is.finite(value)) - sum(keep)
  keep
}


# a smoothing spline with a flexibility cap; GCV where the data support it
#' Smoothing spline with a bounded degrees-of-freedom cap (internal)
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param df_max Upper bound on the fitted degrees of freedom.
#' @returns The fitted smoothing-spline object.
#' @noRd
.aci_spline <- function(x, y, df_max = Inf) {
  df <- min(df_max, length(x) - 2)
  if (is.finite(df) && df >= 3) stats::smooth.spline(x, y, df = df)
  else if (length(x) >= 10)     stats::smooth.spline(x, y)
  else stats::smooth.spline(x, y, df = max(3, min(5, length(x) - 2)))
}


#' Resample irregular observation channels onto the uniform grid the engine needs
#'
#' `channels` is a named list of numeric vectors sharing the times `t`. Each is
#' smoothed and evaluated on a uniform grid of step `dt`, and returned as an
#' `obs_traj`. `keep` (e.g. from `screen_outliers`) subsets the input first, so a
#' dropout screened on one channel removes that DATE from every channel.
#' @param t Observation times.
#' @param channels Named numeric observation channels sharing `t`.
#' @param dt Requested uniform-grid spacing.
#' @param keep Optional logical row selector applied to all channels.
#' @param df_max Maximum smoothing-spline degrees of freedom.
#' @param floor Optional lower floor applied after interpolation.
#'
#' @examples
#' set.seed(1)
#' tt <- sort(runif(30, 0, 10))
#' as_uniform_trajectory(tt, list(ndvi = sin(tt)), dt = 0.5)
#'
#' @export
as_uniform_trajectory <- function(t, channels, dt, keep = NULL, df_max = Inf,
                                  floor = NULL) {
  t <- as.numeric(t)
  if (!is.list(channels) || !length(channels) || is.null(names(channels)) ||
      anyNA(names(channels)) || any(!nzchar(names(channels))) ||
      anyDuplicated(names(channels)))
    aci_abort("aci_error_obs_contract",
              paste(
                "as_uniform_trajectory() needs a list with unique, non-empty",
                "observation-channel names."))
  n_input <- length(t)
  valid_channel <- vapply(channels, function(v)
    is.numeric(v) && length(v) == n_input, logical(1))
  if (any(!valid_channel))
    aci_abort("aci_error_dims", paste(
      "Every observation channel must be numeric and have the same length",
      "as t."))
  if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0)
    aci_abort("aci_error_dims", "dt must be one finite positive number.")
  if (!is.numeric(df_max) || length(df_max) != 1L || is.na(df_max) ||
      df_max < 3)
    aci_abort("aci_error_dims", "df_max must be at least 3 or Inf.")
  if (!is.null(keep)) {
    if (!is.logical(keep) || length(keep) != n_input || anyNA(keep))
      aci_abort("aci_error_dims",
                "keep must be a non-missing logical vector with one value per input time.")
    t <- t[keep]
    channels <- lapply(channels, function(v) v[keep])
  }
  if (length(t) < 4L)
    aci_abort("aci_error_obs_contract",
              sprintf("as_uniform_trajectory(): only %d usable observations after screening.",
                      length(t)))
  values <- unlist(channels, use.names = FALSE)
  if (any(!is.finite(t)) || any(!is.finite(values)))
    aci_abort("aci_error_obs_contract", paste(
      "Times and retained observation-channel values must be finite; use",
      "keep to remove screened dates first."))
  if (any(diff(t) <= 0))
    aci_abort("aci_error_obs_contract",
              "Retained observation times must be unique and strictly increasing.")
  grid <- seq(min(t), max(t), by = dt)
  if (length(grid) < 2L)
    aci_abort("aci_error_obs_contract",
              "dt is too large to produce an observation trajectory over this time span.")
  X <- vapply(channels, function(v)
    as.numeric(stats::predict(.aci_spline(t, as.numeric(v), df_max), grid)$y),
    numeric(length(grid)))
  X <- matrix(X, nrow = length(grid), dimnames = list(NULL, names(channels)))
  if (!is.null(floor)) for (nm in names(floor)) if (nm %in% colnames(X))
    X[, nm] <- pmax(X[, nm], floor[[nm]])
  ob <- observed_trajectory(grid, X, names = names(channels))
  ob$n_obs <- length(t)
  ob
}


#' Deterministic (noise-free) solution of any CGNS model
#'
#' Integrates the drift of a `cgns_model` with the hidden state carried alongside
#' the observed one (the object a calibration compares against data). Generic:
#' it reads the model's own coefficient blocks, so it works for any CGNS, not
#' just the crop model. Returns NULL if the trajectory leaves the finite range.
#' @param model CGNS model.
#' @param times Increasing integration times.
#' @param x0 Initial observed state.
#' @param y0 Initial hidden state.
#'
#' @examples
#' simulate_deterministic(model_dyad(), times = seq(0, 1, by = 0.05),
#'                        x0 = 1, y0 = 2)
#'
#' @export
simulate_deterministic <- function(model, times, x0, y0) {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "simulate_deterministic() needs a cgns_model.")
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
      (length(times) > 1L && any(diff(times) <= 0)))
    aci_abort("aci_error_dims",
              "times must be a non-empty, finite, strictly increasing numeric vector.")
  if (!is.numeric(x0) || length(x0) != model$k || any(!is.finite(x0)) ||
      !is.numeric(y0) || length(y0) != model$l || any(!is.finite(y0)))
    aci_abort("aci_error_dims", paste(
      "x0 and y0 must contain one finite numeric value per observed and",
      "hidden state, respectively."))
  n <- length(times); X <- matrix(NA_real_, n, model$k); Y <- matrix(NA_real_, n, model$l)
  X[1, ] <- x0; Y[1, ] <- y0
  for (j in seq_len(n - 1L)) {
    h  <- times[j + 1L] - times[j]
    co <- eval_coefs(model, times[j], X[j, ])
    X[j + 1L, ] <- X[j, ] + (co$fx + drop(co$Lx %*% Y[j, ])) * h
    Y[j + 1L, ] <- Y[j, ] + (co$fy + drop(co$Ly %*% Y[j, ])) * h
    if (!all(is.finite(c(X[j + 1L, ], Y[j + 1L, ])))) return(NULL)
  }
  list(t = times, x = X, y = Y)
}


#' Calibrate a CGNS model's rate parameters against an observed trajectory
#'
#' Bounded, multi-start least squares of the model's deterministic solution
#' against the observed channels. `build(par)` returns a `cgns_model`; only the
#' names in `bounds` are fitted, everything else in `start` is held fixed, so
#' the cause's amplitude (the pool at sowing, the dated inputs) stays pinned to
#' the management record rather than being absorbed into the fit.
#'
#' Bounds are not cosmetic. Left unconstrained a canopy series will drive the
#' pool-loss rate to zero, and a non-decaying pool violates the finite-memory
#' premise the CIR rests on: attribution windows then reach the start of the
#' record at every reference time. Parameters resting on a bound are reported in
#' `on_bound`, because they are not identified by the data.
#'
#' Process-noise amplitudes are estimated from the increment residuals about the
#' calibrated deterministic path and floored, so the observed-channel Gram stays
#' well conditioned for the Riccati step.
#' @param build Function mapping a named parameter vector to a CGNS model.
#' @param obs Observed trajectory.
#' @param start Named starting parameter values.
#' @param bounds Named lower/upper bounds for fitted parameters.
#' @param y0 Initial hidden state used by deterministic integration.
#' @param weights Optional named observation-channel weights.
#' @param restarts Number of bounded optimization starts.
#' @param sigma_floor Minimum inferred process-noise amplitude.
#' @param maxit Maximum optimizer iterations.
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' build <- function(p) model_dyad(params = list(
#'   d_x = p[["d_x"]], gamma = 2, f_x = 0.5, s_x = 0.5,
#'   d_y = 0.5, f_y = 1, s_y = 1))
#' calibrate_cgns(build, ob, start = c(d_x = 0.4),
#'                bounds = list(d_x = c(0.1, 1)), y0 = 2, restarts = 1, maxit = 5)
#' }
#'
#' @export
calibrate_cgns <- function(build, obs, start, bounds, y0, weights = NULL,
                           restarts = c(1, 2, 0.5), sigma_floor = 0.01,
                           maxit = 4000L) {
  obs <- as_obs(obs)
  nm  <- names(bounds)
  if (!length(nm) || anyDuplicated(nm) || !all(nm %in% names(start)))
    aci_abort("aci_error_model_contract",
              "calibrate_cgns(): every name in `bounds` must appear in `start`.")
  valid_bounds <- vapply(bounds, function(b)
    length(b) == 2L && all(is.finite(b)) && b[1] < b[2], logical(1))
  if (any(!valid_bounds))
    aci_abort("aci_error_model_contract",
              "Every bound must be a finite c(lower, upper) pair with lower < upper.")
  if (length(sigma_floor) != 1L || !is.finite(sigma_floor) || sigma_floor <= 0)
    aci_abort("aci_error_dims", "sigma_floor must be a finite positive scalar.")
  Y  <- obs$x; tt <- obs$t; dt <- obs$dt; k <- ncol(Y)
  channel_names <- colnames(Y) %||% paste0("x", seq_len(k))
  w <- weights %||% rep(1, k)
  if (!is.null(weights) && !is.null(names(weights))) {
    if (anyDuplicated(names(weights)) || !setequal(names(weights), channel_names))
      aci_abort("aci_error_dims",
                "Named weights must name every observed channel exactly once.")
    w <- weights[channel_names]
  }
  if (length(w) != k || any(!is.finite(w)) || any(w < 0))
    aci_abort("aci_error_dims",
              "weights must contain one finite, non-negative value per observed channel.")
  if (!any(w > 0))
    aci_abort("aci_error_dims", "At least one observation-channel weight must be positive.")
  w <- as.numeric(w)
  sc <- apply(Y, 2, stats::sd); sc[!is.finite(sc) | sc == 0] <- 1
  to_p <- function(v, b) stats::qlogis(min(max((v - b[1]) / (b[2] - b[1]), 1e-4), 1 - 1e-4))
  fr_p <- function(p, b) b[1] + (b[2] - b[1]) * stats::plogis(p)
  set_par <- function(lp) { th <- start
    for (i in seq_along(nm)) th[[nm[i]]] <- fr_p(lp[i], bounds[[nm[i]]]); th }
  obj <- function(lp) {
    m <- tryCatch(build(set_par(lp)), error = function(e) NULL)
    if (is.null(m)) return(1e12)
    r <- tryCatch(simulate_deterministic(m, tt, Y[1, ], y0),
                  error = function(e) NULL)
    if (is.null(r)) return(1e12)
    sum(vapply(seq_len(k), function(i)
      w[i] * sum(((r$x[, i] - Y[, i]) / sc[i])^2), numeric(1)))
  }
  best <- NULL
  optim_method <- if (length(nm) == 1L) "BFGS" else "Nelder-Mead"
  for (mult in restarts) {
    p0 <- vapply(nm, function(n) to_p(start[[n]] * mult, bounds[[n]]), numeric(1))
    o  <- stats::optim(p0, obj, method = optim_method,
                       control = list(maxit = maxit, reltol = 1e-11))
    o  <- stats::optim(o$par, obj, method = optim_method,
                       control = list(maxit = maxit, reltol = 1e-13))
    if (is.null(best) || o$value < best$value) best <- o
  }
  par <- set_par(best$par)
  base_model <- build(par)
  det <- simulate_deterministic(base_model, tt, Y[1, ], y0)
  ## increment-residual noise amplitudes, per observed channel
  sig <- vapply(seq_len(k), function(i) {
    drift <- diff(det$x[, i]) / dt
    max(stats::sd(diff(Y[, i]) / dt - drift) * sqrt(dt), sigma_floor)
  }, numeric(1))
  names(sig) <- channel_names

  # Carry the fitted marginal observed-process amplitudes into the returned
  # model.  A generic builder need not expose its diffusion parameters by
  # name, so rescale the rows of both observed diffusion blocks.  This exactly
  # installs `sig` for constant diagonal models (including the crop models),
  # and preserves the builder's state dependence and Wiener-channel structure
  # for more general CGNS models.  Cross-channel correlations are retained but
  # are not estimated by this marginal residual calculation.
  g0 <- cgns_grams(base_model, tt[1], Y[1, ])$gxx
  base_sig <- sqrt(pmax(diag(g0), 0))
  if (length(base_sig) != k || any(!is.finite(base_sig)) || any(base_sig <= 0))
    aci_abort("aci_error_gram", paste(
      "calibrate_cgns() cannot carry residual noise into a model whose",
      "observed diffusion has a zero or non-finite marginal amplitude."))
  row_scale <- unname(sig / base_sig)
  scale_rows <- function(A) sweep(as.matrix(A), 1L, row_scale, `*`)
  fitted_model <- cgns_model(
    Lx = base_model$Lx, fx = base_model$fx,
    Ly = base_model$Ly, fy = base_model$fy,
    Sx1 = function(t, x) scale_rows(base_model$Sx1(t, x)),
    Sx2 = function(t, x) scale_rows(base_model$Sx2(t, x)),
    Sy1 = base_model$Sy1, Sy2 = base_model$Sy2,
    k = base_model$k, l = base_model$l, name = base_model$name,
    meta = base_model$meta)
  fitted_model$meta$calibrated_observed_noise <- sig
  fitted_model$meta$calibration_noise_method <-
    "row-scaled builder diffusion; marginal amplitudes estimated from increment residuals"
  fitted_model$meta$source_status <- paste(
    "Experimental calibration/application extension; outside the supplied",
    "ACI papers and MATLAB reference code.")

  hit <- vapply(nm, function(n) { b <- bounds[[n]]
    abs(par[[n]] - b[1]) < 1e-3 * diff(b) || abs(par[[n]] - b[2]) < 1e-3 * diff(b) },
    logical(1))
  rmse_by_channel <- sqrt(colMeans((det$x - Y)^2))
  names(rmse_by_channel) <- channel_names
  structure(list(par = par, fitted = nm, bounds = bounds, det = det, sigma = sig,
                 on_bound = nm[hit], value = best$value,
                 rmse = sqrt(mean((det$x - Y)^2)),
                 rmse_by_channel = rmse_by_channel,
                 model = fitted_model,
                 meta = list(
                   source_status = paste(
                     "Experimental calibration/application extension; outside",
                     "the supplied ACI papers and MATLAB reference code."),
                   noise_calibration = paste(
                     "Marginal observed-process amplitudes only; cross-channel",
                     "noise correlations and hidden-process noise were not fitted."))),
            class = "cgns_calibration")
}


#' Print a CGNS calibration result
#'
#' @param x A `cgns_calibration` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.cgns_calibration <- function(x, ...) {
  cat(sprintf("<cgns_calibration> %d fitted parameter%s | RMSE %.4g\n",
              length(x$fitted), if (length(x$fitted) == 1) "" else "s", x$rmse))
  cat("  ", paste(sprintf("%s=%.4g", x$fitted, unlist(x$par[x$fitted])),
                  collapse = "  "), "\n", sep = "")
  if (length(x$on_bound))
    cat("  on bound (not identified by the data): ",
        paste(x$on_bound, collapse = ", "), "\n", sep = "")
  cat("  process noise: ", paste(sprintf("%s=%.3g", names(x$sigma), x$sigma),
                                 collapse = "  "), "\n", sep = "")
  if (length(x$rmse_by_channel) > 1L)
    cat("  RMSE by channel: ",
        paste(sprintf("%s=%.4g", names(x$rmse_by_channel), x$rmse_by_channel),
              collapse = "  "), "\n", sep = "")
  invisible(x)
}
