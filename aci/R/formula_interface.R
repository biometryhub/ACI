################################################################################
## formula_interface.R - exported: aci_fit(), generics, plots
## ########################################################################## ##
##
## Contents:
##   * aci_fit() formula front-end and the standard generics:
##       - aci_fit, .fit_prepare, .fit_full, .fit_residual_sds, .fit_template,
##         .fit_partial, print.aci_fit, summary.aci_fit, print.summary.aci_fit,
##         coef.aci_fit, nobs.aci_fit, predict.aci_fit, fitted.aci_fit,
##         residuals.aci_fit, simulate.aci_fit, plot.aci_fit, cir, cir.default,
##         cir.aci_fit
##
##   * base-graphics plot methods for the core classes:
##       - plot.aci_result, plot.da_path_gaussian, plot.cir_result
##
################################################################################


################################################################################
#  aci_fit() formula front-end and the standard generics
################################################################################

#' Fit a latent-cause model from a formula
#'
#' `aci_fit()` provides a modelling front end for a latent-cause formula.
#' Partial-observation learning is experimental and is labelled as such in
#' fitted objects and summaries.
#'
#' @param formula A formula naming the observed variables on the left and the
#'   latent variables on the right.
#' @param data A data frame, or a list of data frames, holding the observed
#'   series.
#' @param time Optional name of the time column; `NULL` assumes a regular grid.
#' @param group Optional name of a grouping column splitting `data` into
#'   segments.
#' @param degree Positive whole number, the polynomial degree of the candidate
#'   library.
#' @param ce_threshold Either `"auto"` or a numeric causation entropy threshold.
#' @param energy_pairs Optional list of index pairs constrained to conserve
#'   energy across equations.
#' @param n_samples Number of posterior path samples drawn when a latent
#'   variable is present.
#' @param hidden_prior List with the prior `decay` and `sd` of the latent
#'   process.
#' @param stepper Either `"explicit"` or `"implicit"`.
#' @param nsub Positive whole number of sub-steps taken per observation.
#' @param seed Non-negative whole number seeding the sampler.
#' @param ... Passed to the partial-observation learner.
#' @returns An object of class `aci_fit`.
#'
#' @seealso [cir()], [learn_model()], [aci()]
#'
#' @examples
#' \donttest{
#' sim <- simulate(model_dyad(), seed = 8, T = 5, dt = 5e-3, burn_in = 1)
#' dat <- data.frame(t = sim$obs$t, x1 = sim$obs$x[, 1], y = sim$hidden[, 1])
#' aci_fit(y ~ x1, data = dat)
#' }
#'
#' @export
aci_fit <- function(formula, data, time = NULL, group = NULL, degree = 2L,
                    ce_threshold = "auto", energy_pairs = NULL,
                    n_samples = 20, hidden_prior = list(decay = 0.7, sd = 1),
                    stepper = c("explicit", "implicit"), nsub = 1L,
                    seed = 1, ...) {
  stepper <- match.arg(stepper)
  cl <- match.call()
  if (!inherits(formula, "formula") || length(formula) != 3L)
    aci_abort("aci_error_model_contract",
      "aci_fit() needs a two-sided formula `latent ~ effect_1 + effect_2`.")
  latent <- all.vars(formula[[2L]]); observed <- all.vars(formula[[3L]])
  if (length(latent) != 1L)
    aci_abort("aci_error_model_contract",
      "aci_fit() needs exactly one latent cause on the left of `~` in v0.")
  if (!length(observed) || latent %in% observed)
    aci_abort("aci_error_model_contract",
      "The right side of the formula must contain at least one observed effect and must not repeat the latent variable.")
  if (!is.null(time) && (length(time) != 1L || !is.character(time) ||
      is.na(time) || !nzchar(time)))
    aci_abort("aci_error_obs_contract", "`time` must be NULL or one non-empty column name.")
  if (!is.null(group) && (length(group) != 1L || !is.character(group) ||
      is.na(group) || !nzchar(group)))
    aci_abort("aci_error_obs_contract", "`group` must be NULL or one non-empty column name.")
  segs <- .fit_prepare(data, observed, latent, time, group)
  has_latent <- vapply(segs, function(d) latent %in% names(d), logical(1))
  if (any(has_latent) && !all(has_latent))
    aci_abort("aci_error_obs_contract",
              "The latent column must be present in every segment or absent from every segment.")
  mode <- if (all(has_latent)) "full" else "partial"
  k <- length(observed); l <- 1L
  lib <- cgns_library(k = k, l = l, degree = degree)
  dt <- segs[[1]]$.dt[1]

  if (mode == "full") {
    fit <- .fit_full(segs, observed, latent, lib, ce_threshold, energy_pairs)
  } else {
    fit <- .fit_partial(segs, observed, latent, lib, ce_threshold,
                        energy_pairs, n_samples, hidden_prior, seed, ...)
  }
  model <- model_from_learned(fit$coefs_x, fit$coefs_y, lib,
                              template = fit$template)
  model$meta$ic_default <- list(x0 = as.numeric(segs[[1]][1, observed]),
                                y0 = fit$y0 %||% 0)
  # representative trajectory: the longest segment
  seg <- segs[[which.max(vapply(segs, nrow, integer(1)))]]
  obs <- observed_trajectory(seg$.t, as.matrix(seg[, observed, drop = FALSE]),
                             names = observed)
  init <- list(mean = rep(0, l), cov = diag(1, l))
  filt <- suppressWarnings(da_filter(model, obs, init = init,
                                     stepper = stepper, nsub = nsub))
  smoo <- da_smooth(model, obs, filter = filt)
  klp <- gaussian_kl_path(smoo, filt)
  structure(list(model = model, formula = formula, call = cl,
                 observed = observed, latent = latent, mode = mode,
                 obs = obs, dt = dt, init = init,
                 paths = list(filter = filt, smoother = smoo),
                 aci = data.frame(t = klp$t, aci = klp$total,
                                  signal = klp$signal,
                                  dispersion = klp$dispersion),
                 report = fit$report, n_segments = length(segs),
                 n_observations = sum(vapply(segs, nrow, integer(1))),
                 cache = new.env(parent = emptyenv())),
            class = "aci_fit")
}


################################################################################
#  data marshalling
################################################################################

#' Split fitting data into uniform segments (internal)
#'
#' @param data A data frame or list of data frames.
#' @param observed Character vector of observed variable names.
#' @param latent Character vector of latent variable names.
#' @param time Optional name of the time column.
#' @param group Optional name of the grouping column.
#' @returns A list of segments, each carrying its observations and step.
#' @noRd
.fit_prepare <- function(data, observed, latent, time, group) {
  if (is.data.frame(data)) {
    if (!is.null(group) && !group %in% names(data))
      aci_abort("aci_error_obs_contract",
                sprintf("Grouping column `%s` is not present in `data`.", group))
    if (!is.null(group) && anyNA(data[[group]]))
      aci_abort("aci_error_obs_contract", "The grouping column must not contain missing values.")
    data <- if (!is.null(group)) split(data, data[[group]], drop = TRUE) else list(data)
  }
  if (!is.list(data) || !length(data))
    aci_abort("aci_error_obs_contract", "`data` must be a data.frame or a list of them.")
  out <- lapply(data, function(d) {
    d <- as.data.frame(d)
    tc <- time %||% intersect(c("t", "time", "date"), names(d))[1]
    if (is.na(tc) || is.null(tc) || !tc %in% names(d))
      aci_abort("aci_error_obs_contract",
        "No time column found; supply `time` or include one of t/time/date.")
    tv <- d[[tc]]
    if (inherits(tv, "Date") || inherits(tv, "POSIXct"))
      tv <- as.numeric(tv) - as.numeric(tv[1])
    tv <- as.numeric(tv)
    miss <- setdiff(observed, names(d))
    if (length(miss))
      aci_abort("aci_error_obs_contract",
        sprintf("Observed effect(s) not in data: %s.", paste(miss, collapse = ", ")))
    if (length(tv) < 3L || any(!is.finite(tv)))
      aci_abort("aci_error_obs_contract",
                "Each segment needs at least three finite time points.")
    needed <- intersect(c(observed, latent), names(d))
    if (any(!vapply(d[needed], is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(d[needed]))))
      aci_abort("aci_error_obs_contract", "Model variables must be finite numeric columns.")
    o <- order(tv); d <- d[o, , drop = FALSE]; tv <- tv[o]
    dtv <- diff(tv)
    if (any(dtv <= 0))
      aci_abort("aci_error_obs_contract", "Times within a segment must be unique and increasing.")
    if (max(abs(dtv - dtv[1])) > 1e-8 * max(abs(dtv[1]), 1e-12))
      aci_abort("aci_error_obs_contract",
                "aci_fit() requires a uniform time grid within every segment.")
    d$.t <- tv - tv[1]; d$.dt <- dtv[1]
    d
  })
  dts <- vapply(out, function(d) d$.dt[1], numeric(1))
  if (max(abs(dts - dts[1])) > 1e-8 * max(abs(dts[1]), 1e-12))
    aci_abort("aci_error_obs_contract",
              "All segments must use the same time step.")
  out
}


################################################################################
#  full-state learning
################################################################################

#' Learn a model from fully observed segments (internal)
#'
#' @param segs List of prepared segments.
#' @param observed Character vector of observed variable names.
#' @param latent Character vector of latent variable names.
#' @param lib A function library; see [cgns_library()].
#' @param ce_threshold Either `"auto"` or a numeric causation entropy threshold.
#' @param energy_pairs Optional list of energy-conserving index pairs.
#' @returns A list with the learned equations and the fitted model.
#' @noRd
.fit_full <- function(segs, observed, latent, lib, ce_threshold, energy_pairs) {
  k <- length(observed)
  Th <- dX <- dY <- NULL
  for (d in segs) {
    Z <- as.matrix(d[, c(observed, latent), drop = FALSE])
    n <- nrow(Z); dt <- d$.dt[1]
    Th <- rbind(Th, eval_library(lib, Z[-n, , drop = FALSE], d$.t[-n]))
    dX <- rbind(dX, (Z[-1, seq_len(k), drop = FALSE] -
                     Z[-n, seq_len(k), drop = FALSE]) / dt)
    dY <- c(dY, diff(Z[, k + 1]) / dt)
  }
  fit_eq <- function(dz, restrict) {
    cand <- if (restrict) lib$names[lib$y_degree <= 1] else lib$names
    ce <- causation_entropy(Th[, cand, drop = FALSE], dz, cand)
    kept <- threshold_structure(ce, ce_threshold)
    if ("1" %in% lib$names) kept <- union("1", kept)
    f <- constrained_mle(Th[, kept, drop = FALSE], dz)
    list(kept = kept, ce = ce, coef = f$coef)
  }
  eqs <- c(lapply(seq_len(k), function(i) fit_eq(dX[, i], TRUE)),
           list(fit_eq(dY, TRUE)))
  names(eqs) <- c(paste0("d", observed), paste0("d", latent))
  if (!is.null(energy_pairs))
    eqs <- .refit_energy_joint(eqs, Th, c(lapply(seq_len(k), function(i) dX[, i]),
                                          list(dY)), lib, energy_pairs, k, 1L,
                                 dt = segs[[1]]$.dt[1])
  sds <- .fit_residual_sds(eqs, Th, dX, dY, segs[[1]]$.dt[1], k)
  list(coefs_x = lapply(seq_len(k), function(i) eqs[[i]]$coef),
       coefs_y = list(eqs[[k + 1]]$coef),
       template = .fit_template(k, sds),
       y0 = segs[[1]][[latent]][1],
       report = list(mode = "full", equations = eqs, sigma = sds))
}


#' Residual standard deviations of learned equations (internal)
#'
#' @param eqs The learned equation coefficients.
#' @param Th Evaluated library matrix.
#' @param dX Observed increments.
#' @param dY Latent increments.
#' @param dt Positive 1-length numeric step.
#' @param k Observed dimension.
#' @returns Numeric vector of residual standard deviations.
#' @noRd
.fit_residual_sds <- function(eqs, Th, dX, dY, dt, k) {
  one <- function(e, dz) {
    r <- dz - drop(Th[, e$kept, drop = FALSE] %*% e$coef)
    stats::sd(r) * sqrt(dt)                # increment-noise amplitude
  }
  c(vapply(seq_len(k), function(i) one(eqs[[i]], dX[, i]), numeric(1)),
    one(eqs[[k + 1]], dY))
}


#' Model template for a learned system (internal)
#'
#' @param k Observed dimension.
#' @param sds Numeric vector of residual standard deviations.
#' @returns A template list consumed by [model_from_learned()].
#' @noRd
.fit_template <- function(k, sds) {
  list(k = k, l = 1L,
       Sx1 = function(t, x) diag(sds[seq_len(k)], k),
       Sy2 = function(t, x) matrix(sds[k + 1], 1, 1))
}


################################################################################
#  partial-observation learning (experimental)
################################################################################

#' Learn a model with a latent variable present (internal)
#'
#' @param segs List of prepared segments.
#' @param observed Character vector of observed variable names.
#' @param latent Character vector of latent variable names.
#' @param lib A function library; see [cgns_library()].
#' @param ce_threshold Either `"auto"` or a numeric causation entropy threshold.
#' @param energy_pairs Optional list of energy-conserving index pairs.
#' @param n_samples Number of posterior path samples drawn.
#' @param hidden_prior List with the prior `decay` and `sd` of the latent
#'   process.
#' @param seed Non-negative whole number seeding the sampler.
#' @param ... Passed to the underlying learner.
#' @returns A list with the learned equations and the fitted model.
#' @noRd
.fit_partial <- function(segs, observed, latent, lib, ce_threshold,
                         energy_pairs, n_samples, hidden_prior, seed, ...) {
  k <- length(observed)
  if (!is.list(hidden_prior))
    aci_abort("aci_error_model_contract", "hidden_prior must be a list.")
  decay <- hidden_prior$decay %||% 0.7
  hidden_sd <- hidden_prior$sd %||% 1
  if (!is.numeric(decay) || length(decay) != 1L || !is.finite(decay) ||
      decay <= 0 || !is.numeric(hidden_sd) || length(hidden_sd) != 1L ||
      !is.finite(hidden_sd) || hidden_sd <= 0)
    aci_abort("aci_error_model_contract",
              "hidden_prior$decay and hidden_prior$sd must be finite positive scalars.")
  if (length(segs) != 1L)
    aci_abort("aci_error_not_implemented",
              "Partial-observation fitting of multiple segments is not implemented; fit segments separately.")
  d <- segs[[which.max(vapply(segs, nrow, integer(1)))]]
  obs <- observed_trajectory(d$.t, as.matrix(d[, observed, drop = FALSE]),
                             names = observed)
  # OU seed for the hidden channel + coupling-free observed fit, then one
  # posterior-sampling learn_model pass around that seed.
  Xo <- as.matrix(d[, observed, drop = FALSE]); n <- nrow(Xo); dt <- d$.dt[1]
  cand0 <- lib$names[lib$y_degree == 0]
  Th0 <- eval_library(lib, cbind(Xo, 0)[-n, , drop = FALSE], d$.t[-n])
  seed_cx <- lapply(seq_len(k), function(i) {
    dz <- diff(Xo[, i]) / dt
    ce <- causation_entropy(Th0[, cand0, drop = FALSE], dz, cand0)
    kept <- threshold_structure(ce, ce_threshold)
    if ("1" %in% lib$names) kept <- union("1", kept)
    cf <- constrained_mle(Th0[, kept, drop = FALSE], dz)$coef
    # seed a small latent-linear coupling so the posterior is informative
    cf2 <- c(cf, stats::setNames(0.1, paste0("y1")))
    cf2
  })
  seed_cy <- list(stats::setNames(-decay, "y1"))
  sdx <- vapply(seq_len(k), function(i)
    stats::sd(diff(Xo[, i])) / sqrt(dt), numeric(1))     # amplitude = inc-sd / sqrt(dt)
  sds <- c(sdx, hidden_sd * sqrt(2 * decay))
  seed_model <- model_from_learned(seed_cx, seed_cy, lib,
                                   template = .fit_template(k, sds))
  lm_ <- suppressWarnings(learn_model(seed_model, obs, lib,
                                      n_samples = n_samples,
                                      ce_threshold = ce_threshold,
                                      energy_pairs = energy_pairs,
                                      enforce_cgns = TRUE, seed = seed, ...))
  eqs <- lm_$equations
  list(coefs_x = lapply(seq_len(k), function(i) eqs[[i]]$coef),
       coefs_y = list(eqs[[k + 1]]$coef),
       template = .fit_template(k, sds),
       y0 = 0,
       report = list(mode = "partial (experimental)", equations = eqs))
}


################################################################################
#  standard generics
################################################################################

#' Print a fitted latent-cause model
#'
#' @param x An `aci_fit` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.aci_fit <- function(x, ...) {
  cat(sprintf("<aci_fit> %s | mode = %s | %d obs steps%s\n",
              deparse1(x$formula), x$mode, length(x$obs$t),
              if (x$n_segments > 1) sprintf(" (%d segments)", x$n_segments) else ""))
  pk <- which.max(x$aci$aci)
  cat(sprintf("  peak ACI %.4g at t = %.4g; see summary(), cir(), plot().\n",
              x$aci$aci[pk], x$aci$t[pk]))
  invisible(x)
}


#' Summarize a fitted latent-cause model
#'
#' @param object An `aci_fit` object.
#' @param verdict `TRUE` to run the surrogate-null test and report its verdict.
#' @param B Number of surrogate replicates used when `verdict` is `TRUE`.
#' @param ... Ignored, for consistency with the generic.
#' @returns An object of class `summary.aci_fit`.
#' @export
summary.aci_fit <- function(object, verdict = FALSE, B = 49, ...) {
  s <- list(fit = object,
            coefficients = lapply(object$report$equations, function(e)
              round(e$coef, 4)),
            peak_aci = max(object$aci$aci),
            mean_aci = mean(object$aci$aci[-seq_len(ceiling(0.1 * nrow(object$aci)))]),
            cir = tryCatch(cir(object), error = function(e) NULL),
            verdict = if (isTRUE(verdict))
              nil_surrogate_test(object$model, object$obs, B = B,
                                 init = object$init) else NULL)
  class(s) <- "summary.aci_fit"
  s
}


#' Print a fitted-model summary
#'
#' @param x A `summary.aci_fit` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.summary.aci_fit <- function(x, ...) {
  print(x$fit)
  cat("\nLearned equations (per-equation kept terms):\n")
  for (nm in names(x$coefficients)) {
    cf <- x$coefficients[[nm]]
    cat(sprintf("  %s/dt = %s\n", nm,
        paste(sprintf("%+.3g %s", cf, names(cf)), collapse = " ")))
  }
  cat(sprintf("\nACI: peak %.4g, post-burn mean %.4g\n", x$peak_aci, x$mean_aci))
  if (!is.null(x$cir)) {
    ok <- is.finite(x$cir$tau)
    if (any(ok)) cat(sprintf("Forward CIR (median over anchors): %.4g time units\n",
                             stats::median(x$cir$tau[ok])))
  }
  if (!is.null(x$verdict)) print(x$verdict)
  if (grepl("partial", x$fit$mode))
    cat("NOTE: partial-observation mode is experimental (see ?learn_model).\n")
  invisible(x)
}


#' Coefficients of a fitted latent-cause model
#'
#' @param object An `aci_fit` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns The learned equation coefficients.
#' @export
coef.aci_fit <- function(object, ...)
  lapply(object$report$equations, `[[`, "coef")


#' Number of observations behind a fitted latent-cause model
#'
#' @param object An `aci_fit` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns 1-length integer, the number of observations used in fitting.
#' @export
nobs.aci_fit <- function(object, ...)
  object$n_observations %||% length(object$obs$t)


#' Reconstructed latent path of a fitted model
#'
#' @param object An `aci_fit` object.
#' @param type Either `"smoother"` or `"filter"`.
#' @param ... Ignored, for consistency with the generic.
#' @returns The requested assimilation path.
#' @export
predict.aci_fit <- function(object, type = c("smoother", "filter"), ...) {
  type <- match.arg(type)
  p <- object$paths[[type]]
  data.frame(t = p$t, mean = p$mean[, 1],
             sd = sqrt(pmax(p$cov[1, 1, ], 0)))
}


#' Fitted latent mean of a fitted model
#'
#' @param object An `aci_fit` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns Numeric matrix of the smoothed latent mean.
#' @export
fitted.aci_fit <- function(object, ...) predict(object)$mean


#' Residuals of a fitted latent-cause model
#'
#' @param object An `aci_fit` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns Numeric matrix of observed-increment residuals.
#' @export
residuals.aci_fit <- function(object, ...) {
  # observed equation increment residuals at the smoother-mean latent path
  ob <- object$obs; sm <- object$paths$smoother; mdl <- object$model
  dt <- ob$dt; N <- length(ob$t) - 1
  r <- matrix(NA_real_, N, mdl$k)
  for (j in seq_len(N)) {
    co <- eval_coefs(mdl, ob$t[j], ob$x[j, ])
    r[j, ] <- (ob$x[j + 1, ] - ob$x[j, ]) / dt -
              (co$fx + drop(co$Lx %*% sm$mean[j, ]))
  }
  colnames(r) <- object$observed
  r
}


#' Simulate from a fitted latent-cause model
#'
#' @param object An `aci_fit` object.
#' @param nsim Positive whole number of realisations to generate.
#' @param seed Optional non-negative whole number seeding the generator.
#' @param T Optional positive total simulated time; `NULL` reuses the fitted
#'   record's span.
#' @param dt Optional positive integration step; `NULL` reuses the fitted step.
#' @param ... Passed to [simulate.stochastic_model()].
#' @returns An `aci_sim` object, or a list of them when `nsim` exceeds one.
#' @export
simulate.aci_fit <- function(object, nsim = 1, seed = NULL, T = NULL,
                             dt = NULL, ...) {
  simulate(object$model, nsim = nsim, seed = seed,
           T = T %||% (object$obs$t[length(object$obs$t)]),
           dt = dt %||% object$dt, ...)
}


#' Plot a fitted latent-cause model
#'
#' @param x An `aci_fit` object.
#' @param which One of `"latent"`, `"metric"` or `"cir"`.
#' @param truth Optional numeric vector of true latent values to overlay on
#'   the `"latent"` plot; it is included in the vertical range so neither the
#'   truth nor the 2 sd band is clipped.
#' @param ... Passed to the underlying base-graphics calls.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.aci_fit <- function(x, which = c("latent", "metric", "cir"),
                         truth = NULL, ...) {
  which <- match.arg(which)
  if (which == "cir") {
    fc <- cir(x, full = TRUE)
    plot(fc$t, fc$tau, type = "l", xlab = "t",
         ylab = expression(tau[f]),
         main = sprintf("Forward CIR: how long %s keeps driving {%s}",
                        x$latent, paste(x$observed, collapse = ", ")), ...)
    return(invisible(x))
  }
  if (which == "latent") {
    p <- predict(x)
    ylim <- range(p$mean + 2 * p$sd, p$mean - 2 * p$sd, truth, na.rm = TRUE)
    plot(p$t, p$mean, type = "l", xlab = "t", ylim = ylim,
         ylab = x$latent, main = sprintf("Reconstructed %s (smoother, +/- 2 sd)",
                                         x$latent), ...)
    graphics::lines(p$t, p$mean + 2 * p$sd, lty = 3, col = "grey40")
    graphics::lines(p$t, p$mean - 2 * p$sd, lty = 3, col = "grey40")
    if (!is.null(truth)) graphics::lines(p$t, truth, col = 2)
  } else {
    plot(x$aci$t, x$aci$aci, type = "l", xlab = "t", ylab = "ACI",
         main = sprintf("ACI(t): %s -> {%s}", x$latent,
                        paste(x$observed, collapse = ", ")), ...)
  }
  invisible(x)
}


################################################################################
# the estimand extractor
################################################################################

#' Causal influence range of a fitted object
#'
#' `cir()` extracts a forward or backward causal influence range from a fitted
#' object.
#'
#' @param object A fitted object, currently an `aci_fit`.
#' @param ... Arguments passed to methods.
#' @returns An object of class `cir_result`, or a data frame when the full table
#'   is requested.
#'
#' @seealso [forward_cir()], [backward_cir()], [aci_fit()]
#'
#' @examples
#' \donttest{
#' sim <- simulate(model_dyad(), seed = 8, T = 5, dt = 5e-3, burn_in = 1)
#' dat <- data.frame(t = sim$obs$t, x1 = sim$obs$x[, 1], y = sim$hidden[, 1])
#' fit <- aci_fit(y ~ x1, data = dat)
#' cir(fit)
#' }
#'
#' @export
cir <- function(object, ...) UseMethod("cir")


#' @describeIn cir Errors for objects with no causal influence range method.
#' @export
cir.default <- function(object, ...)
  aci_abort("aci_error_model_contract",
            "cir() expects an aci_fit object (or use forward_cir()/backward_cir() directly).")


#' @describeIn cir Causal influence range of a fitted latent-cause model.
#' @param direction Either `"forward"` or `"backward"`.
#' @param at Optional numeric anchor times; `NULL` uses the fitted grid.
#' @param tol Positive tolerance for the table's adaptive freeze rule.
#' @param full `TRUE` to return the whole tidy table rather than the range
#'   object.
#' @export
cir.aci_fit <- function(object, direction = c("forward", "backward"),
                        at = NULL, tol = 1e-3, full = FALSE, ...) {
  direction <- match.arg(direction)
  if (direction == "forward") {
    key <- paste0("fwd_table:", format(tol, digits = 17))
    tab <- object$cache[[key]]
    if (is.null(tab)) {
      tab <- lag_table(object$model, object$obs, mode = "forward", tol = tol,
                       init = object$init)
      object$cache[[key]] <- tab
    }
    fc <- forward_cir(tab, ...)
    if (isTRUE(full))
      return(data.frame(t = fc$t, tau = fc$tau, strength = fc$M,
                        direction = "forward"))
    idx <- if (is.null(at)) {
      qs <- stats::quantile(seq_along(fc$t), c(.25, .5, .75))
      unique(round(qs))
    } else vapply(at, function(a) which.min(abs(fc$t - a)), integer(1))
    data.frame(t = fc$t[idx], tau = fc$tau[idx], strength = fc$M[idx],
               direction = "forward")
  } else {
    Ts <- if (is.null(at)) {
      tt <- object$obs$t
      stats::quantile(tt, c(.5, .75, .95))
    } else at
    b <- suppressWarnings(backward_cir(object$model, obs = object$obs,
                                       T = Ts, init = object$init, ...))
    data.frame(t = b$t, tau = b$tau, strength = b$M, direction = "backward")
  }
}


################################################################################
# base-graphics plot methods for the core classes
################################################################################

#' Plot an ACI result
#'
#' @param x An `aci_result` object.
#' @param decompose `TRUE` to draw the signal and dispersion parts alongside the
#'   total.
#' @param ... Passed to the underlying base-graphics calls.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.aci_result <- function(x, decompose = TRUE, ...) {
  plot(x$t, x$aci, type = "l", xlab = "t", ylab = "ACI",
       main = "ACI(t) = KL(smoother || filter)", ...)
  if (decompose && !is.null(x$signal)) {
    graphics::lines(x$t, x$signal, col = 4, lty = 2)
    graphics::lines(x$t, x$dispersion, col = 2, lty = 3)
    graphics::legend("topright", c("total", "signal", "dispersion"),
                     col = c(1, 4, 2), lty = 1:3, bty = "n")
  }
  invisible(x)
}


#' Plot a Gaussian assimilation path
#'
#' @param x A `da_path_gaussian` object.
#' @param component Integer index of the hidden component to draw.
#' @param truth Optional numeric vector of true hidden values to overlay.
#' @param ... Passed to the underlying base-graphics calls.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.da_path_gaussian <- function(x, component = 1, truth = NULL, ...) {
  mu <- x$mean[, component]; sdv <- sqrt(pmax(x$cov[component, component, ], 0))
  ylim <- range(mu + 2 * sdv, mu - 2 * sdv, truth, na.rm = TRUE)
  plot(x$t, mu, type = "l", xlab = "t", ylab = sprintf("y%d", component),
       ylim = ylim, main = sprintf("%s mean +/- 2 sd", x$kind), ...)
  graphics::lines(x$t, mu + 2 * sdv, lty = 3, col = "grey40")
  graphics::lines(x$t, mu - 2 * sdv, lty = 3, col = "grey40")
  if (!is.null(truth)) graphics::lines(x$t, truth, col = 2)
  invisible(x)
}


#' Plot a causal influence range result
#'
#' @param x A `cir_result` object.
#' @param ... Passed to the underlying base-graphics calls.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.cir_result <- function(x, ...) {
  if (length(x$tau) > 1) {
    plot(x$t, x$tau, type = "l", xlab = "t",
         ylab = sprintf("tau_%s", substr(x$direction, 1, 1)),
         main = sprintf("%s CIR (%s)", x$direction, x$method), ...)
  } else {
    graphics::plot.new(); graphics::title(main = sprintf(
      "%s CIR at T = %.4g: tau = %.4g", x$direction, x$t, x$tau))
  }
  invisible(x)
}
