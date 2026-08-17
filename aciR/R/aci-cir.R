# Causal influence range -------------------------------------------------------
#
# The second quantity of the method paper. Where the causal-information metric
# asks HOW MUCH the future of the observed signal says about the unobserved
# state at time t, the causal influence range asks HOW FAR INTO THAT FUTURE one
# must look before the answer stops changing.
#
# The construction is a divergence between two online-smoother posteriors at
# the same instant. One is informed by the whole record, the other only up to
# some later observation. As that later observation moves forward the two
# coincide, and the range is where the gap falls below a threshold.
#
# The reference implementation stores the whole matrix of those divergences
# (one row per time, one column per lagged observation), which is quadratic in
# the window and reaches several gigabytes at the scale of the published
# figure. Nothing downstream needs the matrix, because each row is reduced
# immediately to a handful of scalars. This implementation therefore forms one
# row at a time and discards it, which leaves the memory linear in the window
# while the arithmetic stays what it has to be.

#' Causal influence range
#'
#' Computes the subjective and objective causal influence range of the
#' unobserved component on the observed signal. The causal influence range
#' measures how far forward in the observed record one must look before the
#' estimate of the unobserved state at a given time stops improving, and so
#' complements the causal-information metric, which measures how much it
#' improves in total.
#'
#' At each reported time the function forms the relative entropy between the
#' online-smoother posterior informed by the whole record and the posterior
#' informed only up to each later observation. That sequence decays as the
#' later observation advances. The **subjective** range at a threshold is the
#' elapsed time after which the sequence stays below that threshold; the
#' **objective** range is the threshold-free summary obtained by integrating
#' the sequence and normalising by its peak, which the source paper gives as a
#' computationally efficient underestimate of the range defined by averaging
#' the subjective ranges over all thresholds.
#'
#' The computation is quadratic in the length of `window`, because every
#' reported time is compared against every later observation. Choose the window
#' accordingly. A few thousand steps is comfortable, and reproducing a figure
#' at the scale of the published one is a batch computation rather than an
#' interactive one.
#'
#' A reported time close to the end of the record cannot be resolved, because
#' the observations that would settle it do not exist. Such a time is not
#' unmeasured. Its range is **right-censored**, and the truncated value is a
#' lower bound on the true one. The result therefore returns the value and
#' marks that time `"censored"` in `status`, rather than discarding what the
#' record does support. Only a time with fewer than three later observations,
#' where no quadrature is possible at all, returns `NA`.
#'
#' `objective` and `objective_exact` are different functionals, not two
#' quadratures of one. They coincide when the divergence decreases with lag, by
#' the layer-cake identity. Both this package and the reference measure the
#' range as the **last** time the divergence exceeds a threshold, which is at
#' least the measure of the superlevel set and is strictly larger the moment
#' the sequence is not monotone. A sequence that is not monotone is the common
#' case on a truncated `horizon`. Seeing `objective` below `objective_exact` is
#' that definitional gap, not a numerical defect in either.
#'
#' @param x Numeric vector. The observed signal, one value per time step.
#' @param comp A conditional Gaussian components list; see [aci_components].
#' @param dt Numeric scalar. The integration time step; must be positive.
#' @param filt A list with numeric vectors `mean` and `cov`, as returned by
#'   [aci_filter()]. When `NULL`, the filter is run internally, which requires
#'   `mu0` and `R0`.
#' @param window Integer vector or `NULL`. Indices of the time steps at which
#'   the range is reported. When `NULL`, the whole signal is used, which is
#'   quadratic in its length and is rarely what is wanted for a long record.
#' @param epsilon Numeric vector. Thresholds at which the subjective range is
#'   evaluated, in nats. Defaults to a 129-point logarithmic grid from `1e-6`
#'   to `10^0.5`. The reference implementation spans the same range with 513
#'   points; the cheaper grid is this package's default, and
#'   `10^seq(-6, 0.5, length.out = 513L)` reproduces the reference's
#'   `objective_exact`.
#' @param threshold Numeric scalar. A peak divergence below this value is
#'   treated as no detectable influence, and the objective range is reported as
#'   zero rather than as the ratio of two negligible quantities. Defaults to
#'   `1e-5`.
#' @param margin Numeric scalar in `(0, 1)`. The fraction of each comparison
#'   sequence that must remain unused for the range at that time to count as
#'   resolved. A time whose range consumes more than `1 - margin` of the
#'   sequence it was measured against is marked `status = "censored"`; its
#'   ranges are still returned, as lower bounds. Defaults to `0.1`.
#'
#'   Because a censored time yields a bound rather than a hole, `margin`
#'   governs the flag rather than whether a number is reported at all, so a
#'   slightly wrong value has little consequence. It is this package's device.
#'   The reference guards the same problem with an absolute lookahead chosen
#'   for one figure.
#' @param horizon Integer scalar or `NULL`. How many steps of the record each
#'   reported time may look forward across, counted from the start of the
#'   record rather than from the reported time. `NULL`, the default, uses the
#'   whole record.
#'
#'   The reference implementation truncates this comparison at the end of its
#'   reporting window, which biases both the integral and the range low near
#'   that end; supplying the same value here reproduces its numbers. The
#'   truncation applies only to how far forward the comparison looks, never to
#'   the fully informed posterior it is compared against, which is always taken
#'   over the whole record.
#' @param mu0,R0 Numeric scalars. Initial filtered mean and covariance, used
#'   only when `filt` is `NULL`.
#'
#' @returns An object of class `aci_cir`, a list with the reported `time`, the
#'   `objective` range at each time, the `objective_exact` range obtained by
#'   integrating the subjective ranges over the whole threshold grid rather
#'   than by the efficient approximation, the `subjective` range as a matrix
#'   with one row per threshold and one column per time, the `peak` divergence
#'   at each time, the logical matrix `subjective_censored` marking thresholds
#'   whose range ran past the retained margin, the character `status`
#'   (`"resolved"`, `"censored"`, `"below_threshold"` or `"insufficient"`), the
#'   logical `monotone` marking times whose divergence sequence decreases with
#'   lag (the condition under which `objective` and `objective_exact` are the
#'   same functional), and the logical `saturated`, which is
#'   `status == "censored"`.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*,
#' 36(4), 71. \doi{10.1007/s00332-026-10271-x}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 600, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:250,
#'                mu0 = model$y0, R0 = 0.1)
#' summary(rng$objective)
#'
#' # `horizon` is an index into the record, not a lookahead from each reported
#' # time, so reproducing a published panel means converting from a time. This
#' # is the conversion the reference performs to obtain its `last_idx`.
#' dt <- 0.001
#' time_end <- 0.25
#' lookahead <- 0.05
#' last_idx <- as.integer(round((time_end + lookahead) / dt)) + 1L
#' aci_cir(sim$x, comp, dt = dt, window = 50:250, mu0 = model$y0, R0 = 0.1,
#'         horizon = last_idx)
#'
#' @seealso [aci_online_smoother()], [aci_metric()]
#' @export
aci_cir <- function(x, comp, dt, filt = NULL, window = NULL,
                    epsilon = 10^seq(-6, 0.5, length.out = 129L),
                    threshold = 1e-5, margin = 0.1, horizon = NULL,
                    mu0 = NULL, R0 = NULL) {
  .aci_check_positive(dt, "dt")
  .aci_check_positive(threshold, "threshold")
  is_mv <- .aci_is_mv(comp)
  if (is_mv) {
    x <- .aci_check_signal_mv(x)
    n <- ncol(x)
    if (is.null(filt)) {
      stop(
        paste0(
          "`filt` must be supplied for a vector system: the initial mean and ",
          "covariance have no scalar defaults there."
        ),
        call. = FALSE
      )
    }
    # Checked here rather than left to fail inside the row recursion. A
    # malformed vector posterior used to surface as a subsetting error from
    # three frames down, which tells a caller nothing about which argument was
    # wrong.
    if (!is.list(filt) || !is.matrix(filt$mean) ||
          length(dim(filt$cov)) != 3L) {
      stop(
        paste0(
          "`filt` must be a vector-valued posterior, as returned by ",
          "`aci_filter()` on a vector system: `mean` a matrix and `cov` a ",
          "three-dimensional array."
        ),
        call. = FALSE
      )
    }
    if (ncol(filt$mean) != n || dim(filt$cov)[3L] != n) {
      stop(
        sprintf(
          paste0(
            "`filt` covers %d step(s) but the observed signal has %d. The ",
            "posterior and the signal must be the same length."
          ),
          ncol(filt$mean), n
        ),
        call. = FALSE
      )
    }
    comp <- .aci_check_components_mv(comp, x)
  } else {
    .aci_check_signal(x)
    n <- length(x)
    .aci_check_components(comp, n)
  }
  .aci_check_scalar(margin, "margin")
  if (margin <= 0 || margin >= 1) {
    stop("`margin` must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.null(horizon)) {
    .aci_check_scalar(horizon, "horizon")
    if (horizon != trunc(horizon) || horizon < 1) {
      stop("`horizon` must be a whole number of at least one, or `NULL`.",
           call. = FALSE)
    }
    horizon <- min(as.integer(horizon), n)
  } else {
    horizon <- n
  }
  if (is.null(filt)) {
    .aci_check_scalar(mu0, "mu0")
    .aci_check_positive(R0, "R0")
    filt <- aci_filter(x, comp, dt, mu0 = mu0, R0 = R0)
  }
  if (!is_mv) {
    .aci_check_posterior(filt, "filt", n)
  }
  window <- .aci_check_window(window, n)
  if (!is.numeric(epsilon) || length(epsilon) < 1L || anyNA(epsilon) ||
        any(epsilon <= 0)) {
    stop(
      "`epsilon` must be a vector of positive, non-missing thresholds.",
      call. = FALSE
    )
  }

  aux <- if (is_mv) {
    .aci_online_aux_mv(x, comp, dt, filt)
  } else {
    .aci_online_aux(x, comp, dt, filt)
  }
  eps_sorted <- sort(epsilon)
  n_eps <- length(epsilon)

  objective <- numeric(length(window))
  objective_exact <- rep(NA_real_, length(window))
  peak <- numeric(length(window))
  status <- rep("resolved", length(window))
  monotone <- rep(NA, length(window))
  subjective <- matrix(NA_real_, nrow = n_eps, ncol = length(window))
  subjective_censored <- matrix(FALSE, nrow = n_eps, ncol = length(window))

  for (i in seq_along(window)) {
    j <- window[i]
    re <- if (is_mv) {
      .aci_cir_row_mv(aux, filt, j, n, tol = 1e-18, horizon = horizon)
    } else {
      .aci_cir_row(aux, filt, j, n, horizon = horizon)
    }
    if (is.null(re)) {
      # The final step has no later observation to be compared against.
      status[i] <- "insufficient"
      objective[i] <- NA_real_
      peak[i] <- 0
      next
    }
    if (length(re) < 3L) {
      # Too few later observations to support a range at all. This arises when
      # `horizon` stops the comparison close to the reported time.
      #
      # The reference wraps its quadrature in a try/catch and records zero
      # here, which a reader cannot distinguish from "no detectable
      # influence". The two mean opposite things. This is the one status that
      # genuinely has no number behind it, so it alone returns NA.
      status[i] <- "insufficient"
      objective[i] <- NA_real_
      peak[i] <- max(re)
      next
    }
    peak[i] <- max(re)
    # Whether the divergence decreases with lag decides whether `objective` and
    # `objective_exact` are two quadratures of one functional or two different
    # functionals. The layer-cake identity that makes the efficient integral an
    # underestimate of the threshold average needs a decreasing sequence; the
    # range itself is measured as a LAST exit, which exceeds the measure of the
    # superlevel set the moment the sequence is not monotone.
    monotone[i] <- !is.unsorted(rev(re))

    # ---- Objective range ----------------------------------------------------
    if (peak[i] > threshold) {
      objective[i] <- .aci_simpson(re) * dt / peak[i]
    }

    # ---- Subjective range at every threshold at once ------------------------
    #
    # The range at a threshold is the last position at which the divergence
    # still exceeds it. The suffix maximum of the sequence is non-increasing,
    # and the count of its entries above a threshold is exactly that position,
    # so one sorted search answers the whole grid.
    suffix <- rev(cummax(rev(re)))
    counts <- length(suffix) - findInterval(eps_sorted, rev(suffix))

    # A range that consumes most of the record it was measured against is not
    # a range; it is a statement that the record ended first. The divergence
    # is exactly zero at the final observation, by construction, so a range can
    # never formally reach the end of the row, and testing against the end
    # alone would therefore never fire. The test is against the retained
    # margin.

    # ---- Objective range by its definition ---------------------------------
    #
    # The source paper defines the objective range as the subjective ranges
    # averaged over all thresholds, and gives the integral above as an
    # efficient underestimate of it. Both are reported.
    #
    # This is computed from the UNMASKED counts, deliberately. A quadrature
    # over the threshold grid needs every node; masking one because its range
    # ran past the retained margin would not make the integral conservative,
    # it would make it undefined. The `saturated` flag still marks the time,
    # and is what a caller should test.
    if (peak[i] > threshold) {
      objective_exact[i] <- .aci_simpson(counts * dt, eps_sorted) / peak[i]
    }

    # ---- Resolution by censoring rather than deletion -----------------------
    #
    # A range that consumes most of the record it was measured against is not
    # a resolved range. It is also not an absence of measurement. It is a
    # RIGHT-CENSORED one, and the record does support a statement about it,
    # namely that the range is at least this long. Returning NA here would
    # throw that statement away at the end of the record, which is where a
    # user studying a recent event looks.
    #
    # The value therefore stands and `status` records that it is a bound.
    # `subjective`, `objective` and `objective_exact` are lower bounds wherever
    # the status is "censored", and `subjective_censored` marks the individual
    # thresholds that ran long.
    #
    # Resolution is judged per quantity rather than once for the whole time.
    # The subjective ranges at small thresholds run far longer than the
    # objective range does, and condemning the objective range because the most
    # demanding threshold was unresolved would discard the quantity the method
    # leads with.
    room <- (1 - margin) * length(re)
    subjective_censored[order(epsilon), i] <- counts > room
    subjective[order(epsilon), i] <- counts * dt

    settled <- length(suffix) - findInterval(threshold, rev(suffix))
    status[i] <- if (peak[i] <= threshold) {
      "below_threshold"
    } else if (settled > room) {
      "censored"
    } else {
      "resolved"
    }
  }

  # Retained with its original meaning (a time whose objective range is not
  # resolved) so that callers testing it keep working. The number beside it is
  # a bound rather than a hole.
  saturated <- status == "censored"

  structure(
    list(
      time = (window - 1L) * dt,
      index = window,
      objective = objective,
      objective_exact = objective_exact,
      subjective = subjective,
      subjective_censored = subjective_censored,
      status = status,
      monotone = monotone,
      epsilon = epsilon,
      peak = peak,
      saturated = saturated,
      dt = dt
    ),
    class = "aci_cir"
  )
}

#' One row of the causal-influence-range divergence sequence
#'
#' Forms the relative entropy of the fully informed online-smoother posterior
#' at step `j` from the posterior informed only up to each later observation.
#' The sequence is built from cumulative sums, so the row costs one pass over
#' the remaining record and is discarded by the caller once reduced.
#'
#' @param aux The auxiliary quantities from `.aci_online_aux()`.
#' @param filt A list with numeric vectors `mean` and `cov`.
#' @param j Integer scalar. The reported step.
#' @param n Integer scalar. Length of the observed signal.
#'
#' @returns A numeric vector of divergences, or `NULL` when `j` is the final
#'   step and no later observation exists.
#'
#' @noRd
#' @keywords internal
.aci_cir_row <- function(aux, filt, j, n, horizon = n) {
  if (j >= n) {
    return(NULL)
  }
  k <- seq.int(j, n - 1L)
  offset <- k - j
  # The ordered products for a whole row share the starting step, so the
  # cumulative logarithm is differenced once against a single base point.
  span <- aux$cum_log[k] - aux$cum_log[j]
  sgn <- aux$cum_sign[k] * aux$cum_sign[j]
  d <- sgn * exp(span)
  d[offset == 0L] <- 1

  # Posterior informed up to each later observation, as a running sum.
  mu <- filt$mean[j] + c(0, cumsum(d * aux$innov_mean[k]))
  rr <- filt$cov[j] + c(0, cumsum(d * d * aux$innov_cov[k]))
  mu_end <- mu[length(mu)]
  r_end <- rr[length(rr)]

  # Relative entropy of the fully informed posterior from each partial one,
  # with the dispersion term written so that it stays accurate where the two
  # posteriors nearly agree, which is precisely the tail that sets the range.
  if (any(!is.finite(rr)) || any(rr <= 0)) {
    bad <- which(!is.finite(rr) | rr <= 0)[1L]
    .aci_stop_covariance("online smoother", j + bad - 1L, NA_real_, rr[bad])
  }
  delta <- r_end / rr - 1
  value <- 0.5 * (mu_end - mu)^2 / rr + 0.5 * (delta - log1p(delta))

  # The horizon truncates how far forward the comparison looks, but NOT the
  # fully informed posterior it is compared against. Both `mu_end` and `r_end`
  # above are taken over the whole record either way. That asymmetry is the
  # point.
  # The quantity being measured is how far one must look before the estimate
  # stops changing, and the thing it must stop changing towards is the estimate
  # informed by everything available.
  pmax(value, 0)[seq_len(min(length(value), horizon - j + 1L))]
}

#' Validate a reporting window
#'
#' @param window The window to check, or `NULL` for the whole signal.
#' @param n Integer scalar. Length of the observed signal.
#'
#' @returns An integer vector of indices.
#'
#' @noRd
#' @keywords internal
.aci_check_window <- function(window, n) {
  if (is.null(window)) {
    return(seq_len(n))
  }
  ok <- is.numeric(window) && length(window) >= 1L && !anyNA(window) &&
    all(window == trunc(window)) && all(window >= 1) && all(window <= n)
  if (!ok) {
    stop(
      sprintf(
        paste0(
          "`window` must be whole numbers within the observed signal, which ",
          "has %d steps."
        ),
        n
      ),
      call. = FALSE
    )
  }
  as.integer(window)
}
