# -- causal influence range ----------------------------------------------------
#
# The second quantity of the method paper. Where the causal-information metric
# asks HOW MUCH the future of the observed signal says about the unobserved
# state at time t, the causal influence range asks HOW FAR INTO THAT FUTURE one
# must look before the answer stops changing.
#
# The construction is a divergence between two online-smoother posteriors at
# the same instant: the one informed by the whole record, and the one informed
# only up to some later observation. As that later observation moves forward
# the two coincide, and the range is where the gap falls below a threshold.
#
# The reference implementation stores the whole matrix of those divergences --
# one row per time, one column per lagged observation -- which is quadratic in
# the window and reaches several gigabytes at the scale of the published
# figure. Nothing downstream needs the matrix: each row is reduced immediately
# to a handful of scalars. This implementation therefore forms one row at a
# time and discards it, which leaves the memory linear in the window while the
# arithmetic stays what it has to be.

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
#' accordingly: a few thousand steps is comfortable, and reproducing a
#' figure at the scale of the published one is a batch computation rather than
#' an interactive one.
#'
#' A reported time close to the end of the record cannot be resolved, because
#' the observations that would settle it do not exist. Rather than return the
#' truncated value such a time yields, the result marks it as saturated and
#' returns `NA` for its ranges.
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
#'   evaluated, in nats. Defaults to a logarithmic grid from `1e-6` to
#'   `10^0.5`, the grid of the reference implementation.
#' @param threshold Numeric scalar. A peak divergence below this value is
#'   treated as no detectable influence, and the objective range is reported as
#'   zero rather than as the ratio of two negligible quantities. Defaults to
#'   `1e-5`.
#' @param margin Numeric scalar in `(0, 1)`. The fraction of each comparison
#'   sequence that must remain unused for the range at that time to count as
#'   resolved. A time whose range consumes more than `1 - margin` of the record
#'   available after it is marked saturated and its ranges returned as `NA`.
#'   Defaults to `0.1`.
#' @param mu0,R0 Numeric scalars. Initial filtered mean and covariance, used
#'   only when `filt` is `NULL`.
#'
#' @returns An object of class `aci_cir`, a list with the reported `time`, the
#'   `objective` range at each time, the `subjective` range as a matrix with
#'   one row per threshold and one column per time, the `epsilon` grid, the
#'   `peak` divergence at each time, and the logical `saturated` marking times
#'   the record is too short to resolve.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*.
#' \doi{10.48550/arXiv.2411.05870}
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 600, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:250,
#'                mu0 = model$y0, R0 = 0.1)
#' summary(rng$objective)
#'
#' @seealso [aci_online_smoother()], [aci_metric()]
#' @export
aci_cir <- function(x, comp, dt, filt = NULL, window = NULL,
                    epsilon = 10^seq(-6, 0.5, length.out = 129L),
                    threshold = 1e-5, margin = 0.1, mu0 = NULL, R0 = NULL) {
  .aci_check_signal(x)
  n <- length(x)
  .aci_check_components(comp, n)
  .aci_check_positive(dt, "dt")
  .aci_check_positive(threshold, "threshold")
  .aci_check_scalar(margin, "margin")
  if (margin <= 0 || margin >= 1) {
    stop("`margin` must lie strictly between zero and one.", call. = FALSE)
  }
  if (is.null(filt)) {
    .aci_check_scalar(mu0, "mu0")
    .aci_check_positive(R0, "R0")
    filt <- aci_filter(x, comp, dt, mu0 = mu0, R0 = R0)
  }
  .aci_check_posterior(filt, "filt", n)
  window <- .aci_check_window(window, n)
  if (!is.numeric(epsilon) || length(epsilon) < 1L || anyNA(epsilon) ||
      any(epsilon <= 0)) {
    stop(
      "`epsilon` must be a vector of positive, non-missing thresholds.",
      call. = FALSE
    )
  }

  aux <- .aci_online_aux(x, comp, dt, filt)
  eps_sorted <- sort(epsilon)
  n_eps <- length(epsilon)

  objective <- numeric(length(window))
  peak <- numeric(length(window))
  saturated <- logical(length(window))
  subjective <- matrix(NA_real_, nrow = n_eps, ncol = length(window))

  for (i in seq_along(window)) {
    j <- window[i]
    re <- .aci_cir_row(aux, filt, j, n)
    if (is.null(re)) {
      # The final step has no later observation to be compared against.
      saturated[i] <- TRUE
      peak[i] <- 0
      next
    }
    peak[i] <- max(re)

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
    # never formally reach the end of the row -- testing against the end alone
    # would therefore never fire. The test is against the retained margin.
    #
    # Resolution is judged per quantity rather than once for the whole time.
    # The subjective ranges at small thresholds run far longer than the
    # objective range does, and condemning the objective range because the
    # most demanding threshold was unresolved would discard the quantity the
    # method leads with.
    room <- (1 - margin) * length(re)
    counts[counts > room] <- NA_integer_
    subjective[order(epsilon), i] <- counts * dt

    settled <- length(suffix) - findInterval(threshold, rev(suffix))
    saturated[i] <- settled > room
  }

  objective[saturated] <- NA_real_

  structure(
    list(
      time = (window - 1L) * dt,
      index = window,
      objective = objective,
      subjective = subjective,
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
#' @param aux The auxiliary quantities from [.aci_online_aux()].
#' @param filt A list with numeric vectors `mean` and `cov`.
#' @param j Integer scalar. The reported step.
#' @param n Integer scalar. Length of the observed signal.
#'
#' @returns A numeric vector of divergences, or `NULL` when `j` is the final
#'   step and no later observation exists.
#'
#' @noRd
#' @keywords internal
.aci_cir_row <- function(aux, filt, j, n) {
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
  # posteriors nearly agree -- which is precisely the tail that sets the range.
  if (any(!is.finite(rr)) || any(rr <= 0)) {
    bad <- which(!is.finite(rr) | rr <= 0)[1L]
    .aci_stop_covariance("online smoother", j + bad - 1L, NA_real_, rr[bad])
  }
  delta <- r_end / rr - 1
  value <- 0.5 * (mu_end - mu)^2 / rr + 0.5 * (delta - log1p(delta))
  pmax(value, 0)
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
