# -- methods for the model and result objects ---------------------------------
#
# The methods below are the reporting surface of the package. They exist so
# that inspecting, plotting and exporting a result never requires the user to
# index into the nested representation, and so that the quantities an ACI
# result is sensitive to -- the smallest posterior covariances, the terminal
# identity, the round-off floor of the metric -- are surfaced rather than left
# for the reader to discover.

# -- print methods ------------------------------------------------------------

#' Print a conditional Gaussian nonlinear system model
#'
#' Prints a compact, one-field-per-line summary of an `aci_model` object: its
#' label, the self-drift of the unobserved component, the noise Grammians, the
#' initial state and, when present, the named scalar parameters that define the
#' model.
#'
#' @param x An `aci_model` object; see [aci_cgns_model()].
#' @param ... Ignored, for compatibility with [print()].
#'
#' @returns The model `x`, invisibly.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' print(aci_dyad_model())
#'
#' @seealso [aci_cgns_model()], [aci_dyad_model()]
#' @export
print.aci_model <- function(x, ...) {
  cat("<aci_model> ", x$label, "\n", sep = "")
  cat(sprintf("  unobserved self-drift L_y: %g\n", x$L_y))
  cat(sprintf(
    "  noise Grammians: S_xoS_x = %g, S_yoS_y = %g, S_yoS_x = %g\n",
    x$S_xoS_x, x$S_yoS_y, x$S_yoS_x
  ))
  cat(sprintf("  initial state: x0 = %g, y0 = %g\n", x$x0, x$y0))
  if (!is.null(x$parameters)) {
    pars <- paste(names(x$parameters), unlist(x$parameters),
      sep = " = ", collapse = ", "
    )
    cat("  parameters: ", pars, "\n", sep = "")
  }
  invisible(x)
}

#' Print an assimilative causal inference result
#'
#' Prints a compact summary of an `aci` object: the model it was run for, the
#' number of steps, the integration step and time span, and a five-number
#' summary of the causal-information metric.
#'
#' @param x An `aci` object, as returned by [aci()].
#' @param ... Ignored, for compatibility with [print()].
#'
#' @returns The result `x`, invisibly.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' fit <- aci(aci_simulate(model, n = 2000, seed = 1)$x, model)
#' print(fit)
#'
#' @seealso [aci()], [summary.aci()]
#' @export
print.aci <- function(x, ...) {
  n <- length(x$x)
  cat("<aci> assimilative causal inference\n")
  cat("  model: ", x$model$label, "\n", sep = "")
  cat(sprintf(
    "  steps: %d, dt: %g, time span: [%g, %g]\n", n, x$dt, x$t[1L], x$t[n]
  ))
  cat("  causal-information metric:\n")
  print(summary(x$aci))
  invisible(x)
}

# -- summary ------------------------------------------------------------------

#' Summarise an assimilative causal inference result
#'
#' Summarises an `aci` object into the quantities a reader needs to judge the
#' result rather than only to read it: the distribution of the
#' causal-information metric, where and when it peaks, and the diagnostics that
#' say whether the numerical assumptions held.
#'
#' Three diagnostics are reported alongside the metric. The smallest filtered
#' and smoothed covariances indicate how close the explicit Euler
#' discretisation came to losing positivity; a value near zero means the result
#' is sensitive to `dt`. The terminal residual is the gap between the smoother
#' and the filter at the final step, which is zero analytically and so
#' measures accumulated numerical error. The count of clamped values reports
#' how many steps had a metric at the round-off floor, which is expected
#' wherever the two posteriors agree.
#'
#' @param object An `aci` object, as returned by [aci()].
#' @param ... Ignored, for compatibility with [summary()].
#'
#' @returns An object of class `summary.aci`: a list with the number of steps
#'   `n`, the step `dt`, the time `span`, the model `label`, the metric
#'   five-number summary `metric`, the peak metric value `peak` and its time
#'   `peak_time`, the smallest covariances `min_filter_cov` and
#'   `min_smoother_cov`, the `terminal_residual`, and the number of
#'   round-off-clamped metric values `n_clamped`.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' model <- aci_dyad_model()
#' fit <- aci(aci_simulate(model, n = 5000, seed = 333)$x, model)
#' summary(fit)
#'
#' @seealso [aci()], [plot.aci()], [as.data.frame.aci()]
#' @export
summary.aci <- function(object, ...) {
  n <- length(object$x)
  peak <- which.max(object$aci)
  out <- list(
    n = n,
    dt = object$dt,
    span = c(object$t[1L], object$t[n]),
    label = object$model$label,
    metric = summary(object$aci),
    peak = object$aci[peak],
    peak_time = object$t[peak],
    min_filter_cov = min(object$filter$cov),
    min_smoother_cov = min(object$smoother$cov),
    terminal_residual = abs(object$smoother$mean[n] - object$filter$mean[n]) +
      abs(object$smoother$cov[n] - object$filter$cov[n]),
    n_clamped = sum(object$aci == 0)
  )
  class(out) <- "summary.aci"
  out
}

#' Print a summary of an assimilative causal inference result
#'
#' @param x A `summary.aci` object, as returned by [summary.aci()].
#' @param ... Ignored, for compatibility with [print()].
#'
#' @returns The summary `x`, invisibly.
#'
#' @examples
#' model <- aci_dyad_model()
#' print(summary(aci(aci_simulate(model, n = 2000, seed = 1)$x, model)))
#'
#' @seealso [summary.aci()]
#' @export
print.summary.aci <- function(x, ...) {
  cat("<summary.aci> assimilative causal inference\n")
  cat("  model: ", x$label, "\n", sep = "")
  cat(sprintf(
    "  steps: %d, dt: %g, time span: [%g, %g]\n",
    x$n, x$dt, x$span[1L], x$span[2L]
  ))
  cat("\n  causal-information metric:\n")
  print(x$metric)
  cat(sprintf("  peak %g at time %g\n", x$peak, x$peak_time))
  cat("\n  diagnostics:\n")
  cat(sprintf(
    "    smallest covariance: filter %g, smoother %g\n",
    x$min_filter_cov, x$min_smoother_cov
  ))
  cat(sprintf(
    "    terminal identity residual: %g (zero analytically)\n",
    x$terminal_residual
  ))
  cat(sprintf(
    "    metric at the round-off floor: %d of %d steps\n",
    x$n_clamped, x$n
  ))
  invisible(x)
}

# -- extraction ---------------------------------------------------------------

#' Coerce an assimilative causal inference result to a data frame
#'
#' Flattens the posterior trajectories and the causal-information metric into
#' one tidy row per time step, for plotting, export or joining to other series.
#'
#' @param x An `aci` object, as returned by [aci()].
#' @param row.names Passed to [data.frame()].
#' @param optional Passed to [data.frame()].
#' @param ... Ignored, for compatibility with [as.data.frame()].
#'
#' @returns A data frame with one row per time step and the columns `t`, `x`,
#'   `filter_mean`, `filter_cov`, `smoother_mean`, `smoother_cov` and `aci`.
#'
#' @examples
#' model <- aci_dyad_model()
#' fit <- aci(aci_simulate(model, n = 2000, seed = 1)$x, model)
#' head(as.data.frame(fit))
#'
#' @seealso [aci()], [summary.aci()]
#' @export
as.data.frame.aci <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    t = x$t,
    x = x$x,
    filter_mean = x$filter$mean,
    filter_cov = x$filter$cov,
    smoother_mean = x$smoother$mean,
    smoother_cov = x$smoother$cov,
    aci = x$aci,
    row.names = row.names,
    check.names = !optional
  )
}

# -- plot ---------------------------------------------------------------------

#' Plot an assimilative causal inference result
#'
#' Draws the observed signal and the causal-information metric on a shared time
#' axis, the pair of panels that carries the method's reading: bursts in the
#' observed signal above, and the steps at which the future of that signal
#' sharpens the estimate of the unobserved component below.
#'
#' The method uses base graphics so that plotting costs the package no
#' dependency. The colours are chosen to remain distinguishable in the common
#' forms of colour vision deficiency. For a publication-grade figure, take
#' [as.data.frame()] of the result and draw it with the grammar of graphics;
#' the *Assimilative causal inference on the nonlinear dyad model* vignette
#' shows the recipe.
#'
#' @param x An `aci` object, as returned by [aci()].
#' @param y Ignored, for compatibility with [plot()].
#' @param ... Further graphical parameters passed to the underlying [plot()]
#'   calls.
#'
#' @returns The result `x`, invisibly. Called for the plot it draws.
#'
#' @examples
#' model <- aci_dyad_model()
#' fit <- aci(aci_simulate(model, n = 5000, seed = 333)$x, model)
#' plot(fit)
#'
#' @seealso [aci()], [as.data.frame.aci()]
#' @export
plot.aci <- function(x, y, ...) {
  old <- graphics::par(mfrow = c(2L, 1L), mar = c(4.1, 4.1, 2.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)

  graphics::plot(x$t, x$x,
    type = "l", col = "#0072B2",
    xlab = "", ylab = "observed signal",
    main = x$model$label, ...
  )
  graphics::plot(x$t, x$aci,
    type = "l", col = "#D55E00",
    xlab = "time", ylab = "causal information",
    ...
  )
  invisible(x)
}
