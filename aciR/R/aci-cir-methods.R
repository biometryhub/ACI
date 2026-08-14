# -- reporting surface for the causal influence range --------------------------
#
# The range is returned as a structured object because the quantities that
# matter are not all numbers of the same kind: an objective range, a matrix of
# subjective ranges over a threshold grid, a peak divergence, and -- since a
# range near the end of a record is right-censored rather than missing -- a
# status that says which of those a given time actually earned.
#
# Printing the object without methods dumps a threshold-by-time matrix, which
# buries the three scalars a reader needs.

#' Print a causal influence range
#'
#' @param x An object of class `aci_cir`.
#' @param ... Ignored.
#'
#' @returns `x`, invisibly.
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 900, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
#'                mu0 = model$y0, R0 = 0.1)
#' rng
#'
#' @export
print.aci_cir <- function(x, ...) {
  cat("Causal influence range\n")
  cat(sprintf(
    "  %d reported time%s, %s to %s; %d threshold%s\n",
    length(x$time), if (length(x$time) == 1L) "" else "s",
    format(min(x$time), digits = 4L), format(max(x$time), digits = 4L),
    length(x$epsilon), if (length(x$epsilon) == 1L) "" else "s"
  ))

  counts <- table(factor(
    x$status,
    levels = c("resolved", "censored", "below_threshold", "insufficient")
  ))
  cat("  status: ")
  cat(paste(sprintf("%s %d", names(counts), as.integer(counts)),
            collapse = "  "), "\n", sep = "")

  finite <- is.finite(x$objective)
  if (any(finite)) {
    cat(sprintf(
      "  objective range: median %s (%s to %s)\n",
      format(stats::median(x$objective[finite]), digits = 4L),
      format(min(x$objective[finite]), digits = 4L),
      format(max(x$objective[finite]), digits = 4L)
    ))
  }
  if (any(x$status == "censored")) {
    cat(sprintf(
      "  %d time%s censored -- those ranges are lower bounds\n",
      sum(x$status == "censored"),
      if (sum(x$status == "censored") == 1L) " is" else "s are"
    ))
  }
  invisible(x)
}

#' Summarise a causal influence range
#'
#' Reports the range alongside the two diagnostics a reader needs in order to
#' interpret it: how many reported times are censored rather than resolved, and
#' how often the divergence sequence is monotone. The second matters because
#' `objective` and `objective_exact` are the same functional only when it is.
#'
#' @param object An object of class `aci_cir`.
#' @param ... Ignored.
#'
#' @returns An object of class `summary.aci_cir`.
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 900, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' summary(aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
#'                 mu0 = model$y0, R0 = 0.1))
#'
#' @export
summary.aci_cir <- function(object, ...) {
  finite <- is.finite(object$objective)
  structure(
    list(
      n_time = length(object$time),
      n_epsilon = length(object$epsilon),
      span = range(object$time),
      status = table(factor(
        object$status,
        levels = c("resolved", "censored", "below_threshold", "insufficient")
      )),
      objective = if (any(finite)) {
        stats::quantile(object$objective[finite], c(0, 0.5, 1))
      } else {
        NULL
      },
      exact = if (any(is.finite(object$objective_exact))) {
        stats::median(object$objective_exact[is.finite(object$objective_exact)])
      } else {
        NA_real_
      },
      peak = stats::median(object$peak),
      censored_cells = mean(object$subjective_censored),
      monotone = mean(object$monotone, na.rm = TRUE)
    ),
    class = "summary.aci_cir"
  )
}

#' @rdname summary.aci_cir
#' @param x An object of class `summary.aci_cir`.
#' @export
print.summary.aci_cir <- function(x, ...) {
  cat("Causal influence range -- summary\n\n")
  cat(sprintf("  reported times : %d over [%s, %s]\n", x$n_time,
              format(x$span[1L], digits = 4L),
              format(x$span[2L], digits = 4L)))
  cat(sprintf("  thresholds     : %d\n", x$n_epsilon))
  cat("\n  resolution\n")
  for (nm in names(x$status)) {
    cat(sprintf("    %-16s %d\n", nm, as.integer(x$status[[nm]])))
  }
  cat(sprintf("    %-16s %.1f%% of subjective cells\n", "censored cells",
              100 * x$censored_cells))
  if (!is.null(x$objective)) {
    cat("\n  ranges\n")
    cat(sprintf("    objective        %s  (median %s, max %s)\n",
                format(x$objective[[1L]], digits = 4L),
                format(x$objective[[2L]], digits = 4L),
                format(x$objective[[3L]], digits = 4L)))
    cat(sprintf("    objective_exact  median %s\n",
                format(x$exact, digits = 4L)))
    cat(sprintf("    peak divergence  median %s\n",
                format(x$peak, digits = 4L)))
    cat(sprintf("\n  divergence monotone at %.0f%% of reported times\n",
                100 * x$monotone))
    if (isTRUE(x$monotone < 1)) {
      cat("  Where it is not, `objective` and `objective_exact` are different\n")
      cat("  functionals rather than two quadratures of one: the range is\n")
      cat("  measured as a last exit, which exceeds the superlevel measure.\n")
    }
  }
  invisible(x)
}

#' Coerce a causal influence range to a data frame
#'
#' The per-time quantities only. The subjective matrix is left behind because
#' it is threshold-by-time and does not belong in the same rectangle.
#'
#' @param x An object of class `aci_cir`.
#' @param row.names,optional Passed through to the data frame.
#' @param ... Ignored.
#'
#' @returns A data frame with one row per reported time.
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 900, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' head(as.data.frame(aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
#'                            mu0 = model$y0, R0 = 0.1)))
#'
#' @export
as.data.frame.aci_cir <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    time = x$time,
    index = x$index,
    objective = x$objective,
    objective_exact = x$objective_exact,
    peak = x$peak,
    status = x$status,
    saturated = x$saturated,
    row.names = row.names, stringsAsFactors = FALSE
  )
}

#' Plot a causal influence range
#'
#' Two panels: the objective range against time, with censored times drawn
#' hollow so a lower bound is not mistaken for a measurement, and the peak
#' divergence beneath it.
#'
#' @param x An object of class `aci_cir`.
#' @param y Ignored.
#' @param ... Passed to [graphics::plot()].
#'
#' @returns `x`, invisibly.
#'
#' @examples
#' model <- aci_dyad_model()
#' sim <- aci_simulate(model, n = 900, seed = 1)
#' comp <- aci_dyad_components(sim$x, model$parameters)
#' rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
#'                mu0 = model$y0, R0 = 0.1)
#' old <- graphics::par(mfrow = c(2, 1))
#' plot(rng)
#' graphics::par(old)
#'
#' @export
plot.aci_cir <- function(x, y, ...) {
  censored <- x$status == "censored"
  graphics::plot(
    x$time, x$objective, type = "n",
    xlab = "time", ylab = "objective range",
    main = "Causal influence range", ...
  )
  graphics::lines(x$time, x$objective, col = "grey60")
  graphics::points(x$time[!censored], x$objective[!censored], pch = 16L)
  # Hollow for censored: the marker says "at least this", not "this".
  graphics::points(x$time[censored], x$objective[censored], pch = 1L)

  graphics::plot(
    x$time, x$peak, type = "l",
    xlab = "time", ylab = "peak divergence",
    main = "Peak divergence"
  )
  invisible(x)
}
