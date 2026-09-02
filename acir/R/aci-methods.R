################################################################################
## aci-methods.R - base-graphics plot methods for the core engine classes
## ########################################################################## ##
##
## Contents:
##   * plot.aci_result, plot.da_path_gaussian, plot.cir_result
##
## Moved verbatim from aci 0.0.30 R/formula_interface.R:617-674 (git tree
## 97f6b124); that file is an excluded family in acir 0.1.0, but these three
## methods dispatch on engine classes and have no aci_fit dependency.
##
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
  main_kl <- expression(
    ACI(t) == D[KL] * group("(", smoother ~ "||" ~ filter, ")")
  )
  plot(x$t, x$aci, type = "l", xlab = "t", ylab = "ACI", main = main_kl, ...)
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
  plot(x$t, mu, type = "l", xlab = "t", ylab = bquote(y[.(component)]),
       ylim = ylim, main = bquote(.(x$kind) ~ "mean" %+-% "2 sd"), ...)
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
         ylab = bquote(tau[.(substr(x$direction, 1, 1))]),
         main = sprintf("%s CIR (%s)", x$direction, x$method), ...)
  } else {
    graphics::plot.new()
    graphics::title(main = sprintf("%s CIR at t = %.4g: tau = %.4g",
                                   x$direction, x$t, x$tau))
  }
  invisible(x)
}
