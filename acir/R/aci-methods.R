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
#' @param ... Further arguments to [plot()], which draws the total. A `type`,
#'   `xlab`, `ylab` or `main` given here replaces the default.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.aci_result <- function(x, decompose = TRUE, ...) {
  main_kl <- expression(
    ACI(t) == D[KL] * group("(", smoother ~ "||" ~ filter, ")")
  )
  ## The defaults follow the dots, so a caller's argument of the same name
  ## replaces a default, and every other argument reaches plot() unevaluated
  ## and in its place, as in a direct call.
  draw <- function(..., type = "l", xlab = "t", ylab = "ACI",
                   main = main_kl) {
    plot(x$t, x$aci, type = type, xlab = xlab, ylab = ylab, main = main, ...)
  }
  draw(...)
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
#' @param ... Further arguments to [plot()], which draws the mean. A `type`,
#'   `xlab`, `ylab`, `ylim` or `main` given here replaces the default.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.da_path_gaussian <- function(x, component = 1, truth = NULL, ...) {
  mu <- x$mean[, component]; sdv <- sqrt(pmax(x$cov[component, component, ], 0))
  band <- range(mu + 2 * sdv, mu - 2 * sdv, truth, na.rm = TRUE)
  ## The defaults follow the dots, as in plot.aci_result().
  draw <- function(..., type = "l", xlab = "t",
                   ylab = bquote(y[.(component)]), ylim = band,
                   main = bquote(.(x$kind) ~ "mean" %+-% "2 sd")) {
    plot(x$t, mu, type = type, xlab = xlab, ylab = ylab, ylim = ylim,
         main = main, ...)
  }
  draw(...)
  graphics::lines(x$t, mu + 2 * sdv, lty = 3, col = "grey40")
  graphics::lines(x$t, mu - 2 * sdv, lty = 3, col = "grey40")
  if (!is.null(truth)) graphics::lines(x$t, truth, col = 2)
  invisible(x)
}


#' Plot a causal influence range result
#'
#' @param x A `cir_result` object.
#' @param ... Further arguments to [plot()]. A `type`, `xlab`, `ylab` or
#'   `main` given here replaces the default. A range with a single anchor is
#'   drawn as a title only, so only its `main`, `sub`, `xlab` and `ylab` are
#'   used, by [title()]. A range with more than one anchor whose values are
#'   all `NA`, as when every anchor is masked, needs a `ylim` that is not
#'   `NULL`.
#' @returns `x`, invisibly; called for the plot it draws.
#' @export
plot.cir_result <- function(x, ...) {
  if (length(x$tau) > 1) {
    ## The defaults follow the dots, as in plot.aci_result(). An all-NA
    ## series has no y range, so draw() then needs a ylim that is not NULL,
    ## given by name, partial name or position, and reports whether it drew.
    draw <- function(..., type = "l", xlab = "t",
                     ylab = bquote(tau[.(substr(x$direction, 1, 1))]),
                     main = sprintf("%s CIR (%s)", x$direction, x$method)) {
      if (all(is.na(x$tau))) {
        slot <- match.call(graphics::plot.default,
                           quote(plot(x$t, x$tau, type = type, xlab = xlab,
                                      ylab = ylab, main = main, ...)),
                           envir = environment())$ylim
        if (is.null(slot)) return(FALSE)
        ## The slot holds ..k, the k-th of the dots.
        k <- as.integer(substring(as.character(slot), 3L))
        if (identical(substitute(list(...))[[k + 1L]], quote(expr = )) ||
            is.null(...elt(k)))
          return(FALSE)
      }
      plot(x$t, x$tau, type = type, xlab = xlab, ylab = ylab, main = main,
           ...)
      TRUE
    }
    if (!draw(...))
      aci_abort("aci_error_dims", paste(
        "plot.cir_result: every range value is NA, as when every anchor is",
        "masked, so there is no y range to draw. Supply a ylim to draw the",
        "axes."))
  } else {
    graphics::plot.new()
    ## A single anchor is drawn as a title only. A caller's main, sub, xlab
    ## or ylab replaces the default, and other arguments are not evaluated.
    draw_title <- function(..., main = sprintf("%s CIR at t = %.4g: tau = %.4g",
                                               x$direction, x$t, x$tau),
                           sub = NULL, xlab = NULL, ylab = NULL) {
      graphics::title(main = main, sub = sub, xlab = xlab, ylab = ylab)
    }
    draw_title(...)
  }
  invisible(x)
}
