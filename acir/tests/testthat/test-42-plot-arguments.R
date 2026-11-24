test_that(paste0("plot methods pass caller labels and limits without ",
                 "duplicate arguments"), {
  t <- 0:2
  a <- structure(list(t = t, aci = c(1, 2, 1), signal = NULL),
                 class = "aci_result")
  p <- new_da_path(t, matrix(c(0, 1, 0), 3, 1),
                   array(1, c(1, 1, 3)), "filter")
  cr <- new_cir_result(t, c(1, 2, 1), rep(1, 3), "forward", "exact", 1,
                       cbind(t, t + c(1, 2, 1)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (x in list(a, p, cr))
    expect_identical(plot(x, main = "Custom title", xlab = "Time",
                           ylab = "Value", type = "p", ylim = c(-4, 4)), x)
  original <- utils::getS3method("plot", "default")
  seen <- NULL
  capture <- function(...) {
    seen <<- list(...)
    original(...)
  }
  registerS3method("plot", "default", capture, envir = asNamespace("graphics"))
  on.exit(registerS3method("plot", "default", original,
                          envir = asNamespace("graphics")), add = TRUE)
  for (x in list(a, p, cr)) {
    expect_identical(plot(x), x)
    expect_identical(seen$type, "l")
    expect_identical(seen$xlab, "t")
    expect_identical(plot(x, main = "Custom title", xlab = "Time",
                           ylab = "Value", type = "p", ylim = c(-4, 4)), x)
    expect_identical(seen$main, "Custom title")
    expect_identical(seen$xlab, "Time")
    expect_identical(seen$ylab, "Value")
    expect_identical(seen$type, "p")
    expect_identical(seen$ylim, c(-4, 4))
    expect_identical(plot(x, main = NULL, col = "red"), x)
    expect_true("main" %in% names(seen))
    expect_null(seen$main)
    expect_identical(seen$col, "red")
  }
  ## Defaults that no argument replaced.
  plot(p)
  expect_identical(seen$ylim, c(-2, 3))
  plot(cr)
  expect_identical(seen$main, "forward CIR (exact)")
  expect_identical(seen$ylab, quote(tau["f"]))
})

test_that(paste0("plot methods pass unnamed arguments to plot() by ",
                 "position, as a direct call does"), {
  t <- 0:2
  a <- structure(list(t = t, aci = c(1, 2, 1), signal = NULL),
                 class = "aci_result")
  p <- new_da_path(t, matrix(c(0, 1, 0), 3, 1),
                   array(1, c(1, 1, 3)), "filter")
  cr <- new_cir_result(t, c(1, 2, 1), rep(1, 3), "forward", "exact", 1,
                       cbind(t, t + c(1, 2, 1)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ## The first unnamed argument after a method's own arguments matches xlim,
  ## the first argument of plot.default() that the method leaves unset. With
  ## the default axis style, xlim = c(0, 5) spans -0.2 to 5.2. The default
  ## range of t, from 0 to 2, would span -0.08 to 2.08.
  expect_identical(plot(cr, c(0, 5)), cr)
  expect_equal(graphics::par("usr")[1:2], c(-0.2, 5.2), tolerance = 1e-12)
  expect_identical(plot(a, TRUE, c(0, 5)), a)
  expect_equal(graphics::par("usr")[1:2], c(-0.2, 5.2), tolerance = 1e-12)
  expect_identical(plot(p, 1, NULL, c(0, 5)), p)
  expect_equal(graphics::par("usr")[1:2], c(-0.2, 5.2), tolerance = 1e-12)
  ## A second unnamed argument matches ylim where the method leaves it unset.
  expect_identical(plot(cr, c(0, 5), c(-1, 3)), cr)
  expect_equal(graphics::par("usr")[3:4], c(-1.16, 3.16), tolerance = 1e-12)
  ## An empty argument leaves xlim at its default, so the next one matches
  ## ylim.
  expect_identical(plot(cr, , c(-1, 3)), cr)
  expect_equal(graphics::par("usr")[3:4], c(-1.16, 3.16), tolerance = 1e-12)
  ## A trailing empty argument after the methods' own arguments.
  expect_identical(plot(a, TRUE, col = "red", ), a)
  expect_identical(plot(p, 1, NULL, col = "red", ), p)
  expect_identical(plot(cr, col = "red", ), cr)
})

test_that("panel.first and panel.last draw in the plot the method starts", {
  t <- 0:2
  a <- structure(list(t = t, aci = c(1, 2, 1), signal = NULL),
                 class = "aci_result")
  p <- new_da_path(t, matrix(c(0, 1, 0), 3, 1),
                   array(1, c(1, 1, 3)), "filter")
  cr <- new_cir_result(t, c(1, 2, 1), rep(1, 3), "forward", "exact", 1,
                       cbind(t, t + c(1, 2, 1)))
  ## Operations recorded on a new device by code, closing the device even
  ## when the code fails.
  operations <- function(code) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control("enable")
    force(code)
    vapply(grDevices::recordPlot()[[1]],
           function(op) as.list(op[[2]])[[1]]$name, "")
  }
  for (x in list(a, p, cr)) {
    ## On a device with no plot yet.
    ops <- operations(
      expect_identical(plot(x, panel.first = graphics::abline(h = 0.5),
                            panel.last = graphics::abline(v = 1)), x)
    )
    expect_identical(sum(ops == "C_abline"), 2L)
    ## On a device that already holds a figure, the line belongs to the
    ## second figure, the one the method draws.
    ops <- operations({
      graphics::par(mfrow = c(1, 2))
      plot(0:2, 0:2)
      plot(x, panel.first = graphics::abline(h = 0.5))
    })
    expect_identical(cumsum(ops == "C_plot_new")[ops == "C_abline"], 2L)
  }
})

test_that("plot methods evaluate an argument only when plot() uses it", {
  t <- 0:2
  a <- structure(list(t = t, aci = c(1, 2, 1), signal = NULL),
                 class = "aci_result")
  p <- new_da_path(t, matrix(c(0, 1, 0), 3, 1),
                   array(1, c(1, 1, 3)), "filter")
  cr <- new_cir_result(t, c(1, 2, 1), rep(1, 3), "forward", "exact", 1,
                       cbind(t, t + c(1, 2, 1)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (x in list(a, p, cr))
    expect_identical(plot(x, ann = FALSE, sub = stop("sub was evaluated")),
                     x)
})

test_that("a range whose values are all NA is plotted only with a ylim", {
  t <- 0:2
  cr <- new_cir_result(t, rep(NA_real_, 3), rep(1, 3), "forward", "exact",
                       1, cbind(t, t))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(cr), class = "aci_error_dims")
  expect_error(plot(cr, ylim = NULL), class = "aci_error_dims")
  expect_error(plot(cr, , ), class = "aci_error_dims")
  with_null <- function(x, ylim = NULL) plot(x, ylim = ylim)
  expect_error(with_null(cr), class = "aci_error_dims")
  expect_identical(plot(cr, ylim = c(0, 1)), cr)
  expect_identical(plot(cr, yli = c(0, 1)), cr)
  expect_identical(plot(cr, NULL, c(0, 1)), cr)
})

test_that("a single-anchor CIR forwards title arguments", {
  cr <- new_cir_result(0, 1, 1, "forward", "exact", 1, matrix(c(0, 1), 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(cr, main = "One anchor", xlab = "Time"), cr)
  seen <- NULL
  capture <- function() {
    seen <<- mget(c("main", "sub", "xlab", "ylab"), parent.frame())
  }
  suppressMessages(trace("title", where = asNamespace("graphics"),
                         tracer = substitute(FUN(), list(FUN = capture)),
                         print = FALSE))
  on.exit(suppressMessages(untrace("title", where = asNamespace("graphics"))),
          add = TRUE)
  expect_identical(plot(cr), cr)
  expect_match(seen$main, "forward CIR at t = 0: tau = 1", fixed = TRUE)
  expect_identical(plot(cr, main = "One anchor", xlab = "Time"), cr)
  expect_identical(seen$main, "One anchor")
  expect_identical(seen$xlab, "Time")
  expect_identical(plot(cr, sub = "Sub", ylab = "Value"), cr)
  expect_identical(seen$sub, "Sub")
  expect_identical(seen$ylab, "Value")
  expect_no_warning(plot(cr, main = expression(tau[1]), axes = FALSE,
                         type = "p", ylim = c(0, 2)))
  expect_identical(seen$main, expression(tau[1]))
  ## Other arguments, unnamed or empty ones included, have no effect and are
  ## not evaluated.
  expect_no_warning(plot(cr, c(0, 5), , lwd = 0, cex = 0, pch = NA,
                         xlim = stop("xlim was evaluated"),
                         panel.first = stop("panel.first was evaluated")))
  expect_match(seen$main, "forward CIR at t = 0: tau = 1", fixed = TRUE)
})
