test_that("plot methods pass caller labels and limits without duplicate arguments", {
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
})

test_that("a single-anchor CIR forwards title arguments", {
  cr <- new_cir_result(0, 1, 1, "forward", "exact", 1, matrix(c(0, 1), 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(cr, main = "One anchor", xlab = "Time"), cr)
  seen <- NULL
  capture <- function() {
    seen <<- mget(c("main", "xlab"), parent.frame())
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
  expect_no_warning(plot(cr, main = expression(tau[1]), axes = FALSE,
                         type = "p", ylim = c(0, 2)))
  expect_identical(seen$main, expression(tau[1]))
})
