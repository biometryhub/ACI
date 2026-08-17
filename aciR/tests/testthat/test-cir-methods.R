# The reporting methods for `aci_cir` had no tests at all until the coverage
# gate ran for the first time and reported the file at zero. They are exported,
# they are documented on the package site, and they are what a reader of the
# walkthrough sees, so an error in one of them is user-visible even though it
# touches no arithmetic.
#
# What is worth testing in a print method is not its wording, which will change,
# but the branches that decide what it says at all: whether a bound is announced
# as a bound, whether the singular and plural forms are reachable, and whether
# the quantities survive at the boundary where the record runs out. Asserting
# the exact text would make the tests break on every edit and say nothing about
# correctness, so these assert the branch was taken and the number reached the
# page.

.cir_case <- function(
    which = c("plain", "single", "censored", "below", "end")) {
  which <- match.arg(which)
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 900L, seed = 1L)
  comp <- aci_dyad_components(sim$x, model$parameters)
  args <- list(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
  window <- switch(which,
    plain = 50:150,
    single = 100L,
    censored = seq(600L, 880L, by = 40L),
    below = 50:60,
    end = 896:899
  )
  args$window <- window
  # A threshold above every attainable peak sends every time to
  # "below_threshold", which is the branch where the objective range is
  # reported as zero rather than as a ratio of two negligible quantities.
  if (which == "below") args$threshold <- 1e6
  do.call(aci_cir, args)
}

test_that("print reports resolution, and names a bound as one", {
  plain <- .cir_case("plain")
  out <- utils::capture.output(print(plain))
  expect_match(out[1L], "Causal influence range")
  expect_true(any(grepl("101 reported times", out)))
  expect_true(any(grepl("resolved 101", out)))
  expect_true(any(grepl("objective range: median", out)))
  # Nothing is censored here, so no bound should be announced.
  expect_false(any(grepl("lower bounds", out)))

  censored <- .cir_case("censored")
  out <- utils::capture.output(print(censored))
  expect_true(any(grepl("censored", out)))
  # The point of the line: a censored range is a bound, and the reader is told
  # so rather than left to infer it from a status column.
  expect_true(any(grepl("lower bounds", out)))
})

test_that("print returns its argument invisibly", {
  plain <- .cir_case("plain")
  expect_invisible(print(plain))
  expect_identical(utils::capture.output(x <- print(plain)),
                   utils::capture.output(print(plain)))
  expect_s3_class(x, "aci_cir")
})

test_that("the singular form is reachable for a one-time range", {
  one <- .cir_case("single")
  expect_length(one$time, 1L)
  out <- utils::capture.output(print(one))
  expect_true(any(grepl("1 reported time,", out)))
  expect_false(any(grepl("1 reported times", out)))
})

test_that("summary carries the two diagnostics needed to read the range", {
  plain <- .cir_case("plain")
  s <- summary(plain)
  expect_s3_class(s, "summary.aci_cir")
  expect_identical(s$n_time, length(plain$time))
  expect_identical(s$n_epsilon, length(plain$epsilon))
  expect_equal(s$span, range(plain$time))
  expect_named(s$status)
  # `monotone` is the condition under which the two objective ranges are the
  # same functional, so it has to survive into the summary object rather than
  # existing only inside a print method.
  expect_true(is.numeric(s$monotone))
  expect_gte(s$monotone, 0)
  expect_lte(s$monotone, 1)
  expect_true(is.numeric(s$censored_cells))
})

test_that("summary reports an absent exact range as NA", {
  below <- .cir_case("below")
  expect_true(all(below$status == "below_threshold"))
  expect_true(all(is.na(below$objective_exact)))
  s <- summary(below)
  expect_true(is.na(s$exact))
  # The objective range is still finite here, reported as zero, so the summary
  # must not confuse "no detectable influence" with "not measured".
  expect_false(is.null(s$objective))
})

test_that("print.summary shows the ranges block and the non-monotone warning", {
  plain <- .cir_case("plain")
  s <- summary(plain)
  out <- utils::capture.output(print(s))
  expect_match(out[1L], "Causal influence range, summary")
  expect_true(any(grepl("reported times", out)))
  expect_true(any(grepl("resolution", out)))
  expect_true(any(grepl("objective_exact", out)))
  expect_true(any(grepl("peak divergence", out)))
  expect_true(any(grepl("divergence monotone at", out)))
  # This system is never monotone, so the explanation of why the two ranges
  # differ must appear. Its absence would leave the gap looking like a defect.
  expect_lt(s$monotone, 1)
  expect_true(any(grepl("different", out)))
  expect_true(any(grepl("last exit", out)))
  expect_invisible(print(s))
})

test_that("as.data.frame gives one row per time and keeps the diagnostics", {
  plain <- .cir_case("plain")
  d <- as.data.frame(plain)
  expect_s3_class(d, "data.frame")
  expect_identical(nrow(d), length(plain$time))
  expect_named(
    d,
    c("time", "index", "objective", "objective_exact", "peak", "status",
      "monotone", "saturated")
  )
  expect_type(d$status, "character")
  expect_type(d$monotone, "logical")
  expect_type(d$saturated, "logical")
  expect_equal(d$time, plain$time)
  expect_equal(d$objective, plain$objective)
  # The subjective matrix is threshold-by-time and deliberately left out; a
  # rectangle with one row per time cannot hold it.
  expect_false("subjective" %in% names(d))
})

test_that("as.data.frame honours row.names", {
  one <- .cir_case("single")
  d <- as.data.frame(one, row.names = "only")
  expect_identical(rownames(d), "only")
})

test_that("the data frame survives the end of the record", {
  end <- .cir_case("end")
  d <- as.data.frame(end)
  expect_identical(nrow(d), length(end$time))
  # A time with too few later observations to support any quadrature returns
  # NA rather than a number, and the status says which.
  expect_true(any(end$status == "insufficient"))
  expect_true(any(is.na(d$objective)))
  # Everything else still reports, so one unresolvable time does not discard
  # the rest of the window.
  expect_true(any(is.finite(d$objective)))
})

test_that("plot draws both panels and restores graphics parameters", {
  plain <- .cir_case("plain")
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(path)
  }, add = TRUE)

  before <- graphics::par("mfrow")
  expect_invisible(plot(plain))
  # The method sets a two-panel layout, so it owns restoring it. Without the
  # restore, an interactive session keeps the split for every later plot.
  expect_identical(graphics::par("mfrow"), before)

  # A censored time is drawn hollow so a bound is not read as a measurement,
  # which means the censored branch has to be exercised too.
  expect_invisible(plot(.cir_case("censored")))
})
