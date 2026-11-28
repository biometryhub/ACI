# A floored filter's record follows it into everything built from it ----------
#
# Reusing a stored table carries its record and raises nothing new.

prov_capture <- function(expr) {
  classes <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      classes <<- c(classes, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  list(value = value, n = sum(classes == "aci_warn_regularized"))
}

prov_roundtrip <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(x, f)
  readRDS(f)
}

## Each case floors in its filter. `policy` is what a consuming call needs to
## run at all: the ENSO smoothers refuse the floored filter under "none".
prov_cases <- function() {
  dyad <- aci_dyad_model()
  generic <- aci_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.05, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
  enso <- aci_enso_model()
  list(
    dyad = list(model = dyad,
                obs = as_obs(simulate(dyad, seed = 1, t_end = 1, dt = 0.01)),
                init = list(mean = 0, cov = matrix(20)), policy = "none"),
    generic = list(model = generic,
                   obs = as_obs(simulate(generic, seed = 1, t_end = 0.6,
                                         dt = 0.01,
                                         ic = list(x0 = 0, y0 = 0))),
                   init = list(mean = 0, cov = matrix(5)), policy = "none"),
    enso = list(model = enso,
                obs = as_obs(simulate(enso, seed = 12, t_end = 0.5,
                                      dt = 0.005, burn_in = 0)),
                init = NULL, policy = "floor"))
}

prov_filter <- function(d)
  suppressWarnings(aci_filter(d$model, d$obs, init = d$init,
                              regularize = "floor"))


test_that("the hidden dimensions of the cases are one and two", {
  d <- prov_cases()
  expect_identical(d$dyad$model$l, 1L)
  expect_identical(d$generic$model$l, 1L)
  expect_identical(d$enso$model$l, 2L)
})


test_that("lag_table() on a floored filter carries its record and warns once", {
  for (nm in names(prov_cases())) {
    d <- prov_cases()[[nm]]
    f <- prov_filter(d)
    rf <- f$meta$regularization
    expect_true(rf$fired, info = nm)
    tb <- prov_capture(lag_table(d$model, d$obs, filter = f,
                                 regularize = d$policy))
    expect_identical(tb$n, 1L, info = nm)
    rt <- tb$value$meta$regularization
    expect_identical(rt$policy, "floor", info = nm)
    expect_true(rt$fired, info = nm)
    expect_identical(rt$sites[rt$sites$site %in% rf$sites$site, ], rf$sites,
                     info = nm)
    ## the same record as a table that floors its own filter
    own <- suppressWarnings(lag_table(d$model, d$obs, init = d$init,
                                      regularize = "floor"))
    expect_identical(rt, own$meta$regularization, info = nm)

    ## reuse carries the record and raises nothing new
    a <- prov_capture(aci(d$model, d$obs, table = tb$value))
    expect_identical(a$n, 0L, info = nm)
    expect_identical(a$value$meta$regularization, rt, info = nm)
    expect_match(paste(capture.output(print(a$value)), collapse = "\n"),
                 "regularized:", fixed = TRUE, info = nm)
    cr <- prov_capture(aci_range(tb$value, min_M = 0))
    expect_identical(cr$n, 0L, info = nm)
    expect_identical(cr$value$meta$regularization, rt, info = nm)
  }
})


test_that("aci_online() and aci_smoother() carry a supplied filter's record", {
  for (nm in names(prov_cases())) {
    d <- prov_cases()[[nm]]
    f <- prov_filter(d)
    on <- prov_capture(aci_online(d$model, d$obs, lag = 5, filter = f,
                                  regularize = d$policy))
    expect_identical(on$n, 1L, info = nm)
    own <- suppressWarnings(aci_online(d$model, d$obs, lag = 5, init = d$init,
                                       regularize = "floor"))
    expect_identical(on$value$meta$regularization, own$meta$regularization,
                     info = nm)
    expect_identical(on$value$mean, own$mean, info = nm)
    expect_identical(on$value$cov, own$cov, info = nm)

    sm <- prov_capture(aci_smoother(d$model, d$obs, filter = f,
                                    regularize = "floor"))
    expect_identical(sm$n, 1L, info = nm)
    own <- suppressWarnings(aci_smoother(d$model, d$obs, init = d$init,
                                         regularize = "floor"))
    expect_identical(sm$value$meta$regularization, own$meta$regularization,
                     info = nm)
  }
})


test_that("aci_range() reports a floored aci() result on every keep route", {
  for (nm in names(prov_cases())) {
    d <- prov_cases()[[nm]]
    res <- lapply(c(paths = "paths", none = "none", table = "table"),
                  function(keep)
                    suppressWarnings(aci(d$model, d$obs, init = d$init,
                                         keep = keep, regularize = "floor")))
    ## the table built inside aci() shares its record, counted once
    expect_identical(res$table$table$meta$regularization,
                     res$table$meta$regularization, info = nm)
    streamed <- prov_capture(aci_range(res$paths, min_M = 0))
    rebuilt <- prov_capture(aci_range(res$none, min_M = 0))
    reused <- prov_capture(aci_range(res$table, min_M = 0))
    expect_identical(streamed$n, 1L, info = nm)
    expect_identical(rebuilt$n, 1L, info = nm)
    expect_identical(reused$n, 0L, info = nm)
    expect_true(streamed$value$meta$regularization$fired, info = nm)
    expect_identical(streamed$value$meta$regularization,
                     rebuilt$value$meta$regularization, info = nm)
    expect_identical(reused$value$meta$regularization,
                     res$table$meta$regularization, info = nm)
  }
})


test_that("the record survives saveRDS() of the filter, the table and the result", {
  for (nm in c("dyad", "enso")) {
    d <- prov_cases()[[nm]]
    f <- prov_filter(d)
    ## the model travels with the filter, which checks it by identity
    back <- prov_roundtrip(list(model = d$model, filter = f))
    tb <- prov_capture(lag_table(back$model, d$obs, filter = back$filter,
                                 regularize = d$policy))
    expect_identical(tb$n, 1L, info = nm)
    fresh <- suppressWarnings(lag_table(d$model, d$obs, filter = f,
                                        regularize = d$policy))
    expect_identical(tb$value$meta$regularization,
                     fresh$meta$regularization, info = nm)

    tbr <- prov_roundtrip(fresh)
    a <- prov_capture(aci(tbr$meta$source_model, d$obs, table = tbr))
    expect_identical(a$n, 0L, info = nm)
    expect_identical(a$value$meta$regularization, fresh$meta$regularization,
                     info = nm)

    res <- suppressWarnings(aci(d$model, d$obs, init = d$init, keep = "paths",
                                regularize = "floor"))
    live <- prov_capture(aci_range(res, min_M = 0))
    reloaded <- prov_capture(aci_range(prov_roundtrip(res), min_M = 0))
    expect_identical(live$n, 1L, info = nm)
    expect_identical(reloaded$n, 1L, info = nm)
    expect_identical(reloaded$value$meta$regularization,
                     live$value$meta$regularization, info = nm)
    expect_identical(reloaded$value$tau, live$value$tau, info = nm)
  }
})


test_that("a filter that took no floor changes nothing", {
  d <- prov_cases()$dyad
  ini <- list(mean = 0, cov = matrix(0.5))
  f <- aci_filter(d$model, d$obs, init = ini, regularize = "floor")
  expect_false(f$meta$regularization$fired)
  tb <- prov_capture(lag_table(d$model, d$obs, filter = f))
  expect_identical(tb$n, 0L)
  expect_identical(tb$value, lag_table(d$model, d$obs, init = ini))
  on <- prov_capture(aci_online(d$model, d$obs, lag = 5, filter = f))
  expect_identical(on$n, 0L)
  expect_identical(on$value$meta$regularization$policy, "none")
  expect_false(on$value$meta$regularization$fired)
  a <- aci(d$model, d$obs, init = ini, keep = "paths")
  cr <- prov_capture(aci_range(a, min_M = 0))
  expect_identical(cr$n, 0L)
  expect_identical(cr$value$meta$regularization$policy, "none")
  expect_false(cr$value$meta$regularization$fired)
})
