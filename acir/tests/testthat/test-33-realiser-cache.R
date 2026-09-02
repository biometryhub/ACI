# The generic route's realisation cache: identical arrays on a hit, a fresh
# realisation when a captured parameter moves, a bounded store, an off switch,
# and no entry for the library route.

local_generic_model <- function(env) {
  aci_model(
    k = 1L, l = 1L,
    Lx = function(t, x) matrix(env$gam * x, 1L, 1L),
    fx = function(t, x) -0.5 * x,
    Ly = function(t, x) matrix(-1, 1L, 1L),
    fy = function(t, x) -env$gam * x^2,
    Sx1 = function(t, x) matrix(0.5, 1L, 1L),
    Sx2 = function(t, x) matrix(0, 1L, 1L),
    Sy1 = function(t, x) matrix(0, 1L, 1L),
    Sy2 = function(t, x) matrix(0.5, 1L, 1L)
  )
}

local_record <- function(n, seed = 1L) {
  as_obs(simulate(aci_dyad_model(), seed = seed, t_end = n * 0.01, dt = 0.01,
                  burn_in = 0))
}

init0 <- list(mean = 0, cov = matrix(1, 1L, 1L))

cold <- function(expr) {
  old <- options(aci.realiser_cache = FALSE)
  on.exit(options(old), add = TRUE)
  force(expr)
}

test_that("a second verb on the same model and record reuses the realisation", {
  .realise_cache_clear()
  env <- new.env()
  env$gam <- 1
  m <- local_generic_model(env)
  obs <- local_record(60L)
  f_cold <- cold(aci_filter(m, obs, init = init0))
  s_cold <- cold(aci_smoother(m, obs, filter = f_cold))
  f_hot <- aci_filter(m, obs, init = init0)
  expect_length(.realise_cache$entries, 1L)
  s_hot <- aci_smoother(m, obs, filter = f_hot)
  expect_length(.realise_cache$entries, 1L)
  expect_identical(f_hot$mean, f_cold$mean)
  expect_identical(f_hot$cov, f_cold$cov)
  expect_identical(s_hot$mean, s_cold$mean)
  expect_identical(s_hot$cov, s_cold$cov)
})

test_that("a captured parameter that moves forces a fresh realisation", {
  .realise_cache_clear()
  env <- new.env()
  env$gam <- 1
  m <- local_generic_model(env)
  obs <- local_record(60L)
  f1 <- aci_filter(m, obs, init = init0)
  env$gam <- 2
  f2 <- aci_filter(m, obs, init = init0)
  f2_cold <- cold(aci_filter(m, obs, init = init0))
  expect_false(identical(f2$mean, f1$mean))
  expect_identical(f2$mean, f2_cold$mean)
  expect_identical(f2$cov, f2_cold$cov)
  ## the stale entry was replaced, not kept beside the fresh one
  expect_length(.realise_cache$entries, 1L)
})

test_that("distinct records get distinct entries and the store is bounded", {
  .realise_cache_clear()
  env <- new.env()
  env$gam <- 1
  m <- local_generic_model(env)
  lengths <- 11:15
  for (n in lengths) aci_filter(m, local_record(n), init = init0)
  held <- vapply(.realise_cache$entries, function(e) length(e$t), integer(1))
  expect_length(held, .realise_cache_size)
  expect_false((11L + 1L) %in% held)
  expect_true(all((12:15 + 1L) %in% held))
})

test_that("the option switches the cache off; the library route adds nothing", {
  .realise_cache_clear()
  env <- new.env()
  env$gam <- 1
  m <- local_generic_model(env)
  obs <- local_record(30L)
  cold(aci_filter(m, obs, init = init0))
  expect_length(.realise_cache$entries, 0L)
  aci_filter(aci_dyad_model(), obs, init = init0)
  expect_length(.realise_cache$entries, 0L)
})
