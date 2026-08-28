test_that("Z1: validators catch planted contract violations", {
  expect_error(aci_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0, 1, 1),      # degenerate gxx
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1),
    class = "aci_error_gram")
  # non-affine hidden drift smuggled through the derived g is caught
  m <- aci_dyad_model()
  m$g <- function(t, x, y) y^2
  expect_error(validate_cgns(m), class = "aci_error_model_contract")
})
test_that("Z1b: construction is RNG-neutral and affine conversion is strict", {
  set.seed(999)
  invisible(stats::runif(1)); invisible(aci_dyad_model()); actual <- stats::runif(1)
  set.seed(999)
  invisible(stats::runif(1)); expected <- stats::runif(1)
  expect_equal(actual, expected)

  expect_error(aci_model_from_affine(
    f_full = function(t, x, y) y^2,
    g_full = function(t, x, y) -y,
    Sx = function(t, x) matrix(1, 1, 1),
    Sy_hidden = function(t, x) matrix(1, 1, 1), k = 1, l = 1),
    class = "aci_error_model_contract")
})
test_that("Z2: simulate reproducibility", {
  m <- aci_dyad_model()
  s1 <- simulate(m, seed = 7, T = 1, dt = 1e-3)
  s2 <- simulate(m, seed = 7, T = 1, dt = 1e-3)
  expect_identical(s1$obs$x, s2$obs$x); expect_identical(s1$hidden, s2$hidden)
  expect_error(simulate(m, seed = 7, T = 1, dt = 1e-3, typo = TRUE),
               class = "aci_error_dims")
  expect_error(stochastic_model(
    f = function(t, x, y) y, g = function(t, x, y) -y,
    Sx = function(t, x) matrix(1, 1, 1),
    Sy = function(t, x, y) matrix(1, 1, 1), k = 1, l = 1,
    vectorized_members = "yes"), class = "aci_error_model_contract")
})
test_that("Z2b: seeded simulation contains the caller's generator state", {
  m <- aci_dyad_model()

  ## the seeded path is unchanged: set.seed(seed) still governs the draws
  set.seed(101)
  s1 <- simulate(m, seed = 7, T = 0.2, dt = 1e-2)
  set.seed(202)                                  # a different caller state
  s2 <- simulate(m, seed = 7, T = 0.2, dt = 1e-2)
  expect_identical(s1$obs$x, s2$obs$x)
  expect_identical(s1$hidden, s2$hidden)
  expect_identical(s1$noise$W, s2$noise$W)

  ## the caller's stream is where it was left
  set.seed(303)
  before <- .Random.seed
  invisible(simulate(m, seed = 7, T = 0.2, dt = 1e-2))
  expect_identical(.Random.seed, before)

  ## containment does not defeat the caller's own stream: the draw that follows
  ## a seeded call is the draw that would have followed no call at all
  set.seed(303); invisible(simulate(m, seed = 7, T = 0.2, dt = 1e-2))
  with_sim <- stats::runif(3)
  set.seed(303)
  expect_identical(with_sim, stats::runif(3))

  ## nsim > 1 draws nsim distinct realisations and still contains the state
  set.seed(404)
  before <- .Random.seed
  many <- simulate(m, nsim = 3, seed = 5, T = 0.2, dt = 1e-2)
  expect_length(many, 3L)
  expect_false(identical(many[[1L]]$obs$x, many[[2L]]$obs$x))
  expect_identical(.Random.seed, before)

  ## unseeded calls are untouched: they draw from, and advance, the caller
  set.seed(505)
  u1 <- simulate(m, T = 0.2, dt = 1e-2)
  after_u <- .Random.seed
  expect_false(identical(after_u, {set.seed(505); .Random.seed}))
  set.seed(505)
  u2 <- simulate(m, T = 0.2, dt = 1e-2)
  expect_identical(u1$obs$x, u2$obs$x)
  expect_identical(.Random.seed, after_u)

  ## a rejected call leaves the stream untouched, whether the argument is
  ## checked in this frame (nsim) or downstream in .simulate_one (T)
  set.seed(606)
  before <- .Random.seed
  expect_error(simulate(m, seed = 3, nsim = 0, T = 0.2, dt = 1e-2),
               class = "aci_error_dims")
  expect_identical(.Random.seed, before)
  expect_error(simulate(m, seed = 3, T = -1, dt = 1e-2), class = "aci_error_dims")
  expect_identical(.Random.seed, before)
  expect_error(simulate(m, seed = 3, T = 0.2, dt = 1e-2,
                        ic = list(x0 = c(0, 0), y0 = 0)), class = "aci_error_dims")
  expect_identical(.Random.seed, before)
})
test_that("Z2c: seeded simulation works with no pre-existing .Random.seed", {
  ## a session that has never drawn has no .Random.seed; reproduce that state
  ## here and put the caller's own back afterwards
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  keep <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", keep, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  if (had) rm(".Random.seed", envir = globalenv())
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  m <- aci_dyad_model()
  s <- simulate(m, seed = 7, T = 0.2, dt = 1e-2)

  ## the absent state was materialised rather than left missing, and the draws
  ## are exactly the ones set.seed(7) produces
  expect_true(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(7)
  expect_identical(s$obs$x, simulate(m, T = 0.2, dt = 1e-2)$obs$x)
})
test_that("Z3: dyad energy conservation of the quadratic pair", {
  m <- aci_dyad_model()
  set.seed(1)
  for (i in 1:20) {
    x <- rnorm(1); y <- rnorm(1)
    quad_x <- m$meta$params$gamma * x * y * x      # (gamma x y) contribution to d(x^2)/2
    quad_y <- -m$meta$params$gamma * x^2 * y       # (-gamma x^2) contribution to d(y^2)/2
    expect_lt(abs(quad_x + quad_y), 1e-10)
  }
})
test_that("Z6-lite: zoo models construct and validate", {
  expect_s3_class(aci_dyad_model(), "cgns_model")
  expect_s3_class(aci_predprey_model(), "cgns_model")
  expect_s3_class(aci_enso_model(), "cgns_model")            # (hW, tau) hidden split IS CGNS
  expect_error(aci_enso_model(hidden = c("I")), class = "aci_error_model_contract")
})
test_that("Z8: obs validators", {
  expect_error(observed_trajectory(c(0, .1, .3), matrix(0, 3, 1)),
               class = "aci_error_obs_contract")
  expect_error(observed_trajectory(c(0, .1, .2), matrix(c(0, NA, 0), 3, 1)),
               class = "aci_error_obs_contract")
  expect_error(as_obs(matrix(0, 3, 1), dt = c(0.1, 0.2)),
               class = "aci_error_obs_contract")
})

test_that("predator-prey golden model validates and assimilates both ways", {
  for (h in c("prey", "predator")) {
    m <- aci_predprey_model(hidden = h)
    s <- simulate(m, seed = 42, T = 3, dt = 5e-3)
    f <- suppressWarnings(aci_filter(m, s$obs,
                                     init = list(mean = 4, cov = diag(1, 1))))
    expect_true(all(is.finite(f$mean)))
    expect_true(all(f$cov[1, 1, ] > 0))
  }
})

test_that("enso6 golden variant and u-hidden partitions stay CGNS", {
  for (v in "aci_code") {
    m <- aci_enso_model(hidden = c("u", "hW", "tau"), variant = v)
    expect_s3_class(m, "cgns_model")
    expect_true(all(is.finite(simulate(m, seed = 12, T = 1, dt = 5e-3)$obs$x)))
  }
  expect_equal(aci_enso_model(variant = "aci_code")$meta$params$factor, 0.65)
  expect_error(aci_enso_model(hidden = "TC"), class = "aci_error_model_contract")
})

test_that("a simulation carries the model's observed-channel names", {
  ## The names come from meta$vars$observed, which the built-in constructors
  ## already record, and they are what aci_conditional() resolves block and target
  ## names against.
  for (cs in list(list(m = aci_dyad_model(), nm = "x"),
                  list(m = aci_predprey_model(hidden = "prey"),
                       nm = "predator"),
                  list(m = aci_predprey_model(hidden = "predator"),
                       nm = "prey"),
                  list(m = aci_enso_model(hidden = "hW"),
                       nm = c("u", "TC", "TE", "tau", "I")),
                  list(m = aci_enso_model(hidden = c("u", "hW", "tau")),
                       nm = c("TC", "TE", "I")))) {
    s <- simulate(cs$m, seed = 6, T = 0.2, dt = 0.005)
    expect_identical(colnames(s$obs$x), cs$nm)
    expect_identical(colnames(as_obs(s)$x), cs$nm)
    expect_identical(cs$nm, cs$m$meta$vars$observed)
    ## and a long-form coercion now labels the channels rather than x1, x2, ...
    expect_identical(unique(as.data.frame(s$obs)$var), cs$nm)
  }
  ## a model with no usable meta$vars is left exactly as it was
  m <- aci_dyad_model(); m$meta$vars <- NULL
  expect_null(colnames(simulate(m, seed = 6, T = 0.2, dt = 0.005)$obs$x))
  m$meta$vars <- list(observed = c("a", "b"))          # wrong length for k = 1
  expect_null(colnames(simulate(m, seed = 6, T = 0.2, dt = 0.005)$obs$x))
  m$meta$vars <- list(observed = "")                   # not a usable name
  expect_null(colnames(simulate(m, seed = 6, T = 0.2, dt = 0.005)$obs$x))
})


test_that("simulate() -> as_obs() -> conditional aci() works in three lines", {
  ## The path the README example had to work around: aci_conditional() resolves
  ## names against colnames(obs$x), which a simulation now supplies.
  me <- aci_enso_model(hidden = "hW")
  sim <- simulate(me, seed = 2, T = 0.5, dt = 0.005)
  ob <- as_obs(sim)
  a <- aci(me, ob,
           conditional = aci_conditional(target = "TC", method = "mask"),
           init = list(mean = 0, cov = 0.1))
  expect_s3_class(a, "aci_result")
  expect_true(all(is.finite(a$aci)))
  ## it is the same run as the hand-built trajectory it replaces, to the bit
  hand <- observed_trajectory(sim$obs$t, sim$obs$x,
                              names = me$meta$vars$observed)
  b <- aci(me, hand,
           conditional = aci_conditional(target = "TC", method = "mask"),
           init = list(mean = 0, cov = 0.1))
  expect_identical(a$aci, b$aci)
  ## and naming the complement as given is the same estimand
  cmp <- setdiff(me$meta$vars$observed, "TC")
  d <- aci(me, ob, conditional = aci_conditional(given = cmp, method = "mask"),
           init = list(mean = 0, cov = 0.1))
  expect_identical(a$aci, d$aci)
  ## the names are labels only: stripping them changes no number
  bare <- observed_trajectory(ob$t, unname(ob$x))
  e <- aci(me, bare,
           conditional = aci_conditional(target = 2L, method = "mask"),
           init = list(mean = 0, cov = 0.1))
  expect_identical(a$aci, e$aci)
  expect_identical(a$signal, e$signal)
  expect_identical(a$dispersion, e$dispersion)
})
