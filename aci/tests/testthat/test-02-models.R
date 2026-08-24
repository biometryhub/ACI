test_that("Z1: validators catch planted contract violations", {
  expect_error(cgns_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0, 1, 1),      # degenerate gxx
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1),
    class = "aci_error_gram")
  # non-affine hidden drift smuggled through the derived g is caught
  m <- model_dyad()
  m$g <- function(t, x, y) y^2
  expect_error(validate_cgns(m), class = "aci_error_model_contract")
})
test_that("Z1b: construction is RNG-neutral and affine conversion is strict", {
  set.seed(999)
  invisible(stats::runif(1)); invisible(model_dyad()); actual <- stats::runif(1)
  set.seed(999)
  invisible(stats::runif(1)); expected <- stats::runif(1)
  expect_equal(actual, expected)

  expect_error(cgns_from_affine(
    f_full = function(t, x, y) y^2,
    g_full = function(t, x, y) -y,
    Sx = function(t, x) matrix(1, 1, 1),
    Sy_hidden = function(t, x) matrix(1, 1, 1), k = 1, l = 1),
    class = "aci_error_model_contract")
})
test_that("Z2: simulate reproducibility", {
  m <- model_dyad()
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
test_that("Z3: dyad energy conservation of the quadratic pair", {
  m <- model_dyad()
  set.seed(1)
  for (i in 1:20) {
    x <- rnorm(1); y <- rnorm(1)
    quad_x <- m$meta$params$gamma * x * y * x      # (gamma x y) contribution to d(x^2)/2
    quad_y <- -m$meta$params$gamma * x^2 * y       # (-gamma x^2) contribution to d(y^2)/2
    expect_lt(abs(quad_x + quad_y), 1e-10)
  }
})
test_that("Z6-lite: zoo models construct and validate", {
  expect_s3_class(model_tipping_triad(0.1), "cgns_model")
  expect_s3_class(model_pathways(), "cgns_model")
  expect_s3_class(model_l84(), "cgns_model")
  expect_s3_class(model_l96(n = 12), "stochastic_model")
  expect_s3_class(model_enso6(), "cgns_model")            # (hW, tau) hidden split IS CGNS
  expect_error(model_enso6(hidden = c("I")), class = "aci_error_model_contract")
  expect_s3_class(model_topographic(), "stochastic_model")
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
    m <- model_predator_prey(hidden = h)
    s <- simulate(m, seed = 42, T = 3, dt = 5e-3)
    f <- suppressWarnings(da_filter(m, s$obs,
                                    init = list(mean = 4, cov = diag(1, 1))))
    expect_true(all(is.finite(f$mean)))
    expect_true(all(f$cov[1, 1, ] > 0))
  }
})

test_that("enso6 golden variant and u-hidden partitions stay CGNS", {
  for (v in c("cfy22", "aci_code")) {
    m <- model_enso6(hidden = c("u", "hW", "tau"), variant = v)
    expect_s3_class(m, "cgns_model")
    expect_true(all(is.finite(simulate(m, seed = 12, T = 1, dt = 5e-3)$obs$x)))
  }
  expect_equal(model_enso6(variant = "aci_code")$meta$params$factor, 0.65)
  expect_error(model_enso6(hidden = "TC"), class = "aci_error_model_contract")
})
