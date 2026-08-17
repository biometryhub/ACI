# The six-dimensional model's construction, and the influence range's
# requirements when the system is vector-valued.
#
# `aci_enso_model()` had no test at all: it is a description rather than a
# computation, so nothing else reached it, and a description that nothing
# reads is exactly where a wrong component name survives. Its own accessors
# are asserted here against the names the rest of the package uses.

# ---- the ENSO description ----------------------------------------------------

test_that("the ENSO model describes the split it implements", {
  m <- aci_enso_model()
  expect_named(
    m, c("label", "parameters", "observed", "unobserved", "components")
  )
  # These names appear in the documentation and on the package site, so they
  # are part of the interface rather than an internal label.
  expect_identical(m$observed, c("T_C", "T_E", "I"))
  expect_identical(m$unobserved, c("u", "h_W", "tau"))
  expect_type(m$label, "character")
  expect_type(m$components, "closure")
})

test_that("the components closure matches the function", {
  m <- aci_enso_model()
  T_C <- rep(0.01, 8L)
  T_E <- rep(0.01, 8L)
  I <- rep(1.5, 8L)
  via_closure <- m$components(T_C, T_E, I)
  direct <- aci_enso_components(T_C, T_E, I, p = m$parameters)
  # A convenience wrapper that quietly used different parameters from the
  # model it hangs off would be worse than no wrapper.
  expect_equal(via_closure, direct)
})

test_that("the model refuses a parameter list it did not come from", {
  expect_error(aci_enso_model("nope"), "parameter list")
  expect_error(aci_enso_model(list(a = 1)), "parameter list")
})

test_that("the seasonal scaling must be positive", {
  expect_error(aci_enso_parameters(factor = -1), "`factor` must be positive")
  expect_error(aci_enso_parameters(factor = 0), "`factor` must be positive")
})

test_that("the three observed paths are checked for type, length and count", {
  expect_error(
    aci_enso_components(c("a", "b"), 1:2, 1:2),
    "complete, finite numeric vector"
  )
  expect_error(
    aci_enso_components(c(0.01, NA), c(0.01, 0.01), c(1.5, 1.5)),
    "complete, finite numeric vector"
  )
  expect_error(
    aci_enso_components(1:3, 1:2, 1:2), "must have the same length"
  )
  # One observation is not a record: the smoother has nothing to run backwards
  # over.
  expect_error(aci_enso_components(1, 1, 1), "At least two observations")
})

test_that("observation times must match the paths", {
  T_C <- rep(0.01, 5L)
  expect_error(
    aci_enso_components(T_C, T_C, rep(1.5, 5L), time = 1:3),
    "as long as the observed paths"
  )
  # Supplied times are honoured; omitted times are built from dt. The
  # coefficients are modulated by the season, so this is not cosmetic.
  from_dt <- aci_enso_components(T_C, T_C, rep(1.5, 5L), dt = 0.005)
  explicit <- aci_enso_components(
    T_C, T_C, rep(1.5, 5L), time = (seq_len(5L) - 1L) * 0.005, dt = 0.005
  )
  expect_equal(from_dt, explicit)
})

# ---- the influence range on a vector system ----------------------------------

.vec_case <- function(n = 30L, n_x = 2L) {
  list(
    comp = list(
      L_x = array(rep(diag(-1, n_x, n_x), n), c(n_x, n_x, n)),
      f_x = matrix(0, n_x, n),
      L_y = array(rep(diag(-1, n_x, n_x), n), c(n_x, n_x, n)),
      f_y = matrix(0, n_x, n),
      S_xoS_x = diag(0.1, n_x),
      S_yoS_y = diag(0.1, n_x),
      S_yoS_x = matrix(0, n_x, n_x)
    ),
    x = matrix(stats::rnorm(n_x * n), n_x, n)
  )
}

test_that("a vector system requires an explicit filter posterior", {
  v <- .vec_case()
  # mu0 and R0 have scalar defaults that cannot be guessed in more than one
  # dimension, so the range refuses rather than inventing an initial condition.
  expect_error(
    aci_cir(v$x, v$comp, dt = 0.01, window = 5:10),
    "must be supplied for a vector system"
  )
})

test_that("the supplied posterior must be vector-valued and the right length", {
  v <- .vec_case()
  scalar_shaped <- list(mean = rnorm(30L), cov = rnorm(30L)^2)
  expect_error(
    aci_cir(v$x, v$comp, dt = 0.01, window = 5:10, filt = scalar_shaped),
    "vector-valued posterior"
  )

  # A posterior from a shorter record. Its own components list has to be the
  # matching length, or the filter call fails before the range is reached and
  # the guard under test is never exercised.
  shorter <- .vec_case(n = 20L)
  short <- aci_filter(shorter$x, shorter$comp, dt = 0.01,
                      mu0 = rep(0, 2L), R0 = diag(0.1, 2L))
  expect_error(
    aci_cir(v$x, v$comp, dt = 0.01, window = 5:10, filt = short),
    "same length"
  )
})
