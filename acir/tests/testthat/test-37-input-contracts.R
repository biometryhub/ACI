# Input-shape contracts, checked before the arithmetic ------------------------
#
# Argument-shape contracts that were reached after the arithmetic they are a
# precondition for: shared noise-channel widths, and the nsub integer bound.

quietly_ic <- function(expr) withCallingHandlers(
  expr, warning = function(w) invokeRestart("muffleWarning"))

capture_conditions <- function(expr) {
  warnings <- character(0)
  err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) { err <<- e; NULL }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, error = err, warnings = warnings)
}


test_that(paste0("mismatched shared noise-channel widths report the ",
                 "contract they violate"), {
  ## Sx1 is 2 by 3 and Sy1 is 2 by 2: the first shared channel pair disagrees.
  e <- expect_error(
    aci_model(Lx = function(t, x) diag(2), fx = function(t, x) c(0, 0),
              Ly = function(t, x) diag(-0.5, 2), fy = function(t, x) c(0, 0),
              Sx1 = function(t, x) matrix(1, 2, 3),
              Sy1 = function(t, x) diag(2),
              Sy2 = function(t, x) diag(2), k = 2, l = 2),
    class = "aci_error_model_contract")
  expect_match(conditionMessage(e), "matching column counts")
  expect_false(grepl("non-conformable", conditionMessage(e), fixed = TRUE))

  ## The mirrored case on the second shared channel pair.
  e <- expect_error(
    aci_model(Lx = function(t, x) diag(2), fx = function(t, x) c(0, 0),
              Ly = function(t, x) diag(-0.5, 2), fy = function(t, x) c(0, 0),
              Sx1 = function(t, x) diag(2), Sy1 = function(t, x) diag(2),
              Sx2 = function(t, x) matrix(0.1, 2, 3),
              Sy2 = function(t, x) matrix(0.2, 2, 2), k = 2, l = 2),
    class = "aci_error_model_contract")
  expect_match(conditionMessage(e), "matching column counts")
  expect_false(grepl("non-conformable", conditionMessage(e), fixed = TRUE))
})


test_that(paste0("admissible rectangular shared channels still build ",
                 "and still filter"), {
  dense <- aci_model(
    Lx = function(t, x) matrix(c(0.2, 0.1, -0.05, 0.3), 2, 2),
    fx = function(t, x) c(0.1, -0.2),
    Ly = function(t, x) matrix(c(-0.4, 0.05, 0.02, -0.6), 2, 2),
    fy = function(t, x) c(0, 0.1),
    Sx1 = function(t, x) matrix(c(0.5, 0.1, 0.2, 0.4, -0.1, 0.3), 2, 3),
    Sy1 = function(t, x) matrix(c(0.2, 0.05, -0.1, 0.3, 0.15, 0.1), 2, 3),
    Sx2 = function(t, x) matrix(c(0.3, 0.1), 2, 1),
    Sy2 = function(t, x) matrix(c(0.25, 0.35), 2, 1), k = 2, l = 2)
  expect_s3_class(dense, "cgns_model")
  expect_true(dense$meta$correlated_noise)
  obd <- observed_trajectory(c(0, 0.01, 0.02),
                             cbind(c(0, 0.1, 0.2), c(0, -0.05, -0.1)))
  f <- quietly_ic(aci_filter(dense, obd,
                             init = list(mean = c(0, 0), cov = diag(2))))
  expect_equal(f$meta$loglik, 2.5830017959972604, tolerance = 1e-12)
  expect_true(all(is.finite(f$mean)))
})


test_that("a singly defective model still reports its own defect", {
  base_args <- list(
    fx = function(t, x) c(0, 0), Ly = function(t, x) diag(-0.5, 2),
    fy = function(t, x) c(0, 0), Sx1 = function(t, x) diag(2),
    Sy2 = function(t, x) diag(2), k = 2, l = 2)
  e <- expect_error(
    do.call(aci_model, c(list(Lx = function(t, x) diag(3)), base_args)),
    class = "aci_error_model_contract")
  expect_match(conditionMessage(e), "Lx\\(t,x\\) must be k x l")
  e <- expect_error(
    do.call(aci_model, c(list(Lx = function(t, x) diag(2)),
                         modifyList(base_args,
                                    list(fx = function(t, x) 0)))),
    class = "aci_error_model_contract")
  expect_match(conditionMessage(e), "fx / fy dims wrong")
  e <- expect_error(
    do.call(aci_model, c(list(Lx = function(t, x) diag(2)),
                         modifyList(base_args,
                                    list(Sx1 = function(t, x)
                                      matrix(1, 3, 2))))),
    class = "aci_error_model_contract")
  expect_match(conditionMessage(e), "incorrect row dimensions")
})


test_that(paste0("nsub outside the integer range is refused by the ",
                 "package, not by a coercion"), {
  scalar_model <- aci_dyad_model()
  scalar_obs <- as_obs(simulate(scalar_model, seed = 1, t_end = 0.2,
                                dt = 0.01))
  scalar_init <- list(mean = 0, cov = matrix(1, 1, 1))
  matrix_model <- aci_enso_model(hidden = c("u", "hW", "tau"))
  matrix_obs <- as_obs(simulate(matrix_model, seed = 1, t_end = 0.2,
                                dt = 0.01))
  matrix_init <- list(mean = rep(0, 3), cov = diag(0.1, 3))
  bad <- list(2^40, 3e9, .Machine$integer.max + 1, NA_real_, NaN, Inf, -Inf,
              -1, 0, 2.5, "2", c(1, 2), numeric(0))
  for (v in bad) {
    label <- paste(format(v), collapse = ",")
    for (route in c("scalar", "matrix")) {
      res <- capture_conditions(if (route == "scalar")
        aci_filter(scalar_model, scalar_obs, init = scalar_init, nsub = v) else
        aci_filter(matrix_model, matrix_obs, init = matrix_init, nsub = v))
      info <- paste(route, label)
      expect_s3_class(res$error, "aci_error_dims")
      expect_identical(conditionMessage(res$error),
                       "nsub must be a positive integer.", info = info)
      ## The base coercion warning is gone.
      expect_identical(res$warnings, character(0), info = info)
    }
  }
})


test_that("every supported nsub keeps its result", {
  scalar_model <- aci_dyad_model()
  scalar_obs <- as_obs(simulate(scalar_model, seed = 1, t_end = 2, dt = 0.01))
  scalar_init <- list(mean = 0, cov = matrix(1, 1, 1))
  reference <- lapply(list(1L, 2L, 4L, 20L), function(n)
    quietly_ic(aci_filter(scalar_model, scalar_obs, init = scalar_init,
                          nsub = n)))
  ## Values captured on the pre-change source at the second grid point and at
  ## the end of the record.
  expect_equal(vapply(reference, function(f) as.numeric(f$mean[2L]),
                      numeric(1)),
               c(0.059418475703066916, 0.054380863674442545,
                 0.052292507303202486, 0.050779444791134977),
               tolerance = 1e-12)
  expect_equal(vapply(reference,
                      function(f) as.numeric(f$mean[nrow(f$mean)]),
                      numeric(1)),
               c(-0.37163529584105448, -0.37021584185275108,
                 -0.3695166995368388, -0.36896220614621961),
               tolerance = 1e-12)
  ## The double 1 and the logical TRUE are the integer 1, as before.
  one <- reference[[1L]]
  for (v in list(1.0, TRUE)) {
    f <- quietly_ic(aci_filter(scalar_model, scalar_obs, init = scalar_init,
                               nsub = v))
    expect_identical(as.numeric(f$mean), as.numeric(one$mean),
                     info = format(v))
    expect_identical(f$meta$nsub, 1L, info = format(v))
  }
  matrix_model <- aci_enso_model(hidden = c("u", "hW", "tau"))
  matrix_obs <- as_obs(simulate(matrix_model, seed = 1, t_end = 0.2,
                                dt = 0.01))
  matrix_init <- list(mean = rep(0, 3), cov = diag(0.1, 3))
  fm1 <- quietly_ic(aci_filter(matrix_model, matrix_obs, init = matrix_init,
                               nsub = 1L, stepper = "implicit"))
  for (v in list(1.0, TRUE)) {
    fm <- quietly_ic(aci_filter(matrix_model, matrix_obs, init = matrix_init,
                                nsub = v, stepper = "implicit"))
    expect_identical(fm$mean, fm1$mean, info = format(v))
    expect_identical(fm$cov, fm1$cov, info = format(v))
  }
  fm4 <- quietly_ic(aci_filter(matrix_model, matrix_obs, init = matrix_init,
                               nsub = 4L, stepper = "implicit"))
  expect_identical(fm4$meta$nsub, 4L)
  expect_false(identical(fm4$mean, fm1$mean))
})
