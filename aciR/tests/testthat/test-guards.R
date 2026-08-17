# The guards that reject a malformed model before it reaches the arithmetic.
#
# These matter more than their line count suggests. Every one of them exists
# because the alternative is not an error but a plausible-looking number: a
# missing parameter silently recycled, a direction string quietly ignored, a
# noise covariance of zero inverted into an infinity that propagates. A guard
# that has never been executed is a guard nobody has checked, and the coverage
# gate reported this whole family as unexecuted.

# ---- the predator-prey construction ------------------------------------------

.pp_params <- function() aci_predprey_model()$parameters
.pp_signal <- function() {
  read.csv(
    system.file("extdata", "predprey_signal.csv", package = "aciR")
  )$prey
}

test_that("the parameter list is rejected when it is not one", {
  expect_error(
    aci_predprey_components(.pp_signal(), "not a list", "predator_to_prey"),
    "named list of predator-prey parameters"
  )
})

test_that("a missing parameter is named, singular and plural", {
  one <- .pp_params()
  one$alpha <- NULL
  expect_error(
    aci_predprey_components(.pp_signal(), one, "predator_to_prey"),
    "missing the parameter `alpha`"
  )

  two <- .pp_params()
  two$alpha <- NULL
  two$beta <- NULL
  expect_error(
    aci_predprey_components(.pp_signal(), two, "predator_to_prey"),
    "missing the parameters"
  )
})

test_that("the causal direction must be one of the two the model expresses", {
  expect_error(
    aci_predprey_components(.pp_signal(), .pp_params(), "sideways"),
    "must be either"
  )
  # The message says which process is observed in each case, because the
  # direction names the causal claim and not the observed series, and that is
  # the pairing most easily inverted.
  expect_error(
    aci_predprey_components(.pp_signal(), .pp_params(), "sideways"),
    "prey is observed and the predator is latent"
  )
})

test_that("the observed process must carry noise, either direction", {
  # The observed process's noise covariance is what the filter inverts. Zero
  # there is not a degenerate case to handle, it is an unanswerable question,
  # and it is refused per direction because which series is observed changes.
  expect_error(
    aci_predprey_model(direction = "predator_to_prey", sigma_y = 0),
    "`sigma_y` must be non-zero"
  )
  expect_error(
    aci_predprey_model(direction = "prey_to_predator", sigma_x = 0),
    "`sigma_x` must be non-zero"
  )
  # The other noise may be zero: it belongs to the latent process, which is
  # never inverted.
  expect_no_error(
    aci_predprey_model(direction = "predator_to_prey", sigma_x = 0)
  )
  expect_no_error(
    aci_predprey_model(direction = "prey_to_predator", sigma_y = 0)
  )
})

# ---- the conditional construction --------------------------------------------

.cond_comp <- function(n = 20L, n_x = 2L) {
  list(
    L_x = array(rep(diag(-1, n_x, n_x), n), c(n_x, n_x, n)),
    f_x = matrix(0, n_x, n),
    L_y = array(rep(diag(-1, n_x, n_x), n), c(n_x, n_x, n)),
    f_y = matrix(0, n_x, n),
    S_xoS_x = diag(0.1, n_x),
    S_yoS_y = diag(0.1, n_x),
    S_yoS_x = matrix(0, n_x, n_x)
  )
}

test_that("the conditional question needs a matrix Grammian", {
  comp <- .cond_comp()
  comp$S_xoS_x <- 0.1
  expect_error(aci_conditional(comp, 1L), "must be a matrix")
})

test_that("a character target is matched against the Grammian's row names", {
  comp <- .cond_comp()
  dimnames(comp$S_xoS_x) <- list(c("first", "second"), c("first", "second"))
  by_name <- aci_conditional(comp, "first")
  by_index <- aci_conditional(comp, 1L)
  # Naming a channel and indexing it must give the same object, or the
  # documented character path is a second implementation of the first.
  expect_equal(by_name$S_xoS_x_inv, by_index$S_xoS_x_inv)
})

test_that("a character target naming nothing is refused", {
  comp <- .cond_comp()
  dimnames(comp$S_xoS_x) <- list(c("first", "second"), c("first", "second"))
  expect_error(aci_conditional(comp, "third"), "names no observed component")

  # Without row names there is nothing to match against, which is the state
  # aci_enso_components() leaves the Grammian in. Recorded in
  # walkthrough_debt.yaml; the guard is what makes it a clear error rather
  # than a silently wrong channel.
  unnamed <- .cond_comp()
  expect_error(aci_conditional(unnamed, "first"), "names no observed component")
})

test_that("a target with no observation noise is refused", {
  comp <- .cond_comp()
  comp$S_xoS_x <- diag(c(0, 0.1))
  # The reference implementation adds a small artificial noise here. That is a
  # modelling decision, so this package reports it rather than making it
  # silently on the user's behalf.
  expect_error(aci_conditional(comp, 1L), "not positive definite at step")
})

test_that("the scalar schema admits no conditional question", {
  # One observed component cannot be asked what it says given the others,
  # because there are no others.
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 60L, seed = 2L)
  comp <- aci_dyad_components(sim$x, model$parameters)
  expect_error(aci_conditional(comp, 1L))
})
