test_that("benchmark dyad default initial condition matches the MATLAB source", {
  dyad <- aci_dyad_model()
  ic <- dyad$meta$ic_default
  expect_equal(ic$x0, 1)
  expect_equal(ic$y0, 2)
  expect_equal(dyad$meta$vars, list(observed = "x", hidden = "y"))
  expect_equal(dyad$meta$source_status, "paper + MATLAB checked")
  expect_error(aci_dyad_model(params = c(as.list(dyad$meta$params), typo = 1)),
               class = "aci_error_model_contract")
})

test_that("ENSO aci_code fixed-state coefficients match the MATLAB equations", {
  t <- 0.8
  state <- c(u = 0.2, hW = -0.1, TC = 0.3, TE = -0.2, tau = 0.4, I = 1.7)
  m <- aci_enso_model(hidden = "u", variant = "aci_code")
  x <- state[m$meta$vars$observed]
  y <- state[m$meta$vars$hidden]

  got <- numeric(6)
  names(got) <- m$meta$vars$all
  got[m$meta$vars$observed] <- m$f(t, x, y)
  got[m$meta$vars$hidden] <- m$g(t, x, y)

  factor <- 0.65
  gamma <- 0.75 * factor
  r <- 0.25 * factor
  alpha2 <- 0.125 * factor
  alpha1 <- alpha2 / 2 * factor
  b0 <- 2.5
  mu <- 0.5
  rC <- gamma * b0 * mu / 2
  rE <- 3 * gamma * b0 * mu / 2
  zetaC <- gamma * b0 * mu / 2
  zetaE <- gamma * b0 * mu / 2
  c1 <- (25 * (state[["TC"]] + 0.75 / 7.5)^2 + 0.9) *
    (1 + 0.3 * sin(t * 2 * pi / 6 - pi / 6)) * factor
  c2 <- 1.4 * factor *
    (1 + 0.3 * sin(t * 2 * pi / 6 + 2 * pi / 6) +
       0.25 * sin(2 * t * 2 * pi / 6 + 2 * pi / 6))
  beta <- (1 + (1 - state[["I"]] / 5)) * 0.15 * sqrt(factor)
  expected <- c(
    u = -r * state[["u"]] - alpha1 * b0 * mu *
      (state[["TC"]] + state[["TE"]]) / 2 - 0.2 * beta * state[["tau"]],
    hW = -r * state[["hW"]] - alpha2 * b0 * mu *
      (state[["TC"]] + state[["TE"]]) / 2 - 0.4 * beta * state[["tau"]],
    TC = (rC - c1) * state[["TC"]] + zetaC * state[["TE"]] +
      gamma * state[["hW"]] + state[["I"]] / 5 * factor * state[["u"]] +
      0.03 * factor + 0.8 * beta * state[["tau"]],
    TE = (rE - c2) * state[["TE"]] - zetaE * state[["TC"]] +
      gamma * state[["hW"]] + beta * state[["tau"]],
    tau = -2 * state[["tau"]],
    I = -(2 / 60) * (state[["I"]] - 2)
  )
  expect_equal(got, expected, tolerance = 2e-12)

  grams <- cgns_grams(m, t, x)
  sd_expected <- c(
    hW = 0.02 * sqrt(factor),
    TC = 0.04 * sqrt(factor),
    TE = sqrt(5) * 1e-2 * sqrt(factor),
    tau = 0.9 * (tanh(7.5 * state[["TC"]]) + 1) *
      (1 + 0.3 * cos(t * 2 * pi / 6 + 2 * pi / 6)),
    I = sqrt((2 / 60) * state[["I"]] * (4 - state[["I"]]) +
               1e-3 * (2 / 60))
  )
  expect_equal(sqrt(diag(grams$gxx)), unname(sd_expected), tolerance = 2e-12)
  expect_equal(sqrt(drop(grams$gyy)), 0.04 * sqrt(factor), tolerance = 2e-12)

  expect_false(m$meta$matlab_simulator_parity)
  expect_match(m$meta$simulation_convention, "Euler-Maruyama", fixed = TRUE)
  expect_equal(m$meta$numerical_regularization$I_variance_floor,
               1e-3 * (2 / 60))
  aci_ic <- c(u = 6.9136e-04, hW = -0.0028, TC = 0.0039, TE = 0.0051,
              tau = -0.0256, I = 1.5841)
  expect_equal(m$meta$ic_default$x0, aci_ic[m$meta$vars$observed])
  expect_equal(m$meta$ic_default$y0, aci_ic[m$meta$vars$hidden])
  expect_error(aci_enso_model(hidden = "TC", variant = "aci_code"),
               class = "aci_error_model_contract")
})

test_that("ENSO partitions are unique and source queries are explicit", {
  expect_error(aci_enso_model(hidden = c("u", "u")),
               class = "aci_error_model_contract")
  u <- aci_enso_model(hidden = "u", variant = "aci_code")
  expect_equal(u$meta$source_partition, "u")
  # The recorded estimand is the reference script's, not the set of observed
  # channels the hidden variable happens to reach.
  # ENSO_model_cond_ACI_u_unobs.m:1205 leaves only S_xoS_x_inv(1,1,:)
  # uncommented, so the target is T_C alone and everything else observed is
  # conditioned upon.  Recording the reachable set instead named a complement
  # that is structurally inert: see test-24.
  expect_equal(u$meta$vars$observed[u$meta$target_obs_idx], "TC")
  expect_equal(u$meta$vars$observed[u$meta$conditioning_obs_idx],
               c("hW", "TE", "tau", "I"))
  expect_equal(u$meta$causal_link, "(u) -> (TC) | (hW,TE,tau,I)")
  expect_match(u$meta$estimand_provenance, "ENSO_model_cond_ACI_u_unobs.m:1205",
               fixed = TRUE)
  expect_equal(aci_enso_model()$meta$source_partition, "package_only_partition")
})
