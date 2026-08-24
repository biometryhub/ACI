test_that("benchmark dyad default initial condition matches the MATLAB source", {
  dyad <- model_dyad()
  ic <- dyad$meta$ic_default
  expect_equal(ic$x0, 1)
  expect_equal(ic$y0, 2)
  expect_equal(dyad$meta$vars, list(observed = "x", hidden = "y"))
  expect_equal(dyad$meta$source_status, "paper + MATLAB checked")
  reverse_ic <- model_dyad(observe = "y")$meta$ic_default
  expect_equal(reverse_ic$x0, 2)
  expect_equal(reverse_ic$y0, 1)
  expect_match(model_dyad(observe = "y")$meta$source_status,
               "package extension")
  expect_equal(model_dyad(variant = "p3", observe = "y")$meta$source_status,
               "paper + MATLAB checked (published EnKBS dyad experiment)")
  expect_match(model_dyad(variant = "p3")$meta$provenance,
               "EnKBS causal-inference")
  expect_error(model_dyad(params = c(as.list(dyad$meta$params), typo = 1)),
               class = "aci_error_model_contract")
})

test_that("tipping-triad partitions reproduce both conditional FBCIR questions", {
  eps <- 0.1
  state <- c(x = 0.7, y = -0.2, gamma = 1.3)
  expected <- c(
    x = state[["x"]] - (1 / 3) * state[["x"]]^3 - 4 * state[["y"]],
    y = state[["gamma"]] * state[["x"]] - 0.8 -
      (0.2 / eps) * state[["y"]],
    gamma = -0.5 * (state[["gamma"]] - 1)
  )

  joint <- model_tipping_triad(eps = eps, partition = "joint")
  expect_equal(unname(joint$f(0, state["x"], state[c("y", "gamma")])),
               unname(expected["x"]))
  expect_equal(unname(joint$g(0, state["x"], state[c("y", "gamma")])),
               unname(expected[c("y", "gamma")]))
  expect_equal(joint$meta$vars$hidden, c("y", "gamma"))
  expect_match(joint$meta$source_status, "paper checked")
  expect_length(joint$meta$source_files, 0)
  expect_equal(unname(joint$meta$ic_default$x0), 0)
  expect_equal(unname(joint$meta$ic_default$y0), c(0, 0))

  y_cond <- model_tipping_triad(eps = eps, partition = "y_given_gamma")
  expect_equal(unname(y_cond$f(0, state[c("x", "gamma")], state["y"])),
               unname(expected[c("x", "gamma")]))
  expect_equal(unname(y_cond$g(0, state[c("x", "gamma")], state["y"])),
               unname(expected["y"]))
  expect_equal(y_cond$meta$conditioning,
               list(target = "x", conditioned_on = "gamma"))
  expect_equal(y_cond$meta$target_obs_idx, 1L)
  expect_equal(y_cond$meta$conditioning_obs_idx, 2L)
  expect_equal(y_cond$meta$source_files,
               "FBCIR_code-main/climate_tipping_y_bifurcation_driven.m")
  expect_equal(diag(y_cond$Sx1(0, state[c("x", "gamma")])), c(0.2, 2))
  expect_equal(drop(y_cond$Sy2(0, state[c("x", "gamma")])), 0.3 / sqrt(eps))
  expect_equal(unname(y_cond$meta$ic_default$x0), c(0, 0))

  gamma_cond <- model_tipping_triad(eps = eps, partition = "gamma_given_y")
  expect_equal(unname(gamma_cond$f(0, state[c("x", "y")], state["gamma"])),
               unname(expected[c("x", "y")]))
  expect_equal(unname(gamma_cond$g(0, state[c("x", "y")], state["gamma"])),
               unname(expected["gamma"]))
  expect_equal(gamma_cond$meta$conditioning,
               list(target = "y", conditioned_on = "x"))
  expect_equal(gamma_cond$meta$target_obs_idx, 2L)
  expect_equal(gamma_cond$meta$conditioning_obs_idx, 1L)
  expect_equal(gamma_cond$meta$source_status, "paper + MATLAB checked")
  expect_equal(diag(gamma_cond$Sx1(0, state[c("x", "y")])),
               c(0.2, 0.3 / sqrt(eps)))
  expect_equal(drop(gamma_cond$Sy2(0, state[c("x", "y")])), 2)
  expect_equal(unname(gamma_cond$meta$ic_default$x0), c(0, 0))
  expect_error(model_tipping_triad(params = c(
    as.list(joint$meta$params[names(joint$meta$params) != "eps"]), typo = 1)),
    class = "aci_error_model_contract")
})

test_that("pathways and spectral-topographic provenance is explicit", {
  pathways <- model_pathways()
  expect_equal(pathways$meta$vars,
               list(observed = "u", hidden = c("gamma", "b")))
  expect_match(pathways$meta$provenance, "moser2026extremes equations")
  expect_match(pathways$meta$source_status, "paper checked")
  expect_error(model_pathways(params = c(as.list(pathways$meta$params), typo = 1)),
               class = "aci_error_model_contract")

  topo <- model_topographic()
  expect_equal(topo$meta$vars$observed, "V")
  expect_length(topo$meta$vars$hidden, topo$l)
  expect_match(topo$meta$source_status, "open discrepancy")
  expect_match(topo$meta$preset_caveat, "sigma_psi/sqrt(2)", fixed = TRUE)
  topo_input <- topo$meta$params[names(topo$meta$params) != "beta"]
  expect_error(model_topographic(params = c(topo_input, beta = 1)),
               class = "aci_error_model_contract")
})

test_that("Lorenz-84 separates the P3 and seasonal FBCIR presets", {
  yz <- c(0.4, -0.3)
  p3 <- model_l84()
  fb <- model_l84(variant = "fbcir")

  expect_equal(drop(p3$fy(0, yz)), 0.25 * 8 - sum(yz^2))
  expect_equal(drop(p3$fy(73 / 2, yz)), drop(p3$fy(0, yz)), tolerance = 1e-14)
  expect_equal(diag(p3$Sx1(0, yz)), rep(0.1, 2))
  expect_equal(drop(p3$Sy2(0, yz)), 0.1)

  expect_equal(drop(fb$fy(0, yz)), 0.25 * (8 + 3) - sum(yz^2))
  expect_equal(drop(fb$fy(73 / 2, yz)), 0.25 * (8 - 3) - sum(yz^2),
               tolerance = 1e-14)
  expect_equal(diag(fb$Sx1(0, yz)), rep(0.2, 2))
  expect_equal(drop(fb$Sy2(0, yz)), 0.2)
  expect_equal(fb$meta$forcing$F_period, 73)
  expect_equal(unname(fb$meta$ic_default$x0), c(0, 1))
  expect_equal(unname(fb$meta$ic_default$y0), 1)
  expect_equal(model_l84(variant = "fbcir", target = "z")$meta$causal_link,
               "x -> z | y")
})

test_that("Lorenz-96 P3 defaults and explicit legacy preset are distinct", {
  p3 <- model_l96(n = 8)
  expect_equal(p3$meta$coords$obs, c(2L, 4L, 6L, 8L))
  expect_equal(p3$meta$coords$hidden, c(1L, 3L, 5L, 7L))
  expect_equal(diag(p3$Sx(0, rep(0, 4))), rep(sqrt(0.1), 4))
  expect_equal(diag(p3$Sy(0, rep(0, 4), rep(0, 4))), rep(sqrt(5), 4))
  expect_match(p3$meta$source_status, "non-benchmark configuration")
  expect_equal(p3$meta$vars$observed, paste0("x", c(2, 4, 6, 8)))

  legacy <- model_l96(n = 8, preset = "legacy")
  expect_equal(legacy$meta$coords$obs, c(1L, 3L, 5L, 7L))
  expect_equal(diag(legacy$Sx(0, rep(0, 4))), rep(1, 4))
  expect_equal(diag(legacy$Sy(0, rep(0, 4), rep(0, 4))), rep(1, 4))
  expect_match(legacy$meta$source_status, "package extension")
  expect_match(legacy$meta$preset_caveat, "not jiang2026enkbs", fixed = TRUE)

  exact <- model_l96()
  expect_match(exact$meta$source_status, "paper checked")
  expect_null(exact$meta$preset_caveat)

  explicit <- model_l96(n = 8, observe = c(1, 4), sigma = 0.7)
  expect_equal(explicit$meta$coords$obs, c(1L, 4L))
  expect_equal(diag(explicit$Sx(0, rep(0, 2))), rep(0.7, 2))
  expect_equal(diag(explicit$Sy(0, rep(0, 2), rep(0, 6))), rep(0.7, 6))

  expect_error(model_l96(n = 3), class = "aci_error_model_contract")
  expect_error(model_l96(n = 8.5), class = "aci_error_model_contract")
  expect_error(model_l96(n = NA_real_), class = "aci_error_model_contract")
  expect_error(model_l96(n = 8, F = c(8, 9)), class = "aci_error_model_contract")
  expect_error(model_l96(n = 8, F = Inf), class = "aci_error_model_contract")
  expect_error(model_l96(n = 8, observe = c(1, 2.5)),
               class = "aci_error_model_contract")
  expect_error(model_l96(n = 8, observe = c(0, 2)),
               class = "aci_error_model_contract")
  expect_error(model_l96(n = 8, observe = "even"),
               class = "aci_error_model_contract")
  expect_error(model_l96(n = 8, sigma_observed = "small"),
               class = "aci_error_model_contract")
})

test_that("ENSO aci_code fixed-state coefficients match the MATLAB equations", {
  t <- 0.8
  state <- c(u = 0.2, hW = -0.1, TC = 0.3, TE = -0.2, tau = 0.4, I = 1.7)
  m <- model_enso6(hidden = "u", variant = "aci_code")
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

  cfy <- model_enso6(hidden = "u", variant = "cfy22")
  cfy_ic <- c(u = 0, hW = 0, TC = 0.1, TE = 0.1, tau = 0, I = 2)
  expect_equal(cfy$meta$ic_default$x0, cfy_ic[cfy$meta$vars$observed])
  expect_equal(cfy$meta$ic_default$y0, cfy_ic[cfy$meta$vars$hidden])
  expect_error(model_enso6(hidden = "TC", variant = "aci_code"),
               class = "aci_error_model_contract")
})

test_that("ENSO partitions are unique and source queries are explicit", {
  expect_error(model_enso6(hidden = c("u", "u")),
               class = "aci_error_model_contract")
  u <- model_enso6(hidden = "u", variant = "aci_code")
  expect_equal(u$meta$source_partition, "u")
  expect_equal(u$meta$vars$observed[u$meta$target_obs_idx],
               c("TC", "TE", "I"))
  expect_equal(u$meta$vars$observed[u$meta$conditioning_obs_idx],
               c("hW", "tau"))
  expect_equal(model_enso6()$meta$source_partition, "package_only_partition")
})
