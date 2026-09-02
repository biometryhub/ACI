## acir reserve file
## Origin: aci/tests/testthat/test-11-model-parity.R:21-155
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Cross-family parity blocks (tipping-triad F, pathways/topographic P, L84 E/F, L96 E); filed with the EnKBS majority, see DISPOSITIONS.md.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: tipping-triad (F), pathways + spectral topographic (P), Lorenz-84 (E/F) and Lorenz-96 (E) parity blocks.

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
