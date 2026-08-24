test_that("X1: detect_events recovers planted peaks; separation merges", {
  t <- seq(0, 20, by = 0.01)
  x <- sin(2 * pi * t / 5) * 0.1
  x[t > 4.9 & t < 5.1] <- x[t > 4.9 & t < 5.1] + 3      # spike at t = 5
  x[t > 5.2 & t < 5.35] <- x[t > 5.2 & t < 5.35] + 2.5  # satellite at 5.27
  x[t > 12 & t < 12.2] <- x[t > 12 & t < 12.2] - 3      # negative spike
  ob <- observed_trajectory(t, matrix(x, ncol = 1))
  ev <- detect_events(ob, q = 0.99, two_sided = TRUE, min_separation = 1,
                      event_unit = "threshold_run")
  expect_equal(nrow(ev), 2)
  expect_lt(abs(ev$t_star[1] - 5), 0.15)
  expect_lt(abs(ev$t_star[2] - 12.1), 0.15)
  expect_equal(ev$sign, c(1, -1))
  ev2 <- detect_events(ob, q = 0.99, two_sided = TRUE, min_separation = 0,
                       event_unit = "threshold_run")
  expect_gte(nrow(ev2), 3)                               # satellite kept
})

test_that("event thresholds are validated and can use the event-peak dyad", {
  t <- 0:11
  x <- c(0.5, 1, 0.5, -0.5, -2, -0.5,
         0.5, 3, 0.5, -0.5, -4, -0.5)
  ob <- observed_trajectory(t, matrix(x, ncol = 1))

  expect_error(detect_events(ob, threshold = "quantile", q = 0.5,
                             two_sided = TRUE),
               "q > 0.5", class = "aci_error_dims")
  expect_error(detect_events(ob, threshold = "absolute", value = -1),
               "positive magnitude", class = "aci_error_dims")
  expect_error(detect_events(ob, threshold = "event_peak_quantile", q = 0.5,
                             event_unit = "threshold_run"),
               "sign_excursion", class = "aci_error_dims")

  # The quantile is taken over one peak per sign excursion: amplitudes
  # (1, 2, 3, 4), not over all twelve raw observations.
  ev <- detect_events(ob, threshold = "event_peak_quantile", q = 0.5,
                      two_sided = TRUE, event_unit = "sign_excursion")
  expect_equal(ev$peak_value, c(3, -4))
  expect_equal(ev$peak_index, c(8, 11))
  expect_equal(unname(attr(ev, "threshold")), c(2.5, -2.5))
})

test_that("X2: pathways pipeline -- events, onset, App-A features, clustering", {
  m <- model_pathways()
  s <- simulate(m, seed = 12, T = 40, dt = 5e-3, burn_in = 4)
  ev <- detect_events(s$obs, q = 0.94, two_sided = TRUE,
                      min_separation = 2, burn_in = 2)
  expect_gt(nrow(ev), 5)
  a <- suppressWarnings(aci(m, s$obs))
  evo <- event_onset(a, ev, kappa = 0.8, T_pre = 1.5)
  expect_true(all(is.na(evo$onset) |
                  (evo$onset >= evo$t_star - 1.5 & evo$onset <= evo$t_star)))
  complete <- ev$t_star >= min(s$obs$t) + 1.5 &
              ev$t_star <= max(s$obs$t)
  ev_features <- ev[complete, , drop = FALSE]
  expect_gt(nrow(ev_features), 0)
  F1 <- features_pathways(s, ev_features, source = "hidden")
  expect_equal(dim(F1), c(nrow(ev_features), 18))
  expect_true(all(is.finite(F1)))
  expect_true(all(abs(F1[, "rel_damping"] + F1[, "rel_forcing"] - 1) < 1e-12))
  F2 <- features_pathways(s, ev_features, source = "smoother", paths = a$paths)
  expect_equal(dim(F2), dim(F1))
  if (nrow(ev_features) >= 6) {
    cl <- classify_events(F1, k = min(3, nrow(ev_features) - 1))
    expect_equal(length(cl$cluster), nrow(ev_features))
  }
})

test_that("X3: influence windows and sensitive directions are well-formed", {
  m <- model_pathways()
  s <- simulate(m, seed = 4, T = 5, dt = 5e-3, burn_in = 2)
  ev <- detect_events(s$obs, q = 0.9, min_separation = 2, burn_in = 1)
  ev <- ev[seq_len(min(2, nrow(ev))), ]
  a <- suppressWarnings(aci(m, s$obs))
  tab <- lag_table(m, s$obs, mode = "forward", tol = 2e-3,
                   filter = a$paths$filter)
  infl <- suppressWarnings(event_influence(m, s$obs, ev, forward_table = tab))
  expect_true(all(is.na(infl$tau_forward) | infl$tau_forward >= 0))
  expect_true(all(is.na(infl$tau_backward) | infl$tau_backward > 0))
  expect_true(all(is.na(infl$tau_backward) | infl$window_lo < infl$t_star))
  # This test exercises the direction output; the separate cohort test below
  # covers incomplete-window dropping explicitly.
  es <- event_stats(a$paths, ev, w_pre = 0, w_post = 0)
  sd1 <- sensitive_directions(es, at = 0)
  expect_equal(dim(sd1$covariance_directions), c(2, 2))
  expect_true(all(diff(sd1$covariance_score) <= 1e-12))
  expect_equal(sum(sd1$covariance_directions[, 1]^2), 1, tolerance = 1e-9)
})

test_that("event influence rejects a lag table from different observations", {
  m <- model_dyad()
  s1 <- simulate(m, seed = 31, T = 1, dt = 0.01)
  s2 <- simulate(m, seed = 32, T = 1, dt = 0.01)
  tab <- suppressWarnings(lag_table(m, s1$obs, mode = "forward",
                                    tol = 1e-2))
  ev <- data.frame(t_star = s2$obs$t[51])
  expect_error(event_influence(m, s2$obs, ev, forward_table = tab),
               "different observation values", class = "aci_error_dims")
})

test_that("event influence has explicit, strict lag and CIR controls", {
  m <- model_dyad()
  s <- simulate(m, seed = 33, T = 0.5, dt = 0.01)
  ev <- data.frame(t_star = s$obs$t[31])
  tab <- suppressWarnings(lag_table(m, s$obs, mode = "forward", tol = 1e-2))
  expect_error(event_influence(m, s$obs, ev, forward_table = tab, typo = TRUE),
               class = "aci_error_dims")
  expect_error(event_influence(m, s$obs, ev, forward_table = tab, tol = 1e-2),
               "cannot be supplied", class = "aci_error_dims")
  got <- suppressWarnings(event_influence(m, s$obs, ev, tol = 1e-2,
                                          min_M = 0))
  expect_equal(nrow(got), 1L)
})

test_that("event_stats uses one fixed complete event cohort", {
  tt <- 0:4
  mu <- cbind(tt, 10 + 2 * tt)
  cv <- array(0, c(2, 2, length(tt)))
  for (j in seq_along(tt)) cv[, , j] <- diag(c(1, 2))
  sm <- new_da_path(tt, mu, cv, "smoother")
  ev <- data.frame(t_star = c(0, 2, 4))

  es <- event_stats(sm, ev, w_pre = 1, w_post = 1)
  expect_equal(es$n_events, 1L)
  expect_equal(es$n_events_dropped, 2L)
  expect_equal(es$event_indices_used, 2L)
  expect_equal(es$rel_t, -1:1)
  expect_equal(unname(es$conditional_mean),
               unname(mu[2:4, , drop = FALSE]))
})

test_that("event mixtures include posterior covariance and mean direction", {
  tt <- 0:3
  mu <- cbind(c(0, -1, 1, 0), c(0, 0, 2, 0))
  cv <- array(0, c(2, 2, 4))
  for (j in 1:4) cv[, , j] <- diag(c(4, 1))
  sm <- new_da_path(tt, mu, cv, "smoother")
  ev <- structure(data.frame(t_star = c(1, 2)),
                  class = c("event_set", "data.frame"))
  es <- event_stats(sm, ev, w_pre = 0, w_post = 0)
  # First coordinate: within-posterior variance 4 plus variance of means 1.
  expect_equal(es$conditional_cov[1, 1, 1], 5)
  sd <- sensitive_directions(es, at = 0)
  expect_equal(sum(sd$mean_direction^2), 1, tolerance = 1e-12)
  expect_gt(sd$mean_direction[2], 0)
})

test_that("full sensitive direction maximizes the full projected KL", {
  R0 <- diag(c(2, 0.5))
  RE <- matrix(c(5, 0.8, 0.8, 1), 2, 2)
  mu0 <- c(0, 0)
  muE <- c(0.7, -1.3)
  es <- structure(list(
    rel_t = 0,
    climatology_mean = mu0,
    climatology_cov = R0,
    conditional_mean = matrix(muE, nrow = 1),
    conditional_cov = array(RE, c(2, 2, 1))
  ), class = "aci_event_stats")

  sd <- sensitive_directions(es, at = 0)
  theta <- seq(0, pi, length.out = 20001)
  dirs <- rbind(cos(theta), sin(theta))
  brute <- vapply(seq_along(theta), function(j)
    projected_kl(dirs[, j], mu0, R0, muE, RE), numeric(1))
  best <- dirs[, which.max(brute)]

  expect_equal(sd$full_projected_kl, max(brute), tolerance = 2e-7)
  expect_gt(abs(drop(crossprod(sd$full_direction, best))), 1 - 1e-6)
  expect_equal(drop(projected_kl(sd$full_direction, mu0, R0, muE, RE)),
               sd$full_projected_kl, tolerance = 1e-12)
  expect_equal(drop(sd$directions), sd$full_direction)
})

test_that("feature extraction validates complete windows and event source grid", {
  m <- model_pathways()
  s <- simulate(m, seed = 17, T = 4, dt = 0.01, burn_in = 1)
  boundary <- data.frame(peak_index = 1L, t_star = s$obs$t[1],
                         peak_value = s$obs$x[1, 1], sign = 1,
                         duration = 0)
  expect_error(features_pathways(s, boundary),
               "complete feature pre/post window", class = "aci_error_dims")

  j <- which.min(abs(s$obs$t - 2))
  valid <- data.frame(peak_index = j, t_star = s$obs$t[j],
                      peak_value = s$obs$x[j, 1],
                      sign = sign(s$obs$x[j, 1]), duration = 0)
  foreign <- valid
  foreign$peak_index <- foreign$peak_index + 1L
  expect_error(features_pathways(s, foreign),
               "do not match this simulation grid", class = "aci_error_dims")
  expect_equal(dim(features_pathways(s, valid)), c(1, 18))

  no_truth <- s; no_truth$hidden <- NULL
  expect_error(features_pathways(no_truth, valid, source = "hidden"),
               class = "aci_error_dims")
  a <- suppressWarnings(aci(m, s$obs))
  expect_equal(dim(features_pathways(s, valid, source = "smoother",
                                     paths = a$paths)), c(1, 18))
  other <- simulate(m, seed = 18, T = 4, dt = 0.01, burn_in = 1)
  foreign_paths <- suppressWarnings(aci(m, other$obs))$paths
  expect_error(features_pathways(s, valid, source = "smoother",
                                 paths = foreign_paths),
               class = "aci_error_dims")
})

test_that("event diagnostics validate time spans and clustering controls", {
  sm <- new_da_path(0:2, matrix(0, 3, 1), array(1, c(1, 1, 3)), "smoother")
  expect_error(event_stats(sm, data.frame(t_star = 4), w_pre = 0, w_post = 0),
               class = "aci_error_dims")
  ar <- structure(list(t = 0:2, aci = c(0, 1, 0)), class = "aci_result")
  expect_error(event_onset(ar, data.frame(t_star = 4)), class = "aci_error_dims")
  X <- cbind(a = c(0, 1, 2), b = c(2, 0, 1))
  expect_error(classify_events(X, k = 2, standardize = NA),
               class = "aci_error_dims")
  expect_error(classify_events(X, k = 2, seed = 1.5),
               class = "aci_error_dims")
})

test_that("X4: topographic App-B features compute on a short honest run", {
  m <- model_topographic()
  s <- simulate(m, seed = 2, T = 20, dt = 5e-3, burn_in = 4)
  ev <- detect_events(s$obs, q = 0.9, min_separation = 2, burn_in = 2)
  complete <- ev$t_star >= min(s$obs$t) + 1 &
              ev$t_star <= max(s$obs$t) - 1
  ev_features <- ev[complete, , drop = FALSE]
  skip_if(nrow(ev_features) < 2, "too few complete-window events in short run")
  F3 <- features_topographic(s, ev_features, source = "hidden")
  expect_equal(dim(F3), c(nrow(ev_features), 29))
  expect_true(all(is.finite(F3[, c("V_peak", "E_tot_peak", "R2")])))
  expect_true(all(F3[, "E_tot_pre_max"] >= F3[, "E_tot_pre_mean"] - 1e-12))
})
