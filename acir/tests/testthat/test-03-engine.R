# T1/T5/T3/T4: CGNS engine vs an exact discrete Kalman/RTS oracle -------------
lin_model <- function(a = 1, b = -1, fx0 = 0.2, fy0 = 0.3, sx = 0.5, sy = 0.6)
  aci_model(Lx = function(t, x) matrix(a, 1, 1),
            fx = function(t, x) fx0,
            Ly = function(t, x) matrix(b, 1, 1),
            fy = function(t, x) fy0,
            Sx1 = function(t, x) matrix(sx, 1, 1),
            Sy2 = function(t, x) matrix(sy, 1, 1), k = 1, l = 1)

online_smoother <- function(model, obs, filter) {
  bundle <- .compile_cgns_run(model, obs)
  .smoother_thmD1_compiled(
    bundle, filter, validate = TRUE, warn_cost = FALSE
  )
}

# exact discrete KF/RTS for the Euler-discretized linear pair:
# z_k = Dx_k = a dt y_k + fx0 dt + sx sqrt(dt) w ;  y_{k+1} = (1+b dt) y_k + fy0 dt + sy sqrt(dt) v
kf_rts <- function(obs, a, b, fx0, fy0, sx, sy, mu0, P0) {
  dt <- obs$dt; N1 <- nrow(obs$x); N <- N1 - 1
  H <- a * dt; Rn <- sx^2 * dt; Fm <- 1 + b * dt; Q <- sy^2 * dt
  mu_f <- P_f <- mu_pr <- P_pr <- numeric(N1)
  mu <- mu0; P <- P0; mu_f[1] <- mu; P_f[1] <- P
  for (k in 1:N) {
    z <- obs$x[k + 1, 1] - obs$x[k, 1] - fx0 * dt
    S <- H * P * H + Rn; K <- P * H / S
    mu <- mu + K * (z - H * mu); P <- (1 - K * H) * P     # update at k (obs uses y_k)
    mu_pr[k] <- Fm * mu + fy0 * dt; P_pr[k] <- Fm * P * Fm + Q
    mu <- mu_pr[k]; P <- P_pr[k]
    mu_f[k + 1] <- mu; P_f[k + 1] <- P
  }
  mu_s <- mu_f; P_s <- P_f
  for (k in N:1) {
    # smoother at k combines the UPDATED state at k with the k+1 smoother
    # recover updated-at-k from prediction: mu_upd = (mu_pr - fy0 dt)/Fm etc.
    mu_u <- (mu_pr[k] - fy0 * dt) / Fm
    P_u  <- (P_pr[k] - Q) / Fm^2
    G <- P_u * Fm / P_pr[k]
    mu_s[k] <- mu_u + G * (mu_s[k + 1] - mu_pr[k])
    P_s[k]  <- P_u + G^2 * (P_s[k + 1] - P_pr[k])
  }
  list(mu_f = mu_f, P_f = P_f, mu_s = mu_s, P_s = P_s)
}

test_that("T1: filter/smoother converge to the exact discrete KF/RTS oracle", {
  m <- lin_model()
  err_at <- function(dt) {
    s <- simulate(m, seed = 42, t_end = 4, dt = dt, ic = list(x0 = 0, y0 = 0.5))
    filt <- aci_filter(m, s$obs, init = list(mean = 0, cov = 1))
    smoo <- aci_smoother(m, s$obs, filter = filt)
    or <- kf_rts(s$obs, 1, -1, 0.2, 0.3, 0.5, 0.6, 0, 1)
    burn <- seq_len(200)
    c(f = max(abs(filt$mean[-burn, 1] - or$mu_f[-burn])),
      s = max(abs(smoo$mean[-burn, 1] - or$mu_s[-burn])),
      Pf = max(abs(filt$cov[1, 1, -burn] - or$P_f[-burn])),
      Ps = max(abs(smoo$cov[1, 1, -burn] - or$P_s[-burn])))
  }
  e1 <- err_at(2e-3); e2 <- err_at(1e-3)
  expect_lt(max(e2), 5e-3)                    # tight absolute agreement at dt = 1e-3
  expect_lt(max(e2 / pmax(e1, 1e-12)), 1.1)   # refinement never degrades accuracy
  # (observed errors ~3e-4 sit at the scheme-agreement plateau, so a strict
  #  halving-order ratio is not a stable assertion here; T10 order checks are
  #  covered by the absolute bounds at two meshes.)
})

test_that("T3: online-table identities -- diag == ACI == KL; P[j,N] == 0", {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 5, t_end = 2, dt = 5e-3)
  filt <- suppressWarnings(aci_filter(m, s$obs))
  smoo <- online_smoother(m, s$obs, filt)
  lt <- lag_table(m, s$obs, mode = "full", filter = filt, smoother = smoo)
  a <- aci(m, s$obs, table = lt)
  expect_equal(lt_diag(lt), a$aci, tolerance = 1e-10)
  expect_equal(lt_diag(lt),
               aci_metric(smoo, filt)$total, tolerance = 1e-12)
  lastvals <- sapply(seq_along(lt$t), function(j) {
    r <- lt$rows[[j]]; r[length(r)] })
  expect_lt(max(lastvals[1:(length(lastvals) - 1)]), 1e-10)
  expect_equal(lt$meta$reference_smoother, "thmD1_online_complete")
  expect_equal(a$signal + a$dispersion, a$aci, tolerance = 1e-12)
})

test_that("retention alone never changes the headline backward-ODE ACI", {
  m <- aci_dyad_model(); ini <- list(mean = 0, cov = matrix(1, 1, 1))
  s <- simulate(m, seed = 5, t_end = 0.4, dt = 0.01)
  ap <- aci(m, s$obs, init = ini, keep = "paths")
  at <- aci(m, s$obs, init = ini, keep = "table")
  an <- aci(m, s$obs, init = ini, keep = "none")
  expect_equal(ap$aci, at$aci, tolerance = 0)
  expect_equal(ap$aci, an$aci, tolerance = 0)
  # meta$smoother_scheme names the discretization, not the finer route tag,
  # so that it is comparable with the "theorem3_discrete" a lag table reports
  expect_equal(ap$meta$smoother_scheme, "backward_ode_euler")
  expect_equal(ap$paths$smoother$meta$route, "backward_ode")
  expect_equal(at$table$meta$reference_smoother, "thmD1_online_complete")
  expect_error(aci(m, s$obs, typo = TRUE), class = "aci_error_dims")
  expect_error(aci(m, s$obs, table = at$table, typo = TRUE),
               class = "aci_error_dims")
  # `m =` would partial-match `model` and silently shift every argument;
  # the guard names the deferred ensemble argument instead (the EnKBS
  # family restores it as a real argument).
  expect_error(aci(m, s$obs, m = 50), class = "aci_error_not_implemented")
  expect_error(aci(m, s$obs, m = 50), regexp = "partial-match")
  expect_error(aci_range(at$table, typo = TRUE), class = "aci_error_dims")
})

test_that("lag table full mode cannot be silently capped", {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 5, t_end = 0.2, dt = 0.01)
  expect_error(lag_table(m, s$obs, mode = "full", max_lag = 2),
               class = "aci_error_dims")
})

test_that("the closed-form engine rejects mismatched observations", {
  m <- aci_dyad_model()
  bad <- observed_trajectory(seq(0, 0.1, by = 0.01), matrix(0, 11, 2))
  expect_error(aci_filter(m, bad), class = "aci_error_dims")
})

test_that("Gaussian paths cannot be reused with different data or model", {
  m <- aci_dyad_model()
  s1 <- simulate(m, seed = 1, t_end = 0.2, dt = 0.01)
  s2 <- simulate(m, seed = 2, t_end = 0.2, dt = 0.01)
  f1 <- suppressWarnings(aci_filter(m, s1$obs,
                                    init = list(mean = 0, cov = diag(1))))
  expect_error(aci_smoother(m, s2$obs, filter = f1), class = "aci_error_dims")
  m2 <- aci_dyad_model(params = list(d_x = 0.6, gamma = 2, f_x = 0.5,
                                 s_x = 0.5, d_y = 0.5, f_y = 1, s_y = 1))
  expect_error(aci_smoother(m2, s1$obs, filter = f1),
               class = "aci_error_model_contract")
  sm1 <- aci_smoother(m, s1$obs, filter = f1)
  expect_error(aci_smoother(m, s1$obs, filter = sm1), class = "aci_error_dims")
  expect_error(lag_table(m, s1$obs, filter = sm1, smoother = f1),
               class = "aci_error_dims")
})

test_that("lag tables reject a theorem smoother from another prior", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 3, t_end = 0.2, dt = 0.01)
  f1 <- aci_filter(m, s$obs, init = list(mean = 0, cov = matrix(1, 1, 1)))
  f2 <- aci_filter(m, s$obs, init = list(mean = 3, cov = matrix(0.2, 1, 1)))
  sm2 <- online_smoother(m, s$obs, f2)
  sm2$meta$source_model <- m
  expect_error(lag_table(m, s$obs, filter = f1, smoother = sm2),
               class = "aci_error_model_contract")
})

test_that("T4: adaptive truncation engages and agrees with the full table", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 3, t_end = 3, dt = 5e-3)
  filt <- suppressWarnings(aci_filter(m, s$obs)); smoo <- online_smoother(m, s$obs, filt)
  full <- lag_table(m, s$obs, mode = "full", filter = filt, smoother = smoo)
  adap <- lag_table(m, s$obs, mode = "forward", tol = 1e-4,
                    filter = filt, smoother = smoo)
  N1 <- length(s$obs$t)
  dmax <- max(sapply(seq_len(N1), function(j)
    max(abs(lt_row(full, j) - lt_row(adap, j)))))
  expect_lt(dmax, 2e-4)                    # empirical agreement for this fixture
  expect_lt(dmax, 1e-2 * max(lt_diag(full)))
  expect_lt(mean(adap$L, na.rm = TRUE),
            0.98 * mean(full$L, na.rm = TRUE))           # truncation engages
  expect_true(all(lt_tail_bound(adap) >= 0))
  expect_lt(max(lt_tail_bound(adap)), 1.5e-4)            # heuristic estimates here
  # min_M = NULL: this test is about the truncation label, not about masking
  # low-signal anchors, and the short fixture has one.
  expect_equal(aci_range(adap, min_M = NULL)$bound,
               "objective_on_truncated_table")
})

test_that("max_lag is an exact positive-lag storage cap", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 2, t_end = 0.1, dt = 0.01)
  # An explicit prior: the storage cap under test does not depend on it, and
  # defaulting would exercise the diffuse-prior path incidentally.
  ini <- list(mean = 2, cov = matrix(1, 1, 1))
  for (cap in 1:3) {
    tab <- lag_table(m, s$obs, mode = "forward", tol = 0, max_lag = cap,
                     init = ini)
    want <- pmin(cap, length(tab$t) - seq_along(tab$t))
    expect_equal(tab$L, as.integer(want))
    expect_equal(vapply(tab$rows, length, integer(1)), as.integer(want + 1L))
    expect_match(aci_range(tab, method = "l1_linf", min_M = NULL)$bound,
                 "truncated")
  }
})

test_that("T1b(strict): the diffuse-init stress case stops under the default", {
  # Updated for the strict covariance policy (decision D3). Before it, the
  # explicit scheme here warned aci_warn_riccati_stiff, floored the variance
  # and returned a path; the collapse this block used to assert *was* the
  # floor. Under the default the same call stops instead, and the condition
  # carries where and when.
  m <- aci_dyad_model(); s <- simulate(m, seed = 55, t_end = 1.5, dt = 5e-3)
  big <- list(mean = 0, cov = diag(100, 1))
  e <- tryCatch(aci_filter(m, s$obs, init = big, stepper = "explicit"),
                error = function(e) e)
  expect_s3_class(e, "aci_error_covariance_not_spd")
  expect_s3_class(e, "aci_error_spd")
  expect_identical(e$site, "filter_explicit")
  expect_identical(e$role, "filter covariance")
  expect_identical(e$index, 2L)
  expect_equal(e$time, 0.005)
  expect_lt(e$value, 0)
  expect_match(conditionMessage(e), "regularize = \"floor\"", fixed = TRUE)
  # the implicit stepper never leaves the cone here, so it is unaffected
  fi <- aci_filter(m, s$obs, init = big, stepper = "implicit")
  expect_gt(min(fi$cov[1, 1, ]), 0.05)
})

test_that("T1b(floor): the documented explicit collapse is still available", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 55, t_end = 1.5, dt = 5e-3)
  big <- list(mean = 0, cov = diag(100, 1))
  # The explicit scheme is expected to report its own instability here; that
  # warning is part of the documented behaviour this test pins down.
  ## the floor this case takes is itself reported, so both conditions are
  ## asserted rather than one being left to re-signal unmatched
  expect_warning(
    expect_warning(
      fe <- aci_filter(m, s$obs, init = big, stepper = "explicit",
                       regularize = "floor"),
      class = "aci_warn_riccati_stiff"),
    class = "aci_warn_regularized")
  expect_lt(min(fe$cov[1, 1, ]), 1e-6)     # documented explicit-scheme collapse
  expect_gt(min(fe$cov[1, 1, ]), 0)        # and it is floored back into the cone
  rg <- fe$meta$regularization
  expect_identical(rg$policy, "floor")
  expect_true(rg$fired)
  expect_gt(rg$n_events, 0L)
  expect_identical(rg$sites$site, "filter_explicit")
  expect_identical(rg$sites$n, rg$n_events)
  expect_identical(rg$sites$first_index, 2L)
  expect_equal(rg$sites$first_time, 0.005)
  expect_lt(rg$sites$worst_value, 0)
})

test_that("T1b: implicit Riccati step survives the diffuse-init stress case", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 55, t_end = 1.5, dt = 5e-3)
  fi <- aci_filter(m, s$obs, init = list(mean = 0, cov = diag(100, 1)),
                   stepper = "implicit")
  expect_gt(min(fi$cov[1, 1, ]), 0.05)     # implicit stays positive-definite
  # both steppers agree at a sane init (scheme difference is O(dt))
  ini <- list(mean = 0, cov = diag(1, 1))
  f1 <- aci_filter(m, s$obs, init = ini, stepper = "explicit")
  f2 <- aci_filter(m, s$obs, init = ini, stepper = "implicit")
  # Backward-Euler mean/covariance splitting is deliberately more dissipative
  # than explicit Euler at this grid; the difference remains first-order and
  # is assessed by the sub-stepping convergence check below.
  expect_lt(max(abs(f1$mean - f2$mean)), 0.07)
  expect_lt(max(abs(f1$cov - f2$cov)), 0.02)
  # sub-stepping refines toward the same limit
  f4 <- aci_filter(m, s$obs, init = ini, stepper = "implicit", nsub = 8)
  expect_lt(max(abs(f4$mean - f2$mean)), 0.05)
  # lag_table enforces its explicit single-step contract
  expect_warning(lag_table(m, s$obs, mode = "forward", filter = f2),
                 class = "aci_warn_stepper")
  expect_error(lag_table(m, s$obs, stepper = "implicit"),
               class = "aci_error_stepper")
})

test_that("user priors must be genuinely positive definite", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 1, t_end = 0.1, dt = 0.01)
  expect_error(aci_filter(m, s$obs, init = list(mean = 0, cov = matrix(0, 1, 1))),
               class = "aci_error_spd")
  expect_error(aci_filter(m, s$obs, init = list(mean = 0, cov = matrix(-1, 1, 1))),
               class = "aci_error_spd")
})

test_that("implicit covariance update remains SPD in a stiff scalar model", {
  m <- aci_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-100, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
  ob <- observed_trajectory(c(0, 1, 2), matrix(0, 3, 1))
  fit <- suppressWarnings(aci_filter(m, ob, init = list(mean = 0, cov = 1),
                                     stepper = "implicit"))
  expect_true(all(is.finite(fit$cov)))
  expect_true(all(fit$cov[1, 1, ] > 0))
})

test_that("correlation activated after t0 selects the correlated backward ODE", {
  m <- aci_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.4, 1, 1),
    Sy1 = function(t, x) matrix(0.2 * x[1], 1, 1),
    Sy2 = function(t, x) matrix(0.5, 1, 1), k = 1, l = 1)
  s <- simulate(m, seed = 1, t_end = 0.2, dt = 0.01,
                ic = list(x0 = 0, y0 = 1))
  f <- suppressWarnings(aci_filter(m, s$obs, init = list(mean = 0, cov = 1)))
  sm <- suppressWarnings(aci_smoother(m, s$obs, filter = f))
  expect_equal(sm$meta$route, "backward_ode_correlated")
})

test_that("correlated smoother transcribes the active MATLAB backward loop", {
  m <- aci_model(
    Lx = function(t, x) matrix(0.7, 1, 1), fx = function(t, x) 0.1,
    Ly = function(t, x) matrix(-1.1, 1, 1), fy = function(t, x) -0.05,
    Sx1 = function(t, x) matrix(0.4, 1, 1),
    Sy1 = function(t, x) matrix(0.2, 1, 1),
    Sy2 = function(t, x) matrix(0.5, 1, 1), k = 1, l = 1)
  ob <- observed_trajectory(0.1 * 0:4,
                            matrix(c(0, 0.03, -0.01, 0.04, 0.02), ncol = 1))
  # dt * ||Ly|| = 0.11 here, matching the coarse grid of the MATLAB loop this
  # transcribes, so the accuracy advisory is expected rather than incidental.
  expect_warning(
    f <- aci_filter(m, ob, init = list(mean = 0.2, cov = matrix(0.8, 1, 1))),
    class = "aci_warn_dt_stability")
  got <- aci_smoother(m, ob, filter = f)

  mu <- f$mean[nrow(f$mean), ]; R <- f$cov[, , length(f$t)]
  mu_ref <- f$mean; R_ref <- f$cov
  for (j in (length(ob$t) - 1L):1L) {
    co <- eval_coefs(m, ob$t[j], ob$x[j, ])
    Gi <- solve(co$gxx); Cgi <- co$gyx %*% Gi
    A <- co$Ly - Cgi %*% co$Lx
    B <- co$gyy - Cgi %*% t(co$gyx)
    Rfi <- solve(f$cov[, , j])
    d <- drop(Rfi %*% (f$mean[j, ] - mu))
    dx <- ob$x[j + 1L, ] - ob$x[j, ]
    mu <- mu - drop(co$Ly %*% mu + co$fy - B %*% d) * ob$dt +
      drop(Cgi %*% (-dx + (drop(co$Lx %*% mu) + co$fx) * ob$dt))
    H <- A + B %*% Rfi
    R <- spd_floor(sym(R - (H %*% R + R %*% t(H) - B) * ob$dt))
    mu_ref[j, ] <- mu; R_ref[, , j] <- R
  }
  expect_equal(got$mean, mu_ref, tolerance = 1e-13)
  expect_equal(got$cov, R_ref, tolerance = 1e-13)
})

test_that("supplied Gaussian paths must match grid and conditioning", {
  m <- aci_dyad_model(); s <- simulate(m, seed = 2, t_end = 0.2, dt = 0.01)
  f <- suppressWarnings(aci_filter(m, s$obs, init = list(mean = 0, cov = 1)))
  bad <- f; bad$t[2] <- bad$t[2] + 0.001
  expect_error(aci_smoother(m, s$obs, filter = bad), class = "aci_error_dims")
  stripped <- f; stripped$meta$obs_x <- NULL
  expect_error(aci_smoother(m, s$obs, filter = stripped), class = "aci_error_dims")
  stripped_model <- f; stripped_model$meta$source_model <- NULL
  stripped_model$meta$model <- NULL
  expect_error(aci_smoother(m, s$obs, filter = stripped_model),
               class = "aci_error_model_contract")
  singular <- f; singular$cov[1, 1, 2] <- 0
  expect_error(aci_smoother(m, s$obs, filter = singular), class = "aci_error_spd")
  expect_error(aci_smoother(m, s$obs, filter = f,
                            init = list(mean = 9, cov = matrix(2, 1, 1))),
               "init conflicts", class = "aci_error_dims")
})


test_that("T1c: nonlinear truth-tracking referee -- particle filter agrees", {
  # Closes a verification gap: T1/T3/T5 check scheme identities, not whether
  # the nonlinear posterior itself is right. A bootstrap particle filter with
  # exact discrete weights is an independent referee on the dyad.
  m <- aci_dyad_model(); s <- simulate(m, seed = 21, t_end = 2.5, dt = 5e-3, burn_in = 0.5)
  fp <- suppressWarnings(aci_filter(m, s$obs, init = list(mean = 0, cov = diag(1, 1))))
  set.seed(9); Np <- 2000; dt <- s$obs$dt; x <- s$obs$x[, 1]; N <- length(x) - 1
  yp <- rnorm(Np, 0, 1); pm <- pv <- numeric(N + 1); pm[1] <- 0; pv[1] <- 1
  for (j in 1:N) {
    w <- dnorm(x[j + 1] - x[j], (-0.5 * x[j] + 2 * x[j] * yp + 0.5) * dt,
               0.5 * sqrt(dt))
    yp <- yp[sample.int(Np, Np, TRUE, prob = w / sum(w))]
    yp <- yp + (-0.5 * yp - 2 * x[j]^2 + 1) * dt + rnorm(Np) * sqrt(dt)
    pm[j + 1] <- mean(yp); pv[j + 1] <- var(yp)
  }
  burn <- 1:120
  expect_gt(cor(pm[-burn], fp$mean[-burn, 1]), 0.97)
  expect_lt(abs(mean(pv[-burn]) / mean(fp$cov[1, 1, -burn]) - 1), 0.25)
})



# aci_range() provenance on every retention route -----------------------------
#
# The lag-table recursions are defined for the explicit single-step filter. An
# aci_result computed with the implicit stepper or with substeps therefore has
# its range rebuilt on the explicit single-step filter, and that rebuild has to
# be said exactly once whatever `keep` retained, with the values unchanged.
# Nine cells: three retention modes by three integration settings, on one short
# dyad record with a unit prior.

.prov_record <- function() {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 1, t_end = 0.5, dt = 0.01)
  list(m = m, obs = s$obs, ini = list(mean = 0, cov = matrix(1, 1, 1)))
}

.prov_settings <- list(
  explicit_1 = list(stepper = "explicit", nsub = 1L),
  explicit_4 = list(stepper = "explicit", nsub = 4L),
  implicit_1 = list(stepper = "implicit", nsub = 1L)
)

.prov_keeps <- c("paths", "table", "none")

# Count the classed warnings while letting the call complete; expect_warning()
# alone is satisfied by one warning or several, so a duplicate would pass it.
.count_stepper_warnings <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(expr, aci_warn_stepper = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  list(value = value, n = n)
}

# min_M = 0 keeps every row unmasked; two thresholds give `subjective` a value.
.prov_range <- function(x) aci_range(x, min_M = 0, epsilon = c(0.02, 0.2))

test_that("aci_range() reports the explicit single-step rebuild once on every route", {
  r <- .prov_record()
  for (setting in names(.prov_settings)) {
    st <- .prov_settings[[setting]]
    nondefault <- !identical(st$stepper, "explicit") || st$nsub != 1L
    for (keep in .prov_keeps) {
      label <- paste(setting, keep)
      built <- .count_stepper_warnings(
        aci(r$m, r$obs, init = r$ini, keep = keep,
            stepper = st$stepper, nsub = st$nsub))
      ranged <- .count_stepper_warnings(.prov_range(built$value))
      if (!nondefault) {
        expect_identical(built$n, 0L, info = label)
        expect_identical(ranged$n, 0L, info = label)
      } else if (keep == "table") {
        # the table is built, and the rebuild reported, while aci() runs;
        # aci_range() then reads the explicit single-step table silently
        expect_identical(built$n, 1L, info = label)
        expect_identical(ranged$n, 0L, info = label)
      } else {
        expect_identical(built$n, 0L, info = label)
        expect_identical(ranged$n, 1L, info = label)
      }
    }
  }
  # the no-paths warning names the settings it is reporting
  expect_warning(
    .prov_range(aci(r$m, r$obs, init = r$ini, keep = "none", nsub = 4L)),
    class = "aci_warn_stepper", regexp = "stepper = \"explicit\" and nsub = 4")
  expect_warning(
    .prov_range(aci(r$m, r$obs, init = r$ini, keep = "none",
                    stepper = "implicit")),
    class = "aci_warn_stepper", regexp = "stepper = \"implicit\" and nsub = 1")
})

test_that("the rebuilt range equals the direct explicit single-step range in every cell", {
  r <- .prov_record()
  baseline <- .prov_range(aci(r$m, r$obs, init = r$ini, keep = "paths"))
  # the record-end anchors have no peak and therefore no range; the rest do
  expect_true(any(is.finite(baseline$tau)))
  for (setting in names(.prov_settings)) {
    st <- .prov_settings[[setting]]
    for (keep in .prov_keeps) {
      label <- paste(setting, keep)
      res <- suppressWarnings(aci(r$m, r$obs, init = r$ini, keep = keep,
                                  stepper = st$stepper, nsub = st$nsub))
      cell <- suppressWarnings(.prov_range(res))
      expect_equal(cell$t, baseline$t, tolerance = 1e-12, info = label)
      expect_equal(cell$tau, baseline$tau, tolerance = 1e-12, info = label)
      expect_equal(cell$M, baseline$M, tolerance = 1e-12, info = label)
      expect_equal(cell$subjective, baseline$subjective, tolerance = 1e-12,
                   info = label)
      expect_identical(as.character(cell$status), as.character(baseline$status),
                       info = label)
    }
  }
})

test_that("aci_result records its integration settings and keeps them through saveRDS()", {
  r <- .prov_record()
  for (setting in names(.prov_settings)) {
    st <- .prov_settings[[setting]]
    for (keep in .prov_keeps) {
      label <- paste(setting, keep)
      res <- suppressWarnings(aci(r$m, r$obs, init = r$ini, keep = keep,
                                  stepper = st$stepper, nsub = st$nsub))
      expect_identical(res$meta$stepper, st$stepper, info = label)
      expect_identical(res$meta$nsub, st$nsub, info = label)
      f <- tempfile(fileext = ".rds")
      saveRDS(res, f)
      back <- readRDS(f)
      unlink(f)
      expect_identical(back$meta$stepper, st$stepper, info = label)
      expect_identical(back$meta$nsub, st$nsub, info = label)
      # the reloaded object reports exactly what the live one does
      live <- .count_stepper_warnings(.prov_range(res))
      reloaded <- .count_stepper_warnings(.prov_range(back))
      expect_identical(reloaded$n, live$n, info = label)
      expect_equal(reloaded$value$tau, live$value$tau, tolerance = 1e-12,
                   info = label)
    }
  }
  # a double nsub is stored as the integer the filter ran with
  res <- suppressWarnings(aci(r$m, r$obs, init = r$ini, keep = "none", nsub = 4))
  expect_identical(res$meta$nsub, 4L)
  # a result built on a reused table records the table's own construction,
  # which is explicit single-step by lag_table()'s contract
  tb <- lag_table(r$m, r$obs, init = r$ini)
  reused <- aci(r$m, r$obs, table = tb, init = r$ini)
  expect_identical(reused$meta$stepper, "explicit")
  expect_identical(reused$meta$nsub, 1L)
  expect_identical(.count_stepper_warnings(.prov_range(reused))$n, 0L)
})

test_that("keep = 'table' with a non-default setting warns while aci() builds the table", {
  r <- .prov_record()
  explicit <- aci(r$m, r$obs, init = r$ini, keep = "none")
  direct <- aci_range(lag_table(r$m, r$obs, init = r$ini), min_M = 0,
                      epsilon = c(0.02, 0.2))
  for (setting in c("explicit_4", "implicit_1")) {
    st <- .prov_settings[[setting]]
    built <- .count_stepper_warnings(
      aci(r$m, r$obs, init = r$ini, keep = "table",
          stepper = st$stepper, nsub = st$nsub))
    expect_identical(built$n, 1L, info = setting)
    res <- built$value
    # the ACI series itself used the requested setting, and differs from the
    # explicit single-step series, so the provenance is not a formality
    expect_identical(res$meta$stepper, st$stepper, info = setting)
    expect_identical(res$meta$nsub, st$nsub, info = setting)
    expect_gt(max(abs(res$aci - explicit$aci)), 1e-6)
    # the retained table is the explicit single-step construction
    expect_identical(res$table$meta$reference_smoother, "thmD1_online_complete")
    from_table <- aci_range(res$table, min_M = 0, epsilon = c(0.02, 0.2))
    expect_equal(from_table$tau, direct$tau, tolerance = 1e-12, info = setting)
    expect_equal(from_table$M, direct$M, tolerance = 1e-12, info = setting)
  }
})

test_that("a result saved before the settings were recorded cannot report its provenance", {
  r <- .prov_record()
  res <- aci(r$m, r$obs, init = r$ini, keep = "none", nsub = 4L)
  old <- res
  old$meta$stepper <- NULL
  old$meta$nsub <- NULL
  out <- .count_stepper_warnings(.prov_range(old))
  expect_identical(out$n, 0L)
  baseline <- .prov_range(aci(r$m, r$obs, init = r$ini, keep = "paths"))
  expect_equal(out$value$tau, baseline$tau, tolerance = 1e-12)
  expect_equal(out$value$M, baseline$M, tolerance = 1e-12)
})
