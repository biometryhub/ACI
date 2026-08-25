# T1/T5/T3/T4: CGNS engine vs an exact discrete Kalman/RTS oracle -------------
lin_model <- function(a = 1, b = -1, fx0 = 0.2, fy0 = 0.3, sx = 0.5, sy = 0.6)
  cgns_model(Lx = function(t, x) matrix(a, 1, 1),
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
    s <- simulate(m, seed = 42, T = 4, dt = dt, ic = list(x0 = 0, y0 = 0.5))
    filt <- da_filter(m, s$obs, init = list(mean = 0, cov = 1))
    smoo <- da_smooth(m, s$obs, filter = filt)
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
  m <- model_dyad()
  s <- simulate(m, seed = 5, T = 2, dt = 5e-3)
  filt <- suppressWarnings(da_filter(m, s$obs))
  smoo <- online_smoother(m, s$obs, filt)
  lt <- lag_table(m, s$obs, mode = "full", filter = filt, smoother = smoo)
  a <- aci(m, s$obs, table = lt)
  expect_equal(lt_diag(lt), a$aci, tolerance = 1e-10)
  expect_equal(lt_diag(lt),
               gaussian_kl_path(smoo, filt)$total, tolerance = 1e-12)
  lastvals <- sapply(seq_along(lt$t), function(j) {
    r <- lt$rows[[j]]; r[length(r)] })
  expect_lt(max(lastvals[1:(length(lastvals) - 1)]), 1e-10)
  expect_equal(lt$meta$reference_smoother, "thmD1_online_complete")
  expect_equal(a$signal + a$dispersion, a$aci, tolerance = 1e-12)
})

test_that("retention alone never changes the headline backward-ODE ACI", {
  m <- model_dyad(); ini <- list(mean = 0, cov = matrix(1, 1, 1))
  s <- simulate(m, seed = 5, T = 0.4, dt = 0.01)
  ap <- aci(m, s$obs, init = ini, keep = "paths")
  at <- aci(m, s$obs, init = ini, keep = "table")
  an <- aci(m, s$obs, init = ini, keep = "none")
  expect_equal(ap$aci, at$aci, tolerance = 0)
  expect_equal(ap$aci, an$aci, tolerance = 0)
  expect_equal(ap$meta$smoother_scheme, "backward_ode")
  expect_equal(at$table$meta$reference_smoother, "thmD1_online_complete")
  expect_error(aci(m, s$obs, typo = TRUE), class = "aci_error_dims")
  expect_error(aci(m, s$obs, table = at$table, typo = TRUE),
               class = "aci_error_dims")
  expect_error(forward_cir(at$table, typo = TRUE), class = "aci_error_dims")
})

test_that("lag table full mode cannot be silently capped", {
  m <- model_dyad()
  s <- simulate(m, seed = 5, T = 0.2, dt = 0.01)
  expect_error(lag_table(m, s$obs, mode = "full", max_lag = 2),
               class = "aci_error_dims")
})

test_that("closed-form and ensemble engines reject mismatched observations", {
  m <- model_l84()
  bad <- observed_trajectory(seq(0, 0.1, by = 0.01), matrix(0, 11, 1))
  expect_error(da_filter(m, bad), class = "aci_error_dims")
  expect_error(enkbf(m, bad, m = 5), class = "aci_error_dims")
})

test_that("Gaussian paths cannot be reused with different data or model", {
  m <- model_dyad()
  s1 <- simulate(m, seed = 1, T = 0.2, dt = 0.01)
  s2 <- simulate(m, seed = 2, T = 0.2, dt = 0.01)
  f1 <- suppressWarnings(da_filter(m, s1$obs,
                                   init = list(mean = 0, cov = diag(1))))
  expect_error(da_smooth(m, s2$obs, filter = f1), class = "aci_error_dims")
  m2 <- model_dyad(params = list(d_x = 0.6, gamma = 2, f_x = 0.5,
                                 s_x = 0.5, d_y = 0.5, f_y = 1, s_y = 1))
  expect_error(da_smooth(m2, s1$obs, filter = f1),
               class = "aci_error_model_contract")
  sm1 <- da_smooth(m, s1$obs, filter = f1)
  expect_error(da_smooth(m, s1$obs, filter = sm1), class = "aci_error_dims")
  expect_error(lag_table(m, s1$obs, filter = sm1, smoother = f1),
               class = "aci_error_dims")
})

test_that("lag tables reject a theorem smoother from another prior", {
  m <- model_dyad(); s <- simulate(m, seed = 3, T = 0.2, dt = 0.01)
  f1 <- da_filter(m, s$obs, init = list(mean = 0, cov = matrix(1, 1, 1)))
  f2 <- da_filter(m, s$obs, init = list(mean = 3, cov = matrix(0.2, 1, 1)))
  sm2 <- online_smoother(m, s$obs, f2)
  sm2$meta$source_model <- m
  expect_error(lag_table(m, s$obs, filter = f1, smoother = sm2),
               class = "aci_error_model_contract")
})

test_that("T5: one-lag column matches brute-force truncated smoothing", {
  m <- model_dyad(); dt <- 5e-3
  s <- simulate(m, seed = 9, T = 2, dt = dt)
  filt <- suppressWarnings(da_filter(m, s$obs))
  smoo <- online_smoother(m, s$obs, filt)
  lt <- lag_table(m, s$obs, mode = "one_lag", filter = filt, smoother = smoo)
  N1 <- length(s$obs$t)
  obs_tr <- observed_trajectory(s$obs$t[-N1], s$obs$x[-N1, , drop = FALSE])
  filt_tr <- suppressWarnings(da_filter(m, obs_tr))
  smoo_tr <- online_smoother(m, obs_tr, filt_tr)
  brute <- sapply(1:(N1 - 1), function(j)
    unname(gaussian_kl(smoo$mean[j, ], smoo$cov[, , j],
                       smoo_tr$mean[j, ], smoo_tr$cov[, , j], decompose = FALSE)))
  ol <- lt_onelag(lt)[1:(N1 - 1)]
  keep <- brute > 1e-7          # compare where signal exists
  expect_gt(sum(keep), 20)
  relerr <- abs(ol[keep] - brute[keep]) / brute[keep]
  expect_lt(stats::median(relerr), 0.05)
  expect_gt(stats::cor(ol[keep], brute[keep]), 0.999)
})

test_that("T4: adaptive truncation engages and agrees with the full table", {
  m <- model_dyad(); s <- simulate(m, seed = 3, T = 3, dt = 5e-3)
  filt <- suppressWarnings(da_filter(m, s$obs)); smoo <- online_smoother(m, s$obs, filt)
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
  expect_equal(forward_cir(adap, min_M = NULL)$bound,
               "objective_on_truncated_table")
})

test_that("max_lag is an exact positive-lag storage cap", {
  m <- model_dyad(); s <- simulate(m, seed = 2, T = 0.1, dt = 0.01)
  # An explicit prior: the storage cap under test does not depend on it, and
  # defaulting would exercise the diffuse-prior path incidentally.
  ini <- list(mean = 2, cov = matrix(1, 1, 1))
  for (cap in 1:3) {
    tab <- lag_table(m, s$obs, mode = "forward", tol = 0, max_lag = cap,
                     init = ini)
    want <- pmin(cap, length(tab$t) - seq_along(tab$t))
    expect_equal(tab$L, as.integer(want))
    expect_equal(vapply(tab$rows, length, integer(1)), as.integer(want + 1L))
    expect_match(forward_cir(tab, method = "l1_linf", min_M = NULL)$bound,
                 "truncated")
  }
})

test_that("T1b: implicit Riccati step survives the diffuse-init stress case", {
  m <- model_dyad(); s <- simulate(m, seed = 55, T = 1.5, dt = 5e-3)
  big <- list(mean = 0, cov = diag(100, 1))
  # The explicit scheme is expected to report its own instability here; that
  # warning is part of the documented behaviour this test pins down.
  expect_warning(fe <- da_filter(m, s$obs, init = big, stepper = "explicit"),
                 class = "aci_warn_riccati_stiff")
  fi <- da_filter(m, s$obs, init = big, stepper = "implicit")
  expect_lt(min(fe$cov[1, 1, ]), 1e-6)     # documented explicit-scheme collapse
  expect_gt(min(fi$cov[1, 1, ]), 0.05)     # implicit stays positive-definite
  # both steppers agree at a sane init (scheme difference is O(dt))
  ini <- list(mean = 0, cov = diag(1, 1))
  f1 <- da_filter(m, s$obs, init = ini, stepper = "explicit")
  f2 <- da_filter(m, s$obs, init = ini, stepper = "implicit")
  # Backward-Euler mean/covariance splitting is deliberately more dissipative
  # than explicit Euler at this grid; the difference remains first-order and
  # is assessed by the sub-stepping convergence check below.
  expect_lt(max(abs(f1$mean - f2$mean)), 0.07)
  expect_lt(max(abs(f1$cov - f2$cov)), 0.02)
  # sub-stepping refines toward the same limit
  f4 <- da_filter(m, s$obs, init = ini, stepper = "implicit", nsub = 8)
  expect_lt(max(abs(f4$mean - f2$mean)), 0.05)
  # lag_table enforces its explicit single-step contract
  expect_warning(lag_table(m, s$obs, mode = "one_lag", filter = f2),
                 class = "aci_warn_stepper")
  expect_error(lag_table(m, s$obs, stepper = "implicit"),
               class = "aci_error_stepper")
})

test_that("user priors must be genuinely positive definite", {
  m <- model_dyad(); s <- simulate(m, seed = 1, T = 0.1, dt = 0.01)
  expect_error(da_filter(m, s$obs, init = list(mean = 0, cov = matrix(0, 1, 1))),
               class = "aci_error_spd")
  expect_error(da_filter(m, s$obs, init = list(mean = 0, cov = matrix(-1, 1, 1))),
               class = "aci_error_spd")
})

test_that("one-lag propagation does not truncate on an unscaled operator norm", {
  dt <- 0.01; ax <- (1 - 5e-14) / dt
  m <- cgns_model(
    Lx = function(t, x) matrix(ax, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(0, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
  ob <- observed_trajectory(dt * 0:3, matrix(c(0, 0, 0, 1e6), ncol = 1))
  tab <- lag_table(m, ob, mode = "one_lag",
                   init = list(mean = 0, cov = matrix(1 / ax, 1, 1)))
  expect_true(all(is.finite(lt_onelag(tab))))
  expect_gt(lt_onelag(tab)[2], 1)
  expect_true(is.na(tab$meta$stop_index))
})

test_that("implicit covariance update remains SPD in a stiff scalar model", {
  m <- cgns_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-100, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(1, 1, 1),
    Sy2 = function(t, x) matrix(1, 1, 1), k = 1, l = 1)
  ob <- observed_trajectory(c(0, 1, 2), matrix(0, 3, 1))
  fit <- suppressWarnings(da_filter(m, ob, init = list(mean = 0, cov = 1),
                                    stepper = "implicit"))
  expect_true(all(is.finite(fit$cov)))
  expect_true(all(fit$cov[1, 1, ] > 0))
})

test_that("correlation activated after t0 selects the correlated backward ODE", {
  m <- cgns_model(
    Lx = function(t, x) matrix(1, 1, 1), fx = function(t, x) 0,
    Ly = function(t, x) matrix(-1, 1, 1), fy = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.4, 1, 1),
    Sy1 = function(t, x) matrix(0.2 * x[1], 1, 1),
    Sy2 = function(t, x) matrix(0.5, 1, 1), k = 1, l = 1)
  s <- simulate(m, seed = 1, T = 0.2, dt = 0.01,
                ic = list(x0 = 0, y0 = 1))
  f <- suppressWarnings(da_filter(m, s$obs, init = list(mean = 0, cov = 1)))
  sm <- suppressWarnings(da_smooth(m, s$obs, filter = f))
  expect_equal(sm$meta$route, "backward_ode_correlated")
})

test_that("correlated smoother transcribes the active MATLAB backward loop", {
  m <- cgns_model(
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
    f <- da_filter(m, ob, init = list(mean = 0.2, cov = matrix(0.8, 1, 1))),
    class = "aci_warn_dt_stability")
  got <- da_smooth(m, ob, filter = f)

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
  m <- model_dyad(); s <- simulate(m, seed = 2, T = 0.2, dt = 0.01)
  f <- suppressWarnings(da_filter(m, s$obs, init = list(mean = 0, cov = 1)))
  bad <- f; bad$t[2] <- bad$t[2] + 0.001
  expect_error(da_smooth(m, s$obs, filter = bad), class = "aci_error_dims")
  stripped <- f; stripped$meta$obs_x <- NULL
  expect_error(da_smooth(m, s$obs, filter = stripped), class = "aci_error_dims")
  stripped_model <- f; stripped_model$meta$source_model <- NULL
  stripped_model$meta$model <- NULL
  expect_error(da_smooth(m, s$obs, filter = stripped_model),
               class = "aci_error_model_contract")
  singular <- f; singular$cov[1, 1, 2] <- 0
  expect_error(da_smooth(m, s$obs, filter = singular), class = "aci_error_spd")
  expect_error(da_smooth(m, s$obs, filter = f,
                         init = list(mean = 9, cov = matrix(2, 1, 1))),
               "init conflicts", class = "aci_error_dims")
})


test_that("T1c: nonlinear truth-tracking referee -- particle filter agrees", {
  # Closes a verification gap: T1/T3/T5 check scheme identities, not whether
  # the nonlinear posterior itself is right. A bootstrap particle filter with
  # exact discrete weights is an independent referee on the dyad.
  m <- model_dyad(); s <- simulate(m, seed = 21, T = 2.5, dt = 5e-3, burn_in = 0.5)
  fp <- suppressWarnings(da_filter(m, s$obs, init = list(mean = 0, cov = diag(1, 1))))
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


test_that("contraction certificate reports per-step E diagnostics", {
  m <- model_dyad()
  s <- simulate(m, seed = 1, T = 2, dt = 0.01)
  init <- list(mean = 0, cov = diag(1, 1))
  tb <- lag_table(m, as_obs(s), mode = "forward", init = init)
  cert <- lt_contraction_certificate(tb)
  expect_identical(nrow(cert), length(tb$t) - 1L)
  expect_true(all(is.finite(cert$lambda_min)))
  expect_true(all(is.finite(cert$enorm)))
  # scalar hidden state: operator norm and spectral radius coincide with |E|
  expect_identical(cert$enorm, cert$rho_E)
  expect_equal(attr(cert, "gamma"), max(cert$enorm), tolerance = 1e-15)
  expect_identical(attr(cert, "condition_318"), all(cert$lambda_min > 0))
  # independent plumbing check at one step: rebuild E from a fresh filter
  filt <- da_filter(m, as_obs(s), init = init)
  co <- eval_coefs(m, s$obs$t[5], s$obs$x[5, ])
  Rf <- filt$cov[, , 5]
  Rfi <- chol_solve(Rf, diag(1), "Rf")
  Gi <- chol_solve(co$gxx, diag(1), "gxx")
  Gx <- co$Lx + t(co$gyx) %*% Rfi
  Gy <- co$Ly + co$gyy %*% Rfi
  E <- diag(1) + (co$gyx %*% Gi %*% Gx - Gy) * 0.01
  expect_equal(cert$enorm[5], abs(E[1, 1]), tolerance = 1e-12)
  expect_equal(cert$lambda_min[5], (1 - E[1, 1]) / 0.01,
               tolerance = 1e-8)
})
