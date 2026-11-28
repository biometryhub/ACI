# Sub-steps and steppers, against closed-form references ----------------------
#
# dx = (ax x + cx + Lx y) dt + sx dW1, dy = (fy + Ly y) dt + sy dW2, with the
# references solved in closed form over each interval without package code.

ss_par <- list(ax = -0.5, cx = 0.3, Lx = 1.2, Ly = -0.8, fy = 0.4,
               sx = 0.5, sy = 1.0)
ss_par2 <- list(ax = 0.2, cx = -0.1, Lx = 0.7, Ly = -1.5, fy = -0.2,
                sx = 0.8, sy = 0.6)

ss_record <- function(dt, t_end, seed, k = 1L) {
  n <- round(t_end / dt)
  set.seed(seed)
  x <- vapply(seq_len(k), function(i)
    cumsum(c(0, stats::rnorm(n, 0, 0.5 * sqrt(dt)))), numeric(n + 1L))
  observed_trajectory((0:n) * dt, matrix(x, n + 1L, k))
}

ss_library_model <- function(P)
  aci_linear_model(lambda_x = P$Lx, lambda_y = P$Ly,
                   fx = function(t, x) P$ax * x + P$cx, fy = P$fy,
                   sigma_x = P$sx, sigma_y = P$sy)

ss_closure_model <- function(P)
  aci_model(Lx = function(t, x) matrix(P$Lx),
            fx = function(t, x) P$ax * x + P$cx,
            Ly = function(t, x) matrix(P$Ly),
            fy = function(t, x) P$fy,
            Sx1 = function(t, x) matrix(P$sx),
            Sy2 = function(t, x) matrix(P$sy), k = 1, l = 1)

# Two uncoupled copies, one per component: k = l = 2 runs the matrix kernels.
ss_pair_model <- function(P, Q)
  aci_linear_model(lambda_x = c(P$Lx, 0, 0, Q$Lx),
                   lambda_y = c(P$Ly, 0, 0, Q$Ly),
                   fx = function(t, x) c(P$ax * x[1] + P$cx,
                                         Q$ax * x[2] + Q$cx),
                   fy = c(P$fy, Q$fy), sigma_x = c(P$sx, Q$sx),
                   sigma_y = c(P$sy, Q$sy), k = 2, l = 2)

# Exact filter: Riccati R' = a + 2 b R - c R^2 solved through R = W' / (c W),
# W'' - 2 b W' - a c W = 0, and the mean through the same W.
ss_exact_filter <- function(x, t, P, mu0, R0) {
  a <- P$sy^2; b <- P$Ly; c <- P$Lx^2 / P$sx^2
  n <- length(t); dt <- t[2L] - t[1L]
  M <- numeric(n); V <- numeric(n); M[1L] <- mu0; V[1L] <- R0
  for (j in seq_len(n - 1L)) {
    r <- (x[j + 1L] - x[j]) / dt - (P$ax * x[j] + P$cx)
    if (c == 0) {
      M[j + 1L] <- exp(b * dt) * M[j] + P$fy * expm1(b * dt) / b
      V[j + 1L] <- exp(2 * b * dt) * V[j] + a * expm1(2 * b * dt) / (2 * b)
      next
    }
    s <- sqrt(b^2 + a * c); lp <- b + s; lm <- b - s
    A <- (c * V[j] - lm) / (2 * s); B <- 1 - A
    W <- A * exp(lp * dt) + B * exp(lm * dt)
    Wd <- A * lp * exp(lp * dt) + B * lm * exp(lm * dt)
    I0 <- (A * expm1(s * dt) - B * expm1(-s * dt)) / s
    I1 <- (A * lp * expm1(s * dt) - B * lm * expm1(-s * dt)) / s
    M[j + 1L] <- exp(b * dt) / W *
      (M[j] + P$fy * I0 + P$Lx * r / P$sx^2 / c * I1)
    V[j + 1L] <- Wd / (c * W)
  }
  list(mean = M, cov = V)
}

# The smoother's backward equations with the filter mean fm and variance fv
# held at the start of each interval, solved exactly on that interval.
ss_frozen_smoother <- function(fm, fv, dt, P) {
  a <- P$sy^2; b <- P$Ly; n <- length(fm)
  M <- fm; V <- fv
  for (j in (n - 1L):1L) {
    H <- b + a / fv[j]
    e1 <- exp(-H * dt)
    M[j] <- e1 * M[j + 1L] + (-P$fy + a / fv[j] * fm[j]) * (-expm1(-H * dt)) / H
    V[j] <- e1^2 * V[j + 1L] + a * (-expm1(-2 * H * dt)) / (2 * H)
  }
  list(mean = M, cov = V)
}

ss_gap <- function(path, ref, i = 1L)
  max(abs(path$mean[, i] - ref$mean), abs(path$cov[i, i, ] - ref$cov))

# Successive error ratios over nsub = 8, 16, 32, 64; 2 is first order.
ss_ratios <- function(err) err[-length(err)] / err[-1L]

ss_nsubs <- c(8L, 16L, 32L, 64L)

ss_quiet <- function(expr) withCallingHandlers(
  expr, warning = function(w) invokeRestart("muffleWarning"))


test_that(paste0("the sub-stepped filter approaches the exact filter at ",
                 "first order in dt / nsub, on every route and stepper"), {
  ob <- ss_record(0.04, 4, 1)
  ex <- ss_exact_filter(ob$x[, 1], ob$t, ss_par, 0.2, 0.7)
  ini <- list(mean = 0.2, cov = matrix(0.7))
  models <- list(library = ss_library_model(ss_par),
                 closure = ss_closure_model(ss_par))
  for (route in names(models)) for (stepper in c("explicit", "implicit")) {
    err <- vapply(ss_nsubs, function(ns) ss_gap(ss_quiet(
      aci_filter(models[[route]], ob, init = ini, stepper = stepper,
                 nsub = ns)), ex), numeric(1))
    r <- ss_ratios(err)
    expect_true(all(r > 1.8 & r < 2.2),
                info = paste(route, stepper,
                             paste(signif(r, 4), collapse = " ")))
  }

  ob2 <- ss_record(0.04, 4, 3, k = 2L)
  ini2 <- list(mean = c(0.2, -0.1), cov = diag(c(0.7, 0.3)))
  ex2 <- list(ss_exact_filter(ob2$x[, 1], ob2$t, ss_par, 0.2, 0.7),
              ss_exact_filter(ob2$x[, 2], ob2$t, ss_par2, -0.1, 0.3))
  for (stepper in c("explicit", "implicit")) {
    fs <- lapply(ss_nsubs, function(ns) ss_quiet(
      aci_filter(ss_pair_model(ss_par, ss_par2), ob2, init = ini2,
                 stepper = stepper, nsub = ns)))
    for (i in 1:2) {
      r <- ss_ratios(vapply(fs, ss_gap, numeric(1), ref = ex2[[i]], i = i))
      expect_true(all(r > 1.8 & r < 2.2),
                  info = paste("matrix", stepper, i,
                               paste(signif(r, 4), collapse = " ")))
    }
  }
})


test_that(paste0("the sub-stepped smoother approaches its backward ",
                 "equations with the filter moments held per interval, ",
                 "on both kernels"), {
  ob <- ss_record(0.04, 4, 1)
  ini <- list(mean = 0.2, cov = matrix(0.7))
  m <- ss_library_model(ss_par)
  for (stepper in c("explicit", "implicit")) {
    # The smoother takes nsub from the filter; the reference is computed from
    # that same filter.
    err <- vapply(ss_nsubs, function(ns) {
      f <- ss_quiet(aci_filter(m, ob, init = ini, stepper = stepper,
                               nsub = ns))
      ss_gap(aci_smoother(m, ob, filter = f),
             ss_frozen_smoother(f$mean[, 1], f$cov[1, 1, ], 0.04, ss_par))
    }, numeric(1))
    r <- ss_ratios(err)
    expect_true(all(r > 1.8 & r < 2.2),
                info = paste("scalar", stepper,
                             paste(signif(r, 4), collapse = " ")))
  }

  ob2 <- ss_record(0.04, 4, 3, k = 2L)
  ini2 <- list(mean = c(0.2, -0.1), cov = diag(c(0.7, 0.3)))
  m2 <- ss_pair_model(ss_par, ss_par2)
  pars <- list(ss_par, ss_par2)
  for (stepper in c("explicit", "implicit")) {
    err <- vapply(ss_nsubs, function(ns) {
      f <- ss_quiet(aci_filter(m2, ob2, init = ini2, stepper = stepper,
                               nsub = ns))
      s <- aci_smoother(m2, ob2, filter = f)
      max(vapply(1:2, function(i) ss_gap(s, ss_frozen_smoother(
        f$mean[, i], f$cov[i, i, ], 0.04, pars[[i]]), i = i), numeric(1)))
    }, numeric(1))
    r <- ss_ratios(err)
    expect_true(all(r > 1.8 & r < 2.2),
                info = paste("matrix", stepper,
                             paste(signif(r, 4), collapse = " ")))
  }
})


test_that(paste0("a larger nsub does not bring the smoother or the ACI ",
                 "closer to the exact values, and a smaller dt does"), {
  # With Lx = 0 the record carries no information about the hidden state:
  # the exact smoother is the exact filter, and the exact ACI is 0.
  P0 <- modifyList(ss_par, list(Lx = 0))
  m <- ss_library_model(P0)
  ini <- list(mean = 0.2, cov = matrix(0.7))
  run <- function(dt, ns) {
    ob <- ss_record(dt, 4, 1)
    ex <- ss_exact_filter(ob$x[, 1], ob$t, P0, 0.2, 0.7)
    a <- ss_quiet(aci(m, ob, init = ini, nsub = ns))
    s <- a$paths$smoother
    c(mean = max(abs(s$mean[, 1] - ex$mean)),
      var = max(abs(s$cov[1, 1, ] - ex$cov)), aci = max(a$aci))
  }
  g1 <- run(0.04, 1L)
  g64 <- run(0.04, 64L)
  expect_true(all(g64 >= g1), info = paste(signif(c(g1, g64), 4),
                                            collapse = " "))
  expect_true(all(g1 > 0))
  h <- run(0.02, 64L)[["var"]] / g64[["var"]]
  expect_true(h > 0.4 && h < 0.6, info = signif(h, 4))
})


test_that(paste0("a backward-pass covariance failure under the implicit ",
                 "stepper keeps its class and does not advise that stepper"), {
  Ps <- list(ax = 0, cx = 0, Lx = 10, Ly = -1, fy = 0, sx = 0.1, sy = 1)
  ob <- ss_record(0.1, 2, 2)
  ob2 <- ss_record(0.1, 2, 2, k = 2L)
  calls <- list(
    scalar = function() aci(ss_library_model(Ps), ob,
                            init = list(mean = 0, cov = matrix(1)),
                            stepper = "implicit"),
    matrix = function() aci(ss_pair_model(Ps, Ps), ob2,
                            init = list(mean = c(0, 0), cov = diag(2)),
                            stepper = "implicit"))
  for (route in names(calls)) {
    e <- expect_error(calls[[route]](),
                      class = "aci_error_covariance_not_spd", info = route)
    expect_s3_class(e, "aci_error_spd")
    expect_identical(e$site, "smoother_backward", info = route)
    msg <- conditionMessage(e)
    expect_false(grepl("or use stepper = \"implicit\"", msg, fixed = TRUE),
                 info = route)
    expect_match(msg, "backward pass is explicit under either stepper",
                 fixed = TRUE, info = route)
  }
})
