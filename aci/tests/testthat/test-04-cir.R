# The layer-cake form is integrated with composite Simpson, while the brute
# force below sums subjective durations that count cells as (max index - 1).
test_that("C3: exact layer-cake tracks brute-force eps quadrature", {
  set.seed(11); dt <- 0.01
  for (i in 1:20) {
    p <- abs(rnorm(60)) * exp(-seq(0, 3, length.out = 60)) ; p[60] <- 0
    M <- max(p)
    eps <- seq(0, M, length.out = 20001)[-1]
    taus <- sapply(eps, function(e) { idx <- which(p > e)
      if (!length(idx)) 0 else dt * (max(idx) - 1) })
    brute <- mean(taus)
    ex <- unname(.fwd_lengths(p, dt, "exact")["tau"])
    expect_equal(ex, brute, tolerance = 0.06)   # measured max 0.0407
    # Backward: the computable discrepancy grid ends at T-dt.  The
    # layer-cake duration is measured backwards from that physical endpoint;
    # there is no extra, unobserved cell at T.
    P <- cumsum(abs(rnorm(60))) / 10; cc <- abs(P - P[1]); Mb <- max(cc)
    epsb <- seq(0, Mb, length.out = 20001)[-1]
    lag_end <- dt * (length(P) - 1)
    tb <- sapply(epsb, function(e)
      lag_end - dt * (max(which(cc <= e)) - 1))
    ex_b <- unname(.bwd_lengths(P, dt, "exact")["tau"])
    expect_equal(ex_b, mean(tb), tolerance = 0.06)  # measured max 0.0187
  }
})

test_that("backward exact CIR does not invent a terminal +M cell", {
  dt <- 0.1
  P <- c(0, 1, 0.2)

  # The exact form integrates the suffix minima mm = (0, 0.2, 0.2) with
  # composite Simpson, giving (0 + 4*0.2 + 0.2)/3 = 1/3 and a range of dt/3.
  # The ratio form follows the FBCIR code's active line, the plain Appendix G
  # L1 sum of cc = (0, 1, 0.2), giving 1.2 and a range of 1.2*dt.
  got <- unname(.bwd_lengths(P, dt, "exact")["tau"])
  expect_equal(got, dt / 3, tolerance = 1e-15)
  expect_equal(unname(.bwd_lengths(P, dt, "l1_linf")["tau"]),
               1.2 * dt, tolerance = 1e-15)
  # Retaining the terminal zero sentinel that lt_onelag() drops would send
  # every suffix minimum to zero and collapse the range.
  expect_false(isTRUE(all.equal(got, 0)))
})
test_that("C4/C5: l1_linf bound directions and equality under monotonicity", {
  set.seed(12); dt <- 0.02
  for (i in 1:50) {
    p <- abs(rnorm(40)); p[40] <- 0
    f_ex <- .fwd_lengths(p, dt, "exact")["tau"]; f_ap <- .fwd_lengths(p, dt, "l1_linf")["tau"]
    expect_lte(f_ap, f_ex + 1e-12)
    P <- cumsum(abs(rnorm(40)))
    b_ex <- .bwd_lengths(P, dt, "exact")["tau"]; b_ap <- .bwd_lengths(P, dt, "l1_linf")["tau"]
    expect_gte(b_ap, b_ex - 1e-12)
  }
  p_mono <- sort(abs(rnorm(40)), decreasing = TRUE); p_mono[40] <- 0
  gap <- .fwd_lengths(p_mono, dt, "l1_linf")["tau"] - .fwd_lengths(p_mono, dt, "exact")["tau"]
  expect_equal(unname(gap), 0, tolerance = 1e-12)
  # Backward, the two forms use different quadratures (the ratio keeps the
  # FBCIR code's plain L1 sum, the exact form Simpson), so under monotone
  # decay they agree to the quadrature gap rather than exactly.
  c_mono <- cumsum(abs(rnorm(40)))          # nondecreasing bias-corrected metric
  b_ex <- unname(.bwd_lengths(c_mono, dt, "exact")["tau"])
  b_ap <- unname(.bwd_lengths(c_mono, dt, "l1_linf")["tau"])
  expect_lt(abs(b_ap - b_ex) / b_ex, 0.05)  # measured max 0.0318
})
test_that("C6-lite: dyad CIR pipeline runs and produces sane windows", {
  m <- model_dyad()
  s <- simulate(m, seed = 21, T = 3, dt = 5e-3)
  cp <- suppressWarnings(cir_pair(m, s$obs))
  f <- cp$forward
  expect_true(all(is.na(f$tau) | (f$tau >= 0 & f$tau <= max(f$t) + 1e-9)))
  expect_true(f$bound %in% c("layer_cake_objective",
                             "objective_on_truncated_table"))
  b <- cp$backward
  expect_true(is.na(b$tau) || (b$tau > 0 && b$tau <= max(cp$aci$t)))
  expect_equal(unname(b$t), max(cp$aci$t))
  expect_equal(unname(b$interval[, 2]), max(cp$aci$t) - b$dt)
  expect_equal(b$meta$lagged_grid_end, max(cp$aci$t) - b$dt)
  expect_identical(b$meta$endpoint_convention, "T_minus_dt")
  expect_identical(b$bound,
                   "layer_cake_on_O(dt)_T_minus_dt_grid")
})

test_that("backward CIR excludes the terminal one-lag sentinel", {
  m <- model_dyad()
  s <- simulate(m, seed = 9, T = 1, dt = 0.01)
  lt <- suppressWarnings(lag_table(m, s$obs, mode = "one_lag"))
  P <- lt_onelag(lt)
  expect_length(P, length(lt$t) - 1L)
  expect_equal(backward_cir(lt, method = "exact")$tau,
               unname(.bwd_lengths(P, lt$dt, "exact")["tau"]))
  expect_gt(backward_cir(lt, method = "exact")$tau, 0.1)
  be <- backward_cir(lt, method = "exact")
  ba <- backward_cir(lt, method = "l1_linf")
  lag_end <- max(lt$t) - lt$dt
  expect_equal(unname(be$t), max(lt$t))
  expect_equal(unname(be$interval[, 2]), lag_end)
  expect_equal(be$meta$reference_time, max(lt$t))
  expect_equal(be$meta$lagged_grid_end, lag_end)
  expect_identical(be$meta$endpoint_convention, "T_minus_dt")
  expect_identical(be$bound,
                   "layer_cake_on_O(dt)_T_minus_dt_grid")
  expect_identical(ba$bound,
                   "upper_ratio_on_O(dt)_T_minus_dt_grid")
  sb <- backward_cir(lt, method = "exact", eps = c(0, 1e-3))$subjective
  expect_length(sb, 2L)
  expect_true(all(sb >= 0 & sb <= lag_end))
  expect_equal(unname(backward_cir(lt, method = "exact",
                                   eps = 1e100)$subjective), 0)
})

test_that("masked_value = 'zero' reports the published zero convention", {
  # andreou2026cir Rmks B.4/C.4 and the FBCIR scripts set a low-strength
  # length to 0; the default keeps NA so masking stays visible.
  m <- model_dyad()
  s <- simulate(m, seed = 21, T = 0.5, dt = 5e-3)
  tab <- suppressWarnings(lag_table(m, s$obs, mode = "forward"))
  f_na <- suppressWarnings(forward_cir(tab, min_M = 1e6))
  f_z  <- suppressWarnings(forward_cir(tab, min_M = 1e6,
                                       masked_value = "zero"))
  masked <- is.na(f_na$tau) & is.finite(f_na$M) & f_na$M > 1e-14
  expect_gt(sum(masked), 10)
  expect_true(all(f_z$tau[masked] == 0))
  lt <- suppressWarnings(lag_table(m, s$obs, mode = "one_lag"))
  b_z <- suppressWarnings(backward_cir(lt, min_M = 1e6,
                                       masked_value = "zero"))
  expect_identical(unname(b_z$tau), 0)
  b_na <- suppressWarnings(backward_cir(lt, min_M = 1e6))
  expect_true(is.na(b_na$tau))
})

test_that("quadrature options give both published ratio conventions", {
  m <- model_dyad()
  s <- simulate(m, seed = 9, T = 1, dt = 0.01)
  tab <- suppressWarnings(lag_table(m, s$obs, mode = "forward", tol = 0))
  f_simp <- suppressWarnings(forward_cir(tab, method = "l1_linf",
                                         min_M = NULL))
  f_sum <- suppressWarnings(forward_cir(tab, method = "l1_linf",
                                        min_M = NULL, quadrature = "sum"))
  expect_identical(f_simp$meta$quadrature, "simpson")
  expect_identical(f_sum$meta$quadrature, "sum")
  # The sum variant is the literal andreou2026cir eq. G.8 ratio.
  j <- 10L
  p <- lt_row(tab, j, pad = "zero")
  expect_equal(unname(f_sum$tau[j]), tab$dt * sum(p) / max(p),
               tolerance = 1e-12)
  expect_equal(unname(f_simp$tau[j]), tab$dt * .simpson(p) / max(p),
               tolerance = 1e-12)
  lt <- suppressWarnings(lag_table(m, s$obs, mode = "one_lag"))
  b_sum <- backward_cir(lt, method = "l1_linf")
  b_simp <- backward_cir(lt, method = "l1_linf", quadrature = "simpson")
  cc <- abs(lt_onelag(lt) - lt_onelag(lt)[1])
  expect_equal(unname(b_sum$tau), lt$dt * sum(cc) / max(cc),
               tolerance = 1e-12)
  expect_equal(unname(b_simp$tau), lt$dt * .simpson(cc) / max(cc),
               tolerance = 1e-12)
  expect_identical(b_sum$meta$quadrature, "sum")
  expect_identical(b_simp$meta$quadrature, "simpson")
})


test_that("backward validity gate reports the eq. 21 baseline comparison", {
  m <- model_dyad()
  s <- simulate(m, seed = 4, T = 2, dt = 0.01)
  init <- list(mean = 0, cov = diag(1, 1))
  lt <- lag_table(m, as_obs(s), mode = "one_lag", init = init)
  P <- lt_onelag(lt)
  b <- backward_cir(lt)
  # baseline is the raw t = 0 deficit P^0_T, terminal the raw value at T - dt;
  # the gate compares them without touching the centred metric.
  expect_identical(b$meta$baseline, unname(P[1]))
  expect_identical(b$meta$terminal, unname(P[length(P)]))
  expect_identical(b$meta$above_baseline, unname(P[length(P)] > P[1]))
  b3 <- backward_cir(m, obs = as_obs(s), T = c(1, 1.5, 2), init = init)
  ab <- vapply(b3$meta$per_reference, function(mm) mm$above_baseline,
               logical(1))
  expect_length(ab, 3L)
  ct <- suppressWarnings(cir_table(m, as_obs(s), init = init,
                                   direction = "both", at = c(1, 1.5)))
  expect_true(all(is.na(ct$above_baseline[ct$direction == "forward"])))
  expect_identical(ct$above_baseline[ct$direction == "backward"], ab[1:2])
})
