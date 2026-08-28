# Ledger entry C2b: method = "exact" is now the layer-cake sum itself, which
# is the threshold average exactly rather than a quadrature of it, and it
# averages the counting read-out (max index * dt).  Both are asserted: the
# closed form to the bit, and the brute-force threshold average, whose only
# residual is its own finite eps grid.
test_that("C3: exact layer cake is the brute-force eps average", {
  set.seed(11); dt <- 0.01
  worst_count <- 0; worst_lag <- 0
  for (i in 1:20) {
    p <- abs(rnorm(60)) * exp(-seq(0, 3, length.out = 60)) ; p[60] <- 0
    M <- max(p)
    eps <- seq(0, M, length.out = 20001)[-1]
    last <- sapply(eps, function(e) { idx <- which(p > e)
      if (!length(idx)) 0 else max(idx) })
    ex <- unname(.fwd_lengths(p, dt, "exact")["tau"])
    # the closed form: no quadrature stands between these two
    expect_equal(ex, dt * sum(rev(cummax(rev(p)))) / M, tolerance = 0)
    worst_count <- max(worst_count, abs(ex - mean(dt * last)))
    # the eq. G.7 lag-time average sits exactly one grid step below
    worst_lag <- max(worst_lag, abs(ex - dt - mean(dt * pmax(last - 1, 0))))
  }
  expect_lt(worst_count, 1e-4)                  # measured 2.03e-05
  expect_lt(worst_lag, 1e-4)                    # measured 1.98e-05
})

test_that("C4/C5: l1_linf bound directions and equality under monotonicity", {
  set.seed(12); dt <- 0.02
  for (i in 1:50) {
    p <- abs(rnorm(40)); p[40] <- 0
    f_ex <- .fwd_lengths(p, dt, "exact")["tau"]; f_ap <- .fwd_lengths(p, dt, "l1_linf")["tau"]
    expect_lte(f_ap, f_ex + 1e-12)
    P <- cumsum(abs(rnorm(40)))
  }
  p_mono <- sort(abs(rnorm(40)), decreasing = TRUE); p_mono[40] <- 0
  # Ledger entry C2b: under monotone decay the suffix maximum is the row
  # itself, so the two functionals coincide term by term.  The remaining
  # difference is quadrature alone, and it vanishes against the L1 grid sum
  # that the exact form is: this is now an identity, not a tolerance.
  gap <- .fwd_lengths(p_mono, dt, "l1_linf", quadrature = "sum")["tau"] -
         .fwd_lengths(p_mono, dt, "exact")["tau"]
  expect_identical(unname(gap), 0)
  # against Simpson the same pair differs only by the quadrature rule
  gap_simpson <- .fwd_lengths(p_mono, dt, "l1_linf")["tau"] -
                 .fwd_lengths(p_mono, dt, "exact")["tau"]
  expect_gt(abs(unname(gap_simpson)), 1e-12)
})
test_that("C6-lite: dyad CIR pipeline runs and produces sane windows", {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 21, T = 3, dt = 5e-3)
  a <- suppressWarnings(aci(m, s$obs, keep = "table"))
  f <- suppressWarnings(aci_range(a$table))
  expect_true(all(is.na(f$tau) | (f$tau >= 0 & f$tau <= max(f$t) + 1e-9)))
  expect_true(f$bound %in% c("layer_cake_objective",
                             "objective_on_truncated_table"))
})

test_that("masked_value = 'zero' reports the published zero convention", {
  # andreou2026cir Rmks B.4/C.4 and the FBCIR scripts set a low-strength
  # length to 0; the default keeps NA so masking stays visible.
  m <- aci_dyad_model()
  s <- simulate(m, seed = 21, T = 0.5, dt = 5e-3)
  tab <- suppressWarnings(lag_table(m, s$obs, mode = "forward"))
  f_na <- suppressWarnings(aci_range(tab, min_M = 1e6))
  f_z  <- suppressWarnings(aci_range(tab, min_M = 1e6,
                                     masked_value = "zero"))
  masked <- is.na(f_na$tau) & is.finite(f_na$M) & f_na$M > 1e-14
  expect_gt(sum(masked), 10)
  expect_true(all(f_z$tau[masked] == 0))
})

test_that("quadrature options give both published ratio conventions", {
  m <- aci_dyad_model()
  s <- simulate(m, seed = 9, T = 1, dt = 0.01)
  tab <- suppressWarnings(lag_table(m, s$obs, mode = "forward", tol = 0))
  f_simp <- suppressWarnings(aci_range(tab, method = "l1_linf",
                                       min_M = NULL))
  f_sum <- suppressWarnings(aci_range(tab, method = "l1_linf",
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
})
