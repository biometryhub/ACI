# Golden alignment against the published P1 reference implementation
# (marandmath/ACI_code, dyad_interaction_model.m), ported verbatim in
# helper-golden-p1.R and driven with identical observation series.

golden_setup <- function() {
  prm <- list(d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
              s_x = 0.5, s_y = 1)
  m <- model_dyad()                       # package defaults == golden params
  s <- simulate(m, seed = 7, T = 0.8, dt = 1e-3, burn_in = 0)
  x <- s$obs$x[, 1]; mu0 <- s$hidden[1, 1]
  g <- golden_p1_moments(x, s$obs$dt, prm, mu0 = mu0, R0 = 0.1)
  list(m = m, s = s, x = x, g = g,
       ini = list(mean = mu0, cov = matrix(0.1, 1, 1)))
}

test_that("G1: P1 golden filter, smoother, and ACI agree to machine precision", {
  gs <- golden_setup(); g <- gs$g
  fp <- da_filter(gs$m, gs$s$obs, init = gs$ini)
  sp <- da_smooth(gs$m, gs$s$obs, filter = fp)
  # filter: identical scheme -- machine precision
  expect_lt(max(abs(fp$mean[, 1] - g$fm)), 1e-12)
  expect_lt(max(abs(fp$cov[1, 1, ] - g$fc)), 1e-12)
  # Both implementations use the interval's left/source endpoint.
  expect_lt(max(abs(sp$mean[, 1] - g$sm)), 1e-12)
  expect_lt(max(abs(sp$cov[1, 1, ] - g$sc)), 1e-12)
  a <- aci(gs$m, gs$s$obs, init = gs$ini)
  expect_lt(max(abs(a$aci - g$aci)) / max(g$aci), 1e-12)
  # metric formula identity on the golden moments themselves
  j <- c(100, 400, 700)
  kl <- t(vapply(j, function(i)
    gaussian_kl(g$sm[i], matrix(g$sc[i], 1, 1),
                g$fm[i], matrix(g$fc[i], 1, 1)), numeric(3)))
  expect_lt(max(abs(kl[, 1] - (g$signal[j] + g$dispersion[j]))), 1e-12)
  expect_lt(max(abs(kl[, 2] - g$signal[j])), 1e-12)      # same decomposition
  expect_lt(max(abs(kl[, 3] - g$dispersion[j])), 1e-12)
})

test_that("G2: golden CIR definition matches exact; ratio differs only by quadrature", {
  gs <- golden_setup(); g <- gs$g
  anchors <- seq(60, 700, by = 10)
  gc0 <- golden_p1_cir(g, gs$x, anchors = anchors)
  tab <- lag_table(gs$m, gs$s$obs, mode = "forward", tol = 0, init = gs$ini)
  # Masking stays on: the anchor filter below relies on low-signal rows being
  # NA, so the advisory is asserted rather than disabled.
  expect_warning(fc_ <- forward_cir(tab), class = "aci_warn_low_signal")
  expect_warning(l1 <- forward_cir(tab, method = "l1_linf"),
                 class = "aci_warn_low_signal")
  ok <- anchors[gc0$maxRE[anchors] > 1e-5 &
                is.finite(fc_$tau[anchors]) &
                is.finite(gc0$defn_obj[anchors])]
  expect_gt(length(ok), 40)
  rel_defn <- abs(fc_$tau[ok] - gc0$defn_obj[ok]) /
              pmax(gc0$defn_obj[ok], 1e-12)
  # The residual here is the golden's own route, not a quadrature difference:
  # it integrates the subjective ranges over a 513-point threshold grid, while
  # the package evaluates the layer cake directly in the time domain.
  expect_lt(stats::median(rel_defn), 0.015)          # measured 0.00703
  # The golden helper transcribes simps.m independently, Vandermonde solve
  # included, so this is a parity check against the reference algorithm rather
  # than against a second copy of the package's own rule.  The package's
  # pre-0.0.21 trapezoid close sits at median 3.66e-07 / max 5.45e-05 here.
  rel_appr <- abs(l1$tau[ok] - gc0$approx_obj[ok]) /
              pmax(gc0$approx_obj[ok], 1e-12)
  expect_lt(stats::median(rel_appr), 1e-12)          # measured 3.78e-15
  l1_trap <- suppressWarnings(
    forward_cir(tab, method = "l1_linf", simpson_close = "trapezoid"))
  rel_trap <- abs(l1_trap$tau[ok] - gc0$approx_obj[ok]) /
              pmax(gc0$approx_obj[ok], 1e-12)
  # the legacy close is reachable and is the one that misses the reference
  expect_gt(stats::median(rel_trap), 1e-9)           # measured 3.66e-07
  expect_identical(l1$meta$simpson_close, "quadratic")
  # The golden's own ordering: its Simpson ratio underestimates the
  # definition (their documented approximation), up to grid dust
  expect_true(mean(gc0$approx_obj[ok] <= gc0$defn_obj[ok] + 2e-3) > 0.9)
})
