# Empirical probe for the 2026-08-14 aciR review.
# Reconstructs the unexplained 4.58e-09 CIR residual and tests
# neighbouring hypotheses. Not a package test; a reviewer notebook.

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools needed for load_all")
  }
})

# Run from the repository root, or set ACIR_ROOT. Absolute paths are not
# hardcoded: a probe that runs on one machine only is not a probe.
root <- Sys.getenv("ACIR_ROOT", unset = normalizePath(".", mustWork = TRUE))
pkg <- file.path(root, "aciR")
devtools::load_all(pkg, quiet = TRUE, export_all = TRUE)

out_dir <- file.path(root, "design", "artefacts")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat_line <- function(...) cat(paste0(..., "\n"))

# -- MATLAB simps.m transcribed -----------------------------------------------
#
# Damien Garcia, 2009. Even-length closure is a quadratic interpolant through
# the last three abscissae, integrated over the last interval only. That is
# not the Simpson 3/8 panel aciR uses.

matlab_simps <- function(y, x = NULL) {
  y <- as.numeric(y)
  m <- length(y)
  if (is.null(x)) {
    x <- as.numeric(seq_len(m))
  } else {
    x <- as.numeric(x)
  }
  if (m < 3L) {
    return(sum(diff(x) * (y[-1L] + y[-m]) / 2))
  }
  dx <- diff(x)
  dx1 <- dx[-length(dx)]
  dx2 <- dx[-1L]
  alpha <- (dx1 + dx2) / dx1 / 6
  a0 <- alpha * (2 * dx1 - dx2)
  a1 <- alpha * (dx1 + dx2)^2 / dx2
  a2 <- alpha * dx1 / dx2 * (2 * dx2 - dx1)
  left <- seq.int(1L, length(a0), by = 2L)
  z <- sum(a0[left] * y[left] + a1[left] * y[left + 1L] +
             a2[left] * y[left + 2L])
  if (m %% 2L == 0L) {
    xv <- x[(m - 2L):m]
    yv <- y[(m - 2L):m]
    V <- cbind(xv^2, xv, rep(1, 3L))
    C <- solve(V, yv)
    z <- z + C[1L] * (x[m]^3 - x[m - 1L]^3) / 3 +
      C[2L] * (x[m]^2 - x[m - 1L]^2) / 2 +
      C[3L] * (x[m] - x[m - 1L])
  }
  z
}

# -- 1. Closures on a test function -------------------------------------------

cat_line("=== 1. Simpson closures on exp(sin(x)) over [0, 2] ===")
f <- function(t) exp(sin(t))
exact <- stats::integrate(f, 0, 2, rel.tol = 1e-12)$value
for (n in c(5L, 6L, 7L, 8L, 11L, 12L, 21L, 22L, 129L, 130L)) {
  x <- seq(0, 2, length.out = n)
  y <- f(x)
  a <- aciR:::.aci_simpson(y, x)
  b <- matlab_simps(y, x)
  cat_line(sprintf(
    "  n=%3d intervals=%3d  aciR_err=% .3e  matlab_err=% .3e  |aciR-matlab|=% .3e",
    n, n - 1L, abs(a - exact), abs(b - exact), abs(a - b)
  ))
}

# -- 2. Closures when the last sample is exactly zero -------------------------

cat_line("=== 2. Last sample forced to 0 (CIR tail identity) ===")
for (n in c(6L, 7L, 20L, 21L, 100L, 101L, 1051L, 1052L)) {
  x <- seq_len(n)
  y <- exp(-seq_len(n) / (n / 5))
  y[n] <- 0
  a <- aciR:::.aci_simpson(y)
  b <- matlab_simps(y)
  cat_line(sprintf(
    "  n=%4d  aciR=% .12f  matlab=% .12f  diff=% .3e  last3=(%.3e, %.3e, 0)",
    n, a, b, abs(a - b), y[n - 2L], y[n - 1L]
  ))
}

# -- 3. Reconstruct the 4.58e-09 comparison -----------------------------------

cat_line("=== 3. arbitrary_cross_noise CIR at the reference horizon ===")
ds_path <- file.path(root, "oracle", "parity", "datasets",
                     "arbitrary_cross_noise")
meta <- as.list(read.dcf(file.path(ds_path, "meta.dcf"))[1L, ])
num <- c("N", "dt", "mu0", "R0", "L_y", "Sx_1", "Sx_2", "Sy_1", "Sy_2",
         "CIRStart", "CIREnd")
for (nm in intersect(num, names(meta))) meta[[nm]] <- as.numeric(meta[[nm]])
arr <- utils::read.csv(file.path(ds_path, "arrays.csv"))
s_yox <- meta$Sy_1 * meta$Sx_1 + meta$Sy_2 * meta$Sx_2
comp <- list(
  L_x = arr$L_x, L_y = meta$L_y, f_x = arr$f_x, f_y = arr$f_y,
  S_xoS_x = meta$Sx_1^2 + meta$Sx_2^2,
  S_yoS_y = meta$Sy_1^2 + meta$Sy_2^2,
  S_yoS_x = s_yox, S_xoS_y = s_yox
)
x <- arr$x
dt <- meta$dt
n <- length(x)
filt <- aci_filter(x, comp, dt, mu0 = meta$mu0, R0 = meta$R0)
aux <- aciR:::.aci_online_aux(x, comp, dt, filt)

# Reference last_idx: time_end_plot + lookahead_tolerance = 2 + 0.6
lookahead <- 0.6
first_idx <- as.integer(round(meta$CIRStart / dt) + 1L)
last_idx <- as.integer(round((meta$CIREnd + lookahead) / dt) + 1L)
plot_end <- as.integer(round(meta$CIREnd / dt) + 1L)
cat_line(sprintf(
  "  n=%d dt=%g first_idx=%d last_idx=%d plot_end=%d horizon_len=%d",
  n, dt, first_idx, last_idx, plot_end, last_idx - first_idx + 1L
))

mat <- utils::read.csv(file.path(
  root, "oracle", "parity", "reports",
  "matlab_arbitrary_cross_noise_cir.csv"
))

# Default 129-point grid is not what the reference used. The 4.58e-09
# residual is on the APPROXIMATE objective, which does not use epsilon.
window <- first_idx:last_idx
reported <- window <= plot_end

# Compute one row at a representative even-length and odd-length time
# and also the full comparison against the MATLAB CSV.

row_stats <- function(j, horizon) {
  re <- aciR:::.aci_cir_row(aux, filt, j, n, horizon = horizon)
  peak <- max(re)
  obj_aci <- aciR:::.aci_simpson(re) * dt / peak
  obj_mat <- matlab_simps(re) * dt / peak
  # Naive (MATLAB) dispersion vs log1p form, on this row's partials.
  # Rebuild the two RE sequences from the same posteriors.
  k <- seq.int(j, n - 1L)
  offset <- k - j
  span <- aux$cum_log[k] - aux$cum_log[j]
  sgn <- aux$cum_sign[k] * aux$cum_sign[j]
  d <- sgn * exp(span)
  d[offset == 0L] <- 1
  mu <- filt$mean[j] + c(0, cumsum(d * aux$innov_mean[k]))
  rr <- filt$cov[j] + c(0, cumsum(d * d * aux$innov_cov[k]))
  mu_end <- mu[length(mu)]
  r_end <- rr[length(rr)]
  keep <- seq_len(min(length(mu), horizon - j + 1L))
  mu <- mu[keep]
  rr <- rr[keep]
  delta <- r_end / rr - 1
  re_log1p <- pmax(0.5 * (mu_end - mu)^2 / rr + 0.5 * (delta - log1p(delta)), 0)
  cov_ratio <- r_end / rr
  re_log <- pmax(
    0.5 * (mu_end - mu)^2 / rr + 0.5 * (cov_ratio - 1 - log(cov_ratio)),
    0
  )
  list(
    j = j, n_re = length(re), peak = peak,
    obj_aci = obj_aci, obj_mat = obj_mat,
    quad_diff = abs(obj_aci - obj_mat),
    re_form_max = max(abs(re_log1p - re_log)),
    re_tail = re[length(re)],
    re_penult = re[length(re) - 1L],
    monotone = all(diff(re) <= 1e-15)
  )
}

s_even <- row_stats(first_idx, last_idx)       # length last_idx-first+1 = 1051 odd
s_odd  <- row_stats(first_idx + 1L, last_idx)  # 1050 even
cat_line(sprintf(
  "  j=%d n_re=%d (odd intervals=%s)  quad_diff=% .3e  re_form=% .3e  mono=%s",
  s_even$j, s_even$n_re, (s_even$n_re - 1L) %% 2L == 1L,
  s_even$quad_diff, s_even$re_form_max, s_even$monotone
))
cat_line(sprintf(
  "  j=%d n_re=%d (odd intervals=%s)  quad_diff=% .3e  re_form=% .3e  mono=%s",
  s_odd$j, s_odd$n_re, (s_odd$n_re - 1L) %% 2L == 1L,
  s_odd$quad_diff, s_odd$re_form_max, s_odd$monotone
))

# Full sweep against MATLAB CSV, with aciR's own Simpson and with matlab_simps
# applied to aciR's RE row. Also record saturation at default margin.
t0 <- proc.time()[[3L]]
n_win <- length(window)
obj_aci <- numeric(n_win)
obj_msimps <- numeric(n_win)
peak_aci <- numeric(n_win)
n_re <- integer(n_win)
quad_diff <- numeric(n_win)
mono <- logical(n_win)
settled_frac <- numeric(n_win)
for (i in seq_along(window)) {
  j <- window[i]
  re <- aciR:::.aci_cir_row(aux, filt, j, n, horizon = last_idx)
  n_re[i] <- length(re)
  peak_aci[i] <- max(re)
  if (peak_aci[i] > 1e-5 && length(re) >= 3L) {
    obj_aci[i] <- aciR:::.aci_simpson(re) * dt / peak_aci[i]
    obj_msimps[i] <- matlab_simps(re) * dt / peak_aci[i]
  }
  quad_diff[i] <- abs(obj_aci[i] - obj_msimps[i])
  mono[i] <- all(diff(re) <= 1e-14)
  suffix <- rev(cummax(rev(re)))
  settled <- length(suffix) - findInterval(1e-5, rev(suffix))
  settled_frac[i] <- settled / length(re)
}
cat_line(sprintf("  row sweep %.2fs", proc.time()[[3L]] - t0))

# Align MATLAB CSV on index
stopifnot(identical(mat$index, window))
d_peak <- abs(peak_aci - mat$peak)
d_obj_native <- abs(obj_aci - mat$objective)
d_obj_msimps <- abs(obj_msimps - mat$objective)

rep <- reported
cat_line("  reported region (index <= plot_end):")
cat_line(sprintf("    n=%d  max|peak|=% .3e", sum(rep), max(d_peak[rep])))
cat_line(sprintf(
  "    max|obj aciR-matlabCSV|=% .3e  (n_gt_1e-12=%d  n_gt_1e-10=%d  n_gt_1e-8=%d)",
  max(d_obj_native[rep]),
  sum(d_obj_native[rep] > 1e-12),
  sum(d_obj_native[rep] > 1e-10),
  sum(d_obj_native[rep] > 1e-8)
))
cat_line(sprintf(
  "    max|obj matlab_simps(aciR_re)-matlabCSV|=% .3e",
  max(d_obj_msimps[rep])
))
cat_line(sprintf(
  "    max|aciR_simpson - matlab_simps| on aciR rows =% .3e",
  max(quad_diff[rep])
))
cat_line(sprintf(
  "    even-n_re rows: %d  odd-n_re rows: %d",
  sum(n_re[rep] %% 2L == 0L), sum(n_re[rep] %% 2L == 1L)
))
cat_line(sprintf(
  "    max quad_diff even=% .3e  odd=% .3e",
  max(quad_diff[rep & n_re %% 2L == 0L]),
  max(quad_diff[rep & n_re %% 2L == 1L])
))
cat_line(sprintf(
  "    non-monotone RE rows in reported region: %d / %d",
  sum(!mono[rep]), sum(rep)
))
cat_line(sprintf(
  "    settled_frac (RE last-exit / row length): min=%.3f median=%.3f max=%.3f",
  min(settled_frac[rep]), stats::median(settled_frac[rep]), max(settled_frac[rep])
))
cat_line(sprintf(
  "    fraction of reported times with settled_frac > 0.9: %.3f",
  mean(settled_frac[rep] > 0.9)
))

# Where is the residual largest?
i_worst <- which.max(d_obj_native[rep])
j_worst <- window[rep][i_worst]
cat_line(sprintf(
  "    worst objective residual at j=%d n_re=%d aciR=%.12f matlab=%.12f diff=%.3e",
  j_worst, n_re[window == j_worst],
  obj_aci[window == j_worst],
  mat$objective[mat$index == j_worst],
  d_obj_native[rep][i_worst]
))
cat_line(sprintf(
  "    same row, matlab_simps(aciR_re)=%.12f  residual to CSV=%.3e",
  obj_msimps[window == j_worst],
  abs(obj_msimps[window == j_worst] - mat$objective[mat$index == j_worst])
))

# Default margin vs reference lookahead
rng_default <- aci_cir(
  x, comp, dt, filt = filt, window = window[rep],
  margin = 0.1, horizon = last_idx
)
rng_tiny <- aci_cir(
  x, comp, dt, filt = filt, window = window[rep],
  margin = 0.001, horizon = last_idx
)
cat_line(sprintf(
  "    default margin=0.1 saturated %d / %d reported times",
  sum(rng_default$saturated), length(rng_default$saturated)
))
cat_line(sprintf(
  "    margin=0.001 saturated %d / %d",
  sum(rng_tiny$saturated), length(rng_tiny$saturated)
))

# -- 4. 3/8 on a logarithmic epsilon grid -------------------------------------

cat_line("=== 4. 3/8 closure on a logarithmic abscissa ===")
# The exact objective integrates subjective ranges over a log-spaced epsilon
# grid. The 3/8 panel assumes equal spacing. 129 points -> 128 intervals
# (even) so 1/3 only; 128 points would fire 3/8.
for (n_eps in c(8L, 9L, 64L, 65L, 128L, 129L, 512L, 513L)) {
  eps <- 10^seq(-6, 0.5, length.out = n_eps)
  # Fake a decaying last-exit curve ~ -log10(eps)
  y <- pmax(-log10(eps), 0)
  a <- aciR:::.aci_simpson(y, eps)
  b <- matlab_simps(y, eps)
  # Reference "truth": integrate the interpolant by adaptive quadrature of
  # a monotone interpolant through the nodes (PCHIP-like is overkill);
  # use a very fine trapezoid on the same function as a check of spacing.
  eps_fine <- 10^seq(-6, 0.5, length.out = 20001L)
  y_fine <- pmax(-log10(eps_fine), 0)
  trap_fine <- sum(diff(eps_fine) * (y_fine[-1L] + y_fine[-length(y_fine)]) / 2)
  cat_line(sprintf(
    "  n_eps=%3d intervals=%3d  aciR-matlab=% .3e  aciR-fine_trap=% .3e  matlab-fine=% .3e",
    n_eps, n_eps - 1L, abs(a - b), abs(a - trap_fine), abs(b - trap_fine)
  ))
}

# -- 5. Vector smoother truncation --------------------------------------------

cat_line("=== 5. Vector online-smoother truncation ===")
set.seed(19L)
n_mv <- 400L
dt_mv <- 0.01
x_mv <- rbind(
  cumsum(stats::rnorm(n_mv + 1L, sd = 0.05)),
  cumsum(stats::rnorm(n_mv + 1L, sd = 0.05))
)
comp_mv <- list(
  L_x = matrix(c(-0.6, 0.2, 0.1, -0.5), 2L, 2L),
  L_y = matrix(c(-0.9, 0.15, 0.05, -0.7), 2L, 2L),
  f_x = matrix(0.05, 2L, n_mv + 1L),
  f_y = matrix(0.10, 2L, n_mv + 1L),
  S_xoS_x = diag(c(0.16, 0.09)),
  S_yoS_y = diag(c(0.25, 0.12)),
  S_yoS_x = matrix(0, 2L, 2L)
)
filt_mv <- aci_filter(x_mv, comp_mv, dt_mv, mu0 = c(0, 0), R0 = diag(2L) * 0.1)
aux_mv <- aciR:::.aci_online_aux_mv(
  x_mv, aciR:::.aci_check_components_mv(comp_mv, x_mv), dt_mv, filt_mv
)

spec_rad <- function(m) {
  ev <- eigen(m, only.values = TRUE)$values
  max(Mod(ev))
}
sr <- vapply(seq_len(n_mv + 1L), function(j) spec_rad(aux_mv$E_j[, , j]),
             numeric(1))
cat_line(sprintf(
  "  spectral radius of E_j: min=%.6f median=%.6f max=%.6f  n_ge_1=%d",
  min(sr), stats::median(sr), max(sr), sum(sr >= 1)
))

# Sensitivity of full-lag online smoother to tol
tols <- c(1e-8, 1e-12, 1e-18, 1e-30)
refs <- lapply(tols, function(tol) {
  aci_online_smoother(x_mv, comp_mv, dt_mv, filt_mv, lag = Inf, tol = tol)
})
cat_line("  tol sensitivity (vs tol=1e-30):")
for (i in seq_along(tols)) {
  dmu <- max(abs(refs[[i]]$mean - refs[[length(tols)]]$mean))
  dcov <- max(abs(refs[[i]]$cov - refs[[length(tols)]]$cov))
  cat_line(sprintf(
    "    tol=% .0e  lag_eff=%d  max|dmu|=% .3e  max|dcov|=% .3e",
    tols[i], refs[[i]]$lag_effective, dmu, dcov
  ))
}

# Does CIR-row truncation change mu_end relative to a no-truncation walk?
# Force a huge tol so it truncates immediately after the first tiny product.
re_tight <- aciR:::.aci_cir_row_mv(aux_mv, filt_mv, 30L, n_mv + 1L,
                                   tol = 1e-18, horizon = n_mv + 1L)
re_loose <- aciR:::.aci_cir_row_mv(aux_mv, filt_mv, 30L, n_mv + 1L,
                                   tol = 1e-8, horizon = n_mv + 1L)
re_none  <- aciR:::.aci_cir_row_mv(aux_mv, filt_mv, 30L, n_mv + 1L,
                                   tol = 0, horizon = n_mv + 1L)
cat_line(sprintf(
  "  CIR row j=30: tight-vs-none max|dRE|=% .3e  loose-vs-none=% .3e  n=%d",
  max(abs(re_tight - re_none)), max(abs(re_loose - re_none)), length(re_none)
))

# max(abs(d)) vs Frobenius vs spectral radius along one accumulation
j0 <- 30L
d <- diag(2L)
norms <- matrix(NA_real_, n_mv + 1L - j0, 3L)
colnames(norms) <- c("maxabs", "frob", "spec")
for (k in seq.int(j0, n_mv)) {
  d <- d %*% aux_mv$E_j[, , k]
  i <- k - j0 + 1L
  norms[i, 1L] <- max(abs(d))
  norms[i, 2L] <- sqrt(sum(d^2))
  norms[i, 3L] <- spec_rad(d)
}
cut_maxabs <- which(norms[, 1L] < 1e-18)[1L]
cut_frob   <- which(norms[, 2L] < 1e-18)[1L]
cut_spec   <- which(norms[, 3L] < 1e-18)[1L]
cat_line(sprintf(
  "  product decay from j=30: first <1e-18 at offset maxabs=%s frob=%s spec=%s (of %d)",
  cut_maxabs, cut_frob, cut_spec, nrow(norms)
))

# -- 6. Horizon prefix + saturation shape -------------------------------------

cat_line("=== 6. Horizon prefix and saturation geometry ===")
# Confirm truncated row is a prefix; measure how much objective changes
# when horizon grows, as a check that the target posterior is fixed.
j_mid <- first_idx + 100L
full <- aciR:::.aci_cir_row(aux, filt, j_mid, n, horizon = n)
cut  <- aciR:::.aci_cir_row(aux, filt, j_mid, n, horizon = last_idx)
cat_line(sprintf(
  "  prefix identical: %s  |full|=%d |cut|=%d",
  isTRUE(all.equal(cut, full[seq_along(cut)], tolerance = 0)),
  length(full), length(cut)
))
# If the target moved with the horizon, the prefix would differ at the
# first entry (RE of full-informed vs filter). Record that first entry
# across horizons.
horizons <- c(j_mid + 50L, j_mid + 200L, last_idx, n)
first_re <- vapply(horizons, function(h) {
  aciR:::.aci_cir_row(aux, filt, j_mid, n, horizon = h)[1L]
}, numeric(1))
cat_line(sprintf(
  "  RE[1] (filter vs target) across horizons %s: %s",
  paste(horizons, collapse = ","),
  paste(sprintf("%.6e", first_re), collapse = ", ")
))
cat_line(sprintf(
  "  max |RE[1] - RE[1]_full| = %.3e  (zero iff target ignores horizon)",
  max(abs(first_re - first_re[length(first_re)]))
))

# -- 7. log1p vs log on a CIR row ---------------------------------------------

cat_line("=== 7. Dispersion form: log1p vs log, on the worst row ===")
re_aci <- aciR:::.aci_cir_row(aux, filt, j_worst, n, horizon = last_idx)
# Rebuild MATLAB-style RE from the same posteriors
j <- j_worst
k <- seq.int(j, n - 1L)
offset <- k - j
span <- aux$cum_log[k] - aux$cum_log[j]
sgn <- aux$cum_sign[k] * aux$cum_sign[j]
d <- sgn * exp(span)
d[offset == 0L] <- 1
mu <- filt$mean[j] + c(0, cumsum(d * aux$innov_mean[k]))
rr <- filt$cov[j] + c(0, cumsum(d * d * aux$innov_cov[k]))
mu_end <- mu[length(mu)]
r_end <- rr[length(rr)]
keep <- seq_len(min(length(mu), last_idx - j + 1L))
mu <- mu[keep]
rr <- rr[keep]
delta <- r_end / rr - 1
re_log1p <- pmax(0.5 * (mu_end - mu)^2 / rr + 0.5 * (delta - log1p(delta)), 0)
cov_ratio <- r_end / rr
re_log <- pmax(
  0.5 * (mu_end - mu)^2 / rr + 0.5 * (cov_ratio - 1 - log(cov_ratio)), 0
)
cat_line(sprintf(
  "  max|log1p - log|=% .3e  at index %d (cov_ratio=%.6f)",
  max(abs(re_log1p - re_log)),
  which.max(abs(re_log1p - re_log)),
  cov_ratio[which.max(abs(re_log1p - re_log))]
))
cat_line(sprintf(
  "  simpson(log1p)*dt/peak=% .16f",
  aciR:::.aci_simpson(re_log1p) * dt / max(re_log1p)
))
cat_line(sprintf(
  "  simpson(log)*dt/peak   =% .16f",
  aciR:::.aci_simpson(re_log) * dt / max(re_log)
))
cat_line(sprintf(
  "  matlab_simps(log)*dt/peak=% .16f",
  matlab_simps(re_log) * dt / max(re_log)
))

# -- 8. Timing: scalar filter / CIR / vector online ---------------------------

cat_line("=== 8. Timings (one core) ===")
time_one <- function(expr, nrep = 3L) {
  expr <- substitute(expr)
  env <- parent.frame()
  t <- replicate(nrep, {
    t0 <- proc.time()[[3L]]
    eval(expr, env)
    proc.time()[[3L]] - t0
  })
  sprintf("median %.3fs (nrep=%d)", stats::median(t), nrep)
}
cat_line("  dyad n=2001 filter+smoother+metric: ", local({
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 2001, seed = 1)
  time_one(aci(sim$x, model), 3L)
}))
cat_line("  CIR window 751 x horizon 1301: ",
         time_one(aci_cir(x, comp, dt, filt = filt,
                          window = window[rep], horizon = last_idx,
                          margin = 0.001), 1L))
cat_line("  vector online n=401 dim=2 full lag: ",
         time_one(aci_online_smoother(x_mv, comp_mv, dt_mv, filt_mv,
                                      lag = Inf), 3L))

# -- 9. MV test helper field name ---------------------------------------------

cat_line("=== 9. MV aux field names ===")
cat_line("  names(.aci_online_aux_mv) = ", paste(names(aux_mv), collapse = ", "))
cat_line("  names(.aci_online_aux)    = ", paste(names(aux), collapse = ", "))

# -- 10. Cholesky vs det route on a near-singular slice -----------------------

cat_line("=== 10. Strict PD vs pinv-style near-singular observation noise ===")
# A 2-d observation noise with condition number ~ 1e12. aciR should refuse
# the inverse; a pinv route would silently continue.
s_ill <- matrix(c(1, 1 - 1e-12, 1 - 1e-12, 1), 2L, 2L)
cat_line(sprintf(
  "  cond(S_xx) ~ %.3e  chol_ok=%s",
  kappa(s_ill, exact = TRUE),
  !is.null(tryCatch(chol(s_ill), error = function(e) NULL))
))
comp_ill <- list(
  L_x = diag(2L), L_y = -diag(2L),
  f_x = c(0, 0), f_y = c(0, 0),
  S_xoS_x = s_ill, S_yoS_y = diag(2L),
  S_yoS_x = matrix(0, 2L, 2L)
)
x_ill <- matrix(rnorm(20), 2L, 10L)
ill <- tryCatch(
  aci_filter(x_ill, comp_ill, 0.01, mu0 = c(0, 0), R0 = diag(2L)),
  error = function(e) conditionMessage(e)
)
cat_line("  aci_filter on ill-conditioned S_xx: ",
         if (is.character(ill)) paste("REFUSED:", substr(ill, 1L, 160L))
         else "ACCEPTED")

# Save a small table for the report
tab <- data.frame(
  index = window[rep],
  n_re = n_re[rep],
  peak_aci = peak_aci[rep],
  peak_matlab = mat$peak[rep],
  obj_aci = obj_aci[rep],
  obj_msimps = obj_msimps[rep],
  obj_matlab = mat$objective[rep],
  d_obj_native = d_obj_native[rep],
  d_obj_msimps = d_obj_msimps[rep],
  quad_diff = quad_diff[rep],
  settled_frac = settled_frac[rep],
  monotone = mono[rep]
)
utils::write.csv(
  tab,
  file.path(out_dir, "2026-08-14_cir_residual_rows.csv"),
  row.names = FALSE
)
cat_line("=== wrote 2026-08-14_cir_residual_rows.csv ===")
cat_line("DONE")
