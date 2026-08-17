# validate_core.R, the aciR dyad ACI core (T3-T5), graded against the live oracle
# `oracle.aci_matlab_reference` (T8). Proves the research->oracle->spec->implement->
# grade loop closes GREEN: aciR's OWN filter/smoother/metric on the MATLAB input
# signal must reproduce the MATLAB reference outputs to numerical tolerance (IOP).
# These functions become aciR/R/ (the method core) after this passes.

# the ACI method core (grounded from substrate.md §7; dyad CGNS) ---------------

aci_dyad_components <- function(x, p) {
  list(L_x = p$gamma * x,
       f_x = p$F_x - p$d_x * x,
       L_y = -p$d_y,                       # scalar (constant in time)
       f_y = p$F_y - p$gamma * x^2,
       S_xoS_x = p$sigma_x^2,
       S_yoS_y = p$sigma_y^2,
       S_yoS_x = 0, S_xoS_y = 0)           # dyad: no obs/latent noise cross-term
}

aci_filter <- function(x, comp, dt, mu0, R0) {   # forward CGNS filter (G1.1)
  n <- length(x); m <- numeric(n); v <- numeric(n)
  m[1] <- mu0; v[1] <- R0
  inv <- 1 / comp$S_xoS_x; mu <- mu0; R <- R0
  for (j in 2:n) {
    dx  <- x[j] - x[j - 1]
    aux <- comp$S_yoS_x + R * comp$L_x[j - 1]
    mu  <- mu + (comp$L_y * mu + comp$f_y[j - 1]) * dt +
           aux * inv * (dx - (comp$L_x[j - 1] * mu + comp$f_x[j - 1]) * dt)
    R   <- R + (2 * comp$L_y * R + comp$S_yoS_y - aux * inv * aux) * dt
    m[j] <- mu; v[j] <- R
  }
  list(mean = m, cov = v)
}

aci_smoother <- function(x, comp, dt, filt) {    # backward CGNS smoother (G1.2)
  n <- length(x); m <- numeric(n); v <- numeric(n)
  m[n] <- filt$mean[n]; v[n] <- filt$cov[n]
  inv <- 1 / comp$S_xoS_x; muT <- m[n]; RT <- v[n]
  for (j in (n - 1):1) {
    dx  <- x[j + 1] - x[j]
    A_j <- comp$L_y - comp$S_yoS_x * inv * comp$L_x[j]
    B_j <- comp$S_yoS_y - comp$S_yoS_x * inv * comp$S_xoS_y
    muT <- muT - (comp$L_y * muT + comp$f_y[j] -
                  B_j / filt$cov[j] * (filt$mean[j] - muT)) * dt +
           comp$S_yoS_x * inv * (-dx + (comp$L_x[j] * muT + comp$f_x[j]) * dt)
    RT  <- RT - (2 * (A_j + B_j / filt$cov[j]) * RT - B_j) * dt
    m[j] <- muT; v[j] <- RT
  }
  list(mean = m, cov = v)
}

aci_metric <- function(filt, smooth) {           # ACI relative-entropy metric (G1.3)
  signal     <- 0.5 * (smooth$mean - filt$mean)^2 / filt$cov
  cov_ratio  <- smooth$cov / filt$cov
  dispersion <- 0.5 * (-log(cov_ratio) + cov_ratio - 1)
  signal + dispersion
}

# T8: grade against oracle.aci_matlab_reference --------------------------------

sig <- read.csv("dyad_signal_x.csv", header = FALSE)   # cols: t, x (no header)
x   <- sig$V2
p   <- list(d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
            sigma_x = 0.5, sigma_y = 1)
dt  <- 0.001
comp <- aci_dyad_components(x, p)
filt   <- aci_filter(x, comp, dt, mu0 = p$F_y / p$d_y, R0 = 0.1)
smooth <- aci_smoother(x, comp, dt, filt)
aci    <- aci_metric(filt, smooth)

ref <- read.csv("dyad_reference.csv")                  # MATLAB expected (301 pts)
idx <- seq(1, length(x), by = 100)                     # matches MATLAB 1:100:N+1
stopifnot(length(idx) == nrow(ref))

err <- c(
  filter_mean   = max(abs(filt$mean[idx]   - ref$filter_mean)),
  filter_cov    = max(abs(filt$cov[idx]    - ref$filter_cov)),
  smoother_mean = max(abs(smooth$mean[idx] - ref$smoother_mean)),
  smoother_cov  = max(abs(smooth$cov[idx]  - ref$smoother_cov)),
  ACI_metric    = max(abs(aci[idx]         - ref$ACI_metric))
)
cat("max abs error vs MATLAB oracle (301 sampled points):\n")
for (k in names(err)) cat(sprintf("  %-14s %.3e\n", k, err[k]))
tol <- 1e-6
if (max(err) < tol) {
  cat(sprintf("\nORACLE-VALIDATED: GREEN, aciR core matches the MATLAB reference to < %.0e\n", tol))
  cat(sprintf("  (checksum sum(ACI)=%.10f vs registry 9933.0774195694)\n", sum(aci)))
} else {
  cat("\nMISMATCH, the reimplementation diverges from the oracle; investigate before the spec proceeds.\n")
  quit(status = 1)
}
