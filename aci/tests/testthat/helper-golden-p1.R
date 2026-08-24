# Faithful R port of the P1 golden reference (marandmath/ACI_code,
# dyad_interaction_model.m; MIT). Implements the MATLAB loops verbatim for
# the uncorrelated-noise scalar-hidden case so package output can be compared
# against the published reference on identical observation series.
#
# Convention notes (documented deltas vs this package):
#  * filter: identical scheme (coefficients at the interval start, explicit
#    Euler, gain (Syx + R Lx)/Sxx) -- expected to agree to machine precision.
#  * smoother: both the golden and the package use the left/source endpoint of
#    each backward Euler interval; they should agree to machine precision.
#  * ACI metric: identical signal/dispersion split of the Gaussian KL.
#  * "definition" objective CIR: eps-integral of subjective CIRs -- the
#    layer-cake identity this package's exact estimator computes directly.
#  * approx objective CIR: Simpson L1/Linf ratio. The package's `l1_linf`
#    implements the same continuum ratio with positive-lag rectangle cells.

golden_p1_moments <- function(x, dt, prm, mu0, R0 = 0.1) {
  N1 <- length(x)
  Lx <- prm$gamma * x
  fx <- prm$F_x - prm$d_x * x
  fy <- prm$F_y - prm$gamma * x^2
  Ly <- -prm$d_y
  Sxx <- prm$s_x^2; Syy <- prm$s_y^2   # S_yx = 0 in the golden dyad
  fm <- fc <- numeric(N1)
  fm[1] <- mu0; fc[1] <- R0
  for (j in 2:N1) {
    dx  <- x[j] - x[j - 1]
    aux <- fc[j - 1] * Lx[j - 1]
    fm[j] <- fm[j - 1] + (Ly * fm[j - 1] + fy[j - 1]) * dt +
             aux / Sxx * (dx - (Lx[j - 1] * fm[j - 1] + fx[j - 1]) * dt)
    fc[j] <- fc[j - 1] + (2 * Ly * fc[j - 1] + Syy - aux^2 / Sxx) * dt
  }
  sm <- sc <- numeric(N1)
  sm[N1] <- fm[N1]; sc[N1] <- fc[N1]
  for (j in (N1 - 1):1) {              # target-index convention (golden)
    B <- Syy                           # A_j = Ly, B_j = Syy when S_yx = 0
    sm[j] <- sm[j + 1] - (Ly * sm[j + 1] + fy[j] -
                          B / fc[j] * (fm[j] - sm[j + 1])) * dt
    sc[j] <- sc[j + 1] - (2 * (Ly + B / fc[j]) * sc[j + 1] - B) * dt
  }
  sig <- 0.5 * (sm - fm)^2 / fc
  cr  <- sc / fc
  dis <- 0.5 * (-log(cr) + cr - 1)
  list(fm = fm, fc = fc, sm = sm, sc = sc,
       aci = sig + dis, signal = sig, dispersion = dis,
       Lx = Lx, fx = fx, fy = fy, Ly = Ly, Sxx = Sxx, Syy = Syy, dt = dt)
}

# Composite Simpson's 1/3 rule, uniform grid, spacing 1 (matches the golden
# usage simps(y) * dt); trapezoid on a trailing odd interval.
# Independent transcription of ACI_code-main/simps.m on a unit grid: composite
# Simpson over the leading odd-length block, then, when the point count is
# even, MATLAB's `C = vander(x(end-2:end)) \ y(end-2:end)` quadratic fit
# through the last three points, integrated over the final interval as
# C1 (x3^3 - x2^3)/3 + C2 (x3^2 - x2^2)/2 + C3 (x3 - x2).  The Vandermonde
# solve is kept explicit so this stays an independent check of the package's
# closed form rather than a copy of it.
.simps_u <- function(y) {
  n <- length(y)
  if (n < 3) return(if (n == 2) 0.5 * sum(y) else 0)
  m <- if (n %% 2 == 1) n else n - 1
  i <- seq(1, m - 2, by = 2)
  s <- sum(y[i] + 4 * y[i + 1] + y[i + 2]) / 3
  if (n %% 2 == 0) {
    x <- c(n - 2, n - 1, n)
    C <- solve(cbind(x^2, x, 1), y[c(n - 2, n - 1, n)])
    x2 <- x[2]; x3 <- x[3]
    s <- s + C[1] * (x3^3 - x2^3) / 3 + C[2] * (x3^2 - x2^2) / 2 +
      C[3] * (x3 - x2)
  }
  s
}

# Composite Simpson on a possibly non-uniform grid (parabola per interval
# pair; trapezoid on a trailing odd interval) -- the golden's simps(x, y).
.simps_xy <- function(x, y) {
  n <- length(x); s <- 0; i <- 1
  while (i + 2 <= n) {
    h1 <- x[i + 1] - x[i]; h2 <- x[i + 2] - x[i + 1]; h <- h1 + h2
    s <- s + h / 6 * ((2 - h2 / h1) * y[i] +
                      h^2 / (h1 * h2) * y[i + 1] +
                      (2 - h1 / h2) * y[i + 2])
    i <- i + 2
  }
  if (i + 1 == n) s <- s + 0.5 * (x[n] - x[n - 1]) * (y[n - 1] + y[n])
  s
}

# The golden online fixed-lag smoother (full lag), delta matrix, and the two
# objective-CIR readouts, ported verbatim. Returns, for j in 1..N1:
#   RE[j, n] = delta(T' = t_n; t_j) for n >= j,
#   approx_obj (Simpson L1/Linf ratio) and defn_obj (eps-integral of the
#   subjective CIRs), both on the full window.
golden_p1_cir <- function(g, x, anchors = NULL) {
  N1 <- length(x); dt <- g$dt
  fm <- g$fm; fc <- g$fc; Lx <- g$Lx; fx <- g$fx; fy <- g$fy
  Ly <- g$Ly; Syy <- g$Syy; Sxxi <- 1 / g$Sxx
  Ev <- Fv <- numeric(N1)
  for (j in 1:N1) {
    Gx  <- Lx[j]                        # + S_xy / fc = 0
    Gy  <- Ly + Syy / fc[j]
    Cjj <- 1 - Gy * dt
    H   <- (2 * Ly * fc[j] + Syy) / fc[j]
    K   <- Sxxi * Gx
    Ev[j] <- Cjj                        # + S_yx K dt = 0
    Fv[j] <- -fc[j] * (K + (Gx * K * fc[j] * K - H * K + Ly * K) * dt -
                       Lx[j] * (Sxxi + K * fc[j] * K * dt))
  }
  OM <- OC <- matrix(NA_real_, N1, N1)  # [n, j], n >= j
  OM[1, 1] <- fm[1]; OC[1, 1] <- fc[1]
  step_one <- function(n) {             # p_n(y^{n-1}) from p_n-lineage
    av <- fm[n - 1] -
          Ev[n - 1] * ((1 + Ly * dt) * fm[n - 1] + fy[n - 1] * dt) +
          Fv[n - 1] * (x[n] - x[n - 1] -
                       (Lx[n - 1] * fm[n - 1] + fx[n - 1]) * dt)
    am <- fc[n - 1] - Ev[n - 1] * (1 + Ly * dt) * fc[n - 1] -
          Fv[n - 1] * Lx[n - 1] * fc[n - 1] * dt
    c(Ev[n - 1] * fm[n] + av, Ev[n - 1]^2 * fc[n] + am)
  }
  OM[2, 2] <- fm[2]; OC[2, 2] <- fc[2]
  s1 <- step_one(2); OM[2, 1] <- s1[1]; OC[2, 1] <- s1[2]
  U <- matrix(NA_real_, N1, N1)         # update products, [n-2 lineage, j]
  for (n in 3:N1) {
    OM[n, n] <- fm[n]; OC[n, n] <- fc[n]
    s1 <- step_one(n); OM[n, n - 1] <- s1[1]; OC[n, n - 1] <- s1[2]
    U[n - 2, n - 1] <- 1
    U[n - 2, n - 2] <- Ev[n - 2]
    if (n - 3 >= 1)
      U[n - 2, 1:(n - 3)] <- U[n - 3, 1:(n - 3)] * Ev[n - 2]
    inM <- OM[n, n - 1] - fm[n - 1]
    inC <- OC[n, n - 1] - fc[n - 1]
    jj <- 1:(n - 2)
    OM[n, jj] <- OM[n - 1, jj] + U[n - 2, jj] * inM
    OC[n, jj] <- OC[n - 1, jj] + U[n - 2, jj]^2 * inC
  }
  RE <- matrix(NA_real_, N1, N1)        # [j, n]
  for (j in 1:N1) {
    n <- j:N1
    cr <- OC[N1, j] / OC[n, j]
    RE[j, n] <- 0.5 * (OM[N1, j] - OM[n, j])^2 / OC[n, j] +
                0.5 * (cr - 1 - log(cr))
  }
  maxRE <- apply(RE, 1, max, na.rm = TRUE)
  approx_obj <- vapply(1:N1, function(j) {
    v <- RE[j, j:N1]
    if (maxRE[j] > 1e-5) .simps_u(v) * dt / maxRE[j] else 0
  }, numeric(1))
  if (is.null(anchors)) anchors <- 1:N1
  eps_ord <- rev(seq(-6, 0.5, length.out = 513))
  subj <- matrix(0, length(eps_ord), length(anchors))
  for (ai in seq_along(anchors)) {
    v <- RE[anchors[ai], anchors[ai]:N1]
    sm_ <- rev(cummax(rev(v)))          # suffix maxima: last-exceed is where
    for (ei in seq_along(eps_ord)) {    # the suffix max still clears eps
      eps <- 10^eps_ord[ei]
      idx <- which(sm_ > eps)
      subj[ei, ai] <- if (length(idx)) max(idx) * dt else 0
    }
  }
  eps_asc <- 10^rev(eps_ord)
  defn_obj <- rep(NA_real_, N1)
  for (ai in seq_along(anchors))
    defn_obj[anchors[ai]] <- .simps_xy(eps_asc, rev(subj[, ai])) /
                             maxRE[anchors[ai]]
  list(RE = RE, maxRE = maxRE, approx_obj = approx_obj,
       defn_obj = defn_obj, subj = subj, eps_ord = eps_ord,
       anchors = anchors, E = Ev, F = Fv)
}
