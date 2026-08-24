# Faithful R port of the published EnKBS dyad experiment (jiang2026enkbs;
# EnKBS-main/dyad/utov.m: v observed, u hidden, so u -> v is the direction
# under study). The port mirrors the MATLAB loops statement by statement and
# is driven by supplied noise, so the package's ensemble engine can be graded
# against it at machine precision on identical inputs. MATLAB and R random
# streams differ, which is why the grade drives both sides with one set of
# increments rather than replaying the published run itself.
#
# One deliberate difference of the package is folded in here: the reference
# adds its observation perturbation (+ sigma_v dW), while the package follows
# the paper's innovation (dy - dy_i) and subtracts it. The two are identical
# in distribution; the test passes the package the negated draws so that both
# sides see the same realised innovations.

# Truth trajectory of the dyad by the reference's Euler scheme.
golden_p3_truth <- function(n, dt, p, z_u, z_v, u0 = 0.9, v0 = -0.2) {
  u <- numeric(n); v <- numeric(n)
  u[1] <- u0; v[1] <- v0
  for (i in 2:n) {
    u[i] <- u[i - 1] + ((-p$d_u + p$c * v[i - 1]) * u[i - 1] + p$F_u) * dt +
      sqrt(dt) * z_u[i - 1] * p$sigma_u
    v[i] <- v[i - 1] + (-p$d_v * v[i - 1] - p$c * u[i - 1]^2 + p$F_v) * dt +
      sqrt(dt) * z_v[i - 1] * p$sigma_v
  }
  list(u = u, v = v)
}

# Forward EnKBF: returns the m x n matrix of filter members.
golden_p3_enkbf <- function(u0_ens, v_truth, z_u_ens, z_v_obs, p, dt) {
  n <- length(v_truth); m <- length(u0_ens)
  filt <- matrix(NA_real_, m, n)
  u_ens <- u0_ens
  filt[, 1] <- u_ens
  for (i in 2:n) {
    v_prev <- v_truth[i - 1]
    d_y <- v_truth[i] - v_prev
    f_vals <- (-p$d_u + p$c * v_prev) * u_ens + p$F_u
    h_vals <- -p$d_v * v_prev - p$c * u_ens^2 + p$F_v
    du <- u_ens - mean(u_ens)
    dh <- h_vals - mean(h_vals)
    p_uh <- sum(du * dh) / (m - 1)
    gain <- p_uh / p$sigma_v^2
    u_pred <- u_ens + f_vals * dt + p$sigma_u * sqrt(dt) * z_u_ens[, i - 1]
    innovation <- d_y - h_vals * dt + p$sigma_v * sqrt(dt) * z_v_obs[, i - 1]
    u_ens <- u_pred + gain * innovation
    filt[, i] <- u_ens
  }
  filt
}

# One backward EnKBS pass from horizon `k`, reusing the forward increments.
# Returns the m x k matrix of smoother members over steps 1..k.
golden_p3_enkbs <- function(filt, v_truth, z_u_ens, p, dt, k = ncol(filt)) {
  m <- nrow(filt)
  smoo <- matrix(NA_real_, m, k)
  u_s <- filt[, k]
  smoo[, k] <- u_s
  if (k > 1) for (i in (k - 1):1) {
    v_right <- v_truth[i + 1]
    p_right <- stats::var(filt[, i + 1])
    f_s <- (-p$d_u + p$c * v_right) * u_s + p$F_u
    u_s <- u_s - f_s * dt - p$sigma_u * sqrt(dt) * z_u_ens[, i] -
      (p$sigma_u^2 / p_right) * (u_s - filt[, i + 1]) * dt
    smoo[, i] <- u_s
  }
  smoo
}

# The lagged triangle: one backward pass per horizon, storing the smoother
# mean and variance at (step, horizon), then the relative entropy of the
# complete-record posterior from each partial one.
golden_p3_delta <- function(filt, v_truth, z_u_ens, p, dt) {
  n <- ncol(filt)
  tri_mean <- matrix(NA_real_, n, n)
  tri_var <- matrix(NA_real_, n, n)
  for (k in seq_len(n)) {
    smoo <- golden_p3_enkbs(filt, v_truth, z_u_ens, p, dt, k)
    tri_mean[seq_len(k), k] <- colMeans(smoo)
    tri_var[seq_len(k), k] <- apply(smoo, 2, stats::var)
  }
  delta <- matrix(NA_real_, n, n)
  for (k in seq_len(n)) {
    j <- seq_len(k)
    ratio <- tri_var[j, n] / tri_var[j, k]
    delta[j, k] <- 0.5 * ((tri_mean[j, n] - tri_mean[j, k])^2 / tri_var[j, k] +
                            ratio - 1 - log(ratio))
  }
  list(mean = tri_mean, var = tri_var, delta = delta)
}

# The reference's approximate objective forward CIR from the delta table:
# Simpson ratio per anchor, zeroed below the strength threshold.
golden_p3_cir <- function(delta, dt, threshold = 1e-5) {
  n <- nrow(delta)
  vapply(seq_len(n), function(r) {
    re <- delta[r, r:n]
    if (max(re) <= threshold) return(0)
    .simps_u(re) * dt / max(re)
  }, numeric(1))
}
