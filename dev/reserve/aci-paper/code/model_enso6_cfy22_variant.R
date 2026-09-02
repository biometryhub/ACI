## acir reserve file
## Origin: aci/R/benchmark_models.R:513-654
## Source package: aci 0.0.30, git tree 97f6b124
## Category: aci-paper
## Intended release: 0.1.x, TBD with the supervisor/collaborators
## Reason: Whole aci model_enso6 body; the default 'cfy22' arm is a chen2022enso transcription, not the ACI_code MATLAB.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The whole aci 0.0.30 model_enso6 body, retained so the chen2022enso 'cfy22' transcription arms (parameter preset, c1f/c2f/beta0/zonal advection, coefficient provenance and initial condition) can be re-applied verbatim. The acir mainline keeps only variant 'aci_code'.

model_enso6 <- function(hidden = c("hW", "tau"), sigma_E = NULL,
                        lambda = 2/60, params = NULL,
                        variant = c("cfy22", "aci_code")) {
  variant <- match.arg(variant)
  vars <- c("u", "hW", "TC", "TE", "tau", "I")
  if (!is.character(hidden) || !length(hidden) || anyNA(hidden) ||
      any(!nzchar(hidden)) || anyDuplicated(hidden))
    aci_abort("aci_error_model_contract",
              "hidden must be a non-empty vector of unique variable names.")
  if (variant == "aci_code") {
    # Fixed-state drift and diffusion coefficients from ACI_code
    # ENSO_model_cond_ACI_*_unobs.m.  The generic R simulator remains
    # Euler-Maruyama, whereas the MATLAB scripts use a mixture of Euler and
    # Milstein updates; coefficient agreement therefore does not imply
    # pathwise or numerical-scheme parity.
    fct <- 0.65
    p <- params %||% list(gamma = 0.75 * fct, r = 0.25 * fct,
                          a1 = 0.0625 * fct^2, a2 = 0.125 * fct,
                          b0 = 2.5, mu = 0.5, Cu = 0.03 * fct, d_tau = 2,
                          s_u = 0.04 * sqrt(fct), s_h = 0.02 * sqrt(fct),
                          s_C = 0.04 * sqrt(fct))
    p$factor <- fct
    sigma_E <- sigma_E %||% (sqrt(5) * 1e-2 * sqrt(fct))
  } else {
    p <- params %||% list(gamma = 0.75, r = 0.25, a1 = 0.0625, a2 = 0.125,
                          b0 = 2.5, mu = 0.5, Cu = 0.03, d_tau = 2,
                          s_u = 0.04, s_h = 0.02, s_C = 0.04)
    sigma_E <- sigma_E %||% 0.02
  }
  hid <- match(hidden, vars)
  if (any(is.na(hid))) aci_abort("aci_error_model_contract",
    sprintf("hidden must be a subset of {%s}.", paste(vars, collapse = ", ")))
  if (any(vars[hid] %in% c("TC", "I")))
    aci_abort("aci_error_model_contract", paste(
      "TC and I must remain observed: c1(TC) TC is cubic in TC, sigma_tau",
      "depends on TC via tanh, and sigma_I(I) is multiplicative in I."))
  obs_i <- setdiff(seq_along(vars), hid)
  if (variant == "aci_code") {
    c1f <- function(TC, t)
      (25 * (TC + 0.75 / 7.5)^2 + 0.9) *
      (1 + 0.3 * sin(2 * pi * t / 6 - pi / 6)) * p$factor
    c2f <- function(t)
      1.4 * p$factor * (1 + 0.3 * sin(2 * pi * t / 6 + 2 * pi / 6) +
                          0.25 * sin(2 * (2 * pi * t / 6) + 2 * pi / 6))
    beta0 <- function(I) (1 + (1 - I / 5)) * 0.15 * sqrt(p$factor)
    zonal_advection <- function(I) I / 5 * p$factor
  } else {
    c1f <- function(TC, t)
      (25 * (TC + 0.1)^2 + 0.9) *
      (1 + 0.3 * sin(2 * pi * t / 6 - pi / 3))
    c2f <- function(t)
      1.4 * (1 + 0.2 * sin(2 * pi * t / 6 + pi / 3) +
               0.15 * sin(2 * pi * t / 3 + pi / 3))
    beta0 <- function(I) 0.15 * (2 - I / 5)
    zonal_advection <- function(I) 0.2 * I
  }
  stau <- function(TC, t) 0.9 * (tanh(7.5 * TC) + 1) * (1 + 0.3 * cos(2 * pi * t / 6 + pi / 3))
  # sigma_I from the uniform-stationary construction, with a small variance
  # floor so Gx stays SPD at the natural boundaries I in {0, 4}.  ACI_code
  # instead pseudo-inverts the observed Gramian (pinv, "for stability
  # concerns", ENSO_model_cond_ACI_h_W_unobs.m:1196-1197) and gives the I
  # channel zero precision where sigma_I vanishes (:1206-1208).  Inert for
  # assimilation: the I row of Lx is zero and gyx is zero, so this entry
  # multiplies zero in the gain.
  sI <- function(I) { Ic <- pmin(pmax(I, 0), 4)
    sqrt(pmax(lambda * Ic * (4 - Ic), 0) + 1e-3 * lambda) }
  drift_full <- function(t, s) {
    u <- s[1]; hW <- s[2]; TC <- s[3]; TE <- s[4]; tau <- s[5]; I <- s[6]
    bE <- beta0(I)
    c(-p$r * u  - (p$a1 * p$b0 * p$mu / 2) * (TC + TE) + (-0.2 * bE) * tau,
      -p$r * hW - (p$a2 * p$b0 * p$mu / 2) * (TC + TE) + (-0.4 * bE) * tau,
      (p$gamma * p$b0 * p$mu / 2 - c1f(TC, t)) * TC + (p$gamma * p$b0 * p$mu / 2) * TE +
        p$gamma * hW + zonal_advection(I) * u + p$Cu + (0.8 * bE) * tau,
      p$gamma * hW + (1.5 * p$gamma * p$b0 * p$mu - c2f(t)) * TE -
        (p$gamma * p$b0 * p$mu / 2) * TC + bE * tau,
      -p$d_tau * tau,
      -lambda * (I - 2))
  }
  sd_full <- function(t, s) c(p$s_u, p$s_h, p$s_C, sigma_E, stau(s[3], t), sI(s[6]))
  asm <- function(x, y) { s <- numeric(6); s[obs_i] <- x; s[hid] <- y; s }
  f_split <- function(t, x, y) drift_full(t, asm(x, y))[obs_i]
  g_split <- function(t, x, y) drift_full(t, asm(x, y))[hid]
  Sx <- function(t, x) { s <- numeric(6); s[obs_i] <- x
    diag(sd_full(t, s)[obs_i], length(obs_i)) }
  Sy_h <- function(t, x) { s <- numeric(6); s[obs_i] <- x
    diag(sd_full(t, s)[hid], length(hid), length(hid)) }
  m <- cgns_from_affine(f_split, g_split, Sx, Sy_h,
                        k = length(obs_i), l = length(hid),
                        name = sprintf("ENSO6[%s] (%s hidden)", variant,
                                       paste(vars[hid], collapse = ",")))
  m$meta$vars <- list(all = vars, hidden = vars[hid], observed = vars[obs_i])
  m$meta$params <- c(p, sigma_E = sigma_E, lambda = lambda)
  m$meta$variant <- variant
  matlab_hidden_sets <- list(
    u = "u", hW = "hW", tau = "tau", joint = c("u", "hW", "tau"))
  source_key <- names(Filter(function(z) identical(sort(z), sort(vars[hid])),
                             matlab_hidden_sets))
  m$meta$source_partition <- if (length(source_key)) source_key else NULL
  if (length(source_key)) {
    target_names <- if (identical(source_key, "joint")) c("TC", "TE", "I")
                    else c("TC", "TE", "I")
    target_names <- intersect(target_names, vars[obs_i])
    conditioning_names <- setdiff(vars[obs_i], target_names)
    m$meta$target_obs_idx <- match(target_names, vars[obs_i])
    m$meta$conditioning_obs_idx <- match(conditioning_names, vars[obs_i])
    m$meta$causal_link <- sprintf("(%s) -> (%s)%s",
      paste(vars[hid], collapse = ","), paste(target_names, collapse = ","),
      if (length(conditioning_names)) sprintf(" | (%s)", paste(conditioning_names, collapse = ",")) else "")
  } else {
    m$meta$source_partition <- "package_only_partition"
  }
  m$meta$coefficient_provenance <- if (variant == "aci_code")
    "ACI_code-main/ENSO_model_cond_ACI_*_unobs.m fixed-state coefficients" else
    "chen2022enso transcription"
  m$meta$numerical_regularization <- list(
    I_variance_floor = 1e-3 * lambda,
    reason = paste("The exact sigma_I vanishes at I=0 and I=4, whereas the",
                   "closed-form assimilation implementation requires a",
                   "non-singular observed-noise Gramian. The ACI_code scripts",
                   "instead pseudo-invert that Gramian and give the I channel",
                   "zero precision where sigma_I vanishes. Neither convention",
                   "reaches the hidden posterior here: the I row of Lx and the",
                   "noise cross-Gramian are both zero, so the I precision",
                   "multiplies zero in the filter gain. The floor does perturb",
                   "simulate()'s I path, which matlab_simulator_parity = FALSE",
                   "already records."))
  m$meta$simulation_convention <- paste(
    "simulate() uses Euler-Maruyama. The supplied MATLAB ENSO scripts use",
    "Euler updates for the interannual variables and Milstein-style updates",
    "for I and tau; no pathwise or numerical-scheme parity is claimed.")
  m$meta$matlab_simulator_parity <- FALSE
  m$meta$unsupported_partitions <- list(
    TC = paste("The MATLAB TC-hidden analysis uses a zeroth-order conditional",
               "approximation c1(t,TC) -> c1(t,0), not the exact six-state",
               "CGNS split constructed here."))
  ic_all <- if (variant == "aci_code")
    c(u = 6.9136e-04, hW = -0.0028, TC = 0.0039, TE = 0.0051,
      tau = -0.0256, I = 1.5841) else
    c(u = 0, hW = 0, TC = 0.1, TE = 0.1, tau = 0, I = 2)
  m$meta$ic_default <- list(x0 = ic_all[vars[obs_i]], y0 = ic_all[vars[hid]])
  m
}
