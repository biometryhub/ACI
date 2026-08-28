################################################################################
## aci-model-library.R - benchmark-model constructors
## ########################################################################## ##
##
## Contents:
##   * benchmark models from 'MATLAB Codebase for "Assimilative Causal Inference"':
##       - aci_dyad_model, aci_enso_model, aci_predprey_model
##   * the TC-hidden zeroth-order ENSO branch:
##       - .enso6_prescribed_grid, .enso6_tc_zeroth_order
##
################################################################################


################################################################################
# benchmark models with transcription provenance
################################################################################

#' Complete and validate a scalar parameter list (internal)
#'
#' @param params Optional named list of parameter overrides; `NULL` uses the
#'   defaults unchanged.
#' @param defaults Named list of default parameter values.
#' @param constructor 1-length character naming the constructor in error
#'   messages.
#' @param positive Character vector of parameter names required to be positive.
#' @returns The completed named list, ordered as `defaults`.
#' @noRd
.complete_scalar_params <- function(params, defaults, constructor,
                                    positive = character()) {
  if (is.null(params)) return(defaults)
  if (!is.list(params) || is.null(names(params)) || any(!nzchar(names(params))) ||
      anyDuplicated(names(params)))
    aci_abort("aci_error_model_contract",
              sprintf("%s params must be a uniquely named list.", constructor))
  unknown <- setdiff(names(params), names(defaults))
  missing <- setdiff(names(defaults), names(params))
  if (length(unknown))
    aci_abort("aci_error_model_contract",
              sprintf("Unknown %s parameter(s): %s.", constructor,
                      paste(unknown, collapse = ", ")))
  if (length(missing))
    aci_abort("aci_error_model_contract",
              sprintf("Missing %s parameter(s): %s.", constructor,
                      paste(missing, collapse = ", ")))
  p <- params[names(defaults)]
  good_scalar <- vapply(p, function(z)
    is.numeric(z) && length(z) == 1L && is.finite(z), logical(1))
  if (!all(good_scalar))
    aci_abort("aci_error_model_contract",
              sprintf("All %s parameters must be finite numeric scalars.",
                      constructor))
  if (length(positive) && any(unlist(p[positive], use.names = FALSE) <= 0))
    aci_abort("aci_error_model_contract",
              sprintf("%s parameter(s) must be positive: %s.", constructor,
                      paste(positive, collapse = ", ")))
  p
}


#' Nonlinear dyad (andreou2026aci eq. 1-2)
#'
#' @param variant Paper-specific parameter preset.
#' @param observe Which dyad component is treated as observed.
#' @param params Optional complete parameter list overriding the preset.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications* **17**, 1854. \doi{10.1038/s41467-026-68568-0}
#' @examples
#' aci_dyad_model()
#'
#' @export
aci_dyad_model <- function(variant = "p1", observe = "x", params = NULL) {
  variant <- match.arg(variant); observe <- match.arg(observe)
  defaults <- switch(variant,
    p1 = list(d_x = 0.5, gamma = 2,   f_x = 0.5, s_x = 0.5, d_y = 0.5, f_y = 1,   s_y = 1))
  p <- .complete_scalar_params(params, defaults, "dyad",
                               c("d_x", "gamma", "s_x", "d_y", "s_y"))
  if (observe == "x") {
    # The observed-x dyad has a package batch realiser.  Its coefficient
    # closures share a locked environment so the attached realiser descriptor
    # cannot silently outlive mutation of captured constructor parameters.
    coefficient_env <- list2env(list(p = p), parent = baseenv())
    lockEnvironment(coefficient_env, bindings = TRUE)
    coefficient_functions <- list(
      Lx = function(t, x) matrix(p$gamma * x, 1, 1),
      fx = function(t, x) -p$d_x * x + p$f_x,
      Ly = function(t, x) matrix(-p$d_y, 1, 1),
      fy = function(t, x) -p$gamma * x^2 + p$f_y,
      Sx1 = function(t, x) matrix(p$s_x, 1, 1),
      Sx2 = function(t, x) matrix(0, 1, 1),
      Sy1 = function(t, x) matrix(0, 1, 1),
      Sy2 = function(t, x) matrix(p$s_y, 1, 1)
    )
    coefficient_functions <- lapply(coefficient_functions, function(fun) {
      environment(fun) <- coefficient_env
      fun
    })
    m <- aci_model(
      Lx = coefficient_functions$Lx,
      fx = coefficient_functions$fx,
      Ly = coefficient_functions$Ly,
      fy = coefficient_functions$fy,
      Sx1 = coefficient_functions$Sx1,
      Sx2 = coefficient_functions$Sx2,
      Sy1 = coefficient_functions$Sy1,
      Sy2 = coefficient_functions$Sy2,
      k = 1, l = 1, name = sprintf("dyad[%s] y->x", variant))
  }
  m$meta$energy_conserving <- TRUE
  m$meta$params <- p
  m$meta$vars <- list(observed = "x", hidden = "y")
  m$meta$provenance <- switch(variant,
    p1 = paste("andreou2026aci Sections 3.1 and SI.4.1;",
               "ACI_code-main/dyad_interaction_model.m"))
  m$meta$source_status <- "paper + MATLAB checked"
  m$meta$anti_damping_threshold <- p$d_x / p$gamma
  # The MATLAB reference starts each component at its uncoupled forced
  # equilibrium (F_x / d_x, F_y / d_y).  In particular, the default
  # andreou2026aci configuration starts at (1, 2), not (1, 0).
  m$meta$ic_default <- list(x0 = p$f_x / p$d_x, y0 = p$f_y / p$d_y)
  m <- .attach_cgns_realizer(
    m, "dyad_observed_x_v1", list(params = p)
  )
  m
}


#' Six-variable stochastic conceptual ENSO model (andreou2026aci SI via
#' chen2022enso eqs. 1a-1f).
#'
#' Constructs the six-variable ENSO system split into an observed and a hidden
#' component at `hidden`. `I` must remain observed, and so must `TC` unless the
#' zeroth-order approximation below is requested. The coefficients are the
#' fixed-state ones of the ACI_code conditional ENSO scripts.
#'
#' # T_C hidden (zeroth-order)
#'
#' The six-state ENSO system is not conditionally Gaussian with `T_C` hidden:
#' the damping `c_1(t, T_C) T_C` is cubic in `T_C`. `ACI_code`'s
#' `ENSO_model_cond_ACI_T_C_unobs.m` restores conditional linearity by a
#' zeroth-order Taylor expansion of `c_1` about the climatology `T_C = 0`,
#' replacing the state-dependent damping with the time-only series
#' `r_C - c_1(t, 0)`. `approximation = "zeroth_order_c1"` builds that inference
#' model: observed `(T_E, I)`, hidden `T_C`, with `u`, `h_W` and `tau` entering
#' as prescribed forcings from their observed series, supplied through
#' `prescribed`. This is an approximation of the system, not a re-split of it -
#' the simulator keeps the full nonlinear `c_1(t, T_C)`, and the filter and
#' smoother moments are those of the approximating model. `simulate()` is
#' refused on the result for that reason; generate a path from
#' `aci_enso_model(hidden = c("u", "hW", "tau"))` and build this model from it.
#'
#' The prescribed series are looked up by index on the grid they were supplied
#' on, never interpolated, and assimilation refuses an observation grid that is
#' not that grid.
#'
#' # Observation set for the tau partition
#'
#' `hidden = "tau"` has two estimands, and they are not the same causal
#' quantity. `observations = "reduced"`, the default for that partition, is the
#' reference script's: the observed process is `(T_C, T_E, I)`, and `u` and
#' `h_W` enter the target drifts as prescribed known time series rather than as
#' assimilated channels (`ENSO_model_cond_ACI_tau_unobs.m:1020-1039`,
#' `:1136-1141`). `observations = "full"` assimilates all five observed
#' channels `(u, h_W, T_C, T_E, I)` and so uses strictly more information: both
#' prescribed drifts carry `tau` (`Lx[u] = -0.0407`, `Lx[h_W] = -0.0814`), and
#' prescribing them reproduces their effect on the `T_C`/`T_E` drift but not
#' their own innovations.
#'
#' The script asserts the two agree. On a 4001-point source-derived path they
#' do not: filter means differ by up to 0.247, the ACI series by up to 0.776 -
#' about three times its own mean level, with Pearson correlation 0.905 -
#' while the time-averaged ACI agrees to within 0.5%. The reduction preserves
#' the average level and distorts the time-resolved curve, which is the
#' reported quantity. Any fidelity claim must name which observation set it
#' reproduces. `meta$observations` records which one a model carries.
#'
#' @param hidden Character vector naming hidden ENSO variables.
#' @param sigma_E Eastern-Pacific temperature noise amplitude.
#' @param lambda Decay rate for the diversity index.
#' @param params Optional complete parameter list overriding the preset.
#' @param variant Parameter and coefficient convention.
#' @param observations Observation set the estimand is defined on, `"reduced"`
#'   or `"full"`. `"reduced"` is defined for `hidden = "tau"` and
#'   `hidden = "TC"`, where it is the default and, for `"TC"`, the only value;
#'   every other partition is `"full"`. For `hidden = "tau"` observations are
#'   still supplied on all five observed channels either way, and `"reduced"`
#'   prescribes `u` and `h_W` from them instead of assimilating them.
#' @param approximation Either `"exact"`, the conditionally linear split of the
#'   six-state drift, or `"zeroth_order_c1"`, the `T_C`-hidden substitution
#'   `c1(t, TC) -> c1(t, 0)` described above. `"zeroth_order_c1"` is defined
#'   only for `hidden = "TC"`, and `hidden = "TC"` requires it.
#' @param prescribed For `hidden = "TC"` only, and then required: a data frame
#'   or named list carrying `t`, `u`, `hW` and `tau` as equal-length numeric
#'   series on one strictly increasing uniform time grid. Other elements, such
#'   as the remaining channels of a simulated path, are ignored.
#' @param matlab_defect_compat For `hidden = "TC"` only. The reference script's
#'   assimilation forcing `f_y` (`ENSO_model_cond_ACI_T_C_unobs.m:1053,:1151`)
#'   omits the thermocline term `gamma_C * h_W`, which the same script's
#'   simulator drift (`:1124`) includes and which the sibling `ACI_code`
#'   scripts carry in the corresponding `T_C` coefficient rows. `h_W` is
#'   prescribed and observed here, so the term is available, and `acir`
#'   includes it. Set `TRUE` to reproduce the published script verbatim. On a
#'   14-model-year window the two differ by up to `0.105` in the filter mean
#'   and `2.75` in ACI, and the time-integrated ACI roughly doubles; filter and
#'   smoother covariances are identical, because the term enters only the mean
#'   equations.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications* **17**, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' Chen, N., Fang, X. and Yu, J.-Y. (2022). A multiscale model for El Nino
#' complexity. *npj Climate and Atmospheric Science* **5**, 16. arXiv:2104.07174.
#' @examples
#' aci_enso_model()
#'
#' # The T_C-hidden partition is not self-contained: it needs the u, h_W and
#' # tau series it treats as prescribed forcings.
#' sim <- simulate(aci_enso_model(hidden = c("u", "hW", "tau")),
#'                 seed = 1, T = 1, dt = 0.005)
#' path <- data.frame(t = sim$obs$t, u = sim$hidden[, 1],
#'                    hW = sim$hidden[, 2], tau = sim$hidden[, 3])
#' aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1",
#'             prescribed = path)
#'
#' @export
aci_enso_model <- function(hidden = c("hW", "tau"), sigma_E = NULL,
                           lambda = 2/60, params = NULL,
                           variant = "aci_code",
                           observations = c("reduced", "full"),
                           approximation = c("exact", "zeroth_order_c1"),
                           prescribed = NULL, matlab_defect_compat = FALSE) {
  variant <- match.arg(variant)
  approximation <- match.arg(approximation)
  observations_given <- !missing(observations)
  if (observations_given) observations <- match.arg(observations)
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
  }
  hid <- match(hidden, vars)
  if (any(is.na(hid))) aci_abort("aci_error_model_contract",
    sprintf("hidden must be a subset of {%s}.", paste(vars, collapse = ", ")))
  ## TC alone is admissible under the zeroth-order substitution and nothing
  ## else is: no MATLAB script and no closed form covers a hidden set that
  ## pairs TC with another variable, and I is multiplicative in its own noise.
  tc_only <- identical(vars[hid], "TC")
  if ("I" %in% vars[hid] || ("TC" %in% vars[hid] && !tc_only))
    aci_abort("aci_error_model_contract", paste(
      "TC and I must remain observed: c1(TC) TC is cubic in TC, sigma_tau",
      "depends on TC via tanh, and sigma_I(I) is multiplicative in I."))
  if (tc_only && !identical(approximation, "zeroth_order_c1"))
    aci_abort("aci_error_model_contract", paste(
      "hidden = 'TC' is conditionally linear only under the zeroth-order",
      "substitution c1(t, TC) -> c1(t, 0); pass",
      "approximation = 'zeroth_order_c1'."))
  if (!tc_only && identical(approximation, "zeroth_order_c1"))
    aci_abort("aci_error_model_contract", paste(
      "approximation = 'zeroth_order_c1' is the TC-hidden construction of",
      "ENSO_model_cond_ACI_T_C_unobs.m; the other partitions split the drift",
      "exactly and have nothing to approximate."))
  if (length(matlab_defect_compat) != 1L || !is.logical(matlab_defect_compat) ||
      is.na(matlab_defect_compat))
    aci_abort("aci_error_model_contract",
              "matlab_defect_compat must be TRUE or FALSE.")
  if (isTRUE(matlab_defect_compat) && !tc_only)
    aci_abort("aci_error_model_contract", paste(
      "matlab_defect_compat reproduces one source defect of",
      "ENSO_model_cond_ACI_T_C_unobs.m:1053,:1151 and applies only to",
      "hidden = 'TC'."))
  if (!is.null(prescribed) && !tc_only)
    aci_abort("aci_error_model_contract", paste(
      "prescribed supplies the deterministic u, hW and tau forcings of the",
      "TC-hidden reduced model and applies only to hidden = 'TC'. The tau",
      "partition's own reduction reads its prescribed channels from the",
      "supplied observations; see observations = 'reduced'."))
  ## The reduced observation set is the tau and TC scripts', and only theirs.
  tau_only <- identical(vars[hid], "tau")
  if (!observations_given)
    observations <- if (tau_only || tc_only) "reduced" else "full"
  if (identical(observations, "reduced") && !tau_only && !tc_only)
    aci_abort("aci_error_model_contract", paste(
      "observations = 'reduced' is defined only for hidden = 'tau' and",
      "hidden = 'TC', where ACI_code assimilates a subset of the observed",
      "channels and prescribes the rest",
      "(ENSO_model_cond_ACI_tau_unobs.m:1020-1039, :1136-1141;",
      "ENSO_model_cond_ACI_T_C_unobs.m:999-1021)."))
  if (identical(observations, "full") && tc_only)
    aci_abort("aci_error_model_contract", paste(
      "hidden = 'TC' has no five-channel estimand: the zeroth-order",
      "construction assimilates (TE, I) and prescribes u, hW and tau",
      "(ENSO_model_cond_ACI_T_C_unobs.m:999-1021)."))
  if (tc_only)
    return(.enso6_tc_zeroth_order(p, sigma_E, lambda, variant, vars,
                                  prescribed, matlab_defect_compat))
  obs_i <- setdiff(seq_along(vars), hid)
  # Every ENSO partition has a package batch realiser, so the coefficient
  # closures and the whole-path expressions it evaluates share one locked
  # environment: the attached realiser descriptor cannot outlive mutation of
  # the captured constructor parameters, and the two routes read the same
  # bindings or neither of them runs.
  ce <- list2env(list(p = p, sigma_E = sigma_E, lambda = lambda,
                      hid = hid, obs_i = obs_i, variant = variant),
                 parent = baseenv())
  local({
    if (variant == "aci_code") {
      c1f <- function(TC, t)
        (25 * (TC + 0.75 / 7.5)^2 + 0.9) *
        (1 + 0.3 * sin(2 * pi * t / 6 - pi / 6)) * p$factor
      c2f <- function(t)
        1.4 * p$factor * (1 + 0.3 * sin(2 * pi * t / 6 + 2 * pi / 6) +
                            0.25 * sin(2 * (2 * pi * t / 6) + 2 * pi / 6))
      beta0 <- function(I) (1 + (1 - I / 5)) * 0.15 * sqrt(p$factor)
      zonal_advection <- function(I) I / 5 * p$factor
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
    # drift_full() and sd_full() hold elementwise, so `t` and each `s[[i]]` may
    # be one point, as the per-point coefficient contract supplies them, or a
    # whole observation path, as the batch realiser does.  The expressions are
    # the same either way; that is what makes the two routes bit-identical.
    drift_full <- function(t, s) {
      u <- s[[1]]; hW <- s[[2]]; TC <- s[[3]]; TE <- s[[4]]; tau <- s[[5]]; I <- s[[6]]
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
    sd_full <- function(t, s) { n <- length(t)
      c(rep_len(p$s_u, n), rep_len(p$s_h, n), rep_len(p$s_C, n),
        rep_len(sigma_E, n), stau(s[[3]], t), sI(s[[6]])) }
    asm <- function(x, y) { s <- numeric(6); s[obs_i] <- x; s[hid] <- y; s }
    f_split <- function(t, x, y) drift_full(t, asm(x, y))[obs_i]
    g_split <- function(t, x, y) drift_full(t, asm(x, y))[hid]
    Sx <- function(t, x) { s <- numeric(6); s[obs_i] <- x
      diag(sd_full(t, s)[obs_i], length(obs_i)) }
    Sy_h <- function(t, x) { s <- numeric(6); s[obs_i] <- x
      diag(sd_full(t, s)[hid], length(hid), length(hid)) }
  }, envir = ce)
  lockEnvironment(ce, bindings = TRUE)
  m <- aci_model_from_affine(ce$f_split, ce$g_split, ce$Sx, ce$Sy_h,
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
    # The target is whichever observed channel the reference script leaves
    # uncommented in its masked observation precision, and the conditioning
    # set is everything else observed.  It is NOT the set of observed
    # channels the hidden variable happens to reach: the u and h_W scripts
    # mask down to T_C alone (u_unobs.m:1205, h_W_unobs.m:1202), leaving
    # T_E and I conditioned upon, while the tau and joint scripts leave the
    # whole block uncommented and are unconditional (tau_unobs.m:1205-1211,
    # u_h_W_tau_unobs.m:1218-1226).  Recording the reachable set instead
    # named a complement that is structurally inert - masking {u, tau} on
    # the h_W partition is identical() to not masking at all, because those
    # rows of Lx are zero - so a caller who built the mask from this field
    # got the unconditional run under a conditional name.
    target_names <- switch(source_key,
      u = "TC", hW = "TC", joint = vars[obs_i],
      tau = if (identical(observations, "reduced")) c("TC", "TE", "I")
            else vars[obs_i])
    conditioning_names <- setdiff(vars[obs_i], target_names)
    m$meta$target_obs_idx <- match(target_names, vars[obs_i])
    m$meta$conditioning_obs_idx <- match(conditioning_names, vars[obs_i])
    m$meta$causal_link <- sprintf("(%s) -> (%s)%s",
      paste(vars[hid], collapse = ","), paste(target_names, collapse = ","),
      if (length(conditioning_names)) sprintf(" | (%s)", paste(conditioning_names, collapse = ",")) else "")
    m$meta$estimand_provenance <- switch(source_key,
      u  = "ENSO_model_cond_ACI_u_unobs.m:1205 (S_xoS_x_inv(1,1,:) = 1/sigma_C^2)",
      hW = "ENSO_model_cond_ACI_h_W_unobs.m:1202 (S_xoS_x_inv(1,1,:) = 1/sigma_C^2)",
      tau = paste("ENSO_model_cond_ACI_tau_unobs.m:1205-1211 (every masking",
                  "line commented; the run is unconditional)"),
      joint = paste("ENSO_model_cond_ACI_u_h_W_tau_unobs.m:1218-1226 (every",
                    "masking line commented; the run is unconditional)"))
    m$meta$conditioning_note <- paste(
      "target_obs_idx and conditioning_obs_idx record the reference script's",
      "estimand: which observed channels carry the causal question, and which",
      "are conditioned upon, whether by masking (u, hW) or by prescription",
      "(tau, observations = 'reduced'). They are a record, not a default.",
      "Assimilation conditions only on an aci_conditional() specification",
      "the caller supplies or the constructor declares in",
      "meta$estimand_nontarget.")
  } else {
    m$meta$source_partition <- "package_only_partition"
  }
  m$meta$coefficient_provenance <-
    "ACI_code-main/ENSO_model_cond_ACI_*_unobs.m fixed-state coefficients"
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
  ic_all <- c(u = 6.9136e-04, hW = -0.0028, TC = 0.0039, TE = 0.0051,
              tau = -0.0256, I = 1.5841)
  m$meta$ic_default <- list(x0 = ic_all[vars[obs_i]], y0 = ic_all[vars[hid]])
  m$meta$observations <- observations
  if (identical(observations, "reduced")) {
    ## Observed order for this partition is (u, hW, TC, TE, I); u and hW are
    ## observed columns 1 and 2.  Positional indices, because the observation
    ## contract is positional everywhere else in the package and a name-based
    ## default would refuse an otherwise valid unnamed observation matrix.
    m$meta$estimand_nontarget <- aci_conditional(given = c(1L, 2L),
                                                 method = "reduce")
    m$meta$estimand_note <- paste(
      "Assimilation observes (TC, TE, I) and prescribes u and hW from the",
      "supplied observation columns, the ACI_code tau-script estimand. Use",
      "observations = 'full' for the five-channel estimand; the two differ by",
      "up to 0.247 in the filter mean and 0.776 in ACI on a 4001-point path.")
  }
  # Replaces the descriptor aci_model_from_affine() attached, and carries its
  # two drifts forward: a model whose batch specification does not check out
  # realises through the affine grid route, and one whose coefficient functions
  # have been replaced authenticates as nothing and realises generically.
  m <- .attach_cgns_realizer(m, "enso6_aci_code_v1", list(
    f_full = ce$f_split, g_full = ce$g_split,
    enso6 = list(drift = ce$drift_full, sd = ce$sd_full,
                 obs_i = as.integer(obs_i), hid = as.integer(hid),
                 n_state = length(vars))))
  m
}


#' Validate and store a prescribed forcing grid (internal)
#'
#' `ENSO_model_cond_ACI_T_C_unobs.m:999-1021` assimilates only `(T_E, I)` and
#' drives `u`, `h_W` and `tau` from their observed series. Those series are a
#' constructor input here, stored on their own grid and looked up by index.
#'
#' @param prescribed The constructor's `prescribed` argument.
#' @returns A list with `t0`, `dt`, `n` and the three numeric series.
#' @noRd
.enso6_prescribed_grid <- function(prescribed) {
  need <- c("t", "u", "hW", "tau")
  if (is.null(prescribed))
    aci_abort("aci_error_model_contract", paste(
      "hidden = 'TC' needs prescribed = list(t =, u =, hW =, tau =): unlike",
      "the other partitions this model is not self-contained, because",
      "ENSO_model_cond_ACI_T_C_unobs.m:999-1021 assimilates only (TE, I) and",
      "drives u, hW and tau from their observed series."))
  if (!is.list(prescribed) || is.null(names(prescribed)) ||
      !all(need %in% names(prescribed)))
    aci_abort("aci_error_model_contract",
              "prescribed must be a named list or data frame with t, u, hW and tau.")
  z <- lapply(prescribed[need], function(v) as.numeric(unlist(v, use.names = FALSE)))
  n <- lengths(z)
  if (any(n != n[[1L]]) || n[[1L]] < 2L ||
      !all(vapply(z, function(v) all(is.finite(v)), logical(1))))
    aci_abort("aci_error_model_contract", paste(
      "prescribed t, u, hW and tau must be finite numeric series of one",
      "common length of at least two."))
  step <- diff(z$t)
  if (any(step <= 0) || max(abs(step - step[1L])) > 1e-9 * step[1L])
    aci_abort("aci_error_model_contract", paste(
      "prescribed t must be a strictly increasing uniform grid: the forcings",
      "are looked up by index and never interpolated."))
  list(t0 = z$t[1L], dt = step[1L], n = as.integer(n[[1L]]),
       u = z$u, hW = z$hW, tau = z$tau)
}


#' The T_C-hidden zeroth-order ENSO inference model (internal)
#'
#' `ENSO_model_cond_ACI_T_C_unobs.m` does not build a six-state CGNS with `T_C`
#' hidden; it builds a two-channel reduced inference model whose coefficients
#' are the six-state drift rows for `T_E`, `I` and `T_C`, split about
#' `T_C = 0` with the cubic damping frozen there. Every expression below is the
#' corresponding line of `aci_enso_model()`'s own `drift_full()` and `sd_full()`
#' with `c1f(TC, t)` replaced by `c1f(0, t)`, so the two constructions share
#' one algebra rather than two transcriptions of it.
#'
#' @param p Completed parameter list.
#' @param sigma_E,lambda Constructor scalars.
#' @param variant Parameter and coefficient convention.
#' @param vars The six state names, in order.
#' @param prescribed The constructor's `prescribed` argument.
#' @param matlab_defect_compat `TRUE` to omit `gamma_C * h_W` from `f_y`.
#' @returns A `cgns_model` with `k = 2`, `l = 1`.
#' @noRd
.enso6_tc_zeroth_order <- function(p, sigma_E, lambda, variant, vars,
                                   prescribed, matlab_defect_compat) {
  # Named for the constructor arguments they carry, so the closures below refer
  # only to names this function binds.
  prescribed <- .enso6_prescribed_grid(prescribed)
  matlab_defect_compat <- isTRUE(matlab_defect_compat)
  # One locked environment, as the other partitions have: the coefficient
  # closures cannot outlive mutation of the captured parameters or forcings.
  ce <- list2env(list(p = p, sigma_E = sigma_E, lambda = lambda,
                      prescribed = prescribed,
                      matlab_defect_compat = matlab_defect_compat),
                 parent = baseenv())
  local({
    c1f <- function(TC, t)
      (25 * (TC + 0.75 / 7.5)^2 + 0.9) *
      (1 + 0.3 * sin(2 * pi * t / 6 - pi / 6)) * p$factor
    c2f <- function(t)
      1.4 * p$factor * (1 + 0.3 * sin(2 * pi * t / 6 + 2 * pi / 6) +
                          0.25 * sin(2 * (2 * pi * t / 6) + 2 * pi / 6))
    beta0 <- function(I) (1 + (1 - I / 5)) * 0.15 * sqrt(p$factor)
    sI <- function(I) { Ic <- pmin(pmax(I, 0), 4)
      sqrt(pmax(lambda * Ic * (4 - Ic), 0) + 1e-3 * lambda) }
    ## Prescribed forcings by nearest index, clamped at the ends, exactly as
    ## aci_conditional_reduce()'s prescribed-forcing lookup is. Nothing is
    ## interpolated. The grid is compared against the observation grid once per
    ## run, in .compile_cgns_run(), so a mismatched grid is refused there
    ## rather than silently resolved to a neighbouring step here.
    at <- function(t) {
      j <- as.integer(round((t - prescribed$t0) / prescribed$dt)) + 1L
      min(max(j, 1L), prescribed$n)
    }
    ## zeta_C = zeta_E = r_C = gamma b_0 mu / 2 in the ACI_code parameters.
    zeta <- p$gamma * p$b0 * p$mu / 2
    Lx <- function(t, x) matrix(c(-zeta, 0), 2L, 1L)
    fx <- function(t, x) { j <- at(t)
      c(p$gamma * prescribed$hW[j] +
          (1.5 * p$gamma * p$b0 * p$mu - c2f(t)) * x[1L] +
          beta0(x[2L]) * prescribed$tau[j],
        -lambda * (x[2L] - 2)) }
    Ly <- function(t, x) matrix(zeta - c1f(0, t), 1L, 1L)
    fy <- function(t, x) { j <- at(t)
      zeta * x[1L] + (x[2L] / 5) * p$factor * prescribed$u[j] + p$Cu +
        (0.8 * beta0(x[2L])) * prescribed$tau[j] +
        (if (matlab_defect_compat) 0 else p$gamma * prescribed$hW[j]) }
    Sx1 <- function(t, x) diag(c(sigma_E, sI(x[2L])), 2L, 2L)
    Sy2 <- function(t, x) matrix(p$s_C, 1L, 1L)
  }, envir = ce)
  lockEnvironment(ce, bindings = TRUE)
  m <- aci_model(Lx = ce$Lx, fx = ce$fx, Ly = ce$Ly, fy = ce$fy,
                 Sx1 = ce$Sx1, Sy2 = ce$Sy2, k = 2L, l = 1L,
                 name = sprintf("ENSO6[%s] (TC hidden, zeroth-order c1)",
                                 variant))
  m$meta$vars <- list(all = vars, hidden = "TC", observed = c("TE", "I"),
                      prescribed = c("u", "hW", "tau"))
  m$meta$params <- c(p, sigma_E = sigma_E, lambda = lambda)
  m$meta$variant <- variant
  m$meta$source_partition <- "TC"
  m$meta$approximation <- "zeroth_order_c1"
  m$meta$approximation_note <- paste(
    "L_y is the time-only series r_C - c1(t, 0), a zeroth-order Taylor",
    "expansion of the cubic damping about the climatology TC = 0",
    "(ENSO_model_cond_ACI_T_C_unobs.m:1052, :1150). The filter and smoother",
    "moments are those of this approximating model, not of the six-state",
    "system, whose simulator keeps the full nonlinear c1(t, TC) (:1105-1106,",
    ":1124).")
  m$meta$target_obs_idx <- 1L
  m$meta$conditioning_obs_idx <- 2L
  m$meta$causal_link <- "(TC) -> (TE) | (u,hW,tau,I)"
  m$meta$estimand_provenance <-
    "ENSO_model_cond_ACI_T_C_unobs.m:1209 (S_xoS_x_inv(1,1,:) = 1/sigma_E^2)"
  m$meta$conditioning_note <- paste(
    "target_obs_idx and conditioning_obs_idx record the reference script's",
    "estimand: which observed channels carry the causal question, and which",
    "are conditioned upon, whether by masking (I) or by prescription",
    "(u, hW, tau). They are a record, not a default. Assimilation conditions",
    "only on an aci_conditional() specification the caller supplies or the",
    "constructor declares in meta$estimand_nontarget.")
  ## Positional, as the tau declaration is: TE is observed column 1, and a
  ## name-based default would refuse an otherwise valid unnamed observation
  ## matrix.  The mask is exactly inert on this partition - the I row of Lx is
  ## zero and gyx is zero, so the I precision multiplies zero in the gain - so
  ## declaring it costs nothing and makes the shipped estimand the script's.
  m$meta$estimand_nontarget <- aci_conditional(target = 1L, method = "mask")
  m$meta$estimand_note <- paste(
    "Assimilation observes (TE, I), targets TE and prescribes u, hW and tau,",
    "the ACI_code T_C-script estimand. The TE-target mask is a structural",
    "no-op here: masked and unmasked runs agree bitwise on every moment.")
  m$meta$observations <- "reduced"
  m$meta$prescribed_grid <- list(t0 = prescribed$t0, dt = prescribed$dt,
                                 n = prescribed$n,
                                 channels = c("u", "hW", "tau"))
  m$meta$matlab_defect_compat <- list(
    active = matlab_defect_compat,
    defects = list(fy_gamma_hW_omitted = list(
      source = "ACI_code-main/ENSO_model_cond_ACI_T_C_unobs.m:1053,:1151",
      detail = paste("f_y omits gamma_C * h_W, which the simulator drift at",
                     ":1124 includes and which the sibling scripts carry in",
                     "the corresponding T_C coefficient rows",
                     "(tau_unobs.m:1138; u_h_W_tau_unobs.m:1016,:1131).",
                     "h_W is prescribed and observed in this script."),
      ruling = paste("Intended behaviour by default; the divergence is",
                     "measured and recorded."))))
  m$meta$coefficient_provenance <- paste(
    "ACI_code-main/ENSO_model_cond_ACI_T_C_unobs.m fixed-state coefficients,",
    "split about TC = 0")
  m$meta$coefficient_phase <- "state_time"
  m$meta$coefficient_phase_note <- paste(
    "Seasonal coefficients are evaluated at state time throughout, as every",
    "other acir constructor does. The MATLAB stored arrays hold element 1 at",
    "state time and elements 2..N+1 one step ahead (:1046/:1052 against",
    ":1140/:1150). On a 4001-point path the two conventions differ by up to",
    "9.3e-04 in the filter mean and 1.4e-02 in ACI; no ACI_code result",
    "depends on reproducing the one-step-ahead phase, so acir ships one",
    "convention and no switch.")
  m$meta$numerical_regularization <- list(
    I_variance_floor = 1e-3 * lambda,
    reason = paste("The exact sigma_I vanishes at I=0 and I=4, whereas the",
                   "closed-form assimilation implementation requires a",
                   "non-singular observed-noise Gramian; without the floor",
                   "this two-channel model does not even validate, because",
                   "gxx is degenerate at the I=0 probe point. ACI_code",
                   "instead pseudo-inverts that Gramian and gives the I",
                   "channel zero precision where sigma_I vanishes",
                   "(:1204, :1212-1213). Inert for assimilation here: the I",
                   "row of Lx and the noise cross-Gramian are both zero, and",
                   "the TE-target mask zeroes the I precision anyway."))
  m$meta$simulation_convention <- paste(
    "simulate() is refused on this model. Its hidden drift is the linearised",
    "series r_C - c1(t, 0), so integrating it would produce a path that is",
    "not an ENSO path.")
  m$meta$matlab_simulator_parity <- FALSE
  m$meta$simulate_supported <- FALSE
  m$meta$simulator_model <-
    "aci_enso_model(hidden = c(\"u\", \"hW\", \"tau\"), variant = \"aci_code\")"
  ic_all <- c(u = 6.9136e-04, hW = -0.0028, TC = 0.0039, TE = 0.0051,
              tau = -0.0256, I = 1.5841)
  m$meta$ic_default <- list(x0 = ic_all[c("TE", "I")], y0 = ic_all["TC"])
  # No batch realiser: this partition is a distinct two-channel construction
  # with its own coefficient closures, and it compiles through the generic
  # one-pass route.  Attaching nothing is what keeps the other partitions'
  # descriptor authentication meaningful.
  m
}


#' Noisy predator-prey benchmark model
#'
#' andreou2026aci supplementary model (SI.4.2; ACI_code noisy_predator_prey)
#' Construct either causal partition of the stochastic Lotka-Volterra example
#' used in andreou2026aci. The two partitions should be compared separately
#' because the supplied MATLAB file contains sequential direction-specific
#' blocks.
#'
#' @param hidden Either `"prey"` or `"predator"`, naming the hidden component.
#' @param params Optional named list overriding `alpha`, `beta`, `gamma`,
#'   `delta`, `s_x` and `s_y`.
#' @returns An object of class `cgns_model`.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications* **17**, 1854. \doi{10.1038/s41467-026-68568-0}
#' @seealso [aci_dyad_model()], [aci()]
#'
#' @examples
#' aci_predprey_model(hidden = "prey")
#'
#' @export
aci_predprey_model <- function(hidden = c("prey", "predator"),
                               params = list()) {
  hidden <- match.arg(hidden)
  if (!is.list(params))
    aci_abort("aci_error_model_contract", "params must be a list.")
  p <- utils::modifyList(list(alpha = 0.4, beta = 0.1, gamma = 1.1,
                              delta = 0.4, s_x = 0.3, s_y = 0.3), params)
  if (hidden == "prey") {          # observe predator x; hidden prey y
    m <- aci_model(
      Lx = function(t, x) matrix(p$beta * x[1], 1, 1),
      fx = function(t, x) -p$alpha * x[1],
      Ly = function(t, x) matrix(p$gamma - p$delta * x[1], 1, 1),
      fy = function(t, x) 0,
      Sx1 = function(t, x) matrix(p$s_x, 1, 1),
      Sy2 = function(t, x) matrix(p$s_y, 1, 1),
      k = 1, l = 1, name = "predator_prey[prey hidden]")
    m$meta$ic_default <- list(x0 = 4, y0 = 4)
  } else {                         # observe prey y; hidden predator x
    m <- aci_model(
      Lx = function(t, x) matrix(-p$delta * x[1], 1, 1),
      fx = function(t, x) p$gamma * x[1],
      Ly = function(t, x) matrix(p$beta * x[1] - p$alpha, 1, 1),
      fy = function(t, x) 0,
      Sx1 = function(t, x) matrix(p$s_y, 1, 1),
      Sy2 = function(t, x) matrix(p$s_x, 1, 1),
      k = 1, l = 1, name = "predator_prey[predator hidden]")
    m$meta$ic_default <- list(x0 = 4, y0 = 4)
  }
  if (any(!is.finite(unlist(p, use.names = FALSE))) ||
      any(c(p$s_x, p$s_y) <= 0))
    aci_abort("aci_error_model_contract",
              "Predator-prey parameters must be finite and noise amplitudes positive.")
  m$meta$params <- p
  m$meta$vars <- if (hidden == "prey")
    list(observed = "predator", hidden = "prey") else
    list(observed = "prey", hidden = "predator")
  m$meta$provenance <- "andreou2026aci SI.4.2; ACI_code noisy_predator_prey_model.m"
  m
}
