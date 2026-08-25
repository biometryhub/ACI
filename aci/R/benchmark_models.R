################################################################################
## benchmark_models.R - benchmark-model constructors
## ########################################################################## ##
##
## Contents:
##   * benchmark models from 'MATLAB Codebase for "Assimilative Causal Inference"':
##       - model_dyad, model_tipping_triad, model_pathways, model_l84, model_l96,
##         model_enso6, model_topographic, model_predator_prey
##
##   * benchmark models from 'MATLAB Codebase for "Bridging Prediction and
##     Attribution: Identifying Forward and Backward Causal Influence Ranges
##     Using Assimilative Causal Inference":
##       - .fbcir_multiscale_parameters, model_multiscale_fbcir,
##         .fbcir_layered_topographic_parameters, model_topographic_layered_fbcir
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


#' Nonlinear dyad (andreou2026aci eq. 1-2; jiang2026enkbs s3.2.2;
#' moser2026extremes eq. 4.1-4.2)
#'
#' @param variant Paper-specific parameter preset.
#' @param observe Which dyad component is treated as observed.
#' @param params Optional complete parameter list overriding the preset.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications* **17**, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#'
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' model_dyad()
#' model_dyad(variant = "p4", observe = "x")
#'
#' @export
model_dyad <- function(variant = c("p1", "p3", "p4"),
                       observe = c("x", "y"), params = NULL) {
  variant <- match.arg(variant); observe <- match.arg(observe)
  defaults <- switch(variant,
    p1 = list(d_x = 0.5, gamma = 2,   f_x = 0.5, s_x = 0.5, d_y = 0.5, f_y = 1,   s_y = 1),
    p3 = list(d_x = 0.5, gamma = 2,   f_x = 1,   s_x = 0.5, d_y = 0.5, f_y = 0.8, s_y = 1),
    p4 = list(d_x = 0.8, gamma = 1.2, f_x = 1,   s_x = 0.5, d_y = 0.8, f_y = 1,   s_y = 2))
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
    m <- cgns_model(
      Lx = coefficient_functions$Lx,
      fx = coefficient_functions$fx,
      Ly = coefficient_functions$Ly,
      fy = coefficient_functions$fy,
      Sx1 = coefficient_functions$Sx1,
      Sx2 = coefficient_functions$Sx2,
      Sy1 = coefficient_functions$Sy1,
      Sy2 = coefficient_functions$Sy2,
      k = 1, l = 1, name = sprintf("dyad[%s] y->x", variant))
  } else {
    m <- stochastic_model(
      f = function(t, x, y) -p$d_y * x - p$gamma * y^2 + p$f_y,
      g = function(t, x, y) (-p$d_x + p$gamma * x) * y + p$f_x,
      Sx = function(t, x) matrix(p$s_y, 1, 1),
      Sy = function(t, x, y) matrix(p$s_x, 1, 1),
      k = 1, l = 1, vectorized_members = TRUE,
      name = sprintf("dyad[%s] x->y (non-injective u^2)", variant))
  }
  m$meta$energy_conserving <- TRUE
  m$meta$params <- p
  m$meta$vars <- if (observe == "x")
    list(observed = "x", hidden = "y") else
    list(observed = "y", hidden = "x")
  m$meta$provenance <- switch(variant,
    p1 = paste("andreou2026aci Sections 3.1 and SI.4.1;",
               "ACI_code-main/dyad_interaction_model.m"),
    p3 = "jiang2026enkbs Section 3.2 bidirectional nonlinear-dyad EnKBS causal-inference benchmark",
    p4 = "moser2026extremes equations (4.1)-(4.2)")
  m$meta$source_status <- if (observe == "y" && variant == "p3") {
    "paper + MATLAB checked (published EnKBS dyad experiment)"
  } else if (observe == "y") {
    paste("package extension (source equations; reverse x -> y partition",
          "is not a supplied paper/MATLAB inference benchmark)")
  } else if (variant == "p1") {
    "paper + MATLAB checked"
  } else if (variant == "p3") {
    paste("paper checked; the published EnKBS MATLAB studies this direction",
          "with the ensemble engine")
  } else {
    "paper checked; no corresponding MATLAB supplied"
  }
  m$meta$partition_caveat <- if (observe == "y")
    paste("The reverse x -> y partition is not CGNS because the hidden x",
          "enters through x^2; this constructor is for ensemble methods.",
          if (variant == "p3")
            "jiang2026enkbs Section 3.2 studies this direction with EnKBS." else "") else NULL
  m$meta$anti_damping_threshold <- p$d_x / p$gamma
  # The MATLAB reference starts each component at its uncoupled forced
  # equilibrium (F_x / d_x, F_y / d_y).  In particular, the default
  # andreou2026aci configuration starts at (1, 2), not (1, 0).
  m$meta$ic_default <- if (observe == "x")
    list(x0 = p$f_x / p$d_x, y0 = p$f_y / p$d_y) else
    list(x0 = p$f_y / p$d_y, y0 = p$f_x / p$d_x)
  if (observe == "x")
    m <- .attach_cgns_realizer(
      m, "dyad_observed_x_v1", list(params = p)
    )
  m
}


#' Tipping triad (andreou2026cir Section 3.1, eq. 24a-c)
#'
#' @param eps Fast/slow time-scale ratio. The paper studies both regimes: the
#'   default `0.01` gives the bistable slow-transition regime, while the
#'   supplied FBCIR scripts run the intermittent regime at `0.1`.
#' @param params Optional complete parameter list.
#' @param partition Joint or conditional causal partition.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#' @examples
#' model_tipping_triad()
#'
#' @export
model_tipping_triad <- function(eps = 0.01, params = NULL,
                                partition = c("joint", "y_given_gamma",
                                              "gamma_given_y")) {
  partition <- match.arg(partition)
  defaults <- list(d_x = 1/3, alpha = 4, sigma_x = 0.2, d_y = 0.2,
                   beta = -0.8, sigma_y = 0.3, d_gamma = 0.5,
                   gamma_bar = 1, sigma_gamma = 2)
  p <- .complete_scalar_params(params, defaults, "tipping-triad",
    c("d_x", "alpha", "sigma_x", "d_y", "sigma_y", "d_gamma",
      "sigma_gamma"))
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    aci_abort("aci_error_model_contract", "eps must be one positive finite number.")

  if (partition == "joint") {
    # Observed x; hidden state ordered as (y, gamma).
    m <- cgns_model(
      Lx = function(t, x) matrix(c(-p$alpha, 0), 1, 2),
      fx = function(t, x) x - p$d_x * x^3,
      Ly = function(t, x) matrix(c(-p$d_y / eps, 0, x, -p$d_gamma), 2, 2),
      fy = function(t, x) c(p$beta, p$d_gamma * p$gamma_bar),
      Sx1 = function(t, x) matrix(p$sigma_x, 1, 1),
      Sy2 = function(t, x) diag(c(p$sigma_y / sqrt(eps), p$sigma_gamma)),
      k = 1, l = 2, name = sprintf("tipping_triad[joint,eps=%g]", eps))
    m$meta$vars <- list(observed = "x", hidden = c("y", "gamma"))
    m$meta$conditioning <- NULL
    m$meta$target_obs_idx <- 1L
    m$meta$conditioning_obs_idx <- integer()
    m$meta$causal_link <- "(y,gamma) -> x"
    m$meta$ic_default <- list(x0 = c(x = 0), y0 = c(y = 0, gamma = 0))
  } else if (partition == "y_given_gamma") {
    # Conditional FBCIR question y(t) -> x | gamma.  The observed state is
    # ordered as (x, gamma), so gamma remains in the dynamics and can be masked
    # from the likelihood (or supplied as a prescribed forcing) by the caller.
    m <- cgns_model(
      Lx = function(t, x) matrix(c(-p$alpha, 0), 2, 1),
      fx = function(t, x) c(x[1] - p$d_x * x[1]^3,
                             -p$d_gamma * (x[2] - p$gamma_bar)),
      Ly = function(t, x) matrix(-p$d_y / eps, 1, 1),
      fy = function(t, x) p$beta + x[2] * x[1],
      Sx1 = function(t, x) diag(c(p$sigma_x, p$sigma_gamma)),
      Sy2 = function(t, x) matrix(p$sigma_y / sqrt(eps), 1, 1),
      k = 2, l = 1, name = sprintf("tipping_triad[y->x|gamma,eps=%g]", eps))
    m$meta$vars <- list(observed = c("x", "gamma"), hidden = "y")
    m$meta$conditioning <- list(target = "x", conditioned_on = "gamma")
    m$meta$target_obs_idx <- 1L
    m$meta$conditioning_obs_idx <- 2L
    m$meta$causal_link <- "y -> x | gamma"
    m$meta$ic_default <- list(x0 = c(x = 0, gamma = 0), y0 = c(y = 0))
  } else {
    # Conditional FBCIR question gamma(t) -> y | x.  The observed state is
    # ordered as (x, y); x is the resolved conditioning channel.
    m <- cgns_model(
      Lx = function(t, x) matrix(c(0, x[1]), 2, 1),
      fx = function(t, x) c(x[1] - p$d_x * x[1]^3 - p$alpha * x[2],
                             p$beta - (p$d_y / eps) * x[2]),
      Ly = function(t, x) matrix(-p$d_gamma, 1, 1),
      fy = function(t, x) p$d_gamma * p$gamma_bar,
      Sx1 = function(t, x) diag(c(p$sigma_x, p$sigma_y / sqrt(eps))),
      Sy2 = function(t, x) matrix(p$sigma_gamma, 1, 1),
      k = 2, l = 1, name = sprintf("tipping_triad[gamma->y|x,eps=%g]", eps))
    m$meta$vars <- list(observed = c("x", "y"), hidden = "gamma")
    m$meta$conditioning <- list(target = "y", conditioned_on = "x")
    m$meta$target_obs_idx <- 2L
    m$meta$conditioning_obs_idx <- 1L
    m$meta$causal_link <- "gamma -> y | x"
    m$meta$ic_default <- list(x0 = c(x = 0, y = 0), y0 = c(gamma = 0))
  }
  m$meta$params <- c(p, eps = eps)
  m$meta$partition <- partition
  m$meta$provenance <- "andreou2026cir Section 3.1, equations (24a)-(24c)"
  m$meta$source_files <- switch(partition,
    joint = character(),
    y_given_gamma = "FBCIR_code-main/climate_tipping_y_bifurcation_driven.m",
    gamma_given_y = "FBCIR_code-main/climate_tipping_gamma_noise_induced.m")
  m$meta$source_status <- if (partition == "joint")
    "paper checked; joint partition is not a supplied MATLAB causal query" else
    "paper + MATLAB checked"
  m$meta$partition_caveat <- if (partition == "joint") {
    paste("The joint posterior (y,gamma) | x is the paper's CGNS form, but",
          "the supplied scripts analyse the two conditional causal links.")
  } else {
    paste("The conditioning channel remains in the observed state so its",
          "dynamics are retained; conditional ACI must mask its likelihood",
          "channel or prescribe its realised path, as in the supplied script.")
  }
  m
}


#' Hidden damping + forcing pathways model (moser2026extremes eq. 4.3-4.4)
#'
#' @param params Optional complete parameter list.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' model_pathways()
#'
#' @export
model_pathways <- function(params = NULL) {
  defaults <- list(d_u = 0.8, d_g = 0.8, d_b = 1, c = 1.2,
                   F_u = 1, F_g = 1, s_u = 0.5, s_g = 2, s_b = 2.5)
  p <- .complete_scalar_params(params, defaults, "pathways",
    c("d_u", "d_g", "d_b", "c", "s_u", "s_g", "s_b"))
  m <- cgns_model(
    Lx = function(t, x) matrix(c(p$c * x, 1), 1, 2),
    fx = function(t, x) -p$d_u * x + p$F_u,
    Ly = function(t, x) diag(c(-p$d_g, -p$d_b)),
    fy = function(t, x) c(-p$c * x^2 + p$F_g, 0),
    Sx1 = function(t, x) matrix(p$s_u, 1, 1),
    Sy2 = function(t, x) diag(c(p$s_g, p$s_b)),
    k = 1, l = 2, name = "pathways (u | gamma, b)")
  m$meta$params <- p
  m$meta$vars <- list(observed = "u", hidden = c("gamma", "b"))
  m$meta$provenance <- "moser2026extremes equations (4.3)-(4.4)"
  m$meta$source_status <- "paper checked; no corresponding MATLAB supplied"
  m$meta$anti_damping_threshold <- p$d_u / p$c
  m$meta$onset_defaults <- list(kappa = 0.8, T_pre = 1.5)
  m$meta$ic_default <- list(x0 = 1, y0 = c(0, 0))
  m
}


#' Stochastic Lorenz-84; observed (y, z), hidden x -> CGNS
#'
#' @param params Optional complete parameter list overriding the preset.
#' @param variant jiang2026enkbs constant-forcing or andreou2026cir/FBCIR
#'   seasonal preset.
#' @param target In the FBCIR preset, observed target (`"y"` or `"z"`); the
#'   other observed component is recorded as the conditioning channel.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#' @examples
#' model_l84(variant = "fbcir")
#'
#' @export
model_l84 <- function(params = NULL, variant = c("p3", "fbcir"),
                      target = c("y", "z")) {
  variant <- match.arg(variant)
  target <- match.arg(target)
  p <- params %||% switch(variant,
    p3 = list(a = 0.25, b = 4, F0 = 8, F_amp = 0,
              F_period = 73, g = 1, sig = 0.1),
    fbcir = list(a = 0.25, b = 4, F0 = 8, F_amp = 3,
                 F_period = 73, g = 1, sig = 0.2))
  # Accept the historical custom-parameter spelling Ff while representing the
  # FBCIR seasonal forcing explicitly as F0 + F_amp cos(2 pi t/F_period).
  forcing <- function(t) {
    if (!is.null(p$Ff)) {
      if (is.function(p$Ff)) p$Ff(t) else p$Ff
    } else {
      p$F0 + p$F_amp * cos(2 * pi * t / p$F_period)
    }
  }
  m <- cgns_model(
    Lx = function(t, x) matrix(c(-p$b * x[2] + x[1], p$b * x[1] + x[2]), 2, 1),
    fx = function(t, x) c(p$g - x[1], -x[2]),
    Ly = function(t, x) matrix(-p$a, 1, 1),
    fy = function(t, x) p$a * forcing(t) - x[1]^2 - x[2]^2,
    Sx1 = function(t, x) diag(p$sig, 2),
    Sy2 = function(t, x) matrix(p$sig, 1, 1),
    k = 2, l = 1, name = sprintf("L84[%s] (y,z | x)", variant))
  m$meta$params <- p
  m$meta$variant <- variant
  is_seasonal <- is.null(p$Ff) && !is.null(p$F_amp) && p$F_amp != 0
  m$meta$forcing <- list(type = if (is_seasonal)
    "seasonal_cosine" else "constant_or_custom",
    F0 = p$F0 %||% NA_real_, F_amp = p$F_amp %||% NA_real_,
    F_period = p$F_period %||% NA_real_)
  m$meta$provenance <- if (variant == "fbcir")
    "FBCIR_code-main/lorenz84.m" else "jiang2026enkbs model discovery benchmark"
  m$meta$vars <- list(observed = c("y", "z"), hidden = "x")
  if (variant == "fbcir") {
    m$meta$target_obs_idx <- match(target, c("y", "z"))
    m$meta$conditioning_obs_idx <- setdiff(1:2, m$meta$target_obs_idx)
    m$meta$causal_link <- sprintf("x -> %s | %s", target,
                                  setdiff(c("y", "z"), target))
  }
  # Physical MATLAB ordering is (x, y, z) = (1, 0, 1); this constructor's
  # observed vector is (y, z) and hidden scalar is x.
  m$meta$ic_default <- list(x0 = c(y = 0, z = 1), y0 = c(x = 1))
  m
}


#' Stochastic Lorenz-96 (jiang2026enkbs eq. 16)
#'
#' @param n Number of cyclic Lorenz-96 variables.
#' @param F Constant forcing.
#' @param sigma Legacy common noise amplitude override.
#' @param observe Observed component indices; `NULL` uses the selected preset.
#' @param preset jiang2026enkbs or legacy observation/noise convention.
#' @param sigma_observed Observed component noise amplitude.
#' @param sigma_hidden Hidden component noise amplitude.
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#' @examples
#' model_l96(n = 8)
#'
#' @export
model_l96 <- function(n = 40, F = 8, sigma = NULL, observe = NULL,
                      preset = c("p3", "legacy"),
                      sigma_observed = NULL, sigma_hidden = NULL) {
  preset <- match.arg(preset)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n != as.integer(n) || n < 4L)
    aci_abort("aci_error_model_contract", "n must be one integer of at least 4.")
  n <- as.integer(n)
  if (!is.numeric(F) || length(F) != 1L || !is.finite(F))
    aci_abort("aci_error_model_contract", "F must be one finite number.")
  if (is.null(observe))
    observe <- if (preset == "p3") seq(2, n, by = 2) else seq(1, n, by = 2)
  if (!is.numeric(observe) || !length(observe) || any(!is.finite(observe)) ||
      any(observe != as.integer(observe)) || any(observe < 1L | observe > n))
    aci_abort("aci_error_model_contract",
              "observe must contain integer component indices in 1:n.")
  if (!is.null(sigma)) {
    if (is.null(sigma_observed)) sigma_observed <- sigma
    if (is.null(sigma_hidden)) sigma_hidden <- sigma
  }
  sigma_observed <- sigma_observed %||% if (preset == "p3") sqrt(0.1) else 1
  sigma_hidden <- sigma_hidden %||% if (preset == "p3") sqrt(5) else 1
  if (!is.numeric(sigma_observed) || !is.numeric(sigma_hidden) ||
      length(sigma_observed) != 1L || length(sigma_hidden) != 1L ||
      any(!is.finite(c(sigma_observed, sigma_hidden))) ||
      sigma_observed < 0 || sigma_hidden < 0)
    aci_abort("aci_error_model_contract",
              "sigma_observed and sigma_hidden must be non-negative finite scalars.")
  obs_i <- sort(unique(as.integer(observe))); hid_i <- setdiff(seq_len(n), obs_i)
  if (!length(hid_i)) aci_abort("aci_error_model_contract", "Need at least one hidden component.")
  ip1 <- c(2:n, 1); im1 <- c(n, 1:(n - 1)); im2 <- c(n - 1, n, 1:(n - 2))
  drift_full <- function(Z) {
    if (is.matrix(Z)) (Z[ip1, , drop = FALSE] - Z[im2, , drop = FALSE]) *
                        Z[im1, , drop = FALSE] - Z + F
    else (Z[ip1] - Z[im2]) * Z[im1] - Z + F
  }
  assemble <- function(x, y) {
    if (is.matrix(y)) { Z <- matrix(0, n, ncol(y)); Z[obs_i, ] <- x; Z[hid_i, ] <- y }
    else { Z <- numeric(n); Z[obs_i] <- x; Z[hid_i] <- y }
    Z
  }
  m <- stochastic_model(
    f = function(t, x, y) { D <- drift_full(assemble(x, y))
      if (is.matrix(D)) D[obs_i, , drop = FALSE] else D[obs_i] },
    g = function(t, x, y) { D <- drift_full(assemble(x, y))
      if (is.matrix(D)) D[hid_i, , drop = FALSE] else D[hid_i] },
    Sx = function(t, x) diag(sigma_observed, length(obs_i)),
    Sy = function(t, x, y) diag(sigma_hidden, length(hid_i)),
    k = length(obs_i), l = length(hid_i), vectorized_members = TRUE,
    name = sprintf("L96[%s,n=%d,F=%g]", preset, n, F))
  m$meta$coords <- list(hidden = hid_i, obs = obs_i, period = n)
  m$meta$params <- list(n = n, F = F, sigma_observed = sigma_observed,
                        sigma_hidden = sigma_hidden, preset = preset)
  variable_names <- paste0("x", seq_len(n))
  m$meta$vars <- list(all = variable_names, observed = variable_names[obs_i],
                      hidden = variable_names[hid_i])
  exact_p3_preset <- preset == "p3" && n == 40L && F == 8 &&
    identical(obs_i, seq(2L, n, by = 2L)) &&
    isTRUE(all.equal(sigma_observed, sqrt(0.1))) &&
    isTRUE(all.equal(sigma_hidden, sqrt(5)))
  m$meta$provenance <- if (preset == "p3")
    "jiang2026enkbs equation (16), Lorenz-96 benchmark" else
    "Package legacy preset using the Lorenz-96 drift"
  m$meta$source_status <- if (exact_p3_preset) {
    "paper checked; parameters and scheme match the published EnKBS MATLAB"
  } else if (preset == "p3") {
    "package extension (jiang2026enkbs equations with a non-benchmark configuration)"
  } else {
    "package extension; legacy preset is not in the four ACI papers or supplied MATLAB"
  }
  m$meta$preset_caveat <- if (preset == "legacy") {
    paste("The legacy preset observes odd indices and assigns unit noise to",
          "both partitions; it is retained for compatibility and is not jiang2026enkbs.")
  } else if (!exact_p3_preset) {
    paste("The exact jiang2026enkbs benchmark uses n=40, F=8, even observed indices,",
          "observed variance 0.1, and hidden variance 5.")
  } else NULL
  # Deterministic perturbations keep construction RNG-neutral while retaining
  # the small, non-uniform perturbation used to leave the unstable equilibrium.
  z0 <- F + 0.01 * sin(seq_len(n) * sqrt(2))
  m$meta$ic_default <- list(x0 = z0[obs_i], y0 = z0[hid_i])
  m
}


#' Six-variable stochastic conceptual ENSO model (andreou2026aci SI via
#' chen2022enso eqs. 1a-1f).
#'
#' @param hidden Character vector naming hidden ENSO variables.
#' @param sigma_E Eastern-Pacific temperature noise amplitude.
#' @param lambda Decay rate for the diversity index.
#' @param params Optional complete parameter list overriding the preset.
#' @param variant Parameter and coefficient convention.
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications* **17**, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' Chen, N., Fang, X. and Yu, J.-Y. (2022). A multiscale model for El Nino
#' complexity. *npj Climate and Atmospheric Science* **5**, 16. arXiv:2104.07174.
#' @examples
#' model_enso6()
#'
#' @export
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


#' Spectral barotropic topographic flow (moser2026extremes eqs. 4.5-4.13); V
#' observed.
#' @param noise_convention Complex-noise component convention.
#' @param params Optional complete parameter list.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' model_topographic()
#'
#' @export
model_topographic <- function(noise_convention = c("split", "full"), params = NULL) {
  noise_convention <- match.arg(noise_convention)
  defaults <- list(gam = 0.18, nu = 0.02, nuV = 0.01,
                   s_psi = 0.03, s_V = 0.015)
  p <- .complete_scalar_params(params, defaults, "topographic",
    c("gam", "nu", "nuV", "s_psi", "s_V"))
  p$beta <- 5 * sqrt(p$gam)
  Kp <- rbind(c(0,2), c(1,2), c(0,1), c(1,1), c(2,1),
              c(1,0), c(2,0), c(1,-1), c(2,-1), c(1,-2))
  nK <- nrow(Kp); fullK <- rbind(Kp, -Kp)
  key <- function(k) paste(k[1], k[2], sep = "_")
  hpos <- c("1_1" = 0, "2_1" = -0.75, "1_0" = 0.15, "2_0" = 1.55,
            "1_-1" = 0.90, "2_-1" = 1.65) * p$gam
  hval <- function(k) { v <- hpos[key(k)]; if (is.na(v)) v <- hpos[key(-k)]
    if (is.na(v)) 0 else unname(v) }
  idx <- stats::setNames(seq_len(nK), apply(Kp, 1, key))
  phi_of <- function(z, k) {
    i <- idx[key(k)]; if (!is.na(i)) return(z[i])
    i <- idx[key(-k)]; if (!is.na(i)) return(Conj(z[i]))
    0 + 0i
  }
  in_K <- function(q) !is.na(idx[key(q)]) || !is.na(idx[key(-q)])
  triads <- lapply(seq_len(nK), function(i) {
    k <- Kp[i, ]; out <- list()
    for (r in seq_len(2 * nK)) {
      mv <- fullK[r, ]; q <- k - mv
      if (in_K(q)) {
        kperp <- c(-k[2], k[1]); mperp <- c(-mv[2], mv[1])
        out[[length(out) + 1]] <- list(m = mv, q = q,
          C = -sum((kperp - mperp) * mv), m2 = sum(mv^2), hm = hval(mv))
      }
    }
    out
  })
  k2 <- rowSums(Kp^2)
  drift_c <- function(t, V, z) {
    dz <- complex(nK)
    for (i in seq_len(nK)) {
      k <- Kp[i, ]; acc <- 0 + 0i
      for (tr in triads[[i]])
        acc <- acc + (tr$C / k2[i]) * phi_of(z, tr$q) *
                     (-tr$m2 * phi_of(z, tr$m) + tr$hm)
      dz[i] <- acc + 1i * k[1] * (p$beta / k2[i] - V) * z[i] +
               1i * k[1] * hval(k) / k2[i] * V - p$nu * z[i]
    }
    dV <- -p$nuV * V
    for (r in seq_len(2 * nK)) {
      kv <- fullK[r, ]
      dV <- dV + kv[1] * Im(hval(kv) * Conj(phi_of(z, kv)))
    }
    list(dz = dz, dV = dV)
  }
  sig_c <- if (noise_convention == "split") p$s_psi / sqrt(2) else p$s_psi
  m <- stochastic_model(
    f = function(t, x, y) { z <- complex(real = y[1:nK], imaginary = y[nK + 1:nK])
      drift_c(t, x[1], z)$dV },
    g = function(t, x, y) { z <- complex(real = y[1:nK], imaginary = y[nK + 1:nK])
      d <- drift_c(t, x[1], z)$dz; c(Re(d), Im(d)) },
    Sx = function(t, x) matrix(p$s_V, 1, 1),
    Sy = function(t, x, y) diag(sig_c, 2 * nK),
    k = 1, l = 2 * nK, name = "topographic (V | modes)")
  m$meta$params <- p
  mode_labels <- apply(Kp, 1, function(k) sprintf("%d_%d", k[1], k[2]))
  m$meta$vars <- list(observed = "V",
    hidden = c(paste0("Re_phi_", mode_labels),
               paste0("Im_phi_", mode_labels)))
  m$meta$provenance <- "moser2026extremes equations (4.5)-(4.13) and Appendix B"
  m$meta$source_status <- paste(
    "paper checked; no corresponding MATLAB supplied;",
    "complex-noise component convention remains an open discrepancy")
  m$meta$preset_caveat <- paste(
    "noise_convention='split' assigns sigma_psi/sqrt(2) to each real",
    "component; 'full' assigns sigma_psi to each. moser2026extremes does not disambiguate",
    "these complex-Wiener conventions.")
  m$meta$modes <- Kp
  m$meta$mode_families <- list(upper = c("1_1", "2_1"), zonal = c("1_0", "2_0"),
                               lower = c("1_-1", "2_-1"), merid = c("0_1", "0_2"))
  m$meta$noise_convention <- noise_convention
  m$meta$ic_default <- list(x0 = 0.5, y0 = rep(0.01, 2 * nK))
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
#' @seealso [model_dyad()], [aci()]
#'
#' @examples
#' model_predator_prey(hidden = "prey")
#'
#' @export
model_predator_prey <- function(hidden = c("prey", "predator"),
                                params = list()) {
  hidden <- match.arg(hidden)
  if (!is.list(params))
    aci_abort("aci_error_model_contract", "params must be a list.")
  p <- utils::modifyList(list(alpha = 0.4, beta = 0.1, gamma = 1.1,
                              delta = 0.4, s_x = 0.3, s_y = 0.3), params)
  if (hidden == "prey") {          # observe predator x; hidden prey y
    m <- cgns_model(
      Lx = function(t, x) matrix(p$beta * x[1], 1, 1),
      fx = function(t, x) -p$alpha * x[1],
      Ly = function(t, x) matrix(p$gamma - p$delta * x[1], 1, 1),
      fy = function(t, x) 0,
      Sx1 = function(t, x) matrix(p$s_x, 1, 1),
      Sy2 = function(t, x) matrix(p$s_y, 1, 1),
      k = 1, l = 1, name = "predator_prey[prey hidden]")
    m$meta$ic_default <- list(x0 = 4, y0 = 4)
  } else {                         # observe prey y; hidden predator x
    m <- cgns_model(
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


################################################################################
## FBCIR benchmark models
################################################################################
##
## Note: The multiscale constructor is a direct coefficient transcription of
## Eq. (25) in the FBCIR paper and the three multiscale_*.m scripts. The layered
## topographic constructor is deliberately separate from model_topographic(): it
## comes from FBCIR_code-main/topographic.m and is an extra repository case study,
## not a model analysed in the paper.
##
################################################################################

#' Parameters of the four-state multiscale FBCIR model (internal)
#'
#' @param epsilon Positive scale-separation parameter.
#' @param params Optional named list of parameter overrides.
#' @returns The completed named parameter list.
#' @noRd
.fbcir_multiscale_parameters <- function(epsilon, params) {
  defaults <- list(
    a1 = 1, c1 = 1 / 3, M = 0.5, M1 = 0.5, c2 = 0.4, M2 = -1.5,
    gamma1 = 0.5, N_cross = 4, gamma2 = 1.2,
    I = matrix(c(0.6, 0, 0, 2), 2, 2, byrow = TRUE),
    L = matrix(c(1, 0, 0, 1.5), 2, 2, byrow = TRUE),
    fx1 = 0, fx2_mean = 4, fx2_amplitude = 0.5, fx2_period = 36,
    fy1 = 1, fy2 = -1,
    sigma_x1 = 0.15, sigma_x2 = 0.3, sigma_y1 = 1, sigma_y2 = 2
  )
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0)
    aci_abort("aci_error_model_contract",
              "epsilon must be one positive finite number.")
  if (is.null(params)) params <- list()
  if (!is.list(params))
    aci_abort("aci_error_model_contract", "params must be a list or NULL.")
  unknown <- setdiff(names(params), names(defaults))
  if (length(unknown))
    aci_abort("aci_error_model_contract",
              sprintf("Unknown multiscale parameter(s): %s.", paste(unknown, collapse = ", ")))
  p <- utils::modifyList(defaults, params)
  scalar_names <- setdiff(names(defaults), c("I", "L"))
  scalar_values <- unlist(p[scalar_names], use.names = FALSE)
  if (length(scalar_values) != length(scalar_names) ||
      !is.numeric(scalar_values) || any(!is.finite(scalar_values)))
    aci_abort("aci_error_model_contract",
              "All scalar multiscale parameters must be finite numbers.")
  if (!is.matrix(p$I) || !is.numeric(p$I) || !identical(dim(p$I), c(2L, 2L)) ||
      any(!is.finite(p$I)) ||
      !is.matrix(p$L) || !is.numeric(p$L) || !identical(dim(p$L), c(2L, 2L)) ||
      any(!is.finite(p$L)))
    aci_abort("aci_error_model_contract", "I and L must be finite 2 by 2 matrices.")
  if (p$fx2_period <= 0 || p$gamma1 <= 0 || p$gamma2 <= 0 ||
      any(c(p$sigma_x1, p$sigma_x2, p$sigma_y1, p$sigma_y2) <= 0))
    aci_abort("aci_error_model_contract",
              "The forcing period, damping rates, and noise amplitudes must be positive.")
  p$epsilon <- epsilon
  p
}

#' Four-dimensional multiscale atmospheric FBCIR benchmark
#'
#' Constructs the joint CGNS split `(y1,y2) -> (x1,x2)` or either conditional
#' marginal split `y1 -> x1 | (x2,y2)` and `y2 -> x2 | (x1,y1)`.  In the
#' marginal splits the conditioning variables remain in the observed state;
#' select the target channel recorded in `meta$target_obs_idx` when applying a
#' conditional likelihood mask.
#'
#' @param epsilon Positive fast/slow timescale parameter.  The paper uses 0.1.
#' @param partition One of `"joint"`, `"y1"`, or `"y2"`.
#' @param params Optional named parameter overrides.
#' @return A `cgns_model` with the original shared Wiener channels preserved.
#'
#' @examples
#' model_multiscale_fbcir(epsilon = 0.1)
#'
#' @export
model_multiscale_fbcir <- function(epsilon = 0.1,
                                   partition = c("joint", "y1", "y2"),
                                   params = NULL) {
  partition <- match.arg(partition)
  p <- .fbcir_multiscale_parameters(epsilon, params)
  forcing_x2 <- function(t)
    p$fx2_mean + p$fx2_amplitude * sin(2 * pi * t / p$fx2_period)
  cam <- function(x1, x2) matrix(c(
    p$sigma_y1 * (p$L[1, 1] - p$I[1, 1] * x1) / p$gamma1,
    p$sigma_y2 * (p$L[1, 2] - p$I[1, 2] * x1) / p$gamma2,
    p$sigma_y1 * (p$L[2, 1] - p$I[2, 1] * x2) / p$gamma1,
    p$sigma_y2 * (p$L[2, 2] - p$I[2, 2] * x2) / p$gamma2
  ), 2, 2, byrow = TRUE)
  slow_forcing <- function(t, x1, x2) c(
    p$a1 * x1 - p$c1 * x1^3 -
      x2 * (p$M + p$M1 * x1 + p$M2 * x2) + p$fx1,
    -p$c2 * x2 + x1 * (p$M + p$M1 * x1 + p$M2 * x2) + forcing_x2(t)
  )

  if (partition == "joint") {
    m <- cgns_model(
      Lx = function(t, x) matrix(c(
        p$I[1, 1] * x[1] + p$L[1, 1],
        p$I[1, 2] * x[1] + p$L[1, 2],
        p$I[2, 1] * x[2] + p$L[2, 1],
        p$I[2, 2] * x[2] + p$L[2, 2]
      ), 2, 2, byrow = TRUE),
      fx = function(t, x) slow_forcing(t, x[1], x[2]),
      Ly = function(t, x) matrix(c(
        -p$gamma1 / p$epsilon, p$N_cross,
        -p$N_cross, -p$gamma2 / p$epsilon
      ), 2, 2, byrow = TRUE),
      fy = function(t, x) c(
        -p$L[1, 1] * x[1] - p$L[2, 1] * x[2] -
          p$I[1, 1] * x[1]^2 - p$I[2, 1] * x[2]^2 + p$fy1,
        -p$L[1, 2] * x[1] - p$L[2, 2] * x[2] -
          p$I[1, 2] * x[1]^2 - p$I[2, 2] * x[2]^2 + p$fy2
      ),
      Sx1 = function(t, x) diag(c(p$sigma_x1, p$sigma_x2)),
      Sx2 = function(t, x) cam(x[1], x[2]),
      Sy1 = function(t, x) matrix(0, 2, 2),
      Sy2 = function(t, x) diag(c(p$sigma_y1, p$sigma_y2) / sqrt(p$epsilon)),
      k = 2, l = 2, name = "multiscale_fbcir[(y1,y2)->(x1,x2)]"
    )
    m$meta$vars <- list(observed = c("x1", "x2"), hidden = c("y1", "y2"))
    m$meta$target_obs_idx <- 1:2
    m$meta$conditioning_obs_idx <- integer()
    m$meta$causal_link <- "(y1,y2) -> (x1,x2)"
    m$meta$ic_default <- list(
      x0 = c(x1 = 0, x2 = forcing_x2(0) / p$c2),
      y0 = c(y1 = p$fy1 / (p$gamma1 / p$epsilon),
             y2 = p$fy2 / (p$gamma2 / p$epsilon))
    )
  } else if (partition == "y1") {
    # Observed ordering (x1, x2, y2); hidden y1.  W1 is ordered
    # (Wx1, Wx2, Wy2), while W2 is Wy1, exactly as in the MATLAB script.
    m <- cgns_model(
      Lx = function(t, x) matrix(c(
        p$I[1, 1] * x[1] + p$L[1, 1],
        p$I[2, 1] * x[2] + p$L[2, 1],
        -p$N_cross
      ), 3, 1),
      fx = function(t, x) {
        fslow <- slow_forcing(t, x[1], x[2])
        c(fslow[1] + (p$I[1, 2] * x[1] + p$L[1, 2]) * x[3],
          fslow[2] + (p$I[2, 2] * x[2] + p$L[2, 2]) * x[3],
          -p$gamma2 * x[3] / p$epsilon - p$L[1, 2] * x[1] -
            p$L[2, 2] * x[2] - p$I[1, 2] * x[1]^2 -
            p$I[2, 2] * x[2]^2 + p$fy2)
      },
      Ly = function(t, x) matrix(-p$gamma1 / p$epsilon, 1, 1),
      fy = function(t, x)
        -p$L[1, 1] * x[1] - p$L[2, 1] * x[2] + p$N_cross * x[3] -
        p$I[1, 1] * x[1]^2 - p$I[2, 1] * x[2]^2 + p$fy1,
      Sx1 = function(t, x) {
        C <- cam(x[1], x[2])
        matrix(c(p$sigma_x1, 0, C[1, 2],
                 0, p$sigma_x2, C[2, 2],
                 0, 0, p$sigma_y2 / sqrt(p$epsilon)),
               3, 3, byrow = TRUE)
      },
      Sx2 = function(t, x) {
        C <- cam(x[1], x[2])
        matrix(c(C[1, 1], C[2, 1], 0), 3, 1)
      },
      Sy1 = function(t, x) matrix(0, 1, 3),
      Sy2 = function(t, x) matrix(p$sigma_y1 / sqrt(p$epsilon), 1, 1),
      k = 3, l = 1, name = "multiscale_fbcir[y1->x1|(x2,y2)]"
    )
    m$meta$vars <- list(observed = c("x1", "x2", "y2"), hidden = "y1")
    m$meta$target_obs_idx <- 1L
    m$meta$conditioning_obs_idx <- c(2L, 3L)
    m$meta$causal_link <- "y1 -> x1 | (x2,y2)"
    m$meta$ic_default <- list(
      x0 = c(x1 = 0, x2 = forcing_x2(0) / p$c2,
             y2 = p$fy2 / (p$gamma2 / p$epsilon)),
      y0 = c(y1 = p$fy1 / (p$gamma1 / p$epsilon))
    )
  } else {
    # Observed ordering (x1, x2, y1); hidden y2.  W1 is ordered
    # (Wx1, Wx2, Wy1), while W2 is Wy2, exactly as in the MATLAB script.
    m <- cgns_model(
      Lx = function(t, x) matrix(c(
        p$I[1, 2] * x[1] + p$L[1, 2],
        p$I[2, 2] * x[2] + p$L[2, 2],
        p$N_cross
      ), 3, 1),
      fx = function(t, x) {
        fslow <- slow_forcing(t, x[1], x[2])
        c(fslow[1] + (p$I[1, 1] * x[1] + p$L[1, 1]) * x[3],
          fslow[2] + (p$I[2, 1] * x[2] + p$L[2, 1]) * x[3],
          -p$gamma1 * x[3] / p$epsilon - p$L[1, 1] * x[1] -
            p$L[2, 1] * x[2] - p$I[1, 1] * x[1]^2 -
            p$I[2, 1] * x[2]^2 + p$fy1)
      },
      Ly = function(t, x) matrix(-p$gamma2 / p$epsilon, 1, 1),
      fy = function(t, x)
        -p$L[1, 2] * x[1] - p$L[2, 2] * x[2] - p$N_cross * x[3] -
        p$I[1, 2] * x[1]^2 - p$I[2, 2] * x[2]^2 + p$fy2,
      Sx1 = function(t, x) {
        C <- cam(x[1], x[2])
        matrix(c(p$sigma_x1, 0, C[1, 1],
                 0, p$sigma_x2, C[2, 1],
                 0, 0, p$sigma_y1 / sqrt(p$epsilon)),
               3, 3, byrow = TRUE)
      },
      Sx2 = function(t, x) {
        C <- cam(x[1], x[2])
        matrix(c(C[1, 2], C[2, 2], 0), 3, 1)
      },
      Sy1 = function(t, x) matrix(0, 1, 3),
      Sy2 = function(t, x) matrix(p$sigma_y2 / sqrt(p$epsilon), 1, 1),
      k = 3, l = 1, name = "multiscale_fbcir[y2->x2|(x1,y1)]"
    )
    m$meta$vars <- list(observed = c("x1", "x2", "y1"), hidden = "y2")
    m$meta$target_obs_idx <- 2L
    m$meta$conditioning_obs_idx <- c(1L, 3L)
    m$meta$causal_link <- "y2 -> x2 | (x1,y1)"
    m$meta$ic_default <- list(
      x0 = c(x1 = 0, x2 = forcing_x2(0) / p$c2,
             y1 = p$fy1 / (p$gamma1 / p$epsilon)),
      y0 = c(y2 = p$fy2 / (p$gamma2 / p$epsilon))
    )
  }

  m$meta$params <- p
  m$meta$partition <- partition
  m$meta$provenance <- "FBCIR_paper_equation_25_and_repository_MATLAB"
  m$meta$source_files <- switch(partition,
                                joint = "FBCIR_code-main/multiscale_joint_y_1_y_2_cause.m",
                                y1 = "FBCIR_code-main/multiscale_marginal_y_1_cause.m",
                                y2 = "FBCIR_code-main/multiscale_marginal_y_2_cause.m")
  m$meta$shared_noise_preserved <- TRUE
  m$meta$matlab_pathwise_rng_parity <- FALSE
  m$meta$simulation_convention <- paste(
    "Euler-Maruyama coefficients match the source. R and MATLAB draw shared",
    "Wiener channels in different RNG orders, so equal seeds are not pathwise comparable.")
  m
}

#' Parameters of the seven-state layered topographic model (internal)
#'
#' @param params Optional named list of parameter overrides.
#' @returns The completed named parameter list.
#' @noRd
.fbcir_layered_topographic_parameters <- function(params) {
  defaults <- list(beta = 2, d = 1 / 40, p = 1,
                   topographic_strength = 2)
  if (is.null(params)) params <- list()
  if (!is.list(params))
    aci_abort("aci_error_model_contract", "params must be a list or NULL.")
  unknown <- setdiff(names(params), names(defaults))
  if (length(unknown))
    aci_abort("aci_error_model_contract",
              sprintf("Unknown layered-topographic parameter(s): %s.",
                      paste(unknown, collapse = ", ")))
  p <- utils::modifyList(defaults, params)
  vals <- unlist(p, use.names = FALSE)
  if (length(vals) != length(defaults) || !is.numeric(vals) || any(!is.finite(vals)) ||
      p$d <= 0 || p$topographic_strength <= 0)
    aci_abort("aci_error_model_contract",
              "Layered-topographic parameters must be finite, with d and topographic_strength positive.")
  p
}

#' Seven-dimensional layered topographic FBCIR repository example
#'
#' This is the additional model in `FBCIR_code-main/topographic.m`, not the
#' spectral topographic model from the ACI extreme-events paper.  The observed
#' state is `(v1,...,v6)` and the hidden state is the mean flow `u`.
#'
#' @param target_mode Streamfunction mode (1, 2, or 3) treated as the target
#'   pair by conditional ACI; this changes metadata, not the physical model.
#' @param params Optional named parameter overrides.
#' @return A `cgns_model`.
#'
#' @examples
#' model_topographic_layered_fbcir()
#'
#' @export
model_topographic_layered_fbcir <- function(target_mode = 1L, params = NULL) {
  if (!is.numeric(target_mode) || length(target_mode) != 1L ||
      !is.finite(target_mode) || target_mode != as.integer(target_mode) ||
      !(target_mode %in% 1:3))
    aci_abort("aci_error_model_contract", "target_mode must be one of 1, 2, or 3.")
  target_mode <- as.integer(target_mode)
  p <- .fbcir_layered_topographic_parameters(params)
  modes <- 3L
  H <- p$topographic_strength * sqrt(2) / 2^(seq_len(modes))
  omega <- H / sqrt(2)
  sigma_v <- c(rep(1 / (20 * sqrt(2)), 4),
               rep(1 / (20 * 3^p$p * sqrt(2)), 2))
  sigma_u <- 1 / (20 * sqrt(2))

  m <- cgns_model(
    Lx = function(t, x) matrix(c(
      x[2] - omega[1],
      -x[1],
      2 * x[4] - omega[2] / 2,
      -2 * x[3],
      3 * x[6] - omega[3] / (sqrt(2) * 3^(p$p + 1)),
      -3 * x[5]
    ), 6, 1),
    fx = function(t, x) c(
      -p$d * x[1] - p$beta * x[2],
      -p$d * x[2] + p$beta * x[1],
      -p$d * x[3] - p$beta * x[4] / 2,
      -p$d * x[4] + p$beta * x[3] / 2,
      -p$d * x[5] - p$beta * x[6] / 3,
      -p$d * x[6] + p$beta * x[5] / 3
    ),
    Ly = function(t, x) matrix(-p$d, 1, 1),
    fy = function(t, x)
      omega[1] * x[1] + 2 * omega[2] * x[3] +
      omega[3] * x[5] / (sqrt(2) * 3^(p$p - 1)),
    # topographic.m first defines sigma_v and then sets S_x=diag(sigma_v)/sqrt(2)
    Sx1 = function(t, x) diag(sigma_v / sqrt(2), 6),
    Sy2 = function(t, x) matrix(sigma_u, 1, 1),
    k = 6, l = 1,
    name = sprintf("topographic_layered_fbcir[u->psi%d|other_modes]", target_mode)
  )
  target <- (2L * target_mode - 1L):(2L * target_mode)
  m$meta$params <- c(p, list(H = H, omega = omega,
                             sigma_v = sigma_v, sigma_u = sigma_u))
  m$meta$vars <- list(observed = paste0("v", 1:6), hidden = "u")
  m$meta$target_mode <- target_mode
  m$meta$target_obs_idx <- target
  m$meta$conditioning_obs_idx <- setdiff(1:6, target)
  m$meta$causal_link <- sprintf("u -> psi%d | other streamfunction modes", target_mode)
  m$meta$ic_default <- list(x0 = stats::setNames(rep(0, 6), paste0("v", 1:6)),
                            y0 = c(u = 0))
  m$meta$provenance <- "repository_case_study_named_in_andreou2026cir_code_availability"
  m$meta$source_files <- "FBCIR_code-main/topographic.m"
  m$meta$paper_scope <- FALSE
  m$meta$matlab_pathwise_rng_parity <- FALSE
  m$meta$simulation_convention <- paste(
    "Euler-Maruyama coefficients match the source. R and MATLAB draw observed",
    "and hidden Wiener channels in different RNG orders, so equal seeds are not pathwise comparable.")
  m
}
