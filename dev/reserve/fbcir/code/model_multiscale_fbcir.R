## acir reserve file
## Origin: aci/R/benchmark_models.R:816-1050
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: FBCIR paper eq. (25) and the three multiscale_*.m scripts; includes the FBCIR section banner and .fbcir_multiscale_parameters.
## Verbatim copy from the aci 0.0.30 sources; not modified.

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
