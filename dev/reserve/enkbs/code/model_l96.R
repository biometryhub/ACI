## acir reserve file
## Origin: aci/R/benchmark_models.R:391-491
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: jiang2026enkbs eq. (16) Lorenz-96 benchmark.
## Verbatim copy from the aci 0.0.30 sources; not modified.

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
