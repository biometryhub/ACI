## acir reserve file
## Origin: aci/R/ensemble.R
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: The whole EnKBF/EnKBS ensemble engine; jiang2026enkbs, EnKBS-main.
## Verbatim copy from the aci 0.0.30 sources; not modified.

################################################################################
## ensemble.R - the derivative-free engine: EnKBF/EnKBS with noise reuse
## ########################################################################## ##
##
## Contents:
##   * ensemble Kalman-Bucy filter (perturbed simulated observations):
##       - new_ens_path, print.da_path_ensemble, as_gaussian, .eval_members, enkbf,
##         da_filter.stochastic_model
##
##   * ensemble Kalman-Bucy smoother (forward-noise reuse, hash-enforced):
##       - enkbs, da_smooth.stochastic_model
##
##   * andreou2026cir/jiang2026enkbs forward finite-lag posterior family:
##       - ensemble_lag_table
##
################################################################################


################################################################################
# EnKBF: ensemble Kalman-Bucy filter (perturbed simulated observations)
################################################################################

#' Construct an ensemble assimilation path (internal)
#'
#' @param t Numeric vector of times.
#' @param members Numeric array of ensemble members, hidden dimension by member
#'   by time.
#' @param kind 1-length character, `"filter"` or `"smoother"`.
#' @param meta Optional named list of metadata carried on the object.
#' @returns An object of class `da_path_ensemble`.
#' @noRd
new_ens_path <- function(t, members, kind, meta = list()) {
  structure(list(t = t, members = members, kind = kind, meta = meta),
            class = c("da_path_ensemble", "da_path"))
}


#' Print an ensemble assimilation path
#'
#' @param x A `da_path_ensemble` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.da_path_ensemble <- function(x, ...) {
  d <- dim(x$members)
  cat(sprintf("<da_path_ensemble> kind = %s, l = %d, N+1 = %d, m = %d\n",
              x$kind, d[1], d[2], d[3])); invisible(x)
}


#' Gaussian moments of an ensemble path
#'
#' Reduces an ensemble path to its per-time sample mean and covariance, giving
#' an object the closed-form diagnostics can consume.
#'
#' @param path A `da_path_ensemble` object.
#' @returns A `da_path_gaussian` object carrying the ensemble moments.
#'
#' @seealso [enkbf()], [cross_validate()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' fr <- enkbf(m, ob, m = 20, seed = 1)
#' as_gaussian(fr$path)
#'
#' @export
as_gaussian <- function(path) {
  if (!inherits(path, "da_path_ensemble"))
    aci_abort("aci_error_dims", "path must be a da_path_ensemble object.")
  d <- dim(path$members)
  if (length(d) != 3L || any(d < 1L))
    aci_abort("aci_error_dims", "Ensemble members must be an l x (N+1) x m array.")
  if (!is.numeric(path$t) || length(path$t) != d[2] ||
      any(!is.finite(path$t)))
    aci_abort("aci_error_dims",
              "The ensemble path time grid must contain one finite value per member-array time index.")
  l <- d[1]; N1 <- d[2]
  if (d[3] <= l)
    aci_abort("aci_error_ensemble_rank",
              "Moment-matched Gaussian paths require ensemble size m > l.")
  MU <- matrix(NA_real_, N1, l); CV <- array(NA_real_, c(l, l, N1))
  for (j in seq_len(N1)) {
    Y <- path$members[, j, , drop = FALSE][, 1, ]
    Y <- matrix(Y, nrow = l)
    MU[j, ] <- rowMeans(Y)
    Rj <- sym(tcrossprod(Y - MU[j, ]) / (d[3] - 1))
    if (is.null(tryCatch(chol(Rj), error = function(e) NULL)))
      aci_abort("aci_error_ensemble_rank", sprintf(
        "Ensemble covariance is rank deficient at index %d; increase m or use a separately justified projection.",
        j))
    CV[, , j] <- Rj
  }
  p <- new_da_path(path$t, MU, CV, path$kind, path$meta)
  p$meta$moment_matched_from_m <- d[3]
  p
}


#' Evaluate a drift or diffusion across ensemble members (internal)
#'
#' @param fun Function of `(t, x, y)` to evaluate.
#' @param t 1-length numeric time.
#' @param x Numeric vector of observations at time `t`.
#' @param Y Numeric matrix of hidden members, one column per member.
#' @param vectorized `TRUE` when `fun` accepts the whole member matrix.
#' @param out_rows Number of rows the result must carry.
#' @returns A numeric matrix with `out_rows` rows and one column per member.
#' @noRd
.eval_members <- function(fun, t, x, Y, vectorized, out_rows) {
  v <- if (vectorized) fun(t, x, Y) else
    vapply(seq_len(ncol(Y)), function(i) fun(t, x, Y[, i]),
           numeric(out_rows))
  if (!is.numeric(v) || length(v) != out_rows * ncol(Y) ||
      any(!is.finite(v)))
    aci_abort("aci_error_model_contract", sprintf(
      "A member drift must return a finite %d by %d result.",
      out_rows, ncol(Y)))
  matrix(v, nrow = out_rows)
}


#' Hash the forward Wiener increments (internal)
#'
#' Used to enforce that the backward pass reuses the increments of the forward
#' run it is paired with.
#'
#' @param noise The stored Wiener increments of a forward run.
#' @returns 1-length character hash.
#' @noRd
.noise_fingerprint <- function(noise) {
  tf <- tempfile(pattern = "aci-noise-", fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(list(B = noise$B, W = noise$W), tf, version = 3)
  unname(tools::md5sum(tf))
}


#' Check the ensemble engine's model contract (internal)
#'
#' Rejects models whose observation and signal noise share Wiener channels, and
#' models whose signal diffusion depends on the hidden member.
#'
#' @param model A `stochastic_model` object.
#' @param ob An `obs_traj` object.
#' @param Y Optional matrix of hidden members used to probe the diffusion.
#' @returns Invisibly `TRUE`; called for its error conditions.
#' @noRd
.validate_ensemble_contract <- function(model, ob, Y = NULL) {
  # The EnKBF/EnKBS derivation supplied with the package assumes independent
  # observation and signal Wiener processes. The closed-form CGNS engine is
  # the supported route for shared/correlated channels.
  if (inherits(model, "cgns_model")) {
    corr <- any(vapply(seq_along(ob$t), function(j) {
      .has_cross_noise(eval_coefs(model, ob$t[j], ob$x[j, ]))
    }, logical(1)))
    if (corr)
      aci_abort("aci_error_ensemble_noise_contract", paste(
        "The ensemble engine implements the independent-noise EnKBF/EnKBS",
        "from the supplied paper and cannot be used with shared/correlated",
        "CGNS noise. Use engine = 'cgns'."))
  }
  if (!is.null(Y)) {
    # The implemented paper recursion evaluates signal diffusion as Sigma(t,x),
    # not separately at each hidden member. Refuse hidden-dependent diffusion
    # rather than silently substituting the ensemble mean.
    y0 <- Y[, 1]
    ref <- as.matrix(model$Sy(ob$t[1], ob$x[1, ], y0))
    probes <- unique(c(1L, ncol(Y), max(1L, ceiling(ncol(Y) / 2))))
    for (i in probes) {
      now <- as.matrix(model$Sy(ob$t[1], ob$x[1, ], Y[, i]))
      if (!identical(dim(now), dim(ref)) ||
          max(abs(now - ref)) > 1e-12 * max(1, max(abs(ref))))
        aci_abort("aci_error_ensemble_noise_contract", paste(
          "The supplied EnKBF/EnKBS recursion requires Sy(t, x) to be",
          "independent of the hidden ensemble member."))
    }
  }
  invisible(TRUE)
}


#' Signal diffusion required to be member-independent (internal)
#'
#' @param model A `stochastic_model` object.
#' @param t 1-length numeric time.
#' @param x Numeric vector of observations at time `t`.
#' @param Y Numeric matrix of hidden members.
#' @returns The signal diffusion matrix at `(t, x)`.
#' @noRd
.ensemble_signal_diffusion <- function(model, t, x, Y) {
  ref <- as.matrix(model$Sy(t, x, Y[, 1]))
  for (i in seq_len(ncol(Y))) {
    now <- as.matrix(model$Sy(t, x, Y[, i]))
    if (!identical(dim(now), dim(ref)) ||
        max(abs(now - ref)) > 1e-12 * max(1, max(abs(ref))))
      aci_abort("aci_error_ensemble_noise_contract", paste(
        "The supplied EnKBF/EnKBS recursion requires Sy(t, x) to be",
        "independent of the hidden ensemble member."))
  }
  ref
}


#' Ensemble Kalman-Bucy filter
#'
#' Ensemble Kalman-Bucy filtering following the independent-noise jiang2026enkbs
#' recursion. Shared/correlated observation and signal noise is not accepted by
#' this engine; use the closed-form CGNS engine for that case.
#'
#' @param model A `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param m Ensemble size; a positive whole number.
#' @param ic_sampler Optional function of the ensemble size drawing the initial
#'   members, or a list with the initial `mean` and `cov` from which members
#'   are drawn; `NULL` uses the model's default.
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param inflation Positive multiplicative variance inflation factor.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param seed Optional non-negative whole number seeding the generator.
#' @param noise Optional stored Wiener increments to drive the run.
#' @param keep_members `TRUE` to retain the member paths.
#' @param ... Must be empty; unused arguments are an error.
#' @returns A `da_path_ensemble` object of kind `"filter"`, carrying the Wiener
#'   increments [enkbs()] requires.
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#' @seealso [enkbs()], [as_gaussian()], [localization_spec()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' enkbf(m, ob, m = 20, seed = 1)
#'
#' @export
enkbf <- function(model, obs, m, ic_sampler = NULL, localization = NULL,
                  inflation = 1, nontarget = NULL, seed = NULL, noise = NULL,
                  keep_members = TRUE, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to enkbf().")
  obs <- as_obs(obs)
  if (!inherits(model, "stochastic_model"))
    aci_abort("aci_error_model_contract", "enkbf() requires a stochastic_model.")
  if (obs$k != model$k)
    aci_abort("aci_error_dims", "Observation dimension does not match the stochastic model.")
  rs <- .resolve_nontarget(model, obs, nontarget)
  md <- rs$model; ob <- rs$obs
  N1 <- length(ob$t); N <- N1 - 1L; dt <- ob$dt
  k <- md$k; l <- md$l; sq <- sqrt(dt)
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) ||
      m != floor(m) || m < 2L)
    aci_abort("aci_error_ensemble_rank",
              "m must be an integer of at least 2.")
  m <- as.integer(m)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)))
    aci_abort("aci_error_dims", "seed must be NULL or one non-negative integer.")
  if (length(inflation) != 1L || !is.finite(inflation) || inflation < 1)
    aci_abort("aci_error_dims", "inflation must be a finite scalar >= 1.")
  if (!isTRUE(keep_members))
    aci_abort("aci_error_not_implemented",
              "keep_members = FALSE is unsupported because EnKBS requires the forward members.")
  if (!is.null(localization) && !is.list(localization))
    aci_abort("aci_error_dims", "localization must be NULL or a localization_spec-like list.")
  if (!is.null(seed)) set.seed(seed)
  Y <- if (is.function(ic_sampler)) ic_sampler(m) else {
    ic <- ic_sampler %||% list(mean = md$meta$ic_default$y0 %||% rep(0, l),
                               cov = diag(1, l))
    if (!is.list(ic) || length(ic$mean %||% numeric()) != l ||
        !identical(dim(as.matrix(ic$cov)), c(l, l)) ||
        any(!is.finite(c(ic$mean, ic$cov))))
      aci_abort("aci_error_dims", "ic_sampler list must contain finite mean (l) and cov (l x l).")
    ch <- safe_chol(as.matrix(ic$cov))
    as.numeric(ic$mean) + t(ch) %*% matrix(stats::rnorm(l * m), l, m)
  }
  if (length(Y) != l * m || any(!is.finite(Y)))
    aci_abort("aci_error_dims", "ic_sampler must return l x m finite initial members.")
  Y <- matrix(Y, l, m)
  .validate_ensemble_contract(md, ob, Y)
  dyn <- ncol(as.matrix(md$Sy(ob$t[1], ob$x[1, ], Y[, 1])))
  dxn <- ncol(as.matrix(md$Sx(ob$t[1], ob$x[1, ])))
  if (is.null(noise)) {
    B <- array(stats::rnorm(dyn * N * m), c(dyn, N, m))
    W <- array(stats::rnorm(dxn * N * m), c(dxn, N, m))
    noise <- structure(list(B = B, W = W, seed = seed), class = "noise_store")
  } else {
    if (!inherits(noise, "noise_store"))
      aci_abort("aci_error_noise_mismatch", "noise must be a noise_store from enkbf().")
    B <- noise$B; W <- noise$W
  }
  if (!identical(dim(B), c(dyn, N, m)) || !identical(dim(W), c(dxn, N, m)))
    aci_abort("aci_error_noise_mismatch",
              "Supplied ensemble noise has dimensions inconsistent with the model, grid, or m.")
  if (any(!is.finite(c(B, W))))
    aci_abort("aci_error_noise_mismatch", "Ensemble noise must be finite.")
  noise$fingerprint <- .noise_fingerprint(noise)
  C1 <- localization$C1
  if (!is.null(C1) && !identical(dim(C1), c(l, k)))
    aci_abort("aci_error_dims", "localization$C1 must have dimensions l x k.")
  if (!is.null(C1) && (!is.numeric(C1) || any(!is.finite(C1))))
    aci_abort("aci_error_dims", "localization$C1 must contain finite numeric values.")
  Yf <- array(NA_real_, c(l, N1, m)); Yf[, 1, ] <- Y
  for (kk in seq_len(N)) {
    t <- ob$t[kk]; xk <- ob$x[kk, ]; Dx <- ob$x[kk + 1, ] - xk
    Fy <- .eval_members(md$f, t, xk, Y, md$vectorized_members, k)     # k x m
    Sxk <- as.matrix(md$Sx(t, xk)); Gx <- Sxk %*% t(Sxk)
    Gi <- rs$ginv(Gx)
    dxi <- dt * Fy + sq * (Sxk %*% matrix(W[, kk, ], dxn, m))          # perturbed sim obs
    Ya <- Y - rowMeans(Y); Fa <- Fy - rowMeans(Fy)
    Pyf <- tcrossprod(Ya, Fa) / (m - 1)
    if (!is.null(C1)) Pyf <- C1 * Pyf
    K <- Pyf %*% Gi
    gv <- .eval_members(md$g, t, xk, Y, md$vectorized_members, l)
    Syk <- .ensemble_signal_diffusion(md, t, xk, Y)
    Y <- Y + dt * gv + sq * (Syk %*% matrix(B[, kk, ], dyn, m)) +
         K %*% (Dx - dxi)
    if (inflation > 1) Y <- apply_inflation(Y, inflation)
    if (any(!is.finite(Y)))
      aci_abort("aci_error_ensemble_divergence", sprintf(
        "EnKBF diverged at step %d; try modest inflation, smaller localization radius, or larger m (jiang2026enkbs Table 2 pattern).", kk))
    Yf[, kk + 1, ] <- Y
  }
  path <- new_ens_path(ob$t, Yf, "filter",
                       meta = list(m = m, inflation = inflation,
                                   localization = !is.null(localization),
                                   localization_spec = localization,
                                   nontarget = rs$tag,
                                   noise_fingerprint = noise$fingerprint,
                                   obs_x = ob$x,
                                   resolved_model = md,
                                   source_model = model,
                                   source_obs_x = obs$x,
                                   source_nontarget = nontarget,
                                   engine = "ensemble"))
  list(path = path, noise = noise, model = md)
}


#' @describeIn da_filter Ensemble filter for a general stochastic model.
#' @param m Ensemble size; a positive whole number.
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param inflation Positive multiplicative variance inflation factor.
#' @param seed Optional non-negative whole number seeding the generator.
#' @export
da_filter.stochastic_model <- function(model, obs, m = 100, localization = NULL,
                                       inflation = 1, seed = NULL,
                                       nontarget = NULL, ...) {
  fr <- enkbf(model, obs, m = m, localization = localization, inflation = inflation,
              seed = seed, nontarget = nontarget, ...)
  fr$path$meta$noise_store <- fr$noise
  fr$path
}


################################################################################
# EnKBS: ensemble Kalman-Bucy smoother (forward-noise reuse, hash-enforced)
################################################################################

#' Ensemble Kalman-Bucy smoother
#'
#' Ensemble Kalman-Bucy smoothing following the independent-noise jiang2026enkbs
#' recursion. `enkbs()` must receive the exact forward Wiener increments
#' returned by [enkbf()].
#'
#' @param model A `stochastic_model` object.
#' @param filter A `da_path_ensemble` object of kind `"filter"`.
#' @param noise The Wiener increments of the forward run, which must be the
#'   ones `filter` was driven with.
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param keep_members `TRUE` to retain the member paths.
#' @param ... Must be empty; unused arguments are an error.
#' @returns A `da_path_ensemble` object of kind `"smoother"`.
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#' @seealso [enkbf()], [ensemble_lag_table()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' fr <- enkbf(m, ob, m = 20, seed = 1)
#' enkbs(m, fr$path, fr$noise)
#'
#' @export
enkbs <- function(model, filter, noise, localization = NULL,
                  keep_members = TRUE, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to enkbs().")
  if (!inherits(filter, "da_path_ensemble") || !inherits(noise, "noise_store"))
    aci_abort("aci_error_noise_mismatch",
              "filter and noise must come from the same enkbf() result.")
  if (!identical(filter$kind, "filter"))
    aci_abort("aci_error_noise_mismatch",
              "enkbs() requires a forward ensemble filter path, not a smoother or other path kind.")
  if (!isTRUE(keep_members))
    aci_abort("aci_error_not_implemented",
              "keep_members = FALSE is unsupported because the smoother result is an ensemble path.")
  if (!identical(filter$meta$noise_fingerprint, .noise_fingerprint(noise)))
    aci_abort("aci_error_noise_mismatch",
      "The noise store does not match this filter run; EnKBS must reuse the forward Wiener paths (jiang2026enkbs App. A). Re-run enkbf().")
  d <- dim(filter$members); l <- d[1]; N1 <- d[2]; m <- d[3]
  dt <- filter$t[2] - filter$t[1]; sq <- sqrt(dt)
  stored_localization <- filter$meta$localization_spec
  if (is.null(localization) && !is.null(stored_localization))
    localization <- stored_localization
  if (!is.null(localization) && !is.list(localization))
    aci_abort("aci_error_dims", "localization must be NULL or a localization_spec-like list.")
  if (!is.null(stored_localization) && !identical(localization, stored_localization))
    aci_abort("aci_error_dims",
              "EnKBS must reuse the same localization specification as EnKBF.")
  C2 <- localization$C2
  if (!is.null(C2) && !identical(dim(C2), c(l, l)))
    aci_abort("aci_error_dims", "localization$C2 must have dimensions l x l.")
  if (!is.null(C2) && (!is.numeric(C2) || any(!is.finite(C2))))
    aci_abort("aci_error_dims", "localization$C2 must contain finite numeric values.")
  if (!is.null(C2) && max(abs(C2 - t(C2))) >
      1e-12 * max(1, max(abs(C2))))
    aci_abort("aci_error_dims", "localization$C2 must be symmetric.")
  if (m <= l) {
    c2_ch <- if (is.null(C2)) NULL else
      tryCatch(chol(sym(C2)), error = function(e) NULL)
    if (is.null(c2_ch) || !is.finite(rcond(C2)) || rcond(C2) < 1e-12)
      aci_abort("aci_error_ensemble_rank",
                "EnKBS with m <= l requires a numerically positive-definite C2 localization taper.")
  }
  Ys <- array(NA_real_, c(l, N1, m))
  Ys[, N1, ] <- filter$members[, N1, ]                 # invariant #3
  Ycur <- matrix(Ys[, N1, ], l, m)
  obx <- filter$meta$obs_x
  if (is.null(obx)) aci_abort("aci_error_dims", "Filter path lacks obs_x metadata; re-run enkbf().")
  source_model <- filter$meta$source_model
  resolved_model <- filter$meta$resolved_model
  if (!is.null(source_model) && !identical(model, source_model) &&
      !identical(model, resolved_model))
    aci_abort("aci_error_model_contract",
              "The supplied model does not match the ensemble filter run.")
  model <- resolved_model %||% model
  if (!inherits(model, "stochastic_model") || ncol(as.matrix(obx)) != model$k)
    aci_abort("aci_error_model_contract", "Stored filter metadata is incompatible with the stochastic model.")
  ob_contract <- observed_trajectory(filter$t, obx)
  .validate_ensemble_contract(model, ob_contract, Ycur)
  dyn <- ncol(as.matrix(model$Sy(filter$t[1], obx[1, ], Ycur[, 1])))
  if (!identical(dim(noise$B), c(dyn, N1 - 1L, m)))
    aci_abort("aci_error_noise_mismatch", "Stored signal noise has incompatible dimensions.")
  for (kk in (N1 - 1):1) {
    t1 <- filter$t[kk + 1]
    x1 <- obx[kk + 1, ]
    Yf1 <- matrix(filter$members[, kk + 1, ], l, m)
    Pf <- tcrossprod(Yf1 - rowMeans(Yf1)) / (m - 1)
    # jiang2026enkbs uses sample covariance, 1/(m-1), in both the filter and
    # smoother, as does its published MATLAB (EnKBS repository).
    if (!is.null(C2)) Pf <- C2 * Pf
    Pf <- sym(Pf)
    gv <- .eval_members(model$g, t1, x1, Ycur, model$vectorized_members, l)
    Syk <- .ensemble_signal_diffusion(model, t1, x1, Ycur)
    Gy <- Syk %*% t(Syk)
    # Localization can make the Schur product C2 .* Pf positive definite even
    # when the raw rank is at most m-1 < l (the high-dimensional
    # jiang2026enkbs regime).
    # Without it, report the rank limitation instead of silently regularizing.
    # This inverse is part of the EnKBS equation.  `chol_solve()` deliberately
    # jitters generic covariance inputs, but doing that here would silently
    # change the smoother and let rank-deficient tapers bypass the m <= l
    # contract.  Require an unregularized, numerically nonsingular factor.
    chPf <- tryCatch(chol(Pf), error = function(e) NULL)
    if (is.null(chPf))
      aci_abort("aci_error_ensemble_rank",
        "Localized/filter covariance is not invertible; provide a positive-definite C2 taper or increase m (jiang2026enkbs s2.4).")
    rhs <- forwardsolve(t(chPf), Ycur - Yf1)
    corr <- Gy %*% backsolve(chPf, rhs)
    Ycur <- Ycur - dt * gv - sq * (Syk %*% matrix(noise$B[, kk, ], dyn, m)) -
            dt * corr
    if (any(!is.finite(Ycur)))
      aci_abort("aci_error_ensemble_divergence",
                sprintf("EnKBS diverged at backward step %d.", kk))
    Ys[, kk, ] <- Ycur
  }
  new_ens_path(filter$t, Ys, "smoother",
               meta = c(filter$meta[c("m", "nontarget")], list(engine = "ensemble")))
}


#' Truncate an ensemble run (causal forward EnKBF) to an observation horizon (internal)
#'
#' @param run A stored forward ensemble run.
#' @param n Integer number of time points to retain.
#' @returns The run truncated to its first `n` time points.
#' @noRd
.truncate_ensemble_run <- function(run, n) {
  if (!is.list(run) || !inherits(run$path, "da_path_ensemble") ||
      !inherits(run$noise, "noise_store"))
    aci_abort("aci_error_noise_mismatch",
              "Internal ensemble truncation requires one intact enkbf() result.")
  N1 <- length(run$path$t)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n != floor(n) || n < 2L || n > N1)
    aci_abort("aci_error_dims", "Ensemble truncation index must be between 2 and N+1.")
  n <- as.integer(n)
  ix <- seq_len(n)
  jx <- seq_len(n - 1L)
  noise <- structure(
    list(B = run$noise$B[, jx, , drop = FALSE],
         W = run$noise$W[, jx, , drop = FALSE],
         seed = run$noise$seed),
    class = "noise_store")
  noise$fingerprint <- .noise_fingerprint(noise)
  path <- run$path
  path$t <- path$t[ix]
  path$members <- path$members[, ix, , drop = FALSE]
  path$meta$obs_x <- path$meta$obs_x[ix, , drop = FALSE]
  if (!is.null(path$meta$source_obs_x))
    path$meta$source_obs_x <- path$meta$source_obs_x[ix, , drop = FALSE]
  path$meta$noise_fingerprint <- noise$fingerprint
  list(path = path, noise = noise, model = run$model)
}


#' Build the ensemble lag table using an ensemble run (internal)
#'
#' Uses EnKBS as the posterior engine. The causal EnKBF prefix is reusable, so
#' only the backward pass is repeated here.
#'
#' @param model A `stochastic_model` object.
#' @param obs An observed trajectory.
#' @param run A stored forward ensemble run.
#' @param smoother The complete-record ensemble smoother path.
#' @param nontarget Optional non-target tag recorded on the table.
#' @param localization Optional localization specification.
#' @returns An object of class `lag_table` with an ensemble engine tag.
#' @noRd
.ensemble_lag_table_from_run <- function(model, obs, run, smoother,
                                         nontarget = NULL,
                                         localization = NULL) {
  obs <- as_obs(obs)
  if (!inherits(model, "stochastic_model") || obs$k != model$k ||
      !inherits(run$path, "da_path_ensemble") ||
      !inherits(smoother, "da_path_ensemble"))
    aci_abort("aci_error_model_contract",
              "An ensemble lag table needs matching model, observations, filter, and smoother.")
  d <- dim(run$path$members)
  if (d[3] <= d[1])
    aci_abort("aci_error_ensemble_rank", paste(
      "Ensemble CIR uses full-dimensional Gaussian moment KL and therefore",
      "requires ensemble size m > l."))
  if (!identical(dim(smoother$members), d) ||
      length(run$path$t) != length(obs$t) ||
      max(abs(run$path$t - obs$t)) > 1e-10 * max(1, max(abs(obs$t))))
    aci_abort("aci_error_dims",
              "The ensemble filter/smoother paths do not match the requested observation grid.")
  N1 <- length(obs$t)
  if (N1 > 500L)
    aci_warn("aci_warn_cost", paste(
      "A full ensemble lag table repeats an EnKBS backward pass for every",
      "observation horizon and has O(N^2) time-point work and storage."))

  filt_g <- as_gaussian(run$path)
  full_g <- as_gaussian(smoother)
  dec <- gaussian_kl_path(full_g, filt_g, decompose = TRUE)
  rows <- lapply(seq_len(N1), function(j) rep(NA_real_, N1 - j + 1L))
  for (j in seq_len(N1)) {
    rows[[j]][1] <- dec$total[j]
    rows[[j]][length(rows[[j]])] <- 0
  }

  # At n = j, the lagged smoother endpoint is the filter, already stored in
  # the first cell. At n = N, it is the complete smoother, so the final cell is
  # exactly zero. Only strict intermediate horizons require another pass.
  if (N1 > 2L) for (n in 2L:(N1 - 1L)) {
    prefix <- .truncate_ensemble_run(run, n)
    sm_n <- enkbs(prefix$model, prefix$path, prefix$noise,
                  localization = localization)
    sm_n_g <- as_gaussian(sm_n)
    if (n > 1L) for (j in seq_len(n - 1L))
      rows[[j]][n - j + 1L] <- unname(gaussian_kl(
        full_g$mean[j, ], full_g$cov[, , j],
        sm_n_g$mean[j, ], sm_n_g$cov[, , j], decompose = FALSE))
  }
  if (any(!is.finite(unlist(rows, use.names = FALSE))))
    aci_abort("aci_error_internal",
              "The ensemble lag table construction left a non-finite cell.")

  structure(
    list(t = obs$t, dt = obs$dt, mode = "full",
         diag = dec$total, rows = rows, L = N1 - seq_len(N1),
         diag_signal = dec$signal, diag_dispersion = dec$dispersion,
         tailbnd = rep(0, N1), onelag = NULL,
         meta = list(
           engine = "ensemble", nontarget = nontarget,
           tol = 0, window = Inf, max_lag = Inf, init = NULL,
           source_model = model, source_obs_x = obs$x,
           reference_smoother = "enkbs_full_horizon",
           m = d[3], localization = localization,
           repeated_backward_passes = max(N1 - 2L, 0L),
           source_status = paste(
             "andreou2026cir forward lagged-posterior definition evaluated with",
             "the jiang2026enkbs EnKBS; graded against a literal port of the",
             "published jiang2026enkbs dyad experiment on identical increments."))),
    class = "lag_table")
}


#' Forward CIR ensemble lag table
#'
#' Ensemble lag table for forward causal influence ranges.
#' Evaluates complete-versus-lagged smoother divergence using EnKBF/EnKBS
#' posterior approximation. `ensemble_lag_table()` constructs andreou2026cir's
#' forward finite-lag posterior family by repeating the jiang2026enkbs backward
#' pass at every observation horizon. It has quadratic time-point work/storage
#' and uses full-dimensional Gaussian moment KL, so `m` must exceed the hidden
#' dimension.
#'
#' This route is graded at machine precision against a literal R port of the
#' published jiang2026enkbs dyad experiment, driven by identical increments
#' (see `tests/testthat/helper-golden-p3.R`).
#'
#' @param model A `stochastic_model` object.
#' @param obs An observed trajectory, or anything [as_obs()] accepts.
#' @param m Ensemble size; must exceed the hidden dimension.
#' @param seed Optional non-negative whole number seeding the generator.
#' @param nontarget Optional `nontarget_spec`; see [nontarget()].
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param inflation Positive multiplicative variance inflation factor.
#' @param ic_sampler Optional function of the ensemble size drawing the
#'   initial members, or a list with the initial `mean` and `cov`; `NULL`
#'   uses the model's default.
#' @param ... Must be empty; unused arguments are an error.
#' @returns An object of class `lag_table` built by the ensemble engine.
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
#' @seealso [enkbs()], [forward_cir()]
#'
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' ensemble_lag_table(m, ob, m = 20, seed = 1)
#' }
#'
#' @export
ensemble_lag_table <- function(model, obs, m = 100, seed = NULL,
                               nontarget = NULL, localization = NULL,
                               inflation = 1, ic_sampler = NULL, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to ensemble_lag_table().")
  obs <- as_obs(obs)
  if (!inherits(model, "stochastic_model") || obs$k != model$k)
    aci_abort("aci_error_model_contract",
              "ensemble_lag_table() needs a matching stochastic model and observations.")
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) ||
      m != floor(m) || m <= model$l)
    aci_abort("aci_error_ensemble_rank",
              "Ensemble CIR requires an integer ensemble size m > l.")
  run <- enkbf(model, obs, m = m, seed = seed, nontarget = nontarget,
               localization = localization, inflation = inflation,
               ic_sampler = ic_sampler)
  sm <- enkbs(run$model, run$path, run$noise, localization = localization)
  aci_warn("aci_warn_ensemble_kl", sprintf(
    "Forward CIR uses Gaussian ensemble moments (m = %d); assess Monte Carlo stability across seeds and ensemble sizes.",
    as.integer(m)))
  .ensemble_lag_table_from_run(model, obs, run, sm,
                               nontarget = nontarget,
                               localization = localization)
}


#' @describeIn da_smooth Ensemble smoother for a general stochastic model.
#' @param m Ensemble size; a positive whole number.
#' @param localization Optional localization specification; see
#'   [localization_spec()].
#' @param inflation Positive multiplicative variance inflation factor.
#' @param seed Optional non-negative whole number seeding the generator.
#' @param ic_sampler Optional function of the ensemble size drawing the
#'   initial members, or a list with the initial `mean` and `cov`; `NULL`
#'   uses the model's default.
#' @param noise Optional stored Wiener increments of the forward run.
#' @export
da_smooth.stochastic_model <- function(model, obs, filter = NULL, m = 100,
                                       localization = NULL, inflation = 1,
                                       seed = NULL, nontarget = NULL,
                                       ic_sampler = NULL, noise = NULL, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to da_smooth().")
  ob_requested <- as_obs(obs)
  if (!inherits(model, "stochastic_model") || ob_requested$k != model$k)
    aci_abort("aci_error_model_contract",
              "model must match the requested observation dimension.")
  if (is.null(filter)) {
    filter <- enkbf(model, ob_requested, m = m, localization = localization,
                    inflation = inflation, seed = seed, nontarget = nontarget,
                    ic_sampler = ic_sampler, noise = noise)
  }
  path_to_check <- if (inherits(filter, "da_path_ensemble")) filter else filter$path
  if (!inherits(path_to_check, "da_path_ensemble"))
    aci_abort("aci_error_noise_mismatch", "filter must be an ensemble path or full enkbf() result.")
  same_grid <- length(path_to_check$t) == length(ob_requested$t) &&
    max(abs(path_to_check$t - ob_requested$t)) <=
      1e-10 * max(1, max(abs(ob_requested$t)))
  same_obs <- identical(dim(path_to_check$meta$source_obs_x), dim(ob_requested$x)) &&
    all(abs(path_to_check$meta$source_obs_x - ob_requested$x) <=
          1e-12 * pmax(1, abs(ob_requested$x)))
  if (!same_grid || !isTRUE(same_obs))
    aci_abort("aci_error_dims",
              "The supplied ensemble filter was computed from different observations.")
  if (!identical(path_to_check$meta$source_nontarget %||% NULL,
                 nontarget %||% NULL))
    aci_abort("aci_error_nontarget",
              "The supplied ensemble filter used a different nontarget specification.")
  if (!identical(path_to_check$meta$source_model, model))
    aci_abort("aci_error_model_contract",
              "The supplied ensemble filter was computed with a different model object.")
  if (inherits(filter, "da_path_ensemble")) {
    noise <- filter$meta$noise_store
    if (is.null(noise))
      aci_abort("aci_error_noise_mismatch", paste(
        "This ensemble filter path does not carry its forward noise store;",
        "pass the full enkbf() result or a path returned by da_filter()."))
    return(enkbs(filter$meta$resolved_model %||% model, filter, noise,
                 localization = localization))
  }
  if (!is.list(filter) || is.null(filter$noise) || is.null(filter$path))
    aci_abort("aci_error_noise_mismatch",
      "Pass the full enkbf() result (path + noise) so the backward pass can reuse the forward Wiener paths.")
  enkbs(filter$model %||% model, filter$path, filter$noise, localization = localization)
}
