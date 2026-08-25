################################################################################
## compiled_conditioning.R - one-pass conditioning for compiled CGNS paths
################################################################################


#' Realise a complete CGNS coefficient path once (internal)
#'
#' The original, unreduced model is evaluated exactly once at every observed
#' grid point.  Conditioning is applied to these realised arrays afterwards,
#' so prescribed forcing does not repeat coefficient closure evaluation through
#' a reduced wrapper model.
#'
#' @param model A `cgns_model`.
#' @param obs An `obs_traj` with the model's full observed dimension.
#' @returns Realised coefficient arrays without an observation precision path.
#' @noRd
.realise_cgns_grid_once <- function(model, obs) {
  N1 <- length(obs$t)
  k <- model$k
  l <- model$l
  out <- list(
    Lx = array(NA_real_, c(k, l, N1)),
    fx = matrix(NA_real_, N1, k),
    Ly = array(NA_real_, c(l, l, N1)),
    fy = matrix(NA_real_, N1, l),
    gxx = array(NA_real_, c(k, k, N1)),
    gyy = array(NA_real_, c(l, l, N1)),
    gyx = array(NA_real_, c(l, k, N1))
  )
  path_cross_noise <- FALSE

  for (j in seq_len(N1)) {
    co <- eval_coefs(model, obs$t[j], obs$x[j, ])
    if (!identical(dim(co$Lx), c(k, l)) || length(co$fx) != k ||
        !identical(dim(co$Ly), c(l, l)) || length(co$fy) != l ||
        !identical(dim(co$gxx), c(k, k)) ||
        !identical(dim(co$gyy), c(l, l)) ||
        !identical(dim(co$gyx), c(l, k)) ||
        any(!is.finite(c(
          co$Lx, co$fx, co$Ly, co$fy, co$gxx, co$gyy, co$gyx
        ))))
      aci_abort("aci_error_model_contract", sprintf(
        "CGNS coefficients have incompatible dimensions or values at index %d.",
        j
      ))

    out$Lx[, , j] <- co$Lx
    out$fx[j, ] <- co$fx
    out$Ly[, , j] <- co$Ly
    out$fy[j, ] <- co$fy
    out$gxx[, , j] <- co$gxx
    out$gyy[, , j] <- co$gyy
    out$gyx[, , j] <- co$gyx
    path_cross_noise <- path_cross_noise || .has_cross_noise(co)
  }

  attr(out, "path_cross_noise") <- path_cross_noise
  out
}


#' Realise a constructor-authenticated affine CGNS path (internal)
#'
#' Models created by [cgns_from_affine()] define their coefficient matrices by
#' evaluating the supplied affine drifts at zero and at each hidden basis
#' vector.  The ordinary coefficient closures repeat the zero-state evaluation
#' for `Lx`/`fx` and `Ly`/`fy`.  This batch realiser performs the same ordered
#' evaluations once per grid point and writes the neutral dense contract.
#'
#' @param model An authenticated `cgns_from_affine()` model.
#' @param obs An `obs_traj` matching the model.
#' @returns Realised coefficient arrays without an observation precision path.
#' @noRd
.realise_affine_cgns_grid_once <- function(model, obs) {
  descriptor <- .cgns_realizer_descriptor(model)
  if (is.null(descriptor) ||
      !identical(descriptor$id, "affine_model_v1") ||
      !is.function(descriptor$spec$f_full) ||
      !is.function(descriptor$spec$g_full))
    aci_abort(
      "aci_error_compiled_contract",
      "The affine batch realiser requires an unmodified cgns_from_affine() model."
    )

  N1 <- length(obs$t)
  k <- model$k
  l <- model$l
  zero <- numeric(l)
  out <- list(
    Lx = array(NA_real_, c(k, l, N1)),
    fx = matrix(NA_real_, N1, k),
    Ly = array(NA_real_, c(l, l, N1)),
    fy = matrix(NA_real_, N1, l),
    gxx = array(NA_real_, c(k, k, N1)),
    gyy = array(NA_real_, c(l, l, N1)),
    gyx = array(NA_real_, c(l, k, N1))
  )
  path_cross_noise <- FALSE

  for (j in seq_len(N1)) {
    t <- obs$t[j]
    x <- obs$x[j, ]
    base_f <- as.numeric(descriptor$spec$f_full(t, x, zero))
    base_g <- as.numeric(descriptor$spec$g_full(t, x, zero))
    if (length(base_f) != k || length(base_g) != l)
      aci_abort("aci_error_model_contract", sprintf(
        "CGNS affine drifts have incompatible dimensions at index %d.", j
      ))
    Lx <- vapply(seq_len(l), function(i) {
      as.numeric(descriptor$spec$f_full(t, x, replace(zero, i, 1))) - base_f
    }, numeric(k))
    Ly <- vapply(seq_len(l), function(i) {
      as.numeric(descriptor$spec$g_full(t, x, replace(zero, i, 1))) - base_g
    }, numeric(l))
    Lx <- matrix(Lx, k, l)
    Ly <- matrix(Ly, l, l)

    Sx1 <- as.matrix(model$Sx1(t, x))
    Sx2 <- as.matrix(model$Sx2(t, x))
    Sy1 <- as.matrix(model$Sy1(t, x))
    Sy2 <- as.matrix(model$Sy2(t, x))
    if (nrow(Sx1) != k || nrow(Sx2) != k ||
        nrow(Sy1) != l || nrow(Sy2) != l ||
        ncol(Sx1) != ncol(Sy1) || ncol(Sx2) != ncol(Sy2))
      aci_abort("aci_error_model_contract", sprintf(
        "CGNS diffusion matrices have incompatible dimensions at index %d.", j
      ))
    gxx <- Sx1 %*% t(Sx1) + Sx2 %*% t(Sx2)
    gyy <- Sy1 %*% t(Sy1) + Sy2 %*% t(Sy2)
    gyx <- Sy1 %*% t(Sx1) + Sy2 %*% t(Sx2)
    if (any(!is.finite(c(Lx, base_f, Ly, base_g, gxx, gyy, gyx))))
      aci_abort("aci_error_model_contract", sprintf(
        "CGNS coefficients contain non-finite values at index %d.", j
      ))

    out$Lx[, , j] <- Lx
    out$fx[j, ] <- base_f
    out$Ly[, , j] <- Ly
    out$fy[j, ] <- base_g
    out$gxx[, , j] <- gxx
    out$gyy[, , j] <- gyy
    out$gyx[, , j] <- gyx
    path_cross_noise <- path_cross_noise ||
      .has_cross_noise(list(gxx = gxx, gyy = gyy, gyx = gyx))
  }

  attr(out, "path_cross_noise") <- path_cross_noise
  out
}


#' Compile an authenticated affine CGNS path (internal)
#'
#' @param model An authenticated `cgns_from_affine()` model.
#' @param obs Anything accepted by [as_obs()].
#' @param nontarget Optional `nontarget_spec`.
#' @returns A private `compiled_cgns` bundle.
#' @noRd
.compile_affine_cgns <- function(model, obs, nontarget = NULL) {
  source_obs <- as_obs(obs)
  if (!inherits(model, "cgns_model") || source_obs$k != model$k)
    aci_abort(
      "aci_error_dims",
      "Observation dimension does not match the affine CGNS model."
    )
  full <- .realise_affine_cgns_grid_once(model, source_obs)
  .compile_cgns_complete(
    model, source_obs, nontarget = nontarget, full = full,
    realization = "affine_batch"
  )
}


#' Build a realised observation-precision path (internal)
#'
#' @param gxx Realised observation Gram arrays.
#' @param N Number of observation intervals.
#' @param target Optional target indices for masked-innovation conditioning.
#' @returns An array with one precision matrix per interval start.
#' @noRd
.compiled_precision_path <- function(gxx, N, target = NULL) {
  k <- dim(gxx)[1L]
  precision <- array(NA_real_, c(k, k, N))
  for (j in seq_len(N)) {
    gram <- matrix(gxx[, , j], k, k)
    precision[, , j] <- if (is.null(target)) {
      chol_solve(gram, diag(k), "gxx")
    } else {
      masked_ginv(gram, target)
    }
  }
  precision
}


#' Construct the prescribed-forcing model identity without reevaluation
#' (internal)
#'
#' This is the structural part of [reduce_nontarget()].  Its coefficient
#' closures preserve that function's public model semantics, but construction
#' does not probe them: the full source model has already been validated and
#' its coefficient path has already been realised and checked.
#'
#' @param model Original full CGNS model.
#' @param obs Original observation path.
#' @param target Target observed indices.
#' @param prescribed Prescribed observed indices.
#' @param correlated Whether the realised reduced path has cross noise.
#' @returns A reduced `cgns_model` used as the bundle's resolved identity.
#' @noRd
.compiled_prescribed_model <- function(model, obs, target, prescribed,
                                       correlated) {
  t0 <- obs$t[1L]
  dt <- obs$dt
  XB <- obs$x[, prescribed, drop = FALSE]
  Nrow <- nrow(XB)
  lookupB <- function(t) {
    j <- as.integer(round((t - t0) / dt)) + 1L
    XB[min(max(j, 1L), Nrow), ]
  }
  assemble <- function(t, xA) {
    full <- numeric(length(target) + length(prescribed))
    full[target] <- xA
    full[prescribed] <- lookupB(t)
    full
  }
  mk <- function(fun, rows = NULL) {
    force(fun)
    force(rows)
    function(t, xA) {
      value <- fun(t, assemble(t, xA))
      if (is.null(rows)) return(value)
      if (is.matrix(value)) value[rows, , drop = FALSE] else value[rows]
    }
  }

  Lx <- mk(model$Lx, target)
  fx <- mk(model$fx, target)
  Ly <- mk(model$Ly)
  fy <- mk(model$fy)
  Sx1 <- mk(model$Sx1, target)
  Sx2 <- mk(model$Sx2, target)
  Sy1 <- mk(model$Sy1)
  Sy2 <- mk(model$Sy2)
  reduced <- list(
    Lx = Lx, fx = fx, Ly = Ly, fy = fy,
    Sx1 = Sx1, Sx2 = Sx2, Sy1 = Sy1, Sy2 = Sy2,
    k = as.integer(length(target)), l = model$l,
    name = paste0(model$name, "|reduced"),
    meta = utils::modifyList(model$meta, list(
      correlated_noise = isTRUE(correlated),
      nontarget_reduction = list(
        target = target,
        prescribed = prescribed
      )
    )),
    f = function(t, x, y) drop(fx(t, x) + as.matrix(Lx(t, x)) %*% y),
    g = function(t, x, y) drop(fy(t, x) + as.matrix(Ly(t, x)) %*% y),
    Sx = function(t, x) cbind(as.matrix(Sx1(t, x)), as.matrix(Sx2(t, x))),
    Sy = function(t, x, y) cbind(as.matrix(Sy1(t, x)), as.matrix(Sy2(t, x))),
    vectorized_members = FALSE
  )
  class(reduced) <- c("cgns_model", "stochastic_model")

  if (!is.null(model$meta$ic_default)) {
    ic <- model$meta$ic_default
    if (is.list(ic) && length(ic$x0) == model$k &&
        length(ic$y0) == model$l && all(is.finite(c(ic$x0, ic$y0)))) {
      reduced$meta$ic_default <- list(
        x0 = as.numeric(ic$x0)[target],
        y0 = as.numeric(ic$y0)
      )
    } else {
      reduced$meta$ic_default <- NULL
    }
  }
  reduced
}


#' Restrict full realised arrays to prescribed-forcing targets (internal)
#'
#' @param coefficients Full realised coefficient arrays.
#' @param target Target observed indices.
#' @param N1 Number of grid points.
#' @param l Hidden dimension.
#' @returns Coefficient arrays with the exact reduced dimensions.
#' @noRd
.compiled_prescribed_arrays <- function(coefficients, target, N1, l) {
  kA <- length(target)
  list(
    Lx = array(coefficients$Lx[target, , , drop = FALSE], c(kA, l, N1)),
    fx = matrix(coefficients$fx[, target, drop = FALSE], N1, kA),
    Ly = array(coefficients$Ly, c(l, l, N1)),
    fy = matrix(coefficients$fy, N1, l),
    gxx = array(
      coefficients$gxx[target, target, , drop = FALSE], c(kA, kA, N1)
    ),
    gyy = array(coefficients$gyy, c(l, l, N1)),
    gyx = array(coefficients$gyx[, target, , drop = FALSE], c(l, kA, N1))
  )
}


#' Detect cross noise on a realised coefficient path (internal)
#'
#' Uses the same scale-free tolerance as `.has_cross_noise()`.
#'
#' @param coefficients Realised coefficient arrays.
#' @returns One logical value.
#' @noRd
.compiled_path_has_cross_noise <- function(coefficients) {
  N1 <- dim(coefficients$gxx)[3L]
  k <- dim(coefficients$gxx)[1L]
  l <- dim(coefficients$gyy)[1L]
  any(vapply(seq_len(N1), function(j) {
    .has_cross_noise(list(
      gxx = matrix(coefficients$gxx[, , j], k, k),
      gyy = matrix(coefficients$gyy[, , j], l, l),
      gyx = matrix(coefficients$gyx[, , j], l, k)
    ))
  }, logical(1)))
}


#' Compile a complete conditioned CGNS path in one pass (internal)
#'
#' This generic compiler realises the original model exactly once per observed
#' grid point, then transforms the realised arrays for unconditioned,
#' masked-innovation (`"inflate"`), or prescribed-forcing execution.  It does
#' not alter public routing.
#'
#' @param model A `cgns_model`.
#' @param obs Anything accepted by [as_obs()].
#' @param nontarget Optional `nontarget_spec`.
#' @returns A private `compiled_cgns` bundle.
#' @noRd
.compile_cgns_complete <- function(
    model, obs, nontarget = NULL, full = NULL,
    realization = "generic_closure_one_pass") {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "compiled_cgns requires a cgns_model.")
  source_obs <- as_obs(obs)
  if (source_obs$k != model$k)
    aci_abort(
      "aci_error_dims",
      "Observation dimension does not match the CGNS model."
    )
  if (!is.null(nontarget) && !inherits(nontarget, "nontarget_spec"))
    aci_abort(
      "aci_error_nontarget",
      "nontarget must be created by nontarget()."
    )

  N1 <- length(source_obs$t)
  N <- N1 - 1L
  if (is.null(full)) full <- .realise_cgns_grid_once(model, source_obs)

  if (is.null(nontarget)) {
    coefficients <- full
    coefficients$gxx_weight <- .compiled_precision_path(full$gxx, N)
    rs <- list(
      model = model,
      obs = source_obs,
      likelihood_idx = seq_len(model$k),
      tag = NULL
    )
    correlated <- isTRUE(model$meta$correlated_noise) ||
      isTRUE(attr(full, "path_cross_noise"))
  } else {
    ix <- .nt_indices(nontarget, source_obs)
    if (identical(nontarget$strategy, "prescribed_forcing")) {
      bad_cross <- any(vapply(seq_len(N1), function(j) {
        gram <- full$gxx[, , j]
        max(abs(gram[ix$A, ix$B, drop = FALSE])) >
          1e-12 * max(abs(gram), 1e-300)
      }, logical(1)))
      if (bad_cross)
        aci_abort(
          "aci_error_nontarget_crossnoise",
          paste(
            "gxx has a nonzero A-B cross-block; use",
            "nontarget(strategy = 'inflate') (SPEC-01 s6, pending SI",
            "equivalence transcription)."
          )
        )

      coefficients <- .compiled_prescribed_arrays(
        full, ix$A, N1, model$l
      )
      coefficients$gxx_weight <- .compiled_precision_path(
        coefficients$gxx, N
      )
      correlated <- .compiled_path_has_cross_noise(coefficients)
      reduced_model <- .compiled_prescribed_model(
        model, source_obs, ix$A, ix$B, correlated
      )
      reduced_obs <- observed_trajectory(
        source_obs$t, source_obs$x[, ix$A, drop = FALSE]
      )
      rs <- list(
        model = reduced_model,
        obs = reduced_obs,
        likelihood_idx = seq_len(length(ix$A)),
        tag = nontarget
      )
    } else {
      coefficients <- full
      coefficients$gxx_weight <- .compiled_precision_path(
        full$gxx, N, target = ix$A
      )
      rs <- list(
        model = model,
        obs = source_obs,
        likelihood_idx = ix$A,
        tag = nontarget
      )
      correlated <- isTRUE(model$meta$correlated_noise) ||
        isTRUE(attr(full, "path_cross_noise"))
    }
  }

  # Attributes on the temporary coefficient list are not part of the neutral
  # compiled contract.
  attr(coefficients, "path_cross_noise") <- NULL
  .new_compiled_cgns(
    rs = rs,
    source_model = model,
    source_obs = source_obs,
    coefficients = coefficients,
    correlated_noise = correlated,
    realization = realization
  )
}
