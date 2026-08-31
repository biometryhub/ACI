################################################################################
## aci-conditional.R - one-pass conditioning for compiled CGNS paths
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
#' Models created by [aci_model_from_affine()] define their coefficient
#' matrices by evaluating the supplied affine drifts at zero and at each
#' hidden basis vector.  The ordinary coefficient closures repeat the
#' zero-state evaluation for `Lx`/`fx` and `Ly`/`fy`.  This batch realiser
#' performs the same ordered evaluations once per grid point and writes the
#' neutral dense contract.
#'
#' The drift evaluations, their order, and the Gram arithmetic are exactly
#' those of the earlier per-step form; what is removed is the scaffolding
#' around them.  Every array is filled through contiguous writes into one flat
#' numeric vector per quantity and given its dimensions once at the end, the
#' basis vectors and the coefficient closures are bound once for the whole
#' grid, the affine differencing is applied to a whole coefficient block at a
#' time rather than column by column, and the cross-noise test is a single
#' whole-path reduction over the realised arrays instead of a per-step call.
#'
#' @param model An authenticated `aci_model_from_affine()` model.
#' @param obs An `obs_traj` matching the model.
#' @param batch `FALSE` realises point by point even when the constructor
#'   declared a whole-path realisation of the same expressions. The two must
#'   agree bit for bit, and that is what the argument exists to check.
#' @returns Realised coefficient arrays without an observation precision path.
#' @noRd
.realise_affine_cgns_grid_once <- function(model, obs, batch = TRUE) {
  descriptor <- .cgns_realizer_descriptor(model)
  if (is.null(descriptor) ||
      !(identical(descriptor$id, "affine_model_v1") ||
        identical(descriptor$id, "enso6_aci_code_v1")) ||
      !is.function(descriptor$spec$f_full) ||
      !is.function(descriptor$spec$g_full))
    aci_abort(
      "aci_error_compiled_contract",
      paste("The affine batch realiser requires an unmodified",
            "aci_model_from_affine() model.")
    )

  ## A constructor may also declare a whole-path realisation of the same
  ## expressions.  Everything else, a declaration that does not check out
  ## included, is realised one grid point at a time below.
  e6 <- descriptor$spec$enso6
  if (isTRUE(batch) && identical(descriptor$id, "enso6_aci_code_v1") &&
      is.list(e6) && is.function(e6$drift) && is.function(e6$sd) &&
      is.integer(e6$obs_i) && is.integer(e6$hid) &&
      is.integer(e6$n_state) && length(e6$n_state) == 1L &&
      identical(length(e6$obs_i), model$k) &&
      identical(length(e6$hid), model$l) &&
      identical(sort(c(e6$obs_i, e6$hid)), seq_len(e6$n_state)) &&
      identical(obs$k, model$k))
    return(.realise_enso6_cgns_grid_once(e6, obs, model$k, model$l))

  f_full <- descriptor$spec$f_full
  g_full <- descriptor$spec$g_full
  Sx1_fun <- model$Sx1
  Sx2_fun <- model$Sx2
  Sy1_fun <- model$Sy1
  Sy2_fun <- model$Sy2
  tgrid <- obs$t
  xgrid <- obs$x
  N1 <- length(tgrid)
  k <- model$k
  l <- model$l
  zero <- numeric(l)
  ## One basis vector per hidden component, reused over the whole grid: a
  ## coefficient function cannot modify its caller's argument, so re-forming
  ## these per step buys nothing.
  basis <- lapply(seq_len(l), function(i) replace(zero, i, 1))
  dim_msg <- "CGNS affine drifts have incompatible dimensions at index %d."

  kl <- k * l
  ll <- l * l
  kk <- k * k
  lk <- l * k
  vLx <- numeric(kl * N1)
  vfx <- numeric(k * N1)
  vLy <- numeric(ll * N1)
  vfy <- numeric(l * N1)
  vgxx <- numeric(kk * N1)
  vgyy <- numeric(ll * N1)
  vgyx <- numeric(lk * N1)
  block_x <- numeric(kl)
  block_y <- numeric(ll)

  for (j in seq_len(N1)) {
    t <- tgrid[j]
    x <- xgrid[j, ]
    base_f <- as.numeric(f_full(t, x, zero))
    base_g <- as.numeric(g_full(t, x, zero))
    if (length(base_f) != k || length(base_g) != l)
      aci_abort("aci_error_model_contract", sprintf(dim_msg, j))
    at <- 0L
    for (i in seq_len(l)) {
      col <- as.numeric(f_full(t, x, basis[[i]]))
      if (length(col) != k)
        aci_abort("aci_error_model_contract", sprintf(dim_msg, j))
      block_x[(at + 1L):(at + k)] <- col
      at <- at + k
    }
    at <- 0L
    for (i in seq_len(l)) {
      col <- as.numeric(g_full(t, x, basis[[i]]))
      if (length(col) != l)
        aci_abort("aci_error_model_contract", sprintf(dim_msg, j))
      block_y[(at + 1L):(at + l)] <- col
      at <- at + l
    }
    ## Each column of the block is differenced against the same zero-state
    ## drift, so the whole block recycles that vector in one subtraction.
    Lx <- block_x - base_f
    Ly <- block_y - base_g
    dim(Lx) <- c(k, l)
    dim(Ly) <- c(l, l)

    Sx1 <- Sx1_fun(t, x); if (!is.matrix(Sx1)) Sx1 <- as.matrix(Sx1)
    Sx2 <- Sx2_fun(t, x); if (!is.matrix(Sx2)) Sx2 <- as.matrix(Sx2)
    Sy1 <- Sy1_fun(t, x); if (!is.matrix(Sy1)) Sy1 <- as.matrix(Sy1)
    Sy2 <- Sy2_fun(t, x); if (!is.matrix(Sy2)) Sy2 <- as.matrix(Sy2)
    if (nrow(Sx1) != k || nrow(Sx2) != k ||
        nrow(Sy1) != l || nrow(Sy2) != l ||
        ncol(Sx1) != ncol(Sy1) || ncol(Sx2) != ncol(Sy2))
      aci_abort("aci_error_model_contract", sprintf(
        "CGNS diffusion matrices have incompatible dimensions at index %d.", j
      ))
    tSx1 <- t.default(Sx1)
    tSx2 <- t.default(Sx2)
    gxx <- Sx1 %*% tSx1 + Sx2 %*% tSx2
    gyy <- Sy1 %*% t.default(Sy1) + Sy2 %*% t.default(Sy2)
    gyx <- Sy1 %*% tSx1 + Sy2 %*% tSx2
    if (any(!is.finite(c(Lx, base_f, Ly, base_g, gxx, gyy, gyx))))
      aci_abort("aci_error_model_contract", sprintf(
        "CGNS coefficients contain non-finite values at index %d.", j
      ))

    at <- (j - 1L) * kl; vLx[(at + 1L):(at + kl)] <- Lx
    at <- (j - 1L) * k;  vfx[(at + 1L):(at + k)] <- base_f
    at <- (j - 1L) * ll; vLy[(at + 1L):(at + ll)] <- Ly
    at <- (j - 1L) * l;  vfy[(at + 1L):(at + l)] <- base_g
    at <- (j - 1L) * kk; vgxx[(at + 1L):(at + kk)] <- gxx
    at <- (j - 1L) * ll; vgyy[(at + 1L):(at + ll)] <- gyy
    at <- (j - 1L) * lk; vgyx[(at + 1L):(at + lk)] <- gyx
  }

  dim(vLx) <- c(k, l, N1)
  dim(vLy) <- c(l, l, N1)
  dim(vgxx) <- c(k, k, N1)
  dim(vgyy) <- c(l, l, N1)
  dim(vgyx) <- c(l, k, N1)
  ## fx and fy are filled one grid point at a time, which is contiguous in the
  ## k by N1 orientation and strided in the N1 by k one the contract wants.
  dim(vfx) <- c(k, N1)
  dim(vfy) <- c(l, N1)
  out <- list(
    Lx = vLx,
    fx = t.default(vfx),
    Ly = vLy,
    fy = t.default(vfy),
    gxx = vgxx,
    gyy = vgyy,
    gyx = vgyx
  )
  attr(out, "path_cross_noise") <- .compiled_path_has_cross_noise(out)
  out
}


#' Realise a constructor-authenticated ENSO6 path in one whole-path pass
#' (internal)
#'
#' [aci_enso_model()] writes its drift and its noise amplitudes as expressions
#' that hold elementwise, so the whole observation grid can be pushed through
#' them at
#' once.  This realiser evaluates exactly what the per-point coefficient
#' closures evaluate, in the same order and with the same associations: the
#' affine columns are still the drift at each hidden basis vector differenced
#' against the drift at the zero hidden state, and the Gram arrays are still
#' what the constructor's diagonal diffusion blocks and its auto-generated zero
#' cross-channels multiply out to, which is the squared amplitudes on the
#' diagonal and exact zeros elsewhere.  What is removed is the per-grid-point
#' call, not any arithmetic.
#'
#' @param spec The `enso6` element of an authenticated realiser descriptor.
#' @param obs An `obs_traj` matching the model.
#' @param k Observed dimension.
#' @param l Hidden dimension.
#' @returns Realised coefficient arrays without an observation precision path.
#' @noRd
.realise_enso6_cgns_grid_once <- function(spec, obs, k, l) {
  tgrid <- obs$t
  xgrid <- obs$x
  N1 <- length(tgrid)
  obs_i <- spec$obs_i
  hid <- spec$hid
  ns <- spec$n_state
  ## The full state at every grid point with the hidden components at zero:
  ## what `asm(x, numeric(l))` builds per point, and what the diffusion blocks
  ## see, since they are given the observed components only.
  state <- rep(list(numeric(N1)), ns)
  for (i in seq_len(k)) state[[obs_i[i]]] <- xgrid[, i]
  drift_grid <- function(s) {
    d <- as.numeric(spec$drift(tgrid, s))
    if (length(d) != ns * N1)
      aci_abort("aci_error_model_contract",
                "The ENSO6 drift did not return one value per state and grid point.")
    dim(d) <- c(N1, ns)
    d
  }

  base <- drift_grid(state)
  fx <- base[, obs_i, drop = FALSE]
  fy <- base[, hid, drop = FALSE]
  Lx <- numeric(k * l * N1); dim(Lx) <- c(k, l, N1)
  Ly <- numeric(l * l * N1); dim(Ly) <- c(l, l, N1)
  unit <- rep(1, N1)
  for (i in seq_len(l)) {
    lifted <- state
    lifted[[hid[i]]] <- unit
    col <- drift_grid(lifted)
    Lx[, i, ] <- t.default(col[, obs_i, drop = FALSE] - fx)
    Ly[, i, ] <- t.default(col[, hid, drop = FALSE] - fy)
  }

  sd <- as.numeric(spec$sd(tgrid, state))
  if (length(sd) != ns * N1)
    aci_abort("aci_error_model_contract",
              "The ENSO6 diffusion did not return one amplitude per state and grid point.")
  dim(sd) <- c(N1, ns)
  gxx <- numeric(k * k * N1); dim(gxx) <- c(k, k, N1)
  for (i in seq_len(k)) { a <- sd[, obs_i[i]]; gxx[i, i, ] <- a * a }
  gyy <- numeric(l * l * N1); dim(gyy) <- c(l, l, N1)
  for (i in seq_len(l)) { a <- sd[, hid[i]]; gyy[i, i, ] <- a * a }
  gyx <- numeric(l * k * N1); dim(gyx) <- c(l, k, N1)

  if (!all(is.finite(Lx), is.finite(fx), is.finite(Ly), is.finite(fy),
           is.finite(gxx), is.finite(gyy))) {
    ## Exceptional, so the index may be found the slow way: the message names
    ## the grid point the per-point realiser would have stopped at.
    j <- 1L
    while (j <= N1 && all(is.finite(Lx[, , j]), is.finite(fx[j, ]),
                          is.finite(Ly[, , j]), is.finite(fy[j, ]),
                          is.finite(gxx[, , j]), is.finite(gyy[, , j]))) j <- j + 1L
    aci_abort("aci_error_model_contract", sprintf(
      "CGNS coefficients contain non-finite values at index %d.", j
    ))
  }

  out <- list(Lx = Lx, fx = fx, Ly = Ly, fy = fy,
              gxx = gxx, gyy = gyy, gyx = gyx)
  attr(out, "path_cross_noise") <- .compiled_path_has_cross_noise(out)
  out
}


#' Compile an authenticated affine CGNS path (internal)
#'
#' @param model An authenticated `aci_model_from_affine()` model.
#' @param obs Anything accepted by [as_obs()].
#' @param conditional Optional `aci_conditional_spec`.
#' @returns A private `compiled_cgns` bundle.
#' @noRd
.compile_affine_cgns <- function(model, obs, conditional = NULL) {
  source_obs <- as_obs(obs)
  if (!inherits(model, "cgns_model") || source_obs$k != model$k)
    aci_abort(
      "aci_error_dims",
      "Observation dimension does not match the affine CGNS model."
    )
  full <- .realise_affine_cgns_grid_once(model, source_obs)
  .compile_cgns_complete(
    model, source_obs, conditional = conditional, full = full,
    realization = "affine_batch"
  )
}


#' Build a realised observation-precision path (internal)
#'
#' This is the single point at which conditional ACI's masked observation
#' precision is realised, and therefore the only place the first-slice
#' convention can be applied. See [aci_conditional()] for what `first_step`
#' selects.
#'
#' @param gxx Realised observation Gram arrays.
#' @param N Number of observation intervals.
#' @param target Optional target indices for masked-innovation conditioning.
#' @param first_step First-slice convention for the masked branch, `"uniform"`
#'   or `"matlab"`; ignored when `target` is `NULL`, which has no mask.
#' @returns An array with one precision matrix per interval start.
#' @noRd
.compiled_precision_path <- function(gxx, N, target = NULL,
                                     first_step = c("uniform", "matlab")) {
  first_step <- match.arg(first_step)
  k <- dim(gxx)[1L]
  kk <- k * k
  if (!is.null(target)) {
    ## `masked_ginv()` rejects these per slice; the path is one target, so the
    ## rejection is hoisted rather than repeated, and the fast pair below is
    ## then free to subset without re-checking.  Same class, same message.
    if (!is.numeric(target) || any(!is.finite(target)) ||
        any(target != floor(target)) || any(target < 1L) || any(target > k) ||
        anyDuplicated(target))
      aci_abort("aci_error_dims", "idxA contains invalid Gram-matrix indices.")
    target <- as.integer(target)
  }

  ## Two solver pairs computing the same thing.  The guarded pair IS the
  ## definition.  The fast pair is `chol_solve()` and `masked_ginv()` with
  ## their loop-invariant re-entry checks lifted out; sym(), chol() and, on the
  ## masked side, forwardsolve() and backsolve() are theirs unchanged, so a
  ## path whose Grams all factor gets the same bits without paying 4N
  ## validations for it.
  ##
  ## The unmasked inverse is taken with `chol2inv()` on that same factor rather
  ## than by two triangular solves against an identity: one `dpotri` instead of
  ## a `dtrsm` pair, and exactly symmetric where the triangular route is not.
  ## On a dense Gram the two routes differ in the last bit or two, which the
  ## tests budget; on a diagonal Gram, which is every Gram in the ACI_code
  ## scope, they agree to the bit EXCEPT that `dpotri` writes `-0.0` into the
  ## off-diagonal zeros.  `+ 0` normalises those back to `+0.0` and is a no-op
  ## on every other value, so the realised arrays stay bitwise what they were.
  ## `identical()` and a max-absolute-difference both miss a sign-of-zero flip,
  ## so it has to be removed at the source rather than tested for.
  rhs_full <- diag(k)
  fast_full <- function(gram) chol2inv(chol.default(sym(gram))) + 0
  safe_full <- function(gram) chol_solve(gram, rhs_full, "gxx")
  fast_part <- if (length(target) > 0L) {
    rhs_part <- diag(length(target))
    function(gram) {
      ch <- chol.default(sym(gram[target, target, drop = FALSE]))
      M <- matrix(0, k, k)
      M[target, target] <- backsolve(ch, forwardsolve(t.default(ch), rhs_part))
      M
    }
  } else function(gram) matrix(0, k, k)
  safe_part <- function(gram) masked_ginv(gram, target)

  fill <- function(inv_full, inv_part) {
    precision <- numeric(kk * N)
    inv <- if (is.null(target)) inv_full else inv_part
    from <- 1L
    if (!is.null(target) && identical(first_step, "matlab") && N >= 1L) {
      ## ACI_code fills S_xoS_x_inv(:,:,1) with the full inverse of the first
      ## Gram and then overwrites only the target block, which leaves the
      ## first slice unmasked; every later slice is written by the target-only
      ## assignment alone (ENSO_model_cond_ACI_h_W_unobs.m:1197 and :1202,
      ## against the commented-out :1250 refresh).
      gram <- gxx[, , 1L]
      dim(gram) <- c(k, k)
      precision[seq_len(kk)] <- inv_full(gram)
      from <- 2L
    }
    if (from <= N) for (j in from:N) {
      gram <- gxx[, , j]
      dim(gram) <- c(k, k)
      at <- (j - 1L) * kk
      precision[(at + 1L):(at + kk)] <- inv(gram)
    }
    precision
  }

  ## One guard for the whole path rather than one per slice: this loop carries
  ## no state between slices, so a single failure anywhere can simply refill it
  ## through the guarded pair, which reproduces the jitter ladder and the
  ## classed errors exactly.  The backward smoother recursion in
  ## `aci-kernels-matrix.R` guards per slice instead, because there a restart
  ## would be a restart of the recursion.
  precision <- tryCatch(fill(fast_full, fast_part), error = function(e) NULL)
  ## The refill also runs when the fast pair "succeeded" through a LAPACK
  ## that completes on non-finite input without signalling: a poisoned fill
  ## is a failed fill, and the guarded pair owns the classed errors.
  if (is.null(precision) || !all(is.finite(precision)))
    precision <- fill(safe_full, safe_part)
  dim(precision) <- c(k, k, N)
  precision
}


#' Construct the prescribed-forcing model identity without reevaluation
#' (internal)
#'
#' This is the structural part of [aci_conditional_reduce()].  Its coefficient
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
      ## Already reduced: the declaration must not survive onto the result,
      ## where it would ask for the same reduction a second time.
      estimand_nontarget = NULL,
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
#' Uses the same scale-free tolerance as `.has_cross_noise()`, applied to the
#' same per-slice maxima; the reduction is written over the whole path rather
#' than one slice at a time.
#'
#' @param coefficients Realised coefficient arrays.
#' @returns One logical value.
#' @noRd
.compiled_path_has_cross_noise <- function(coefficients) {
  N1 <- dim(coefficients$gxx)[3L]
  ## The per-slice maximum is taken by folding the slice's entries together
  ## with pmax, which is one vectorised pass per entry rather than one max()
  ## call per grid point.  Selection is exact, so this is the same maximum.
  slice_max <- function(a) {
    m <- matrix(a, length(a) / N1, N1)
    got <- abs(m[1L, ])
    for (i in seq_len(nrow(m))[-1L]) got <- pmax(got, abs(m[i, ]))
    got
  }
  cross <- slice_max(coefficients$gyx)
  scale <- sqrt(slice_max(coefficients$gxx) * slice_max(coefficients$gyy))
  tol <- 100 * .Machine$double.eps
  any(is.finite(cross) & is.finite(scale) &
      cross > tol * pmax(scale, .Machine$double.xmin))
}


#' Compile a complete conditioned CGNS path in one pass (internal)
#'
#' This generic compiler realises the original model exactly once per observed
#' grid point, then transforms the realised arrays for unconditioned,
#' masked-innovation (`"mask"`), or prescribed-forcing (`"reduce"`)
#' execution.  It does
#' not alter public routing.
#'
#' @param model A `cgns_model`.
#' @param obs Anything accepted by [as_obs()].
#' @param conditional Optional `aci_conditional_spec`.
#' @returns A private `compiled_cgns` bundle.
#' @noRd
.compile_cgns_complete <- function(
    model, obs, conditional = NULL, full = NULL,
    realization = "generic_closure_one_pass") {
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "compiled_cgns requires a cgns_model.")
  source_obs <- as_obs(obs)
  if (source_obs$k != model$k)
    aci_abort(
      "aci_error_dims",
      "Observation dimension does not match the CGNS model."
    )
  if (!is.null(conditional) && !inherits(conditional, "aci_conditional_spec"))
    aci_abort(
      "aci_error_nontarget",
      "conditional must be created by aci_conditional()."
    )

  N1 <- length(source_obs$t)
  N <- N1 - 1L
  if (is.null(full)) full <- .realise_cgns_grid_once(model, source_obs)

  if (is.null(conditional)) {
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
    ix <- .nt_indices(conditional, source_obs)
    if (identical(conditional$method, "reduce")) {
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
            "aci_conditional(method = 'mask') (SPEC-01 s6, pending SI",
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
        tag = conditional
      )
    } else {
      coefficients <- full
      coefficients$gxx_weight <- .compiled_precision_path(
        full$gxx, N, target = ix$A,
        first_step = conditional$first_step %||% "uniform"
      )
      rs <- list(
        model = model,
        obs = source_obs,
        likelihood_idx = ix$A,
        tag = conditional
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
