################################################################################
## model_classes.R - model classes and simulation
## ########################################################################## ##
##
## Contents:
##   * the general stochastic_model class and validators:
##       - stochastic_model, validate_stochastic, print.stochastic_model
##
##   * the CGNS class: coefficients, Grams, validation, affine constructor:
##       - cgns_model, cgns_grams, eval_coefs, validate_cgns,
##         conditionally_linear_model, cgns_from_affine
##
##   * Euler-Maruyama simulation and reusable noise stores:
##       - simulate.stochastic_model, .simulate_one
##
################################################################################


################################################################################
# general stochastic_model class and validators
################################################################################

#' Validate and coerce state dimensions (internal)
#'
#' @param k Observed dimension; a positive whole number.
#' @param l Hidden dimension; a positive whole number.
#' @returns Named integer vector with elements `k` and `l`.
#' @noRd
.state_dims <- function(k, l) {
  valid_one <- function(z)
    is.numeric(z) && length(z) == 1L && is.finite(z) &&
    z >= 1 && z <= .Machine$integer.max && z == floor(z)
  if (!valid_one(k) || !valid_one(l))
    aci_abort("aci_error_model_contract", "k and l must be positive integers.")
  c(k = as.integer(k), l = as.integer(l))
}


#' General stochastic model
#'
#' Constructs the general (not necessarily conditional-Gaussian) stochastic
#' model consumed by the ensemble engine. The observed and hidden drifts and
#' diffusions are supplied as functions and validated on a small probe set
#' before the object is returned.
#'
#' @param f Function giving the observed drift at `(t, x, y)`.
#' @param g Function giving the hidden drift at `(t, x, y)`.
#' @param Sx Function giving the observed diffusion coefficient at `(t, x)`.
#' @param Sy Function giving the hidden diffusion coefficient at `(t, x, y)`.
#' @param k Observed dimension; a positive whole number.
#' @param l Hidden dimension; a positive whole number.
#' @param vectorized_members `TRUE` when the drifts accept a matrix of ensemble
#'   members rather than a single hidden state.
#' @param name Optional 1-length character label for the model.
#' @param meta Optional named list of metadata carried on the object.
#' @returns An object of class `stochastic_model`.
#'
#' @seealso [cgns_model()], [enkbf()]
#'
#' @examples
#' stochastic_model(
#'   f = function(t, x, y) -x + y,
#'   g = function(t, x, y) -y,
#'   Sx = function(t, x) matrix(0.3, 1, 1),
#'   Sy = function(t, x, y) matrix(0.3, 1, 1),
#'   k = 1, l = 1)
#'
#' @export
stochastic_model <- function(f, g, Sx, Sy, k, l,
                             vectorized_members = FALSE, name = NULL,
                             meta = list()) {
  dims <- .state_dims(k, l)
  if (!all(vapply(list(f, g, Sx, Sy), is.function, logical(1))))
    aci_abort("aci_error_model_contract", "f, g, Sx, and Sy must be functions.")
  if (!is.list(meta))
    aci_abort("aci_error_model_contract", "meta must be a list.")
  if (!is.logical(vectorized_members) || length(vectorized_members) != 1L ||
      is.na(vectorized_members))
    aci_abort("aci_error_model_contract",
              "vectorized_members must be TRUE or FALSE.")
  m <- structure(list(f = f, g = g, Sx = Sx, Sy = Sy, k = dims[["k"]],
                      l = dims[["l"]],
                      vectorized_members = isTRUE(vectorized_members),
                      name = name %||% "stochastic_model", meta = meta),
                 class = "stochastic_model")
  validate_stochastic(m)
}


#' Validate a general stochastic model (internal)
#'
#' Probes the supplied drift and diffusion functions at several states and
#' times and checks that their dimensions match the declared `k` and `l`.
#'
#' @param m A `stochastic_model` object.
#' @param n_probe Number of probe evaluations used to check dimensions.
#' @returns The validated model `m`.
#' @noRd
validate_stochastic <- function(m, n_probe = 5) {
  dims <- .state_dims(m$k, m$l)
  m$k <- dims[["k"]]; m$l <- dims[["l"]]
  if (!all(vapply(m[c("f", "g", "Sx", "Sy")], is.function, logical(1))))
    aci_abort("aci_error_model_contract", "f, g, Sx, and Sy must be functions.")
  if (length(n_probe) != 1L || !is.numeric(n_probe) || !is.finite(n_probe) ||
      n_probe < 1L || n_probe != as.integer(n_probe))
    aci_abort("aci_error_model_contract", "n_probe must be a positive integer.")
  dx_width <- dy_width <- NULL
  for (i in seq_len(n_probe)) {
    # Deterministic probes validate constructors without changing the caller's
    # random-number stream (model creation must be referentially transparent).
    t <- (i - 1) * sqrt(2)
    x <- sin(seq_len(m$k) * (i + 0.25))
    y <- cos(seq_len(m$l) * (i + 0.5))
    fv <- m$f(t, x, y); gv <- m$g(t, x, y)
    if (!is.numeric(fv) || length(fv) != m$k)
      aci_abort("aci_error_model_contract", "f must return a numeric vector of length k.")
    if (!is.numeric(gv) || length(gv) != m$l)
      aci_abort("aci_error_model_contract", "g must return a numeric vector of length l.")
    Sx1 <- as.matrix(m$Sx(t, x)); Sy1 <- as.matrix(m$Sy(t, x, y))
    if (!is.numeric(Sx1) || !is.numeric(Sy1) ||
        any(!is.finite(c(fv, gv, Sx1, Sy1))))
      aci_abort("aci_error_model_contract", "Model coefficients must be finite numeric values.")
    if (nrow(Sx1) != m$k) aci_abort("aci_error_model_contract", "Sx must have k rows.")
    if (nrow(Sy1) != m$l) aci_abort("aci_error_model_contract", "Sy must have l rows.")
    if (is.null(dx_width)) {
      dx_width <- ncol(Sx1); dy_width <- ncol(Sy1)
    } else if (ncol(Sx1) != dx_width || ncol(Sy1) != dy_width) {
      aci_abort("aci_error_model_contract",
                "Diffusion channel counts must remain constant over state and time.")
    }
    if (!dx_width || !dy_width)
      aci_abort("aci_error_model_contract", "Sx and Sy must each expose at least one Wiener channel.")
    # SPEC-04 contract: Gx SPD; Sx has no y argument by signature already.
    Gx <- Sx1 %*% t(Sx1)
    rc <- rcond(Gx)
    if (!is.finite(rc) || rc < 1e-12)
      aci_abort("aci_error_gram",
        "Gx = Sx Sx' is degenerate; the engines require nondegenerate observation noise (jiang2026enkbs eq. 1b / andreou2026aci SI regularity).")
  }
  m
}


#' Print a stochastic model
#'
#' @param x A `stochastic_model` or `cgns_model` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.stochastic_model <- function(x, ...) {
  cat(sprintf("<%s> '%s': k = %d observed, l = %d hidden%s\n",
              class(x)[1], x$name, x$k, x$l,
              if (isTRUE(x$meta$correlated_noise)) " (correlated noise)" else ""))
  invisible(x)
}


################################################################################
# CGNS class: coefficients, Grams, validation, affine constructor
################################################################################

#' Conditional-Gaussian nonlinear system
#'
#' Constructs a conditional-Gaussian nonlinear system, in which the observed
#' drift is affine in the hidden state and the hidden drift is linear in it, so
#' that the conditional statistics of the hidden component are Gaussian in
#' closed form. The noise is described by up to two shared Wiener channels, so
#' that correlated observation and signal noise can be represented.
#'
#' @param Lx Function of `(t, x)` giving the coupling of the hidden state into
#'   the observed drift, a `k` by `l` matrix.
#' @param fx Function of `(t, x)` giving the remaining observed drift.
#' @param Ly Function of `(t, x)` giving the hidden self-drift, an `l` by `l`
#'   matrix.
#' @param fy Function of `(t, x)` giving the remaining hidden drift.
#' @param Sx1 Function of `(t, x)` giving the observed diffusion on the first
#'   Wiener channel.
#' @param Sx2 Optional function of `(t, x)` giving the observed diffusion on the
#'   second Wiener channel; `NULL` is a zero block matched to `Sy2`.
#' @param Sy1 Optional function of `(t, x)` giving the hidden diffusion on the
#'   first Wiener channel; `NULL` is a zero block matched to `Sx1`.
#' @param Sy2 Function of `(t, x)` giving the hidden diffusion on the second
#'   Wiener channel.
#' @param k Observed dimension; a positive whole number.
#' @param l Hidden dimension; a positive whole number.
#' @param name Optional 1-length character label for the model.
#' @param meta Optional named list of metadata carried on the object.
#' @returns An object of class `cgns_model`, which also inherits from
#'   `stochastic_model`.
#'
#' @seealso [cgns_from_affine()], [da_filter()], [aci()]
#'
#' @examples
#' cgns_model(
#'   Lx = function(t, x) matrix(1, 1, 1),
#'   fx = function(t, x) -0.5 * x,
#'   Ly = function(t, x) matrix(-0.5, 1, 1),
#'   fy = function(t, x) 0,
#'   Sx1 = function(t, x) matrix(0.5, 1, 1),
#'   Sy2 = function(t, x) matrix(1, 1, 1),
#'   k = 1, l = 1)
#'
#' @export
cgns_model <- function(Lx, fx, Ly, fy, Sx1, Sx2 = NULL, Sy1 = NULL, Sy2,
                       k, l, name = NULL, meta = list()) {
  dims <- .state_dims(k, l); k <- dims[["k"]]; l <- dims[["l"]]
  if (!all(vapply(list(Lx, fx, Ly, fy, Sx1, Sy2), is.function, logical(1))) ||
      (!is.null(Sx2) && !is.function(Sx2)) ||
      (!is.null(Sy1) && !is.function(Sy1)))
    aci_abort("aci_error_model_contract",
              "CGNS drifts and diffusion coefficients must be functions.")
  if (!is.list(meta))
    aci_abort("aci_error_model_contract", "meta must be a list.")
  # zero cross-channels must match the PARTNER channel width (shared Wiener
  # channels W1, W2): Sx2 pairs with Sy2, Sy1 pairs with Sx1.
  if (is.null(Sx2)) Sx2 <- function(t, x) matrix(0, k, ncol(as.matrix(Sy2(t, x))))
  if (is.null(Sy1)) Sy1 <- function(t, x) matrix(0, l, ncol(as.matrix(Sx1(t, x))))
  obj <- list(Lx = Lx, fx = fx, Ly = Ly, fy = fy,
              Sx1 = Sx1, Sx2 = Sx2, Sy1 = Sy1, Sy2 = Sy2,
              k = k, l = l, name = name %||% "cgns_model", meta = meta)
  # derived parent drifts so every CGNS runs on the ensemble engine (SPEC-05 s2.2)
  obj$f <- function(t, x, y) drop(fx(t, x) + as.matrix(Lx(t, x)) %*% y)
  obj$g <- function(t, x, y) drop(fy(t, x) + as.matrix(Ly(t, x)) %*% y)
  obj$Sx <- function(t, x) cbind(as.matrix(Sx1(t, x)), as.matrix(obj$Sx2(t, x)))
  obj$Sy <- function(t, x, y) cbind(as.matrix(obj$Sy1(t, x)), as.matrix(Sy2(t, x)))
  obj$vectorized_members <- FALSE
  m <- structure(obj, class = c("cgns_model", "stochastic_model"))
  validate_cgns(m)
}


#' Noise Gram products of a CGNS at one point (internal)
#'
#' @param m A `cgns_model` object.
#' @param t 1-length numeric time.
#' @param x Numeric vector of observations at time `t`.
#' @returns A list with the Gram matrices `gxx` (`k` by `k`), `gyy` (`l` by `l`)
#'   and `gyx` (`l` by `k`).
#' @noRd
cgns_grams <- function(m, t, x) {
  Sx1 <- as.matrix(m$Sx1(t, x)); Sx2 <- as.matrix(m$Sx2(t, x))
  Sy1 <- as.matrix(m$Sy1(t, x)); Sy2 <- as.matrix(m$Sy2(t, x))
  list(gxx = Sx1 %*% t(Sx1) + Sx2 %*% t(Sx2),
       gyy = Sy1 %*% t(Sy1) + Sy2 %*% t(Sy2),
       gyx = Sy1 %*% t(Sx1) + Sy2 %*% t(Sx2))
}


#' Scale-free test for non-zero noise cross-covariance (internal)
#'
#' A scale-free test for a non-zero instantaneous observation/signal noise
#' cross-covariance. Cauchy-Schwarz gives ||gyx|| <= sqrt(||gxx||||gyy||), so
#' comparing to that natural scale avoids changing the routing when a model's
#' physical units are rescaled.
#'
#' @param co Coefficient list from `eval_coefs()`.
#' @param tol Relative tolerance against the Cauchy-Schwarz scale.
#' @returns 1-length logical, `TRUE` when the cross-covariance is non-zero.
#' @noRd
.has_cross_noise <- function(co, tol = 100 * .Machine$double.eps) {
  cross <- max(abs(co$gyx))
  scale <- sqrt(max(abs(co$gxx)) * max(abs(co$gyy)))
  is.finite(cross) && is.finite(scale) &&
    cross > tol * max(scale, .Machine$double.xmin)
}


#' Evaluate all CGNS coefficients at one point (internal)
#'
#' @param m A `cgns_model` object.
#' @param t 1-length numeric time.
#' @param x Numeric vector of observations at time `t`.
#' @returns A list with `Lx`, `fx`, `Ly`, `fy` and the Gram matrices `gxx`,
#'   `gyy` and `gyx`.
#' @noRd
eval_coefs <- function(m, t, x) {
  g <- cgns_grams(m, t, x)
  list(Lx = as.matrix(m$Lx(t, x)), fx = as.numeric(m$fx(t, x)),
       Ly = as.matrix(m$Ly(t, x)), fy = as.numeric(m$fy(t, x)),
       gxx = g$gxx, gyy = g$gyy, gyx = g$gyx)
}


#' Validate a conditional-Gaussian model (internal)
#'
#' Probes the coefficient functions, checks that the observed drift is affine
#' and the hidden drift linear in the hidden state, and records whether the
#' realised coefficients carry a non-zero noise cross-covariance.
#'
#' @param m A `cgns_model` object.
#' @param n_probe Number of probe evaluations used to check the contract.
#' @returns The validated model `m`, with `meta$correlated_noise` set.
#' @noRd
validate_cgns <- function(m, n_probe = 5) {
  dims <- .state_dims(m$k, m$l)
  m$k <- dims[["k"]]; m$l <- dims[["l"]]
  if (length(n_probe) != 1L || !is.numeric(n_probe) || !is.finite(n_probe) ||
      n_probe < 1L || n_probe != as.integer(n_probe))
    aci_abort("aci_error_model_contract", "n_probe must be a positive integer.")
  corr <- FALSE
  widths <- NULL
  for (i in seq_len(n_probe)) {
    t <- (i - 1) * sqrt(2)
    x <- sin(seq_len(m$k) * (i + 0.25))
    co <- eval_coefs(m, t, x)
    if (!all(dim(co$Lx) == c(m$k, m$l)))
      aci_abort("aci_error_model_contract", "Lx(t,x) must be k x l.")
    if (!all(dim(co$Ly) == c(m$l, m$l)))
      aci_abort("aci_error_model_contract", "Ly(t,x) must be l x l.")
    if (length(co$fx) != m$k || length(co$fy) != m$l)
      aci_abort("aci_error_model_contract", "fx / fy dims wrong.")
    Sx1 <- as.matrix(m$Sx1(t, x)); Sx2 <- as.matrix(m$Sx2(t, x))
    Sy1 <- as.matrix(m$Sy1(t, x)); Sy2 <- as.matrix(m$Sy2(t, x))
    if (nrow(Sx1) != m$k || nrow(Sx2) != m$k ||
        nrow(Sy1) != m$l || nrow(Sy2) != m$l)
      aci_abort("aci_error_model_contract", "CGNS diffusion matrices have incorrect row dimensions.")
    if (ncol(Sx1) != ncol(Sy1) || ncol(Sx2) != ncol(Sy2))
      aci_abort("aci_error_model_contract",
                "Shared CGNS channel pairs Sx1/Sy1 and Sx2/Sy2 must have matching column counts.")
    if (any(!is.finite(c(co$Lx, co$fx, co$Ly, co$fy, Sx1, Sx2, Sy1, Sy2))))
      aci_abort("aci_error_model_contract", "Non-finite CGNS coefficient values at a probe point.")
    if (rcond(co$gxx) < 1e-12)
      aci_abort("aci_error_gram",
        "gxx degenerate at probe point (SPEC-01 s5.0).")
    now_widths <- c(ncol(Sx1), ncol(Sx2))
    if (is.null(widths)) widths <- now_widths else if (!identical(now_widths, widths))
      aci_abort("aci_error_model_contract",
                "CGNS shared-channel counts must remain constant over state and time.")
    # Affinity probe on the model drifts (defense-in-depth, SPEC-05 s2.2).
    ya <- sin(seq_len(m$l) * (i + 0.7))
    yc <- cos(seq_len(m$l) * (i + 1.3)); yb <- (ya + yc) / 2
    sc <- max(1, sum(abs(m$f(t, x, ya))) + sum(abs(m$f(t, x, yc))))
    if (sum(abs(m$f(t, x, ya) + m$f(t, x, yc) - 2 * m$f(t, x, yb))) > 1e-8 * sc)
      aci_abort("aci_error_model_contract", "Observed drift is not affine in the hidden state.")
    sg <- max(1, sum(abs(m$g(t, x, ya))) + sum(abs(m$g(t, x, yc))))
    if (sum(abs(m$g(t, x, ya) + m$g(t, x, yc) - 2 * m$g(t, x, yb))) > 1e-8 * sg)
      aci_abort("aci_error_model_contract", "Hidden drift is not affine in the hidden state.")
    if (.has_cross_noise(co)) corr <- TRUE
  }
  m$meta$correlated_noise <- corr
  m
}


#' Conditionally linear model from constant coefficients
#'
#' Convenience constructor for a conditional-Gaussian system whose coupling and
#' self-drift are constant, supplied either as scalars or as full coefficient
#' matrices.
#'
#' @param lambda_x Coupling of the hidden state into the observed drift; a
#'   finite scalar or a `k` by `l` matrix.
#' @param lambda_y Hidden self-drift; a finite scalar or an `l` by `l` matrix.
#' @param fx Function of `(t, x)` giving the remaining observed drift.
#' @param fy Function of `(t, x)` giving the remaining hidden drift.
#' @param sigma_x Observed noise amplitude.
#' @param sigma_y Hidden noise amplitude.
#' @param k Observed dimension; a positive whole number.
#' @param l Hidden dimension; a positive whole number.
#' @returns An object of class `cgns_model`.
#'
#' @seealso [cgns_model()]
#'
#' @examples
#' conditionally_linear_model(
#'   lambda_x = 1, lambda_y = -0.5,
#'   fx = function(t, x) -0.5 * x, fy = 0,
#'   sigma_x = 0.5, sigma_y = 1)
#'
#' @export
conditionally_linear_model <- function(lambda_x, lambda_y, fx, fy,
                                       sigma_x, sigma_y, k = 1, l = 1) {
  dims <- .state_dims(k, l); k <- dims[["k"]]; l <- dims[["l"]]
  valid_numeric <- function(z, lengths)
    is.numeric(z) && length(z) %in% lengths && all(is.finite(z))
  if (!valid_numeric(lambda_x, c(1L, k * l)) ||
      !valid_numeric(lambda_y, c(1L, l * l)))
    aci_abort("aci_error_model_contract",
              "lambda_x and lambda_y must be finite scalars or full coefficient matrices.")
  if (!is.function(fx) && !valid_numeric(fx, c(1L, k)))
    aci_abort("aci_error_model_contract", "fx must be a function or a finite scalar/vector of length k.")
  if (!is.function(fy) && !valid_numeric(fy, c(1L, l)))
    aci_abort("aci_error_model_contract", "fy must be a function or a finite scalar/vector of length l.")
  if (!valid_numeric(sigma_x, c(1L, k)) || any(sigma_x <= 0) ||
      !valid_numeric(sigma_y, c(1L, l)) || any(sigma_y < 0))
    aci_abort("aci_error_model_contract",
              "sigma_x must be positive and sigma_y non-negative, each scalar or state-sized.")
  expand_linear <- function(z, nr, nc, label) {
    if (length(z) == 1L) {
      if (nr != nc && nr * nc > 1L)
        aci_abort("aci_error_model_contract", sprintf(
          "Scalar %s is ambiguous for a %d x %d map; supply all coefficients.",
          label, nr, nc))
      if (nr == nc) return(diag(as.numeric(z), nr))
    }
    matrix(z, nr, nc)
  }
  Lx0 <- expand_linear(lambda_x, k, l, "lambda_x")
  Ly0 <- expand_linear(lambda_y, l, l, "lambda_y")
  m <- cgns_model(
    Lx = function(t, x) Lx0,
    fx = if (is.function(fx)) fx else function(t, x) rep_len(fx, k),
    Ly = function(t, x) Ly0,
    fy = if (is.function(fy)) fy else function(t, x) rep_len(fy, l),
    Sx1 = function(t, x) diag(rep_len(sigma_x, k), k),
    Sy2 = function(t, x) diag(rep_len(sigma_y, l), l),
    k = k, l = l, name = "conditionally_linear")
  class(m) <- c("conditionally_linear_model", class(m))
  m
}


#' Build cgns_model  from a general affine-in-hidden drift pair by exact
#'
#' Recovers the conditional-Gaussian coefficients from drift functions written
#' in terms of the full state, by exact affine differencing. The supplied
#' functions are checked before they are replaced by the reconstructed
#' representation, so a drift that is not affine in the hidden state is
#' rejected rather than silently linearised.
#'
#' @param f_full Function of `(t, x, y)` giving the full observed drift.
#' @param g_full Function of `(t, x, y)` giving the full hidden drift.
#' @param Sx Function of `(t, x)` giving the observed diffusion on the first
#'   Wiener channel.
#' @param Sy_hidden Function of `(t, x)` giving the hidden diffusion on its own
#'   Wiener channel.
#' @param k Observed dimension; a positive whole number.
#' @param l Hidden dimension; a positive whole number.
#' @param name Optional 1-length character label for the model.
#' @param Sx2 Optional function of `(t, x)` giving the observed diffusion on the
#'   second Wiener channel.
#' @param Sy_shared Optional function of `(t, x)` giving the hidden diffusion on
#'   the shared Wiener channel, which introduces correlated noise.
#' @param meta Optional named list of metadata carried on the object.
#' @returns An object of class `cgns_model`.
#'
#' @seealso [cgns_model()]
#'
#' @examples
#' cgns_from_affine(
#'   f_full = function(t, x, y) -0.5 * x + y,
#'   g_full = function(t, x, y) -0.5 * y,
#'   Sx = function(t, x) matrix(0.5, 1, 1),
#'   Sy_hidden = function(t, x) matrix(1, 1, 1),
#'   k = 1, l = 1)
#'
#' @export
cgns_from_affine <- function(f_full, g_full, Sx, Sy_hidden, k, l, name = NULL,
                             Sx2 = NULL, Sy_shared = NULL, meta = list()) {
  dims <- .state_dims(k, l); k <- dims[["k"]]; l <- dims[["l"]]
  if (!all(vapply(list(f_full, g_full, Sx, Sy_hidden), is.function, logical(1))) ||
      (!is.null(Sx2) && !is.function(Sx2)) ||
      (!is.null(Sy_shared) && !is.function(Sy_shared)))
    aci_abort("aci_error_model_contract",
              "Affine drifts and diffusion coefficients must be functions.")
  # Check the functions supplied by the caller, before replacing them by the
  # finite-difference affine representation.  Validating only the reconstructed
  # functions would silently linearise nonlinear inputs such as y^2.
  for (i in seq_len(5L)) {
    t <- (i - 1) * sqrt(2)
    x <- sin(seq_len(k) * (i + 0.25))
    ya <- sin(seq_len(l) * (i + 0.7))
    yc <- cos(seq_len(l) * (i + 1.3)); yb <- (ya + yc) / 2
    for (item in list(list(fun = f_full, n = k, label = "Observed"),
                      list(fun = g_full, n = l, label = "Hidden"))) {
      fa <- as.numeric(item$fun(t, x, ya)); fb <- as.numeric(item$fun(t, x, yb))
      fc <- as.numeric(item$fun(t, x, yc))
      if (length(fa) != item$n || length(fb) != item$n || length(fc) != item$n)
        aci_abort("aci_error_model_contract",
                  sprintf("%s drift has the wrong dimension.", item$label))
      scale <- max(1, sum(abs(fa)) + sum(abs(fc)))
      if (sum(abs(fa + fc - 2 * fb)) > 1e-8 * scale)
        aci_abort("aci_error_model_contract",
                  sprintf("%s drift is not affine in the hidden state.", item$label))
    }
  }
  Lx <- function(t, x) {
    base <- f_full(t, x, rep(0, l))
    vapply(seq_len(l), function(i) f_full(t, x, replace(rep(0, l), i, 1)) - base,
           numeric(k))
  }
  fx <- function(t, x) f_full(t, x, rep(0, l))
  Ly <- function(t, x) {
    base <- g_full(t, x, rep(0, l))
    vapply(seq_len(l), function(i) g_full(t, x, replace(rep(0, l), i, 1)) - base,
           numeric(l))
  }
  fy <- function(t, x) g_full(t, x, rep(0, l))
  LxM <- function(t, x) matrix(Lx(t, x), k, l)
  LyM <- function(t, x) matrix(Ly(t, x), l, l)
  cgns_model(Lx = LxM, fx = fx, Ly = LyM, fy = fy,
             Sx1 = Sx, Sx2 = Sx2, Sy1 = Sy_shared, Sy2 = Sy_hidden,
             k = k, l = l, name = name, meta = meta)
}


################################################################################
# Euler-Maruyama simulation and reusable noise stores
################################################################################

#' Simulate a stochastic or conditional-Gaussian model
#'
#' Generates truth-twin realisations by Euler-Maruyama integration, optionally
#' retaining the hidden path and the driving Wiener increments so that a later
#' ensemble smoother can reuse them.
#'
#' @param object A `stochastic_model` or `cgns_model` object.
#' @param nsim Positive whole number of realisations to generate.
#' @param seed Optional non-negative whole number seeding the generator.
#' @param T Positive 1-length numeric total simulated time, excluding burn-in.
#' @param dt Positive 1-length numeric integration step.
#' @param ic Optional list with elements `x0` and `y0` giving the initial state;
#'   `NULL` uses the model's default initial condition.
#' @param burn_in Non-negative 1-length numeric time discarded before recording.
#' @param keep_hidden `TRUE` to retain the hidden path.
#' @param keep_noise `TRUE` to retain the driving Wiener increments.
#' @param ... Must be empty; unused arguments are an error.
#' @returns An object of class `aci_sim` when `nsim` is one, otherwise a list of
#'   such objects.
#'
#' @seealso [stochastic_model()], [enkbf()]
#' @export
simulate.stochastic_model <- function(object, nsim = 1, seed = NULL,
                                      T, dt, ic = NULL, burn_in = 0,
                                      keep_hidden = TRUE, keep_noise = TRUE, ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to simulate().")
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim))
    aci_abort("aci_error_dims", "nsim must be a positive integer.")
  nsim <- as.integer(nsim)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)))
    aci_abort("aci_error_dims", "seed must be NULL or one non-negative integer.")
  if (length(keep_hidden) != 1L || is.na(keep_hidden) || !is.logical(keep_hidden) ||
      length(keep_noise) != 1L || is.na(keep_noise) || !is.logical(keep_noise))
    aci_abort("aci_error_dims", "keep_hidden and keep_noise must be TRUE or FALSE.")
  if (!is.null(seed)) set.seed(seed)
  one <- function() .simulate_one(object, T = T, dt = dt, ic = ic,
                                  burn_in = burn_in, keep_hidden = keep_hidden,
                                  keep_noise = keep_noise)
  if (nsim == 1) one() else lapply(seq_len(nsim), function(i) one())
}


#' @exportS3Method stats::simulate
simulate.cgns_model <- simulate.stochastic_model


#' Generate one Euler-Maruyama realisation (internal)
#'
#' @param m A `stochastic_model` or `cgns_model` object.
#' @param T Positive 1-length numeric total simulated time.
#' @param dt Positive 1-length numeric integration step.
#' @param ic Optional list with elements `x0` and `y0`.
#' @param burn_in Non-negative 1-length numeric time discarded before recording.
#' @param keep_hidden `TRUE` to retain the hidden path.
#' @param keep_noise `TRUE` to retain the driving Wiener increments.
#' @returns An object of class `aci_sim`.
#' @noRd
.simulate_one <- function(m, T, dt, ic, burn_in, keep_hidden, keep_noise) {
  if (length(T) != 1L || length(dt) != 1L || !is.finite(T) || !is.finite(dt) ||
      T <= 0 || dt <= 0)
    aci_abort("aci_error_dims", "T and dt must be finite positive scalars.")
  if (length(burn_in) != 1L || !is.finite(burn_in) || burn_in < 0)
    aci_abort("aci_error_dims", "burn_in must be a finite non-negative scalar.")
  N <- round(T / dt); Nb <- round(burn_in / dt)
  if (N < 1L || abs(N * dt - T) > 1e-8 * max(T, dt))
    aci_abort("aci_error_dims", "T must be an integer multiple of dt.")
  if (abs(Nb * dt - burn_in) > 1e-8 * max(burn_in, dt))
    aci_abort("aci_error_dims", "burn_in must be an integer multiple of dt.")
  N <- as.integer(N); Nb <- as.integer(Nb); Ntot <- Nb + N
  ic <- ic %||% m$meta$ic_default
  if (is.null(ic)) ic <- list()
  if (!is.list(ic)) aci_abort("aci_error_dims", "ic must be a list with x0 and y0.")
  x <- if (!is.null(ic$x0)) as.numeric(ic$x0) else rep(0, m$k)
  y <- if (!is.null(ic$y0)) as.numeric(ic$y0) else rep(0, m$l)
  if (length(x) != m$k || length(y) != m$l || any(!is.finite(c(x, y))))
    aci_abort("aci_error_dims", "Initial states x0 and y0 have incompatible dimensions or values.")
  shared <- inherits(m, "cgns_model")   # cgns shares channels W1 (dx_n1), W2
  X <- matrix(NA_real_, Ntot + 1, m$k); Y <- matrix(NA_real_, Ntot + 1, m$l)
  X[1, ] <- x; Y[1, ] <- y
  if (shared) {
    d1 <- ncol(as.matrix(m$Sx1(0, x)))
    d2 <- ncol(as.matrix(m$Sx2(0, x)))
    W <- array(stats::rnorm(d1 * Ntot), c(d1, Ntot))
    B <- array(stats::rnorm(d2 * Ntot), c(d2, Ntot))
  } else {
    dx_n <- ncol(as.matrix(m$Sx(0, x)))
    dy_n <- ncol(as.matrix(m$Sy(0, x, y)))
    W <- array(stats::rnorm(dx_n * Ntot), c(dx_n, Ntot))
    B <- array(stats::rnorm(dy_n * Ntot), c(dy_n, Ntot))
  }
  sq <- sqrt(dt)
  for (j in seq_len(Ntot)) {
    # Run burn-in on [-burn_in, 0), so the retained trajectory's t = 0 clock
    # agrees with every time-dependent coefficient used to generate it.
    t <- (j - 1 - Nb) * dt
    if (shared) {
      Sx1 <- as.matrix(m$Sx1(t, x)); Sx2 <- as.matrix(m$Sx2(t, x))
      Sy1 <- as.matrix(m$Sy1(t, x)); Sy2 <- as.matrix(m$Sy2(t, x))
      w1 <- W[, j]; w2 <- B[, j]
      xn <- x + m$f(t, x, y) * dt + sq * drop(Sx1 %*% w1 + Sx2 %*% w2)
      yn <- y + m$g(t, x, y) * dt + sq * drop(Sy1 %*% w1 + Sy2 %*% w2)
    } else {
      xn <- x + m$f(t, x, y) * dt + sq * drop(as.matrix(m$Sx(t, x)) %*% W[, j])
      yn <- y + m$g(t, x, y) * dt + sq * drop(as.matrix(m$Sy(t, x, y)) %*% B[, j])
    }
    if (any(!is.finite(c(xn, yn))))
      aci_abort("aci_error_sim_divergence",
                sprintf("Simulation diverged at step %d (t = %.4f).", j, t),
                step = j, last_state = list(x = x, y = y))
    x <- xn; y <- yn; X[j + 1, ] <- x; Y[j + 1, ] <- y
  }
  keep <- (Nb + 1):(Ntot + 1)
  obs <- observed_trajectory(dt * (seq_along(keep) - 1), X[keep, , drop = FALSE])
  structure(list(
    obs = obs,
    hidden = if (keep_hidden) Y[keep, , drop = FALSE] else NULL,
    noise = if (keep_noise) structure(list(
      B = B[, (Nb + 1):Ntot, drop = FALSE], W = W[, (Nb + 1):Ntot, drop = FALSE]),
      class = "noise_store") else NULL,
    model = m, meta = list(perfect_model = TRUE, scheme = "euler_maruyama",
                           dt = dt, burn_in = burn_in,
                           burn_interval = c(-burn_in, 0))),
    class = "aci_sim")
}
