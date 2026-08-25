################################################################################
## discovery.R - model discovery: entropy, sampling, constrained fits
## ########################################################################## ##
##
## Contents:
##   * polynomial CGNS dictionaries:
##       - cgns_library, print.function_library, eval_library
##
##   * causation entropy and structure thresholding:
##       - causation_entropy, threshold_structure
##
##   * exact FFBS posterior sampling for CGNS:
##       - sample_paths
##
##   * constrained MLE, cross-equation energy constraints, learn_model:
##       - constrained_mle, model_from_learned, learn_model, print.learned_model,
##         .refit_energy_joint
##
################################################################################


################################################################################
# polynomial CGNS dictionaries
################################################################################

#' Polynomial dictionary over the full state z = (x, y) with CGNS tags.
#'
#' @param k Number of observed variables.
#' @param l Number of hidden variables.
#' @param degree Maximum polynomial degree.
#' @param intercept Whether to include the constant term.
#'
#' @examples
#' cgns_library(k = 1, l = 1, degree = 2)
#'
#' @export
cgns_library <- function(k, l, degree = 2, intercept = TRUE) {
  dims <- .state_dims(k, l); k <- dims[["k"]]; l <- dims[["l"]]
  if (length(intercept) != 1L || is.na(intercept) || !is.logical(intercept))
    aci_abort("aci_error_dims", "intercept must be TRUE or FALSE.")
  n <- k + l
  if (!is.numeric(degree) || length(degree) != 1L || !is.finite(degree) || degree < 0 ||
      degree != as.integer(degree))
    aci_abort("aci_error_dims", "degree must be a non-negative integer.")
  degree <- as.integer(degree)
  nm_z <- c(paste0("x", seq_len(k)), if (l) paste0("y", seq_len(l)))
  terms <- list(); names_ <- character(); ydeg <- integer()
  if (intercept) { terms[["1"]] <- function(Z, t) rep(1, nrow(Z))
    names_ <- "1"; ydeg <- 0L }
  exponent_vectors <- function(nvar, total) {
    rec <- function(pos, left, prefix) {
      if (pos == nvar) return(list(c(prefix, left)))
      unlist(lapply(left:0, function(a) rec(pos + 1L, left - a,
                                            c(prefix, a))), recursive = FALSE)
    }
    rec(1L, total, integer())
  }
  if (degree > 0L) for (d in seq_len(degree)) for (e in exponent_vectors(n, d)) {
    pieces <- vapply(which(e > 0L), function(i) {
      if (e[i] == 1L) nm_z[i] else paste0(nm_z[i], "^", e[i])
    }, character(1))
    nm <- paste(pieces, collapse = "*")
    f <- local({ ee <- e; function(Z, t) {
      ans <- rep(1, nrow(Z))
      for (ii in which(ee > 0L)) ans <- ans * Z[, ii]^ee[ii]
      ans
    }})
    terms[[nm]] <- f; names_ <- c(names_, nm)
    ydeg <- c(ydeg, if (l > 0L) sum(e[(k + 1L):n]) else 0L)
  }
  if (!length(terms))
    aci_abort("aci_error_dims", "The requested function library contains no terms.")
  structure(list(funs = terms, names = names_, y_degree = ydeg, k = k, l = l),
            class = "function_library")
}


#' Print a function library
#'
#' @param x A `function_library` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.function_library <- function(x, ...) {
  cat(sprintf("<function_library> %d terms over k=%d observed + l=%d hidden\n",
              length(x$funs), x$k, x$l)); invisible(x)
}


#' Evaluate a function library
#'
#' Evaluates every candidate term of a library at the supplied states and
#' times.
#'
#' @param lib A `function_library` object; see [cgns_library()].
#' @param Z Numeric matrix of states, one row per time.
#' @param t Numeric vector of times, one per row of `Z`.
#' @returns A numeric matrix with one column per library term.
#'
#' @seealso [cgns_library()], [learn_model()]
#'
#' @examples
#' lib <- cgns_library(k = 1, l = 1, degree = 2)
#' head(eval_library(lib, Z = rbind(c(0.5, 0.2), c(0.6, 0.1)), t = c(0, 0.1)))
#'
#' @export
eval_library <- function(lib, Z, t) {
  if (!inherits(lib, "function_library"))
    aci_abort("aci_error_dims", "lib must be a function_library.")
  Z <- as.matrix(Z)
  if (!is.numeric(Z) || ncol(Z) != lib$k + lib$l || !nrow(Z) ||
      any(!is.finite(Z)))
    aci_abort("aci_error_dims", "Z must be a finite matrix with k + l columns.")
  if (length(t) != 1L && length(t) != nrow(Z))
    aci_abort("aci_error_dims", "t must be scalar or have one value per row of Z.")
  if (!is.numeric(t) || any(!is.finite(t)))
    aci_abort("aci_error_dims", "t must contain finite numeric values.")
  M <- vapply(lib$funs, function(f) {
    ans <- f(Z, t)
    if (!is.numeric(ans) || length(ans) != nrow(Z) || any(!is.finite(ans)))
      aci_abort("aci_error_model_contract",
                "Every library function must return one finite numeric value per row.")
    ans
  }, numeric(nrow(Z)))
  # vapply() returns a bare vector when Z has a single row, which colnames()
  # cannot take; the reshape keeps the result a matrix for every record length.
  M <- matrix(M, nrow = nrow(Z))
  colnames(M) <- lib$names
  M
}


################################################################################
# causation entropy and structure thresholding
################################################################################

#' Gaussian causation entropy
#'
#' CE of each candidate column into `target`, conditioned on the remaining
#' candidates (jiang2026enkbs's Gaussian CE): 0.5 log( RSS(rest) / RSS(rest +
#' candidate) ).
#'
#' Convention: the target here is the increment, following the derivative
#' formulation of the framework jiang2026enkbs cites. The published EnKBS
#' discovery script instead conditions on the next state via covariance
#' log-determinants, which inflates each equation's own state term and yields
#' different absolute CE values, so thresholds do not transfer between the two
#' conventions.
#' @param Theta Numeric design matrix.
#' @param target Numeric response vector.
#' @param candidates Candidate column names in `Theta`.
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#'
#' Sun, J. and Bollt, E. M. (2014). Causation entropy identifies indirect
#' influences, dominance of neighbors and anticipatory couplings. *Physica D*
#' **267**, 49-57. \doi{10.1016/j.physd.2013.07.001}
#'
#' Elinger, J. and Rogers, J. (2021). Causation entropy method for covariate
#' selection in dynamic models. *2021 American Control Conference (ACC)*,
#' 2842-2847. \doi{10.23919/ACC50511.2021.9483371}
#' @examples
#' set.seed(1)
#' Th <- cbind(a = rnorm(50), b = rnorm(50))
#' causation_entropy(Th, target = Th[, "a"] * 2 + rnorm(50, sd = 0.1))
#'
#' @export
causation_entropy <- function(Theta, target, candidates = colnames(Theta)) {
  Theta <- as.matrix(Theta); target <- as.numeric(target)
  if (!is.numeric(Theta) || !nrow(Theta) || nrow(Theta) != length(target) ||
      any(!is.finite(c(Theta, target))) || is.null(colnames(Theta)) ||
      anyDuplicated(colnames(Theta)) || !is.character(candidates) ||
      !length(candidates) || anyNA(candidates) || anyDuplicated(candidates) ||
      any(!candidates %in% colnames(Theta)))
    aci_abort("aci_error_dims", paste(
      "Theta must be a finite, uniquely named design matrix aligned with",
      "target, and candidates must name its columns."))
  # No implicit intercept: the dictionary's own constant term (if any) is a
  # first-class candidate, so conditioning sets are exactly `cols`.
  fitR <- function(cols) {
    if (!length(cols)) return(sum(target^2))
    q <- qr(Theta[, cols, drop = FALSE]); sum(qr.resid(q, target)^2)
  }
  # `candidates` defines both the tested set and its conditioning set.  Public
  # callers may pass a wider design matrix containing nuisance columns; those
  # must not silently change the requested candidate-only CE calculation.
  all_ <- candidates
  vapply(candidates, function(cn) {
    rest <- setdiff(all_, cn)
    r0 <- fitR(rest); r1 <- fitR(all_)
    max(0.5 * log(max(r0, 1e-300) / max(r1, 1e-300)), 0)
  }, numeric(1))
}


#' Apply structure threshold to causation entropies
#'
#' Selects the candidate terms whose causation entropy exceeds a threshold,
#' retaining at least a minimum number of them.
#'
#' @param ce Named numeric vector of causation entropies.
#' @param threshold Either `"auto"` or a numeric threshold. `"auto"` is a
#'   relative floor, `1e-2` of the largest CE; the jiang2026enkbs experiments
#'   use an absolute threshold instead (`1e-3` in the paper), which can be
#'   reproduced by passing that number directly.
#' @param value Optional numeric threshold used when `threshold` is `"auto"`
#'   and an explicit value is preferred.
#' @param keep_min Minimum number of terms retained.
#' @returns Character vector of the retained term names.
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#'
#' @seealso [causation_entropy()], [eval_library()]
#'
#' @examples
#' threshold_structure(c(a = 0.9, b = 0.001))
#'
#' @export
threshold_structure <- function(ce, threshold = "auto",
                                value = NULL, keep_min = 1) {
  if (!is.numeric(ce) || is.null(names(ce)) || anyNA(names(ce)) ||
      any(!nzchar(names(ce))) || anyDuplicated(names(ce)) || !length(ce) ||
      any(!is.finite(ce)) || any(ce < 0) ||
      !is.numeric(keep_min) || length(keep_min) != 1L ||
      !is.finite(keep_min) || keep_min < 0 || keep_min != as.integer(keep_min) ||
      keep_min > length(ce))
    aci_abort("aci_error_dims",
              "ce must be a finite non-negative named vector and keep_min a valid integer.")
  if (is.numeric(threshold)) {
    if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0)
      aci_abort("aci_error_dims", "A numeric threshold must be one finite non-negative value.")
    mode <- "value"; value <- threshold
  } else {
    if (!is.character(threshold) || length(threshold) != 1L ||
        is.na(threshold) || !threshold %in% c("auto", "gap", "value"))
      aci_abort("aci_error_dims", "threshold must be 'auto', 'gap', 'value', or a non-negative number.")
    mode <- threshold
  }
  if (mode == "auto") {
    thr <- 1e-2 * max(ce)                  # robust relative floor (default)
  } else if (mode == "gap") {
    srt <- sort(ce[ce > max(ce) * 1e-8], decreasing = TRUE)
    if (length(srt) > 1L) {
      gap <- which.max(-diff(log(srt)))
      thr <- sqrt(srt[gap] * srt[gap + 1L])
    } else {
      thr <- if (length(srt)) srt / 2 else 0
    }
  } else {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0)
      aci_abort("aci_error_dims", "value must be supplied as one finite non-negative threshold.")
    thr <- value
  }
  kept <- names(ce)[ce > thr]
  if (length(kept) < keep_min) kept <- names(sort(ce, decreasing = TRUE))[seq_len(keep_min)]
  structure(kept, threshold = thr, ce = ce)
}


################################################################################
# Euler-consistent FFBS posterior sampling for CGNS
################################################################################

#' Sample hidden trajectories conditioned on an observed path
#'
#' Method `"ffbs"` uses the one-step Euler CGNS transition conditional on each
#' observed increment. Method `"enkbs"` returns EnKBS smoother members for a
#' compatible stochastic model.
#'
#' The FFBS pass conditions on a filter it computes itself when `filter` is not
#' supplied. The transition below is the single-step Euler CGNS transition, so
#' FFBS requires the matching explicit, one-step filter. It does not silently mix
#' this transition with an implicit or sub-stepped filter.
#'
#' @param model Model used for assimilation.
#' @param obs Observed trajectory.
#' @param n_samples Number of hidden paths.
#' @param method FFBS or EnKBS sampler.
#' @param seed Reproducibility seed.
#' @param filter Optional compatible precomputed filter.
#' @param nontarget Optional conditional ACI masking specification.
#' @param init Optional Gaussian initial distribution. For FFBS this is passed
#'   to the closed-form filter; for EnKBS it initializes the ensemble sampler.
#' @param stepper Filter covariance stepper used by FFBS.
#' @param nsub Number of filter substeps.
#' @param ... Additional filter or ensemble arguments.
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' sample_paths(m, ob, n_samples = 3, seed = 1)
#'
#' @export
sample_paths <- function(model, obs, n_samples = 20,
                         method = c("ffbs", "enkbs"), seed = NULL,
                         filter = NULL, nontarget = NULL,
                         init = NULL, stepper = c("explicit", "implicit"),
                         nsub = 1L, ...) {
  method <- match.arg(method); stepper <- match.arg(stepper)
  dots <- list(...)
  obs <- as_obs(obs)
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples < 1L ||
      n_samples != as.integer(n_samples))
    aci_abort("aci_error_dims", "n_samples must be a positive integer.")
  n_samples <- as.integer(n_samples)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)))
    aci_abort("aci_error_dims", "seed must be NULL or one non-negative integer.")
  if (!is.numeric(nsub) || length(nsub) != 1L || !is.finite(nsub) ||
      nsub < 1L || nsub != floor(nsub))
    aci_abort("aci_error_dims", "nsub must be a positive integer.")
  nsub <- as.integer(nsub)
  if (!is.null(seed)) set.seed(seed)
  if (method == "enkbs") {
    sm <- if (is.null(filter)) {
      if (!is.null(init) && "ic_sampler" %in% names(dots))
        aci_abort("aci_error_dims",
                  "Supply only one of init and ic_sampler for ensemble path sampling.")
      args <- c(list(model = model, obs = obs, m = n_samples, seed = seed,
                     nontarget = nontarget), dots)
      if (!is.null(init)) args$ic_sampler <- init
      fr <- do.call(enkbf, args)
      enkbs(model, fr$path, fr$noise)
    } else {
      if (!is.null(init))
        aci_abort("aci_error_dims",
                  "init cannot be combined with a supplied ensemble filter.")
      do.call(da_smooth.stochastic_model,
              c(list(model = model, obs = obs, filter = filter,
                     nontarget = nontarget), dots))
    }
    return(sm$members)
  }
  if (!inherits(model, "cgns_model"))
    aci_abort("aci_error_model_contract", "method='ffbs' requires a cgns_model.")
  if (length(dots))
    aci_abort("aci_error_dims", "Unused arguments were supplied to FFBS path sampling.")
  source_model <- model
  bundle <- .compile_cgns_run(model, obs, nontarget)
  m <- bundle$model; ob <- bundle$obs
  if (stepper != "explicit" || nsub != 1L)
    aci_abort("aci_error_stepper",
              "FFBS requires stepper = 'explicit' and nsub = 1 to match its Euler transition.")
  if (!is.null(filter) &&
      (!identical(filter$meta$stepper %||% "explicit", "explicit") ||
       (filter$meta$nsub %||% 1L) != 1L))
    aci_abort("aci_error_stepper", "The supplied FFBS filter is not explicit single-step.")
  if (!is.null(filter))
    .validate_gaussian_path(filter, ob, m$l, "filter", bundle$nontarget, model = m,
                            source_model = source_model)
  if (!is.null(filter) && !is.null(init) &&
      !.same_gaussian_init(init, filter$meta$init, m$l))
    aci_abort("aci_error_dims", "init conflicts with the prior stored on the supplied filter.")
  filt <- filter %||% .cgns_filter_compiled(
    bundle, init = init, stepper = stepper, nsub = nsub, validate = FALSE
  )
  N1 <- length(ob$t); N <- N1 - 1L; dt <- ob$dt; l <- m$l
  Y <- array(NA_real_, c(l, N1, n_samples))
  chN <- safe_chol(filt$cov[, , N1])
  Y[, N1, ] <- filt$mean[N1, ] + t(chN) %*% matrix(stats::rnorm(l * n_samples), l)
  for (j in N:1) {
    co <- .compiled_co(bundle, j)
    Gi <- .compiled_ginv(bundle, j)
    A  <- diag(l) + (co$Ly - co$gyx %*% Gi %*% co$Lx) * dt
    cc <- co$fy * dt + drop(co$gyx %*% Gi %*%
            (ob$x[j + 1, ] - ob$x[j, ] - co$fx * dt))
    Qc <- sym(co$gyy - co$gyx %*% Gi %*% t(co$gyx)) * dt
    Rf <- filt$cov[, , j]; muf <- filt$mean[j, ]
    Rfi <- chol_solve(Rf, diag(l), "Rf")
    Qci <- chol_solve(spd_floor(Qc), diag(l), "Qc")
    Sig <- chol_solve(sym(Rfi + t(A) %*% Qci %*% A), diag(l), "FFBS")
    ch  <- safe_chol(spd_floor(Sig))
    base <- drop(Rfi %*% muf)
    Ynext <- matrix(Y[, j + 1, ], l)
    Mn <- Sig %*% (base + t(A) %*% Qci %*% (Ynext - cc))
    Y[, j, ] <- Mn + t(ch) %*% matrix(stats::rnorm(l * n_samples), l)
  }
  Y
}


################################################################################
# constrained MLE, cross-equation energy constraints, learn_model
################################################################################

#' constrained least squares
#'
#' Minimize ||dZ - Theta xi||^2 (+ ridge) subject to A xi = b (KKT system).
#'
#' @param Theta Numeric design matrix.
#' @param dZ Increment-rate response.
#' @param constraints Optional list containing equality matrix `A` and vector `b`.
#' @param ridge Non-negative ridge penalty.
#'
#' @examples
#' set.seed(1)
#' Th <- cbind(a = rnorm(50), b = rnorm(50))
#' constrained_mle(Th, dZ = Th[, "a"] * 2 + rnorm(50, sd = 0.1))
#'
#' @export
constrained_mle <- function(Theta, dZ, constraints = NULL, ridge = 0) {
  Theta <- as.matrix(Theta); dZ <- as.numeric(dZ)
  if (!is.numeric(Theta) || nrow(Theta) < 2L || ncol(Theta) < 1L ||
      nrow(Theta) != length(dZ) || any(!is.finite(c(Theta, dZ))) ||
      is.null(colnames(Theta)) || any(!nzchar(colnames(Theta))) ||
      anyDuplicated(colnames(Theta)))
    aci_abort("aci_error_dims",
              "Theta must be a finite, uniquely named design matrix with at least two rows and aligned dZ.")
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0)
    aci_abort("aci_error_dims", "ridge must be one finite non-negative number.")
  p <- ncol(Theta)
  G <- crossprod(Theta) + diag(ridge, p)
  h <- crossprod(Theta, dZ)
  if (is.null(constraints)) {
    xi <- drop(chol_solve(spd_floor(G), h, "normal equations"))
  } else {
    if (!is.list(constraints) || is.null(constraints$A) || is.null(constraints$b))
      aci_abort("aci_error_dims", "constraints must be a list containing A and b.")
    A <- as.matrix(constraints$A); b <- as.numeric(constraints$b)
    if (!is.numeric(A) || nrow(A) < 1L || ncol(A) != p ||
        length(b) != nrow(A) || any(!is.finite(c(A, b))))
      aci_abort("aci_error_dims", "Constraint A and b dimensions or values are invalid.")
    q <- nrow(A)
    K <- rbind(cbind(G, t(A)), cbind(A, matrix(0, q, q)))
    sol <- tryCatch(solve(K, c(h, b)), error = function(e)
      aci_abort("aci_error_model_contract",
                sprintf("Constrained normal equations are singular: %s", conditionMessage(e))))
    xi <- sol[seq_len(p)]
  }
  names(xi) <- colnames(Theta)
  res <- dZ - drop(Theta %*% xi)
  structure(list(coef = xi, rss = sum(res^2), n = length(dZ),
                 sigma2_dt = stats::var(res)), class = "cmle_fit")
}


#' Rebuild cgns_model from learned polynomial coefficients
#'
#' All kept terms must have hidden-degree <= 1, and noise is carried over from the template.
#'
#' @param coefs_x Learned observed equation coefficients.
#' @param coefs_y Learned hidden equation coefficients.
#' @param library Function library specification.
#' @param template Model supplying dimensions and diffusion coefficients.
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' lib <- cgns_library(k = 1, l = 1, degree = 2)
#' lm1 <- learn_model(m, ob, lib, n_samples = 3, seed = 1)
#' lm1$model
#'
#' @export
model_from_learned <- function(coefs_x, coefs_y, library, template) {
  if (!inherits(library, "function_library"))
    aci_abort("aci_error_model_contract", "library must be a function_library.")
  if (!is.list(template) || is.null(template$k) || is.null(template$l))
    aci_abort("aci_error_model_contract", "template must provide k, l, and diffusion functions.")
  dims <- .state_dims(template$k, template$l)
  k <- dims[["k"]]; l <- dims[["l"]]
  if (library$k != k || library$l != l)
    aci_abort("aci_error_model_contract", "library dimensions do not match the template.")
  if (!is.function(template$Sx1) || !is.function(template$Sy2) ||
      (!is.null(template$Sx2) && !is.function(template$Sx2)) ||
      (!is.null(template$Sy1) && !is.function(template$Sy1)))
    aci_abort("aci_error_model_contract", "template diffusion entries must be functions.")
  if (!is.null(template$meta) && !is.list(template$meta))
    aci_abort("aci_error_model_contract", "template meta must be a list when supplied.")
  validate_coefs <- function(xs, rows, label) {
    if (!is.list(xs) || length(xs) != rows)
      aci_abort("aci_error_model_contract",
                sprintf("%s must be a list with one coefficient vector per equation.", label))
    lapply(xs, function(cf) {
      if (!is.numeric(cf) || any(!is.finite(cf)) ||
          (length(cf) && (is.null(names(cf)) || any(!nzchar(names(cf))) ||
                          anyDuplicated(names(cf)) ||
                          any(!names(cf) %in% library$names))))
        aci_abort("aci_error_model_contract",
                  sprintf("%s contains invalid or unknown named coefficients.", label))
      if (length(cf) && any(library$y_degree[match(names(cf), library$names)] > 1L))
        aci_abort("aci_error_model_contract",
                  "A learned CGNS drift cannot contain terms nonlinear in the hidden state.")
      cf
    })
  }
  coefs_x <- validate_coefs(coefs_x, k, "coefs_x")
  coefs_y <- validate_coefs(coefs_y, l, "coefs_y")
  mk_drift <- function(coef_list, rows) function(t, x, y) {
    Z <- matrix(c(x, y), 1)
    vapply(seq_len(rows), function(i) {
      cf <- coef_list[[i]]
      if (!length(cf)) return(0)
      sum(cf * vapply(names(cf), function(nm) library$funs[[nm]](Z, t), numeric(1)))
    }, numeric(1))
  }
  cgns_from_affine(mk_drift(coefs_x, k), mk_drift(coefs_y, l),
                   Sx = template$Sx1, Sx2 = template$Sx2,
                   Sy_shared = template$Sy1, Sy_hidden = template$Sy2,
                   k = k, l = l, name = "learned",
                   meta = utils::modifyList(template$meta %||% list(),
                     list(provenance = "learned_model_package_implementation")))
}


#' Model discovery from observations + posterior hidden samples
#'
#' jiang2026enkbs s5 loop (sample -> CE prune -> constrained MLE), iterated as
#' a stochastic-EM.
#' n_iter = 1 (default) is a single stochastic-EM M-step: structure recovery
#' is reliable, but hidden equation coefficients carry the posterior path
#' h-transform drift bias (~20-40% on quadratic terms in the dyad demo).
#' n_iter >= 2 re-samples under the learned model (structure frozen after
#' the first pass); this loop is EXPERIMENTAL - without Robbins-Monro
#' damping it can amplify rather than remove the bias. Validation required against MATLAB Codebase.
#'
#' @param model Seed model.
#' @param obs Observed trajectory.
#' @param library Function library specification.
#' @param n_samples Number of posterior hidden paths.
#' @param method Posterior sampling method.
#' @param ce_threshold Causation entropy selection threshold.
#' @param constraints Optional equality constraints.
#' @param energy_pairs Optional energy-conserving equation pairs.
#' @param enforce_cgns Restrict selected terms to the CGNS class.
#' @param n_iter Number of stochastic-EM iterations.
#' @param seed Reproducibility seed.
#' @param burn_frac Fraction discarded before fitting.
#' @param ... Additional arguments passed to path sampling.
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#'
#' Chen, N. and Zhang, Y. (2023). A causality-based learning approach for
#' discovering the underlying dynamics of complex systems from partial
#' observations with stochastic parameterization. *Physica D* **449**, 133743.
#' \doi{10.1016/j.physd.2023.133743}
#' @examples
#' \donttest{
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' lib <- cgns_library(k = 1, l = 1, degree = 2)
#' learn_model(m, ob, lib, n_samples = 3, seed = 1)
#' }
#'
#' @export
learn_model <- function(model, obs, library, n_samples = 20,
                        method = c("ffbs", "enkbs"),
                        ce_threshold = "auto", constraints = NULL,
                        energy_pairs = NULL, enforce_cgns = FALSE,
                        n_iter = 1, seed = 1, burn_frac = 0.1, ...) {
  method <- match.arg(method)
  obs <- as_obs(obs)
  if (!inherits(model, "stochastic_model") || obs$k != model$k)
    aci_abort("aci_error_model_contract", "model and observation dimensions are incompatible.")
  if (!inherits(library, "function_library") || library$k != model$k ||
      library$l != model$l)
    aci_abort("aci_error_model_contract", "library dimensions must match model k and l.")
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) ||
      n_samples < 1L || n_samples != floor(n_samples))
    aci_abort("aci_error_dims", "n_samples must be a positive integer.")
  n_samples <- as.integer(n_samples)
  if (!is.numeric(n_iter) || length(n_iter) != 1L || !is.finite(n_iter) ||
      n_iter < 1L || n_iter != floor(n_iter))
    aci_abort("aci_error_dims", "n_iter must be a positive integer.")
  n_iter <- as.integer(n_iter)
  if (!is.numeric(burn_frac) || length(burn_frac) != 1L ||
      !is.finite(burn_frac) || burn_frac < 0 || burn_frac >= 1)
    aci_abort("aci_error_dims", "burn_frac must be in [0, 1).")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    aci_abort("aci_error_dims", "seed must be one non-negative integer.")
  if (seed + 97 * (n_iter - 1) > .Machine$integer.max)
    aci_abort("aci_error_dims", "seed and n_iter produce a seed outside R's integer range.")
  if (length(enforce_cgns) != 1L || is.na(enforce_cgns) || !is.logical(enforce_cgns))
    aci_abort("aci_error_dims", "enforce_cgns must be TRUE or FALSE.")
  if (!is.null(constraints) && !is.list(constraints))
    aci_abort("aci_error_dims", "constraints must be NULL or a named list by equation.")
  eq_names <- c(paste0("x", seq_len(model$k)), paste0("y", seq_len(model$l)))
  if (!is.null(constraints) &&
      (is.null(names(constraints)) || anyNA(names(constraints)) ||
       any(!nzchar(names(constraints))) || anyDuplicated(names(constraints)) ||
       any(!names(constraints) %in% eq_names)))
    aci_abort("aci_error_dims", sprintf(
      "constraints must be uniquely named by equation: %s.",
      paste(eq_names, collapse = ", ")))
  if (!is.null(constraints) && !is.null(energy_pairs))
    aci_abort("aci_error_not_implemented",
              "Combining per-equation constraints with energy_pairs is not implemented.")
  keep_start <- ceiling(burn_frac * length(obs$t)) + 1L
  keep_end <- length(obs$t) - 1L
  if (keep_start > keep_end)
    aci_abort("aci_error_dims", "burn_frac leaves no increments for model learning.")
  if (n_iter > 1) enforce_cgns <- TRUE
  sampler_model <- model
  frozen_kept <- NULL          # structure is identified once, then held fixed
  for (it in seq_len(n_iter)) {
  Y <- sample_paths(sampler_model, obs, n_samples = n_samples, method = method,
                    seed = seed + 97 * (it - 1), ...)
  N1 <- length(obs$t); dt <- obs$dt
  keep_t <- seq.int(keep_start, keep_end)
  k <- model$k; l <- model$l
  ThetaL <- list(); dXL <- vector("list", k); dYL <- vector("list", l)
  for (s in seq_len(n_samples)) {
    Z <- cbind(obs$x, t(matrix(Y[, , s], l)))[keep_t, , drop = FALSE]
    Th <- eval_library(library, Z, obs$t[keep_t])
    ThetaL[[s]] <- Th
    for (i in seq_len(k)) dXL[[i]] <- c(dXL[[i]], diff(obs$x[, i])[keep_t] / dt)
    for (i in seq_len(l)) dYL[[i]] <- c(dYL[[i]], diff(Y[i, , s])[keep_t] / dt)
  }
  Theta <- do.call(rbind, ThetaL)
  fit_one <- function(dz, eqname) {
    cand <- library$names
    if (enforce_cgns) cand <- cand[library$y_degree <= 1]
    ce <- causation_entropy(Theta[, cand, drop = FALSE], dz, cand)
    kept <- if (!is.null(frozen_kept)) frozen_kept[[eqname]]
            else threshold_structure(ce, ce_threshold)
    if ("1" %in% library$names) kept <- union("1", kept)
    cons <- constraints[[eqname]]
    if (!is.null(cons)) {
      keep2 <- union(kept, cons$involves %||% character())
      kept <- intersect(library$names, keep2)
    }
    fit <- constrained_mle(Theta[, kept, drop = FALSE], dz, constraints = cons)
    list(kept = kept, ce = ce, coef = fit$coef,
         sigma_hat = sqrt(max(fit$sigma2_dt * dt, 0)),
         fit = fit)
  }
  out <- c(lapply(seq_len(k), function(i) fit_one(dXL[[i]], paste0("x", i))),
           lapply(seq_len(l), function(i) fit_one(dYL[[i]], paste0("y", i))))
  names(out) <- c(paste0("dx", seq_len(k)), paste0("dy", seq_len(l)))
  if (!is.null(energy_pairs))
    out <- .refit_energy_joint(out, Theta, c(dXL, dYL), library,
                               energy_pairs, k, l, dt = dt)
  if (it == 1 && n_iter > 1)
    frozen_kept <- stats::setNames(lapply(out, `[[`, "kept"),
      c(paste0("x", seq_len(k)), paste0("y", seq_len(l))))
  if (it < n_iter) {
    cx <- lapply(seq_len(k), function(i) out[[i]]$coef)
    cy <- lapply(seq_len(l), function(i) out[[k + i]]$coef)
    sampler_model <- tryCatch(
      model_from_learned(cx, cy, library, template = model),
      error = function(e) {
        aci_warn("aci_warn_truncation",
                 sprintf("iteration %d: learned model rejected (%s); keeping previous sampler.",
                         it, conditionMessage(e)))
        sampler_model })
  }
  }
  structure(list(equations = out, library = library,
                 meta = list(n_samples = n_samples, method = method,
                             n_iter = n_iter)),
            class = "learned_model")
}


#' Print a learned model
#'
#' @param x A `learned_model` object.
#' @param ... Ignored, for consistency with the generic.
#' @returns `x`, invisibly; called for the printed output.
#' @export
print.learned_model <- function(x, ...) {
  cat("<learned_model>\n")
  for (nm in names(x$equations)) {
    e <- x$equations[[nm]]
    cf <- e$coef[abs(e$coef) > 1e-10]
    cat(sprintf("  %s/dt = %s\n", nm,
        paste(sprintf("%+.3g %s", cf, names(cf)), collapse = " ")))
  }
  invisible(x)
}


#' Joint cross-equation refit under energy-conservation equalities
#'
#' Imposes coef("xi*ym" in dx_i) + coef("xi^2" in dy_m) = 0,
#' the quadratic-pair conservation of the dyad/triad/pathways class. The
#' touched equations are stacked into one block-diagonal KKT system; the
#' involved terms are unioned into the kept sets if CE pruned them.
#'
#' @param out Per-equation fit objects.
#' @param Theta Design matrix.
#' @param targets Per-equation response vectors.
#' @param library Function library specification.
#' @param pairs Observed/hidden energy-pair constraints.
#' @param k Number of observed equations.
#' @param l Number of hidden equations.
#' @returns The refitted per-equation coefficient list.
#' @noRd
.refit_energy_joint <- function(out, Theta, targets, library, pairs, k, l,
                                dt = NULL) {
  if (!is.null(names(pairs)) && all(c("obs", "hid") %in% names(pairs)))
    pairs <- list(pairs)                    # a single bare pair was passed
  if (!is.list(pairs) || !length(pairs))
    aci_abort("aci_error_model_contract", "energy_pairs must be a non-empty list of obs/hid index pairs.")
  pairs <- lapply(pairs, function(p) {
    if (is.null(names(p)) || !all(c("obs", "hid") %in% names(p)))
      aci_abort("aci_error_model_contract", "Each energy pair must name obs and hid indices.")
    i <- p[["obs"]]; m <- p[["hid"]]
    if (!is.numeric(i) || length(i) != 1L || !is.finite(i) || i != floor(i) ||
        i < 1L || i > k || !is.numeric(m) || length(m) != 1L ||
        !is.finite(m) || m != floor(m) || m < 1L || m > l)
      aci_abort("aci_error_model_contract", "Energy-pair obs/hid indices are out of range.")
    c(obs = as.integer(i), hid = as.integer(m))
  })
  kept <- lapply(out, `[[`, "kept")
  for (p in pairs) {
    i <- as.integer(p[["obs"]]); m <- as.integer(p[["hid"]])
    txy <- sprintf("x%d*y%d", i, m); tx2 <- sprintf("x%d^2", i)
    if (!all(c(txy, tx2) %in% library$names))
      aci_abort("aci_error_model_contract",
                sprintf("energy pair (%d, %d): library lacks %s / %s.", i, m, txy, tx2))
    kept[[i]] <- union(kept[[i]], txy)
    kept[[k + m]] <- union(kept[[k + m]], tx2)
  }
  nE <- length(out); pj <- vapply(kept, length, integer(1)); off <- cumsum(c(0, pj))
  nr <- nrow(Theta)
  TJ <- matrix(0, nE * nr, sum(pj)); yJ <- numeric(nE * nr)
  joint_names <- unlist(lapply(seq_len(nE), function(e)
    paste0(names(out)[e] %||% paste0("eq", e), "::", kept[[e]])),
    use.names = FALSE)
  colnames(TJ) <- make.unique(joint_names)
  for (e in seq_len(nE)) {
    rows <- (e - 1) * nr + seq_len(nr)
    TJ[rows, off[e] + seq_len(pj[e])] <- Theta[, kept[[e]], drop = FALSE]
    yJ[rows] <- targets[[e]]
  }
  A <- matrix(0, length(pairs), sum(pj)); b <- rep(0, length(pairs))
  for (r in seq_along(pairs)) {
    i <- as.integer(pairs[[r]][["obs"]]); m <- as.integer(pairs[[r]][["hid"]])
    A[r, off[i] + match(sprintf("x%d*y%d", i, m), kept[[i]])] <- 1
    A[r, off[k + m] + match(sprintf("x%d^2", i), kept[[k + m]])] <- 1
  }
  fit <- constrained_mle(TJ, yJ, constraints = list(A = A, b = b))
  for (e in seq_len(nE)) {
    cf <- fit$coef[off[e] + seq_len(pj[e])]; names(cf) <- kept[[e]]
    res <- targets[[e]] - drop(Theta[, kept[[e]], drop = FALSE] %*% cf)
    ef <- structure(list(coef = cf, rss = sum(res^2), n = length(res),
                         sigma2_dt = stats::var(res)), class = "cmle_fit")
    out[[e]]$kept <- kept[[e]]; out[[e]]$coef <- cf
    out[[e]]$fit <- ef
    if (!is.null(dt)) out[[e]]$sigma_hat <- sqrt(max(ef$sigma2_dt * dt, 0))
    out[[e]]$energy_constrained <- TRUE
  }
  out
}
