## acir reserve file
## Origin: aci/R/assimilation.R:613-683
## Source package: aci 0.0.30, git tree 97f6b124
## Category: aci-paper
## Intended release: 0.1.x, TBD with the supervisor/collaborators
## Reason: andreou2026smoother eqs. 3.18-3.19; paper-derived with no ACI_code MATLAB backing.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: andreou2026smoother eqs. 3.18-3.19; paper-derived with no ACI_code MATLAB backing (plan section 3.3).

#' Per-step contraction certificate for the Theorem 3 update
#'
#' Evaluates the checkable per-step contraction condition of the online
#' smoother update matrices (andreou2026smoother eqs. 3.18-3.19): when the
#' Hermitian part of `(I - E^j) / dt` is positive definite at every step,
#' the spectral radius of every `E^j` falls below one for a small enough
#' step. Reported per step: `lambda_min`, the smallest eigenvalue of that
#' Hermitian part (positive at the step means condition 3.18 holds there);
#' `enorm`, the operator 2-norm of `E^j`; and `rho_E`, its spectral radius.
#'
#' The condition is stated for the continuous-time generator and holds for
#' `dt` sufficiently small, so a positive `lambda_min` alongside an `enorm`
#' at or above one flags a step-size margin rather than a broken model. The
#' spectral radius is not submultiplicative, so no bound on the accumulated
#' update products follows from `rho_E < 1` alone; [lt_tail_bound()] remains
#' a heuristic estimate.
#'
#' @param x A `lag_table` object.
#' @returns A data frame with one row per grid step, carrying `j`, `t`,
#'   `lambda_min`, `enorm` and `rho_E`, with attributes `gamma` (the largest
#'   `enorm`) and `condition_318` (`TRUE` when every `lambda_min` is
#'   positive).
#'
#' @references
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
#' **36**(4), 71. \doi{10.1007/s00332-026-10271-x}
#'
#' @seealso [lag_table()], [lt_tail_bound()], [truncation_profile()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' tb <- lag_table(m, as_obs(sim), mode = "forward",
#'                 init = list(mean = 0, cov = diag(1, 1)))
#' cert <- lt_contraction_certificate(tb)
#' attr(cert, "gamma"); attr(cert, "condition_318")
#'
#' @export
lt_contraction_certificate <- function(x) {
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  model <- x$meta$source_model
  if (is.null(model) || is.null(x$meta$source_obs_x))
    aci_abort("aci_error_dims",
              "This lag table carries no model/observation handles.")
  obs <- observed_trajectory(x$t, x$meta$source_obs_x)
  bundle <- .compile_cgns_run(model, obs, x$meta$nontarget)
  filt <- .cgns_filter_compiled(
    bundle, x$meta$init, stepper = "explicit", nsub = 1L,
    validate = FALSE
  )
  l <- bundle$l; dt <- x$dt; N1 <- length(x$t)
  lam <- en <- rh <- numeric(N1 - 1L)
  for (j in seq_len(N1 - 1L)) {
    co <- .compiled_co(bundle, j)
    aux <- .thmD1_aux_compiled(bundle, j, filt$cov[, , j], co = co)
    Hm  <- (diag(l) - aux$E) / dt
    Hm  <- (Hm + t(Hm)) / 2
    lam[j] <- min(eigen(Hm, symmetric = TRUE, only.values = TRUE)$values)
    en[j]  <- if (l == 1L) abs(aux$E[1, 1]) else norm(aux$E, "2")
    rh[j]  <- if (l == 1L) abs(aux$E[1, 1]) else
      max(abs(eigen(aux$E, only.values = TRUE)$values))
  }
  out <- data.frame(j = seq_len(N1 - 1L), t = x$t[seq_len(N1 - 1L)],
                    lambda_min = lam, enorm = en, rho_E = rh)
  attr(out, "gamma") <- max(en)
  attr(out, "condition_318") <- all(lam > 0)
  out
}
