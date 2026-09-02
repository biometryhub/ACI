## acir reserve file
## Origin: aci/R/causal_metrics.R:780-979
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Backward CIR is FBCIR-only (FBCIR_code-main/climate_tipping_y_bifurcation_driven.m eq. 21); absent from ACI_code.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: backward_cir generic plus the .lag_table and .aci_result methods, including the ensemble guard at 907-911.

#' Backward causal influence range
#'
#' Summarizes the duration of influence on the discrete time grid, backward from
#' a reference time. A backward result is anchored at `T`, but its computable
#' lagged grid and physical interval end at `T - dt`; its bound label records
#' this convention. The backward `l1_linf` ratio is evaluated with the plain
#' Appendix G L1 sum, following the FBCIR code's active line; the exact form
#' integrates its suffix minima with composite Simpson.
#'
#' A backward result also reports the andreou2026cir Section 2.3.4 validity
#' gate in its `meta`: `baseline` (the initial-time information deficit, the
#' one-lag metric at t = 0), `terminal` (the metric at the `T - dt` end of the
#' grid), and `above_baseline` (their comparison, eq. 21 and Remark C.2).
#' A large `M` with `above_baseline = FALSE` means the attribution evidence
#' never rises above the initial-time baseline; treat such a range with
#' caution rather than as a long attribution window.
#'
#' @param x A `lag_table`, `aci_result` or `cgns_model` object.
#' @param ... Arguments passed to methods.
#' @returns An object of class `cir_result` with `direction` `"backward"`.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' @seealso [forward_cir()], [lag_table()], [cir_pair()]
#'
#' @examples
#' m <- model_dyad()
#' sim <- simulate(m, seed = 1, T = 2, dt = 0.01)
#' ob <- as_obs(sim)
#' tb <- lag_table(m, ob, mode = "one_lag")
#' backward_cir(tb)
#'
#' @export
backward_cir <- function(x, ...) UseMethod("backward_cir")


#' @describeIn backward_cir Backward range from a one-lag table.
#' @param method Either `"exact"` or `"l1_linf"`.
#' @param eps Optional numeric vector of finite non-negative thresholds at which
#'   the subjective range is also reported.
#' @param min_M Either `"auto"` or one finite non-negative number; a range whose
#'   strength falls below it is masked.
#' @param masked_value What a masked range reports: `"na"` (the default) keeps
#'   the masking visible as `NA`; `"zero"` follows andreou2026cir Remarks B.4
#'   and C.4 and the FBCIR scripts, which set the length to 0.
#' @param quadrature Quadrature for the `l1_linf` ratio: `"sum"` (the default)
#'   is the literal andreou2026cir eq. G.14 L1 grid-function sum, the FBCIR
#'   scripts' active line; `"simpson"` is the smoother alternative those
#'   scripts retain in comments. The choice is recorded in the result's
#'   `meta$quadrature`; the `exact` form always integrates with Simpson.
#' @param simpson_close Closing rule for the leftover interval when a grid has
#'   an even number of points: `"quadratic"` (the default) fits a quadratic
#'   through the last three points, following `simps.m` in the ACI reference
#'   code; `"trapezoid"` is the package's pre-0.0.21 rule and reproduces
#'   results reported by earlier versions. Odd-length grids are unaffected.
#' @export
backward_cir.lag_table <- function(x, method = c("exact", "l1_linf"),
                                   eps = NULL, min_M = "auto",
                                   masked_value = c("na", "zero"),
                                   quadrature = c("sum", "simpson"),
                                   simpson_close = c("quadratic", "trapezoid"),
                                   ...) {
  if (length(list(...)))
    aci_abort("aci_error_dims", "Unused arguments were supplied to backward_cir().")
  masked_value <- match.arg(masked_value)
  quadrature <- match.arg(quadrature)
  simpson_close <- match.arg(simpson_close)
  method <- match.arg(method)
  if (!inherits(x, "lag_table"))
    aci_abort("aci_error_dims", "x must be a lag_table.")
  P <- lt_onelag(x)
  r <- .bwd_lengths(P, x$dt, method, quadrature, simpson_close)
  ## andreou2026cir Section 2.3.4 validity gate: causal intensity at the
  ## T - dt end of the grid versus the initial-time information deficit
  ## P^0_T (eq. 21; Remark C.2).  The gate is a diagnostic beside tau and M,
  ## not a change to the centred metric.
  baseline <- unname(P[1L])
  terminal <- unname(P[length(P)])
  Tend <- x$t[length(x$t)]
  lag_end <- Tend - x$dt
  subj <- NULL
  if (!is.null(eps)) {
    if (!is.numeric(eps) || !length(eps) || any(!is.finite(eps)) || any(eps < 0))
      aci_abort("aci_error_dims", "eps must contain finite non-negative thresholds.")
    cc <- abs(P - P[1]); pt <- x$t[seq_along(P)]
    subj <- vapply(eps, function(e) {
      hit <- which(cc <= e)
      if (!length(hit)) NA_real_ else lag_end - max(pt[hit])
    }, numeric(1))
    names(subj) <- format(eps)
  }
  if (!is.null(min_M) && is.finite(r["M"])) {
    mm <- .cir_min_strength(min_M)
    if (r["M"] < mm) {
      r["tau"] <- if (masked_value == "zero") 0 else NA_real_
      aci_warn("aci_warn_low_signal", sprintf(
        "Backward CIR masked (M < %.3g); interpret it jointly with ACI.", mm))
    }
  }
  new_cir_result(Tend, unname(r["tau"]), unname(r["M"]), "backward", method,
                 x$dt, interval = cbind(lag_end - unname(r["tau"]), lag_end),
                 subjective = subj,
                 meta = list(quadrature = quadrature,
                             simpson_close = simpson_close,
                             reference_time = Tend,
                             lagged_grid_end = lag_end,
                             endpoint_convention = "T_minus_dt",
                             dt_caveat = "one-lag limit approximated to O(dt) (andreou2026cir eq. G.12)",
                             baseline = baseline,
                             terminal = terminal,
                             above_baseline = terminal > baseline,
                             table_nontarget = x$meta$nontarget),
                 bound = if (method == "exact")
                   "layer_cake_on_O(dt)_T_minus_dt_grid" else
                   "upper_ratio_on_O(dt)_T_minus_dt_grid")
}


#' @describeIn backward_cir Backward range from an ACI result, at one or more
#'   reference times.
#' @param T Either `"end"` or numeric reference times within the observed grid.
#' @export
backward_cir.aci_result <- function(x, T = "end", method = "exact", ...) {
  if (identical(x$meta$engine, "ensemble") ||
      identical(x$table$meta$engine %||% NULL, "ensemble"))
    aci_abort("aci_error_not_implemented", paste(
      "jiang2026enkbs implements the EnKBS forward CIR but identifies ensemble backward",
      "CIR as future work; this package does not invent that extension."))
  mdl <- x$handles$model; ob <- x$handles$obs; nt <- x$handles$nontarget
  Ts <- if (identical(T, "end")) ob$t[length(ob$t)] else as.numeric(T)
  if (!length(Ts) || any(!is.finite(Ts)) ||
      any(Ts < ob$t[1] - 1e-12) || any(Ts > ob$t[length(ob$t)] + 1e-12))
    aci_abort("aci_error_dims",
              "T must contain finite reference times within the observation grid span.")
  if (length(Ts) > 20)
    aci_warn("aci_warn_truncation",
             "backward_cir over >20 reference times is O(|T| N l^3) in v0.")
  prefix_lengths <- vapply(Ts, function(Tv)
    sum(ob$t <= Tv + 1e-12), integer(1))
  if (any(prefix_lengths < 5L))
    aci_abort("aci_error_dims", "Reference time too early.")
  # Coefficients after the latest requested reference time are not part of any
  # backward-CIR estimand. Compile the maximal required prefix once so multiple
  # references share work without evaluating irrelevant future coefficients.
  max_N1 <- max(prefix_lengths)
  prefix_obs <- if (max_N1 == length(ob$t)) ob else observed_trajectory(
    ob$t[seq_len(max_N1)], ob$x[seq_len(max_N1), , drop = FALSE],
    names = colnames(ob$x)
  )
  bundle <- .compile_cgns_run(mdl, prefix_obs, nt)
  full_filter <- x$paths$filter %||% NULL
  can_reuse <- !is.null(full_filter) &&
    nrow(full_filter$mean) >= max_N1 &&
    identical(full_filter$meta$stepper %||% "explicit", "explicit") &&
    (full_filter$meta$nsub %||% 1L) == 1L
  if (can_reuse) {
    full_filter <- .slice_compiled_filter(full_filter, bundle)
    compatible <- tryCatch({
      .validate_gaussian_path(
        full_filter, bundle$obs, bundle$l, "filter", bundle$nontarget,
        model = bundle$model, source_model = bundle$source_model
      )
      TRUE
    }, error = function(e) FALSE)
    can_reuse <- isTRUE(compatible)
  }
  if (!can_reuse)
    full_filter <- .cgns_filter_compiled(
      bundle, init = x$handles$init, stepper = "explicit", nsub = 1L,
      validate = FALSE
    )
  res <- lapply(prefix_lengths, function(N1) {
    prefix <- .slice_compiled_cgns(bundle, N1)
    prefix_filter <- .slice_compiled_filter(full_filter, prefix)
    tb <- .lag_table_compiled(
      prefix, mode = "one_lag", filter = prefix_filter,
      init = prefix_filter$meta$init, validate = FALSE
    )
    backward_cir(tb, method = method, ...)
  })
  actual_T <- vapply(res, function(z) as.numeric(z$t)[1], numeric(1))
  if (any(abs(actual_T - Ts) > 1e-10))
    aci_warn("aci_warn_grid_snap",
             "Some requested backward-CIR reference times were snapped down to the observation grid.")
  if (length(res) == 1) return(res[[1]])
  subj <- lapply(res, `[[`, "subjective")
  if (all(vapply(subj, is.null, logical(1)))) subj <- NULL
  common_bound <- unique(vapply(res, function(z) z$bound, character(1)))
  new_cir_result(actual_T, sapply(res, `[[`, "tau"), sapply(res, `[[`, "M"),
                 "backward", method, ob$dt,
                 interval = do.call(rbind, lapply(res, `[[`, "interval")),
                 subjective = subj,
                 meta = list(per_reference = lapply(res, `[[`, "meta")),
                 bound = if (length(common_bound) == 1L) common_bound else
                   "mixed_reference_semantics")
}
