# -- S4b: predator-prey, aciR against the reference ----------------------------
#
# The dyad half of S4 builds a dataset and drives the extracted kernels on it.
# Predator-prey cannot be graded that way, because `noisy_predator_prey_model.m`
# must not be run top to bottom (F6): its two causal directions assign the same
# names, so a whole-script run pairs direction one's smoother with direction
# two's filter. The reference values graded here are therefore the ones the
# COMPOSED runner produced -- one direction's declared blocks, in file order --
# and gate G1 has already shown, on 14 outputs per direction at a maximum
# absolute difference of 0, that the extracted kernels reproduce them exactly.
#
# What was missing until now is only the other side of the comparison: aciR had
# never been run on this system at all. This file runs it and compares quantity
# by quantity, on the reference's own conventions.
#
# -- grading like with like ---------------------------------------------------
#
# The causal influence range is the quantity with the designed differences, and
# getting the setup wrong here once produced a 4.58e-09 "disagreement" that was
# not a disagreement at all -- two different questions differenced. The three
# arguments that make the questions the same are:
#
#   * `horizon = last_idx`. The reference compares each reporting time against
#     later observations only as far as its own last_idx, where aciR's default
#     runs to the end of the record. Different horizons, different integrals.
#   * `margin` stood down. aciR flags a time whose range consumes most of the
#     sequence it was measured against; the reference has no such flag, so the
#     flag must not be allowed to remove times from the comparison.
#   * the reference's own threshold grid, 513 logarithmic points, against
#     aciR's default of 129.
#
# aciR's default horizon is reported too, separately and with an infinite
# tolerance, so the size of the design difference is on the record rather than
# hidden by the recipe that removes it.

# -- tolerances, declared before the first run ---------------------------------
#
# Declared here, in the file, before any predator-prey number had been measured,
# and taken from the classes the dyad comparison already uses rather than chosen
# for this system: 1e-12 for the linear recursions, 1e-11 for anything that
# passes through a logarithm, a division by a covariance or a quadrature. A
# tolerance picked after seeing the difference is a description, not a test.
#
# `cir_objective_exact` is new to this comparison -- the dyad's driver never
# compared it -- and is declared at 1e-11 for the same reason its neighbours
# are: it is a quadrature over a 513-point grid divided by a peak.
.predprey_tolerances <- c(
  filter_mean = 1e-12,
  filter_cov = 1e-12,
  smoother_mean = 1e-12,
  smoother_cov = 1e-12,
  ACI_metric = 1e-11,
  E_j = 1e-12,
  F_j = 1e-12,
  online_mean = 1e-11,
  online_cov = 1e-11,
  cir_peak = 1e-11,
  cir_objective = 1e-11,
  cir_subjective = 1e-11,
  cir_objective_exact = 1e-11
)

# -- export --------------------------------------------------------------------

#' Export one captured predator-prey workspace for comparison.
#'
#' @param direction Integer scalar, 1 or 2.
#' @param workspace_root Directory holding the captured profiles.
#' @param bundle_root Directory to write the dataset bundle into.
#' @param report_dir Directory to write the reference side into.
#'
#' @returns The bundle directory, invisibly.
export_predprey <- function(direction,
                            workspace_root = file.path("oracle", "parity",
                                                       "workspaces"),
                            bundle_root = file.path("oracle", "parity",
                                                    "datasets"),
                            report_dir = file.path("oracle", "parity",
                                                   "reports")) {
  profile <- sprintf("predprey_dir%d", direction)
  workspace <- file.path(workspace_root, profile,
                         "noisy_predator_prey_model_workspace.mat")
  if (!file.exists(workspace)) {
    stop(
      sprintf("no captured workspace at %s; run the capture first.", workspace),
      call. = FALSE
    )
  }
  bundle <- file.path(bundle_root, profile)
  run_matlab(
    sprintf("export_predprey(%s, %d, %s, %s);",
            as_matlab_path(workspace), as.integer(direction),
            as_matlab_path(bundle), as_matlab_path(report_dir)),
    label = sprintf("export_%s", profile)
  )
  invisible(bundle)
}

#' Read a predator-prey bundle.
#'
#' The dyad bundle carries the noise as four feedback matrices and a scalar
#' latent self-drift. Neither holds here: under the direction-one mirror the
#' four feedback names change role as well as value, and both directions give
#' the latent process a self-drift that varies with the observed state. The
#' bundle therefore carries the Grammians the reference itself formed, and a
#' per-step `L_y`, which `aci_cir()` and the core both accept.
#'
#' @param path Directory holding `meta.dcf` and `arrays.csv`.
#'
#' @returns A list with `meta`, `arrays` and the aciR components list.
read_predprey_bundle <- function(path) {
  meta <- as.list(read_manifest(file.path(path, "meta.dcf"))[1L, ])
  numeric_fields <- c("N", "dt", "mu0", "R0", "S_xoS_x", "S_yoS_y", "S_yoS_x",
                      "S_xoS_y", "CIRStart", "CIREnd", "FirstIdx", "LastIdx",
                      "PlotLen", "PlotEndIdx", "Lookahead", "EpsilonResolution",
                      "FixedLag", "DefnLen", "Direction")
  for (field in numeric_fields) {
    if (!field %in% names(meta)) {
      stop(sprintf("bundle %s declares no `%s`.", path, field), call. = FALSE)
    }
    meta[[field]] <- as.numeric(meta[[field]])
  }
  arrays <- utils::read.csv(file.path(path, "arrays.csv"),
                            stringsAsFactors = FALSE)
  if (nrow(arrays) != meta$N + 1L) {
    stop(
      sprintf("arrays.csv has %d rows; meta declares N = %d, so %d expected.",
              nrow(arrays), meta$N, meta$N + 1L),
      call. = FALSE
    )
  }
  comp <- list(
    L_x = arrays$L_x, L_y = arrays$L_y,
    f_x = arrays$f_x, f_y = arrays$f_y,
    S_xoS_x = meta$S_xoS_x, S_yoS_y = meta$S_yoS_y,
    S_yoS_x = meta$S_yoS_x, S_xoS_y = meta$S_xoS_y
  )
  list(meta = meta, arrays = arrays, comp = comp)
}

# -- the comparison ------------------------------------------------------------

#' Run aciR on one predator-prey direction and compare with the reference.
#'
#' @param direction Integer scalar, 1 or 2.
#' @param tolerances Named numeric vector of relative tolerances, declared in
#'   advance at the top of this file.
#' @param bundle_root Directory holding the exported bundles.
#' @param report_dir Directory holding the exported reference side.
#'
#' @returns A data frame with one row per compared quantity.
compare_predprey <- function(direction,
                             tolerances = .predprey_tolerances,
                             bundle_root = file.path("oracle", "parity",
                                                     "datasets"),
                             report_dir = file.path("oracle", "parity",
                                                    "reports")) {
  name <- sprintf("predprey_dir%d", direction)
  bundle <- read_predprey_bundle(file.path(bundle_root, name))
  meta <- bundle$meta
  comp <- bundle$comp
  x <- bundle$arrays$x
  stem <- file.path(report_dir, sprintf("matlab_%s", name))

  reference <- utils::read.csv(sprintf("%s.csv", stem),
                               stringsAsFactors = FALSE)

  filt <- aciR::aci_filter(x, comp, meta$dt, mu0 = meta$mu0, R0 = meta$R0)
  smooth <- aciR::aci_smoother(x, comp, meta$dt, filt)
  metric <- aciR::aci_metric(filt, smooth)
  aux <- getFromNamespace(".aci_online_aux", "aciR")(x, comp, meta$dt, filt)

  ours <- list(
    filter_mean = filt$mean, filter_cov = filt$cov,
    smoother_mean = smooth$mean, smoother_cov = smooth$cov,
    ACI_metric = metric, E_j = aux$E, F_j = aux$F
  )
  rows <- lapply(names(ours), function(quantity) {
    .compare_vectors(quantity, ours[[quantity]], reference[[quantity]],
                     tolerances[[quantity]], name)
  })

  # ---- the online smoother ---------------------------------------------------
  #
  # The reference sets fixed_lag = N+1 unconditionally, so the full path is the
  # only lag it defines and the only one there is anything to compare against.
  # aciR's `Inf` is that lag.
  lag_path <- sprintf("%s_online_lag%d.csv", stem, as.integer(meta$FixedLag))
  if (!file.exists(lag_path)) {
    stop(sprintf("no online-smoother export at %s.", lag_path), call. = FALSE)
  }
  theirs_lag <- utils::read.csv(lag_path, stringsAsFactors = FALSE)
  online <- aciR::aci_online_smoother(x, comp, meta$dt, filt, lag = Inf)
  rows <- c(rows, list(
    .compare_vectors("online_mean(lag=Inf)", online$mean,
                     theirs_lag$online_mean, tolerances[["online_mean"]], name),
    .compare_vectors("online_cov(lag=Inf)", online$cov,
                     theirs_lag$online_cov, tolerances[["online_cov"]], name)
  ))

  rows <- c(rows, .compare_predprey_cir(stem, meta, comp, filt, x, tolerances))

  result <- do.call(rbind, rows)
  .assert_graded(result, name)
  result
}

# -- causal influence range ----------------------------------------------------

.compare_predprey_cir <- function(stem, meta, comp, filt, x, tolerances) {
  theirs <- utils::read.csv(sprintf("%s_cir.csv", stem),
                            stringsAsFactors = FALSE)
  epsilon <- utils::read.csv(sprintf("%s_cir_epsilon.csv", stem),
                             header = FALSE)[[1L]]
  their_subjective <- as.matrix(utils::read.csv(
    sprintf("%s_cir_subjective.csv", stem), header = FALSE))
  their_defn <- utils::read.csv(sprintf("%s_cir_defn.csv", stem),
                                header = FALSE)[[1L]]
  name <- meta$Name

  # ---- the region the reference stands behind --------------------------------
  #
  # The reference computes over first_idx..last_idx but plots only to
  # time_end_plot; the trailing `lookahead_tolerance` worth of reporting times
  # is a guard it computes and discards. Comparing over the guard compares
  # against numbers the reference does not stand behind. Its own exact objective
  # range is defined over exactly the plotted region, which is the cross-check
  # asserted below rather than assumed.
  plot_end <- as.integer(meta$PlotEndIdx)
  reported <- theirs$index <= plot_end
  last_idx <- as.integer(meta$LastIdx)
  if (sum(reported) != length(their_defn)) {
    stop(
      sprintf(
        paste0(
          "the plotted region holds %d times but the reference's exact ",
          "objective range covers %d. One of the two is not what it is ",
          "believed to be, and nothing is graded until they agree."
        ),
        sum(reported), length(their_defn)
      ),
      call. = FALSE
    )
  }
  if (sum(reported) == 0L) {
    stop("the plotted region is empty; there would be nothing to compare.",
         call. = FALSE)
  }
  if (nrow(their_subjective) != length(epsilon) ||
        ncol(their_subjective) != nrow(theirs)) {
    stop(
      sprintf(
        "subjective_CIR is %d x %d; expected %d thresholds by %d times.",
        nrow(their_subjective), ncol(their_subjective), length(epsilon),
        nrow(theirs)
      ),
      call. = FALSE
    )
  }

  # ---- the reference's own conventions ---------------------------------------
  same <- aciR::aci_cir(
    x, comp, meta$dt, filt = filt, window = theirs$index, epsilon = epsilon,
    margin = 1e-9, horizon = last_idx
  )
  ours_sub <- same$subjective[order(epsilon), reported, drop = FALSE]
  their_sub <- their_subjective[order(epsilon), reported, drop = FALSE]
  both_sub <- is.finite(ours_sub) & is.finite(their_sub)

  # aciR withholds the exact objective range where the peak divergence is below
  # its detection threshold; the reference divides by the tiny peak anyway. Only
  # the times both resolve are compared, and how many were dropped is reported
  # rather than left implicit.
  ours_defn <- same$objective_exact[reported]
  both_defn <- is.finite(ours_defn) & is.finite(their_defn)

  out <- list(
    .compare_vectors("cir_peak [reference conventions]", same$peak[reported],
                     theirs$peak[reported], tolerances[["cir_peak"]], name),
    .compare_vectors("cir_objective [reference conventions]",
                     same$objective[reported], theirs$objective[reported],
                     tolerances[["cir_objective"]], name),
    .compare_vectors("cir_subjective [reference conventions]",
                     ours_sub[both_sub], their_sub[both_sub],
                     tolerances[["cir_subjective"]], name),
    .compare_vectors("cir_objective_exact [reference conventions]",
                     ours_defn[both_defn], their_defn[both_defn],
                     tolerances[["cir_objective_exact"]], name)
  )

  # ---- aciR's default horizon: a designed difference, measured not hidden ----
  default <- aciR::aci_cir(
    x, comp, meta$dt, filt = filt, window = theirs$index, epsilon = epsilon
  )
  out <- c(out, list(
    .compare_vectors("cir_objective [aciR default -- by design]",
                     default$objective[reported], theirs$objective[reported],
                     Inf, name)
  ))

  message(sprintf(
    paste0(
      "   CIR %s: %d times, %d reported; at reference conventions %d ",
      "censored, %d of %d subjective cells and %d of %d exact-objective ",
      "times compared"
    ),
    name, length(same$peak), sum(reported),
    sum(same$status == "censored"), sum(both_sub), length(both_sub),
    sum(both_defn), length(both_defn)
  ))
  out
}

# -- driver --------------------------------------------------------------------

#' Export and compare both causal directions, and write the report.
#'
#' Source `tools/manifest.R`, `tools/matlab.R`, `tools/compare.R` and this file
#' from the repository root, then call it. `export = FALSE` re-grades from the
#' exports already on disk without paying for MATLAB again.
#'
#' @param directions Integer vector of causal directions to run.
#' @param export Logical. Re-export from the captured workspaces first.
#' @param out Path to write the report to.
#'
#' @returns The combined data frame, invisibly.
run_predprey_parity <- function(directions = c(1L, 2L), export = TRUE,
                                out = file.path("oracle", "parity", "reports",
                                                "parity_predprey.csv")) {
  rows <- lapply(directions, function(direction) {
    if (isTRUE(export)) {
      export_predprey(direction)
    }
    message(sprintf("-- comparing predator-prey direction %d", direction))
    compare_predprey(direction)
  })
  result <- do.call(rbind, rows)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, out, row.names = FALSE)
  message(sprintf("wrote %d comparisons to %s", nrow(result), out))
  invisible(result)
}

# -- the guard against a vacuous pass ------------------------------------------
#
# `max(numeric(0))` is -Inf, which sails under any threshold, and this harness
# once reported `ok` for having compared nothing. `.compare_vectors()` names an
# empty comparison rather than passing it; this asserts that no row of the
# result is empty, mismatched or over budget, so a driver cannot walk past a
# failure it printed.
.assert_graded <- function(result, name) {
  if (is.null(result) || nrow(result) == 0L) {
    stop(sprintf("%s produced no comparisons at all.", name), call. = FALSE)
  }
  empty <- is.na(result$n) | result$n == 0L
  if (any(empty)) {
    stop(
      sprintf("%s compared nothing for: %s", name,
              paste(result$quantity[empty], collapse = ", ")),
      call. = FALSE
    )
  }
  bad <- result$verdict != "ok"
  if (any(bad)) {
    message(sprintf(
      "%s: %d quantit%s did not return `ok` -- %s", name, sum(bad),
      if (sum(bad) == 1L) "y" else "ies",
      paste(sprintf("%s (%s)", result$quantity[bad], result$verdict[bad]),
            collapse = ", ")
    ))
  }
  invisible(TRUE)
}
