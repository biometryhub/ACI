# -- S5: re-grade the committed fixtures against the authors' own code ---------
#
# The fixtures under aciR/inst/extdata were produced by harnesses in oracle/
# that transcribe the reference's algebra by hand. Those harnesses and the R
# implementation share an author, so a transcription error common to both is
# invisible to every test the package owns -- the oracle would agree with the
# code because it came from the same reading of the paper, not because either
# is right.
#
# This step removes that shared step. The reference's setup, signal generation,
# filter and smoother are all O(N), so they can be run at the PUBLISHED N =
# 30000 through the byte-verified extracts, and the fixture compared against
# values the authors' code produced. Anything that disagrees is a finding about
# the fixture, not about the reference.

#' Regenerate the dyad reference quantities at publication settings.
#'
#' @param out_path Where MATLAB should write the subsampled quantities.
#' @param stride Integer. Subsampling stride, matching the committed fixture.
#'
#' @returns The path written.
regrade_dyad_reference <- function(
    out_path = file.path("oracle", "parity", "reports",
                         "regrade_dyad_publication.csv"),
    stride = 100L) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  run_matlab(
    c(
      paste0(
        "[x, y, L_x, f_x, f_y, L_y, Sx_1, Sx_2, Sy_1, Sy_2, N, dt, T, gamma, ",
        "d_x, d_y, F_x, F_y, sigma_x, sigma_y] = ref_dyad_setup_and_signals();"
      ),
      paste0(
        "[filter_mean, filter_cov, S_xoS_x, S_xoS_x_inv, S_yoS_y, S_yoS_x, ",
        "S_xoS_y] = ref_filter_scalar(x, y, L_x, f_x, f_y, L_y, Sx_1, Sx_2, ",
        "Sy_1, Sy_2, N, dt);"
      ),
      paste0(
        "[smoother_mean, smoother_cov, E_j_matrices, F_j_matrices] = ",
        "ref_smoother_scalar(x, filter_mean, filter_cov, L_x, f_x, f_y, L_y, ",
        "S_xoS_x_inv, S_yoS_y, S_yoS_x, S_xoS_y, N, dt);"
      ),
      paste0(
        "[ACI_metric, sig, cvr, disp] = ref_aci_metric_scalar(filter_mean, ",
        "filter_cov, smoother_mean, smoother_cov);"
      ),
      sprintf("idx = 1:%d:N+1; t = (idx-1)*dt;", stride),
      paste0(
        "writetable(table(t(:), x(idx)', y(idx)', L_x(idx)', f_x(idx)', ",
        "f_y(idx)', filter_mean(idx)', filter_cov(idx)', smoother_mean(idx)', ",
        "smoother_cov(idx)', ACI_metric(idx)', 'VariableNames', ",
        "{'t','x','y','L_x','f_x','f_y','filter_mean','filter_cov',",
        "'smoother_mean','smoother_cov','ACI_metric'}), ",
        as_matlab_path(out_path), ");"
      ),
      "fprintf('regrade wrote %d rows at N = %d\\n', numel(idx), N);"
    ),
    label = "regrade_dyad"
  )
  out_path
}

#' Compare a committed fixture against the re-graded values.
#'
#' @param fixture_path Path to the committed fixture CSV.
#' @param regrade_path Path to the re-graded CSV.
#' @param columns Character vector of columns to compare.
#'
#' @returns A data frame with one row per column.
compare_to_fixture <- function(fixture_path, regrade_path, columns) {
  fixture <- utils::read.csv(fixture_path, stringsAsFactors = FALSE)
  regrade <- utils::read.csv(regrade_path, stringsAsFactors = FALSE)
  n <- min(nrow(fixture), nrow(regrade))

  rows <- lapply(columns, function(column) {
    if (!column %in% names(fixture) || !column %in% names(regrade)) {
      return(data.frame(
        quantity = column, n = NA_integer_, max_abs = NA_real_,
        max_rel = NA_real_, verdict = "ABSENT", stringsAsFactors = FALSE
      ))
    }
    a <- fixture[[column]][seq_len(n)]
    b <- regrade[[column]][seq_len(n)]
    absolute <- abs(a - b)
    scale <- pmax(abs(a), abs(b))
    relative <- ifelse(scale > 0, absolute / scale, 0)
    data.frame(
      quantity = column, n = n,
      max_abs = max(absolute), max_rel = max(relative),
      verdict = if (max(absolute) == 0) {
        "exact"
      } else if (max(relative) < 1e-12) {
        "machine"
      } else {
        "DIVERGENT"
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
