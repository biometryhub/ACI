# gate G1: capture and validate a reference workspace --------------------------
#
# Runs each declared (script, profile) through MATLAB and checks the result is
# not degenerate before anything is graded against it. The check is necessary:
# a reduced problem size silently emptied one of the reference's own outputs
# the first time this was run, and a pairing graded against an empty vector
# passes without comparing anything.

#' Capture every declared workspace and assert it is usable.
#'
#' @param manifest_path Path to `manifest/knobs.dcf`.
#' @param reference_dir Directory holding the unmodified reference scripts.
#' @param out_dir Directory holding the reduced scripts and workspaces.
#'
#' @returns A data frame with one row per (script, profile).
capture_all <- function(manifest_path,
                        reference_dir = "matlab_reference",
                        out_dir = file.path("tools", "oracle", "parity", "workspaces"),
                        profiles = NULL) {
  records <- read_manifest(manifest_path)
  if (!is.null(profiles)) {
    records <- records[records[, "Profile"] %in% profiles, , drop = FALSE]
  }
  results <- list()

  for (i in seq_len(nrow(records))) {
    script <- manifest_field(records, i, "Script")
    profile <- manifest_field(records, i, "Profile")
    stem <- sub("\\.m$", "", script)
    reduced <- file.path(out_dir, profile, script)
    workspace <- file.path(out_dir, profile, sprintf("%s_workspace.mat", stem))

    conditions <- character(0)
    if ("Requires" %in% colnames(records) && !is.na(records[i, "Requires"])) {
      raw <- gsub("[\r\n]+", " ", records[i, "Requires"])
      conditions <- trimws(strsplit(raw, ";", fixed = TRUE)[[1L]])
      conditions <- conditions[nzchar(conditions)]
    }

    message(sprintf("-- capturing %s / %s", script, profile))
    output <- run_matlab(
      c(
        sprintf(
          "capture_workspace(%s, %s, %s);",
          as_matlab_path(reduced), as_matlab_path(reference_dir),
          as_matlab_path(workspace)
        ),
        sprintf(
          "check_profile(%s, %s);",
          as_matlab_path(workspace), as_matlab_cell(conditions)
        )
      ),
      label = sprintf("capture_%s_%s", stem, profile)
    )
    cat(paste(output, collapse = "\n"), "\n")

    results[[length(results) + 1L]] <- data.frame(
      script = script, profile = profile, workspace = workspace,
      conditions = length(conditions), stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}
