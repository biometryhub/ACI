# gate G1: extracted function against the script's own workspace ---------------
#
# The extraction is mechanical, but "mechanical" is a claim about a program,
# and programs have bugs. G1 settles it empirically: call each extracted
# function with the inputs the script itself used, and require every output to
# equal the value the script itself produced. Not to a tolerance, but equal. A
# hoist that changed the arithmetic cannot survive this, and a hoist that
# changed nothing cannot fail it.

#' Build and run the G1 comparison for one captured workspace.
#'
#' @param manifest_path Path to `manifest/extracts.dcf`.
#' @param workspace Path to a captured `.mat` workspace.
#' @param profile Character scalar, used to name the report.
#' @param only Optional character vector of function names to restrict to.
#'
#' @returns A data frame with one row per compared output.
run_gate_g1 <- function(manifest_path, workspace, profile, only = NULL) {
  records <- read_manifest(manifest_path)
  if (!is.null(only)) {
    records <- records[records[, "Function"] %in% only, , drop = FALSE]
  }
  report <- file.path("tools", "oracle", "parity", "reports",
                      sprintf("g1_%s.csv", profile))
  dir.create(dirname(report), recursive = TRUE, showWarnings = FALSE)

  code <- c(
    sprintf("w = load(%s);", as_matlab_path(workspace)),
    "fun = {}; quantity = {}; same = []; maxabs = []; info = {};"
  )

  for (i in seq_len(nrow(records))) {
    name <- records[i, "Function"]
    inputs <- .split_list(records[i, "Inputs"])
    outputs <- .split_list(records[i, "Outputs"])
    locals <- sprintf("g1_%s", outputs)

    code <- c(
      code,
      sprintf("%% ---- %s ----", name),
      sprintf(
        "[%s] = %s(%s);",
        paste(locals, collapse = ", "), name,
        paste(sprintf("w.%s", inputs), collapse = ", ")
      )
    )
    for (j in seq_along(outputs)) {
      code <- c(code, sprintf(
        paste0(
          "[s, m, d] = bitcompare(%s, w.%s); fun{end+1} = '%s'; ",
          "quantity{end+1} = '%s'; same(end+1) = s; maxabs(end+1) = m; ",
          "info{end+1} = d;"
        ),
        locals[j], outputs[j], name, outputs[j]
      ))
    }
  }

  code <- c(
    code,
    paste0(
      "writetable(table(fun(:), quantity(:), logical(same(:)), maxabs(:), ",
      "info(:), 'VariableNames', {'fun','quantity','same','maxabs','info'}), ",
      as_matlab_path(report), ");"
    ),
    sprintf("fprintf('G1 wrote %%d comparisons\\n', numel(fun));")
  )

  run_matlab(code, label = sprintf("g1_%s", profile))
  result <- utils::read.csv(report, stringsAsFactors = FALSE)
  result$profile <- profile
  result
}
