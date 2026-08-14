# -- hoisting reference kernels into callable functions ------------------------
#
# The reference is seven top-to-bottom scripts. To run its algebra on a dataset
# other than the one the script generates, the computational passages have to
# become functions -- and the moment we retype algebra we have re-implemented
# it, which is precisely the failure the harness exists to detect.
#
# So nothing is retyped. Each generated function is a signature, a byte-exact
# slice of the reference between two markers, and `end`. `check_extraction()`
# re-reads the slice from the reference and compares bytes, so a generated file
# that drifted from its source cannot pass. The generated files are not
# committed: they are derived from a repository we do not vendor, and
# regenerating them is one command.

.slice_lines <- function(spec) {
  parts <- trimws(strsplit(spec, ",", fixed = TRUE)[[1L]])
  bounds <- lapply(parts, function(p) {
    ends <- as.integer(trimws(strsplit(p, "-", fixed = TRUE)[[1L]]))
    if (length(ends) != 2L || anyNA(ends) || ends[2L] < ends[1L]) {
      stop(sprintf("malformed line range `%s`.", p), call. = FALSE)
    }
    seq.int(ends[1L], ends[2L])
  })
  unlist(bounds)
}

.split_list <- function(spec) {
  if (is.na(spec)) {
    return(character(0))
  }
  out <- trimws(strsplit(gsub("[\r\n]+", " ", spec), ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

.begin_marker <- function(script, spec) {
  sprintf("%% >>> BEGIN VERBATIM %s:%s", script, spec)
}

.end_marker <- function(script, spec) {
  sprintf("%% <<< END VERBATIM %s:%s", script, spec)
}

# -- generation ----------------------------------------------------------------

#' Generate one extracted MATLAB function.
#'
#' @param record A one-row manifest slice, as a named character vector.
#' @param reference_dir Directory holding the unmodified reference scripts.
#' @param out_dir Directory to write the generated function into.
#'
#' @returns The path written.
extract_function <- function(record, reference_dir, out_dir) {
  script <- record[["Script"]]
  spec <- record[["Lines"]]
  name <- record[["Function"]]
  inputs <- .split_list(record[["Inputs"]])
  outputs <- .split_list(record[["Outputs"]])

  source_file <- .read_lines_preserving(file.path(reference_dir, script))
  index <- .slice_lines(spec)
  if (max(index) > length(source_file$lines)) {
    stop(
      sprintf("%s: range %s exceeds %s (%d lines).", name, spec, script,
              length(source_file$lines)),
      call. = FALSE
    )
  }
  body <- source_file$lines[index]

  # ---- consumption assertions ------------------------------------------------
  #
  # An input the slice touches only at one index is safe to supply as a stub on
  # an arbitrary dataset -- the reference's filter reads the true latent path
  # `y` only at y(1), as its initial mean. Asserting that is what lets the
  # dataset contract pass a one-element stand-in without hand-waving.
  consumed <- .split_list(record[["ConsumedOnly"]])
  # Comments mention the variables they describe -- "% ... the latent variable
  # y." -- so the check runs on code only. This drops anything after a `%`,
  # which would also truncate a `%` inside a string literal; none of the
  # extracted slices contain one, and the extraction check would fail loudly if
  # that ever changed.
  code_only <- sub("%.*$", "", body)
  for (rule in consumed) {
    variable <- trimws(sub("=.*$", "", rule))
    allowed <- trimws(sub("^[^=]*=", "", rule))
    occurrences <- gregexpr(sprintf("\\b%s\\b\\s*(\\([^)]*\\))?", variable),
                            code_only)
    hits <- unlist(regmatches(code_only, occurrences))
    hits <- hits[nzchar(hits)]
    illegal <- hits[hits != sprintf("%s(%s)", variable, allowed)]
    if (length(illegal) > 0L) {
      stop(
        sprintf(
          paste0(
            "%s declares `%s` consumed only as %s(%s), but the slice also ",
            "uses it as: %s"
          ),
          name, variable, variable, allowed,
          paste(unique(illegal), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  signature <- if (length(outputs) == 0L) {
    sprintf("function %s(%s)", name, paste(inputs, collapse = ", "))
  } else {
    sprintf(
      "function [%s] = %s(%s)",
      paste(outputs, collapse = ", "), name, paste(inputs, collapse = ", ")
    )
  }

  note <- if (is.na(record[["Note"]])) "" else gsub("[\r\n]+", " ",
                                                    record[["Note"]])
  header <- c(
    signature,
    sprintf("%%%s  Extracted from %s.", toupper(name), script),
    "%",
    strwrap(note, width = 76, prefix = "%   "),
    "%",
    sprintf("%%   Source: %s lines %s.", script, spec),
    "%",
    "%   The block between the VERBATIM markers below is reproduced byte for",
    "%   byte from the reference implementation of Andreou, Chen and Bollt,",
    "%   Copyright (c) 2025 Marios Andreou, under the MIT Licence; see",
    "%   matlab_reference/LICENSE. Nothing in that block was retyped, and",
    "%   tools/extract.R re-checks it against the source on every run.",
    "%",
    "%   GENERATED FILE -- edit manifest/extracts.dcf, not this.",
    ""
  )

  lines <- c(
    header,
    .begin_marker(script, spec),
    body,
    .end_marker(script, spec),
    "",
    "end"
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, sprintf("%s.m", name))
  # Written with the reference's own terminator so the verbatim block compares
  # byte-for-byte without normalisation.
  .write_lines_preserving(lines, path, source_file$eol, TRUE)
  path
}

#' Generate every declared extraction.
#'
#' @param manifest_path Path to `manifest/extracts.dcf`.
#' @param reference_dir Directory holding the reference scripts.
#' @param out_dir Directory to write generated functions into.
#'
#' @returns A data frame of generated functions.
extract_all <- function(manifest_path,
                        reference_dir = "matlab_reference",
                        out_dir = file.path("oracle", "parity", "matlab",
                                            "extracted")) {
  records <- read_manifest(manifest_path)
  out <- vector("list", nrow(records))
  for (i in seq_len(nrow(records))) {
    record <- records[i, ]
    path <- extract_function(record, reference_dir, out_dir)
    out[[i]] <- data.frame(
      fun = record[["Function"]], script = record[["Script"]],
      lines = record[["Lines"]], path = path, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# -- verification --------------------------------------------------------------

#' Check every generated function still matches its source slice byte for byte.
#'
#' @param manifest_path Path to `manifest/extracts.dcf`.
#' @param reference_dir Directory holding the reference scripts.
#' @param out_dir Directory holding the generated functions.
#'
#' @returns A data frame with a `verdict` column; errors if any file drifted.
check_extraction <- function(manifest_path,
                             reference_dir = "matlab_reference",
                             out_dir = file.path("oracle", "parity", "matlab",
                                                 "extracted")) {
  records <- read_manifest(manifest_path)
  rows <- vector("list", nrow(records))

  for (i in seq_len(nrow(records))) {
    name <- records[i, "Function"]
    script <- records[i, "Script"]
    spec <- records[i, "Lines"]
    path <- file.path(out_dir, sprintf("%s.m", name))

    if (!file.exists(path)) {
      stop(sprintf("%s has not been generated.", name), call. = FALSE)
    }
    generated <- .read_lines_preserving(path)$lines
    begin <- which(generated == .begin_marker(script, spec))
    end <- which(generated == .end_marker(script, spec))
    if (length(begin) != 1L || length(end) != 1L || end <= begin) {
      stop(sprintf("%s: verbatim markers missing or malformed.", name),
           call. = FALSE)
    }
    block <- generated[seq.int(begin + 1L, end - 1L)]

    source_lines <- .read_lines_preserving(file.path(reference_dir, script))$lines
    expected <- source_lines[.slice_lines(spec)]

    identical_block <- identical(block, expected)
    # Nothing outside the block may contain executable reference algebra: the
    # only lines permitted there are the signature, comments, blanks and `end`.
    outside <- generated[-seq.int(begin, end)]
    stray <- outside[!grepl("^\\s*(%|$)", outside) &
                       !grepl("^function\\b", outside) &
                       !grepl("^end\\s*$", outside)]

    rows[[i]] <- data.frame(
      fun = name, source = sprintf("%s:%s", script, spec),
      lines_verbatim = length(block),
      verdict = if (!identical_block) {
        "DRIFTED"
      } else if (length(stray) > 0L) {
        "STRAY CODE"
      } else {
        "ok"
      },
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  bad <- result[result$verdict != "ok", , drop = FALSE]
  if (nrow(bad) > 0L) {
    stop(
      sprintf(
        "extraction check failed for: %s",
        paste(sprintf("%s (%s)", bad$fun, bad$verdict), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  result
}
