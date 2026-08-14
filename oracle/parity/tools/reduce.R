# -- minimal-substitution reduction of a reference script ----------------------
#
# The reference scripts fix their problem size in top-level literal assignments
# and then run to completion, allocating objects that do not fit on a
# workstation at the published settings. To obtain a gold workspace we must run
# them smaller, and the honest way to do that is a substitution small enough to
# be proved harmless rather than asserted to be.
#
# `reduce_script()` replaces declared literal assignments and then checks its
# own work: exactly the declared lines changed, each of them a bare assignment
# on both sides, and the file still the same length so that the line ranges the
# extraction manifest depends on still address the same code.

# -- byte-preserving line access -----------------------------------------------
#
# The reference scripts are CRLF-terminated with an unterminated final line.
# Reading with `readLines()` and writing with `writeLines()` silently rewrites
# every terminator, which turns a three-line substitution into a whole-file
# diff and destroys the reviewer's ability to see what actually changed. The
# terminator and the trailing-newline state are therefore carried through.

.read_lines_preserving <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  bytes <- readBin(connection, "raw", file.info(path)$size)
  text <- rawToChar(bytes)
  Encoding(text) <- "UTF-8"
  eol <- if (grepl("\r\n", text, fixed = TRUE)) "\r\n" else "\n"
  trailing <- endsWith(text, eol)
  if (trailing) {
    text <- substr(text, 1L, nchar(text) - nchar(eol))
  }
  list(
    lines = strsplit(text, eol, fixed = TRUE)[[1L]],
    eol = eol,
    trailing = trailing
  )
}

.write_lines_preserving <- function(lines, path, eol, trailing) {
  text <- paste(lines, collapse = eol)
  if (trailing) {
    text <- paste0(text, eol)
  }
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  writeBin(charToRaw(enc2utf8(text)), connection)
}

# -- helpers -------------------------------------------------------------------

# Split a `name=value; name=value` knob field into a named character vector.
.parse_knobs <- function(spec) {
  parts <- trimws(strsplit(spec, ";", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  keys <- sub("=.*$", "", parts)
  values <- sub("^[^=]*=", "", parts)
  stats::setNames(trimws(values), trimws(keys))
}

# A bare top-level assignment of a numeric literal, and nothing else.
.assignment_pattern <- function(name) {
  sprintf("^%s\\s*=\\s*[-+0-9.eE]+\\s*;\\s*$", name)
}

# -- reduction -----------------------------------------------------------------

#' Apply the declared knobs to a reference script.
#'
#' @param source_path Path to the unmodified reference script.
#' @param knobs Named character vector of replacement literals.
#' @param target_path Where to write the reduced script.
#'
#' @returns A data frame of the substitutions actually made, one row per knob.
reduce_script <- function(source_path, knobs, target_path) {
  source_file <- .read_lines_preserving(source_path)
  original <- source_file$lines
  reduced <- original
  hits <- integer(length(knobs))

  for (i in seq_along(knobs)) {
    name <- names(knobs)[i]
    matched <- grep(.assignment_pattern(name), original)
    if (length(matched) != 1L) {
      stop(
        sprintf(
          paste0(
            "knob `%s` must match exactly one bare assignment in %s; matched ",
            "%d line(s). A knob that matches zero lines is a stale manifest; ",
            "a knob that matches several cannot be applied unambiguously."
          ),
          name, basename(source_path), length(matched)
        ),
        call. = FALSE
      )
    }
    hits[i] <- matched
    reduced[matched] <- sprintf("%s = %s;", name, knobs[[i]])
  }

  # ---- self-check: the substitution reached nothing but the declared lines ----
  if (length(reduced) != length(original)) {
    stop("reduction changed the line count; line ranges would shift.",
         call. = FALSE)
  }
  changed <- which(reduced != original)
  if (!identical(sort(changed), sort(hits))) {
    stop(
      sprintf(
        "reduction changed lines %s but declared %s.",
        paste(sort(changed), collapse = ", "), paste(sort(hits), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  for (line in changed) {
    ok_before <- grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*=\\s*[-+0-9.eE]+\\s*;\\s*$",
                       original[line])
    ok_after <- grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*=\\s*[-+0-9.eE]+\\s*;\\s*$",
                      reduced[line])
    if (!ok_before || !ok_after) {
      stop(
        sprintf("line %d is not a bare literal assignment on both sides.", line),
        call. = FALSE
      )
    }
  }

  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  .write_lines_preserving(
    reduced, target_path, source_file$eol, source_file$trailing
  )

  data.frame(
    knob = names(knobs),
    line = hits,
    before = original[hits],
    after = reduced[hits],
    stringsAsFactors = FALSE
  )
}

#' Read the knob manifest and reduce every declared (script, profile).
#'
#' @param manifest_path Path to `manifest/knobs.dcf`.
#' @param reference_dir Directory holding the unmodified reference scripts.
#' @param out_dir Directory to write reduced scripts into, one per profile.
#'
#' @returns A data frame of every substitution made.
reduce_all <- function(manifest_path, reference_dir, out_dir) {
  records <- read_manifest(manifest_path)
  made <- list()
  for (i in seq_len(nrow(records))) {
    script <- records[i, "Script"]
    profile <- records[i, "Profile"]
    knobs <- .parse_knobs(records[i, "Knobs"])
    target <- file.path(out_dir, profile, script)
    result <- reduce_script(file.path(reference_dir, script), knobs, target)
    result$script <- script
    result$profile <- profile
    made[[length(made) + 1L]] <- result
  }
  do.call(rbind, made)
}
