# -- manifest reading ----------------------------------------------------------
#
# The manifests are Debian-control files: a flat, diffable, dependency-free
# key/value format that both R and a ten-line MATLAB parser can read. R's
# `read.dcf()` has no comment convention, so comments are stripped here before
# parsing. They earn their place -- a manifest a reviewer cannot read cold is a
# manifest that will drift from what it describes.

#' Read a DCF manifest, stripping `#` comment lines.
#'
#' @param path Path to the manifest.
#'
#' @returns A character matrix, one row per record.
read_manifest <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  connection <- textConnection(lines)
  on.exit(close(connection))
  read.dcf(connection)
}

#' Fetch one field from a manifest row, erroring when it is absent.
#'
#' @param records A manifest matrix.
#' @param i Integer scalar. Row index.
#' @param field Character scalar. Field name.
#'
#' @returns A character scalar.
manifest_field <- function(records, i, field) {
  if (!field %in% colnames(records)) {
    stop(sprintf("manifest has no `%s` field.", field), call. = FALSE)
  }
  value <- records[i, field]
  if (is.na(value)) {
    stop(sprintf("record %d has no `%s`.", i, field), call. = FALSE)
  }
  value
}
