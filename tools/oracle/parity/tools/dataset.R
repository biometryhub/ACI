# the dataset contract ---------------------------------------------------------
#
# What both implementations have to agree on before they can be compared: the
# observed signal, the conditional Gaussian components, the step, and the
# initial condition. The format is one DCF of scalars and one CSV of columns,
# chosen because R reads both with `read.dcf`/`read.csv` and MATLAB with a
# short parser and `readtable`, so neither side gains a package dependency for
# a development harness.
#
# One asymmetry is baked into the contract and is stated here rather than
# hidden. The reference parameterises noise by the feedback matrices Sx_1,
# Sx_2, Sy_1, Sy_2 and forms the Grammians itself; aciR takes the Grammians
# directly. A dataset therefore carries the feedbacks, and the R side derives
# what aciR wants. The two are equivalent here, but only the reference's
# direction is lossy: not every admissible Grammian pair is expressible as two
# feedback columns of this shape, so aciR's surface is the wider one.

#' Write a dataset bundle.
#'
#' @param path Directory to create.
#' @param meta Named list of scalars and strings for `meta.dcf`.
#' @param arrays Data frame of per-step columns, `n + 1` rows.
#'
#' @returns The directory path, invisibly.
write_dataset <- function(path, meta, arrays) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  fields <- vapply(meta, function(v) as.character(v), character(1L))
  write.dcf(as.data.frame(as.list(fields), stringsAsFactors = FALSE),
            file.path(path, "meta.dcf"))
  utils::write.csv(arrays, file.path(path, "arrays.csv"), row.names = FALSE)
  invisible(path)
}

#' Read a dataset bundle.
#'
#' @param path Directory holding `meta.dcf` and `arrays.csv`.
#'
#' @returns A list with `meta` (named list) and `arrays` (data frame).
read_dataset <- function(path) {
  raw <- read_manifest(file.path(path, "meta.dcf"))
  meta <- as.list(raw[1L, ])
  numeric_fields <- c("N", "dt", "mu0", "R0", "L_y", "Sx_1", "Sx_2", "Sy_1",
                      "Sy_2")
  for (field in intersect(numeric_fields, names(meta))) {
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
  list(meta = meta, arrays = arrays)
}

#' Build the aciR components list from a dataset.
#'
#' @param dataset A list from [read_dataset()].
#'
#' @returns A conditional Gaussian components list.
dataset_components <- function(dataset) {
  meta <- dataset$meta
  arrays <- dataset$arrays
  s_yox <- meta$Sy_1 * meta$Sx_1 + meta$Sy_2 * meta$Sx_2
  list(
    L_x = arrays$L_x,
    L_y = meta$L_y,
    f_x = arrays$f_x,
    f_y = arrays$f_y,
    S_xoS_x = meta$Sx_1^2 + meta$Sx_2^2,
    S_yoS_y = meta$Sy_1^2 + meta$Sy_2^2,
    S_yoS_x = s_yox,
    S_xoS_y = s_yox
  )
}
