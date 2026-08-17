# Fixture provenance gate ------------------------------------------------------
#
# The MATLAB oracle harnesses in this directory write their fixtures with bare
# filenames, so the files land wherever the script was run from. The copies the
# package actually ships live in `aciR/inst/extdata/`. Two directories therefore
# hold what is meant to be one set of bytes, and nothing reads the copies here,
# which is precisely the arrangement in which they can disagree for months
# without anyone noticing.
#
# `aciR/tests/testthat/test-oracle-manifest.R` pins the shipped copies against
# the hashes recorded in the manifest, so a change to `inst/extdata/` is caught.
# What that test cannot see is this directory: `tools/` is excluded from the
# build, so under `R CMD check` the path does not exist. A testthat test
# reaching out here would skip on every continuous-integration run, and the
# manifest test's own header records why that is not acceptable: a validation
# gate that skips is not a gate.
#
# So the gate lives here instead, as a script run from the repository root where
# both directories are present, and is wired into the lint workflow. Run it
# after regenerating any fixture:
#
#   Rscript tools/oracle/check_fixture_provenance.R
#
# It exits 0 when the two directories agree and 1 otherwise, so continuous
# integration fails on drift rather than reporting it into a log nobody reads.

check_fixture_provenance <- function(
    oracle_dir = file.path("tools", "oracle"),
    shipped_dir = file.path("aciR", "inst", "extdata")) {
  for (d in c(oracle_dir, shipped_dir)) {
    if (!dir.exists(d)) {
      stop("directory not found: ", d,
           ". Run this from the repository root.", call. = FALSE)
    }
  }

  # Only the top level of the oracle directory: `parity/` below it holds the
  # side-by-side harness, whose outputs are a different artefact and are not
  # shipped.
  oracle <- list.files(oracle_dir, pattern = "\\.csv$", full.names = FALSE)
  shipped <- list.files(shipped_dir, pattern = "\\.csv$", full.names = FALSE)
  common <- intersect(oracle, shipped)

  if (length(common) == 0L) {
    stop("no fixture appears in both directories, which means one of the two ",
         "paths is wrong rather than that everything agrees.", call. = FALSE)
  }

  same <- vapply(common, function(f) {
    a <- unname(tools::md5sum(file.path(oracle_dir, f)))
    b <- unname(tools::md5sum(file.path(shipped_dir, f)))
    identical(a, b)
  }, logical(1L))

  # A fixture present in the oracle directory but absent from the shipped one is
  # the failure this gate exists to catch: a harness was run and its output was
  # never copied across. The reverse is unremarkable, because several shipped
  # fixtures are derived rather than written directly by a harness.
  orphaned <- setdiff(oracle, shipped)

  list(
    common = common,
    same = same,
    differing = common[!same],
    orphaned = orphaned,
    ok = all(same) && length(orphaned) == 0L
  )
}

report_fixture_provenance <- function(res) {
  cat(sprintf("Fixtures in both directories : %d\n", length(res$common)))
  cat(sprintf("Byte-identical               : %d\n", sum(res$same)))
  cat(sprintf("Differing                    : %d\n", length(res$differing)))
  cat(sprintf("Written but never shipped    : %d\n", length(res$orphaned)))

  if (length(res$differing) > 0L) {
    cat("\nThese differ between tools/oracle/ and aciR/inst/extdata/:\n")
    cat(paste0("  ", res$differing, collapse = "\n"), "\n")
    cat("\nThe shipped copy is the authority. If the harness output is the\n",
        "correct one, copy it across and refresh the hashes in\n",
        "aciR/inst/extdata/oracle-manifest.yml.\n", sep = "")
  }
  if (length(res$orphaned) > 0L) {
    cat("\nWritten by a harness but not shipped:\n")
    cat(paste0("  ", res$orphaned, collapse = "\n"), "\n")
  }
  invisible(res)
}

if (sys.nframe() == 0L) {
  res <- check_fixture_provenance()
  report_fixture_provenance(res)
  if (!res$ok) {
    cat("\nFAIL: the two fixture directories disagree.\n")
    quit(status = 1L)
  }
  cat("\nOK: every shared fixture is byte-identical.\n")
}
