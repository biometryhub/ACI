# The shipped evidence register, and the coverage it is required to have -------
#
# inst/evidence/register.csv is one row per checked feature of the public
# surface: what is checked, how, against what, at which tolerance class, and
# with which hash-pinned fixture behind it.  This file is the gate that keeps
# it honest.  Three things fail the build:
#
#   * an exported verb with no row at all, so that a new export cannot ship
#     without its evidence being stated;
#   * a row marked `checked` whose fixture is absent from the tree or whose
#     bytes no longer match the sha256 the row carries;
#   * a register that does not parse with the exact column set below.
#
# The `state` column separates fixture-graded rows from the rest, and nothing
# more.  `checked` means the feature is graded numerically against the named
# hash-pinned fixture.  `behavioural_only` means no such fixture stands behind
# that row; the `check_method` column still records what does, which for many
# of those rows is an exact relation rather than a behavioural assertion.  The
# two columns are read together or not at all.
#
# The register's sha256 values are the oracle manifests' own.  The last
# assertion below requires every cited fixture to be manifest-pinned, so a row
# cannot introduce a fixture that no manifest is responsible for.

.register_columns <- c("export", "feature", "check_method", "against",
                       "tolerance_class", "state", "fixture_path", "sha256")

.register_methods <- c("authors_source_comparison", "source_derived_comparison",
                       "independent_transcription", "exact_relation",
                       "second_implementation_agreement", "behavioural")

.register_tolerances <- c("exact", "machine_1e-12", "numerical_1e-6",
                          "numerical_1e-8", "numerical_1e-10", "none")

.register_states <- c("checked", "behavioural_only")

.register_prefix <- "tests/testthat/fixtures/oracles/"

.register_path <- function() {
  path <- system.file("evidence", "register.csv", package = "acir")
  expect_true(nzchar(path) && file.exists(path),
              info = "inst/evidence/register.csv is not installed")
  path
}

.register_read <- function() {
  utils::read.csv(.register_path(), stringsAsFactors = FALSE,
                  colClasses = "character")
}

# The register stores package-root-relative paths; a test runs from
# tests/testthat.  Resolve rather than re-spell, so the stored path stays the
# one a reader of the file can follow.
.register_resolve <- function(fixture_path) {
  expect_true(startsWith(fixture_path, .register_prefix),
              info = sprintf("fixture_path outside the oracle directory: %s",
                             fixture_path))
  testthat::test_path("fixtures", "oracles",
                      substring(fixture_path, nchar(.register_prefix) + 1L))
}


test_that("the evidence register parses with its exact column set", {
  reg <- .register_read()
  expect_identical(names(reg), .register_columns)
  expect_gt(nrow(reg), 0L)

  # A controlled vocabulary in every graded column.  A free-text state or
  # method would make the coverage assertions below unenforceable.
  expect_true(all(reg$check_method %in% .register_methods),
              info = paste(setdiff(reg$check_method, .register_methods),
                           collapse = ", "))
  expect_true(all(reg$tolerance_class %in% .register_tolerances),
              info = paste(setdiff(reg$tolerance_class, .register_tolerances),
                           collapse = ", "))
  expect_true(all(reg$state %in% .register_states))

  # No empty cells outside the two that are empty by construction.
  for (column in c("export", "feature", "check_method", "against",
                   "tolerance_class", "state"))
    expect_true(all(nzchar(reg[[column]])),
                info = sprintf("empty %s cell", column))

  # One row per feature, not one row repeated.
  expect_identical(anyDuplicated(paste(reg$export, reg$feature, sep = " | ")),
                   0L)
})


test_that("every exported verb appears in the evidence register", {
  reg <- .register_read()
  exports <- getNamespaceExports("acir")
  # S3 methods are registered, not exported, so this is the caller-visible
  # surface and nothing else.
  missing <- setdiff(exports, reg$export)
  expect_identical(
    missing, character(0L),
    info = sprintf("exported with no register row: %s",
                   paste(missing, collapse = ", "))
  )

  # The reverse direction catches a row left behind by a rename, which would
  # otherwise satisfy the assertion above while documenting nothing.
  stale <- setdiff(reg$export, exports)
  expect_identical(
    stale, character(0L),
    info = sprintf("register rows naming no export: %s",
                   paste(stale, collapse = ", "))
  )
})


test_that("every checked row names a fixture whose bytes still match", {
  sha256 <- tryCatch(get("sha256sum", envir = asNamespace("tools")),
                     error = function(e) NULL)
  reg <- .register_read()
  checked <- reg[reg$state == "checked", , drop = FALSE]
  expect_gt(nrow(checked), 0L)

  # Both cells are mandatory on a checked row.  A checked row with no fixture
  # would pass every assertion below by having nothing to check.
  expect_true(all(nzchar(checked$fixture_path)))
  expect_true(all(nchar(checked$sha256) == 64L))
  expect_true(all(grepl("^[0-9a-f]{64}$", checked$sha256)))

  for (i in seq_len(nrow(checked))) {
    path <- .register_resolve(checked$fixture_path[i])
    expect_true(
      file.exists(path),
      info = sprintf("register cites a fixture that is not in the tree: %s",
                     checked$fixture_path[i])
    )
    if (is.function(sha256))
      expect_identical(
        unname(sha256(path)), checked$sha256[i],
        info = sprintf("%s no longer matches the sha256 the register carries",
                       checked$fixture_path[i])
      )
  }
  if (!is.function(sha256))
    skip("tools::sha256sum() is not available in this R; existence only")
})


test_that("a behavioural-only row claims no fixture", {
  reg <- .register_read()
  other <- reg[reg$state != "checked", , drop = FALSE]
  expect_gt(nrow(other), 0L)
  expect_true(all(!nzchar(other$fixture_path)))
  expect_true(all(!nzchar(other$sha256)))
})


test_that("every fixture the register cites is pinned by a manifest", {
  reg <- .register_read()
  cited <- unique(reg$fixture_path[nzchar(reg$fixture_path)])
  pinned <- unlist(lapply(
    c("oracle-manifest.yml", "oracle-manifest-partitions.yml"),
    function(mf) .oracle_manifest_declared_csvs(
      testthat::test_path("fixtures", "oracles", mf)
    )
  ))
  expect_true(
    all(basename(cited) %in% pinned),
    info = sprintf("cited but unpinned: %s",
                   paste(setdiff(basename(cited), pinned), collapse = ", "))
  )
})
