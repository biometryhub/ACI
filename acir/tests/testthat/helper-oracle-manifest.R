# One manifest reader for both oracle manifests.
#
# fixtures/oracles carries two provenance manifests, not one:
#
#   oracle-manifest.yml             the authors-source and independent-MATLAB
#                                   fixtures migrated from aciR;
#   oracle-manifest-partitions.yml  the source-derived scalar ENSO partitions.
#
# They are kept apart because they describe different evidence classes and
# must never be cited as one another, and each is read by its own test file
# (test-19 and test-28). The parsing lives here so both files use ONE reader:
# a per-file copy is how the two drifted apart in the first place. The
# original reader matched `^[A-Za-z_]+\\.csv:$`, which excludes digits. Every
# name in oracle-manifest.yml happens to be digit-free, so the omission was
# invisible until the `enso6_*` names arrived, at which point that reader
# would have skipped all six SILENTLY: no error, no failing expectation, just
# no byte pinning. The class below admits digits, and
# `.oracle_manifest_declared_csvs()` exists so that a future name class the
# strict reader cannot see fails loudly instead.

# Strict reader: the block headers this reader is willing to pin, and the hash
# it pulls from each block.
.oracle_manifest_hashes <- function(path, algorithm = "md5") {
  lines <- trimws(readLines(path, warn = FALSE))
  entries <- grep("^[A-Za-z0-9_]+\\.csv:$", lines)
  hashes <- character(0L)
  for (i in entries) {
    block <- lines[seq.int(i + 1L, min(i + 3L, length(lines)))]
    hit <- grep(sprintf("^%s: [0-9a-f]+$", algorithm), block, value = TRUE)
    expect_length(hit, 1L)
    hashes[sub(":$", "", lines[i])] <- sub(
      sprintf("^%s: ", algorithm), "", hit[1L]
    )
  }
  hashes
}

# Permissive scan: every line that looks like a `<something>.csv:` block
# header, whatever characters the name uses. Deliberately NOT the strict
# pattern -- the point is to disagree with it whenever the strict pattern is
# too narrow, so the tests can assert the two agree.
.oracle_manifest_declared_csvs <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  sub(":$", "", lines[grep("^[^#[:space:]]+\\.csv:$", lines)])
}
