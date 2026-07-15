# The fixtures are the package's evidence, so their identity is part of the
# contract. This test replaces an earlier one that reached outside the package
# for the same files and skipped when it could not find them: a validation gate
# that skips is not a gate. Everything here runs from the installed package,
# needs no package beyond base R, and never skips.
#
# The manifest carries two hashes per fixture. The SHA-256 is the audit trail,
# reproducible by anyone with `shasum -a 256`; the MD5 is what this gate checks,
# because `tools::md5sum()` is base R at every version the package supports
# whereas `tools::sha256sum()` arrived only in R 4.5.0. Collision resistance is
# not the property wanted here -- detecting accidental drift in a committed
# fixture is -- and the two hashes are asserted to describe the same file by the
# refresh procedure recorded in the manifest.

.aci_manifest_path <- function() {
  path <- system.file("extdata", "oracle-manifest.yml", package = "aciR")
  testthat::expect_true(
    file.exists(path) && nzchar(path),
    info = "the oracle manifest must ship in inst/extdata"
  )
  path
}

# Rather than take a YAML dependency for eight scalars, read the hash lines
# directly: the manifest is the authority. A parser that silently matched
# nothing would make this test vacuous, so the count is asserted separately.
.aci_manifest_hashes <- function(algorithm = "md5") {
  lines <- trimws(readLines(.aci_manifest_path(), warn = FALSE))
  entries <- grep("^[A-Za-z_]+\\.csv:$", lines)
  hashes <- character(0L)
  for (i in entries) {
    block <- lines[seq(i + 1L, min(i + 3L, length(lines)))]
    hit <- grep(sprintf("^%s: [0-9a-f]+$", algorithm), block, value = TRUE)
    testthat::expect_length(hit, 1L)
    hashes[sub(":$", "", lines[i])] <- sub(
      sprintf("^%s: ", algorithm), "", hit[1L]
    )
  }
  hashes
}

test_that("the manifest records both hashes for every shipped fixture", {
  expected <- c(
    "dyad_signal_x.csv", "dyad_reference.csv",
    "cross_signal_x.csv", "cross_reference.csv"
  )
  md5 <- .aci_manifest_hashes("md5")
  sha <- .aci_manifest_hashes("sha256")

  expect_setequal(names(md5), expected)
  expect_setequal(names(sha), expected)
  expect_true(all(nchar(md5) == 32L))
  expect_true(all(nchar(sha) == 64L))
})

test_that("every shipped fixture matches its recorded hash", {
  hashes <- .aci_manifest_hashes("md5")
  for (name in names(hashes)) {
    path <- system.file("extdata", name, package = "aciR")
    expect_true(file.exists(path) && nzchar(path))
    expect_identical(
      unname(tools::md5sum(path)),
      unname(hashes[[name]]),
      info = sprintf(
        paste0(
          "%s does not match the hash recorded in oracle-manifest.yml. A ",
          "fixture must not change without a deliberate refresh; see the ",
          "refresh policy in the manifest."
        ),
        name
      )
    )
  }
})

test_that("the manifest states what each oracle does and does not grade", {
  # The scope of an oracle is itself a claim. The dyad fixture pins the noise
  # cross-covariance at zero, and a reader who took it as validating the whole
  # core would be wrong; the manifest must say so in as many words.
  text <- paste(readLines(.aci_manifest_path(), warn = FALSE), collapse = " ")
  expect_match(text, "does_not_grade")
  expect_match(text, "cross-covariance")
  expect_match(text, "analytic_oracles")
})
