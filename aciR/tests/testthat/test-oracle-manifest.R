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
# not the property wanted here (detecting accidental drift in a committed
# fixture is), and the two hashes are asserted to describe the same file by
# the refresh procedure recorded in the manifest.

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

# Read a block of `key: value` numbers from the manifest, e.g. the recorded
# errors under observed_max_abs_error.
#
# The lines are read raw, not trimmed: indentation is what distinguishes a
# member of the block from its sibling at the parent level, and trimming it
# away made this reader swallow `measured_on: 2026-07-15` from the enclosing
# block, a date that a permissive number pattern then accepted, because a
# hyphen is also a minus sign. The pattern below matches a number in decimal
# or scientific form and nothing else.
.aci_manifest_numbers <- function(block, indent = 4L) {
  lines <- readLines(.aci_manifest_path(), warn = FALSE)
  start <- grep(sprintf("^\\s*%s:\\s*$", block), lines)
  testthat::expect_length(start, 1L)

  pattern <- sprintf(
    "^ {%d}([a-z_]+): (-?[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?)\\s*$", indent
  )
  out <- numeric(0L)
  for (i in seq(start + 1L, length(lines))) {
    hit <- regmatches(lines[i], regexec(pattern, lines[i]))[[1L]]
    if (length(hit) < 3L) {
      break
    }
    out[hit[2L]] <- as.numeric(hit[3L])
  }
  testthat::expect_false(anyNA(out))
  out
}

# Recompute an oracle's maximum absolute error from the shipped fixtures, the
# same way the oracle tests do.
.aci_oracle_error <- function(id) {
  signal <- read.csv(
    system.file("extdata", sprintf("%s_signal_x.csv", id), package = "aciR"),
    header = FALSE
  )
  x <- signal$V2
  ref <- read.csv(
    system.file("extdata", sprintf("%s_reference.csv", id), package = "aciR")
  )

  comp <- if (id == "dyad") {
    aci_dyad_components(x, aci_dyad_model()$parameters)
  } else {
    list(
      L_x = 2 * x, f_x = 0.5 - 0.5 * x, L_y = -0.5, f_y = 1 - 2 * x^2,
      S_xoS_x = 0.6^2 + 0.3^2, S_yoS_y = 0.5^2 + 0.8^2,
      S_yoS_x = 0.5 * 0.6 + 0.8 * 0.3, S_xoS_y = 0.5 * 0.6 + 0.8 * 0.3
    )
  }
  filt <- aci_filter(x, comp, dt = 0.001, mu0 = 2, R0 = 0.1)
  smooth <- aci_smoother(x, comp, dt = 0.001, filt)
  metric <- aci_metric(filt, smooth)
  idx <- seq(1, length(x), by = 100)

  max(
    abs(filt$mean[idx] - ref$filter_mean),
    abs(filt$cov[idx] - ref$filter_cov),
    abs(smooth$mean[idx] - ref$smoother_mean),
    abs(smooth$cov[idx] - ref$smoother_cov),
    abs(metric[idx] - ref$ACI_metric)
  )
}

test_that("the manifest records both hashes for every shipped fixture", {
  # Enumerated rather than discovered, so that adding a fixture is a deliberate
  # act that updates this list. The causal-influence-range entries carry no
  # signal of their own: that harness reads the dyad signal already pinned
  # here, which is why only reference files appear for it.
  expected <- c(
    "dyad_signal_x.csv", "dyad_reference.csv",
    "cross_signal_x.csv", "cross_reference.csv",
    "cir_online_reference.csv", "cir_range_reference.csv",
    "predprey_signal.csv",
    "predprey_reference_predator_to_prey.csv",
    "predprey_reference_prey_to_predator.csv",
    "mv_signal.csv", "mv_reference.csv",
    "cir_cross_online_reference.csv", "cir_cross_range_reference.csv",
    "enso_signal.csv", "enso_reference.csv",
    "mv_online_reference.csv"
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

test_that("the manifest's recorded errors are the errors actually observed", {
  # The manifest asserts a measured quantity, the agreement between this
  # package and its oracles. A number like that means something only if it
  # came from a run, and a hand-typed one is indistinguishable from a
  # remembered or invented one until something checks it. This test is that
  # something: it recomputes both errors and holds the manifest to them.
  #
  # The tolerance is a factor of ten, which is loose against the floating-point
  # variation a different platform or BLAS produces (observed: well under 2x)
  # and tight against the failure it exists to catch, a plausible-looking
  # value that was never measured.
  recorded <- .aci_manifest_numbers("observed_max_abs_error")
  expect_setequal(names(recorded), c("dyad", "cross"))

  measured <- c(
    dyad = .aci_oracle_error("dyad"),
    cross = .aci_oracle_error("cross")
  )

  for (id in names(recorded)) {
    expect_lt(measured[[id]], 1e-6)
    ratio <- measured[[id]] / recorded[[id]]
    expect_true(
      ratio > 0.1 && ratio < 10,
      info = sprintf(
        paste0(
          "oracle-manifest.yml records observed_max_abs_error for `%s` as %g, ",
          "but the measured error is %g. Either the numerics moved, or the ",
          "recorded value was not produced by a run. Record what the run says."
        ),
        id, recorded[[id]], measured[[id]]
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
