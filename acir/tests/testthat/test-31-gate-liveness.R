# Induced violations for the gates that had no recorded failing run -----------
#
# inst/evidence/gate_liveness.md is the register: one row per gate, where it
# is enforced, the induced violation that proves it fails when it should, and
# where that run is recorded.  Most of the package's gates are already proved
# live by a test that induces the violation as part of its own subject: the
# strict covariance policy has its three firing probes, the trusted-filter
# token has its mutation and forgery blocks, the contract validators have
# their planted-violation blocks, and the conditional split has its rejected
# specifications.  This file exists for the gates that had no such run.
#
# The blocks below induce a violation and assert the gate reports it.  They
# are deliberately cheap: one flipped byte, one tampered copy of a register,
# one unstable model, one refused argument, and one perturbation of an oracle
# column.  None of them touch the shipped fixtures or the shipped register;
# every tamper is made on a copy, in memory or under tempdir().

.gl_oracle_dir <- testthat::test_path("fixtures", "oracles")
.gl_oracle_file <- function(name) file.path(.gl_oracle_dir, name)

.gl_register_columns <- c("export", "feature", "check_method", "against",
                          "tolerance_class", "state", "fixture_path", "sha256")


# Gate: fixture byte-pinning ---------------------------------------------------

test_that("byte-pinning fails on one flipped byte of a fixture copy", {
  src <- .gl_oracle_file("dyad_reference.csv")
  pinned <- .oracle_manifest_hashes(
    .gl_oracle_file("oracle-manifest.yml"), "md5"
  )[["dyad_reference.csv"]]

  # The gate passes on the shipped bytes.  Without this line the induced
  # failure below would also be satisfied by a manifest that pins nothing.
  expect_identical(unname(tools::md5sum(src)), pinned)

  tampered <- file.path(tempdir(), "gl-dyad_reference.csv")
  on.exit(unlink(tampered), add = TRUE)
  bytes <- readBin(src, "raw", file.size(src))
  # One decimal digit in the interior of the file, moved by one.  The file
  # length does not change, so nothing but the hash can see it.
  j <- which(bytes >= as.raw(0x30) & bytes <= as.raw(0x38))[100L]
  bytes[j] <- as.raw(as.integer(bytes[j]) + 1L)
  writeBin(bytes, tampered)

  expect_identical(file.size(tampered), file.size(src))
  expect_false(identical(unname(tools::md5sum(tampered)), pinned))
  sha256 <- tryCatch(get("sha256sum", envir = asNamespace("tools")),
                     error = function(e) NULL)
  if (is.function(sha256))
    expect_false(identical(unname(sha256(tampered)), unname(sha256(src))))
  else
    skip("tools::sha256sum() is not available in this R; md5 tamper check only")
})


test_that("the shipped manifest reader sees the names the narrow one went blind to", {
  # The induced violation here is the historical reader, applied to today's
  # manifest.  `^[A-Za-z_]+\\.csv:$` excludes digits, so it returns none of the
  # six `enso6_partition_*` entries; the shipped scan returns all six.  A
  # future name class the strict reader cannot see fails the same way.
  path <- .gl_oracle_file("oracle-manifest-partitions.yml")
  lines <- trimws(readLines(path, warn = FALSE))
  narrow <- sub(":$", "", lines[grep("^[A-Za-z_]+\\.csv:$", lines)])
  expect_length(grep("^enso6_partition_", narrow, value = TRUE), 0L)

  shipped <- .oracle_manifest_declared_csvs(path)
  expect_length(grep("^enso6_partition_", shipped, value = TRUE), 6L)
  expect_length(shipped, 14L)
})


# Gate: the evidence register's coverage test ---------------------------------

test_that("register coverage fails on a dropped row, a moved hash and a renamed column", {
  path <- system.file("evidence", "register.csv", package = "acir")
  expect_true(nzchar(path))
  reg <- utils::read.csv(path, stringsAsFactors = FALSE,
                         colClasses = "character")
  exports <- getNamespaceExports("acir")

  # (a) coverage.  The shipped register names every export; a register with
  # one verb's rows removed does not.
  expect_identical(setdiff(exports, reg$export), character(0L))
  dropped <- reg[reg$export != "safe_chol", , drop = FALSE]
  expect_identical(setdiff(exports, dropped$export), "safe_chol")

  # (b) the fixture hash.  The shipped sha256 matches the file; a hash with one
  # hex digit moved does not, and a checked row naming a fixture that is not in
  # the tree fails on existence before the hash is ever computed.
  row <- reg[reg$state == "checked", , drop = FALSE][1L, ]
  fixture <- .gl_oracle_file(basename(row$fixture_path))
  expect_true(file.exists(fixture))
  sha256 <- tryCatch(get("sha256sum", envir = asNamespace("tools")),
                     error = function(e) NULL)
  if (is.function(sha256)) {
    expect_identical(unname(sha256(fixture)), row$sha256)
    moved <- sub("^.", if (startsWith(row$sha256, "a")) "b" else "a", row$sha256)
    expect_false(identical(unname(sha256(fixture)), moved))
  }
  expect_false(file.exists(.gl_oracle_file("no_such_reference.csv")))

  # (c) the column set, induced through a file so the reader is exercised and
  # not only the comparison.
  expect_identical(names(reg), .gl_register_columns)
  renamed <- file.path(tempdir(), "gl-register.csv")
  on.exit(unlink(renamed), add = TRUE)
  out <- reg
  names(out)[3L] <- "method"
  utils::write.csv(out, renamed, row.names = FALSE)
  back <- utils::read.csv(renamed, stringsAsFactors = FALSE,
                          colClasses = "character")
  expect_false(identical(names(back), .gl_register_columns))
  if (!is.function(sha256))
    skip("tools::sha256sum() is not available in this R; hash gates not exercised")
})


# Gate: the oracle tolerance ---------------------------------------------------

test_that("the dyad oracle tolerance rejects a perturbation in every graded column", {
  # Non-vacuity, induced column by column.  A grade that silently stopped
  # comparing one of its five series would still pass its own test; it cannot
  # pass this one, because each column is moved on its own and the comparison
  # is required to see it.  The perturbation is 1e-5, ten times the 1e-6
  # tolerance and far below any real regression.
  signal <- utils::read.csv(.gl_oracle_file("dyad_signal_x.csv"),
                            header = FALSE)
  ref <- utils::read.csv(.gl_oracle_file("dyad_reference.csv"))
  obs <- observed_trajectory(signal$V1, matrix(signal$V2, ncol = 1L))
  bundle <- .compile_cgns_run(aci_dyad_model(), obs)
  filter <- .cgns_filter_compiled(
    bundle, init = list(mean = 2, cov = matrix(0.1, 1L, 1L)),
    stepper = "explicit", nsub = 1L
  )
  smoother <- .cgns_smoother_compiled(bundle, filter)
  metric <- .gaussian_kl_path_compiled(bundle, smoother, filter,
                                       decompose = TRUE)
  idx <- seq.int(1L, nrow(signal), by = 100L)

  graded <- list(
    filter_mean    = filter$mean[idx, 1L],
    filter_cov     = filter$cov[1L, 1L, idx],
    smoother_mean  = smoother$mean[idx, 1L],
    smoother_cov   = smoother$cov[1L, 1L, idx],
    ACI_metric     = metric$total[idx]
  )
  error <- function(reference)
    max(vapply(names(graded),
               function(nm) max(abs(graded[[nm]] - reference[[nm]])),
               numeric(1L)))

  expect_lt(error(ref), 1e-6)
  for (nm in names(graded)) {
    perturbed <- ref
    perturbed[[nm]] <- ref[[nm]] + 1e-5
    expect_gt(error(perturbed), 1e-6, label = nm)
  }
})


# Gate: the simulator divergence guard -----------------------------------------

test_that("the divergence guard stops an unstable simulation and names the step", {
  m <- aci_model(
    Lx  = function(t, x) matrix(0, 1L, 1L),
    fx  = function(t, x) 0 * x,
    Ly  = function(t, x) matrix(1e6, 1L, 1L),
    fy  = function(t, x) 0,
    Sx1 = function(t, x) matrix(0.5, 1L, 1L),
    Sy2 = function(t, x) matrix(1, 1L, 1L),
    k = 1L, l = 1L, name = "gate-liveness-divergent"
  )
  # The hidden self-drift is 1e6, so on a 1e-3 grid the Euler factor is 1001
  # per step and the state leaves the double range inside 110 steps.
  err <- tryCatch(stats::simulate(m, seed = 1L, T = 0.2, dt = 1e-3,
                                  burn_in = 0),
                  error = function(e) e)
  expect_s3_class(err, "aci_error_sim_divergence")
  expect_true(is.integer(err$step) || is.numeric(err$step))
  expect_gt(err$step, 1L)
  expect_match(conditionMessage(err), "diverged at step")

  # The same model on a grid fine enough to stay finite returns normally, so
  # the abort is the guard firing and not the constructor refusing the model.
  ok <- stats::simulate(m, seed = 1L, T = 2e-4, dt = 1e-6, burn_in = 0)
  expect_true(all(is.finite(ok$obs$x)))
})


# Gate: the non-CGNS model routes ----------------------------------------------

test_that("a non-CGNS model is refused by every assimilation verb", {
  # The ensemble engine is out of scope in this release, and each entry point
  # says so with a classed condition rather than dispatching to a CGNS route.
  # aci_online()'s refusal is induced in test-27-online-lag.R; the other three
  # had no induced run before this block.
  m <- structure(list(k = 1L, l = 1L), class = "stochastic_model")
  ob <- observed_trajectory(seq(0, 0.1, by = 0.01),
                            matrix(seq(0, 0.1, by = 0.01), ncol = 1L))
  expect_error(aci_filter(m, ob), class = "aci_error_not_implemented")
  expect_error(aci_smoother(m, ob), class = "aci_error_not_implemented")
  expect_error(aci_online(m, ob, lag = 1), class = "aci_error_not_implemented")
  expect_error(lag_table(m, ob, mode = "forward"),
               class = "aci_error_not_implemented")
})


# Gate: the aci_range direction gate -------------------------------------------

test_that("the backward range is refused on both aci_range methods", {
  m <- aci_dyad_model()
  s <- stats::simulate(m, seed = 1L, T = 0.3, dt = 0.01, burn_in = 0)
  ob <- as_obs(s)
  init <- list(mean = 2, cov = matrix(0.1, 1L, 1L))
  table <- lag_table(m, ob, mode = "forward", init = init)
  result <- aci(m, ob, keep = "table", init = init)

  expect_error(aci_range(table, direction = "backward"),
               class = "aci_error_not_implemented")
  expect_error(aci_range(result, direction = "backward"),
               class = "aci_error_not_implemented")
  expect_match(
    conditionMessage(tryCatch(aci_range(table, direction = "backward"),
                              error = function(e) e)),
    "not in this release"
  )

  # The gate is the direction, not the verb: forward still returns, and says
  # which direction it returned.
  expect_identical(suppressWarnings(aci_range(table))$direction, "forward")
  expect_identical(suppressWarnings(aci_range(result))$direction, "forward")
  expect_error(aci_range(table, direction = "sideways"))
})
