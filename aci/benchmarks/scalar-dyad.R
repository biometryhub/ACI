#!/usr/bin/env Rscript

# Development benchmark for the production compiled-CGNS architecture.
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}
has_arg <- function(name) paste0("--", name) %in% args
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot determine benchmark script path.")
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             mustWork = TRUE)
package_dir <- normalizePath(
  arg_value("package", dirname(dirname(script_file))), mustWork = TRUE)
output_dir <- normalizePath(
  arg_value("output", file.path(dirname(script_file), "results",
                                 "production-compiled-cgns")),
  mustWork = FALSE)
allocation_stage <- arg_value("allocation-stage", NULL)
quick <- has_arg("quick")
no_alloc <- has_arg("no-alloc")
no_profile <- has_arg("no-profile")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

needed <- c("pkgload", "microbenchmark", "digest")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("Install development benchmark package(s): ",
       paste(missing, collapse = ", "))
pkgload::load_all(package_dir, attach = FALSE, export_all = TRUE,
                  helpers = FALSE, quiet = TRUE)
aci_ns <- asNamespace("aci")
aci_fun <- function(name) get(name, envir = aci_ns, inherits = FALSE)

# These are the same front door and kernels used by the public CGNS workflows.
compile_run <- aci_fun(".compile_cgns_run")
compile_generic <- aci_fun(".compile_cgns_complete")
filter_compiled <- aci_fun(".cgns_filter_compiled")
smoother_compiled <- aci_fun(".cgns_smoother_compiled")
metric_compiled <- aci_fun(".gaussian_kl_path_compiled")
complete_compiled <- aci_fun(".aci_cgns_compiled")
lag_compiled <- aci_fun(".lag_table_compiled")
cir_compiled <- aci_fun(".forward_cir_compiled")
model_dyad <- aci_fun("model_dyad")
simulate_model <- aci_fun("simulate.stochastic_model")
as_obs <- aci_fun("as_obs")
observed_trajectory <- aci_fun("observed_trajectory")
cgns_from_affine <- aci_fun("cgns_from_affine")
cgns_model <- aci_fun("cgns_model")
nontarget <- aci_fun("nontarget")
da_filter <- aci_fun("da_filter")
da_smooth <- aci_fun("da_smooth")
aci_public <- aci_fun("aci")
lag_public <- aci_fun("lag_table")
forward_cir <- aci_fun("forward_cir")

hash_double <- function(x) digest::digest(
  writeBin(as.double(x), raw(), size = 8L, endian = "little"),
  algo = "sha256", serialize = FALSE)

# Frozen 3,001-point seed-1 dyad record retained from the 0.0.21 shootout.
dyad_model <- model_dyad()
dyad_sim <- simulate_model(
  dyad_model, seed = 1, T = 3, dt = 0.001, burn_in = 0)
dyad_obs <- as_obs(dyad_sim)
dyad_init <- list(mean = 2, cov = matrix(0.1, 1L, 1L))
dyad_hash <- hash_double(dyad_obs$x[, 1L])
expected_dyad_hash <-
  "64d44c33341708e49c6946b779757e02d2c37033e567ee40a1bcd0784c5cbed1"
stopifnot(length(dyad_obs$t) == 3001L,
          identical(dyad_hash, expected_dyad_hash))

# Fixed affine k=2,l=2 input: authenticated batch realization plus matrix kernels.
matrix_model <- cgns_from_affine(
  f_full = function(t, x, y) c(
    -0.25 * x[1L] + 0.08 + (0.65 + 0.04 * x[1L]) * y[1L] -
      0.12 * y[2L],
    -0.18 * x[2L] - 0.03 + 0.08 * y[1L] +
      (0.45 + 0.03 * x[2L]) * y[2L]),
  g_full = function(t, x, y) c(
    0.10 * sin(t) - 0.06 * x[1L] - 0.70 * y[1L] - 0.04 * y[2L],
    0.12 * cos(t) + 0.04 * x[2L] + 0.09 * y[1L] - 0.55 * y[2L]),
  Sx = function(t, x) matrix(c(0.65, 0.05, 0, 0.55), 2L, 2L),
  Sy_hidden = function(t, x) matrix(c(0.72, 0.03, 0, 0.64), 2L, 2L),
  Sy_shared = function(t, x) matrix(c(0.05, 0.01, 0, 0.04), 2L, 2L),
  k = 2L, l = 2L, name = "benchmark-affine-matrix")
matrix_t <- seq(0, 2, by = 0.005)
matrix_x <- cbind(0.25 + 0.08 * sin(4 * matrix_t),
                  -0.15 + 0.06 * cos(3 * matrix_t))
matrix_obs <- observed_trajectory(matrix_t, matrix_x)
matrix_init <- list(
  mean = c(0.15, -0.08),
  cov = matrix(c(0.42, 0.025, 0.025, 0.31), 2L, 2L))

# Fixed generic closure input: masked-innovation conditioning.
conditioning_model <- cgns_model(
  Lx = function(t, x)
    matrix(c(0.7 + 0.1 * x[1L], -0.25 + 0.05 * x[2L]), 2L, 1L),
  fx = function(t, x)
    c(-0.3 * x[1L] + 0.1 * sin(t), -0.2 * x[2L] - 0.05),
  Ly = function(t, x) matrix(-0.8 + 0.02 * x[1L], 1L, 1L),
  fy = function(t, x) 0.15 * cos(t) - 0.1 * x[1L],
  Sx1 = function(t, x) diag(c(0.7, 0.6), 2L),
  Sx2 = function(t, x) matrix(0, 2L, 1L),
  Sy1 = function(t, x) matrix(c(0.04, 0), 1L, 2L),
  Sy2 = function(t, x) matrix(0.8, 1L, 1L),
  k = 2L, l = 1L, name = "benchmark-conditioned-generic")
conditioning_t <- seq(0, 2, by = 0.005)
conditioning_x <- cbind(
  target = 0.3 + 0.1 * sin(3 * conditioning_t),
  nontarget = -0.2 + 0.08 * cos(2 * conditioning_t))
conditioning_obs <- observed_trajectory(conditioning_t, conditioning_x)
conditioning_nt <- nontarget(2L, "inflate")
conditioning_init <- list(mean = 0.2, cov = matrix(0.4, 1L, 1L))

# Bounded range scenario: 201 points and a maximum of 25 retained lags.
range_idx <- seq_len(201L)
range_obs <- observed_trajectory(
  dyad_obs$t[range_idx], dyad_obs$x[range_idx, , drop = FALSE])
range_init <- dyad_init

input_manifest <- data.frame(
  scenario = c("dyad_3001", "affine_matrix", "conditioned_generic",
               "bounded_range"),
  points = c(length(dyad_obs$t), length(matrix_t), length(conditioning_t),
             length(range_idx)),
  k = c(1L, 2L, 2L, 1L), l = c(1L, 2L, 1L, 1L),
  dt = c(dyad_obs$dt, matrix_obs$dt, conditioning_obs$dt, range_obs$dt),
  observation_sha256 = c(
    dyad_hash, hash_double(matrix_x), hash_double(conditioning_x),
    hash_double(range_obs$x)),
  construction = c(
    "model_dyad; simulate seed=1,T=3,dt=0.001,burn_in=0",
    "fixed analytic observations; cgns_from_affine",
    "fixed analytic observations; generic closures; inflate channel 2",
    "first 201 points of dyad_3001; max_lag=25 for lag table"),
  stringsAsFactors = FALSE)
utils::write.csv(input_manifest, file.path(output_dir, "inputs.csv"),
                 row.names = FALSE)

# Trusted warm inputs are prepared outside timed expressions.
dyad_bundle <- compile_run(dyad_model, dyad_obs)
dyad_generic_bundle <- compile_generic(dyad_model, dyad_obs)
dyad_filter <- filter_compiled(dyad_bundle, dyad_init, validate = FALSE)
dyad_smoother <- smoother_compiled(dyad_bundle, dyad_filter, validate = FALSE)
dyad_private <- complete_compiled(dyad_bundle, dyad_init)
dyad_public <- aci_public(dyad_model, dyad_obs, init = dyad_init)

matrix_bundle <- compile_run(matrix_model, matrix_obs)
matrix_filter <- filter_compiled(matrix_bundle, matrix_init, validate = FALSE)
matrix_smoother <- smoother_compiled(
  matrix_bundle, matrix_filter, validate = FALSE)
matrix_private <- complete_compiled(matrix_bundle, matrix_init)
matrix_public <- aci_public(matrix_model, matrix_obs, init = matrix_init)

conditioning_bundle <- compile_run(
  conditioning_model, conditioning_obs, conditioning_nt)
conditioning_private <- complete_compiled(
  conditioning_bundle, conditioning_init)
conditioning_public <- aci_public(
  conditioning_model, conditioning_obs, init = conditioning_init,
  nontarget = conditioning_nt)

range_bundle <- compile_run(dyad_model, range_obs)
range_filter <- filter_compiled(range_bundle, range_init, validate = FALSE)
range_result <- complete_compiled(range_bundle, range_init)
range_lag_private <- lag_compiled(
  range_bundle, mode = "forward", max_lag = 25L,
  filter = range_filter, validate = FALSE)
range_lag_public <- lag_public(
  dyad_model, range_obs, mode = "forward", max_lag = 25L,
  filter = range_filter, init = range_init)
range_cir_private <- cir_compiled(
  range_bundle, filter = range_filter, method = "exact", min_M = 0)
range_cir_public <- forward_cir(range_result, method = "exact", min_M = 0)

max_error <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || any(is.na(a) != is.na(b))) return(Inf)
  keep <- !is.na(a)
  if (any(is.infinite(a[keep]) != is.infinite(b[keep])) ||
      any(a[keep][is.infinite(a[keep])] != b[keep][is.infinite(b[keep])]))
    return(Inf)
  d <- abs(a[keep & is.finite(a)] - b[keep & is.finite(b)])
  if (!length(d)) 0 else max(d)
}
bundle_field_error <- function(a, b) max(vapply(
  names(a$coefficients),
  function(field) max_error(a$coefficients[[field]], b$coefficients[[field]]),
  numeric(1L)))
result_parity <- function(prefix, a, b, tol = 1e-12) data.frame(
  quantity = paste0(prefix, c("_filter_mean", "_filter_cov",
                              "_filter_loglik", "_smoother_mean",
                              "_smoother_cov", "_metric")),
  max_abs_error = c(
    max_error(a$paths$filter$mean, b$paths$filter$mean),
    max_error(a$paths$filter$cov, b$paths$filter$cov),
    max_error(a$paths$filter$meta$loglik, b$paths$filter$meta$loglik),
    max_error(a$paths$smoother$mean, b$paths$smoother$mean),
    max_error(a$paths$smoother$cov, b$paths$smoother$cov),
    max_error(a$aci, b$aci)),
  tolerance = tol, stringsAsFactors = FALSE)
parity <- rbind(
  data.frame(
    quantity = c("dyad_production_vs_generic_coefficients",
                 "range_lag_diagonal", "range_lag_tail_bound",
                 "range_cir_tau", "range_cir_strength"),
    max_abs_error = c(
      bundle_field_error(dyad_bundle, dyad_generic_bundle),
      max_error(range_lag_private$diag, range_lag_public$diag),
      max_error(range_lag_private$tailbnd, range_lag_public$tailbnd),
      max_error(range_cir_private$tau, range_cir_public$tau),
      max_error(range_cir_private$M, range_cir_public$M)),
    tolerance = c(0, rep(1e-12, 4L)), stringsAsFactors = FALSE),
  result_parity("dyad", dyad_private, dyad_public),
  result_parity("matrix", matrix_private, matrix_public, 1e-11),
  result_parity("conditioning", conditioning_private, conditioning_public,
                1e-11))
parity$pass <- parity$max_abs_error <= parity$tolerance
utils::write.csv(parity, file.path(output_dir, "parity.csv"), row.names = FALSE)
if (!all(parity$pass)) stop("Production/public parity failed; timings rejected.")
stopifnot(identical(dyad_bundle$realization, "dyad_directed"),
          identical(matrix_bundle$realization, "affine_batch"),
          identical(conditioning_bundle$realization,
                    "generic_closure_one_pass"))

expressions <- list()
stage_rows <- list()
reps <- function(full, fast = max(3L, ceiling(full / 10)))
  if (quick) as.integer(fast) else as.integer(full)
add_stage <- function(name, scenario, layer, fun, repetitions, contract,
                      compilation = FALSE, likelihood = FALSE,
                      validation = FALSE, construction = FALSE,
                      retained_lag_rows = FALSE) {
  expressions[[name]] <<- fun
  stage_rows[[name]] <<- data.frame(
    stage = name, scenario = scenario, layer = layer,
    repetitions = repetitions, contract = contract,
    includes_compilation = compilation,
    includes_predictive_likelihood = likelihood,
    includes_validation = validation,
    includes_result_construction = construction,
    retains_lag_rows = retained_lag_rows,
    stringsAsFactors = FALSE)
}

add_stage("dyad_compile_production", "dyad_3001", "compilation",
          function() compile_run(dyad_model, dyad_obs), reps(51L),
          "authenticated production front door", TRUE, FALSE, TRUE, TRUE)
add_stage("dyad_compile_generic_fallback", "dyad_3001", "compilation",
          function() compile_generic(dyad_model, dyad_obs), reps(11L),
          "one-pass generic closure fallback", TRUE, FALSE, TRUE, TRUE)
add_stage("dyad_filter_warm", "dyad_3001", "warm_kernel",
          function() filter_compiled(dyad_bundle, dyad_init, validate = FALSE),
          reps(51L), "precompiled filter; likelihood and path construction",
          FALSE, TRUE, FALSE, TRUE)
add_stage("dyad_smoother_warm", "dyad_3001", "warm_kernel",
          function() smoother_compiled(
            dyad_bundle, dyad_filter, validate = FALSE), reps(51L),
          "precompiled smoother with trusted filter", FALSE, FALSE, FALSE, TRUE)
add_stage("dyad_metric_warm", "dyad_3001", "warm_kernel",
          function() metric_compiled(
            dyad_bundle, dyad_smoother, dyad_filter, validate = FALSE), reps(51L),
          "precompiled decomposed KL path", FALSE, FALSE, FALSE, TRUE)
add_stage("dyad_complete_bundle", "dyad_3001", "bundle_complete",
          function() complete_compiled(dyad_bundle, dyad_init), reps(31L),
          "complete CGNS execution over existing bundle", FALSE, TRUE, TRUE,
          TRUE)
add_stage("dyad_public_filter", "dyad_3001", "public_workflow",
          function() da_filter(dyad_model, dyad_obs, init = dyad_init), reps(31L),
          "public filter; compiles and validates", TRUE, TRUE, TRUE, TRUE)
add_stage("dyad_public_smoother_supplied", "dyad_3001", "public_workflow",
          function() da_smooth(dyad_model, dyad_obs, filter = dyad_filter),
          reps(31L), "public smoother; compiles and validates supplied filter",
          TRUE, FALSE, TRUE, TRUE)
add_stage("dyad_public_aci", "dyad_3001", "public_workflow",
          function() aci_public(dyad_model, dyad_obs, init = dyad_init), reps(31L),
          "complete public ACI contract", TRUE, TRUE, TRUE, TRUE)

add_stage("matrix_compile_production", "affine_matrix", "compilation",
          function() compile_run(matrix_model, matrix_obs), reps(21L),
          "authenticated affine batch realiser", TRUE, FALSE, TRUE, TRUE)
add_stage("matrix_filter_warm", "affine_matrix", "warm_kernel",
          function() filter_compiled(
            matrix_bundle, matrix_init, validate = FALSE), reps(21L),
          "precompiled k=2,l=2 filter", FALSE, TRUE, FALSE, TRUE)
add_stage("matrix_smoother_warm", "affine_matrix", "warm_kernel",
          function() smoother_compiled(
            matrix_bundle, matrix_filter, validate = FALSE), reps(21L),
          "precompiled k=2,l=2 smoother", FALSE, FALSE, FALSE, TRUE)
add_stage("matrix_metric_warm", "affine_matrix", "warm_kernel",
          function() metric_compiled(
            matrix_bundle, matrix_smoother, matrix_filter, validate = FALSE),
          reps(21L), "precompiled k=2,l=2 decomposed KL path", FALSE, FALSE,
          FALSE, TRUE)
add_stage("matrix_complete_bundle", "affine_matrix", "bundle_complete",
          function() complete_compiled(matrix_bundle, matrix_init), reps(11L),
          "complete matrix execution over existing bundle", FALSE, TRUE, TRUE,
          TRUE)
add_stage("matrix_public_aci", "affine_matrix", "public_workflow",
          function() aci_public(matrix_model, matrix_obs, init = matrix_init),
          reps(11L), "complete public affine-matrix ACI", TRUE, TRUE, TRUE, TRUE)

add_stage("conditioning_compile_production", "conditioned_generic",
          "compilation", function() compile_run(
            conditioning_model, conditioning_obs, conditioning_nt), reps(11L),
          "generic compile plus inflate conditioning", TRUE, FALSE, TRUE, TRUE)
add_stage("conditioning_complete_bundle", "conditioned_generic",
          "bundle_complete", function() complete_compiled(
            conditioning_bundle, conditioning_init), reps(21L),
          "conditioned execution over existing bundle", FALSE, TRUE, TRUE, TRUE)
add_stage("conditioning_public_aci", "conditioned_generic", "public_workflow",
          function() aci_public(
            conditioning_model, conditioning_obs, init = conditioning_init,
            nontarget = conditioning_nt), reps(11L),
          "complete public conditioned ACI", TRUE, TRUE, TRUE, TRUE)

add_stage("range_lag_bundle", "bounded_range", "warm_lag_cir",
          function() lag_compiled(
            range_bundle, mode = "forward", max_lag = 25L,
            filter = range_filter, validate = FALSE), reps(11L),
          "precompiled forward lag table; max_lag=25", FALSE, FALSE, TRUE,
          TRUE, TRUE)
add_stage("range_public_lag", "bounded_range", "public_workflow",
          function() lag_public(
            dyad_model, range_obs, mode = "forward", max_lag = 25L,
            filter = range_filter, init = range_init), reps(11L),
          "public bounded forward lag table", TRUE, FALSE, TRUE, TRUE, TRUE)
add_stage("range_cir_streamed_bundle", "bounded_range", "warm_lag_cir",
          function() cir_compiled(
            range_bundle, filter = range_filter, method = "exact", min_M = 0),
          reps(11L), "precompiled streaming CIR; no retained triangle", FALSE,
          FALSE, TRUE, TRUE, FALSE)
add_stage("range_public_cir", "bounded_range", "public_workflow",
          function() forward_cir(
            range_result, method = "exact", min_M = 0), reps(11L),
          "public CIR from ACI result without retained table", TRUE, FALSE,
          TRUE, TRUE, FALSE)

stages <- do.call(rbind, stage_rows)
row.names(stages) <- NULL
utils::write.csv(stages, file.path(output_dir, "stages.csv"), row.names = FALSE)

# Child mode profiles one stage after identical setup and four JIT warmups.
if (!is.null(allocation_stage)) {
  if (!allocation_stage %in% names(expressions))
    stop("Unknown allocation stage: ", allocation_stage)
  for (i in seq_len(4L)) invisible(expressions[[allocation_stage]]())
  gc()
  log_path <- file.path(output_dir, paste0("alloc-", allocation_stage, ".log"))
  Rprofmem(log_path, threshold = 0)
  invisible(expressions[[allocation_stage]]())
  Rprofmem(NULL)
  quit(save = "no", status = 0L)
}

timing_rows <- vector("list", length(expressions))
names(timing_rows) <- names(expressions)
for (stage in names(expressions)) {
  fun <- expressions[[stage]]
  invisible(fun())
  gc()
  measured <- microbenchmark::microbenchmark(
    fun(), times = stages$repetitions[match(stage, stages$stage)],
    unit = "ns", control = list(warmup = 2L))
  timing_rows[[stage]] <- data.frame(
    stage = stage, iteration = seq_along(measured$time),
    nanoseconds = as.numeric(measured$time),
    seconds = as.numeric(measured$time) / 1e9,
    stringsAsFactors = FALSE)
}
timings_raw <- do.call(rbind, timing_rows)
row.names(timings_raw) <- NULL
utils::write.csv(timings_raw, file.path(output_dir, "timings-raw.csv"),
                 row.names = FALSE)
summary_rows <- lapply(split(timings_raw$seconds, timings_raw$stage), function(z)
  data.frame(n = length(z), min = min(z), q1 = unname(quantile(z, 0.25)),
             median = median(z), q3 = unname(quantile(z, 0.75)),
             max = max(z), iqr = IQR(z)))
timing_summary <- do.call(rbind, summary_rows)
timing_summary$stage <- row.names(timing_summary)
row.names(timing_summary) <- NULL
timing_summary <- merge(
  stages[, c("stage", "scenario", "layer")], timing_summary, by = "stage",
  sort = FALSE)
utils::write.csv(timing_summary,
                 file.path(output_dir, "timings-summary.csv"), row.names = FALSE)

allocation_summary <- data.frame()
if (!no_alloc) {
  child_script <- file.path(tempdir(), "aci-production-benchmark-child.R")
  if (!file.copy(script_file, child_script, overwrite = TRUE))
    stop("Could not prepare the allocation child script.")
  common <- c(
    "--vanilla", child_script,
    shQuote(paste0("--package=", package_dir)),
    shQuote(paste0("--output=", output_dir)), "--no-profile")
  if (quick) common <- c(common, "--quick")
  for (stage in names(expressions)) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c(common, shQuote(paste0("--allocation-stage=", stage))),
                      stdout = FALSE, stderr = FALSE)
    if (!identical(status, 0L))
      stop("Allocation child failed for stage: ", stage)
  }
  allocation_summary <- do.call(rbind, lapply(names(expressions), function(stage) {
    lines <- readLines(
      file.path(output_dir, paste0("alloc-", stage, ".log")), warn = FALSE)
    sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    data.frame(
      stage = stage, total_events = length(lines),
      sized_events = sum(is.finite(sizes)),
      sized_bytes = sum(sizes[is.finite(sizes)]),
      new_pages = sum(startsWith(lines, "new page:")),
      stringsAsFactors = FALSE)
  }))
  allocation_summary <- merge(
    stages[, c("stage", "scenario", "layer")], allocation_summary,
    by = "stage", sort = FALSE)
}
utils::write.csv(allocation_summary,
                 file.path(output_dir, "allocations-summary.csv"),
                 row.names = FALSE)

profile_repetitions <- c(
  dyad_compile_production = 30L,
  dyad_complete_bundle = 10L,
  dyad_public_aci = 10L,
  matrix_public_aci = 2L,
  conditioning_public_aci = 2L,
  range_public_lag = 2L,
  range_public_cir = 2L)
if (!no_profile) {
  for (stage in names(profile_repetitions)) {
    path <- file.path(output_dir, paste0("profile-", stage, ".Rprof"))
    gc()
    Rprof(path, interval = 0.001, memory.profiling = TRUE)
    for (i in seq_len(if (quick) 1L else profile_repetitions[[stage]]))
      invisible(expressions[[stage]]())
    Rprof(NULL)
    profile <- tryCatch(
      summaryRprof(path, memory = "both"),
      error = function(e) summaryRprof(path, memory = "none"))
    utils::write.csv(
      profile$by.total,
      file.path(output_dir, paste0("profile-", stage, "-by-total.csv")))
    utils::write.csv(
      profile$by.self,
      file.path(output_dir, paste0("profile-", stage, "-by-self.csv")))
  }
}

git_identity <- function(path) {
  root <- tryCatch(suppressWarnings(system2(
    "git", c("-C", shQuote(path), "rev-parse", "--show-toplevel"),
    stdout = TRUE, stderr = FALSE)), error = function(e) character())
  if (!length(root) || !is.null(attr(root, "status")))
    return(list(root = "absent", head = "absent", dirty = "not_applicable"))
  head <- tryCatch(suppressWarnings(system2(
    "git", c("-C", shQuote(path), "rev-parse", "HEAD"), stdout = TRUE,
    stderr = FALSE)), error = function(e) character())
  status <- tryCatch(suppressWarnings(system2(
    "git", c("-C", shQuote(path), "status", "--porcelain=v1",
             "--untracked-files=all"),
    stdout = TRUE, stderr = FALSE)),
    error = function(e) structure(character(), status = 1L))
  list(root = root[[1L]],
       head = if (length(head)) head[[1L]] else "unavailable",
       dirty = if (!is.null(attr(status, "status"))) "unknown" else
         as.character(length(status) > 0L))
}
source_manifest <- function(path) {
  files <- c(file.path(path, c("DESCRIPTION", "NAMESPACE")),
             sort(list.files(file.path(path, "R"), pattern = "[.]R$",
                             full.names = TRUE)),
             file.path(path, "benchmarks", c("scalar-dyad.R", "README.md",
                                             "DECISION.md")))
  files <- files[file.exists(files)]
  data.frame(
    source = substring(files, nchar(path) + 2L),
    sha256 = vapply(files, digest::digest, character(1L),
                    algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE)
}
source_hashes <- source_manifest(package_dir)
source_hash_path <- file.path(output_dir, "source-hashes.csv")
utils::write.csv(source_hashes, source_hash_path, row.names = FALSE)
description <- read.dcf(file.path(package_dir, "DESCRIPTION"))
git <- git_identity(package_dir)
environment <- c(
  paste0("timestamp=", format(Sys.time(), usetz = TRUE)),
  paste0("package_dir=", package_dir),
  paste0("package=", description[1L, "Package"]),
  paste0("package_version=", description[1L, "Version"]),
  paste0("git_root=", git$root), paste0("git_head=", git$head),
  paste0("git_dirty=", git$dirty),
  paste0("source_hash_manifest_sha256=", digest::digest(
    source_hash_path, algo = "sha256", file = TRUE)),
  paste0("R.version=", R.version.string),
  paste0("platform=", R.version$platform),
  paste0("OS=", Sys.info()[["sysname"]], " ", Sys.info()[["release"]]),
  paste0("machine=", Sys.info()[["machine"]]),
  paste0("locale=", Sys.getlocale()),
  paste0("timezone=", Sys.timezone()),
  paste0("BLAS=", unname(extSoftVersion()["BLAS"])),
  paste0("LAPACK=", unname(extSoftVersion()["LAPACK"])),
  paste0("La_library=", paste(La_library(), collapse = " ")),
  paste0("La_version=", La_version()),
  paste0("microbenchmark=", utils::packageVersion("microbenchmark")),
  paste0("pkgload=", utils::packageVersion("pkgload")),
  paste0("digest=", utils::packageVersion("digest")),
  paste0("RNGkind=", paste(RNGkind(), collapse = "/")),
  paste0("quick=", quick), paste0("allocations_enabled=", !no_alloc),
  paste0("profiles_enabled=", !no_profile),
  paste0("dyad_realization=", dyad_bundle$realization),
  paste0("matrix_realization=", matrix_bundle$realization),
  paste0("conditioning_realization=", conditioning_bundle$realization),
  "timing_policy=one untimed expression warmup; gc per stage; microbenchmark timer warmup 2; raw observations retained",
  "allocation_policy=fresh R process per stage; four warmups; gc; one Rprofmem expression; sized bytes and unsized new-page events retained",
  "comparison_policy=no headline ratio; compilation, warm kernels, bundle completion, and public workflows are separate contracts",
  "aci_0.0.21_policy=not loaded into this process; use archived results or a separate R process")
writeLines(environment, file.path(output_dir, "environment.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))

cat("Production compiled-CGNS benchmark complete\n")
cat("Output:", output_dir, "\n")
print(timing_summary, row.names = FALSE)
