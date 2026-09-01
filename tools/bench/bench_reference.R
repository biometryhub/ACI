#!/usr/bin/env Rscript
# =============================================================================
# bench_reference.R -- the Section 8 performance stages, timed and compared
# =============================================================================
#
# Times acir on the stages the core-engine specification budgets (Section 8):
# the filter, the smoother, the online smoother at all lags, the influence
# range at all reporting anchors, on the authors' dyad record (N = 3000, the
# parity harness's dyad_reference_head), and the filter on a simulated
# climate record (six-variable ENSO, three hidden, N = 20000). One extra
# stage times an R-level loop: every other stage is reported as a ratio to
# it, so a baseline recorded on one machine compares with a run on another
# without absolute seconds meaning anything.
#
# Usage, from the repository root:
#
#   Rscript tools/bench/bench_reference.R                       # time, print
#   Rscript tools/bench/bench_reference.R --out bench.csv        # and save
#   Rscript tools/bench/bench_reference.R --baseline tools/bench/baseline.csv
#                                                               # and compare
#   Rscript tools/bench/bench_reference.R --baseline ... --gate  # fail on drift
#
# Comparison rule: a stage whose ratio exceeds the baseline's by more than
# 25 percent, or whose seconds exceed its Section 8 budget by more than 25
# percent, prints a GitHub `::warning::` line. With --gate the script exits
# 1 on any warning; without it the exit status is always 0. The budgets are
# gates from PR-10 onward, warnings before.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag) { i <- match(flag, args); if (is.na(i) || i == length(args)) NULL else args[[i + 1L]] }
gate <- "--gate" %in% args
self <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
root <- normalizePath(file.path(dirname(if (length(self)) self else "tools/bench/x"), "..", ".."))
suppressMessages(pkgload::load_all(file.path(root, "acir"), quiet = TRUE, export_all = FALSE))
source(file.path(root, "tools", "oracle", "parity", "tools", "manifest.R"))
source(file.path(root, "tools", "oracle", "parity", "tools", "dataset.R"))

budgets <- c(
  "filter, library model" = 0.02, "filter, generic model" = 0.02,
  "smoother, library model" = 0.02, "smoother, generic model" = 0.02,
  "online smoother, all lags" = 0.25, "online smoother, lag 25" = NA,
  "influence range, reporting anchors" = 0.40,
  "climate filter, 20000 steps" = 0.7, "climate smoother, 20000 steps" = NA
)

median_time <- function(f, reps) {
  secs <- numeric(reps)
  for (r in seq_len(reps)) { st <- proc.time()[["elapsed"]]; f(); secs[r] <- proc.time()[["elapsed"]] - st }
  stats::median(secs)
}
rows <- list()
record <- function(stage, seconds, reps) rows[[length(rows) + 1L]] <<- data.frame(stage = stage, seconds = seconds, reps = reps)

# -- calibration: an R-level loop, the cost every kernel here is bound by ----
calibration <- median_time(function() { s <- 0; for (i in seq_len(8e6)) s <- s + i * 1e-7; s }, 3L)

# -- the authors' dyad record ---------------------------------------------------
ds <- read_dataset(file.path(root, "tools", "oracle", "parity", "datasets", "dyad_reference_head"))
meta <- ds$meta; ar <- ds$arrays; dt <- meta$dt; N1 <- nrow(ar)
idx <- function(t) min(max(as.integer(round(t / dt)) + 1L, 1L), N1)
generic <- acir::aci_model(
  Lx = function(t, x) matrix(ar$L_x[idx(t)], 1L, 1L), fx = function(t, x) ar$f_x[idx(t)],
  Ly = function(t, x) matrix(meta$L_y, 1L, 1L), fy = function(t, x) ar$f_y[idx(t)],
  Sx1 = function(t, x) matrix(meta$Sx_1, 1L, 1L), Sx2 = function(t, x) matrix(meta$Sx_2, 1L, 1L),
  Sy1 = function(t, x) matrix(meta$Sy_1, 1L, 1L), Sy2 = function(t, x) matrix(meta$Sy_2, 1L, 1L),
  k = 1L, l = 1L, name = "dyad, generic route")
library_model <- acir::aci_dyad_model()
ob <- acir::observed_trajectory(seq(0, by = dt, length.out = N1), ar$x, names = "x")
init <- list(mean = meta$mu0, cov = matrix(meta$R0, 1L, 1L))

for (route in c("library", "generic")) {
  m <- if (route == "library") library_model else generic
  f <- NULL
  record(sprintf("filter, %s model", route), median_time(function() f <<- acir::aci_filter(m, ob, init = init), 3L), 3L)
  record(sprintf("smoother, %s model", route), median_time(function() acir::aci_smoother(m, ob, filter = f), 3L), 3L)
}
f <- acir::aci_filter(library_model, ob, init = init)
record("online smoother, all lags", median_time(function() acir::aci_online(library_model, ob, lag = Inf, filter = f), 3L), 3L)
record("online smoother, lag 25", median_time(function() acir::aci_online(library_model, ob, lag = 25L, filter = f), 3L), 3L)
res <- acir::aci(library_model, ob, init = init)
anchors <- seq.int(round(as.numeric(meta$CIRStart) / dt) + 1L, round(as.numeric(meta$CIREnd) / dt) + 1L)
record("influence range, reporting anchors",
       median_time(function() suppressWarnings(acir::aci_range(res, method = "l1_linf", quadrature = "simpson", anchors = anchors)), 3L), 3L)

# -- a climate record of the reference's size, simulated by the package --------
enso <- acir::aci_enso_model(hidden = c("u", "hW", "tau"))
sim <- stats::simulate(enso, seed = 1L, T = 100, dt = 0.005)
ob_c <- acir::as_obs(sim); init_c <- list(mean = rep(0, enso$l), cov = diag(0.1, enso$l))
fc <- NULL
record("climate filter, 20000 steps", median_time(function() fc <<- acir::aci_filter(enso, ob_c, init = init_c), 3L), 3L)
record("climate smoother, 20000 steps", median_time(function() acir::aci_smoother(enso, ob_c, filter = fc), 3L), 3L)

# -- table ----------------------------------------------------------------------
tab <- do.call(rbind, rows)
tab$ratio <- tab$seconds / calibration
tab$budget <- unname(budgets[tab$stage])
tab$calibration <- calibration
tab$machine <- paste(Sys.info()[["sysname"]], Sys.info()[["machine"]])
tab$r <- paste(R.version$major, R.version$minor, sep = ".")
tab$acir <- as.character(utils::packageVersion("acir"))
tab$commit <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
tab$date <- format(Sys.Date())
out <- opt("--out"); if (!is.null(out)) utils::write.csv(tab, out, row.names = FALSE)

warnings_out <- character(0)
warn <- function(msg) { warnings_out <<- c(warnings_out, msg); cat("::warning::", msg, "\n", sep = "") }
base <- opt("--baseline")
if (!is.null(base)) {
  b <- utils::read.csv(base, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(tab))) {
    j <- match(tab$stage[i], b$stage)
    if (is.na(j)) next
    drift <- tab$ratio[i] / b$ratio[j] - 1
    if (is.finite(drift) && drift > 0.25)
      warn(sprintf("%s: +%.0f%% against the baseline (ratio %.1f vs %.1f, %s)", tab$stage[i], 100 * drift, tab$ratio[i], b$ratio[j], b$machine[j]))
  }
}
# The Section 8 budgets are seconds on the machine the baseline was recorded
# on. With a baseline in hand they are applied as ratios to that machine's
# calibration, so a slower runner is judged on the code, not on its clock;
# without one they are applied as the seconds they were written as.
base_cal <- if (!is.null(base)) b$calibration[1] else NA_real_
for (i in seq_len(nrow(tab))) {
  bud <- tab$budget[i]
  if (!is.finite(bud)) next
  if (is.finite(base_cal)) {
    if (tab$ratio[i] > 1.25 * bud / base_cal)
      warn(sprintf("%s: ratio %.2f against a Section 8 budget of %.2f s on the baseline machine (ratio %.2f, +25%% allowance); %.3f s here",
                   tab$stage[i], tab$ratio[i], bud, bud / base_cal, tab$seconds[i]))
  } else if (tab$seconds[i] > 1.25 * bud) {
    warn(sprintf("%s: %.3f s against a Section 8 budget of %.2f s (+25%% allowance)", tab$stage[i], tab$seconds[i], bud))
  }
}
cat(sprintf("\ncalibration loop %.3f s on %s, R %s, acir %s, %s\n", calibration, tab$machine[1], tab$r[1], tab$acir[1], tab$commit[1]))
print(tab[, c("stage", "seconds", "reps", "ratio", "budget")], row.names = FALSE, right = FALSE, digits = 4)
cat(sprintf("\n%d warning(s)%s\n", length(warnings_out), if (gate && length(warnings_out)) ", gate is on: failing" else ""))
if (gate && length(warnings_out)) quit(status = 1L)
