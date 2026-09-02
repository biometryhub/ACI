## =========================================================================
## make-partition-fixtures.R
##
## Regenerates the scalar-ENSO-partition oracle fixtures graded by
## tests/testthat/test-28-partition-oracles.R and pinned by
## tests/testthat/fixtures/oracles/oracle-manifest-partitions.yml.
##
## Nothing here is hand-copied: the driving signal is re-simulated from its
## seed, every coefficient array and every filter/smoother/ACI series is
## recomputed, and the pinned SHA-256 of the signal is checked before anything
## else runs.
##
## The producer of the pinned reference values is `aci` 0.0.30, deliberately
## NOT the package under test: the fixtures grade acir against a second
## implementation's numbers rather than against its own.  Run with
## PRODUCER = "acir" to see whether this package reproduces the pinned bytes,
## which is the same question test-28 asks through the public API.
##
##   Rscript tools/fixtures/make-partition-fixtures.R <outdir> \
##           [producer] [producer_lib]
##
## <outdir> is required and must be OUTSIDE the package.  Fixtures are never
## regenerated in place: write to a scratch directory, diff the SHA-256
## against oracle-manifest-partitions.yml, and replace the shipped files only
## as a deliberate maintainer step with new hashes and a NEWS entry.
##
## Optional regeneration proof.  The original run also re-derived the C3
## comparison outputs in their original layout and compared their SHA-256
## against the C3 originals.  Those originals live in the comparison workspace,
## not in this repository, so the proof step is opt-in:
##
##   ACIR_FIXTURE_C3_DIR=/path/to/comparison/c3 Rscript ... <outdir>
##
## Without it the step is skipped and reported as skipped.  The proof is
## recorded under `refresh: regeneration_proof:` in the manifest and does not
## need to re-run to reproduce the shipped bytes.
## =========================================================================

args <- commandArgs(trailingOnly = TRUE)
if( !length(args) )
  stop("usage: Rscript make-partition-fixtures.R <outdir> [producer] [producer_lib]")
OUTDIR   <- normalizePath(args[[1L]], mustWork = FALSE)
PRODUCER <- if( length(args) >= 2L ) args[[2L]] else "aci"
if( length(args) >= 3L ) .libPaths(c(args[[3L]], .libPaths()))
C3       <- Sys.getenv("ACIR_FIXTURE_C3_DIR", "")

FIX <- file.path(OUTDIR, "fixtures")
VER <- file.path(OUTDIR, "verify")
dir.create(FIX, showWarnings = FALSE, recursive = TRUE)
dir.create(VER, showWarnings = FALSE, recursive = TRUE)

suppressMessages(requireNamespace(PRODUCER, quietly = TRUE))
P  <- asNamespace(PRODUCER)
gg <- function(nm) get(nm, envir = P)

sha256 <- function(f) {
  z <- system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE)
  sub(" .*$", "", z)
}
md5 <- function(f) unname(tools::md5sum(f))

cat(sprintf("producer: %s %s   out: %s\n", PRODUCER,
            as.character(utils::packageVersion(PRODUCER)), OUTDIR))
cat(sprintf("regeneration proof against C3: %s\n\n",
            if( nzchar(C3) ) C3 else "SKIPPED (ACIR_FIXTURE_C3_DIR unset)"))

## -------------------------------------------------------------------------
## 1. The driving signal, re-simulated from the pinned generator settings.
##    aci_enso_model(variant = "aci_code", hidden = c("u","hW","tau")) then
##    simulate(seed = 42, T = 20, dt = 0.005).
## -------------------------------------------------------------------------
SIGNAL_SHA <- "77412008549a2980af8bc843c3548a673a4c2df952a40aea12ef2d30d8dec250"

sig_file <- file.path(FIX, "enso6_partition_signal.csv")
em <- gg("aci_enso_model")(variant = "aci_code", hidden = c("u", "hW", "tau"))
es <- simulate(em, seed = 42, T = 20, dt = 0.005)
sx <- as.matrix(es$obs$x); sy <- as.matrix(es$hidden)
colnames(sx) <- c("TC", "TE", "I"); colnames(sy) <- c("u", "hW", "tau")
utils::write.csv(cbind(t = es$obs$t, sx, sy), sig_file, row.names = FALSE)
cat(sprintf("signal regenerated: sha256 %s  match_pinned=%s\n",
            substr(sha256(sig_file), 1L, 16L),
            identical(sha256(sig_file), SIGNAL_SHA)))
stopifnot(identical(sha256(sig_file), SIGNAL_SHA))

path <- utils::read.csv(sig_file)
N1   <- nrow(path)
IDX  <- seq(1L, N1, by = 20L)          ## 201 sampled indices, house idiom
stopifnot(N1 == 4001L, length(IDX) == 201L, IDX[1L] == 1L, IDX[201L] == N1)

## -------------------------------------------------------------------------
## 2. Helpers
## -------------------------------------------------------------------------
## Column block for a realised coefficient bundle, in the C3 layout. gxx is
## stored by its diagonal only; the off-diagonal is verified exactly zero on
## every run and the check is recorded rather than the zeros.
coef_block <- function(b, ov) {
  co <- b$coefficients; k <- b$k
  d  <- data.frame(t = b$t)
  for( i in seq_len(k) ) d[[paste0("Lx_", ov[i])]]       <- co$Lx[i, 1L, ]
  for( i in seq_len(k) ) d[[paste0("fx_", ov[i])]]       <- co$fx[, i]
  d$Ly <- co$Ly[1L, 1L, ]
  d$fy <- co$fy[, 1L]
  for( i in seq_len(k) ) d[[paste0("gxx_diag_", ov[i])]] <- co$gxx[i, i, ]
  d$gyy <- co$gyy[1L, 1L, ]
  for( i in seq_len(k) ) d[[paste0("gyx_", ov[i])]]      <- co$gyx[1L, i, ]
  return(d)
}

## Structural invariants pinned alongside the numbers.
invariants <- function(b) {
  co  <- b$coefficients
  off <- max(abs(co$gxx - array(apply(co$gxx, 3L, function(z) diag(diag(z))),
                                dim(co$gxx))))
  gi  <- max(abs(vapply(seq_len(b$N),
                        function(j) max(abs(co$gxx_weight[, , j] %*% co$gxx[, , j] -
                                            diag(b$k))),
                        numeric(1L))))
  return(list(gxx_offdiag_max = off,
              gyx_abs_max     = max(abs(co$gyx)),
              gram_inverse_max_residual = gi,
              gxx_weight_slices = dim(co$gxx_weight)[3L],
              coefficient_slices = b$N1))
}

## Per-series full-record summary. Closes the gap left by sampling: a defect
## at any of the 3800 unsampled steps moves at least one of these.
summarise <- function(part, d) {
  nm <- setdiff(names(d), "t")
  do.call(rbind, lapply(nm, function(k) data.frame(
    partition = part, series = k, n = nrow(d),
    min = min(d[[k]]), max = max(d[[k]]),
    mean = mean(d[[k]]), sum_abs = sum(abs(d[[k]])))))
}

## Bit-for-bit regeneration proof against a C3 original. Skipped, and reported
## as skipped, when the comparison workspace is not on hand.
verify <- function(df, name) {
  if( !nzchar(C3) ) { cat(sprintf("  regen %-34s SKIPPED\n", name)); return(NA) }
  f <- file.path(VER, name)
  utils::write.csv(df, f, row.names = FALSE)
  o <- file.path(C3, name)
  ok <- file.exists(o) && identical(sha256(f), sha256(o))
  cat(sprintf("  regen %-34s sha256 %s  c3_identical=%s\n",
              name, substr(sha256(f), 1L, 16L), ok))
  return(ok)
}

## Round-trip proof for the output series. The C3 outputs files interleave the
## aciR arm, so their bytes cannot be rebuilt here. Instead push the freshly
## computed doubles through the SAME write.csv formatting C3 used, read them
## back, and require bit-for-bit equality with the C3 columns. Equality proves
## the regenerated doubles agree with C3's to the 15 significant digits the
## file preserves; any residual in-memory difference is the file's
## quantisation, not a computation difference.
roundtrip <- function(regen, c3file, cmp, tag) {
  if( !nzchar(C3) ) { cat(sprintf("  roundtrip %-14s SKIPPED\n", tag)); return(NA) }
  o3 <- utils::read.csv(file.path(C3, c3file))
  c3cols <- stats::setNames(o3[cmp], names(cmp))
  f <- file.path(VER, sprintf("roundtrip_%s.csv", tag))
  utils::write.csv(regen, f, row.names = FALSE)
  rb <- utils::read.csv(f)
  bit <- vapply(names(regen), function(k) identical(rb[[k]], c3cols[[k]]),
                logical(1L))
  ulp <- vapply(names(regen), function(k) {
    v <- abs(c3cols[[k]]); v <- max(v[is.finite(v)])
    if( v == 0 ) return(Inf)
    max(abs(regen[[k]] - c3cols[[k]])) / (v * 1e-15)
  }, numeric(1L))
  cat(sprintf("  roundtrip %-14s bit_identical=%s  worst in-memory gap = %.2f x (max|v| * 1e-15)\n",
              tag, all(bit), max(ulp[is.finite(ulp)])))
  if( !all(bit) ) cat("    columns not bit-identical:",
                      paste(names(bit)[!bit], collapse = ", "), "\n")
  return(all(bit))
}

## The reference-series column names the C3 outputs files use.
CMP_HIDDEN <- c(filter_mean = "aci_filter_mean", filter_cov = "aci_filter_cov",
                smooth_mean = "aci_smooth_mean", smooth_cov = "aci_smooth_cov",
                aci = "aci_aci", truth = "truth")
CMP_RED3   <- c(filter_mean = "red3_filter_mean", filter_cov = "red3_filter_cov",
                smooth_mean = "red3_smooth_mean", smooth_cov = "red3_smooth_cov",
                aci = "red3_aci", truth = "truth")

## The shipped fixture layout: sampled indices, coefficients and reference
## series in one file.
write_fixture <- function(cb, ref, file) {
  fx <- cbind(data.frame(index = IDX), cb[IDX, , drop = FALSE],
              ref[IDX, setdiff(names(ref), "t"), drop = FALSE])
  rownames(fx) <- NULL
  names(fx)[names(fx) %in% names(cb)[-1L]] <- paste0("coef_", names(cb)[-1L])
  names(fx)[names(fx) %in% c("filter_mean", "filter_cov", "smooth_mean",
                             "smooth_cov", "aci")] <-
    paste0("ref_", c("filter_mean", "filter_cov", "smooth_mean",
                     "smooth_cov", "aci"))
  utils::write.csv(fx, file.path(FIX, file), row.names = FALSE)
  return(invisible(NULL))
}

## The observation set is named on every construction.  `acir` defaults
## hidden = "tau" to the reduced three-channel estimand (C2c) while `aci`
## 0.0.30 has no such argument, so the argument is passed only when the
## producer accepts it.  Leaving it to the default would make the two
## producers realise DIFFERENT estimands from the same call.
build_model <- function(hidden, observations) {
  f <- gg("aci_enso_model")
  if( "observations" %in% names(formals(f)) )
    return(f(variant = "aci_code", hidden = hidden,
             observations = observations))
  if( !identical(observations, "full") )
    stop("this producer's aci_enso_model() has no `observations` argument, so ",
         "only the full observation set can be built directly")
  return(f(variant = "aci_code", hidden = hidden))
}

## -------------------------------------------------------------------------
## 3. The three scalar-hidden partitions (C3 step 1-2 recomputed)
## -------------------------------------------------------------------------
inv <- list(); summ <- list(); regen_ok <- c()

for( h in c("u", "hW", "tau") ) {
  m    <- build_model(h, "full")
  ov   <- m$meta$vars$observed
  ob   <- gg("observed_trajectory")(path$t, as.matrix(path[, ov, drop = FALSE]))
  init <- list(mean = as.numeric(path[[h]][1L]), cov = matrix(0.1, 1L, 1L))

  b  <- get(".compile_cgns_run", envir = P)(m, ob, NULL)
  fA <- gg("aci_filter")(m, ob, init = init)
  sA <- gg("aci_smoother")(m, ob, filter = fA, init = init)
  rA <- gg("aci")(m, ob, init = init, decompose = TRUE)

  cb <- coef_block(b, ov)
  inv[[h]] <- invariants(b)

  ## --- regeneration proof, C3 layout, full length -------------------------
  regen_ok[sprintf("components_%s_hidden.csv", h)] <-
    verify(cb, sprintf("components_%s_hidden.csv", h))

  ref <- data.frame(
    t             = b$t,
    filter_mean   = as.numeric(fA$mean),
    filter_cov    = as.numeric(fA$cov),
    smooth_mean   = as.numeric(sA$mean),
    smooth_cov    = as.numeric(sA$cov),
    aci           = as.numeric(rA$aci),
    truth         = path[[h]])

  regen_ok[sprintf("outputs_%s_hidden.csv", h)] <-
    roundtrip(ref[names(CMP_HIDDEN)], sprintf("outputs_%s_hidden.csv", h),
              CMP_HIDDEN, sprintf("%s_hidden", h))

  summ[[h]] <- rbind(summarise(h, cb), summarise(h, ref))

  write_fixture(cb, ref, sprintf("enso6_partition_%s_reference.csv", h))
}

## -------------------------------------------------------------------------
## 4. The tau reduced-3-dim partition (C3 step 3 recomputed)
## -------------------------------------------------------------------------
## Built as the pinned bytes were built: the five-channel model with the
## reduction supplied explicitly.  On acir this is bit-identical to
## aci_enso_model(hidden = "tau", observations = "reduced"), which carries the
## same specification as a declaration; test-28 asserts that identity.
m    <- build_model("tau", "full")
ov   <- m$meta$vars$observed
ob   <- gg("observed_trajectory")(path$t, as.matrix(path[, ov, drop = FALSE]))
init <- list(mean = path$tau[1L], cov = matrix(0.1, 1L, 1L))

spec  <- gg("aci_conditional")(given = c("u", "hW"), method = "reduce")
specI <- gg("aci_conditional")(given = c("u", "hW"), method = "mask")
fB <- gg("aci_filter")(m, ob, init = init, conditional = spec)
sB <- gg("aci_smoother")(m, ob, filter = fB, init = init, conditional = spec)
rB <- gg("aci")(m, ob, init = init, conditional = spec)
fC <- gg("aci_filter")(m, ob, init = init, conditional = specI)
sC <- gg("aci_smoother")(m, ob, filter = fC, init = init, conditional = specI)
rC <- gg("aci")(m, ob, init = init, conditional = specI)

bB  <- get(".compile_cgns_run", envir = P)(m, ob, spec)
ovr <- colnames(bB$x)
cbr <- coef_block(bB, ovr)
inv[["tau_reduced3"]] <- invariants(bB)
inv[["tau_reduced3"]]$mask_identical_to_reduce <-
  identical(as.numeric(fB$mean), as.numeric(fC$mean)) &&
  identical(as.numeric(fB$cov),  as.numeric(fC$cov))  &&
  identical(as.numeric(sB$mean), as.numeric(sC$mean)) &&
  identical(as.numeric(sB$cov),  as.numeric(sC$cov))  &&
  identical(as.numeric(rB$aci),  as.numeric(rC$aci))
inv[["tau_reduced3"]]$reduced_obs <- paste(ovr, collapse = ",")

regen_ok["components_tau_reduced3dim.csv"] <-
  verify(cbr, "components_tau_reduced3dim.csv")

refr <- data.frame(
  t           = bB$t,
  filter_mean = as.numeric(fB$mean),
  filter_cov  = as.numeric(fB$cov),
  smooth_mean = as.numeric(sB$mean),
  smooth_cov  = as.numeric(sB$cov),
  aci         = as.numeric(rB$aci),
  truth       = path$tau)

regen_ok["outputs_tau_reduced.csv[red3]"] <-
  roundtrip(refr[names(CMP_RED3)], "outputs_tau_reduced.csv", CMP_RED3,
            "tau_reduced3")

summ[["tau_reduced3"]] <- rbind(summarise("tau_reduced3", cbr),
                                summarise("tau_reduced3", refr))

write_fixture(cbr, refr, "enso6_partition_tau_reduced3_reference.csv")

## -------------------------------------------------------------------------
## 5. Full-record summary fixture and the invariants
## -------------------------------------------------------------------------
S <- do.call(rbind, summ); rownames(S) <- NULL
utils::write.csv(S, file.path(FIX, "enso6_partition_fullpath_summary.csv"),
                 row.names = FALSE)

cat("\ninvariants:\n")
for( k in names(inv) ) cat(sprintf("  %-13s %s\n", k,
  paste(sprintf("%s=%s", names(inv[[k]]), unlist(inv[[k]])), collapse = "  ")))

cat("\nregeneration vs C3:\n")
for( k in names(regen_ok) ) cat(sprintf("  %-40s %s\n", k, regen_ok[[k]]))

cat("\nshipped fixture hashes:\n")
ff <- sort(list.files(FIX, pattern = "\\.csv$", full.names = TRUE))
for( f in ff ) cat(sprintf("  %-46s %8.1f KB  sha256 %s  md5 %s\n",
                           basename(f), file.size(f) / 1024, sha256(f), md5(f)))

saveRDS(list(invariants = inv, regen_ok = regen_ok, idx = IDX),
        file.path(OUTDIR, "generation-meta.rds"))
cat("\nDONE\n")
