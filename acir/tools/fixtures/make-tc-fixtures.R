## =========================================================================
## make-tc-fixtures.R
##
## Regenerates the T_C-hidden zeroth-order oracle fixtures graded by
## tests/testthat/test-29-tc-zeroth-order.R and pinned by
## tests/testthat/fixtures/oracles/oracle-manifest-partitions.yml.
##
## THE ORACLE IS THE TRANSCRIPTION BELOW, NOT A PACKAGE.  Everything under
## "independent transcription" is a statement-by-statement R rendering of the
## filter (:1198-1257), smoother (:1263-1324) and ACI (:1394-1397) sections of
## ACI_code-main/ENSO_model_cond_ACI_T_C_unobs.m, written for wave-3
## preparation task P1 without reference to acir's kernels.  MATLAB 1-based
## indexing is kept verbatim, loops are written as loops, and the masked
## observation-precision array is carried as a stored 2x2 array.  No MATLAB was
## executed anywhere in this work: the evidence class is SOURCE-DERIVED from an
## independent transcription, one step stronger than package-to-package
## agreement and one step weaker than an authors-source fixture.
##
## The package is used for exactly one thing: simulating the two driving paths.
## A path is an input to a filter, not an oracle for one.
##
## COEFFICIENT PHASE.  The transcription can evaluate the seasonal
## coefficients on either convention (see tc_phase() below).  The SHIPPED
## fixtures use "state_time", because that is the convention every acir
## constructor uses.  The MATLAB stored arrays hold element 1 at state time and
## elements 2..N+1 one step ahead (:1046/:1052 against :1140/:1150); the
## coefficient fixture carries both columns so the documented difference stays
## visible, and the manifest records its measured size.
##
##   Rscript tools/fixtures/make-tc-fixtures.R <outdir> [producer] [producer_lib]
##
## <outdir> is required and must be OUTSIDE the package.  Fixtures are never
## regenerated in place: write to a scratch directory, diff the SHA-256 against
## the manifest, and replace the shipped files only as a deliberate maintainer
## step with new hashes and a NEWS entry.
## =========================================================================

args <- commandArgs(trailingOnly = TRUE)
if( !length(args) )
  stop("usage: Rscript make-tc-fixtures.R <outdir> [producer] [producer_lib]")
OUTDIR   <- normalizePath(args[[1L]], mustWork = FALSE)
PRODUCER <- if( length(args) >= 2L ) args[[2L]] else "acir"
if( length(args) >= 3L ) .libPaths(c(args[[3L]], .libPaths()))

FIX <- file.path(OUTDIR, "fixtures")
dir.create(FIX, showWarnings = FALSE, recursive = TRUE)

suppressMessages(requireNamespace(PRODUCER, quietly = TRUE))
P <- asNamespace(PRODUCER)

sha256 <- function(f)
  sub(" .*$", "", system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE))

## The pinned driving signal for case A: the shared 4001-point path that the
## scalar-partition fixtures are also graded on.  Located relative to this
## script so the generator does not depend on the working directory.
SIGNAL_SHA <- "77412008549a2980af8bc843c3548a673a4c2df952a40aea12ef2d30d8dec250"
here <- tryCatch(dirname(normalizePath(sys.frame(1L)$ofile)),
                 error = function(e) NA_character_)
pkg_root <- if( is.na(here) ) getwd() else normalizePath(file.path(here, "..", ".."))
SIGNAL <- file.path(pkg_root, "tests", "testthat", "fixtures", "oracles",
                    "enso6_partition_signal.csv")


## -------------------------------------------------------------------------
## Parameters (ENSO_model_cond_ACI_T_C_unobs.m, line references inline)
## -------------------------------------------------------------------------
tc_params <- function() {
  fct     <- 0.65                          # :913  factor
  b_0     <- 2.5                           # :920
  mu      <- 0.5                           # :922
  alpha_2 <- 0.125 * fct                   # :927
  alpha_1 <- alpha_2 / 2 * fct             # :928
  gamma_C <- 0.75 * fct                    # :938
  gamma_E <- 0.75 * fct                    # :940
  list(factor = fct, b_0 = b_0, mu = mu, alpha_1 = alpha_1, alpha_2 = alpha_2,
       r = 0.25 * fct,                     # :936
       gamma_C = gamma_C, gamma_E = gamma_E,
       r_C = gamma_C * b_0 * mu / 2,       # :942
       r_E = 3 * gamma_E * b_0 * mu / 2,   # :944
       zeta_C = gamma_C * b_0 * mu / 2,    # :946
       zeta_E = gamma_E * b_0 * mu / 2,    # :948
       C_u = 0.03 * fct,                   # :951
       lambda = 2 / 60,                    # :955
       m = 2,                              # :957
       sigma_C = 0.04 * sqrt(fct),         # :962
       sigma_E = sqrt(5) * 1e-2 * sqrt(fct), # :973
       dt = 0.005,                         # :894
       k_dt = 100)                         # :903  k_dt = 0.5/dt
}

## MATLAB stored-coefficient phase.  L_y(1) and f_x(:,1) are written at
## :1052/:1046 with the literal argument 0*dt, the state time of index 1.
## Inside the loop, :1150/:1141 use j*dt at array index j, whose state time is
## (j-1)*dt, so every element after the first carries a one-step-ahead phase.
tc_phase <- function(j, dt, convention) {
  if( identical(convention, "state_time") ) return((j - 1) * dt)
  ifelse(j == 1L, 0, j * dt)
}

## The full coefficient set of the reduced inference model.
##   observed x = (T_E, I)                (proxy reduction, :999-1021)
##   hidden   y = T_C
##   L_x = [-zeta_E; 0], constant         (:1024)
##   L_y(t)  = r_C - c_1(t, 0)            (:1052, :1150)   zeroth-order Taylor
##   f_y(t)  from prescribed u, h_W, tau  (:1053, :1151)
##   S_x(t)  = diag(sigma_E, sigma_I(I))  (:1046, :1144)
##   S_y     = sigma_C                    (:1037)
tc_coefficients <- function(path, p, gamma_hW, phase) {
  dt <- p$dt
  N1 <- nrow(path)
  ph <- tc_phase(seq_len(N1), dt, phase)
  u <- path$u; hW <- path$hW; TE <- path$TE; tau <- path$tau; I <- path$I

  beta0 <- (1 + (1 - I / 5)) * 0.15 * sqrt(p$factor)
  ## c_1(t, 0): the zeroth-order substitution, 25*(0 + 0.75/7.5)^2 + 0.9.
  c1_0 <- (25 * (0 + 0.75 / 7.5)^2 + 0.9) *
          (1 + 0.3 * sin(ph * 2 * pi / 6 - pi / 6)) * p$factor
  c2 <- 1.4 * p$factor * (1 + 0.3 * sin(ph * 2 * pi / 6 + 2 * pi / 6) +
                              0.25 * sin(2 * ph * 2 * pi / 6 + 2 * pi / 6))

  f_y <- p$zeta_C * TE + I / 5 * p$factor * u + p$C_u + (0.8 * beta0) * tau
  if( identical(gamma_hW, "intended") ) f_y <- f_y + p$gamma_C * hW

  f_x <- cbind(TE = (p$r_E - c2) * TE + p$gamma_E * hW + beta0 * tau,
               I  = -p$lambda * (I - p$m))
  ## S_x(:,:,j) = [sigma_E 0; 0 sqrt(lambda*(4-I)*I)]   :1144
  sigma_I <- sqrt(pmax(p$lambda * (4 - I) * I, 0))
  ## The variance floor acir applies so the observed Gram stays non-singular
  ## at the natural boundaries I in {0, 4}.  Transcribed here from
  ## aci_enso_model()'s documented convention, not read from the package.
  sigma_I_floored <- local({
    Ic <- pmin(pmax(I, 0), 4)
    sqrt(pmax(p$lambda * Ic * (4 - Ic), 0) + 1e-3 * p$lambda)
  })
  list(t = path$t, N1 = N1, dt = dt, L_x = c(-p$zeta_E, 0),
       L_y = p$r_C - c1_0, f_y = f_y, f_x = f_x,
       sigma_E = p$sigma_E, sigma_I = sigma_I,
       sigma_I_floored = sigma_I_floored,
       S_yoS_y = p$sigma_C^2,               # :1216
       x_TE = path$TE, x_I = path$I)
}

## Observed-noise Gram inverse.  The target is T_E (:1209), so the masked
## precision is diag(1/sigma_E^2, 0) at every slice: the uniform convention.
tc_gram_inverse <- function(co) {
  G <- array(0, c(2, 2, co$N1))
  G[1, 1, ] <- 1 / co$sigma_E^2
  G
}


## -------------------------------------------------------------------------
## Independent transcription: filter, smoother, ACI
## -------------------------------------------------------------------------
matlab_filter <- function(co, Ginv, filter_mean1, filter_cov1) {
  N1 <- co$N1; dt <- co$dt
  L_x <- matrix(co$L_x, nrow = 2L, ncol = 1L)             # :1024
  S_yoS_y <- co$S_yoS_y                                   # :1216
  filter_mean <- numeric(N1); filter_cov <- numeric(N1)
  filter_mean[1L] <- filter_mean1                         # :1223
  filter_cov[1L]  <- filter_cov1                          # :1229
  x <- cbind(co$x_TE, co$x_I)
  for( j in 2:N1 ) {
    ## :1239  dx = [T_E(j)-T_E(j-1); I(j)-I(j-1)]
    dx <- matrix(c(x[j, 1L] - x[j - 1L, 1L], x[j, 2L] - x[j - 1L, 2L]), nrow = 2L)
    Gi <- matrix(Ginv[, , j - 1L], 2L, 2L)
    ## :1247-1249
    innov <- dx - (L_x * filter_mean[j - 1L] +
                     matrix(co$f_x[j - 1L, ], nrow = 2L)) * dt
    filter_mean[j] <- filter_mean[j - 1L] +
      (co$L_y[j - 1L] * filter_mean[j - 1L] + co$f_y[j - 1L]) * dt +
      as.numeric(filter_cov[j - 1L] * t(L_x) %*% Gi %*% innov)
    ## :1250-1252
    filter_cov[j] <- filter_cov[j - 1L] +
      (co$L_y[j - 1L] * filter_cov[j - 1L] +
         filter_cov[j - 1L] * co$L_y[j - 1L] + S_yoS_y) * dt -
      as.numeric(filter_cov[j - 1L] * t(L_x) %*% Gi %*% L_x *
                   filter_cov[j - 1L]) * dt
  }
  list(mean = filter_mean, cov = filter_cov)
}

matlab_smoother <- function(co, Ginv, filt) {
  N1 <- co$N1; dt <- co$dt
  L_x <- matrix(co$L_x, nrow = 2L, ncol = 1L)
  S_yoS_y <- co$S_yoS_y
  smoother_mean <- numeric(N1); smoother_cov <- numeric(N1)
  smoother_mean[N1] <- filt$mean[N1]                      # :1270
  smoother_cov[N1]  <- filt$cov[N1]                       # :1271
  muT <- smoother_mean[N1]; RT <- smoother_cov[N1]
  for( j in (N1 - 1L):1L ) {
    A_j <- co$L_y[j]                                      # :1306
    B_j <- S_yoS_y                                        # :1307
    ## :1315
    mu <- muT - (co$L_y[j] * muT + co$f_y[j] -
                   B_j / filt$cov[j] * (filt$mean[j] - muT)) * dt
    ## :1316
    R <- RT - ((A_j + B_j / filt$cov[j]) * RT +
                 RT * (A_j + B_j / filt$cov[j]) - B_j) * dt
    smoother_mean[j] <- mu; smoother_cov[j] <- R
    muT <- mu; RT <- R
  }
  list(mean = smoother_mean, cov = smoother_cov)
}

## :1394-1397
matlab_aci <- function(filt, smoo) {
  signal <- 0.5 * (smoo$mean - filt$mean)^2 / filt$cov
  ratio  <- smoo$cov / filt$cov
  dispersion <- 0.5 * (-log(ratio) + ratio - 1)
  list(signal = signal, dispersion = dispersion, total = signal + dispersion)
}

matlab_run <- function(path, p, gamma_hW, phase) {
  co <- tc_coefficients(path, p, gamma_hW, phase)
  Ginv <- tc_gram_inverse(co)
  filt <- matlab_filter(co, Ginv, filter_mean1 = path$TC[1L], filter_cov1 = 0.1)
  smoo <- matlab_smoother(co, Ginv, filt)
  list(co = co, filter = filt, smoother = smoo, aci = matlab_aci(filt, smoo))
}


## -------------------------------------------------------------------------
## 1. Driving paths
## -------------------------------------------------------------------------
p <- tc_params()

stopifnot(file.exists(SIGNAL))
cat("case A signal sha256:", sha256(SIGNAL), "\n")
if( !identical(sha256(SIGNAL), SIGNAL_SHA) )
  stop("case A signal does not match its pinned SHA-256")
sigA <- utils::read.csv(SIGNAL)
pathA <- data.frame(t = sigA$t, u = sigA$u, hW = sigA$hW, TC = sigA$TC,
                    TE = sigA$TE, tau = sigA$tau, I = sigA$I)

## Case B: a longer path, so that a 14-model-year analysis window with the
## script's two-year buffers actually exists.  Case A's 4000 steps are 40 model
## months and hold no such window.
simulate_six <- function(seed, T_end) {
  m <- get("aci_enso_model", envir = P)(variant = "aci_code",
                                     hidden = c("u", "hW", "tau"))
  s <- stats::simulate(m, seed = seed, T = T_end, dt = p$dt)
  x <- as.matrix(s$obs$x); y <- as.matrix(s$hidden)
  data.frame(t = s$obs$t, u = y[, 1L], hW = y[, 2L], TC = x[, 1L],
             TE = x[, 2L], tau = y[, 3L], I = x[, 3L])
}
pathB <- simulate_six(4242L, 110)

## Case A regenerates from its own seed too, which is what lets the signal be
## pinned by bytes rather than by a seed.
regenA <- simulate_six(42L, 20)
cat("case A regenerates from seed 42:",
    max(abs(as.matrix(regenA) - as.matrix(pathA))), "\n")

## MATLAB window arithmetic (:1351-1354) with sim_year_start = 3 and
## ACI_period_years = 14, 1-based into the full record.
idxB <- local({
  a <- 3 * 12 * p$k_dt + 1
  a:(a + 14 * 12 * p$k_dt)
})
stopifnot(max(idxB) <= nrow(pathB))
pathB_window <- pathB[idxB, ]
tmp <- file.path(OUTDIR, "tc_path_caseB_window.csv")
utils::write.csv(pathB_window, tmp, row.names = FALSE)
cat("case B window sha256:", sha256(tmp), " rows:", nrow(pathB_window),
    " t:", min(pathB_window$t), "-", max(pathB_window$t), "\n")


## -------------------------------------------------------------------------
## 2. Coefficient fixture (case A, both phase conventions)
## -------------------------------------------------------------------------
coM <- tc_coefficients(pathA, p, "intended", "matlab")
coS <- tc_coefficients(pathA, p, "intended", "state_time")
coL <- tc_coefficients(pathA, p, "literal",  "state_time")
utils::write.csv(
  data.frame(t = coS$t,
             L_y_state_time      = coS$L_y,
             L_y_matlab_phase    = coM$L_y,
             f_x_TE_state_time   = coS$f_x[, "TE"],
             f_x_TE_matlab_phase = coM$f_x[, "TE"],
             f_y_intended        = coS$f_y,
             f_y_literal         = coL$f_y,
             f_x_I               = coS$f_x[, "I"],
             sigma_I_matlab      = coS$sigma_I,
             sigma_I_floored     = coS$sigma_I_floored),
  file.path(FIX, "tc_coefficients_caseA.csv"), row.names = FALSE)


## -------------------------------------------------------------------------
## 3. Filter / smoother / ACI, both D1 variants, on the shipped phase
## -------------------------------------------------------------------------
dump <- function(run, idx, file)
  utils::write.csv(
    data.frame(t = run$co$t[idx],
               filter_mean    = run$filter$mean[idx],
               filter_cov     = run$filter$cov[idx],
               smoother_mean  = run$smoother$mean[idx],
               smoother_cov   = run$smoother$cov[idx],
               aci_signal     = run$aci$signal[idx],
               aci_dispersion = run$aci$dispersion[idx],
               aci_total      = run$aci$total[idx]),
    file.path(FIX, file), row.names = FALSE)

rAI <- matlab_run(pathA, p, "intended", "state_time")
rAL <- matlab_run(pathA, p, "literal",  "state_time")
rBI <- matlab_run(pathB, p, "intended", "state_time")
rBL <- matlab_run(pathB, p, "literal",  "state_time")

dump(rAI, seq_len(nrow(pathA)), "tc_outputs_caseA_intended.csv")
dump(rAL, seq_len(nrow(pathA)), "tc_outputs_caseA_literal.csv")
## The case B window is 16801 rows.  Pin a monthly subsample (k_dt = 100) plus
## the first 200 consecutive steps, which is what a step-level regression
## needs.  Both are read-outs of a filter that ran from t = 0 over the whole
## 22001-point record, exactly as the script's own analysis window is.
monthly <- idxB[seq(1L, length(idxB), by = p$k_dt)]
headwin <- idxB[seq_len(200L)]
dump(rBI, monthly, "tc_outputs_caseB_window_monthly_intended.csv")
dump(rBL, monthly, "tc_outputs_caseB_window_monthly_literal.csv")
dump(rBI, headwin, "tc_outputs_caseB_window_head_intended.csv")
dump(rBL, headwin, "tc_outputs_caseB_window_head_literal.csv")


## -------------------------------------------------------------------------
## 4. The D1 divergence, pinned as a measurement rather than as a default
## -------------------------------------------------------------------------
## The reference script's f_y omits gamma_C * h_W (:1053, :1151).  acir
## includes it and reproduces the omission under matlab_defect_compat = TRUE.
## The covariances are bit-identical between the two, because f_y enters
## neither the Riccati equation nor the backward covariance ODE; that is the
## sharpest available guard that the flag touches only the mean channel.
divergence <- function(case, a, b, idx) {
  q <- list(filter_mean = c("filter", "mean"), filter_cov = c("filter", "cov"),
            smoother_mean = c("smoother", "mean"),
            smoother_cov = c("smoother", "cov"),
            aci_signal = c("aci", "signal"),
            aci_dispersion = c("aci", "dispersion"),
            aci_total = c("aci", "total"))
  do.call(rbind, lapply(names(q), function(nm) {
    k <- q[[nm]]
    va <- a[[k[1L]]][[k[2L]]][idx]; vb <- b[[k[1L]]][[k[2L]]][idx]
    data.frame(case = case, quantity = nm, n = length(idx),
               identical = identical(va, vb),
               max_abs = max(abs(va - vb)),
               rmsd = sqrt(mean((va - vb)^2)))
  }))
}
## Trapezoidal time integral of the ACI series, the script's reported figure.
trapz <- function(t, v) sum(diff(t) * (utils::head(v, -1L) + v[-1L]) / 2)
D1 <- rbind(divergence("A whole record", rAI, rAL, seq_len(nrow(pathA))),
            divergence("B 14-year window", rBI, rBL, idxB))
iA <- seq_len(nrow(pathA))
D1 <- rbind(D1, data.frame(
  case = c("A whole record", "B 14-year window"),
  quantity = "aci_total_integrated",
  n = c(length(iA), length(idxB)),
  identical = FALSE,
  ## max_abs holds the intended total, rmsd the literal one: two numbers with
  ## no natural home in a difference table, kept here so the headline D1
  ## figure is pinned in the same file as the per-quantity differences.
  max_abs = c(trapz(rAI$co$t[iA], rAI$aci$total[iA]),
              trapz(rBI$co$t[idxB], rBI$aci$total[idxB])),
  rmsd = c(trapz(rAL$co$t[iA], rAL$aci$total[iA]),
           trapz(rBL$co$t[idxB], rBL$aci$total[idxB]))))
utils::write.csv(D1, file.path(FIX, "tc_d1_divergence.csv"), row.names = FALSE)
cat("\nD1 divergence (max_abs / rmsd; the integrated rows hold",
    "intended and literal totals):\n")
print(D1, row.names = FALSE, digits = 6)


## -------------------------------------------------------------------------
## 5. Hashes
## -------------------------------------------------------------------------
cat("\nshipped fixture hashes:\n")
for( f in sort(list.files(FIX, pattern = "\\.csv$", full.names = TRUE)) )
  cat(sprintf("  %-46s %8.1f KB  sha256 %s  md5 %s\n", basename(f),
              file.size(f) / 1024, sha256(f), unname(tools::md5sum(f))))
cat("\nDONE\n")
