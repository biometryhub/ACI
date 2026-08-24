################################################################################
## aci-package.R - maintainer's map of the package architecture
## ########################################################################## ##
##
## Layered design (each layer depends only on those above it):
##
##   utils.R          conditions, SPD-safe linear algebra, obs contract
##   model_classes.R  stochastic_model / cgns_model classes + simulation
##   benchmark_models.R benchmark constructors with per-object provenance
##   benchmark_models_fbcir.R andreou2026cir multiscale model + distinct
##                    repository-only layered topographic benchmark
##   assimilation.R   closed-form CGNS filter/smoother; Theorem 3 lag
##                    table; conditional ACI strategies; localization
##   ensemble.R       EnKBF / EnKBS and forward ensemble lag table
##                    (derivative-free; noise-reuse invariant)
##   causal_metrics.R Gaussian KL + aci() metric + forward/backward CIRs;
##                    tidy CIR read-out
##   discovery.R      dictionaries, causation entropy, FFBS sampling,
##                    constrained/energy-joint MLE, learn_model()
##   validation_diagnostics.R nil checks, cross-engine validation, LR test,
##                    one-call check, and perfect-model OSSE
##   extremes.R       moser2026extremes events, onsets, composites, features,
##                    clustering
##   formula_interface.R aci_fit() formula front-end, generics, plot methods
##   applied_workflows.R reusable observation preprocessing and CGNS calibration
##
## Things to know before editing:
##   * The Theorem 3 lag table is exact for the EXPLICIT single-step
##     filter/smoother discretization; paths carry meta$stepper / meta$nsub
##     and lag_table() recomputes (with aci_warn_stepper) on mismatch.
##   * EnKBS must reuse the forward Wiener increments; enforced by
##     meta$noise_hash (aci_error_noise_mismatch).
##   * All errors/warnings are classed ("aci_error_*", "aci_warn_*").
##   * Golden alignment vs the published reference codes lives in
##     tests/testthat/test-10-golden.R + helper-golden-p1.R.
##
## Tests: tests/testthat/test-01..14; grouped runs stay under CI timeouts
## (01-05 | 06,08,10 | 07 | 09). tests/testthat.R runs everything.
##
################################################################################

#' @keywords internal
"_PACKAGE"

#' Imported generics and helpers
#'
#' The generics extended by this package's S3 methods, and the base helpers used
#' throughout, are imported here so that the generated namespace carries them.
#'
#' @importFrom graphics plot
#' @importFrom stats coef fitted nobs predict residuals simulate rnorm var sd
#' @name aci-imports
#' @keywords internal
NULL
