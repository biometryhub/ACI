################################################################################
## acir-package.R - maintainer's map of the package architecture
## ########################################################################## ##
##
## Layered design (each layer depends only on those above it):
##
##   aci-utils.R            conditions, SPD-safe linear algebra, obs contract
##   aci-model.R            stochastic_model / cgns_model classes + simulation
##   aci-model-library.R    benchmark constructors with per-object provenance
##   aci-assimilation.R     closed-form CGNS filter/smoother; Theorem 3 lag
##                          table; conditional ACI strategies
##   aci-kernels-scalar.R   compiled scalar CGNS kernels and the dyad fast path
##   aci-kernels-matrix.R   compiled matrix CGNS kernels
##   aci-conditional.R      compiled conditional-ACI construction
##   aci-online-smoother.R  compiled Theorem 3 lag-table core
##   aci-cir.R              compiled streaming forward-CIR reducer
##   aci-core.R             Gaussian KL + aci() metric + forward CIR
##   aci-methods.R          base-graphics methods for the engine result classes
##
## Things to know before editing:
##   * The Theorem 3 lag table is exact for the EXPLICIT single-step
##     filter/smoother discretization; paths carry meta$stepper / meta$nsub
##     and lag_table() recomputes (with aci_warn_stepper) on mismatch.
##   * All errors/warnings are classed ("aci_error_*", "aci_warn_*").
##   * Golden alignment vs the published reference code lives in
##     tests/testthat/test-10-golden.R + helper-golden-p1.R.
##
## Scope: this release implements the closed-form ACI_code surface only. The
## ensemble (EnKBS) engine, model discovery, the extreme-event family, the
## formula front end, the applied workflows, the validation-diagnostics suite
## and backward CIR are held in reserve/ and are not part of the mainline.
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
#' @importFrom stats simulate rnorm
#' @name acir-imports
#' @keywords internal
NULL
