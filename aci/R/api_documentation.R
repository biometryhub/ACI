################################################################################
## This file contains documentation-only roxygen blocks.  It keeps related
## low-level functions on a small number of discoverable help topics and adds
## argument documentation to functions whose implementation comments carry the
## scientific detail.  There is deliberately no executable code here.
################################################################################


#' Model and observation contracts
#'
#' Constructors and coercion helpers for uniformly sampled observations,
#' general stochastic models, and conditional-Gaussian nonlinear systems.
#' `cgns_from_affine()` verifies that the supplied drifts are affine in the
#' hidden state before constructing a CGNS representation.
#'
#' @name model_contracts
#' @seealso [da_filter()], [stats::simulate()]
NULL

#' Data assimilation and finite-lag API
#'
#' `da_filter()` and `da_smooth()` reconstruct hidden states. `lag_table()`
#' stores the finite-lag divergences used by the CIR estimators. The `lt_*`
#' helpers access a table without depending on its storage representation.
#' `nontarget()` describes conditional ACI masking; `reduce_nontarget()` is the
#' prescribed-forcing reduction used when requested. The historical
#' `lt_tail_bound()` name is retained for compatibility, but its value is a
#' heuristic tail estimate, not a certified mathematical error bound. A lag
#' table uses the complete online smoother of andreou2026cir Theorem 3
#' (Appendix G.1) as its reference. That
#' reference costs O(N) time-point work; table construction then costs work
#' proportional to the retained lag cells, with O(N^2) cells for a full table
#' in the worst case.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother with
#' closed-form solutions and information-theoretic lag selection for
#' conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
#' **36**(4), 71. \doi{10.1007/s00332-026-10271-x}
#'
#' @name assimilation_api
#' @seealso [aci()], [forward_cir()]
NULL

#' Assimilative causal metrics and influence ranges
#'
#' Gaussian relative entropy is oriented as smoother relative to filter. A
#' normal `aci()` call uses the supplied-code backward-ODE headline smoother,
#' including its correlated-noise correction, independently of `keep`.
#' `lag_table()` and `aci(table = ...)` instead use the complete online
#' Theorem 3 smoother; their finite-grid diagonal can therefore differ from
#' headline ACI. `forward_cir()` and `backward_cir()` summarize the duration of
#' influence on the discrete time grid. A finite adaptive table is labelled
#' `objective_on_truncated_table`; its `tail_bound` field is a heuristic tail
#' estimate and must not be interpreted as a certified error bound. The
#' `l1_linf` estimators are ratios: the forward ratio and both exact forms are
#' integrated with composite Simpson, following the ACI reference code, while
#' the backward ratio uses the plain Appendix G.3 L1 grid-function sum
#' (andreou2026cir eq. G.14), following the FBCIR code's active line. A
#' backward result is anchored at `T`, but its computable lagged grid and
#' physical interval end at `T - dt`; its
#' bound label records this convention.
#' For a general stochastic model, `aci(engine = "ensemble", keep = "table")`
#' evaluates the andreou2026cir forward lagged-posterior family using repeated
#' jiang2026enkbs EnKBS passes; `forward_cir()` then uses the same estimators.
#' Ensemble backward CIR is not exposed because jiang2026enkbs identifies it as
#' future work.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#' @name causal_metrics
#' @seealso [lag_table()], [cir()]
NULL

#' Ensemble assimilation API
#'
#' Ensemble Kalman-Bucy filtering and smoothing following the independent-noise
#' jiang2026enkbs recursion. `enkbs()` must receive the exact forward Wiener
#' increments returned by `enkbf()`. Shared/correlated observation and signal
#' noise is not accepted by this engine; use the closed-form CGNS engine for
#' that case. `ensemble_lag_table()` constructs andreou2026cir's forward
#' finite-lag posterior family by repeating the jiang2026enkbs backward pass at
#' every observation horizon. It has quadratic time-point work/storage and uses
#' full-dimensional Gaussian moment KL, so `m` must exceed the hidden dimension.
#' This route is graded at machine precision against a literal R port of the
#' published jiang2026enkbs dyad experiment, driven by identical increments
#' (see `tests/testthat/helper-golden-p3.R`).
#'
#' @references
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#'
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#' @name ensemble_api
#' @seealso [da_filter()], [localization_spec()]
NULL

#' Formula interface
#'
#' `aci_fit()` provides a modelling front end for a latent-cause formula.
#' `cir()` extracts a forward or backward causal influence range from a fitted
#' object. Partial-observation learning is experimental and is labelled as such
#' in fitted objects and summaries.
#'
#' @name formula_interface
#' @seealso [learn_model()], [aci()]
NULL

#' Validation diagnostics
#'
#' `nil_causality_check()` is a structural zero-coupling check.
#' `cross_validate()` compares ensemble moments with a compatible closed-form
#' CGNS result. These are diagnostics, not substitutes for checking model
#' assumptions or source provenance.
#'
#' @name validation_diagnostics
#' @seealso [nil_surrogate_test()]
NULL

#' Model discovery helpers
#'
#' Evaluate a function library, apply a structure threshold, or request the
#' default feature names used by the event classifiers.
#'
#' @name discovery_helpers
#' @seealso [cgns_library()], [learn_model()]
NULL


################################################################################
# References
################################################################################

#' References
#'
#' The published sources this package implements, and the shorthand keys used
#' to cite them in help pages and source comments. A key stands in for the work
#' itself, so `andreou2026aci eq. 7` means equation 7 of the first entry below,
#' and `moser2026extremes Appendix A` means Appendix A of the fourth.
#'
#' The `andreou2026aci`, `andreou2026cir`, and `jiang2026enkbs` papers have
#' associated MATLAB code, which the R implementation uses wherever possible.
#' Material attributed to `moser2026extremes` is implemented from the equations
#' alone and carries no claim of replaying any published reference run.
#'
#' \describe{
#'   \item{`andreou2026aci`}{Andreou, M., Chen, N. and Bollt, E. (2026).
#'     Assimilative causal inference. *Nature Communications* **17**, 1854.
#'     \doi{10.1038/s41467-026-68568-0} Code: <https://github.com/marandmath/ACI_code>.}
#'   \item{`andreou2026cir`}{Andreou, M. and Chen, N. (2026). Bridging
#'     prediction and attribution: identifying forward and backward causal
#'     influence ranges using assimilative causal inference. arXiv:2510.21889v2,
#'     4 August 2026. \doi{10.48550/arXiv.2510.21889}. Code: <https://github.com/marandmath/FBCIR_code>.}
#'   \item{`jiang2026enkbs`}{Jiang, Andreou, Reich and Chen (2026). A
#'     continuous-time ensemble Kalman-Bucy smoother for causal inference and
#'     model discovery. arXiv:2604.25157.
#'     \doi{10.48550/arXiv.2604.25157}. Code: <https://github.com/jiangzh67/EnKBS>.}
#'   \item{`moser2026extremes`}{Moser, Chen and Andreou (2026). Mechanisms and
#'     pathways of extreme events in partially-observed stochastic dynamical
#'     systems. arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}}
#'   \item{`andreou2026smoother`}{Andreou, M., Chen, N. and Li, Y. (2026). An
#'     adaptive online smoother with closed-form solutions and
#'     information-theoretic lag selection for conditional Gaussian nonlinear
#'     systems. *Journal of Nonlinear Science* **36**(4), 71. arXiv:2411.05870.
#'     \doi{10.1007/s00332-026-10271-x}.}
#'   \item{`chen2022enso`}{Chen, N., Fang, X. and Yu, J.-Y. (2022). A
#'     multiscale model for El Nino complexity. *npj Climate and Atmospheric
#'     Science* **5**, 16. arXiv:2104.07174.}
#' }
#'
#' @name aci_references
NULL
