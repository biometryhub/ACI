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
#' `aci_model_from_affine()` verifies that the supplied drifts are affine in
#' the hidden state before constructing a CGNS representation.
#'
#' @name model_contracts
#' @seealso [aci_filter()], [stats::simulate()]
NULL

#' Data assimilation and finite-lag API
#'
#' `aci_filter()`, `aci_smoother()` and `aci_online()` reconstruct hidden
#' states: from the record up to each time, from the whole record, and from
#' the record up to a fixed number of steps ahead of each time. `lag_table()`
#' stores the finite-lag divergences used by the CIR estimators. The `lt_*`
#' helpers access a table without depending on its storage representation.
#' `aci_conditional()` describes the conditional question, and
#' `aci_conditional_reduce()` carries out the model reduction its
#' `method = "reduce"` asks for. The historical `lt_tail_bound()` name is
#' retained for compatibility, but its value is a
#' heuristic tail estimate, not a certified mathematical error bound. A lag
#' table uses the complete online smoother of andreou2026cir Theorem 3
#' (Appendix G.1) as its reference. That
#' reference costs O(N) time-point work; table construction then costs work
#' proportional to the retained lag cells, with O(N^2) cells for a full table
#' in the worst case. `aci_online()` costs O(N) whatever the lag.
#'
#' @section Scheme:
#' Two discretizations of the same continuous-time smoothing problem are in
#' use, and `meta$scheme` on a path says which one produced it.
#' `"backward_ode_euler"` is [aci_smoother()]: the continuous backward smoothing
#' equations integrated with an Euler step. `"theorem3_discrete"` is
#' [aci_online()] and the lag table's reference smoother: the exact conditional
#' law of the hidden state given the observed increments on the sampling grid,
#' under the explicit single-step discretization. They agree only to first
#' order in the step, so [aci_online()] at `lag = Inf` does not reproduce
#' [aci_smoother()], and the gap grows with the length of the record rather than
#' settling to a constant. See the Scheme section of [aci_online()] for the
#' measured size on the packaged ENSO partition. `aci()` reports the scheme its
#' own result was built under in `meta$smoother_scheme`.
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
#' @seealso [aci()], [aci_range()]
NULL

#' Assimilative causal metrics and influence ranges
#'
#' Gaussian relative entropy is oriented as smoother relative to filter. A
#' normal `aci()` call uses the supplied-code backward-ODE headline smoother,
#' including its correlated-noise correction, independently of `keep`.
#' `lag_table()` and `aci(table = ...)` instead use the complete online
#' Theorem 3 smoother; their finite-grid diagonal can therefore differ from
#' headline ACI. `aci_range()` summarizes the duration of influence on the
#' discrete time grid. A finite adaptive table is labelled
#' `objective_on_truncated_table`; its `tail_bound` field is a heuristic tail
#' estimate and must not be interpreted as a certified error bound. The
#' `l1_linf` estimator is a ratio: the forward ratio and the exact form are
#' integrated with composite Simpson, following the ACI reference code.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#' @name causal_metrics
#' @seealso [lag_table()], [aci_range()]
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
