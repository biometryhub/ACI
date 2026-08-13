# -- conditional assimilative causal inference ---------------------------------
#
# The metric as it stands asks what the observed process, taken whole, says
# about the unobserved one. Conditional ACI asks a sharper question: what does
# ONE observed component say about the unobserved state, given that the others
# are also being watched?
#
# The construction is the one the reference implementation uses, and it is
# simpler than it first looks. The filter weights each observed component by
# the inverse of the observation-noise covariance. Inflating a component's
# observational uncertainty without bound sends its weight to zero, which is
# the same as declining to condition on it. So the conditional question is
# asked by handing the filter an inverse Grammian that is zero everywhere
# except on the target block -- an object that is deliberately not the inverse
# of any covariance, which is why the components schema admits
# `S_xoS_x_inv` directly.
#
# What this changes is the ESTIMAND, not the arithmetic. The result is a
# statement about a different quantity, and the assumptions vignette says so;
# it is not a variance-reduction trick applied to the same one.

#' Condition the causal question on a subset of the observed components
#'
#' Rewrites a vector-valued components list so that the causal-information
#' metric computed from it measures what the *target* observed components say
#' about the unobserved state, rather than what the whole observed process
#' says.
#'
#' The non-target components are not removed. They continue to drive the
#' system, and their own dynamics are unchanged; what changes is that the
#' filter stops treating them as observations to be assimilated. This is what
#' distinguishes the conditional question from simply running the method on a
#' shorter signal, in which the non-target components would not be present at
#' all.
#'
#' The construction inflates the observational uncertainty of the non-target
#' components without bound, which sends their weight in the filter to zero. It
#' is implemented by supplying the filter with an inverse noise Grammian
#' supported only on the target block. That object is not the inverse of a
#' covariance, and it is not required to be.
#'
#' A component whose observation noise is exactly zero cannot be a target: the
#' construction needs its noise covariance to be invertible on the target
#' block. The reference implementation meets this by adding a small artificial
#' noise to such a component, which is a modelling decision and is reported
#' here as an error rather than made silently.
#'
#' @param comp A vector-valued conditional Gaussian components list; see
#'   [aci_components]. The scalar schema has only one observed component and so
#'   admits no conditional question.
#' @param target Integer or character vector. Which observed components the
#'   causal question is asked about. Characters are matched against the row
#'   names of the observation-noise Grammian.
#'
#' @returns A components list with `S_xoS_x_inv` supported on the target block,
#'   ready for [aci_filter()], [aci_smoother()] and [aci_metric()].
#'
#' @references
#' Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
#' *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
#'
#' @examples
#' # A two-component observed process; ask what only the first says about the
#' # unobserved state.
#' comp <- list(
#'   L_x = diag(2), f_x = c(0, 0), L_y = -diag(2), f_y = c(0, 0),
#'   S_xoS_x = diag(c(0.5, 0.8)), S_yoS_y = diag(2),
#'   S_yoS_x = matrix(0, 2, 2)
#' )
#' conditioned <- aci_conditional(comp, target = 1)
#' conditioned$S_xoS_x_inv
#'
#' @seealso [aci_components], [aci_filter()], [aci_metric()]
#' @export
aci_conditional <- function(comp, target) {
  if (!.aci_is_mv(comp)) {
    stop(
      paste0(
        "`comp` must be a vector-valued components list: a conditional ",
        "question needs more than one observed component, and the scalar ",
        "schema has exactly one."
      ),
      call. = FALSE
    )
  }
  s_xx <- .aci_slice(comp$S_xoS_x, 1L)
  if (!is.matrix(s_xx)) {
    stop("`comp$S_xoS_x` must be a matrix.", call. = FALSE)
  }
  n_x <- nrow(s_xx)

  if (is.character(target)) {
    matched <- match(target, rownames(s_xx))
    if (anyNA(matched)) {
      stop(
        sprintf(
          "`target` names no observed component: %s.",
          paste(sQuote(target[is.na(matched)]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    target <- matched
  }
  ok <- is.numeric(target) && length(target) >= 1L && !anyNA(target) &&
    all(target == trunc(target)) && all(target >= 1) && all(target <= n_x) &&
    !anyDuplicated(target)
  if (!ok) {
    stop(
      sprintf(
        paste0(
          "`target` must be distinct whole numbers between 1 and %d, ",
          "indexing the observed components."
        ),
        n_x
      ),
      call. = FALSE
    )
  }
  if (length(target) == n_x) {
    stop(
      paste0(
        "`target` names every observed component, so there is nothing to ",
        "condition on. Use the components list unchanged for the ",
        "unconditional question."
      ),
      call. = FALSE
    )
  }

  target <- as.integer(sort(target))

  # The masked weight inherits the time-variation of the covariance it is
  # derived from. Building it from the first step alone would freeze a weight
  # that should move -- the same silent failure the vector validator once had,
  # and it would bite exactly where it is hardest to notice, since the
  # conditional result has no unconditional counterpart to be compared against.
  varying <- length(dim(comp$S_xoS_x)) == 3L
  n <- if (varying) dim(comp$S_xoS_x)[3L] else 1L
  steps <- seq_len(n)
  inv <- if (varying) {
    array(0, c(n_x, n_x, n))
  } else {
    matrix(0, n_x, n_x, dimnames = dimnames(s_xx))
  }

  for (j in steps) {
    block <- .aci_slice(comp$S_xoS_x, j)[target, target, drop = FALSE]
    factor <- .aci_chol(block)
    if (is.null(factor)) {
      stop(
        sprintf(
          paste0(
            "The observation-noise covariance of the target components is ",
            "not positive definite at step %d, so their contribution to the ",
            "filter cannot be weighted. A component with no observation ",
            "noise of its own cannot be a target until some is modelled for ",
            "it."
          ),
          j
        ),
        call. = FALSE
      )
    }
    if (varying) {
      inv[target, target, j] <- chol2inv(factor)
    } else {
      inv[target, target] <- chol2inv(factor)
    }
  }

  comp$S_xoS_x_inv <- inv
  comp$conditional_target <- target
  comp
}
