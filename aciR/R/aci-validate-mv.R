# -- contracts for the vector-valued components schema -------------------------
#
# The scalar contract boundary in aci-validate.R checks lengths and signs. The
# vector one has to check conformability as well: a components list is a set of
# matrices that must compose, and a mismatch produces either a low-level error
# from the first multiplication or, worse, a silently recycled one.
#
# The dimensions are inferred from the noise Grammians rather than declared,
# because those are the two entries whose shape is unambiguous: S_xoS_x is
# square in the observed dimension and S_yoS_y is square in the unobserved one.
# Everything else is checked against them.

#' Is this a vector-valued components list?
#'
#' The discriminator is the shape of the latent-noise covariance: a matrix
#' means the unobserved component is a vector, a bare number means it is
#' scalar. A one-by-one matrix therefore routes to the vector path, which is
#' deliberate -- someone who wrote their system in matrices gets the matrix
#' recursions, and the two agree.
#'
#' @param comp A components list.
#'
#' @returns `TRUE` for the vector schema, `FALSE` for the scalar one.
#'
#' @noRd
#' @keywords internal
.aci_is_mv <- function(comp) {
  is.list(comp) &&
    (is.matrix(comp$S_yoS_y) || is.matrix(comp$S_xoS_x) ||
       is.matrix(comp$L_y) || length(dim(comp$L_y)) == 3L)
}

#' Check a coefficient's shape against the state dimensions
#'
#' @param value The coefficient.
#' @param rows,cols Integer scalars. Required dimensions.
#' @param n Integer scalar. The number of time steps.
#' @param label Character. The entry name, used in the error message.
#' @param name Character. The components argument name.
#'
#' @returns `TRUE`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_mv_coef <- function(value, rows, cols, n, label, name) {
  d <- dim(value)
  ok <- is.numeric(value) && !is.null(d) &&
    ((length(d) == 2L && d[1L] == rows && d[2L] == cols) ||
       (length(d) == 3L && d[1L] == rows && d[2L] == cols && d[3L] == n))
  if (!ok) {
    stop(
      sprintf(
        paste0(
          "`%s$%s` must be a %d by %d numeric matrix, constant in time, or a ",
          "%d by %d by %d array with time in the last margin; it is %s."
        ),
        name, label, rows, cols, rows, cols, n, .aci_describe(value)
      ),
      call. = FALSE
    )
  }
  if (!all(is.finite(value))) {
    stop(
      sprintf("`%s$%s` must be finite throughout.", name, label),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check a drift term's shape against the state dimensions
#'
#' @param value The drift term.
#' @param rows Integer scalar. The state dimension.
#' @param n Integer scalar. The number of time steps.
#' @param label,name Character. Entry and argument names.
#'
#' @returns `TRUE`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_mv_drift <- function(value, rows, n, label, name) {
  ok <- is.numeric(value) &&
    ((is.null(dim(value)) && length(value) == rows) ||
       (is.matrix(value) && nrow(value) == rows && ncol(value) == n))
  if (!ok) {
    stop(
      sprintf(
        paste0(
          "`%s$%s` must be a numeric vector of length %d, constant in time, ",
          "or a %d by %d matrix with time in columns; it is %s."
        ),
        name, label, rows, rows, n, .aci_describe(value)
      ),
      call. = FALSE
    )
  }
  if (!all(is.finite(value))) {
    stop(
      sprintf("`%s$%s` must be finite throughout.", name, label),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate a vector-valued components list and derive its inverse Grammian
#'
#' Checks conformability, symmetry and admissibility, then returns the list
#' with `S_xoS_x_inv` filled in. The inverse is computed once here rather than
#' at every step of the recursion, and it may also be SUPPLIED by the caller:
#' the conditional causal-inference construction works by handing in an inverse
#' that is deliberately not the inverse of anything, so that the non-target
#' observations carry no weight.
#'
#' @param comp The components list.
#' @param x The observed signal, a matrix with one row per observed component.
#' @param name The argument name, used in error messages.
#'
#' @returns The components list, with `S_xoS_x_inv` present.
#'
#' @noRd
#' @keywords internal
.aci_check_components_mv <- function(comp, x, name = "comp") {
  required <- c("L_x", "f_x", "L_y", "f_y", "S_xoS_x", "S_yoS_y", "S_yoS_x")
  absent <- setdiff(required, names(comp))
  if (length(absent) > 0L) {
    stop(
      sprintf(
        "`%s` is missing the component%s %s.",
        name, if (length(absent) > 1L) "s" else "",
        paste(sprintf("`%s`", absent), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  n <- ncol(x)
  n_x <- nrow(x)
  n_y <- nrow(as.matrix(.aci_slice(comp$S_yoS_y, 1L)))

  s_xx <- .aci_slice(comp$S_xoS_x, 1L)
  if (!is.matrix(s_xx) || nrow(s_xx) != n_x || ncol(s_xx) != n_x) {
    stop(
      sprintf(
        paste0(
          "`%s$S_xoS_x` must be a %d by %d matrix, square in the observed ",
          "dimension implied by `x`; it is %s."
        ),
        name, n_x, n_x, .aci_describe(s_xx)
      ),
      call. = FALSE
    )
  }

  .aci_check_mv_coef(comp$L_x, n_x, n_y, n, "L_x", name)
  .aci_check_mv_coef(comp$L_y, n_y, n_y, n, "L_y", name)
  .aci_check_mv_coef(comp$S_yoS_y, n_y, n_y, n, "S_yoS_y", name)
  .aci_check_mv_coef(comp$S_yoS_x, n_y, n_x, n, "S_yoS_x", name)
  .aci_check_mv_drift(comp$f_x, n_x, n, "f_x", name)
  .aci_check_mv_drift(comp$f_y, n_y, n, "f_y", name)

  # The latent-noise covariance is a covariance, so it must be symmetric. An
  # asymmetric one is a transposition error rather than an exotic model, and
  # the recursion would propagate it silently.
  for (entry in c("S_yoS_y", "S_xoS_x")) {
    m <- .aci_slice(comp[[entry]], 1L)
    if (max(abs(m - t(m))) > 1e-10 * max(1, max(abs(m)))) {
      stop(
        sprintf(
          "`%s$%s` must be symmetric; it is not.", name, entry
        ),
        call. = FALSE
      )
    }
  }

  # The observation-noise covariance is inverted by the filter, so it must be
  # positive definite and not merely non-negative -- at EVERY step it is
  # supplied for, not merely the first.
  #
  # Inverting only the first slice of a time-varying covariance is a silent
  # wrong answer rather than an error: the recursion runs, the result looks
  # plausible, and the variation the caller supplied is discarded. This code
  # did exactly that until it was probed by perturbing a late slice and
  # observing the output not move at all.
  if (is.null(comp$S_xoS_x_inv)) {
    varying <- length(dim(comp$S_xoS_x)) == 3L
    steps <- if (varying) seq_len(n) else 1L
    inverse <- if (varying) array(0, c(n_x, n_x, n)) else NULL
    for (j in steps) {
      slice <- .aci_slice(comp$S_xoS_x, j)
      factor <- .aci_chol(slice)
      if (is.null(factor)) {
        stop(
          sprintf(
            paste0(
              "`%s$S_xoS_x` must be symmetric positive definite at every ",
              "step: the filter inverts it. It is not at step %d. Supply ",
              "`%s$S_xoS_x_inv` directly if a weighting that is not the ",
              "inverse of a covariance is intended, as the conditional ",
              "construction requires."
            ),
            name, j, name
          ),
          call. = FALSE
        )
      }
      if (varying) {
        inverse[, , j] <- chol2inv(factor)
      } else {
        inverse <- chol2inv(factor)
      }
    }
    comp$S_xoS_x_inv <- inverse
  } else {
    d <- dim(comp$S_xoS_x_inv)
    ok <- is.numeric(comp$S_xoS_x_inv) && !is.null(d) &&
      d[1L] == n_x && d[2L] == n_x &&
      (length(d) == 2L || (length(d) == 3L && d[3L] == n)) &&
      all(is.finite(comp$S_xoS_x_inv))
    if (!ok) {
      stop(
        sprintf(
          paste0(
            "`%s$S_xoS_x_inv` must be a finite %d by %d matrix, constant in ",
            "time, or a %d by %d by %d array with time in the last margin; ",
            "it is %s."
          ),
          name, n_x, n_x, n_x, n_x, n, .aci_describe(comp$S_xoS_x_inv)
        ),
        call. = FALSE
      )
    }
  }
  comp
}

#' Validate an observed signal for the vector path
#'
#' @param x The observed signal.
#' @param name The argument name.
#'
#' @returns `x` as a matrix with one row per observed component.
#'
#' @noRd
#' @keywords internal
.aci_check_signal_mv <- function(x, name = "x") {
  if (is.numeric(x) && is.null(dim(x))) {
    x <- matrix(x, nrow = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) < 2L) {
    stop(
      sprintf(
        paste0(
          "`%s` must be a numeric matrix with one row per observed component ",
          "and at least two columns, one per time step; it is %s."
        ),
        name, .aci_describe(x)
      ),
      call. = FALSE
    )
  }
  bad <- which(!is.finite(x))
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "`%s` must be complete and finite; it holds %s at position %d.",
        name, .aci_bad_kind(x[bad[1L]]), bad[1L]
      ),
      call. = FALSE
    )
  }
  x
}
