# aci-validate.R -- Shared argument validation for the aciR contract boundary.
#
# Every exported function routes its checks through the helpers below, so the
# general conditional Gaussian surface and the flagship dyad path enforce one
# contract and raise one class of message. The mathematical admissibility of a
# noise covariance -- a positive observation-noise covariance, a non-negative
# latent-noise covariance and a positive semidefinite joint covariance -- is
# stated once, in .aci_check_noise_covariance(), and is reused by the model
# constructor and by the components validator alike.
#
# The tolerances are scaled to the magnitude of the quantities being compared
# rather than fixed, because the noise Grammians of a physical system carry
# the units of that system and a fixed absolute tolerance would be strict on a
# small system and permissive on a large one.

# -- tolerance and description helpers ----------------------------------------

#' Round-off tolerance scaled to the magnitude of the compared quantities
#'
#' @param ... Numeric scalars whose magnitudes set the scale.
#'
#' @returns A positive numeric scalar.
#'
#' @noRd
#' @keywords internal
.aci_tol <- function(...) {
  sqrt(.Machine$double.eps) * max(1, abs(c(...)))
}

#' Describe a value for an error message
#'
#' Renders the type and shape of an offending value so a failed check can name
#' what it received rather than only what it wanted.
#'
#' @param value Any object.
#'
#' @returns A character string such as `"a character vector of length 3"`.
#'
#' @noRd
#' @keywords internal
.aci_describe <- function(value) {
  if (!is.null(dim(value))) {
    return(sprintf(
      "an object of class %s with dimension %s",
      class(value)[1L], paste(dim(value), collapse = " x ")
    ))
  }
  if (is.function(value)) {
    return("a function")
  }
  sprintf("%s of length %d", class(value)[1L], length(value))
}

#' Name the kind of a non-finite value
#'
#' @param value A length-one numeric.
#'
#' @returns One of `"NA"`, `"NaN"` or `"an infinite value"`.
#'
#' @noRd
#' @keywords internal
.aci_bad_kind <- function(value) {
  if (is.nan(value)) {
    return("NaN")
  }
  if (is.na(value)) {
    return("NA")
  }
  "an infinite value"
}

# -- scalar and signal checks -------------------------------------------------

#' Require a single finite numeric value
#'
#' @param value The value to check.
#' @param name The argument name, used in the error message.
#'
#' @returns `value`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(
      sprintf("`%s` must be a single finite numeric value.", name),
      call. = FALSE
    )
  }
  invisible(value)
}

#' Require a positive scalar
#'
#' @param value The value to check.
#' @param name The argument name, used in the error message.
#'
#' @returns `value`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_positive <- function(value, name) {
  .aci_check_scalar(value, name)
  if (value <= 0) {
    stop(sprintf("`%s` must be positive.", name), call. = FALSE)
  }
  invisible(value)
}

#' Validate an observed signal
#'
#' The numerical core assumes a plain numeric vector of at least two complete,
#' finite observations on a regular time grid. This check enforces everything
#' but the regularity of the grid, which only [aci()] can see.
#'
#' @param x The observed signal.
#' @param name The argument name, used in the error message.
#'
#' @returns `x`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_signal <- function(x, name = "x") {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(
      sprintf(
        "`%s` must be a plain numeric vector; it is %s.",
        name, .aci_describe(x)
      ),
      call. = FALSE
    )
  }
  if (length(x) < 2L) {
    stop(
      sprintf("`%s` must be a numeric vector of at least two observations.",
        name
      ),
      call. = FALSE
    )
  }
  bad <- which(!is.finite(x))
  if (length(bad) > 0L) {
    stop(
      sprintf(
        paste0(
          "`%s` must be complete and finite; it holds %s at index %d. ",
          "The filter has no missing-observation contract: interpolate or ",
          "subset to a complete, regularly sampled span before calling."
        ),
        name, .aci_bad_kind(x[bad[1L]]), bad[1L]
      ),
      call. = FALSE
    )
  }
  invisible(x)
}

#' Validate a time vector and derive its step
#'
#' The numerical core integrates a fixed step, so an observed time vector is
#' only admissible when it is strictly increasing and equally spaced. An
#' irregular grid is rejected rather than approximated: the closed-form
#' recursions have no contract for one.
#'
#' @param time The time vector.
#' @param n The number of observations it must cover.
#' @param name The argument name, used in the error message.
#'
#' @returns The derived time step, a positive numeric scalar.
#'
#' @noRd
#' @keywords internal
.aci_check_time <- function(time, n, name = "time") {
  if (!is.numeric(time) || !is.null(dim(time))) {
    stop(
      sprintf("`%s` must be a plain numeric vector; it is %s.",
        name, .aci_describe(time)
      ),
      call. = FALSE
    )
  }
  if (length(time) != n) {
    stop(
      sprintf(
        "`%s` must have one value per observation: expected length %d, got %d.",
        name, n, length(time)
      ),
      call. = FALSE
    )
  }
  bad <- which(!is.finite(time))
  if (length(bad) > 0L) {
    stop(
      sprintf("`%s` must be finite; it holds %s at index %d.",
        name, .aci_bad_kind(time[bad[1L]]), bad[1L]
      ),
      call. = FALSE
    )
  }
  steps <- diff(time)
  bad <- which(steps <= 0)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "`%s` must be strictly increasing; it is not at index %d.",
        name, bad[1L] + 1L
      ),
      call. = FALSE
    )
  }
  dt <- mean(steps)
  if (max(abs(steps - dt)) > .aci_tol(dt)) {
    stop(
      sprintf(
        paste0(
          "`%s` must be equally spaced: the recursions integrate a fixed ",
          "step, and the spacing ranges from %g to %g. Resample onto a ",
          "regular grid before calling."
        ),
        name, min(steps), max(steps)
      ),
      call. = FALSE
    )
  }
  dt
}

#' Validate a model label
#'
#' @param label The label.
#'
#' @returns `label`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_label <- function(label) {
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop(
      "`label` must be a single non-missing character string.",
      call. = FALSE
    )
  }
  invisible(label)
}

#' Validate a named list of model parameters
#'
#' @param parameters `NULL`, or a named list of finite numeric scalars.
#'
#' @returns `parameters`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_parameters <- function(parameters) {
  if (is.null(parameters)) {
    return(invisible(NULL))
  }
  scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }
  named <- !is.null(names(parameters)) && all(nzchar(names(parameters)))
  all_scalar <- all(vapply(parameters, scalar, logical(1L)))
  if (!is.list(parameters) || !named || !all_scalar) {
    stop(
      paste0(
        "`parameters` must be NULL or a named list of finite numeric ",
        "scalars; it records the values that define the model."
      ),
      call. = FALSE
    )
  }
  invisible(parameters)
}

# -- coefficient evaluation ---------------------------------------------------

#' Evaluate a model coefficient at the observed signal and validate the result
#'
#' Coefficient functions are supplied by the user, so their return value is
#' part of the package's input contract rather than of its internals. A
#' function that returns a scalar is rejected rather than recycled: silent
#' recycling would leave the filter reading a constant where the model
#' promised a per-step coefficient, and a genuinely constant coefficient is
#' already expressible as a numeric scalar at construction.
#'
#' @param fn The coefficient function.
#' @param name The coefficient name, used in the error message.
#' @param x The observed signal.
#'
#' @returns The evaluated coefficient, a finite numeric vector of
#'   `length(x)`.
#'
#' @noRd
#' @keywords internal
.aci_eval_coef <- function(fn, name, x) {
  value <- fn(x)
  n <- length(x)
  if (!is.numeric(value) || !is.null(dim(value))) {
    stop(
      sprintf(
        "`%s` must return a numeric vector of length %d; it returned %s.",
        name, n, .aci_describe(value)
      ),
      call. = FALSE
    )
  }
  if (length(value) != n) {
    stop(
      sprintf(
        paste0(
          "`%s` must return one value per observation: expected length %d, ",
          "got %d. Coefficient functions are not recycled; supply a ",
          "vectorised function of the observed signal, or a numeric scalar ",
          "for a coefficient that is constant in time."
        ),
        name, n, length(value)
      ),
      call. = FALSE
    )
  }
  bad <- which(!is.finite(value))
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "`%s` must return finite values; it returned %s at index %d.",
        name, .aci_bad_kind(value[bad[1L]]), bad[1L]
      ),
      call. = FALSE
    )
  }
  value
}

#' Coerce a coefficient argument to a vectorised function
#'
#' A numeric scalar is the documented shorthand for a coefficient that is
#' constant in time, and is widened here to the function form the core
#' consumes.
#'
#' @param value A function or a finite numeric scalar.
#' @param name The argument name, used in the error message.
#'
#' @returns A function of the observed signal.
#'
#' @noRd
#' @keywords internal
.aci_as_coef <- function(value, name) {
  if (is.function(value)) {
    return(value)
  }
  if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
    force(value)
    return(function(x) rep_len(value, length(x)))
  }
  stop(
    sprintf(
      paste0(
        "`%s` must be a vectorised function of the observed signal or a ",
        "finite numeric scalar; it is %s."
      ),
      name, .aci_describe(value)
    ),
    call. = FALSE
  )
}

# -- noise covariance ---------------------------------------------------------

#' Validate a scalar noise covariance system
#'
#' Enforces the admissibility of the joint noise covariance of the observed
#' and unobserved processes. The observation-noise covariance must be strictly
#' positive because the filter inverts it; the latent-noise covariance must be
#' non-negative; and the joint covariance must be positive semidefinite, which
#' for the scalar system reduces to a non-negative determinant.
#'
#' A singular system, whose determinant is zero, is admissible: it describes
#' perfectly correlated noise sources. It is not rejected here because the
#' per-step covariance checks in the filter and the smoother are the right
#' place to catch the run it may destabilise.
#'
#' @param S_xoS_x The observation-noise covariance.
#' @param S_yoS_y The latent-noise covariance.
#' @param S_yoS_x The latent-to-observation noise cross-covariance.
#'
#' @returns `TRUE`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_noise_covariance <- function(S_xoS_x, S_yoS_y, S_yoS_x) {
  .aci_check_scalar(S_xoS_x, "S_xoS_x")
  .aci_check_scalar(S_yoS_y, "S_yoS_y")
  .aci_check_scalar(S_yoS_x, "S_yoS_x")
  if (S_xoS_x <= 0) {
    stop(
      "`S_xoS_x` must be positive; it is the observation-noise covariance.",
      call. = FALSE
    )
  }
  if (S_yoS_y < 0) {
    stop(
      "`S_yoS_y` must be non-negative; it is the latent-noise covariance.",
      call. = FALSE
    )
  }
  determinant <- S_xoS_x * S_yoS_y - S_yoS_x^2
  if (determinant < -.aci_tol(S_xoS_x, S_yoS_y, S_yoS_x)) {
    stop(
      sprintf(
        paste0(
          "The joint noise covariance is not positive semidefinite: ",
          "`S_xoS_x` * `S_yoS_y` - `S_yoS_x`^2 is %g. The noise ",
          "cross-covariance cannot exceed sqrt(`S_xoS_x` * `S_yoS_y`) = %g ",
          "in magnitude."
        ),
        determinant, sqrt(S_xoS_x * S_yoS_y)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# -- components and posteriors ------------------------------------------------

#' The names of a conditional Gaussian components list
#'
#' @noRd
#' @keywords internal
.aci_component_names <- c(
  "L_x", "f_x", "L_y", "f_y",
  "S_xoS_x", "S_yoS_y", "S_yoS_x", "S_xoS_y"
)

#' Validate a conditional Gaussian components list
#'
#' The components list is the package's expert-level extension surface: a user
#' may hand-build one for a system aciR provides no constructor for. It is
#' therefore validated as rigorously as a model object.
#'
#' @param comp The components list.
#' @param n The number of time steps the per-step coefficients must cover.
#' @param name The argument name, used in the error message.
#'
#' @returns `TRUE`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_components <- function(comp, n, name = "comp") {
  if (!is.list(comp)) {
    stop(
      sprintf(
        paste0(
          "`%s` must be a conditional Gaussian components list; it is %s. ",
          "See `?aci_components` for the expected shape."
        ),
        name, .aci_describe(comp)
      ),
      call. = FALSE
    )
  }
  absent <- setdiff(.aci_component_names, names(comp))
  if (length(absent) > 0L) {
    stop(
      sprintf(
        paste0(
          "`%s` is missing the component%s %s. See `?aci_components` for ",
          "the expected shape."
        ),
        name,
        if (length(absent) > 1L) "s" else "",
        paste(sprintf("`%s`", absent), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Per-step coefficients carry one value per observation; the self-drift of
  # the unobserved component is constant in time for the systems this core
  # integrates.
  for (entry in c("L_x", "f_x", "f_y")) {
    value <- comp[[entry]]
    if (!is.numeric(value) || !is.null(dim(value)) || length(value) != n) {
      stop(
        sprintf(
          paste0(
            "`%s$%s` must be a numeric vector of length %d, one value per ",
            "observation; it is %s."
          ),
          name, entry, n, .aci_describe(value)
        ),
        call. = FALSE
      )
    }
    bad <- which(!is.finite(value))
    if (length(bad) > 0L) {
      stop(
        sprintf(
          "`%s$%s` must be finite; it holds %s at index %d.",
          name, entry, .aci_bad_kind(value[bad[1L]]), bad[1L]
        ),
        call. = FALSE
      )
    }
  }
  .aci_check_scalar(comp$L_y, sprintf("%s$L_y", name))
  .aci_check_noise_covariance(comp$S_xoS_x, comp$S_yoS_y, comp$S_yoS_x)
  .aci_check_scalar(comp$S_xoS_y, sprintf("%s$S_xoS_y", name))
  asymmetry <- abs(comp$S_xoS_y - comp$S_yoS_x)
  if (asymmetry > .aci_tol(comp$S_yoS_x, comp$S_xoS_y)) {
    stop(
      sprintf(
        paste0(
          "`%s$S_xoS_y` must equal `%s$S_yoS_x`: the noise cross-covariance ",
          "of a scalar system is symmetric, but they are %g and %g."
        ),
        name, name, comp$S_xoS_y, comp$S_yoS_x
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate a posterior structure
#'
#' @param post The posterior, a list of `mean` and `cov`.
#' @param name The argument name, used in the error message.
#' @param n Optional. The number of time steps the posterior must cover.
#'
#' @returns `TRUE`, invisibly.
#'
#' @noRd
#' @keywords internal
.aci_check_posterior <- function(post, name, n = NULL) {
  if (!is.list(post) || !all(c("mean", "cov") %in% names(post))) {
    stop(
      sprintf(
        paste0(
          "`%s` must be a list with numeric `mean` and `cov`, as returned by ",
          "`aci_filter()`; it is %s."
        ),
        name, .aci_describe(post)
      ),
      call. = FALSE
    )
  }
  paired <- is.numeric(post$mean) && is.numeric(post$cov) &&
    length(post$mean) == length(post$cov)
  if (!paired) {
    stop(
      sprintf(
        paste0(
          "`%s$mean` and `%s$cov` must be numeric vectors of equal length; ",
          "they are %s and %s."
        ),
        name, name, .aci_describe(post$mean), .aci_describe(post$cov)
      ),
      call. = FALSE
    )
  }
  if (!is.null(n) && length(post$mean) != n) {
    stop(
      sprintf(
        "`%s` must cover %d time steps; it covers %d.",
        name, n, length(post$mean)
      ),
      call. = FALSE
    )
  }
  bad <- which(!is.finite(post$mean))
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "`%s$mean` must be finite; it holds %s at index %d.",
        name, .aci_bad_kind(post$mean[bad[1L]]), bad[1L]
      ),
      call. = FALSE
    )
  }
  bad <- which(!is.finite(post$cov) | post$cov <= 0)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        paste0(
          "`%s$cov` must be finite and strictly positive; it holds %s at ",
          "index %d. A non-positive posterior covariance has no ",
          "relative-entropy interpretation."
        ),
        name, format(post$cov[bad[1L]]), bad[1L]
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# -- runtime diagnostics ------------------------------------------------------

#' Report a covariance that left the admissible domain during a recursion
#'
#' Called only from the failure branch of the filter and smoother loops, so
#' the per-step cost of the guard is one comparison.
#'
#' @param algorithm Either `"filter"` or `"smoother"`.
#' @param index The step index at which the covariance failed.
#' @param time The time corresponding to `index`.
#' @param value The offending covariance.
#'
#' @returns Never returns; raises an error.
#'
#' @noRd
#' @keywords internal
.aci_stop_covariance <- function(algorithm, index, time, value) {
  stop(
    sprintf(
      paste0(
        "The %s covariance must stay finite and strictly positive; it became ",
        "%s at index %d (time %g). Reduce `dt`, or check the model's noise ",
        "covariance: an explicit Euler step too large for the system can ",
        "drive the covariance out of its domain even when the model is ",
        "admissible."
      ),
      algorithm, format(value), index, time
    ),
    call. = FALSE
  )
}
