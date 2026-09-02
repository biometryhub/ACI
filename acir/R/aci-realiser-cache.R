################################################################################
## aci-realiser-cache.R - one realisation of a closure model per grid
################################################################################

# A model supplied as coefficient closures is realised on the observed grid by
# .realise_cgns_grid_once(), one closure call per grid point per coefficient.
# The library models avoid that through a sealed batch realiser; the generic
# route paid it on every verb, so filter, smoother, online smoother and range
# on one model and one record realised the same arrays four times. This file
# keeps the last few realisations and hands one back when the same model meets
# the same grid again.
#
# The key is the model object and the grid. identical() on the model compares
# the closures' formals and bodies and requires the same enclosing
# environments, so two models built from the same expressions in different
# calls do not collide. What identical() cannot see is a value the closures
# read from their environment that moved after the realisation was stored,
# so a hit is confirmed by re-evaluating the coefficients at three grid
# points, the first, the middle and the last, and comparing bit for bit with
# the stored arrays there. Storing an entry evaluates no closure, so the
# generic route keeps its contract of one evaluation per grid point. A
# captured value that changed everywhere fails the probe; one that changed
# only away from all three points would not, and the option below switches
# the cache off for a model of that kind.

.realise_cache <- new.env(parent = emptyenv())
.realise_cache$entries <- list()
.realise_cache_size <- 4L


#' The probe indices of a grid (internal)
#'
#' @param n1 Number of grid points.
#' @returns The first, middle and last index, without repeats.
#' @noRd
.realise_probe_index <- function(n1) unique(c(1L, (n1 + 1L) %/% 2L, n1))


#' Evaluate the coefficients at the probe points of a grid (internal)
#'
#' @param model A `cgns_model`.
#' @param obs An `obs_traj`.
#' @returns One plain numeric vector per probe point, the coefficients in the
#'   order `Lx`, `fx`, `Ly`, `fy`, `gxx`, `gyy`, `gyx`, column-major.
#' @noRd
.realise_probe_eval <- function(model, obs) {
  lapply(.realise_probe_index(length(obs$t)), function(j) {
    co <- eval_coefs(model, obs$t[j], obs$x[j, ])
    c(as.vector(co$Lx), as.vector(co$fx), as.vector(co$Ly), as.vector(co$fy),
      as.vector(co$gxx), as.vector(co$gyy), as.vector(co$gyx))
  })
}


#' Read the probe points back from realised arrays (internal)
#'
#' The same vectors as `.realise_probe_eval()` would produce, taken from the
#' arrays a realisation stored, so that storing an entry evaluates no closure.
#'
#' @param full Realised coefficient arrays.
#' @param n1 Number of grid points.
#' @returns One plain numeric vector per probe point.
#' @noRd
.realise_probe_stored <- function(full, n1) {
  lapply(.realise_probe_index(n1), function(j) {
    c(as.vector(full$Lx[, , j]), as.vector(full$fx[j, ]),
      as.vector(full$Ly[, , j]), as.vector(full$fy[j, ]),
      as.vector(full$gxx[, , j]), as.vector(full$gyy[, , j]),
      as.vector(full$gyx[, , j]))
  })
}


#' Empty the realisation cache (internal)
#'
#' @returns `NULL`, invisibly.
#' @noRd
.realise_cache_clear <- function() {
  .realise_cache$entries <- list()
  invisible(NULL)
}


#' Realise a closure model on a grid, reusing a stored realisation (internal)
#'
#' Returns the arrays `.realise_cgns_grid_once()` would return for this model
#' and grid. When the same model has met the same grid before and its
#' coefficients still evaluate to the stored values at the probe points, the
#' stored arrays are returned without re-evaluating the closures. Otherwise
#' the grid is realised afresh and the result stored, displacing the least
#' recently used entry once `.realise_cache_size` are held. The option
#' `aci.realiser_cache` set to `FALSE` bypasses the cache entirely.
#'
#' @param model A `cgns_model`.
#' @param obs An `obs_traj` with the model's observed dimension.
#' @returns Realised coefficient arrays, as from `.realise_cgns_grid_once()`,
#'   carrying the unconditioned precision path as the attribute `gxx_weight`
#'   whenever the cache is in use.
#' @noRd
.realise_cgns_grid_cached <- function(model, obs) {
  if (!isTRUE(getOption("aci.realiser_cache", TRUE)))
    return(.realise_cgns_grid_once(model, obs))
  entries <- .realise_cache$entries
  for (i in seq_along(entries)) {
    e <- entries[[i]]
    same <- identical(e$model, model) && identical(e$t, obs$t) &&
      identical(e$x, obs$x)
    if (!same) next
    if (identical(.realise_probe_eval(model, obs), e$probes)) {
      .realise_cache$entries <- c(list(e), entries[-i])
      return(e$full)
    }
    entries <- entries[-i]
    break
  }
  full <- .realise_cgns_grid_once(model, obs)
  ## The unconditioned precision path is a pure function of these arrays, so
  ## it is realised once with them; the conditioned routes derive their own.
  attr(full, "gxx_weight") <- .compiled_precision_path(
    full$gxx, length(obs$t) - 1L
  )
  entry <- list(model = model, t = obs$t, x = obs$x,
                probes = .realise_probe_stored(full, length(obs$t)),
                full = full)
  kept <- c(list(entry), entries)
  .realise_cache$entries <- utils::head(kept, .realise_cache_size)
  full
}
