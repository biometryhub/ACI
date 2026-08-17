# the vector dataset contract --------------------------------------------------
#
# The scalar bundle carries one column per per-step quantity. The vector bundle
# has to carry arrays of matrices, and how those reach CSV is a real choice.
# Three options, all viable:
#
#   (a) One row per time step, the slice flattened COLUMN-MAJOR (chosen).
#       MATLAB's A(:,:,j)(:) and R's as.vector(A[,,j]) are both column-major, so
#       neither side transposes and a transposition bug has nowhere to hide.
#   (b) Long format: one row per (slice, i, j, value). Self-describing, but
#       three times the bytes and both sides must reassemble, which is more
#       code on precisely the axis where the mistakes happen.
#   (c) One file per matrix entry, A11.csv .. A33.csv. Trivially readable, but
#       nine files per component and the naming becomes the contract.
#
# Noise is carried as the feedback matrices S_x and S_y rather than as the
# Grammians, matching the scalar bundle, and deliberately, because forming the
# Grammians and inverting them is itself a place the two implementations part
# company: the reference pseudo-inverts, aciR factorises.

#' Read a vector dataset bundle.
#'
#' @param path Directory holding `meta.dcf` and the component CSVs.
#'
#' @returns A list with `meta`, the observed signal `x`, and the components.
read_vector_dataset <- function(path) {
  meta <- as.list(read_manifest(file.path(path, "meta.dcf"))[1L, ])
  for (field in c("N", "dt", "k", "l", "R0")) {
    meta[[field]] <- as.numeric(meta[[field]])
  }
  n <- meta$N + 1L
  k <- as.integer(meta$k)
  l <- as.integer(meta$l)

  flat <- function(name) {
    as.matrix(utils::read.csv(file.path(path, paste0(name, ".csv")),
                              header = FALSE))
  }
  # Row j holds the slice column-major, so t() then array() rebuilds it without
  # any transposition on either side.
  cube <- function(name, r, c) array(t(flat(name)), c(r, c, n))

  list(
    meta = meta,
    x = t(flat("x")),
    mu0 = as.numeric(flat("mu0")),
    R0 = meta$R0 * diag(l),
    S_x = cube("S_x", k, k),
    S_y = cube("S_y", l, l),
    L_x = cube("L_x", k, l),
    L_y = cube("L_y", l, l),
    f_x = t(flat("f_x")),
    f_y = t(flat("f_y"))
  )
}

#' Build the aciR vector components from a bundle.
#'
#' @param d A list from [read_vector_dataset()].
#'
#' @returns A conditional Gaussian components list for the vector path.
vector_components <- function(d) {
  n <- dim(d$L_x)[3L]
  gram <- function(A, B) {
    out <- array(0, c(dim(A)[1L], dim(B)[1L], n))
    for (j in seq_len(n)) out[, , j] <- A[, , j] %*% t(B[, , j])
    out
  }
  # Whether the observed and unobserved processes are driven by the SAME Wiener
  # increments or by disjoint ones decides the cross-Grammian entirely, and it
  # cannot be read off the two blocks: both are 3x3 either way. Declared in the
  # bundle, and absent means refuse. Reading the ENSO blocks as shared gives a
  # cross-Grammian of 0.38 and drives the filter covariance to NA by step 1872,
  # where that model's own script states the cross terms are absent.
  coupling <- d$meta$NoiseCoupling
  if (is.null(coupling) || is.na(coupling)) {
    stop(
      paste0(
        "bundle does not declare `NoiseCoupling`. It must be `disjoint` (the ",
        "observed and unobserved processes are driven by separate Wiener ",
        "increments, so the noise cross-Grammian is exactly zero) or `shared` ",
        "(the same increments, so it is S_y S_x'). This cannot be inferred ",
        "from the blocks themselves."
      ),
      call. = FALSE
    )
  }
  s_yox <- switch(
    coupling,
    disjoint = array(0, c(dim(d$S_y)[1L], dim(d$S_x)[1L], n)),
    shared = gram(d$S_y, d$S_x),
    stop(sprintf("unknown NoiseCoupling `%s`.", coupling), call. = FALSE)
  )

  list(
    L_x = d$L_x, L_y = d$L_y, f_x = d$f_x, f_y = d$f_y,
    S_xoS_x = gram(d$S_x, d$S_x),
    S_yoS_y = gram(d$S_y, d$S_y),
    S_yoS_x = s_yox
  )
}
