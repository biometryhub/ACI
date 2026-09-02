# Oracle models shared by more than one test file.
#
# The mv model backs the independent-MATLAB matrix oracle in
# test-19-compiled-oracles.R and the matrix arm of the public online smoother
# in test-27-online-lag.R. It lives here rather than in either file so the two
# grade the same coefficients rather than two transcriptions of them.

.compiled_oracle_mv_model <- function() {
  noise <- matrix(
    c(
      0.60, 0.10, 0.25, 0.05,
      0.20, 0.50, 0.10, 0.30,
      0.15, 0.05, 0.70, 0.10,
      0.05, 0.10, 0.15, 0.55
    ),
    4L, 4L
  )
  sx <- noise[1:2, , drop = FALSE]
  sy <- noise[3:4, , drop = FALSE]
  aci_model(
    Lx = function(t, x) matrix(
      c(0.8 + 0.3 * x[1L], 0.1 * x[2L],
        0.2 * sin(t), 0.6 - 0.2 * x[1L]),
      2L, 2L
    ),
    fx = function(t, x) c(0.4 - 0.5 * x[1L], -0.3 * x[2L] + 0.2),
    Ly = function(t, x) matrix(
      c(-1.2 + 0.1 * x[1L], 0.2,
        0.3, -0.9 - 0.1 * x[2L]),
      2L, 2L
    ),
    fy = function(t, x) c(0.5 - 0.2 * x[1L]^2,
                           0.1 - 0.15 * x[2L]),
    Sx1 = function(t, x) sx[, 1:2, drop = FALSE],
    Sx2 = function(t, x) sx[, 3:4, drop = FALSE],
    Sy1 = function(t, x) sy[, 1:2, drop = FALSE],
    Sy2 = function(t, x) sy[, 3:4, drop = FALSE],
    k = 2L, l = 2L, name = "independent-matrix-oracle"
  )
}
