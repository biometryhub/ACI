# Shared fixture for the public provenance grid and stored/streamed CIR checks.
.compiled_cir_matrix_setup <- function() {
  model <- aci_model(
    Lx = function(t, x) matrix(c(
      0.65 + 0.04 * x[1], -0.12,
      0.08, 0.45 + 0.03 * x[2]
    ), 2, 2, byrow = TRUE),
    fx = function(t, x) c(-0.25 * x[1] + 0.08,
                           -0.18 * x[2] - 0.03),
    Ly = function(t, x) matrix(c(-0.7, -0.04, 0.09, -0.55),
                                2, 2, byrow = TRUE),
    fy = function(t, x) c(0.1 * sin(t) - 0.06 * x[1],
                           0.12 * cos(t) + 0.04 * x[2]),
    Sx1 = function(t, x) matrix(c(0.65, 0.05, 0, 0.55), 2, 2),
    Sy1 = function(t, x) matrix(c(0.05, 0.01, 0, 0.04), 2, 2),
    Sy2 = function(t, x) matrix(c(0.72, 0.03, 0, 0.64), 2, 2),
    k = 2, l = 2, name = "compiled-cir-matrix"
  )
  t <- seq(0, 0.12, by = 0.01)
  obs <- observed_trajectory(
    t, cbind(0.25 + 0.08 * sin(4 * t), -0.15 + 0.06 * cos(3 * t))
  )
  list(
    model = model, obs = obs,
    init = list(
      mean = c(0.15, -0.08),
      cov = matrix(c(0.42, 0.025, 0.025, 0.31), 2, 2)
    )
  )
}

# Coefficients use only base functions. Avoid serializing the test runner's
# environment when the public retention checks save and reload this model.
environment(.compiled_cir_matrix_setup) <- asNamespace("acir")
