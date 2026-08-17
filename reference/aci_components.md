# Conditional Gaussian components list

Several functions in aciR share a `comp` argument, a named list holding
the per-step coefficients of a conditional Gaussian nonlinear system
(CGNS). For an observed process `x` and an unobserved (conditionally
Gaussian) process `y`, the coefficients are those of the pair of
stochastic differential equations governing `x` and `y`. The list has
the following entries.

## Details

- `L_x`:

  Numeric vector, length `n`. The coupling of the unobserved component
  into the drift of the observed process at each time step.

- `f_x`:

  Numeric vector, length `n`. The remaining drift of the observed
  process at each time step.

- `L_y`:

  Numeric scalar, or numeric vector of length `n`. The linear self-drift
  of the unobserved component. A scalar is a self-drift constant in
  time, as in the supplied dyad model; a vector carries one value per
  observation, as the supplied predator-prey model needs. It may depend
  on the observed signal but never on the unobserved component, which is
  what keeps the system conditionally Gaussian.

- `f_y`:

  Numeric vector, length `n`. The remaining drift of the unobserved
  process at each time step.

- `S_xoS_x`:

  Numeric scalar. The observation-noise covariance, the product of the
  observed-process noise coefficient with its transpose. Must be
  strictly positive.

- `S_yoS_y`:

  Numeric scalar. The latent-noise covariance of the unobserved process.
  Must be non-negative.

- `S_yoS_x`:

  Numeric scalar. The latent-to-observation noise cross-covariance.

- `S_xoS_y`:

  Numeric scalar. The observation-to-latent noise cross-covariance, the
  transpose of `S_yoS_x`. For the scalar systems this package integrates
  the two are equal, and a components list in which they disagree is
  rejected.

The joint noise covariance must be positive semidefinite, which for a
scalar system means `S_xoS_x * S_yoS_y - S_yoS_x^2` is non-negative. A
components list that violates any of these conditions is rejected before
the recursion starts.

This schema is the package's extension surface for advanced use. Build a
components list directly to run the core on a conditional Gaussian
system for which aciR supplies no constructor. See
[`aci_dyad_components()`](https://biometryhub.github.io/ACI/reference/aci_dyad_components.md)
for a worked example, and
[`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md)
for the higher-level alternative.
