# Model and observation contracts

Constructors and coercion helpers for uniformly sampled observations,
general stochastic models, and conditional-Gaussian nonlinear systems.
[`aci_model_from_affine()`](https://biometryhub.github.io/ACI/reference/aci_model_from_affine.md)
verifies that the supplied drifts are affine in the hidden state before
constructing a CGNS representation.

## See also

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html)
