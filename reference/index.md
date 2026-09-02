# Package index

## Assimilation and the causal measure

- [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) :
  Assimilative causal inference
- [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
  : Data assimilation filter
- [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
  : Data assimilation smoother
- [`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
  : Gaussian relative entropy along a pair of paths
- [`aci_metric_pair()`](https://biometryhub.github.io/ACI/reference/aci_metric_pair.md)
  : Gaussian relative entropy
- [`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
  : Fixed-lag online data assimilation

## Influence range and conditioning

- [`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
  : Causal influence range
- [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
  : Conditional ACI specification
- [`aci_conditional_reduce()`](https://biometryhub.github.io/ACI/reference/aci_conditional_reduce.md)
  : Reduce a model by prescribing the conditioning channels
- [`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
  : Finite-lag divergence table
- [`lt_diag()`](https://biometryhub.github.io/ACI/reference/lt_diag.md)
  : Diagonal of a lag table
- [`lt_row()`](https://biometryhub.github.io/ACI/reference/lt_row.md) :
  One row of a lag table
- [`lt_tail_bound()`](https://biometryhub.github.io/ACI/reference/lt_tail_bound.md)
  : Heuristic tail estimate of a lag table

## Models and simulation

- [`aci_model()`](https://biometryhub.github.io/ACI/reference/aci_model.md)
  : Conditional-Gaussian nonlinear system
- [`aci_model_from_affine()`](https://biometryhub.github.io/ACI/reference/aci_model_from_affine.md)
  : CGNS model from an affine-in-hidden drift pair
- [`aci_linear_model()`](https://biometryhub.github.io/ACI/reference/aci_linear_model.md)
  : Conditionally linear model from constant coefficients
- [`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)
  : Nonlinear dyad (andreou2026aci eq. 1-2)
- [`aci_predprey_model()`](https://biometryhub.github.io/ACI/reference/aci_predprey_model.md)
  : Noisy predator-prey benchmark model
- [`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md)
  : Six-variable stochastic conceptual ENSO model (andreou2026aci SI via
  chen2022enso eqs. 1a-1f).
- [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  : Simulate a path from a model
- [`as_obs()`](https://biometryhub.github.io/ACI/reference/as_obs.md) :
  Coerce to an observed trajectory
- [`observed_trajectory()`](https://biometryhub.github.io/ACI/reference/observed_trajectory.md)
  : Observed trajectory on a uniform time grid

## Numerical policies

- [`safe_chol()`](https://biometryhub.github.io/ACI/reference/safe_chol.md)
  : Cholesky with escalating jitter ladder
- [`spd_floor()`](https://biometryhub.github.io/ACI/reference/spd_floor.md)
  : Project to SPD by eigenvalue flooring (used after Euler covariance
  updates)

## Methods

Print, plot, coercion and simulation methods for the package’s classes.

- [`print(`*`<aci_result>`*`)`](https://biometryhub.github.io/ACI/reference/print.aci_result.md)
  : Print aci() result
- [`print(`*`<cir_result>`*`)`](https://biometryhub.github.io/ACI/reference/print.cir_result.md)
  : Print the causal influence range
- [`print(`*`<da_path_gaussian>`*`)`](https://biometryhub.github.io/ACI/reference/print.da_path_gaussian.md)
  : Print a Gaussian assimilation path
- [`print(`*`<lag_table>`*`)`](https://biometryhub.github.io/ACI/reference/print.lag_table.md)
  : Print a lag table
- [`print(`*`<obs_traj>`*`)`](https://biometryhub.github.io/ACI/reference/print.obs_traj.md)
  : Print an observed trajectory
- [`print(`*`<aci_conditional_spec>`*`)`](https://biometryhub.github.io/ACI/reference/print.aci_conditional_spec.md)
  : Print a conditional ACI specification
- [`print(`*`<stochastic_model>`*`)`](https://biometryhub.github.io/ACI/reference/print.stochastic_model.md)
  : Print a stochastic model
- [`plot(`*`<aci_result>`*`)`](https://biometryhub.github.io/ACI/reference/plot.aci_result.md)
  : Plot an ACI result
- [`plot(`*`<cir_result>`*`)`](https://biometryhub.github.io/ACI/reference/plot.cir_result.md)
  : Plot a causal influence range result
- [`plot(`*`<da_path_gaussian>`*`)`](https://biometryhub.github.io/ACI/reference/plot.da_path_gaussian.md)
  : Plot a Gaussian assimilation path
- [`as.data.frame(`*`<aci_result>`*`)`](https://biometryhub.github.io/ACI/reference/as.data.frame.aci_result.md)
  : aci() result to data.frame
- [`as.data.frame(`*`<da_path_gaussian>`*`)`](https://biometryhub.github.io/ACI/reference/as.data.frame.da_path_gaussian.md)
  : Coerce a Gaussian assimilation path to a data frame
- [`as.data.frame(`*`<lag_table>`*`)`](https://biometryhub.github.io/ACI/reference/as.data.frame.lag_table.md)
  : Coerce a lag table to a data frame
- [`as.data.frame(`*`<obs_traj>`*`)`](https://biometryhub.github.io/ACI/reference/as.data.frame.obs_traj.md)
  : Coerce an observed trajectory to a data frame
- [`simulate(`*`<stochastic_model>`*`)`](https://biometryhub.github.io/ACI/reference/simulate.stochastic_model.md)
  : Simulate a stochastic or conditional-Gaussian model

## Topics

The model contract, the assimilation surface, the causal measures, and
the references.

- [`model_contracts`](https://biometryhub.github.io/ACI/reference/model_contracts.md)
  : Model and observation contracts
- [`assimilation_api`](https://biometryhub.github.io/ACI/reference/assimilation_api.md)
  : Data assimilation and finite-lag API
- [`causal_metrics`](https://biometryhub.github.io/ACI/reference/causal_metrics.md)
  : Assimilative causal metrics and influence ranges
- [`aci_references`](https://biometryhub.github.io/ACI/reference/aci_references.md)
  : References
