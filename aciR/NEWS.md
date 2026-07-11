# aciR 0.0.0.9000

* Initial development scaffold.
* Adds the conditional Gaussian nonlinear system (CGNS) core: a forward
  filter (`aci_filter()`), a backward smoother (`aci_smoother()`) and the
  relative-entropy causal-information metric (`aci_metric()`), each taking a
  general CGNS components list so they are not tied to a single model.
* Adds `aci_dyad_components()` and `aci_simulate_dyad()` for the nonlinear
  dyad model with intermittent extreme events used as the worked example.
* Adds a model layer over the core. `aci_cgns_model()` is the general
  constructor for a CGNS `aci_model` object, `aci_dyad_model()` builds the
  flagship nonlinear dyad model, `aci_simulate()` draws an Euler-Maruyama
  realisation of a model, and `aci()` is the single user-facing entry point
  that runs the filter, smoother and metric and returns an `aci` object.
  `print()` methods are provided for both classes.

## Roadmap

* A noisy predator-prey (Lotka-Volterra) model constructor. Its unobserved
  process has a time-varying self-drift, whereas the current numerical core
  integrates a time-invariant self-drift; the constructor waits on a core that
  admits a time-varying self-drift and on its own independent-oracle fixture.
* The causal-influence-range (`aci_cir()`): the subjective and objective
  influence-range lengths from the fixed-lag online smoother. This is an
  order-N-squared computation and waits on its own independent-oracle fixture.
