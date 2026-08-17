# Package index

## Model layer

High-level entry points. Build a conditional Gaussian nonlinear system
as a model object, simulate a realisation of it, and run the whole
assimilative causal inference workflow in a single call.

- [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) : Run
  assimilative causal inference
- [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md)
  : Conditional Gaussian nonlinear system model
- [`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)
  : Nonlinear dyad model
- [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  : Simulate a conditional Gaussian nonlinear system

## Conditional Gaussian core

The numerical core. The forward filter, the backward smoother and the
causal-information metric, each written against a general components
list so they apply to any conditional Gaussian nonlinear system.

- [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
  : Forward conditional Gaussian filter
- [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
  : Backward conditional Gaussian smoother
- [`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
  : Assimilative causal-information metric

## Causal influence range

How far ahead in the observed record the estimate of the unobserved
state keeps improving, and the fixed-lag online smoother it is built
from. These answer a different question from the causal-information
metric, which measures how much the record says rather than how long it
takes to say it.

- [`aci_online_smoother()`](https://biometryhub.github.io/ACI/reference/aci_online_smoother.md)
  : Fixed-lag online conditional Gaussian smoother
- [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  : Causal influence range

## Vector states and conditional questions

The numerical core takes vector-valued states with matrix coefficients,
and the conditional construction asks what one observed component says
about the unobserved state given that the others are also watched.

- [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
  : Condition the causal question on a subset of the observed components

## The stochastic ENSO model

The six-dimensional El Nino Southern Oscillation system of the paper’s
case study. The largest system the package expresses, and the only one
whose noise covariances vary in time.

- [`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md)
  : Stochastic ENSO model
- [`aci_enso_components()`](https://biometryhub.github.io/ACI/reference/aci_enso_components.md)
  : Conditional Gaussian components of the stochastic ENSO model
- [`aci_enso_parameters()`](https://biometryhub.github.io/ACI/reference/aci_enso_parameters.md)
  : Parameters of the stochastic ENSO model

## The noisy predator-prey model

A stochastic Lotka-Volterra pair, studied in either causal direction. It
is the model whose latent self-drift is set by the observed state, and
so the one that exercises a self-drift varying in time.

- [`aci_predprey_model()`](https://biometryhub.github.io/ACI/reference/aci_predprey_model.md)
  : Noisy predator-prey model
- [`aci_predprey_components()`](https://biometryhub.github.io/ACI/reference/aci_predprey_components.md)
  : Conditional Gaussian components of the noisy predator-prey model

## The nonlinear dyad model

The worked-example dyad model with intermittent extreme events, at the
components level, and the components schema it illustrates.

- [`aci_components`](https://biometryhub.github.io/ACI/reference/aci_components.md)
  : Conditional Gaussian components list
- [`aci_dyad_components()`](https://biometryhub.github.io/ACI/reference/aci_dyad_components.md)
  : Conditional Gaussian components of the nonlinear dyad model

## Methods

Reporting methods for the model and result objects. Inspecting, plotting
or exporting a result should never require indexing into it.

- [`print(`*`<aci>`*`)`](https://biometryhub.github.io/ACI/reference/print.aci.md)
  : Print an assimilative causal inference result
- [`print(`*`<aci_model>`*`)`](https://biometryhub.github.io/ACI/reference/print.aci_model.md)
  : Print a conditional Gaussian nonlinear system model
- [`summary(`*`<aci>`*`)`](https://biometryhub.github.io/ACI/reference/summary.aci.md)
  : Summarise an assimilative causal inference result
- [`print(`*`<summary.aci>`*`)`](https://biometryhub.github.io/ACI/reference/print.summary.aci.md)
  : Print a summary of an assimilative causal inference result
- [`as.data.frame(`*`<aci>`*`)`](https://biometryhub.github.io/ACI/reference/as.data.frame.aci.md)
  : Coerce an assimilative causal inference result to a data frame
- [`plot(`*`<aci>`*`)`](https://biometryhub.github.io/ACI/reference/plot.aci.md)
  : Plot an assimilative causal inference result
- [`print(`*`<aci_cir>`*`)`](https://biometryhub.github.io/ACI/reference/print.aci_cir.md)
  : Print a causal influence range
- [`summary(`*`<aci_cir>`*`)`](https://biometryhub.github.io/ACI/reference/summary.aci_cir.md)
  [`print(`*`<summary.aci_cir>`*`)`](https://biometryhub.github.io/ACI/reference/summary.aci_cir.md)
  : Summarise a causal influence range
- [`as.data.frame(`*`<aci_cir>`*`)`](https://biometryhub.github.io/ACI/reference/as.data.frame.aci_cir.md)
  : Coerce a causal influence range to a data frame
- [`plot(`*`<aci_cir>`*`)`](https://biometryhub.github.io/ACI/reference/plot.aci_cir.md)
  : Plot a causal influence range

## Package

- [`aciR`](https://biometryhub.github.io/ACI/reference/aciR-package.md)
  [`aciR-package`](https://biometryhub.github.io/ACI/reference/aciR-package.md)
  : aciR: Assimilative Causal Inference for Conditional Gaussian
  Nonlinear Systems
