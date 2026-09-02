# acir

[![R-CMD-check](https://github.com/biometryhub/ACI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/biometryhub/ACI/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://biometryhub.github.io/ACI/API_STABILITY.md)
[![Licence:
MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](https://biometryhub.github.io/ACI/LICENSE.md)

acir is a closed-form assimilative causal inference engine for
conditional Gaussian nonlinear systems, consolidating the `aci` and
`aciR` packages into a single implementation.

Its outputs are verified number by number against the authors’ published
MATLAB reference implementation, and the evidence documents accompany
this repository.

## Status

Version 0.1.0 is the parity milestone: every quantity the reference
MATLAB implementation computes is reproduced to the tolerance recorded
in the evidence register, and the performance table of the specification
is met. Maintainer: Aidan Moller. Authors: Aidan Moller and Max
Moldovan.

The numerical core is fixed at parity with the reference implementation;
the public interface may still change before 1.0, and every such change
is announced in `NEWS.md`. Details in `API_STABILITY.md`.

## Installation

From a source checkout of this repository (the package lives in the
`acir/` subdirectory):

``` sh
R CMD INSTALL acir
```

or directly from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("biometryhub/ACI", subdir = "acir")
```

## A first run

``` r

library(acir)
model <- aci_dyad_model()
sim <- aci_simulate(model, t_end = 5, dt = 0.005, seed = 1)
fit <- aci(model, sim)
fit
plot(fit)
```

## Licence

MIT; see `LICENSE` and `inst/COPYRIGHTS` for the notices of the
reference implementations the package is verified against.
