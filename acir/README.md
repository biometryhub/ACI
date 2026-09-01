# acir

acir is a closed-form assimilative causal inference engine for conditional
Gaussian nonlinear systems, consolidating the `aci` and `aciR` packages into a
single implementation.

Its outputs are verified number by number against the authors' published MATLAB
reference implementation, and the evidence documents accompany this repository.

## Status

acir is under active joint development. Authorship, citation and the maintainer
field are under joint review, and nothing in `DESCRIPTION` pre-empts the outcome
of that review.

## Installation

From a source checkout of this repository (the package lives in the
`acir/` subdirectory):

```sh
R CMD INSTALL acir
```

or directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("biometryhub/ACI", subdir = "acir")
```

## A first run

```r
library(acir)
model <- aci_dyad_model()
sim <- aci_simulate(model, T = 5, dt = 0.005, seed = 1)
fit <- aci(model, sim)
fit
plot(fit)
```

## Licence

MIT; see `LICENSE` and `inst/COPYRIGHTS` for the notices of the reference
implementations the package is verified against.
