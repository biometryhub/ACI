# aciR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`aciR` implements assimilative causal inference (ACI) for conditional Gaussian
nonlinear systems. A conditional Gaussian nonlinear system (CGNS) pairs an
observed signal with an unobserved component whose statistics, conditional on
the observed path, are Gaussian and available in closed form. ACI measures,
step by step, how strongly the observed signal carries information about that
unobserved component -- and so where the two are causally coupled.

The workflow rests on data assimilation. A forward filter estimates the hidden
component from the observed path seen so far, a backward smoother estimates it
from the whole observed path, and the causal-information metric -- the relative
entropy of the smoother posterior from the filter posterior -- reads off how
much the future of the observed signal sharpens the estimate at each step.

The method reimplemented here is due to Andreou, Chen and Bollt (2026); the
numerical core is validated against the authors' reference implementation.

## Installation

Install the development version from a local checkout with:

```r
# install.packages("remotes")
remotes::install_local("aciR")
```

## Usage

The flagship example is the nonlinear dyad model with intermittent extreme
events. Build the model, simulate a realisation, then run the whole ACI
workflow in a single call to `aci()`.

```r
library(aciR)

# Build the dyad model and simulate a realisation of it
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 5000, seed = 333)

# Run the filter, the smoother and the causal-information metric
fit <- aci(sim$x, model)
fit
#> <aci> assimilative causal inference
#>   model: nonlinear dyad model with intermittent extreme events
#>   steps: 5000, dt: 0.001, time span: [0, 4.999]
#>   causal-information metric:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#>  0.0000  0.0645  0.1209  0.2518  0.2962  2.2935
```

The returned `aci` object holds the `filter` and `smoother` posteriors of the
hidden component -- each a list of `mean` and `cov` -- and the
causal-information metric in `aci`, one value per time step, ready for
inspection or plotting.

```r
summary(fit$aci)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#>  0.0000  0.0645  0.1209  0.2518  0.2962  2.2935
```

The three core functions -- `aci_filter()`, `aci_smoother()` and
`aci_metric()` -- work on any conditional Gaussian nonlinear system through a
general components list, so they are not tied to the dyad model.
`aci_cgns_model()` builds a model object for a system of one's own. The
vignette *Assimilative causal inference on the nonlinear dyad model* is a full
walkthrough.

## Citation

If you use `aciR` in published work, please cite both the package and the
method paper. See `citation("aciR")` for the current entries.

The method is:

> Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
> *Nature Communications*, 17, 1854.
> <https://doi.org/10.1038/s41467-026-68568-0>

## Licence

MIT (c) 2026 Max Moldovan. See the `LICENSE` file.
