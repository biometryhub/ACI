# acir

[![R-CMD-check](https://github.com/biometryhub/ACI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/biometryhub/ACI/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://biometryhub.github.io/ACI/API_STABILITY.md)
[![Licence:
MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](https://biometryhub.github.io/ACI/LICENSE.md)

acir is a closed-form assimilative causal inference engine for
conditional Gaussian nonlinear systems, consolidating the `aci` and
`aciR` packages into a single implementation.

Its outputs are graded against hash-pinned reference records. The
filter, the smoother and the total ACI metric are compared with the
authors’ published MATLAB outputs for the dyad and predator-prey
systems, and the fixed-lag online smoother and the causal influence
range with those outputs for the dyad record only; the paths those
scripts do not exercise (the noise cross-covariance terms, the matrix
online smoother, the general auxiliary matrices) are compared with
independent MATLAB or R transcriptions of the published equations, with
analytic identities, or with a source-derived run.
`inst/evidence/register.csv` names the comparator, the source class and
the tolerance for each checked feature.

## Status

Version 0.1.0 is the parity milestone for the graded surface. The
evidence register carries 68 rows: 28 are graded against a hash-pinned
fixture and 40 record an exact relation, a behavioural check or an
in-test independent transcription with no fixture behind them. Ten of
the fixture-backed rows cite a fixture the method authors produced: nine
compare a computed output against it, and the tenth records that
[`observed_trajectory()`](https://biometryhub.github.io/ACI/reference/observed_trajectory.md)
ingests the authors’ pinned input signal and carries it onto the seven
dyad grades; the two predator-prey grades run on the authors’ pinned
predator-prey record, which has no register row of its own. The
remaining fixture-backed rows compare against source-derived,
independently transcribed or analytic references, each named in
`inst/evidence/register.csv` with its tolerance. The stages of the
specification’s performance table are timed against the committed
baseline by `tools/bench/bench_reference.R`.

Maintainer: Aidan Moller. Authors: Aidan Moller and Max Moldovan.

The numerical core is fixed: a change that moves a graded number beyond
round-off is a change of method, not a release. The public interface may
still change before 1.0, and every such change is announced in
`NEWS.md`. Details in `API_STABILITY.md`.

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
