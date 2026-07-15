# aciR

Assimilative causal inference (ACI) for conditional Gaussian nonlinear systems
in R, reimplementing the method of Andreou, Chen and Bollt (2026),
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0).

This repository holds the R package and the independent-oracle harnesses that
validate it. The package documentation is in **[`aciR/`](aciR/)**.

```r
# install.packages("remotes")
remotes::install_github("max578/aciR", subdir = "aciR")
```

## Layout

| Path | Contents |
|---|---|
| [`aciR/`](aciR/) | The R package. Start at its [README](aciR/README.md). |
| [`oracle/`](oracle/) | MATLAB harnesses that generate the validation fixtures, and the fixtures themselves. |

The package is a subdirectory rather than the repository root so that the
oracle harnesses can live beside it: the fixtures they produce are the
package's evidence, and `aciR/inst/extdata/oracle-manifest.yml` records the
provenance of each one, including the harness that generated it.

## Validation

The numerical core is graded against fixtures produced by an independent MATLAB
harness reproducing the deterministic core of the authors' reference
implementation ([marandmath/ACI_code](https://github.com/marandmath/ACI_code),
MIT). aciR runs its own filter, smoother and metric on the reference signal and
must reproduce all five output series to a maximum absolute error below `1e-6`;
the observed error is `4.6e-14`. The gate never skips.

Because every scalar model in that reference sets the noise cross-covariance to
zero, a second harness and a set of analytic Kalman-Bucy identities grade the
terms the reference cannot reach. See the package vignette *Validation and the
independent oracle* for the design and, for each oracle, what it does and does
not cover.

## Licence

MIT (c) 2026 Max Moldovan. The MATLAB reference implementation is separately
licensed by its authors (MIT) and is not redistributed here.
