# ACI

Assimilative causal inference for conditional Gaussian nonlinear systems,
implementing the method of Andreou, Chen and Bollt (2026),
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0).

The method asks how much the future of an observed signal tells us about the
present state of something we cannot observe. For a conditional Gaussian
nonlinear system, the forward filter and the backward smoother are both
available in closed form, and the relative entropy between them, evaluated at
each instant, is the causal-information metric. A second quantity, the causal
influence range, measures how far ahead one must look before that answer stops
changing.

This repository holds the R implementation together with the evidence that it
computes the method correctly.

## The package

| Path | Contents |
|---|---|
| [`acir/`](acir/) | The R package `acir`, the closed-form engine: filter, smoother, causal measure, fixed-lag online smoother, influence range, conditional questions, and the benchmark systems of the reference MATLAB code. |

```r
# install.packages("remotes")
remotes::install_github("biometryhub/ACI", subdir = "acir")
```

`acir` consolidates two earlier packages by the same authors: `aciR` (Max
Moldovan), whose numerical core was graded against the authors' MATLAB to
round-off, and `aci` (Aidan Moller), which contributed the model interface,
the one-time evaluation of model coefficients, the influence-range machinery
and the benchmark systems of the later papers. Both are retained in the
history at the tag `parents-final`, and the merged package is compared with
each of them on every shared quantity as part of its checks. The package is
under active joint development; authorship and citation metadata are
interim until that review closes.

## How the implementation is graded

A reimplementation checked only against fixtures its own author transcribed can
demonstrate self-consistency and nothing more. If the author misread an
equation, the fixture encodes the same misreading and passes.

Every numerical claim here is therefore graded against a source that did not
also produce the code under test. The authors publish their work as MATLAB
scripts rather than as a callable library, so their computational passages are
hoisted into callable functions as byte-exact slices, never retyped, and a
checker fails on a single byte of drift. Both sides then run on the same data
and are compared quantity by quantity. The package ships the register of what
is checked, how, and against what (`acir/inst/evidence/register.csv`), and a
test fails the build if an exported function has no row or a row names a
fixture whose bytes have moved.

For scale, an experiment comparing independently developed implementations of
the same algorithms on identical input found agreement degrading from six
significant figures to one (Hatton 1997,
[doi:10.1109/99.609829](https://doi.org/10.1109/99.609829)).

## Layout

| Path | Contents |
|---|---|
| [`acir/`](acir/) | The R package. |
| [`tools/oracle/`](tools/oracle/) | The MATLAB harnesses that generate the validation fixtures, and the byte-exact parity harness that grades against the authors' own code (`ACI_code`, and the backward influence range of `FBCIR_code`). |
| [`tools/design/`](tools/design/) | Design records, review rounds and audits, dated and kept as written. |
| [`tools/ledger/`](tools/ledger/) | Decision records, each with the alternatives considered and the cost the choice carries forward. |
| [`dev/`](dev/) | Development records of the consolidation, and the staged material of the later paper families that is not yet in the package. |

`tools/` and `dev/` are excluded from the R build, so none of it enters the
package tarball. They are kept in the repository because a claim and its
evidence should travel together.

## Documents

The project website is <https://biometryhub.github.io/ACI/>, built from
`acir/` by continuous integration.

## Contributing

Reports of numerical disagreement are the most useful contribution, and a
report that names the system, the parameters and the observed difference can
be acted on directly. Open an issue at
<https://github.com/biometryhub/ACI/issues>.

## Citation

`citation("acir")` after installation. The method is the authors' and should
be cited alongside any use of this software:

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
*Nature Communications*, 17, 1854.

## Security

See [`SECURITY.md`](SECURITY.md) for how to report a vulnerability privately.
Reports of numerical error are handled as defects rather than as
vulnerabilities.

## Licence

MIT. See [`LICENSE`](LICENSE) for the repository and
[`acir/inst/COPYRIGHTS`](acir/inst/COPYRIGHTS) for the notices of the
reference implementations the package is verified against.
