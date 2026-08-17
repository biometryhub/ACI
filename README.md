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

This repository holds implementations of the method together with the evidence
that they compute it correctly. The R implementation is the first;
the layout leaves room for others beside it.

## Implementations

| Path | Language | Status |
|---|---|---|
| [`aciR/`](aciR/) | R | Version 0.2.2. 778 tests, `R CMD check --as-cran` clean. |

```r
# install.packages("remotes")
remotes::install_github("biometryhub/ACI", subdir = "aciR")
```

Start at the package's own [README](aciR/README.md), then the article
*Assimilative causal inference on the nonlinear dyad model*.

## How the implementations are graded

A reimplementation checked only against fixtures its own author transcribed can
demonstrate self-consistency and nothing more. If the author misread an
equation, the fixture encodes the same misreading and passes.

Every numerical claim here is therefore graded against a source that did not
also produce the code under test. The authors publish their work as MATLAB
scripts rather than as a callable library, so their computational passages are
hoisted into callable functions as byte-exact slices, never retyped, and a
checker fails on a single byte of drift. Both sides then run on the same data
and are compared quantity by quantity.

| | |
|---|---|
| Extract outputs reproducing the reference's own workspace | 156 of 156, difference exactly 0 |
| Verbatim reference lines, none retyped | 2,159 across 33 extracts |
| R implementation against the reference, scalar core | 26 of 26 quantities, worst 1.5e-14 |
| Predator-prey, both causal directions, filter through influence range | 26 graded quantities, worst 1.5e-12 |

For scale, an experiment comparing independently developed implementations of
the same algorithms on identical input found agreement degrading from six
significant figures to one (Hatton 1997,
[doi:10.1109/99.609829](https://doi.org/10.1109/99.609829)).

## Documents

The project website is <https://biometryhub.github.io/ACI/>. It carries the
function reference, the three articles, and two ledgers that are published
nowhere else.

| Document | Answers |
|---|---|
| [Parity ledger](https://biometryhub.github.io/ACI/ledgers/parity_ledger.html) | Whether the numbers are right. The reference and the implementation side by side, quantity by quantity. |
| [Development ledger](https://biometryhub.github.io/ACI/ledgers/development_ledger.html) | How the software came to be. Design decisions, the alternatives set aside, what the defects taught, and what remains open. |

Both are self-contained pages, and their sources are in
[`aciR/pkgdown/assets/ledgers/`](aciR/pkgdown/assets/ledgers/). They live there
rather than in a top-level `docs/` directory because they are HTML, which
GitHub serves as source rather than as a page, so the rendered site is the only
place they can actually be read.

## Layout

| Path | Contents |
|---|---|
| [`aciR/`](aciR/) | The R package, and under `pkgdown/` the sources of the project website. |
| [`tools/oracle/`](tools/oracle/) | The MATLAB harnesses that generate the validation fixtures, and the byte-exact parity harness that grades against the authors' own code. |
| [`tools/design/`](tools/design/) | Design records, review rounds and audits, dated and kept as written. |
| [`tools/ledger/`](tools/ledger/) | Decision records. Twelve, each with the alternatives considered and the cost the choice carries forward. |

`tools/` is excluded from the R build, so none of it enters the package
tarball. It is kept in the repository because a claim and its evidence should
travel together.

## Contributing

See [`aciR/CONTRIBUTING.md`](aciR/CONTRIBUTING.md). Reports of numerical
disagreement are the most useful contribution, and a report that names the
system, the parameters and the observed difference can be acted on directly.

## Citation

`citation("aciR")` after installation, or
[`aciR/inst/CITATION`](aciR/inst/CITATION). The method is the authors' and
should be cited alongside any use of this software.

## Security

See [`SECURITY.md`](SECURITY.md) for how to report a vulnerability privately.
Reports of numerical error are handled as defects rather than as
vulnerabilities, and [`aciR/SUPPORT.md`](aciR/SUPPORT.md) says what to include.

## Licence

MIT, Copyright (c) 2026 Max Moldovan. See [`LICENSE`](LICENSE). The licence
covers the whole repository: the package, the harnesses in `tools/`, and the
published documents.

The reference implementation
([marandmath/ACI_code](https://github.com/marandmath/ACI_code)) is separately
licensed MIT, Copyright (c) 2025 Marios Andreou, and is not redistributed here.
Excerpts quoted in the parity ledger carry that notice beside the code they
quote, and the copyright holder is named in the package metadata and in both
licence files.
