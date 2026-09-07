# Getting help with acir

## Reading first

Three articles ship with the package, and between them they answer most
questions.

- *Assimilative causal inference (ACI) in R* walks the whole workflow on
  a worked example.
- *The closed-form ACI engine* documents the machinery underneath: how
  models are written, what the filter, smoother and online-smoother
  arguments do, and which conventions each reported number is computed
  under.
- *Reproducing the reference MATLAB codebase* records which quantities
  are graded against the authors’ implementation, where the package
  departs from it deliberately, and how far each correspondence has been
  checked.

`vignette(package = "acir")` lists them. Function documentation is
available through
[`?aci`](https://biometryhub.github.io/ACI/reference/aci.md),
[`?aci_range`](https://biometryhub.github.io/ACI/reference/aci_range.md)
and the rest.

## Asking a question

Open an issue at <https://github.com/biometryhub/ACI/issues>.

A question about a result is easiest to answer with the model, the
parameters and the observed output. A short reproducible example is
ideal, and
[`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
takes a `seed`, so a simulated case can be made reproducible in one
line.

## Reporting a numerical disagreement

These are the most valuable reports this package receives, and they are
treated as findings rather than as support requests. The issue tracker
offers a template for them.

Please include the system, the parameters, the quantity that disagrees,
and the size of the disagreement. If you have compared against another
implementation, say which one and at what settings. Conventions differ
between implementations in ways that look like disagreements and are
not, so naming the settings matters.

`tools/oracle/parity/` in the repository holds the harness used to grade
this package against the authors’ reference, and
`inst/evidence/register.csv` records the tolerance each graded quantity
meets, so a disagreement can be read against the tolerance it is
supposed to sit inside.

## Reporting a defect

Also the issue tracker. A defect report is most useful with the output
of [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), the call
that produced the problem, and what you expected instead.
