# Getting help with aciR

## Reading first

Three articles ship with the package, and between them they answer most
questions.

- *Assimilative causal inference on the nonlinear dyad model* walks the
  whole workflow on a worked example.
- *Assumptions and interpretation* states the conditions under which the
  method applies, and what a peak in the metric does and does not
  support.
- *Validation and the independent oracle* records what each oracle
  grades and, for each, what it does not.

`vignette(package = "aciR")` lists them. Function documentation is
available through
[`?aci`](https://biometryhub.github.io/ACI/reference/aci.md),
[`?aci_cir`](https://biometryhub.github.io/ACI/reference/aci_cir.md) and
the rest.

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
treated as findings rather than as support requests.

Please include the system, the parameters, the quantity that disagrees,
and the size of the disagreement. If you have compared against another
implementation, say which one and at what settings. Conventions differ
between implementations in ways that look like disagreements and are
not, so naming the settings matters.

`tools/oracle/parity/` in the repository holds the harness used to grade
this package against the authors’ reference, and its report records the
conventions under which each comparison was made.

## Reporting a defect

Also the issue tracker. A defect report is most useful with the output
of [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), the call
that produced the problem, and what you expected instead.
