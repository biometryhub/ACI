# Intent: refuse inputs that are ignored without a condition

Author: the `acir` authors. Status: draft 2026-11-24.

## Problem

Two inputs are accepted and then ignored, so a call runs with settings the
user did not ask for.

- `aci(table = )` uses the table's own construction and ignores `stepper`,
  `nsub` and `regularize`. It also accepts `nsub = "x"`. On the dyad model
  with `t_end = 2`, `dt = 0.01`, a table call with `nsub = 4` reports a peak
  ACI of 3.7170, while the streamed call with `nsub = 4` gives 3.7116.
  `lag_table()` refuses the same `stepper` and `nsub` values.
- A list `init` with a name other than `mean` and `cov` is ignored.
  `init = list(m = 5, cov = matrix(1))` runs with prior mean 0 and no
  condition. This holds for the scalar and the matrix kernels and for every
  verb that takes `init`.

A third case is not part of this change: assigning a coefficient into a
model (`m$Lx <- ...`) changes filtering and not simulation, while the help
for `aci_model()` says both stay on the original coefficients. See the open
questions.

## Proposed outcome

- `aci(table = )` with an explicitly supplied `stepper`, `nsub` or
  `regularize` stops with `aci_error_dims` and says to give these to
  `lag_table()`.
- An `init` list with an unknown or empty name stops with `aci_error_dims`,
  naming the unknown entries and the accepted ones. The accepted names are
  `mean` and `cov`, and `source`, which the package writes into the prior it
  stores in `meta$init`, so a stored prior can be passed back.
- Every call accepted before returns the same numbers.

## Affected code and users

`acir/R/aci-core.R` (`aci()`), `acir/R/aci-assimilation.R`
(`.same_gaussian_init()` and a new name check), `acir/R/aci-kernels-scalar.R`
and `acir/R/aci-kernels-matrix.R` (the two prior checks), their help, a new
test file and `acir/NEWS.md`. Users who pass these arguments with a reused
table, or who misspell an entry of `init`.

## Constraints

No change to any value on an accepted input, to any fixture or tolerance, or
to any existing test.

## Open questions

- Coefficient assignment into a model. A `$<-` method that refuses it makes
  six existing tests error, because those tests assign `fx` or `g` into a
  model on purpose to exercise the realiser identity and the validators
  (`tests/testthat/test-02-models.R`, `test-15-compiled-scalar.R`,
  `test-17-compiled-conditioning.R`, `test-25-enso6-batch-realiser.R`). The
  authors decide between changing those tests to build the altered model
  another way, or correcting the help paragraph to describe what assignment
  does.
