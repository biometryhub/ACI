# Intent: finish a scalar row that the floor policy splits into many passes

Author: the `acir` authors. Status: draft 2026-11-25.

## Problem

Under `regularize = "floor"`, `lag_table()` and `aci_range()` on an `aci()`
result stop with the base R error "missing value where TRUE/FALSE needed"
when the floor fires at several lags of one anchor, for a model with one
hidden component. On `aci_dyad_model()` simulated with seed 1 to
`t_end = 1.5`, thinned to `dt = 0.1` (16 points) and started from a prior
variance of 0.1, the row of the third anchor has 13 cells and the variance is
not positive at cells 2 to 13. Both calls stop at its eleventh cell.

`.cir_scalar_row()` forms a row in chunks and multiplies the chunk length by
four after every pass. A cell that needs the covariance policy ends its chunk
and is formed alone, so each such cell takes one pass, and each stretch of
cells between two of them takes one more. The length starts at 2,048 and the
tenth multiplication gives 2^31, which an integer cannot hold. A row that
needs an eleventh pass stops. A row that needs exactly ten returns, with the
base R warning "NAs produced by integer overflow": a linear model with 13
points at `dt = 0.12` does this.

## Proposed outcome

Both calls return under `regularize = "floor"` however many passes a row
takes. They raise one `aci_warn_regularized`, record every event in
`meta$regularization`, and raise no base R warning. A floored row equals the
cell-by-cell recursion under the same policy. Every row that was returned
before is bit-identical.

## Affected code and users

`.cir_scalar_row()` in `acir/R/aci-cir-rows.R`, reached from `lag_table()`
and from `aci_range()` on an `aci()` result, for models with one hidden
component. Users who call either with `regularize = "floor"` on a record
whose step is too coarse for the explicit scheme.

## Constraints

No number changes on any input that returned before. The running sums are
carried across chunk boundaries, so the boundaries stay where they are.
`regularize = "none"` still stops at the first such cell with
`aci_error_spd`. Models with several hidden components use a separate
recursion and are not touched.

## Open questions

None.
