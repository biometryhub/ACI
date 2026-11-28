# Intent: label a range row cut by `max_lag` as censored when it has not decayed

Author: the `acir` authors. Status: draft 2026-11-24.

## Problem

`aci_range()` on a table built with `lag_table(max_lag = )` labels most rows
`"resolved"`. On `aci_dyad_model()` simulated with seed 1 to `t_end = 3`
(301 anchors), `max_lag = 5` gives 290 resolved rows against 28 on the full
table. Anchor 100 reads tau 0.037, resolved; the full table gives 0.557,
censored. Its last kept divergence is 0.0154, far above the strength floor of
1e-5.

The status check runs on the row after it is padded with zeros to the end of
the record. The padding puts the last divergence above the floor in the first
half of the row, so the check reports that the row decayed.

## Proposed outcome

A row cut by `max_lag` whose last kept divergence is above the strength floor
has status `"censored"`. A row that reaches the end of the record already
gets `"censored"` in that case, since its last exceedance is then its last
cell. Every other row keeps its current status, and every range, strength,
tail estimate and subjective range is unchanged.

## Affected code and users

`.cir_status()` and `aci_range.lag_table()` in `acir/R/aci-core.R`. Users who
build a table with a finite `max_lag` and read `status` from its range,
directly or through `aci(table = )`.

## Constraints

No number changes. Status changes only on rows cut by `max_lag`. No new
status level. Tables saved by earlier releases keep working.

## Open questions

None.
