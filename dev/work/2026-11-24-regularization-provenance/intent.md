# Intent: carry a floored filter's record into what is built from it

Author: the `acir` authors. Status: draft 2026-11-24.

## Problem

`aci_filter(regularize = "floor")` records its floor events in
`meta$regularization` and raises `aci_warn_regularized`. A call that is given
that filter starts a new, empty record. On the dyad model with prior variance
20, `dt = 0.01` and one second of record, the filter takes one floor.
`lag_table(m, ob, filter = f)` then raises nothing and records policy
`"none"` with no event, and `aci(m, ob, table = tb)` returns a peak ACI of
5.758e+24 whose record says nothing fired. `aci_online(filter = f)` does the
same. `aci_range()` on `aci(regularize = "floor", keep = "paths")` reuses the
retained filter and reports zero events, where the same result with
`keep = "none"` reports the floor. The 0.2.0 notes say a floored result
always says so.

## Proposed outcome

A call that consumes a supplied filter starts its record from the filter's
record. A lag table, online path or smoother built on a floored filter, a
range computed from a retained floored filter, and everything later built
from them, carry policy `"floor"` and every floor event of the filter, and
the call raises one `aci_warn_regularized`. The record survives `saveRDS()`
and `readRDS()`. No number changes, and a call whose filter took no floor
behaves exactly as before.

## Affected code and users

`acir/R/aci-utils.R` (the recorder), `acir/R/aci-online-smoother.R`
(`.lag_table_compiled()`), `acir/R/aci-cir.R` (`.forward_cir_compiled()`)
and `acir/R/aci-assimilation.R` (`aci_smoother()`, `aci_online()`). Users who
pass a floored filter to `lag_table()`, `aci_online()` or `aci_smoother()`,
and users who call `aci_range()` on a floored `aci()` result that kept its
paths.

## Constraints

No change to any returned number, interface, fixture or tolerance. A call's
own policy still governs the recursions it runs. An existing test that
changes is a stop.

## Open questions

Whether reusing a floored table in `aci(table = )` or `aci_range()` should
warn again. The current tests assert that reuse carries the record and its
printed line and raises nothing new; this change keeps that contract. The
authors decide.
