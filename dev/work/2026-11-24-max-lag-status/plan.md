# Plan: status of rows cut by `max_lag` (from intent.md 2026-11-24)

Status: draft 2026-11-24.

## Files that change

- `acir/R/aci-core.R`: `.cir_status()` takes the number of cells kept by a
  `max_lag` cut and returns `"censored"` when the last kept cell is above the
  floor; `aci_range.lag_table()` passes that number for rows shorter than the
  record and longer than `meta$max_lag`.
- `acir/tests/testthat/test-47-max-lag-status.R`: new.
- `acir/NEWS.md`: one bullet under 0.2.2, Bug fixes.

## Order of work

1. Confirm the defect on `ba3559b` and write the tests; they fail there.
2. Change `.cir_status()` and its caller.
3. Run the new tests, the full suite and the bit-identity comparison of
   tables and ranges against `ba3559b`.

## Risks

- A row cut by the `tol` rule could change status. That rule cuts only after
  `window` cells whose divergence plus tail estimate is below `tol`, so its
  last kept cell is below `tol`. Only rows whose kept length exceeds
  `max_lag` are treated as cut by it.
- A table saved before `meta$max_lag` existed: read as `Inf`, so no row
  changes.
- The streamed route on an `aci()` result builds its rows with no lag limit
  and is not changed.

## Proof

The new test file: the reported case, the full table unchanged, every number
bit-identical to the same table under the old rule, a cut row below the floor
staying resolved, one and two hidden components, `aci(table = )`, a
`saveRDS()`/`readRDS()` round trip, and a `tol` cut with no `max_lag`
unchanged. The full suite and the series identity comparison.

## Departures

- `acir/vignettes/vignette-2-advanced.Rmd` also changes. Its list of statuses
  said that censoring is judged inside the row and never against a horizon,
  which the new rule for a row cut by `max_lag` makes incomplete. The
  `censored` entry and the paragraph after the list now state the rule for a
  cut row. No code chunk changes.
