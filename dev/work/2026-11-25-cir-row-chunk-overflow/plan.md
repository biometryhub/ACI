# Plan: cap the chunk growth of a scalar row (from intent.md 2026-11-25)

Status: draft 2026-11-25.

## Files that change

- `acir/R/aci-cir-rows.R`: `.cir_scalar_row()` stops multiplying the chunk
  length once it covers the row.
- `acir/tests/testthat/test-51-cir-row-chunk.R`: new.
- `acir/NEWS.md`: one bullet under 0.2.2, Bug fixes.

## Order of work

1. Reproduce the stop and the leaked warning on `ba3559b` and write the
   tests; they fail there.
2. Cap the chunk length at the row length `C`, in integer arithmetic.
3. Run the new tests, the full suite, and the bit-identity comparison of
   tables and ranges against `ba3559b`. The reference records are shorter
   than the first chunk, so also compare rows of 100,000 cells from the old
   and the new kernel with `identical()`.

## Risks

- A chunk boundary moves and a carried sum changes in its last bit. A chunk
  of `C` cells or more always ends at the last cell of the row, so the cap
  moves no boundary on a row that completed before. Step 3 measures it.
- The cap itself overflows. It compares the length with `C %/% 4L` and
  multiplies only at or below that.
- A long run of floored cells costs one pass over the rest of the row per
  cell. That is the cost of the existing per-cell rule and is not changed.

## Proof

The new test file: a row with a run of 40 floored cells against a
cell-by-cell recursion that floors each cell, with equal event counts, with
and without `max_lag`; `lag_table()` and `aci_range()` on an `aci()` result
for each `keep`, on the dyad record of the intent and on a model built from
closures; a row that takes exactly ten passes raising no base R warning;
`regularize = "none"` still stopping. The full suite and the two identity
comparisons of step 3.

## Departures

Empty at approval.
