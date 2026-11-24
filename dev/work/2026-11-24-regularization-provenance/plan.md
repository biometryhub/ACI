# Plan: carry a floored filter's record (from intent.md 2026-11-24)

Status: draft 2026-11-24.

## Files that change

- `acir/R/aci-utils.R`: new internal `.aci_reg_carry(rec, reg)`, which adds
  a fired frozen record to a recorder: its sites, its event count and policy
  `"floor"`. A record that did not fire, or is `NULL`, leaves the recorder
  untouched. The recorder's `strict` flag is not changed.
- `acir/R/aci-online-smoother.R`: `.lag_table_compiled()` carries the
  supplied filter's record when it keeps that filter and created its own
  recorder. A recorder passed in by `aci(keep = "table")` already holds the
  filter's events, so it is not carried twice.
- `acir/R/aci-cir.R`: `.forward_cir_compiled()` does the same for the
  retained filter that `aci_range()` passes from an `aci_result`.
- `acir/R/aci-assimilation.R`: `aci_smoother()` and `aci_online()` carry a
  supplied filter's record. Help for `regularize` on `lag_table()`,
  `aci_smoother()` and `aci_online()` says so.
- `acir/man/`: the three pages, regenerated.
- `acir/tests/testthat/test-46-regularization-provenance.R`: new.
- `acir/NEWS.md`: one bullet under 0.2.2 Bug fixes.

## Order of work

1. Write the test file and confirm it fails on `ba3559b`.
2. Add the helper and the four call sites.
3. Regenerate the help and keep only the changed pages.
4. Run the full suite, the base run of the new file and the prose guard.

## Risks

- Double counting when a recorder is shared. Carrying only where the
  recorder is created in the same call rules it out; the test checks the
  event count on `aci(keep = "table")` and on each consuming route.
- A changed number. The test checks that a table built on a supplied
  unfloored filter is identical to one that builds its own filter, and that
  objects built on a floored filter differ from their base form only in
  `meta$regularization`. A digest of every numeric field of each object, on
  the base and on the branch, is compared outside the package.
- An existing test that asserts silence. The reuse routes keep their
  contract, so none should change.

## Proof

`test-46-regularization-provenance.R`: the dyad (one hidden component), the
ENSO library model with two hidden components, and a generic `aci_model()`
from closures; each consuming route; a fresh and a `saveRDS()`/`readRDS()`
filter, table and result; exact warning counts; bit-identity of the
unfloored objects. The full suite.

## Departures

Empty at approval.
