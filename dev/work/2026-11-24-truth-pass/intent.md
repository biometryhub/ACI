# Intent: make false statements in the help, vignette, README and evidence true

Author: the `acir` authors. Status: draft 2026-11-24.

## Problem

Several statements in the documentation and the evidence files are false
when checked against the 0.2.2 code:

- `vignette-2-advanced.Rmd` says `mode = "full"` keeps negative-lag cells.
  A full table stores the same rows as `mode = "forward", tol = 0`, with
  positive lags only.
- The same vignette says `tail_bound` overstates the error under a hard
  `max_lag` cap. Its own table gives a tail estimate of 0.39 against an
  error of 0.63.
- `?as.data.frame.lag_table` says the columns are the anchor time and the
  lag. They are two grid indices, `j` and `n`, and the divergence `P`.
- `?as.data.frame.da_path_gaussian` says the data frame carries the
  marginal variance. The column is `sd`, its square root.
- `?lt_row` says the row runs at increasing positive lag. Its first element
  is lag 0.
- `?lag_table` offers `stepper = "implicit"` and sub-steps, which it
  refuses, and calls `tol` positive, where 0 is accepted and used.
- `?aci_linear_model` says `fx` and `fy` are functions. Constants are
  accepted.
- `?aci_metric_pair` gives precision figures measured at `R_q = 1` without
  saying that other scales are worse: over `R_q` from 1e-3 to 1e3 the
  relative error at `delta = 1e-07` reaches 0.11, against 0.012 at 1.
- `?aci_range` and the error for `direction = "backward"` point to
  `reserve/fbcir/`, and `simpson_close` to a "pre-0.0.21 rule". Neither
  exists for a user of the package.
- `README.md` says the filter, smoother and metric are compared with "the
  authors' published MATLAB outputs". The fixtures are outputs of MATLAB
  harnesses in `tools/oracle/` adapted from the authors' scripts.
- `test-03-engine.R` compares `lt_diag(lt)` with `aci(table = lt)$aci`,
  which returns the same vector by construction, and three register rows
  cite that comparison as agreement with the ACI series.
- The oracle manifest cites `test-identities.R` for the analytic oracles
  it calls the primary grounding of the cross-noise terms. No such file
  exists here, and no acir test holds the Riccati fixed point.

## Proposed outcome

Each statement says what the code does, with the smallest change of words.
No number, test result, fixture or tolerance moves.

## Affected code and users

Roxygen in `acir/R/aci-assimilation.R`, `aci-core.R` and `aci-model.R`,
with the regenerated pages; `vignettes/vignette-2-advanced.Rmd`;
`acir/README.md`; one comment in `tests/testthat/test-03-engine.R`; three
rows of `inst/evidence/register.csv`, reworded in place; the annotation
block of `tests/testthat/fixtures/oracles/oracle-manifest.yml`. Readers of
the help, the vignette and the README, and reviewers of the evidence.

## Constraints

No change to code paths or numbers. The one reworded error message keeps
its class and the phrase a test matches. The manifest's original text below
its annotation stays as it is, because it is the hash record for the
fixtures. No register row is added or removed.

## Open questions

The README wording of the reference claim is for the co-author to confirm.
