# Plan: correct false statements (from intent.md 2026-11-24)

Status: draft 2026-11-24.

## Files that change

- `acir/R/aci-assimilation.R`: roxygen of `lag_table()` (`stepper`, `nsub`,
  `tol`), `lt_row()`, `as.data.frame.lag_table()` and the data-frame method
  for `da_path_gaussian`.
- `acir/R/aci-core.R`: roxygen of `aci_metric_pair()` and `aci_range()`;
  the message of the `direction = "backward"` error, keeping "not in this
  release".
- `acir/R/aci-model.R`: roxygen of `aci_linear_model()` for `fx` and `fy`.
- `acir/man/`: the pages regenerated from the roxygen above.
- `acir/vignettes/vignette-2-advanced.Rmd`: the `tail_bound` sentence and
  the `mode = "full"` sentence.
- `acir/README.md`: the sentence on what the filter, smoother and metric are
  compared with.
- `acir/tests/testthat/test-03-engine.R`: a comment labelling the
  `aci(table = )` comparison as an identity. The assertion is unchanged.
- `acir/inst/evidence/register.csv`: the `lag_table`, `lt_diag` and
  `aci_metric` rows that cite that comparison, reworded in place.
- `acir/tests/testthat/fixtures/oracles/oracle-manifest.yml`: the
  annotation block states which acir test holds each analytic oracle, and
  that none holds the Riccati fixed point.
- `acir/NEWS.md`: one Documentation bullet under 0.2.2.

## Order of work

1. Reproduce each false statement on `ba3559b`, one script per item.
2. Edit the roxygen, regenerate `man/`, keep only the pages that changed.
3. Edit the vignette and render it against the loaded source.
4. Edit the README, the test comment, the register rows and the manifest.
5. Run the full suite, the prose guard and the register tests.

## Risks

- A reworded message could break a test that matches its text. The only
  such test matches "not in this release", which stays.
- Regenerating `man/` with a newer roxygen rewrites unrelated files. Those
  are restored from the base.
- A register row could be weakened by accident. Each reworded row claims
  only what its test checks, and the before and after are listed in the
  pull request.

## Proof

The full suite with 0 failures and 7,818 expectations, `test-19`, `test-30`
and `test-31`, which read the manifest and the register, a render of
`vignette-2-advanced.Rmd`, and the reproduction scripts and their output.

## Departures

Empty at approval.

### 2026-11-24, NEWS

- **Two bullets.** The documentation corrections and the evidence
  corrections are listed separately under 0.2.2, because the second changes
  what the register claims.
