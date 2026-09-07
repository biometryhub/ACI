# Plan: the performance work package (from intent.md 2026-09-01)

Status: approved 2026-09-01 by Aidan Moller, on v0.3 of the plan, after his
review of v0.2. v0.1 was written on 2026-08-31.

Evidence files named below as `comparison/...` are the authors' working
records and are not in this repository; the numbers they support are restated
in the pull-request descriptions and in `tools/bench/baseline.csv`.

## Order of work

Each step is measured by the one before it, so the order is not negotiable.

| Step | Change | Closes | Lane |
|:---|:---|:---|:---|
| PR-4 | Performance tracking: a benchmark script over the Section 8 stages on the reference records, a committed baseline, a CI job that records and warns at +25% | the instrument | Max |
| PR-5 | Scalar range rows vectorised: the cumulative-log row construction of `aciR` ported into the streamed range and the scalar lag table, cumulative sums blocked every 512 steps; on the scalar path the adaptive truncation becomes a cut by index with the tail estimate kept; the 79 exact-equality assertions restated to tolerance and the resummation declared; the worst case added to the tests | mechanism 1 on the scalar path, the 270x; the scalar half of mechanism 2 | Aidan's `R/`, from a draft branch with the ported kernels and tests |
| PR-6 | Adaptive freeze made cheap on the matrix path: no bound recursions when nothing can freeze, no numerical core touched | mechanism 2 | Aidan |
| PR-7 | Realiser cache on the generic-model route, keyed on the specification's fingerprint | mechanism 3 | Aidan |
| PR-8 | Batch axis for matrix systems: all active anchors' update matrices advanced per step as one array operation | mechanism 1 where logarithms do not apply | joint |
| PR-9 | Climate-model filter to its budget: per-step allocations removed from the matrix kernel | the last budget | Aidan |
| PR-10 | Budgets become gates: the warnings of PR-4 fail the build, the baseline re-committed | regression-proofing | Max |

The 0.1.0 tag sits after PR-5.

## Files that change

| Step | Files |
|:---|:---|
| PR-4 | `tools/bench/bench_reference.R` (new), `tools/bench/baseline.csv` (new), `.github/workflows/bench.yaml` (new) |
| PR-5 | `acir/R/aci-cir-rows.R` (new), `acir/R/aci-cir.R`, `acir/R/aci-online-smoother.R`, `acir/tests/testthat/test-32-scalar-rows.R` (new), the four test files holding the 79 assertions, `acir/inst/evidence/register.csv`, `acir/NEWS.md` |
| PR-6 to PR-9 | named in each step's own plan when it is written; the areas are the matrix-path freeze, the realiser, the matrix kernel |
| PR-10 | `.github/workflows/bench.yaml`, `tools/bench/baseline.csv` |

## Risks

- The resummation moves bits. Counted: 79 exact-equality assertions across
  test-04 (4), test-18 (19), test-20 (8) and test-23 (48). They are restated
  as tolerance assertions and the change is declared.
- Cumulative logarithms cancel on long records. Measured; blocked sums
  restarted every 512 steps err by 3.6e-15 where the plain form errs by
  1.6e-12 at N = 20,000. The blocked form is used and the worst case is
  tested.
- The CI runner is slower than the development machine, so seconds-budgets
  written on one fail on the other. Timing is relative to the runner's own
  baseline.
- A check that passes on the machine it was written on may not pass
  elsewhere; every step runs the full six-platform matrix before it merges.

## Proof

1. **Numbers.** Every fixture in `acir/inst/evidence/register.csv` at its
   tolerance, and the package's stored outputs on the reference records
   before and after the change within 1e-12. A change that moves a number
   beyond round-off is not a performance change and is refused as one.
2. **Time.** The PR-4 table on the CI runner against the committed baseline,
   relative to the runner. Warn at +25% until PR-10, fail after.
3. **Hygiene.** `R CMD check --as-cran` at zero errors and zero warnings,
   coverage not below the floor, the evidence register updated where a check
   method changed.
4. **Provenance.** Every ported kernel names its source (`aciR` 0.2.3, file
   and function) in the roxygen and in `NEWS.md`, against the tag
   `parents-final`.

## Changes from v0.2, after review

The review of v0.2 (2026-09-01) raised four points. Each was checked against
the tree before being taken.

1. The vignette was not the reason to gate the tag: the introductory vignette
   subsamples to 41 anchors on a 1,001-point record and builds in about three
   seconds in CI. Taken: the tag is placed on the budget argument.
2. PR-5 moves bits. Taken: the 79 assertions are restated as tolerance
   assertions inside PR-5, and the resummation is declared.
3. Cumulative logs need a cancellation test. Taken: the blocked form is in
   PR-5's scope and the worst case is in its test set, with the measurements
   recorded in `spec.md`.
4. PR-5 and PR-7 overlapped on the scalar path. Taken: the scalar cut is
   folded into PR-5, and the matrix-path freeze fix moves ahead of the
   realiser cache.

## What approval means

Approving this plan approves the order 4 to 10 as revised, the four gates,
the PR-5 lane arrangement, and the placement of the 0.1.0 tag after PR-5.
Nothing in it changes a number, an interface, or the specification's scope
rule.

## Departures

Recorded 2026-09-02.

- **PR-5b added**, stacked on PR-5: the online smoother's auxiliaries on the
  scalar path formed as vector arithmetic by the same construction, so that
  the lag table and the online smoother share one row builder. Not in the
  approved order; carried as its own pull request with the same four gates.
- **A CI chore added**: every push to `acir-package` runs the full check
  matrix, and the site is published from `main` only. Not in the approved
  order; no numerical code touched.

## Status

| Step | Pull request | State |
|:---|:---|:---|
| PR-4 | #6 | merged 2026-09-01 |
| PR-5 | #7 | open, checks green |
| PR-5b | #8 | open, stacked on #7, checks green |
| CI chore | #9 | open, checks green |
| PR-6 to PR-10 | | not started |
