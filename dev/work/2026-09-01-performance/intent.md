# Intent: bring every Section 8 stage of `acir` inside its performance budget

Author: Max Moldovan. Status: accepted 2026-09-01, with the plan that carried
it, reviewed by Aidan Moller.

Consolidated on 2026-09-02 from the performance plan v0.1 (2026-08-31), where
this intent was first recorded, and from the review of its v0.2.

## Problem

`acir` was set up to take the best of its two parents: the scope and model
interface of `aci`, the speed and evidence of `aciR`. The first half is done:
the package agrees with the authors' MATLAB to 1e-14 on the filter, smoother,
metric, online smoother and influence range. The second half is not. Measured
on 2026-08-31 on the authors' dyad record (N = 3,000):

| Stage | Budget | `acir` | Verdict |
|:---|---:|---:|:---|
| influence range, all anchors | 0.40 s | 106.8 s | fails by about 270x; `aciR` does the same task in 0.28 s with the same numbers to 1.4e-14 |
| filter, model supplied as closures | 0.02 s | 0.142 s | passes only on the library route (0.022 s) |
| smoother, model supplied as closures | 0.02 s | 0.110 s | same (0.009 s on the library route) |
| online smoother, all lags | 0.25 s | 0.264 s | passes, two to five times faster than MATLAB |
| climate model, 20,001 steps, filter | 0.7 s | not measured | open |

## Proposed outcome

Every Section 8 stage of the specification inside its budget on the reference
records, in base R with no compiled inner loop, with no graded number moved
beyond round-off and no interface changed. A release that fails its own range
budget by two orders of magnitude is not the package that was set out to be
built, so the 0.1.0 tag waits for the range.

## Affected code and users

The scalar path of the influence range and the lag table (`R/aci-cir.R`,
`R/aci-online-smoother.R`); the adaptive freeze on the matrix path; the
generic-model realiser; the benchmark under `tools/bench/` and its workflow;
the placement of the 0.1.0 tag. Users notice only speed.

## Constraints

The scope rule (closed forms only). Every fixture in the evidence register at
its tolerance, and the stored outputs on the reference records within 1e-12
before and after. `R CMD check --as-cran` at zero errors and zero warnings.
Every ported kernel names its source in `aciR` 0.2.3 at the tag
`parents-final`. `R/` is in Aidan's lane; ported kernels arrive as a draft
branch.

## Open questions at acceptance

Whether rows built from cumulative logarithms cancel on the 20,001-step
climate record. How many exact-equality assertions the change of summation
order moves. Whether the tag should wait for the range at all. All three were
settled in the review of the plan and are recorded in `plan.md`.
