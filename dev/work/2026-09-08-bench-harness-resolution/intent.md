# Intent: make the benchmark harness measure the cost of the stages it times

Author: the `acir` authors. Status: draft 2026-09-08.

## Problem

The Time gate has been reporting a breach it cannot support. `bench.yaml`
warns whenever a stage drifts more than 25 percent from `tools/bench/
baseline.csv` after dividing by an R-loop calibration, and since 2026-09-07 the
stage "online smoother, all lags" has warned on every run: 0.005 s against a
baseline of 0.003 s, which is +25 to +33 percent normalised. Two runs made with
`--gate` exited non-zero on it. The stage's Section 8 budget is 0.25 s.

Measured on 2026-09-08 on the same commit (`1e3d0df`), the breach is in the
instrument.

- `median_time()` (`tools/bench/bench_reference.R:47-51`) times three calls
  with `proc.time()[["elapsed"]]` and takes their median.
- `proc.time()`'s elapsed tick on the measuring machine is 1 ms.
  `Sys.time()`'s is 9.54e-7 s, so the coarseness is the function chosen, not
  the platform. Every seconds value in `baseline.csv` and in every run is an
  exact multiple of 1 ms.
- The stage's first two calls are warm-up. Call by call the cost is 4.0 ms
  (32.7 ms in a bare process), then 42.8 ms (47.3 ms), then flat at 0.84 ms
  from call 3 onward. Call 2 is the maximum and call 3 the minimum, so the
  median of three is call 1: the harness reports a cold call, quantised.
- The warm steady-state cost is 0.828 ms (mean of 20 batches of 50 calls;
  batch median 0.820 ms, IQR 0.040 ms, min 0.800 ms). Per-call at high
  resolution the median is 0.793 ms. That is a ratio of 0.0085 against the
  calibration, 78 percent below the baseline ratio of 0.0395, and about 300
  times inside the stage's budget.

Five ticks against three ticks is a two-tick difference between two cold calls
under two R versions (4.5.2 and 4.6.1); at this magnitude one tick is 20 to 33
percent, which is the whole of the reported drift. There is no code regression
in this stage. Any other stage costing less than a few milliseconds carries the
same defect silently, whether or not it has warned yet.

The measurement, its commands and its `sessionInfo()` output are the release
evidence for 2026-09-08; the numbers above are restated here so this record
stands on its own.

## Proposed outcome

A harness whose reported figure for a stage is that stage's steady-state cost,
resolved well enough that a 25 percent threshold means a 25 percent change in
the work done.

1. `median_time()` discards a warm-up before it measures. Two calls are enough
   on the stage measured here; the number is a property of the harness, not of
   the stage, and is stated in the script.
2. Repetitions are raised until the reported statistic is stable at the
   threshold that reads it. A stage costing 0.8 ms cannot be resolved to 25
   percent by three draws of a 1 ms clock.
3. The clock is finer than 1 ms. `Sys.time()` resolves to 9.54e-7 s on the
   measuring machine and is already available; whatever is chosen, the script
   records the tick it measured alongside the timings, so a future reader can
   see whether the resolution supported the comparison.
4. `tools/bench/baseline.csv` is re-recorded under the repaired harness, on a
   named machine and R version, and reviewed as a change to committed evidence
   rather than regenerated in passing. Until it is, no comparison against it
   means anything.
5. Only then does `--gate` go on and the Time gate become blocking, which is
   the state `dev/PROCESS.md` and `bench.yaml` both describe as pending.

Checkable: on the repaired harness, three consecutive runs of the whole table
on one machine agree within the gate's own allowance, and "online smoother, all
lags" reports a figure of the order of 0.8 ms rather than 5 ms.

## Affected code and users

`tools/bench/bench_reference.R` (the timing function, the repetition counts,
the recorded clock and the warning format at `:116`), `tools/bench/
baseline.csv` (re-recorded), `.github/workflows/bench.yaml` (`--gate`, last),
`dev/PROCESS.md` and `CONTRIBUTING.md` (the Time gate's status, once it
changes). No package source is touched and no user-visible behaviour changes:
`tools/` is excluded from the build and no test reads either file.

Carried into this package rather than taken in the release round: the warning
at `bench_reference.R:116` formats both ratios with `%.1f` and so renders this
stage as "ratio 0.1 vs 0.0", where 0.0 reads as zero. `%.4g` renders 0.05155
and 0.03947. It changes a warning string only.

## Constraints

- The baseline is committed evidence. It is re-recorded once, under review,
  with the machine, R version and commit recorded in the row, and never reset
  to clear a warning.
- The 25 percent allowance and the Section 8 budgets are not widened. The
  defect is in the measurement; widening the threshold would hide it.
- No package source and no graded number is touched by this work package.
- `--gate` goes on after the baseline is re-recorded, not before, and not in
  the same change.

## Open questions

- How many warm-up calls and repetitions the stages need. The 0.8 ms stage is
  the hardest case measured so far; the 0.66 s climate filter may need none.
  This is decided by measurement in the specification, not chosen here.
- Whether the baseline should be recorded on the CI runner rather than on an
  author's machine. The calibration loop was introduced to divide the runner's
  speed out, and whether it does so well enough to compare across machines is
  the question a re-recorded baseline has to answer.
- Whether stages whose cost is far inside their budget should be timed at all,
  or reported without a drift comparison.

## What this work package does not do

It does not change `tools/bench/` in the 2026-09-08 initial-release correction
round. Both files were left exactly as committed there, verified by sha256
before and after, and the correction round recorded the Time gate as advisory
with this measurement as the reason.
