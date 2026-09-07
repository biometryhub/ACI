# Specification: the performance work package

From `intent.md` (2026-09-01). Status: signed off 2026-09-01 by Aidan Moller
as part of the plan's approval.

## Requirements

1. Every Section 8 stage of the core-engine specification within its budget on
   the reference records, measured on the CI runner relative to a committed
   baseline.
2. No graded number moved beyond round-off: every register fixture at its
   tolerance class, and the stored outputs on the reference records within
   1e-12 of their values before the change.
3. No interface changed and no change to the scope rule.
4. Base R, no compiled inner loop.
5. Every ported kernel names its source file, function and tag.

## Design

Three mechanisms account for the shortfall, and each has a fix in code the
project already owns.

1. **The range walks the divergence table one cell at a time.** An R loop
   evaluates one relative entropy per (anchor, lag) cell, about 26
   microseconds over some four million cells. On the scalar path each
   anchor's whole row is one vectorised expression in the cumulative
   logarithms of the update factors, so the same work runs as vector
   arithmetic. The cumulative sums are blocked, restarted every 512 steps,
   so that differencing them never cancels. The adaptive truncation becomes a
   cut of the row by index, and the tail estimate is kept.
2. **The adaptive freeze costs more than it saves when nothing freezes.** On
   the dyad record it retained 1,499 of 1,500 lags and took 142 s against
   107 s for the full table. On the matrix path the bound recursions are
   skipped when nothing can freeze.
3. **The generic-model route realises coefficients on every call.** A model
   supplied as closures is evaluated at every grid point each time a verb
   runs; the library route caches its realisation and meets the budget. A
   cache on the generic route is keyed on the fingerprint the specification
   already defines.

Where logarithms do not apply, on matrix systems, all active anchors' update
matrices are advanced per step as one array operation; the measured spikes
were 2.9x to 7.7x for matrix systems and 27x for a thousand scalar records.

## Acceptance

The four gates every performance pull request carries: numbers (requirement
2), time (the benchmark against the baseline, warn at +25% until the gates
become failures), hygiene (`R CMD check --as-cran` at 0/0, coverage not below
the floor, the register updated where a check method changed), provenance
(requirement 5). The scalar-row change adds the worst case to its test set: a
long record, a short window, near-equal covariances.

## Flagged concerns, resolved at review

- The change of summation order moves bits: 79 exact-equality assertions
  across four test files pin route-to-route agreement. Restated as tolerance
  assertions inside the change, with the resummation declared in the
  pull request.
- Cumulative logarithms cancel on long records: measured on synthetic update
  factors in (0.90, 0.9995), differenced cumulative logs err by 1.8e-13 at
  N = 3,000, 1.6e-12 at N = 20,000 and 9.0e-12 at N = 100,000, so the climate
  record would breach the 1e-12 gate. Blocked sums restarted every 512 steps
  err by 3.6e-15. The blocked form is in scope.
- Two changes overlapped on the scalar path. The scalar cut is folded into the
  row change, and the matrix-path freeze fix moves ahead of the cache since
  it touches no numerical core.
