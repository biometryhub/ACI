# PR #29: warning provenance and numerical comparisons

## Scope and authority

The author approved completing the supervisor's four follow-ups on PR #29
(a679cbbb5f47868be8661a9eb9479a384a960e58, acir 0.2.1). Prepare changes,
tests and evidence without committing, pushing or posting. The usual
commit-before-code convention is represented by saved plans and source
identities because no commit is authorised. Preserve the current paper and
its package snapshot. No engine, API, dependency, default or fixture changes.

The original retained-table warning fix and its verification were completed
first. The author then approved the wording and matrix-test recommendations
and asked for a better round-off policy. This plan records the final scope.
The initial plan and prior evidence remain preserved outside this checkout.

## Implementation

- Check aci_result integration metadata before returning a retained table.
  The old missing-filter condition also excluded this route. Use the retained
  filter's settings for older table results lacking result metadata.
- Warn once per aci_range() call for non-default settings. The message must
  distinguish reuse of the already explicit table from recomputation. Leave
  the existing streamed-path warning in place without duplication.
- Strengthen the scalar grid's live, serialized and repeated-call counts.
  Share the existing correlated matrix fixture between tests 03 and 20 and
  extend that model to the complete public settings/retention grid.
- Require test 28's five expected complete finite series and exact matching
  lengths before numerical comparisons. Retain its strict maximum absolute
  1e-12 limit. This is a provisional regression gate, not a proved error bound.
- Enable optional per-series measurement output in CI and retain successful
  check artifacts. Include maximum gaps, scales, index/time, fixture identity
  and R/BLAS details. Do not fabricate unavailable Ubuntu measurements.
- Qualify the range help by common inputs/settings and numerical precision.
  Replace platform-independent bit-equality claims in the advanced vignette
  and rationale. Keep the historical measured values and their limited scope.
  Align NEWS and the two existing evidence descriptions.

## Verification

The strengthened warning regression failed on the original code. Reuse that
preserved failure and compare all 18 saved scalar/matrix numerical snapshots
with this final version. Run focused tests and the full suite, fixture
provenance, source/Rd consistency, package build and R CMD check --as-cran in
an isolated directory using existing dependencies. Temporary installation by
the standard build/check tools is confined to that directory. No dependency
installation, reference-code execution, new study or remote CI dispatch.

The matrix fixture initially captured the test runner's environment, causing
serialization warnings in source-loaded tests. Its factory now uses the acir
namespace as its enclosing environment. Equations and numerical results are
unchanged. Retain the initial failure log as evidence.

Cross-platform CI, lint and coverage remain the normal PR gates. Local tests
cannot certify an unrun Ubuntu job. Choosing a tighter tolerance requires
actual Ubuntu measurements and a stated error rationale, not another familiar
constant or automatic adjustment after a failure.
