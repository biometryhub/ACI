# Design record

Why the package is shaped the way it is. These documents are **historical**:
they record decisions and their reasoning at a point in time, and they are not
maintained against the code. Where they disagree with the package, the package
is right and these are a record of what was believed when they were written.

For what the package currently claims, read, in order of authority:

1. the code and its tests;
2. `aciR/inst/extdata/oracle-manifest.yml` — what each oracle grades, and what
   it does not;
3. `aciR/API_STABILITY.md` — what a user may rely on;
4. `aciR/NEWS.md` — what changed and why.

## Contents

| Document | What it is | Status |
|---|---|---|
| [`2026-07-15_critical_review.md`](2026-07-15_critical_review.md) | An external critical review of 0.0.0.9000: findings F1–F12, a package-standard scorecard, and a proposed release sequence. | Addressed. All P0 and P1 findings closed in 0.1.0. |
| [`2026-07-15_uplift_spec.md`](2026-07-15_uplift_spec.md) | The response: an adjudication of the review, twelve fixed design decisions (D1–D12), five work packages, a twenty-row test matrix, and a definition of done for 0.1.0. | Executed. Shipped as 0.1.0. |

## Why they are kept

Three reasons, in increasing order of importance.

**They answer "why is it like this?"** The spec fixes twelve decisions that the
code cannot explain on its own — why `aci_simulate_dyad()` was deleted rather
than deprecated, why the noise cross-covariance transpose is derived rather
than supplied, why a singular joint covariance is admissible while an
indefinite one is not, why the seed is contained. A maintainer who reopens one
of these should reopen it knowing what it cost to close.

**They record a disagreement that mattered.** The spec does not accept the
review verbatim. It amends it in six places and adds six findings the review
missed, and the amendments are argued rather than asserted. One of them was
load-bearing: the review proposed asserting that the causal-information metric
is exactly zero when the observed signal carries no information, which is true
in continuous time and false under the package's explicit Euler scheme. Written
as proposed, the test would have failed or been loosened into meaninglessness.
The version that shipped asserts the residual falls at second order in the step
— a stronger claim, and a true one. A reviewer's diagnosis can be right while
their prescription is wrong, and the record of that distinction is worth more
than either document alone.

**They contain the finding neither document started with.** The review asked
for the oracle's provenance to be written down. Writing it down forced the
question of what the oracle actually covers, and the answer was that the
flagship fixture pins the noise cross-covariance at zero and therefore never
grades the terms carrying it — while the package exposed those terms publicly.
The finding came from the act of stating the scope, not from either party's
list. That is the most useful thing in this directory and it is not in the
review.

## Reading them fairly

The review is critical by design and its verdict — "not ready for a public
scientific release" — was accurate for the version it examined. It was also
substantially correct: every finding whose evidence could be re-checked against
the source was confirmed. Neither document is an indictment or a defence. They
are the working-out.
