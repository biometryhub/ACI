# Intent: correct the release candidate's claims, guards and process statements

Author: the `acir` authors. Status: draft 2026-09-08.

## Problem

The release candidate at commit `1e3d0df` was audited against its own claims.
The engine's numbers held: the graded surface reproduces its comparators, the
suite is green, and `R CMD check --as-cran` is clean. What did not hold is the
account the package gives of itself, in four kinds.

**Claims that outrun their evidence.** `acir/README.md:11-19`,
`API_STABILITY.md:5-10` and `inst/validation/README.md:3-6` say, in three
different wordings, that the graded quantities are verified against the method
authors' MATLAB implementation. The register itself does not say that: of its
65 rows, 28 are graded against a hash-pinned fixture and 10 of those against
the authors' own outputs; the rest compare against a source-derived run, an
independent transcription of the published equations, or an analytic identity.
The register is right and the three summaries are wrong. Two register rows and
two test descriptions carry the same overreach at a smaller scale: row 2
attributes the ACI signal and dispersion parts to a fixture that holds the
total alone, and `test-29` labels an independent transcription as
source-derived. The full-record partition summary is described as grading the
whole path when it is four order-insensitive reductions over 4,001 steps with
position fixed only at 201 sampled indices; a swap of two unsampled indices
moves the series by 0.0458 and leaves every sample and every reduction
identical.

**Mathematical wording that says more than the derivation supports.** The
Theorem 3 auxiliaries are described as exact for the discretization when what
is derived is a scheme-specific closed form; the causal influence range's
description of its quadrature does not match the finite sum the code computes
(`aci-core.R:799-802`, `aci-documentation.R:78-79`); the introductory vignette
states the causal reading of the metric more strongly than the assumptions
carry; and the tail and range descriptions in the register, the specification
and the roxygen do not agree with each other.

**Guards that are missing at the boundary.** The hidden-covariance integration
path solves per-slice Gram systems that the observation-Gram check upstream has
already passed, so an ill-conditioned block reaches `fill()` and returns a
number rather than an error. The filter and smoother kernels store means and
accumulate log-likelihoods without a post-loop finiteness test, so a diverged
run returns `NaN` in place of stopping. Regularization events are recorded in
`meta$regularization` and are silent at run time, so a user who never inspects
the metadata cannot tell a regularized result from a clean one. Channel widths,
`nsub` bounds and one internal token argument are unvalidated.

**Process statements that contradict the repository.** `dev/PROCESS.md:36-38`
lists Time as a gate every pull request carries; `bench.yaml` runs without
`--gate`, and two other statements in the same repository already say so.
`PROCESS.md:47-49` claims every listed local command exits non-zero on failure;
the lint command it lists exits zero. `PROCESS.md:74` records host branch
protection as verified when this repository holds no record of the host
configuration. The declared R floor is 4.0.0 and the only floor ever built is
4.1; `inst/CITATION` carries a transcribed version that DESCRIPTION has already
moved past. The one standing benchmark alarm is a measurement artefact: the
harness's three-call median on a 1 ms clock reports a cold call for a stage
that costs 0.83 ms warm.

None of this is a defect in the numbers. All of it is a defect in what a reader
is told, or in what the package does when it is handed something it cannot
compute.

## Proposed outcome

A release candidate whose every claim is one a reader can check from the
repository, and whose failure modes are refusals rather than quiet numbers.
Stated so it can be checked:

1. No sentence in the package says a quantity is verified against a source
   that its register row does not name. Every summary sentence either names
   the comparator or points at the register row that does.
2. Every documented mathematical property is the one the code computes, at the
   scope of the discretization or the quadrature it actually uses.
3. An ill-conditioned Gram on the hidden-covariance path, a non-finite mean and
   a non-finite log-likelihood are classed errors naming the quantity, the
   index and the time. A run that regularizes says so once, at run time.
4. Stated and enforced policy agree. A check described as blocking fails its
   job; a check that reports is described as reporting, with the reason.
5. The declared R floor is the floor that is built, and the installed citation
   reports the installed version.
6. The standing benchmark alarm has a recorded cause and a recorded repair,
   and neither the baseline nor the allowance is touched to clear it.

Everything the audit found is either fixed, or recorded as a named limitation
with what is missing. Nothing unrun is reported as a pass.

## Affected code and users

Users of the numerical core notice nothing: no graded number moves, no default
changes, no interface is added or removed. Two behavioural changes are visible
and are the point of the work.

- New error classes on paths that previously returned a number:
  `aci_error_gram_path` and `aci_error_nonfinite`. Code that relied on
  receiving `NaN` from a diverged run now receives a condition.
- A new warning class `aci_warn_regularized`, raised once per run when the
  regularization recorder fired. Code running under `options(warn = 2)` with a
  firing floor now errors where it previously returned. This is a
  compatibility-sensitive signal change and is presented for review as such.

Code touched: the conditioning and realiser paths (`aci-conditional.R`,
`aci-realiser-cache.R`), both kernels (`aci-kernels-scalar.R`,
`aci-kernels-matrix.R`), `aci-utils.R`, `aci-core.R`, `aci-assimilation.R`,
`aci-online-smoother.R`, `aci-cir.R`, `aci-cir-rows.R`, `aci-model.R`,
`aci-model-library.R`, `aci-documentation.R`; both vignettes; the evidence
register, `README.md`, `API_STABILITY.md`, `inst/validation/README.md`,
`inst/CITATION`, `DESCRIPTION`; four test files and one oracle manifest; the
fixture producer under `tools/`; `dev/PROCESS.md`, `CONTRIBUTING.md` and the
header comments of two workflows.

## Constraints

- No graded number moves. Every register fixture at its tolerance class, and
  the package's stored outputs on the reference records within 1e-12 of their
  values before the change.
- No fixture bytes, no tolerance, no `eps = 1e-12`, no prior default, no scheme
  label and no graded numerical value changes.
- No new export, no new public argument, no new S3 method, no new dependency.
- The scope rule holds: closed forms only.
- `acir/NEWS.md`'s 0.1.0 entry is history. `acir-v0.1.0` is a tag, so the text
  under it is left verbatim and corrected by a 0.1.0.9000 entry above it.
- `tools/bench/bench_reference.R` and `tools/bench/baseline.csv` are frozen for
  this round. The baseline is committed evidence and its re-recording is a
  separate decision.
- Existing authorship, scientific citations, licences and the source-port
  provenance statements are preserved exactly.

## Open questions

- Whether `aci_warn_regularized` should ship in this release at all, given that
  it can turn a previously returning call into an error under
  `options(warn = 2)`. The proposal is in the uncommitted diff for that
  decision to be made on.
- Whether to raise the declared R floor to 4.1.0, or to keep 4.0.0 and add a
  passing R 4.0 job to the check matrix. No R 4.0 build has ever run, so
  neither compatibility nor incompatibility at 4.0 is demonstrated.
- Whether an authorised run of the pinned fixture producer against a library
  holding `aci` 0.0.30 is available. Without it the producer repair is a
  static, unexecuted adapter and reproduction remains a documented limitation.
- When the benchmark harness is repaired and its baseline re-recorded, which
  decides when the Time gate becomes blocking.
- Whether the host's branch-protection settings will be queried before release
  and the "not verified" wording replaced with a captured result.
