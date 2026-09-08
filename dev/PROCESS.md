# How `acir` is developed

This document describes the development process of the repository: the
stages a change passes through, the artefact each stage commits, the gate
that closes it, and who decides. It follows the stage-and-artefact discipline
of the AI-native software development lifecycle playbook (Claxton, 2026),
adapted to a two-author scientific package whose deliverable is a numerical
result graded against the method authors' reference implementation.

Development uses AI coding assistance. The process below is what makes that
assistance auditable: what was asked is committed before what was built, what
was built is graded against a source its author did not write, and every
change is approved by an author who did not make it. `dev/acir-process-rationale.md`
records how the package was assembled before this process was written down;
this process applies from 2026-09-02.

## The stages

| Stage | Artefact | Where it lives | Gate | Who decides |
|:---|:---|:---|:---|:---|
| Plan | `intent.md`: problem, proposed outcome, affected code, constraints, open questions | `dev/work/<date>-<slug>/intent.md` | accepted by the other author when its pull request merges | the authors |
| Design | `spec.md`: requirements, design, flagged concerns | same directory; the package-level specification is `dev/specification.md` | flagged concerns resolved before engineering starts; the other author signs off | the authors |
| Build | `plan.md`: files that change, order of work, risks, proof | same directory, committed before code; a departure is recorded in the same commit as the change that causes it | a reader unfamiliar with the change could implement it from the plan alone | the author of the change |
| Test | the diff, its tests, the evidence register | `acir/tests/`, `acir/inst/evidence/register.csv`, `.github/workflows/` | the four gates below, locally before the pull request and again in CI | CI, deterministically |
| Deploy | review findings ranked by severity | the pull-request thread, under `REVIEW.md` | one approving review from an author who did not write the change; branch protection is the intended enforcement, see the correspondence table below | the reviewing author |
| Maintain | incident record | a new `intent.md` | a bench breach or a numerical disagreement re-enters at Plan | the authors |

## The four gates

Every pull request into `acir-package` carries these, and CI runs them again.
Numbers, the check and coverage arms of Hygiene, and fixture provenance fail
their jobs; Time and lint report and do not.

1. **Numbers.** Every fixture in `acir/inst/evidence/register.csv` at its
   tolerance class, and the package's stored outputs on the reference records
   within 1e-12 of their values before the change. A change that moves a
   number beyond round-off is a change of method and is reviewed as one.
2. **Time.** The Section 8 stages timed against the committed baseline
   (`tools/bench/baseline.csv`) on the CI runner, relative to the runner
   (`.github/workflows/bench.yaml`). Advisory for this release: the workflow
   runs without `--gate`, so drift is reported as a warning and does not fail
   the job. The reason is in the measurement, not in the policy. The harness
   times three calls and takes their median, on a clock whose tick is 1 ms
   here; the first two calls of a stage are warm-up, so for a stage costing
   less than a millisecond the median is a quantised cold call rather than
   the steady-state cost. The standing warning on "online smoother, all
   lags" is a two-tick difference between two such cold calls: measured warm
   on 2026-09-08 the stage costs 0.83 ms against a budget of 0.25 s. Time
   becomes blocking when the harness resolves the costs it measures; that
   change and its re-recorded baseline are
   `dev/work/2026-09-08-bench-harness-resolution/intent.md`.
3. **Hygiene.** `R CMD check --as-cran` at zero errors and zero warnings on
   the six configurations of `.github/workflows/R-CMD-check.yaml`, coverage
   not below the floor (`.github/workflows/test-coverage.yaml`, blocking),
   fixture provenance verified (`tools/oracle/check_fixture_provenance.R`,
   blocking), and lint reported (`.github/workflows/lint.yaml`, advisory
   until the style pass recorded in
   `tools/design/2026-08-17_air_formatter_decision.md`).
4. **Provenance.** Every ported kernel names its source file and function and
   the tag it was taken from, in the roxygen and in `acir/NEWS.md`.

## Verifying a change before opening a pull request

Run from the repository root, with the package installed from source. Every
command except the lint one exits non-zero on failure; `lintr::lint_package()`
prints its findings and exits zero unless `LINTR_ERROR_ON_LINT` is set, which
is the form CI uses. Each command's output is quoted in the pull request under
*How it was checked*.

```sh
Rscript -e 'devtools::test("acir")'
Rscript -e 'lintr::lint_package("acir")'  # reports (LINTR_ERROR_ON_LINT fails)
R CMD build acir && R CMD check --as-cran acir_*.tar.gz
Rscript tools/oracle/check_fixture_provenance.R
```

If a test fails, the fix is to the code. A change to a test, to a fixture or
to a tolerance is a change to the evidence, is listed in the pull request,
and is a review finding in its own right.

## Correspondence to the playbook

| Play in the playbook | Practice here | Status |
|:---|:---|:---|
| Capture as `intent.md` | `dev/work/<date>-<slug>/intent.md` | in use |
| Requirements and design as `spec.md` | `dev/work/<date>-<slug>/spec.md`; package-level `dev/specification.md` | in use |
| Plan before code, `plan.md` | `dev/work/<date>-<slug>/plan.md`, committed before the change | in use |
| Institutional knowledge as agent configuration (`CLAUDE.md`, skills, hooks) | not part of this repository; each author's tooling configuration stays with that author | not used in the repository |
| Parallel sessions in worktrees | one branch per work package, one writer per branch | in use |
| Feedback loop before reporting done | the four gates, run locally and in CI | in use |
| Continuous evaluation of agent configuration | not applicable to the repository | not used |
| Review loop under a written policy | `REVIEW.md`; approval by a human author is required | in use |
| Hooks as approval gates | intended branch protection on `acir-package` and `main`: pull request required, one approving review, the required checks green, no force push | intended; the host configuration is not recorded in this repository and has not been verified here |
| An agent inside CI/CD | not used; CI is deterministic | not used |
| Control bands in maintenance | `bench.yaml` compares every run with `baseline.csv`; a breach is triaged by the authors and re-enters as an intent | in use, with a manual response |
| Recurring security scans | not used; the package opens no network connection (`SECURITY.md`) | not used |
| On-call automation | not applicable | not used |

## What is measured

Read from the pull-request history when needed, not automated: the share of
pull requests whose checks pass on the first run, the time from a pull request
opening to its first review, and whether the merged diff matches the committed
plan. The bench workflow records the timing of every run as an artefact.

## First recorded instance

`dev/work/2026-09-01-performance/` is the performance work package: the
intent and specification were recorded inside its plan at the time and are
separated here, and the plan is the version approved on 2026-09-01 after
review. Later work packages start from `dev/work/_template/`.

## Reference

Claxton, L. (2026, August 21). *The AI-native SDLC playbook*. Anthropic.
https://claude.com/blog/the-ai-native-sdlc-playbook
