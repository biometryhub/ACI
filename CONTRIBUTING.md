# Contributing

`acir` is developed by its two authors under the process in `dev/PROCESS.md`.
Contributions from outside are welcome and pass through the same stages.

## Proposing a change

Open an issue, or write an `intent.md` from `dev/work/_template/` in a new
directory `dev/work/<date>-<slug>/` and open a pull request containing only
it. An accepted intent is one whose pull request merges.

## Building a change

Write `plan.md` before code: the files that change, the order of work, the
risks, and the proof that the change works. Commit it. Then make the change,
and if the implementation departs from the plan, update the plan in the same
commit.

## Before opening a pull request

Run the four gates in `dev/PROCESS.md` and quote their output in the pull
request under *How it was checked*. A change to the numerical core needs the
oracle tests, not only a green suite, because those are graded against a
source this repository did not author.

## Review

Every pull request is reviewed under `REVIEW.md` and needs one approving
review from an author who did not write it. Branch protection on
`acir-package` and `main` enforces this and the required checks.

## Reporting

Bugs and numerical disagreements have issue templates under
`.github/ISSUE_TEMPLATE/`. Security reports go through `SECURITY.md`.
