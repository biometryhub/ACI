# Review policy

This applies to every pull request into `acir-package` and `main`, and to
every reviewer, human or automated. A review is a list of findings ranked by
severity. Approval is a separate act, by an author who did not write the
change, and branch protection requires it before a merge.

## Three passes

1. **Numbers.** Does any graded quantity move? Check the evidence register
   (`acir/inst/evidence/register.csv`), the fixtures it names and their
   tolerance classes. A test, fixture or tolerance that changed is a finding
   whether or not the pull request declares it.
2. **Code.** Logic errors and edge cases: short records, one-dimensional
   against multi-dimensional states, the freeze paths, the generic-model
   route against the library route. The provenance of any ported kernel.
3. **Compliance.** Does the change match its `plan.md`, with departures
   recorded in the same commit? Is there a `NEWS.md` entry for anything a
   user would notice, and a register row for any new export?

## Blocking against nit

Blocking: a number outside its tolerance; a test weakened or removed; a gate
bypassed; a fixture changed without provenance; a change to the numerical
core without a plan. Everything about style, naming or wording is a nit.
Report at most five nits and give the rest as a count.

## Not reported

Anything CI already enforces (lint, check, coverage, bench), generated files
(`man/`, `NAMESPACE`), and the parents' code at the tag `parents-final`.
