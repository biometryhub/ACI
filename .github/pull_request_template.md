<!--
Link the issue this closes, so that merging closes it:

  Fixes #<issue-number>
-->

## What this changes

<!-- One or two sentences. The why matters more than the what, which the diff
     already shows. -->

## How it was checked

<!-- Which of these ran, and what they said. A change to the numerical core
     needs more than a green test suite: it needs the oracle tests, because
     those are the ones graded against a source this repository did not
     author. -->

- [ ] `devtools::test()`
- [ ] `lintr::lint_package()`
- [ ] `R CMD check --as-cran`
- [ ] `Rscript tools/oracle/check_fixture_provenance.R`, if any fixture changed
- [ ] `NEWS.md` entry, for anything a user would notice

## Anything left open

<!-- A known gap stated here is a contribution. A known gap discovered later
     is a defect. -->
