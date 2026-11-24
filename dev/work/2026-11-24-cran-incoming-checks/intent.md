# Intent: get the 0.2.2 tarball through CRAN's incoming checks

Author: Max Moldovan. Status: draft.

## Problem

`R CMD check --as-cran` on the 0.2.2 tarball is clean on macOS with R 4.6.1
(0 errors, 0 warnings, 0 notes, 109 s). Several things would still be raised
when the tarball reaches CRAN.

- `README.md` links to `API_STABILITY.md` and `LICENSE.md`. Both are in
  `.Rbuildignore`, so in the tarball the two links point at files that are not
  there, and CRAN reports them as invalid file URIs.
- The Description names MATLAB without the single quotes that *Writing R
  Extensions* asks for around the names of external software, and uses CGNS
  without spelling it out.
- `Suggests` allows testthat 3.0.0, but the tests use expectations that
  arrived in 3.1.5.
- Tests take 68 s of the 109 s. Two blocks of `test-32-scalar-rows.R`, both on
  a 100,000-step record, take 24 s between them. CRAN's Windows check has a
  10-minute limit and runs several times slower than a local check.

## Proposed outcome

The same tarball checks clean with CRAN's remote incoming checks switched on,
and the local check takes under 90 s. No number in the evidence register moves, and no test is removed.

## Affected code and users

`acir/README.md`, `acir/DESCRIPTION` with its two generated copies
(`codemeta.json`, `man/acir-package.Rd`),
`acir/tests/testthat/test-32-scalar-rows.R`, and the two workflows that run
the tests. Users see no change in behaviour.

## Constraints

No change under `acir/R/`. The two slow
test blocks still run on every platform in the CI matrix and in the coverage
job; only CRAN's machines skip them.

## Open questions

- The Title is 84 characters on two lines. *Writing R Extensions* asks for no
  continuation lines, and CRAN listings cut titles at 65 characters. Whether
  to shorten it is for both authors, since the article uses the same title.
- `cran-comments.md` is the maintainer's to write once win-builder and R-hub
  have run on the final commit. It is added to `.Rbuildignore` here.
