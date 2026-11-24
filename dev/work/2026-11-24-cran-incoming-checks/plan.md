# Plan: CRAN incoming checks for 0.2.2 (from intent.md 2026-11-24)

Status: draft, not yet approved.

## Files that change

- `acir/README.md`, changed: the lifecycle and licence badges link to the two
  files on GitHub by full URL.
- `acir/DESCRIPTION`, changed: `'MATLAB'` in single quotes in the Description,
  and CGNS spelt out as conditional Gaussian nonlinear systems where it first
  appears. `testthat (>= 3.1.5)` in Suggests, because the tests call
  `expect_no_error()`, `expect_no_warning()` and `expect_no_condition()`
  16 times and testthat's NEWS introduces them in 3.1.5.
- `acir/codemeta.json`, `acir/man/acir-package.Rd`, changed: the same
  Description text, kept identical to DESCRIPTION.
- `acir/tests/testthat/test-32-scalar-rows.R`, changed: `skip_on_cran()` at
  the top of the two blocks that build a 100,000-step record.
- `.github/workflows/R-CMD-check.yaml`, `.github/workflows/test-coverage.yaml`,
  changed: `NOT_CRAN: 'true'` in `env`, so both jobs run the two blocks.
- `acir/.Rbuildignore`, changed: `^cran-comments\.md$`.
- `acir/inst/WORDLIST`, changed: the 13 words the en-GB spelling check
  reports in the help, the vignettes, `NEWS.md` and `README.md` (`nats`,
  `Ornstein`, `Uhlenbeck`, `nondegenerate` and nine more), added with
  `spelling::update_wordlist()`. No text is changed.

## Order of work

1. Record the check time and the test time on the
   unchanged tarball.
2. Make the edits above.
3. Rebuild and repeat step 1.
4. Run the suite with `NOT_CRAN` unset and with `NOT_CRAN=true` and compare
   the counts.

## Risks

- Register row 29 (the cancellation of a differenced cumulative logarithm)
  rests on one of the two skipped blocks. On CRAN's machines that row is not
  exercised. It is exercised on all six CI configurations, in the coverage
  job and in any local run with `NOT_CRAN=true`.
- If `NOT_CRAN` were not set in the coverage job, the two blocks would be
  skipped there and coverage would fall. Step 4 shows the variable does what
  the workflows rely on.

## Proof

- Check time and test time before and after.
- Expectation counts: the full suite with `NOT_CRAN=true` is unchanged; with
  `NOT_CRAN` unset it is lower by the expectations of the two blocks, which
  are reported as skipped.
- The list of URIs in the built package contains no relative file link.
- `spelling::spell_check_package()` reports no word.

## Departures

None.
