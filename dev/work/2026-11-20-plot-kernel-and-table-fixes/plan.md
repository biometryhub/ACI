# Plan: plot arguments, kernel and table reuse (from intent.md 2026-11-20)

Status: draft 2026-11-20.

## Files that change

- `acir/R/aci-methods.R`: the three plot methods and their help for `...`.
- `acir/R/aci-kernels-scalar.R`: the scalar KL kernel, and the new internal
  helpers `.two_sum()`, `.two_product()` and `.exact_sum_sign()`.
- `acir/R/aci-assimilation.R`: `.validate_lag_table_source()`.
- `acir/R/aci-core.R`: two roxygen paragraphs wrapped at 80 characters.
- `acir/vignettes/vignette-2-advanced.Rmd`: one paragraph reworded to fit 80
  characters.
- `acir/man/`: the three plot pages, `aci_metric.Rd` and `aci_metric_pair.Rd`,
  regenerated from the roxygen.
- `acir/tests/testthat/`: new tests in test-24, test-42, test-43 and test-45,
  and the single-anchor test in test-42 also checks `sub` and `ylab`.
- `acir/NEWS.md`: the 0.2.1 section.

## Order of work

1. Write the tests and confirm that they fail on `17487dc`.
2. Change the plot methods, the kernel and the table check.
3. Regenerate the help, and keep only the five pages above.
4. Run the four gates.

## Risks

- A plot call that worked before PR #30 could draw differently. The methods
  are compared with `b0f7141` call by call.
- The kernel could change an ordinary result. Its results are compared with
  `17487dc` bit for bit, and the overflow decision with exact rational
  arithmetic.

## Proof

The new tests named above, the full suite with the MATLAB parity data,
`R CMD check --as-cran` and the fixture provenance check.

## Departures

Empty at approval.
