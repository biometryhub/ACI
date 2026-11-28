# Plan: refuse ignored table settings and unknown prior entries (from intent.md 2026-11-24)

Status: draft 2026-11-24.

## Files that change

- `acir/R/aci-core.R`: `aci()` refuses `stepper`, `nsub` or `regularize`
  supplied with `table`, tested with `missing()`; the `@param table` and
  `@param init` text.
- `acir/R/aci-assimilation.R`: new internal `.check_init_names()`, called
  from `.same_gaussian_init()` on the supplied prior; the `@param init` text
  of `aci_smoother()`, `aci_online()` and `lag_table()`.
- `acir/R/aci-kernels-scalar.R`: `.scalar_filter_init()` calls
  `.check_init_names()` first.
- `acir/R/aci-kernels-matrix.R`: `.compiled_filter_init()` calls
  `.check_init_names()` first.
- `acir/R/aci-*.R`: any other `@param init` line, so each names the accepted
  entries.
- `acir/man/`: the pages regenerated from that roxygen.
- `acir/tests/testthat/test-49-ignored-inputs.R`: new.
- `acir/NEWS.md`: two bullets under `# acir 0.2.2`, `## Bug fixes`.

## Order of work

1. Write the tests and confirm that they fail on `ba3559b`.
2. Add `.check_init_names()` and its three calls.
3. Add the table refusal in `aci()`.
4. Update the help text, regenerate, and keep only the pages that changed
   for this work.
5. Run the full suite and the reproduction script on the branch.

## Risks

- A prior the package stores and passes back (`meta$init`, which carries
  `source`) could be refused. `source` is accepted, and the `aci_range()`
  rebuild and `keep = "table"` routes, which pass a stored prior, are in the
  existing suite.
- An existing call could pass one of the three arguments with a table. A
  search of `R/`, `tests/`, `vignettes/` and `man/` found none.
- A number could move. The check runs before any arithmetic and returns
  nothing; the suite's stored outputs and fixtures are unchanged.

## Proof

`test-49-ignored-inputs.R` fails on `ba3559b` and passes here; the full suite
passes unchanged.

## Departures

Empty at approval.

### 2026-11-24, help text

- Only the five exported `@param init` entries changed (`aci()`,
  `aci_filter()`, `aci_smoother()`, `aci_online()`, `lag_table()`). The
  internal `@noRd` entries describe arguments that receive a prior already
  checked by the verbs.

### 2026-11-24, `lag_table()` settings

- The help of `lag_table()` said `stepper` may be `"implicit"` and `nsub`
  any positive whole number. It refuses both with `aci_error_stepper`. The
  two entries now say so, because the new error from `aci(table = )` sends
  the reader to `lag_table()` for these settings.

### 2026-11-25, assignment into a model

- The open question in `intent.md` is settled by correcting the help of
  `aci_model()`. Assignment is still accepted, and the six existing tests
  that assign `fx` or `g` are unchanged.
- The help said that assignment leaves simulation and filtering on the
  original coefficients. Measured on closure models with one and two hidden
  components and on the library dyad: an assigned drift coefficient is used
  by `aci_filter()`, `aci_smoother()`, `aci_online()`, `lag_table()` and
  `aci()`, with results bit-identical to a model rebuilt with that
  coefficient, and `simulate()` keeps the original. An assigned diffusion
  coefficient is used by all of them.
- Three tests in `test-49-ignored-inputs.R` pin this. They pass on
  `ba3559b` as well, because no behaviour changed.
- A third NEWS bullet follows the two already planned.

### 2026-11-25, check time

- The two assignment tests on closure models take about three seconds
  each and are skipped on CRAN. The library-model test of the same
  behaviour still runs there, and all three run with `NOT_CRAN=true`.
