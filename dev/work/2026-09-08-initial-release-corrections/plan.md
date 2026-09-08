# Plan: the initial-release correction round (from intent.md 2026-09-08)

Status: draft 2026-09-08.

Base: `1e3d0df`, clean. The pre-change reference build is an installed copy of
that commit, kept aside so the 1e-12 comparison is against the tree as it was
and not against a re-run of the candidate.

## Order of work

Each step is measured by the one before it where that applies, and the shared
files fix the order more than the subject matter does. `acir/NEWS.md`,
`acir/inst/evidence/register.csv`, `acir/R/aci-utils.R`,
`acir/R/aci-kernels-scalar.R`, `acir/R/aci-kernels-matrix.R` and
`acir/vignettes/vignette-2-advanced.Rmd` each carry changes from several steps,
so each is brought to its final state once rather than in overlapping pieces.

| Step | Change | Closes |
|:---|:---|:---|
| 1 | Capture the pre-change state: install the untouched tree to a separate library, record the stored outputs on the reference records, the register fixtures at their tolerances, `sha256` of every fixture and of `tools/bench/baseline.csv`, and one strict benchmark run | the 1e-12 and no-byte-moved comparisons have something to compare against |
| 2 | Guards on the numerical core: the hidden-covariance Gram test and its error class; the four mean-store and the log-likelihood finiteness tests; `.aci_stop_cov()`'s two added sentences | ADV-01, ADV-02 |
| 3 | Boundary checks: channel widths, `nsub` bounds on both blocks, `.calc_tau()`'s token assertion, the `aci_metric()` and `aci_metric_pair()` help note | SBC-01, SBC-02, SBC-03, SBC-04 |
| 4 | Regularization visibility: the recorder latch, `.aci_reg_report()`, the eleven attachment sites, the four conditional print lines, the two existing tests restated to assert both classes | NUM-02 |
| 5 | The multi-hidden ENSO example and the diffuse-init rewording, in the roxygen, the advanced vignette and a pinning test | NUM-01, NUM-01b |
| 6 | Scientific wording: the Theorem 3 scope, the causal reading, the tail and range descriptions, the range's quadrature, and the `c(1, .5, 0)` identity test | MATH-01, MATH-02, MATH-03, RANGE-01, RANGE-02 |
| 7 | Evidence wording: the three parity summaries, register rows 2 and 66, the two `test-29` descriptions, `test-28` and the partition manifest, the probe-limitation and seed and observation notes | E01, E02, E05, ADV-03, RNG-01, OBS-01 |
| 8 | The fixture producer's static adapter and version check | E03 |
| 9 | Process and release metadata: the four-gate text, the local commands, the branch-protection row, the two workflow headers, the R floor, `inst/CITATION` | ENG-01, ENG-02 |
| 10 | The timing record: the measurement, its conclusion, and the harness work package | the benchmark alarm |
| 11 | `acir/NEWS.md` written once from the collected entries; `man/` regenerated with roxygen2 8.1.0 and the three tooling-only hunks reverted | the release record |
| 12 | Validation: the full suite, fixture provenance, `R CMD build` and `check --as-cran`, the 1e-12 comparison against step 1, the strict benchmark | acceptance |

Steps 2 to 9 are independent in subject and dependent in files. Step 2 edits
`aci-utils.R` before step 4 does; step 4 edits both kernels after step 2's
finiteness tests are in them; steps 4 and 5 share the advanced vignette; steps
6 and 7 share the register. Step 11 is last because `NEWS.md` collects what
steps 2 to 9 produce.

## Files that change

Package sources.

| File | Change | Item |
|:---|:---|:---|
| `acir/R/aci-conditional.R` | new internal `.check_gram_path()`; new internal `tgrid` argument on `.compiled_precision_path()`; the check called before `fill()`; three call sites pass `source_obs$t` | ADV-01 |
| `acir/R/aci-realiser-cache.R` | pass `obs$t` to `.compiled_precision_path()` | ADV-01 |
| `acir/R/aci-kernels-scalar.R` | the same `rcond` test on the constant dyad Gram; post-loop finite-mean guards in the filter and smoother kernels; post-loop finite-loglik guard; `nsub` bounds; freeze to report at three sites; the diffuse-init message | ADV-01, ADV-02, SBC-02, NUM-02, NUM-01b |
| `acir/R/aci-kernels-matrix.R` | post-loop finite-mean guards in both kernels; post-loop finite-loglik guard; `nsub` bounds; freeze to report at three sites; the diffuse-init message | ADV-02, SBC-02, NUM-02, NUM-01b |
| `acir/R/aci-utils.R` | `.aci_stop_nonfinite()`; `.aci_stop_cov()`'s two added sentences and the non-finite branch of `.cov_guard_scalar()`; the recorder latch and `.aci_reg_report()`; the `as_obs` and `observed_trajectory` notes. `safe_chol()` and `spd_floor()` unchanged | ADV-01, ADV-02, NUM-02, ADV-03, OBS-01 |
| `acir/R/aci-core.R` | `.calc_tau()`'s token assertion; the `aci_metric` and `aci_metric_pair` help note; freeze to report; two conditional print lines; the `aci()` description; the `aci_range` topic and its quadrature | SBC-03, SBC-04, NUM-02, MATH-02, MATH-03, RANGE-01 |
| `acir/R/aci-assimilation.R` | freeze to report; two conditional print lines; the `regularize` parameter blocks; the Theorem 3 scope wording | NUM-02, MATH-01 |
| `acir/R/aci-online-smoother.R` | freeze to report at three sites; the auxiliaries' scope wording; one comment | NUM-02, MATH-01, MATH-03 |
| `acir/R/aci-cir.R` | freeze to report; the tail and range description | NUM-02, MATH-03 |
| `acir/R/aci-cir-rows.R` | one comment on the row construction's scope | MATH-03 |
| `acir/R/aci-model.R` | channel-width validation inside `validate_cgns()`; the probe-limitation note and its example; the `seed` parameter | SBC-01, ADV-03, RNG-01 |
| `acir/R/aci-model-library.R` | the stated-prior multi-hidden example and one sentence naming the bare default | NUM-01 |
| `acir/R/aci-documentation.R` | the assimilation scheme section; the `causal_metrics` description, its tail wording and its quadrature | MATH-01, MATH-02, MATH-03, RANGE-01 |
| `acir/vignettes/vignette-1-intro.Rmd` | six sentence-level edits and one clause; the Theorem 3 sentence. Introductory prose: each diff presented before it is applied | MATH-01, MATH-02, NUM-01b |
| `acir/vignettes/vignette-2-advanced.Rmd` | the multi-hidden ENSO subsection; the regularization paragraph; the scheme and range prose; the observations section | NUM-01, NUM-02, MATH-01, RANGE-01, OBS-01 |
| `acir/README.md` | the three parity sentences | E01 |
| `acir/API_STABILITY.md` | the first bullet | E01 |
| `acir/inst/validation/README.md` | the opening paragraph | E01 |
| `acir/inst/evidence/register.csv` | row 2 free-text columns; row 36 feature and against text; row 66 feature text; one new behavioural row for the regularization contract. Rows 16, 33, 39, 44, 64 and 65 unchanged; no sha256, tolerance or column list moves | E02, MATH-03, ADV-03, NUM-02 |
| `acir/inst/CITATION` | `note` and `textVersion` built from `meta$Version`; the article bibentry untouched | ENG-02 |
| `acir/DESCRIPTION` | `Depends: R (>= 4.1.0)`. Nothing else | ENG-02 |
| `acir/NEWS.md` | one 0.1.0.9000 entry collecting every bullet below, including the breaking-changes bullet for the R floor. The 0.1.0 entry is under a tag and stays verbatim | all |
| `acir/man/*.Rd` | regenerated with roxygen2 8.1.0; the `RoxygenNote`, `NAMESPACE` layout and `acir-package.Rd` copyright hunks reverted | mechanical |

Tests and fixtures.

| File | Change | Item |
|:---|:---|:---|
| `acir/tests/testthat/test-04-cir.R` | one new block on `c(1, .5, 0)` at tolerance 0 | RANGE-02 |
| `acir/tests/testthat/test-15-compiled-scalar.R` | two class vectors gain `aci_warn_regularized`; one explicit assertion pins it | NUM-02 |
| `acir/tests/testthat/test-03-engine.R` | the stiff-Riccati expectation wrapped so both warnings are asserted | NUM-02 |
| `acir/tests/testthat/test-28-partition-oracles.R` | one comment sentence. No assertion, index set, reduction or tolerance | E05 |
| `acir/tests/testthat/test-29-tc-zeroth-order.R` | two `test_that()` description strings | E02 |
| `acir/tests/testthat/fixtures/oracles/oracle-manifest-partitions.yml` | the `exists_because`, `grades` and `does_not_grade` prose. No sha256, file name or column list | E05 |
| new test files | the new error and warning classes, the multi-hidden ENSO example pinned, the Gram and finiteness paths asserted by class and by message field | ADV-01, ADV-02, NUM-01, NUM-02 |
| every fixture | unchanged, verified by sha256 before and after | constraint |

Tooling, process and workflows.

| File | Change | Item |
|:---|:---|:---|
| `acir/tools/fixtures/make-partition-fixtures.R` | the producer version check, the verb table, `make_spec()`, `with_cond()` and the horizon adapter. `build_model()` unchanged | E03 |
| `dev/PROCESS.md` | the four-gate lead; Time recorded as advisory with the measured cause; Hygiene's four arms separated into blocking and reporting; the local lint command described as it behaves; the branch-protection row recorded as intended and unverified | ENG-01 |
| `CONTRIBUTING.md` | the two reporting gates named; branch protection recorded as intended | ENG-01 |
| `.github/workflows/bench.yaml` | header comment only. No step, flag or environment variable | ENG-01 |
| `.github/workflows/lint.yaml` | header comments only; `continue-on-error` and `LINTR_ERROR_ON_LINT` unchanged | ENG-01 |
| `dev/work/2026-09-08-bench-harness-resolution/intent.md` | new: the breach, its measurement, the proposed harness change and its baseline decision | the benchmark alarm |
| `dev/specification.md` | dated additive annotations on the tail and range statements; the existing lines verbatim | MATH-03 |
| `tools/bench/bench_reference.R`, `tools/bench/baseline.csv` | unchanged, frozen, verified by sha256 | constraint |
| `acir/CITATION.cff` | one comment line naming this file and `acir/codemeta.json` as 0.1.0 release records regenerated at release. No field changed | ENG-02 |
| `acir/codemeta.json` | unchanged: JSON carries no comment, so the note sits in `CITATION.cff` and here. Both are release records at 0.1.0 and are regenerated when the release version and the R floor are final | ENG-02 |
| `dev/work/2026-09-08-initial-release-corrections/quantity-to-comparator.md` | new: the quantity-to-comparator table behind the parity rewording, one entry per checked register row | E02 |

## Repository-only changes

Three changes in this round reach the repository but not the installed
package: `acir/tools/` is `.Rbuildignore`d, and `dev/` and `CONTRIBUTING.md`
sit outside `acir/`, so none of them is in the source tarball. They were
drafted as `NEWS.md` bullets and are recorded here instead, on 2026-09-08,
because `NEWS.md` is the changelog of the installed package.

- **The partition fixture producer is pinned to the build it declares.**
  `acir/tools/fixtures/make-partition-fixtures.R` resolves the pinned
  producer's own verb and argument names, so the script matches the `aci`
  0.0.30 build it declares, and refuses any other producer version. The
  adapter is static: the producer has not been run from here, and regeneration
  of the shipped fixtures has not been reproduced. The manifest's
  `regeneration_proof` carries the same dated qualification.
- **The development documents separate the checks that fail from the checks
  that report.** `dev/PROCESS.md` and `CONTRIBUTING.md` now say that the
  package check, the coverage floor and the fixture-provenance gate fail their
  jobs, while the timing comparison and lint report, each with its reason
  recorded in the workflow that runs it. Branch protection is described as the
  intended enforcement mechanism rather than a verified one, because this
  repository holds no record of the host configuration.
- **The benchmark's standing warning on the full-lag online smoother has a
  recorded cause.** The measurement, the cause and the proposed harness repair
  are in `dev/work/2026-09-08-bench-harness-resolution/intent.md`. Neither
  `tools/bench/baseline.csv` nor the +25 percent allowance was touched.

## Risks

- **The regularization warning breaks a caller.** Under `options(warn = 2)` a
  run that reaches a firing floor now errors. Mitigated by raising once per
  recorder rather than per event, by a classed condition a caller can suppress
  by class, and by leaving the decision to ship it explicit rather than silent.
  Two existing tests observe regularization warnings and are restated to assert
  both classes so the change is pinned, not tolerated.
- **The Gram test costs time on the two heaviest filter stages.** Measured at
  +12 to +15 percent against a +25 percent allowance, on audited ratios of 6.97
  against a baseline of 7.72. If the strict benchmark warns on either stage,
  the Gershgorin pre-screen replaces the `rcond` test and step 12 is re-run.
- **The Gram test rejects a matrix `fill()` would have solved.** The whole
  realised slice is tested at every interval start, which is stricter than the
  sub-block a masked route inverts. Accepted deliberately: the constructor's
  contract is on the whole Gram, so the whole-slice test is the contract the
  model was admitted under, and it makes a run's strictness independent of the
  route and of `options(aci.realiser_cache)`. The constructor's own `rcond`
  function and threshold are used, so the two cannot disagree about what is
  singular, and the new tests exercise the unmasked, masked, `"matlab"`
  first-slice and empty-target calls.
- **A post-loop finiteness test changes where a failure is reported.** A run
  that used to return `NaN` now stops, naming the first non-finite index rather
  than the step at which the divergence began. That is the trade taken for not
  paying a per-step test on the hottest loops.
- **The wording edits change meaning where none was intended.** The
  introductory vignette is human-written prose; its diff is presented sentence
  by sentence. The specification is annotated additively and its existing lines
  are left verbatim.
- **Several items share `NEWS.md`, the register and both kernels.** Each of
  those files is brought to its final state once, and the register's untouched
  rows are enumerated so a stray edit is visible in review.
- **The producer adapter cannot be run here.** It closes as a documented
  limitation. The risk is that it is later read as demonstrated reproduction;
  mitigated by saying so in the script, the intent and the NEWS entry.
- **The benchmark alarm reappears on every run.** It is expected, its cause is
  recorded, and neither the baseline nor the allowance is touched to clear it.

## Proof

1. **Numbers.** Every fixture in `acir/inst/evidence/register.csv` at its
   tolerance class. The package's stored outputs on the reference records
   compared with the pre-change installed build of `1e3d0df` at 1e-12,
   quantity by quantity, in two separate R processes so the two libraries are
   never loaded together. `sha256` of every fixture in
   `acir/tests/testthat/fixtures/` and `acir/inst/` identical before and after.
2. **Classes.** Each new condition asserted by class and by the fields of its
   message: `aci_error_gram_path` on the unconditioned and both conditional
   routes, under both cache settings, and on the compiled dyad's constant
   Gram; `aci_error_nonfinite` at each of the four
   mean-store sites and each likelihood accumulator; `aci_warn_regularized`
   once per recorder, with `test-15` and `test-03` restated to assert it
   alongside the classes they already assert.
3. **Examples and vignettes.** The multi-hidden ENSO example runs from the
   roxygen and from the advanced vignette and is pinned by a test; both
   vignettes build under `R CMD check`; the `c(1, .5, 0)` range identity holds
   at tolerance 0.
4. **Time.** `tools/bench/bench_reference.R --baseline tools/bench/baseline.csv
   --gate` on the candidate. The known "online smoother, all lags" artefact is
   the only permitted warning; a warning on `filter, library model` or
   `climate filter, 20000 steps` sends the Gram test back to the pre-screen.
   `sha256` of `baseline.csv` and `bench_reference.R` unchanged.
5. **Hygiene.** `R CMD build acir` and `R CMD check --as-cran` at zero errors
   and zero warnings, with examples and vignettes; `Rscript
   tools/oracle/check_fixture_provenance.R`; the full suite under
   `devtools::test("acir")`; coverage not below the floor.
6. **Provenance.** Every source-port statement, scientific citation, licence
   line and authorship line byte-identical to before.
   `utils::readCitationFile()` with `meta` supplied prints both bibentries and
   reports the DESCRIPTION version.
7. **Evidence.** The quantity-to-comparator table behind the parity rewording
   is `quantity-to-comparator.md` beside this plan, one entry per checked
   register row. The measurement behind the timing conclusion is restated in
   `dev/work/2026-09-08-bench-harness-resolution/intent.md`. The pre-change and
   post-change captures are release evidence for 2026-09-08 and are cited from
   the pull request rather than restated in it.

## Departures

Recorded in the same commit as the change that causes them. Empty at approval.

### 2026-09-08, numerical contract

- **ADV-01 test assertion.** The plan's `new_tests` bullet asks the test to
  assert that the covariance message does not contain `regularize = "floor"`.
  The plan's own proposed message text requires the message to say that
  substeps, the implicit stepper and `regularize = "floor"` do not repair the
  fault. The message specification was followed and the assertion written as
  the intent behind the test bullet: the message must not contain
  "raise nsub" and must not contain the recommending phrase "call with
  regularize" that `.aci_stop_cov()` uses.
- **ADV-01 `tgrid` argument.** `.compiled_precision_path()` takes `tgrid` with
  a `NULL` default and runs the realised-Gram check only when it is supplied,
  rather than unconditionally. `tests/testthat/test-25-enso6-batch-realiser.R`
  calls that function positionally on a hand-built array whose slice 7 is
  deliberately singular, to exercise the jitter-ladder refill; an
  unconditional check would refuse that array before the refill runs.
  All four compile-time call sites supply `tgrid`, so every public route is
  covered, and `tgrid` is placed after `first_step` so existing positional
  calls are untouched. The reason is stated in the function's roxygen.
- **ADV-02 helper location.** `.aci_stop_nonfinite()` is defined in
  `acir/R/aci-kernels-scalar.R`, not `acir/R/aci-utils.R` as the plan places
  it, so that it sits beside the guards that raise it. Both kernel files use
  it; the package namespace makes the location irrelevant to behaviour.
- **SBC-01 closure evaluations.** The diffusion evaluations and their shape
  tests were hoisted above `eval_coefs()`, but the evaluated matrices are not
  passed into `eval_coefs()`/`cgns_grams()` as the plan suggests. Leaving the
  count at what it was avoids changing how many times a user's coefficient
  closure is called during validation, which is observable and which no test
  pins. The reported error and its class are unaffected.
- **ADV-01 masked-route prediction.** The plan predicts that `method = "mask"`
  with a well-conditioned target block and a singular non-target block still
  succeeds. It does not, and the check does not mirror `fill()`'s block
  selection. Mirroring it made the refusal depend on where the precision path
  was realised: `.compile_cgns_complete()` obtains its arrays from
  `.realise_cgns_grid_cached()`, which realises the unconditioned whole-slice
  precision path eagerly, so with the cache in use the whole slice was tested
  and with `options(aci.realiser_cache = FALSE)` only the target sub-block
  was. The constructor's contract is on the whole Gram, and a principal
  sub-block of a well-conditioned symmetric positive-definite matrix is itself
  well-conditioned, so the whole slice is tested at every interval start on
  every route. The property is pinned at the unit level on
  `.compiled_precision_path()` and at route level under both
  `options(aci.realiser_cache = TRUE)` and `FALSE`.
- **New test file placement.** The SBC-01 and SBC-02 assertions are in
  `tests/testthat/test-37-input-contracts.R` rather than appended to
  `test-02-models.R`, `test-15-compiled-scalar.R` and
  `test-16-compiled-matrix.R` as the plan proposes.

### 2026-09-08, regularization policy and ENSO example

- **SBC-04 absolute-error figure.** The plan says to state that the absolute
  error stays below 1.3e-17 nats. Re-measured on the loaded tree, the maximum
  over the documented sweep is 2.88e-17 at `delta = 1e-07`. The help therefore
  says "below 3e-17 nats". The plan's figure would have been a false bound.
- **MATH-03 closing sentence.** The plan's draft reads "It is a diagnostic
  under the retained record, not a guarantee about observations beyond it."
  What the truncation drops is lag cells, not observations, so the sentence in
  `aci-core.R`, `aci-assimilation.R` and `aci-documentation.R` reads "not a
  guarantee about the cells the truncation dropped".
- **One helper beyond the plan.** `.aci_reg_cat(reg)` in `aci-utils.R`. The
  plan describes the same one-liner added to four print methods in two files;
  a single internal helper keeps them identical by construction. No export, no
  argument, and a no-op when `fired` is not `TRUE`.
- **Two documentation spans beyond the task detail.** The `aci_range()` topic
  description at the Simpson-default sentence (RANGE-01 item 4) and the
  tail-estimate sentence in the same block (MATH-03 item 5), both in
  `aci-core.R`. The verified-correct `count` versus `lag_time` span was left
  unchanged.
- **Test file name.** The plan names `test-36-regularization-visibility.R`;
  the file created is `tests/testthat/test-38-regularization-policy.R`,
  because `test-36` and `test-39` are in use by other work packages.

### 2026-09-08, scientific wording

- **Test file name.** The ENSO startup pins are in
  `tests/testthat/test-39-enso-example.R`, not the proposed
  `test-35-enso-startup.R`, because `test-35` is the observation-Gram file.
- **Refinement pins retained.** All four refinement pins are kept; the file
  runs in 1.21 s. The plan's stated-prior explicit `nsub` 1/2/4/20 sequence
  and its per-index Gram nondegeneracy sweep are shown in the `@examples`
  block and quoted in the help from the existing audit evidence rather than
  asserted, so no new evidence artefact was created.
- **Three plan edits applied outside the itemised list.** MATH-01 item (4),
  the second half of MATH-03 item (5), and the MATH-02 units and
  supplied-likelihood sentences, all in `aci-documentation.R`, because no
  other work package could reach them and they would otherwise appear nowhere
  in the help.
- **`\donttest{}` on the ENSO example checks.** The refinement,
  automatic-prior and prior-sensitivity runs measure 2.33 s together, above
  the ~1 s guidance; the unwrapped part is the stated configuration and the
  refusal only. The prior-sensitivity 8x refusal is stated in a comment rather
  than executed, to avoid a second `tryCatch` in the block.

### 2026-09-08, vignettes

- **`vignette-1-intro.Rmd` V5 grammar repair.** The plan places the nats and
  units sentence between the KL display and the sentence beginning "and it
  splits into". That sentence grammatically continued the sentence the display
  interrupted, so the insertion forces "and it splits into" to become "It
  splits into". This is the only change to existing human wording beyond the
  seven itemised edits, and is reversible by placing V5 after the
  decomposition sentence instead.
- **NUM-01b's vignette-1 clause not added.** A clause at
  `vignette-1-intro.Rmd:101-103` about a diffuse prior destabilising the
  explicit step was considered. This round's scope for that file is the six
  MATH-02 edits plus the MATH-01 sentence, and NUM-01b is the warning message
  alone, so the clause is out of scope. The statement is made in vignette-2's
  new ENSO subsection instead.
- **OBS-01's vignette-2 edit applied although not itemised.** RNG-01 and OBS-01
  are in scope as planned, and the reworded refusal message the chunk prints
  had already landed, so the chunk's expected output would otherwise be stale.
  Two hunks, revertible in isolation.
- **ADV-03 vignette demonstration.** None added; the item's documentation
  lands in the roxygen only.

### 2026-09-08, evidence

- **Dash typography in the E01 replacement text.** The plan's three
  replacement sentences use dash asides. The public documents deliberately
  carry none, following
  `tools/design/2026-08-17_register_and_typography_sweep.md`; the three asides
  are rendered as parentheses. Every other word is verbatim.
- **`README.md` paragraph break.** A break was inserted before "Maintainer:
  Aidan Moller." because the Status paragraph grew from three lines to eight
  under the replacement. Text unchanged.
- **`README.md` trailing sentence kept.** The plan's replacement stops at
  "announced in `NEWS.md`."; the sentence that follows in the file, "Details
  in `API_STABILITY.md`.", was kept verbatim.
- **E05 counterexample figure.** The plan's 0.045789292793430025 under a swap
  of one-based indices 2 and 3 did not reproduce. Measured on each partition's
  own 4001-step series, the move is 0.0395 (u), 0.0330 (hW), 0.0718 (tau) and
  0.0349 (tau_reduced3). The comment names the `u` partition and states
  0.0395, which reproduces from the shipped signal fixture.
- **`make-partition-fixtures.R` banner and `else` placement.** The banner
  prints `PRODUCER_VERSION` instead of calling `utils::packageVersion()` a
  second time, a consequence of the version gate; the printed output is
  unchanged. The plan's `VERB` assignment places `else` at the start of a
  line, which is a parse error at top level in R, and is written with braces.

### 2026-09-08, engineering

- **`dev/PROCESS.md` line 25.** The stages table's Deploy row carried the same
  unverified branch-protection claim corrected at line 74 and was corrected
  with it.
- **`lint.yaml` second comment.** The plan asked only for the comment at the
  old lines 36-40 to be extended. The comment at the old 41-43 asserted that a
  new diagnostic fails the job, which `continue-on-error` contradicts, and was
  corrected too. No step, flag or environment variable changed.
- **Inline lint comment shortened.** The plan's text put the line at 98
  characters; it reads "# reports (LINTR_ERROR_ON_LINT fails)" at exactly 80.
  The full form is in the prose immediately above the block.

### 2026-09-08, shared files and alignment

- **`.cov_guard_scalar()` message.** The plan names `.aci_stop_cov()` only. The
  same observation-Gram qualification was applied to the non-finite branch of
  `.cov_guard_scalar()` as well, because a user who reaches that branch has
  likewise already had the observation Gram accepted. The class, the fields and
  the remedies are unchanged.
- **Condition classes named on `aci_filter()` only.** The plan places the
  `aci_error_gram_path` and `aci_error_nonfinite` sentences on `aci_filter()`,
  and on `aci_smoother()`, `aci_online()`, `lag_table()` and `aci()` if their
  condition paragraphs list `aci_error_covariance_not_spd`. Those four
  delegate to `aci_filter()` rather than naming the class, so the sentences
  were added to `aci_filter()` only, at the end of the `@param regularize`
  paragraph rather than mid-paragraph, so that the "none" and "floor" policies
  are described without interruption.
- **`aci-documentation.R` tail sentence aligned.** Changed from the plan's
  literal "not a guarantee about observations beyond it" to "not a guarantee
  about the cells the truncation dropped", matching `aci-core.R` and
  `aci-assimilation.R`.
- **`test-29-tc-zeroth-order.R` labels.** All four `[source-derived]`
  `test_that()` descriptions became `[independent transcription]`, not only
  the two the plan itemises: the register rows behind the other two (47 and
  48) are also `independent_transcription`. The file's opening comment was
  aligned with them. Descriptions and one comment only; no assertion,
  tolerance or fixture touched.
- **Three register rows added, and the README counts updated.** Rows were added
  for the realised-observation-Gram contract, the non-finite recursion guards
  and the regularization-visibility contract. `README.md`
  said "65 rows: 28 ... and 37"; it now says "68 rows: 28 ... and 40", which
  is what `inst/evidence/register.csv` now holds. The ten
  authors-source-comparison rows are unchanged.
- **`test-33-scalar-aux.R` floor call wrapped.** The kernel call at line 133
  now fires `aci_warn_regularized` and surfaced as an unhandled warning in the
  suite. It is wrapped in `expect_warning(..., class = "aci_warn_regularized")`,
  which adds an assertion and weakens none.
- **`tools/bench/bench_reference.R` and the release-time metadata left
  frozen.** The `%.4g` warning-format change is deferred to
  `dev/work/2026-09-08-bench-harness-resolution/`, and `CITATION.cff` and
  `codemeta.json` are release-time regeneration items; neither was touched.
- **Register row 61 keeps its machine columns.** Its `against` text now says
  the check is ingestion and identity of the authors' pinned input record,
  but `tolerance_class = numerical_1e-6` and `state = checked` are unchanged
  because the register's machine columns are frozen in this work package.
  Reclassifying that row (an ingestion check has no numerical tolerance) is a
  release-time register decision, together with whether the predator-prey
  input record receives a row of its own.

### 2026-09-08, review corrections

- **Whole-Gram test on the `method = "reduce"` route.** The plan placed the
  test at the top of `.compiled_precision_path()` and left each route to pass
  the arrays it inverts. On the reduce route those arrays are the reduced
  ones, so a degeneracy on a dropped non-target channel was invisible to the
  test whenever `options(aci.realiser_cache = FALSE)` also removed the eager
  check in the realisation cache, and the record was accepted. The reduce
  branch of `.compile_cgns_complete()` now calls `.check_gram_path()` on the
  full realised Gram before the reduction and before the cross-block test,
  which is the order the cache produces when it is in use. The later test on
  the reduced arrays is kept: it is a principal sub-block of a slice that has
  already passed, so it cannot fail, and removing it would make the route
  depend on the caller for its own contract. `test-35` pins all six
  combinations of `NULL`, `"mask"` and `"reduce"` against both cache states.
- **The comparator table placed in the work package.** The round's evidence
  plan kept `quantity-to-comparator.md` outside the repository. It is the
  proof behind the parity rewording that Proof item 7 cites, and a proof a
  reader of this repository cannot reach is not a proof of anything they can
  check, so the table was placed beside this plan and added to the file table
  above.
- **`dev/specification.md` range annotations.** The plan promised dated
  additive annotations on the tail and range statements. The first round
  annotated the tail statement only. The forward-range capability row and the
  *Closed forms only* policy bullet now carry their own dated annotations,
  saying that `method = "exact"` is a finite threshold sum with no time-axis
  quadrature and that the `l1_linf` ratio uses composite Simpson by default,
  and a one-line dated pointer sits under the *What each source contributes*
  table so the tail annotation is reachable from both places the round names.
  Every existing line is verbatim.
- **Three `NEWS.md` bullets moved here.** The partition fixture producer, the
  development documents and the benchmark cause are changes to files that no
  installed package carries, so they were moved out of the changelog into
  *Repository-only changes* above.
