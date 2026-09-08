# Specification: the initial-release correction round

From `intent.md` (2026-09-08). Status: draft 2026-09-08.

## Requirements

Numbered, each checkable from the tree or from a named test.

1. **No unsupported parity claim.** No sentence in `acir/README.md`,
   `API_STABILITY.md` or `inst/validation/README.md` states that the graded
   surface is verified against the method authors' MATLAB implementation
   without naming the comparator classes the register records. The 0.1.0 NEWS
   entry stays verbatim and is corrected by a 0.1.0.9000 entry.
2. **Every register row says what its fixture holds.** Row 2 is narrowed to
   the total scalar ACI, in its free-text columns only; row 39 is retained as
   it stands; the two `test-29` descriptions name an independent transcription
   rather than a source-derived run. No sha256, tolerance, column list or
   assertion changes.
3. **The partition summary's coverage is stated as it is.** `test-28` and
   `oracle-manifest-partitions.yml` say the summary bounds magnitude over all
   4,001 steps and fixes position only at the 201 sampled indices. No sample,
   reduction, tolerance or fixture byte moves.
4. **The fixture producer resolves the producer it declares.** The verb and
   argument names of `aci` 0.0.30 are mapped explicitly, the producer's version
   is checked against the pinned one, and the `acir` branch of every mapping is
   the identity. The adapter is static; reproduction is not claimed.
5. **The Theorem 3 auxiliaries are described at their scope.** Both
   "exact for that discretization" warning strings become
   "derived for that discretization", and the roxygen, the specification's
   annotations and the vignette prose agree with each other and with the code.
6. **The causal reading is stated at the strength the assumptions carry.** Six
   sentence-level edits in `vignette-1-intro.Rmd` and the matching help
   wording; nothing else in that vignette.
7. **The influence range's quadrature is described correctly.** The exact form
   is the finite sum the code computes, not Simpson's rule.
   `aci-core.R:844-850` is verified correct and unchanged; a test on
   `c(1, .5, 0)` pins the identity at tolerance 0.
8. **An ill-conditioned hidden-covariance Gram is an error.** A per-slice
   `rcond(block) < 1e-12` test on exactly the blocks `fill()` selects, plus the
   same test on the constant Gram of the compiled dyad, raising
   `c("aci_error_gram_path", "aci_error_gram")` with the index, the time and
   the measured `rcond`.
9. **A non-finite result is an error.** Post-loop vector finiteness at all four
   mean-store sites and `is.finite()` on each likelihood accumulator, raising
   `aci_error_nonfinite` naming the quantity, the index and the time.
10. **A regularized run says so once.** A recorder that fired raises
    `aci_warn_regularized` (parent `aci_warning`) once per recorder, reporting
    the first event and pointing at `meta$regularization` for the total, at the
    eleven attachment sites only. The print methods of `aci_result`,
    `lag_table`, `cir_result` and `da_path_gaussian` carry one line, only when
    the recorder fired.
11. **The smaller boundary checks hold.** Channel widths, `nsub` bounds on both
    the scalar and the matrix blocks, and `.calc_tau()`'s token argument are
    validated; `aci_metric()` and `aci_metric_pair()` carry the help note.
12. **The multi-hidden ENSO example starts.** A stated-prior example on
    `aci_enso_model(hidden = c("u", "hW", "tau"))` runs from the roxygen and
    from the vignette, with refinement and prior sensitivity shown as checks
    and the bare constructor's default stated. The diffuse-init warning no
    longer implies that a burn-in window repairs a first-update failure.
13. **Stated and enforced policy agree.** Every check `dev/PROCESS.md` and
    `CONTRIBUTING.md` describe as blocking fails its job; Time and lint are
    described as reporting, with the reason; the local lint command is
    described as it behaves; host branch protection is recorded as intended
    and unverified. No workflow step, flag or environment variable changes.
14. **The benchmark alarm has a recorded cause.** The measurement is release
    evidence and the harness repair is a separate work package. Neither
    `tools/bench/bench_reference.R` nor `tools/bench/baseline.csv` changes,
    verified by sha256 before and after.
15. **The declared floor is the built floor and the citation is current.**
    `Depends: R (>= 4.1.0)`; `inst/CITATION` reads `meta$Version`; both
    bibentries still print. `CITATION.cff` and `codemeta.json` stay at 0.1.0
    as release records and are listed for regeneration at release.
16. **Nothing else moves.** Requirement set aside: no graded number beyond
    round-off, no fixture byte, no tolerance, no default, no export, no public
    argument, no S3 method, no dependency.

## Design

### Claims (requirements 1-4)

The register is the authority and the summaries are rewritten to defer to it.
Each of the three blanket sentences is replaced by one that names the
comparator classes and points at the row. The 0.1.0 NEWS text is under a tag
and is therefore history: it is corrected forward by a 0.1.0.9000 entry rather
than edited, which is the only form of correction available for released text.

The partition summary's gap is real and cannot be closed cheaply. There is no
full expected series in the fixture set: the four reference files are 201
sampled rows each and the summary is four order-insensitive reductions, so a
positional gate needs either new pinned bytes or a producer run. Both are out
of scope here, so the wording is corrected and the limitation stays named. The
alternative considered and rejected was extending the summary with an
order-sensitive statistic, which changes the fixture's sha256.

The producer adapter is written by reading the parent's `NAMESPACE` and sources
only. The parent is not executed. A static adapter does not establish
reproduction; only an authorised run of the pinned producer, compared with the
six shipped fixtures, does. That is stated in the script and in the intent.

### Wording (requirements 5-7)

Three separate errors with one shape: a property is asserted at a wider scope
than it was derived at. The auxiliaries are exact for the scheme they are
derived for, not for the continuous system; the range's ratio is a finite sum
over the grid, not a Simpson quadrature; the metric's causal reading depends on
the assimilation assumptions being met. Each is narrowed in place, in the
roxygen, the vignette and the register together, so the four statements of each
property agree. `dev/specification.md` is annotated additively and dated, its
existing lines left verbatim, because it is the record of what was specified.

`vignette-1-intro.Rmd` is human-written introductory prose. Its edits are
sentence-level, presented one at a time before they are applied, and nothing
else in that file changes.

### Guards (requirements 8-11)

The Gram check tests the whole realised Gram at every interval start, slices
1..N, whatever the route and whether or not the realisation cache is in use.
The constructor's contract is on the whole Gram, and a principal sub-block of a
well-conditioned symmetric positive-definite matrix is itself well-conditioned,
so the whole-slice test is sufficient for every block a route inverts and, being
independent of that selection, makes a run's strictness independent of
`options(aci.realiser_cache)`. It uses the constructor's own
`rcond` function and its threshold, so the two paths cannot disagree about what
is singular. It is placed at the top of `.compiled_precision_path()`, which
takes a new internal `tgrid` argument passed from all four call sites, so the
error can name the time as well as the index. `safe_chol()` is untouched:
the observation Gram is a different matrix with a different acceptance, and
`.aci_stop_cov()` gains two sentences saying so, since a user who reaches the
new error has already had that Gram accepted.

The alternative considered was a Gershgorin pre-screen, which is cheaper but
rejects a superset. It is held in reserve for the case where the strict
benchmark warns on the affected filter stages; the measured cost of the
`rcond` test on the two heaviest stages is +12 to +15 percent, and the climate
filter's audited ratio is 6.97 against a baseline of 7.72, so the candidate
sits inside the +25 percent allowance.

The finiteness tests are post-loop, not per-step: a per-step test on the
hottest loops would be paid on every step of every clean run for a condition
that ends the run once. Post-loop detection loses the exact step but names the
first non-finite index, which is what a user needs.

The regularization signal is a latch on the recorder, not a new argument and
not a change to any floor. `.aci_reg_report()` raises once and then delegates
to the existing `.aci_reg_freeze()`, substituted at the eleven attachment sites
only, so the guards themselves are untouched and no eps, clamp or cutoff moves.
Two alternatives were rejected: warning at the guard, which fires per event on
a stiff record, and printing unconditionally, which reports on clean runs.

### Process (requirements 13-15)

The unenforced Time gate is resolved by text, not by adding `--gate`. Adding it
would turn a job red permanently on a warning whose cause is the harness's own
statistic, and the only ways to clear it - resetting the baseline or widening
the allowance - are exactly the two that must not be taken. Lint is left
advisory for the same shape of reason: `lintr` cannot be run in the environment
where this decision is being made, so the tree's lint state is unknown, and a
gate should be declared blocking only when someone has seen it pass.

The R floor is raised rather than defended, because no R 4.0 build has ever run
and nothing in the package uses syntax newer than 4.0. Raising it makes the
declaration equal to the only floor that is built. Keeping 4.0.0 would require
adding an R 4.0 job and having it pass, which has not been done and might fail
for dependency-availability reasons that say nothing about `acir`.

## Acceptance

The four gates in `dev/PROCESS.md`, plus what is specific to this round.

- **Numbers.** Every fixture in `acir/inst/evidence/register.csv` at its
  tolerance class, and the package's stored outputs on the reference records
  within 1e-12 of their values before the change, compared against the
  untouched pre-change build rather than against a re-run of the candidate.
- **Time.** The strict benchmark on the candidate. The only permitted warning
  is the known "online smoother, all lags" measurement artefact; a warning on
  either of the two filter stages carrying the new Gram test is a result, not
  noise, and sends the Gershgorin pre-screen back into scope.
- **Hygiene.** `R CMD check --as-cran` at zero errors and zero warnings,
  including the vignettes and the examples; fixture provenance verified;
  coverage not below the floor.
- **Provenance.** Every source-port statement, citation, licence and
  authorship line is byte-identical to before the change.
- **Specific to this round.** Every new condition class is asserted by a named
  test, not merely tolerated. Every fixture sha256 in the tree is unchanged.
  `tools/bench/baseline.csv` and `bench_reference.R` are unchanged.

## Flagged concerns

- **`aci_warn_regularized` is a compatibility-sensitive signal change.** Code
  under `options(warn = 2)` that reaches a firing floor now errors where it
  previously returned a value. This cannot be made invisible and still be a
  signal. It is proposed in the uncommitted diff, with the two existing tests
  that observe regularization warnings updated to assert both classes, and is
  left for the authors' decision rather than assumed.
- **Two new error classes replace returned numbers.** `aci_error_gram_path`
  and `aci_error_nonfinite` stop runs that previously returned a number, which
  in the non-finite case was `NaN`. The judgement taken is that returning
  `NaN` from a diverged filter is worse than stopping; it is recorded here so
  the judgement is reviewable.
- **The new Gram test costs measurable time.** +12 to +15 percent on the two
  heaviest filter stages. Accepted against the +25 percent allowance on the
  audited ratios; if the strict benchmark warns, the cheaper pre-screen
  replaces it and the acceptance is re-run.
- **`NEWS.md` and `inst/evidence/register.csv` are each touched by several
  items in this round.** Their entries are collected and applied once, so each
  file has a single coherent state rather than a sequence of partial edits.
- **The producer repair is unexecuted.** Without an authorised run against a
  library holding `aci` 0.0.30, item 4 closes as a documented evidence
  limitation and not as a pass. That must be stated wherever regeneration is
  described.
- **Host branch protection is disputed.** The supplied feedback reports that
  the host has no branch protection; no host response and no remote CI record
  was captured here. Neither side of the fact is asserted; the text records
  the intent and says it is unverified.
