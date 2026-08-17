# aciR 0.1.0 revision and uplift specification

**Date:** 2026-07-15
**Baseline:** `aciR` 0.0.0.9000 at commit `254523c` (main)
**Inputs:** `aciR_package_critical_review_2026-07-15.md` (the review); full
independent re-inspection of `aciR/R`, `aciR/tests`, metadata, and the
`tools/oracle/` harness performed 2026-07-15.
**Target:** a defensible 0.1.0 research preview.

---

## 1. Verdict on the review

The review is high quality: evidence-based, file-and-line cited, correctly
prioritised (admissibility before UX before release engineering), and its
mathematics checks out. Every finding whose evidence could be re-verified by
independent code inspection was confirmed. Its strongest judgement -- keep the
MATLAB-oracle test as an immutable, never-skipping release gate -- is adopted
unchanged.

### 1.1 Verification status of the review's claims

| Claim | Re-verified? | Result |
|---|---|---|
| F1: `S_yoS_y = -1` and asymmetric `S_xoS_y` construct successfully | Yes, code trace of `aci_cgns_model()` (`aci-model.R:118-130`) | Confirmed: only `S_xoS_x > 0` is checked |
| F2: scalar-returning coefficient accepted; `NA` observations accepted | Yes, code trace of `.aci_as_coef()` and `aci()` | Confirmed: `is.function()` is the only coefficient check; `aci()` never checks finiteness of `x` |
| F3: no per-step covariance guards; `aci_metric()` domain unchecked | Yes, `aci-core.R:88-96, 148-157, 202-207` | Confirmed |
| F3: proposed `log1p` KL form is equivalent | Yes, algebra | Confirmed: `0.5*(-log(r) + r - 1) == 0.5*((r-1) - log1p(r-1))` exactly |
| F4: `aci_filter(1, ...)` reaches a replacement error | Yes, code trace | Confirmed: `for (j in 2:n)` with `n = 1` iterates `c(2, 1)`; at `j = 1`, `x[0]` yields `numeric(0)` and `m[1] <- numeric(0)` errors |
| F5: metadata gaps (URL, BugReports, cph, R floor, CI, `url: ~`) | Yes, read `DESCRIPTION`, `_pkgdown.yml`, tree | Confirmed; no `.github/`, no git remote |
| F6: duplicate dyad simulation, no `seed` in legacy path | Yes, read `dyad-model.R` | Confirmed |
| F7: `test-oracle-grade.R` duplicates the packaged oracle test behind a skip | Yes, read both test files | Confirmed; `_snaps/` is empty (no snapshot tests) |
| F9: README overstates ("validated against the authors' reference implementation") | Yes, README vs `tools/oracle/aci_oracle_dyad.m` header | Confirmed: the oracle is an independent harness reproducing the deterministic core of `marandmath/ACI_code`, not an unmodified upstream run |
| Fixture SHA-256 identity | Yes, re-hashed both copies | Confirmed, values match the review exactly |
| Check/coverage/lint numbers (0E/0W/2N, 99.57%, 23 lints) | Not re-run | Accepted on the review's recorded evidence; re-established by WP5 CI anyway |

### 1.2 Amendments to the review

Six places where the review is amended rather than adopted verbatim:

- **AM1 -- Use the pre-release breaking window.** The review proposes removing
  `S_xoS_y` "in the next permitted breaking release" and building lifecycle
  badges, deprecation warnings and a removal timetable for the legacy dyad
  helpers. The package is 0.0.0.9000 with no users and no remote: the permitted
  breaking window is *now*. Deprecation machinery for never-released symbols is
  process overhead. Consolidate by deletion before 0.1.0 (D2, D3).
- **AM2 -- Fix the loop idiom, not just the boundary.** F4's replacement errors
  are symptoms of the descending-range footgun (`for (j in 2:n)`,
  `for (j in (n - 1):1)`), which silently iterates backwards when `n < 2`.
  The spec mandates guarded lengths plus `seq_len()`-family loop bounds
  (r_style invariant 12) so the whole failure class is closed structurally.
- **AM3 -- Identity tests must be discretisation-aware.** F7.1
  (zero-information identity) holds exactly in continuous time but only to
  O(`dt`) under the Euler scheme: the filter integrates the prior ODE forward
  and the smoother re-integrates it backward, and forward-then-backward Euler
  do not invert exactly. Asserting `ACI == 0` would fail spuriously or force a
  meaninglessly loose tolerance. Frame these as convergence tests: error bound
  scaled to `dt`, and error shrinking ~linearly as `dt` halves. The terminal
  identity (F7.2) *is* exact by construction and can be asserted exactly.
- **AM4 -- Helper naming follows the package convention.** The review proposes
  `.check_observed_signal()` etc.; the package already uses an `.aci_` prefix
  (`.aci_check_scalar`). New helpers are `.aci_check_*` / `.aci_eval_*` (D6).
- **AM5 -- The singular-covariance boundary needs a policy, not a phrase.**
  F1's acceptance criterion "valid singular ... cases behave as documented" is
  under-specified. A singular joint covariance (`det == 0`) makes the
  conditional latent noise `B_j = 0`, which admits filter covariance decaying
  toward zero, which the smoother then divides by. Policy in D5: construction
  admits `det >= -tol` (mathematical admissibility), the per-step F3 guards own
  the runtime failure, and the singular boundary is documented and tested.
- **AM6 -- `plot.aci()` stays base-graphics, but the vignette shows ggplot2.**
  The review's base-graphics recommendation is right for the S3 method (zero
  new hard dependencies); ggplot2 is already in Suggests and the vignette
  should carry the publication-grade recipe.

### 1.3 Additions the review missed

Ordered by importance:

- **A1 (P0, scientific) -- The noise cross-covariance path is
  oracle-ungrounded.** The MATLAB harness sets `Sy_1 = 0; Sx_2 = 0`, so
  `S_yoS_x = 0` in the only fixture: the dyad oracle never exercises the
  `S_yoS_x != 0` branches of the filter (`aux <- comp$S_yoS_x + R * L_x`),
  smoother (`A_j`, `B_j`, and the `S_yoS_x * inv * (-dx + ...)` transport
  term). Yet `aci_cgns_model(S_yoS_x = 0.2, ...)` + `aci()` happily computes
  on that path today. Under the Independent Oracle Principle this path is
  `[unverified]`: live public code computing scientific output with no
  independent grounding. 0.1.0 must either ground it (WP3.3: analytic
  constant-coefficient Kalman-Bucy comparison with correlated noise, plus an
  optional second MATLAB fixture) or refuse it (error "correlated noise is not
  yet validated") until grounded. Grounding is the better path and is
  specified below.
- **A2 (P1, performance + contract) -- `aci_simulate()` draws `rnorm(1)` twice
  per step inside an R loop.** The increments are exogenous: draw
  `stats::rnorm(n - 1L)` twice before the loop and index. Order-of-magnitude
  constant-factor win at zero risk, but it changes the seed-to-path mapping --
  so it must land *before* the RNG contract freezes at 0.1.0 (same release as
  F11/D4, with a NEWS entry).
- **A3 (P2, canon) -- `Language: en-GB` should be `en-AU`** (r_style
  invariant 13).
- **A4 (P2, canon) -- Inline qualification.** `rnorm` is used bare with a
  scattered `importFrom(stats, rnorm)`; canon is `stats::rnorm()` at the call
  site (r_style invariant 10).
- **A5 (P2, project hygiene) -- No `PROJECT_LOG.md` and no L0 claim register**
  in the dev repo. The package asserts external-authority facts (paper
  citation/DOI, upstream repo URL and licence, MATLAB version): these belong in
  a light claim register; the oracle manifest (WP3.5) doubles as the grounding
  record for the numerical claims. `PROJECT_LOG.md` created 2026-07-15.
- **A6 (P3, docs wording) -- DESCRIPTION grammar**: "The method reimplements
  the approach of ..." should read "Reimplements the method of ...".

---

## 2. Design decisions

Each decision is fixed here so implementation does not re-litigate it.

**D1 -- Public surface at 0.1.0.**
Exports: `aci()`, `aci_cgns_model()`, `aci_dyad_model()`, `aci_simulate()`,
`aci_filter()`, `aci_smoother()`, `aci_metric()`, `aci_dyad_components()`,
plus S3 methods `print.aci_model`, `print.aci`, `summary.aci`,
`print.summary.aci`, `as.data.frame.aci`, `plot.aci`.
The components-level trio (`aci_filter`/`aci_smoother`/`aci_metric`) **is** the
supported expert surface (the docs already promise it), so it gets the full
validation contract (F4) and a versioned components schema documented on the
`aci_components` help page. `aci_dyad_components()` stays as the worked
example of that schema.

**D2 -- `aci_simulate_dyad()` is removed** (not wrapped, not deprecated).
Rationale: duplicate implementation with weaker validation and no `seed`
(F6), unreleased symbol (AM1), fully covered by
`aci_simulate(aci_dyad_model(...), ...)`. All examples and tests that use it
switch to the model path. NEWS records the removal.

**D3 -- `S_xoS_y` is removed as a constructor argument.** For the scalar
systems this package supports, the cross-Grammian is symmetric by
construction; the model stores `S_xoS_y = S_yoS_x` internally so the
components schema (and MATLAB naming parity) is unchanged. The components
validator requires the two entries to be equal (tolerance
`sqrt(.Machine$double.eps)`-scaled) for hand-built `comp` lists.

**D4 -- RNG contract: contained seed, `stats::simulate()` precedent, no new
dependency.** In `aci_simulate()`: if `seed` is `NULL`, consume the global RNG
stream as usual; if supplied, save `.Random.seed` (creating it first if
absent), `set.seed(seed)`, and restore the prior state via `on.exit()`.
Documented on the help page and tested (same seed twice -> identical path;
global RNG state unchanged after a seeded call). No `withr` import.

**D5 -- Admissibility boundary.** Constructor (and components validator)
require: `S_xoS_x > 0`; `S_yoS_y >= 0`;
`det = S_xoS_x * S_yoS_y - S_yoS_x^2 >= -tol` with
`tol = sqrt(.Machine$double.eps) * max(1, S_xoS_x, S_yoS_y)`. Singular systems
(`det` within `[−tol, tol]`, or `S_yoS_y = 0`) construct successfully -- they
are mathematically admissible -- but the documentation states that they can
drive the filter covariance toward zero, in which case the per-step guards
(D7) stop the run with a diagnostic. One test exercises this boundary.

**D6 -- Validator naming and placement.** New internal helpers, all in a new
`R/aci-validate.R`, all `@noRd` + `@keywords internal`, all erroring with
`stop(..., call. = FALSE)` and backticked argument names (no `cli`, no
`checkmate`):
`.aci_check_signal(x)`, `.aci_eval_coef(fn, name, x)`,
`.aci_check_noise_covariance(S_xoS_x, S_yoS_y, S_yoS_x)`,
`.aci_check_components(comp, n)`, `.aci_check_posterior(post, name, n)`.
`.aci_check_scalar()` moves in beside them.

**D7 -- Per-step guards fail fast, never clip.** After each covariance update
in filter and smoother: `if (!is.finite(R) || R <= 0)` stop with algorithm
name, step index, time `(j - 1) * dt` (filter) or `(j - 1) * dt` at index `j`
(smoother), offending value, and the guidance "reduce `dt` or check the
model's noise covariance". Inline comparisons (not a helper call) to keep the
30k-step loop cheap. No silent clipping of substantive negatives anywhere; the
only clamp in the package is the metric's round-off clamp (D8).

**D8 -- Metric hardening.** `aci_metric()` validates both posteriors
(strictly positive covariances, equal lengths, finite means), computes the
dispersion with the cancellation-resistant form
`ratio_delta <- smooth$cov / filt$cov - 1;`
`dispersion <- 0.5 * (ratio_delta - log1p(ratio_delta))`,
clamps values in `(-1e-10, 0)` to zero (documented round-off tolerance,
matching the existing test bound), and errors on anything `<= -1e-10`
(defensive: cannot occur with valid inputs).

**D9 -- Optional `time` argument on `aci()`.** `aci(x, model, dt, ...,
time = NULL)`: when `time` is supplied it must be numeric, the length of `x`,
strictly increasing and equally spaced (spacing agreeing within
`sqrt(.Machine$double.eps)`-scaled tolerance); `dt` is then derived and must
not also be supplied inconsistently. Irregular grids are rejected with a
message stating the algorithm's regular-sampling requirement. The regular,
complete-observation contract moves into the `aci()` docs and the assumptions
vignette prominently (F11).

**D10 -- Result methods stay thin and pedagogical** (F10):
`summary.aci()` returns a `summary.aci` object (not printed text) holding: n,
`dt`, time span, metric five-number summary + peak value/time, minimum filter
and smoother covariance (stability indicators), terminal-identity residual,
and the count of round-off-clamped metric values. `as.data.frame.aci()`
returns columns `t, x, filter_mean, filter_cov, smoother_mean, smoother_cov,
aci`. `plot.aci()` uses base graphics: two stacked panels (observed signal;
ACI trace), colourblind-safe defaults, `...` passed through. ggplot2 recipe
lives in the vignette.

**D11 -- Oracle gate is immutable.** `test-oracle-dyad.R` is untouched except
for additions; the `1e-6` maximum-absolute-error gate must pass after every
WP. `test-oracle-grade.R` (workspace-dependent duplicate) is deleted and
replaced by a fixture SHA-256 identity test against the manifest (WP3.5).

**D12 -- Style tooling.** `.lintr` at package root keeps default linters but
configures `object_name_linter` with a narrow additional regex admitting the
documented paper-aligned notation (`L_x`, `f_y`, `S_xoS_x`, `R0`, `mu0`,
`A_j`, `B_j`, `muT`, `RT`, ...), with a comment naming the source equations.
No whole-file exclusions. `inst/WORDLIST` covers ACI, CGNS, Andreou, Bollt,
Chen, Kullback, Leibler, Grammian(s), dyad, smoother, SDE(s), etc. New
unexplained lint output fails CI.

---

## 3. Work packages

Dependency order: WP1 -> WP2 -> WP3 -> WP4 -> WP5. WP2 and WP3 can overlap
after WP1 lands. Every WP ends with `devtools::test()` green including the
untouched oracle gate, run from a vanilla session.

### WP1 -- Correctness boundary (release blocker; review Phase 1 + A1 gate, AM2)

Files: `R/aci-validate.R` (new), `R/aci-core.R`, `R/aci-model.R`,
`R/dyad-model.R`, tests.

1. Implement the five validators of D6 with the error catalogue of §6.
2. `aci_cgns_model()`: call `.aci_check_noise_covariance()`; validate `label`
   (single non-NA character) and `parameters` (NULL or named list of finite
   numeric scalars); remove `S_xoS_y` argument (D3).
3. `aci()`: `.aci_check_signal(x)`; evaluate each coefficient through
   `.aci_eval_coef()` (rejects wrong type, wrong length -- no implicit
   recycling -- and non-finite values, naming the coefficient and the received
   type/length); assemble and `.aci_check_components()`.
4. `aci_filter()` / `aci_smoother()`: `.aci_check_signal()`,
   `.aci_check_components(comp, n)`, scalar checks on `dt > 0`, `mu0`,
   `R0 > 0`; `.aci_check_posterior(filt, "filt", n)` in the smoother; loop
   bounds rewritten via `seq_len()` arithmetic so no descending range can
   occur (AM2); per-step guards per D7.
5. `aci_metric()`: D8.
6. `aci_simulate()`: vectorise increment draws (A2); inline
   `stats::rnorm()` qualification (A4).
7. Tests: covariance admissibility (negative `S_yoS_y`, negative determinant,
   singular boundary per D5); coefficient contract (scalar return, wrong
   length, wrong type, `NA`/`Inf` return, non-numeric); signal contract
   (`NA`, `Inf`, matrix, length 1); per-step guard triggering (a valid model
   with `dt` large enough that Euler loses positivity -- e.g. the dyad model
   with `dt = 1`); posterior contract at `aci_metric()` and `aci_smoother()`;
   error snapshots for the public messages (populates `_snaps/`).

**Gate:** oracle max-abs error < 1e-6 unchanged; every malformed input in the
test matrix fails before recursion starts; `aci_filter(1, ...)` and friends
fail with the boundary message, not a replacement error.

### WP2 -- API consolidation (review Phase 2, AM1)

1. Delete `aci_simulate_dyad()`; migrate its uses in examples
   (`aci_filter`, `aci_smoother`, `aci_metric`, `aci_dyad_components`) and
   `test-aci-core.R` to `aci_simulate(aci_dyad_model(...), ...)` (D2).
2. `aci_dyad_components()` gains parameter validation (named list, the seven
   scalars, each finite; `sigma_x != 0`) and keeps its schema-example role.
3. Write `API_STABILITY.md`: lifecycle stage per export (model layer:
   maturing; components layer: maturing-expert; everything documented
   "no breaking changes without a minor-version NEWS entry" post-0.1.0),
   components-schema version statement, and the policy that the oracle gate
   is a release invariant.
4. NEWS.md restructure per r_style invariant 13 (`## Breaking changes` with
   before/after blocks for D2/D3, `## New features`, `## Minor improvements
   and fixes`) under a `# aciR 0.1.0` heading written at integration time
   (version bump deferred to the release commit).

**Gate:** exactly one implementation per operation; `R CMD check` clean; all
examples run on the consolidated surface.

### WP3 -- Scientific validation and documentation (review Phases 2-3, A1, AM3)

1. **Zero-information convergence test** (F7.1 + AM3): `L_x = 0`,
   `S_yoS_x = 0`; assert `max(aci)` below a `dt`-scaled bound and that the
   bound tightens ~linearly when `dt` halves.
2. **Terminal identity test** (F7.2): exact equality at index `n`; terminal
   metric exactly 0.
3. **Correlated-noise grounding** (A1, P0): constant-coefficient linear CGNS
   with `S_yoS_x != 0` compared against the analytic stationary Kalman-Bucy
   solution (Riccati fixed point for the filter; the analogous backward fixed
   point for the smoother) -- convergence in `dt`, plus the exact stationary
   values at tight tolerance for small `dt`. Optional (Max decision O2): a
   second MATLAB fixture from a correlated-noise variant of the harness, same
   manifest treatment. Until (3) lands, the cross-noise path is documented
   `[unverified]`; with (3), the README/vignette may claim the general core is
   validated on both independent and correlated noise.
4. **Discretisation stress test**: dyad model, `dt` sequence, document the
   convergence order; plus the D7 failure test at coarse `dt`.
5. **Oracle provenance manifest** (F9): `inst/extdata/oracle-manifest.yml`
   per the template in §7; a test re-hashes the shipped fixtures against the
   manifest (replaces `test-oracle-grade.R`, D11); a maintainer refresh script
   (`tools/oracle/refresh_fixtures.R` guidance) that regenerates to a temp dir,
   diffs hashes and requires explicit copy.
6. **Assumptions and interpretation vignette** (F8): estimand and
   directionality; what is conditioned on; model-correctness, known-parameter
   and identifiability assumptions; regular complete sampling and fixed `dt`;
   scalar-only scope; sensitivity to `mu0`/`R0`/`dt`; ACI evidence under the
   model vs causal identification from observational data; what an ACI peak
   does and does not support; recommended sensitivity reporting.
7. **Validation vignette**: oracle design (harness, seed, checksums), fixture
   hashes, maximum observed errors for the release, the correlated-noise
   grounding, and limitations. README rewording (F9): "validated against an
   independent MATLAB harness that reproduces the deterministic core of the
   authors' reference implementation
   (`github.com/marandmath/ACI_code`, MIT)".

**Gate:** every scientific assumption and every public input class maps to at
least one test (traceability table in the validation vignette); no
`[unverified]` external claim remains in shipped prose.

### WP4 -- Result UX and contracts (review Phase 3)

1. `summary.aci()` + `print.summary.aci()`, `as.data.frame.aci()`,
   `plot.aci()` per D10, each with runnable examples and tests (including
   snapshot of the summary print).
2. `time` argument on `aci()` per D9, with irregular-grid rejection tests.
3. RNG contract per D4, with state-restoration test; NEWS entry for the
   changed seed-to-path mapping (A2).
4. Benchmark script (dev repo `bench/bench_aci.R`, Rbuildignored): 30,001 and
   3e5 steps, simulate + fit timings, recorded in `PROJECT_LOG.md`; decides
   whether any further optimisation is warranted (expected: no).

**Gate:** inspect/plot/export workflows need no internal-list indexing; `aci`
object layout unchanged.

### WP5 -- Release engineering (review Phase 4; blocked on decision O1)

1. DESCRIPTION: `URL`, `BugReports`, `cph` role, `Depends: R (>= 4.1.0)`
   (provisional -- confirmed by CI oldrel), `Language: en-AU` (A3),
   Description rewording (A6).
2. CI: r-lib/actions -- check matrix (Linux devel/release/oldrel, macOS
   release, Windows release), coverage (covr -> Codecov), pkgdown build+deploy,
   lint job honouring `.lintr` (D12).
3. `_pkgdown.yml` url; README install instructions for the chosen home;
   `CITATION.cff` (cffr) and `codemeta.json` (codemetar), both hand-reviewed;
   `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `cran-comments.md`.
4. `.lintr` + `inst/WORDLIST` (D12); spelling check wired into CI.
5. Release checklist: vanilla-session `R CMD build` + `R CMD check --as-cran`
   including PDF manual; tarball hygiene review; sensitive-info grep guard
   (r_style hard rule) over `R/ tests/ vignettes/ inst/`.
6. Version bump to 0.1.0 and the NEWS heading land in the integration commit
   only (parallel-reconciliation invariant), tagged on green CI confirmed via
   `gh run view <id> --json conclusion`.

**Gate:** green CI on the declared matrix; URLs/version/attribution agree
across DESCRIPTION, README, CITATION, CFF, codemeta, pkgdown; check clean
apart from justified incoming notes recorded in `cran-comments.md`.

---

## 4. Test plan matrix (new tests only; oracle tests unchanged)

| # | Contract | Case | Expectation | File |
|---|---|---|---|---|
| T1 | Covariance admissibility | `S_yoS_y = -1` | E03 error at construction | test-aci-model.R |
| T2 | Covariance admissibility | `S_yoS_x^2 > S_xoS_x * S_yoS_y` | E04 error | test-aci-model.R |
| T3 | Singular boundary | `S_yoS_y = 0`, `S_yoS_x = 0` | constructs; run either completes or stops via E15 | test-aci-model.R |
| T4 | Label/parameters | `label = NA`, unnamed `parameters` | E05/E06 | test-aci-model.R |
| T5 | Coefficient return | `function(x) 1` on length-300 signal | E08 names coefficient, expected/received length | test-aci-model.R |
| T6 | Coefficient return | function returning character / `NaN` | E08/E09 | test-aci-model.R |
| T7 | Signal | `NA`, `Inf`, matrix, length-1 into `aci()` and `aci_filter()` | E10-E12 before recursion | test-validate.R |
| T8 | Components schema | missing entry, wrong length, `S_xoS_y != S_yoS_x` | E13 | test-validate.R |
| T9 | Posterior | mismatched lengths, zero/negative cov into `aci_metric()`/`aci_smoother()` | E14 | test-validate.R |
| T10 | Per-step guard | dyad model, `dt = 1` | E15 with algorithm, index, time | test-aci-core.R |
| T11 | Zero-information | `L_x = 0`, `S_yoS_x = 0`, two `dt` values | `max(aci)` under dt-scaled bound; shrinks with `dt` | test-identities.R |
| T12 | Terminal identity | any valid fit | exact filter/smoother equality and zero metric at `n` | test-identities.R |
| T13 | Correlated noise (A1) | constant-coefficient CGNS, `S_yoS_x != 0` | matches analytic stationary Kalman-Bucy; converges in `dt` | test-identities.R |
| T14 | Convergence order | dyad, `dt` halving | error ratio ~2 (first-order) | test-identities.R |
| T15 | Metric near ratio 1 | `smooth$cov == filt$cov` exactly | dispersion exactly 0, no negatives | test-aci-core.R |
| T16 | RNG contract | seeded call | identical repeat; global `.Random.seed` restored | test-aci-model.R |
| T17 | Time grid | irregular `time`; inconsistent `dt`+`time` | E18 | test-aci-model.R |
| T18 | Error snapshots | all public messages | snapshot-locked | test-errors.R |
| T19 | Fixture identity | shipped CSVs vs manifest hashes | SHA-256 equal | test-oracle-manifest.R |
| T20 | Methods | summary/as.data.frame/plot | shapes, snapshot of summary print, plot returns invisibly | test-methods.R |

---

## 5. Definition of done for 0.1.0 (amended from the review)

The review's DoD is adopted with three additions (marked +):

- All P0/P1 findings closed or documented with maintainer-approved rationale.
- No invalid covariance, coefficient, signal, components or posterior contract
  can reach numerical recursion.
- Oracle maximum absolute error < 1e-6, unchanged gate.
- Zero-information, terminal-identity, invalid-domain and discretisation tests
  pass on all supported R/OS combinations.
- **+** The correlated-noise path is independently grounded (T13) or made
  unreachable -- no `[unverified]` public computation path.
- **+** `aci_simulate()` seed contract is contained (D4) and the vectorised
  path change is NEWS-recorded (A2).
- Every export has a declared lifecycle stage (`API_STABILITY.md`).
- Assumptions and interpretation limits visible from README and pkgdown.
- `R CMD check --as-cran` clean (PDF manual built) apart from justified
  incoming notes; run in CI and a local vanilla session.
- URLs, version and attribution agree across all metadata surfaces.
- **+** Sensitive-info grep guard clean; `Language: en-AU`; lint and spelling
  jobs quiet under the committed `.lintr`/WORDLIST policy.

---

## 6. Error-message catalogue (locked by T18 snapshots)

Final wording set at implementation; identifiers and required content fixed
here. All errors use `stop(..., call. = FALSE)` with backticked arguments.

| ID | Trigger | Required content |
|---|---|---|
| E01 | non-scalar/non-finite scalar arg | arg name; "single finite numeric value" (existing wording kept) |
| E02 | `S_xoS_x <= 0` | existing wording kept |
| E03 | `S_yoS_y < 0` | arg name; "non-negative"; role (latent-noise covariance) |
| E04 | determinant `< -tol` | the computed determinant; "not positive semidefinite" |
| E05 | bad `label` | "single non-missing character string" |
| E06 | bad `parameters` | "NULL or a named list of finite numeric scalars" |
| E07 | coefficient wrong type at construction | existing wording kept |
| E08 | coefficient wrong return type/length | coefficient name; expected length; received type and length; "no recycling" |
| E09 | coefficient non-finite return | coefficient name; first offending index |
| E10 | `x` wrong type/class | "plain numeric vector" |
| E11 | `length(x) < 2` | existing wording kept |
| E12 | non-finite `x` | offending value class (`NA`/`NaN`/`Inf`) and first index |
| E13 | malformed `comp` | missing/misshapen entry names; expected schema pointer to `?aci_components` |
| E14 | malformed posterior | argument name; equal-length and strictly-positive-cov requirements; first offending index |
| E15 | per-step covariance failure | algorithm (`filter`/`smoother`); step index; time; value; "reduce `dt` or check the model's noise covariance" |
| E16 | metric `<= -1e-10` | internal-inconsistency wording; ask to report a bug |
| E17 | bad `n`/`dt`/`R0` | existing wordings kept |
| E18 | bad `time` | "strictly increasing and equally spaced"; regular-sampling requirement |

## 7. Oracle manifest template (`inst/extdata/oracle-manifest.yml`)

```yaml
oracle: dyad_interaction_model deterministic core
upstream:
  repository: https://github.com/marandmath/ACI_code
  commit: "<immutable upstream SHA>"          # to record (O2)
  file: dyad_interaction_model.m
  licence: MIT
harness:
  file: oracle/aci_oracle_dyad.m
  matlab_release: "<R2024b or R2025b as run>" # to record (O2)
  command: matlab -batch "run('aci_oracle_dyad.m')"
  seed: 333
  n_steps: 30001
  dt: 0.001
  parameters: {d_x: 0.5, d_y: 0.5, gamma: 2, F_x: 0.5, F_y: 1,
               sigma_x: 0.5, sigma_y: 1}
  sampled_indices: seq(1, 30001, by = 100)    # 301 rows
fixtures:
  dyad_signal_x.csv:
    sha256: f2b4d5284e19151432d1e1216619c15e00183704fde4d4ba6e378618a648d602
  dyad_reference.csv:
    sha256: 6e7cebd1ccf5b0ac99e9e4ca78f3772c978c0706ae60e976650ae1350cc80090
validation:
  package_version: 0.1.0
  max_abs_error: "<recorded at release>"
  tolerance_gate: 1.0e-6
```

## 8. Out of scope for 0.1.0 (roadmap, unchanged from NEWS)

Time-varying `L_y` core + noisy predator-prey constructor (needs its own
oracle fixture); causal-influence-range `aci_cir()` (O(N^2), own fixture);
vector-state CGNS; Rcpp acceleration (revisit only if WP4 benchmark says so);
parameter estimation.

## 9. Open decisions (Max)

- **O1 (blocks WP5 only):** repository home. Personal GitHub vs AAGI-AUS org.
  Note: an AAGI-AUS home activates the AAGI canon, including the
  GRDC-Manager licence consultation gate *before* the repo goes public --
  MIT is already declared, so this must be resolved before pushing, not after.
  Personal home: no gate, proceed with MIT.
- **O2 (optional, WP3.3):** run MATLAB once more for (a) the upstream commit
  SHA + MATLAB release fields of the manifest, and (b) optionally a second
  correlated-noise fixture. The analytic T13 test proceeds regardless.
- **O3 (confirm):** removal of `aci_simulate_dyad()` (D2) -- recommended, but
  it deletes a public symbol from the current draft.
- **O4 (confirm):** R floor `>= 4.1.0` once CI oldrel is green.
