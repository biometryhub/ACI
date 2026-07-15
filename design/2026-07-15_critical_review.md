# Critical and constructive review of `aciR`

**Review date:** 15 July 2026  
**Package version:** 0.0.0.9000  
**Review scope:** package structure, numerical core, model and simulation APIs,
tests, independent-oracle evidence, documentation, metadata, CI/release readiness,
maintainability, usability and likely scientific failure modes.

## Executive verdict

`aciR` has a small, legible and unusually well-grounded numerical core. Its
strongest asset is the committed MATLAB-oracle fixture and the mandatory tests
that compare all five principal output series at 301 points to an absolute
tolerance of `1e-6`. The package builds, installs, runs its examples, passes its
tests and rebuilds its vignette. Measured line coverage is 99.57%.

The package is nevertheless **not ready for a public scientific release**.
This is not because the validated dyad calculation is visibly wrong. The main
risk is that the general public API accepts mathematically inadmissible models
and malformed coefficient functions, then returns `NA`, `NaN`, non-positive
covariances or low-level replacement errors instead of rejecting the inputs.
Release metadata, CI, API-stability policy and scientific limitations also need
to be completed.

The recommended release sequence is:

1. Make invalid models impossible to construct.
2. Add shared validation and fail-fast numerical diagnostics to every exported
   low-level function.
3. Add tests for mathematical identities, invalid domains and discretisation
   sensitivity, while retaining the oracle gate unchanged.
4. Consolidate the duplicate dyad APIs and publish an API-stability policy.
5. Complete metadata, CI, pkgdown and release checks.

With those changes, the package would have a credible route to a strong 0.1.0
research preview.

## What was inspected and run

The review examined all files under `aciR/R`, `aciR/tests`, `aciR/vignettes`,
`aciR/man`, package metadata, the root `oracle` harness and both copies of the
oracle fixtures.

Commands and outcomes:

| Check | Outcome |
|---|---|
| `R --vanilla CMD build aciR` | Passed; produced the source tarball. |
| `R --vanilla CMD check --as-cran --no-manual aciR_0.0.0.9000.tar.gz` | 0 errors, 0 warnings, 2 notes. |
| Examples | Passed. |
| Tests | Passed; 18 `test_that()` blocks. |
| Vignette build and rebuild | Passed. |
| `covr::package_coverage()` | 99.57%; only `R/aci-model.R:346` was uncovered. |
| `goodpractice::gp()` | No findings. |
| `lintr::lint_package()` | 23 diagnostics: 18 mathematical-name style findings, 2 indentation findings and 3 object-usage findings. |
| `spelling::spell_check_package()` | Domain vocabulary and proper names require a package wordlist; no clear prose defect was established. |
| Sensitive-path scan | No personal filesystem paths or working-note filenames found in the package. |
| Fixture identity | Packaged and root oracle CSVs have identical SHA-256 hashes. |

The two check notes were:

- CRAN incoming URL/ORCID checks could not run because the review environment
  had no network access; the development version also contains a large version
  component (`9000`).
- The current time could not be verified by the environment.

The ordinary `R CMD build` command, without `--vanilla`, was disrupted by a
user-level R startup profile that attempted to create a parallel socket cluster.
That is not a package defect, but it is a reproducibility warning: release and
CI commands should run in clean sessions.

## Major strengths worth preserving

### 1. Independent-oracle evidence is a real differentiator

`aciR/tests/testthat/test-oracle-dyad.R` always runs from installed package
fixtures. It compares filter mean, filter covariance, smoother mean, smoother
covariance and the ACI metric against MATLAB-derived reference output at 301
indices. It also applies an explicit maximum absolute-error gate. This is much
stronger than tests that merely assert output shape or self-consistency.

The root MATLAB harness is readable and records the reference seed, model,
equations, sampling plan and deterministic checksums. The packaged and root
fixtures are byte-identical:

- `dyad_reference.csv`: `6e7cebd1ccf5b0ac99e9e4ca78f3772c978c0706ae60e976650ae1350cc80090`
- `dyad_signal_x.csv`: `f2b4d5284e19151432d1e1216619c15e00183704fde4d4ba6e378618a648d602`

Keep this as a non-skipping release gate.

### 2. The package architecture is compact and understandable

The code separates:

- the general filter/smoother/metric core;
- model construction and high-level inference;
- the worked dyad model; and
- package documentation.

There are only nine exported functions plus two registered print methods, one
hard dependency (`stats`) and no native-code or system-library burden. The
public surface is comfortably below the A1 statistical-method budget.

### 3. The high-level workflow is coherent

`aci_dyad_model()` -> `aci_simulate()` -> `aci()` is easy to understand. The
`aci` object retains the model, observed signal, time, posterior trajectories,
metric and time step. The README and vignette demonstrate the same path.

### 4. Documentation is already better than a typical development scaffold

Every exported function and print method has an example. There are no
`\dontrun{}` or `\donttest{}` escapes. The vignette gives the governing SDEs,
causal question, simulation, inference and plots. The method paper has a DOI
and a dedicated citation entry.

### 5. Dependency discipline is excellent

The package depends only on `stats`, while documentation and testing tools are
correctly placed in `Suggests`. There is no evidence of unnecessary framework
adoption or dependency-driven API complexity.

## Prioritised findings

| ID | Priority | Finding | Risk | Effort |
|---|---:|---|---|---:|
| F1 | P0 | Invalid covariance systems can be constructed. | Incorrect or non-finite scientific output. | Small |
| F2 | P0 | General coefficient and observed-signal contracts are not validated. | Silent recycling, `NA` propagation and opaque errors. | Medium |
| F3 | P0 | Filter/smoother covariance failures are not detected at the failing step. | Invalid KL divergence and misleading downstream results. | Medium |
| F4 | P1 | Low-level exported functions have inconsistent or absent validation. | Public API is fragile outside the flagship path. | Medium |
| F5 | P1 | Release metadata and CI are incomplete. | Package is not discoverable, supportable or release-gated. | Medium |
| F6 | P1 | API stability is undeclared and dyad functionality is duplicated. | Divergence and avoidable compatibility burden. | Medium |
| F7 | P1 | Tests optimise line coverage more than domain coverage. | High coverage can mask untested scientific contracts. | Medium |
| F8 | P1 | Assumptions and causal-interpretation limits need a dedicated treatment. | Users may over-interpret model-dependent information as causal proof. | Medium |
| F9 | P2 | Oracle provenance is strong but not packaged as a reproducible evidence record. | Harder independent audit and future fixture refresh. | Small |
| F10 | P2 | Result UX lacks summary, plotting and tidy extraction methods. | Unnecessary user-side object manipulation. | Medium |
| F11 | P2 | Time-grid and RNG side effects need clearer contracts. | Reproducibility surprises. | Small |
| F12 | P3 | Style tooling needs project-specific configuration. | CI noise and inconsistent maintenance practice. | Small |

## Detailed findings and actionable fixes

### F1 — Reject mathematically inadmissible noise covariance systems

**Evidence:** `aciR/R/aci-model.R:118-145` validates scalar finiteness and only
requires `S_xoS_x > 0`. It does not require:

- `S_yoS_y >= 0`;
- `S_yoS_x == S_xoS_y` for the scalar covariance matrix; or
- positive semidefiniteness of the joint noise covariance.

The following invalid model currently constructs successfully and simulation
then produces `NaN`:

```r
bad <- aci_cgns_model(
  L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
  S_xoS_x = 1, S_yoS_y = -1
)
aci_simulate(bad, n = 3, seed = 1)
```

**Fix:** introduce one internal covariance validator used by the constructor
and any low-level component validator. For the scalar case, enforce:

```r
S_xoS_x > 0
S_yoS_y >= 0
isTRUE(all.equal(S_yoS_x, S_xoS_y, tolerance = tol))
S_xoS_x * S_yoS_y - S_yoS_x^2 >= -tol
```

If only symmetric covariance is supported, remove `S_xoS_y` as an independent
constructor argument in the next permitted breaking release and derive it from
`S_yoS_x`. Until then, reject disagreement explicitly.

Also validate `label` as one non-missing character value and document the
expected shape of `parameters`.

**Acceptance criteria:** construction fails with specific, tested messages for
negative marginal variance, asymmetric cross-covariance and a negative
determinant; valid singular and positive-definite cases behave as documented.

### F2 — Validate coefficient functions and the observed signal before inference

**Evidence:** `.aci_as_coef()` checks that an input is a function, but not what
the function returns. `aci()` calls each function on all of `x` and assumes a
finite numeric vector of exactly `length(x)`. A scalar-returning function is
accepted and later creates missing indexed coefficients. `aci()` also accepts
non-finite observations.

Observed behaviour:

- a coefficient function `function(x) 1` is accepted and eventually produces
  `NA` filter entries;
- `aci(c(1, NA, 3), aci_dyad_model())` returns an `aci` object containing only
  missing metric values instead of rejecting the input.

**Fix:** add shared helpers such as `.check_observed_signal()` and
`.evaluate_coefficient()`.

For `x`, require an atomic, unclassed or explicitly supported numeric vector,
length at least two, all finite, and no missing values. For every evaluated
coefficient, require numeric type, exact length `n`, and finite values. Do not
allow implicit scalar recycling from user-supplied functions; constant scalars
are already supported directly by the constructor.

**Acceptance criteria:** every malformed input fails before the filter starts,
and the message names the coefficient, expected length and received type/length.

### F3 — Fail at the first non-positive or non-finite covariance

**Evidence:** `aciR/R/aci-core.R:88-96` and `aciR/R/aci-core.R:148-157` update
Euler approximations without checking that the covariance remains finite and
strictly positive. `aci_metric()` then divides by filter covariance and logs a
covariance ratio without checking its domain.

This can arise from an invalid model, a time step that is too large, or a valid
continuous-time system whose explicit Euler discretisation becomes unstable.
Returning a metric in this state would be scientifically unsafe.

**Fix:** after each filter and smoother update, check the covariance. On
failure, stop with:

- the algorithm (`filter` or `smoother`);
- index and corresponding time;
- failed value; and
- guidance to inspect model covariance and reduce `dt`.

Do not silently clip substantive negative covariance to a small positive value.
If tiny negative values attributable to floating-point error are tolerated,
make the tolerance explicit and test it.

For the KL dispersion term, use a cancellation-resistant form near a ratio of
one:

```r
ratio_delta <- smooth$cov / filt$cov - 1
dispersion <- 0.5 * (ratio_delta - log1p(ratio_delta))
```

Validate both posterior structures, equal lengths, finite means and strictly
positive covariances first. Clamp only tiny negative final values within a
documented round-off tolerance.

**Acceptance criteria:** unstable runs fail at the first bad index; valid runs
retain oracle agreement within `1e-6`; near-equal covariance tests do not
produce artificial negative KL values.

### F4 — Give every exported core function a complete contract

**Evidence:** `aci_filter()`, `aci_smoother()`, `aci_metric()`,
`aci_dyad_components()` and `aci_simulate_dyad()` are public but largely trust
their inputs. For example, `aci_filter(1, ...)` reaches an opaque replacement
error, and `aci_simulate_dyad(n = 1)` also fails inside its loop rather than at
the boundary.

**Fix:** decide which layer is public:

- If the low-level functions are a supported extension surface, validate them
  rigorously and document a versioned components schema.
- If they are implementation details, stop exporting them before 0.1.0 and
  expose one carefully designed extension mechanism instead.

The current documentation promises a general CGNS components interface, so the
first option is the less surprising path. Implement one `.check_components()`
helper that validates names, scalar/vector shapes, finiteness, covariance and
length alignment. Validate `filt` and `smooth` as named posterior structures.

**Acceptance criteria:** direct calls and high-level calls enforce the same
mathematical contract and produce the same class of error messages.

### F5 — Complete release metadata and continuous integration

**Evidence:** `DESCRIPTION` has no `URL`, `BugReports` or tested R version floor.
The author has `aut` and `cre` roles but not `cph`, although the licence names
the same person as copyright holder. `_pkgdown.yml` has `url: ~`. There are no
GitHub Actions workflows. `CITATION.cff`, `codemeta.json`, `cran-comments.md`,
`CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` and `API_STABILITY.md` are absent.

There is currently no git remote, so a canonical URL cannot be inferred safely.

**Fix:** once the repository home is chosen:

1. Add canonical repository and site URLs to `DESCRIPTION`, `_pkgdown.yml`,
   README, citation metadata and pkgdown configuration.
2. Add `BugReports: <repository>/issues`.
3. Add `cph` to the appropriate `Authors@R` entry.
4. Choose an R floor based on a tested compatibility target, not the current
   development machine.
5. Generate and then hand-review `CITATION.cff` and `codemeta.json`.
6. Add R CMD check CI across Linux, macOS and Windows and R old-release,
   release and devel; include coverage and pkgdown jobs.
7. Run release checks from a clean vanilla session and record justified notes
   in `cran-comments.md`.

**Acceptance criteria:** CI is green on the declared matrix; repository, site,
issues, version and attribution surfaces agree; `R CMD check --as-cran` has no
package-caused errors, warnings or notes.

### F6 — Consolidate the dyad API and declare stability before users depend on it

**Evidence:** two overlapping dyad paths exist:

- `aci_dyad_model()` + `aci_simulate()`; and
- `aci_dyad_components()` + `aci_simulate_dyad()`.

The second simulation independently reimplements the first, has weaker
validation, returns a list instead of a data frame and has no `seed` argument.
This creates two places for equations and defaults to drift.

The README labels the package experimental, but there is no statement of which
symbols are stable or how changes will be handled.

**Fix:** make the model-object path canonical. Rewrite legacy helpers as thin
wrappers over it so there is one implementation. Before 0.1.0, decide whether
the component helpers are:

- a supported expert API;
- experimental; or
- internal.

Publish `API_STABILITY.md` with lifecycle stages, compatibility promises and a
deprecation timetable. If retaining the old helpers, add lifecycle badges and
snapshot their warnings before eventual removal.

**Acceptance criteria:** one numerical implementation per operation; identical
results under matched seed and parameters; documented lifecycle status for
every export.

### F7 — Convert high line coverage into high domain coverage

**Evidence:** 99.57% line coverage is excellent but several invalid domains and
scientific identities are untested. The two oracle files also duplicate much of
the same comparison; `test-oracle-grade.R` skips when the root workspace is not
available, while `test-oracle-dyad.R` already provides the stronger installed
package gate.

**Fix:** retain the packaged oracle test and add the following focused tests:

1. **Zero-information identity:** with zero observation coupling and zero noise
   cross-covariance, filter and smoother should agree and ACI should be zero.
2. **Terminal identity:** smoother equals filter at the final index and terminal
   ACI is zero.
3. **Closed-form/simple-system comparison:** validate a constant-coefficient
   linear Gaussian case against an analytical or established Kalman result.
4. **Covariance admissibility:** constructor rejects all invalid covariance
   combinations.
5. **Coefficient contract:** wrong type, wrong length, `NA`, `Inf` and
   non-vectorised functions fail clearly.
6. **Time-step stress:** document or test convergence as `dt` is reduced and
   failure when explicit Euler loses positivity.
7. **Posterior contract:** mismatched lengths and non-positive covariance fail
   before metric calculation.
8. **Error snapshots:** lock the public error interface.
9. **RNG behaviour:** deterministic same-seed output and declared handling of
   global RNG state.

Replace the workspace-skipping oracle test with either a fixture-hash check or
a deliberately separate developer script outside the package test suite.

**Acceptance criteria:** coverage remains high, but the test plan explicitly
maps each scientific assumption and public input class to at least one test.

### F8 — State the scientific assumptions and limits of causal interpretation

**Evidence:** the vignette explains the mechanism clearly but makes strong
statements such as large ACI marking strong causal coupling. The software result
is conditional on the supplied CGNS being an adequate model, regular sampling,
known parameters, the scalar latent structure, Gaussian conditional laws and
the numerical discretisation.

**Fix:** add an *Assumptions and interpretation* vignette or prominent section
covering at least:

- what causal estimand is represented;
- directionality and what is conditioned on;
- model correctness and identifiability assumptions;
- regular, complete observations and fixed `dt`;
- known versus estimated model coefficients;
- scalar-only current scope;
- sensitivity to initial mean/covariance and time step;
- distinction between ACI evidence under the model and causal identification
  from observational data generally; and
- diagnostics or sensitivity analyses users should report.

Add a validation vignette that presents the oracle design, maximum errors,
hashes and limitations. This would turn the strongest engineering asset into a
visible scientific asset.

**Acceptance criteria:** a reader can state when the method is valid, when the
software should refuse to run and what claims an ACI peak does and does not
support.

### F9 — Make oracle provenance durable and independently auditable

**Evidence:** the root MATLAB file records the source repository in a comment,
but the package carries only the CSVs. It does not ship a machine-readable
manifest containing source revision, MATLAB version, generation command,
fixture hashes, parameterisation and licence/provenance statement.

The README says the core is validated against “the authors' reference
implementation”. The current evidence more precisely shows agreement with an
independent MATLAB harness reproducing the deterministic core grounded in the
authors' reference code. The distinction should be explicit unless the fixtures
were produced directly by running an unmodified upstream implementation.

**Fix:** add `inst/extdata/oracle-manifest.yml` or a validation vignette with:

- upstream repository URL and immutable commit;
- upstream licence;
- exact upstream file/function mapped;
- MATLAB release and platform;
- generation command;
- parameters, seed and sampled indices;
- SHA-256 hashes; and
- maximum observed errors for the released package version.

Prefer a refresh script that generates to a temporary directory, compares
hashes, and requires an explicit maintainer step to replace fixtures.

**Acceptance criteria:** another researcher can reconstruct why each fixture is
trusted and determine whether a changed fixture is intentional.

### F10 — Add scientific result methods without enlarging the core unnecessarily

**Evidence:** `print.aci()` prints a useful compact view, but there is no
`summary.aci`, `plot.aci`, `as.data.frame.aci` or direct diagnostic method.
Users must know and manipulate the nested representation.

**Fix:** consider a small method layer:

- `summary.aci()` returning structured metric, covariance and stability
  summaries rather than only printing `summary(x$aci)`;
- `as.data.frame.aci()` returning time, observed signal, filtered/smoothed means
  and covariances, and ACI;
- `plot.aci()` using base graphics to avoid a hard dependency, or a guarded
  ggplot2 method; and
- a diagnostic helper reporting minimum covariance, terminal identity and
  non-finite values.

Keep these methods pedagogical: explain the result and flag assumption-sensitive
quantities rather than merely dumping fields.

**Acceptance criteria:** the common inspect, plot and export workflows require
no direct indexing into internal lists, while the underlying object remains
simple and stable.

### F11 — Clarify the time-grid and random-number contracts

**Evidence:** `aci()` creates time from `dt` and assumes equally spaced complete
observations, but the docs do not foreground this limitation. There is no path
for an observed time vector or a check that externally supplied data use the
declared step.

`aci_simulate(seed = ...)` calls `set.seed()` directly and leaves the global RNG
state changed. `aci_simulate_dyad()` instead requires callers to seed globally.

**Fix:** document regular sampling prominently. Consider accepting `time = NULL`
at the high level, validating strict monotonicity and equal spacing, and deriving
`dt` when supplied. Reject irregular grids until the algorithm explicitly
supports them.

Choose one RNG contract. A contained seed (for example with `withr::with_seed()`)
is friendlier in simulation APIs, but a documented global-state contract is also
defensible if avoiding another import is more important. Apply the same contract
to both simulation entry points.

**Acceptance criteria:** time and RNG behaviour are explicit, consistent and
covered by tests.

### F12 — Configure style checks for mathematical software

**Evidence:** `lintr` reports 23 diagnostics. Most object-name findings concern
paper-aligned notation such as `L_x`, `S_xoS_x`, `R0`, `A_j` and `B_j`; renaming
all of these mechanically could reduce traceability to the method equations.
There are also two genuine indentation findings and three object-usage findings
that appear related to package load order rather than missing functions.

**Fix:** add a project `.lintr` that:

- preserves deliberate mathematical notation through narrow exclusions;
- keeps meaningful indentation, object-usage, sequence and scalar checks;
- treats new unexplained diagnostics as CI failures; and
- avoids a blanket exclusion of whole files.

Add `inst/WORDLIST` for ACI, CGNS, author names, mathematical terms and domain
vocabulary. Apply `styler` only after reviewing diffs so equation-aligned code
remains readable.

**Acceptance criteria:** lint and spelling jobs are quiet for documented domain
exceptions and fail on new actionable issues.

## Package-standard scorecard

This manual scorecard uses the statistical-method (A1) archetype. It is stricter
than the automated audit because it includes the style/release invariant and
manually evaluates the stochastic legacy API.

| Invariant | Status | Evidence and interpretation |
|---|---|---|
| I1 Metadata | **Partial** | Licence and ORCID exist; URL, BugReports, cph role and tested R floor are missing. |
| I2 Managed namespace | **Pass** | Roxygen-managed, explicit exports and S3 registrations. |
| I3 Lean hard dependencies | **Pass** | One justified import. |
| I4 Runnable examples | **Pass** | Examples on all exports/methods; no run-suppression escapes. |
| I5 Automated tests | **Pass** | testthat edition 3, 18 blocks, mandatory oracle tests, 99.57% coverage. |
| I6 CI matrix | **Fail** | No workflows. |
| I7 Structured NEWS | **Pass** | Current development entry and roadmap exist. |
| I8 Long-form documentation | **Pass** | Method paper, package citation and worked vignette; assumptions/validation vignettes remain desirable. |
| I9 API stability | **Fail** | Experimental badge only; no strategy or deprecation policy. |
| I10 Reverse-dependency process | **N/A now** | No released reverse dependencies; add a release checklist before publication. |
| I11 Reproducibility | **Partial** | High-level simulation has `seed`; legacy simulation does not; global RNG side effect is undeclared. |
| I12 Clean as-CRAN check | **Partial** | Package execution passed; two environment/development notes remain and manual PDF was not checked. |
| I13 Code/style/release hygiene | **Partial** | No sensitive leaks; 23 lint diagnostics and no committed lint policy. |

Strict current score: **6 Pass / 12 applicable = 50%**. This score reflects
release completeness, not the quality of the validated dyad calculation.

## A1 statistical-method uplift scorecard

| A1 uplift | Status | Evidence / next step |
|---|---|---|
| A1.1 Long-form method documentation | **Met** | Peer-reviewed method paper and worked vignette. Add assumptions and validation articles. |
| A1.2 Independent numerical validation | **Met** | Mandatory MATLAB-oracle fixture tests at `1e-6`. Add analytic/simple-system tests. |
| A1.3 Stochastic reproducibility | **Partial** | `aci_simulate()` has `seed`; align `aci_simulate_dyad()` and define RNG-state behaviour. |
| A1.4 Pedagogical print/summary methods | **Partial** | Print methods exist; add structured `summary.aci()` and diagnostics. |
| A1.5 API stability | **Missing** | Add stability policy and lifecycle treatment. |

## Proposed implementation plan

### Phase 1 — Correctness boundary (release blocker)

- Add `.check_noise_covariance()`, `.check_observed_signal()`,
  `.evaluate_coefficient()`, `.check_components()` and `.check_posterior()`.
- Enforce covariance symmetry and positive semidefiniteness in the constructor.
- Validate coefficient results at the `aci()` boundary.
- Add per-step finite/positive covariance checks to filter and smoother.
- Stabilise the KL dispersion calculation near covariance ratio one.
- Add focused failure and mathematical-identity tests.

**Gate:** all old oracle expectations still pass to `1e-6`; invalid domains fail
before returning a result.

### Phase 2 — API consolidation and scientific contract

- Make the model-object workflow canonical.
- Turn legacy dyad functions into validated thin wrappers or deprecate them.
- Add `API_STABILITY.md` and lifecycle badges.
- Add assumptions, interpretation and validation documentation.
- Add oracle provenance manifest and fixture hash test.

**Gate:** one implementation per operation; public extension contract and causal
claim boundary are documented.

### Phase 3 — User experience and diagnostics

- Add `summary.aci()`, `as.data.frame.aci()` and, if useful, `plot.aci()`.
- Add a regular-time-grid contract and optional validated `time` argument.
- Decide and test RNG-state behaviour.
- Add a small benchmark representative of 30,001 and larger trajectories.

**Gate:** common analysis and reporting workflows do not depend on internal list
structure; benchmark establishes whether optimisation is necessary.

### Phase 4 — Release engineering

- Choose canonical repository and pkgdown URLs.
- Complete DESCRIPTION, CFF/codemeta, contribution, conduct and support files.
- Add 3-OS / old-release-release-devel CI, coverage and pkgdown workflows.
- Add `.lintr`, `inst/WORDLIST` and release checklist.
- Run full `R CMD check --as-cran`, including PDF manual, in clean CI and a local
  vanilla session.

**Gate:** green CI, coherent metadata, clean source tarball and no unexplained
check output.

## Suggested 0.1.0 definition of done

- All P0 and P1 findings above are closed or explicitly documented with a
  maintainer-approved rationale.
- No invalid covariance or coefficient contract can reach numerical recursion.
- Oracle maximum absolute error remains below `1e-6`.
- Zero-information, terminal-identity, invalid-domain and discretisation tests
  pass on all supported R/OS combinations.
- Every export has a declared lifecycle stage.
- Assumptions and interpretation limits are visible from README and pkgdown.
- `R CMD check --as-cran` is clean apart from explicitly environmental incoming
  checks, and the PDF manual has been built.
- Package URLs, version and attribution agree across all metadata surfaces.

## Final assessment

The package is a promising research-software foundation, not a throwaway
prototype. Its numerical core is concise, the main workflow is elegant, and the
oracle evidence is worth building around. The immediate task is to turn the
implicit mathematical preconditions into enforced software contracts. Once that
boundary is hardened and the release surface is completed, `aciR` can credibly
present itself as a small, auditable and scientifically responsible
implementation of assimilative causal inference rather than merely a faithful
dyad reproduction.
