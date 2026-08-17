# aciR — coding and package-quality review

**Review date:** 16 July 2026  
**Package version:** 0.1.0  
**Reviewed at commit:** `73196fd` (main; the state every finding below refers to)  
**Scope:** package architecture, numerical core, contracts and validation, tests and oracles, metadata/CI/release posture, maintainability, usability.  
**Companion document:** [scientific and vignette review](2026-07-16_scientific_vignette_review.md). Finding IDs are prefixed `PQ-` here and `SV-` there so cross-references stay unambiguous.

---

## Executive verdict

`aciR` 0.1.0 is an unusually strong research-preview R package. It is small (~1.8k lines of R including a large shared validation module), dependency-light (`Imports`: `graphics`, `stats` only), and organised so that the scientific core is separable from the user-facing workflow. The standout strength is **independent-oracle grading**: committed MATLAB fixtures, a never-skipping absolute-error gate of `1e-6`, observed agreement near `1e-14`, and an explicit manifest of what each oracle does *and does not* grade.

Relative to the 15 July critical review of 0.0.0.9000, the 0.1.0 release has closed essentially every P0/P1 software-correctness finding: inadmissible noise covariances are rejected, coefficient and signal contracts are enforced, filter/smoother covariances fail closed, dyad APIs were consolidated, `summary` / `plot` / `as.data.frame` exist, and API lifecycle is declared.

**Overall grade (coding / package engineering): A− / excellent research preview.**

| Dimension | Grade | One-line judgment |
|---|---|---|
| Architecture & API shape | A | Clean core / model / methods split; one coherent path |
| Numerical correctness evidence | A+ | Independent oracles + analytic identities + hash-pinned fixtures |
| Input contracts & fail-fast behaviour | A | Shared validators; messages are part of the interface |
| Test design | A | Domain tests dominate; coverage is a byproduct, not the goal |
| Documentation (Rd / README / NEWS) | A− | Thorough, honest about limits; README numbers are generated |
| Metadata, CI, release hygiene | A− | Full r-lib matrix incl. declared R floor; private-URL NOTE only |
| Performance & scalability | B+ | Pure-R Euler loops fine for paper-scale n; no vector path yet |
| Completeness vs method paper | B | CIR, vector state, time-varying L_y still roadmap — correctly so |
| Production / CRAN readiness | B+ | Ready as private research preview; public CRAN later |

**Bottom line:** ship 0.1.0 as a research preview with confidence in the *validated scalar CGNS core*. Do not oversell completeness relative to the full ACI paper (CIR and multivariate/conditional ACI are still out of scope). Remaining software work is mostly polish, discoverability, and careful extension under the same oracle discipline.

---

## What was inspected

| Artefact | Notes |
|---|---|
| `aciR/R/*.R` | All six source files (~1.8k LOC) |
| `aciR/tests/testthat/*` | Nine test files (~1.4k LOC); no skips |
| `aciR/vignettes/*` | Three vignettes (covered more deeply in the companion review) |
| `aciR/man/*`, `NAMESPACE`, `DESCRIPTION` | Exports, metadata, examples |
| `inst/extdata/*`, `oracle-manifest.yml` | Fixture provenance |
| Root `tools/oracle/`, `matlab_reference/` | MATLAB harnesses and authors’ reference |
| `aciR.Rcheck/00check.log` | Prior `--as-cran` run: 0 errors, 0 warnings, 3 notes |
| `.github/workflows/*` | R-CMD-check, coverage, lint, pkgdown |
| `API_STABILITY.md`, `NEWS.md`, `cran-comments.md` | Release policy |
| Method paper PDF | Cross-check of equations and claims (companion review) |

Local probe (source load via `devtools::load_all`): flagship dyad path, diagnostic summary, correlated-noise simulation rejection, and wall time for `n = 30001` (median 0.036 s over five runs on `aarch64-apple-darwin20`, R 4.5.2; artefact: `tools/design/artefacts/2026-07-16_probe.txt`).

---

## Package architecture

### Strengths

1. **Layered design.**  
   - Numerical core: `aci_filter` / `aci_smoother` / `aci_metric` on a general components list (`aci-core.R`).  
   - Model layer: `aci_cgns_model` / `aci_dyad_model` / `aci_simulate` / `aci` (`aci-model.R`).  
   - Worked dyad components: `aci_dyad_components` (`dyad-model.R`).  
   - Shared contract boundary: `aci-validate.R`.  
   - Reporting surface: `aci-methods.R`.  

   This is the right split for an A1 statistical-method package: experts can drop to the core; most users never leave `aci()`.

2. **Small public surface.** Nine functional exports plus S3 methods. No framework bloat, no compiled code, no heavy Suggests on the critical path.

3. **Paper-aligned naming.** Symbols (`L_x`, `S_xoS_x`, `A_j`, `B_j`) match the governing equations and the MATLAB reference. A narrow `.lintr` exception preserves that correspondence without disabling style checks wholesale.

4. **Single dyad implementation after 0.1.0.** Removal of `aci_simulate_dyad()` and derivation of `S_xoS_y` from `S_yoS_x` remove two silent-drift hazards that the previous review flagged.

### Minor architectural notes

| ID | Severity | Note |
|---|---|---|
| PQ-A1 | Low | Two routes to dyad components remain (`aci_dyad_components()` vs model coefficients). They are equality-tested; keep that test forever. |
| PQ-A2 | Low | `aci_components` is documented as an internal name with `@keywords internal` but is a public contract in `API_STABILITY.md`. That tension is intentional (schema doc, not a function); ensure pkgdown still surfaces it. |
| PQ-A3 | Info | Pure-R scalar Euler loops are clear and oracle-matched. If vector CGNS arrives, resist rewriting the scalar path “for performance” without a second oracle. |

---

## Numerical core and software correctness

### Strengths

1. **Filter / smoother match the MATLAB reference term-for-term** (scalar reduction of the matrix form `L_y R + R L_y` → `2 L_y R`). Oracle agreement at ~1e−14 over a 30,001-step trajectory is as strong as floating-point evidence gets for this class of method.

2. **Fail-closed covariance guards.** Non-positive or non-finite filter/smoother covariances stop at the offending index with algorithm, time, and value named. Relative entropy is never scored on an invalid posterior.

3. **Metric numerics.** Dispersion uses `ratio_delta - log1p(ratio_delta)` so that ratios near 1 do not cancel to exact zero; round-off negatives clamp inside `1e-10`; larger negatives error. Domain tests cover identical posteriors, near-unit ratios, and full-trajectory non-negativity.

4. **RNG containment.** Seeded `aci_simulate()` restores `.Random.seed`; unseeded use is ordinary. The seed→path map is correctly declared unstable across versions.

5. **Time-grid contract.** Strictly increasing, equally spaced `time` only; irregular grids are rejected rather than silently approximated.

### Residual software risks

| ID | Severity | Finding | Recommendation |
|---|---|---|---|
| PQ-S1 | Medium | **`n_clamped` conflates structural zeros with round-off clamps.** `summary.aci` counts `sum(aci == 0)`. The terminal metric is *exactly* zero by construction, so the count always includes at least that step even when no clamping occurred. The printed line says “metric at the round-off floor”, which overclaims. | Count only values that were clamped (`value < 0` before clamp), or rename the diagnostic to “metric exactly zero (includes terminal identity)”. |
| PQ-S2 | Medium | **`aci_simulate()` refuses correlated noise**, while filter/smoother support `S_yoS_x ≠ 0` and the cross oracle grades that path. Users who want to *generate* cross-noise data must leave the package. | Document as a known gap; optionally add Euler–Maruyama with a 2×2 Cholesky factor of the joint noise covariance, oracle-tested against `aci_oracle_cross.m`’s simulation block. |
| PQ-S3 | Low | **No intermediate progress / interrupt contract** for very long series. Fine at paper scale (`n ~ 3e4` runs in ~0.04 s; probe artefact above); less friendly if someone runs `n ~ 1e7`. | Not blocking; consider compiled or vectorised paths only after vector CGNS is designed. |
| PQ-S4 | Low | **`plot.aci` is Experimental** and base-graphics only. Restores `par`; colour choices are CVD-aware. Acceptable. | Keep Experimental until layout stabilises; vignette already shows ggplot recipe. |
| PQ-S5 | Low | **Coefficient functions are free-form.** Validation of return type/length/finiteness is excellent; there is no check that they are *vectorised* over a vector `x` at construction time (failure is deferred to first evaluation). | Optional: smoke-eval at construction with a tiny probe vector when `n` is known — probably overkill. |
| PQ-S6 | Info | **Integer handling of `n` in `aci_simulate`** uses a careful whole-number check; good. Elsewhere length comparisons use `==` on integers correctly. | — |

### Mathematical-software alignment (brief)

The discrete filter update

```text
mu ← mu + (L_y mu + f_y) dt + aux / Sxx * (dx − (L_x mu + f_x) dt)
R  ← R  + (2 L_y R + Syy − aux²/Sxx) dt
```

and the backward smoother (with `A_j`, `B_j` and the cross-noise transport term) match `matlab_reference/dyad_interaction_model.m` and the CGNS continuous-time structure used by the paper. The scalar KL form matches `KL(N_s ‖ N_f)` for one-dimensional Gaussians. **No coding-side transcription defect was found in the core.**

---

## Contracts, errors, and API stability

### Strengths

- Single validation module with magnitude-scaled tolerances (`.aci_tol`) — correct for dimensioned Grammians.
- Noise PSD enforced once and reused by model constructor and components validator.
- Coefficient evaluation rejects silent recycling of scalars returned from functions.
- Incomplete / non-finite signals are refused with a missing-observation message that does not pretend to impute.
- Error text is snapshot-tested (`test-errors.R` + `_snaps/errors.md`) and declared part of the interface in `API_STABILITY.md`.
- Lifecycle table: everything Maturing except `plot.aci` (Experimental); nothing falsely marked Stable. Honest for a 0.1.0 preview.

### Residual issues

| ID | Severity | Finding |
|---|---|---|
| PQ-C1 | Low | `API_STABILITY.md` is excellent but lives only in the package root; first-time users may miss it. Link from README “Scope and interpretation” is partial; a one-line pointer under Installation or Citation would help. |
| PQ-C2 | Low | Deprecation path for the 0.1.0 breaking changes is documented in NEWS with before/after code — good. Future removals should follow the same template. |
| PQ-C3 | Info | External function qualification (`stats::`, `graphics::`) without Imports-import of symbols is correct and CRAN-friendly. |

---

## Testing and validation evidence

### What is exemplary

| Layer | Evidence |
|---|---|
| Authors’ MATLAB oracle (dyad) | 301 indices × 5 series; max abs error gate `1e-6`; measured ~`4.57e-14`; high-level `aci()` path also graded |
| Cross-noise fixture | Independent harness (no upstream scalar model with `S_yoS_x ≠ 0`); measured ~`1.38e-14` |
| Analytic oracles | Riccati fixed point with non-zero cross-noise; smoother fixed point; terminal identity; zero-information *and* O(dt²) residual rate |
| Contract tests | Large `test-validate.R`; model, core, methods, errors |
| Provenance | `oracle-manifest.yml` + hash checks in `test-oracle-manifest.R`; refresh policy forbids in-place regeneration |
| Non-skipping policy | Tests use `system.file(..., package = "aciR")` only — works under R CMD check and after install |

This is the right evidence hierarchy for a method-reproduction package. **Self-consistency alone is never treated as validation.**

### Gaps in the test suite (software lens)

| ID | Severity | Gap |
|---|---|---|
| PQ-T1 | Medium | **No property-based / fuzz tests** on random admissible constant-coefficient systems (e.g. random PSD noise, random slopes) beyond the hand-chosen Riccati case. Low cost, high confidence for regressions. |
| PQ-T2 | Low | **No benchmark / regression timing test** (even a soft upper bound). Optional. |
| PQ-T3 | Low | **`n_clamped` semantics untested** against a case that actually hits the clamp path vs terminal zero. |
| PQ-T4 | Info | Coverage measured at 100.00% (`covr::package_coverage()`; artefact: `tools/design/artefacts/2026-07-16_coverage.txt`); the *important* fact is that ungraded *paths* (cross-noise) were fixed by new oracles rather than by chasing line coverage. Keep that discipline. |

---

## Metadata, CI, and release posture

### Strengths

- `DESCRIPTION`: Authors@R with ORCID and `cph`, SPDX-style MIT, URL/BugReports, R ≥ 4.1.0, Encoding UTF-8, Language en-AU, roxygen markdown.
- Citation dual-entry (package + method paper DOI).
- CI matrix: macOS/Windows release; Ubuntu devel/release/oldrel-1/**4.1** (declared floor built, not just claimed); `error-on: warning`.
- Companion workflows: test-coverage, lint, pkgdown.
- `NEWS.md` for 0.1.0 is a model of breaking-change documentation.
- `cran-comments.md` correctly records that 0.1.0 is **not** a CRAN submission and adjudicates notes in advance.
- `CODE_OF_CONDUCT`, `CONTRIBUTING`, `CITATION.cff`, `codemeta.json`, `_pkgdown.yml` present.

### Issues

| ID | Severity | Finding | Recommendation |
|---|---|---|---|
| PQ-M1 | Medium | **GitHub URLs return 404** to anonymous checkers because the repo is private. R CMD check NOTE and README badges claim a public Actions URL that outsiders cannot see. | Either make the repo public before advertising badges widely, or qualify install instructions (“private; needs access”) and accept the NOTE until then. Documented honestly in `cran-comments.md`. |
| PQ-M2 | Low | **GitHub Pages blocked** for private repos (project log). pkgdown builds but is unserved. | Public repo, paid plan, or stop implying a live site. |
| PQ-M3 | Low | **No revdep story** — fine for a new package with zero reverse depends. |
| PQ-M4 | Info | HTML Tidy environmental NOTE is machine-local; CI covers manuals. |
| PQ-M5 | Info | Version `0.1.0` is appropriate for research preview; avoid claiming 1.0.0 until CIR and/or external use exist (API_STABILITY already says this). |

---

## Usability and developer experience

### Strengths

- Canonical path is three calls: model → simulate → `aci()`.
- `print` / `summary` surface scientific diagnostics (min cov, terminal residual), not only pretty numbers.
- `as.data.frame()` enables tidy plotting without digging into nested lists.
- Rd examples are runnable (`\donttest` free); README is generated from `README.Rmd` so printed numbers cannot silently rot.
- Expert path (`comp` → filter → smoother → metric) is documented and vignette-demonstrated.

### Friction points

| ID | Severity | Note |
|---|---|---|
| PQ-U1 | Low | Users coming from potential-outcomes “causal inference” may still misread the package name despite excellent assumptions vignette. README already warns; keep that warning above the fold. |
| PQ-U2 | Low | Dual y-scale plotting of metric vs signal is left to the user; base `plot.aci` does the right two-panel split. |
| PQ-U3 | Info | Installation requires `subdir = "aciR"` because of monorepo layout — correctly documented. |

---

## Comparison to the 0.0.0.9000 critical review

| Prior finding | 0.1.0 status |
|---|---|
| F1 Inadmissible covariance constructible | **Fixed** |
| F2 Coefficient / signal contracts | **Fixed** |
| F3 Covariance failure not detected | **Fixed** |
| F4 Inconsistent low-level validation | **Fixed** |
| F5 Incomplete metadata / CI | **Fixed** (private URL remaining) |
| F6 Undeclared API + dyad duplication | **Fixed** (`API_STABILITY.md`, simulate consolidation) |
| F7 Coverage over domain tests | **Fixed** (identities + cross oracle) |
| F8 Causal interpretation limits | **Fixed** (assumptions vignette) |
| F9 Oracle provenance packaging | **Fixed** (manifest + hash tests) |
| F10 Missing summary / plot / tidy | **Fixed** |
| F11 Time grid / RNG | **Fixed** |
| F12 Style tooling | **Fixed** (`.lintr`) |

The package has executed its own uplift plan. Remaining work is not “make invalid input safe” — it is “extend the method without breaking the evidence model.”

---

## Prioritised recommendations (coding / package)

### Do soon (quality of life / correctness polish)

1. **Fix `n_clamped` semantics** (PQ-S1) — small, honest diagnostics matter for publication-facing `summary()`.
2. **Clarify private-repo discoverability** (PQ-M1) — badges and install story for external collaborators.
3. **Either implement or clearly “wontfix for 0.x” correlated-noise simulation** (PQ-S2) — closes the asymmetry between what the core can filter and what the package can generate.

### Do before any 1.0 / CRAN push

4. Keep oracle non-skipping + hash gates as release blockers (already policy).
5. Public repository (or drop public URL claims).
6. Live pkgdown or remove dead site expectations.
7. One external user dry-run of the assumptions vignette before promoting APIs to Stable.

### Explicitly defer (correct)

8. Time-varying `L_y`, CIR, vector states — only with new fixtures (NEWS roadmap).
9. Compiled speedups without paired oracles.
10. Interventional causal claims or automatic model checking — out of method scope.

---

## Suggested scoring rubric (for future audits)

| Criterion | 0.1.0 score |
|---|---|
| I1 DESCRIPTION metadata | Pass |
| I2 Managed NAMESPACE / explicit exports | Pass |
| I3 Lean dependencies | Pass |
| I4 Runnable examples | Pass |
| I5 Automated tests | Pass (strong) |
| I6 CI multi-OS / multi-R | Pass (incl. floor) |
| I7 NEWS structured | Pass |
| I8 Long-form docs | Pass (3 vignettes + pkgdown) |
| I9 API stability mechanism | Pass (declared Maturing) |
| I10 Revdep | N/A (new package) |
| I11 Reproducibility floor | Pass (seed plumbing; fixtures) |
| I12 Clean R CMD check | Pass (notes adjudicated) |

---

## Closing judgment

From a **software-engineering and package-soundness** perspective, `aciR` 0.1.0 is release-quality as a research preview. The numerical core is not merely tested — it is **independently graded**, with scope limitations stated in machine-readable form. Validation, error contracts, and documentation discipline are well above the typical academic-method package of similar size.

The main residual risks are **scope communication** (what the package does not yet implement relative to the full Nature Communications method) and a few **diagnostic/UX nits**, not latent numerical defect in the graded paths. Treat the oracle gate as sacred when extending the package.

---

**Post-review action (2026-07-16, same day).** PQ-S1 and PQ-T3 were fixed immediately after this review: `aci()` now records the true clamp count (`n_clamped`), decoupled from the terminal exact zero, with deterministic counting tests; `summary()` reports that count and refuses a pre-fix object rather than miscounting. SV-V1, SV-V4, SV-O2 and SV-R5 from the companion review were addressed in the same pass. A follow-up meta-review also found and closed a licence-attribution gap both reviews missed: `Authors@R`, `LICENSE` and `LICENSE.md` did not carry the upstream MATLAB implementation's copyright (Marios Andreou, MIT), which the oracle manifest and README already recorded. See `NEWS.md` (development version) and the commits following `73196fd`.

*See also:* [scientific and vignette review](2026-07-16_scientific_vignette_review.md).
