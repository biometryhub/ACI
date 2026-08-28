# DISPOSITIONS - the acir 0.1.0 audit ledger

One line per item: `item | origin (file:lines) | category | action | why | evidence ref`.

Everything is measured against **`aci` 0.0.30, git tree `97f6b124`**.

**Evidence refs.** Each row names the working record that backs it: **SEP** =
the separability report (the verified extraction map and its fifteen cut
points), **SCO** = the scope-partition report (constructors, tests, fixtures
and assets in both source packages), **PLAN** = the 0.1.0 plan of record,
**MIG** = the migration report, **VER** = the verification report. Those five
records are retained with the project and are not published in this directory;
this register is the published summary of what they conclude.

**The extraction was faithful.** No numerical, algorithmic or policy change was
made in it. Every `action = mainline` line below means "carried verbatim from
`aci` 0.0.30 except the enumerated cuts and their one-line adaptations", and
the acceptance bar, 116 quantities bit-identical to `aci` 0.0.30 on the
retained surface, was met (VER §3).

**Currency.** Sections 1 to 8 record the extraction and do not move. Sections 9
and 10 record the state at the close of the extraction; later work has since
closed several of their rows, and each such row says so inline.

---

## 1. Mainline - R sources

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `R/utils.R` (331) | `aci/R/utils.R` | ACI_code | mainline | byte-identical; zero excluded-file references | VER §3b |
| `R/compiled_scalar.R` (575) | `aci/R/compiled_scalar.R` | ACI_code | mainline | byte-identical; the `dyad_observed_x_v1` fast path | VER §3b, SEP §7 |
| `R/compiled_matrix.R` (420) | `aci/R/compiled_matrix.R` | ACI_code | mainline | byte-identical | VER §3b |
| `R/compiled_conditioning.R` (434) | `aci/R/compiled_conditioning.R` | ACI_code | mainline | byte-identical; conditional ACI for the `h_W`/`T_C`/`u` scripts | VER §3b, SCO B1 |
| `R/model_classes.R` (699) | `aci/R/model_classes.R` (705) | ACI_code | mainline | roxygen-only delta: `@export`→`@noRd` on `stochastic_model`, two `[enkbf()]` `@seealso`, ensemble prose | SEP §5(c), VER §3b |
| `R/assimilation.R` (769) | `aci/R/assimilation.R` (1025) | ACI_code | mainline | filter/smoother/lag-table core; five blocks excised (rows 3.x below) | SEP §7 cuts 4-6, MIG |
| `R/causal_metrics.R` (613) | `aci/R/causal_metrics.R` (1162) | ACI_code | mainline | ACI metric + forward CIR; ensemble arm and backward family excised | SEP §7 cuts 1-3, §3 |
| `R/compiled_lag.R` (398) | `aci/R/compiled_lag.R` (444) | ACI_code | mainline | Theorem-3 lag core; `one_lag` and unreachable `smoother_only` removed | SEP §3.4, §4 note |
| `R/compiled_cir.R` (290) | `aci/R/compiled_cir.R` (374) | ACI_code | mainline | forward CIR reducer; the two backward slice helpers removed | SEP §3.5 |
| `R/benchmark_models.R` (333) | `aci/R/benchmark_models.R` (1151) | ACI_code | mainline | `.complete_scalar_params`, `model_dyad(p1,x)`, `model_enso6(aci_code)`, `model_predator_prey` only | SCO A1, A5, A6 |
| `R/api_documentation.R` (116) | `aci/R/api_documentation.R` (188) | ACI_code | mainline | four doc blocks deleted, ensemble/backward prose trimmed | SEP §7 cut 8 |
| `R/plots.R` (73) | `aci/R/formula_interface.R:613-674` | ACI_code | mainline | three engine `plot` methods were stranded in an excluded file; moved verbatim rather than lost | SEP §6.2, cut 7 |
| `R/acir-package.R` (47) | `aci/R/aci-package.R` (54) | ACI_code | mainline | layer map rewritten, imports narrowed; comment-only | SEP §7 cut 9 |

## 2. Mainline - package metadata, docs, tests, assets

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `NAMESPACE` (22 exports, 25 S3, 3 importFrom) | `aci/NAMESPACE` (75/49) | ACI_code | mainline | SEP §5(a) keep list minus `stochastic_model`; hand-written, reproduced line-for-line by roxygen2 7.3.3 | SEP §5, VER §3b |
| `man/` (43 `.Rd`) | `aci/man/` (119) | ACI_code | mainline | retained topics only; independently regenerated, 0 of 43 differ | SCO C, VER §3b |
| `DESCRIPTION` | `aci/DESCRIPTION` | ACI_code | mainline | re-scoped Title/Description/Authors; `Imports` narrowed to `stats, graphics, utils` (no unused entry, `R CMD check` OK). **Awaiting author sign-off** | MIG, VER §RESIDUALS 7 |
| `LICENSE`, `LICENSE.md` | `aci/LICENSE*` | ACI_code | mainline | adapted: EnKBS holder and the `'FBCIR_code'` mention dropped. **Awaiting author sign-off** | MIG |
| `inst/CITATION` | `aci/inst/CITATION` | ACI_code | mainline | re-scoped to the ACI paper. **Awaiting author sign-off** | SCO C, MIG |
| `tests/testthat/` 15 test files + 2 helpers + `testthat.R` | `aci/tests/testthat/` (20 + 3) | ACI_code | mainline | engine-only set; `helper-load.R`, `helper-golden-p1.R`, test-10/15/16/17/19 byte-identical | SCO B1, B3, VER §3b |
| `tests/testthat/fixtures/oracles/` 16 CSV + `oracle-manifest.yml` | same in `aci` | ACI_code | mainline | all 17 sha256-identical; **shared evidence with `aciR`, not a second oracle** | SCO B2, MIG |
| `benchmarks/` (5 entries) | `aci/benchmarks/` | ACI_code | mainline (dev asset, `.Rbuildignore`d) | `scalar-dyad.R` benchmarks the `model_dyad` scalar path; copied unadjusted, so it still calls `aci::` and will not run against `acir` unretargeted | SCO C, MIG, VER §RESIDUALS 8 |
| vignettes | `aci/vignettes/` (2) | ACI_code / X | **not migrated** | vignette-1 is built on the excluded `aci_fit()` front-end; vignette-2 is mixed and needs a rewrite. `VignetteBuilder` omitted, `knitr`/`rmarkdown` dropped from Suggests | SCO C, PLAN 3.6 |

## 3. Reserved - `reserve/fbcir/` (FBCIR_code, `andreou2026cir`; release 0.2.0 or 0.3.0)

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `backward_cir` generic + `.lag_table`/`.aci_result` methods | `causal_metrics.R:780-979` | fbcir | reserved-to-fbcir | backward CIR appears only in `FBCIR_code-main/climate_tipping_y_bifurcation_driven.m` (eq. 21); ACI_code has no backward CIR | SEP §3, SCO preamble |
| `backward_cir.cgns_model` | `causal_metrics.R:1024-1038` | fbcir | reserved-to-fbcir | dies with the generic | SEP §3 |
| `.bwd_lengths` | `causal_metrics.R:553-573` | fbcir | reserved-to-fbcir | backward-only reducer, `andreou2026cir` eq. 13 | SEP §4 |
| `cir_pair` | `causal_metrics.R:982-1021` | fbcir | reserved-to-fbcir | runs both directions; a forward-only survivor is a one-line wrapper | SEP §5(c) |
| `lt_onelag` | `assimilation.R:486-513` | fbcir | reserved-to-fbcir | only consumers are backward CIR and the one-lag arm of `as.data.frame.lag_table` | SEP §5(c) |
| `as.data.frame.lag_table` one-lag arm | `assimilation.R:709-712` | fbcir | reserved-to-fbcir | unreachable without `mode = "one_lag"` | SEP §3.7 |
| `.lagtable_core_compiled` `mode == "one_lag"` block | `compiled_lag.R:126-156` | fbcir | reserved-to-fbcir | self-contained early return feeding backward CIR only | SEP §3.4 |
| `.slice_compiled_cgns`, `.slice_compiled_filter` | `compiled_cir.R:293-374` | fbcir | reserved-to-fbcir | exist solely for `backward_cir.aci_result` prefix sharing | SEP §3.5 |
| `model_tipping_triad` | `benchmark_models.R:179-282` | fbcir | reserved-to-fbcir | `andreou2026cir` eqs. (24a)-(24c); `climate_tipping_*.m` | SCO A2 |
| `model_multiscale_fbcir` + `.fbcir_multiscale_parameters` | `benchmark_models.R:816-1050` | fbcir | reserved-to-fbcir | FBCIR eq. (25) and the three `multiscale_*.m` | SCO A2 |
| `model_topographic_layered_fbcir` + its parameter helper | `benchmark_models.R:1052-1151` | fbcir | reserved-to-fbcir | `FBCIR_code-main/topographic.m` repository case study | SCO A2 |
| `model_l84(variant = "fbcir")` | `benchmark_models.R:346-350,378-383` | fbcir | reserved-to-enkbs (pointer) | inseparable arm of a constructor whose default is EnKBS; filed once at `reserve/enkbs/code/model_l84.R` | SCO A5 |
| backward test blocks in test-03/04/18/20 | `test-03-engine.R:140-159`, `test-04-cir.R` (6 blocks), `test-18-compiled-lag.R:136`, `test-20-compiled-cir.R:144-183` | fbcir | reserved-to-fbcir | assert the backward family | SCO B1(c) |
| `test-13-fbcir-models.R` | whole file | fbcir | reserved-to-fbcir | drives the two FBCIR constructors only | SCO B1(b) |
| `FBCIR-MULTISCALE.md`, `generate_fbcir_multiscale_oracles.m` | `tests/testthat/fixtures/oracles/` | fbcir | reserved-to-fbcir | referenced by no test and no R file; the generator needs MATLAB and has never been run | SCO B2 |
| `fbcir.patch` | generated | fbcir | reserved-to-fbcir | applies cleanly to `main`; installs; passes. Regenerated since: the shipped patch is 1031 lines at **pass 6397 / fail 0 / error 0**, logged in `reserve/fbcir/NOTES.md` | `reserve/fbcir/NOTES.md` |

## 4. Reserved - `reserve/enkbs/` (EnKBS_code, `jiang2026enkbs`; release 0.2.0 or 0.3.0)

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `R/ensemble.R` (746) | whole file | enkbs | reserved-to-enkbs | EnKBF/EnKBS engine; the only excluded file the engine calls into | SEP §1 |
| Gaspari-Cohn localization + inflation | `assimilation.R:902-1025` | enkbs | reserved-to-enkbs | excluded-family code living in a kept engine file; sole consumers `ensemble.R` and test-06 | SEP §1 reverse-located |
| `aci()` `engine == "ensemble"` arm | `causal_metrics.R:383-415` | enkbs | reserved-to-enkbs | the single code-level engine→excluded dependency in the package | SEP §1, §7 cut 1 |
| `forward_cir.aci_result` ensemble guard | `causal_metrics.R:766-769` | enkbs | reserved-to-enkbs | unreachable once the ensemble engine is out | SEP §7 cut 3 |
| `R/discovery.R` (759) | whole file | enkbs | reserved-to-enkbs | model-discovery family; no engine file calls into it | SEP §1 verified negatives |
| `model_l96` | `benchmark_models.R:391-491` | enkbs | reserved-to-enkbs | `jiang2026enkbs` eq. (16); `preset = "legacy"` is package-only but inseparable | SCO A3 |
| `model_l84` (both variants) | `benchmark_models.R:321-388` | enkbs | reserved-to-enkbs | default variant `p3` is EnKBS; **no ACI_code arm exists**, so the constructor leaves acir entirely | SCO A5 |
| `model_dyad(variant = "p3")` | `benchmark_models.R:93,144,153-155` | enkbs | reserved-to-enkbs | published EnKBS dyad experiment; needed by test-14 | SCO A3 |
| `model_dyad` whole body (carrier file) | `benchmark_models.R:88-176` | enkbs | reserved-to-enkbs | the p3/p4/`observe="y"` arms interleave with the retained p1 arm and are not re-appliable in isolation | MIG deviation 3 |
| `test-06-ensemble.R`, `test-14-golden-enkbs.R`, `helper-golden-p3.R`, `test-08-discovery.R` | whole files | enkbs | reserved-to-enkbs | family test surface; test-14 carries the machine-precision dyad golden grade | SCO B1(b), B3 |
| `test-03-engine.R:106-111` ensemble-mismatch block | mixed block | enkbs | reserved-to-enkbs | used `model_l84` + `enkbf`; the mainline block is retargeted to `model_dyad`, closed-form route only | SCO B1(c), MIG deviation 9 |
| `test-11-model-parity.R:8-16` p3 / observe-y block | mixed block | enkbs | reserved-to-enkbs | asserts the p3 preset and the reverse partition | SCO B1(c) |
| `test-02-models.R:50-58` (Z6-lite) | mixed block | enkbs (cross-family) | reserved-to-enkbs | asserts F, E and P constructors in one block; filed once with the EnKBS majority, cross-referenced from `fbcir/tests/CROSS-FAMILY.md` and `paper-extremes/tests/CROSS-FAMILY.md` | SCO B1(c) |
| `test-11-model-parity.R:21-155` | mixed blocks | enkbs (cross-family) | reserved-to-enkbs | same; tipping-triad (F), pathways/topographic (P), L84 (E/F), L96 (E) | SCO B1(c) |
| `enkbs.patch` | generated | enkbs | reserved-to-enkbs | applies cleanly to `main`; installs; passes; two test blocks left commented with their cross-category dependency named. Regenerated since: the shipped patch is 2610 lines at **pass 6483 / fail 0 / error 0**, logged in `reserve/enkbs/NOTES.md` | `reserve/enkbs/NOTES.md` |

## 5. Reserved - `reserve/aci-paper/` (ACI-adjacent papers, no ACI_code MATLAB)

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `lt_contraction_certificate` | `assimilation.R:613-683` | aci-paper | reserved-to-aci-paper | `andreou2026smoother` eqs. 3.18-3.19; paper-derived, no ACI_code backing | SEP §5(c) |
| `model_enso6(variant = "cfy22")` (carrier: whole body) | `benchmark_models.R:513-654` | aci-paper | reserved-to-aci-paper | `aci`'s **default** variant is a `chen2022enso` transcription, not the MATLAB; acir keeps only `variant = "aci_code"`, so the default changes | SCO A1 |
| `test-03-engine.R:353-378` | block | aci-paper | reserved-to-aci-paper | tests the contraction certificate | SCO B1(c) |
| `test-02-models.R:79-87` cfy22 arm, `test-11-model-parity.R:223-226` | blocks | aci-paper | reserved-to-aci-paper | cfy22 construction and initial conditions | SCO B1(c) |

## 6. Reserved - `reserve/paper-extremes/` (`moser2026extremes`, no MATLAB)

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `R/extremes.R` (875) | whole file | paper-extremes | reserved-to-paper-extremes | fourth paper, no supplied MATLAB; comes after the three MATLAB-backed families | SCO A4, PLAN release map |
| `model_pathways` | `benchmark_models.R:285-318` | paper-extremes | reserved-to-paper-extremes | eqs. (4.3)-(4.4); "paper checked; no corresponding MATLAB supplied" | SCO A4 |
| `model_topographic` | `benchmark_models.R:657-750` | paper-extremes | reserved-to-paper-extremes | eqs. (4.5)-(4.13); **not** the FBCIR layered model | SCO A4 correction |
| `model_dyad(variant = "p4")` | `benchmark_models.R:94,145,157` | paper-extremes | reserved-to-enkbs (pointer) | inseparable arm; carrier body filed once under enkbs | SCO A5 |
| `test-07-extremes.R` | whole file | paper-extremes | reserved-to-paper-extremes | family test surface | SCO B1(b) |

## 7. Reserved / dropped - `reserve/extensions/` (package-only, no paper, no MATLAB)

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `R/formula_interface.R:1-610` | `aci/R/formula_interface.R` | extensions | reserved-to-extensions | `aci_fit()` front-end + 12 generics; no paper or MATLAB backing. The three engine `plot` methods at `:613-674` are mainline, not reserved | SCO B1(b), SEP §6.2 |
| `R/applied_workflows.R` (384) | whole file | extensions | reserved-to-extensions | applied on-ramps; lightest coupling of the six excluded files | SEP §2 |
| `R/validation_diagnostics.R` (763) | whole file | extensions | reserved-to-extensions | validation suite; package-only infrastructure | SCO B1(c) |
| `test-09-interface.R`, `test-12-extensions.R` | whole files | extensions | reserved-to-extensions | their test surface | SCO B1(b),(c) |
| `model_dyad(observe = "y")` | `benchmark_models.R:127-135` | extensions | reserved-to-enkbs (pointer) | `aci` labels it "package extension … not a supplied paper/MATLAB inference benchmark" (`:148-150`); inseparable arm | SCO A5 |
| `kl_increment` | `causal_metrics.R:112-131` | extensions | **dropped (DEAD)** | zero callers in `R/`, `tests/` or `vignettes/` at 0.0.30 - only its own definition and its `.Rd`. Text kept in `extensions/code/` for the audit trail only | SEP §5(c), obs. 1; PLAN C.2 |
| `cir_table` | `causal_metrics.R:1041-1162` | extensions | **dropped** | self-disclaimed in the object it returns (`:1154-1157`): "outside the supplied papers and MATLAB reference code" - the package's own evidence against a fidelity-scoped 0.1.0 | SEP §5(c), obs. 5 |
| `empirical_kl` | `causal_metrics.R:250-287` | extensions | **dropped (stub)** | its `estimator = "knn"` is an explicit not-implemented stub (`:270-271`, "v0.2 stub (SPEC-02)"); no in-package consumer; ensemble-flavoured condition class | SEP §5(c), obs. 4 |
| `truncation_profile` | `assimilation.R:585-610` | extensions | **dropped (no callers)** | zero callers in `aci/R/` and zero in `aci/tests/`; can return with adaptive-lag work | SEP §5(c), obs. 2 |
| `projected_kl` | `causal_metrics.R:210-247` | extensions | **dropped** | only in-package consumer is the excluded extremes family (`extremes.R:498,503,518,581`); effectively an extremes export sitting in an engine file. Filed here with its test, filed by origin like everything else; **it is a hard dependency of `paper-extremes/code/extremes.R`** and must return with that family | SEP §5(c), SCO B1(c) |
| `.lagtable_core_compiled` `mode == "smoother_only"` | `compiled_lag.R:105,158,220,259,270,288` | extensions | **dropped (unreachable)** | `.lag_table_compiled:316` never passes it; `"smoother_only"` appears nowhere else in `aci/R/`. Recorded as a diff, not as re-appliable code; removal verified behaviour-preserving | SEP §4 note, VER §3b(2) |
| `test-01-kl.R:23-31,48-52,66-70` blocks | mixed | extensions | reserved-to-extensions | `projected_kl` and `empirical_kl` assertions | SCO B1(c) |

## 8. Items left verbatim in the mainline although dead or stale

Recorded so the next pass does not mistake them for oversights. Under the
extraction rule these are *not* touched: the acceptance bar is bit-identity
with `aci` 0.0.30.

| item | origin | category | action | why | evidence |
|---|---|---|---|---|---|
| `.calc_tau`'s `direction` switch | `causal_metrics.R:509-515` | ACI_code | mainline (dead arm kept) | two dead lines; SEP marks them harmless and touching them could only add risk | SEP §3.2 |
| `new_cir_result`'s `bound` ternary | `causal_metrics.R:613-615` | ACI_code | mainline (dead arm kept) | same | SEP §3.3 |
| `lt_row`'s "mode = 'one_lag'" error string | `acir/R/assimilation.R:518` | ACI_code | mainline (fixed, session 16) | parenthetical dropped, guard kept: "Table has no stored rows." is accurate in both mainline and family contexts | VER §RESIDUALS 4 |
| `lag_table$onelag` field, always `NULL` | `assimilation.R` | ACI_code | mainline (kept, by finding) | session 16 established it is the fbcir re-entry point (the family's one-lag mode produces it, `lt_onelag()` reads it); a removal broke the applied family patch and was reverted, the field now carries a comment naming the staging role | VER §RESIDUALS 5 |
| `aci_references` six-key `\describe` list | `api_documentation.R` | ACI_code | mainline (kept verbatim) | its keys are still used in retained source comments and error strings (`jiang2026enkbs` at `model_classes.R:138`, `chen2022enso` in `benchmark_models.R`); re-scoping is an editorial call, not an enumerated cut. **Flagged for author decision** | MIG deviation 10 |
| `aci(mod, ob, m = 50)` message quality | `causal_metrics.R` | ACI_code | mainline (fixed, session 16) | staged-absence guard on the literal name `m` via `sys.call()`, `aci_error_not_implemented`, no argument added; the EnKBS family removes the guard when it restores `m` | VER §RESIDUALS 3 |
| `meta$source_model` closure environments | all constructors | ACI_code | mainline (structural) | closure bodies identical, `environment()` differs `namespace:aci` vs `namespace:acir`; unavoidable, non-numerical | VER §3 supplementary |

## 9. Decision defaults in force (PLAN Revision 3, section E)

**Read with the extraction rule.** D2, D5 and D8 were realised in the
extraction itself. D1, D3, D4, D6 and D7 were decisions taken at the same time
whose implementation was deliberately held back from the extraction, so that
each could land individually against the verified baseline. All five have since
been implemented; the "not yet implemented" entries below are the status at the
close of the extraction, kept as the historical record and annotated inline.

| id | decision | status in this tree |
|---|---|---|
| D1 | intended-behaviour + documented divergence for the two MATLAB source defects | adopted; not implemented at extraction close, mainline was verbatim `aci`. **Implemented since** |
| D2 | backward CIR | **resolved and realised**: excised from the mainline, filed in `reserve/fbcir/` with `fbcir.patch` for immediate re-application (supersedes the earlier keep-unexported recommendation) |
| D3 | strict covariance policy by default, explicit opt-in flooring | adopted; not implemented at extraction close. **Implemented since** |
| D4 | MATLAB/`aciR` count convention as the default subjective read-out, `aci` eq. G.7 convention available and named | adopted; not implemented at extraction close, the mainline then shipped `aci`'s convention only. **Implemented since** |
| D5 | API naming | `aci`-style naming was in force at extraction close, the mainline being `aci`'s surface by construction. **Superseded since**: the public surface was renamed to the `aciR`-style `aci_*` interface |
| D6 | uniform masking default + a MATLAB-compatibility first-step profile | adopted; not implemented at extraction close. **Implemented since.** Note the trap: the conditional scripts' first filter step uses an unmasked inverse before later slices are masked; that is a convention, not a kernel failure |
| D7 | evidence additions: the scalar-partition arrays become source-derived fixtures | adopted; not implemented at extraction close, acir then inheriting `aci`'s 16 CSV + manifest unchanged. **Implemented since** |
| D8 | name and attribution | **Decision: the package name is `acir`, final.** Rationale: renaming at the birth of a package is cheap and renaming after a CRAN release is not, so the name was fixed before any release. The CRAN case-insensitivity consequence against the existing `aciR` stands and is to be raised with the supervisor. DESCRIPTION / LICENSE / CITATION credit lines are drafted and **await author sign-off** |

**Note on D8.** The same case-insensitivity that affects CRAN also affected the build: on a case-insensitive filesystem a directory named `acir` beside the existing `aciR` tree is the same directory. The nested package root and the separate install library used during construction are a local workaround for that collision only, and carry no meaning for the published package.

## 10. Open, unclassifiable, or not closed by the extraction

| item | status |
|---|---|
| `T_C`-hidden zeroth-order ACI_code approximation | **Closed since.** At extraction close it was filed nowhere, because it existed in neither source package: `aci` hard-rejects `TC` as hidden (`benchmark_models.R:545-548`, rationale `:644-647`) and `aciR` computes the full `c_1` (`enso-model.R:196`). It was implemented afterwards as a fresh port of the reference script, and was the last open ACI_code completeness item |
| `reserve/fbcir-paper/`, `reserve/enkbs-paper/` | **Empty by finding, not by omission.** Every paper-only fragment in those two families (`model_tipping_triad(partition="joint")`, `model_l96(preset="legacy")`) is an inseparable arm of a MATLAB-backed constructor and travels with it. See each directory's README |
| `reserve/aci_code-future/` | **Category retired.** It was empty of code by finding, every ACI_code-scoped block in `aci` 0.0.30 having been retained in the mainline, and the `T_C` gap it recorded is closed (row above). No such directory ships |
| author-output parity for FBCIR | **open**: needs MATLAB to run `generate_fbcir_multiscale_oracles.m`, or an independent transcription of it. Neither was available during construction |
| EnKBS verification breadth (L96, `v -> u`, L84) | **open**: later work. The dyad golden grade (G3/G4) does travel with `enkbs.patch` |
| `fbcir.patch` + `enkbs.patch` stacking | **does not stack.** Each applies cleanly to a clean `main`; both edit `NAMESPACE`, `R/assimilation.R`, `R/causal_metrics.R` and overlapping tests, and `git apply` (plain and `--3way`) refuses the second. Manual merge required |
| `man/` for both patches | neither patch touches `man/`; the preview branches install with a short help index. `roxygen2::roxygenise()` regenerates |
| package path and install library | The nested package root and the separate preview install library used during construction were forced by the case-insensitive filesystem, which collapses `acir` and `aciR` into one directory. Both are local build arrangements, not properties of the package, and they go away once the two trees no longer sit side by side |
