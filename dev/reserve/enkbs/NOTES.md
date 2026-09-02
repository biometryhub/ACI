# reserve/enkbs - verification status and regeneration rule

Family: `EnKBS-main/`, paper `jiang2026enkbs`. Everything here was excised
from `aci` 0.0.30 (git tree `97f6b124`) when acir 0.1.0 was cut.

## Regeneration log

| regenerated | against `main` | branch tip | patch |
|---|---|---|---|
| 2026-08-28, session 16, current | session-16 update | `0268fe2` | 2636 lines |
| 2026-08-28, post-file-and-argument rename | `c7b8720` | `9698d39` | 2610 lines |
| 2026-08-28, post-surface-rename | `3db97dd` | `0b8da79` | 2610 lines |
| 2026-08-28, pre-rename | `bded686` | `b67532e` | 2620 lines |
| at the 0.1.0 extraction | `8fba4bc` | `6d1b5a7` | 2625 lines |

**Three mainline commits killed the previous patch.** Measured with
`git apply --check` on `c7b8720` before regenerating:

```
error: R/assimilation.R: No such file or directory
error: R/benchmark_models.R: No such file or directory
error: R/causal_metrics.R: No such file or directory
error: R/model_classes.R: No such file or directory
```

The mainline renamed twelve of its thirteen `R/` sources to the `aci-*`
convention, so every edited target path in the patch named a file that no
longer exists. A second commit renamed `nontarget =` to `conditional =` and
`eps =` to `epsilon =`; a third was metadata only. Note that the file rename
alone is a harder failure than the surface rename before it: a surface rename
fails on context lines, which invites forcing, while a file rename fails on the
header and leaves nothing to force. Branch tips are recorded so the branches
are recreatable; the branches themselves are deleted after verification, and
the patch reconstructs them exactly.

**Two session-16 mainline changes killed the c7b8720 patch** (both land in
this branch's single session-16 update commit). Measured with
`git apply --check` against the session-16 tree before regenerating: the
predator-prey metadata change moved the `aci-model-library.R` context, and
the staged-absence guard moved the top of `aci()` in `aci-core.R`. The regeneration also adds adaptation E (below):
the family deletes that guard, whose whole purpose is this family's absence.

## What `enkbs.patch` restores

Generated `git diff main..preview/enkbs`, 2610 lines, over eleven files
(six of them new).

| area | restored |
|---|---|
| `R/aci-ensemble.R` (new) | the whole EnKBF/EnKBS engine (746 lines): `enkbf`, `enkbs`, `as_gaussian`, `ensemble_lag_table`, `.ensemble_lag_table_from_run`, `aci_filter.stochastic_model`, `aci_smoother.stochastic_model`, `print.da_path_ensemble` |
| `R/aci-discovery.R` (new) | the whole model-discovery family (759 lines): `cgns_library`, `eval_library`, `causation_entropy`, `threshold_structure`, `sample_paths`, `constrained_mle`, `model_from_learned`, `learn_model` |
| `R/aci-assimilation.R` | the Gaspari-Cohn localization and inflation section (`gaspari_cohn`, `localization_spec`, `apply_inflation`); `lag_table`'s ensemble error text; **removal** of the two not-implemented stubs, which the ensemble engine supersedes |
| `R/aci-core.R` | `aci()`'s `engine = "ensemble"` arm with `m`, `seed`, `localization`, `inflation`, `ic_sampler`, the engine-conditional metric and lag-table calls, and the `meta$m` / `smoother_scheme` fields; the ensemble guard at the head of `aci_range.aci_result` |
| `R/aci-model.R` | `@export` on `stochastic_model` |
| `R/aci-model-library.R` | `aci_dyad_model`'s full body (`p1`/`p3`/`p4`, `observe = "x"`/`"y"`); `model_l96` |
| `NAMESPACE` | 20 lines |
| tests (new) | `test-06-ensemble.R`, `test-08-discovery.R`, `test-14-golden-enkbs.R`, `helper-golden-p3.R` |

The discovery source and `test-08` are inside the patch because
`test-06-ensemble.R` calls `sample_paths()`; discovery has no dependency
outside the EnKBS family, so pulling it in restores real coverage rather than
deleting an assertion.

## Two blocks deliberately NOT restored (both marked in the test source)

1. The ensemble guard on the backward range in `test-06-ensemble.R`. The
   backward range is FBCIR-scope (`reserve/fbcir/`). Commented out with the
   dependency named. Apply `fbcir.patch` as well and uncomment to recover it.
   **The surface rename adds a second reason**: `aci_range()` now has a
   `direction = "backward"` release gate that raises
   `aci_error_not_implemented` on every result, ensemble or not, so restoring
   the line on this branch alone would assert nothing about the ensemble
   engine while appearing to pass.
2. `test-06-ensemble.R`'s last two blocks - `nil_causality_check` and
   `nil_surrogate_test`. These are `validation_diagnostics.R`, i.e.
   package-only **extensions**, not EnKBS science; restoring them inside an
   EnKBS patch would cross the filing categories. Commented out with the
   dependency named (`reserve/extensions/code/validation_diagnostics.R`).

Also not restored: any `man/` topic. The patch does not touch `man/`, so the
branch installs with a short help index; `roxygen2::roxygenise()` regenerates
it. `model_l84` (`code/model_l84.R`) and the cross-family test blocks in
`tests/` are outside the patch and return with the constructor work.

## The three mechanical maps

The reserve `.R` files are untouched: they stay the verbatim `aci` 0.0.30
excision. Three renames the mainline has since made are applied **on the way
in**, using the same rules and the same guards the mainline passes used, and
every substitution is counted so that drift shows up as a number rather than
as a load error.

**1. The public-surface map.** 63 substitutions across six reserve files:
`model_dyad`→`aci_dyad_model` 29, `da_smooth`→`aci_smoother` 11,
`da_filter`→`aci_filter` 7, `cgns_model(`→`aci_model(` 5,
`forward_cir`→`aci_range` 4, `nontarget_spec`→`aci_conditional_spec` 2,
`nontarget(`→`aci_conditional(` 2, `gaussian_kl_path`→`aci_metric` 1,
`gaussian_kl`→`aci_metric_pair` 1,
`cgns_from_affine`→`aci_model_from_affine` 1. The rules' `(?<![\w.])` guard
means the compiled internals that kept their names -
`.gaussian_kl_path_compiled`, `.forward_cir_compiled`,
`.da_filter_authenticated` - are never touched, and the quoted class label
`"cgns_model"` survives every occurrence.

**2. The file-name map**, applied to the patch's target paths and to any
in-body file reference: `causal_metrics.R`→`aci-core.R`,
`assimilation.R`→`aci-assimilation.R`, `model_classes.R`→`aci-model.R`,
`benchmark_models.R`→`aci-model-library.R`, and the rest of the twelve.
**0 in-body substitutions** for this family: the only file references the
reserve sources carry are their `## Origin:` provenance headers, which name
the frozen `aci` 0.0.30 tree rather than a sibling acir file, and which
re-application strips before insertion. Test file names were not renamed by
the mainline and are unchanged here.

**3. The argument map**, applied to the family code being inserted *and* to
the mainline anchors, so that the anchors match the current tree and the
inserted bodies call current names. 42 substitutions:
`nontarget`→`conditional` 41, `$nontarget`→`$conditional` 1. No `eps`
substitution: the range surface is FBCIR-scope, and this family never spells
that argument. The mainline's own lookbehind guards are reused verbatim, so
the internals (`.resolve_nontarget`), the condition classes
(`aci_error_nontarget[_crossnoise]`) and the model metadata
(`meta$estimand_nontarget`, `meta$nontarget_reduction`) are all out of reach.
Checked after the fact: the produced patch adds no line naming `nontarget`, a
bare `eps`, or an old file path.

The maps are applied surface-first, then arguments: otherwise the argument
rule would turn `nontarget(` into `conditional(` instead of letting the
surface rule turn it into `aci_conditional(`.

## Adaptations beyond the mechanical maps

Everything in this section is a change to something other than a name, and
each one is marked at its point of change during re-application.

- **The two new `R/` files are named to the current convention.** The reserve
  carries them as `ensemble.R` and `discovery.R`, the `aci` 0.0.30 spelling.
  Every other `R/` source in acir is now `aci-*`, so they enter as
  `R/aci-ensemble.R` and `R/aci-discovery.R`. **This is a judgement call, not
  a mechanical map**: the mainline file rename covered only files that already
  existed, and these two do not. It is safe for the reason the mainline rename
  recorded - `DESCRIPTION` has no `Collate` field, no source carries
  `@include`, and the package defines no S4 class or load-time cross-file
  dependency, so `R/` load order does not matter. The file **contents** are
  untouched by this choice.
- **Three anchors that moved with the surface rename, re-pointed.** The two
  not-implemented stubs are `aci_filter.stochastic_model` and
  `aci_smoother.stochastic_model` under `@describeIn`, so the removal anchor is
  spelled the new way; the ensemble guard that sat at the head of
  `forward_cir.aci_result` now goes in after `.aci_range_direction(direction)`
  in `aci_range.aci_result`; and the model-library body swap anchors on
  `aci_dyad_model <- function(variant = "p1", ...)` with a second guard that
  aborts if the renamed reserve body does not begin `aci_dyad_model <-
  function(`, so a rename miss cannot write a file defining the retired
  constructor.
- **One call inside the `aci()` splice.** Its fallback metric call is
  `aci_metric(smoo, filt, decompose = decompose)`, the current name for
  `gaussian_kl_path()`. The compiled CGNS branch still calls
  `.gaussian_kl_path_compiled()`, which kept its name.
- **Line width is left as the reserve wrote it.** The reserve bodies carry
  132 added `R/` lines over 80 columns from the 0.0.30 excision. **Exactly one
  of them crossed 80 because of the argument map**: a continuation line in the
  ensemble engine reading `inflation = inflation, seed = seed, conditional =
  conditional,`, 78 columns before and 82 after. It is left verbatim. The
  mainline would have re-wrapped it; re-wrapping a family body here is a body
  edit, and the point of this directory is that the bodies are not edited.
  Recorded so it is a known deviation rather than an unnoticed one.

Three adaptations from earlier regenerations are kept as they were:

1. **The reserve is read from outside the package tree.**
2. **The ensemble arm is spliced into the current `aci()`, not swapped in.**
   Seven targeted edits - the roxygen model/engine line, the five ensemble
   `@param` blocks, the signature, the `engine = "auto"` resolution, the
   ensemble arm plus the engine-conditional metric path, the
   engine-conditional lag table and recorded init, and `meta$m` with the
   smoother-scheme fallback. The mainline's `regularize` policy, `.aci_reg_new`
   / `.aci_reg_freeze`, `loglik` switch and reused-table `scheme` fallback all
   survive. **The ensemble code inside the splice is aci 0.0.30 verbatim apart
   from the counted maps.**
3. **`meta$smoother_scheme`** reads
   `smoo$meta$scheme %||% if (engine == "ensemble") "enkbs" else "unspecified"`,
   because the mainline renamed the CGNS smoother's tag from `route` to
   `scheme` and the ensemble smoother (`as_gaussian`) carries neither.

**`aci_dyad_model` was checked, not assumed.** The reserve copy
(`code/model_dyad_p3_p4_and_observe_y.R`) is a strict superset of the current
mainline body - it already carries the locked coefficient environment,
`.attach_cgns_realizer` and `source_status` - so the wholesale body swap loses
nothing. Note that this makes the file's own "Verbatim excision ... not
modified" header inaccurate: it carries post-0.0.30 material. Re-application
guards the swap with a four-symbol check that aborts rather than writing if a
future reserve copy falls behind the mainline.

**S16 ADAPTATION E.** The session-16 mainline guards
`aci()` against a literal `m =` in the call: on the mainline, `m`
partial-matches `model` and silently shifts every argument, so the guard
raises `aci_error_not_implemented` naming this family's deferred
ensemble-size argument. This family restores `m` as a real argument, so
re-application deletes the guard and its two `test-03` liveness assertions -
staged-absence scaffolding of exactly the same class as the out-of-scope
engine error the family replaces.

## Verification status

- **Applies:** `git apply --check --directory=acir dev/reserve/enkbs/enkbs.patch`
  from the ACI-project root, against the session-16 update -
  clean.
- **Installs:** `R CMD INSTALL` into a separate preview library:
  `* DONE (acir)`, no warnings.
- **Tests (session 16; mainline suite is 6410 / 0 / 0 / 1):** run in
  isolation, `test-06-ensemble.R`, `test-08-discovery.R` and
  `test-14-golden-enkbs.R` with `helper-golden-p3.R` give **247 / 0 / 0 / 0**,
  the whole family delta, G3 and G4 included - identical to every previous
  regeneration, so the session-16 mainline changes cost the family nothing.
  The full suite on branch `preview/enkbs` reads **6648 pass, 3 fail, 1
  error, 1 skip** (the same `ACI_ORACLE_PARITY_ROOT` skip). All four deltas
  are mainline evidence gates shipped after the previous regeneration,
  asserting the mainline state this patch legitimately changes:
  `test-30-evidence-register.R:96` and `test-31-gate-liveness.R:80/:82`
  (register coverage - the patch adds exports `gaspari_cohn`,
  `constrained_mle`, `stochastic_model`, `ensemble_lag_table`, `as_gaussian`
  with no register rows yet) and `test-31-gate-liveness.R:191` (the mainline's
  non-CGNS refusal, which the ensemble engine exists to lift). Integrating
  the family closes all four by adding the register rows and updating those
  gate tests; they are recorded here, not worked around.
- **Evidence grade - mixed, and this is the important line.**
  - `test-14-golden-enkbs.R` (G3/G4) is a genuine **golden grade against the
    published EnKBS dyad experiment**: `helper-golden-p3.R` is a faithful R
    port of `EnKBS-main/dyad/utov.m`, and G3/G4 grade EnKBF, EnKBS and the
    ensemble lag table / forward range against it. That evidence travels with
    this patch.
  - **Breadth beyond the dyad remains open.** The Lorenz-96 (`model_l96`),
    `v -> u` and Lorenz-84 cases have behavioural tests only; there is no
    authors-source fixture for them in this repository. Closing that breadth
    is later work, recorded rather than glossed. This patch ships
    implementation plus existing tests, not completed verification.
  - The discovery source's coverage is likewise behavioural (`test-08`), not
    graded against EnKBS-main output.

## Stacking with fbcir.patch

They do not stack - re-measured on the session-16 update, in both orders. Each applies
cleanly to a clean `main`; neither applies on top of the other:

| order | plain `git apply --check` refuses on | `--3way --check` conflicts on |
|---|---|---|
| enkbs, then fbcir | `NAMESPACE`, `R/aci-core.R` | `NAMESPACE`, `R/aci-core.R` |
| fbcir, then enkbs | `NAMESPACE`, `R/aci-core.R` | `NAMESPACE`, `R/aci-core.R` |

The overlap is smaller than it was: `R/aci-assimilation.R` now merges cleanly
under `--3way` in both orders, where the equivalent file conflicted before.
`NAMESPACE` and `R/aci-core.R` remain genuinely shared. Applying both still
means a manual merge, or regenerating one patch against a branch that already
carries the other. Recorded, not worked around.

## Regeneration rule

Patches rot as the mainline moves, and the last four regenerations were each
forced by a mainline rename rather than by a change of algorithm. The standing
rule this has settled into:

**Regenerate both patches after the last commit that touches the package
surface, immediately before any push. `git apply --check` against the current
`main` is the only test of a patch's validity - a patch that was clean last
week is evidence of nothing.**

Regeneration is maintainer-side work; the patch in this directory is the
deliverable. The shape of it:

```sh
git checkout -b preview/enkbs main
# re-apply this directory's sources on top, per the maps and adaptations above
R CMD INSTALL . --library=<preview lib> && <run test_dir>
git commit -am "preview/enkbs: ..." && git checkout main
git diff main..preview/enkbs > <this directory>/enkbs.patch
git apply --check <this directory>/enkbs.patch
git branch -D preview/enkbs      # record the tip above first
```

**The acceptance test for a candidate patch, end to end.** No scripts and no
external state are needed:

```sh
git apply --check <this directory>/enkbs.patch   # must be silent
git checkout -b check/enkbs main
git apply <this directory>/enkbs.patch
R CMD INSTALL . --library=<a scratch library>
Rscript -e '.libPaths(c("<a scratch library>", .libPaths()));
            library(testthat); library(acir);
            print(test_dir("tests/testthat", package = "acir",
                           filter = "^(06|08|14)-",
                           reporter = "summary"))'
```

Expect `* DONE (acir)` from the install and **247 pass, 0 fail, 0 error,
0 skip** from those three files. Run the same command without `filter =` for
the whole suite; expect **6483 / 0 / 0 / 1 skip**, the skip being the optional
`ACI_ORACLE_PARITY_ROOT` MATLAB oracle tree.

**Start from the current adaptations, not the earlier ones.** Every earlier
re-application is stale, and each in a different way: the extraction-time one
read the reserve from inside the package and swapped `aci()` wholesale; the
pre-rename one has dead anchors and inserts bodies that call retired names;
the post-surface-rename one writes to `R/` paths that no longer exist.

Re-application is anchor-based and fails loudly on a moved anchor. **The
largest risk here is still `aci()`, but it is no longer an anchor risk.** The
whole-block swap was replaced with seven targeted edits, so a mainline change
to any one of those seven sites fails loudly and locally instead of silently
reverting the rest of `aci()`. Re-merge only the site that moved. The same
caution applies to the model library: `aci_dyad_model`'s body is still a
whole-body swap, guarded but not anchored, so diff the reserve copy against
the mainline body before trusting it. **A rename is a third failure mode**,
invisible to both: an anchor can still match while the body it inserts calls a
name the mainline no longer exports. **And a file rename is a fourth**, which
is louder than all of them - the patch does not apply at all, because its
headers name files that are gone. The counted maps above exist so that the
third class shows up as a number and the fourth as a mapped path, not as a
load error or a dead patch.
