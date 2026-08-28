# reserve/fbcir — verification status and regeneration rule

Family: `FBCIR_code-main/`, paper `andreou2026cir`. Everything here was
excised from `aci` 0.0.30 (git tree `97f6b124`) when acir 0.1.0 was cut.

## Regeneration log

| regenerated | against `main` | branch tip | patch |
|---|---|---|---|
| 2026-08-28, current, post-file-and-argument rename | `c7b8720` | `8bee983` | 1031 lines |
| 2026-08-28, post-surface-rename | `3db97dd` | `3bb4f05` | 1031 lines |
| 2026-08-28, pre-rename | `bded686` | `8f75f9f` | 958 lines |
| at the 0.1.0 extraction | `8fba4bc` | `31035a3` | 941 lines |

**Three mainline commits killed the previous patch, and the failure was total
rather than partial.** Measured with `git apply --check` on `c7b8720` before
regenerating:

```
error: R/assimilation.R: No such file or directory
error: R/causal_metrics.R: No such file or directory
error: R/compiled_cir.R: No such file or directory
error: R/compiled_lag.R: No such file or directory
error: patch failed: tests/testthat/test-05-nontarget.R:35
```

The mainline renamed twelve of its thirteen `R/` sources to the `aci-*`
convention, so **every target path in the patch named a file that no longer
exists**; a second commit renamed `nontarget =` to `conditional =` and `eps =`
to `epsilon =`, which broke the surviving test-file hunks on their context
lines; a third was metadata only. Note the difference from the previous
rotting: a *surface* rename fails on context lines, and a reader might be
tempted to force it. A *file* rename fails on the header, and there is nothing
to force. Branch tips are recorded so the branches are recreatable; the
branches themselves are deleted after verification, and the patch reconstructs
them exactly.

## The framing change the surface rename forced

The rename did not only turn `forward_cir()` into `aci_range()`. It gave
`aci_range()` a `direction = c("forward", "backward")` argument and a release
gate, `.aci_range_direction()`, that raised `aci_error_not_implemented` on
`"backward"` and named this reserve directory in the message. The backward
range therefore no longer re-enters as a separate exported `backward_cir()`
generic: **it re-enters as `aci_range(direction = "backward")`**, which is the
surface the rename already promised, and the gate is what the patch opens.

`backward_cir` is not a name in the patched package any more. The two family
reducers keep their bodies verbatim and become internals,
`.aci_range_backward_table()` and `.aci_range_backward_result()`; the
convenience method that dispatched straight from a model becomes the exported
`aci_range.cgns_model()`. `cir_pair()` and `lt_onelag()` keep their names.

## What `fbcir.patch` restores

Generated `git diff main..preview/fbcir`, 1031 lines, over ten files.

| area | restored |
|---|---|
| `R/aci-core.R` | `.bwd_lengths`; `.aci_range_backward_table` and `.aci_range_backward_result` (the reserve's `backward_cir.lag_table` / `backward_cir.aci_result` bodies); `aci_range.cgns_model`; `cir_pair`; the backward branch of `aci_range.lag_table` and `aci_range.aci_result`; `.aci_range_direction` reduced from a gate to a resolver |
| `R/aci-assimilation.R` | `lt_onelag`; the `is.null(x$rows)` arm of `as.data.frame.lag_table`; `mode = "one_lag"` in `lag_table` |
| `R/aci-online-smoother.R` | the `mode == "one_lag"` early return in `.lagtable_core_compiled`; `mode = "one_lag"` in `.lag_table_compiled` |
| `R/aci-cir.R` | `.slice_compiled_cgns`, `.slice_compiled_filter` |
| `NAMESPACE` | `S3method(aci_range, cgns_model)`, `export(cir_pair)`, `export(lt_onelag)` |
| tests | the backward blocks spliced into the current test-04-cir.R and test-20-compiled-cir.R (**not** those files reverted to their `aci` form; see the retained earlier adaptations below); the T5 one-lag and one-lag-propagation blocks in test-03-engine.R; the `one_lag` rows in test-18-compiled-lag.R; the `one_lag` route in test-05-nontarget.R |

**Not restored, deliberately:** `cir_table`. It carries its own provenance
disclaimer in the `aci` 0.0.30 source and is a separate disposition (dropped,
`reserve/extensions/code/cir_table.R`). The three assertions at the end of
`test-04-cir.R` that used it are replaced by a comment saying so. Restoring
backward CIR *with* the cir_table read-out means applying this patch and then
re-adding `reserve/extensions/code/cir_table.R` by hand.

Also not restored: any `man/` topic. The patch does not touch `man/`, so the
branch installs but its help index is short of `cir_pair` and `lt_onelag`, and
`aci_range.Rd` still carries the not-in-this-release paragraph.
`roxygen2::roxygenise()` over the restored roxygen regenerates all of it; the
acir NAMESPACE is hand-written, so nothing else depends on that run.

`code/model_tipping_triad.R`, `code/model_multiscale_fbcir.R` and
`code/model_topographic_layered_fbcir.R` are **not** in `fbcir.patch`; they
are the constructor half of the family and come back with
`tests/test-13-fbcir-models.R`, which the patch also does not carry. That is a
scope choice, not a failure: the patch covers the part of the family that sits
inside files acir still ships.

## The three mechanical maps

The reserve `.R` files are untouched: they stay the verbatim `aci` 0.0.30
excision. Three renames the mainline has since made are applied **on the way
in**, using the same rules and the same guards the mainline passes used, and
every substitution is counted so that drift shows up as a number rather than
as a load error.

**1. The public-surface map.** 18 substitutions across six reserve files:
`model_dyad`→`aci_dyad_model` 7, `forward_cir`→`aci_range` 3,
`cgns_model(`→`aci_model(` 3, `da_filter`→`aci_filter` 2,
`gaussian_kl`→`aci_metric_pair` 1, `nontarget_spec`→`aci_conditional_spec` 1,
`nontarget(`→`aci_conditional(` 1. The `(?<![\w.])` guard means the compiled
internals that kept their names are never touched.

**2. The file-name map**, applied to the patch's target paths and to any
in-body file reference: `causal_metrics.R`→`aci-core.R`,
`assimilation.R`→`aci-assimilation.R`, `compiled_lag.R`→
`aci-online-smoother.R`, `compiled_cir.R`→`aci-cir.R`, and the rest of the
twelve. **0 in-body substitutions** for this family: the only file references
the reserve sources carry are their `## Origin:` provenance headers, which
name the frozen `aci` 0.0.30 tree rather than a sibling acir file, and which
re-application strips before insertion. Test file names were not renamed by
the mainline and are unchanged here.

**3. The argument map**, applied to the family code being inserted *and* to
the mainline anchors, so that the anchors match the current tree and the
inserted bodies call current names. 72 substitutions:
`nontarget`→`conditional` 10, `$nontarget`→`$conditional` 6,
`table_nontarget`→`table_conditional` 1, `eps_grid`→`epsilon_grid` 2,
`eps`→`epsilon` 53. The mainline's own lookbehind guards are reused verbatim,
so the internals (`.resolve_nontarget`, `.cir_eps_grid`), the condition
classes (`aci_error_nontarget[_crossnoise]`), the model metadata
(`meta$estimand_nontarget`, `meta$nontarget_reduction`) and the **values**
`"matlab_eps_grid"` and `"eps_grid_objective"` are all out of reach. Checked
after the fact: the produced patch adds no line naming `nontarget`, a bare
`eps`, or an old file path.

The maps are applied surface-first, then arguments: otherwise the argument
rule would turn `nontarget(` into `conditional(` instead of letting the
surface rule turn it into `aci_conditional(`.

## Adaptations beyond the mechanical maps

Everything in this section is a change to something other than a name, and
each one is marked at its point of change during re-application.

1. **`.aci_range_direction()` becomes a resolver.** The release gate's abort
   on `"backward"` is removed and both call sites now bind its value instead
   of discarding it. This is the gate the family exists to open.
2. **`aci_range()`'s roxygen** loses the not-in-this-release paragraph and
   gains the reserve's backward-range prose, plus a backward `@examples` line.
   **One sentence of that prose is corrected.** The reserve says the exact
   form "integrates its suffix minima with composite Simpson"; the change that
   fixed the forward CIR objective falsified that when it made `.calc_tau`'s
   non-`l1_linf` branch the exact layer-cake sum. An earlier regeneration
   corrected the test for that change (adaptation 4 below); the sentence is
   corrected for the same reason. The family code is still not touched.
3. **`aci_range.lag_table()` routes the backward direction** to
   `.aci_range_backward_table()` before its own `match.arg` block. The two
   directions are different functionals with different parameter sets, so
   they are not folded into one body: backward defaults to the eq. G.14 L1
   sum rather than Simpson, and has no `anchors`, `convention` or
   `epsilon_grid`. Supplying any of those three with `direction = "backward"`,
   or `quadrature = "matlab_eps_grid"`, is now an `aci_error_dims`; the
   reserve had no equivalent, because its backward entry point had no such
   formals. `aci_range.aci_result()` gets the same two-line route to
   `.aci_range_backward_result()`.
4. **The generic is retired, the bodies are not.** `backward_cir.lag_table`
   and `backward_cir.aci_result` are spliced in by slicing the reserve file at
   its signature lines, so their bodies arrive verbatim. Four things move:
   the two definition lines; the signature's continuation indent (the new name
   is three characters longer, and the one line that would then pass 80
   columns is wrapped); `"Unused arguments were supplied to backward_cir()."`
   → `... aci_range()."`; and the reference-count warning text, which said
   `"backward_cir over >20 reference times ..."` and now says
   `"A backward range over >20 reference times ..."`. The internal recursive
   call `backward_cir(tb, ...)` becomes a direct
   `.aci_range_backward_table(tb, ...)` — the family's own S3 dispatch is not
   preserved as an internal generic, because with exactly two reducers and a
   known argument class it bought nothing and would have needed its own
   registration. The reducers' `@param` documentation is condensed into
   `@noRd` blocks; the substantive prose moved to `aci_range()` (adaptation 2)
   and the backward `quadrature` semantics into `aci_range.lag_table()`'s
   `@param quadrature`.
5. **`backward_cir.cgns_model` becomes `aci_range.cgns_model`**, an exported
   `@describeIn` method. Its body is the reserve's. One behaviour is new:
   `aci_range(model, direction = "forward", ...)` raises
   `aci_error_not_implemented` naming the lag-table and result routes, where
   aci 0.0.30 gave R's "no applicable method" — there was no
   `forward_cir.cgns_model` there to dispatch to. **`cir_pair()`** keeps its
   name and body; only its two range calls move to the `aci_range` spelling,
   and the `structure()` call is re-wrapped over three lines to stay inside 80
   columns.
6. **Every `backward_cir(...)` call site in the reserve test blocks** becomes
   an `aci_range(..., direction = "backward")` call site — 12 in
   `tests/test-04-backward-cir.R`, 4 in `tests/test-20-backward-cir.R`. They
   are written out in full, with the continuation lines realigned, rather than
   regex-substituted, so the re-wrapping is visible at the point of change.
   Two of them also become keyword calls: `backward_cir(model, obs, T = ...)`
   passed `obs` positionally, which still resolves after the named
   `direction`, but is written `obs = obs` rather than left to argument-order
   reasoning.
7. **A duplicated one-lag propagation block is dropped.** An earlier
   regeneration inserted `test_that("one-lag propagation does not truncate on
   an unscaled operator norm")` twice: once from
   `tests/test-03-one-lag.R`, which carries both staged blocks, and again by
   lifting the same block out of `aci/tests/testthat/test-03-engine.R`. The two
   texts were checked and are byte-identical, so the second insertion was a
   duplicate; its 3 assertions are the whole difference between the 6400 of
   that regeneration and the 6397 here. Dropping it also removed the last read
   of the frozen `aci` 0.0.30 tree, which neither family's re-application now
   needs.
8. **`lt_onelag()`'s `@seealso`** named the retired `backward_cir()` topic and
   now names `aci_range()`.
9. **Line width is left as the reserve wrote it.** The reserve bodies carry
   13 added `R/` lines over 80 columns from the 0.0.30 excision. The argument
   map lengthened three of them further; **none of the three was inside 80
   columns before the map**, so no line crossed the boundary here, and none
   was re-wrapped. Re-wrapping a family body is a body edit, and the point of
   this directory is that the bodies are not edited.

Four adaptations from earlier regenerations are kept as they were: the reserve
is read from outside the package tree; `test-04-cir.R` and
`test-20-compiled-cir.R` are edited in place, never overwritten with `aci`'s
copies; the backward brute-force fragment is restored as its own block
replaying aci 0.0.30's seed-11 draw order; and the terminal +M cell
expectation stays at `0.4 * dt` with the closed-form assertion at
`tolerance = 0`.

## Verification status

- **Applies:** `git apply --check reserve/fbcir/fbcir.patch` on a clean `main`
  (`c7b8720`) — clean.
- **Installs:** `R CMD INSTALL` into a separate preview library:
  `* DONE (acir)`, no warnings.
- **Tests:** `testthat::test_dir()` on branch `preview/fbcir` —
  **pass = 6397, fail = 0, error = 0, warning = 0, skip = 1**. The skip is the
  optional external MATLAB oracle tree gated on `ACI_ORACLE_PARITY_ROOT`.
  Mainline `main` is pass = 6236, so the patch reinstates **161 assertions**,
  in `test-03`, `test-04`, `test-05`, `test-18` and `test-20`. Those five
  files run in isolation give **537 / 0 / 0 / 0**. Both figures are identical
  to the previous regeneration's, so the file and argument renames cost the
  family nothing.
- **Evidence grade of what is restored: behavioural and analytic only.** The
  backward-CIR tests are identity checks (`.bwd_lengths` against a brute-force
  layer-cake grid, the terminal-sentinel convention, the two quadrature
  conventions, prefix-sharing agreement between shared and separately built
  tables) and package-to-package consistency. **There is no authors-source
  FBCIR fixture in this repository.** `fixtures/FBCIR-MULTISCALE.md` and
  `fixtures/generate_fbcir_multiscale_oracles.m` are the unrun generator and
  its status note; the generator's own header states that the bytes it would
  produce are authors-source *execution* fixtures, not author-provided
  outputs, and that authors-source numerical parity remains open. **Running
  it needs MATLAB, which was not available during construction.**
  Author-output parity for the FBCIR family is therefore an open 0.2.0/0.3.0
  item, not something this patch closes.

## Stacking with enkbs.patch

`fbcir.patch` and `enkbs.patch` each apply cleanly to a clean `main`, but
**they do not stack** — re-measured on `c7b8720`, in both orders:

| order | plain `git apply --check` refuses on | `--3way --check` conflicts on |
|---|---|---|
| fbcir, then enkbs | `NAMESPACE`, `R/aci-core.R` | `NAMESPACE`, `R/aci-core.R` |
| enkbs, then fbcir | `NAMESPACE`, `R/aci-core.R` | `NAMESPACE`, `R/aci-core.R` |

The overlap is smaller than it was: `R/aci-assimilation.R` now merges cleanly
under `--3way` in both orders, where the equivalent file conflicted before.
`NAMESPACE` and `R/aci-core.R` remain genuinely shared. Applying both still
means a manual merge, or regenerating one patch against a branch that already
carries the other. Recorded, not worked around.

Note the one real cross-dependency in the other direction:
`reserve/enkbs/tests/test-06-ensemble.R` asserts the ensemble guard on the
backward range, so `enkbs.patch` comments that line out and names this patch.
Under the renamed surface it has a second reason to stay out: the release
gate answers `direction = "backward"` before the ensemble guard can, on every
result, so the assertion would pass without testing the ensemble at all.

## Regeneration rule

Patches rot as the mainline moves, and the last four regenerations were each
forced by a mainline rename rather than by a change of algorithm. The standing
rule this has settled into:

**Regenerate both patches after the last commit that touches the package
surface, immediately before any push. `git apply --check` against the current
`main` is the only test of a patch's validity — a patch that was clean last
week is evidence of nothing.**

Regeneration is maintainer-side work; the patch in this directory is the
deliverable. The shape of it:

```sh
git checkout -b preview/fbcir main
# re-apply this directory's sources on top, per the maps and adaptations above
R CMD INSTALL . --library=<preview lib> && <run test_dir>
git commit -am "preview/fbcir: ..." && git checkout main
git diff main..preview/fbcir > <this directory>/fbcir.patch
git apply --check <this directory>/fbcir.patch
git branch -D preview/fbcir      # record the tip above first
```

**The acceptance test for a candidate patch, end to end.** No scripts and no
external state are needed:

```sh
git apply --check <this directory>/fbcir.patch   # must be silent
git checkout -b check/fbcir main
git apply <this directory>/fbcir.patch
R CMD INSTALL . --library=<a scratch library>
Rscript -e '.libPaths(c("<a scratch library>", .libPaths()));
            library(testthat); library(acir);
            print(test_dir("tests/testthat", package = "acir",
                           filter = "^(03|04|05|18|20)-",
                           reporter = "summary"))'
```

Expect `* DONE (acir)` from the install and **537 pass, 0 fail, 0 error,
0 skip** from those five files. Run the same command without `filter =` for
the whole suite; expect **6397 / 0 / 0 / 1 skip**, the skip being the optional
`ACI_ORACLE_PARITY_ROOT` MATLAB oracle tree.

**Start from the current adaptations, not the earlier ones.** Every earlier
re-application is stale, and each in a different way: the extraction-time one
read the reserve from inside the package and copied whole test files; the
pre-rename one has a dead `forward_cir.aci_result` anchor, inserts bodies that
call retired names, and carries the duplicated one-lag block; the
post-surface-rename one writes to `R/` paths that no longer exist.

Re-application is anchor-based, so it fails loudly on a moved anchor rather
than producing a silent half-patch. If an anchor is gone, the mainline changed
under the reserve: fix the anchor, do not force the patch. **Anchors are not
the only failure mode.** Anywhere a whole `aci` file or function body is
copied, diff it against the current acir one first: on the extraction-time
tree those copies were no-ops, and on any later tree they silently revert
mainline work. Three such places were found and converted to in-place edits.
**A rename is a third failure mode**, invisible to both: an anchor can still
match while the body it inserts calls a name the mainline no longer exports.
**And a file rename is a fourth**, which is louder than all of them — the
patch does not apply at all, because its headers name files that are gone.
The counted maps above exist so that the third class shows up as a number and
the fourth as a mapped path, not as a load error or a dead patch.
