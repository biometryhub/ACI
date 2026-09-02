# `reserve/` - what acir 0.1.0 left out, and how to get it back

This directory is the other half of the acir 0.1.0 extraction. The mainline is
the closed-form `ACI_code-main/` surface only; everything that was cut is
filed here by **origin**, with a one-line audit entry in
[`DISPOSITIONS.md`](DISPOSITIONS.md).

`reserve/` lives outside the package tree entirely (`dev/reserve/`, beside
`acir/`, since the repo-hygiene commit `3106967`), so it costs the built
package nothing; no build-ignore entry is involved.

Everything here came from **`aci` 0.0.30, git tree `97f6b124`**. Every file
carries a header naming its origin path, that tree, its category, its intended
release, and one line of reason.

---

## Categories

| directory | what lives here | intended release |
|---|---|---|
| `fbcir/` | `FBCIR_code-main`-derived: the backward-CIR family, the FBCIR constructors, the unrun fixture generator | 0.2.0 or 0.3.0 |
| `fbcir-paper/` | FBCIR-paper-backed with no MATLAB | *empty by finding* - see its README |
| `enkbs/` | `EnKBS-main`-derived: the ensemble engine, localization, discovery, the L84/L96 constructors | 0.2.0 or 0.3.0 |
| `enkbs-paper/` | EnKBS-paper-backed with no MATLAB | *empty by finding* - see its README |
| `aci_code-future/` | ACI_code-derived but deferred from 0.1.0 | *category retired*: it was empty of code by finding, and the `T_C`-hidden gap it recorded is closed, so no directory ships |
| `aci-paper/` | backed by an ACI-adjacent paper with no ACI_code MATLAB (`lt_contraction_certificate`, the `cfy22` ENSO variant) | 0.1.x, TBD |
| `paper-extremes/` | `moser2026extremes` (fourth paper, no MATLAB) | after the three MATLAB-backed families |
| `extensions/` | package-only infrastructure with neither paper nor MATLAB backing, plus the five outright drops | unscheduled |

The release order for FBCIR vs EnKBS is not settled here.

Within a family: `code/` holds R sources, `tests/` holds test files and
excised test blocks, `fixtures/` holds non-R assets, `NOTES.md` records
verification status, and `<family>.patch` is the re-application patch.

## How the patches work

`fbcir.patch` and `enkbs.patch` are `git diff main..preview/<family>` against
the acir mainline. Each was built on a real branch, installed, and run through
the full test suite before the diff was taken.

In this repository the package is the `acir/` subdirectory, so from the
ACI-project root (the patches are package-root-relative):

```sh
git checkout -b restore/fbcir main
git apply --directory=acir dev/reserve/fbcir/fbcir.patch
R CMD INSTALL acir --library=<a scratch library>
```

Both `--directory=acir` apply-checks were run and returned clean at the
session-16 update (this branch's single commit after the `acir-import` tag).

Both were last verified against the session-16 update; the
figures below are from that verification, and each family's `NOTES.md`
carries the full log.

| | applies to clean mainline | installs | family isolation | full suite on the preview branch |
|---|---|---|---|---|
| `fbcir.patch` (1031 lines, regenerated at `c7b8720`) | clean at the session-16 update | `* DONE (acir)` | **539 / 0 / 0 / 0** | 6562 pass, 3 fail, 1 error, 1 skip |
| `enkbs.patch` (2636 lines, regenerated at the session-16 update) | clean at the session-16 update | `* DONE (acir)` | **247 / 0 / 0 / 0** | 6648 pass, 3 fail, 1 error, 1 skip |

Mainline is 6410 pass / 0 fail. The full-suite deltas are not family defects:
they are `inst/evidence/` register-coverage and staged-absence gates shipped
after the last patch regeneration, asserting the MAINLINE state - a family
patch legitimately adds exports the register does not yet row, and removes
refusals those gates assert. Integrating a family closes them by adding
register rows and updating the named gate tests; each family's `NOTES.md`
lists them test by test. The single skip in every run is the optional
external MATLAB oracle tree gated on `ACI_ORACLE_PARITY_ROOT`.

**Three limits, on the record:**

1. **The two patches do not stack.** Both edit `NAMESPACE` and `R/aci-core.R`;
   `git apply`, plain and `--3way`, refuses the second in either order.
   Applying both means a manual merge, or regenerating the second against a
   branch already carrying the first.
2. **Neither patch touches `man/`.** The preview branches install with a help
   index short of the restored topics. `roxygen2::roxygenise()` fixes that;
   the acir `NAMESPACE` is hand-written, so nothing else depends on it.
3. **A patch ships implementation plus existing tests, not completed
   verification.** `enkbs.patch` carries the machine-precision dyad golden
   grade (G3/G4); L96 / `v -> u` / L84 breadth is open. `fbcir.patch` carries
   behavioural and analytic evidence only - author-output parity needs the
   MATLAB generator run, and there is no MATLAB here. Each family's
   `NOTES.md` states this per item.

## How the patches are regenerated

**Patches rot as the mainline moves. Every acir release regenerates both** -
a stated release-checklist step, and the practical form of the
architectural-fit guardrail: keep the engine surface they touch stable, or
knowingly pay the rebase.

Regeneration is maintainer-side work and the patch is the deliverable here.
It is done by branching from the mainline, re-applying the reserve sources,
installing, running the suite, and re-taking the diff; the re-application is
**anchor-based**, meaning every edit asserts that its anchor text occurs
exactly once and aborts loudly otherwise. A failed anchor means the mainline
moved under the reserve: fix the anchor, never force the patch. Each family's
`NOTES.md` states the regeneration procedure and the failure modes to watch.

## Cross-references you will trip over

Three excised items are single arms of functions that also carry retained or
other-family code, so they cannot be filed twice. Each is filed once with its
primary family, and a `POINTER.md` sits in the other category:

- `model_dyad`'s `p3` (EnKBS), `p4` (extremes) and `observe = "y"` (extension)
  arms - the whole body is at `reserve/enkbs/code/model_dyad_p3_p4_and_observe_y.R`.
- `model_l84`'s `fbcir` variant - the whole constructor is at
  `reserve/enkbs/code/model_l84.R` (its *default* variant is EnKBS).

Two excised **test blocks** assert three families at once
(`test-02-models.R:50-58`, `test-11-model-parity.R:21-155`); both are filed
under `reserve/enkbs/tests/`, with `CROSS-FAMILY.md` notes in
`reserve/fbcir/tests/` and `reserve/paper-extremes/tests/`.

One genuine dependency crosses categories in the other direction:
`paper-extremes/code/extremes.R` needs the backward range (`reserve/fbcir/`)
and `projected_kl` (`reserve/extensions/`), and
`reserve/enkbs/tests/test-06-ensemble.R` asserts the backward range's ensemble
guard. `enkbs.patch` comments that one line out and names the patch it needs.
Note the surface: the reserve sources carry the `aci` 0.0.30 spelling
`backward_cir()`, but on the current mainline the backward direction re-enters
as `aci_range(direction = "backward")`, and that is what `fbcir.patch`
restores.

## Five things here are dropped, not reserved

`kl_increment` (dead), `cir_table` (self-disclaimed provenance),
`empirical_kl` (not-implemented stub), `truncation_profile` (no callers), and
the unreachable `mode == "smoother_only"` branch. Their text is kept under
`extensions/` **for the audit trail only** - no release is planned for them.
`DISPOSITIONS.md` §7 gives the reason for each.
