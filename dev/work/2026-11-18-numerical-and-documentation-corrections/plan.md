# Implementation and evidence plan

The author approved adoption and staging on 18 November 2026. The external
PLAN.md and before/ hashes were saved before implementation. No commit is
authorised, so the commit-before-code convention is represented by those
saved plans and source identities, as in the 17 November retained-table work.

1. Preserve the maintained checkout and approved default source as separate
   snapshots. Run the approved baseline suite with the saved oracle enabled.
2. Restrict the new production correction to the scalar KL kernel. Test safe
   evaluation orders and near-overflow classification before selecting one.
   Add focused public-route regression tests and a short help explanation.
3. Carry forward the approved core/lag/online/range documentation, generated
   help, plotting fixes, vignette window count and regression tests. Add NEWS
   and these work records. Do not include the optional near-equality branch.
4. Run exact-input accuracy, unchanged-arithmetic, complete-suite and saved
   source comparisons. Run timing separately from other heavy computations.
   Run fixture provenance and an isolated local package build/as-cran check
   using existing dependencies only.
5. Present the complete maintained-to-final change and the incremental signal
   correction in the existing inline colour review. Preserve stable item IDs
   and verify exact source reconstruction.
6. Verify the maintained source and index still match the saved baseline.
   Apply tested files, stage explicit paths and compare index bytes with the
   tested candidate. Do not commit or push. Record the deferred alternative
   and the remaining verification limits in the workspace handover.

The main risks are spurious infinity close to the output range limit and
extra work on the ordinary numerical path. Retain failed prototype evidence.
The ordinary path's arithmetic must not change. The optional series touches
the same kernel and must be ported onto this fix if subsequently accepted,
rather than restoring its older complete kernel over the signal correction.

Evidence lives in the external output named in intent.md, including the full
report, source hashes, commands, results and exact colour review. Temporary
installation by R CMD check is confined to its isolated hygiene directory.
No dependency installation, reference-code execution or remote CI run is
authorised. Cross-platform CI and missing local tools remain explicit limits.
