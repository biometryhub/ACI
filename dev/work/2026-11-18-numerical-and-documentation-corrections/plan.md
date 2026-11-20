# Implementation and evidence plan

1. Preserve the base source as a separate snapshot. Run the baseline suite
   with the saved oracle enabled.
2. Restrict the new production correction to the scalar KL kernel. Test safe
   evaluation orders and near-overflow classification before selecting one.
   Add focused public-route regression tests and a short help explanation.
3. Include the core/lag/online/range documentation, generated help, plotting
   fixes, vignette window count and regression tests. Add NEWS and these work
   records. Do not include the optional near-equality branch.
4. Run exact-input accuracy, unchanged-arithmetic, complete-suite and saved
   source comparisons. Run timing separately from other heavy computations.
   Run fixture provenance and an isolated local package build/as-cran check
   using existing dependencies only.

The main risks are spurious infinity close to the output range limit and
extra work on the ordinary numerical path. The ordinary path's arithmetic must
not change. The optional series touches the same kernel and must be ported
onto this fix if it is accepted later.

Temporary installation by R CMD check is confined to its isolated hygiene
directory. No dependency was installed and no reference code was run.
Cross-platform CI and missing local tools remain explicit limits.
