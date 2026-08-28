# reserve/extensions — package-only, no paper and no MATLAB backing

Infrastructure `aci` grew that is not traceable to any of the four papers or
any of the three MATLAB codebases. Nothing here is scheduled.

- `code/formula_interface.R` — the `aci_fit()` front-end and its twelve
  generics, `aci/R/formula_interface.R:1-610`. The three engine `plot`
  methods that were physically stranded in that file (`:613-674`) are **not**
  here: they are mainline acir, at `R/plots.R`.
- `code/applied_workflows.R`, `code/validation_diagnostics.R` — the applied
  on-ramps and the validation suite.
- `code/cir_table.R`, `code/kl_increment.R`, `code/empirical_kl.R`,
  `code/projected_kl.R`, `code/truncation_profile.R` — five exported symbols
  recorded as **dropped** in DISPOSITIONS.md. They are kept here as text so
  the audit can see exactly what left the surface, not because a release is
  planned for them.
- `code/lagtable_smoother_only_mode.diff` — the unreachable
  `mode == "smoother_only"` branch, recorded as a diff rather than as
  re-appliable code.
- `code/model_dyad_observe_y.POINTER.md` — the reverse-partition branch.
- `tests/` — the matching test files and blocks.
