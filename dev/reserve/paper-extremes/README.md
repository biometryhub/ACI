# reserve/paper-extremes — moser2026extremes family, no MATLAB

The fourth paper's extreme-event machinery. No supplied MATLAB codebase
corresponds to it, so under the release map it comes after all three
MATLAB-backed families.

- `code/extremes.R` — the whole family (`detect_events`, `event_influence`,
  `sensitive_directions`, `features_*`, `classify_events`).
- `code/model_pathways.R`, `code/model_topographic.R` — its two constructors.
  Note `model_topographic` is the moser2026extremes spectral barotropic
  model, **not** FBCIR's `model_topographic_layered_fbcir`; `aci` says so at
  `benchmark_models.R:820-824` and `:1079-1081`.
- `code/model_dyad_p4_arm.POINTER.md` — the `p4` preset.
- `tests/test-07-extremes.R`, `tests/CROSS-FAMILY.md`.

Dependency to carry when this family returns: `extremes.R` is the only
in-package consumer of `projected_kl` (`aci/R/extremes.R:498,503,518,581`)
and it also calls the backward range (`:301`, spelled `backward_cir()` in the
`aci` 0.0.30 source) and the engine internal `.cir_min_strength`.
`projected_kl` is filed in `reserve/extensions/`; the backward range in
`reserve/fbcir/`, where it re-enters as `aci_range(direction = "backward")`.
