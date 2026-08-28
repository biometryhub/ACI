# reserve/aci-paper — ACI-paper-backed extensions, no ACI_code MATLAB

Blocks whose provenance is one of the ACI-adjacent papers rather than
`ACI_code-main/`. They are closed-form and CGNS-native, so they are the most
likely 0.1.x additions, but they are not part of the 0.1.0 fidelity claim.

- `code/lt_contraction_certificate.R` — `andreou2026smoother` eqs. 3.18-3.19.
- `code/model_enso6_cfy22_variant.R` — the `chen2022enso` transcription that
  is `aci`'s *default* `model_enso6` variant. The acir mainline keeps only
  `variant = "aci_code"`, so the default changed; this is a deliberate
  scope decision, recorded in DISPOSITIONS.md.
- `tests/` — the three excised test blocks that drive them.

No patch is prepared for this category: patches are scoped to the two
MATLAB-backed families.
