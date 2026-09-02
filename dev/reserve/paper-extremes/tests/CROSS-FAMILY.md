# Cross-family test blocks that also cover the extremes family

Two excised blocks assert `model_pathways` / `model_topographic`
(moser2026extremes) alongside FBCIR and EnKBS constructors. Both are filed
once, under EnKBS:

- `reserve/enkbs/tests/test-02-zoo-models.R` - its `model_pathways()` and
  `model_topographic()` lines.
- `reserve/enkbs/tests/test-11-other-families.R` - its pathways and
  spectral-topographic parity blocks.

`projected_kl`, whose only in-package consumer is this family
(`aci/R/extremes.R:498,503,518,581`), is filed under `reserve/extensions/`
with its test, filed by origin like everything else here. See DISPOSITIONS.md.
