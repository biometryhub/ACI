# reserve/fbcir-paper - FBCIR-paper-backed, no MATLAB

**Currently empty.** No block excised from `aci` 0.0.30 for acir 0.1.0 is
FBCIR-paper-backed *without* corresponding `FBCIR_code-main/` MATLAB. The one
near-miss is `model_tipping_triad(partition = "joint")`, which `aci` marks
paper-only at `benchmark_models.R:270-271` while the other two partitions map
to `climate_tipping_*.m`; because the arm is inseparable from the
MATLAB-backed constructor, the whole file sits in `reserve/fbcir/code/`.

The directory exists so the category is visible in the audit rather than
silently absent.
