# Pointer: `model_l84(variant = "fbcir")`

The FBCIR arm of `model_l84` is not a separable block. `model_l84`
(`aci/R/benchmark_models.R:321-388`, tree `97f6b124`) is one function whose
two presets differ only in the parameter row at `:346-350` (`F_amp` 0 vs 3,
`sig` 0.1 vs 0.2) plus the FBCIR conditional-metadata block at `:378-383`.
Its **default** variant is `"p3"` (jiang2026enkbs), so the whole constructor
is filed under EnKBS:

    reserve/enkbs/code/model_l84.R

There is no ACI_code-scoped arm; the constructor is absent from the acir
0.1.0 mainline entirely. Whichever of FBCIR or EnKBS lands first restores the
file; the second release then only needs its own preset row verified.
