# Pointer: `model_dyad(observe = "y")`

The reverse `x -> y` partition is the `else` branch at
`aci/R/benchmark_models.R:127-135` (tree `97f6b124`), which builds a bare
`stochastic_model` rather than a `cgns_model`. `aci` labels it itself at
`:148-150`: "package extension ... the reverse x -> y partition is not a
supplied paper/MATLAB inference benchmark". It is therefore an **extensions**
item, not an ACI_code one.

It is not separable from the constructor: the complete original `model_dyad`
body is filed once at

    reserve/enkbs/code/model_dyad_p3_p4_and_observe_y.R

Restoring this branch alone means taking `:127-135` plus the `observe == "y"`
arms of the metadata switches at `:159-163` and `:168-170` from that file.
Note `model_dyad(variant = "p3", observe = "y")` is EnKBS-scope, not an
extension: it is the published EnKBS dyad experiment, and `test-14` needs it.
