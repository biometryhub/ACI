# Pointer: `model_dyad(variant = "p4")`

The `p4` preset (`moser2026extremes` eqs. (4.1)-(4.2),
`aci/R/benchmark_models.R:94,145,157`, tree `97f6b124`) is a single `switch`
arm interleaved with the retained `p1` arm and four metadata switches, so it
is not re-appliable in isolation. The complete original `model_dyad` body,
carrying `p1`, `p3`, `p4` and both `observe` branches, is filed once at:

    reserve/enkbs/code/model_dyad_p3_p4_and_observe_y.R

filed there because `p3` is the arm with a MATLAB-backed family behind it.
Restoring `p4` means taking the `p4` rows of `:92`, `:145`, `:157`, `:161`
and `:169` from that file, not applying the file wholesale.
