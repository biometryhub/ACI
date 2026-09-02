# Data assimilation and finite-lag API

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
and
[`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
reconstruct hidden states: from the record up to each time, from the
whole record, and from the record up to a fixed number of steps ahead of
each time.
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
stores the finite-lag divergences used by the CIR estimators. The `lt_*`
helpers access a table without depending on its storage representation.
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
describes the conditional question, and
[`aci_conditional_reduce()`](https://biometryhub.github.io/ACI/reference/aci_conditional_reduce.md)
carries out the model reduction its `method = "reduce"` asks for. The
historical
[`lt_tail_bound()`](https://biometryhub.github.io/ACI/reference/lt_tail_bound.md)
name is retained for compatibility, but its value is a heuristic tail
estimate, not a certified mathematical error bound. A lag table uses the
complete online smoother of andreou2026cir Theorem 3 (Appendix G.1) as
its reference. That reference costs O(N) time-point work; table
construction then costs work proportional to the retained lag cells,
with O(N^2) cells for a full table in the worst case.
[`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
costs O(N) whatever the lag.

## Scheme

Two discretizations of the same continuous-time smoothing problem are in
use, and `meta$scheme` on a path says which one produced it.
`"backward_ode_euler"` is
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md):
the continuous backward smoothing equations integrated with an Euler
step. `"theorem3_discrete"` is
[`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
and the lag table's reference smoother: the exact conditional law of the
hidden state given the observed increments on the sampling grid, under
the explicit single-step discretization. They agree only to first order
in the step, so
[`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
at `lag = Inf` does not reproduce
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
and the gap grows with the length of the record rather than settling to
a constant. See the Scheme section of
[`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
for the measured size on the packaged ENSO partition.
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) reports
the scheme its own result was built under in `meta$smoother_scheme`.

## References

Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
identifying forward and backward causal influence ranges using
assimilative causal inference. arXiv:2510.21889v2, 4 August 2026.
[doi:10.48550/arXiv.2510.21889](https://doi.org/10.48550/arXiv.2510.21889)

Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother
with closed-form solutions and information-theoretic lag selection for
conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
**36**(4), 71.
[doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x)

## See also

[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
