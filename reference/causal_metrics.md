# Assimilative causal metrics and influence ranges

Gaussian relative entropy is oriented as smoother relative to filter. A
normal [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)
call uses the supplied-code backward-ODE headline smoother, including
its correlated-noise correction, independently of `keep`.
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
and `aci(table = ...)` instead use the complete online Theorem 3
smoother; their finite-grid diagonal can therefore differ from headline
ACI.
[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
summarizes the duration of influence on the discrete time grid. A finite
adaptive table is labelled `objective_on_truncated_table`; its
`tail_bound` field is a heuristic tail estimate and must not be
interpreted as a certified error bound. The `l1_linf` estimator is a
ratio: the forward ratio and the exact form are integrated with
composite Simpson, following the ACI reference code.

## References

Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
identifying forward and backward causal influence ranges using
assimilative causal inference. arXiv:2510.21889v2, 4 August 2026.
[doi:10.48550/arXiv.2510.21889](https://doi.org/10.48550/arXiv.2510.21889)

## See also

[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md),
[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
