# Assimilative causal metrics and influence ranges

Every quantity on this page is conditional on the supplied model, prior
and observed record. The model supplies the likelihood and is not
estimated from the record, so a positive value is influence under the
dynamics that were supplied: it is not an empirically identified causal
effect, an intervention effect or a significance statement. Gaussian
relative entropy is oriented as smoother relative to filter, and its
per-time value is a pointwise information gain in nats; summed over the
record with the step it is a time-integrated value in nats times model
time. The value carries the prior and the integration step as well as
the model; [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)
gives the structurally independent null's return at two steps and at a
step coarse against the prior, each with the prior variance, horizon,
scheme and regularization status it was measured under. A normal
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) call uses
the supplied-code backward-ODE headline smoother, including its
correlated-noise correction, independently of `keep`.
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
and `aci(table = ...)` instead use the complete online Theorem 3
smoother; their finite-grid diagonal can therefore differ from headline
ACI.
[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
summarizes the duration of influence on the discrete time grid. A finite
adaptive table is labelled `objective_on_truncated_table`; its
`tail_bound` field is a heuristic tail estimate and must not be
interpreted as a certified error bound. It is a diagnostic under the
retained record, not a guarantee about the cells the truncation dropped.
The `l1_linf` estimator is a ratio, integrated with composite Simpson by
default, following the ACI reference code; `quadrature = "sum"` uses the
L1 grid-function sum instead. The `exact` objective is a finite
threshold sum with no time-axis quadrature, so it is unaffected by that
choice.

## References

Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
identifying forward and backward causal influence ranges using
assimilative causal inference. arXiv:2510.21889v2, 4 August 2026.
[doi:10.48550/arXiv.2510.21889](https://doi.org/10.48550/arXiv.2510.21889)

## See also

[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md),
[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
