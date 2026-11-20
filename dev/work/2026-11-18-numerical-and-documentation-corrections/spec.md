# Behaviour and verification requirements

- Use the existing scalar Gaussian signal `(mu_q - mu_p)^2 / (2 * R_q)`.
  Change evaluation only when the existing signal is infinite. Retain finite
  signal values bit for bit, and leave the dispersion results unchanged by
  this additional correction.
- Recover finite signals in the square-overflow case, overflow of the full
  mean difference and the tested subnormal-variance cases. Preserve
  infinity when the result exceeds the rounding range of double precision.
  Test the overflow boundary explicitly. This is not a claim of universal
  correctly rounded Gaussian KL evaluation.
- Preserve return names, dimensions, decomposition and errors. The scope is
  the scalar path calculation. `aci_metric_pair()` and multivariate arithmetic
  retain their existing numerical limits.
- Exclude the optional near-equality series and its tests from this batch.
- Run the complete suite with saved MATLAB comparisons enabled. Require the
  existing moment tolerance below 1e-12 and component tolerance below 1e-11.
  Require before/after differences at most 1e-12 and identical moment bytes.
  MATLAB itself is not run.
- Check exact represented-input high-precision references separately from
  source agreement. Preserve and report adverse boundary results.
- Compare helper cost and all nine benchmark stages against the base commit.
  Repository stage timing is advisory. Investigate a
  new helper slowdown over 25% on previously valid workloads. Timings of
  formerly incorrect overflow cases are reported separately.
- Run available fixture provenance, package build and local as-cran checks.
  Record unavailable lint, coverage and other CI configurations explicitly.

No new export, fixture, tolerance or dependency is introduced.
