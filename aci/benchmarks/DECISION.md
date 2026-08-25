# Production benchmark protocol

## Scope

The package now routes closed-form CGNS execution through:

> function-valued model -> production compiled coefficient bundle -> scalar or
> matrix filter, smoother, metric, lag, and CIR kernels

The benchmark therefore measures that architecture directly. It does not retain
or call the pre-0.0.30 implementation. Historical result directories remain
untouched as evidence of the earlier package state.

## Evidence gates

A reportable run must satisfy all of the following:

1. the frozen 3,001-step dyad input hash is exact;
2. production dyad compilation and the generic one-pass fallback realize equal
   coefficients;
3. precompiled complete scalar, affine-matrix, and conditioned runs agree with
   the current public API for filter moments, predictive log-likelihood,
   smoother moments, and ACI;
4. bounded private/public lag and streaming-CIR results agree;
5. raw observations, allocations, profiles, stage contracts, source hashes, and
   environment details are retained; and
6. compilation, warm execution, complete-bundle execution, and public workflows
   are reported separately.

These parity checks show that timed private stages implement the same current
production work as the public calls. MATLAB/source-fixture correctness belongs in
the package test suite; it is not replaced by a benchmark self-comparison.

## Comparison policy

There is no single package-wide ratio. A matched numerical-kernel comparison may
compare warm kernels with warm kernels. A user-facing comparison must compare
complete public workflows and state their different contracts. Allocation totals
must report both sized bytes and unsized `new page` events.

aci 0.0.21 is either represented by the preserved raw 2026-08-25 run or executed
in a separate vanilla R process from `aci-prev`. It must not share a namespace
with the candidate package, and current retired internals must not be recreated
for benchmarking.

Optional aciR or MATLAB timings are separate comparators with separate contracts
and environments. They are not silently merged with current aci timings. In
particular, predictive likelihood, validation, provenance, result construction,
conditioning, and numerical-policy differences must be stated alongside any
comparison.
