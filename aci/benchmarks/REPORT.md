# aci 0.0.30 compiled-CGNS benchmark report

## Result

The production architecture is proven on the frozen 3,001-point scalar dyad:

> function-valued model -> compiled coefficient bundle -> scalar kernels

All 23 benchmark parity gates passed exactly. The full package suite passed
2,552 expectations with no failures or warnings (one opt-in external-fixture
skip), and the built source package completed `R CMD check --no-manual` with
`Status: OK`.

## Same-machine dyad comparison

Medians are followed by the interquartile interval in milliseconds. Every row
uses the same frozen observation vector (SHA-256
`64d44c33341708e49c6946b779757e02d2c37033e567ee40a1bcd0784c5cbed1`)
and the same R, machine, BLAS and package process setup recorded with the raw
results.

| Operation | aciR 0.2.3 | aci 0.0.21 | aci 0.0.30 | aci improvement |
|---|---:|---:|---:|---:|
| Filter call | 1.517 [1.478-1.557] | 255.209 [254.219-257.608] | 2.144 [2.057-2.188] | 119.0x faster |
| Smoother call with supplied filter | 1.888 [1.861-1.918] | 375.906 [373.994-377.766] | 30.404 [30.109-30.929] | 12.4x faster |
| Complete ACI call | 3.533 [3.468-3.598] | 752.409 [749.317-757.243] | 3.041 [2.963-3.190] | 247.4x faster |

The improvement column is the like-for-like change between aci releases. The
aciR figures are orientation, not speedup denominators: its filter does not
compute aci's predictive log-likelihood, and its complete result does not have
aci's identical validation, decomposition, provenance and construction work.
The supplied-filter smoother exposes the remaining public-boundary cost in
aci: the trusted warm smoother kernel itself is 0.686 ms, while the public call
recompiles coefficients and validates every supplied covariance.

## aci 0.0.30 stage evidence

| Scenario and contract | Median [Q1-Q3] ms |
|---|---:|
| Dyad production compilation | 0.125 [0.118-0.137] |
| Dyad generic closure compilation | 136.417 [135.169-138.187] |
| Dyad warm filter, including likelihood/path | 1.964 [1.926-2.046] |
| Dyad warm smoother | 0.686 [0.665-0.721] |
| Dyad warm decomposed metric | 0.218 [0.208-0.229] |
| Dyad complete execution over bundle | 3.076 [2.992-3.311] |
| Affine matrix production compilation (`k=2,l=2`, 401 points) | 19.447 [19.272-19.673] |
| Affine matrix complete execution over bundle | 53.590 [53.247-53.897] |
| Affine matrix complete public ACI | 73.273 [72.675-73.334] |
| Generic conditioned compilation (`k=2,l=1`, 401 points) | 20.011 [19.559-20.318] |
| Generic conditioned complete public ACI | 58.678 [58.466-59.541] |
| Bounded forward lag table (201 points, maximum 25 lags) | 141.089 [140.180-142.408] |
| Streaming forward CIR (201 points, no retained triangle) | 463.692 [458.377-466.348] |

The current public dyad ACI allocation record contains 89 sized events,
2,009,784 sized bytes and no unsized `new page` events. The 0.0.21 complete
call recorded 3,373 events, including 3,345 unsized `new page` events and
877,664 sized bytes. Sized bytes must not be added to unsized pages, so these
are reported separately rather than collapsed into a misleading allocation
total.

## Cause of the improvement

aci 0.0.21 repeatedly called coefficient closures and executed scalar problems
through 1-by-1 matrix algebra. aci 0.0.30 instead:

- realises deterministic coefficient functions once per observation grid;
- authenticates directed dyad and affine-model batch realisers while retaining
  a generic one-pass closure fallback;
- dispatches explicit scalar filter, likelihood, smoother and KL arithmetic;
- dispatches the corresponding matrix, implicit, conditioning, lag and online
  routes over the same compiled contract;
- streams forward-CIR row reductions and reuses a maximal compiled prefix for
  backward CIR; and
- no longer installs the retired per-step filter, smoother or lag kernels.

No public numerical default, equation, indexing convention, covariance floor
or CIR quadrature policy was changed to obtain these timings. Valid-path
warning classes and result metadata are preserved. Coefficient closures are now
explicitly required to be deterministic over the realised grid, so an invalid
terminal coefficient can be rejected during compilation earlier than it was by
the interval-by-interval 0.0.21 filter.

## Evidence and limitation

The test suite includes authors-source scalar dyad, predator-prey and CIR
fixtures; independent correlated scalar/matrix and CIR transcriptions; a
source-derived ENSO path; direct likelihood/conditioning equations; frozen
0.0.21 implicit-stepper regression values; and structural/public-route tests.
The supplied FBCIR archive contains source but no author-produced numerical
outputs, and the installed MATLAB cannot obtain a licence. FBCIR author-output
parity is therefore not claimed; the package's equations and coefficient
routes are tested, and a hash-pinned authors-source fixture generator is
included for a future licensed run.

Raw current evidence is under
`benchmarks/results/production-compiled-cgns-20260826/`. The untouched 0.0.21
and aciR measurements are under
`benchmarks/results/scalar-dyad-20260825-final/`. See `README.md` and
`DECISION.md` for the reproducibility and comparison protocol.
