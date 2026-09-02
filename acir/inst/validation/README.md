# Numerical validation of acir

acir's numerical validation does not use a scenario registry. Every graded
quantity is checked against outputs of the method authors' own MATLAB
programs, hoisted into callable form as byte-exact extracts and pinned by
hash, so the oracle is a source this package's authors did not write.

The evidence is in three places in this package:

- `inst/evidence/register.csv`: one row per checked feature of the public
  surface, naming the check, what it is checked against, the tolerance class
  and the hash-pinned fixture behind it. A test fails the build if an
  exported function has no row or a row names a fixture whose bytes moved.
- `tests/testthat/fixtures/oracles/`: the reference outputs and their
  manifests, compared file by file in the test suite.
- `inst/evidence/gate_liveness.md`: every gate shown to fail once on a
  deliberate violation, with the test that performs it.

Edge-of-domain behaviour is tested directly: short records, a scalar against
a multi-dimensional hidden state, the covariance policy at the edge of the
positive-definite cone, cancellation on long records, and the generic against
the library model route.

## Waiver

The scenario-registry and execution-grid formats are waived for this release.

- Reason: the byte-pinned reference oracle above is the stronger evidence for
  a reimplementation, and every capability of the public surface is executed
  by the test suite from the built package under `R CMD check`.
- Version: 0.1.0.
- Owner: Max Moldovan.
