# Numerical validation of acir

acir's numerical validation does not use a scenario registry. Each graded
quantity is checked against a hash-pinned reference record whose class its
register row names. Where the method authors' own MATLAB programs compute the
quantity, that program's output is the comparator, hoisted into callable form
as a byte-exact extract. Where no upstream counterpart exists (the noise
cross-covariance terms, the matrix online smoother, the general auxiliary
matrices), the comparator is an independent MATLAB or R transcription of the
published equations, or an analytic identity solved in the test file. The
register separates the two, and a reader must not read one as the other.

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
