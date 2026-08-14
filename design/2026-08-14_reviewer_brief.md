# aciR — external reviewer brief

**Version** 0.1.0.9000 · **Date** 2026-08-14 · **Status** pre-CRAN, private repo

Intended for an external reviewer working from the source tree. Written to be
pasted whole. Read §6 before circulating: one class of finding in here concerns
a third party's published code and its disclosure is not yet decided.

---

## 1. What the package is

An R implementation of Assimilative Causal Inference for Conditional Gaussian
Nonlinear Systems — Andreou, Chen and Bollt, *Nature Communications* 17:1854
(`10.1038/s41467-026-68568-0`), with the fixed-lag online smoother of Andreou,
Chen and Li (`10.48550/arXiv.2411.05870`).

It computes, for a partially observed system, how much the future of the
observed signal says about the unobserved state (the **ACI metric**) and how far
into that future one must look before the answer stops changing (the **causal
influence range**). Scalar and vector states, both causal directions, and
conditional ACI.

**Base R plus nothing.** No compiled code, no `Imports` beyond base. 652 tests.
`R CMD check --as-cran`: **0 errors, 0 warnings, 1 NOTE** (the NOTE is three
404s from the declared URLs, which resolve once the repository is public).

## 2. What has already been verified — do not re-derive this

The upstream MATLAB is not a library; it is seven top-to-bottom scripts,
12,339 lines, of which only two files are callable and both are third-party.
A harness under `oracle/parity/` hoists the authors' computational passages into
callable functions as **byte-exact slices** between markers — 1,681 verbatim
lines across 25 extracts, nothing retyped, with a checker that fails on one byte
of drift.

| Gate | Result |
|---|---|
| Each extract vs the script's own workspace | **124 / 124, difference exactly 0** |
| Committed dyad fixture vs authors' code, published N = 30000 | **exactly 0**, 8 quantities |
| Predator-prey fixtures vs authors' code, published N = 12000 | **0** (dir 1), 1.5e-13 (dir 2) |
| aciR vs reference, scalar core, two datasets | 26 / 26, worst 1.53e-14 |
| aciR vs reference, ENSO vector core | 5 / 5, worst 5.73e-14 |

Reproduce: `oracle/parity/tools/`, driven from R; needs MATLAB and a clone of
`github.com/marandmath/ACI_code` at `matlab_reference/`.

## 3. Where to spend review effort, in order

**(a) The unexplained 4.58e-09.** With the reference's comparison horizon, the
subjective range agrees to 2.22e-16 and the exact objective range to 1.59e-14,
but `objective` — the efficient approximation, `simpson(re) * dt / peak` — agrees
only to 4.58e-09. The obvious suspect, a differing odd-length Simpson closure,
was tested and **eliminated**: the closures do differ (up to 4.3e-03 on a test
function) but the CIR integrand decays to ~0 at its endpoint, so the closure term
is multiplied by nothing. No replacement hypothesis. This is the single best
target in the package.

**(b) The saturation policy.** aciR marks a reported time unresolved when its
range consumes more than `1 - margin` of what it examined; the reference uses an
absolute lookahead guard. At the reference's horizon aciR marks **every** time
saturated where the reference reports numbers. Is `margin` the right instrument,
is 0.1 the right default, and is a fraction the right shape?

**(c) `.aci_cir_row()` and the horizon.** The comparison sequence is truncated
but the fully informed posterior it is compared against is not. That asymmetry is
deliberate and load-bearing. `tests/testthat/test-cir-horizon.R` asserts the
truncated row is a prefix of the untruncated one. Attack the argument, not just
the code.

**(d) The vector online smoother**, `R/aci-online-smoother-mv.R`. Matrices do not
commute, so the scalar version's cumulative-logarithm reconstruction has no
analogue; the ordered products are accumulated and truncated on a spectral-radius
argument. Check the truncation tolerance is doing what the comment claims.

**(e) Numerical routes that differ by design.** The multivariate relative entropy
uses a Cholesky factorisation where the reference uses `det` and a right
division; the observation-noise inverse uses Cholesky where the reference uses
`pinv`. Measured differences 6.15e-14 and 1.14e-13 respectively. aciR is stricter
— it requires positive-definiteness at every slice and will decline systems the
reference analyses. Is that the right trade?

## 4. Where NOT to spend effort

The filter, smoother and ACI metric on scalar and vector states are graded
against the authors' code at machine precision on multiple datasets, including
one built specifically to carry a non-zero noise cross-covariance that the
published fixtures cannot exercise. Re-deriving the CGNS algebra is unlikely to
pay.

## 5. Known gaps, stated so they are not "found"

- Four ENSO scripts are captured and validated but their extracts are **not
  gated**; each needs an input list derived individually.
- aciR has **not** been compared against the reference on predator-prey — the
  extracts and gate are complete, the aciR-side comparison is not run.
- `aci_conditional()` has an executable upstream counterpart (see §6) but is not
  yet graded against it.
- The ENSO simulation uses a Milstein correction whose derivation is unclear;
  documented at `design/2026-08-13_milstein_anomaly.md`, unresolved.

## 6. Handling note — read before circulating

Two findings concern **the reference implementation's own behaviour**, not this
package:

- `noisy_predator_prey_model.m` computes direction one's smoother from direction
  two's filter if run top to bottom, silently. Its own section banners instruct
  otherwise, so it is a hazard for automation rather than a defect in use.
- Of the five ENSO scripts, exactly one ships with conditional ACI **enabled** —
  a line distinguished from surrounding guidance comments by leading whitespace
  rather than a comment marker — so it answers a different question from its four
  siblings.

Both are reported here because a reviewer of this package needs them to judge the
harness. Neither has been raised with the authors, and whether to do so is the
maintainer's decision alone. Treat as internal until told otherwise.
