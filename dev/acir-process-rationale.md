# How `acir` 0.1.0 was built, and why

**Audience.** The supervisor and collaborators, and anyone who later asks why a
particular line of `aci` or `aciR` did or did not survive into the joint
package. It reads end to end without the working documents, but every
consequential claim carries a pointer to the record that holds the measurement.

**Status of the numbers.** Every figure in this document was produced by code
run in this workspace. An earlier draft carried fifteen provisional markers on
figures that later work could move; all fifteen are resolved here, and no
figure is now provisional. Two kinds of figure are closed and will not move:
those from the comparison of the two frozen packages (`aci` 0.0.30 against
`aciR` 0.2.3), and those from a change entry, each of which was measured
against a specific commit that is in the history.

**The state this describes.** The package tree at commit
`3db97dd` (the surface rename, on top of the 41-commit merge history). Full
test suite 6236 pass, 0 fail, 0 error, 1 skip. `R CMD check --no-manual` on
the built tarball, with vignettes, tests and examples all run: `Status: OK`,
zero warnings, zero notes. 24 exports, 13 R source files (7451 lines), 22
test files (6733 lines), 45 help topics, 30 oracle fixture CSVs under two
manifests, and three vignettes. After this document was finalised the public
surface was renamed to the `aciR`-style `aci_*` interface, so exported names
below read as they now ship, while `aci` 0.0.30 and `aciR` 0.2.3 names stay
historical.

**How to read this.** Changes carry a short id throughout; the key below gives
each a plain name, and section 4 treats them in that order, 4.1 to 4.18,
with 4.19 collecting candidates that were measured and not adopted. Each
is one commit with its own evidence gate, and the ids fall in three rounds of
adoption, C1, C2 and C3. Section 8 names and resolves decisions D1 to D8.
Measurements are cited by the study or entry that made them, and those records
are in the evidence archive at the end.

| Id | Change |
|---|---|
| C1a | Cancellation-resistant scalar KL arithmetic |
| C1b | Optional predictive log-likelihood |
| C1c | RNG containment on seeded simulation |
| C1d | Matrix recursion arithmetic |
| C1e | Vectorised affine coefficient realisation |
| C2a | Trusted supplied-filter smoother route |
| C2b | Forward CIR window, exact objective, named read-outs |
| C2c | Target by name, first-slice convention, two tau estimands |
| C2d | Authenticated whole-path ENSO realisation |
| C3a | Explicit inverses, and the signed-zero finding |
| C3b | Scalar log-likelihood kernel; memory-table withdrawal |
| C3c | Exact positivity test, scalar hidden state |
| C3d | Strict covariance policy, with opt-in flooring |
| C3e | Simulated trajectories carry channel names |
| C3f | `aci_online()`, the fixed-lag online smoother |
| C3g | Source-derived fixtures, scalar ENSO partitions |
| C3h | The `T_C`-hidden zeroth-order partition |
| C3i | Two vignettes, and the closing check |

---

## 1. Why one package, and what "closed form" decides

### 1.1 The starting position

Two R packages implemented overlapping parts of the same science. `aci`
(0.0.30) is a large package with a compiled conditional-Gaussian kernel
architecture, an ensemble engine, a model-discovery loop, an extremes family,
a formula front end, and applied workflows. `aciR` (0.2.3) is a smaller
package, closed form throughout, with a different components schema and a
different set of numerical policies. Both are graded against the same MATLAB
reference codebase, `ACI_code-main/`.

The supervisor's direction settled two things at once: the joint artefact is a
new package rather than an integration into either existing tree, and the
engine package should carry only closed-form solutions, taken from `aci` where
possible.

### 1.2 What "closed form" means operationally

The framework rests on the conditional Gaussian nonlinear system (CGNS) class:
observed `x` and hidden `y` coupled so that every coefficient may depend
nonlinearly on `(t, x)`, the hidden state enters both drifts only linearly, and
no coefficient depends on `y`. Conditioned on the observed path, the hidden
process is then exactly linear Gaussian, and every quantity the method needs is
either an exact finite-dimensional recursion or an algebraic formula in the
resulting means and covariances: the filter, the backward smoother, the
fixed-lag online smoother, the ACI metric as a Gaussian Kullback-Leibler
divergence, and the forward causal influence range. Deterministic given the
data, exact up to time discretisation, one pass over the grid.

The criterion therefore draws a scientific-class boundary, not a quality
boundary. It admits exactly the estimators that evaluate CGNS conditional
statistics analytically, or algebraic functionals of them. It excludes the
ensemble family (EnKBF and EnKBS approximate the same conditional distributions
by Monte Carlo for systems outside the CGNS class) and, by extension, the
applied layers built on either engine.

Boundary cases were classified individually rather than by family label:

| Case | Classification | Where it landed in `acir` |
|---|---|---|
| ENSO with `T_C` hidden | The exact system leaves the CGNS class. The MATLAB authors restore membership by the zeroth-order substitution `c1(t, T_C) -> c1(t, 0)` in the assimilation model only, then solve exactly. Closed-form solution of an approximated model. Inside the criterion. | Implemented in 0.1.0 (C3h) |
| Conditional ACI (precision masking, prescribed forcing) | Closed form. Inside. | Mainline; API completed in C2c |
| Backward CIR and adaptive lag (the FBCIR family) | Closed-form mathematics. Excluded from 0.1.0 by the release scope rule, not by the criterion. Returns in a later release. | `reserve/fbcir/` with a patch |
| Simulation (Euler-Maruyama, Milstein) | Data generation, not a solution. Retained because the package needs it. | Mainline |
| FFBS sampling | Uses the closed-form smoother but is a Monte Carlo routine serving discovery. Travels with discovery. | `reserve/enkbs/` |
| `aci`'s implicit stepper | A discretisation policy for the same exact equations. Closed form, no MATLAB counterpart. Retained. | Mainline |

### 1.3 Was the extraction actually possible

This was verified rather than assumed, before any code moved; the separability
report holds the full map.

In the whole of `aci` 0.0.30 there are exactly four code statements in the
engine set that call a symbol defined in an excluded file. All four sit inside
one function, in the `engine == "ensemble"` branch of `aci()`
(`aci/R/causal_metrics.R:393, 396, 397, 412`). The five compiled kernel files
have zero references to excluded files in either direction, and nothing in the
excluded families reaches compiled internals. There is no shared mutable state,
no cross-family cache, no registry, no load hooks, and no `Collate` ordering to
unwind. One block of excluded-family code was found living inside a retained
file (Gaspari-Cohn localization and inflation, `aci/R/assimilation.R:903-1025`,
consumed only by the ensemble engine) and one group of retained methods was
found living inside an excluded file (three `plot` methods for engine classes,
`aci/R/formula_interface.R:617-674`); both were relocated rather than lost.

The report enumerates fifteen cut points, each with its one-line adaptation.
The extraction itself was then verified rather than trusted: 116 quantities
bit-identical to `aci` 0.0.30 on the retained surface, including all 2 003 001
stored lag-table cells, with `NAMESPACE` and `man/` independently regenerated
and matching line for line.

### 1.4 What the criterion does not decide

The criterion does not choose between the two packages. `aciR` is entirely
closed form by construction: nothing in its sixteen R files belongs to the
excluded families, and it cannot build a non-CGNS model at all. So
"take them from `aci`" settles the engine chassis, not the merge. `aci`
supplies the compiled kernel architecture and the model interface; `aciR`
supplies the better half of several numerical policies and result semantics.
Which half is which is decided per component, on measurement, in section 4.

One consequence is worth stating explicitly: CRAN treats package names
case-insensitively for conflict purposes, so `acir` and `aciR` cannot both
exist on CRAN. Shipping `acir` forecloses a separate narrow `aciR` submission.
That is a feature of the one-package decision rather than a cost, but better
said now than discovered later. The package name `acir` is settled and final
(settled by the package author, 28 August 2026).

---

## 2. How the evidence was gathered

### 2.1 The three rules that shaped the method

**Correctness is a gate, efficiency is a rank.** An implementation enters the
comparison only if it is faithful to the reference equations and passes its
evidence gates. Efficiency then ranks implementations that are already correct.
A faster kernel that loses cancellation resistance or a verified numerical
policy is not more efficient; it is wrong faster. No pre-fixed pass or fail
performance threshold exists anywhere in this merge.

**Contracts are mapped before numbers are compared.** Two functions with the
same name in the two packages routinely compute different things. Every
comparison states the contract on each side first, states the normalisation
applied, and separates agreement claims from convention claims.

**No package-wide ratios.** Performance is reported per component and per stage
contract: component realisation, warm numerical kernel, and complete public
workflow, separately. A single package-wide ratio would have been misleading in
both directions on this pair, and section 3.2 shows how.

### 2.2 The instrument

Four numerical comparisons and one timing study were run on the two frozen
packages, all on shared deterministic inputs, both loaded from a scratch
library, and with nothing outside the scratch tree created or modified.

| Study | Scope |
|---|---|
| Joint ENSO comparison | Joint ENSO: observed `(T_C, T_E, I)`, hidden `(u, h_W, tau)`. Filter, smoother, ACI, fixed-lag online moments, CIR divergence rows, CIR summaries. Two records, 4001 steps each. |
| Dyad online/CIR comparison | Dyad fixed-lag online smoother and forward CIR, against the pinned fixture and a 3001-point path. Includes an independent R transcription of `ACI_code-main/simps.m`. |
| Scalar-partition comparison | The three scalar ENSO partitions (`u`, `h_W`, `tau` hidden), kernel comparison on identical realised coefficient arrays. Includes the tau estimand test. |
| Conditional comparison | Conditional (masked) ACI: `aci`'s inflate strategy against `aciR`'s `aci_conditional`, plus the MATLAB first-step convention quantified. |
| Timing study | Per-component stage-contract timings and memory, forty stage definitions, contracts recorded per row. |

Four features of the instrument are what make the results usable.

**A sanity gate before any comparison.** The joint ENSO comparison first
reproduced the pinned `enso_reference.csv` from `enso_signal.csv` through both
packages, on all eleven recorded columns at 201 sampled indices: worst case
8.39e-14 (`aci`) and 8.51e-14 (`aciR`). That gate shows neither package has
drifted from the pinned path. It is explicitly not a second oracle, because the
fixture set is a byte-identical copy shared by both packages.

**Pre-registered normalisations.** Each report lists the transformations
applied before comparison (mean layout, hidden-component order, subjective
read-out offset, time-grid handling, KL orientation) and, where a normalisation
could flatter one side, carries the over-corrected control in the output CSV.
The joint ENSO subjective read-out is the clearest case: applying the
one-grid-step offset unconditionally leaves a spurious 5.00e-03 residual on
8152 of 51471 cells, so the clamped form is the correct alignment and the
unclamped form is recorded beside it as the control.

**Built-in controls.** The `h_W` case against the tau case (section 6.1), and
the timing study's dyad rows against the previously recorded benchmark rows,
are controls rather than results. That cross-check found eleven dyad rows all
within +15.4% and -6.3% of the earlier record, nine of them inside 5%, with the
largest departure on an `aciR` row while the `aci` rows it is paired against
moved at most 1.6%. The instrument does not flatter `aci`.

**Timing method.** `microbenchmark::get_nanotime()` at roughly 41 ns measured
granularity, three untimed warm-ups per stage, one garbage collection per
group, one timed call per repetition, and stages interleaved inside the
repetition loop so that machine drift is shared across packages rather than
accumulating against one of them. Warm inputs are built outside the timed
region. Routes are confirmed at run time rather than assumed.

### 2.3 The gates each adopted change had to pass

After the comparison, each change into the new package is one reviewable commit
with a fixed gate structure, recorded in the per-change engineering log:

1. **A regression capture.** A fixed set of quantities (dyad filter, smoother,
   ACI, log-likelihood; the ENSO joint case; a scalar ENSO partition; a
   conditional case; forward CIR at four anchors; the lag-table diagonal;
   seeded, unseeded and rejected simulation) computed before and after, and
   compared value by value. The standard is bit identity unless a tolerance is
   explicitly budgeted and justified.
2. **The full test suite**, with the new assertions counted separately.
3. **A same-process A/B timing**, where the pre-change code is reconstructed
   from the previous commit, byte-compiled, bound into the namespace, and
   proved to reproduce the pre-change values before any timing is taken. This
   removes between-build drift from the measurement.
4. **A tradeoff paragraph**, recording in plain terms what the change costs and
   what was deliberately not done.
5. **A verdict**, including where a stated target was missed.

Two consequences of this structure show that the gates are real rather than
decorative. The vectorised realiser was adopted with its 3x target openly
missed at 1.53x, with the decomposition showing why the premise of the target
was wrong. And the batch realiser's 30 ms target was recorded as met only after
scaling for a machine measured to be running 1.7 to 2.2x slow, with the
sentence "I did not measure it under 30 ms on a comparable machine, and I am
not claiming that I did" left in the log rather than smoothed away.

**The gate itself was upgraded mid-project, the single most important process
event in the build.** Through the C1 and C2 entries the regression comparator
used `identical()` and `max(abs(a - b))` on each recorded leaf. In the C3
entries, the explicit-inverse change (C3a) found a change that both instruments
reported as zero and that was nevertheless a real change to the bytes the
package produces: LAPACK's `dpotri` writes `-0.0` into the off-diagonal zeros
of a diagonal matrix's inverse where the triangular-solve route writes `+0.0`,
and `0 == -0` in R at the default settings of both instruments. The comparator
was rewritten to compare the serialized bytes of every leaf, so that a
sign-of-zero flip counts as a difference, and every C3 entry from C3a onward is
graded that way: 263 leaves and 1 475 087 values in C3a, rising to 330 leaves
and 1 777 953 values at the close of the round. Section 4.10 gives the finding
in full. The gate caught what the previous gate could not see, and the fix was
to change the gate rather than to argue that the difference did not matter.

### 2.4 How the work was sequenced

Everything after the extraction is one commit per change, in three sequential
rounds, made one at a time against the tree so that no two changes are ever in
flight in the same file.

| Stage | Commits | Entries | Suite at close |
|---|---|---|---|
| Faithful extraction | `293cc8e`..`43556ad` (17) | none; graded by the verification report | 1802 pass / 0 fail / 1 skip |
| First adoption round | `6d06eeb`..`8c381de` (6) | C1a to C1e | 2109 |
| Second adoption round | `c6f9933`..`0bec004` (7) | C2a to C2d | 3436 |
| Third adoption round | `227d67d`..`5e5ee5e` (10) | C3a to C3i | 6236 |
| Repository hygiene | `3106967` (1) | recorded in this document, section 7.1 | 6236 |

Preparation work ran alongside the second round in a scratch tree that never
touched the package: two prototype packages, the `T_C` port and the partition
fixtures, and six evidence studies covering the covariance policy, bundle
memory, inverse routes, the online smoother, the process and metadata
documents, and the vignette drafts. Nothing from that work entered the package
without its own entry and its own gate, and in two cases the C3 integration
contradicted the preparation and the integration won: the bundle-memory study's
recommendation to drop the scalar flat views was rejected on runtime cost, and
the process-document study's filing-register row count was wrong and is
corrected in section 7.2.

---

## 3. The verdict

### 3.1 Numerically, this is one algorithm in two packagings

Wherever the schemes and the contracts match, the two packages compute the same
numbers to machine precision. This is the single most important result of the
comparison, because it removes numerical quality as a selection criterion
almost everywhere.

| Quantity | Agreement, `aci` against `aciR` | Source |
|---|---|---|
| Joint ENSO filter and smoother moments, full 3-vector and full 3x3, 4001 steps, two records | <= 2.4e-15 | Joint ENSO |
| Joint ENSO ACI total | <= 6.3e-15 | Joint ENSO |
| Fixed-lag prefix-smoother moments, 36 `(j, n)` pairs | <= 1.7e-16 | Joint ENSO |
| CIR divergence-row cells, 62 `(j, lag)` pairs | <= 4.0e-15 | Joint ENSO |
| CIR row peak and objective range, 399 times | <= 1.3e-14 | Joint ENSO |
| Subjective range, 51471 threshold-time cells, after the documented one-grid-step normalisation | <= 2.2e-16 | Joint ENSO |
| Dyad online-smoother moments, all 75 fixture pairs | <= 4.2e-15 | Dyad online/CIR |
| The three scalar ENSO partitions, filter, smoother and ACI | <= 2.1e-15, filter covariance bit-identical | Scalar partitions |
| Conditional operator: the masked Gram itself | exactly 0 over all 4000 slices, both cases | Conditional |
| Conditional posteriors and ACI | <= 8.9e-16 | Conditional |

There is no numerical basis to prefer either kernel anywhere that was
measurable. Selection therefore falls entirely to efficiency and interface, per
component.

Two qualifications belong with that claim. First, the residual differences that
do exist trace to three documented arithmetic choices, none of which changes the
estimand: `aci` applies a symmetrise-and-floor step per iteration where `aciR`
symmetrises and guards the Cholesky; the smoother forms an explicit inverse on
one side and a solve on the other; and the metric dispatches to a scalar kernel
on one side and a matrix kernel on the other for hidden dimension one. Second,
agreement between the two packages is not by itself evidence of MATLAB
fidelity, because the fixture set that grades both is one shared oracle, copied
byte for byte. That limitation is stated in each report rather than papered
over.

### 3.2 Efficiency is mixed, and contract-sensitive

The timing result is split, which is why the stage-contract method mattered.
Everything in this subsection is from the comparison of the two frozen
packages; section 3.4 records where the merged package finished.

**The scalar dyad favours `aci` where it counts.** The complete public workflow
is 3.089 ms against `aciR`'s 3.554 ms while doing strictly more work:
compilation, predictive likelihood, signal and dispersion decomposition, and
provenance. Its warm scalar smoother at 0.684 ms is 3.18x faster than `aciR`'s
public smoother. `aciR` wins the filter pair (1.591 against 2.173 ms) but that
pair is mismatched, because `aci`'s filter also computes a likelihood.

**The matrix path reverses it, on matched contracts.** Warm ENSO smoother:
`aciR` 96.1 ms against `aci` 183.0 ms, a factor of 1.90, with no likelihood on
either side. Total-only matrix KL: 81.2 against 137.3 ms, a factor of 1.69.
Complete over already-realised inputs: 259.9 against 525.5 ms, a factor of
2.02. The gap is not explained by the likelihood, because the smoother pair
carries no likelihood term at all.

**The single largest finding is coefficient realisation.** `aci`'s ENSO
compilation measured 540.1 ms against 13.1 ms for the `aciR` analogue that
produces the same per-step arrays plus the observation-Gram inverse: a factor
of 41.3 on a matched contract, and roughly half of `aci`'s 1066 ms public ENSO
call.

**Forward CIR route cost.** `aciR`'s windowed row-at-a-time computation
answered the fifteen fixture anchors in about 0.0 s and the whole 3001-point
record in 0.5 s, where `aci`'s routes took 37.6 s and 94.9 s respectively.
`aci` computes every anchor, has no window argument, and needs
`options(aci.default_tol = 0)` to reproduce the fixture at all.

**Memory appeared to track the same story, and did not.** The timing study
reported the ENSO compiled bundle at 3.04 MB in `aci` (0.50 MB of it
provenance) against 1.56 MB for `aciR`'s components plus the Gram inverse, and
the dyad at 0.90 MB against 0.07 MB. That measurement was later falsified: it
used `utils::object.size()`, which counts a shared object once for every place
it appears. The correction is in section 4.11, and the work item it generated
was withdrawn rather than implemented.

A package-wide ratio computed from this table would have said "`aciR` is
roughly twice as fast" and would have been wrong about the scalar complete
workflow, wrong about the warm scalar smoother, and silent about the two
findings that actually mattered, which are the realisation row and the CIR
route cost.

### 3.3 Therefore: per-component selection

No package won. Each component was assigned on its own evidence, and every
assignment was then implemented under its own gate.

| Component | Selected | Deciding evidence | Entry |
|---|---|---|---|
| Model interface, constructors, grid and contract validation | `aci` | Time-grid safety (joint ENSO comparison, finding 6), realiser identity, and `model_enso6` is the only route to the scalar-partition coefficients | Extraction |
| Coefficient realisation, scalar (authenticated dyad) | `aci` | 0.117 ms; complete-workflow win | Extraction |
| Coefficient realisation, matrix and affine | Rebuild along `aciR`'s vectorised approach inside `aci`'s realiser architecture | 41.3x on a matched contract; identical arrays out | C1e, then C2d |
| Scalar filter and smoother kernels | `aci` compiled | Dyad complete 3.089 against 3.554 ms doing more; warm smoother 0.684 ms | Extraction, C3b, C3c |
| Scalar KL arithmetic | `aciR`'s `log1p` form, into the compiled kernel | Accuracy, not speed; see 4.1 | C1a |
| Matrix filter, smoother and KL kernels | `aciR`'s recursion organisation, with the hidden-dimension-1 defect fixed and the likelihood made optional | 1.90x, 1.69x and 2.02x on matched contracts; machine-precision agreement makes the swap safe | C1b, C1d, C3a |
| Fixed-lag online smoother | `aciR`'s public lag surface and scheme labels | Efficiency, plus `aci` has no public accessor at all | C3f |
| Forward CIR | `aciR`'s windowed row-at-a-time engine, plus `aci`'s persisted lag table as an optional diagnostic, with a merged API | 0.5 s against 94.9 s; the naming hazards of section 6.3 | C2b |
| Conditioning | Shared construction, since the two are bit-identical; merged API surface | Masked Gram identical to zero | C2c |
| Result objects | `aci` provenance plus `aciR` statuses | C2b statuses; the timing study's memory row was withdrawn, see 4.11 | C2b, C3b |
| Public plumbing | Trusted supplied-filter smoother route, compile-once reuse of the Gram inverse, optional likelihood, CIR window | Timing study, dyad online/CIR comparison | C1b, C2a, C2b, C3c |
| ENSO tau estimand | Reduced form as the fidelity default, full form as the documented alternative | Scalar-partition comparison, section 5 of that report; see 6.1 | C2c |
| `T_C` zeroth-order partition | New port from MATLAB | Coverage; neither package has it | C3h |

### 3.4 Where the merged package finished

Every efficiency row the comparison recorded against `aci` was worked, and the
closing positions are below. All are from a change entry's own measurement,
with that entry's own caveats about machine load carried across.

| Comparison finding | Closing position in `acir` | Entry |
|---|---|---|
| Matrix contracts 1.69x to 2.02x against `aci` | Reversed in the same process: warm smoother 80.79 against `aciR`'s 101.06 ms, matrix KL 69.13 against 84.23, warm filter without likelihood 46.43 against 70.23 | C1d |
| ENSO realisation 41.3x against `aci` | The realiser itself 212.6x faster (775.47 to 3.65 ms in one process); the ENSO compile is now 1.9x the `aciR` analogue, against 27x at the start of that entry | C1e, C2d |
| Forward CIR 37.6 s for fifteen anchors | 719 ms streamed and 713 ms through the public route, 57.9x, reproducing the pinned MATLAB peak to 1.03e-14 | C2b |
| Public supplied-filter smoother, a factor of 44 over its own kernel | 27.975 to 0.833 ms on the dyad, 33.6x, of which 0.680 ms is the kernel and 0.016 ms the authentication | C2a |
| Mandatory in-loop predictive likelihood | Optional; `loglik = FALSE` cuts the dyad warm filter 54.6% and the ENSO warm filter 30.6%. The scalar likelihood kernel itself is a further 2.44x | C1b, C3b |
| No public accessor for fixed-lag online moments | `aci_online()` exported, O(N) in the record and flat in the lag; the oracle workload drops 2.391 s to 0.209 s | C3f |
| Bundle memory, 0.90 MB against 0.07 MB on the dyad | Withdrawn. The measurement was an `object.size()` artefact; the provenance list costs 440 bytes | C3b |

Two comparison rows were not closed and are recorded as such: `lag_table()`
remains quadratic and has no window, so the fast CIR path and the storable one
are different objects (C2b tradeoffs); and the compiled bundle still carries a
flattened scalar view that costs 216 968 B on the dyad, kept because removing
it costs the scalar kernels 7 to 10% of their runtime (C3b).

---

## 4. What was adopted from where, and why

This section is sourced from the per-change engineering log, which holds the
full evidence for each entry. Everything below is one commit with its own gate.
Eighteen entries, in the order they landed.

### 4.1 Cancellation-resistant scalar KL arithmetic (C1a, from `aciR`)

The Gaussian KL dispersion term was computed as
`0.5 * (R_p/R_q - 1 + log(R_q) - log(R_p))` and is now computed as
`0.5 * (delta - log1p(delta))` with `delta = R_p/R_q - 1`. `aciR` already wrote
it this way.

This is an accuracy change, not a speed change, and it matters exactly where an
ACI run spends its time: near equal covariances. Measured against a
four-term Maclaurin series that is exact to double precision at these
magnitudes, the old form's relative error at a variance ratio within 1e-8 of
one was 5.00e-01, meaning it returned half the true value; the new form's was
2.60e-09, an improvement of 1.9e8. The improvement runs from about 2.0e4 at
`delta = 1e-4` to 5.9e5 at `delta = 1e-10`. The old form's best case is
`R_q = 1`, where one logarithm vanishes exactly; away from that, at
`R_q = 3.7`, the same deltas give direct-form relative errors of 1.00, 4.4e4
and 4.4e8.

Two things were learned that the plan had not anticipated. The reach is wider
than the kernel's name suggests: every hidden-dimension-1 model dispatches to
this scalar kernel, including all three scalar ENSO partitions and conditional
ACI on a single hidden variable, not only systems with one observed and one
hidden channel. And the regression gate's own premise was wrong in one place,
which the log records rather than hides: the conditional `h_W` case was
expected to be bit-identical because it is a matrix system, and it was not,
because with one hidden variable it runs on the scalar kernel. The 8.88e-16
move is the expected last-bit change.

Gate: 39 of 50 recorded quantities bit-identical, maximum absolute movement
8.88e-16 against a 1e-13 bound. Timing improved slightly (one fewer square
root, one fewer logarithm).

Left open and recorded: the matrix KL keeps the direct form and the same
weakness. That is a separate change needing its own evidence, and it was not
taken in 0.1.0.

### 4.2 Optional predictive log-likelihood (C1b, a contract change from the timing study)

`aci`'s filter computed a predictive log-likelihood inside the moment loop,
unconditionally. ACI never reads it. It is consumed only by the validation
diagnostics, which are not in this package. It is also real work: on the ENSO
system it is 60.5 ms of a 1052 ms call, and on the dyad it is most of the
filter.

A `loglik` argument was threaded through the whole filter chain, defaulting to
`TRUE` so that no existing call changes behaviour. The decision to default to
`TRUE` rather than `FALSE` was made on evidence rather than taste: a grep found
that `tests/testthat/test-16-compiled-matrix.R:169-170` reads the likelihood
from a default `aci()` call through `keep = "paths"`. That is a tested
exposure, so the saving is opt-in.

With `loglik = FALSE`: the dyad warm filter drops 54.6% (1.879 to 0.853 ms),
the ENSO warm filter 30.6% (197.8 to 137.3 ms), and the complete dyad `aci()`
call 36.2%. The flag itself is free when on: all four measured rows move by
less than 0.8%, inside the quartile spread.

Gate: 59 of 59 quantities bit-identical at default settings; 29 of 29
bit-identical against the pre-change reference in the `loglik = FALSE` arm.

Recorded as still available and not taken: the smoother and lag-table entry
points still build their internal filter with the likelihood on and then
discard it, which is a free saving of the same order. Still not taken at 0.1.0.

### 4.3 RNG containment on seeded simulation (C1c, from `aciR`)

`simulate()` with a `seed` now saves and restores the caller's `.Random.seed`,
materialising it first if the session has never drawn. `aciR` already had this
contract. Before the change, all six measured seeded arms left the caller's
stream advanced, and a call rejected for an invalid argument also advanced it,
because seeding happened before validation.

Gate: 64 of 64 quantities bit-identical, including the seeded paths themselves.
Overhead on a 3001-point simulation is below the noise floor of the
measurement.

### 4.4 Matrix recursion arithmetic (C1d, organisation from `aciR`, adopted at bit identity)

This is the change that closed the matrix-path efficiency gap. What was adopted
is `aciR`'s organisation of the recursions, not its code: per-step slicing off
bound components rather than rebuilding a list of matrix copies at every step,
one factorisation against the filtered covariance rather than a re-validated
solve, and a Cholesky-only metric. The bundle interface, the conditioning
consumption, the likelihood flag, substepping, the implicit stepper and the
covariance policy are unchanged.

Profiling identified the three targets precisely: on the smoother, the
per-step coefficient reader was 25.0% of total time, the covariance solve
26.7%, and the floor 20.7%; on the KL path, the strict Cholesky wrapper was
58.8%, of which a single coercion call was 19% self time.

| Row | `aciR` 0.2.3 | `aci` 0.0.30 before (timing study) | `acir` after | Speedup |
|---|---:|---:|---:|---:|
| ENSO warm smoother | 101.06 ms [a] | 182.98 ms | 80.79 ms | 2.27x |
| ENSO warm matrix KL, total only | 84.23 ms [a] | 137.28 ms | 69.13 ms | 1.99x |
| ENSO warm filter, likelihood on | 70.23 ms [b] | 200.98 ms | 93.83 ms | 2.14x |

`aciR` 0.2.3 column: this entry's own same-process re-measurement, not quoted
from the timing study. [a] Matched contract, no likelihood on either side:
`aci_smoother()` and `aci_metric()` with `S_xoS_x_inv` supplied. [b]
`aci_filter()` with `S_xoS_x_inv` supplied, and the contract is mismatched:
`aciR` computes no predictive likelihood anywhere, so its filter is timed
without one against a row that has it on.

Measured in the same process against `aciR` 0.2.3 rather than quoted from the
timing study, the three contracts that `aciR` previously won are now reversed:
warm smoother 80.79 against 101.06 ms, matrix KL 69.13 against 84.23 ms, warm
filter without likelihood 46.43 against 70.23 ms. The `aciR` rows measure 4 to
5% slower today than in the timing study, so this machine is slightly slower,
not faster, and the gains are if anything understated.

Gate: 86 of 86 quantities bit-identical, maximum absolute difference 0. A
1e-12 allowance was available and unused.

Recorded as the largest remaining win and deliberately not taken at the time:
the smoother still forms an explicit inverse of the filtered covariance, where
`aciR` avoids it. **That claim was later withdrawn on measurement**; see
section 4.10.

### 4.5 Vectorised affine coefficient realisation (C1e, target missed, honestly)

The 41.3x realisation row was the merge's largest single efficiency item. The
realiser was rewritten to bind closures and the grid once for the whole path,
fill each realised quantity by contiguous writes into one flat vector, apply
the affine differencing to a whole block at a time, and reduce the cross-noise
test to a single whole-path pass.

| Row | `aciR` 0.2.3 | `acir` before | `acir` after | Speedup |
|---|---:|---:|---:|---:|
| ENSO compile, the graded row | 12.56 ms [a] | 525.09 ms | 343.13 ms | 1.53x |
| ENSO affine realisation stage | 0.47 ms [b] | 474.47 ms | 297.10 ms | 1.60x |
| Complete public ENSO `aci()` | 297.26 ms [c] | 749.76 ms | 572.18 ms | 1.31x |

`aciR` 0.2.3 column. [a] and [b] are this entry's own in-session `aciR`
controls, measured in the same post-change run as the after column:
`aci_enso_components()` plus the `S_xoS_x_inv` construction for [a], and
`aci_enso_components()` alone for [b], which builds no Gram inverse. Both are
whole-path array construction with no per-point coefficient contract to honour,
which is the difference this entry diagnoses below. [c] The timing study, a
different entry: `aciR` has no single public complete call for a vector system,
so its row is a composed pipeline with no likelihood, no decomposition and no
provenance. This entry's own before column puts the machine within 3% of the
timing study's.

The stated target was 3x and was missed. The log records why, because the
diagnosis matters more than the number. After the change, the realiser's own
scaffolding is 11.7% of the compile and the model's per-point coefficient
closures are 74.8%. A whole-path vectorised transcription of the *same* ENSO
expressions over the same grid costs 0.18 ms. The remaining gap against `aciR`
is therefore not realisation waste at all: it is the per-point coefficient
contract, under which the constructor takes functions of one `(t, x)` point and
the affine construction must call them a fixed number of times per grid point
whatever the realiser does. `aciR`'s ENSO model is written as whole-path array
construction and has no such contract to honour. Closing the rest of that row
is a constructor contract change, not a realisation change. That diagnosis is
what C2d then acted on, and it is the reason this entry's honest miss was
productive rather than merely honest.

Two hypotheses were tested and killed here, which is worth recording so they
are not raised again. Provenance deep copies are not a factor: bundle assembly
including the whole provenance block is 0.05 ms of a 343 ms compile, 0.015%.
And halving the drift evaluations, which is available in principle because the
ENSO model computes one joint drift and subsets it twice, was measured at a
ceiling of about 2.0x rather than 3x and needs either a public API addition or
exactly the mutable model state the contract excludes. (C2d then took it for
free through a different mechanism; see 4.9.)

Gate: 564 of 564 recorded quantities bit-identical, including every realised
coefficient array, across the joint ENSO model, all three scalar partitions,
both conditioning strategies, three generic affine models including one with
correlated noise, a generic-closure control, and the directed dyad control.

Recorded tradeoff: the auto-generated zero cross-channels now memoise their
partner's channel width, which is state in a place that previously had none. It
does not change a result, and it makes one existing check vacuous and another
stricter; both directions are pinned by new tests.

The complete public ENSO `aci()` row moved again twice after this entry, and
the two later measurements are recorded in their own entries rather than
chained into a single cumulative figure: C2d measured 1104.33 to 494.00 ms
(2.24x) and 1209.48 to 512.05 ms (2.36x) in two cross-process A/B runs, and C3a
measured 282.64 to 259.64 ms (1.09x) within one process on a differently loaded
machine. The absolute milliseconds across those three entries are not
comparable to each other, and the log says so at each point.

### 4.6 Trusted supplied-filter smoother route (C2a, from the timing study)

The worst public-boundary row in the timing study was `aci`'s supplied-filter
smoother: 30.38 ms against a 0.684 ms kernel on the dyad, a factor of 44, and
769.3 against 183.0 ms on ENSO. Decomposed, 97.7% of the dyad call was
re-validating the filter path that the same namespace had just produced.

`aci_filter()` now seals the path it returns with a token, and `aci_smoother()`
skips the per-step re-derivation when the path still authenticates against the
run being smoothed. Every other skipped check has a live counterpart that is
re-established rather than trusted: grid, dimensions, observation provenance,
non-target specification and model provenance are all checked against the
authoritative live copy. The one genuinely retired check is the per-step
Cholesky, which the filter kernels guarantee by construction because they floor
every covariance. Mean finiteness is the one precondition the kernels do not
guarantee, so a path whose moment sums are not finite refuses to be sealed and
goes through full validation.

| Row | `aciR` 0.2.3 | `acir` before | `acir` after | Speedup |
|---|---:|---:|---:|---:|
| Dyad `aci_smoother()`, supplied filter | 2.178 ms [a] | 27.975 ms | 0.833 ms | 33.6x |
| ENSO joint, supplied filter | 122.589 ms [b] | 472.661 ms | 421.097 ms | 1.12x |

`aciR` 0.2.3 column: the timing study, same machine, contract as stated, and a
different entry from the before and after columns, which are this entry's own
measurements. [a] `aci_smoother()` with a supplied filter, the pair the timing
study records as matched against `aci`'s public smoother. [b] The same call on
ENSO, where `aciR`'s validation rebuilds the Gram inverse its filter has
already built. `aciR` returns bare lists on both rows, with no path metadata or
provenance.

Of the 0.833 ms, 0.680 ms is the smoother kernel itself and the authentication
is 0.016 ms. The ENSO row moves only 1.12x because validation was never its
problem: it is dominated by a recompile, which is the same per-point coefficient
contract cost identified in 4.5 and closed in 4.9.

Gate: 174 of 178 quantities bit-identical; the four that differ are object-weight
rows, not numerical results. Nineteen rejection conditions were captured as
class plus message and all nineteen are unchanged, including the two that only
the validation can raise. Authentication failure was proved for thirteen
distinct tampering routes, each falling back to the unchanged validated path.

Two design deviations were made against this entry's plan and both are
evidence-backed in the log: the token is a plain data list rather than an
environment, and the compiled bundle is deliberately not carried inside it.

Recorded costs, in the log's own words rather than softened: a sealed filter
path serialises 35.9% larger on the dyad and 55.0% larger on ENSO, though it
costs no memory because the token holds references. Those bytes land on users
who archive paths for inspection, and a path that has been through
serialisation cannot be supplied to the smoother anyway, before or after this
change. Carrying the compiled bundle in the token would remove the ENSO
recompile entirely, was designed for, and was rejected on measurement: it would
grow every filter path from 682 KB to about 3.4 MB, always, to speed up one
call pattern.

### 4.7 Forward CIR: the reporting window, the exact objective, and the named read-outs (C2b)

This entry resolves the CIR route cost and all three naming hazards of section
6.3, in two commits (`2da162c` engine, `577bc21` estimands).

**The window.** `aci_range()` computed every anchor of the record and the
caller discarded all but the ones asked for. It now takes `anchors`, forms each
requested row once, reads every quantity off it, and skips every per-interval
primitive before the earliest requested anchor.

| Row | `aciR` 0.2.3 | `acir` before | `acir` after | Ratio |
|---|---:|---:|---:|---:|
| 15 fixture anchors, streamed | 0.0 s [a] | 41 621.6 ms | 719.2 ms | 57.9x |
| 15 fixture anchors, public route | 0.0 s [a] | not available | 713.1 ms | |
| objective reduction, one 2001-cell row, exact | . [b] | 0.042 ms | 0.030 ms | 1.40x |
| all 2001 anchors, streamed (control) | . [c] | 41 974.9 ms | 43 710.0 ms | 0.96x |

`aciR` 0.2.3 column. [a] The dyad online/CIR comparison's route cost for the
same fifteen fixture anchors (interface finding 1), recorded to one decimal in
seconds and reported there as 0.0 s. It is a different entry from the before
and after columns, which are this entry's own measurements, so the comparison
is indicative of the route cost rather than a paired timing. `aciR` has one
windowed public route, and that single measurement stands against both `acir`
routes for this workload. [b] No `aciR` timing exists for a single-row
reduction. [c] `aciR`'s whole-record figure in the same comparison finding, 0.5
s, was measured over the 3001-point record rather than this row's 2001-point
one, so it is not the same workload and is not carried across.

The graded row is 57.9x against a 1 s target, and the public route reaches the
same 713 ms while reproducing the pinned MATLAB peak to 1.03e-14 and the
objective to 2.47e-15. Where the remaining 713 ms goes was measured: the
Theorem 3 smoother over the whole record 130.9 ms, the per-interval primitives
from the window start 148.8 ms, the fifteen rows themselves about 399 ms, and
compile and filter 1.8 ms.

**The exact objective.** The objective range is defined as the subjective range
averaged over every threshold. Both packages approximated that average and
neither had to: the count of cells above a threshold, integrated over the
threshold, is identically the sum of the row's suffix maximum, so the average
is `dt * sum(q) / M` with no quadrature error at all. The exact form is also
the cheaper one. Every value moves up by close to half a grid step, which is
the quadrature error the old form carried: over all 2000 finite anchors of the
oracle record the maximum absolute move is 6.603e-04, every one inside the 2e-3
gate, against `dt/2 = 5.0e-04`. The reference quadrature it replaces is
retained under its own name and graded against the independent `simps.m`
transcription at 3.43e-14 relative over 65 anchors.

**The estimand map**, which is the durable output of this entry and the answer
to hazards one and two of section 6.3:

| Quantity | `acir` | `aciR` 0.2.3 | `aci` 0.0.30 | Reference MATLAB |
|---|---|---|---|---|
| efficient ratio `dt * simps(row) / M` | `method = "l1_linf"`, `quadrature = "simpson"` | `objective` | `method = "l1_linf"` | `approx_objective_CIR` |
| eq. G.8 L1 grid sum `dt * sum(row) / M` | `method = "l1_linf"`, `quadrature = "sum"` | absent | same | FBCIR active lines |
| definitional threshold average, exact | `method = "exact"` (default) | absent | absent | absent, the limit both approximate |
| definitional average, threshold-grid quadrature | `method = "exact"`, `quadrature = "matlab_eps_grid"` | `objective_exact` | absent | `defn_objective_CIR` |
| definitional average, time-axis Simpson | retired | absent | `method = "exact"` | absent |
| subjective range, `index * dt` | `convention = "count"` (default) | `subjective` | absent | `subj_CIR_idx * dt` |
| subjective range, eq. G.7 lag time | `convention = "lag_time"` | absent | the only read-out | absent |
| reporting thresholds | `eps` | `epsilon` (overloaded) | `eps` | `eps_ord_values` |
| threshold quadrature nodes | `eps_grid` | `epsilon` (overloaded) | absent | `10.^flip(eps_ord_values)` |
| reporting window | `anchors` | `window` | absent | `first_idx:last_idx` |
| resolution status | `status` | `status` (with `margin`) | absent | absent |

Hazard three, the overloaded `epsilon`, is settled structurally: `eps_grid` is
refused by every mode that has no nodes to place, and `eps` is never read as a
quadrature grid.

**The status vocabulary** was adopted from `aciR` but its record-end criterion
was not. `aciR`'s `margin` is a caller-set fraction of the row. The shipped rule
is row-internal: it compares the record left after the last exceedance against
the exceedance itself, on the reasoning that a record that did not outlast the
influence it measures did not resolve it. On the fifteen fixture anchors this
gives 13 censored and 2 resolved, where `aciR` at `margin = 0.1` gives 11 and 4,
and the two orderings by anchor agree. Over all 2001 anchors: 1836 censored,
163 resolved, 2 insufficient.

Gates: suite 2368 pass / 0 fail / 0 error / 1 skip, from 2200. 197 of 268
recorded leaves bit-identical, with every one of the other 71 enumerated. Four
existing test files carry an intentionally changed output; no tolerance was
loosened anywhere and two were tightened to exact identities. The subjective
read-out change was verified on the grid index rather than on the difference:
`new == dt * k` and `old == dt * max(k - 1, 0)` hold bit-identically on all 8004
oracle cells, all 8004 cross-oracle cells and all 516 standard-set cells.

Recorded and not shipped: `lag_table()` has no window and remains quadratic, so
the fast path and the storable path are different objects and a caller who
wants both pays twice; the backward branch of the tau reducer is unreachable in
this release and therefore untested; `aciR`'s `margin`, `horizon`, `monotone`,
`subjective_censored` and `saturated` were not adopted.

### 4.8 Target by name, the MATLAB first-slice convention, and the two tau estimands (C2c)

Three API items, one of them with a numerical default change behind it.

**Target by name.** `aci_conditional()` gains `target`, mutually exclusive with
`given`, resolved through one shared resolver that derives the complement and
refuses a side that would leave the other empty. The MATLAB comments name
targets (`ENSO_model_cond_ACI_h_W_unobs.m:1199-1202`, `h_W(t) -> T_C |
(u,T_E,tau,I)`), so `target = "TC"` reads with the source where the complement
read against it. Graded as what it is: bit-identical on every route and every
quantity.

**The first-slice convention.** `aci_conditional()` gains `first_step`,
implemented at the single point where the masked precision is realised.
`"matlab"` writes slice 1 as the full Gram inverse and masks from slice 2. It is
refused with `method = "reduce"`, which has no masked precision path to apply
it to. The conditional comparison had recorded that `aci` had no injection
point for this at all, while `aciR` could express it only because it accepts a
caller-supplied inverse array; neither is an API. The `"matlab"` arm reproduces
that comparison's independent measurement to the digits it gave: step-2 filter
mean -1.079777623214551e-01 against -1.0798e-01, peak ACI 1.561914587488461e-01
against 1.5619e-01, filter covariance 3.65625e-02 against 3.6563e-02, and
time-integrated ACI moving 2.014773 to 1.986729, -1.3919%, against -1.392%.
Where that comparison said it must be inert, on the `u` partition, it is
asserted as an exact identity rather than a tolerance.

**The tau estimands.** `aci_enso_model()` gains `observations`, and `"reduced"`
is now the default for `hidden = "tau"`. This is the one intentional numerical
default change in the entry, and it moves because the previous default was not
what the reference script computes. Measured on the shared path, old against
new: filter mean up to 0.247021532 absolute (RMSD 0.0649, both series cross
zero so the relative column is uninformative), filter covariance up to
0.008491542 (posterior variance inflates 7.0% on average), and ACI up to
0.776486104, which is 3.06 times the mean ACI level, while the mean ACI itself
moves -0.51%. `likelihood_idx` moves from `1:5` to `1:3` and the compiled
observed dimension from 5 to 3, because the run now assimilates three channels.
`observations = "full"` restores the previous behaviour exactly.

The scalar-partition and conditional comparison figures those arms reproduce
were originally measured through `aci` and `aciR`; the numbers above were
measured in `acir`, on the same path, independently, and agree with the
scalar-partition stored columns to 5.1e-15, which is that file's write
precision.

Gates: suite 2484 pass / 0 fail / 0 error / 1 skip, from 2368, with all 116 new
assertions in one new test file and no existing test changed. 107 of 141
recorded leaves bit-identical; the 3 that changed are the one intentional
default. Overhead below 2% on every conditional row against an untouched control
that itself moved 0.38%.

Recorded and not shipped: a declared estimand cannot be composed with a
caller's specification, which is the conditional comparison's interface finding
A3 and is refused rather than silently resolved; `"reduced"` exists only for
`hidden = "tau"`, because that comparison verified the reduction is exactly
equivalent to masking for `h_W` and `u` and a name with no estimand behind it
is worse than no name; and the `T_C` script was still unrepresentable at this
point, recorded in `meta$unsupported_partitions` until C3h closed it.

### 4.9 Authenticated whole-path realisation for the ENSO benchmark (C2d)

This is the vectorised realiser's own conclusion (C1e), applied through its own
pattern. The coefficient closures of `aci_enso_model()` are built inside one
locked environment, a batch realiser is attached to the constructor as an
authenticated descriptor, and the whole-path route is taken only when the
declaration checks out structurally against the model's own dimensions and
state layout. Anything else realises point by point exactly as before.

The Gram-inverse path is in the same entry because it is what the graded row
becomes once the realiser leaves: realisation falls to 3.6 ms and the per-slice
inverse is then 94% of the compile, so a batch realiser on its own would have
reported an 8x improvement on a row that was still 94% something else.

| Row (one process, interleaved) | `aciR` 0.2.3 | `acir`, C1e route | `acir`, C2d route | Speedup |
|---|---:|---:|---:|---:|
| ENSO affine realisation | . [b] | 775.47 ms | 3.65 ms | 212.6x |
| ENSO compile, the graded row | 24.79 ms [a] | 813.84 ms | 47.26 ms | 17.2x |
| ENSO conditioned compile, mask, `hW` | . [b] | 573.42 ms | 56.39 ms | 10.2x |
| ENSO Gram inverses, plain | . [b] | 92.32 ms | 44.62 ms | 2.07x |
| ENSO Gram inverses, masked | . [b] | 126.82 ms | 52.49 ms | 2.42x |

`aciR` 0.2.3 column. [a] This entry's in-process `aciR` control,
`aci_enso_components()` plus the `S_xoS_x_inv` construction, timed in the same
process as both `acir` routes; it is the analogue the 1.9x claim below is made
against. The machine was loaded for the whole session, so the in-process ratio
is the reading and the absolute figure is not comparable across entries. [b] No
`aciR` control was measured for these rows in this entry: the one `aciR`
analogue here bundles realisation and the Gram inverse into a single call and
does not decompose into a realisation stage, a conditioned compile, or a
Gram-inverse stage.

The machine this ran on was measurably slow: an unrelated process held a core
at 99% for the whole session and every untouched control ran 1.7 to 2.2x slower
than the same row in C1e. The graded comparison is therefore made within one
process, where both routes exist in the shipped build and every pair is
asserted `identical()` before it is timed. The 30 ms target was measured at
47.26 ms and the log states plainly that the target was not met on a comparable
machine. What is machine-independent, because it is measured in the same
process: the ENSO compile is now 1.9x the cost of the `aciR` analogue that
produces the same per-step arrays plus the Gram inverse, against 27x at the
start of the entry. Cross-process A/B against an installed pre-change build
agrees and is recorded twice: complete public ENSO `aci()` 2.24x and 2.36x,
conditional 2.23x and 2.31x, the tau reduced compile 8.19x and 8.63x.

Gate: 1046 of 1048 recorded leaves bit-identical, with no tolerance used
anywhere. The two exceptions are a private descriptor id and one hole being
closed. Suite 3436 pass / 0 fail / 0 error / 1 skip, from 2484, with all 952 new
assertions in one new test file, no existing test changed and no tolerance
loosened.

**The locked-environment hole, which the harness found by accident.** On the
pre-change build, evaluating `p$r <- 99` in the environment behind
`aci_enso_model()`'s coefficient closures succeeded, and the first version of
the capture script silently poisoned its own reference model with `r = 99` for
every subsequent probe. Locking the environment closes that, and closing it is a
precondition for attaching a batch declaration at all: the declaration and the
per-point closures must read the same immutable parameters or the two routes can
diverge without any code changing. The tradeoff is recorded: anyone who was
reaching into `environment(model$fx)` to retune a parameter without rebuilding
the model can no longer do so, and the supported route,
`aci_enso_model(params = )`, is unchanged.

Also recorded: the batch declaration is `aci_enso_model()`-specific and does not
generalise, because it rests on two constructor facts (elementwise drift and
noise expressions, diagonal diffusion blocks) that are not checked at runtime.
That is exactly why the declaration sits behind an authenticated descriptor and
why the gate is `identical()` against the generic realiser rather than a
tolerance. A second model wanting this pays the same price: its own declaration,
its own identity gate. `chol2inv()` was measured here and deliberately not
taken, because on this model it is bit-identical only by virtue of a diagonal
Gram; it landed in C3a with its own budget and its own dense-Gram evidence.

### 4.10 The explicit inverses, and the signed-zero finding (C3a)

The two explicit inverses the package forms per step are now taken with one
`chol2inv()` on a factor it already has, instead of two triangular solves
against an identity right-hand side. At the dimensions this package works at
nothing is FLOP-bound: the cost is R call overhead and argument coercion, and
`backsolve()` and `forwardsolve()` coerce their arguments where `chol2inv()`
does not.

| Row (one process, interleaved) | `aciR` 0.2.3 | `acir` before | `acir` after | Speedup |
|---|---:|---:|---:|---:|
| precision path, ENSO joint, k = 3, 4000 slices | . [c] | 25.48 ms | 14.99 ms | 1.70x |
| precision path, ENSO `hW`, k = 5, 4000 slices | . [c] | 29.05 ms | 17.91 ms | 1.62x |
| warm matrix smoother, ENSO joint, l = 3 | 96.11 ms [a] | 84.88 ms | 72.83 ms | 1.17x |
| ENSO compile | 13.06 ms [a] | 28.70 ms | 17.68 ms | 1.62x |
| complete public ENSO `aci()` | 297.26 ms [b] | 282.64 ms | 259.64 ms | 1.09x |
| complete public ENSO `hW` `aci()`, l = 1 | . [c] | 209.35 ms | 185.61 ms | 1.13x |

`aciR` 0.2.3 column: the timing study, same machine, contract as stated. No
`aciR` control was measured inside this entry, and section 10 records that C3a
ran at a load average between 2.8 and 4.2, so the column carries the contract
rather than pairing with this row's own before and after. [a] `aci_smoother()`
with `S_xoS_x_inv` supplied, matched contract with no likelihood on either
side; and `aci_enso_components()` plus the `S_xoS_x_inv` construction, which
the timing study records as the analogue of compilation. [b] The composed
pipeline, since `aciR` has no single public complete call for a vector system,
with no likelihood, no decomposition and no provenance. [c] No `aciR` analogue
was measured: no entry times an isolated Gram-inverse or precision stage on the
`aciR` side, and none times an `aciR` scalar ENSO partition.

**The signed-zero finding.** The preparatory study concluded that on a diagonal
Gram, `chol2inv()` and two triangular solves against an identity produce the
same bits, and C2d had recorded `identical()` TRUE with maximum absolute
difference 0 on this build. Both statements are true and both miss the thing
that matters. LAPACK's `dpotri` writes `-0.0` into the off-diagonal zeros of a
diagonal Gram's inverse where the triangular route writes `+0.0`. `identical()`
at its default and `max(abs(a - b))` each report `0 == -0`, so neither
instrument could see it. The serialized bytes can, and they do: the realised
observation-Gram weight array of every ENSO partition with more than one
observed channel changed its bytes while comparing equal by every numeric test.
This is precisely the hazard C1e's own tradeoff note had named for a different
optimisation, in the abstract, two adoption rounds earlier.

Two things follow, and both shipped. The arithmetic carries an explicit `+ 0`,
which normalises `-0.0` to `+0.0`, is a no-op on every other double, and costs
0.260 ms of a 26.76 ms row, 1.0%. It is paid to keep the realised coefficient
arrays bitwise what they were. And the regression comparator was rewritten to
compare serialized bytes, package-wide, for this entry and every entry after
it. The log names the risk plainly: a reader who deletes the `+ 0` will see
`identical()` TRUE, `max(abs(.))` zero and a green suite except for one
byte-level assertion, which is the only thing standing between the package and
silently changing the bytes of every realised Gram weight array.

Gate, on the bytes: 263 leaves, 1 475 087 values. 230 of 263 bitwise identical,
including the entire ACI_code scope, every realised bundle over 22 compiles,
every filter moment and log-likelihood, every hidden-dimension-1 smoother, the
lag table and both forward-CIR methods. The precision-path change moves nothing
at all inside the scope. The 33 leaves that moved are all downstream of the
matrix smoother at hidden dimension above one, or on a synthetic non-diagonal
Gram built to exercise the route; the worst is 3.60e-15 relative against budgets
of 1e-14 absolute and 1e-13 relative, and eleven orders inside the 1e-6 fixture
gate. Suite 3900 pass / 0 fail / 0 error / 1 skip, from 3436, with 464 net new
assertions. Twenty-one existing assertions were rewritten, and none was
loosened without an exact route pin replacing it.

**A withdrawal.** The matrix-recursion entry's tradeoff note (C1d) had recorded
`aciR`'s `t(solve(R_f, t(B_j)))` organisation, which never forms the explicit
inverse, as the remaining smoother win the package was giving up. That claim is
withdrawn. Built and measured as a variant, it runs 0.75x to 0.77x, 24 to 25
percent slower, on all four smoother rows in two independent interleaved runs,
and 1.00x on the complete public `aci()`. The mechanism is in the per-call
table: the variant pays 4.63 microseconds for the triangular solve against the
transposed coupling block plus 7.92 for the vector solve, where the current
kernel pays 4.07 for one inverse it then uses twice, and `forwardsolve()` on a
vector right-hand side is the single most expensive call in the table, more
expensive than the same solve against a full matrix. It is the more accurate
route, by a factor of 1.3 to 2.5 on the solve residual, which is the textbook
result; but the accuracy it buys is 1e-15 on ACI, no better than what was
adopted, and it costs a quarter of the graded row. On this arithmetic, at these
shapes, in R, avoiding the explicit inverse is not an optimisation. Not
adopted, and no longer listed as an open win.

Recorded and not taken: the masked branch still inverts the target sub-block
with two triangular solves, so every conditioned run stays bit-identical and the
masked reference remains the definition it is checked against without a
tolerance.

### 4.11 The scalar log-likelihood kernel, and the memory-table withdrawal (C3b)

One line. `rate`, the observed drift, the coupling row and the Gram are read out
of the bundle's flat scalar view once, before the loop, instead of four list
lookups at every step. This is the hoist the filter and smoother kernels already
did; the likelihood kernel was the one that did not. The arithmetic, its
parenthesisation, the order of accumulation, the finiteness check, the floor call
and the classed error are character for character what they were.

| Row (one process, interleaved) | `aciR` 0.2.3 | `acir` before | `acir` after | Speedup |
|---|---:|---:|---:|---:|
| scalar log-likelihood kernel, dyad N = 3000 | . [b] | 1.100 ms | 0.451 ms | 2.44x |
| public dyad `aci_filter()`, N1 = 3001 | 1.591 ms [a] | 2.310 ms | 1.686 ms | 1.37x |

`aciR` 0.2.3 column: the timing study, same machine, contract as stated. No
`aciR` control was measured inside this entry, and section 10 records that C3b
ran at a load average between 2.8 and 4.2, so the column carries the contract
rather than pairing with this row's own before and after. [a] `aciR`'s
`aci_filter()`, public and validating, returning bare mean and covariance
vectors; the pair is mismatched because it computes no predictive likelihood
where `acir`'s does. [b] `aciR` computes no likelihood anywhere, so there is no
analogue of the kernel this row times.

The reach is exactly the dyad, because the flat scalar view is materialised only
for one observed and one hidden channel, so the ENSO partitions are a control on
the injection rather than a second case. Gate: 263 of 263 leaves bitwise
identical on the bytes. Suite 3900, unchanged, with no test added and none
changed, because the kernel's output is bit-identical and any new assertion
would be a second copy of one that already passes.

**The withdrawal of the timing study's memory table.** It reported the compiled
bundle at 3.04 MB on the ENSO benchmark with 0.50 MB of it provenance, and the
plan carried that into a work item: store provenance as identities rather than
deep copies. **That row is withdrawn and the interpretation behind it is
corrected.** It measured with `utils::object.size()`, which walks a structure
and counts a shared object once for every place it appears. Re-measured with an
allocator-aware tool that counts each object once, and confirmed by address
identity, the bundle's model, source model and their two provenance mirrors are
one object, the four observation slots are one object, and the time and
observation vectors are shared. `object.size()` over-reports the bundle by 391
704 B (41.5%) on the dyad and 932 632 B (28.9%) on ENSO. **The provenance list
costs 440 bytes**, on both cases. There are no provenance deep copies to
remove, because there are no provenance deep copies.

The identity-based replacement the work item implied was built and measured
anyway: replacing the five mirrored slots with recomputable fingerprints reports
a 224 KB `object.size` saving and is a 1960-byte real increase, roughly doubles
the cost of the provenance check, and would add a hashing dependency to a package
whose `Imports` is `stats, graphics, utils`. It also cannot work in principle,
because the five identity comparisons in the bundle validator are the only reads
of the provenance block and they compare pointers to objects some twenty other
sites need the bundle to keep anyway. Not adopted, and no longer an open work
item. The one real duplication the study did find, the flattened scalar view at
216 968 B and 39.3% of the dyad bundle, is a different row with a different
tradeoff: it costs the scalar kernels 7 to 10% to remove, so it was neither
adopted nor withdrawn, and it stays.

### 4.12 An exact positivity test for one-dimensional hidden states (C3c)

When the hidden dimension is 1, the Gaussian-path validator's per-slice strict
Cholesky loop is replaced by one vectorised finite-and-positive test over the
whole covariance path. Hidden dimension above one walks the same loop as before,
in the same order, with the same messages.

This is an equivalence, not an approximation, and the log proves it term by
term against the strict routine's own body: the dimension checks are
established by a check that precedes the loop, the symmetry test can never
fire because a 1x1 matrix equals its transpose exactly, symmetrising is the
identity, and the Cholesky succeeds if and only if the single value is finite
and strictly positive. That last line is the only substantive one, and it is
tested in both directions over sixteen edge values including `-0`, subnormals
and both infinities, plus a 450-value random sweep spanning eight orders of
magnitude, with the accept/reject decision required to equal `is.finite(v) &
v > 0` on every one. The error contract is part of the equivalence: the
branch reports the first failing index, with the same condition class and the
same message string, and the tests pin the class, the message and the index
at the first, an interior and the last slice.

| Row (one process, interleaved) | `aciR` 0.2.3 | `acir` before | `acir` after | Speedup |
|---|---:|---:|---:|---:|
| validator alone, dyad l = 1, N1 = 3001 | . [b] | 28.003 ms | 0.074 ms | 379.9x |
| dyad `aci_smoother(force_validate = TRUE)` | 2.178 ms [a] | 29.718 ms | 0.944 ms | 31.5x |
| dyad `lag_table(forward)`, N1 = 2001 | . [b] | 49 415 ms | 49 550 ms | 1.00x |
| control: ENSO joint `aci_smoother(force_validate = TRUE)`, l = 3 | 122.589 ms [a] | 122.435 ms | 122.136 ms | 1.00x |

`aciR` 0.2.3 column: the timing study, same machine, contract as stated, and a
different entry from the before and after columns, so the column carries the
contract rather than pairing with this row's own figures. [a] `aci_smoother()`
with a supplied filter, which validates on every public call and has no
`validate = FALSE` route; the timing study records the pair as matched against
`aci`'s public smoother, and `aciR` returns a bare list. [b] No `aciR`
analogue: it exposes no separate path validator, and it has no lag table.

The two flat rows are the honest reading of the predecessor entry's claim that
this would be "useful to `lag_table()` and `aci_range()` too". The change does
reach them, but one `lag_table()` on that record is 49 seconds, so 28 ms of
validation is 0.06% and disappears into the spread. The reach is real; the
benefit is confined to the row where validation is most of the call.

Gate: 263 of 263 leaves bitwise identical, including the five recorded classed
rejections captured as data. Suite 4407, from 3900, with 507 new assertions and
no existing test changed.

### 4.13 A strict covariance policy, with recorded opt-in flooring (C3d, decision D3)

Every implicit covariance floor inside a state recursion, a metric input or the
likelihood is now governed by one argument. `regularize = "none"` is the default
and stops the run with a classed condition naming the site, the grid index, the
time and the offending value; `regularize = "floor"` is the previous behaviour
and records every floor it takes in the result's metadata. Eleven governed sites
across the matrix and scalar filters, the smoothers, the metric inputs, the
Theorem 3 auxiliaries, the lag-table core and the streaming CIR reducer; four
public signatures gained the argument.

**The evidence for changing the default is a fire count, not an argument.** An
inventory of 24 implicit floor sites counted how often they act: **29 events in
4 371 482 floor-capable invocations across 47 runs**, and every one of the 29
was inside one of three deliberately broken stress probes. No authors-source,
source-derived or independent-transcription fixture (the evidence classes
defined in section 10) fires a floor anywhere. The
floor is therefore not smoothing rounding error in normal use. Where it does
act it can fabricate a number far from what the recursion produced: the decisive
measured example is a smoother variance of -1.07e13 replaced by 1e-12, in a run
that completed and returned an ACI path, with the only signal a warning about a
different quantity. The package was also inconsistent with itself, since the
public KL and the lag table's diagonal were already strict while the same
table's off-diagonal cells floored, and inconsistent with `aciR`, which has been
strict throughout. 0.1.0 is where the two stop disagreeing.

**A boundary was drawn and recorded rather than acted on.** The conditioning-Gram
inverses are deliberately not governed. They invert an observed-noise Gram, a
model input whose positive definiteness the model validator already contracts,
so the remedy for a failure there is to fix the model, not to reduce the time
step: a different error with a different message. The same inventory measured
**0 jitter-ladder entries in 226 413 invocations** at those sites. A test block
named for that boundary is what fails if it is ever moved by accident.

Gate, on the bytes: 147 of 147 leaves, 553 417 values, bitwise unchanged under
the default. The three firing probes were re-run with `regularize = "floor"` on
the new build and with no argument at all on the old one, and their ACI paths,
covariances, log-likelihood and warning classes are bitwise equal, which
re-establishes the prototype's proof on the integrated code rather than on a
shim. Suite 4407 to 5151, with 744 new assertions. The two test blocks the
inventory predicted would change behaviour are the only two that did, and they
failed on the first build at exactly the value predicted, -700.495, at exactly
the site named.

Timing: no measurable cost. Worst ratio 1.004 against a 2% budget, and two rows
faster (the implicit filter 0.953x, the lag table 0.973x) because both reach the
factorisation helper with a recorder, which skips a per-call finiteness scan on
the accepting path. Two structural details paid for that flat result and both
were found by measurement rather than inspection: the empty site frame is built
once at load, because constructing a zero-row data frame per kernel costs about
30 microseconds, which is 2% of a millisecond-scale dyad filter; and the scalar
guards keep two `if` statements rather than one `||`, because joining them reads
better and costs 1.2% of the log-likelihood kernel.

Recorded tradeoffs: it is a breaking change, and the release note says so.
Anyone running the explicit stepper on a coarse grid now gets an error where
they got a warning and a path; the migration is one argument, the error message
names it, and a session option restores the old default. Strict also throws a
partial result away: on one probe the first 138 points were usable and strict
returns none of them, with the condition carrying the index and time so a caller
can truncate and re-run. A partial-result contract would need a new return shape
and was not attempted.

### 4.14 Simulated trajectories carry the model's channel names (C3e)

`simulate()` now names the observed columns of the trajectory it returns, from
the metadata the built-in constructors already record, and `as_obs()` on a
simulation carries them through. The reason is a usability failure found while
writing the README: `aci_conditional()` resolves block and target names against
the observation column names, a simulation had none, and the single most likely
first workflow, simulate then run a conditional `aci()`, failed with "Named
blocks need named obs columns". The workaround was a three-argument call into a
lower-level constructor that a new user will not guess, and the README had to
use two different idioms in two adjacent examples to work around it.

The rule is "when available": a model whose recorded names are missing, the
wrong length, empty or duplicated gets unnamed columns exactly as before. The
names are ignored, never repaired, and no vocabulary is invented for a model
that declares none.

Gate: five models run twice, once on the named simulation and once on the same
matrix with the names stripped, compared on the bytes across filter, smoother,
ACI and its decomposition, the lag table and forward CIR: 60 of 60 leaves
bitwise equal. Suite 5151 to 5181. The visible consequences are recorded: a data
frame built from a simulated trajectory now labels the channels by name, and a
simulation object saved by the old build is no longer `identical()` to a new
one, because the observation matrix gained a dimension-names attribute. Every
number in them is the same, and the gate is what says so.

### 4.15 `aci_online()`, the public fixed-lag online smoother (C3f)

The one substantive assimilation surface in either source package that `acir`
had no answer to. The estimate at index `j` conditions on the observed record
through index `j + lag`, saturating at the end of the record. Nothing here
restates a numerical expression: the auxiliaries, the one-lag statistics and the
complete backward sweep are the existing Theorem 3 internals. What is new is a
two-stack aggregation queue that delivers a fixed lag at every anchor in time
linear in the record and independent of the lag.

Four design decisions, each with a measurement behind it. `lag` has no default,
because `aciR`'s default of infinity is the one value most likely to be mistaken
for the full smoother, which it is not. `aciR`'s tolerance argument is not
carried, because its purpose is to stop a product that has already underflowed
and with a lag-independent cost early stopping buys nothing: measured across
tolerances from 1e-18 to 1e-300 on a 2001-step record, no value and no effective
lag changed. The effective lag is a per-anchor integer vector rather than a
single scalar, because that is what tells the caller which estimates saturated
at the end of the record. And the returned path has kind `"online"`, not
`"smoother"`, which is what makes handing it to the lag table fail at the
existing kind check with a message that names the problem instead of reaching a
provenance comparison and failing there misleadingly. No new guard was needed;
the kind is the guard.

Accuracy of `acir`'s `aci_online()`, against every oracle available:

| Check | Quantity | Result |
|---|---|---|
| authors-source dyad reference, 75 pairs | mean / covariance | 5.33e-15 / 5.00e-16 |
| independent MATLAB cross reference, 75 pairs | mean / covariance | 6.22e-15 / 4.44e-16 |
| independent MATLAB multivariate reference, 45 pairs, l = 2 | mean / covariance | 6.66e-15 / 5.27e-16 |
| `lag = 0` against `aci_filter()` | mean, covariance | bitwise |
| `lag = Inf` against the complete backward sweep | mean, covariance | bitwise |
| window against direct accumulator, four lag scales | mean / covariance | 4.44e-15 / 8.88e-16 |
| `aciR::aci_online_smoother`, 9 lags | mean / covariance | 8.44e-15 / 1.33e-15 |
| off-by-one: against the prefix ending at `j + lag`, 18 pairs | absolute | 4.44e-16 |
| off-by-one: against the nearest wrong neighbours | absolute | 6.85e-03 and 3.37e-03 |

The last two rows are the pin that matters: the value agrees with the intended
prefix some thirteen orders more closely than with either neighbour, so the test
asserts the boundary rather than the magnitude.

Cost, and the claim that the cost is flat in the lag:

| Case | lag | `aciR` 0.2.3 | `acir` window | `acir` per-anchor accumulator |
|---|---:|---:|---:|---:|
| dyad l = 1, N+1 = 2001 | 1 | 1 ms [a] | 8.1 ms | 16 ms |
| | 100 | 8 ms [a] | 21.1 ms | 942 ms |
| | 2000 | 121 ms [a] | 22.4 ms | 9674 ms |
| ENSO l = 3, N+1 = 4001 | 1 | . [b] | 16.2 ms | 35 ms |
| | 100 | . [b] | 43.6 ms | 2023 ms |
| | 4000 | . [b] | 51.2 ms | 40 978 ms |

`aciR` 0.2.3 column. [a] `aciR::aci_online_smoother()` at a fixed lag over
every anchor of the same 2001-point dyad record, from the online-smoother
prototype study, recorded in seconds to three decimals as 0.001, 0.008 and
0.121. That study ran on a more loaded machine than the rows beside it, and its
own note is that absolute numbers are not comparable across the two runs while
ratios are, so the column carries the contract rather than pairing with the
window column. [b] No analogue: `aciR`'s route is scalar only, a vectorised
log-product over one hidden channel, and no `aciR` row was measured for a
three-dimensional hidden state.

A 2000-fold increase in the lag costs the window a factor of 2.8 and the
per-anchor route a factor of 605. On the workload this replaces, the package's
own oracle harness, 2.391 s becomes 0.209 s, a ratio of 0.087, against 0.088
measured independently in the prototype on a differently loaded machine.

**Scheme labels, the second half of the entry.** A smoother path now records
which of the two discretisations produced it, in a field with exactly two
values: the backward Euler ODE integration, and the Theorem 3 discrete
recursion. The finer implementation tag is unchanged and stays separate. Before
this, the metric reported one string from one vocabulary for a direct call and a
different string from a different vocabulary for a reused lag table, for a field
whose whole purpose is to say why the two can disagree. They disagree because
they are different discretisations of the same continuous-time object, and the
field now says so. The numbers are what they were.

**The scheme gap is not a constant**, and this is why the documentation refuses
`aciR`'s sentence that at infinite lag the online smoother "reproduces" the
smoother. Measured on ENSO at `dt = 0.005`, the ACI gap between the two schemes
is 0.104 at 401 points, 0.103 at 1001, 0.156 at 2001 and 0.483 at 4001, so it is
a growing quantity and must never be quoted as a constant. `aciR`'s own test file
asserts a first-order convergence rate, not an equality; its documentation says
otherwise, and the test is right.

Gate: 101 of 101 pre-existing leaves, 217 451 values, bitwise identical. Eleven
leaves moved and all eleven are label strings, each enumerated. Suite 5181 to
5279, with the one existing assertion that had to change moving to the new value
rather than being relaxed. A separate positivity survey found 0 non-positive-definite
slices in 30 014 accumulated covariances across three models and every lag from
1 to full.

Recorded and not shipped: the sparse `(j, lag)` grid route stays internal,
because `aciR`'s public shape is the fixed lag at every anchor and a second
public surface has no source-package precedent; the prototype's third route was
dropped as slower than the window at every lag measured; and the lag table still
recomputes exactly what the new auxiliary pass produces, which is de-duplication
in the hottest loop in the package and therefore a candidate for its own change
entry with its own gate, deferred past 0.1.0.

### 4.16 Source-derived oracle fixtures for the four scalar ENSO partitions (C3g, decision D7)

Test assets and test code only; no file under `R/` was touched and no exported
behaviour changed. Six new files: the driving path, four partition references
at 201 sampled indices, and a full-path summary reducing every graded series
over all 4001 steps. Before this entry the only ENSO fixture in the package was
the joint three-hidden pair, so the three scalar-hidden partitions were
ungraded and with them two dispatch shapes no other fixture reaches: the matrix
branch at five observed and one hidden channel, and the non-target reduction at
three observed and one hidden with its cross-Gram admissibility gate. C2c had
moved the tau default onto the reduced estimand on the strength of a
measurement that lived in the comparison scripts and in the log; nothing in the
test suite held the two tau estimands apart. It does now, and the size of the
gap between them is pinned so that a later change cannot quietly close it and
be read as a fix.

**The evidence class is stated in the manifest, and it is not the strongest
one.** These are source-derived, not authors-source. No MATLAB output was read,
generated or compared for any of them. Specifically ungraded, and named as such
in the manifest: the MATLAB stochastic realisation, since `simulate()` is
Euler-Maruyama against a mixed Euler/Milstein scheme; `aciR`'s coefficient
realisation, since `aciR` has no constructor for any of these partitions, so
the recorded 2.11e-15 corroboration ran `aciR`'s kernels on `aci`-realised
arrays through an adapter and does not close `aciR`'s own open grade; the
upstream first-slice conditional convention, since these fixtures grade the
unconditional observation; and the noise cross-covariance, which this model
sets to exactly zero.

Every number the preparation run had recorded was stale by four commits and all
of them were re-measured against the current tree before the copy. Nothing moved:
the sampled maximum absolute error stayed at 5.329071e-15 on each of the four
arms, the scaled full-record maximum at 4.729359e-15, the coefficient arrays
bit-identical across eight arrays, the filter and smoother moments bit-identical
with maximum difference exactly 0, and the driving signal regenerated from its
seed to the same SHA-256.

**A silent-failure hole was found and closed twice over.** The existing manifest
reader matched fixture names with a character class that excludes digits. All
sixteen existing fixture names are digit-free, so the omission had never
mattered; all six new names begin with a digit-bearing prefix. Measured: against
the new manifest the old class sees 0 of 6 entries, and in the counterfactual
where the two manifests had been merged it sees 16 of 22 while the pre-existing
set-equality assertion against the hand-maintained name list **would have
passed**, because both sides agree about a set of names neither of them pins.
Six fixtures would have shipped with no byte pinning, silently. Two safeguards
ship: the manifests stay separate, which is the structural fix, because they
describe different evidence classes; and one shared reader is used by both, with
a strict parser whose class admits digits and a deliberately permissive scan that
finds block headers by shape. The tests assert the two agree on both manifests
for both hash algorithms, that the new manifest returns exactly six entries all
bearing a digit, and that between them the two manifests account for every CSV
in the directory with no name pinned twice. The permissive scan is not
decoration: it is the only assertion in the pair that does not inherit the strict
reader's blind spots.

Gate: suite 5680 pass / 0 fail / 0 error / 1 skip, from 5279, with 401 new
expectations, no existing expectation changed and no tolerance loosened. 257
recorded leaves, 1 721 842 values, all bitwise identical. The grading gate is
1e-8 and the worst measured agreement is 5.329071e-15, a factor of 1.9e6 of
headroom, with the reasoning for that particular gate recorded: two decades
tighter than the rest of the directory because these fixtures grade a
same-project producer with no cross-language printing slack, and not tighter
than 1e-8 because 4001 sequential steps under a different BLAS could move the
accumulated sums with no code change, while a real regression would move a
coefficient far further. Verified from its new build-ignored home, the generator
regenerates all six shipped files byte for byte when driven by `aci` 0.0.30, and
differs in exactly one column by 9.992e-16 when driven by `acir`, which is the
documented round-off in the cancellation-resistant KL arithmetic of 4.1 and is
why the producer is deliberately not the package under test.

Recorded and not shipped: no grading register was rebuilt, although these
fixtures would supply four rows to one; sampling at every twentieth step is
closed rather than accepted, by the full-path summary that reduces every series
over all 4001 steps; and the fixture payload grows 28% gzipped.

### 4.17 The `T_C`-hidden zeroth-order partition, and the estimand the metadata was naming (C3h)

The last coverage gap in the 0.1.0 scientific scope. `ACI_code` has five ENSO
scripts and this is the fifth, the only one whose split is not exact: with the
central Pacific temperature hidden, its damping term is cubic in the hidden
state, and the reference script restores conditional linearity by a zeroth-order
Taylor expansion of the damping coefficient about zero. What it assimilates is a
two-channel reduced model with three of the six states supplied as prescribed
forcings.

The constructor admits `TC` as hidden only as the exact singleton set and only
with the approximation named explicitly; every other set containing it keeps the
pre-existing exact-split error. It does not route through the affine
construction machinery, because that machinery splits one drift and here the
drift used for inference is not the drift used for simulation. The prescribed
series are looked up by nearest index and clamped at the ends, never
interpolated, and the check that the observation grid is the forcing grid lives
where the model meets the observations rather than in the closure, because the
model validator probes at times that are off any regular grid and a strict
closure would make the constructor unusable. `simulate()` refuses on this model
and names the joint constructor instead, because falling through to the generic
path would integrate the linearised system and return a plausible-looking path
of a system nobody wrote down.

**The construction is the exact partitions' own algebra with one substitution**,
and that is asserted rather than asserted-of: the tests compare the drift
coefficients against the exact three-hidden model's drift evaluated at zero, at
five times over three states, so there is one algebra in the package and not two
transcriptions of it.

Evidence, and its class. The shipped chain was compared in one process against
an independent transcription of the reference script written for this purpose:
worst 6.328271e-15 over all seven graded series across four arms, with filter and
smoother moments at most 1.943e-16. Eight fixtures ship, graded at 1e-12 with a
worst per-file maximum of 9.325873e-15, a factor of 107 inside the gate. The
evidence class is source-derived from an independent transcription: one step
stronger than package-to-package agreement, one step weaker than an
authors-source fixture. No MATLAB was executed.

**The omitted thermocline term (decision D1), measured.** The reference
script's assimilation drift omits a thermocline term that its own simulator
applies, that the prescribed and observed thermocline depth makes available,
and that the sibling scripts carry in the corresponding rows. `acir` includes
the term, and `matlab_defect_compat = TRUE` reproduces the published script.
The two arms differ only in the mean channel: the filter and smoother
covariances are `identical()` between them, maximum difference exactly 0,
asserted as an identity. On the script's own fourteen-year analysis window the
time-integrated ACI is **13.603754 intended against 27.516198 literal,
+102.3%**, with a filter-mean maximum absolute difference of 0.105225129 and an
ACI maximum absolute difference of 2.775189422. On the whole shorter record the
same comparison is +26.3%. This is what makes that a measured decision rather
than a stylistic one: the omission is not a small bias, it doubles the reported
quantity on the window the script itself analyses.

**The metadata defect, found in passing and fixed here.** The ENSO constructor's
recorded conditioning fields named the complement of the channels the hidden
variable reaches, not the estimand the reference scripts compute. That recorded
complement is structurally inert: two of its rows in the observation-coupling
array are identically zero, so masking it is `identical()` to not masking at all,
with maximum difference exactly 0. A caller who built a conditioning
specification from the field therefore got the unconditional run under a
conditional name. The scripts mask down to a single target channel
(`ENSO_model_cond_ACI_h_W_unobs.m:1202`, `ENSO_model_cond_ACI_u_unobs.m:1205`),
and masking the script's complement instead moves the peak ACI from 2.6900 to
1.5657, a maximum absolute difference of 1.3241. The fields now record the
script's estimand with the source line cited on the object, and a note on the
object states that the fields are a record and not a default. Declaring them as
a default would change the numerical default of two partitions, which was not
this entry's scope and is named as a real divergence rather than a tidy-up.

Gates: suite **6236 pass / 0 fail / 0 error / 1 skip**, from 5680, and the +556
accounts exactly across five files. 285 of 285 comparable regression leaves
bitwise identical over 1 721 925 values; the 20 leaves that moved are all
character metadata, enumerated. Re-run after a check fix and an added example,
330 of 330 bitwise identical with nothing moved. `R CMD check`: `Status: OK`.

Recorded and not shipped: the shipped fixtures are the state-time regeneration
of the prototype's, not its bytes, because the reference stores its coefficient
arrays one step ahead after the first element while `acir` evaluates every
constructor's coefficients at state time. Shipping a second convention knob on
one constructor was judged not worth its surface, so the generator carries the
transcription unchanged and runs it on the shipped convention; the gap is
9.3e-04 on one filter mean and 2.57e-02 on one ACI, small but not round-off,
the coefficient fixture keeps both phase columns so the difference stays
visible, and the prototype's original artifacts remain in the evidence archive
as the record. Every prototype figure reproduces up to exactly that phase delta
and nothing else. Also recorded: this partition cannot grade the first-slice
convention, because both the mask and the convention are exactly inert here for
the same structural reason, asserted as a negative regression rather than left
implicit; `h_W` remains the only live conditional test in the package. And the
package's variance floor on one noise amplitude is carried where the MATLAB
expression has none, because without it the two-channel Gram is singular at the
validator's probe and the constructor aborts; the deviation is at most 3.7e-04,
is asserted as an exact identity, and is inert because that channel's coupling
row is zero.

### 4.18 Two vignettes, and the closing check (C3i)

Documentation and packaging only. Two vignettes, 414 lines and 23 chunks, and
899 lines and 43 chunks. Their drafts had been written and fully executed
against a commit that predates all of the C3 entries and two of the C2 entries,
so every one of nine adaptations below was forced by a surface that had moved
and each was re-derived against the current tree rather than transcribed:
target by name, the corrected estimand metadata, the tau default, the
first-slice switch, `aci_online()`, the covariance policy, the `T_C` partition,
the scheme labels, and the simulated channel names. The drafts carried exactly
one unevaluated chunk, showing a then-unlanded call; there is now no
unevaluated chunk in either file. **All 66 chunks run, and every number in the
two vignettes is computed at build time on the vignette's own record.**

One correction was made in the tree rather than argued: this entry's plan said
the CIR window had landed and asked for the window prose to become an executed
chunk. What landed was the reporting window; the lookahead window the draft
described is a different argument and is not in this build, which was verified
at HEAD by inspecting the function's formals and finding the internal window
hard-coded as a freeze rule. So the reporting window became the executed chunk,
measured live, and the closing subsection was retitled to say that the
lookahead cap is not in this build, keeping the naming trap it warns about.

Gate: 330 of 330 comparable leaves, 1 777 953 values, bitwise identical on the
bytes, which is what a documentation change should look like and is checked
rather than assumed. Suite 6236, unchanged, adding no test and changing none.
`R CMD check --no-manual` on the built tarball with vignettes, tests and
examples all run: **`Status: OK`, zero warnings, zero notes**, with all five
vignette-specific checks reading OK. The tarball ships the sources and the
rendered output and nothing from the reserve tree, the benchmarks, the fixture
generator home or the development records.

The build ignore file was audited rather than edited, including one entry that
suppresses the prebuilt vignette index: it was re-tested with vignettes actually
present, and installing the tarball still lists both vignettes with their correct
titles because the index is rebuilt at install time.

Vignette authorship was subsequently ruled on by the package author: the
introductory vignette's prose in `aci` is the author's own, multiply screened,
and the `acir` introduction must be based on that text with minimal prose
changes. This entry therefore stands as the integration of the drafts, and the
drafts are demoted to input for a rework of the introduction on the author's
own prose, presented as diffs. That rework is in progress at the time of
writing and is listed in section 9.

### 4.19 Measured and not taken

Recorded here so they are not proposed again without new evidence. Every one
was built or costed, not merely considered.

| Candidate | Why not | Entry |
|---|---|---|
| `aciR`'s inverse-free smoother organisation | Measured 24 to 25% slower at these shapes; the accuracy it buys is no better than what was adopted | C3a withdrawal |
| Provenance stored as identities rather than copies | The premise was an `object.size()` artefact; the provenance list costs 440 bytes and the replacement measured as a 1960-byte increase | C3b withdrawal |
| Dropping the flattened scalar bundle view | Real, 216 968 B and 39.3% of the dyad bundle, but costs the scalar kernels 7 to 10% | C3b |
| Carrying the compiled bundle in the filter token | Would remove the ENSO recompile, and would grow every filter path from 682 KB to about 3.4 MB, always | C2a |
| Halving the drift evaluations through memoisation | Ceiling measured at 2.0x and needs either a public API addition or the mutable model state the contract excludes; later taken for free by the batch declaration | C1e, then C2d |
| `chol2inv()` on the masked conditioning branch | Would move the masked assertions off exact identity for a row the timing does not name | C3a |
| `aciR`'s `margin`, `horizon`, `monotone`, `subjective_censored`, `saturated` | The plan named the status vocabulary, not those; each needs its own evidence and its own entry | C2b |
| `aciR`'s online-smoother tolerance argument | With a lag-independent cost, early stopping buys nothing; measured across 282 orders of magnitude with no value changing | C3f |
| A public sparse `(j, lag)` online API | No source-package precedent; the demand is in the package's own tests, where it is available | C3f |
| Composing a declared estimand with a caller's specification | Overriding would silently swap the estimand; accepting would need composition semantics nobody has specified | C2c |
| A batch declaration for the predator-prey model, or for generic affine models | The declaration rests on constructor facts that are not checked at runtime; a second model pays the same price, its own declaration and its own identity gate | C2d |
| Storing the largest fixture compressed | Would work, and no other fixture in the directory is compressed | C3g |

---

## 5. Defects found, and what was done about them

All of these were found by the comparison or by the merge, rather than looked
for. Each is recorded with the file and line that exhibits it.

### 5.1 In `aciR` 0.2.3

**The vector path cannot run any hidden-dimension-1 system.** `.aci_slice`
(`aciR/R/aci-core-mv.R:56`) uses a default-drop subscript, so an array with a
singleton non-time margin returns a bare vector rather than a matrix. At hidden
dimension one this hits the observation-coupling array, and the filter fails
with a non-conformable-arguments error from the linear algebra rather than a
contract message. The schema validator accepts the arrays first, so the failure
surfaces late.

This is not a corner case for this project. All three scalar ENSO partitions
and all three MATLAB conditional scripts are exactly this shape: a scalar
hidden state with a matrix observed process. Nothing in shipped `aciR`
exercises it, because its own vector path was only ever run at three-by-three.
The fix is one line, and it was validated as a bit-identical no-op on `aciR`'s
own shipped three-by-three path before being used for any comparison, with
`identical()` true for filter, smoother and metric.

Disposition: mandatory in the merged package. The scalar-partition family
cannot exist without it, and four of the shipped partition fixtures grade
exactly that shape.

**Silent time-grid desynchronisation.** `aciR` takes `dt` as a separate scalar
from the components. Components built on the fixture time grid and then
filtered with a different `dt` are accepted with no warning and produce a wrong
answer (filter mean 0.0288869 against the correct 0.0266807 at one measured
index). Components built on a genuinely non-uniform time grid and filtered with
one scalar `dt` are also accepted. `aci` cannot desynchronise, because it
derives `dt` from the trajectory, and it refuses a non-uniform grid outright.

Disposition: `aci`'s behaviour is adopted. This is the strongest single
interface argument in `aci`'s favour found anywhere in the comparison, and it
is why the merged components carry their grid.

**Docstring drift.** `aci_online_smoother`'s documentation says the infinite-lag
boundary "reproduces `aci_smoother()`" and calls the boundaries exact
identities, while `aciR/tests/testthat/test-online-smoother.R:108-112` asserts
the opposite and pins a first-order convergence rate. The test is right; the
measured gap on the ENSO record is 1.89e-02 in the mean, and section 4.15 shows
it growing with the record rather than converging. Disposition: the sentence is
not carried forward, and `aci_online()`'s documentation states the gap and its
growth instead.

### 5.2 In `aci` 0.0.30

These are the plumbing findings from the timing study and the dyad online/CIR
comparison. None is a numerical error; all are public-contract costs, and all
are now closed one way or another.

- The public supplied-filter smoother recompiled the bundle and re-validated
  every supplied covariance, a factor of 13.95 on the matched public pair.
  **Fixed** (4.6, and the recompile half in 4.9).
- Forward CIR had no window argument and computed every anchor, which is why
  fifteen fixture anchors cost 37.6 s. **Fixed** (4.7): 719 ms streamed, 713 ms
  through the public route.
- The predictive likelihood was mandatory inside the filter loop. **Fixed**
  (4.2), with the kernel itself a further 2.44x (4.11).
- The compiled bundle is dense and appeared to store provenance as deep copies
  (0.90 MB against 0.07 MB on the dyad). **Withdrawn** (4.11): the measurement
  was an `object.size()` artefact, the provenance list costs 440 bytes, and
  there are no deep copies to remove. Provenance was already not a compile-time
  cost either; that hypothesis was tested and killed in 4.5.
- There was no injection point for the MATLAB first-slice convention: the
  precision path was built internally, so the convention could not be expressed
  at all. **Fixed** (4.8): `first_step = c("uniform", "matlab")`.
- Reproducing the pinned CIR fixture required setting a global option to zero,
  because at the shipped default the routine truncates rows adaptively. **Still
  true**, and now visible: the option is named in the vignette and in the
  regression harness, and the shipped status vocabulary reports censoring
  explicitly rather than leaving it to a truncation default.
- There was no public accessor for fixed-lag online moments, so those
  quantities could only be compared through the KL functional (joint ENSO
  comparison). **Fixed** (4.15): `aci_online()` is exported and graded against
  three oracles.

### 5.3 In the MATLAB reference source

Three findings, all verified directly against the source in this workspace.

**`f_y_s` is used but never defined, in the tau script only.**
`ENSO_model_cond_ACI_tau_unobs.m` reads `f_y_s` at lines 1455 and 1484 inside
the online-smoother recursion. There is no assignment of `f_y_s` anywhere in
that file; the other four ENSO scripts all assign it from the windowed drift
(for example `ENSO_model_cond_ACI_h_W_unobs.m:1376`). The script therefore
errors if run standalone. The intended value is zero, because in that
partition the hidden drift has no additive term.

**The `T_C` assimilation drift omits a term its own simulator has.**
In `ENSO_model_cond_ACI_T_C_unobs.m`, the simulator advances `T_C` with a
`gamma_C * h_W` term at line 1124, while the assimilation model's hidden drift
at line 1151 omits it. The same script's line 1150 shows the deliberate
zeroth-order substitution that makes the system conditionally Gaussian, with
the quadratic coefficient evaluated at `T_C = 0`; the omission is separate from
that substitution. **The size of the omission is now measured** (4.17): on the
script's own fourteen-year analysis window the time-integrated ACI is 13.603754
with the term and 27.516198 without it, +102.3%, and the effect is confined to
the mean channel, with both arms' covariances bit-identical.

**Global versus windowed indexing in the same recursion.** In the same script,
the analysis window is extracted at lines 1373 to 1389, and every companion
quantity is used in its windowed form, but the hidden-drift coefficient is read
with a global index at lines 1457 and 1486. This was verified benign rather
than assumed benign: the window offset is `838 * 12 * k_dt = 1005600` indices
and the seasonal factor has index period `6 / dt = 1200`, so the offset is
exactly 838 whole seasonal cycles. Measured over 200001 indices, the
mis-indexed coefficient differs from the intended one by at most 5.76e-13
absolute, which is 9.60e-13 of the factor's own range, and is pure
floating-point argument reduction, reproduced in R.

**Disposition for all three (decision D1):** implement the intended behaviour
and document the divergence. The alternative of bug-compatible reproduction was
considered and rejected, because the first defect makes the script
non-executable rather than differently-valued, and the second is an omission
against the script's own simulator rather than a modelling choice. As shipped,
the published script remains reproducible under `matlab_defect_compat = TRUE`,
so the divergence is documented, measured and switchable rather than merely
asserted. Consulting the authors remains available and is listed in the open
items.

### 5.4 In `acir` itself, found during the merge

Four, all found by an instrument rather than by inspection, and all closed.

**A mutable coefficient environment behind a sealed constructor** (4.9). On the
pre-change build, assigning into the environment behind `aci_enso_model()`'s
coefficient closures succeeded, and the capture harness silently poisoned its own
reference model before the hole was written down as a test. Closed by locking the
environment, which is also a precondition for attaching a batch realiser.

**A byte-level divergence invisible to every numeric gate** (4.10). The
signed-zero finding, and the comparator upgrade it forced. This is the merge's
clearest demonstration that a gate can pass while a package's output changes.

**A fixture manifest reader that would have pinned nothing** (4.16). Its name
class excluded digits, all existing names were digit-free, and the pre-existing
set-equality assertion would have passed on a merged manifest while six fixtures
shipped unpinned. Closed structurally, by keeping the manifests separate, and by
fixing the reader with a counterfactual proof recorded.

**Constructor metadata naming an inert conditioning set** (4.17). The recorded
complement was structurally inert, so a caller who built a conditioning
specification from it got the unconditional run under a conditional name. Fixed,
with the source line cited on the object and the difference measured on both
sides.

---

## 6. The three estimand findings, in plain language

These are the findings most likely to be asked about, because in each case two
things that look like the same quantity are not. All three are now expressed in
the shipped API rather than only in a report.

### 6.1 The tau "operational shortcut" is not an identity

**What the MATLAB says.** The tau-hidden script observes three channels
`(T_C, T_E, I)` rather than the full five, treating the other two hidden-family
channels as known forcing, and states that the two forms yield the same results
for conditional ACI analysis.

**What was tested.** Both forms were run on the same path with the same realised
coefficients, using a validated instrument. The equivalence claim was tested
rather than assumed.

**What was found.** For `h_W` the claim is exactly true: every compared
quantity, filter mean and covariance, smoother mean and covariance, and ACI,
differs by exactly zero across 4001 steps. For tau it is false. Time-resolved
ACI differs by up to 0.776, which is 3.1 times its own mean level, with a
Pearson correlation of 0.905, while the time-averaged ACI agrees to within
0.51%. Filter accuracy against the true hidden path degrades 3.3% and the
posterior variance inflates 7.0%.

That pairing is the whole finding: the quantity the method reports is the
time-resolved curve, and the quantity that agrees is the time average. A
summary statistic that looks fine conceals a curve that does not.

**Why the two cases differ.** The reduction is exact when the prescribed
channels carry no information about the hidden variable, which happens when
their coupling rows are zero. For `h_W` they are. For tau they are not: both
prescribed drifts contain tau. Prescribing them reproduces their effect on the
retained channels' drift but discards their innovations, and the innovations
are where the tau information is. This is information loss, not
re-parameterisation.

**Why this is a strong result and where its limits are.** The `h_W` case is a
control measurement with a known answer, run through the same instrument, and
it returns exact zeros. That is what makes the tau result a refutation by
counterexample rather than a suspicious number. The magnitudes are specific to
this path; the structural conclusion is not.

**Disposition, as shipped.** The two forms are different estimands and the
package ships both, through `aci_enso_model(hidden = "tau", observations = )`,
with the reduced three-channel form as the fidelity default because that is what
the MATLAB script computes, and the full-observation form as the documented
alternative that reproduces the previous default exactly. The default moved in
C2c and the divergence was reproduced independently in `acir` on the same path,
agreeing with the original measurement to 1e-15. Since C3g the two estimands are
also held apart by pinned bytes, with the size of the gap between them recorded
so that a later change cannot quietly close it and be read as a fix. A fidelity
claim about the tau partition must say which observation set it reproduces.

**Related exact result.** On this model family the two conditioning strategies,
prescribed forcing and precision masking, are the same operator: results are
bit-identical. That holds because the observation-noise Gram is exactly
diagonal and the noise cross-term is exactly zero, so zeroing the non-target
rows of the inverse gives exactly the gain that deleting those channels gives.
The equivalence is model-specific rather than general, so both routes are kept
and the reduction keeps its cross-noise admissibility guard. It was re-measured
on the `T_C` partition too, at 2.776e-16 on ACI.

### 6.2 The first-slice conditional convention is a verification requirement

**What differs.** Both R packages mask every slice of the observation precision.
The MATLAB conditional scripts leave the first slice unmasked, then mask
everything after it.

**Why it is not cosmetic.** On the `h_W` path, measured by running both
conventions through the same package so that nothing else varies: the step-2
filter mean moves by 0.108 (a 75% relative change against a value of 0.144),
the filter covariance by 3.66e-02, the peak ACI by 0.574, and the
time-integrated ACI by -1.39%. The difference decays roughly one decade per
decade of steps and never reaches round-off inside a 4001-step record: 4000 of
4001 steps remain above 1e-6.

**Why it is that large.** Decomposing the information in the first update by
channel shows that 76% of it comes from the `T_E` channel, which the mask
removes at every later step. The unmasked first slice is therefore not a small
initialisation detail; it is one step of a materially different estimator.

**Where it does not matter.** For the `u` partition the effect is exactly zero,
because the observation-coupling array has exactly one non-zero row there, so
masking the others cannot change the gain. That also means the `u` case does
not exercise the mask at all. Anyone grading conditional fidelity on the `u`
script is grading nothing; `h_W` is the real test. The `T_C` partition added in
C3h is inert for the same structural reason, which was asserted there as a
negative regression rather than left implicit, so `h_W` remains the only live
conditional test in the package.

**Why it matters for scale.** In the MATLAB scripts' own million-step analysis
windows the transient is dead long before it can affect an integrated
quantity. On the record lengths that are affordable to run repeatedly in
testing, it dominates the early record. Both facts are true at once.

**Disposition (decision D6), as shipped.** Uniform masking is the default, with
`first_step = "matlab"` reproducing the unmasked first slice. It is implemented
at the single point where the masked precision is realised, and it is refused
with prescribed forcing, which has no masked precision path to apply it to. The
`"matlab"` arm reproduces the original independent measurement to the digits
that report gave (4.8), and is asserted as an exact identity where it must be
inert. The switch exists so that future MATLAB-derived conditional fixtures can
be matched, since no conditional fixture exists yet in either package and
`aciR`'s own ledger records its conditional grade as open.

### 6.3 The forward CIR functional map, and its naming hazards

Forward CIR reduces a divergence row to a small number of summary lengths.
Three separate hazards were quantified.

**Hazard one: the same word names different functionals.** `aciR`'s `objective`
is the same functional as `aci`'s `l1_linf` method, agreeing to 1.7e-15. It is
not `aci`'s *default* method, which is named `exact`. Mapping the two by name
pairs quantities that differ by up to a factor of 1.9 on this data, and all
fifteen fixture anchors are non-monotone, which is exactly the regime where the
two functionals separate. `aciR`'s own documentation states the precondition
under which they coincide; that precondition does not hold at any fixture
anchor.

**Hazard two: the same word names different quadratures.** `aciR`'s
`objective_exact` and `aci`'s `exact` are the same estimand integrated along
different axes: `aci` integrates in the lag variable against the suffix-maximum
envelope, `aciR` integrates in the threshold variable on the supplied grid.
Refining the threshold grid moves `aciR` toward the analytic limit, from
1.43e-02 at 129 points to 7.78e-04 at 4097 points, which is the signature of a
quadrature gap rather than an estimand gap. But it cannot close the gap to
`aci`'s `exact`, because that is a different axis: the irreducible residual is
about 5.2e-04 absolute and 8.3e-03 relative.

**Hazard three: one argument doing two jobs.** `aciR` uses a single `epsilon`
argument both as the reporting thresholds for the subjective read-out and as
the quadrature nodes for `objective_exact`. Honouring a pre-registered
four-threshold normalisation therefore makes its objective go negative
(-0.28677 against a correct value near 0.19150), because a four-node
logarithmic grid closed by a quadratic is not a usable quadrature. This is not
a defect in either package's arithmetic; it is a defect in the argument's
contract.

**One more convention, verified to the last bit.** The subjective read-out
differs by exactly one grid step between the packages: `aciR` follows the MATLAB
reference and counts cells, `aci` reports the lag time of the last exceedance.
Across 11996 cells the difference from one grid step is at most 3.34e-16, and
both report zero where nothing exceeds, so the alignment is a one-step offset
with a zero clamp. Both packages document their own choice correctly.

**Disposition, as shipped (4.7).** The merged API names three functionals
separately, splits the two epsilon roles into `eps` and `eps_grid` with each
refused where it has no meaning, and carries the read-out convention as
`convention = c("count", "lag_time")` with the MATLAB and `aciR` count
convention as the default (decision D4). The definitional objective is computed
by the exact layer-cake sum, which is both exact and cheaper than either
package's current quadrature, with the threshold-grid quadrature retained solely
as the MATLAB compatibility mode and graded against an independent transcription
of `simps.m` at 3.43e-14. The read-out change was verified on the grid index
rather than on the difference and holds bit-identically on more than sixteen
thousand cells, which makes the pinned fixture comparison direct for the first
time: both CIR oracles now match the pinned columns with no normalisation at all.

One correction belongs on the record, because an earlier framing overstated it:
computing the objective exactly does not change agreement with the existing
fixtures, which pin the peak, the efficient objective and the subjective
read-out only. It is exactness against the definition, not a fixture-score
improvement. It is, however, a real change to a published-looking number: any
value computed with the previous `method = "exact"` moves by close to half a
grid step and will not reproduce. Both are tabulated in the log for exactly
that reason.

**A related scheme split, present in both packages equally.** Both carry two
discretisations of the fully-informed posterior, a continuous backward equation
and a discrete recursion, and their public surfaces are crossed: what is public
in one is internal in the other. The two agree only to first order in the time
step, and the measured gap is the *same number* in both packages to all printed
digits (0.104276 in the ACI metric on the ENSO record). Any figure that crosses
the two schemes is a scheme difference, not a package difference, and is
labelled as such throughout the reports. The merged API names the scheme
explicitly, in a field with exactly two values (4.15), and the documentation
records that the gap grows with the record rather than converging: 0.104 at 401
points, 0.103 at 1001, 0.156 at 2001 and 0.483 at 4001.

---

## 7. What is deferred, and where it lives

Nothing was deleted from either source package: `aci/` and `aciR/` remain
untouched as evidence trees. Within the new package, everything excluded from
the 0.1.0 mainline is filed rather than discarded.

### 7.1 The filing system, and where it now lives

`reserve/` is organised by origin, not by convenience:

| Directory | Contents |
|---|---|
| `reserve/fbcir/` | FBCIR-derived code and tests, with `fbcir.patch` |
| `reserve/fbcir-paper/` | FBCIR-paper-backed with no MATLAB. Empty by finding |
| `reserve/enkbs/` | EnKBS-derived code and tests, with `enkbs.patch` |
| `reserve/enkbs-paper/` | EnKBS-paper-backed with no MATLAB. Empty by finding |
| `reserve/aci_code-future/` | ACI_code-derived but deferred. Category retired: empty of code by finding, and the gap it recorded is closed |
| `reserve/aci-paper/` | ACI-paper-backed extensions with no MATLAB backing |
| `reserve/paper-extremes/` | The fourth paper's family, no MATLAB |
| `reserve/extensions/` | Package-only infrastructure with neither paper nor MATLAB backing |

Two of those categories were added beyond the original plan, because the
partition exposed them: the extremes family belongs to a fourth paper with no
supplied MATLAB, and the formula interface and validation diagnostics have
neither paper nor MATLAB backing. Filing them under an existing family label
would have misrepresented their provenance.

Three of those categories hold no code, and the distinction matters: they are
empty *by finding*, not by omission. Every paper-only fragment in the FBCIR and
EnKBS families turned out to be an inseparable arm of a MATLAB-backed
constructor and travels with it; and every ACI_code-scoped block in `aci`
0.0.30 was retained in the mainline, so there was nothing to defer. The two
paper categories keep a README recording the finding. The third,
`reserve/aci_code-future/`, was where the one open completeness item was
tracked, the `T_C`-hidden zeroth-order case that existed in neither parent
package. **That item is now closed** (4.17), so the category is retired and no
such directory ships.

Tests, documentation blocks and fixtures follow their code into the reserve
directory rather than being stranded. 73 files in all.

**Where these records live, and why that changed.** The per-change engineering
log and the reserve tree were originally inside the package tree, build-ignored
so that they never reached the tarball. **Decision: development and process
artifacts do not live inside the package repository at all**, executed as
commit `3106967`. The rationale is that a build-ignore keeps material out of
the tarball but not out of the published history, and these records are project
material rather than package material. Both were moved to a sibling directory
outside the package; removed from the git index; backed by `.gitignore` entries
under both names so that later work cannot silently re-create them inside the
tree; and their build-ignore entries removed, leaving a `.Rbuildignore` whose
remaining entries need no explanation beyond the benchmarks directory, which
matches the practice in `aci` itself. The resulting diff against the close of
the C3 entries
is 77 files, every one of them a removal of a build-ignored development record
plus two build-ignore lines and one new ignore file. Nothing under `R/`,
`tests/`, `man/`, `vignettes/`, `DESCRIPTION` or `NAMESPACE` changed, which is
why the suite tally and the check status recorded at that close carry forward
to HEAD unchanged. Past commits still contain the files; a clean initial commit
at publication resolves the history.

### 7.2 The filing register

The filing register carries one line per item: what it is, where it came from
with file and line references, its category, the action taken (mainline,
reserved to a named directory, or dropped), why, and the evidence reference.

Counted programmatically from the file, the register holds ten tables and 97
data rows: 81 item rows (74 dispositions of migrated or excluded material and
7 recording things deliberately left verbatim) plus 16 rows across the
decision register and the open-items list. The action column of the 81: 20
mainline verbatim, 8 mainline with a recorded caveat, 1 not migrated, 46
reserved to a named family directory (including 3 pointer rows for
inseparable constructor arms), and 6 dropped.

Six rows carry a **dropped** action, five of them dropped outright:

| Item, from `aci` 0.0.30 | Evidence for dropping |
|---|---|
| `kl_increment` | Dead: zero callers in the source, tests or vignettes of `aci` 0.0.30 |
| `cir_table` | Self-disclaimed in the object it returns: the package's own text places it outside the supplied papers and MATLAB reference code |
| `empirical_kl` | Its nearest-neighbour estimator is an explicit not-implemented stub, and it has no in-package consumer |
| `truncation_profile` | Zero callers in source and zero in tests; can return with adaptive-lag work |
| An unreachable mode in the lag-table core | The only caller never passes it, and the string appears nowhere else |

The sixth, `projected_kl`, is marked dropped from the engine but is filed with
its test rather than discarded, because it is a hard dependency of the extremes
family and must return with it.

A separate section records items left verbatim in the mainline *although* dead
or stale, so that the next pass does not mistake them for oversights. Under the
extraction rule the acceptance bar was bit identity with `aci` 0.0.30, so
cosmetic cleanups were deliberately not made.

**One caveat on the register's currency.** Tables 9 and 10 describe the state
at the close of the extraction. Table 9's status column still reads
"adopted; not yet implemented" for D1, D3, D4, D6 and D7, all five of which are
implemented since, and table 10's first row records the `T_C` item as filed
nowhere because it existed in neither package, which C3h closed. Both tables
carry those corrections inline. The per-change engineering log is the current
authority on all six; the filing register is the authority on what was
migrated, reserved or dropped.

### 7.3 The patch mechanism, and what it does and does not deliver

`fbcir.patch` (1031 lines) and `enkbs.patch` (2610 lines) were regenerated
against the renamed tree at HEAD `3db97dd`, and both apply cleanly to a clean
mainline, install and pass. In isolation the FBCIR family
tests give 537 pass / 0 fail and the EnKBS family 247 / 0, the golden dyad
grading included; on their preview branches the full suite reads
6397 / 0 / 0 / 1 and 6483 / 0 / 0 / 1, which is +161 and +247 assertions over
the mainline's 6236. Non-stacking was re-verified on this tree rather than
carried over. The backward direction is no longer a verb of its own: it
re-enters as `aci_range(direction = "backward")`.

Three caveats are on the record rather than discovered later:

1. **The patches do not stack.** Each applies cleanly on its own; both edit the
   same four files and overlapping tests, and applying the second after the
   first is refused, plain and three-way. A manual merge is required when both
   families return.
2. **Patches rot as the mainline evolves.** Each release regenerates them. That
   is a stated release-checklist step, and it is the practical form of the
   architectural-fit guardrail: either keep the engine surface they touch
   stable, or knowingly pay the rebase. Eighteen adoption entries into the
   engine surface is exactly the situation that rule anticipates.
3. **A patch ships implementation plus existing tests, not completed
   verification.** The EnKBS patch carries a machine-precision dyad grade; the
   breadth of its verification across the other systems remains later work. The
   FBCIR patch's parity against the authors' own outputs still needs either
   MATLAB to run the supplied generator or an independent transcription of it,
   and neither was available during construction.

---

## 8. Decisions D1 to D8

Each was raised as a decision rather than resolved silently. All are recorded
with the rationale for the default, what the alternative would have cost, and
where the implementation landed.

**D1. The MATLAB source defects.** *Adopted: implement the intended behaviour
and document the divergence. Implemented in C3h.* The alternatives were
bug-compatible reproduction and consulting the authors first. Bug compatibility
is not meaningful for the first defect, since the script does not run at all in
that state; and for the second, reproducing an omission against the script's own
simulator would embed an error rather than a convention. As shipped, the
divergence is documented, the published script remains reproducible under
`matlab_defect_compat = TRUE`, and the size of the divergence is measured rather
than asserted: +102.3% on the time-integrated ACI over the script's own analysis
window, confined to the mean channel. Consulting the authors remains open and is
listed in section 9.

**D2. Backward CIR.** *Resolved by the filing system: excised from the mainline,
filed in `reserve/fbcir/` with a patch for immediate re-application.* The
earlier recommendation was to migrate it but leave it unexported, which is
lower churn. The filing system supersedes that: it keeps the public surface
scoped to the reference codebase, keeps the code recoverable in one command,
and avoids carrying roughly 300 lines of unreachable code through a release
whose acceptance bar is fidelity. One consequence is recorded in C2b: the
backward branch of the tau reducer is unreachable in this release and therefore
untested, and it is to be settled when backward CIR returns.

**D3. Covariance policy.** *Adopted: strict by default, with explicit opt-in
flooring. Implemented in C3d.* `aciR` refuses an invalid covariance; `aci`
floored it implicitly inside the recursion. The measured evidence is that the
floor is numerically inert on the paths tested (29 fires in 4 371 482
floor-capable invocations, all 29 inside deliberately broken stress probes), so
the decision is about what the package tells the user rather than about the
numbers. A silent floor converts a modelling error into a plausible answer, and
the decisive measured case is a smoother variance of -1.07e13 floored to 1e-12 in
a run that completed and returned an ACI path. It is a breaking change, it is in
the release notes, the migration is one argument, and a session option restores
the old behaviour exactly, bitwise.

**D4. Subjective read-out convention.** *Adopted: the MATLAB and `aciR` count
convention as the default, with the paper's lag-time convention available and
named. Implemented in C2b.* Both are correct; they answer slightly different
questions and sit exactly one grid step apart, verified to the last bit and then
verified again on the grid index across more than sixteen thousand cells. The
default follows the reference implementation because 0.1.0 is a fidelity
release, and the effect is that both CIR oracles now match the pinned columns
with no normalisation at all. The alternative is one named argument away and
retains its old values exactly.

**D5. API naming.** *Originally adopted: `aci`-style names; superseded by the
surface alignment.* Renaming at the birth of a new package is cheap and
renaming after a CRAN release is not, so this was raised explicitly rather than
defaulted into, and the original decision was not to spend effort renaming (new
surface followed the same convention, `da_online()` beside `da_filter()` and
`da_smooth()`). **Decision, superseding that one: the public surface is the
`aciR`-style `aci_*` interface.** The rationale is one consistent prefix
across the joint package, fixed before any release fixes the names. It landed
as one commit with outputs bitwise identical and the suite and check unchanged.
The names in this paragraph are the historical ones.

**D6. Conditional first-step convention.** *Adopted: uniform masking as the
default, plus a MATLAB-compatibility profile. Implemented in C2c.* Section 6.2
is the evidence. The compatibility profile is required by the completeness
criterion, because the conditional scripts are part of the reference codebase
and their convention is part of what they compute. It is refused with prescribed
forcing, which has no masked precision path to apply it to, and it is asserted as
exactly inert on the two partitions where it must be.

**D7. Evidence additions.** *Resolved by the completeness criterion: the
comparison's scalar-partition arrays become pinned source-derived fixtures.
Implemented in C3g.* The alternative was to record the pin gaps as accepted.
Since the arrays were already produced and verified during the comparison, and
are lossless for this model family because the relevant Gram is exactly diagonal
and the cross term exactly zero, pinning them costs almost nothing and closes
three open coverage rows. As shipped they are four arms rather than three, since
the two tau estimands are pinned separately, and the evidence class is stated as
source-derived rather than authors-source in the manifest and in the tests.

**D8. Name and attribution.** *Decision: the package name is `acir`, final.*
The rationale is that the name should be fixed before any release, since
renaming later is expensive. The CRAN case-insensitivity consequence against
the existing `aciR` stands and should be raised with the supervisor. The
collision also bit locally during construction: the working filesystem is
case-insensitive, so a directory named `acir` beside the existing `aciR` tree
*is* that directory. The package therefore sat one level down and installed to
a separate library, purely as a local workaround while both trees coexisted,
with no consequence for the published package. DESCRIPTION, LICENSE and
CITATION credit lines are drafted and await sign-off; the drafts are in the
accompanying metadata documents, which are not part of this published
record.

---

## 9. Open items

Recorded plainly rather than closed by assertion.

**Scientific and verification**

1. **No conditional fixture exists in either package.** `aciR`'s own ledger
   records its conditional grade as open, and `aci` has no conditional oracle
   either. The first-slice switch exists so that a future MATLAB-derived
   conditional fixture can be matched, but the fixture itself does not yet
   exist. Note that a fixture built on the `u` script would not close this,
   because that case is structurally inert, and neither would one built on the
   `T_C` script, which C3h showed is inert for the same reason. `h_W` is the only
   live conditional test in the package.
2. **The shared-oracle limitation.** The sixteen fixtures inherited from the
   parent packages are one oracle, copied byte for byte into both. Agreement
   between `aci` and `aciR` is useful and was measured, but it does not by itself
   establish fidelity to the MATLAB. The fourteen fixtures added in this build
   improve the position but do not close it: six are source-derived from the
   package's own constructors, and eight are source-derived from an independent
   transcription of one reference script. Independent grading of the remainder
   would need either MATLAB runs or further independent transcription.
3. **Simulator parity is a separate claim from inference parity.** Neither R
   package reproduces the mixed MATLAB simulation scheme path for path, and
   exact stochastic-path parity is impossible across the two random number
   generators in any case. Distributional parity is what is claimed, the ENSO
   constructor records `matlab_simulator_parity = FALSE` in its own metadata,
   and the package documentation says so.
4. **`aciR`'s scalar ENSO coefficient-realisation grade remains open** on its own
   ledger. The partition fixtures do not close it: `aciR` has no constructor for
   any of these partitions, so the corroboration recorded in the manifest ran
   `aciR`'s kernels on `aci`-realised arrays through an adapter and an in-process
   shim.
5. **FBCIR author-output parity** needs MATLAB to run the supplied generator, or
   an independent transcription of it. Not available in this workspace.
6. **EnKBS verification breadth** across the other benchmark systems is later
   work; only the dyad grade travels with the patch today.

**Build and release work still to do**

7. **The two reserve patches were regenerated and verified** against the
   renamed tree, with current figures in section 7.3; they need regenerating
   again only if the mainline surface moves.
8. **The filing register's last two tables record the state at the close of
   the extraction**: five decisions recorded there as "not yet implemented" are
   implemented, and the `T_C` row records an item that is now closed. Both are
   annotated inline in the register (section 7.2).
9. **The introductory vignette is being reworked** on the author's own
   human-written prose, which is the authorship policy for the vignettes; the
   version in the tree at HEAD is the integrated draft, fully executed and
   check-clean, not the final text.
10. **The version string is `0.0.0.9000`**, not `0.1.0`. Stamping the release
    version, and the release notes that go with it, is a release-time step that
    has not been taken.
11. **A stale header in the older fixture manifest** cites two test file names
    that did not survive the migration. Flagged three times across the C3
    entries and left untouched each time as out of scope; it needs a pass of
    its own.

**Decisions for the supervisor or the package author**

12. **Whether to consult the reference authors** about the two source defects
    before publishing a divergence note.
13. **The CRAN name collision** between `acir` and `aciR`, and whether a
    separate narrow `aciR` submission was ever wanted.
14. **Where the ensemble family eventually ships**: inside the same package
    with a documented engine boundary, or as a companion package depending on
    this one. The second is the literal reading of "engine package". Nothing in
    the next two releases depends on the answer.
15. **DESCRIPTION Title, Description, Authors@R, LICENSE and CITATION** await
    sign-off, including one contributor-role question that both files depend on
    and one optional courtesy notice for the third-party quadrature routine
    whose rule was independently reimplemented. Drafts and the specific
    judgement calls are in the accompanying metadata documents, which are not
    part of this published record.
16. **Whether the package's references topic is re-scoped.** It lists six
    papers, including two whose families are not in this release. The filing
    register records this as deliberate rather than an oversight, because those
    keys are still referenced by retained source comments and error strings.
    Re-scoping is an editorial call.

**Known residuals, cosmetic**

17. A stale error string on a branch that is now unreachable, an always-null
    field kept for behavioural fidelity with `aci`, and a partial-matching
    positional argument message. All are recorded in the filing register as
    deliberate, not as oversights, because the extraction's acceptance bar was
    bit identity. They are a later editorial pass. One further item of the same
    kind was retired during the C3 entries rather than left: the constructor
    metadata field that named an inert conditioning set was a correctness
    issue, not a cosmetic one, and it was fixed.

---

## 10. Provenance: what is measured, what is quoted, what is unverified

**Measured in this workspace.** Every number in sections 3, 4, 5.3, 5.4, 6 and
7.2, and every agreement figure in the comparison tables. The scripts, result
files and per-entry harnesses are in the evidence archive described at the end
of this document, one directory per study and one per change entry.

**How the numerical gates were graded, and when that changed.** The C1 and C2
entries were graded with `identical()` and absolute-difference comparisons on
each recorded leaf. From C3a onward every entry is graded on the serialized
bytes of every leaf, because a real change was found that both earlier
instruments reported as zero (section 4.10). Where a C1 or C2 entry claims bit
identity, that claim rests on the earlier instrument; the later byte-level
captures re-established the whole standard set (the fixed regression capture
described in section 2.3) at the bit and found it
unchanged, at 147, 257, 263, 305 and 330 leaves in successive entries, so
nothing in the earlier rounds is known to have slipped past the weaker gate.
That is a re-establishment, not a proof about the past, and it is stated as
such.

**Independent transcription.** Two, in different places. `ACI_code-main/simps.m`
was transcribed into R independently and used as a third reference for the
quadrature check, including its even-length closure rule; both packages reproduce
it at 5.7e-16 and 1.1e-15 on a row that does not need the closure, and 5.6e-15 and
2.6e-15 on a row that does. And `ENSO_model_cond_ACI_T_C_unobs.m` was
transcribed into R directly from the MATLAB source, before the partition was
implemented and independently of it, and the shipped chain is graded against
that transcription at 6.328271e-15 in one process on identical paths.

**Read but not run.** The MATLAB source. No MATLAB was executed at any
point. Every MATLAB claim in this document is a source-reading
claim, with file and line given, and the one quantitative MATLAB claim that is
not a transcription result (the benign indexing in section 5.3) was verified by
reproducing the arithmetic in R rather than by running the script.

**Evidence classes are kept distinct.** Authors-source fixtures,
source-derived harnesses, independent transcriptions, analytic identities,
package-to-package agreement, and behavioural tests are different kinds of
evidence and are not treated as interchangeable anywhere in this process. In
particular: package-to-package agreement is not MATLAB fidelity; a shared
fixture is not a second oracle; and a passing test is not a numerical grade.
The two fixture packages added in this build each state their own class in their
manifest, and the weaker of the two states in its own text what it does not
grade.

**Non-comparable items were recorded, not forced.** Each comparison report
carries a section listing what could not be compared and why: the signal and
dispersion decomposition, which one package does not expose for vector systems;
the discrete-scheme fixed-lag moments, reachable only through a functional; the
status vocabulary, which had no counterpart on the other side; the predictive
likelihood, which one package does not compute at all; and the endpoint
contracts, which differ by design. Where two APIs implement different
estimands, the semantics were mapped and the gap recorded rather than scored.
Three of those five gaps are closed in the merged package, by the public online
accessor, the status vocabulary and the optional likelihood.

**Reproducibility caveat on timings.** Single machine, single session, one path
per case. R 4.5.2 on aarch64-apple-darwin20, macOS 26.5.2, with the Accelerate
BLAS. Timings are indicative and rank implementations; they are not acceptance
thresholds, and no threshold exists in this merge. Machine load varied
substantially across the rounds and three entries say so explicitly and quote
ratios rather than absolute milliseconds for that reason: C2d measured its
untouched controls running 1.7 to 2.2x slower than the same rows in C1e, and
C3a and C3b ran at load averages between 2.8 and 4.2. Absolute milliseconds are
not comparable across entries; the within-entry, within-process ratios are, and
every such ratio was taken with both arms in one process and each pair asserted
equal before it was timed.

**The state this document describes.** The package tree at commit
`3106967`, read live rather than pinned. The per-change engineering log and the
filing register, read at the same point. The final counts in the header were
taken from that tree by direct enumeration, and the filing-register row count
in section 7.2 was regenerated from the file by a parser rather than carried
from either earlier figure. The suite tally (6236 pass / 0 fail / 0 error / 1
skip; the skip is a guarded external-oracle parity test) and the `R CMD check`
status (`Status: OK`, zero warnings, zero notes, with vignettes, tests and
examples run) were recorded at commit `5e5ee5e`; the diff from there to HEAD is
77 files, all of them removals of build-ignored development records plus two
build-ignore lines and one new ignore file, with nothing under `R/`, `tests/`,
`man/`, `vignettes/`, `DESCRIPTION` or `NAMESPACE` touched. That diff was
checked rather than assumed, which is why those two figures are quoted at HEAD.

---

### Evidence archive

The complete measurement records are retained with the project, outside the
package and outside this published set: every script, comparator, timing
harness and result file for the four numerical comparisons, the timing study,
the six preparation studies, and each of the eighteen change entries. Four
documents sit with them: the plan of record, the verified extraction map and
its fifteen cut points (the separability report), the verified scope partition
over constructors, tests, fixtures and assets in both packages (the
scope-partition report), and the build-state and key-issues register this
document is drawn from. With them sits the per-change engineering log, one
entry per adopted or rejected change with its gate, timing, tradeoffs and
verdict, eighteen entries, C1a to C3i. They are available on request; this
document is written to be read without them.

Two of those records are published. The filing register, one line per item with
origin, action, reason and evidence, ten tables and 97 data rows, ships as
`reserve/DISPOSITIONS.md`; and the reserve tree itself holds the filed code,
tests and fixtures for the deferred families, with a README per category and a
patch per MATLAB-backed family. Inside the package, `NEWS.md` is the
user-facing record of what changed, per entry.
