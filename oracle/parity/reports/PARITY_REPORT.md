# MATLAB / aciR parity report

Generated from `oracle/parity/`, 2026-08-14. Scope: the dyad path in full,
the predator-prey filter/smoother/metric in both causal directions, and the
ENSO vector filter/smoother.
MATLAB R2025b (25.2.0.3042426 Update 1), macOS arm64, 32 GB.

Reproduce with `tools/run_parity.R`; every number below comes from a file in
`reports/`.

---

## 1. What is being compared

The reference publishes **seven top-to-bottom scripts and no function
library** -- 12,339 lines, of which only `simps.m` and `legendUnq.m` are
callable, and both are third-party FileExchange code rather than the authors'.
There was nothing to place beside `aci_filter()`.

Seven kernels were therefore hoisted out of `dyad_interaction_model.m` into
callable functions. **343 lines of verbatim reference code, none retyped.**
Each generated file is a signature, a byte-exact slice between `BEGIN VERBATIM`
/ `END VERBATIM` markers, and `end`; `check_extraction()` re-reads the source
and fails on one byte of drift.

| Extract | Reference lines | aciR counterpart |
|---|---|---|
| `ref_dyad_setup_and_signals` | 48--125 | `aci_dyad_model`, `aci_dyad_components` |
| `ref_filter_scalar` | 131--172 | `aci_filter` |
| `ref_smoother_scalar` | 178--241 | `aci_smoother` **and** `.aci_online_aux` |
| `ref_aci_metric_scalar` | 305--308 | `aci_metric` |
| `ref_online_smoother_scalar` | 326--431 | `aci_online_smoother` |
| `ref_cir_scalar` | 433--546 | `aci_cir` |
| `ref_objective_cir_exact` | 549--561 | **absent in aciR** |

The pairing is not one-to-one: the reference computes the backward smoother and
the online smoother's auxiliaries `E_j`, `F_j` in a single loop, where aciR
factors the auxiliaries out.

## 2. Gate G1 -- do the extracts reproduce the script?

Each extract is called with the inputs the script itself used, and every output
must **equal** -- not approximate -- the value the script itself produced.

> **124 of 124 outputs across five profiles, maximum absolute difference 0.**
> 58 for the dyad, 58 for predator-prey (both causal directions, core and
> range), 8 for the ENSO vector filter and smoother -- over **1,681 verbatim
> reference lines in 25 extracts**.

## 3. The fixture re-grade

The committed fixture `aciR/inst/extdata/dyad_reference.csv` was produced by a
hand transcription in `oracle/aci_oracle_dyad.m`. That transcription and the R
implementation share an author, so an error common to both would be invisible
to every test the package owns.

The setup, signal, filter and smoother are all O(N), so this was settled at the
**published N = 30000** through the byte-verified extracts:

> **All 8 quantities over 301 rows, maximum absolute difference exactly 0.**
> The transcription introduced no error.

The same was then done for both predator-prey fixtures at the published
N = 12000, which matters more because that script cannot be run top to bottom
(see F6):

| Fixture | Direction | Agreement |
|---|---|---|
| `predprey_reference_predator_to_prey.csv` | dir1, x(t) -> y | **exactly 0** |
| `predprey_reference_prey_to_predator.csv` | dir2, y(t) -> x | 1.5e-13 |

dir1 is bit-identical; dir2 agrees only to round-off, so the transcription
reproduced dir1's expression tree exactly and dir2's only up to reassociation
over 12,000 steps. Benign -- and the informative part is the scale: had the
transcription conflated the two directions, the hazard that script's structure
creates, the disagreement would be of order one rather than 1e-13.

## 4. Core parity on arbitrary datasets

Two bundles. `dyad_reference_head` is the published signal truncated to 3001
steps. `arbitrary_cross_noise` has nothing to do with any model in the
reference -- a random-walk observation with trigonometric components and, most
importantly, a **non-zero noise cross-covariance** (`S_yoS_x = 0.5`), which the
published fixture cannot exercise because its `Sx_2` and `Sy_1` are both zero.

Criterion: `|a - b| <= atol + rtol*|b|`, `atol = 1e-12`.

> **26 of 26 quantities pass on both datasets.** Worst absolute difference
> 1.53e-14; tightest headroom 143x. Filter, smoother, ACI metric, `E_j`, `F_j`,
> and the online smoother at lags 0, 25 and full.

No row returned `SUSPECT-EXACT`, which is itself the evidence that both sides
genuinely execute: two implementations in two languages agreeing bit-for-bit
would mean one of them had not run.

## 5. Findings

### F1 -- The CIR horizon differs, and the package's own CIR oracle shared aciR's convention

**Severity: high. This is the finding the exercise was built to produce.**

For each reporting time the reference compares partially-informed posteriors
over `obs = j … last_idx` -- the reporting window plus its `lookahead_tolerance`
guard. aciR runs to the **end of the record**. aciR therefore sees more of the
tail and reports a larger objective integral and a longer subjective range.

Measured on `arbitrary_cross_noise`, restricted to the region the reference
actually reports:

| Quantity | Reported region | Reference guard region |
|---|---|---|
| `peak` | **9.28e-15** | 1.74e-02 |
| `objective` | 4.27e-03 | -- |
| `subjective` | 1.16e+00 | -- |

The peak agrees to machine precision; the integral and the range do not, and
the discrepancy is the right scale for a horizon gap of 700 steps.

On `dyad_reference_head` the record is short relative to the influence range
and the divergence is larger still (`peak` 6.57e-02) -- and aciR marks **every**
reported time saturated, returning `NA`, while the reference returns a
truncated number anyway. On this evidence aciR's behaviour is the more
defensible: the reference's value there is an artefact of where its matrix
happened to end.

**The part that matters for the package.** `oracle/aci_oracle_cir.m:162` loops
`for obs = j:(N+1)` -- to the end of the record, which is *aciR's* convention,
not the reference's. Every CIR grading the package owns was therefore taken
against a harness that shared the convention under test, and this divergence
was structurally invisible to it. The agreements previously reported for the
CIR (2.67e-15, 7.08e-15) are real, but they graded aciR against itself on this
axis.

**CLOSED.** `aci_cir()` takes `horizon`, defaulting to `NULL` = the whole
record, so existing behaviour is unchanged. The fully informed posterior
is still taken over the whole record even when the comparison sequence is
truncated, which is the part easiest to get wrong and is covered by a test
asserting the truncated row is a prefix of the untruncated one.

With the reference's horizon the subjective range agrees to **2.22e-16**,
against 1.16 time units at the default.

**The 4.58e-09 residual is resolved, and my elimination of the Simpson
hypothesis was wrong.** An external review re-tested it properly and the
closure *is* the cause. My argument had been that the CIR integrand decays to
~0 at its endpoint so the closure term is multiplied by nothing. That holds for
a comparison run to the end of the record, where the last entry is the
divergence of a posterior from itself and is exactly zero. At the **reference
horizon** the last entry is the divergence of a still-partial posterior from
the fully informed one, and it is not small. Measured at `j = 408`: full-record
row ends at 0.0000e+00, truncated row ends at 9.2541e-06, which is 2.7e-04 of
the peak.

The error was not the hypothesis. It was substituting an argument for a
measurement at the point where the argument was load-bearing.

`.aci_simpson()` now closes an odd interval count the way the reference does --
the quadratic through the last three samples, integrated over the last interval
-- derived rather than transcribed, and asserted to reduce to
`h/12 * (-y0 + 8*y1 + 5*y2)` at equal spacing. `objective` went from **4.58e-09
to 3.57e-15**.

The same change removed a latent defect: the old 3/8 panel assumed equal
spacing, and `objective_exact` integrates a logarithmic grid. Neither the
129-point default nor the reference's 513 triggered it, both having an even
interval count, but any even-length `epsilon` would have.

A by-product of that measurement, against the exact integral: aciR's 3/8
closure is consistently about three times more accurate than `simps.m`'s for
odd interval counts at n >= 6 (errors 0.00479 vs 0.00913 at n = 6, 0.00029 vs
0.00100 at n = 12).

### F2 -- aciR does not compute the exact objective CIR

The reference computes the objective CIR two ways: the efficient underestimate
(inside `ref_cir_scalar`) and the **definition** at line 561,
`defn_objective_CIR`, integrating the subjective ranges over the threshold grid
and normalising by the peak. aciR implements the underestimate and documents
the definition, but does not compute it.

This **corrects a claim I made in the previous session** that capability parity
with the reference was complete. It was not.

**CLOSED.** `aci_cir()` returns `objective_exact` alongside `objective`. It agrees with the reference's `defn_objective_CIR` to
**1.59e-14** once the horizons match.

Cheaper than the earlier estimate in this report, which was wrong on a checkable
point: `.aci_simpson()` already accepted an abscissa. It was never unit-spacing
only.

It is computed from the **unmasked** subjective ranges, deliberately: a
quadrature over the threshold grid needs every node, and masking one because
its range ran long would make the integral undefined rather than conservative.
The `saturated` flag still marks the time.

### F3 -- The reference's exact objective CIR silently evaluates to empty on a short window

`subjective_CIR(:, 1:end-lookahead_tolerance/dt)` produces an **empty range**
rather than an error when the reporting window is shorter than the lookahead,
so `defn_objective_CIR` becomes `1x0` with no diagnostic. Encountered while
building the reduced profiles, not at the reference's own settings.

Not a defect at publication settings, and not ours to fix. Guarded here by the
`Requires` conditions in `manifest/knobs.dcf`, checked after every capture. The
same class of failure was found and fixed **inside this harness** -- see F4.

### F4 -- Harness defect, fixed: an empty comparison passed vacuously

`max(numeric(0))` is `-Inf`, which sails under any threshold. A comparison in
which no element was resolved by both sides reported `ok` for having compared
nothing. Fixed: an empty comparison is now the verdict
`EMPTY -- nothing compared` and is a failure. It fired on the first run
afterwards, on `dyad_reference_head`'s objective CIR.

### F6 -- `noisy_predator_prey_model.m` cannot be run top to bottom

**Severity: high, and it is a property of the reference, not of aciR.**

The script carries fourteen blocks headed "RUN THIS CODE SECTION TO STUDY THE
CAUSAL RELATIONSHIP", seven per causal direction, and every one assigns the
same names:

```
225  filter_mean = dir1     ->   276  filter_mean = dir2   (dir1 destroyed, unused)
321  smoother dir1   <- consumes filter_mean, which is now dir2's
601  ACI metric dir1 <- built on that
```

So a top-to-bottom run computes direction one's smoother and ACI metric from
direction two's filter. Nothing raises; the numbers look plausible. The split
reaches into the model setup, where `f_x` is a scalar in one direction and an
array in the other under the same name.

This is not a defect when the script is used as documented -- its own banners
tell the reader to run one section at a time. It is a hazard for anyone
automating it, and the dyad's linear structure makes a whole-script capture
look like the general method when it is not.

**Handled**, not merely reported: each direction is captured from a runner
composed of that direction's declared line ranges, in file order, and every
profile declares a discriminating condition (`numel(f_x) == 1` holds only in
direction one) checked before anything is graded.

### F7 -- Conditional ACI is enabled in exactly one of the five ENSO scripts

**This finding has been corrected twice and the current form is the measured
one.** First logged as "pending extraction", then as "no executable counterpart
upstream" on the evidence of `u_h_W_tau_unobs.m` alone. Reading all five settles
it:

| Script | Conditional ACI |
|---|---|
| `ENSO_..._h_W_unobs.m` | **ENABLED** -- line 1202 is live |
| `ENSO_..._T_C_unobs.m` | commented, full ACI |
| `ENSO_..._tau_unobs.m` | commented, full ACI |
| `ENSO_..._u_unobs.m` | commented, full ACI |
| `ENSO_..._u_h_W_tau_unobs.m` | commented, full ACI |

In `h_W_unobs.m`, `S_xoS_x_inv(1, 1, :) = 1/sigma_C^2;` carries leading
whitespace rather than a comment marker, and sits among five comment lines that
read as guidance. It computes `h_W(t) -> T_C | (u,T_E,tau,I)` where its four
siblings compute unconditional ACI.

This is a **reproducibility hazard in the reference**, not in aciR: a user
running all five expecting comparable outputs gets four answers to one question
and one answer to another, distinguished by one character of indentation. It
also means `aci_conditional()` does have an executable counterpart to be graded
against -- which the earlier version of this finding denied.

The original wording follows, for the record.

`ENSO_model_cond_ACI_u_h_W_tau_unobs.m` lines 1214--1226 are **commented
guidance**, not code. They tell the reader which single `S_xoS_x_inv` line to
substitute by hand for each of three conditional analyses:

```matlab
% USE THIS DEFINITION OF S_xoS_x_inv INSTEAD FOR CONDITIONAL ACI.
%   S_xoS_x_inv(1, 1, :) = 1/sigma_C^2;      % (u,h_W,tau) -> T_C | (T_E,I)
%   S_xoS_x_inv(2, 2, :) = 1/sigma_E^2;      % ... -> T_E | (T_C,I)
%   sigma_I = sqrt(lambda .* (4-I) .* I);    % ... -> I   | (T_C,T_E)
```

So there is nothing to extract and nothing to grade against. `aci_conditional()`
is a generalisation of a documented manual edit, and its parity status is
**absent-upstream**, not pending. It remains graded only by aciR's own tests.

### F8 -- The reference pseudo-inverts where aciR refuses

The ENSO filter forms `S_xoS_x_inv(:, :, j) = pinv(S_x * S_x')`, a Moore-Penrose
pseudoinverse the authors chose "for stability concerns". aciR uses a Cholesky
inverse and its validator requires positive-definiteness at **every** slice.

On a singular slice the reference returns a pseudoinverse and continues; aciR
errors. This is not hypothetical -- the reference's own conditional-ACI guidance
computes `sigma_I = sqrt(lambda .* (4-I) .* I)`, which vanishes at `I = 0` and
`I = 4`, and then patches the result with
`S_xoS_x_inv(3, 3, sigma_I == 0) = 0`.

aciR is the stricter of the two and rejects `I` outside `(0, 4)` at construction.
That is a defensible choice and arguably the safer one, but it is a **behavioural
divergence**, not an implementation detail: a system the reference will analyse,
aciR will decline. Documented, not changed -- silently pseudo-inverting a
singular observation-noise Grammian is exactly the kind of quiet degradation
this package exists to avoid.

### F5 -- Interface differences, documented not repaired

- The reference hardcodes the initial covariance as `0.1` inside the filter
  range. `R0` is not reachable; the driver **refuses** a dataset declaring any
  other value rather than grading at the wrong one. aciR's `R0` argument is an
  extension.
- The reference parameterises noise by feedback matrices and forms the
  Grammians itself; aciR takes Grammians directly. aciR's surface is the wider
  one -- not every admissible Grammian pair is two feedback columns of that
  shape.
- The reference fixes 513 thresholds inside the CIR range; aciR defaults to
  129. The driver passes the reference's grid explicitly.

## 6. Not yet covered

Honest scope statement. The machinery is general; only the dyad has been
carried through it.

| Reference script | Status |
|---|---|
| `dyad_interaction_model.m` | **complete** -- 7 extracts, G1 58/58, parity run |
| `noisy_predator_prey_model.m` (1,332 lines) | **complete, both directions** -- setup, filter, smoother, metric, online smoother, CIR, exact objective range. G1 58/58, both fixtures re-graded |
| `ENSO_model_cond_ACI_u_h_W_tau_unobs.m` (2,103 lines) | **complete for the core** -- G1 8/8, aciR vs reference 5/5. ACI/CIR sections not extracted |
| 4 remaining `ENSO_model_cond_ACI_*.m` | **workspaces captured and validated**; extracts NOT verified -- see below |

Consequently `aci_conditional` has no parity row and never will have one from
this reference (F7); `aci_enso_*` and `aci_predprey_*` have rows for the core
but not for the range.

The ENSO comparison of **aciR against the reference** is not yet run, only the
extracts against the script (G1). It needs a vector form of the dataset
contract, which is scalar-only today -- and it cannot be done as a fixture
re-grade, because `oracle/aci_oracle_enso.m` deliberately generates its signal
by plain Euler-Maruyama rather than the reference's Milstein scheme. The right
test is to drive aciR with the reference's own captured signal. They are graded only by the existing transcription-based
fixtures -- which, per F1, is exactly the grading whose blind spot this
exercise demonstrated. Extending the manifest to those scripts is the next
step; the generator, gate and driver need no changes.


## 7. Open items

**The four scalar-latent ENSO scripts are captured but not gated.** Their
workspaces run and pass their validity conditions -- including the
`size(filter_mean, 1) == 1` discriminator that distinguishes them from the
vector configuration -- but the extracted filter and smoother do not yet call,
because each script's block needs its own input list derived individually
rather than by the shared template. `h_W_unobs.m` needs `sigma_C` for the
reason in F7; the others each want something different. The eight unverified
records were **removed from the manifest** rather than left in place claiming a
grading they do not have.

**The 4.58e-09 `objective` residual is unexplained.** The Simpson-closure
hypothesis was tested and eliminated. No replacement hypothesis is offered
here.

**aciR has not been compared against the reference on predator-prey.** The
extracts and G1 are complete for both directions, but the aciR-side comparison
was run only for the dyad (scalar) and ENSO (vector).


## 8. Second-round results (2026-08-14, after external review)

An external review of the tree corrected one finding of mine and raised eight
more. Its central correction is recorded in F1 above. The parity numbers that
changed:

| Comparison | Before | After |
|---|---|---|
| `cir_objective`, reference conventions | 4.58e-09 | **1.37e-14** |
| `cir_subjective`, reference conventions | 2.22e-16 | 2.22e-16 |
| `cir_peak`, reference conventions | 9.28e-15 | 9.28e-15 |
| Scalar core, both datasets | 26 / 26 | 26 / 26 |

**The CIR is graded for the first time.** Before `horizon` existed the harness
compared aciR's whole-record answer with the reference's truncated one and
reported the gap as a disagreement; those are two different questions. The
harness now grades on the reference's own conventions -- `horizon = last_idx`,
`margin` stood down, the reference's threshold grid -- and reports aciR's
default separately, with its magnitude, as a designed difference rather than a
failure.

Also changed, from the same review:

* A time whose range is not resolved returns a **right-censored lower bound**
  rather than `NA`, with a four-valued `status`. On the arbitrary dataset at
  the reference horizon this is 1049 finite bounds where there were 1049 holes.
* Every surface that still said the online smoother and CIR were "scalar-only"
  -- `DESCRIPTION`, `README`, the assumptions vignette, two oracle-manifest
  notes -- has been corrected. `DESCRIPTION` was CRAN-visible and denied a
  capability that shipped in `e4cdc41`.
* `aci_cir()` gained `print`, `summary`, `as.data.frame` and `plot`; a vector
  `filt` is now validated; a `$` partial match that let a test pass for the
  wrong reason is fixed and made un-reintroducible by
  `warnPartialMatchDollar` in the test setup.

Declined, with reasons, in `design/2026-08-14_reviewer_OUTPUT.md` order: S7
(parallel CIR -- 0.082 s, and `mclapply` is Unix-only), S3 (`pinv` escape --
the reviewer's own probe showed the *filter* fails at step 2 on the motivating
case, so the escape would not have rescued it), and the S6 algorithm change,
whose `O(n)` prefix-sum form requires `exp(+/- cum_log)` reaching `e^{+/-60}`
over a published-scale record and needs blocked rescaling to be safe.
