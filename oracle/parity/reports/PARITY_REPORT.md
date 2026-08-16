# MATLAB / aciR parity report

Generated from `oracle/parity/`, 2026-08-14, extended 2026-08-16. Scope: the
dyad path in full, predator-prey in full in both causal directions (filter,
smoother, metric, online smoother and causal influence range, reference side
and aciR side), the ENSO vector filter/smoother, and the reference side of all
four scalar-latent ENSO configurations.
MATLAB R2025b (25.2.0.3042426 Update 1), macOS arm64, 32 GB.

Reproduce with `tools/run_parity.R`; every number below comes from a file in
`reports/`.

---

## 1. What is being compared

The reference publishes **seven top-to-bottom scripts and no function
library** -- 12,339 lines, of which three files are callable. `simps.m` and
`legendUnq.m` are third-party FileExchange code; `progress_bar.m` is the
authors' own, and is a console progress indicator. None of the three carries
any of the method, so there was nothing to place beside `aci_filter()`.

> Corrected 2026-08-16. This paragraph previously said that only two files were
> callable and that both were third-party. `progress_bar.m` is a third, it is
> the authors' own, and all six model scripts call it. The conclusion is
> unchanged, because it does no numerical work.

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

> **156 of 156 outputs across nine profiles, maximum absolute difference 0.**
> 58 for the dyad, 58 for predator-prey (both causal directions, core and
> range), 8 for the ENSO vector filter and smoother, and 32 for the four
> scalar-latent ENSO configurations -- over **2,159 verbatim reference lines in
> 33 extracts**.

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

### 4b. Predator-prey, both causal directions (2026-08-16)

The predator-prey extracts and gate G1 were complete in both directions, but
aciR had never been run on that system at all. It has now been, on the
reference's own conventions, and the numbers are in `reports/parity_predprey.csv`.

The reference values are the captured workspaces' own -- the composed runner's,
since this script must not be run top to bottom (F6) -- so nothing was
recomputed to make the comparison. Tolerances were declared at the head of
`tools/compare_predprey.R` before the first run and are the dyad's classes:
`1e-12` for the linear recursions, `1e-11` for anything passing through a
logarithm, a division by a covariance or a quadrature.

> **28 of 28 quantities pass across the two directions.** Worst absolute
> difference among the graded rows 1.50e-12; tightest headroom **18.8x**, on
> direction one's ACI metric. No row returned `SUSPECT-EXACT` and no comparison
> was empty. Two of the 28 rows are the untested `aciR default -- by design`
> measurements, reported at infinite tolerance and excluded from those two
> figures.

| Quantity | dir1 max_abs (headroom) | dir2 max_abs (headroom) |
|---|---|---|
| `filter_mean` | 6.57e-14 (41x) | 6.39e-14 (179x) |
| `filter_cov` | 5.41e-16 (2037x) | 6.66e-15 (362x) |
| `smoother_mean` | 1.15e-14 (341x) | 5.51e-14 (208x) |
| `smoother_cov` | 8.33e-17 (12775x) | 4.88e-15 (438x) |
| `ACI_metric` | 6.56e-13 (**19x**) | 3.69e-14 (231x) |
| `E_j` / `F_j` | 5.55e-16 / 5.84e-16 | 5.11e-15 / 2.52e-15 |
| `online_mean` / `online_cov`, lag `Inf` | 3.20e-14 / 1.21e-15 | 1.10e-13 / 3.28e-14 |
| `cir_peak` | 1.21e-13 (46x) | 1.44e-13 (58x) |
| `cir_objective` | 3.30e-13 (27x) | 5.73e-13 (24x) |
| `cir_subjective`, 616,113 cells | 8.88e-16 (46387x) | 1.78e-15 (45824x) |
| `cir_objective_exact` | 3.30e-13 (42x) | 1.50e-12 (32x) |

`cir_objective_exact` -- the reference's `defn_objective_CIR`, F2's quantity --
is compared here for the first time inside the driver rather than by hand. The
subjective ranges agree to less than one part in `dt`, so the two
implementations select the *same index* at all 616,113 threshold-time cells in
each direction; the residual is round-off in `counts * dt`.

The CIR was set up on the reference's conventions throughout: `window` its own
`first_idx:last_idx`, `horizon = last_idx`, `margin` stood down to `1e-9`, and
its own 513-point threshold grid. The comparison is restricted to the region
the reference plots, and that region is not assumed -- the driver refuses to
grade unless its length equals that of the reference's own
`defn_objective_CIR`, which is defined over exactly that region.

**The horizon divergence of F1 is invisible here, and that is informative.**
On the dyad, running aciR's default whole-record horizon instead of the
reference's shifted `objective` by 4.27e-03 and 0.567. On predator-prey
direction one it shifts it by **3.30e-13** -- the same number as the
reference-conventions row, so the design difference is below round-off -- and on
direction two by 1.81e-04. The two runs are genuinely different objects: over
the reference's discarded guard region they differ by 0.246, and across the
full 2,201-time window 589 times are censored under the reference's horizon
against 438 under aciR's. Inside the reported region neither censors any.

The cause is measured, not inferred. At direction one's last reported time,
`j = 1201`, the divergence sequence peaks at 1.0019 and has fallen to 6.39e-13
by the reference's horizon, with the entire discarded 800-step tail no larger
than that. The reference's lookahead guard is 5 time units here against the
dyad's 0.6, and the influence range is short relative to it, so the truncation
cuts nothing that is still alive. F1's discrepancy is therefore a function of
how much live tail a horizon removes, not a standing bias -- which is the
opposite of what a single-system measurement would have suggested.

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
`h/12 * (-y0 + 8*y1 + 5*y2)` at equal spacing. `objective` went from **4.58e-09 to 1.37e-14**, the maximum over the full
reported region of both graded datasets.

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

### F7 -- Conditional ACI is enabled in THREE of the five ENSO scripts

**Severity: high. This finding has now been corrected three times, and only the
current form was measured rather than read.** It was first logged as "pending
extraction", then as "no executable counterpart upstream" on the evidence of
`u_h_W_tau_unobs.m` alone, then as "enabled in exactly one" on a reading of all
five. That third form was still wrong, and the reason it was wrong is the
reason it kept being wrong: the question was answered by looking at source
lines, and the distinguishing feature of those lines is a comment marker.

The current answer comes from the captured workspaces instead -- from the
observation-noise inverse each script's own filter actually consumed. Two
structures are possible and they are far apart: with the per-step pseudoinverse
live, `S_xoS_x_inv` varies from slice to slice and is dense; with it commented
out and one diagonal entry overwritten across all slices, every slice from the
second onward is the same matrix with exactly **one** non-zero entry. Measured
in `reports/enso_conditional_structure.csv`:

| Script | shape | non-zeros, slice 2 | slices 2..N identical | Verdict |
|---|---|---|---|---|
| `ENSO_..._h_W_unobs.m` | 3x3 | 1 | yes | **CONDITIONAL** -- `h_W(t) -> T_C \| (u,T_E,tau,I)` |
| `ENSO_..._T_C_unobs.m` | 2x2 | 1 | yes | **CONDITIONAL** -- `T_C(t) -> T_E \| (u,h_W,tau,I)` |
| `ENSO_..._u_unobs.m` | 3x3 | 1 | yes | **CONDITIONAL** -- `u(t) -> T_C \| (h_W,T_E,tau,I)` |
| `ENSO_..._tau_unobs.m` | 3x3 | 3 | no (spread 144.5) | unconditional |
| `ENSO_..._u_h_W_tau_unobs.m` | 3x3 | 3 | no (spread 144.5) | unconditional |

So three of the five ship with conditional ACI enabled, not one. The previous
version of this table was wrong about `T_C_unobs.m` and `u_unobs.m`.

The hazard is also worse than the previous wording described. In
`h_W_unobs.m` the live line carries three spaces of leading whitespace, which
is what the earlier reading fastened on -- "distinguished by one character of
indentation". In `T_C_unobs.m:1209` and `u_unobs.m:1205` there is **no
indentation at all**: the line sits in column one, between commented siblings,
under a banner that reads `THE FOLLOWING CODE IS FOR IMPLEMENTING THE
CONDITIONAL ACI FRAMEWORK ... USE THIS DEFINITION OF S_xoS_x_inv INSTEAD`. The
banner says "instead", the line beneath it is live, and there is nothing to
notice.

A user running all five scripts and comparing their ACI outputs gets three
answers to conditional questions and two to unconditional ones, with nothing in
the output to say which is which.

Two smaller reference defects were found in the same pass and are recorded
here rather than separately. `tau_unobs.m` carries the h_W script's guidance
comments verbatim, so its five commented alternatives all name **h_W** as the
latent variable in a script whose latent variable is tau -- harmless, since the
lines are inert, but it is the copy-paste that plausibly explains how three
files came to have a line uncommented. And `T_C_unobs.m`'s live line is
`1/sigma_E^2` where the other two use `1/sigma_C^2`, because its observation is
the two-dimensional `(T_E, I)`; a reader pattern-matching on `sigma_C` would
miss it.

`aci_conditional()` therefore has three executable counterparts upstream, not
one -- but see F9 for the interface difference that has to be handled before
any of them can be graded against it.

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

### F9 -- The reference's conditional filter leaves its first step unmasked, and aciR's does not

**Severity: high for anyone grading `aci_conditional()`. Found 2026-08-16 while
gating the four scalar-latent ENSO scripts, and measured rather than argued.**

The conditional scripts build their observation-noise inverse in two passes. A
pseudoinverse fills slice one -- `S_xoS_x_inv(:, :, 1) = pinv(S_x(:,:,1) *
S_x(:,:,1)')` -- and only then is the single target entry overwritten across
every slice. The per-step pseudoinverse that would refill the later slices is
commented out. **Slice one therefore keeps the pseudoinverse's other entries,
and the filter consumes slice one at `j = 2`.** Every step but the first
assimilates the target channel alone.

aciR's `aci_conditional()` masks the non-target block at **every** step, slice
one included. That is the more defensible reading of the construction -- the
whole point is to decline to condition on the other channels -- but it is not
what the reference computes.

Measured on `h_W_unobs.m` at `T = 100`, by running the *same* byte-verified
extract twice and changing only its input: `S_x(:, :, 1)` is diagonal in this
model, so zeroing its two non-target entries makes the slice-one pseudoinverse
vanish off the target block, which is exactly the weight aciR supplies.

| Quantity | difference at step 2 | worst | above 1e-6 until | above 1e-12 until |
|---|---|---|---|---|
| `filter_mean` | **0.2376** | 0.2376 (step 2) | step 6256 | step 15682 |
| `filter_cov` | 0.0366 | 0.0366 (step 2) | step 1313 | step 6008 |

The reference's own filter mean reaches 0.3246 in absolute value over the whole
record, so the first-step divergence is **73% of the largest value the filtered
state ever takes**. It is a transient -- the filter forgets its initial
condition -- but a slow one: some 6,000 of 20,001 steps pass before it falls
below `1e-6`, and 15,682 before `1e-12`. It is not negligible at any tolerance
this harness uses.

Neither side is wrong. The measured weights at slice one are `1/sigma_C^2` on
the target and up to 3076.9 -- that is `1/sigma_E^2` -- on a channel the
analysis says it is conditioning on. This is recorded as an interface
difference so that the eventual `aci_conditional()` grading is set up as a
like-for-like comparison from the start, rather than reporting a 0.238
disagreement and then having to explain it. The like-for-like route is to hand
aciR the reference's own `S_xoS_x_inv`, which the components schema already
admits directly.

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
| `noisy_predator_prey_model.m` (1,332 lines) | **complete, both directions** -- setup, filter, smoother, metric, online smoother, CIR, exact objective range. G1 58/58, both fixtures re-graded, **aciR vs reference 28/28** (§4b) |
| `ENSO_model_cond_ACI_u_h_W_tau_unobs.m` (2,103 lines) | **complete for the core** -- G1 8/8, aciR vs reference 5/5. ACI/CIR sections not extracted |
| `ENSO_model_cond_ACI_h_W_unobs.m` | **gated** -- filter + smoother extracted, G1 8/8 at maximum absolute difference 0. aciR side not run |
| `ENSO_model_cond_ACI_T_C_unobs.m` | **gated** -- G1 8/8, difference 0. aciR side not run |
| `ENSO_model_cond_ACI_tau_unobs.m` | **gated** -- G1 8/8, difference 0. aciR side not run |
| `ENSO_model_cond_ACI_u_unobs.m` | **gated** -- G1 8/8, difference 0. aciR side not run |

The four scalar-latent ENSO scripts were gated on 2026-08-16, closing the open
item below. Each needed its own input list, and the differences were not
cosmetic: `h_W` consumes `sigma_C` and `sigma_h`, `T_C` consumes `sigma_C` and
`sigma_E` and has a **two**-dimensional observation, `u` consumes `sigma_C` and
`sigma_u`, and `tau` consumes a time-varying `S_y` and no sigma at all. `L_x` is
constant in `h_W` and `T_C` and per-step in `tau` and `u`; `L_y` is the mirror.
No template covers that, which is exactly why the earlier attempt failed.

`aci_conditional` now has three executable counterparts upstream rather than
none (F7, corrected), though grading against them requires handling F9 first.
`aci_enso_*` and `aci_predprey_*` have rows for the core, and `aci_predprey_*`
now has them for the range as well.

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

Three items stood here on 2026-08-15. Two are closed and the third was already
contradicted by this document's own §8; all three are recorded below with what
replaced them.

**CLOSED 2026-08-16 -- the four scalar-latent ENSO scripts are gated.** Eight
extracts were added, `check_extraction()` passes on all 33 records in the
manifest, and G1 returns **32 of 32 outputs equal at a maximum absolute
difference of 0** across `h_W`, `T_C`, `tau` and `u`. The diagnosis in the old
wording was right -- each script needs its own input list -- and the
per-script differences are listed in §6. Gating them is what produced the F7
correction and F9.

**CLOSED 2026-08-16 -- aciR has been compared against the reference on
predator-prey.** 28 of 28 quantities in both causal directions, tightest
headroom 18.8x; see §4b.

**WITHDRAWN -- "the 4.58e-09 `objective` residual is unexplained".** This item
was stale when it was written: §8 records the residual as resolved at
1.37e-14, and F1 gives the cause (the Simpson closure, whose elimination was
the error F1 now documents). Removed rather than left standing, because a
reader meeting the two statements in one document cannot tell which is
current.

### Still open

**aciR has not been compared against the reference on any scalar-latent ENSO
script.** The reference side is complete -- G1 32/32 -- and `aci_conditional()`
now has three upstream counterparts (F7). Two things stand in the way, and both
are now precisely known rather than vague:

1. The vector dataset contract in `tools/dataset_mv.R` assumes a square
   configuration and reads noise as feedback matrices. These scripts are
   three-observed-by-one-latent (two-by-one for `T_C`) with the latent noise
   given as a scalar sigma. aciR's own validator takes that shape without
   complaint -- `L_x` is `n_x` by `n_y`, `L_y` is `n_y` by `n_y` -- so the work
   is in the bundle export, not in the package.
2. F9. A comparison set up naively would report a 0.238 disagreement at step 2
   decaying over some 6,000 steps, and that number is an interface difference,
   not a defect. The like-for-like route is to hand aciR the reference's own
   `S_xoS_x_inv`.

**The ENSO ACI and CIR sections are still outside the captured range.** The
`Sections` in `manifest/knobs.dcf` stop at the smoother for all five ENSO
profiles, so the ACI metric and the causal influence range have no reference
values captured for any of them.

**`reports/parity_predprey.csv` is not in the rendered table.**
`tools/render.R` reads `parity_scalar.csv` only. The predator-prey rows are a
separate file, following `parity_vector.csv`'s precedent, and the renderer has
not been extended to pick them up.


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

---

## Disclosure status — RELEASED 2026-08-14

The maintainer has released the hold. The two observations about the reference
implementation's own behaviour -- that `noisy_predator_prey_model.m` computes
direction one's smoother from direction two's filter when run top to bottom,
and that some of the five ENSO scripts ship with conditional ACI enabled --
may now be published with the rest of the project's findings.

Both are reported as properties of published code, reproducible from the
manifests in `oracle/parity/`, and neither is a defect when the scripts are
used as their own instructions direct. Contacting the authors remains a
separate decision.

**Amended 2026-08-16.** This paragraph said "exactly one of the five ENSO
scripts". It is three, measured from the captured workspaces rather than read
off the sources; see the corrected F7. The scope of the release is unchanged --
the same observation about the same files -- but the count in it was wrong, and
a released statement carrying a wrong number is worse than a held one.
