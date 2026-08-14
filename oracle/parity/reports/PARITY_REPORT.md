# MATLAB / aciR parity report

Generated from `oracle/parity/`, 2026-08-14. Scope: the scalar dyad path.
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

> **58 of 58 outputs, two profiles, maximum absolute difference 0.**

## 3. The fixture re-grade

The committed fixture `aciR/inst/extdata/dyad_reference.csv` was produced by a
hand transcription in `oracle/aci_oracle_dyad.m`. That transcription and the R
implementation share an author, so an error common to both would be invisible
to every test the package owns.

The setup, signal, filter and smoother are all O(N), so this was settled at the
**published N = 30000** through the byte-verified extracts:

> **All 8 quantities over 301 rows, maximum absolute difference exactly 0.**
> The transcription introduced no error.

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

**Recommended repair.** Add a `horizon` argument to `aci_cir()` defaulting to
the whole record, i.e. current behaviour. `.aci_cir_row()` must keep computing
the fully-informed posterior over the *whole* record -- the reference's `{end}`
is the full record even when its sequence is truncated -- and return only the
first `horizon - j + 1` entries. Not implemented here: it changes the signature
of an exported function and needs its own tests and grading. Roughly 15 lines
plus validation, documentation and tests.

### F2 -- aciR does not compute the exact objective CIR

The reference computes the objective CIR two ways: the efficient underestimate
(inside `ref_cir_scalar`) and the **definition** at line 561,
`defn_objective_CIR`, integrating the subjective ranges over the threshold grid
and normalising by the peak. aciR implements the underestimate and documents
the definition, but does not compute it.

This **corrects a claim I made in the previous session** that capability parity
with the reference was complete. It was not.

**Recommended repair.** `.aci_simpson()` currently assumes unit spacing;
the threshold grid is logarithmic, so it needs an abscissa-aware form. Then
`objective_exact[i] = simpson(epsilon_sorted, subjective[order(epsilon), i]) /
peak[i]`. Also not implemented here -- same reason, it adds a field to an
exported return value.

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
| `noisy_predator_prey_model.m` (1,332 lines) | not extracted |
| 5 x `ENSO_model_cond_ACI_*.m` (~2,000 lines each) | not extracted |

Consequently `aci_conditional`, `aci_predprey_*` and `aci_enso_*` have **no
parity row yet**. They are graded only by the existing transcription-based
fixtures -- which, per F1, is exactly the grading whose blind spot this
exercise demonstrated. Extending the manifest to those scripts is the next
step; the generator, gate and driver need no changes.
