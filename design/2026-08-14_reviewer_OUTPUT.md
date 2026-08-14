# aciR — external reviewer report

**Package** 0.1.0.9000 · **Date** 2026-08-14 · **Scope** whole tree, with
empirical probes · **Status** reviewer notes, not a patch

Intended as the written counterpart of
`design/2026-08-14_reviewer_brief.md`. The brief asked for a one-to-one
reading against the authors' MATLAB, for effort on five named targets, and
for constructive suggestions on coding, mathematics, causal method, and
efficiency. Those five targets are answered first. Suggestions follow, ranked.
Evidence for the 4.58e-09 residual is in
`design/artefacts/2026-08-14_cir_residual_rows.csv`, produced by
`design/artefacts/2026-08-14_reviewer_probe.R`.

This report does not re-derive the CGNS filter, smoother, or ACI metric.
Those are already graded against the authors' code at machine precision, and
the brief was right that another derivation would not pay. It also does not
raise the two reference-implementation behaviours recorded in brief §6 with
anyone but the maintainer.

---

## 1. Verdict

The package is a careful, unusually well-evidenced reimplementation of
Assimilative Causal Inference for conditional Gaussian nonlinear systems. The
scalar and vector cores, the online smoother identities, and the CIR
reduction are in genuine agreement with the MATLAB reference wherever both
sides compute the same quantity. The local test suite completed with zero
failures. The remaining disagreements with the authors' scripts are almost
all *designed*, and they are designed for the right reasons.

One disagreement was not designed, only mis-diagnosed. The unexplained
4.58e-09 on the approximate objective CIR is the odd-interval Simpson
closure. It is not a defect in the CIR row, the online smoother, or the
relative-entropy formula. The earlier elimination of that hypothesis tested
the wrong integrand.

The package is not yet a one-to-one replica of every published figure, and it
should not pretend to be. Reproducing those figures needs three knobs used
together: `horizon` set to the reference's `last_idx`, `margin` stood down,
and Garcia's even-length Simpson closure in place of aciR's 3/8 panel. The
first two already exist. The third is a ten-line change and is the single
highest-leverage remaining step if figure-level parity is a goal.

**Recommendation.** Treat the current tree as a research preview that is
numerically trustworthy on the graded surface. Close the documentation drift
(several surfaces still say the online smoother and CIR are scalar-only).
Decide, as a named choice, whether `objective` is allowed to differ from
`simps.m` at the 1e-9 to 1e-7 level, or whether figure-level transcription
wins. Do not reverse the saturation policy; change what it *returns*.

| Dimension | Grade | One line |
|---|---|---|
| Fidelity of graded filter / smoother / metric | A+ | Machine precision on multiple fixtures, including correlated noise |
| CIR arithmetic vs the authors' code | A | Peak and subjective at 1e-15 once the horizon matches; objective residual now identified |
| Saturation and horizon policy | A− | Scientifically the better reading; the default then disagrees with every published CIR panel |
| Vector online smoother | A− | Equations (3.5)–(3.7) in full generality; the truncation story overclaims |
| Numerical strictness (Chol vs `det` / `pinv`) | A | Right default for a package; needs an explicit escape, not a silent `pinv` |
| Causal hygiene / estimand | A | Model-conditional KL, stated as such; the name still invites misreading |
| Package engineering | A | Lean, contracted, oracle-gated; `aci_cir` is the incomplete reporting surface |
| Documentation freshness | B | Several surfaces still deny capabilities that shipped |
| Performance vs the "truncation bounds the work" claim | B | `tol` does not fire at paper scale; full-lag online smoother is quadratic |
| Completeness vs the MATLAB scripts | A− | Conditional ACI and predator-prey CIR are implemented and ungraded against the authors |

---

## 2. What was inspected and how

Read: every file under `aciR/R/`, the CIR / horizon / quadrature / online
smoother / oracle / grading-matrix tests, `DESCRIPTION`, `NAMESPACE`,
`NEWS.md`, `API_STABILITY.md`, the three vignettes, the saturation and
online-smoother cairns, `oracle/parity/reports/PARITY_REPORT.md`,
`matlab_reference/simps.m`, the CIR block of `dyad_interaction_model.m`
(lines 433–561), and `oracle/aci_oracle_cir.m`.

Ran, locally:

* `devtools::test()` on the installed-from-source tree. Zero failures. The
  brief's count of 652 assertions is consistent with the file-by-file
  output; this review does not re-count them.
* A dedicated probe
  (`design/artefacts/2026-08-14_reviewer_probe.R`) that transcribes
  `simps.m` into R, rebuilds the `arbitrary_cross_noise` CIR at the
  reference horizon, and measures the five brief targets. Wall time of the
  CIR sweep over 1,051 times: 0.13 s.

Did not run: MATLAB, the parity harness, or `R CMD check --as-cran`. The
brief already records that check as 0 errors, 0 warnings, 1 NOTE (private
URL 404s). Those were not re-derived.

---

## 3. The five targets

### (a) The unexplained 4.58e-09 — identified

**Finding.** The residual is the even-length closure of composite Simpson.
aciR and `simps.m` implement different rules on an odd number of intervals.
They agree to machine precision on an even number of intervals. The CIR
*row* is not the source.

**Evidence.** On `arbitrary_cross_noise` (`dt = 0.002`, `first_idx = 251`,
`last_idx = 1301`, reported region `index <= 1001`, 751 times):

| Comparison | max abs, reported region |
|---|---|
| `peak` vs MATLAB CSV | 9.28e-15 |
| `objective` via `.aci_simpson` vs MATLAB CSV | **6.97e-07** |
| `objective` via a transcription of `simps.m` applied to *aciR's* RE row vs MATLAB CSV | **1.29e-14** |
| `.aci_simpson(re)` vs `simps.m(re)` on odd-length rows (even intervals) | 6.7e-15 |
| `.aci_simpson(re)` vs `simps.m(re)` on even-length rows (odd intervals) | 6.97e-07 |

The previously reported 4.58e-09 is not a typical scale, it is a specific
row. At `j = 408` (`n_re = 894`, even) the residual is 4.580e-09. The
first quartile of even-length residuals on this dataset is 4.60e-09. Odd
length residuals on the same dataset never exceed 6.7e-15.

The two RE formulae (`delta - log1p(delta)` versus
`cov_ratio - 1 - log(cov_ratio)`) agree to 0 on the worst row. The
dispersion form is not the residual.

**Why the earlier elimination failed.** The brief records that the closures
were compared on a test function (they do differ, by up to ~4e-3) and then
dismissed because "the CIR integrand decays to ~0 at its endpoint, so the
closure term is multiplied by nothing." That argument holds only for a
comparison that runs to the end of the record, where the last entry is the
relative entropy of a posterior against itself and is exactly zero. At the
*reference* horizon the last entry is the relative entropy of a still-partial
posterior against the fully informed one, and it is not small.

On this dataset, in the reported region, the last-exit time of the
`threshold = 1e-5` level set sits at a median fraction 1.00 of the truncated
row (minimum 0.91). Every one of the 751 reported RE sequences is
non-monotone. The integrand has not decayed. Garcia's even-length adjustment
integrates a quadratic through the last *three* samples over the last
interval; `y_n = 0` is not the situation, and even if it were, `y_{n-1}` and
`y_{n-2}` still contribute.

The two closures, written out:

* **MATLAB `simps.m` (Garcia 2009).** Composite Simpson 1/3 on the longest
  odd-length prefix. If `length(y)` is even, a quadratic interpolant through
  the last three abscissae is integrated over the last interval by a
  Vandermonde solve.
* **aciR `.aci_simpson`.** Composite 1/3 on an even interval count. If the
  interval count is odd, 1/3 up to `n - 3` and a Simpson 3/8 panel on the
  last three intervals. The 3/8 panel further assumes equal spacing
  (`aci-quadrature.R`, `.aci_simpson_38`).

On equal-spaced unit abscissae (the approximate objective, after the caller
scales by `dt`) both are third-order and both are legitimate. They are not
the same number. On a decaying *full-record* CIR row the difference is
negligible against the rest of the package's grades. On a *truncated-horizon*
row it is 1e-9 to 1e-7, which is exactly the open residual.

**What to do.** See suggestion S1. If the package's claim is one-to-one
reproduction of the authors' *numbers*, adopt Garcia's closure for
`objective`. If the claim is a better quadrature of the same integrand, keep
the 3/8 panel, document the 7e-7 as a designed difference, and stop calling
the residual unexplained. Either is defensible. Mixing them is not.

A by-product, confirming the brief: on `exp(sin x)` over `[0, 2]`, aciR's
3/8 closure is more accurate than `simps.m` at several even lengths (n = 12:
6.3e-5 vs 1.4e-4; n = 130: 1.3e-9 vs 9.8e-9) and identical at odd lengths.
Accuracy and transcription are different goals.

---

### (b) The saturation policy — keep the idea, change the return

**Finding.** `margin` is the right *kind* of instrument and 0.1 is a
reasonable default for a reusable package. A fraction is a better shape than
the reference's absolute `lookahead_tolerance = 0.6`, which is a dyad-plot
constant dressed as a method. The policy then does what the cairn
(`2026-08-13-cir-saturation-margin`) says it will: at the reference's
horizon, aciR marks every reported time saturated.

Measured, same 751 reported times, `horizon = last_idx = 1301`:

* `margin = 0.1` (default): **751 / 751 saturated**
* `margin = 0.001`: 587 / 751 saturated
* last-exit of `1e-5` as a fraction of the truncated row: min 0.91, median 1

So the brief's observation is not a corner case. Inside the reference's own
comparison window the divergence has not fallen below the reporting
threshold, because the window was chosen for a figure, not for resolution.
The reference still prints a number. aciR prints `NA`. aciR is right that
the number is not a resolved range. It is wrong to throw the number away.

**Attack on the default.** A fraction of *whatever record remains after this
time* is not a property of the system. On a 2,000-step record, 0.1 is 200
steps. On a 200-step window it is 20 steps. The amount of lookahead a range
needs is the range itself, which is circular if encoded as a static
fraction. The reference's 0.6 time units has the opposite problem: it is a
physical scale that happened to suit one figure.

**Attack on returning `NA`.** A time whose last-exit consumes 91–100 % of
the available record is not "no measurement". It is a right-censored
measurement: the range is *at least* as long as the unused-margin cut. Survival
analysis would report a lower bound and a censoring flag. Returning `NA` for
both `objective` and `objective_exact` discards the only information the
record actually contains, and it does so at the end of the record, which is
where a user studying a recent event looks.

Resolution is already judged per quantity, which is correct (a demanding
subjective threshold outruns the objective range, and condemning the
objective range for that would discard the headline number). The objective
range is judged against decay of the divergence to `threshold`, not against
the most demanding epsilon. Keep that.

**What to do.** See S2. Keep `margin`. Keep the flag. Return the truncated
value as a lower bound rather than as `NA`, and make the bound visible
(`objective_censored`, or a `status` column with values `resolved` /
`censored` / `below_threshold`). `margin = 1e-9` remains the switch that
makes the oracle compare like with like.

Is 0.1 the right default? It is a judgement, as the cairn already says. It
should stay documented as this package's device. I would not move it without
a calibration study across the dyad, the predator-prey pair, and ENSO. I
would not replace it with an absolute lookahead.

---

### (c) `.aci_cir_row()` and the horizon — the asymmetry is right

**Finding.** Truncating the comparison sequence and *not* truncating the
fully informed posterior is the quantity the method asked for, and it is
also what the reference computes. The prefix property is the right test.

The argument, stated so it can be attacked:

> The range at time *t* is how far forward one must look before the
> online-smoother posterior at *t* stops changing. "Stops changing" is
> only meaningful against a fixed target. The target is the estimate
> informed by everything available. Cutting the target at the same
> horizon as the comparison asks a different question: how far one must
> look to match a still-partial estimate. That question's rows are not
> prefixes of each other, and a time that would have been resolved
> against the full record can look resolved against a nearby partial
> one.

Measured on `arbitrary_cross_noise`, time `j = 351`, horizons
`{401, 551, 1301, 2001}`: the first entry of the row — relative entropy of
the filter against the target — is `1.990093e-02` at every horizon, to all
reported digits. The truncated row is a byte-identical prefix of the
untruncated one. If the target had moved, the first entry would have moved.

The reference does the same thing, somewhat by accident. In
`dyad_interaction_model.m` the comparison loops `obs = j:last_idx`, but the
fully informed side is `online_fixed_mean{end}(j)`, and `{end}` is the last
cell of a triangle of length `N+1`, the full record. aciR's `horizon`
reproduces that convention on purpose. The test in `test-cir-horizon.R` that
asserts the prefix property is the one that would fail if this were
implemented the tempting way.

**Attacks that do not land.**

* "The user who passed `horizon` asked for a shorter study." They asked for
  a shorter *comparison*, or they asked to reproduce a figure. In the first
  case the scientific target is still the best estimate the record supports.
  In the second case the reference itself uses the full-record target. Either
  reading keeps the asymmetry.
* "Near the end of a long record the extra tail has no influence." True
  once the ordered products have decayed, which is the saturation region.
  The asymmetry is load-bearing exactly where saturation is being decided.
  That is a reason to keep it, not to drop it.

**Attacks that do land, as documentation rather than as a code change.**

* The default `horizon = NULL` (whole record) will not reproduce the
  published CIR figures. `NEWS.md` says this. The Rd page says this. A user
  who copies the dyad vignette and stares at Figure 1 of the paper will
  still not understand why the whiskers differ. The Rd should show the
  three-argument invocation that recovers the figure (`horizon`, `margin`,
  and — once S1 lands — the quadrature).
* `horizon` is counted from the start of the record, not from the reported
  time. That matches the reference's `last_idx`. It is easy to misread as a
  per-time lookahead. The argument name is slightly too short for its
  meaning. `compare_until` or `last_index` would have been clearer. I would
  not rename it now; I would add one sentence to the Rd and an example that
  sets it from a time.

The implementation comment at `aci-cir.R:323-329` is the right length and
the right claim. Do not weaken it.

---

### (d) The vector online smoother — the truncation tolerance is not doing what the comment claims

**Finding.** The recursion is the paper's (3.5)–(3.7) and (3.12)–(3.16) in
full generality, including a non-zero noise cross-covariance that the
reference specialises away. At one dimension it collapses onto the graded
scalar path at every lag, not merely at full lag. The independent MATLAB
transcription of the general form is graded to 1e-6 on a fixture that
exercises lags past 1,000. I did not re-derive the algebra.

The truncation story is overstated.

The comment and the Rd say the products are truncated once they fall below
`tol` relative to the leading term, "which bounds the work at
`O(n * lag_effective)`". The test is `max(abs(d)) < tol` with default
`tol = 1e-18`. That is a max-entry test on the accumulated product, not a
spectral-radius test on the factor, and at paper scale it does not fire.

Measured, two-dimensional fixture, `n = 401`, `dt = 0.01`:

* Spectral radius of each `E_j`: min 0.9920, median 0.9920, max 0.9951.
  None reach 1. The paper's bound holds on this system.
* `lag_effective` at `tol` in `{1e-8, 1e-12, 1e-18, 1e-30}`: **400** every
  time, i.e. the whole record.
* Product norms from `j = 30`: `max(abs)`, Frobenius, and spectral radius
  of the running `D` never fall below 1e-18 in the remaining 371 steps.
* CIR-row values at `tol = 1e-18`, `1e-8`, and `0` are identical.

The same thing happens on the scalar dyad, where the products *can* be
formed in logarithms and so the cost of not truncating is only the number
of offsets visited:

| `n` | `lag_effective` | median `|E|` | time |
|---|---|---|---|
| 500 | 499 | 0.9929 | 0.018 s |
| 2,001 | 2,000 | 0.9956 | 0.097 s |
| 5,001 | 5,000 | 0.9984 | 0.491 s |
| 10,001 | 10,000 | 0.9984 | 1.872 s |

Times scale as `n^2`. `tol` never cuts the loop. At `|E| ≈ 0.9984`, reaching
`1e-18` takes about `18 * log(10) / -log(0.9984) ≈ 26,000` steps. Published
dyad records are 30,001 steps. The tolerance is a safety catch for a
product that has already underflowed, not a bound on the work.

On the vector path the same arithmetic is a loop of matrix products. At
`n = 401`, dimension 2, full lag: 0.255 s. Scaling as `n^2 * n_y^3` puts
the ENSO case (`n ~ 3e4`, `n_y = 3`) at tens of minutes in pure R. That is
consistent with the cairn that moved published-scale reproduction out of
the package. It is not consistent with the comment that truncation has
already solved this.

`max(abs(d))` is a sufficient condition that the next innovation cannot
move the estimate at double precision, and it is cheap. It is not the
spectral-radius argument the comment cites. A matrix can have every entry
above `1e-18` and still be a contraction; the test then continues, which is
conservative (extra work, not a wrong number). The opposite failure —
stopping while a non-normal product is still in a transient-growth regime —
is possible in principle and not visible on the fixtures I ran. I would not
change the test to a spectral radius of `E_j`. I would change the comment,
and I would set the lag from a measured contraction rate (S6).

The CIR-row truncation in `.aci_cir_row_mv` fills the rest of the sequence
with the current estimate and then takes `mu_end` from that fill. If
truncation is valid, the fully informed posterior is unchanged. If it is
not, two errors land together: the target moves, and the tail of the
divergence is written as exact zero. On the fixtures I ran, `tol` does not
fire, so this path is not exercised. A test that *forces* an early
truncation and checks `mu_end` against a `tol = 0` walk would lock the
claim the comment is making.

---

### (e) Numerical routes that differ by design — keep the stricter one

**Finding.** Using Cholesky for the multivariate relative entropy and for
the observation-noise inverse is the right default. The measured
disagreements with the reference (6e-14, 1e-13) are the cost of a better
conditioned evaluation, not a defect. Requiring positive-definiteness at
every slice, and refusing systems the reference will analyse with `pinv`,
is a package-versus-script decision and should stay on the package side.

Reasons, in order:

1. A relative entropy is undefined for a covariance that is not symmetric
   positive definite. Returning a number from `det` / `pinv` on that slice
   is a silent change of estimand. The reference is a set of analysis
   scripts whose author can see the `pinv` and know what it did. A package
   is used by people who will not.
2. The Cholesky evaluation of the multivariate KL (quadratic form by
   triangular solve, trace as a sum of squares, log-determinant as a
   difference of logged diagonals) is the form that stays accurate where
   the two posteriors nearly agree — which is the tail that sets the CIR.
   Forming `R_f^{-1}` and `det(R_f)` loses exactly that tail.
3. The observation-noise inverse is computed once, in the validator, and
   may also be *supplied* (`S_xoS_x_inv`). That is how conditional ACI is
   expressed. A `pinv` default would collapse the distinction between "this
   covariance is singular and I want the Moore-Penrose weight" and "I am
   deliberately zeroing a block". The second is a modelling choice and is
   already a first-class path.

The strictness will refuse some systems the reference analyses. I probed a
two-dimensional observation-noise covariance with condition number ~2e12.
`chol` itself succeeded; the *filter* then drove the filtered covariance
non-positive at step 2 and stopped with the named-step error. That is the
right failure: an explicit Euler step on an almost-singular observation
noise is not a well-posed discrete filter, and `pinv` would have hidden it.

**What to do.** See S3. Keep Cholesky as the only silent path. If a user
needs the reference's behaviour, give them `on_singular = "pinv"` (or
`"warn"`) as an explicit, documented escape that is off by default and that
records itself on the returned object. Do not add it in order to pick up
one more MATLAB fixture.

---

## 4. Additional findings

These are outside the brief's five targets and inside the audit.

### A1. Documentation still says the online smoother and CIR are scalar-only

`DESCRIPTION` (lines 31–32), `README.Rmd`, `vignettes/assumptions.Rmd`
(lines 100–105), and two `does_not_grade` notes in
`inst/extdata/oracle-manifest.yml` all state that the online smoother and
the causal influence range remain scalar-only. The code, the tests, and
`NEWS.md` say otherwise. `API_STABILITY.md` is the one surface that has
been updated.

This is not a wording nit. The assumptions vignette is the document a
careful user reads before trusting a vector result, and it currently tells
them the function they are about to call does not exist.

### A2. The default epsilon grid is not the reference's grid

`aci_cir()` documents `epsilon` as "the grid of the reference
implementation" and then defaults to
`10^seq(-6, 0.5, length.out = 129L)`. The reference uses 513 points
(`dyad_interaction_model.m:441`). The approximate `objective` does not use
the grid, so this does not touch the 4.58e-09. `objective_exact` *is* a
quadrature over that grid. 129 versus 513 is a coarsening, not a
transcription. Either change the default to 513, or document 129 as this
package's cheaper grid and show the 513-point call that matches
`defn_objective_CIR`.

### A3. The 3/8 panel is wrong on a logarithmic abscissa

`.aci_simpson_38` sets `h <- (x[i[4]] - x[i[1]]) / 3` and applies the
equal-spacing 3/8 weights. `objective_exact` integrates subjective ranges
over a logarithmically spaced epsilon grid. The default 129 and the
reference 513 both have an even interval count, so the 3/8 branch does not
run on those grids, which is why `objective_exact` already agrees with
`defn_objective_CIR` to 1.6e-14. An even-length `epsilon` (odd interval
count) fires the branch. On an 8-point log grid the equal-spacing 3/8 and
Garcia's unequal-spacing closure already differ by several units. This is a
latent defect, waiting for a caller to pass `length(epsilon) %% 2 == 0`.

### A4. `aci_cir()` has no reporting methods

`print`, `summary`, `plot`, and `as.data.frame` exist for `aci` and not for
`aci_cir`. The returned object is a list with a class and no behaviour. A
user who types `rng` at the prompt sees a raw dump of a 129-by-n subjective
matrix. Given that saturation, `objective` versus `objective_exact`, and
the peak are the quantities the method leads with, this is the incomplete
half of an otherwise careful reporting surface.

The high-level `aci()` is also still scalar-only. A vector system has to be
assembled by hand (`aci_filter` / `aci_smoother` / `aci_metric`). That is
defensible for an expert surface. It is inconsistent with shipping
`aci_enso_model()` as a user-facing constructor.

### A5. `aci_cir()` does not validate a vector `filt`

The scalar path runs `.aci_check_posterior`. The vector path checks that
`filt` is not `NULL` and then uses it. A malformed vector posterior fails
inside `.aci_cir_row_mv` with a low-level subsetting error rather than with
the contract message the rest of the package raises.

### A6. Partial matching in the vector online-smoother grade

`tests/testthat/test-online-smoother-mv.R` reconstructs the lagged estimate
with `aux$E[, , k]`. `.aci_online_aux_mv` returns `E_j`, not `E`. The test
passes because `$` partial-matches `E` to `E_j`. The scalar auxiliary
*does* return `E`. A later field named `E`, or
`options(warnPartialMatchDollar = TRUE)`, turns a passing grade into a
failure that looks like a numerical regression. Write `aux$E_j`.

### A7. Merged roxygen on two internal helpers

In `R/aci-validate.R` the title and description of `.aci_stop_covariance`
sit immediately above the title of `.aci_expand`, and `.aci_stop_covariance`
itself is left with only `@returns Never returns`. Both are `@noRd`, so
the manual is unaffected. The source is the contract for the next editor.

### A8. Known ungraded surfaces, restated so they are not "found"

As the brief already says, and as I confirmed by reading the tests rather
than by running MATLAB:

* Predator-prey filter / smoother / metric is graded to 1e-6 in both
  directions. Predator-prey CIR is not compared to the authors.
* `aci_conditional()` has tests for the construction (the inverse is
  supported on the target block; a zero-noise target is refused). It is not
  graded against the one ENSO script that has conditional ACI enabled.
* Four of five ENSO scripts are captured and not gated.

These are scope statements, not defects. They become defects if the README
or the validation vignette claims capability parity with the reference on
those surfaces.

### A9. `aci_simulate()` still refuses correlated noise

The filter, the smoother, the metric, the online smoother, and the CIR all
accept a non-zero `S_yoS_x`. The simulator does not. A user who wants a
path from the cross-noise system the rest of the package is proud to grade
has to write their own Euler step. The joint covariance is already checked
for positive-semidefiniteness at construction. Drawing a pair of increments
from its Cholesky factor is the natural widening and is additive.

---

## 5. Package-quality scorecard

A light pass against the estate's R-package invariants, not a `/rpkg audit`
re-run.

| Invariant | Status | Evidence |
|---|---|---|
| DESCRIPTION metadata | Pass | UTF-8, MIT, ORCID on the maintainer, URL + BugReports, `Language: en-AU`, R >= 4.1.0 empirically on CI |
| NAMESPACE | Pass | Explicit exports, no `Depends` beyond R |
| Lean Imports | Pass | `graphics`, `stats` only; ggplot2 / knitr in Suggests |
| Runnable examples | Pass | Every export I opened carries one; `aci_cir` example is short enough for check |
| Tests | Pass | testthat edition 3; local `devtools::test()` green; oracle tests never skip |
| CI matrix | Pass | `.github/workflows/R-CMD-check.yaml`: macOS, Windows, Ubuntu × devel / release / oldrel-1 / 4.1 |
| NEWS | Pass | Breaking / features / fixes shape; before/after on `horizon` |
| Long-form docs | Pass | Three vignettes + pkgdown config |
| API stability | Pass | `API_STABILITY.md`, lifecycle stages, nothing marked Stable at 0.1.0 (honest) |
| `R CMD check` | Pass, on the brief's word | 0 E / 0 W / 1 NOTE (private URLs). Not re-run here |
| Hand-crafted R | Pass, with A6–A7 | Dash banners, delegated validation, `call. = FALSE`, leading-dot internals, paragraphed roxygen |

Andreou is in `Authors@R` as `cph` of the MATLAB reference, without ORCID.
That is the correct role. An ORCID would be a courtesy, not a requirement.

The 16 July quality review (grade A−) described a package that did not yet
have CIR, vector states, or time-varying `L_y`. Those have since landed.
The present tree is a different package, and the completeness grades in
that review should not be cited as current.

---

## 6. Constructive suggestions

Ranked by (severity of the gap they close) × (how much work they are).
Each is a suggestion, not a patch. The first three are the ones I would
want a decision on before any more features.

### S1. Make the approximate-objective quadrature a named choice

**Why.** The 4.58e-09 is this choice, un-named. Figure-level parity with
the authors needs Garcia's closure. A slightly more accurate integral of
the same samples wants the 3/8 panel. The package currently does the
second and is graded against the first.

**What.** Add `quadrature = c("reference", "simpson38")` to `aci_cir()`,
defaulting to `"reference"` if one-to-one numbers are the charter, or to
`"simpson38"` if they are not, with the default written down in
`API_STABILITY.md`. `"reference"` is a transcription of `simps.m` (the
probe already contains one, in R, that recovers the MATLAB CSV to 1e-14).
`"simpson38"` is the current `.aci_simpson`. Report which one was used on
the returned object.

If a second argument feels heavy, do the smaller thing: switch `objective`
to Garcia's closure, keep `.aci_simpson` for `objective_exact` on equal
grids, and write one paragraph in the Rd that says the two rules differ at
the 1e-9 to 1e-7 level on truncated rows and that this package follows the
authors on the headline number.

**Do not** raise the oracle tolerance to 1e-6 and call the residual closed.
That is the move the mixed-tolerance cairn exists to forbid.

### S2. Report a censored lower bound, not `NA`

**Why.** Section 3(b). The flag is right. The missing number is not.

**What.** For a time that fails the margin test, keep `saturated = TRUE`
and return the computed `objective` / subjective values as *lower bounds*.
Surface this in `summary.aci_cir()` as, for example, `>= 0.18 (censored)`.
A user reproducing a figure stands the margin down, as now. A user asking
a scientific question about the end of a record is told the range is at
least this long, which is the only true statement the record supports.

A more ambitious version treats the last-exit time as a right-censored
observation and reports a Kaplan–Meier-style curve of subjective range
versus threshold, with the censoring times marked. That is a paper, not a
patch. The lower bound is the patch.

I would not change the default `margin`. I would add one calibration
figure, in the assumptions vignette or a short methods note, showing
resolved fraction against `margin` on the dyad, the predator-prey pair, and
ENSO. 0.1 may survive that figure. It should not be inherited unseen.

### S3. Keep Cholesky; add an explicit singular escape

**Why.** Section 3(e). The stricter route is the package route. The
reference's `pinv` is sometimes what a user with a rank-deficient
observation-noise model actually wants.

**What.** `on_singular = c("error", "pinv")` on `aci_filter()` / the
validator, default `"error"`. `"pinv"` uses a Moore-Penrose inverse, warns
once, and records `singular_method = "pinv"` on the result. Conditional
ACI continues to go through `S_xoS_x_inv` and does not use this switch.

### S4. Refresh every surface that still says "scalar-only"

`DESCRIPTION`, `README.Rmd` (then rebuild the README),
`vignettes/assumptions.Rmd` § "State dimension", and the two
`does_not_grade` notes in `oracle-manifest.yml`. Replace with: the online
smoother and CIR now dispatch on the components the same way the core
does; the vector online smoother is graded by collapse onto the scalar
path and by an independent MATLAB transcription of the general form; it is
not graded against the authors' ENSO scripts.

While there, fix A2 (129 is not 513) and A5 (validate a vector `filt`).

### S5. Give `aci_cir` a reporting surface, and `aci()` a vector path

`print.aci_cir`, `summary.aci_cir`, `as.data.frame.aci_cir`, and a two-panel
`plot.aci_cir` (objective range against time, with saturated times drawn
differently; subjective range as an image against `log10(epsilon)`). The
`aci` methods are the template. This is the difference between a result
object and a list.

`aci()` should accept a vector model, or a sibling `aci_fit()` should. The
enso constructor without a one-call workflow is a constructor for experts
who will write the three-call sequence anyway.

### S6. Set the online-smoother lag from a measured contraction, not from `1e-18`

**Why.** Section 3(d). `tol` is not a work bound at paper scale.

**What, cheap.** After `.aci_online_aux` / `_mv`, estimate
`rho_hat = quantile(spectral radius or |E_j|, 0.9)` and set
`lag_cap = min(lag, ceiling(log(tol) / log(rho_hat)))`, with a floor so a
noisy `rho_hat >= 1` does not explode. Visit at most `lag_cap` offsets.
Record `rho_hat` and `lag_cap` on the result so a user can see whether
truncation happened.

**What, vector path.** The products cannot be reconstructed from endpoints.
They can be accumulated with periodic restart: every `k` steps, drop a
factor whose running product is below `tol` in Frobenius norm and start a
new block. That is the same bound the comment already claims, implemented
as arithmetic rather than as a hope about `1e-18`.

**What, full-lag scalar path.** At `lag = Inf` the online smoother is the
discrete smoother. The package already knows this and tests the
convergence rate against `aci_smoother()`. For the *value* at full lag,
forming the CIR-style cumulative sum once is `O(n)` after the `O(n)`
auxiliary pass, and does not walk `n` offsets. The offset loop is the
right algorithm for an intermediate lag; it is the wrong algorithm for
the two boundaries, one of which is the default.

I would not introduce compiled code for this. The pure-R cairn is still
right for the package's advertised scale. The above are R-level shape
changes.

### S7. Parallelise the CIR time loop, not the inner row

Each reported time is independent once the auxiliary quantities exist.
`window` of a few thousand is the expensive call, and it is embarrassingly
parallel. `parallel::mclapply` over `window`, with `mc.cores` defaulting to
1 and documented as safe up to the machine, would cut the published-figure
batch by the number of cores without touching the arithmetic. The inner
row is already a `cumsum`. Do not parallelise that.

I would cap this at the 6 cores the brief budgeted, and I would not make
it the default. CRAN check has one core and a memory cap; a default of 1
keeps the check honest.

### S8. Fix the latent 3/8-on-log-grid defect

Three options, in decreasing preference:

1. Use Garcia's unequal-spacing closure everywhere `x` is supplied. Then
   `objective_exact` on any epsilon grid, even or odd, matches the
   reference's `simps(x, y)` and the 3/8 equal-spacing assumption goes
   away.
2. Refuse an even-length `epsilon` (odd interval count) when calling
   `objective_exact`, with a message that names the 3/8 limitation.
3. Implement an unequal-spacing 3/8. More work, no transcription benefit.

(1) plus S1 is one change, not two.

### S9. Grade the two surfaces the brief already listed as open

In the order the brief implied:

1. Predator-prey CIR against the extracts that already exist. This is the
   cheapest remaining authors'-code grade, and it is the first CIR grade
   on a system whose `L_y` moves.
2. `aci_conditional()` against the one ENSO script that has it enabled.
   The harness note in brief §6 matters here: that script is not the same
   question as its four siblings. The grade has to say so.
3. The four ungated ENSO extracts, each with its own input list. These
   are infrastructure, not package arithmetic.

I would not start (3) until (1) and (2) are green. The package's claim is
already stronger than most method-reimplementation packages; these three
close the last "we implemented it but we did not grade it against them"
gaps.

### S10. Simulate correlated noise

`aci_simulate()` draws two independent `rnorm` streams and refuses
`S_yoS_x != 0`. Widen it: factor the 2-by-2 (or `(n_x + n_y)`-by-same)
Grammian once, multiply, and keep the contained-seed contract. The
cross-noise oracle currently has to bring its own signal. That is the
right way to *grade*. It is the wrong way to *demo*.

### S11. Coding hygiene, small and local

* `aux$E_j` in `test-online-smoother-mv.R` (A6).
* Split the roxygen of `.aci_expand` and `.aci_stop_covariance` (A7).
* `options(warnPartialMatchDollar = TRUE)` in `tests/testthat/setup.R`,
  so A6 cannot recur.
* A `print.aci_cir` that does not dump the subjective matrix (S5).
* Drop the "grid of the reference implementation" clause, or make it true
  (A2).

None of these change a number.

---

## 7. Mathematics and causal method

The assumptions vignette is the best scientific page in the package. What
follows is additive to it, not a replacement.

### M1. Last-exit is not the layer-cake integral, and the rows are not monotone

The paper defines the objective range as the average of the subjective
ranges over thresholds, and offers `integral(RE) / peak` as a cheaper
underestimate. Those two coincide when RE is a decreasing function of lag,
by the layer-cake identity
`∫_0^{peak} λ({s : RE(s) > ε}) dε = ∫ RE`.

The implementation, and the reference, do not use the Lebesgue measure of
the superlevel set. They use the *last* time the divergence exceeds the
threshold (`find(..., 1, 'last')`; in aciR, a suffix-maximum plus
`findInterval`). Last-exit is at least as large as the measure of the
superlevel set, and is strictly larger the moment RE is non-monotone.

On the truncated-horizon rows I computed, **751 / 751 were non-monotone**.
That is the regime in which `objective` and `objective_exact` are
different functionals, not two quadratures of the same one. The package
already returns both, which is the right response. The Rd should say that
they part company when the comparison sequence is not decreasing, and that
a truncated horizon is exactly when that happens. A user who sees
`objective < objective_exact` and concludes that Simpson is biased is
misreading a definitional gap as a numerical one.

I would also add one diagnostic to `summary.aci_cir`: the fraction of
reported times at which the RE row is monotone. If that fraction is small,
the efficient underestimate is the wrong headline number.

### M2. ACI is not Pearl, and conditional ACI is not a back-door

The vignette already denies interventional claims. Two sharper sentences
would help the readers who will cite this package from a causal-inference
venue:

* The metric is `KL(p(y_t | x_{0:T}, M) || p(y_t | x_{0:t}, M))`. It is a
  property of a model-conditional posterior pair. It is not
  `p(y | do(x))`, not a Granger test, and not transfer entropy, even
  though all four can be large at the same instants on a CGNS.
* `aci_conditional()` changes the *observation* the filter assimilates,
  not the *system* that generated the path. Zeroing a block of
  `S_xoS_x_inv` is "stop updating on this channel", not "this channel was
  never there" and not "we have adjusted for this channel in the
  do-calculus sense". The non-target components still drive the drift.
  That is the whole point, and it is also the thing a reader coming from
  Pearl will get backwards.

I would put those two sentences in the `aci_conditional` Rd, not only in
the vignette.

### M3. There is still no model criticism

The metric of a wrong model is a smooth, plausible, wrong curve. The
package says this and then offers no residual. Three cheap diagnostics
would make the warning operational:

1. Filter innovations `x_{j} - x_{j-1} - (L_x μ_f + f_x) dt` against the
   predicted observation-noise scale. A CGNS that fits has standardised
   innovations close to white. One that does not will still produce an ACI
   peak.
2. The probability-integral transform of the unobserved component, when
   the unobserved component is available (simulations, the dyad vignette).
   PIT histograms are the natural "does the filter know its own
   uncertainty" check.
3. Sensitivity of the metric to a one-at-a-time perturbation of each
   named parameter, reported as a companion series. The package already
   refuses to pretend the metric has a standard error. A tornado plot is
   the honest substitute.

None of these need a new oracle. They need a vignette section and a
helper that returns the innovation series `aci()` already computed and
threw away.

### M4. Parameter uncertainty is not optional once the coefficients are fitted

The vignette is clear that a metric computed at a point estimate is
conditional on that point. The moment someone estimates `gamma` from the
same record they then score, the ACI peak is a function of the data twice
and its uncertainty is not the Euler discretisation error. I would not
build a bootstrap into 0.1.x. I would add a one-paragraph recipe:
refit or perturb the parameters, rerun `aci()`, look at the band. The
contained `seed` on `aci_simulate` already makes that recipe
reproducible.

### M5. The discrete / continuous smoother split is a feature, not a discrepancy

`aci_online_smoother(..., lag = Inf)` is the discrete smoother.
`aci_smoother()` is an Euler discretisation of the continuous backward
equation. They agree at first order in `dt`, and the package tests the
rate. That split is the right thing to have found. I would surface it as
a diagnostic: `summary.aci` already reports the terminal identity
residual of the *continuous* pair. Reporting
`max |online_full - smoother|` next to it would tell a user whether `dt`
is small enough for the CIR (which is built from the discrete object) to
be read as a property of the continuous system the paper writes down.

### M6. An adaptive epsilon grid is already in the authors' comments

`dyad_interaction_model.m:446-455` sketches a two-piece mesh, denser
above `10^{-2}`, and leaves it commented out. `objective_exact` is a
quadrature of a monotone function of `ε` (subjective range, decreasing).
A geometrically spaced grid wastes nodes where the range has already
collapsed. If S1/S8 put Garcia's closure under `objective_exact`, the
next improvement is the authors' own adaptive mesh, exposed as
`epsilon = "adaptive"`. I would not make it the default until it is
graded against the 513-point geometric grid they actually ran.

### M7. A short comparison to Granger / transfer entropy belongs in the assumptions vignette

Not because those methods are competitors on a CGNS — the closed-form
smoother *is* the point — but because every new reader will ask. One
figure on the dyad, ACI metric against a sliding Granger F and a
binned transfer entropy, with the perfect-model caveat in the caption,
would save the package a decade of "why not just use …" reviews. It
would also make the estimand section concrete: the three series peak at
related times and answer different questions.

---

## 8. Efficiency, collected

The package is fast enough for everything it advertises, and honest about
the published-figure batch being a batch. The remaining gains are shape,
not language.

| Call | What it actually costs | What to do |
|---|---|---|
| `aci()` on the dyad, n = 2,001 | 0.011 s | Nothing |
| `aci_cir` on 751 times × horizon 1,301 | 0.082 s | Nothing for interactive use; S7 for figure batches |
| `aci_online_smoother`, scalar, full lag, n = 10,001 | 1.87 s, quadratic | S6: `O(n)` cumulative sum at `lag = Inf`; contraction-based `lag_cap` otherwise |
| `aci_online_smoother`, vector, n = 401, dim 2, full lag | 0.255 s, quadratic in n | S6; leave ENSO-scale in the companion repository |
| Vector CIR row | One pass of matrix products per time | Already the right shape; the `tol` fill-forward is unexercised (3(d)) |

I would not add Rcpp to this package to make the online smoother faster.
The cairn that reversed that decision is still right: the advertised
scale fits in R, and the published scale does not belong in a CRAN
example. I *would* stop walking `n` offsets to produce a quantity that is
a prefix-sum.

Memory is already the win. Not forming the reference triangle is what
makes CIR shippable. That design is correct and should not be revisited.

---

## 9. What I would not do

* Re-derive the CGNS algebra. The oracles have it.
* Replace `margin` with the reference's 0.6. That constant is a figure
  margin, not a method.
* Apply the horizon to the fully informed posterior. The prefix property
  is the quantity.
* Switch the multivariate KL to `det` and a right-division to pick up
  6e-14 of MATLAB agreement. That is the wrong direction.
* Raise any oracle tolerance to hide S1. The mixed-tolerance cairn is
  the right doctrine.
* Compile the inner loops in this package. Measure first, and measure
  after S6.
* Contact the authors about brief §6. That call is the maintainer's.

---

## 10. Suggested decision list for the maintainer

Three choices, then the mechanical work.

1. **Is `objective` allowed to differ from `simps.m` at 1e-9 to 1e-7 on
   truncated rows?** If no, S1 defaults to `"reference"`. If yes, document
   the residual as designed and close the open item in the parity report.
2. **Does a saturated time return `NA` or a censored lower bound?** I
   recommend the bound (S2). The flag stays.
3. **Is figure-level transcription of the published CIR panels a 0.1.x
   goal?** If yes, the README needs a four-line recipe: `horizon = last_idx`,
   `margin` stood down, 513-point epsilon, Garcia quadrature. If no, the
   default whole-record horizon stays, and the docs should say more
   loudly that the figures will not match.

Then S4 (docs), S5 (reporting), S11 (hygiene), S8 (3/8), S9 (the two
ungraded surfaces), S10 (correlated simulation), S6 (contraction-based
lag), S7 (optional parallel CIR), in that order.

---

## Appendix. Probe commands

```r
Rscript --vanilla design/artefacts/2026-08-14_reviewer_probe.R
```

Writes `design/artefacts/2026-08-14_cir_residual_rows.csv` (751 reported
times: peaks, both quadratures, residuals, last-exit fractions,
monotonicity). The 4.58e-09 row is `index == 408`.

```r
cd aciR && Rscript --vanilla -e 'devtools::test(reporter = "summary")'
```

Local suite, zero failures, 2026-08-14.

MATLAB was not used. The Garcia closure was transcribed from
`matlab_reference/simps.m` and checked by recovering the committed MATLAB
CIR CSV to 1.29e-14 when applied to aciR's own RE rows.
