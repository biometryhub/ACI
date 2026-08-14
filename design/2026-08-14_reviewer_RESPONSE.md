# aciR — response to external review

**Package** 0.1.0.9000 · **Date** 2026-08-14 · **Responding to**
`design/2026-08-14_reviewer_OUTPUT.md`

Every numbered item in the review is dispositioned below: **done**, **declined**
with a reason, or **deferred** with a target. Nothing is left unaddressed and
nothing is claimed that was not measured.

Gates after the work: **718 tests, 0 failures, 0 errors, 0 skips**;
`R CMD check --as-cran` **0 errors, 0 warnings, 1 NOTE** (new submission plus
three 404s from the declared URLs, which resolve when the repository is public).

---

## 1. The correction, accepted

The review's central finding is that my elimination of the odd-interval Simpson
closure was wrong. It was, and the way it was wrong is worth stating plainly
because it is a method error rather than an arithmetic one.

I had compared `.aci_simpson()` against `simps.m` on a test function, found
they differ by up to 4.3e-03, and then **argued** — rather than measured — that
this could not reach the CIR because its integrand decays to ~0 at the endpoint,
so the closure term is multiplied by nothing.

That argument is true of a comparison run to the end of the record, where the
last entry is the divergence of a posterior from itself and is exactly zero. It
is false of a truncated-horizon row, which is the only case in question. Now
measured, at `j = 408` on `arbitrary_cross_noise`:

| Row | length | last value | last / peak |
|---|---|---|---|
| Full record | 1594 | **0.0000e+00** | 0 |
| Truncated at `last_idx = 1301` | 894 | **9.2541e-06** | **2.71e-04** |

The failure was substituting an argument for a measurement at the point where
the argument was load-bearing. That is the same shape as a fixture that
annihilates the term it claims to grade, which is the defect class this package
was built around, and it is the third instance in the project. It is recorded
in `ledger/r-package/2026-08-14-reference-quadrature-closure.cairn.md` under
Forward cost, where the generalisable lesson lives rather than in a changelog.

Thank you for testing it properly.

---

## 2. Maintainer's rulings on the three decisions

From review §10, answered before any work began:

1. **May `objective` differ from `simps.m` at 1e-9 to 1e-7?** **No.** S1
   defaults to the reference rule.
2. **`NA` or a censored lower bound?** **Bound**, flag retained.
3. **Is figure-level transcription a 0.1.x goal?** **Yes.** The README carries
   the recipe.

---

## 3. The five targets

### (a) The 4.58e-09 — closed

`.aci_simpson()` now closes an odd interval count with the quadratic through
the last three samples integrated over the final interval, which is the
reference's rule.

**Derived, not transcribed.** Integrating the Lagrange basis through abscissae
at `t = -h1, 0, h2` over `[0, h2]` gives

```
-h2^3 / (6 h1 (h1+h2)) * y0
  + h2 (h2 + 3 h1) / (6 h1) * y1
  + h2 (2 h2 + 3 h1) / (6 (h1+h2)) * y2
```

with a test asserting it reduces to `h/12 * (-y0 + 8 y1 + 5 y2)` at `h1 == h2`,
and further tests asserting exactness for quadratics on unequal grids including
a logarithmic one. The general form was written rather than the equal-spacing
special case precisely so that S8 closes with it.

| | before | after |
|---|---|---|
| `objective` vs the reference, reference conventions | 4.58e-09 | **3.57e-15** |

The 3/8 panel this replaces is the more accurate rule on an equally spaced
grid — your `exp(sin x)` measurements match mine — and that accuracy is given
up deliberately under ruling 1. The trade is recorded as a decision, not as an
implementation detail.

The quadrature tests changed rather than being deleted. The old assertion was
"Simpson is exact for polynomials up to cubic", which the new rule does not
satisfy on an odd interval count. It is now split: cubic exactness where the
interval count is even, quadratic exactness on any spacing where it is odd, and
the equal-spacing reduction above. Asserting cubic exactness would have been
asserting a rule the package no longer uses.

### (b) Saturation — done as recommended

`margin` kept, default kept, flag kept, return changed.

`aci_cir()` now returns a four-valued `status`: `"resolved"`, `"censored"`,
`"below_threshold"`, `"insufficient"`. A censored time returns the computed
value as a lower bound. `saturated` is retained and equals
`status == "censored"`, so existing callers keep working. `subjective_censored`
is a logical matrix marking the individual thresholds that ran long, because
resolution is judged per quantity.

Only `"insufficient"` — fewer than three later observations, where no
quadrature exists at all — still returns `NA`. That is the one status with
genuinely nothing behind it, and separating it from `"censored"` is what stops
the flag claiming a bound exists when it does not.

Measured on `arbitrary_cross_noise` at the reference horizon: **1049 of 1051
reported times now carry a finite bound where all 1049 were holes.**

`summary.aci_cir()` reports the censored fraction. Not done: the calibration
figure for `margin` across the dyad, predator-prey and ENSO. It is the review
trigger on the cairn rather than a claim that 0.1 is right.

### (c) The horizon asymmetry — endorsed, docs improved

No code change; your reading matches the implementation's intent and the
`test-cir-horizon.R` prefix test is the one that would fail if it were done the
tempting way.

Done from your documentation attacks: the README now carries the
three-argument invocation that recovers the published panels. **Not done:** an
Rd example that sets `horizon` from a time rather than an index. The Rd already
says it is counted from the start of the record; the worked example is not
there.

### (d) The truncation tolerance — comment fixed, algorithm deferred

You are right that the comment claimed a bound the test does not deliver. The
Rd now says so explicitly: `tol` is a safety catch on a product that has
already underflowed, not a work bound; at a per-step contraction of 0.998 it
needs some 26,000 steps to reach `1e-18` and published records are shorter, so
the loop is not cut and the cost is quadratic. It also now states that
`max(abs(d)) < tol` is a sufficient condition rather than the spectral-radius
argument, and that erring this way costs work and never accuracy.

**Deferred: the algorithm (S6).** One hazard your note does not address. The
`O(n)` prefix-sum at `lag = Inf` requires factorising
`D[j,k] = exp(cum_log[k] - cum_log[j])` into `exp(-cum_log[j]) * exp(cum_log[k])`.
With `|E| ~ 0.998` over a 30,000-step record `cum_log` reaches about `-60`, so
those factors are `e^{+/-60}`. The current offset loop keeps the *difference*
bounded and the factorised form does not. The reformulation is right
mathematically and needs blocked rescaling to be numerically safe; that is not
a ten-line change and it should be measured, not assumed.

**Not done:** the test that forces an early truncation and checks `mu_end`
against a `tol = 0` walk. Accepted as correct — the fill-forward path in
`.aci_cir_row_mv` is unexercised on every fixture we have, which is exactly why
it deserves a test. It is the highest-value item on the deferred list.

### (e) Cholesky over `pinv` — endorsed; the escape declined

Kept as the only silent path, for your three reasons.

**S3 declined**, on evidence from your own probe. You reported that at
condition number ~2e12 `chol` *succeeded* and the **filter** then drove the
covariance non-positive at step 2 and stopped with the named-step error. So
`pinv` would not have rescued the motivating case; it would have moved the
failure somewhere less legible. An escape that does not help the case that
motivates it, while re-opening the silent-degradation door the package closed
deliberately, is surface for a user who has not appeared. Reconsider on a real
report.

---

## 4. Additional findings A1–A9

| # | Disposition |
|---|---|
| A1 | **Done.** `DESCRIPTION`, `README.Rmd` (rebuilt), `vignettes/assumptions.Rmd`, both `oracle-manifest.yml` notes. `DESCRIPTION` was CRAN-visible and denied a capability shipped in `e4cdc41` — the most serious item in the review by that measure |
| A2 | **Done, your second option.** 129 kept as this package's cheaper default; the false "grid of the reference implementation" clause removed; the 513-point call that reproduces `defn_objective_CIR` shown in the Rd and in the README recipe |
| A3 | **Done**, as part of (a). The replacement closure is exact for unequal spacing, so the logarithmic-grid defect cannot fire |
| A4 | **Half done.** `print`, `summary`, `as.data.frame` and `plot` for `aci_cir`; censored times are drawn hollow in the plot so a bound is not read as a measurement. **`aci()` still has no vector path** |
| A5 | **Done.** A vector `filt` is validated for shape and length with the contract message the rest of the package raises |
| A6 | **Done**, and made un-reintroducible: `aux$E_j` in the test, plus `options(warnPartialMatchDollar = TRUE)` in a new `tests/testthat/setup.R`. Warnings are errors under testthat, so the class cannot recur silently |
| A7 | **Done.** The `.aci_stop_covariance` header is back above its own function |
| A8 | Scope statements; nothing to do. They remain true and are restated in the parity report |
| A9 | **Not done** — see S10 |

---

## 5. Suggestions S1–S11

| # | Disposition |
|---|---|
| S1 | **Done**, defaulting to the reference rule per ruling 1. The second argument was *not* added: two rules mean two sets of numbers to grade and a user decision with one right answer once ruling 1 is made |
| S2 | **Done** as recommended |
| S3 | **Declined** — §3(e) |
| S4 | **Done** |
| S5 | **Half** — methods yes, `aci()` vector path no |
| S6 | **Comment done, algorithm deferred** with the overflow hazard above |
| S7 | **Declined.** 0.082 s for 751 times, and `mclapply` is Unix-only. Complexity and a platform split for no advertised benefit; your own note says "nothing for interactive use" |
| S8 | **Done**, your option (1), fused with S1 as you suggested |
| S9 | **Not done.** All three remain open |
| S10 | **Not done** |
| S11 | **Done**, all five items |

---

## 6. Mathematics and causal method

| # | Disposition |
|---|---|
| M1 | **Done, including the diagnostic.** The Rd now states that `objective` and `objective_exact` are different functionals wherever the sequence is not monotone, that the range is a last exit rather than a superlevel measure, and that a truncated horizon is when they part. `aci_cir()` returns a logical `monotone` per time and `summary.aci_cir()` reports the fraction. On the dyad at `n = 2500` it reports **0% monotone**, with `objective` 0.056 against `objective_exact` 0.097 — your finding corroborated on a second system |
| M2 | **Done.** An `@section Estimand:` on `aci_conditional`, in the Rd rather than only the vignette, stating the model-conditional KL and that masking a block means "stop updating on this channel", not an adjustment in the do-calculus sense |
| M3 | **Deferred to 0.2.0.** Agreed on merit. New scientific scope, not defect repair, and it wants a vignette section and a helper rather than a patch |
| M4 | **Deferred to 0.2.0**, same reason |
| M5 | **Not done.** Cheap and worth doing; it did not make this pass |
| M6 | **Deferred.** Correctly blocked behind S1/S8, which have now landed, so it is unblocked for 0.2.0 |
| M7 | **Deferred to 0.2.0.** Agreed that every reader will ask |

---

## 7. What changed, measured

| Comparison | Before | After |
|---|---|---|
| `cir_objective`, reference conventions | 4.58e-09 | **1.37e-14** |
| `cir_subjective`, reference conventions | 2.22e-16 | 2.22e-16 |
| `cir_peak`, reference conventions | 9.28e-15 | 9.28e-15 |
| Scalar core, two datasets | 26 / 26 | 26 / 26 |
| Finite objective ranges at the reference horizon | 0 of 1051 | **1049 of 1051** |
| Parity rows overall | 32, six reported as failures | **34, all passing** |

One structural change in the harness deserves a note, because it means the
package's CIR was never graded before this pass. Until `horizon` existed the
harness compared aciR's whole-record answer against the reference's truncated
one and reported the gap as a disagreement. Those are answers to different
questions. The harness now grades on the reference's own conventions and
reports aciR's default separately **with its magnitude**, as a designed
difference rather than a failure. The design choice is on the record; it is not
hidden by the change.

---

## 8. Two questions back

1. **`margin` calibration.** We did not run the study. Your text implies 0.1
   may survive it. Is your expectation that the resolved fraction is flat in
   `margin` over a useful range, or that the right default differs by system?
   The second would argue for making it a per-model default rather than a
   package constant.
2. **The forced-truncation test for `.aci_cir_row_mv`.** Agreed it should
   exist. Since `tol` never fires on any fixture we own, the test has to
   manufacture a contraction fast enough to trigger it. Would you accept a
   synthetic system with `|E| ~ 0.5` for that purpose, or does the test only
   mean something on a system the package would actually be used on?

---

## 9. Handling note, unchanged

The two findings about the reference implementation's own behaviour — the
predator-prey execution hazard, and that exactly one of the five ENSO scripts
ships with conditional ACI enabled — remain **internal**. They have not been
raised with the authors. That decision is the maintainer's and has not been
taken.


---

# Addendum: response to the second-round review

**Date** 2026-08-14 · **Responding to** `design/2026-08-14_reviewer_ROUND2.md`

All seven checkable claims were verified independently before acting. All hold.
Everything in F1 to F7 and both vignette passes is done. Version bumped to
**0.2.0**; the `aci_cir()` return and quadrature changes are breaking, so this
is a minor bump rather than a patch.

## The two corrections, accepted

**A1 was not done.** `aciR.Rmd` and `validation.Rmd` still denied capabilities.
Fixed.

**The grade is 1.37e-14, not 3.57e-15.** Fixed in `NEWS.md`, the cairn, and
this document. Worth recording *why*, so it does not recur: these were two
different measurements, not one number mistyped. Mine used `margin = 0.001`,
which censors more times and so takes a maximum over a **subset**; the harness
uses `margin = 1e-9` over the full reported region. The quoted figure is now
the maximum **and the set it is over** -- 1.37e-14 across 751 reported times on
the cross-noise record, 1.24e-14 across 1001 on the dyad.

## The pattern behind both, fixed as a mechanism

Round 1 caught a `summary.aci_cir()` docstring claiming a diagnostic the code
did not compute. Round 2 caught A1 marked done from four surfaces without a
search over the rest. Same error: **asserting completeness from the surfaces
edited rather than from a search over the surfaces that could carry the claim.**

`tests/testthat/test-retired-claims.R` now fails the build if any of six
retired sentences reappears in `R/`, `man/`, the vignettes, `DESCRIPTION`,
`README.Rmd`, `NEWS.md` or the oracle manifest. Entries are retired
*sentences*, not keywords, because the package is often obliged to discuss a
rule it no longer uses; what it may not do is assert one in the present tense.

It caught `aciR.Rmd` on its first run, which is the item F1 names.

## F1 to F7

| | |
|---|---|
| F1 | Done. `aciR.Rmd` scope note trimmed to one sentence; `validation.Rmd` "What is not validated" replaced by two lists, graded and open, as specified |
| F2 | Done. `plot.aci_cir()` sets `mfrow` with `on.exit`, matching `plot.aci()`; the example's outer `par()` removed, since it was the workaround |
| F3 | Done. `@param margin`, the `.aci_simpson` roxygen, and the file-header cubic claim, which your note did not list but is the same defect |
| F4 | Done. `as.data.frame.aci_cir()` carries `monotone` |
| F5 | Done, with the set named |
| F6 | Done. The Rd example converts a time to `last_idx` exactly as your two lines do |
| F7 | Done, to all four of your constraints |

Vignettes: the KL display formula added to the dyad on-ramp; a
*Reading the causal influence range* subsection added to
*Assumptions and interpretation* covering last-exit, the two objective ranges,
and censoring; the maintainer-procedure paragraph removed from
*Validation and the independent oracle*.

Not taken from §5.1: the optional `plot(rng)` demo in the dyad article. On a
short dyad record the range saturates completely -- at `n = 600`,
`window = 50:250` it is 201 of 201 -- so every point would render hollow and
teach a first-time reader that the range never resolves. It works at
`n = 2500` (11 resolved, 6 censored). If it goes in, it needs the longer
record; it did not go in this pass.

## Your Q1 and Q2 answers

Both adopted. Q1's decisive point was one I had missed: now that a censored
time returns a bound, `margin` governs the flag rather than whether a number is
reported, which makes a slightly wrong 0.1 cheap and removes the urgency from
per-model defaults. The calibration figure stays on the cairn's review trigger,
with your plot specification recorded.

Q2's constraint that the fixture must go through `.aci_online_aux_mv` rather
than a hand-built `E` array is the one that made the test worth having.

## One thing back to you, from F7

Building that fixture surfaced something bearing on S6. You measured spectral
radii in `[0.9920, 0.9951]` on your two-dimensional fixture -- all below one --
and noted that stopping while a non-normal product is still in a
transient-growth regime is "possible in principle and not visible on the
fixtures I ran".

It is visible on this one. The first two factors have spectral radius **1.200
and 1.100**, settling to 0.700 by about step 20. The cause is mundane: the
filter covariance starts at `R0 = 0.1` rather than at its fixed point, so
`S_yoS_y R_f^{-1}` is briefly large enough to flip the sign of the drift term
in `E`. The paper's bound describes the settled filter.

The consequence for S6 is concrete. A `lag_cap` derived from a contraction rate
estimated over all steps -- your `quantile(radius, 0.9)` sketch -- would be
contaminated by the transient on any system started away from its covariance
fixed point, which is every system a user supplies. The current
`max(abs(d)) < tol` test cannot have this failure, because the accumulated
product is the thing that contracts and is the thing tested. That is now
asserted in `test-online-truncation.R` rather than left as a remark.

If S6 is ever built, the rate wants estimating over the settled portion, and
the settling point itself needs detecting.

---

## Disclosure status — RELEASED 2026-08-14

The maintainer has released the hold. The two observations about the reference
implementation's own behaviour -- that `noisy_predator_prey_model.m` computes
direction one's smoother from direction two's filter when run top to bottom,
and that exactly one of the five ENSO scripts ships with conditional ACI
enabled -- may now be published with the rest of the project's findings.

Both are reported as properties of published code, reproducible from the
manifests in `oracle/parity/`, and neither is a defect when the scripts are
used as their own instructions direct. Contacting the authors remains a
separate decision.
