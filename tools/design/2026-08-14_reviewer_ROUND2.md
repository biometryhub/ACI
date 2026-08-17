# aciR — second-round review

**Date** 2026-08-14 · **Responding to** `tools/design/2026-08-14_reviewer_RESPONSE.md`
· **Against** the tree as it stands after that pass

This is a check of the response, answers to the two questions it asked, a
short list of remaining work, and a vignette review. Vignette notes are in
§5, not a separate file.

Local gates I re-ran: `testthat::test_local()` **0 failures**; CIR residual
against the committed MATLAB CSV, reported region, **max 1.37e-14** (was
6.97e-07 / 4.58e-09). I did not re-run `R CMD check --as-cran`.

---

## 1. Verdict on the response

The three rulings were the right ones, and the work matches them. The
Simpson closure is now the reference's rule, derived rather than pasted, and
the residual is gone. Censored ranges return a bound. The comment on `tol`
no longer overclaims. A1 through A7 and S1, S2, S4, S8, S11 are in the
code, not only in the table.

Two things in the response are not quite true, and the vignettes were left
behind.

1. **A1 is not done.** `DESCRIPTION`, `README.Rmd` (mostly),
   `assumptions.Rmd`, and the oracle-manifest notes were updated. The two
   vignettes a user actually opens were not. *Assimilative causal inference
   on the nonlinear dyad model* still says "What aciR still does not
   implement is the vector and conditional case". *Validation and the
   independent oracle* still says the online smoother and CIR "on vector
   states" are "Not implemented, so not graded." That is the same class of
   defect A1 named, on the two pages that outrank `DESCRIPTION` for a
   reader.
2. **The residual is 1.37e-14, not 3.57e-15.** Section 3(a), `NEWS.md`, and
   the cairn all report 3.57e-15. Section 7 reports 1.37e-14. I measure
   **1.366e-14** as the maximum on the reported region of
   `arbitrary_cross_noise`, and **3.75e-16** at the original `j = 408` row.
   Both are machine precision. Quote the maximum. 3.57e-15 is a typical
   or single-row figure presented as the grade.

The rest of the dispositions I accept: S3 declined for a reason that holds,
S6 deferred with a real overflow hazard I had not written down, S7 declined
at this scale, S9 / S10 / M3–M7 parked with honest labels.

---

## 2. Claim-by-claim check

| Response claim | Check |
|---|---|
| Closure is Garcia / Lagrange over the last interval | **Yes.** `.aci_simpson_closure()` matches the weights I derived from `simps.m` on equal and unequal grids. The equal-spacing reduction `h/12 * (-y0 + 8 y1 + 5 y2)` is tested. |
| `objective` residual closed | **Yes, at 1.37e-14 max**, not 3.57e-15. `j = 408` is 3.75e-16. |
| Tests retargeted (cubic even / quadratic odd) | **Yes.** `test-quadrature.R` is the right split. The fourth-order test still uses only odd `n` (even interval counts), which is correct for that claim and should say so in the comment. |
| Four-valued `status`, bound not `NA` | **Yes.** 1049 / 1051 finite on the reference window; 2 `insufficient` (the last two indices, no quadrature). Tests in `test-cir.R` and `test-cir-horizon.R` assert the bound. |
| `saturated` kept as `status == "censored"` | **Yes.** |
| Horizon docs / README recipe | **Yes.** Recipe is the right three arguments. One leftover sentence in the README still opens "Current scope: a scalar observed process" and then walks it back. |
| `tol` Rd corrected | **Yes.** |
| S3 declined | **Accept.** Your reading of my probe is right: `pinv` would not have saved that run. |
| S6 overflow hazard | **Accept**, and it is a better reason to defer than "not this pass". See Q notes below. |
| A5 vector `filt` check | **Yes.** |
| A6 `aux$E_j` + `warnPartialMatchDollar` | **Yes.** `tests/testthat/setup.R` is the right shape. |
| A7 roxygen split | **Yes.** |
| M1 `monotone` + summary diagnostic | **Yes.** I also see 0% monotone on the truncated cross-noise window, same finding. |
| M2 `@section Estimand` on `aci_conditional` | **Yes.** The wording is the wording I wanted. |
| `print` / `summary` / `as.data.frame` / `plot` for `aci_cir` | **Present, two defects** — see §4. |
| A1 vignettes | **Not done** — see §5. |
| 718 tests, 0 failures | **0 failures** re-confirmed. I did not re-count 718. |

Stale surfaces the response's own table does not list:

* `R/aci-quadrature.R` lines 17–23 (the `.aci_simpson` roxygen) still
  describe the Simpson 3/8 panel.
* `R/aci-cir.R` `@param margin` still says ranges are "returned as `NA`".
  `man/aci_cir.Rd` carries that sentence into the installed manual. The
  `@returns` block was updated; the `@param` was not. A reader of `?aci_cir`
  is told both stories.
* `plot.aci_cir()` does not set `mfrow`. `plot.aci()` does. Calling
  `plot(rng)` draws the objective panel and immediately overwrites it with
  the peak panel. The example works only because it sets `par(mfrow)`
  *outside* the method.

---

## 3. The two questions

### Q1. `margin` calibration

The resolved fraction will **not** be flat in `margin`, and the right cut
will **differ by system**. Last-exit of the reporting threshold is a
property of how fast the ordered products decay, which is a property of
`(L_y, S, dt)`, not of the package. A dyad at `dt = 0.001` and ENSO at
`dt = 0.005` should not share a calibrated fraction except by accident.

That does **not** argue for per-model defaults today. Now that a censored
time returns a bound, `margin` mainly controls the *flag*, not whether the
user gets a number. A slightly wrong 0.1 is cheap: the value is still
there, marked. Per-model defaults become worth the surface only after the
calibration figure exists and shows that one number is actively misleading
on a system you ship. Until then keep 0.1, keep it documented as this
package's device, and put the figure on the cairn's review trigger as you
have.

What I would plot, when you run it: for each of dyad, predator-prey (both
directions), and ENSO, the fraction of the reporting window with
`status == "resolved"` against `margin` in `{0.02, 0.05, 0.1, 0.2, 0.4}`,
at the default (whole-record) horizon. If the curves fan out, say so and
leave 0.1. If one system is still fully censored at 0.4, *that* is the
case for a per-model default.

### Q2. Forced-truncation test

**Accept a synthetic system with `|E| ~ 0.5`.** The identity under test is
algebraic: once the running product is negligible, filling forward must
leave `mu_end` equal to a `tol = 0` walk. That identity does not care
whether the contraction came from a published model. A fixture the package
would "actually be used on" never reaches `1e-18` in a record it can hold,
so it cannot exercise the path.

Constraints I would put on the synthetic system, so the test still means
something:

* It must be an admissible CGNS and go through `.aci_online_aux_mv`, not a
  hand-built `E` array. Otherwise you are testing a different function.
* `|E| ~ 0.5` (or spectral radius ~ 0.5) is plenty. At 0.5, `1e-18` is
  reached in ~60 steps.
* Assert two things, not one: (i) `tol = 1e-18` fires (`lag_effective` or
  the fill-forward length is shorter than the record); (ii) `mu_end` and
  the filled tail of the CIR row agree with `tol = 0` to ~1e-14.
* Keep it out of the oracle manifest. It is a path test, not a grade
  against the authors.

Do not wait for a "realistic" fast contraction. That is how this path stays
untested.

---

## 4. Further suggestions

Small enough to do in the same sitting as the vignette edits. None of these
reopen a ruling.

### F1. Finish A1 in the two vignettes

See §5. This is the one that is still a factual error in shipped
documentation.

### F2. `plot.aci_cir()` should own its layout

Copy the `par(mfrow = c(2, 1))` / `on.exit` pattern from `plot.aci()`.
Without it the method's first panel is invisible in interactive use. The
example setting `mfrow` outside is a workaround, not a design.

### F3. Align the CIR manual with the new return

* `@param margin`: delete "returned as `NA`"; say the time is marked
  `status = "censored"` and the value is a lower bound.
* `.aci_simpson` roxygen: describe the last-interval quadratic, not 3/8.
* Rebuild `man/`. The installed `@returns` is current; the installed
  `@param margin` is not.

### F4. `as.data.frame.aci_cir()` should carry `monotone`

It is a per-time column, it is the condition under which the two objective
ranges are the same functional, and it is the one thing `summary()` tells
the reader to look at. Leaving it off the tidy frame means the diagnostic
you just added is available only through the print method.

### F5. Quote 1.37e-14 as the grade

One number, the maximum, in `NEWS.md`, the cairn, and any future parity
table. 3.57e-15 can stay as a parenthetical typical if you want, but it
should not be the headline.

### F6. One Rd example that sets `horizon` from a time

You already noted this as not done. Two lines:

```r
last_idx <- as.integer(round((time_end + lookahead) / dt)) + 1L
aci_cir(..., horizon = last_idx)
```

That is the conversion every reader of the paper will have to write.

### F7. The forced-truncation test (your Q2)

Do it on the synthetic `|E| ~ 0.5` system, with the constraints in §3.
Highest-value remaining test, as you said.

I would still not start S6, S9, S10, or M3–M7 in this sitting. S6 in
particular should wait on a measurement of the blocked-rescaling sketch
against the current offset loop, on `n` in `{5e3, 1e4, 3e4}`. The overflow
argument is enough to keep it off the critical path.

---

## 5. Vignettes

The request was: clearer, more concise, well displayed, formulas present,
not overdone. All three vignettes are already in a good house voice. What
they need is a factual pass and a light structural trim, not a rewrite.

No Unicode math, no em dashes. That side is clean.

### 5.1 *Assimilative causal inference on the nonlinear dyad model*

This is the on-ramp. It should show the workflow, one formula for the
metric, and a short pointer at the range. It should not deny capabilities.

**Must fix**

* Lines 221–222. Delete "What aciR still does not implement is the vector
  and conditional case, which `NEWS.md` and `API_STABILITY.md` record."
  Both are implemented. Referring to those filenames in vignette prose is
  also the wrong register; if a pointer remains, name the documents by
  title.
* The scope note that remains can stay, and should be one sentence: the
  range is `aci_cir()`, it answers a different question, and this article
  does not walk through it.

**Should add, small**

* After "The causal question", one display formula — the same KL that
  *Assumptions and interpretation* already carries. The dyad article
  currently describes the metric only in words. A reader who starts here
  and never opens the second vignette never sees what is being plotted.

$$
\mathrm{ACI}(t)
  = \mathrm{KL}\!\bigl(p(y_t \mid x_{0:T}) \,\|\, p(y_t \mid x_{0:t})\bigr).
$$

  That is enough. Do not add the Gaussian closed form here; it belongs in
  the validation article, where it already is.

* Optionally, six to eight lines at the end of the metric section showing
  `aci_cir()` on a short window and `plot(rng)`, with one sentence that a
  hollow point is a lower bound. Do not build a second worked example
  around it. The current reason for skipping CIR ("a different question")
  is fine if the function at least appears.

**Do not add**

* A derivation of the filter. A second pair of figures. A CIR heatmap. Any
  mention of Simpson, `horizon`, or `margin` beyond the one-sentence
  pointer. Those belong in the README recipe and the Rd.

**Tighten**

* "The causal question" and "The filter and the smoother" overlap. Keep
  the first as the estimand (plus the formula) and let the second be only
  the figure. You can lose about a paragraph without losing a step.

The two ggplot panels are the right display density. Leave them.

### 5.2 *Assumptions and interpretation*

The best of the three. The estimand formula is the right size. The four
things a peak does not support are the right list. The A1 fix here was
actually done.

**Should add, small**

* A short subsection on the range, because `status` and last-exit are now
  interpretive facts a user can get wrong. Four or five sentences, no new
  figure:

  > The causal influence range is a last-exit time, not a duration. The
  > efficient `objective` and the threshold-average `objective_exact` are
  > the same functional only when the comparison sequence decreases with
  > lag; `summary()` reports how often that happens. A time marked
  > `censored` is a lower bound, not a missing value.

* One sentence in "Conditional questions" pointing at `?aci_conditional`
  for the estimand note you just wrote. Do not duplicate the Rd paragraph.

**Do not add**

* Granger / transfer-entropy comparison (M7, correctly deferred).
* A residual / PIT section (M3, same). A pointer that those are not in the
  package is already implied by "there is no goodness-of-fit test".

**Tighten**

* Almost nothing. This article is allowed to be the long one.

### 5.3 *Validation and the independent oracle*

Conceptually the right article: three oracles, a manifest, and an honest
negative space. The negative space is now stale, and the heading "What is
not validated" mixes things that are graded with things that are not.

**Must fix**

Rewrite that section as two lists.

*Still open*

* Simulation paths (R's generator vs MATLAB's). Keep as written.
* `aci_conditional()` against the one ENSO script that has it enabled.
* Predator-prey CIR against the authors (filter / smoother / metric are
  graded; the range is not).
* Four of five ENSO extracts captured and not gated.
* Model adequacy (cannot be graded).

*Graded, with the scope named*

* Time-varying `L_y`, both predator-prey directions, to 1e-6. This bullet
  currently sits under "not validated" and then says it *is* graded.
* Online-smoother and CIR cross-noise, against `cir_cross` (a second
  implementation, not the authors). The current bullet says this path is
  unverified. That has been wrong since the cross-noise CIR fixture
  landed.
* Vector filter / smoother / metric, against the ENSO fixture and the
  block-diagonal Riccati identities.
* Vector online smoother and CIR, by collapse onto the scalar path and by
  the independent MATLAB transcription of (3.5)–(3.7). Not against the
  authors' ENSO scripts.

Delete the sentence "The online smoother and the causal influence range on
vector states. Not implemented, so not graded."

**Display**

* Keep the KL closed form and the Riccati. That is the right amount of
  mathematics for this article: two formulae the oracles actually grade,
  not a derivation.
* The agreement table (dyad 4.57e-14, cross 1.38e-14) is the right shape.
  When you next refresh it, add a CIR row: peak 9.28e-15, objective
  1.37e-14, under the reference conventions. One table, current numbers.
* The provenance section can lose the paragraph that starts "Fixtures are
  never regenerated in place." That is maintainer procedure, not a
  validation claim a reader needs.

**Do not add**

* Extracts of `simps.m`. A second table of every fixture in the manifest.
  Screenshots of MATLAB.

### 5.4 What "not overdone" means here

One display formula per idea, and only if the idea is used. The dyad
article needs the KL. The assumptions article already has it. The
validation article needs the KL closed form and the Riccati. Nobody needs
the Garcia weights, the last-exit suffix-max, or the Euler residual
expansion in a vignette.

Figures: two in the dyad article (already). Zero more unless you add the
optional six-line `plot(rng)` call, which reuses the new method and does
not invent a ggplot. The assumptions and validation articles should stay
unillustrated.

Length: the dyad article can drop a paragraph. The validation "not
validated" section should get shorter by becoming two short lists. The
assumptions article stays.

---

## 6. Suggested order for the leftover work

1. F1 / §5.1 and §5.3 factual fixes (vignettes). Same sitting.
2. F2, F3, F4 (plot layout, Rd, `as.data.frame`). Same sitting.
3. F5 (one number in NEWS / cairn). Same sitting.
4. F6 (Rd `horizon` example) and F7 (synthetic truncation test). Next
   sitting.
5. Everything the response already deferred, in the order it gave.

I would not reopen S3, S6, or S7 on the basis of this pass.

---

## 7. Handling note

Unchanged. The two observations about the reference implementation remain
internal.
