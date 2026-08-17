# aciR — scientific and expositional review

**Review date:** 16 July 2026  
**Package version:** 0.1.0  
**Reviewed at commit:** `73196fd` (main; the state every finding below refers to)  
**Scope:** mathematical fidelity of the implemented method; scientific interpretation and estimand; quality of vignettes as exposition; coverage relative to Andreou, Chen & Bollt (2026, *Nature Communications*); gaps, risks of misreading, and recommendations for scientific communication.  
**Companion document:** [package quality review](2026-07-16_package_quality_review.md). Finding IDs are prefixed `SV-` here and `PQ-` there so cross-references stay unambiguous.  
**Primary sources:** package vignettes (`aciR`, `assumptions`, `validation`), Rd/README, numerical core, MATLAB reference (`dyad_interaction_model.m`), method paper PDF (`literature/s41467-026-68568-0.pdf`).

---

## Executive verdict

Scientifically, `aciR` is a **faithful, carefully scoped reimplementation of the scalar CGNS core of assimilative causal inference**, not a full reproduction of every capability in the method paper. The mathematics of the filter, smoother, and relative-entropy metric agree with the authors’ MATLAB reference at floating-point noise level. The package’s **interpretive discipline is stronger than much of the surrounding literature rhetoric**: the assumptions vignette states a precise model-conditional estimand and explicitly denies interventional and model-free causal claims.

Expository quality of the three vignettes is high: clear prose, honest negative space, and a validation narrative that treats oracle *scope* as part of the scientific claim. The main scientific gaps are **declared** (CIR not implemented; vector/conditional ACI absent; time-varying self-drift roadmapped) rather than silently papered over — which is the right failure mode for a research preview.

**Overall grade (science & exposition): A− for what is claimed; B for completeness vs the full paper.**

| Dimension | Grade | One-line judgment |
|---|---|---|
| Fidelity of CGNS filter/smoother/metric | A+ | Term-level match to reference; ~1e−14 absolute error |
| Estimand statement & causal hygiene | A | Model-conditional KL; interventional claims denied |
| Dyad vignette (worked example) | A− | Complete workflow; perfect-model reading caveats present |
| Assumptions vignette | A | Best scientific artefact in the package |
| Validation vignette | A | Oracle epistemology done properly |
| Coverage of paper’s full method | C+ / roadmap | CIR, conditional ACI, ENSO/real data, vector states missing |
| Risk of overclaim / mis-citation | Medium | Name “causal inference” still invites misreading |
| Pedagogical path for non-specialists | B+ | Strong for DA-literate readers; thin on derivation detail |

---

## The method, as the paper defines it

Andreou, Chen & Bollt (2026) frame causality as an **inverse problem in uncertainty quantification**. Given a dynamical model for \((x, y)\) and a continuous observation of the path of \(x\), they ask whether future information in \(x\) reduces uncertainty about the latent state \(y(t)\). Operationally:

1. **Filter** posterior \(p^f_t = p(y_t \mid x_{0:t})\).  
2. **Smoother** posterior \(p^s_t = p(y_t \mid x_{0:T})\).  
3. **ACI metric** — relative entropy  
   \[
   \mathcal{P}(p^s_t, p^f_t)
   = \mathrm{KL}\!\left(p^s_t \,\|\, p^f_t\right)
   = \int p^s_t \log\frac{p^s_t}{p^f_t}\,dy
   \]
   (paper eq. (7); nonzero values taken as evidence that \(y(t)\) influences future \(x\) under the model).  
4. **Causal influence range (CIR)** — subjective thresholds on lagged smoothers, plus an objective integral over thresholds (paper eqs. (8)–(9); \(O(N^2)\) naively).  
5. **Conditional ACI** with non-target variables via inflated observational uncertainty on \(x_B\) (paper §4.3).

Closed-form recursion exists when the system is a **conditional Gaussian nonlinear system (CGNS)**: drift of \(y\) linear in \(y\), noise independent of \(y\), so filter/smoother posteriors stay Gaussian.

The flagship scalar example is the **nonlinear dyad** with intermittent extremes:

\[
\begin{aligned}
\mathrm{d}x &= (-d_x x + \gamma x y + F_x)\,\mathrm{d}t + \sigma_x\,\mathrm{d}W_x, \\
\mathrm{d}y &= (-d_y y - \gamma x^2 + F_y)\,\mathrm{d}t + \sigma_y\,\mathrm{d}W_y,
\end{aligned}
\]

with parameters \(d_x = d_y = 0.5\), \(\gamma = 2\), \(F_x = 0.5\), \(F_y = 1\), \(\sigma_x = 0.5\), \(\sigma_y = 1\) — exactly the package defaults.

---

## Mathematical fidelity of the package

### What is implemented correctly

| Piece | Assessment |
|---|---|
| Dyad SDEs and default parameters | Match paper and MATLAB reference |
| CGNS coefficient map (`L_x = γx`, `f_x = F_x − d_x x`, `L_y = −d_y`, `f_y = F_y − γx²`) | Correct |
| Euler–Maruyama simulation (independent noise) | Correct structure; generator differs from MATLAB by design |
| Forward conditional Gaussian filter (Euler discretisation) | Matches MATLAB update; scalar \(2 L_y R\) form of matrix Lyapunov term |
| Backward smoother with \(A_j\), \(B_j\), cross-noise transport | Matches MATLAB |
| Scalar Gaussian KL (signal + dispersion) | Algebraically correct form of \(\mathrm{KL}(\mathcal{N}_s\|\mathcal{N}_f)\); numerically improved near unit covariance ratio |
| Terminal identity \(p^s_T = p^f_T\) ⇒ \(\mathrm{ACI}(T) = 0\) | Exact by construction; tested |
| Zero-information limit (no coupling, no cross-noise) | Continuous-time identity recovered at \(O(\mathrm{d}t^2)\) under Euler — correctly tested as a *rate*, not false equality |
| Stationary Riccati fixed point with \(S_{yx} \neq 0\) | Independent of both R and MATLAB; strong algebraic grounding of cross-noise terms |

**Oracle evidence** (scientific, not merely software):

| Oracle | What it grounds | Observed max \|error\| |
|---|---|---|
| Dyad MATLAB fixture (authors’ reference lineage) | Full transient, independent noise, state-dependent \(L_x\) | \(\approx 4.57 \times 10^{-14}\) |
| Cross-noise MATLAB harness | Cross-covariance terms over full transient | \(\approx 1.38 \times 10^{-14}\) |
| Analytic Riccati / smoother fixed points | Cross-noise *exactly* at stationarity | machine precision |

This is exceptional grounding for a scalar method package.

### Discretisation and approximation structure (scientifically important)

The package (and the reference) use **first-order Euler** for the continuous-time filter/smoother ODEs. Consequences, correctly reflected in docs and tests:

1. Continuous-time identities hold only approximately (zero-information residual \(\propto \mathrm{d}t^2\) because KL is quadratic in posterior discrepancies that are themselves \(O(\mathrm{d}t)\)).  
2. Too large \(\mathrm{d}t\) can destroy positivity of the covariance even for an admissible model; the software fails closed rather than inventing a KL.  
3. The method paper assumes continuous, noise-free observation of \(x\); the package enforces a regular complete grid and does not invent discrete-time observation noise. That is faithful to the paper’s simplified setup, not a hidden defect — but it **is** a limitation for real instruments (paper itself flags discrete-time extensions as future work).

**Scientific recommendation:** keep stating “Euler discretisation of continuous-time CGNS posteriors” in any methods section that cites the package; do not imply continuous-time exactness of numerical output.

### Gaussian KL formula check

For one-dimensional Gaussians,

\[
\mathrm{KL}(\mathcal{N}(m_s,R_s)\,\|\,\mathcal{N}(m_f,R_f))
= \tfrac12\Biggl[
  \frac{(m_s-m_f)^2}{R_f}
  + \frac{R_s}{R_f}
  - 1
  - \log\frac{R_s}{R_f}
\Biggr].
\]

Package: signal term \(\tfrac12 (m_s-m_f)^2/R_f\) plus dispersion \(\tfrac12\bigl((R_s/R_f-1) - \log(R_s/R_f)\bigr)\) via `log1p`. **Correct**, and better-conditioned near \(R_s \approx R_f\) than the naive rearrangement.

Direction \(\mathrm{KL}(p^s \| p^f)\) matches the paper’s integrating density \(p^s\). Reversing the arguments would be a different (wrong for ACI) quantity.

---

## Estimand, causality language, and scientific hygiene

### What the package claims (and should claim)

The assumptions vignette states the estimand as

\[
\mathrm{ACI}(t)
= \mathrm{KL}\!\bigl(
  p(y_t \mid x_{0:T},\mathcal{M})
  \;\big\|\;
  p(y_t \mid x_{0:t},\mathcal{M})
\bigr),
\]

i.e. a **property of a model-conditional posterior pair**, in nats, non-negative, asymmetric, and **not** an interventional effect.

This is scientifically cleaner than the paper’s occasionally stronger prose (“\(y\) is identified as the cause of \(x\) at time \(t\)” when relative entropy is positive). The package correctly emphasises:

| Claim a peak supports | Claim a peak does **not** support |
|---|---|
| Under \(\mathcal{M}\), future \(x\) is unusually informative about \(y_t\) | \(\mathcal{M}\) is correct |
| Timing of information flow under a specified mechanism | Intervention would change \(x\) |
| — | Causality discovered from data alone |
| Comparability across time within one model/record | Comparability across models, parameterisations, sampling rates |

**This is the right scientific posture** for an R package that will be cited by users outside data-assimilation communities. It prevents the common category error of treating a model-dependent information functional as a discovery of causal structure.

### Remaining risks of misinterpretation

| Risk | Severity | Mitigation already present | Residual gap |
|---|---|---|---|
| Package name “causal inference” vs potential-outcomes literature | Medium | Assumptions vignette; README scope section | First Google hit / CRAN title still primes the wrong literature |
| Perfect-model demos read as validation of real-world causality | Medium | Dyad vignette ends with explicit caveat; assumptions vignette | No imperfect-model demo yet (paper has preliminary ENSO real-data discussion) |
| ACI peak ≅ “y causes x” without model literacy | Medium | Asymmetry of the metric is stated | No short “not Granger / not TE / not do-calculus” comparison table in the dyad vignette |
| Units (nats) ignored when comparing runs | Low | Assumptions § comparability | Could show a wrong comparison anti-example |
| Parameters treated as known when they were fit | Medium | “Known parameters” assumption explicit | No sensitivity-to-parameters vignette chunk |

### Relation to classical causal frameworks (expositional opportunity)

The package correctly *distances* ACI from interventional identification but does not yet give a compact positioning table. Scientifically useful contrasts for a future short section:

| Framework | Question | ACI |
|---|---|---|
| Granger / transfer entropy | Does past of \(y\) improve prediction of \(x\)? | Inverse: does future of \(x\) reduce uncertainty in \(y\)? |
| Pearl / potential outcomes | What happens under \(\mathrm{do}(\cdot)\)? | Not answered |
| CCM / dynamical systems reconstruction | Shared attractor / state reconstruction | Requires model equations, not just embedding |
| Ensemble information transfer (Liang–Kleeman et al.) | Forward model ensembles | Single path + model, Bayesian DA |

Adding half a page of this positioning would raise the dyad vignette from “good walkthrough” to “good scientific onboarding.”

---

## Vignette-by-vignette review

### 1. `aciR.Rmd` — Assimilative causal inference on the nonlinear dyad model

**Role:** flagship walkthrough.  
**Grade: A−**

**Strengths**

- Opens with a plain-language statement of filter / smoother / metric before any code.  
- Writes the dyad SDEs with parameters matching the paper.  
- States the causal question as *when* coupling is informationally active, not merely *whether* a static edge exists — aligned with the paper’s intermittency motivation.  
- Shows filter vs smoother vs true \(y\) (simulation ground truth) — pedagogically excellent for DA intuition.  
- Surfaces `summary()` diagnostics (min cov, terminal residual).  
- Ends with a clear pointer to the assumptions article and a general components path for non-dyad systems.  
- References the method paper with DOI.

**Expositional weaknesses**

| ID | Severity | Issue | Suggestion |
|---|---|---|---|
| SV-V1 | Medium | **Metric and signal plotted on a shared y-axis** in the ggplot chunk. For the default seed the ranges happen to be comparable (\(x \in [-0.12,1.8]\), ACI \(\in [0,2.46]\)), so the figure “works” — but the design does not generalise and can visually equate incommensurable quantities. | Prefer two panels (as `plot.aci` does) or a secondary axis with an explicit dual-scale warning. |
| SV-V2 | Medium | **No quantitative filter/smoother error summary** (e.g. RMSE of filter mean vs true \(y\), smoother improvement). Readers see a plot but leave without a number. | Add 2–3 lines of RMSE / mean ACI at burst times. |
| SV-V3 | Low | **Physical time span is short** (\(T \approx 5\)) vs paper figures that highlight \(t \in [5,25]\). Pedagogically fine; for visual alignment with the paper, a longer run or a note about the paper’s plotting window would help. | Optional longer seed or window note. |
| SV-V4 | Low | **No mention of CIR** even as “not in this package yet,” except via roadmap elsewhere. A one-sentence “the paper also defines CIR; not implemented here” prevents readers from thinking the vignette *is* the whole method. | One sentence + link to NEWS/API_STABILITY. |
| SV-V5 | Info | ggplot is optional (`requireNamespace`); base path exists via `plot()`. Good dependency hygiene. | — |

**Scientific content density:** appropriate for a first vignette. It does not derive the filter equations (correctly deferred).

---

### 2. `assumptions.Rmd` — Assumptions and interpretation

**Role:** scientific requirements document and anti-misuse charter.  
**Grade: A (best artefact)**

**Strengths**

- Explicit purpose: the other vignette’s “peak = coupling” reading is conditional.  
- Estimand written formally with \(\mathcal{M}\) in the conditioning.  
- Assumptions listed as checkable software contracts vs uncheckable scientific ones:  
  - model correctness (uncheckable; most dangerous)  
  - CGNS structure  
  - known parameters  
  - scalar state  
  - time-invariant \(L_y\)  
  - regular complete sampling  
  - Euler discretisation  
- Live `error = TRUE` examples for missing data and unstable \(\mathrm{d}t\) — teaching by refusal.  
- Sensitivity section (`dt`, `R0`, `mu0`) with concrete reporting advice.  
- Four denied claims (model truth, intervention, discovery from data alone, cross-model comparison).  
- Publication checklist of diagnostics to report.  

This vignette should be **required reading** before any applied citation of the package. As scientific exposition, it is unusually responsible.

**Minor improvements**

| ID | Severity | Issue | Suggestion |
|---|---|---|---|
| SV-A1 | Low | **No formal derivation** of why CGNS ⇒ Gaussian closed form (pointer to SI of the paper would suffice). | One paragraph + SI reference. |
| SV-A2 | Low | **“Units of nats”** stated once; a one-line conversion note to bits (\(/\log 2\)) helps information-theory readers. | Optional. |
| SV-A3 | Medium | **Model misspecification** is named as the chief danger but not illustrated. A short anti-example (wrong \(\gamma\), or fitted-as-true parameters) would make the abstract warning visceral. | Strong candidate for a future vignette chunk or fourth vignette. |
| SV-A4 | Low | Links to validation vignette and paper; could also link `API_STABILITY.md` roadmap for CIR. | Cross-link. |

---

### 3. `validation.Rmd` — Validation and the independent oracle

**Role:** epistemology of numerical evidence.  
**Grade: A**

**Strengths**

- Correct philosophical setup: author-written tests prove self-consistency, not method fidelity.  
- Three-oracle design with **scope**, not just tolerance:  
  1. authors’ dyad fixture,  
  2. cross-noise fixture (necessary because reference scalar models zero \(S_{yx}\)),  
  3. analytic identities independent of both languages.  
- Explicit “line coverage can lie” lesson about annihilated cross-covariance terms — this is research-grade methodological literacy.  
- Zero-information test explained as rate-in-\(\mathrm{d}t\), not false equality.  
- Manifest provenance, dual hashes, byte-for-byte MATLAB reproduction note, non-regeneration policy.  
- Clear “what is not validated” list (simulation paths, time-varying \(L_y\), CIR, vector states, model adequacy).  

**Minor improvements**

| ID | Severity | Issue | Suggestion |
|---|---|---|---|
| SV-O1 | Low | Error table is excellent; could add **which code paths** each row exercises (e.g. `aux`, `A_j`, transport). | Optional column. |
| SV-O2 | Low | Does not show the **algebraic form of the scalar Gaussian KL** used in the metric. | One displayed equation would connect validation to science. |
| SV-O3 | Info | Cross fixture “not authors’-reference grounding” is stated honestly — keep that sentence forever; it is a trust asset. | — |

---

## Coverage relative to the full method paper

| Paper capability | In aciR 0.1.0? | Scientific impact of absence |
|---|---|---|
| Scalar CGNS filter / smoother / ACI metric | **Yes**, oracle-graded | Core scientific claim of the package — solid |
| Nonlinear dyad example | **Yes** | Matches paper § Results opening example |
| Noise cross-covariance | **Yes** in core; simulation no | Filtering path solid; generative path incomplete |
| Causal influence range (subjective/objective) | **No** (roadmap) | **Major** relative to paper’s selling points (Fig. 1 whiskers, objective CIR) |
| Online smoother for CIR | **No** | Blocks CIR without a separate design |
| Noisy predator–prey / time-varying \(L_y\) | **No** | Blocks a second paper example class |
| Conditional ACI / non-target variables | **No** | Blocks ENSO-style multi-variable claims |
| Vector / high-dimensional CGNS | **No** | Scalability claims of the paper are out of package scope |
| Real-world ENSO observational case | **No** (MATLAB repo only) | Limits “applied science” demonstration |
| Discrete-time / gappy observations | **No** (explicitly rejected) | Faithful to paper simplicity; limits practice |
| Ensemble / approximate DA for non-CGNS | **No** | Out of current design (correct) |

**Judgment:** as a package named for the *method*, 0.1.0 implements the **computational heart** of ACI for the scalar CGNS case and the paper’s main pedagogical system. It does **not** yet implement the paper’s second signature contribution (CIR) nor the multi-variable conditional machinery used in the ENSO showcase. The roadmap and assumptions vignette make this transparent; marketing and titles should stay equally transparent.

---

## The mathematics “beyond” the package

This section evaluates the conceptual layer users must bring, and how well the package prepares them.

### What the package teaches well

1. **Filter vs smoother as asymmetric information sets** — central to ACI, well explained.  
2. **Relative entropy as mean + dispersion information**, not Shannon-entropy difference alone — aligned with the paper’s rationale.  
3. **Model-conditionality of every numerical claim.**  
4. **Discretisation as a scientific assumption**, not a software detail.  
5. **Oracle scope as part of scientific validation.**

### What remains under-taught (relative to a methods-literate user)

| Topic | Why it matters | Current state |
|---|---|---|
| Why CGNS stays Gaussian | Justifies closed form vs particle filters | Asserted, not derived |
| Continuous-time observation of \(x\) with process noise vs measurement noise | Users often import Kalman intuition with separate \(R\) measurement noise | Implicit; could confuse |
| Connection of ACI peaks to antidamping threshold \(y > d_x/\gamma\) | Paper’s physical mechanism for dyad extremes | Mentioned lightly as “bursts”; threshold not plotted |
| CIR mathematics (threshold integral) | Paper’s objective influence length | Absent (feature missing) |
| Nil-causality / conditional principles under ACI | Paper SI verification topics | Absent |
| Model error sensitivity | Paper Discussion prioritises this as future work | Named as risk only |

### Physical reading of the dyad (opportunity)

Paper Fig. 1 ties large ACI to **onset/peak of extremes** when \(-\,d_x + \gamma y\) becomes anti-damping, and short CIR after the peak when filter already sees high SNR. The package vignette shows co-located spikes of \(x\) and ACI but does not draw the antidamping threshold \(d_x/\gamma = 0.25\) on \(y\) or discuss the post-peak collapse of ACI. Adding that physical narrative would:

- deepen scientific exposition without new algorithms,  
- align the R demo with the paper’s storytelling,  
- reinforce that ACI peaks are mechanism-timed, not generic “nonlinearity detectors.”

---

## Scientific risks for users and for the package’s reputation

| ID | Severity | Risk | Mitigation |
|---|---|---|---|
| SV-R1 | High | Citing package as establishing real-world causality from a fitted CGNS without model critique | Assumptions vignette already forbids this; enforce in README and any paper templates |
| SV-R2 | Medium | Treating MATLAB–R agreement as proof the *method* is right, not that the *transcription* is | Validation vignette already distinguishes; keep that framing in any methods paragraph |
| SV-R3 | Medium | Publishing results without \(\mathrm{d}t\) / \(R_0\) sensitivity | Checklist in assumptions vignette — make it a `summary()` “report card” later? |
| SV-R4 | Medium | Comparing ACI magnitudes across models or sampling rates | Explicitly denied in assumptions; still easy to do with `as.data.frame` |
| SV-R5 | Medium | Overstating package completeness vs Nature Communications method (especially CIR) | Title/DESCRIPTION say CGNS ACI; add “scalar CGNS metric; CIR not yet” in abstract-like DESCRIPTION if space |
| SV-R6 | Low | Parameter uncertainty ignored after offline estimation | Document; future propagation is research, not a bug |

---

## Recommendations (science & exposition)

### Highest value, low effort

1. **Dyad vignette:** two-panel metric plot; one sentence on CIR-not-implemented; optional antidamping threshold on a \(y\) panel.  
2. **Assumptions vignette:** short misspecification anti-example (wrong \(\gamma\) or swapped roles).  
3. **Validation vignette:** display the scalar Gaussian KL formula once.  
4. **Fix `n_clamped` wording/count** so publication diagnostics do not overstate round-off clamping (see package review PQ-S1).  
5. **DESCRIPTION / README one-liner** on scope: scalar CGNS ACI metric; CIR and vector states roadmap.

### Medium effort, high scientific value

6. **Misspecification / sensitivity vignette** (or section): vary \(\gamma\), \(d_y\), \(\mathrm{d}t\); show that peaks move or inflate under wrong models.  
7. **Positioning paragraph** vs Granger, transfer entropy, CCM, and interventional causality.  
8. **Correlated-noise simulation** so the graded cross path is also a generative path users can explore scientifically.  
9. **Longer dyad window** matching paper Fig. 1 aesthetics for visual cross-check by human readers.

### Large, correctly deferred until oracles exist

10. **Causal influence range** (`aci_cir`) with independent fixture — paper’s second headline deliverable.  
11. **Time-varying \(L_y\)** + predator–prey constructor.  
12. **Conditional ACI** for non-target variables — required for multi-state climate demos.  
13. **Real ENSO walkthrough** using public reanalysis (package or companion data package), with imperfect-model honesty.

---

## Suggested language for methods sections citing aciR

A citation-safe sketch consistent with the package’s own hygiene:

> We compute the assimilative causal-information metric of Andreou, Chen & Bollt (2026) using the R implementation aciR (Moldovan, 2026). For a user-specified conditional Gaussian nonlinear system \(\mathcal{M}\), aciR evaluates the relative entropy between the smoother and filter posteriors of the unobserved state given a regularly sampled complete path of the observed process. Reported peaks are therefore statements about information flow *under \(\mathcal{M}\)*, not tests of model adequacy or interventional causal effects. Numerical agreement of the scalar core with the authors’ MATLAB reference is within \(10^{-13}\) absolute error on packaged fixtures (gate \(10^{-6}\)). The causal influence range of the original method is not computed by aciR 0.1.0.

---

## Closing judgment

From a **scientific and expositional** standpoint, `aciR` 0.1.0 is a high-integrity research preview of **scalar CGNS assimilative causal inference**. The mathematics of the implemented core is solid and independently verified. The assumptions and validation vignettes are exemplary scientific communication: they teach the estimand, the failure modes, and the limits of numerical evidence.

What the package does *not* yet do — CIR, conditional multi-variable ACI, real-data imperfect-model showcases — is the difference between “reproduces the paper’s computational heart” and “reproduces the paper’s full scientific toolkit.” That gap is documented, not hidden. Closing it under the same oracle discipline will determine whether aciR becomes the reference R implementation of the *method* or remains a reference implementation of the *metric core*.

For users and collaborators: trust the numbers on graded paths; treat every causal sentence as model-conditional; read the assumptions vignette before the abstract of any derived paper.

---

**Post-review action (2026-07-16, same day).** SV-V1 (two-panel metric figure), SV-V4 (CIR-not-implemented sentence), SV-O2 (displayed scalar Gaussian KL formula) and SV-R5 (`DESCRIPTION` scope sentence) were addressed immediately after this review, alongside PQ-S1/PQ-T3 from the companion review and a licence-attribution fix identified in a follow-up meta-review. See `NEWS.md` (development version) and the commits following `73196fd`.

*See also:* [package quality review](2026-07-16_package_quality_review.md).
