# aci 0.0.21

## Features

* `backward_cir()` results report the andreou2026cir Section 2.3.4 validity
  gate: `meta$baseline` (the initial-time information deficit), `meta$terminal`
  (the metric at the `T - dt` end of the grid) and `meta$above_baseline`
  (eq. 21; Remark C.2). `print.cir_result()` prints the anchor count above
  baseline, and `cir_table()` backward rows carry an `above_baseline` column.
* `lt_contraction_certificate()` evaluates the checkable per-step contraction
  condition of the online smoother update (andreou2026smoother eqs.
  3.18-3.19), reporting per-step `lambda_min`, `enorm` and `rho_E` with
  `gamma` and `condition_318` attributes. Diagnostic only: `lt_tail_bound()`
  is unchanged and remains a heuristic estimate.
* `forward_cir()` and `backward_cir()` gain `simpson_close`, selecting how
  composite Simpson closes the leftover interval on an even-length grid. The
  new default `"quadratic"` fits a quadratic through the last three points,
  following `simps.m` in the ACI reference code; `"trapezoid"` is the
  pre-0.0.21 rule and reproduces results reported by earlier versions. The two
  differ by `-(y_n - 2 y_{n-1} + y_{n-2}) / 12`, which is zero on a flat or
  linear tail, so odd-length grids and zero-padded truncated rows are
  unaffected and the median change across a table is zero; the largest
  observed shifts are 3.8% (adaptively truncated) and 12.3% (untruncated), on
  rows that end on live, curving values. Against an independent transcription
  of `simps.m`, the new default measures 3.78e-15 median relative error where
  the old rule measured 3.66e-07. This changes reported `l1_linf` ranges: pass
  `simpson_close = "trapezoid"` to reproduce 0.0.20 numbers exactly. The
  filter, smoother and ACI metric are untouched, and their machine-precision
  gradings against the MATLAB oracles are unchanged.
* `plot.aci_fit()` gains a `truth` argument, and its `"latent"` view now sets a
  vertical range covering the reconstruction, its 2 sd band and any supplied
  truth, so none of the three can be clipped off the panel.
* `print.nilcheck()` reports the queried direction, the direct sensitivity, the
  empirical peak and both flags, in place of the raw list.
* The low-signal masking warning from `forward_cir()` names its source in
  author-year form rather than by internal citation key.

## Documentation

* Per-function `@references` with full citations and DOIs replace the pointer
  lines on every paper-bearing help topic; `?aci_references` remains the key
  legend. `gaspari_cohn()` now cites Gaspari and Cohn (1999);
  `causation_entropy()` and `learn_model()` cite the causation entropy
  lineage (Sun and Bollt 2014; Elinger and Rogers 2021; Chen and Zhang 2023);
  `osse_twin()` promotes its informal mention of Arnold and Dey (1986) to a
  full entry. `?gaussian_kl` and `?observed_trajectory` cite the published
  basis for the strict public KL and the noise-free observation scope.
* The vignettes ship as a pair: `vignette-1-intro` (the front door) and
  `vignette-2-advanced` (the engine reference, retitled "The ACI engine: from
  models to extremes"). Every demonstrated claim in the engine vignette was
  checked against the source; the ensemble comparison starts both smoothers
  from the same prior, the conditional ACI example uses the model's recorded
  reference partition, and event detection uses the moser2026extremes dyad
  convention.
* Help-page corrections: the `osse_twin()` title no longer says "fitted", since
  the function takes any `cgns_model`; `ensemble_lag_table()` no longer states
  that no reference MATLAB exists; `detect_events()` documents its returned
  columns, including that `duration` is the half-peak width rather than the run
  length; `aci()`, `enkbf()` and `ensemble_lag_table()` document the list form
  of `ic_sampler`; `?nil_surrogate_test` records that the test is not out of
  sample when the model was learned from the same series; `forward_cir()`
  records that its subjective length follows eq. G.7 while the reference script
  counts cells, a one-step offset.
* README citation keys in prose replaced with author-year citations;
  `URL`/`BugReports` added to `DESCRIPTION`.

# aci 0.0.20

GitHub release.

## Features

* Closed-form conditional-Gaussian (CGNS) filter and backward-ODE smoother,
  with an explicit and a positivity-preserving implicit Riccati stepper,
  sub-stepping, and an accumulated one-step predictive log-likelihood.
* The assimilative causal-information metric `aci()`, oriented as smoother
  relative to filter, with its signal/dispersion decomposition.
* Each range estimator now follows its own reference codebase's active
  quadrature. The forward estimators and the backward exact form integrate
  with composite Simpson over the whole span, following the ACI reference
  code; the backward ratio uses the plain Appendix G L1 sum, following the
  FBCIR code's active line. (The two reference codebases differ on this
  point: FBCIR retains Simpson only as a commented alternative.) The forward
  `l1_linf` ratio now agrees with the golden R port to machine precision
  (measured 3.9e-15 median relative error, previously ~3%). That port shared
  this version's trapezoid close on even-length grids, so the agreement is
  port parity rather than parity with `simps.m` itself; see 0.0.21, which
  closes the remaining gap. The earlier right-endpoint sum over positive-lag cells was
  first order and biased low for a decaying sequence; reported ranges therefore
  change slightly. Under the paper's monotonicity conditions -- the forward
  CIR metric nonincreasing in lag (andreou2026cir eq. 12) and the complete
  backward metric nondecreasing (its Theorem 1) -- the two forward forms
  coincide exactly, the backward pair agrees to the gap between the Simpson
  and plain-sum rules, and the bound orderings hold in both directions. Each
  ratio's alternative quadrature is reachable through the `quadrature`
  argument of `forward_cir()` and `backward_cir()`.
* Forward and backward causal influence ranges (`forward_cir()`,
  `backward_cir()`, `cir_pair()`, `cir_table()`), built on the
  andreou2026cir Theorem 3 online smoother and a divergence table with
  tolerance-based row truncation (`lag_table()`).
* Conditional ACI by two strategies, prescribed forcing and innovation
  masking (`nontarget()`, `reduce_nontarget()`).
* Ensemble Kalman-Bucy filter and smoother (`enkbf()`, `enkbs()`) with
  Gaspari-Cohn localization and multiplicative inflation, and the forward
  ensemble lag table (`ensemble_lag_table()`).
* Model discovery: polynomial CGNS dictionaries, Gaussian causation entropy,
  FFBS posterior path sampling, and constrained/energy-joint fits
  (`learn_model()`).
* The partially implemented moser2026extremes extreme-event workflow: event
  detection, onsets, composites, sensitive directions, feature sets and
  clustering.
* A formula front end (`aci_fit()`) with the standard generics, and
  benchmark model constructors for the dyad, predator-prey, tipping triad,
  ENSO, Lorenz-84/96, topographic and FBCIR multiscale systems.

## Scientific provenance

* The published sources each component implements, and the shorthand keys used
  to cite them throughout the help pages, are listed in `?aci_references`. Read
  it before citing any claim of parity with the published reference
  implementations.
* The ensemble, model discovery and extreme-event layers are equation-checked
  against their source papers (recursions, parameter sets, feature definitions
  and estimand equations). The ensemble engine and forward ensemble range are
  additionally graded at machine precision against a literal R port of the
  published EnKBS dyad experiment (<https://github.com/jiangzh67/EnKBS>, MIT),
  with both sides driven by identical increments; replaying published MATLAB
  runs themselves is not claimed, since the random streams differ. The
  extreme-event paper publishes no code, so that layer remains
  equation-checked only.

## Known limitations

* Observations must lie on a uniform time grid and are assumed effectively
  noise-free; irregular or noisy records must be resampled first with
  `as_uniform_trajectory()`.
* `forward_cir()` and `backward_cir()` default to `method = "exact"`, the
  threshold-averaged last-exit functional. This is not the same quantity as
  the efficient integral reported in the source papers, which is
  `method = "l1_linf"`; the two differ whenever the lagged discrepancy is not
  monotone. Each result records which was used in its `bound` field.
* `lt_tail_bound()` returns a heuristic tail estimate, not a certified error
  bound.
* Ensemble backward CIR is deliberately absent, following the source paper's
  identification of it as future work.
