# acir 0.2.2

## Bug fixes

* `aci(table = )` stops with `aci_error_dims` when `stepper`, `nsub` or
  `regularize` is also supplied. The table already fixes these settings, and
  the call used to ignore them; give them to `lag_table()`. `?lag_table` now
  says that it takes only `stepper = "explicit"` and `nsub = 1`.
* An `init` list with an entry other than `mean` and `cov` stops with
  `aci_error_dims` naming the unknown entries. `init = list(m = 5, cov = 1)`
  used to run with prior mean 0. A longer name that R matched to `cov`, as
  in `covariance =`, used to set the prior covariance and is refused as
  well. The `source` entry of a stored `meta$init` is accepted, so a stored
  prior can be passed back.
* `?aci_model` now says what assigning into an existing model does. An
  assigned drift coefficient (`Lx`, `fx`, `Ly` or `fy`) changes the results
  of `aci_filter()`, `aci_smoother()`, `aci_online()`, `lag_table()` and
  `aci()` but not those of `simulate()`. The help used to say that both stay
  on the original coefficients. No behaviour changed.
* Plot methods for ACI results, Gaussian paths and ranges pass extra
  arguments to `plot()` as 0.2.0 did, including unnamed arguments,
  `panel.first` and `panel.last`, and named labels such as `main` still
  replace the defaults. Plots of long series are fast again.
* The covariance error raised in the smoother's backward pass no longer
  advises `stepper = "implicit"`, which applies to the forward filter only.
  The condition class is unchanged.
* A range with one anchor, drawn as a title only, no longer stops on
  graphical settings such as `lwd = 0`. Plotting a range with more than one
  anchor whose values are all `NA` stops with an `aci_error_dims` condition
  that asks for `ylim`.
* `aci_range()` on a table built with a finite `max_lag` labels a row
  `"censored"` when the row was cut at `max_lag` and its last kept divergence
  is above the strength floor. Such rows were labelled `"resolved"`. Ranges
  and every other value are unchanged, and so is the status of rows the cut
  did not reach.
* `aci()` and `aci_metric()` no longer return `Inf` for some KL signals just
  below the overflow limit, for models and paths with one hidden component.
  Results for ordinary inputs are unchanged.
* `aci_online()` at a finite lag no longer returns a covariance that is not
  positive definite. Under the default policy it stops with
  `aci_error_covariance_not_spd` naming the anchor's index and time, as
  `lag = Inf` and `lag_table()` do; under `regularize = "floor"` it floors
  the covariance and records the event in `meta$regularization`. Results
  that were positive definite are unchanged.
* `aci_online()`, `lag_table()` and the ranges built on them warn with
  `aci_warn_leading_order` when the hidden-noise step ratio exceeds 1:
  one step of hidden noise exceeds the filter variance at the start of an
  interval. A prior variance below one step of hidden noise does this at the
  first interval, and a filter variance that falls that low does it later in
  the record. The smoothed variance could then exceed the filter variance
  with no condition, 0.826 against 0.001 on the packaged dyad at
  `dt = 0.01`. Returned values are unchanged, and the condition carries the
  ratio, index and time. Intervals from the first floored filter covariance
  on are left out, since `aci_warn_regularized` reports the floor.
  `?aci_online` explains this diagnostic and its limits. The update also
  depends on drift, and a ratio at or below 1 does not establish accuracy.
* `aci()` accepts a lag table that `lag_table()` built for a model that
  declares its own estimand, such as `aci_enso_model(hidden = "tau")`.
* `lag_table()`, `aci_online()` and `aci_smoother()` given a filter that
  `regularize = "floor"` floored, and `aci_range()` on a floored `aci()`
  result that kept its paths, now record the filter's floor events and
  policy `"floor"` in `meta$regularization` and raise one
  `aci_warn_regularized`. Before, they started an empty record, so a table
  built on a floored filter, and `aci(table = )` on it, reported no floor.
  The record survives `saveRDS()`. Results are unchanged.

## Documentation

* `?aci_filter` has a new section on `nsub` and `stepper`. A larger `nsub`
  refines the filter integration on the fixed observation grid. It can also
  improve the smoother and ACI, but does not in general remove the error from
  holding filter moments fixed within each interval. `stepper` applies to the
  forward filter only. `?aci`, `?aci_smoother` and `?aci_online` point to the
  section, and the advanced vignette says the same of its sub-stepped
  examples.
* The plot help says which arguments replace the defaults, and when a range
  plot needs `ylim`.

# acir 0.2.1

## Bug fixes

* Scalar `aci_metric()` dispersion retains small positive variance ratios and
  avoids intermediate ratio overflow when the dispersion is representable.
  The scalar signal calculation also recovers finite results when its initial
  squared standardised mean difference overflows. Ordinary finite signal
  calculations retain their previous arithmetic.
* Plot methods for ACI results, Gaussian paths and causal influence ranges
  honour caller-supplied labels without passing duplicate arguments to the
  plotting function. A single-anchor range plot retains its default title.

* **`aci_range()` reports explicit single-step provenance on every retention
  route.** `aci()` records
  `stepper` and `nsub` in `meta`, recording the filter's integration settings
  (`"explicit"` and `1L` for a result built on a reused table). The backward
  smoother inherits `nsub` and remains explicit. When a
  result computed with `stepper = "implicit"` or `nsub > 1` is passed to
  `aci_range()` with no retained paths (`keep = "none"`), the range is, as
  before, computed from a recomputed explicit single-step filter and smoother,
  since the lag-table recursions are defined for that discretization only;
  the rebuild is now reported by one `aci_warn_stepper` warning naming the
  recorded settings, as the `keep = "paths"` route already did off its
  retained filter. With `keep = "table"`, `aci()` warns while building the
  table and every subsequent `aci_range()` call warns once about using that
  explicit single-step table, including after `saveRDS()`/`readRDS()`.
  Older retained-table results use the settings on their retained filter.
  Range values are unchanged on every route. A result saved by an
  earlier release with `keep = "none"` carries no record of its settings, so
  it cannot report the rebuild; `?aci_range` states the behaviour.

## Tests and documentation

* The range status descriptions distinguish length and strength checks from
  possible sensitivity to the end of the record. Lag-storage controls and
  threshold coverage are clarified. The online-smoother example counts
  anchors whose windows reach the end of the record.
* Scalar KL regression tests cover small variance ratios and intermediate
  dispersion and signal overflow.

* The retention and serialization checks now include a correlated
  matrix-observation model. The range help states the common inputs and
  settings needed for comparing routes and allows floating-point differences.
* The conditional reduction/masking comparison requires all five complete,
  finite series with matching names and lengths. Its maximum absolute
  difference limit remains `1e-12`. CI records the per-series differences
  and runtime details. The advanced vignette no longer promises bit-for-bit
  agreement across platforms.

# acir 0.2.0

## Release

* The corrections below are released as 0.2.0. Citation metadata now agrees
  with the package version and its R 4.1.0 requirement.
* The deprecated simulation argument `T` remains accepted with a warning in
  this release. Use `t_end`; no removal version is scheduled here. Supplying
  both arguments remains an error.

## Breaking changes

* **The declared minimum R version is 4.1.0.** It was 4.0.0, and the only
  floor the check matrix has ever built is 4.1. No source change was needed:
  the package uses no syntax newer than R 4.0, so the raise replaces an
  untested claim with one that is built on every run. Installation on R 4.0.x
  is now refused by R.

## Behaviour changes

* **A degenerate observation-noise Gram is refused wherever it occurs, not
  only at the constructor's probe points.** `aci_model()` has always required
  `rcond(gxx) >= 1e-12` at the five points it probes. The same test now runs
  on the whole realised Gram at every interval start of the record, whatever
  the route and whether or not the realisation cache is in use, before any
  precision is formed. A model whose observation Gram becomes exactly zero
  between probe points previously completed with
  a silently jittered precision and no recorded regularisation event; it now
  raises `aci_error_gram_path`, a subclass of `aci_error_gram`, naming the
  grid index, the time and the offending `rcond`. The check runs before the
  covariance policy, so `regularize = "floor"` does not bypass it. Exported
  `safe_chol()` keeps its documented jitter ladder unchanged; what changed is
  that an invalid observation model no longer reaches it. No shipped library
  model is refused on its test or packaged reference records: the suite pins a
  realised `rcond` above 1e-4 on the dyad record and on the two ENSO partition
  records it builds, and the smallest value measured across the packaged
  reference records is 7.7e-05, more than seven orders of magnitude above the
  1e-12 threshold.

* **Overflow in a filter, smoother or likelihood recursion is an error, not a
  result.** A run whose mean or predictive log-likelihood becomes `Inf` or
  `NaN` now raises `aci_error_nonfinite` naming the quantity, the grid index
  and the time, instead of returning a Gaussian path with non-finite means,
  finite positive covariances and no recorded events. Finite results are
  unchanged, and the new class is deliberately distinct from
  `aci_error_covariance_not_spd`: this is invalid arithmetic, not a covariance
  the policy can floor.

* **Regularized results say so.** When `regularize = "floor"` actually takes a
  floor, the call now raises one `aci_warn_regularized` naming the first
  floored site, its grid index and its time, and `print()` on an `aci_result`,
  `lag_table`, `cir_result` or assimilation path adds a line reporting the
  event count. Previously a floored run could complete with no condition of
  any kind: a floor taken only in the backward smoother raised nothing,
  because the stiffness diagnostic exists only in the explicit filter. Every
  number is unchanged, including the eigenvalue floor, and a run that takes no
  floor prints exactly what it printed before. This is a compatibility-relevant
  signal change: code running under `options(warn = 2)` on a record whose floor
  fires now stops where it previously returned. The condition is classed, so
  `suppressWarnings()` or a `withCallingHandlers()` muffle on
  `aci_warn_regularized` restores the previous flow. The help for `regularize`
  now
  states that flooring changes the numerical covariance so the recursion can
  continue and establishes nothing about the accuracy of the reconstruction or
  of the resulting information score; `aci_metric()` carries no record of its
  own and reports whatever its input paths were computed under.

* **Covariance failures no longer offer remedies that cannot apply.** When a
  covariance leaves the positive-definite cone, or a variance becomes
  non-finite on the scalar path, the message now says that the observation
  Gram was accepted against the constructor's relative-conditioning contract
  at every step of this record, so the failure is integration instability on
  a valid observation model, and the positive-definite message adds that
  flooring changes the numerical covariance without establishing that the
  step resolves the dynamics. The `regularize` help states what that contract
  is and is not: the Gram test is a conditioning test on `rcond()`, which is
  invariant under uniform rescaling, and not a noise-floor test. For a
  one-dimensional observed state, small positive Grams can pass; extreme
  floating-point scales can also cause condition estimation to fail.
  The condition classes,
  their `site`, `role`, `index`, `time` and `value` fields and the
  `regularize = "floor"` instruction are unchanged.

* **`nsub` above the integer range is refused by the package, not by a
  coercion.** `nsub = 2^40` previously produced R's own "NAs introduced by
  coercion" warning followed by "missing value where TRUE/FALSE needed", on
  both the scalar and the matrix route; it now raises `aci_error_dims` with
  the usual message. Every accepted value behaves exactly as before.

* **Mismatched shared noise-channel widths report the contract they violate.**
  `Sx1`/`Sy1` or `Sx2`/`Sy2` pairs with different column counts previously
  failed inside a matrix product with "non-conformable arguments"; the widths
  are now checked before that product and the existing
  `aci_error_model_contract` message is what the caller sees. Admissible
  rectangular shared channels are unaffected.

* **The diffuse-prior warning no longer offers a burn-in window as a repair.**
  It now says that the opening steps are prior-dominated and that a prior far
  wider than the hidden state's own scale can destabilise the explicit step
  outright, which is a refusal rather than a window to discard. The class
  `aci_warn_diffuse_init`, the condition it is raised on and every number are
  unchanged.

* **The stepper warning says what the recursions are derived for.** The
  `aci_warn_stepper` text raised by `lag_table()` and by the forward-CIR
  reducer said the Theorem 3 recursions are "exact for" the explicit
  single-step discretization; it now says they are "derived for" it. The
  condition class and the recomputation it announces are unchanged.

* **The observation contract's refusals no longer oversell their remedy.** The
  non-uniform-grid message no longer stops at "resample first": it says that
  interpolation onto a uniform grid does not create independent observations
  and does not remove observation error, and the noise-free refusal adds that
  the model's process diffusion is not sensor variance and cannot stand in for
  it. `observed_trajectory()`'s help now opens with what this release
  assimilates: a complete record of the observed state on a uniform grid,
  taken as noise-free, with every observed channel present at every time. The
  classes and the inputs accepted or refused are unchanged.

* Five condition messages no longer cite internal specification sections or
  development notes: the singular observation-noise Gram, the non-uniform
  time grid, the noise-free observation contract, and the inadmissible
  `method = "reduce"` reduction.

* **`citation("acir")` reports the installed version.** `inst/CITATION` reads
  it from the package metadata that `utils::readCitationFile()` already
  supplies, instead of carrying a transcribed number that went stale at the
  first version bump.

## Documentation

* **Small observation noise is not guaranteed to cause a later failure.**
  `aci_filter()` help now states that implicit integration can return finite
  positive covariances without a warning or regularization event, even when
  the time step does not resolve the covariance dynamics. Explicit integration
  need not fail either; its behaviour depends on the coupling and time step.
  Regression tests cover these returns, valid large reductions in uncertainty
  and changes of units. No calculation, default or acceptance threshold changes;
  a further scale diagnostic remains deferred.

* **A worked ENSO startup, and what a floored covariance is worth.**
  `aci_enso_model()` gains a fully specified worked example for the
  three-state partition `c("u", "hW", "tau")`, and a help section explaining
  why the automatically supplied scalar prior is refused there at the first
  update: it is one scalar for every hidden component, and on this partition
  the component scales differ by three orders of magnitude. The example states
  its prior componentwise as the Ornstein-Uhlenbeck stationary scale read off
  the model's own coefficients and uses the default explicit stepper at
  `nsub = 1`; refinement and prior sensitivity are shown as checks on that
  configuration rather than as alternatives to it. The bare constructor's own
  default partition, `hidden = c("hW", "tau")`, is now stated in the help.
  *The closed-form ACI engine* carries the same partition end to end, and adds
  the contrast that motivates the example: the same record under
  `regularize = "floor"` returns a peak of 6.0e+21 with two floor events
  against 1.83 with none, so a floored number is a diagnostic that the prior
  or the step did not resolve the dynamics, not a reconstruction and not an
  information score.

* **ACI is described as conditional on the supplied model.** The introduction,
  `aci()` and the causal-metrics topic now state that the metric scores
  information gained under the model, prior and observed record supplied to
  it, give the units (nats per time point, nats times model time when
  integrated), distinguish a supplied likelihood from parameter estimation,
  and note that the prior and the integration step move the value as well as
  the model does. On a structurally independent Brownian null with prior
  variance 1, one observed interval of length `dt`, the observed path
  `x(t) = 2t`, the explicit stepper at `nsub = 1` and `regularize = "none"`,
  the backward smoother returns 1.0135e-4 nats at `dt = 0.1` and 2.4419e-8 at
  `dt = 0.0125`; the same null at prior variance 0.1 and `dt = 0.07` returns
  1.4660115 nats, which `nsub = 20` reduces to 0.0140779. All of those runs
  are positive definite and finite and record zero floor events, which does
  not certify that the step resolved the problem. A positive value alone is
  not an empirically identified causal effect, an intervention effect or a
  significance statement.

* **The online smoother's finite-step output is described as an
  approximation.** `aci_online()` composes the published leading-order
  Theorem 3 updates; at a finite step that composition approximates the
  continuous-time conditional law and is not generally the exact posterior of
  an Euler-sampled record. The help of `aci_online()` and `assimilation_api`
  now says so with a worked comparison, `aci_online(lag = Inf)` is described
  as composing those updates over the whole record rather than returning the
  exact posterior, and the vignettes and the package-level comments follow.
  The `lag = 0` boundary is still exact, both scheme labels
  (`"theorem3_discrete"`, `"backward_ode_euler"`) are unchanged, and every
  measured figure is unchanged.

* **The truncation tail quantity is labelled a heuristic estimate
  everywhere.** The shipped value is a norm-based suffix accumulation with a
  1.5 multiplier, not the spectral-radius condition of andreou2026smoother
  eq. 3.19, and it is a diagnostic under the retained record rather than a
  guarantee about the cells the truncation dropped. The evidence register and
  the source comments now agree with the help, and the historical
  specification's attribution is annotated rather than removed. No arithmetic
  changed, including the 1.5 multiplier.

* **The two influence-range functionals are documented with their
  quadrature.** `method = "exact"` is a finite threshold sum with no time-axis
  quadrature, and is exact for its discrete counting functional rather than
  for continuous-time inference; the `l1_linf` ratio is integrated with
  composite Simpson by default, with `quadrature = "sum"` selecting the L1
  grid-function sum. The help said the two functionals coincide wherever the
  divergence decreases with lag. They do so term by term, which makes them
  equal only under `quadrature = "sum"`; under the shipped Simpson default
  they still differ. On the row `c(1, 0.5, 0)` with `M = 1` the objective is
  `1.5 * dt`, the summed ratio `1.5 * dt` and the Simpson ratio `1 * dt`, and
  `tests/testthat/test-04-cir.R` pins all three at tolerance 0. No estimator
  changed.

* **The relative-entropy help records where its arithmetic loses precision.**
  `aci_metric()` uses a cancellation-resistant `log1p` form for a
  one-dimensional hidden state and the trace and log-difference form for two
  or more; `aci_metric_pair()` uses the trace form in every dimension, so the
  two public routes need not agree in the last digits on identical
  one-dimensional inputs. The help now gives the measured relative errors as
  the two covariances converge and records that the absolute error stays below
  3e-17 nats across that range. No arithmetic changed.

* **Models and observed trajectories are documented as immutable.**
  `aci_model()`'s help now says that a model is fixed once constructed, that
  assigning a new coefficient into an existing model leaves simulation and
  filtering on the original coefficients, and shows a parameter being changed
  by rebuilding. The realisation cache's help names its probe points as the
  first, middle and last grid point and states plainly that a captured change
  visible only away from all three is not detected, with
  `options(aci.realiser_cache = FALSE)` as the switch for a model of that
  kind; the evidence register's cache row is worded to match. The 0.1.0
  performance entry below says the probe "catches a parameter the closures
  read from their environment moving between calls"; it carries the earlier,
  broader wording and is left as the historical record. `as_obs()` returns an
  existing `obs_traj` without re-validating it, so an object whose fields were
  assigned into after construction keeps whatever `dt`, `k` or `noise_free` it
  was left with.

* **The `seed` guarantee is stated precisely.** `simulate()`'s help now says
  that an existing `.Random.seed` is restored bit for bit, that when none
  exists one is created before the seeded draw so that there is a state to
  restore and the caller is left holding it, and that an unseeded call
  advances the caller's stream. No seeded path changed.

* **The parity claim is narrowed to what the register records.** `README.md`
  no longer says that every graded quantity reproduces the method authors'
  MATLAB reference; it and `inst/evidence/register.csv` are what an installed
  package carries on this question. Of the evidence register's
  68 rows, 28 are graded against a hash-pinned fixture and 40 have no fixture
  behind them and record an exact relation, a behavioural check or an in-test
  independent transcription. Ten of the fixture-backed rows cite a fixture the
  method authors produced: nine compare a computed output against it, and the
  tenth records that `observed_trajectory()` ingests the authors' pinned input
  signal and carries it onto the seven dyad grades; the two predator-prey
  grades run on the authors' pinned predator-prey record, which has no
  register row of its own. The remaining fixture-backed rows compare against
  a source-derived run, an independent transcription of the published
  equations, or an analytic identity. The 0.1.0 entry below
  carries the earlier, broader wording and is left as the historical record;
  `inst/evidence/register.csv` is the authority for what each feature is
  graded against.

* **The evidence register's `aci` dyad row claims only the total ACI metric.**
  `dyad_reference.csv` carries no signal or dispersion column, and the
  comparison that runs is of the total alone; the decomposition's graded
  evidence is the T_C zeroth-order rows and the `aci_metric_pair` row, all
  independent transcriptions. Three rows are added for the contracts
  introduced here: the realised observation-noise Gram, the non-finite
  recursion guards, and the regularization-visibility signal.

* **The partition full-record summary fixture is described for what it
  bounds.** It bounds magnitude over the whole record and fixes position only
  at the sampled indices. Its four reductions are order-insensitive, so a
  rearrangement confined to unsampled steps is invisible to them; the earlier
  wording said the summary closed the sampling gap. No fixture, sample index,
  reduction or tolerance changed.

* The vignettes show the package's refusals as the classed condition and its
  message, printed as ordinary output, instead of a halted chunk. A genuine
  error in those chunks now fails the vignette build.

# acir 0.1.0

The parity release: every quantity the reference MATLAB implementation
computes is reproduced to the tolerance recorded in
`inst/evidence/register.csv`, and every stage of the specification's
performance table is inside its budget.

## Breaking changes

* **The simulation horizon is `t_end`, not `T`.** `simulate()` and
  `aci_simulate()` take `t_end`; the name `T` shadowed the alias of `TRUE`
  in R and is retired. A call that still passes `T` works, with a warning of class
  `aci_warning_deprecated`, until acir 0.2.0; passing both is an error.


## Bug fixes

* **The built package carries its vignette index.** The build-ignore rule
  `^build($|/)`, inherited from aci, removed the `build/` directory that
  `R CMD build` creates, so an installed package listed no vignettes under
  `vignette()` or `browseVignettes()` and the CRAN incoming check reported a
  missing vignette index. The rule is gone.

## Performance

* **A model supplied as closures is realised once per record, not once per
  verb.** The generic route evaluated every coefficient closure at every grid
  point on each call, so a filter, a smoother, an online smoother and an
  influence range on one model and one record realised the same arrays four
  times. The last four realisations are now kept, keyed on the model object
  and the grid, and a stored one is reused after its coefficients are
  re-evaluated at three grid points and found unchanged, which catches a
  parameter the closures read from their environment moving between calls.
  The arrays are identical either way, so no number moves. The option
  `aci.realiser_cache = FALSE` bypasses the cache (`R/aci-realiser-cache.R`,
  the specification's fingerprint, plan v0.3 PR-7).

* For a scalar hidden state the Theorem 3 auxiliaries (the update factor,
  the gain row and the one-lag increments of every interval) are formed as
  vector expressions over the record (`.online_aux_scalar()`,
  `R/aci-online-scalar.R`) instead of one interval at a time, and feed the
  Theorem 3 smoother, the forward primitives of the range, the online window
  route and the lag table. On the dyad record the smoother falls from 0.17 s
  to 0.004 s and the primitives from 0.15 s to 0.02 s, within one rounding
  of the per-interval kernels (bit-identical where the BLAS's triangular
  solve divides); the covariance policy is reached by the per-interval route
  whenever a variance needs it. Ported from aciR 0.2.3 (`.aci_online_aux()`
  in `R/aci-online-smoother.R`, tag `parents-final`). The matrix path is
  unchanged.

* The forward influence range and the lag table form each anchor's row of
  divergences as one vector expression when the hidden state is scalar
  (`.cir_scalar_row()`, `R/aci-cir-rows.R`), instead of advancing every
  anchor one cell at a time. On the authors' dyad record (N = 3,000, 1,001
  reporting anchors) the range falls from 44 s to under a second and the
  adaptive lag table from 137 s to one second, with every reported quantity
  within 1e-13 of the cell-by-cell recursion and the freeze indices, tail
  estimates and statuses unchanged. Ported from aciR 0.2.3 (`.aci_cir_row()`
  in `R/aci-cir.R`, the cumulative logarithms of `.aci_online_aux()` in
  `R/aci-online-smoother.R`; tag `parents-final`), with one change of
  summation: the logarithms of the update factors are summed from the anchor
  outward in blocks of 512 cells rather than differenced from a record-length
  cumulative sum, which errs by at least the rounding of that sum on any
  platform (6e-13 at 100,000 steps of contracting factors on x86) and by
  1.7e-12 at 20,000 steps and 2.8e-11 at 100,000 where `cumsum` accumulates
  in double (arm64); the blocked form stays below 3e-13. The
  relative entropy of each cell is evaluated in the operations of
  `.kl_fast()`, so a cell is bit-identical to the recursion's given the same
  posterior. The matrix path (`l > 1`) is unchanged.
