# aciR (development version)

## New features

* The fixed-lag online smoother and the causal influence range now run on
  **vector-valued states**. `aci_online_smoother()` and `aci_cir()` dispatch on
  the shape of the components exactly as the core does.

  One thing does not carry over from the scalar implementation, and it is the
  thing that made the scalar one fast. In one dimension the ordered product of
  the per-step auxiliary matrices reduces to a difference of cumulative
  logarithms, so any range is recoverable in constant time. Matrices do not
  commute and there is no such reduction. What survives is the reason the
  reduction was worth having: the products decay geometrically, so the
  accumulation is truncated once its norm falls below tolerance rather than
  reconstructed from endpoint summaries.

  The auxiliary matrices are implemented in full generality, equations (3.5) to
  (3.7) of the source paper. The reference implementation carries only the
  zero-cross-noise specialisation, equation (3.8), which is what its own models
  need; the general form is here because the components schema exposes the
  cross-covariance and would otherwise offer a path the recursion could not
  take.

* `aci_enso_model()`, `aci_enso_components()` and `aci_enso_parameters()` build
  the stochastic ENSO model of the paper's case study: three observed
  variables -- the central- and eastern-Pacific temperature anomalies and the
  Walker circulation -- and three unobserved ones, the zonal current, the
  thermocline depth and the intraseasonal wind burst.

  It is the largest system the package expresses and the only one whose noise
  covariances vary in time. The Walker circulation's observation noise is
  multiplicative in its own state, and the latent noise of the wind burst depends on
  the observed central-Pacific temperature and on the season. A path whose
  Walker circulation reaches either end of its domain is rejected rather than
  regularised: the noise vanishes there, and the filter inverts it.

  `aci_conditional()` composes with it directly, which is the configuration the
  case study is built around.

* `aci_simulate()` gains the **Milstein scheme** for systems whose diffusion
  varies with the state it integrates, alongside the Euler-Maruyama scheme it
  has always used. Pass `scheme = "milstein"` with `sigma_x` and `d_sigma_x`,
  the diffusion and its derivative as functions of the observed state.

  The correction term is not estimated numerically on the caller's behalf: a
  derivative supplied by finite differences would silently degrade the
  convergence order the scheme exists to provide, so it must be given.

  The two schemes coincide exactly, not approximately, when the diffusion is
  constant. When it is not, the tests measure what distinguishes them --
  Euler-Maruyama converges strongly at order one half and Milstein at order
  one, against the closed-form solution of geometric Brownian motion driven by
  the same Wiener path the simulator integrated.

  A caveat worth stating plainly: this is a **simulation** capability. The
  filter still requires a constant observation-noise covariance in the scalar
  schema, so a scalar path generated with a state-dependent diffusion cannot
  yet be assimilated. The vector schema has no such restriction.


* The numerical core takes **vector-valued states**. `aci_filter()`,
  `aci_smoother()` and `aci_metric()` dispatch on the shape of the components
  they are given: a matrix latent-noise covariance selects the vector path, a
  bare number the scalar one. The observed signal may then be a matrix with one
  row per observed component, the coefficients matrices or arrays with time in
  the last margin.

  The scalar path is unchanged and is not routed through the new code. At one
  dimension the two agree bit for bit, which is graded, but the scalar
  recursion is some thirty times faster and is what the package's oldest oracle
  grades; replacing it would spend a validated asset to buy nothing.

  In the vector case the covariance is re-symmetrised at each step, since an
  explicit Euler step breaks symmetry at round-off and the asymmetry compounds;
  positive-definiteness is checked by attempting a Cholesky factorisation
  rather than by a proxy for it; and the relative entropy is evaluated through
  Cholesky factors throughout, never an explicit inverse or a determinant
  formed as a product.

* `aci_conditional()` asks the causal question of one observed component given
  that the others are also observed. The construction is the reference
  implementation's: inflating a component's observational uncertainty without
  bound sends its weight in the filter to zero, which is implemented by
  supplying an inverse noise Grammian supported only on the target block. That
  object is deliberately not the inverse of any covariance, which is why the
  components schema now admits `S_xoS_x_inv` directly.

  This changes the estimand rather than the arithmetic. The result answers a
  different question, and the assumptions article says so.

* The unobserved component's self-drift may now vary in time. `L_y` accepts one
  value per observation in a components list, and a function of the observed
  signal in `aci_cgns_model()`, alongside the constant it accepted before. A
  constant still means what it meant, so components lists and models built
  against the earlier contract are unaffected.

  This admits systems whose latent damping is set by the observed state, which
  is still a conditional Gaussian system -- the coefficient is measurable with
  respect to the observed path. What remains out of scope is a self-drift
  depending on the unobserved component itself, which would leave the
  conditional Gaussian class altogether.

  `model$L_y` is now a coefficient function, like `model$L_x` and `model$f_y`
  before it, rather than a bare number. Code reading it as a number should read
  `model$L_y_constant`, which holds the value when the self-drift is constant
  and `NA` when it is not.

* `aci_predprey_model()` and `aci_predprey_components()` build the noisy
  predator-prey model, a stochastic Lotka-Volterra pair, in either causal
  direction: `"prey_to_predator"` observes the predator and treats the prey as
  latent, and `"predator_to_prey"` is the converse. The two are different
  questions rather than one question and its mirror, and the metric is not
  symmetric between them.

  This is the system that needed a self-drift varying in time: in both
  directions the latent population's growth rate is set by the population being
  watched. It is graded against the reference implementation in both
  directions, over a self-drift that changes sign, so the latent process is
  damped at some times and driven at others.

* `aci_cir()` computes the causal influence range, the second quantity of the
  method paper. Where the causal-information metric measures how much the
  future of the observed signal says about the unobserved state at a given
  time, the range measures how far into that future one must look before the
  answer stops changing. Both the subjective range, at a grid of thresholds,
  and the objective range, the threshold-free summary, are returned, along
  with the peak divergence each is derived from.

  A time close to the end of the record cannot be resolved, because the
  observations that would settle it do not exist. Such a time yields a small
  number that looks like a confident short range, so the result marks it
  saturated and returns `NA` rather than the truncated value. Resolution is
  judged separately for the objective range and for each threshold of the
  subjective one, since the demanding thresholds run far longer than the
  objective range does.

  The reference implementation stores the whole matrix of comparisons, one row
  per time and one column per later observation, which is quadratic in the
  window and reaches several gigabytes at the scale of the published figure.
  Nothing downstream needs the matrix, so `aci_cir()` forms one row at a time
  and reduces it immediately, leaving the memory linear in the window.

* `aci_online_smoother()` runs the fixed-lag forward-in-time online smoother of
  Andreou, Chen and Li (2026) <doi:10.48550/arXiv.2411.05870>. Where
  `aci_smoother()` conditions every step on the whole observed path, the online
  smoother conditions step `j` on the path up to `j + lag`, so it is the
  estimator available while the signal is still arriving. It is the object
  `aci_cir()` is built from.

  The two boundaries of the lag family are the estimators the package already
  had: at `lag = 0` the result is the filter, exactly, and as `lag` grows the
  result approaches the backward smoother. The first of those is asserted as a
  bitwise identity. The second is a first-order limit rather than an identity,
  because the online recursion is the discrete smoother while `aci_smoother()`
  integrates the continuous backward equation; the test asserts the measured
  convergence order, which a wrongly stated recursion would not produce.

  Each new observation updates earlier steps through an ordered product that
  decays geometrically, a property the source paper establishes by bounding the
  spectral radius of each factor. The implementation accumulates those products
  in logarithms and truncates them once they fall below tolerance, so the work
  is bounded by the effective lag and the memory is linear in the signal
  length. The reference algorithm forms the full quadratic triangle instead;
  a literal transcription of it is carried in the tests and the two agree to
  machine precision.

## Minor improvements and fixes

* A time-varying observation-noise covariance is now honoured rather than read
  at its first step. The vector components schema accepted a three-dimensional
  `S_xoS_x` and inverted only its first slice, so a covariance that changed
  over the record produced a result identical, to the last bit, to a constant
  one. Positive-definiteness is now checked at every step rather than only the
  first, and `aci_conditional()` likewise builds a weight that tracks the
  covariance it is derived from.

  Nothing caught this: not the oracles, whose fixtures all held that covariance
  constant; not the contracts, which accepted the array; not line coverage,
  which was total. The package's grading register now tests the general
  property instead -- a coefficient is genuinely consumed only if perturbing it
  at a late step changes the answer.

* The clamp diagnostic in `summary()` no longer counts the terminal step. The
  metric at the final step is exactly zero because the smoother is the filter
  there by construction, but the old count (`sum(aci == 0)`) treated every
  zero as sitting at the round-off floor, so it reported at least one clamped
  step on every run, including runs in which nothing was clamped. `aci()` now
  records the number of values actually clamped -- negative by no more than
  round-off before being set to zero -- in a new `n_clamped` field of the
  result, and `summary()` reports that count. An `aci` object saved by 0.1.0
  has no such field, and `summary()` now says so rather than miscounting;
  re-run `aci()` to refresh it.

* The dyad vignette draws the observed signal and the causal-information
  metric in two panels rather than on one shared axis. The two series are in
  different units -- the signal in its own units, the metric in nats -- and a
  shared axis invited magnitude comparisons across them that carry no meaning.

* The dyad vignette says which of the paper's two quantities it covers and
  which it does not, and the validation vignette displays the scalar Gaussian
  relative-entropy formula the metric evaluates. The scope note was written
  when the causal influence range was unimplemented; it now distinguishes the
  two questions instead, since the range arrived later in this same development
  cycle.

* `DESCRIPTION` states the release's scope in one sentence, and `Authors@R`,
  `LICENSE` and `LICENSE.md` carry the upstream copyright of the reference
  MATLAB implementation (Marios Andreou, MIT) alongside the package's own.
  The attribution was always in the README, the citation entries and the
  oracle manifest; it now also travels with the licence, where a derived
  work's notice belongs.

# aciR 0.1.0

First research preview. The validated dyad calculation is unchanged and still
agrees with the MATLAB oracle to well inside `1e-6`; the work in this release
is to turn the method's implicit mathematical preconditions into enforced
software contracts, and to state what the package's evidence does and does not
cover.

## Breaking changes

* `aci_simulate_dyad()` is removed. It was a second implementation of the same
  Euler-Maruyama integration as `aci_simulate()`, with weaker validation and no
  `seed` argument, and two copies of one set of equations are two places for
  them to drift apart. Use the model object instead.

  ``` r
  # Before
  set.seed(1)
  sim <- aci_simulate_dyad(n = 2000, p = list(d_x = 0.5, d_y = 0.5, gamma = 2,
                                              F_x = 0.5, F_y = 1,
                                              sigma_x = 0.5, sigma_y = 1))
  comp <- aci_dyad_components(sim$x, p)

  # After
  model <- aci_dyad_model()
  sim <- aci_simulate(model, n = 2000, seed = 1)
  comp <- aci_dyad_components(sim$x, model$parameters)
  ```

  `aci_simulate()` returns a data frame rather than a list, and `aci_simulate_dyad()`'s
  parameter list is now reachable as `aci_dyad_model()$parameters`.

* `aci_cgns_model()` no longer takes an `S_xoS_y` argument. The noise
  cross-covariance of a scalar system is its own transpose, so the value is now
  derived from `S_yoS_x` and the two can no longer be set to disagree. The
  `S_xoS_y` entry of the components schema is unchanged, and a hand-built
  components list in which the two disagree is now rejected.

  ``` r
  # Before -- the two could silently disagree
  aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1,
                 S_yoS_y = 1, S_yoS_x = 0.5, S_xoS_y = 0.2)

  # After
  aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1,
                 S_yoS_y = 1, S_yoS_x = 0.5)
  ```

* Inputs that were previously accepted and produced `NA`, `NaN` or a low-level
  replacement error are now rejected with a message naming the offending
  argument and index. This is a breaking change for code that relied on the old
  behaviour, though such code was reading missing values as results. The cases
  are: an inadmissible noise covariance, a coefficient function returning the
  wrong type or length, a non-finite or incomplete observed signal, a malformed
  components list, and a malformed posterior.

* `aci_simulate()` draws its random increments before the recursion rather than
  two at a time inside it. Paths are still reproducible from a `seed`, but the
  mapping from seed to path has changed: a seed that produced a given path in
  0.0.0.9000 produces a different one here. This does not affect the oracle
  grade, which runs on the MATLAB signal rather than on a simulated path.

* `aci_simulate(seed = )` no longer leaves the global random-number state
  changed. The generator state is saved before the draw and restored on exit, so
  a seeded call no longer consumes the caller's stream.

## New features

* An independent-oracle fixture for the noise-cross-covariance path
  (`cross_signal_x.csv`, `cross_reference.csv`) and the harness that generates
  it (`oracle/aci_oracle_cross.m`). Every scalar model in the reference
  implementation sets the noise cross-covariance to zero, so the dyad fixture
  executes the terms that carry it but pins them at zero and never grades them,
  while `aci_cgns_model(S_yoS_x = )` exposes those terms publicly. The new
  fixture grades them over the full transient, and `test-identities.R` grades
  them exactly in the stationary regime against the algebraic Kalman-Bucy fixed
  points, which are independent of both implementations.

* `summary()`, `plot()` and `as.data.frame()` methods for `aci` objects. The
  summary reports the metric distribution and its peak alongside three
  diagnostics -- the smallest posterior covariances, the terminal-identity
  residual, and how many metric values sat at the round-off floor -- so the
  quantities a result is sensitive to are surfaced rather than left to be
  discovered.

* `aci()` accepts an observed `time` grid instead of a `dt`, validates that it
  is strictly increasing and equally spaced, and derives the step from it. An
  irregular grid is rejected: the recursions integrate a fixed step and have no
  contract for one.

* `inst/extdata/oracle-manifest.yml` records, for each fixture, the upstream
  repository and commit, the licence, the MATLAB release, the generating
  command, the parameters and seed, both hashes, the measured error -- and,
  for each, what it grades and what it does not. `test-oracle-manifest.R`
  checks the shipped bytes against the recorded hashes.

* Two vignettes. *Assumptions and interpretation* states the estimand, the
  conditions under which the method is valid, and what an ACI peak does and does
  not support. *Validation and the independent oracle* records the oracle
  design, provenance and observed agreement.

* `API_STABILITY.md` declares a lifecycle stage for every export.

## Minor improvements and fixes

* The filter and smoother check the posterior covariance at every step and stop
  at the first one that is not finite and strictly positive, naming the
  algorithm, the index, the time and the value. An explicit Euler scheme can
  drive the covariance out of its domain when `dt` is too large for the system
  even though the model is admissible, and the relative entropy that scores the
  result has no meaning in that state. Nothing is silently clipped.

* `aci_metric()` evaluates the dispersion term in a cancellation-resistant form.
  Writing the covariance ratio as `1 + d`, the previous expression subtracted
  two nearly equal quantities and lost the value exactly where the metric is
  smallest: at a ratio of `1 + 1e-12` it returned `0` where the correct value is
  `2.5e-25`. Values at the round-off floor are clamped to zero within a
  documented tolerance, and anything more negative is an error rather than a
  result.

* `aci_cgns_model()` rejects a negative latent-noise covariance, an asymmetric
  noise cross-covariance and a joint covariance that is not positive
  semidefinite. Previously `S_yoS_y = -1` constructed successfully and produced
  `NaN` on simulation. A singular joint covariance remains admissible: it
  describes perfectly correlated noise, and the per-step guards own the runtime
  consequences.

* The filter and smoother loop bounds no longer descend when given a single
  observation, which is what turned `aci_filter(1, ...)` into an opaque
  replacement error rather than a boundary check.

* `aci_dyad_components()` validates its parameter list, and `aci_dyad_model()`
  rejects a parameterisation whose initial state is undefined.

* External functions are qualified at the call site rather than imported.

* `DESCRIPTION` gains `URL`, `BugReports`, a `cph` role and a tested R version
  floor; `Language` is `en-AU`.

* A committed `.lintr` admits the paper's equation notation (`L_x`, `S_xoS_x`,
  `A_j` and relatives) by a narrow rule rather than a blanket exclusion, so the
  correspondence to the published equations survives a style pass.

# aciR 0.0.0.9000

* Initial development scaffold.
* Adds the conditional Gaussian nonlinear system (CGNS) core: a forward
  filter (`aci_filter()`), a backward smoother (`aci_smoother()`) and the
  relative-entropy causal-information metric (`aci_metric()`), each taking a
  general CGNS components list so they are not tied to a single model.
* Adds `aci_dyad_components()` and `aci_simulate_dyad()` for the nonlinear
  dyad model with intermittent extreme events used as the worked example.
* Adds a model layer over the core. `aci_cgns_model()` is the general
  constructor for a CGNS `aci_model` object, `aci_dyad_model()` builds the
  flagship nonlinear dyad model, `aci_simulate()` draws an Euler-Maruyama
  realisation of a model, and `aci()` is the single user-facing entry point
  that runs the filter, smoother and metric and returns an `aci` object.
  `print()` methods are provided for both classes.

## Roadmap

* A noisy predator-prey (Lotka-Volterra) model constructor. Its unobserved
  process has a time-varying self-drift, whereas the current numerical core
  integrates a time-invariant self-drift; the constructor waits on a core that
  admits a time-varying self-drift and on its own independent-oracle fixture.
* The causal-influence-range (`aci_cir()`): the subjective and objective
  influence-range lengths from the fixed-lag online smoother. This is an
  order-N-squared computation and waits on its own independent-oracle fixture.
