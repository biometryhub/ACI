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
