# Changelog

## aciR 0.2.3

### Testing

- Coverage rises from 81% to 97%, against the 95% floor the continuous
  integration enforces. The gate had never run before this release,
  because nothing had been pushed with the workflows active, and the
  first run failed.

  The largest gap was the whole family of reporting methods for
  [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md),
  at zero. [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) are all
  exported and documented, and none had a test. They now do, including
  the branches that decide whether a censored range is announced as a
  bound.

  The remainder were the guards that reject a malformed model: a missing
  parameter, an unrecognised causal direction, an observed process with
  no noise for the filter to invert, a target channel with no
  observation noise of its own. Each exists because the alternative is
  not an error but a plausible-looking number, and each is now exercised
  through the public entry point rather than by calling the internal
  checker.

### Documentation

- The reference for the fixed-lag online smoother, which
  [`aci_online_smoother()`](https://biometryhub.github.io/ACI/reference/aci_online_smoother.md)
  and
  [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  both cite, now points at the published article rather than at the
  preprint. It appeared in *Journal of Nonlinear Science* **36**(4), 71,
  [doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x),
  after the preprint was cited here.

- The project website is at <https://biometryhub.github.io/ACI/>. It
  carries the function reference and the three articles, together with
  two ledgers published nowhere else: one recording how this package’s
  numbers were graded against the authors’ own code, quantity by
  quantity, and one recording the design decisions, the alternatives set
  aside, and what remains open.

### Internal

- `DESCRIPTION` lists the website alongside the repository in `URL`,
  which is what lets the documentation link back to itself.

- `CITATION.cff` and `codemeta.json` are regenerated from `DESCRIPTION`.
  Both had been left at an earlier version, and the citation file
  disagreed with itself about which one.

## aciR 0.2.2

### New features

- [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  gains `increments`, which takes the standard normal variates to
  integrate instead of drawing any. It defaults to `NULL`, and every
  existing call keeps drawing its own exactly as before.

  This closes the last ungraded piece of the numerical surface. The
  simulator was the one component with no independent oracle. The reason
  on record was that R and MATLAB draw normal variates by different
  algorithms, so a simulated path cannot coincide with the reference
  path even at a matching seed. That statement is about the two
  generators rather than about the integrator, which is the component a
  grade should reach. Euler-Maruyama is invertible. Subtract the drift
  from a captured transition, divide by the noise coefficient, and out
  comes the variate that produced it. The variates the reference
  actually drew are therefore recoverable from its own captured path,
  and driving the simulator with them reproduces that path rather than
  another draw from the same law.

  Run that way against the reference’s 3001-step dyad capture, the
  integrator reproduces the observed component to **5.1e-15** and the
  unobserved to **1.8e-15**, against a round-off bound of **1.15e-10**
  derived from the arithmetic before the comparison was made. That bound
  is five ulps of the largest state visited, summed over the steps and
  weighted by the amplification each injection still has to pass
  through.

  The change carries two costs. The simulator now has two ways to obtain
  its increments, so the seeded path’s behaviour rests on an asserted
  identity (supplied variates equal to the seeded draws reproduce a
  seeded call bit for bit) rather than on there being only one route
  through the function. The grade also needs both components of the
  reference path on the integration grid, which no packaged oracle
  fixture carried. `inst/extdata/dyad_signal_x.csv` holds the observed
  path alone and `inst/extdata/dyad_reference.csv` holds the pair only
  at every hundredth step, so a 3001-step capture of both joins the test
  fixtures.

  Supplying `increments` together with `seed` is an error. A seed
  governs nothing on a path that draws no random numbers, and an
  argument that did nothing would misdescribe the call.

## aciR 0.2.1

### Minor improvements and fixes

- The vignettes typeset their mathematics as MathML, converted when the
  page is built, rather than fetching MathJax when it is opened.

  Two things were wrong with the previous arrangement. MathJax matches
  the x-height of its own fonts to the surrounding text, and the
  vignette body is a sans-serif face with a large x-height, so every
  displayed equation was scaled up by about a third and sat out of
  proportion with the prose around it. Less visibly, the mathematics was
  fetched from a remote host at read time. With no network the equations
  did not render at all, and the reader saw the TeX source. MathML is
  set by the browser at the surrounding size and needs nothing external.
  `pkgdown` already defaults to it, so the two surfaces now agree.

- `inst/CITATION` derives the package version from the installed
  metadata rather than repeating it as a literal, which is one fewer
  place for a release to leave a stale number behind.

## aciR 0.2.0

### Breaking changes

- [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  closes an odd interval count with the quadratic through the last three
  samples, where it previously used a Simpson 3/8 panel over the last
  three intervals.

  This changes `objective`. The quadratic closure is the rule the
  reference implementation uses, and adopting it takes agreement with
  the published numbers from **4.58e-09 to 1.37e-14** on a truncated
  comparison horizon, taken as the maximum over the full reported region
  of both graded datasets. The 3/8 panel is the slightly more accurate
  rule on an equally spaced grid, measurably so, and it was chosen for
  that reason. We judged fidelity to the reference, on the quantity the
  method leads with, to be the more valuable of the two properties.
  Whether to offer the more accurate closure as an option is an open
  question.

  The change also removes a latent defect. The 3/8 panel assumed equal
  spacing, and `objective_exact` integrates over a logarithmic threshold
  grid. Neither the 129-point default nor the reference’s 513-point grid
  triggered it, since both have an even interval count, but any
  even-length `epsilon` would have. The replacement is exact for unequal
  spacing.

- A reported time whose range is not resolved now returns a
  right-censored lower bound rather than `NA`. A new `status` field
  carries the qualification (`"resolved"`, `"censored"`,
  `"below_threshold"`, `"insufficient"`). `saturated` is retained and
  equals `status == "censored"`.

  Such a time is not unmeasured. The record supports the statement that
  the range is at least this long, and returning `NA` discarded that
  statement at the end of the record, which is where a user studying a
  recent event looks. The one status that still returns `NA` is
  `"insufficient"`, where fewer than three later observations leave no
  quadrature to evaluate.

- [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  returns the logical matrix `subjective_censored`, marking the
  individual thresholds whose range ran past the retained margin.

### New features

- A synthetic fast-contracting fixture exercises the online smoother’s
  truncation branch for the first time. On every system this package
  ships the update product never reaches `tol` inside a record it can
  hold, so the branch was dead on all real fixtures; the new fixture
  reaches it in 121 steps of 400 and checks that truncating leaves the
  answer where an untruncated walk leaves it.

  The fixture also measured the filter’s transient. The source paper
  bounds the spectral radius of each per-step factor below one, and that
  bound describes the settled filter. Here the first two factors have
  radius 1.200 and 1.100, because the covariance starts away from its
  fixed point, before settling at 0.700. A lag bound estimated from a
  contraction rate measured over all steps would be misled by that
  transient. The implementation tests the accumulated product rather
  than the individual factors, which is the form the transient cannot
  mislead.

- A regression guard, `tests/testthat/test-retired-claims.R`, fails the
  build if a retired claim reappears anywhere in `R/`, `man/`, the
  vignettes, `DESCRIPTION`, `README.Rmd`, `NEWS.md` or the oracle
  manifest. Two review rounds each found documentation asserting
  something that had ceased to be true, corrected on the surfaces the
  author had touched and left standing on the others. The check costs
  one search and was not being run, so it is a test rather than a habit.

- [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  gains [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
  returned object was previously a classed list with no behaviour, so
  printing it dumped a threshold-by-time matrix over the three scalars a
  reader needs.

- The README carries a recipe for reproducing the published
  causal-influence-range panels, which needs `horizon`, `margin` stood
  down, and the reference’s 513-point threshold grid.

- [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  gains a `horizon` argument, and returns `objective_exact` alongside
  `objective`.

  Both come from grading the package against the authors’ own MATLAB
  rather than against a transcription of it. The reference truncates its
  comparison sequence at the end of its reporting window; this
  implementation runs to the end of the record, which is the more
  faithful reading of the definition but means the two do not produce
  the same numbers. `horizon` makes the reference’s convention
  available, and with it the subjective range agrees to **2.2e-16**
  where it differed by 1.16 time units at the default.

  `objective_exact` is the objective range by its definition, the
  subjective ranges integrated over the whole threshold grid, as
  distinct from the efficient underestimate `objective` already
  reported. The source paper gives both, and only the second was
  implemented. It agrees with the reference’s `defn_objective_CIR` to
  **1.6e-14** once the horizons match.

  The default is unchanged. `horizon = NULL` uses the whole record, so
  existing results are not affected.

- A reported time whose comparison sequence is too short to support a
  range, which `horizon` can produce near the end of a window, is now
  marked saturated and returns `NA` rather than being integrated over
  two points. The reference records a zero there, which reads as “no
  detectable influence” when it means “not measured”.

- The fixed-lag online smoother and the causal influence range now run
  on vector-valued states.
  [`aci_online_smoother()`](https://biometryhub.github.io/ACI/reference/aci_online_smoother.md)
  and
  [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  dispatch on the shape of the components exactly as the core does.

  One property does not carry over from the scalar implementation, and
  it is the property that made the scalar one fast. In one dimension the
  ordered product of the per-step auxiliary matrices reduces to a
  difference of cumulative logarithms, so any range is recoverable in
  constant time. Matrices do not commute, so there is no such reduction,
  and no matrix analogue of it is known. What survives is the property
  the reduction exploited. The products decay geometrically, so the
  accumulation is truncated once its norm falls below tolerance rather
  than reconstructed from endpoint summaries.

  The auxiliary matrices are implemented in full generality, equations
  (3.5) to (3.7) of the source paper. The reference implementation
  carries only the zero-cross-noise specialisation, equation (3.8),
  which is what its own models need; the general form is here because
  the components schema exposes the cross-covariance and would otherwise
  offer a path the recursion could not take.

- [`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md),
  [`aci_enso_components()`](https://biometryhub.github.io/ACI/reference/aci_enso_components.md)
  and
  [`aci_enso_parameters()`](https://biometryhub.github.io/ACI/reference/aci_enso_parameters.md)
  build the stochastic ENSO model of the paper’s case study. It carries
  three observed variables (the central- and eastern-Pacific temperature
  anomalies and the Walker circulation) and three unobserved ones (the
  zonal current, the thermocline depth and the intraseasonal wind
  burst).

  It is the largest system the package expresses and the only one whose
  noise covariances vary in time. The Walker circulation’s observation
  noise is multiplicative in its own state, and the latent noise of the
  wind burst depends on the observed central-Pacific temperature and on
  the season. A path whose Walker circulation reaches either end of its
  domain is rejected rather than regularised, because the noise vanishes
  there and the filter inverts it.

  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
  composes with it directly, which is the configuration the case study
  is built around.

- [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  gains the Milstein scheme for systems whose diffusion varies with the
  state it integrates, alongside the Euler-Maruyama scheme it has always
  used. Pass `scheme = "milstein"` with `sigma_x` and `d_sigma_x`, the
  diffusion and its derivative as functions of the observed state.

  The correction term is not estimated numerically on the caller’s
  behalf. A derivative supplied by finite differences would silently
  degrade the convergence order the scheme exists to provide, so it must
  be given.

  The two schemes coincide exactly, not approximately, when the
  diffusion is constant. When it is not, the tests measure what
  distinguishes them, namely that Euler-Maruyama converges strongly at
  order one half and Milstein at order one, against the closed-form
  solution of geometric Brownian motion driven by the same Wiener path
  the simulator integrated.

  The capability is confined to simulation. The filter still requires a
  constant observation-noise covariance in the scalar schema, so a
  scalar path generated with a state-dependent diffusion cannot yet be
  assimilated. The vector schema carries no such restriction, and is the
  route available today. The scalar restriction would lift with a scalar
  schema that admits a time-varying observation-noise covariance, as the
  vector schema already does.

- The numerical core takes vector-valued states.
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
  [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
  and
  [`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
  dispatch on the shape of the components they are given. A matrix
  latent-noise covariance selects the vector path, a bare number the
  scalar one. The observed signal may then be a matrix with one row per
  observed component, the coefficients matrices or arrays with time in
  the last margin.

  The scalar path is unchanged and is not routed through the new code.
  At one dimension the two agree bit for bit, which is graded, but the
  scalar recursion is some thirty times faster and is what the package’s
  oldest oracle grades. Keeping both leaves two implementations to
  maintain. Retiring the scalar one would cost that speed and would move
  the oldest oracle onto code it has never graded, which we judged the
  larger loss.

  In the vector case the covariance is re-symmetrised at each step,
  since an explicit Euler step breaks symmetry at round-off and the
  asymmetry compounds; positive-definiteness is checked by attempting a
  Cholesky factorisation rather than by a proxy for it; and the relative
  entropy is evaluated through Cholesky factors throughout, never an
  explicit inverse or a determinant formed as a product.

- [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
  asks the causal question of one observed component given that the
  others are also observed. The construction follows the reference
  implementation. Inflating a component’s observational uncertainty
  without bound sends its weight in the filter to zero, which is
  implemented by supplying an inverse noise Grammian supported only on
  the target block. That object is deliberately not the inverse of any
  covariance, which is why the components schema now admits
  `S_xoS_x_inv` directly.

  This changes the estimand rather than the arithmetic. The result
  answers a different question, and the assumptions article says so.

- The unobserved component’s self-drift may now vary in time. `L_y`
  accepts one value per observation in a components list, and a function
  of the observed signal in
  [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md),
  alongside the constant it accepted before. A constant still means what
  it meant, so components lists and models built against the earlier
  contract are unaffected.

  This admits systems whose latent damping is set by the observed state.
  Such a system is still a conditional Gaussian system, since the
  coefficient is measurable with respect to the observed path. A
  self-drift depending on the unobserved component itself remains out of
  scope. It would leave the conditional Gaussian class, and with it the
  closed-form posterior this package is built on, so no route to
  admitting it is available here.

  `model$L_y` is now a coefficient function, like `model$L_x` and
  `model$f_y` before it, rather than a bare number. Code reading it as a
  number should read `model$L_y_constant`, which holds the value when
  the self-drift is constant and `NA` when it is not.

- [`aci_predprey_model()`](https://biometryhub.github.io/ACI/reference/aci_predprey_model.md)
  and
  [`aci_predprey_components()`](https://biometryhub.github.io/ACI/reference/aci_predprey_components.md)
  build the noisy predator-prey model, a stochastic Lotka-Volterra pair,
  in either causal direction. `"prey_to_predator"` observes the predator
  and treats the prey as latent, and `"predator_to_prey"` is the
  converse. The two are different questions rather than one question and
  its mirror, and the metric is not symmetric between them.

  This is the system that needed a self-drift varying in time. In both
  directions the latent population’s growth rate is set by the
  population being watched. It is graded against the reference
  implementation in both directions, over a self-drift that changes
  sign, so the latent process is damped at some times and driven at
  others.

- [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  computes the causal influence range, the second quantity of the method
  paper. Where the causal-information metric measures how much the
  future of the observed signal says about the unobserved state at a
  given time, the range measures how far into that future one must look
  before the answer stops changing. Both the subjective range, at a grid
  of thresholds, and the objective range, the threshold-free summary,
  are returned, along with the peak divergence each is derived from.

  A time close to the end of the record cannot be resolved, because the
  observations that would settle it do not exist. Such a time yields a
  small number that looks like a confident short range, so the result
  marks it censored and returns the value as a lower bound rather than
  as a resolved range. Only a time with fewer than three later
  observations, where no quadrature exists, returns `NA`. Resolution is
  judged separately for the objective range and for each threshold of
  the subjective one, since the demanding thresholds run far longer than
  the objective range does.

  The reference implementation stores the whole matrix of comparisons,
  one row per time and one column per later observation, which is
  quadratic in the window and reaches several gigabytes at the scale of
  the published figure. Nothing downstream needs the matrix, so
  [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  forms one row at a time and reduces it immediately, leaving the memory
  linear in the window.

- [`aci_online_smoother()`](https://biometryhub.github.io/ACI/reference/aci_online_smoother.md)
  runs the fixed-lag forward-in-time online smoother of Andreou, Chen
  and Li (2026),
  [doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x).
  Where
  [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
  conditions every step on the whole observed path, the online smoother
  conditions step `j` on the path up to `j + lag`, so it is the
  estimator available while the signal is still arriving. It is the
  object
  [`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)
  is built from.

  The two boundaries of the lag family are the estimators the package
  already had. At `lag = 0` the result is the filter, exactly, and as
  `lag` grows the result approaches the backward smoother. The first of
  those is asserted as a bitwise identity. The second is a first-order
  limit rather than an identity, because the online recursion is the
  discrete smoother while
  [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
  integrates the continuous backward equation; the test asserts the
  measured convergence order, which a wrongly stated recursion would not
  produce.

  Each new observation updates earlier steps through an ordered product
  that decays geometrically, a property the source paper establishes by
  bounding the spectral radius of each factor. The implementation
  accumulates those products in logarithms and truncates them once they
  fall below tolerance, so the work is bounded by the effective lag and
  the memory is linear in the signal length. The reference algorithm
  forms the full quadratic triangle instead; a literal transcription of
  it is carried in the tests and the two agree to machine precision.

### Minor improvements and fixes

- A time-varying observation-noise covariance is now honoured rather
  than read at its first step. The vector components schema accepted a
  three-dimensional `S_xoS_x` and inverted only its first slice, so a
  covariance that changed over the record produced a result identical,
  to the last bit, to a constant one. Positive-definiteness is now
  checked at every step rather than only the first, and
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
  likewise builds a weight that tracks the covariance it is derived
  from.

  Nothing caught this. The oracles did not, since their fixtures all
  held that covariance constant. The contracts did not, since they
  accepted the array. Line coverage did not, since it was total. The
  package’s grading register now tests the general property instead,
  that a coefficient is genuinely consumed only if perturbing it at a
  late step changes the answer.

- The clamp diagnostic in
  [`summary()`](https://rdrr.io/r/base/summary.html) no longer counts
  the terminal step. The metric at the final step is exactly zero
  because the smoother is the filter there by construction, but the old
  count (`sum(aci == 0)`) treated every zero as sitting at the round-off
  floor, so it reported at least one clamped step on every run,
  including runs in which nothing was clamped.
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) now
  records in a new `n_clamped` field of the result the number of values
  actually clamped, meaning those negative by no more than round-off
  before being set to zero, and
  [`summary()`](https://rdrr.io/r/base/summary.html) reports that count.
  An `aci` object saved by 0.1.0 has no such field, and
  [`summary()`](https://rdrr.io/r/base/summary.html) now says so rather
  than miscounting; re-run
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) to
  refresh it.

- The dyad vignette draws the observed signal and the causal-information
  metric in two panels rather than on one shared axis. The two series
  are in different units (the signal in its own units, the metric in
  nats), and a shared axis invited magnitude comparisons across them
  that carry no meaning.

- The dyad vignette says which of the paper’s two quantities it covers
  and which it does not, and the validation vignette displays the scalar
  Gaussian relative-entropy formula the metric evaluates. The scope note
  was written when the causal influence range was unimplemented; it now
  distinguishes the two questions instead, since the range arrived later
  in this same development cycle.

- `DESCRIPTION` states the release’s scope in one sentence, and
  `Authors@R`, `LICENSE` and `LICENSE.md` carry the upstream copyright
  of the reference MATLAB implementation (Marios Andreou, MIT) alongside
  the package’s own. The attribution was always in the README, the
  citation entries and the oracle manifest; it now also travels with the
  licence, where a derived work’s notice belongs.

## aciR 0.1.0

First research preview. The validated dyad calculation is unchanged and
still agrees with the MATLAB oracle to well inside `1e-6`; the work in
this release is to turn the method’s implicit mathematical preconditions
into enforced software contracts, and to state what the package’s
evidence does and does not cover.

### Breaking changes

- `aci_simulate_dyad()` is removed. It was a second implementation of
  the same Euler-Maruyama integration as
  [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md),
  with weaker validation and no `seed` argument, and two copies of one
  set of equations are two places for them to drift apart. Use the model
  object instead.

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

  [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  returns a data frame rather than a list, and `aci_simulate_dyad()`’s
  parameter list is now reachable as `aci_dyad_model()$parameters`.

- [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md)
  no longer takes an `S_xoS_y` argument. The noise cross-covariance of a
  scalar system is its own transpose, so the value is now derived from
  `S_yoS_x` and the two can no longer be set to disagree. The `S_xoS_y`
  entry of the components schema is unchanged, and a hand-built
  components list in which the two disagree is now rejected.

  ``` r

  # Before (the two could silently disagree)
  aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1,
                 S_yoS_y = 1, S_yoS_x = 0.5, S_xoS_y = 0.2)

  # After
  aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1,
                 S_yoS_y = 1, S_yoS_x = 0.5)
  ```

- Inputs that were previously accepted and produced `NA`, `NaN` or a
  low-level replacement error are now rejected with a message naming the
  offending argument and index. This is a breaking change for code that
  relied on the old behaviour, though such code was reading missing
  values as results. The cases are an inadmissible noise covariance, a
  coefficient function returning the wrong type or length, a non-finite
  or incomplete observed signal, a malformed components list, and a
  malformed posterior.

- [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  draws its random increments before the recursion rather than two at a
  time inside it. Paths are still reproducible from a `seed`, but the
  mapping from seed to path has changed. A seed that produced a given
  path in 0.0.0.9000 produces a different one here. This does not affect
  the oracle grade, which runs on the MATLAB signal rather than on a
  simulated path.

- `aci_simulate(seed = )` no longer leaves the global random-number
  state changed. The generator state is saved before the draw and
  restored on exit, so a seeded call no longer consumes the caller’s
  stream.

### New features

- An independent-oracle fixture for the noise-cross-covariance path
  (`cross_signal_x.csv`, `cross_reference.csv`) and the harness that
  generates it (`tools/oracle/aci_oracle_cross.m`). Every scalar model
  in the reference implementation sets the noise cross-covariance to
  zero, so the dyad fixture executes the terms that carry it but pins
  them at zero and never grades them, while `aci_cgns_model(S_yoS_x = )`
  exposes those terms publicly. The new fixture grades them over the
  full transient, and `test-identities.R` grades them exactly in the
  stationary regime against the algebraic Kalman-Bucy fixed points,
  which are independent of both implementations.

- [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for `aci` objects. The summary reports the metric distribution and its
  peak alongside three diagnostics (the smallest posterior covariances,
  the terminal-identity residual, and how many metric values sat at the
  round-off floor), so the quantities a result is sensitive to are
  surfaced rather than left to be discovered.

- [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) accepts
  an observed `time` grid instead of a `dt`, validates that it is
  strictly increasing and equally spaced, and derives the step from it.
  An irregular grid is rejected, because the recursions integrate a
  fixed step and have no contract for one.

- `inst/extdata/oracle-manifest.yml` records, for each fixture, the
  upstream repository and commit, the licence, the MATLAB release, the
  generating command, the parameters and seed, both hashes, the measured
  error, and what that fixture grades and what it does not.
  `test-oracle-manifest.R` checks the shipped bytes against the recorded
  hashes.

- Two vignettes. *Assumptions and interpretation* states the estimand,
  the conditions under which the method is valid, and what an ACI peak
  does and does not support. *Validation and the independent oracle*
  records the oracle design, provenance and observed agreement.

- `API_STABILITY.md` declares a lifecycle stage for every export.

### Minor improvements and fixes

- The filter and smoother check the posterior covariance at every step
  and stop at the first one that is not finite and strictly positive,
  naming the algorithm, the index, the time and the value. An explicit
  Euler scheme can drive the covariance out of its domain when `dt` is
  too large for the system even though the model is admissible, and the
  relative entropy that scores the result has no meaning in that state.
  Nothing is silently clipped.

- [`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
  evaluates the dispersion term in a cancellation-resistant form.
  Writing the covariance ratio as `1 + d`, the previous expression
  subtracted two nearly equal quantities and lost the value exactly
  where the metric is smallest. At a ratio of `1 + 1e-12` it returned
  `0` where the correct value is `2.5e-25`. Values at the round-off
  floor are clamped to zero within a documented tolerance, and anything
  more negative is an error rather than a result.

- [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md)
  rejects a negative latent-noise covariance, an asymmetric noise
  cross-covariance and a joint covariance that is not positive
  semidefinite. Previously `S_yoS_y = -1` constructed successfully and
  produced `NaN` on simulation. A singular joint covariance remains
  admissible. It describes perfectly correlated noise, and the per-step
  guards own the runtime consequences.

- The filter and smoother loop bounds no longer descend when given a
  single observation, which is what turned `aci_filter(1, ...)` into an
  opaque replacement error rather than a boundary check.

- [`aci_dyad_components()`](https://biometryhub.github.io/ACI/reference/aci_dyad_components.md)
  validates its parameter list, and
  [`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)
  rejects a parameterisation whose initial state is undefined.

- External functions are qualified at the call site rather than
  imported.

- `DESCRIPTION` gains `URL`, `BugReports`, a `cph` role and a tested R
  version floor; `Language` is `en-AU`.

- A committed `.lintr` admits the paper’s equation notation (`L_x`,
  `S_xoS_x`, `A_j` and relatives) by a narrow rule rather than a blanket
  exclusion, so the correspondence to the published equations survives a
  style pass.

## aciR 0.0.0.9000

- Initial development scaffold.
- Adds the conditional Gaussian nonlinear system (CGNS) core, comprising
  a forward filter
  ([`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)),
  a backward smoother
  ([`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md))
  and the relative-entropy causal-information metric
  ([`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)),
  each taking a general CGNS components list so they are not tied to a
  single model.
- Adds
  [`aci_dyad_components()`](https://biometryhub.github.io/ACI/reference/aci_dyad_components.md)
  and `aci_simulate_dyad()` for the nonlinear dyad model with
  intermittent extreme events used as the worked example.
- Adds a model layer over the core.
  [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md)
  is the general constructor for a CGNS `aci_model` object,
  [`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)
  builds the flagship nonlinear dyad model,
  [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  draws an Euler-Maruyama realisation of a model, and
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) is the
  single user-facing entry point that runs the filter, smoother and
  metric and returns an `aci` object.
  [`print()`](https://rdrr.io/r/base/print.html) methods are provided
  for both classes.

### Roadmap

- A noisy predator-prey (Lotka-Volterra) model constructor. Its
  unobserved process has a time-varying self-drift, whereas the current
  numerical core integrates a time-invariant self-drift; the constructor
  waits on a core that admits a time-varying self-drift and on its own
  independent-oracle fixture.
- The causal-influence-range
  ([`aci_cir()`](https://biometryhub.github.io/ACI/reference/aci_cir.md)),
  the subjective and objective influence-range lengths from the
  fixed-lag online smoother. This is an order-N-squared computation and
  waits on its own independent-oracle fixture.
