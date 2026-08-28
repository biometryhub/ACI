# acir 0.0.0.9000

* **The evidence a feature rests on now ships with the package, and two tests
  keep it honest.** `inst/evidence/register.csv` is one row per checked feature
  of the public surface: what is checked, by which method, against what, at
  which tolerance class, and with which hash-pinned fixture behind it. Every
  exported verb appears, and a row that claims a fixture must name one that is
  in the tree and whose bytes still match the `sha256` the row carries;
  `tests/testthat/test-30-evidence-register.R` fails the build otherwise, so a
  new export cannot ship without its evidence being stated.
  `inst/evidence/gate_liveness.md` is the companion register for the package's
  gates: one row per gate, where it is enforced, the induced violation that
  proves it fails when it should, and where that run is recorded. Most gates
  were already proved live by a test that induces the violation as part of its
  own subject; `tests/testthat/test-31-gate-liveness.R` supplies the induced
  runs for the gates that had none, among them the fixture byte-pinning
  comparison, the simulator divergence guard and the `aci_range()` backward
  direction gate. No R code, no exported object and no numerical output is
  touched.

* **Interim authorship and citation metadata for the public development
  branch.** `Authors@R` now records on its face that the author roles, their
  order and the maintainer field are under joint review, and `inst/CITATION` is
  removed rather than shipped as a settled citation; `citation("acir")` falls
  back to the automatic entry built from `DESCRIPTION` until the review closes.
  The one piece of metadata that is *not* interim ships as a settled fact: a new
  `inst/COPYRIGHTS` carries the third-party notice for the reference MATLAB
  implementation `ACI_code`, reproducing its MIT licence verbatim, naming the
  commit the oracle manifest pins, and recording that `simps.m` was never ported
  into this package. A minimal `README.md` states the same interim status. No R
  code, no exported object and no numerical output is touched.

* **`nontarget = ` is now `conditional = `.** On `aci()`, `aci_filter()`,
  `aci_smoother()`, `aci_online()` and `lag_table()`. The value has been an
  `aci_conditional_spec` built by `aci_conditional()` since the surface rename,
  so `aci(m, ob, nontarget = aci_conditional(...))` read against itself.
  Nothing else about the argument changed: same position, same `NULL` default,
  same object, same numbers. The stored field moves with it, so
  `result$meta$nontarget` is now `result$meta$conditional` and a compiled
  bundle's `$nontarget` is `$conditional`. The condition classes
  `aci_error_nontarget` and `aci_error_nontarget_crossnoise`, and the model
  metadata `meta$estimand_nontarget` and `meta$nontarget_reduction`, keep their
  names: the first pair is a condition surface of its own, and the second is a
  declaration a constructor makes, not an argument a caller passes.

* **`eps = ` is now `epsilon = ` and `eps_grid = ` is now `epsilon_grid = `** on
  `aci_range()`, which is the specification's spelling. The recorded value moves
  with them, so `meta$eps_grid` is now `meta$epsilon_grid`. The **values**
  `quadrature = "matlab_eps_grid"` and `bound = "eps_grid_objective"` are
  unchanged: they name the reference MATLAB script's own `eps_grid` variable,
  which is the point of a compatibility label.

* **The `R/` source files are renamed to the `aciR` file convention.** A pure
  file-layout change: `utils.R` becomes `aci-utils.R`, `model_classes.R`
  becomes `aci-model.R`, `benchmark_models.R` becomes `aci-model-library.R`,
  `assimilation.R` becomes `aci-assimilation.R`, `causal_metrics.R` becomes
  `aci-core.R`, `compiled_scalar.R` and `compiled_matrix.R` become
  `aci-kernels-scalar.R` and `aci-kernels-matrix.R`, `compiled_conditioning.R`
  becomes `aci-conditional.R`, `compiled_lag.R` becomes
  `aci-online-smoother.R`, `compiled_cir.R` becomes `aci-cir.R`, `plots.R`
  becomes `aci-methods.R`, and `api_documentation.R` becomes
  `aci-documentation.R`. `acir-package.R` is unchanged. No code moved between
  files and no code changed: only the file-header banners, the architecture
  map in `acir-package.R` and two cross-file comments name the new paths.
  `DESCRIPTION` has no `Collate` field, so load order is unaffected.

* **The public surface is renamed to the specification's `aci_*` verbs.**
  Nothing had shipped, so the old names are gone outright: there are no
  deprecated aliases and no wrapper duplication. The complete map, with the
  only intentional behavioural change called out under `aci_conditional()`:

  | 0.0.0.9000 (before) | now |
  |---|---|
  | `cgns_model()` | `aci_model()` |
  | `cgns_from_affine()` | `aci_model_from_affine()` |
  | `conditionally_linear_model()` | `aci_linear_model()` |
  | `model_dyad()` | `aci_dyad_model()` |
  | `model_predator_prey()` | `aci_predprey_model()` |
  | `model_enso6()` | `aci_enso_model()` |
  | `da_filter()` | `aci_filter()` |
  | `da_smooth()` | `aci_smoother()` |
  | `da_online()` | `aci_online()` |
  | `gaussian_kl_path()` | `aci_metric()` |
  | `gaussian_kl()` | `aci_metric_pair()` |
  | `forward_cir()` | `aci_range()` |
  | `nontarget()` | `aci_conditional()` |
  | `reduce_nontarget()` | `aci_conditional_reduce()` |
  | class `nontarget_spec` | class `aci_conditional_spec` |
  | class `conditionally_linear_model` | class `aci_linear_model` |

  `aci()`, `lag_table()`, `lt_row()`, `lt_diag()`, `lt_tail_bound()`,
  `as_obs()`, `observed_trajectory()`, `safe_chol()`, `spd_floor()` and the
  `simulate()` method keep their names. Every argument list is unchanged apart
  from `aci_conditional()` and the new `direction` on `aci_range()`; in
  particular the `nontarget =` argument of `aci()`, `aci_filter()`,
  `aci_smoother()`, `aci_online()` and `lag_table()` still has that name, and
  the S3 class labels `cgns_model`, `stochastic_model`, `aci_result`,
  `cir_result`, `da_path_gaussian`, `obs_traj`, `lag_table` and `aci_sim` are
  unchanged.

* **`aci_conditional()` replaces `nontarget()`, and renames its arguments and
  its values.** `blocks =` becomes `given =`, `strategy =` becomes
  `method =`, `"inflate"` becomes `"mask"` and `"prescribed_forcing"` becomes
  `"reduce"`. `target =`, `first_step =`, the target/complement equivalence
  and the `"uniform"` / `"matlab"` first-slice switch are all carried over
  unchanged. **The default moves**: `nontarget()` defaulted to
  `"prescribed_forcing"`, `aci_conditional()` defaults to `"mask"`. This is
  the one behavioural change in the rename, and it is deliberate - masking is
  what the reference conditional scripts do. Two call sites that relied on the
  old default now name `method = "reduce"` explicitly, so no recorded number
  moved.

* **`aci_range()` takes a `direction`.** `direction = "forward"` is the
  default and the only direction in this release; it is exactly what
  `forward_cir()` computed. `direction = "backward"` raises
  `aci_error_not_implemented` naming the development reserve, rather than
  returning a forward answer under a backward label. The backward range is a
  `FBCIR_code-main` feature with no `ACI_code` counterpart.

* **`aci_simulate()` is a new verb** for the existing `simulate()` method, so
  that the public surface reads as one family. It is the same call with the
  same arguments and the same draws, and the S3 method stays registered on
  `stats::simulate()`, so `simulate(model, ...)` keeps working.

* `aci_sample()` (FFBS hidden-path sampling) is in the specification's
  interface but not in this release stage; it is staged, not silently absent.

* **Three vignettes.** `vignette("vignette-1-intro")`, "Assimilative causal
  inference (ACI) in R", is the introduction: what question ACI answers, the
  model as input rather than output, a simulated record, filter, smoother, the
  ACI metric, and the forward causal influence range. Its prose is carried
  over from the introductory vignette of the predecessor `aci` package, with
  the code retargeted to the `acir` API and the sections describing the
  formula interface and the partial-observation learner replaced by scope
  notes, both being outside 0.1.0.
  `vignette("vignette-2-advanced")`, "The closed-form ACI engine", is the
  machinery: writing CGNS models, the strict covariance policy and what it
  records, the steppers, the lag table's storage, accessor and quadrature
  arguments, the three named objective functionals and the two subjective
  read-outs, conditional ACI with `target` and the two masking strategies, and
  the two smoothing discretizations. Every convention argument is presented as
  a choice that changes the number, with the size of the change measured in
  the vignette.
  `vignette("vignette-3-matlab")`, "Reproducing the reference MATLAB
  codebase", is new and collects the fidelity material: what each benchmark
  constructor's provenance claims and what it does not, both
  `model_predator_prey()` partitions, the `T_C` zeroth-order approximation,
  `matlab_defect_compat`, the conditional scripts' first-slice convention, the
  two tau estimands, the fields to read a result's conventions off, and the
  grades of evidence behind the package's fidelity claims. The `T_C`,
  `first_step` and tau material moved here from the advanced vignette rather
  than being duplicated. All 60 chunks across the three vignettes execute;
  nothing is written against a specification. `DESCRIPTION` gains
  `Suggests: knitr, rmarkdown` and `VignetteBuilder: knitr`, which is what
  makes `R CMD build` build them.

* **`model_enso6()` now covers the `T_C`-hidden partition**, the fifth and last
  `ACI_code` ENSO script and the only one whose split is not exact.
  `model_enso6(hidden = "TC", approximation = "zeroth_order_c1", prescribed =)`
  builds the two-channel reduced inference model
  `ENSO_model_cond_ACI_T_C_unobs.m` actually assimilates: observed `(T_E, I)`,
  hidden `T_C`, with `u`, `h_W` and `tau` supplied as prescribed forcings and
  looked up by index on their own grid, never interpolated. The six-state
  system is not conditionally Gaussian with `T_C` hidden - the damping
  `c_1(t, T_C) T_C` is cubic in `T_C` - and the script restores conditional
  linearity by a zeroth-order Taylor expansion of `c_1` about `T_C = 0`. What
  the model returns are the moments of that approximating system, so
  `simulate()` on it is refused and names the joint constructor instead. Every
  other `TC`-containing hidden set keeps the existing exact-split error, and
  `hidden = "I"` is still refused outright.

* **`matlab_defect_compat`, and the defect it names.** The reference script's
  assimilation forcing `f_y` (`ENSO_model_cond_ACI_T_C_unobs.m:1053,:1151`)
  omits the thermocline term `gamma_C * h_W`, which the same script's simulator
  drift at `:1124` includes, which `h_W` being prescribed and observed makes
  available, and which the sibling scripts carry in the corresponding `T_C`
  coefficient rows. `acir` includes the term; `matlab_defect_compat = TRUE`
  reproduces the published script verbatim. The difference is large: on the
  script's own 14-model-year window the filter mean moves by up to `0.105` and
  ACI by `2.75`, and the time-integrated ACI roughly doubles, `13.60` against
  `27.52`. Both filter and smoother covariances are bit-identical between the
  two, because the term enters only the mean equations, and the shipped tests
  assert that as an identity rather than as a tolerance.

* Eight new source-derived fixtures grade the new partition at `1e-12` in
  `tests/testthat/test-29-tc-zeroth-order.R`: the realised coefficients, the
  filter, smoother and ACI series on the shared 4001-point path and on a
  14-model-year window of a 22001-point one, both defect arms, and the defect
  measurement itself. Their oracle is an **independent R transcription** of the
  script's filter, smoother and ACI sections, which ships in
  `tools/fixtures/make-tc-fixtures.R`; the package supplies only the driving
  paths. That is stronger than package-to-package agreement and weaker than an
  authors-source fixture - no MATLAB was executed. The conditional mask is
  **exactly inert** on this partition, masked and unmasked runs agreeing
  bitwise, because the `I` row of `Lx` and the noise cross-Gram are both zero;
  the tests assert those zeros, and `h_W` remains the only live conditional
  test in the package.

* **Fixed: `model_enso6()$meta$conditioning_obs_idx` named the wrong
  complement.** It recorded the observed channels the hidden variable does not
  reach - `{u, tau}` for `hidden = "hW"` - rather than the complement of the
  reference script's target. Masking the recorded set is exactly inert, so a
  caller who built a `nontarget()` from the field got the unconditional run
  under a conditional name. The `u` and `h_W` scripts mask down to `T_C` alone
  (`ENSO_model_cond_ACI_u_unobs.m:1205`,
  `ENSO_model_cond_ACI_h_W_unobs.m:1202`), and the field now says so:
  `causal_link` reads `(hW) -> (TC) | (u,TE,tau,I)`, the complement built from
  the field reproduces the script's arm, and a new
  `meta$estimand_provenance` cites the line it came from. The `tau` and joint
  scripts leave every masking line commented and now record an empty masked
  complement. No numerical default changed: these fields are a record, and
  assimilation still conditions only on a specification the caller supplies or
  the constructor declares.

* Always-on oracle grades for the three scalar-hidden ENSO partitions and the
  reduced three-channel tau observation, in
  `tests/testthat/test-28-partition-oracles.R`. Six new fixtures pin the
  realised coefficient arrays and the filter, smoother and ACI series of
  `hidden = "u"`, `"hW"` and `"tau"` at 201 sampled indices, plus a
  full-record reduction of every series over all 4001 steps so nothing is
  graded on the sampled points alone. The matrix branch at `k = 5, l = 1` and
  the non-target reduction at `k = 3, l = 1` had no fixture before this. These
  are **source-derived**, not authors-source: the driving path is an R
  simulation, no MATLAB output was read for any of them, and they do not grade
  the MATLAB realisation, `aciR`'s coefficient realisation, or the upstream
  scripts' first-slice conditional convention.
  `tests/testthat/fixtures/oracles/oracle-manifest-partitions.yml` is the
  authority for the scope, and `tools/fixtures/make-partition-fixtures.R`
  regenerates every shipped byte.

* The oracle fixture directory now carries two manifests, and the byte-pinning
  test covers both. The manifest reader's name pattern excluded digits, so
  every `enso6_*` entry would have been skipped silently had the two been
  merged; the pattern now admits digits, the reader is shared through
  `tests/testthat/helper-oracle-manifest.R`, and `test-19` asserts that every
  entry of every manifest is visible to it and that between them the manifests
  account for every `.csv` in the directory.

* New `da_online()`: the fixed-lag online smoother of andreou2026smoother. The
  estimate at index `j` conditions on the observed record through index
  `j + lag`, saturating at the end of the record, and the returned path carries
  `meta$lag`, a per-anchor `meta$lag_effective` and `meta$saturated`. `lag = 0`
  returns the filter moments unchanged, value for value; `lag = Inf` returns
  the complete Theorem 3 posterior. `lag` has no default. The cost is O(N)
  whatever the lag: a sliding window over the whole 2001-step dyad record costs
  8.1 ms at lag 1 and 22.4 ms at lag 2000, against 16 ms and 9674 ms for
  accumulating each anchor separately. Graded against the pinned MATLAB online
  fixtures (dyad, cross and mv, agreement to 6.7e-15) and against
  `aciR::aci_online_smoother` over nine lags (8.4e-15).

* **`aci()$meta$smoother_scheme` now names a discretization, not a route.** A
  path records which of the two smoothing discretizations produced it in
  `meta$scheme`: `"backward_ode_euler"` for `da_smooth()`, and
  `"theorem3_discrete"` for `da_online()` and the lag table's reference
  smoother. `meta$smoother_scheme` reads that, so it reports
  `"backward_ode_euler"` where it previously said `"backward_ode"`, and
  `"theorem3_discrete"` where a reused lag table previously reported
  `"thmD1_online_complete"`. The finer implementation tag is unchanged and
  still in `meta$route`. This is a label change only: every number is bitwise
  what it was. The two schemes agree to first order in the step and no further,
  and the gap grows with the length of the record: on the packaged ENSO
  partition the smoothed means differ by up to 1.9e-02 over 401 steps and
  9.6e-02 over 4001, and the ACI values by 0.10 and 0.48 against scales of 1.09
  and 2.35.

* `simulate()` now names the observed columns of the trajectory it returns,
  from the model's own `meta$vars$observed`, and `as_obs()` on a simulation
  carries them through. `nontarget()` resolves block and target names against
  those column names, so `aci(m, as_obs(sim), nontarget = nontarget(target =
  "TC"))` works directly where it previously failed with "Named blocks need
  named obs columns" and had to be routed through an explicit
  `observed_trajectory(..., names = )`. The names are labels: every numeric
  result is bitwise what it was, and a model that declares no usable names
  still gets unnamed columns.

* **Breaking.** Covariances are strict by default. A covariance that leaves the
  positive-definite cone inside a state recursion, a metric input or the
  likelihood now stops the run with a classed `aci_error_covariance_not_spd`
  condition naming the site, the grid index, the time and the offending value,
  where it was previously projected back by `spd_floor()` and the run
  continued. `aci()`, `da_filter()`, `da_smooth()` and `lag_table()` take
  `regularize = "floor"` for the previous behaviour, and
  `options(aci.regularize = "floor")` restores it session-wide. Every floor
  taken under `"floor"` is recorded in the result's `meta$regularization`,
  which is always present and reports `fired = FALSE` with a zero-row site
  table on a clean run. Nothing in the ACI_code scope is affected: across 47
  runs covering every fixture in this workspace, a floor fired 29 times in
  4 371 482 opportunities, and all 29 were inside three deliberately broken
  stress probes. The exported `spd_floor()` and `safe_chol()` are unchanged and
  remain the implementation of the opt-in. `forward_cir()` on an `aci_result`
  reads the policy off the result rather than the option, so re-analysing a
  saved result cannot silently change it. The conditioning-Gram inverses stay
  on `safe_chol()`'s jitter ladder deliberately: they invert a model input
  whose positive definiteness `validate_cgns()` already contracts.

* Validating a supplied one-dimensional hidden-state path no longer factorises
  every 1x1 covariance slice in turn. For a 1x1 matrix a Cholesky factorisation
  succeeds exactly when the value is finite and positive, so the same paths are
  accepted and the same paths are rejected, at the same index and with the same
  error; the tests check that equivalence in both directions rather than assume
  it. Paths with more than one hidden dimension are unchanged. A
  `da_smooth(force_validate = TRUE)` call on the dyad benchmark drops from 29.7
  to 0.94 ms.

* The scalar predictive log-likelihood reads its four coefficient vectors once
  instead of once per grid point, as the filter and smoother kernels already
  did. The accumulation is the same additions in the same order and every
  result is bitwise unchanged. The kernel is 2.4x faster and a `da_filter()`
  call on the dyad benchmark 1.4x.

* The two explicit matrix inverses the compiled engine forms per step are now
  taken with one `chol2inv()` on a factor it already has, rather than two
  triangular solves against an identity: the unmasked observation-precision
  path and the smoother's filtered-covariance inverse. Every ACI_code result is
  bitwise what it was, checked on the serialized bytes rather than with
  `identical()`, which cannot see the sign of a zero. Above one hidden
  dimension the smoother's moments move by at most 1.4e-17 and the ACI, its
  decomposition and both matrix KL contracts by at most 1.8e-15; on a Gram that
  is not diagonal the precision path moves by at most 3.6e-15 relative, out to
  condition number 1e6. The Gram path is 1.7x faster, an ENSO compile 1.6x, the
  warm matrix smoother 1.2x, and a complete `aci()` call on the ENSO benchmark
  1.1x.

* `model_enso6()` now realises its coefficient path over the whole observation
  grid at once, for every partition and both conditioning strategies, instead
  of calling its coefficient closures at each of the 4001 grid points. Nothing
  about the model changes: the arguments, the documentation, the metadata, the
  declared tau estimand and every realised array are what they were, checked
  with `identical()` and never a tolerance, against both the previous realiser
  and the generic one. The compile drops 17x and the realisation step itself
  212x, measured with both routes interleaved in one process; a complete
  `aci()` call on the ENSO benchmark drops 2.2-2.4x. The route is selected from
  the constructor's sealed descriptor, so a model whose coefficient functions
  have been replaced still realises the slow, generic way. `model_enso6()`'s
  captured parameters are now locked: retuning them by reaching into the
  closure environment, which was possible and silent, is an error, and
  `model_enso6(params = )` is the supported route.

* The observation-precision path no longer repeats `chol_solve()`'s and
  `masked_ginv()`'s entry validation at every step of the grid. The
  factorisation and the two triangular solves are unchanged to the bit, an
  invalid target is rejected once with the same error, and a Gram that does not
  factor sends the whole path back through the guarded form, jitter ladder
  included. It is 2.1x faster on the plain branch and 2.4x on the masked one,
  which after the change above is most of what an ENSO compile costs.

* `nontarget()` gains `target`, naming the target observed channels directly
  instead of their complement. The reference scripts write the causal question
  that way - `h_W(t) -> T_C | (u, T_E, tau, I)` - so `target = "TC"` now reads
  with the source where `blocks = c("u","TE","tau","I")` read against it.
  `blocks` is unchanged and either side may be given, never both. The
  documentation states the estimand the two sides split.

* `nontarget()` gains `first_step`, the first-slice convention of the masked
  observation precision. `"uniform"` is the shipped behaviour and stays the
  default; `"matlab"` leaves the first slice as the full Gram inverse, which is
  what the reference scripts compute before their target-only overwrite. It is
  not a round-off-level choice: on a 4001-point ENSO path with `h_W` hidden and
  `T_C` the target it moves the step-2 filter mean by 0.108 - a 75% move on a
  filter mean of 0.144 - the peak ACI by 0.574, and the time-integrated ACI by
  -1.39%, decaying about a decade per decade of steps without reaching
  round-off inside the record. It is exactly inert on the `u`-hidden case,
  where the mask itself is inert.

* `model_enso6(hidden = "tau")` gains `observations`, and its default output
  changes. `"reduced"`, the new default for that partition only, is the
  reference script's estimand: the observed process is `(T_C, T_E, I)` and `u`
  and `h_W` enter as prescribed known time series. `"full"` keeps the
  five-channel construction and is bit-identical to the previous default.
  The script asserts the two agree; on a 4001-point path they do not. Filter
  means differ by up to 0.247, the ACI series by up to 0.776 - about three
  times its own mean level, Pearson 0.905 - while the time-averaged ACI moves
  only 0.51%. They are two estimands, and `meta$observations` records which one
  a model carries. Observations are still supplied on all five channels either
  way.

* `forward_cir()` gains an `anchors` argument selecting which anchor times to
  report, and the streamed engine now forms only those rows and skips the
  per-interval primitives before the earliest of them. Fifteen anchors on a
  2001-point record drop from 41.6 s to 719 ms, and the public route -
  `forward_cir()` on an `aci()` result that kept no lag table - reaches the
  same 713 ms while reproducing the pinned MATLAB peak to 1.03e-14 and
  objective to 2.47e-15. `lag_table()` itself still builds every row, so the
  window pays off on the streamed route rather than on a retained table.

* `method = "exact"` now computes the definitional objective range exactly.
  The objective is the subjective range averaged over every threshold, and
  read off the running maximum that average is a finite sum, `dt * sum(suffix
  max) / M`, with no quadrature error - and it is cheaper than the composite
  Simpson it replaces. Reported ranges move up by close to half a step: at
  most 6.60e-04 anywhere on a 2001-point dyad record, 5.22e-04 at the pinned
  fixture anchors. `method = "l1_linf"` is untouched and bit-identical.

* The reference script's `defn_objective_CIR` is retained under its own name,
  `quadrature = "matlab_eps_grid"`, which integrates the subjective read-out
  over the threshold grid `eps_grid` - by default the reference 513-point grid
  `10^seq(-6, 0.5, length.out = 513)`. It reproduces a verbatim port of that
  script to 3.4e-14. Quadrature nodes are only ever taken from `eps_grid` and
  reporting thresholds only ever from `eps`; supplying `eps_grid` to any other
  mode is an error, so one vector can never be asked to do both jobs.

* The subjective range gains `convention`, defaulting to `"count"`: the
  reference script's `index * dt`. Nonzero cells sit exactly one grid step
  above the previous read-out and cells with no exceedance stay zero, so the
  pinned MATLAB columns now compare directly with no normalisation.
  `convention = "lag_time"` is the andreou2026cir eq. G.7 lag time and
  reproduces the previous values to the bit.

* `cir_result` gains a per-anchor `status` factor: `"resolved"`, `"censored"`,
  `"below_threshold"` or `"insufficient"`. It is metadata - no reported number
  changes. A range is censored when the record left after the last exceedance
  of the strength floor is shorter than the exceedance itself; the test is
  taken from the row, with no margin or horizon asked of the caller.

* `da_smooth()` no longer re-validates a filter path that `da_filter()` has
  just produced for the same run. Such a path is sealed against its own
  moments, and when it still authenticates - same model, same observations,
  same non-target specification, nothing altered since - the smoother skips the
  per-step re-derivation of a Cholesky factor for every filtered covariance the
  filter kernel had already floored. On a 3001-point dyad, `da_smooth()` with a
  supplied filter drops from 27.98 ms to 0.83 ms, of which 0.68 ms is the
  smoother kernel itself; on a 4001-point three-dimensional ENSO system it
  drops from 472.7 ms to 421.1 ms, and on the ENSO scalar partition from
  368.3 ms to 326.7 ms. The smoother it returns is identical to the bit, and so
  is every condition raised for a path that is not accepted. Any other supplied
  path - one from a different model or grid, one that has been edited, one that
  has been through `saveRDS()` - is validated in full exactly as before, and
  the new `force_validate = TRUE` argument to `da_smooth()` asks for that
  unconditionally. Sealing costs `da_filter()` 0.58% on the dyad and 0.009% on
  ENSO, and no memory: the seal holds references to the path's own arrays. It
  does make a filter path 36% (dyad) to 55% (ENSO) larger when written with
  `saveRDS()`, which `object.size()` over-reports further still.

* Compiling a model built by `cgns_from_affine()` is faster, with every
  realised coefficient identical to the bit. The coefficient path is now
  assembled from flat contiguous writes and given its dimensions once, the
  affine differencing is applied to a whole coefficient block at a time, the
  closures and hidden basis vectors are bound once for the whole grid, the
  cross-noise test is a single whole-path reduction, and the observation-Gram
  inverses are built in a tighter loop. Automatically generated zero
  cross-channels also stop re-evaluating their partner coefficient function on
  every call to read its channel width, which speeds up generic-closure
  compilation and simulation for the same reason. On a 4001-point
  three-dimensional ENSO system a complete compile drops from 525 ms to 343 ms,
  the same compile under non-target conditioning from 438 ms to 270 ms, a
  generic-closure compile from 176 ms to 155 ms, and a complete public `aci()`
  call from 750 ms to 572 ms. A partner diffusion block whose channel count
  changes part-way along the observation grid, which the shared-channel
  contract does not allow, is now reported rather than silently followed.

* The compiled matrix filter and smoother and the matrix branch of
  `gaussian_kl_path()` are substantially faster, with every value unchanged to
  the bit. The coefficient arrays are now sliced in place instead of being
  rebuilt into a list at every step, the filtered covariance and the path
  covariances are factorised without repeating validation the caller has
  already done, and the transposes and diagonal extractions are taken off the
  per-step path. On a 4001-point three-dimensional ENSO system the warm
  smoother drops from 199.0 ms to 80.7 ms, the warm matrix KL path from
  145.5 ms to 68.6 ms, the warm filter from 216.7 ms to 93.5 ms with the
  predictive likelihood and from 151.4 ms to 45.8 ms without it, and a complete
  public `aci()` call from 1066 ms to 794 ms. Filter and smoother moments,
  log-likelihoods, ACI and its decomposition, the lag table and the causal
  influence ranges are all bit-identical, on the explicit and implicit
  steppers, with and without substepping, and on conditional non-target
  configurations.

* Seeding a `simulate()` call is now contained. When `seed` is supplied the
  caller's `.Random.seed` is saved and restored on exit, so a reproducible path
  no longer moves the caller's generator forward; a session that has never
  drawn gets its state materialised first, so there is something to put back.
  `set.seed(seed)` still governs the draws, so seeded paths are bit-identical
  to before. Unseeded calls are unchanged: they draw from, and advance, the
  caller's stream. A rejected call leaves the stream where it was. Overhead on
  a 3001-point simulate is below the noise floor.

* `da_filter()` and `aci()` gain a `loglik` argument, `TRUE` by default, so
  existing calls are unchanged. With `loglik = FALSE` the filter skips the
  predictive log-likelihood entirely and leaves `meta$loglik` `NULL`; the
  filter moments, the smoother, and every ACI quantity are bit-identical
  either way. ACI never reads the likelihood, and it is a large share of the
  filter: on a 3001-point dyad the warm filter drops from 1.88 ms to 0.85 ms
  and a complete `aci()` call from 2.85 ms to 1.82 ms, and on a 4001-point
  three-dimensional ENSO system the warm filter drops from 197.9 ms to
  137.3 ms. Turning the flag on costs nothing measurable.

* The scalar Gaussian KL dispersion term is now computed as
  `0.5 * (delta - log1p(delta))` with `delta = R_p/R_q - 1`, replacing the
  algebraically equivalent but cancellation-prone
  `0.5 * (R_p/R_q - 1 + log(R_q) - log(R_p))`. Near `R_p == R_q`, which is
  where a converged ACI run sits, the old form lost the whole term: at a
  variance ratio within 1e-8 of one it returned half the true value, and away
  from `R_q == 1` it was wrong by orders of magnitude. Relative accuracy
  improves by four to eight orders of magnitude across variance ratios from
  1e-4 to 1e-13 of the tie. Every hidden-dimension-1 model is affected,
  including the scalar ENSO partitions and conditional ACI on a single hidden
  variable; the matrix KL is unchanged. Well-conditioned values move by at most
  1e-15, and the term is slightly cheaper (one fewer `sqrt`, one fewer `log`).
