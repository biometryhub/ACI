# Causal influence range

Summarizes the duration of influence on the discrete time grid, forward
from each anchor time. A finite adaptive table is labelled
`objective_on_truncated_table`; its `tail_bound` field is a heuristic
tail estimate and must not be interpreted as a certified error bound.
The `l1_linf` estimators are ratios, integrated with composite Simpson
over the whole span, following the ACI reference code.

## Usage

``` r
aci_range(x, direction = c("forward", "backward"), ...)

# S3 method for class 'lag_table'
aci_range(
  x,
  direction = c("forward", "backward"),
  method = c("exact", "l1_linf"),
  epsilon = NULL,
  min_M = "auto",
  masked_value = c("na", "zero"),
  quadrature = c("simpson", "sum", "matlab_eps_grid"),
  simpson_close = c("quadratic", "trapezoid"),
  anchors = NULL,
  convention = c("count", "lag_time"),
  epsilon_grid = NULL,
  ...
)

# S3 method for class 'aci_result'
aci_range(x, direction = c("forward", "backward"), ...)
```

## Arguments

- x:

  A `lag_table` or `aci_result` object.

- direction:

  `"forward"`, the only direction in this release. See above.

- ...:

  Arguments passed to methods.

- method:

  The objective functional. `"exact"` is the definitional objective
  range, the subjective range averaged over every threshold; reading the
  range off the running maximum makes that average a finite sum,
  `dt * sum(suffix max) / M`, with no quadrature error. `"l1_linf"` is
  the efficient ratio the ACI reference script computes,
  `dt * integral(row) / M`. They are different functionals, not two
  quadratures of one: they coincide only where the divergence decreases
  with lag.

- epsilon:

  Optional numeric vector of finite non-negative thresholds at which the
  subjective range is also reported. Its read-out convention is set by
  `convention`. This vector reports thresholds and nothing else; the
  MATLAB-compatibility quadrature takes its nodes from `epsilon_grid`, a
  separate argument, so that one vector is never asked to do both jobs.

- min_M:

  Either `"auto"` or one finite non-negative number; ranges whose
  strength falls below it are masked. It also sets the strength floor
  the `status` vocabulary is judged against, and `min_M = NULL` turns
  off the masking but not the status, which then falls back to the
  `aci.cir_min_strength` option.

- masked_value:

  What a masked range reports: `"na"` (the default) keeps the masking
  visible as `NA`; `"zero"` follows andreou2026cir Remarks B.4 and C.4
  and the FBCIR scripts, which set the length to 0. The paper's forward
  figures are computed on an untruncated table; the adaptive truncation
  of the default
  [`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
  is a package storage device whose use is recorded in the result's
  `bound` label, and `tol = 0` reproduces the untruncated convention.

- quadrature:

  How a row is reduced to its objective range, recorded in
  `meta$quadrature`. For `method = "l1_linf"`: `"simpson"` (the default)
  follows the ACI and EnKBS codebases' active lines, and `"sum"` is the
  literal andreou2026cir eq. G.8 L1 grid-function sum, which the FBCIR
  scripts' active lines use. For `method = "exact"` the threshold
  average is a finite sum, so there is no time-axis quadrature to choose
  and `"simpson"` and `"sum"` return the same exact value;
  `"matlab_eps_grid"` instead reproduces the reference script's
  `defn_objective_CIR`, a Simpson quadrature of the subjective read-out
  over the threshold grid `epsilon_grid`. `"matlab_eps_grid"` is defined
  for `method = "exact"` only.

- simpson_close:

  Closing rule for the leftover interval when a grid has an even number
  of points: `"quadratic"` (the default) fits a quadratic through the
  last three points, following `simps.m` in the ACI reference code;
  `"trapezoid"` is the package's pre-0.0.21 rule and reproduces results
  reported by earlier versions. Odd-length grids are unaffected.

- anchors:

  Optional whole numbers indexing the anchor times to report; `NULL`
  (the default) reports every anchor time. Only the requested rows are
  formed and reduced, so a handful of anchors on a long record costs a
  fraction of the whole-record computation. The result is subset to the
  anchors asked for, in the order asked for, and records them in
  `meta$anchors`.

- convention:

  The subjective read-out. `"count"` (the default) is the ACI reference
  script's `subj_CIR_idx * dt`, the 1-based index of the last exceedance
  times the step. `"lag_time"` is andreou2026cir eq. G.7, the lag time
  of that exceedance with the first cell at lag 0, one grid step below.
  Cells with no exceedance report 0 under both. The objective range is
  not affected: `"exact"` and `"matlab_eps_grid"` both average the
  counting read-out, which is the convention the definition is written
  in.

- epsilon_grid:

  Threshold quadrature nodes for `quadrature = "matlab_eps_grid"`,
  strictly increasing. `NULL` (the default) uses the reference script's
  513-point grid, `10^seq(-6, 0.5, length.out = 513)`. Supplying it with
  any other `quadrature` is an error, so the reporting thresholds in
  `epsilon` and these quadrature nodes can never be taken from one
  vector.

## Value

An object of class `cir_result` with `direction` `"forward"`.

## Details

Only `direction = "forward"` is in this release. The backward range is a
`FBCIR_code-main` feature: it appears in no ACI_code script, and it is
held in the development reserve (`reserve/fbcir/`) for the release that
brings that family in. `direction = "backward"` raises
`aci_error_not_implemented` rather than returning a forward answer under
a backward label.

## Methods (by class)

- `aci_range(lag_table)`: Forward range from a precomputed lag table.

- `aci_range(aci_result)`: Forward range from an ACI result, building
  the table when it was not retained. The covariance policy is not
  re-resolved here: it is read off the result's own
  `meta$regularization`, so re-analysing a saved result in a new session
  cannot silently change it.

## References

Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
identifying forward and backward causal influence ranges using
assimilative causal inference. arXiv:2510.21889v2, 4 August 2026.
[doi:10.48550/arXiv.2510.21889](https://doi.org/10.48550/arXiv.2510.21889)

## See also

[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
tb <- lag_table(m, ob, mode = "forward")
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
aci_range(tb)
#> Warning: 1 forward CIR values masked (M < 1e-05); interpret CIRs jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
#> <cir_result> forward | method = exact (layer_cake_objective) | masked/NA: 1 of 201
#>   status: resolved 21, censored 178, insufficient 2
```
