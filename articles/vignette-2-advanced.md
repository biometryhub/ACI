# The closed-form ACI engine

*Assimilative causal inference (ACI) in R*
([`vignette("vignette-1-intro")`](https://biometryhub.github.io/ACI/articles/vignette-1-intro.md))
runs the engine end to end. This one is about the machinery underneath:
how models are written, what the filter, smoother and online-smoother
arguments do, what happens when a covariance leaves the
positive-definite cone, which conventions the causal influence range is
computed under, and which of the conditional questions are genuinely
different estimands rather than different estimators.

Every convention below is a choice that changes the number reported. The
package states them on the result objects rather than leaving them to be
inferred, and this vignette is largely a guide to reading those
statements.

Where a choice exists only to match the reference MATLAB codebase, it is
named here and worked through in *Reproducing the reference MATLAB
codebase*
([`vignette("vignette-3-matlab")`](https://biometryhub.github.io/ACI/articles/vignette-3-matlab.md)),
which also covers the benchmark constructors’ provenance and the kinds
of evidence behind the package’s fidelity claims.

## 1. Models and observations

### The CGNS class

ACI’s closed form rests on a structural property: **given the observed
path, the hidden dynamics are linear and Gaussian.** A conditional
Gaussian nonlinear system (CGNS) is written

\mathrm{d}x = \big\[L_x(t,x)\\y + f_x(t,x)\big\]\mathrm{d}t +
\Sigma\_{x,1}\mathrm{d}W_1 + \Sigma\_{x,2}\mathrm{d}W_2, \qquad
\mathrm{d}y = \big\[L_y(t,x)\\y + f_y(t,x)\big\]\mathrm{d}t +
\Sigma\_{y,1}\mathrm{d}W_1 + \Sigma\_{y,2}\mathrm{d}W_2 .

The coefficients may be arbitrarily nonlinear in (t, x). What is
forbidden is nonlinearity in y. Two Wiener channels are carried
separately so that observation and signal noise can be correlated, which
the filter’s gain accounts for.

### Writing a model

[`aci_model()`](https://biometryhub.github.io/ACI/reference/aci_model.md)
takes the blocks directly:

``` r

m2 <- aci_model(
  Lx  = function(t, x) matrix(c(1, 0), 2, 1),
  fx  = function(t, x) -0.5 * x,
  Ly  = function(t, x) matrix(-0.5, 1, 1),
  fy  = function(t, x) 0,
  Sx1 = function(t, x) matrix(c(0.5, 0.3, 0, 0.4), 2, 2),
  Sy2 = function(t, x) matrix(1, 1, 1),
  k = 2, l = 1, name = "correlated observation noise")
m2
#> <cgns_model> 'correlated observation noise': k = 2 observed, l = 1 hidden
```

[`aci_model_from_affine()`](https://biometryhub.github.io/ACI/reference/aci_model_from_affine.md)
takes the full drift instead and does the split for you, after checking
that the split exists. A drift that is not affine in the hidden state is
rejected rather than silently linearised:

``` r

aci_model_from_affine(
  f_full    = function(t, x, y) -0.5 * x + y^2,
  g_full    = function(t, x, y) -0.5 * y,
  Sx        = function(t, x) matrix(0.5, 1, 1),
  Sy_hidden = function(t, x) matrix(1, 1, 1),
  k = 1, l = 1)
#> Error in `aci_model_from_affine()`:
#> ! Observed drift is not affine in the hidden state.
```

``` r

hand_dyad <- aci_model_from_affine(
  f_full    = function(t, x, y) -0.5 * x + 2 * x * y + 0.5,
  g_full    = function(t, x, y) -0.5 * y - 2 * x^2 + 1,
  Sx        = function(t, x) matrix(0.5, 1, 1),
  Sy_hidden = function(t, x) matrix(1, 1, 1),
  k = 1, l = 1, name = "hand-built dyad")
hand_dyad
#> <cgns_model> 'hand-built dyad': k = 1 observed, l = 1 hidden
```

[`aci_linear_model()`](https://biometryhub.github.io/ACI/reference/aci_linear_model.md)
is the shortcut for the linear-Gaussian case, where the two coupling
blocks and the two diffusions are constants rather than functions of (t,
x). It expands scalars to the right shapes and returns an ordinary
`cgns_model` object, so the engine takes it like any other:

``` r

lin <- aci_linear_model(lambda_x = 1, lambda_y = -0.5,
                        fx = 0, fy = 0,
                        sigma_x = 0.5, sigma_y = 1,
                        k = 1, l = 1)
lin
#> <aci_linear_model> 'conditionally_linear': k = 1 observed, l = 1 hidden
c(Lx = lin$Lx(0, 1), Ly = lin$Ly(0, 1))
#>   Lx   Ly 
#>  1.0 -0.5

s_lin <- simulate(lin, seed = 2, t_end = 4, dt = 0.01, burn_in = 1)
max(aci(lin, as_obs(s_lin), init = list(mean = 0, cov = matrix(1, 1, 1)))$aci)
#> [1] 1.506861
```

Coefficient functions are **mathematical coefficients**: for a fixed
(t,x) they must return deterministic values with stable shapes. The
closed-form routes may realise each coefficient once on the observation
grid and reuse that realised path, so random draws or result-changing
mutable state inside a coefficient function are outside the contract.

### Benchmark constructors

Three constructors carry the reference MATLAB codebase’s models:
[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md),
[`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md)
and
[`aci_predprey_model()`](https://biometryhub.github.io/ACI/reference/aci_predprey_model.md).
Each records its provenance and how far the correspondence has been
checked. *Reproducing the reference MATLAB codebase* goes through those
records model by model and exercises the predator-prey partitions.

``` r

sapply(list(dyad          = aci_dyad_model(),
            enso6         = aci_enso_model(hidden = "tau"),
            predator_prey = aci_predprey_model(hidden = "prey")),
       function(mm) c(observed = mm$k, hidden = mm$l))
#>          dyad enso6 predator_prey
#> observed    1     5             1
#> hidden      1     1             1
```

Two of those records matter for the rest of this vignette. First,
[`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md)
covers five hidden partitions, and the fifth, T_C hidden, is an
approximation rather than an exact split: it is built only when the
substitution is named, it assimilates a reduced two-channel observed
process, and it refuses to simulate. Second, coefficient and inference
agreement is the claim made here. Pathwise simulator agreement is a
separate claim: [`simulate()`](https://rdrr.io/r/stats/simulate.html) is
Euler-Maruyama throughout, and `meta$matlab_simulator_parity` records
where the reference does something else.

### Observations

[`observed_trajectory()`](https://biometryhub.github.io/ACI/reference/observed_trajectory.md)
derives `dt` from the supplied times and refuses a non-uniform grid. A
desynchronised time column is an error, not a wrong answer:

``` r

observed_trajectory(c(0, 0.01, 0.03), matrix(c(1, 2, 3), ncol = 1))
#> Error in `observed_trajectory()`:
#> ! v0 requires a uniform time grid (invariant #11); resample first.
```

Naming the columns is what makes the conditional questions of section 4
readable, because
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
accepts channel names.

## 2. The engine: filter, smoother, metric

``` r

m    <- aci_dyad_model()
sim  <- simulate(m, seed = 7, t_end = 8, dt = 0.01, burn_in = 1)
ob   <- observed_trajectory(sim$obs$t, sim$obs$x, names = "x")
init <- list(mean = 0, cov = matrix(1, 1, 1))

filt <- aci_filter(m, ob, init = init)
smoo <- aci_smoother(m, ob, filter = filt)
a    <- aci(m, ob, init = init)
```

### Covariances are strict

A covariance that leaves the positive-definite cone inside a state
recursion, a metric input or the likelihood stops the run. The condition
is classed and names the site, the grid index, the time and the
offending value, rather than reporting a number the recursion did not
produce. Take the same record on a twentyfold coarser grid, where the
explicit Riccati step overshoots:

``` r

idx    <- seq.int(1L, length(ob$t), by = 20L)
ob_bad <- observed_trajectory(ob$t[idx], ob$x[idx, , drop = FALSE],
                              names = "x")
aci_filter(m, ob_bad, init = init)
#> Error in `.aci_stop_cov()`:
#> ! The filter covariance must stay finite and positive definite; it reached -0.1336333 at index 37 (time 7.2), in the explicit Riccati step. Reduce dt, raise nsub, or use stepper = "implicit", which preserves positivity. To keep the previous behaviour, call with regularize = "floor"; every floored step is then recorded in the result's meta$regularization.
```

`regularize = "floor"` is the previous behaviour, available on
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
and
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md),
and `options(aci.regularize = "floor")` sets it session-wide. It does
not make the run correct. It makes what it did visible:

``` r

bad <- aci_filter(m, ob_bad, init = init, regularize = "floor")
#> Warning in .cgns_filter_scalar(bundle, init, nsub, validate = validate, :
#> Explicit Riccati step is unstable (max ||Lx' gxx^-1 Lx R|| dt = 2.44 > 1): the
#> covariance can overshoot, leave the positive-definite cone, and oscillate. Use
#> the positivity-preserving implicit stepper, or reduce dt / increase nsub.
bad$meta$regularization[c("policy", "fired", "n_events")]
#> $policy
#> [1] "floor"
#> 
#> $fired
#> [1] TRUE
#> 
#> $n_events
#> [1] 1
bad$meta$regularization$sites
#>              site              role n first_index first_time worst_value
#> 1 filter_explicit filter covariance 1          37        7.2  -0.1336333
```

The record is always present, so a clean run says so explicitly rather
than by omission:

``` r

a$meta$regularization[c("policy", "fired", "n_events")]
#> $policy
#> [1] "none"
#> 
#> $fired
#> [1] FALSE
#> 
#> $n_events
#> [1] 0
nrow(a$meta$regularization$sites)
#> [1] 0
```

[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
on an `aci_result` reads the policy off the result rather than off the
option, so re-analysing a saved result in a later session cannot
silently change it.
[`spd_floor()`](https://biometryhub.github.io/ACI/reference/spd_floor.md)
and
[`safe_chol()`](https://biometryhub.github.io/ACI/reference/safe_chol.md)
are exported. They are the implementation of the opt-in, not a default.
Called directly,
[`spd_floor()`](https://biometryhub.github.io/ACI/reference/spd_floor.md)
lifts eigenvalues to a floor and
[`safe_chol()`](https://biometryhub.github.io/ACI/reference/safe_chol.md)
factors through a jitter ladder, and gives up rather than returning a
factor of something that was never a covariance:

``` r

spd_floor(matrix(c(1, 0, 0, -1e-14), 2, 2))
#>      [,1]  [,2]
#> [1,]    1 0e+00
#> [2,]    0 1e-12
safe_chol(matrix(c(2, 1, 1, 2), 2, 2))
#>          [,1]      [,2]
#> [1,] 1.414214 0.7071068
#> [2,] 0.000000 1.2247449
try(safe_chol(matrix(-1, 1, 1)))
#> Error in safe_chol(matrix(-1, 1, 1)) : 
#>   Matrix (covariance) is not positive definite even after jitter ladder; min eig = -1.000e+00.
```

### Steppers

The default explicit stepper is the reference scheme, and it is the
convention the lag table of section 3 assumes. A positivity-preserving
implicit stepper is available for stiff models or very diffuse priors,
with `nsub` sub-steps per observation:

``` r

filt_i <- aci_filter(m, ob, init = init, stepper = "implicit", nsub = 5)
a_i    <- aci(m, ob, init = init, stepper = "implicit", nsub = 5)

c(filter_mean_gap = max(abs(filt$mean - filt_i$mean)),
  aci_gap         = max(abs(a$aci - a_i$aci)),
  peak_aci        = max(a$aci))
#> filter_mean_gap         aci_gap        peak_aci 
#>      0.03779429      0.07789620      1.51958662
```

The two steppers are different discretisations, so they give different
numbers, and the difference here is not negligible against the size of
the metric itself. Pick one deliberately and record which.

[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
refuses anything else outright, because the Theorem 3 recursions it
accumulates are exact for the explicit single-step discretisation and
for no other:

``` r

lag_table(m, ob, mode = "forward", init = init, stepper = "implicit")
#> Error in `lag_table()`:
#> ! lag_table requires stepper = 'explicit' and nsub = 1.
```

### The predictive likelihood

The filter chain can accumulate a predictive log-likelihood. ACI never
uses it, so it is switchable:

``` r

filt$meta$loglik
#> [1] 1261.46
is.null(aci_filter(m, ob, init = init, loglik = FALSE)$meta$loglik)
#> [1] TRUE
```

`loglik = FALSE` leaves every ACI quantity unchanged and skips the work.
The saving on the warm filter is recorded in the development record
(section *Evidence*).

### Supplying a filter to the smoother

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
seals the path it returns with a private token, and
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
skips the per-step covariance re-validation of a path that still
authenticates against the run being smoothed. Any other supplied path,
including one that has been through
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html), is validated in full
as before, and `force_validate = TRUE` validates unconditionally. The
smoother result is the same either way:

``` r

identical(aci_smoother(m, ob, filter = filt)$mean,
          aci_smoother(m, ob, filter = filt, force_validate = TRUE)$mean)
#> [1] TRUE
```

The development record has the measurement behind this: on a 3001-point
dyad, validation was 97.7% of
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)’s
time when the documented
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
then
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
pattern was used.

### The metric, and its parts

[`aci_metric_pair()`](https://biometryhub.github.io/ACI/reference/aci_metric_pair.md)
and
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
are exposed directly, so the metric can be recomputed from posteriors
obtained any other way:

``` r

aci_metric_pair(mu_p = 0, R_p = matrix(1, 1, 1),
                mu_q = 0.5, R_q = matrix(2, 1, 1))
#>      total     signal dispersion 
#> 0.15907359 0.06250000 0.09657359
max(abs(aci_metric(smoo, filt)$total - a$aci))
#> [1] 0
```

The signal and dispersion parts are each floored at zero independently
before being summed, so a decomposition never reports a negative part
from round-off.

## 3. The forward causal influence range

### Two routes, and what they cost

The finite-lag divergences can be stored or streamed.

- [`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
  **retains** the rows. Cost is proportional to the retained lag cells,
  which is O(N^2) for a full table, and the object holds them. The table
  can then be reduced repeatedly at no further cost.
- [`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)
  on an `aci_result` **streams**, forming and reducing one row at a time
  and keeping none. Memory is linear.

Both reduce identically. On this record, asking the streamed route for a
subsample of anchors costs a fraction of the whole-record computation:

``` r

anchors <- seq(1L, length(ob$t), by = 8L)
fc <- aci_range(a, anchors = anchors)
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
l1 <- aci_range(a, anchors = anchors, method = "l1_linf")
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
fc
#> <cir_result> forward | method = exact (objective_on_truncated_table) | masked/NA: 1 of 101
#>   status: resolved 55, censored 45, insufficient 1
```

`anchors` selects which anchor times are reported, and the streamed
engine forms only those rows, skipping the per-interval primitives
before the earliest of them. Asking for fifteen anchors instead of every
one is a saving in work, not a subsample taken afterwards, and the
values it returns are the same ones to the bit:

``` r

ob_w <- observed_trajectory(ob$t[1:401], ob$x[1:401, , drop = FALSE],
                            names = "x")
a_w  <- aci(m, ob_w, init = init)
w15  <- seq(1L, 401L, by = 28L)

t_all <- system.time(f_all <- aci_range(a_w))[["elapsed"]]
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
t_w15 <- system.time(f_w15 <- aci_range(a_w, anchors = w15))[["elapsed"]]

c(all_401_anchors = t_all, fifteen_anchors = t_w15, ratio = t_all / t_w15)
#> all_401_anchors fifteen_anchors           ratio 
#>           0.077           0.007          11.000
identical(f_all$tau[w15], f_w15$tau)    # the same values, not a resample
#> [1] TRUE
```

The saving is on the streamed route only.
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
still builds every row, so a retained table plus `aci_range(anchors = )`
pays the whole quadratic cost and saves only the reduction. The
development record has the measurements behind that split.

What `anchors` does not do is cap how far forward each row looks. Row
length is governed separately, at the end of this section.

### Table storage arguments

When the table is retained, three arguments govern how much of it is
kept. They trade storage and time against the tail of each row:

- `tol` is the level below which a row becomes a candidate for freezing.
- `window` is how many consecutive steps a row must stay below `tol`
  before it is actually frozen, which stops a single dip from truncating
  a row that later recovers.
- `max_lag` is a hard cap on the retained positive lag.

Truncation is not free and it is not hidden. Measured against an
untruncated reference on a shorter record:

``` r

ob_s <- observed_trajectory(ob$t[1:301], ob$x[1:301, , drop = FALSE],
                            names = "x")
reduce_tab <- function(...) {
  tb <- lag_table(m, ob_s, mode = "forward", init = init, ...)
  f  <- aci_range(tb, min_M = 0)
  c(mean_retained_lag = round(mean(tb$L), 1),
    Mb                = round(as.numeric(object.size(tb)) / 2^20, 2),
    max_tail_bound    = signif(max(f$tail_bound), 2),
    tau               = list(f$tau))
}
untruncated <- reduce_tab(tol = 0)
compare <- function(...) {
  r <- reduce_tab(...)
  c(unlist(r[1:3]),
    max_tau_error = signif(max(abs(r$tau - untruncated$tau), na.rm = TRUE), 2))
}
rbind(`tol=0 (reference)` = compare(tol = 0),
      `tol=1e-2`          = compare(tol = 1e-2),
      `tol=1e-1`          = compare(tol = 1e-1),
      `max_lag=50`        = compare(tol = 0, max_lag = 50))
#>                   mean_retained_lag   Mb max_tail_bound max_tau_error
#> tol=0 (reference)             150.0 0.45         0.0000        0.0000
#> tol=1e-2                      116.9 0.37         0.0099        0.0012
#> tol=1e-1                       81.4 0.29         0.0990        0.0550
#> max_lag=50                     45.8 0.21         0.3900        0.6300
```

Three accessors read a retained table without reaching into its
internals.
[`lt_diag()`](https://biometryhub.github.io/ACI/reference/lt_diag.md)
returns the zero-lag diagonal,
[`lt_row()`](https://biometryhub.github.io/ACI/reference/lt_row.md)
returns one anchor’s row padded out to full length, and
[`lt_tail_bound()`](https://biometryhub.github.io/ACI/reference/lt_tail_bound.md)
returns the heuristic tail estimate. `lt_row(pad = "na")` is the one to
reach for when the question is *where the table stopped* rather than
what the row was worth, because the zero padding of the default is
indistinguishable from a genuinely decayed tail:

``` r

tb_s <- lag_table(m, ob_s, mode = "forward", init = init, tol = 1e-2)
r10  <- lt_row(tb_s, 10, pad = "na")

c(row_cells_padded_out = length(r10),
  cells_actually_kept  = sum(!is.na(r10)),
  retained_lag_L       = tb_s$L[10])
#> row_cells_padded_out  cells_actually_kept       retained_lag_L 
#>                  292                  192                  191
lt_tail_bound(tb_s, 10)
#> [1] 0.007551115
round(utils::head(lt_row(tb_s, 10), 5), 4)
#> [1] 0.2442 0.2292 0.1824 0.1905 0.2109
```

Truncation is recorded on the reduced result as well: `bound` says
whether the number came off a truncated table, and `tail_bound`
estimates what the frozen tails would have contributed.

``` r

c(bound = fc$bound, l1_bound = l1$bound)
#>                                 bound                              l1_bound 
#>        "objective_on_truncated_table" "lower_ratio_on_truncated_table_only"
max(fc$tail_bound)
#> [1] 9.342417e-09
```

`tail_bound` is an estimate, not a certified error bound. In the table
above it tracks the realised error to within an order of magnitude for
the `tol` rule, and overstates it for the hard `max_lag` cap. Setting
`options(aci.default_tol = 0)` disables truncation and reproduces the
untruncated convention of the published figures. Note that the streamed
route reads `tol` from that option rather than from an argument.

`mode = "full"` retains the negative-lag cells as well, at
correspondingly more work, and requires `max_lag = Inf`.

### Three named functionals

The objective range can be asked for in three ways, and the names matter
more than usual here.

- `method = "exact"` (the default) is the **definitional objective**:
  the subjective range averaged over every threshold. Read off the
  running maximum of the row, that average collapses to a finite sum,
  `dt * sum(suffix max) / M`. There is no quadrature in it, which is why
  `quadrature = "simpson"` and `quadrature = "sum"` return the same
  number for this method.
- `method = "l1_linf"` is the **efficient ratio** the reference script
  computes, `dt * integral(row) / M`.
- `method = "exact", quadrature = "matlab_eps_grid"` is the
  **MATLAB-compatibility mode**, reproducing the reference script’s
  `defn_objective_CIR` by taking a Simpson quadrature of the subjective
  read-out over a threshold grid.

The first two are different functionals, not two quadratures of one.
They coincide only where the divergence decreases with lag. On this
record they differ by a factor of

``` r

round(max(fc$tau / l1$tau, na.rm = TRUE), 3)
#> [1] 3.169
```

This is a naming hazard the comparison of the two predecessor packages
flagged: the word “objective” in one of them names the *ratio*, not the
definitional form. Mapping by name alone pairs quantities that differ by
a factor of up to 1.9 on the shared fixture. Read `method`, `bound` and
`meta$quadrature` off the result rather than trusting a remembered
default.

The first and third are the *same* functional on different axes. The
compatibility mode carries a quadrature error that the exact sum does
not, and refining its grid moves it towards the exact value:

``` r

mg   <- aci_range(a, anchors = anchors, quadrature = "matlab_eps_grid")
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
mg4k <- aci_range(a, anchors = anchors, quadrature = "matlab_eps_grid",
                  epsilon_grid = 10^seq(-6, 0.5, length.out = 4097))
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
mg$meta$quadrature
#> [1] "matlab_eps_grid"
c(default_513_nodes = max(abs(mg$tau   - fc$tau), na.rm = TRUE),
  refined_4097      = max(abs(mg4k$tau - fc$tau), na.rm = TRUE))
#> default_513_nodes      refined_4097 
#>       0.006239067       0.001413819
```

Use it for parity with the reference script, and the default for
anything else.

Two arguments fix the ratio’s quadrature:

- `quadrature`: `"simpson"` (the default) follows the ACI reference
  code’s active line, `"sum"` is the literal L1 grid-function sum used
  by the FBCIR scripts.
- `simpson_close`: the closing rule for the leftover interval when a
  grid has an even number of points. `"quadratic"` reproduces `simps.m`,
  `"trapezoid"` reproduces results reported by earlier package versions.

``` r

sum_rule <- aci_range(a, anchors = anchors, method = "l1_linf",
                      quadrature = "sum")
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
c(max_abs = max(abs(sum_rule$tau - l1$tau), na.rm = TRUE),
  max_rel = max(abs(sum_rule$tau - l1$tau) / l1$tau, na.rm = TRUE))
#>    max_abs    max_rel 
#> 0.00794573 0.26913812
```

A quadrature choice is worth tens of percent here, so it is not a
formatting detail. `simpson_close` only bites on rows with an even
number of points, so on a record whose rows are all odd-length the two
closing rules agree exactly. That is a property of the record, not
evidence that the argument is inert.

### The subjective range, and the two read-out conventions

The \varepsilon-resolved **subjective** range, the lag of the last
exceedance of each threshold, is reached through `epsilon`:

``` r

sc <- aci_range(a, anchors = anchors,
                epsilon = 10^seq(-4, 0, length.out = 30))
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
dim(sc$subjective)
#> [1] 101  30
sc$meta$convention
#> [1] "count"
```

There are two published read-outs, one grid step apart, and `convention`
names them:

- `"count"` (the default) is the ACI reference script’s
  `subj_CIR_idx * dt`, the 1-based index of the last exceedance times
  the step.
- `"lag_time"` is the paper convention, the lag time of that exceedance
  with the first cell at lag 0.

``` r

lag_time <- aci_range(a, anchors = anchors,
                      epsilon = 10^seq(-4, 0, length.out = 30),
                      convention = "lag_time")
#> Warning in .forward_cir_compiled(bundle, filter = x$paths$filter, init =
#> x$handles$init, : 1 forward CIR values masked (M < 1e-05); interpret CIRs
#> jointly with the ACI metric (Andreou & Chen 2026, Remark B.4).
gap <- sc$subjective - lag_time$subjective
c(dt              = sc$dt,
  min_gap         = min(gap),
  max_gap         = max(gap),
  only_0_or_dt    = all(abs(gap) < 1e-12 | abs(gap - sc$dt) < 1e-12))
#>           dt      min_gap      max_gap only_0_or_dt 
#>         0.01         0.00         0.01         1.00
```

The difference is exactly `dt`, or zero where nothing exceeded the
threshold. Neither is an error. The one-step offset was verified to the
unit in last place across 51,471 cells of the shared comparison record
(development record).

`epsilon` and `epsilon_grid` are separate arguments on purpose.
`epsilon` reports thresholds, and `epsilon_grid` supplies quadrature
nodes to the compatibility mode. Asking one vector to do both jobs drove
a predecessor implementation’s objective negative on a sparse grid, so
the package refuses the combination outright:

``` r

aci_range(a, anchors = anchors, epsilon_grid = c(1e-3, 1e-2))
#> Error in `.cir_compat_grid()`:
#> ! epsilon_grid supplies quadrature nodes for quadrature = 'matlab_eps_grid' only; use epsilon for reporting thresholds.
```

### Statuses and masking

``` r

table(fc$status)
#> 
#>        resolved        censored below_threshold    insufficient 
#>              55              45               0               1
sum(is.na(fc$tau))
#> [1] 1
```

`status` classifies each anchor, and the four levels are not four grades
of the same thing: two of them report a range and two of them decline
to.

- `resolved`: the record comfortably outlasts the influence it measured,
  so the reported range is supported by the observations that follow the
  anchor.
- `censored`: it does not. The last exceedance of the strength floor
  sits in the later half of the row, leaving less record after the
  influence than the influence itself occupied, and the reported range
  should be read as a lower bound. The test is scale-free and taken
  entirely from the row, with no horizon or margin supplied by the
  caller.
- `below_threshold`: the row’s peak divergence never cleared the floor,
  so the range is masked rather than reported as an inflated number.
- `insufficient`: fewer than three observations follow the anchor. The
  final anchor is always `insufficient`.

Censoring is judged inside the row, not against a caller-supplied
horizon. The obvious test, “did the exceedance run into the last cell”,
can never fire: the last cell of a complete row is identically zero,
because given the whole record the online estimate is the smoother it is
scored against. The test used instead compares the record remaining
after the last exceedance against the extent of the exceedance itself. A
record that does not outlast the influence it measured did not resolve
it, and its range is a lower bound. The decline in \tau_f at the end of
a record is the same finite-record artefact seen from the other side.

Masking is not smoothing. `min_M` masks anchors whose peak divergence is
too small for the ratio to mean anything, `masked_value = "na"` keeps
that visible, and `masked_value = "zero"` restores the published zero
convention. A masked anchor is a refusal to report, not an estimate of
zero influence.

### Row length

`anchors` caps how many rows are formed. It does not cap how far forward
each of them looks, and on a long record the far tail of a row
contributes nothing but cost. On the streamed route, row length is
governed by the truncation rule through `options(aci.default_tol = )`.
On a retained table it is `max_lag`, which applies to the whole table
rather than per anchor. `lag_table(window = )` means something else, the
consecutive-steps-below-`tol` freeze rule of the previous section.

## 4. Conditional ACI

### The construction

With several observed effects, “does y influence the observations”
becomes “does y influence x_A once x_B is accounted for”.
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
names x_B, the **complement** of the causal question, and offers two
constructions:

- `method = "mask"` keeps x_B in the dynamics but gives its innovations
  no assimilation weight.
- `method = "reduce"` substitutes x_B as known forcing and removes it
  from the observed block entirely.

These are different operators. `method = "reduce"` is admissible only
when the observation-noise Gram has no A-B cross-block, and the package
checks that along the whole path rather than at its first point:

``` r

s2 <- simulate(m2, seed = 1, t_end = 1, dt = 0.01)
o2 <- observed_trajectory(s2$obs$t, s2$obs$x, names = c("xA", "xB"))
aci_conditional_reduce(m2, o2, aci_conditional("xB", "reduce"))
#> Error in `aci_conditional_reduce()`:
#> ! gxx has a nonzero A-B cross-block; use aci_conditional(method = 'mask') (SPEC-01 s6, pending SI equivalence transcription).
```

`method = "mask"` carries no such restriction and works on the same
model:

``` r

init2 <- list(mean = 0, cov = matrix(1, 1, 1))
c(total = max(aci(m2, o2, init = init2)$aci),
  conditional = max(aci(m2, o2, init = init2,
                        conditional = aci_conditional("xB", "mask"))$aci))
#>       total conditional 
#>   0.7190804   0.5099931
```

On models whose observation noise *is* diagonal and uncorrelated with
the hidden noise, the two methods agree. That was verified bit for bit
across the ENSO partitions on the shared comparison record (development
record). On the shorter h_W record built here, which the rest of this
section also uses, the two paths agree to round-off rather than bitwise:
the arithmetic takes a different route to the same operator, not to a
second answer.

``` r

me     <- aci_enso_model(hidden = "hW")
se     <- simulate(me, seed = 12, t_end = 4, dt = 5e-3, burn_in = 1)
obe    <- as_obs(se)
init_e <- list(mean = me$meta$ic_default$y0, cov = matrix(0.1, 1, 1))

a_inf <- aci(me, obe, init = init_e,
             conditional = aci_conditional(target = "TC", method = "mask"))
a_pf  <- aci(me, obe, init = init_e,
             conditional = aci_conditional(target = "TC",
                                           method = "reduce"))
max(abs(a_inf$aci - a_pf$aci))
#> [1] 2.775558e-17
```

That agreement is a property of this model’s noise structure, not a
general identity, and it should not be carried over to a model whose
noise structure has not been checked. On `m2` above, whose observation
noise has an A-B cross-block, `method = "reduce"` is not merely
different. It is refused.

### State the estimand

Each of the reference conditional scripts fixes a specific target. The
h_W script studies h_W \to T_C \mid (u, T_E, \tau, I), with two further
targets present but commented out.
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
takes either side of that split, `target` or `given`, never both, and
the constructor records which estimand the script asked for and where it
read it:

``` r

me$meta$causal_link
#> [1] "(hW) -> (TC) | (u,TE,tau,I)"
me$meta$estimand_provenance
#> [1] "ENSO_model_cond_ACI_h_W_unobs.m:1202 (S_xoS_x_inv(1,1,:) = 1/sigma_C^2)"

by_target     <- aci_conditional(target = "TC", method = "mask")
by_complement <- aci_conditional(given = me$meta$vars$observed[
                             me$meta$conditioning_obs_idx],
                           method = "mask")
by_target
#> <aci_conditional_spec> x_A = {TC}, method = mask
by_complement
#> <aci_conditional_spec> x_B = {u, TE, tau, I}, method = mask

a_total <- aci(me, obe, init = init_e)
a_tc    <- aci(me, obe, init = init_e, conditional = by_target)
a_comp  <- aci(me, obe, init = init_e, conditional = by_complement)

c(total = max(a_total$aci), target_TC = max(a_tc$aci))
#>     total target_TC 
#>  2.689974  1.565727
identical(a_tc$aci, a_comp$aci)         # one split, named from either side
#> [1] TRUE
```

The choice is not cosmetic. Conditioning removes a real part of the
metric here, but only because of *which* channel it removes. The hidden
state reaches the observed drift only through the non-zero rows of L_x:

``` r

Lx_max <- sapply(seq_along(obe$t),
                 function(j) abs(drop(me$Lx(obe$t[j], obe$x[j, ]))))
setNames(round(apply(Lx_max, 1, max), 4), me$meta$vars$observed)
#>      u     TC     TE    tau      I 
#> 0.0000 0.4875 0.4875 0.0000 0.0000
```

Three of the five rows are identically zero, so conditioning on those
three is not a small correction, it is bit-for-bit the unconditional
run:

``` r

a_inert <- aci(me, obe, init = init_e,
               conditional = aci_conditional(given = c("u", "tau", "I"),
                                             method = "mask"))
identical(a_inert$aci, a_total$aci)     # masking {u, tau, I} changes nothing
#> [1] TRUE
max(abs(a_total$aci - a_tc$aci))        # removing TE changes this much
#> [1] 1.324153
```

Whether a conditional mask bites is a property of the model’s coupling
structure, and it is worth checking on any new model before reading a
conditional result as evidence. Note also what the metadata is and is
not: `causal_link` and `conditioning_obs_idx` record the estimand the
source script defines, and nothing is conditioned on until a
specification is supplied to
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) or
declared by the constructor.

### Two conditional conventions covered elsewhere

Two further conditional choices are convention-matching questions rather
than API questions, and *Reproducing the reference MATLAB codebase*
treats them with the rest of the reference-codebase material:
`aci_conditional(first_step = "matlab")`, which reproduces the reference
scripts’ unmasked first slice of the observation-precision array, and
`aci_enso_model(observations = )`, which selects between the \tau-hidden
script’s three-channel estimand and the full five-channel construction.
Both change the numbers reported, and the second changes the estimand
rather than the estimator.

## 5. Two smoothers, two schemes

The package carries two smoothing **discretisations**, and they are not
interchangeable at finite \Delta t.

- [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
  is an Euler discretisation of the **continuous** backward equation.
  This is what
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) uses for
  its headline metric.
- The **discrete** Theorem 3 posterior is what
  [`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
  accumulates and what
  [`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
  returns, and it is what every CIR quantity is scored against.

They agree only to first order in \Delta t. Every object records which
discretisation produced it, in `meta$scheme`, and an `aci_result` reads
that into `meta$smoother_scheme`. The finer implementation tag is
separate, in a path’s `meta$route`, so the label a reader needs and the
label a maintainer needs are not the same field:

``` r

tb <- lag_table(m, ob, mode = "forward", init = init, max_lag = 1)

c(aci          = a$meta$smoother_scheme,
  aci_smoother    = smoo$meta$scheme,
  lag_table    = tb$meta$scheme,
  aci_online    = aci_online(m, ob, lag = Inf, init = init)$meta$scheme)
#>                  aci         aci_smoother            lag_table 
#> "backward_ode_euler" "backward_ode_euler"  "theorem3_discrete" 
#>           aci_online 
#>  "theorem3_discrete"
c(smoother_route = smoo$meta$route,
  table_reference_smoother = tb$meta$reference_smoother)
#>           smoother_route table_reference_smoother 
#>           "backward_ode"  "thmD1_online_complete"
```

`max_lag = 1` is the cheap way to obtain the discrete diagonal without
building a table: it retains only the cells needed for the zero-lag
entry, so the cost is linear in the record rather than quadratic, and
the diagonal is the same one a full table would give.

The gap between the two schemes on this record:

``` r

c(max_abs_gap = max(abs(lt_diag(tb) - a$aci)),
  peak_aci    = max(a$aci))
#> max_abs_gap    peak_aci 
#>   0.1227121   1.5195866
```

This is a real difference between two discretisations of the same
object, not a defect in either. The comparison of the two predecessor
packages measured it on the ENSO record in both independently and got
the same number to all printed digits, 1.04276e-01 on a scale of
1.09324, which is what establishes it as a scheme property rather than a
package difference (development record).

Feeding a table to
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) swaps the
metric onto the discrete scheme deliberately, and the result says so:

``` r

a_discrete <- aci(m, ob, init = init, table = tb)
c(engine = a_discrete$meta$engine,
  scheme = a_discrete$meta$smoother_scheme)
#>              engine              scheme 
#>      "reused_table" "theorem3_discrete"
```

When comparing an ACI curve against a published figure, check which
smoother the figure used before concluding anything from a discrepancy
of this size.

### The fixed-lag online smoother

[`aci_online()`](https://biometryhub.github.io/ACI/reference/aci_online.md)
is the discrete scheme exposed on its own, at any lag. The estimate at
index j conditions on the observed record through j + \mathtt{lag} and
saturates at the end of the record, which the result records per anchor:

``` r

on <- aci_online(m, ob, lag = 20, init = init)
on
#> <da_path_gaussian> kind = online, l = 1, N+1 = 801
on$meta$lag
#> [1] 20
c(saturated_anchors = sum(on$meta$saturated),
  lag_at_the_start  = on$meta$lag_effective[1],
  lag_at_the_end    = utils::tail(on$meta$lag_effective, 1))
#> saturated_anchors  lag_at_the_start    lag_at_the_end 
#>                 0                20                 0
```

The requested lag is not always the lag obtained: the last anchors have
no record left to condition on, so `lag_effective` falls away to zero
and `saturated` marks where it did. Nothing is padded or extrapolated to
hide that.

The two ends are exact, and they are exact against different things.
`lag = 0` returns the filter moments value for value, and `lag = Inf`
returns the complete Theorem 3 posterior, which is the lag table’s own
reference smoother rather than
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md):

``` r

on_full <- aci_online(m, ob, lag = Inf, init = init)

identical(aci_online(m, ob, lag = 0, init = init)$mean, filt$mean)
#> [1] TRUE
c(full_lag_vs_table_diagonal =
    max(abs(aci_metric(on_full, filt)$total - lt_diag(tb))),
  full_lag_vs_da_smooth = max(abs(on_full$mean - smoo$mean)))
#> full_lag_vs_table_diagonal      full_lag_vs_da_smooth 
#>               4.440892e-16               6.333986e-02
```

That last number is the scheme gap again, seen on the means instead of
on the metric. `aci_online(lag = Inf)` does not reproduce
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
and it is not meant to. The gap grows with the length of the record
rather than staying a fixed offset: on the packaged ENSO partition the
smoothed means differ by up to 1.89e-02 over 401 steps and 9.58e-02 over
4001, and the ACI values by 0.104 against a scale of 1.093 and 0.482
against 2.347 (development record).

`lag` has no default, deliberately: the lag is the argument the function
exists for, and a defaulted full-lag result is exactly the one that
would be mistaken for
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).
Cost is linear in the record whatever the lag, so choosing the lag is a
modelling decision and not a budget:

``` r

lags <- c(1, 20, 400, 800)
setNames(sapply(lags, function(L)
  system.time(aci_online(m, ob, lag = L, init = init))[["elapsed"]]),
  paste0("lag_", lags))
#>   lag_1  lag_20 lag_400 lag_800 
#>   0.013   0.025   0.023   0.002
```

An online path carries `kind = "online"`, and that is what keeps it out
of the places a complete smoother is required: `lag_table(smoother = )`
rejects it rather than treating a truncated lag as the whole record.

``` r

lag_table(m, ob, mode = "forward", init = init,
          filter = filt, smoother = on)
#> Error in `.validate_gaussian_path()`:
#> ! smoother has kind 'online', not 'smoother'.
```

## Evidence

The measurements cited above by name are in the repository’s development
record, `dev/acir-process-rationale.md`, which records each adopted
change with the measurement and the accuracy gate behind it. Numbers
computed in this vignette are computed on its own short records and are
not those measurements. The kinds of evidence behind the package’s
fidelity claims, and what each establishes, are set out in *Reproducing
the reference MATLAB codebase*.
