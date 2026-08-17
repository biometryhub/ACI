# Causal influence range

Computes the subjective and objective causal influence range of the
unobserved component on the observed signal. The causal influence range
measures how far forward in the observed record one must look before the
estimate of the unobserved state at a given time stops improving, and so
complements the causal-information metric, which measures how much it
improves in total.

## Usage

``` r
aci_cir(
  x,
  comp,
  dt,
  filt = NULL,
  window = NULL,
  epsilon = 10^seq(-6, 0.5, length.out = 129L),
  threshold = 1e-05,
  margin = 0.1,
  horizon = NULL,
  mu0 = NULL,
  R0 = NULL
)
```

## Arguments

- x:

  Numeric vector. The observed signal, one value per time step.

- comp:

  A conditional Gaussian components list; see
  [aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).

- dt:

  Numeric scalar. The integration time step; must be positive.

- filt:

  A list with numeric vectors `mean` and `cov`, as returned by
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).
  When `NULL`, the filter is run internally, which requires `mu0` and
  `R0`.

- window:

  Integer vector or `NULL`. Indices of the time steps at which the range
  is reported. When `NULL`, the whole signal is used, which is quadratic
  in its length and is rarely what is wanted for a long record.

- epsilon:

  Numeric vector. Thresholds at which the subjective range is evaluated,
  in nats. Defaults to a 129-point logarithmic grid from `1e-6` to
  `10^0.5`. The reference implementation spans the same range with 513
  points; the cheaper grid is this package's default, and
  `10^seq(-6, 0.5, length.out = 513L)` reproduces the reference's
  `objective_exact`.

- threshold:

  Numeric scalar. A peak divergence below this value is treated as no
  detectable influence, and the objective range is reported as zero
  rather than as the ratio of two negligible quantities. Defaults to
  `1e-5`.

- margin:

  Numeric scalar in `(0, 1)`. The fraction of each comparison sequence
  that must remain unused for the range at that time to count as
  resolved. A time whose range consumes more than `1 - margin` of the
  sequence it was measured against is marked `status = "censored"`; its
  ranges are still returned, as lower bounds. Defaults to `0.1`.

  Because a censored time yields a bound rather than a hole, `margin`
  governs the flag rather than whether a number is reported at all, so a
  slightly wrong value has little consequence. It is this package's
  device. The reference guards the same problem with an absolute
  lookahead chosen for one figure.

- horizon:

  Integer scalar or `NULL`. How many steps of the record each reported
  time may look forward across, counted from the start of the record
  rather than from the reported time. `NULL`, the default, uses the
  whole record.

  The reference implementation truncates this comparison at the end of
  its reporting window, which biases both the integral and the range low
  near that end; supplying the same value here reproduces its numbers.
  The truncation applies only to how far forward the comparison looks,
  never to the fully informed posterior it is compared against, which is
  always taken over the whole record.

- mu0, R0:

  Numeric scalars. Initial filtered mean and covariance, used only when
  `filt` is `NULL`.

## Value

An object of class `aci_cir`, a list with the reported `time`, the
`objective` range at each time, the `objective_exact` range obtained by
integrating the subjective ranges over the whole threshold grid rather
than by the efficient approximation, the `subjective` range as a matrix
with one row per threshold and one column per time, the `peak`
divergence at each time, the logical matrix `subjective_censored`
marking thresholds whose range ran past the retained margin, the
character `status` (`"resolved"`, `"censored"`, `"below_threshold"` or
`"insufficient"`), the logical `monotone` marking times whose divergence
sequence decreases with lag (the condition under which `objective` and
`objective_exact` are the same functional), and the logical `saturated`,
which is `status == "censored"`.

## Details

At each reported time the function forms the relative entropy between
the online-smoother posterior informed by the whole record and the
posterior informed only up to each later observation. That sequence
decays as the later observation advances. The **subjective** range at a
threshold is the elapsed time after which the sequence stays below that
threshold; the **objective** range is the threshold-free summary
obtained by integrating the sequence and normalising by its peak, which
the source paper gives as a computationally efficient underestimate of
the range defined by averaging the subjective ranges over all
thresholds.

The computation is quadratic in the length of `window`, because every
reported time is compared against every later observation. Choose the
window accordingly. A few thousand steps is comfortable, and reproducing
a figure at the scale of the published one is a batch computation rather
than an interactive one.

A reported time close to the end of the record cannot be resolved,
because the observations that would settle it do not exist. Such a time
is not unmeasured. Its range is **right-censored**, and the truncated
value is a lower bound on the true one. The result therefore returns the
value and marks that time `"censored"` in `status`, rather than
discarding what the record does support. Only a time with fewer than
three later observations, where no quadrature is possible at all,
returns `NA`.

`objective` and `objective_exact` are different functionals, not two
quadratures of one. They coincide when the divergence decreases with
lag, by the layer-cake identity. Both this package and the reference
measure the range as the **last** time the divergence exceeds a
threshold, which is at least the measure of the superlevel set and is
strictly larger the moment the sequence is not monotone. A sequence that
is not monotone is the common case on a truncated `horizon`. Seeing
`objective` below `objective_exact` is that definitional gap, not a
numerical defect in either.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother
with closed-form solutions and information-theoretic lag selection for
conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*,
36(4), 71.
[doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x)

## See also

[`aci_online_smoother()`](https://biometryhub.github.io/ACI/reference/aci_online_smoother.md),
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 600, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:250,
               mu0 = model$y0, R0 = 0.1)
summary(rng$objective)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02162 0.03528 0.04530 0.04636 0.05477 0.08105 

# `horizon` is an index into the record, not a lookahead from each reported
# time, so reproducing a published panel means converting from a time. This
# is the conversion the reference performs to obtain its `last_idx`.
dt <- 0.001
time_end <- 0.25
lookahead <- 0.05
last_idx <- as.integer(round((time_end + lookahead) / dt)) + 1L
aci_cir(sim$x, comp, dt = dt, window = 50:250, mu0 = model$y0, R0 = 0.1,
        horizon = last_idx)
#> Causal influence range
#>   201 reported times, 0.049 to 0.249; 129 thresholds
#>   status: resolved 0  censored 201  below_threshold 0  insufficient 0
#>   objective range: median 0.03335 (0.01451 to 0.0592)
#>   201 times are censored; those ranges are lower bounds
```
