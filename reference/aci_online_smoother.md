# Fixed-lag online conditional Gaussian smoother

Runs the forward-in-time online smoother of a conditional Gaussian
nonlinear system. The online smoother estimates the unobserved component
at step `j` from the observed path up to step `j + lag`, so it
interpolates between the forward filter, which uses no future
observations, and the backward smoother, which uses all of them.

## Usage

``` r
aci_online_smoother(x, comp, dt, filt, lag = Inf, tol = 1e-18)
```

## Arguments

- x:

  Numeric vector. The observed signal, one value per time step; at least
  two complete, finite observations.

- comp:

  A conditional Gaussian components list; see
  [aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).

- dt:

  Numeric scalar. The integration time step; must be positive.

- filt:

  A list with numeric vectors `mean` and `cov`, as returned by
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).

- lag:

  Numeric scalar. The number of future steps each estimate may condition
  on. Must be a non-negative whole number, or `Inf` for the full path.
  Defaults to `Inf`.

- tol:

  Numeric scalar. Relative magnitude below which the geometric update
  products are truncated; must be positive. Defaults to `1e-18`.

## Value

A list with numeric vectors `mean` and `cov`, the online smoothed mean
and covariance of the unobserved component at each time step, and the
integer `lag_effective`, the longest update range actually accumulated
before truncation.

## Details

The `lag` argument selects a member of that family, and its two
boundaries are exact identities rather than approximations. At `lag = 0`
the estimator reproduces
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md);
at `lag = Inf`, or any lag at least as long as the signal, it reproduces
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).
Intermediate lags are the online estimator proper, and they are what the
causal influence range is built from.

Each new observation updates every retained earlier step through an
ordered product of the per-step auxiliary matrices. That product decays
geometrically, which the source paper establishes by bounding the
spectral radius of each factor below one, so the influence of an
observation on distant earlier steps falls away exponentially. The
implementation exploits this by accumulating the products in logarithms
rather than forming the full `O(n^2)` triangle, which is what keeps the
memory at `O(n)`.

`tol` stops the accumulation once the product falls below it. It is a
safety catch on a product that has already underflowed rather than a
bound on the work. At a typical per-step contraction of `0.998` the
product needs some 26,000 steps to reach the default `1e-18`, and
published records are shorter than that, so on the systems this package
ships the loop is not cut and the cost is quadratic in the record
length. The test applied is `max(abs(d)) < tol` on the accumulated
product, a sufficient condition that the next innovation cannot move the
estimate at double precision, rather than the spectral-radius argument
above. Erring that way costs work, never accuracy.

## References

Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother
with closed-form solutions and information-theoretic lag selection for
conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*,
36(4), 71.
[doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x)

## See also

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 500, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)

# A lag of zero is the filter; a full lag is the backward smoother.
online <- aci_online_smoother(sim$x, comp, dt = 0.001, filt, lag = 50)
str(online)
#> List of 3
#>  $ mean         : num [1:500] 2.05 2.05 2.04 2.04 2.03 ...
#>  $ cov          : num [1:500] 0.0908 0.0913 0.0918 0.0923 0.0929 ...
#>  $ lag_effective: int 50
```
