# Gaussian relative entropy

Relative entropy of one multivariate Gaussian from another, optionally
split into its signal and dispersion parts.

## Usage

``` r
aci_metric_pair(mu_p, R_p, mu_q, R_q, decompose = TRUE)
```

## Arguments

- mu_p:

  Numeric vector, mean of the first distribution.

- R_p:

  Covariance matrix of the first distribution.

- mu_q:

  Numeric vector, mean of the second distribution.

- R_q:

  Covariance matrix of the second distribution.

- decompose:

  `TRUE` to return the signal and dispersion parts alongside the total.

## Value

A named numeric vector with the `total` and, when `decompose` is `TRUE`,
the `signal` and `dispersion` parts.

## Details

Public KL values never apply a covariance ridge; callers wanting
regularisation opt in explicitly with
[`spd_floor()`](https://biometryhub.github.io/ACI/reference/spd_floor.md).
andreou2026cir (Section 2.2, closing paragraph) states that
regularization is expected only in the degenerate limit, which is the
published basis for this function's strictness. The state recursions are
strict by default too, and regularise only when a call asks for it with
`regularize = "floor"`; see
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).

The dispersion part is `O(delta^2)` in `delta = R_p / R_q - 1`, and this
function evaluates it in the trace and log-difference form
`0.5 * (tr(R_q^-1 R_p) - l + log det R_q - log det R_p)` at every
dimension, including `l = 1`. That form loses relative precision to
cancellation as the two covariances approach each other: measured
against the exact series at `R_q = 1`, `R_p = 1 + delta`, its relative
error is 4.5e-05 at `delta = 1e-06`, 1.2e-02 at 1e-07 and a factor of
two from 1e-08 down.
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
instead routes a one-dimensional hidden state through a
cancellation-resistant `log1p` form, so the two public routes disagree
on identical one-dimensional inputs in that regime. The absolute error
stays below 3e-17 nats across that sweep, because the quantity itself is
quadratic in `delta`, so this is a precision boundary and not a
demonstrated accuracy failure on ordinary inputs. Positivity is imposed
by `max(., 0)` on each part; that is a guard, not an accuracy proof.

## See also

[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
aci_metric_pair(mu_p = 0, R_p = matrix(1), mu_q = 1, R_q = matrix(2))
#>      total     signal dispersion 
#> 0.34657359 0.25000000 0.09657359 
```
