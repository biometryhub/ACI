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

## See also

[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
aci_metric_pair(mu_p = 0, R_p = matrix(1), mu_q = 1, R_q = matrix(2))
#>      total     signal dispersion 
#> 0.34657359 0.25000000 0.09657359 
```
