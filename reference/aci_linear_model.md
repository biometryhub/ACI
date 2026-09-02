# Conditionally linear model from constant coefficients

Convenience constructor for a conditional-Gaussian system whose coupling
and self-drift are constant, supplied either as scalars or as full
coefficient matrices.

## Usage

``` r
aci_linear_model(lambda_x, lambda_y, fx, fy, sigma_x, sigma_y, k = 1, l = 1)
```

## Arguments

- lambda_x:

  Coupling of the hidden state into the observed drift; a finite scalar
  or a `k` by `l` matrix.

- lambda_y:

  Hidden self-drift; a finite scalar or an `l` by `l` matrix.

- fx:

  Function of `(t, x)` giving the remaining observed drift.

- fy:

  Function of `(t, x)` giving the remaining hidden drift.

- sigma_x:

  Observed noise amplitude.

- sigma_y:

  Hidden noise amplitude.

- k:

  Observed dimension; a positive whole number.

- l:

  Hidden dimension; a positive whole number.

## Value

An object of class `cgns_model`.

## See also

[`aci_model()`](https://biometryhub.github.io/ACI/reference/aci_model.md)

## Examples

``` r
aci_linear_model(
  lambda_x = 1, lambda_y = -0.5,
  fx = function(t, x) -0.5 * x, fy = 0,
  sigma_x = 0.5, sigma_y = 1)
#> <aci_linear_model> 'conditionally_linear': k = 1 observed, l = 1 hidden
```
