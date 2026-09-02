# Project to SPD by eigenvalue flooring (used after Euler covariance updates)

Returns the matrix unchanged when it already admits a Cholesky factor,
and otherwise floors its eigenvalues at a small positive multiple of the
largest absolute eigenvalue before reassembling it.

## Usage

``` r
spd_floor(R, eps = 1e-12)
```

## Arguments

- R:

  Symmetric matrix.

- eps:

  Minimum eigenvalue.

## Value

A symmetric positive-definite matrix of the same dimension as `R`.

## See also

[`safe_chol()`](https://biometryhub.github.io/ACI/reference/safe_chol.md)

## Examples

``` r
spd_floor(matrix(c(1, 2, 2, 1), 2, 2))
#>      [,1] [,2]
#> [1,]  1.5  1.5
#> [2,]  1.5  1.5
```
