# Cholesky with escalating jitter ladder

Returns the Cholesky factor of a square numeric matrix, adding an
escalating ridge to the diagonal when the unmodified matrix is not
numerically positive definite. The ladder runs from `1e-12` to `1e-6`
relative to the mean diagonal entry.

## Usage

``` r
safe_chol(R, where = "covariance")
```

## Arguments

- R:

  Symmetric covariance matrix.

- where:

  Context label used in numerical error messages.

## Value

Upper-triangular Cholesky factor of `R`, possibly of `R` plus a diagonal
ridge.

## See also

[`spd_floor()`](https://biometryhub.github.io/ACI/reference/spd_floor.md)

## Examples

``` r
safe_chol(matrix(c(2, 0.5, 0.5, 1), 2, 2))
#>          [,1]      [,2]
#> [1,] 1.414214 0.3535534
#> [2,] 0.000000 0.9354143
```
