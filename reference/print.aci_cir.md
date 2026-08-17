# Print a causal influence range

Print a causal influence range

## Usage

``` r
# S3 method for class 'aci_cir'
print(x, ...)
```

## Arguments

- x:

  An object of class `aci_cir`.

- ...:

  Ignored.

## Value

`x`, invisibly.

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 900, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
               mu0 = model$y0, R0 = 0.1)
rng
#> Causal influence range
#>   101 reported times, 0.049 to 0.149; 129 thresholds
#>   status: resolved 101  censored 0  below_threshold 0  insufficient 0
#>   objective range: median 0.04139 (0.02497 to 0.05295)
```
