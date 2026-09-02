# Conditional-Gaussian nonlinear system

Constructs a conditional-Gaussian nonlinear system, in which the
observed drift is affine in the hidden state and the hidden drift is
linear in it, so that the conditional statistics of the hidden component
are Gaussian in closed form. The noise is described by up to two shared
Wiener channels, so that correlated observation and signal noise can be
represented.

## Usage

``` r
aci_model(
  Lx,
  fx,
  Ly,
  fy,
  Sx1,
  Sx2 = NULL,
  Sy1 = NULL,
  Sy2,
  k,
  l,
  name = NULL,
  meta = list()
)
```

## Arguments

- Lx:

  Function of `(t, x)` giving the coupling of the hidden state into the
  observed drift, a `k` by `l` matrix.

- fx:

  Function of `(t, x)` giving the remaining observed drift.

- Ly:

  Function of `(t, x)` giving the hidden self-drift, an `l` by `l`
  matrix.

- fy:

  Function of `(t, x)` giving the remaining hidden drift.

- Sx1:

  Function of `(t, x)` giving the observed diffusion on the first Wiener
  channel.

- Sx2:

  Optional function of `(t, x)` giving the observed diffusion on the
  second Wiener channel; `NULL` is a zero block matched to `Sy2`.

- Sy1:

  Optional function of `(t, x)` giving the hidden diffusion on the first
  Wiener channel; `NULL` is a zero block matched to `Sx1`.

- Sy2:

  Function of `(t, x)` giving the hidden diffusion on the second Wiener
  channel.

- k:

  Observed dimension; a positive whole number.

- l:

  Hidden dimension; a positive whole number.

- name:

  Optional 1-length character label for the model.

- meta:

  Optional named list of metadata carried on the object.

## Value

An object of class `cgns_model`, which also inherits from
`stochastic_model`.

## Details

CGNS coefficient functions are mathematical coefficients: for a fixed
`(t, x)` they must return deterministic values with stable shapes and
diffusion-channel counts. Random-number generation or result-changing
mutable state inside a coefficient function is outside the model
contract. Closed-form execution may realise each coefficient once on the
observation grid and reuse that realised path.

## See also

[`aci_model_from_affine()`](https://biometryhub.github.io/ACI/reference/aci_model_from_affine.md),
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
aci_model(
  Lx = function(t, x) matrix(1, 1, 1),
  fx = function(t, x) -0.5 * x,
  Ly = function(t, x) matrix(-0.5, 1, 1),
  fy = function(t, x) 0,
  Sx1 = function(t, x) matrix(0.5, 1, 1),
  Sy2 = function(t, x) matrix(1, 1, 1),
  k = 1, l = 1)
#> <cgns_model> 'cgns_model': k = 1 observed, l = 1 hidden
```
