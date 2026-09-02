# CGNS model from an affine-in-hidden drift pair

Recovers the conditional-Gaussian coefficients from drift functions
written in terms of the full state, by exact affine differencing. The
supplied functions are checked before they are replaced by the
reconstructed representation, so a drift that is not affine in the
hidden state is rejected rather than silently linearised.

## Usage

``` r
aci_model_from_affine(
  f_full,
  g_full,
  Sx,
  Sy_hidden,
  k,
  l,
  name = NULL,
  Sx2 = NULL,
  Sy_shared = NULL,
  meta = list()
)
```

## Arguments

- f_full:

  Function of `(t, x, y)` giving the full observed drift.

- g_full:

  Function of `(t, x, y)` giving the full hidden drift.

- Sx:

  Function of `(t, x)` giving the observed diffusion on the first Wiener
  channel.

- Sy_hidden:

  Function of `(t, x)` giving the hidden diffusion on its own Wiener
  channel.

- k:

  Observed dimension; a positive whole number.

- l:

  Hidden dimension; a positive whole number.

- name:

  Optional 1-length character label for the model.

- Sx2:

  Optional function of `(t, x)` giving the observed diffusion on the
  second Wiener channel.

- Sy_shared:

  Optional function of `(t, x)` giving the hidden diffusion on the
  shared Wiener channel, which introduces correlated noise.

- meta:

  Optional named list of metadata carried on the object.

## Value

An object of class `cgns_model`.

## See also

[`aci_model()`](https://biometryhub.github.io/ACI/reference/aci_model.md)

## Examples

``` r
aci_model_from_affine(
  f_full = function(t, x, y) -0.5 * x + y,
  g_full = function(t, x, y) -0.5 * y,
  Sx = function(t, x) matrix(0.5, 1, 1),
  Sy_hidden = function(t, x) matrix(1, 1, 1),
  k = 1, l = 1)
#> <cgns_model> 'cgns_model': k = 1 observed, l = 1 hidden
```
