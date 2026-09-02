# Conditional ACI specification

Describes which observed channels carry the causal question and which
are conditioned out of it, and by which method. `aci_conditional()`
describes conditional ACI masking;
[`aci_conditional_reduce()`](https://biometryhub.github.io/ACI/reference/aci_conditional_reduce.md)
is the model reduction `method = "reduce"` asks for.

## Usage

``` r
aci_conditional(
  given = NULL,
  method = c("mask", "reduce"),
  target = NULL,
  first_step = c("uniform", "matlab")
)
```

## Arguments

- given:

  Integer or character vector naming the conditioning observed channels
  `x_B`. Supply this or `target`, not both.

- method:

  Either `"mask"`, which gives the conditioning channels' innovations
  zero weight in the filter, or `"reduce"`, which substitutes them as
  known forcing.

- target:

  Integer or character vector naming the target observed channels `x_A`.
  Supply this or `given`, not both.

- first_step:

  Either `"uniform"`, which masks the observation precision at every
  step, or `"matlab"`, which leaves the first slice unmasked as the
  reference scripts do. `"matlab"` requires `method = "mask"`.

## Value

An object of class `aci_conditional_spec`.

## Estimand

Split the observed process into target and conditioning channels,
`x = (x_A, x_B)`, with hidden process `y`. Conditional ACI is the
estimand `y(t) -> x_A | x_B`: only the target channels `x_A` transfer
information into the hidden posterior, and the conditioning channels
`x_B` are conditioned upon rather than assimilated. `"mask"` realises it
by giving the `x_B` innovations zero weight in the filter gain, the
Riccati term and the online-smoother gain, through an
observation-precision matrix supported only on the `A` block; `"reduce"`
realises it by rewriting the model so `x_B` enters as a known time
series (prescribed forcing). The two coincide when the `A`-`B` noise
cross-block vanishes, which the reduction checks along the whole path.

`target` names `x_A` directly, which is how the reference scripts write
the question (`h_W(t) -> T_C | (u, T_E, tau, I)`,
`ENSO_model_cond_ACI_h_W_unobs.m:1199-1202`). `given` names the
complement `x_B`. Supply exactly one; the other side is derived.

## First-slice convention

The reference scripts fill the first slice of the observation-precision
array with the full Gram inverse before the target-only overwrite, and
mask only the later slices (`ENSO_model_cond_ACI_h_W_unobs.m:1197`
against `:1250`). `first_step` selects between masking every slice,
`"uniform"`, and reproducing that asymmetry, `"matlab"`. It is not a
round-off-level choice: on a 4001-point ENSO path with `h_W` hidden and
`T_C` the target it moves the step-2 filter mean by 0.108 and the peak
ACI by 0.574, and the difference decays through the record rather than
vanishing. It is inert wherever the mask itself is inert. `"matlab"`
applies only to `"mask"`, which is where a masked precision path exists.

## See also

[`aci_conditional_reduce()`](https://biometryhub.github.io/ACI/reference/aci_conditional_reduce.md),
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)

## Examples

``` r
aci_conditional(given = 2, method = "mask")
#> <aci_conditional_spec> x_B = {2}, method = mask
aci_conditional(target = 1, method = "mask")
#> <aci_conditional_spec> x_A = {1}, method = mask
```
