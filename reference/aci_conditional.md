# Condition the causal question on a subset of the observed components

Rewrites a vector-valued components list so that the causal-information
metric computed from it measures what the *target* observed components
say about the unobserved state, rather than what the whole observed
process says.

## Usage

``` r
aci_conditional(comp, target)
```

## Arguments

- comp:

  A vector-valued conditional Gaussian components list; see
  [aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).
  The scalar schema has only one observed component and so admits no
  conditional question.

- target:

  Integer or character vector. Which observed components the causal
  question is asked about. Characters are matched against the row names
  of the observation-noise Grammian.

## Value

A components list with `S_xoS_x_inv` supported on the target block,
ready for
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
and
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md).

## Details

The non-target components are not removed. They continue to drive the
system, and their own dynamics are unchanged; what changes is that the
filter stops treating them as observations to be assimilated. This is
what distinguishes the conditional question from simply running the
method on a shorter signal, in which the non-target components would not
be present at all.

The construction inflates the observational uncertainty of the
non-target components without bound, which sends their weight in the
filter to zero. It is implemented by supplying the filter with an
inverse noise Grammian supported only on the target block. That object
is not the inverse of a covariance, and it is not required to be.

A component whose observation noise is exactly zero cannot be a target,
because the construction needs its noise covariance to be invertible on
the target block. The reference implementation meets this by adding a
small artificial noise to such a component, which is a modelling
decision and is reported here as an error rather than made silently.

## Estimand

The metric is the relative entropy \\KL(p(y_t \mid x\_{0:T}, M) \\ p(y_t
\mid x\_{0:t}, M))\\, a property of a pair of model-conditional
posteriors. It is not \\p(y \mid do(x))\\, not a Granger test, and not
transfer entropy, although all four can be large at the same instants on
a conditional Gaussian system.

Conditioning here changes the **observation the filter assimilates**,
not the system that generated the path. Masking a block of the inverse
observation-noise Grammian means "stop updating on this channel"; it
does not mean the channel was absent, and it is not an adjustment in the
do-calculus sense. The non-target components still drive the drift,
which is the entire point of the construction and the point most easily
misread by a reader arriving from the interventional literature.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md),
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)

## Examples

``` r
# A two-component observed process; ask what only the first says about the
# unobserved state.
comp <- list(
  L_x = diag(2), f_x = c(0, 0), L_y = -diag(2), f_y = c(0, 0),
  S_xoS_x = diag(c(0.5, 0.8)), S_yoS_y = diag(2),
  S_yoS_x = matrix(0, 2, 2)
)
conditioned <- aci_conditional(comp, target = 1)
conditioned$S_xoS_x_inv
#>      [,1] [,2]
#> [1,]    2    0
#> [2,]    0    0
```
