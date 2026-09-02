# Six-variable stochastic conceptual ENSO model (andreou2026aci SI via chen2022enso eqs. 1a-1f).

Constructs the six-variable ENSO system split into an observed and a
hidden component at `hidden`. `I` must remain observed, and so must `TC`
unless the zeroth-order approximation below is requested. The
coefficients are the fixed-state ones of the ACI_code conditional ENSO
scripts.

## Usage

``` r
aci_enso_model(
  hidden = c("hW", "tau"),
  sigma_E = NULL,
  lambda = 2/60,
  params = NULL,
  variant = "aci_code",
  observations = c("reduced", "full"),
  approximation = c("exact", "zeroth_order_c1"),
  prescribed = NULL,
  matlab_defect_compat = FALSE
)
```

## Arguments

- hidden:

  Character vector naming hidden ENSO variables.

- sigma_E:

  Eastern-Pacific temperature noise amplitude.

- lambda:

  Decay rate for the diversity index.

- params:

  Optional complete parameter list overriding the preset.

- variant:

  Parameter and coefficient convention.

- observations:

  Observation set the estimand is defined on, `"reduced"` or `"full"`.
  `"reduced"` is defined for `hidden = "tau"` and `hidden = "TC"`, where
  it is the default and, for `"TC"`, the only value; every other
  partition is `"full"`. For `hidden = "tau"` observations are still
  supplied on all five observed channels either way, and `"reduced"`
  prescribes `u` and `h_W` from them instead of assimilating them.

- approximation:

  Either `"exact"`, the conditionally linear split of the six-state
  drift, or `"zeroth_order_c1"`, the `T_C`-hidden substitution
  `c1(t, TC) -> c1(t, 0)` described above. `"zeroth_order_c1"` is
  defined only for `hidden = "TC"`, and `hidden = "TC"` requires it.

- prescribed:

  For `hidden = "TC"` only, and then required: a data frame or named
  list carrying `t`, `u`, `hW` and `tau` as equal-length numeric series
  on one strictly increasing uniform time grid. Other elements, such as
  the remaining channels of a simulated path, are ignored.

- matlab_defect_compat:

  For `hidden = "TC"` only. The reference script's assimilation forcing
  `f_y` (`ENSO_model_cond_ACI_T_C_unobs.m:1053,:1151`) omits the
  thermocline term `gamma_C * h_W`, which the same script's simulator
  drift (`:1124`) includes and which the sibling `ACI_code` scripts
  carry in the corresponding `T_C` coefficient rows. `h_W` is prescribed
  and observed here, so the term is available, and `acir` includes it.
  Set `TRUE` to reproduce the published script verbatim. On a
  14-model-year window the two differ by up to `0.105` in the filter
  mean and `2.75` in ACI, and the time-integrated ACI roughly doubles;
  filter and smoother covariances are identical, because the term enters
  only the mean equations.

## Value

A `cgns_model` for the chosen ENSO partition, ready for
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and the
assimilation verbs.

## T_C hidden (zeroth-order)

The six-state ENSO system is not conditionally Gaussian with `T_C`
hidden: the damping `c_1(t, T_C) T_C` is cubic in `T_C`. `ACI_code`'s
`ENSO_model_cond_ACI_T_C_unobs.m` restores conditional linearity by a
zeroth-order Taylor expansion of `c_1` about the climatology `T_C = 0`,
replacing the state-dependent damping with the time-only series
`r_C - c_1(t, 0)`. `approximation = "zeroth_order_c1"` builds that
inference model: observed `(T_E, I)`, hidden `T_C`, with `u`, `h_W` and
`tau` entering as prescribed forcings from their observed series,
supplied through `prescribed`. This is an approximation of the system,
not a re-split of it - the simulator keeps the full nonlinear
`c_1(t, T_C)`, and the filter and smoother moments are those of the
approximating model.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) is refused on the
result for that reason; generate a path from
`aci_enso_model(hidden = c("u", "hW", "tau"))` and build this model from
it.

The prescribed series are looked up by index on the grid they were
supplied on, never interpolated, and assimilation refuses an observation
grid that is not that grid.

## Observation set for the tau partition

`hidden = "tau"` has two estimands, and they are not the same causal
quantity. `observations = "reduced"`, the default for that partition, is
the reference script's: the observed process is `(T_C, T_E, I)`, and `u`
and `h_W` enter the target drifts as prescribed known time series rather
than as assimilated channels
(`ENSO_model_cond_ACI_tau_unobs.m:1020-1039`, `:1136-1141`).
`observations = "full"` assimilates all five observed channels
`(u, h_W, T_C, T_E, I)` and so uses strictly more information: both
prescribed drifts carry `tau` (`Lx[u] = -0.0407`, `Lx[h_W] = -0.0814`),
and prescribing them reproduces their effect on the `T_C`/`T_E` drift
but not their own innovations.

The script asserts the two agree. On a 4001-point source-derived path
they do not: filter means differ by up to 0.247, the ACI series by up to
0.776 - about three times its own mean level, with Pearson correlation
0.905 - while the time-averaged ACI agrees to within 0.5%. The reduction
preserves the average level and distorts the time-resolved curve, which
is the reported quantity. Any fidelity claim must name which observation
set it reproduces. `meta$observations` records which one a model
carries.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications* **17**, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

Chen, N., Fang, X. and Yu, J.-Y. (2022). A multiscale model for El Nino
complexity. *npj Climate and Atmospheric Science* **5**, 16.
arXiv:2104.07174.

## Examples

``` r
aci_enso_model()
#> <cgns_model> 'ENSO6[aci_code] (hW,tau hidden)': k = 4 observed, l = 2 hidden

# The T_C-hidden partition is not self-contained: it needs the u, h_W and
# tau series it treats as prescribed forcings.
sim <- simulate(aci_enso_model(hidden = c("u", "hW", "tau")),
                seed = 1, t_end = 1, dt = 0.005)
path <- data.frame(t = sim$obs$t, u = sim$hidden[, 1],
                   hW = sim$hidden[, 2], tau = sim$hidden[, 3])
aci_enso_model(hidden = "TC", approximation = "zeroth_order_c1",
            prescribed = path)
#> <cgns_model> 'ENSO6[aci_code] (TC hidden, zeroth-order c1)': k = 2 observed, l = 1 hidden
```
