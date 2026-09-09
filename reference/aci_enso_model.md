# Six-variable stochastic conceptual ENSO model (andreou2026aci SI via chen2022enso eqs. 1a-1f).

Constructs the six-variable ENSO system split into an observed and a
hidden component at `hidden`. `I` must remain observed, and so must `TC`
unless the zeroth-order approximation below is requested. The
coefficients are the fixed-state ones of the ACI_code conditional ENSO
scripts. The constructor's own default partition is
`hidden = c("hW", "tau")`; a partition whose hidden components differ
widely in stationary variance needs a stated prior, as the next section
and the examples show.

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

## Prior for a multi-state partition

When `init` is omitted the assimilation verbs supply one scalar times
the identity, derived from the mean over hidden components of the
hidden-noise Gram and of the hidden damping. On
`hidden = c("u", "hW", "tau")` that mean is dominated by `tau`: the
componentwise Ornstein-Uhlenbeck stationary scales
`gyy_ii / (2 |Ly_ii|)` are `(0.0032, 0.0008, 0.283697)` while the
automatic scalar is `1.221602`, about 380 times too wide on `u` and
`h_W`. The explicit single-step Riccati update then leaves the
positive-definite cone at the first update: on the record of the
examples below, `aci(m, ob)` aborts `aci_error_covariance_not_spd` at
index 2, time 0.005, smallest eigenvalue -7.149360359. The bare
constructor's two-state default refuses the same way on its own record,
at -8.265707365.

That refusal is hidden-covariance integration instability on a valid
observation model, not a rank-deficient observation Gram: on the same
801 points the observation Gram has smallest eigenvalue 0.000325 and
smallest reciprocal condition number 0.00244 throughout. Because the
failure is at the first update there is no opening window whose
discarding would repair it. State a prior instead - the examples compute
the componentwise stationary scale from the model's own coefficients -
or refine the integration.

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

# Assimilating the three hidden states u, h_W and tau, with the prior
# stated componentwise: gyy_ii / (2 |Ly_ii|) is the Ornstein-Uhlenbeck
# stationary variance of hidden component i, read off the model's own
# coefficients at the first observation. Default explicit stepper,
# nsub = 1.
m   <- aci_enso_model(hidden = c("u", "hW", "tau"))
ob  <- as_obs(simulate(m, seed = 12, t_end = 4, dt = 0.005, burn_in = 0))
d0  <- function(f) diag(as.matrix(f(ob$t[1], ob$x[1, ])))
gyy <- d0(function(t, x) tcrossprod(m$Sy1(t, x)) + tcrossprod(m$Sy2(t, x)))
Ly  <- d0(m$Ly)
ini <- list(mean = m$meta$ic_default$y0, cov = diag(gyy / (2 * abs(Ly))))
diag(ini$cov)                      # 0.0032, 0.0008, 0.283697
#> [1] 0.0032000 0.0008000 0.2836974
a <- aci(m, ob, init = ini)
max(a$aci)                         # 1.830865
#> [1] 1.830865
a$meta$regularization[c("policy", "fired", "n_events")]
#> $policy
#> [1] "none"
#> 
#> $fired
#> [1] FALSE
#> 
#> $n_events
#> [1] 0
#> 

# Omitting init supplies the scalar automatic prior, which is far too wide
# on u and h_W, and the recursion is refused at the first update.
e <- tryCatch(suppressWarnings(aci(m, ob)), error = function(e) e)
class(e)[1]
#> [1] "aci_error_covariance_not_spd"
c(index = e$index, time = e$time, value = e$value)
#>    index     time    value 
#>  2.00000  0.00500 -7.14936 

# \donttest{
# Checks on the configuration above, not alternatives to it. Refining the
# step under the stated prior stays event-free and settles near 1.81. The
# explicit and implicit schemes are not required to agree at a finite
# step; this is a self-consistency check on one record, not an accuracy
# guarantee.
max(aci(m, ob, init = ini, nsub = 20)$aci)                    # 1.813052
#> [1] 1.813052
max(aci(m, ob, init = ini, stepper = "implicit", nsub = 20)$aci)
#> [1] 1.803244

# Refining enough to make the automatic-prior record run at all does not
# recover the same number, which is why the prior is stated rather than
# worked around.
max(suppressWarnings(aci(m, ob, nsub = 20))$aci)              # 4.821106
#> [1] 4.821106
max(suppressWarnings(aci(m, ob, stepper = "implicit"))$aci)   # 4.821235
#> [1] 4.821235

# Sensitivity of the stated prior: 0.5x to 4x run event-free, 8x is
# refused at index 2, and the automatic prior is about 380x.
vapply(c(0.5, 1, 2, 4), function(s) {
  max(aci(m, ob, init = list(mean = ini$mean, cov = s * ini$cov))$aci)
}, numeric(1))
#> [1] 1.814876 1.830865 1.851048 1.971858
# }
```
