# Simulate a stochastic or conditional-Gaussian model

Generates truth-twin realisations by Euler-Maruyama integration,
optionally retaining the hidden path and the driving Wiener increments
so that a later ensemble smoother can reuse them.

## Usage

``` r
# S3 method for class 'stochastic_model'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  t_end,
  dt,
  ic = NULL,
  burn_in = 0,
  keep_hidden = TRUE,
  keep_noise = TRUE,
  ...
)
```

## Arguments

- object:

  A `stochastic_model` or `cgns_model` object.

- nsim:

  Positive whole number of realisations to generate.

- seed:

  Optional non-negative whole number seeding the generator. Seeding is
  contained: the caller's `.Random.seed` is restored when the call
  returns, so a reproducible path leaves the caller's stream where it
  was.

- t_end:

  Positive 1-length numeric total simulated time, excluding burn-in.
  Called `T` before 0.1.0; that name is accepted with a warning until
  0.2.0.

- dt:

  Positive 1-length numeric integration step.

- ic:

  Optional list with elements `x0` and `y0` giving the initial state;
  `NULL` uses the model's default initial condition.

- burn_in:

  Non-negative 1-length numeric time discarded before recording.

- keep_hidden:

  `TRUE` to retain the hidden path.

- keep_noise:

  `TRUE` to retain the driving Wiener increments.

- ...:

  Must be empty; unused arguments are an error.

## Value

An object of class `aci_sim` when `nsim` is one, otherwise a list of
such objects. Its observed trajectory carries the model's own channel
names when `meta$vars$observed` supplies one unique non-empty name per
observed channel, as the built-in constructors do, so that
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)'s
`given` and `target` can be named straight off a simulation. The names
are labels: no numeric result depends on them.

## See also

[`aci_model()`](https://biometryhub.github.io/ACI/reference/aci_model.md),
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
