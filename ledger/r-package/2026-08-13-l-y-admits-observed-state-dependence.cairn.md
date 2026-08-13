---
id: 2026-08-13-l-y-admits-observed-state-dependence
schema_version: 1.4
date: 2026-08-13
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "Latent self-drift widened to admit dependence on the observed state"
tags: [public_api, components_schema, predator_prey, cgns]
triggers:
  - public_api_change
  - adr_class_commitment
reversal_cost: medium
decision_pressure: null
review_due: null
review_trigger: "when the vector core lands: L_y becomes a matrix per step and the same question is asked again in a form where the answer may differ"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-online-smoother-log-reconstruction
---

## Status

Decided and implemented 2026-08-13 as stage R2 of the parity roadmap, to admit
the reference implementation's noisy predator-prey model. Shipped at version
0.1.0.9000. Until this change the numerical core integrated a self-drift
constant in time, which covers the dyad but not the Lotka-Volterra pair, where
the latent population's growth rate is set by the population being watched.

The widening is what the earlier `API_STABILITY.md` had promised would be
additive, and this cairn records that the promise was kept for the components
schema and *not* kept for one model field.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Admit only a length-n vector in the components schema; leave the model layer alone |
| B (CHOSEN) | Admit scalar or vector in the schema, and make `L_y` a coefficient function in the model layer like `L_x`, `f_x` and `f_y` already were |
| C | Add a separate time-varying constructor and leave the existing one untouched |
| D | Expand a scalar to a vector at construction so only one shape ever reaches the core |

## Rationale for rejection

### Option A

Would have left the model layer unable to express the predator-prey system at
all: `aci_simulate()` reads the self-drift once per step from the model, and a
components-only widening never reaches it. The result would have been a schema
that admitted a system the constructors could not build, which is a worse
contract than the constant-only one it replaced.

### Option C

Avoids touching a working constructor, which was its whole appeal, and was
rejected because it duplicates the constructor's validation, its noise
admissibility checks and its printed form for a difference of one argument's
type. The package already carries one lesson about two implementations of one
thing drifting apart: `aci_simulate_dyad()` was removed in 0.1.0 for exactly
that reason.

### Option D

Tempting because it makes the core see one shape. Rejected because expanding at
construction discards the information that the drift *is* constant, which the
printed form uses and which a future reader of a saved model would want. The
expansion happens instead at the top of each recursion, where it is a local
convenience rather than a change to what the object records.

## Implemented option (B)

The components schema takes `L_y` as a scalar or as one value per observation;
a scalar keeps its old meaning exactly, so every components list built against
the earlier contract still means what it meant. `aci_cgns_model()` takes `L_y`
on the same terms as the other coefficients -- a constant or a vectorised
function of the observed signal -- and stores it as a function.

The mathematics is unchanged in kind. A self-drift depending on the observed
path is still conditionally Gaussian, because the coefficient is measurable
with respect to that path. What remains impossible to express, deliberately, is
a coefficient depending on the *unobserved* component, which would leave the
class and the closed-form recursions with no posterior to propagate.

## Forward cost

- **One field's type changed, and that is a real break.** `model$L_y` was a
  number and is now a function. Code reading it as a number breaks loudly
  rather than silently, and `model$L_y_constant` carries the value when there
  is one. The components schema widened additively as promised; the model
  object did not, and saying so plainly is the point of recording it.
- **Three recursions now index a coefficient they used to hold constant** --
  the filter, the backward smoother and the online smoother's auxiliary
  quantities. A fourth, the causal influence range, inherits it through the
  online smoother. Any future recursion must index it too, and the failure
  would be a silently wrong answer rather than an error.
- **The expansion allocates.** A constant self-drift becomes a length-n vector
  at the top of each recursion. A few hundred kilobytes on the longest signals
  this core integrates, bought deliberately for a loop body that reads like the
  equation.
- **The same question returns, harder, with vector states.** `L_y` becomes a
  matrix per step, and whether it may be time-varying interacts with whether
  the ordered products of the online smoother still reduce. Recorded as this
  cairn's review trigger.

## References

### Methodological

- Andreou, M., Chen, N. and Bollt, E. (2026). *Assimilative causal inference*.
  Nature Communications, 17, 1854.
  https://doi.org/10.1038/s41467-026-68568-0 -- the noisy predator-prey model,
  Supplementary Information section 4.2.

### Empirical

- `aciR/R/aci-validate.R` -- the widened `L_y` clause of
  `.aci_check_components()` and the `.aci_expand()` helper.
- `aciR/R/predprey-model.R` -- `aci_predprey_components()` and
  `aci_predprey_model()`, both causal directions.
- `aciR/tests/testthat/test-oracle-predprey.R` -- the grade, and the
  non-degeneracy test asserting the self-drift genuinely moves and changes
  sign before the grade is trusted.
- Measured 2026-08-13 against `oracle/aci_oracle_predprey.m`: maximum absolute
  error 4.281020e-13 for direction predator-to-prey and 6.572520e-14 for
  prey-to-predator, against a gate of 1e-6. Self-drift ranges over
  [-0.3705, 0.8174] and [-1.1583, 0.8386] respectively.

### Operational

- Cairn: `2026-08-13-online-smoother-log-reconstruction.cairn.md` (this
  project) -- the online smoother whose auxiliary quantities this change
  reaches into.
- `aciR/API_STABILITY.md` -- the components contract, which recorded this
  widening as foreseeable before it happened.
