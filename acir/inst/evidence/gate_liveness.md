# Gate liveness

A gate is a check that is supposed to stop something. The rule this file
records is that every gate must be shown to fail once, on a violation induced
on purpose, and that run must be written down. A gate with no such record is
not a gate: it may never have executed, or it may be comparing a quantity
against itself.

One row per gate. `Where enforced` is the code the gate lives in and the test
that exercises it on clean input. `Induced violation` is the deliberate
breakage. `Where recorded` is the block that performs it.

Every induced violation is made on a copy. No shipped fixture, manifest or
register is modified by any test in this package.

## Evidence gates

| Gate | Where enforced | Induced violation | Where recorded |
|---|---|---|---|
| Fixture byte-pinning, authors-source and independent-transcription set | `tests/testthat/fixtures/oracles/oracle-manifest.yml`, compared file by file in `tests/testthat/test-19-compiled-oracles.R` over 16 fixtures | one decimal digit of a copy of `dyad_reference.csv` moved by one, at unchanged file length | `tests/testthat/test-31-gate-liveness.R`, "byte-pinning fails on one flipped byte of a fixture copy" |
| Fixture byte-pinning, source-derived partition set | `tests/testthat/fixtures/oracles/oracle-manifest-partitions.yml`, compared file by file in `tests/testthat/test-28-partition-oracles.R` over 14 fixtures | the same flipped-byte copy; the two manifests are read by one shared reader | `tests/testthat/test-31-gate-liveness.R`, same block |
| Manifest reader coverage: no manifest entry may be invisible to the reader | `tests/testthat/helper-oracle-manifest.R`, asserted in `tests/testthat/test-19-compiled-oracles.R` by requiring the strict and permissive readers to agree | the earlier name class `^[A-Za-z_]+\.csv:$` applied to the current partition manifest, where it returns none of the six `enso6_partition_*` entries the shipped reader returns | `tests/testthat/test-31-gate-liveness.R`, "the shipped manifest reader sees the names the narrow one went blind to" |
| Evidence register coverage: every export has a row, every checked row has a fixture whose bytes match, the column set is exact | `inst/evidence/register.csv`, enforced in `tests/testthat/test-30-evidence-register.R` | three, on copies: one verb's rows removed from the table, one hex digit of a stored sha256 moved, and a column renamed through a written file and read back | `tests/testthat/test-31-gate-liveness.R`, "register coverage fails on a dropped row, a moved hash and a renamed column" |
| Numerical oracle tolerance, dyad grade at 1e-6 | `tests/testthat/test-19-compiled-oracles.R` | each of the five graded columns of `dyad_reference.csv` shifted by 1e-5 in turn, ten times the tolerance, with the comparison required to see every one | `tests/testthat/test-31-gate-liveness.R`, "the dyad oracle tolerance rejects a perturbation in every graded column" |

The remaining numerical grades use the same comparison shape as the dyad grade
against their own references: the online auxiliaries and the matrix paths at
1e-6, the scalar ENSO partitions at 1e-8, the streamed forward range at 1e-10
for the peak and objective and 1e-12 for the subjective read-out, and the T_C
zeroth-order chain at 1e-12. The column-by-column induction above is performed
on the dyad grade only, and is not claimed for the others.

## Policy and authentication gates

| Gate | Where enforced | Induced violation | Where recorded |
|---|---|---|---|
| Strict covariance policy: a covariance that leaves the positive-definite cone stops the run and names the site and index | `.cov_guard()` in `R/aci-utils.R`, called from the scalar and matrix recursions | three unstable configurations that make the guard fire: a dyad grid coarsened 200 times, a scalar prior covariance of 1e-14, and the same prior on a two-hidden matrix model | `tests/testthat/test-26-covariance-policy.R`, blocks S4, S5 and S5b. The same three configurations were run as standalone probes with their firing site, index, time and value recorded in the adoption ledger at entry C3d |
| Non-finite covariance is refused under both policies | `.cov_guard()` and `spd_floor()`'s finiteness contract, both in `R/aci-utils.R` | a covariance matrix carrying `NA` | `tests/testthat/test-26-covariance-policy.R`, "a non-finite covariance stops under both policies" |
| Cholesky factorisation of a non-symmetric-positive-definite matrix aborts rather than returning | `safe_chol()` in `R/aci-utils.R` | an indefinite two-by-two, under policy `none` and policy `floor` | `tests/testthat/test-26-covariance-policy.R`, "the exported regularisers keep their documented behaviour" |
| Batch realiser authentication: the fast affine route is taken only for a model whose constructor identity is intact | `.attach_cgns_realizer()` and the realiser selection in `R/aci-model.R` | a coefficient closure replaced by hand on an otherwise valid model, which then authenticates as nothing and falls back to the generic route | `tests/testthat/test-25-enso6-batch-realiser.R`, "ENSO6 realiser selection uses sealed constructor identity" |
| Trusted-filter authentication: a supplied filter is reused without revalidation only when its token is genuine | the token check in `R/aci-assimilation.R` | four, in separate blocks: a hand-mutated covariance, a hand-mutated model and provenance, a foreign and a forged token, and a filter built for a different observation grid | `tests/testthat/test-22-trusted-filter.R`, four blocks named for those four cases |

## Contract and argument gates

| Gate | Where enforced | Induced violation | Where recorded |
|---|---|---|---|
| CGNS model contract: dimensions, symmetry, callable coefficients, unknown arguments | `aci_model()` and `validate_cgns()` in `R/aci-model.R` | planted violations of each kind, refused with `aci_error_model_contract` or `aci_error_dims` | `tests/testthat/test-02-models.R` |
| Observation contract: numeric, finite, strictly increasing time | `observed_trajectory()` and `as_obs()` in `R/aci-utils.R` | non-numeric, non-finite and non-increasing inputs, refused with `aci_error_obs_contract` | `tests/testthat/test-02-models.R`, with a further case in `tests/testthat/test-15-compiled-scalar.R` |
| Conditional admissibility: a split needs exactly one non-empty side, indices must be whole and distinct, and a model that declares its own estimand refuses a second specification | `aci_conditional()` in `R/aci-assimilation.R` and the conditioning entry points | empty, fractional and duplicate index sets, and a second specification supplied over a declared reduced estimand | `tests/testthat/test-05-nontarget.R`, `tests/testthat/test-24-conditional-estimands.R`, `tests/testthat/test-28-partition-oracles.R` |
| Online lag arguments: stepper, initial-condition dimension | `aci_online()` in `R/aci-assimilation.R` | an implicit-stepper filter and a wrongly sized initial covariance, refused with `aci_error_stepper` and `aci_error_dims` | `tests/testthat/test-27-online-lag.R` |
| `aci_range()` argument gates: unknown argument, out-of-range anchors, even-length threshold grid, incompatible method and quadrature, unknown read-out convention | the argument validation in `R/aci-core.R` | one deliberately invalid value for each | `tests/testthat/test-03-engine.R`, `tests/testthat/test-23-forward-cir-api.R` |
| `aci_range()` direction gate: the backward range is refused rather than answered forward | `.aci_range_direction()` in `R/aci-core.R`, called by both methods | `direction = "backward"` on the `lag_table` method and on the `aci_result` method, and an unknown direction | `tests/testthat/test-31-gate-liveness.R`, "the backward range is refused on both aci_range methods" |
| Non-CGNS model routes are refused rather than dispatched | the `stochastic_model` methods in `R/aci-assimilation.R` | a bare `stochastic_model` object passed to `aci_filter()`, `aci_smoother()`, `aci_online()` and `lag_table()` | `tests/testthat/test-27-online-lag.R` for `aci_online()`; `tests/testthat/test-31-gate-liveness.R`, "a non-CGNS model is refused by every assimilation verb", for the other three |
| Simulator divergence guard: a trajectory that leaves the double range stops and reports the step | the Euler-Maruyama loop in `R/aci-model.R` | a valid CGNS model with a hidden self-drift of 1e6 on a 1e-3 grid, which overflows inside 110 steps; the same model on a 1e-6 grid returns finite output, so the abort is the guard and not the constructor | `tests/testthat/test-31-gate-liveness.R`, "the divergence guard stops an unstable simulation and names the step" |

## What is not listed, and why

`aci_error_internal` is raised in three places by option-gated debug assertions
on the sign of a relative entropy, and in one place by a check on an internal
realiser descriptor. None of them can be reached by any caller-visible input,
so they are internal invariants rather than gates, and no induced run is
claimed for them.

The register above covers gates that reject caller input or shipped evidence.
It does not cover warnings, which report rather than stop, and it does not
cover the masking conventions of `aci_range()`, which return a value.
