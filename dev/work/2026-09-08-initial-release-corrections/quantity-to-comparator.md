# Quantity to comparator, `inst/evidence/register.csv`

Evidence for E02, 2026-09-08. One entry per `state == "checked"` row of
`acir/inst/evidence/register.csv` (28 of the register's rows). Line numbers
are the CSV's own file lines, so line 2 is the first data row. Test line ranges
are the `test_that()` block boundaries in the tree of this round; where the
assertion lives in a shared helper, the helper's line is named as well.

Register line numbers and test line numbers below are as of the register that
carried 65 rows, after the row-2 and row-36 edits described at the end and
before the three behavioural rows this round added at the end of the file. Row
36 is `behavioural_only` and therefore outside this table.

Column list is the fixture's header as shipped.

Tolerance abbreviations: `t19 = .compiled_oracle_tolerance` (1e-6,
`test-19-compiled-oracles.R:24`); `t28 = .partition_oracle_tolerance` (1e-8,
`test-28-partition-oracles.R:53`); `t29 = .tc_oracle_tolerance` (1e-12,
`test-29-tc-zeroth-order.R:57`).

---

## Rows resting on `dyad_reference.csv`

Columns: `t, x, y, filter_mean, filter_cov, smoother_mean, smoother_cov,
ACI_metric` (301 rows).
Comparator helper: `.compiled_oracle_scalar_error()`,
`test-19-compiled-oracles.R:58-66`, five `abs()` terms, one per fixture output
column, at `idx = seq(1, 30001, by = 100)`.

| line | export | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|
| 2 | `aci` | `test-19:188-205`, assertion at `:64` | t19 = 1e-6 | authors_source_comparison | `run$metric$total[idx]` vs `ref$ACI_metric` | MISMATCH before this round; corrected |
| 6 | `aci_filter` | `test-19:188-205`, assertions at `:60-61` | t19 = 1e-6 | authors_source_comparison | filter mean and covariance | OK |
| 11 | `aci_smoother` | `test-19:188-205`, assertions at `:62-63` | t19 = 1e-6 | authors_source_comparison | smoother mean and covariance | OK |
| 41 | `aci_dyad_model` | `test-19:188-205` | t19 = 1e-6 | authors_source_comparison | the four moment columns and `ACI_metric` produced by the constructed model | OK, indirect and disclosed (see note A) |

Row 2, before the edit, read feature "ACI metric with its signal and dispersion
parts, scalar hidden variable" and against "ACI_code dyad_interaction_model.m
stored outputs". The fixture has no signal or dispersion column, and the only
metric term in the helper is `abs(run$metric$total[idx] - ref$ACI_metric)`
(`test-19:64`). `decompose = TRUE` at `test-19:53` computes the parts but
nothing compares them against an authors-source column. The row now reads
"total ACI metric, scalar hidden variable" and its against text names the
column, the restriction to the total, and `test-19-compiled-oracles.R:64`.
The part-level evidence is register rows 4, 10, 14 (independent transcription,
1e-12, whole record) and row 39 (`aci_metric_pair`, independent transcription,
behavioural_only); row 39 is unchanged.

## Rows resting on `cross_reference.csv`

Columns: `t, x, y, filter_mean, filter_cov, smoother_mean, smoother_cov,
ACI_metric` (301 rows). Same helper as above; the model is
`.compiled_oracle_cross_model()` (`test-19:68-81`), with
`S_yo S_x = 0.5 * 0.6 + 0.8 * 0.3 = 0.54`.

| line | export | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|
| 12 | `aci_smoother` | `test-19:208-225` | t19 = 1e-6 | independent_transcription | smoother mean and covariance under non-zero noise cross-covariance | OK |
| 50 | `aci_model` | `test-19:208-225` | t19 = 1e-6 | independent_transcription | the generic model contract carrying all five output columns | OK |

## Rows resting on the predator-prey fixtures

Columns of both files: `index, filter_mean, filter_cov, smoother_mean,
smoother_cov, ACI_metric` (301 rows each). Compared at `ref$index`.

| line | export | fixture | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|---|
| 7 | `aci_filter` | `predprey_reference_predator_to_prey.csv` | `test-19:228-253` | t19 = 1e-6 | authors_source_comparison | filter mean and covariance, predator hidden | OK |
| 43 | `aci_predprey_model` | `predprey_reference_prey_to_predator.csv` | `test-19:228-253` | t19 = 1e-6 | authors_source_comparison | both directions, all five columns | OK |

## Rows resting on `mv_reference.csv`

Columns: `index, fm1, fm2, fc11, fc12, fc22, sm1, sm2, sc11, sc12, sc22,
ACI_metric` (201 rows).

| line | export | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|
| 8 | `aci_filter` | `test-19:256-292` | t19 = 1e-6 | independent_transcription | matrix filter mean (`fm1, fm2`) and covariance (`fc11, fc12, fc22`); `expect_gt(max(abs(ref$fc12)), 1e-3)` proves the cross term is exercised | OK |

## Rows resting on `enso_reference.csv`

Columns: `index, fm1, fm2, fm3, fc11, fc33, fc13, sm1, sm3, sc11, sc33,
ACI_metric` (201 rows).

| line | export | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|
| 44 | `aci_enso_model` | `test-19:295-338` | t19 = 1e-6 | source_derived_comparison | the eleven moment and metric columns of the joint three-hidden run; `matlab_simulator_parity == FALSE` asserted at `:334` | OK, indirect on the word "coefficients" (see note A) |

## Rows resting on the online fixtures

`cir_online_reference.csv` and `cir_cross_online_reference.csv`: columns
`j, n, online_mean, online_cov`, 75 rows each, lags 0 to 1950.
`mv_online_reference.csv`: columns `j, n, om1, om2, oc11, oc12, oc22`, 45 rows,
lags 0 to 3949.

| line | export | fixture | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|---|
| 17 | `aci_online` | `cir_online_reference.csv` | `test-27-online-lag.R:68-82` (also `test-19:341-357`) | 1e-6 literal at `test-27:79-80` | authors_source_comparison | fixed-lag mean and covariance at the 75 pinned `(j, n)` pairs | OK |
| 18 | `aci_online` | `cir_cross_online_reference.csv` | `test-19:360-379` | t19 = 1e-6 | independent_transcription | the same two columns under non-zero cross-noise | OK |
| 19 | `aci_online` | `mv_online_reference.csv` | `test-27-online-lag.R:201-222` (also `test-19:382-399`) | 1e-6 literal at `test-27:217-220` | independent_transcription | matrix fixed-lag mean (`om1, om2`) and covariance (`oc11, oc12, oc22`); `max(n - j) = 3949` matches the row's text | OK |

## Rows resting on the forward-CIR range fixtures

Columns of both files: `j, peak, objective, subj_0.1, subj_0.01, subj_0.001,
subj_0.0001`, 15 rows. Comparator helper `.compiled_oracle_cir()`,
`test-19:402-441`, which also re-runs under `convention = "lag_time"`.

| line | export | fixture | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|---|
| 23 | `aci_range` | `cir_range_reference.csv` | `test-19:443-452` | 1e-10 at `:448-449` | authors_source_comparison | `peak` and `objective` at the 15 pinned anchors, anchors asserted identical | OK |
| 24 | `aci_range` | `cir_range_reference.csv` | `test-19:443-452` | 1e-12 at `:450-451` | authors_source_comparison | the four `subj_*` columns under both read-out conventions | OK |
| 25 | `aci_range` | `cir_cross_range_reference.csv` | `test-19:455-465` | 1e-10 at `:460-461` | independent_transcription | `peak` and `objective` under `S_yo S_x = 0.54` | OK |

## Rows resting on the scalar ENSO partition fixtures

`enso6_partition_u_reference.csv`, `..._hW_...`, `..._tau_...`: 201 rows,
columns `index, t`, twenty-odd `coef_*` columns, then `ref_filter_mean,
ref_filter_cov, ref_smooth_mean, ref_smooth_cov, ref_aci, truth`.
`enso6_partition_tau_reduced3_reference.csv`: the same layout at k = 3.
Comparator `.partition_max_error()`, `test-28:131-137`, at `idx = ref$index`
(201 sampled indices), plus `.partition_summary_rows()` (`test-28:159-176`)
against `enso6_partition_fullpath_summary.csv`.

| line | export | fixture | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|---|
| 3 | `aci` | `enso6_partition_u_reference.csv` | `test-28:281-356` (`hidden = "u"`) | t28 = 1e-8 | source_derived_comparison | `ref_aci`, the total from `aci(..., decompose = TRUE)$aci` | OK |
| 9 | `aci_filter` | `enso6_partition_hW_reference.csv` | `test-28:281-356` (`hidden = "hW"`) | t28 = 1e-8 | source_derived_comparison | `ref_filter_mean`, `ref_filter_cov` | OK |
| 13 | `aci_smoother` | `enso6_partition_tau_reference.csv` | `test-28:281-356` (`hidden = "tau"`) | t28 = 1e-8 | source_derived_comparison | `ref_smooth_mean`, `ref_smooth_cov` | OK |
| 45 | `aci_enso_model` | `enso6_partition_tau_reference.csv` | `test-28:281-356` | t28 = 1e-8 | source_derived_comparison | the `coef_*` block, realised coefficients column by column | OK |
| 46 | `aci_enso_model` | `enso6_partition_tau_reduced3_reference.csv` | `test-28:359-395` | t28 = 1e-8 | source_derived_comparison | the k = 3 `coef_*` block and the five `ref_*` series | OK |
| 56 | `aci_conditional` | `enso6_partition_tau_reduced3_reference.csv` | `test-28:359-395` for the grade, `test-28:398-414` for the identity | t28 = 1e-8, and `identical()` for the specification | source_derived_comparison | the declared reduced run against the pinned bytes, with `aci_conditional(given = c("u","hW"), method = "reduce")` asserted bit-identical to the declaration | OK |

Note on the full-record summary: `enso6_partition_fullpath_summary.csv`
(columns `partition, series, n, min, max, mean, sum_abs`, 108 rows) is pinned by
`oracle-manifest-partitions.yml` and asserted at `test-28:342-346` and
`:383-387`. It is not itself a register row. Its four reductions are
order-insensitive; the wording in both the manifest and `test-28` was corrected
under ADV-E05 to say it bounds magnitude over the record and fixes position only
at the sampled indices.

## Rows resting on the T_C zeroth-order fixtures

`tc_outputs_*` files: columns `t, filter_mean, filter_cov, smoother_mean,
smoother_cov, aci_signal, aci_dispersion, aci_total`. Case A files carry 4001
rows; `..._caseB_window_head_...` carries 200.
`tc_coefficients_caseA.csv`: columns `t, L_y_state_time, L_y_matlab_phase,
f_x_TE_state_time, f_x_TE_matlab_phase, f_y_intended, f_y_literal, f_x_I,
sigma_I_matlab, sigma_I_floored`, 4001 rows.
`tc_d1_divergence.csv`: columns `case, quantity, n, identical, max_abs, rmsd`,
16 rows.
The graded quantity vector is `.tc_graded` (`test-29:122-123`), all seven
output columns at every step.

| line | export | fixture | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|---|
| 4 | `aci` | `tc_outputs_caseA_intended.csv` | `test-29:402-417` (`defect = FALSE` arm) | t29 = 1e-12 | independent_transcription | all seven `.tc_graded` columns at 4001 steps, `aci_signal` and `aci_dispersion` included | OK |
| 10 | `aci_filter` | `tc_outputs_caseB_window_head_intended.csv` | `test-29:419-447` (`part = "head"`, `tag = "intended"`) | t29 = 1e-12 | independent_transcription | `filter_mean`, `filter_cov` over the 200-row window head | OK |
| 14 | `aci_smoother` | `tc_outputs_caseA_literal.csv` | `test-29:402-417` (`defect = TRUE` arm) | t29 = 1e-12 | independent_transcription | `smoother_mean`, `smoother_cov` on the literal-defect arm | OK |
| 47 | `aci_enso_model` | `tc_coefficients_caseA.csv` | `test-29:352-400` | t29 = 1e-12 | independent_transcription | `L_y_state_time`, `f_x_TE_state_time`, `f_x_I`, `f_y_intended`, `f_y_literal` and `sigma_I_floored` compared at 1e-12; `sigma_I_matlab`, `L_y_matlab_phase` and `f_x_TE_matlab_phase` used as the other side of the recorded convention gaps (see note B) | OK |
| 48 | `aci_enso_model` | `tc_d1_divergence.csv` | `test-29:452-512` | `identical()` where `identical == TRUE`, else 1e-6 relative on `max_abs` and `rmsd` | independent_transcription | the recorded divergence between the two D1 arms, per quantity and per case | OK |

## Row resting on `dyad_signal_x.csv`

| line | export | fixture | test block | tolerance | source class | graded quantity | verdict |
|---|---|---|---|---|---|---|---|
| 61 | `observed_trajectory` | `dyad_signal_x.csv` (headerless, two columns `t` and `x`, 30001 rows) | `test-19:188-205` (read at `:189-191`, container built at `:192`) | t19 = 1e-6 on the run it drives | authors_source_comparison | MISMATCH of kind: the fixture is the pinned INPUT signal, not an output of `observed_trajectory()`. Nothing compares a column of this file against a value the export computes; what the row grades is that the container carries `t`, `x` and `dt` onto the graded paths. | MISMATCH, disclosed; no change proposed |

The row's against text already reads "the pinned MATLAB signal every
authors-source grade is driven by", which states the situation accurately, so
the plan proposes no edit. Recorded here so the distinction is inspectable.

---

## Notes

**A. Constructor rows graded through the run, not through a coefficient
column.** Rows 41 (`aci_dyad_model`) and 44 (`aci_enso_model`, joint case) name
"coefficients and Grammians" and "coefficients and vector recursions" as the
feature, but `dyad_reference.csv` and `enso_reference.csv` carry only moment and
metric columns. The grade is indirect: the constructed model is the one the
pinned run is computed from, so a coefficient regression moves the graded
moments. Both rows' against text says so ("the constructor realises the model
the pinned grade runs on"; "a source-derived harness on an Euler-Maruyama
realisation"). `test-19:335-336` does assert one coefficient property directly
(the `gxx[3,3,]` seasonal ratio exceeds 2), but against a threshold, not against
a fixture column. This is weaker than the partition rows 45 and 46, whose
fixtures carry the `coef_*` columns and grade them one by one. No change is
proposed here; the observation is recorded because the feature wording is
stronger than the column list on its own would support.

**B. Convention columns in `tc_coefficients_caseA.csv`.** All ten columns are
used, but not all as direct comparators. `sigma_I_floored` is graded at 1e-12
(`test-29:380`). The alternative-convention columns are used as the other side
of a stated gap rather than as expectations: `sigma_I_matlab` appears in an
exact identity (`sigma_I^2 = sigma_I_matlab^2 + 1e-3 * lambda`,
`test-29:386-387`) and in a two-sided bound (`:389-390`), and
`L_y_matlab_phase` and `f_x_TE_matlab_phase` are asserted to differ from the
shipped state-time columns by more than 1e-6 (`:398-399`), which is what makes
the constructor's phase choice a documented decision rather than an accident.
Nothing in the fixture is dead.

**C. Class distribution over the 28 checked rows.** authors_source_comparison
10 (lines 2, 6, 7, 11, 17, 23, 24, 41, 43, 61); source_derived_comparison 7
(3, 9, 13, 44, 45, 46, 56); independent_transcription 11 (4, 8, 10, 12, 14, 18,
19, 25, 47, 48, 50). The other register rows are `behavioural_only`, with no
fixture behind them: 17 `behavioural`, 15 `exact_relation` and 8
`independent_transcription`, 37 of them before this round's three new rows and
40 after. The 28 checked rows, the 10 authors-source rows and the nine of
those that compare a computed output against the authors' fixture (row 61
being the tenth, and an ingestion check) are the figures the corrected
`README.md` and `NEWS.md` wording states.

**D. Fixtures cited by more than one row.** `dyad_reference.csv` (rows 2, 6,
11, 41), `cross_reference.csv` (12, 50), `cir_range_reference.csv` (23, 24),
`enso6_partition_tau_reference.csv` (13, 45),
`enso6_partition_tau_reduced3_reference.csv` (46, 56). Each pair grades a
different column group of the same file, so the sharing is not double-counting.

## Register edits made in this round

- Line 2: feature and against narrowed to the total ACI metric, naming the
  fixture's lack of a parts column and the comparator line. `check_method`,
  `tolerance_class`, `state`, `fixture_path` and `sha256` are byte-identical to
  before.
- Line 36 (`lt_tail_bound`, behavioural_only): "truncation tail bound" became
  "heuristic truncation tail estimate", and the against text now says the
  1.5e-4 figure is a recorded value checked behaviourally, not a certificate.
  This is MATH-03's register clause. Machine-read columns unchanged.
- Line 39 (`aci_metric_pair`) unchanged, as decided.
