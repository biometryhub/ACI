# Intent: correct the plot, kernel and table-reuse faults

Author: the `acir` authors. Status: draft 2026-11-20.

## Problem

Since PR #30 the plot methods for `aci_result`, `da_path_gaussian` and
`cir_result` build their call from `utils::modifyList()` on `list(...)` and
`do.call()`. Unnamed arguments are dropped and an empty argument fails. Every
argument is evaluated before the plot starts, so `panel.first` and
`panel.last` fail on a device without a plot and draw on the previous figure
otherwise. The data are copied into the call, which slows plots of long
series. A range with one anchor stops on values such as `lwd = 0`.

PR #30's scalar KL kernel no longer recycles arguments of different lengths,
stops with an untyped error when both variances are infinite, returns `Inf`
for some signals just below the overflow threshold, and does extra work when
nearly every signal overflows.

Two faults predate PR #30. `aci(table = )` refuses a lag table that
`lag_table()` built for a model that declares its own estimand, and plotting
a range whose values are all `NA` fails with a base graphics error.

## Proposed outcome

Plot calls that worked before PR #30 draw the same plot again, and named
labels still replace the defaults. The kernel recycles its arguments, returns
`NaN` for two infinite variances and decides overflow exactly at the rounding
midpoint. `aci(table = )` accepts the table, and the all-`NA` range stops with
a classed condition that asks for `ylim`.

## Affected code and users

`acir/R/aci-methods.R`, `acir/R/aci-kernels-scalar.R` and
`acir/R/aci-assimilation.R`, with their help and tests. Users who pass extra
arguments to `plot()`, and users who reuse lag tables for a model that
declares its own estimand.

## Constraints

No change to any finite result on ordinary inputs, to the exported interface,
or to any fixture or tolerance.

## Open questions

None.
