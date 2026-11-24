# Specification: plot arguments, scalar KL kernel and table reuse

From `intent.md` (2026-11-20). Status: draft 2026-11-20.

## Requirements

1. Each plot method passes its dots to `plot()` unevaluated and in order, as
   a direct call would. A caller's `type`, `xlab`, `ylab` or `main`, and
   `ylim` for a Gaussian path, replaces the method's default.
2. A range with one anchor passes only `main`, `sub`, `xlab` and `ylab` to
   `title()` and evaluates no other argument.
3. A range with more than one anchor whose values are all `NA` stops with
   `aci_error_dims` unless the caller gives a `ylim` that is not `NULL`, by
   name, partial name or position.
4. The scalar KL kernel recycles its arguments, returns a `NaN` dispersion
   for two infinite variances, and returns an infinite signal exactly when the
   exact signal reaches the overflow midpoint `2^1024 - 2^970`. Its results on
   ordinary inputs are unchanged bit for bit.
5. When the caller of `aci(table = )` gives no `conditional`, the table must
   carry the specification the model declares, if it declares one.

## Design

Each plot method draws through a local function whose defaults follow the
dots. R's argument matching then lets a caller's argument of the same name
replace a default and passes every other argument on as a promise. A range
with one anchor uses a second local function that takes only the title text.
The all-`NA` check matches the call against `plot.default()` with
`match.call()`, which finds `ylim` however it was given, and evaluates only
that argument, to see whether it is `NULL`.

The kernel recycles its arguments before the branches that index single
elements. It decides overflow by the exact sign of a sum of error-free terms
(TwoSum, TwoProduct and Shewchuk's Grow-Expansion), and an argument half that
underflows to zero breaks an exact tie by its sign. Elements whose
standardised difference exceeds `2^514` skip that work, because their signal
cannot round back into the double range.

`.validate_lag_table_source()` uses `.model_estimand_spec()` to supply the
model's declared specification when the caller gives none.

## Acceptance

The four gates of `dev/PROCESS.md`. The new tests in test-24, test-42,
test-43 and test-45 fail on `17487dc` and pass after the change.

## Flagged concerns

None.
