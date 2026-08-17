---
id: 2026-08-14-section-runner-for-multi-direction-scripts
schema_version: 1.4
date: 2026-08-14
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: infrastructure
project: aciR
status: accepted
title: "Reference scripts that collide on variable names are composed section by section, not run whole"
tags: [parity_harness, independent_oracle, provenance, near_miss]
triggers:
  - adr_class_commitment
  - near_miss_caught
  - methodology_change_downstream
reversal_cost: low
decision_pressure: null
review_due: null
review_trigger: "if an upstream release separates the two directions into distinct variable names, the composition can collapse back to a whole-script run"
supersedes: []
superseded_by: null
related:
  - 2026-08-14-minimal-substitution-capture-gate
  - 2026-08-14-extraction-by-byte-verified-hoist
---

## Status

Decided and implemented 2026-08-14, on discovering that
`noisy_predator_prey_model.m` **cannot be run top to bottom without silently
producing wrong results**.

The script carries fourteen blocks headed "RUN THIS CODE SECTION TO STUDY THE
CAUSAL RELATIONSHIP", seven for each causal direction, and every block assigns
the same names. `filter_mean` is written for direction one at line 225 and
overwritten by direction two at line 276 before anything consumes it, so the
direction-one smoother at line 321 reads direction TWO's filter. The
direction-one ACI metric at line 601 is then built on that. The split reaches
into the model setup: `f_x` is a scalar in direction one (line 97) and an array
in direction two (line 126), same name.

None of this raises. A whole-script run completes and produces plausible
numbers for a direction-one quantity that was computed from direction two's
posterior.

`tools/reduce.R` therefore composes a **runner** from declared line ranges of
the reference, in file order -- the same lines a user executes following the
script's own instructions. Nothing is retyped or reordered; only selected. Each
profile declares a discriminating condition checked after capture:
`numel(f_x) == 1` holds only in direction one, `numel(f_y) == 1` only in
direction two, so a runner that ever picked up the wrong blocks fails before
anything is graded against it.

Composed runners are 454 and 453 lines from a 1,332-line script. Gate G1 then
passed **30 of 30 outputs at difference 0** across the two directions.

## Alternatives considered

**1. Run the script whole, as the dyad is run**, and take whatever the
workspace holds.

**2. Capture twice with a `clear` between directions**, or reorder execution by
editing the file.

**3. Drive the two directions by MATLAB code sections** (`%%` cells) via the
editor's section-execution API.

**4. Grade the predator-prey extracts against the existing transcribed
harness** rather than a captured workspace.

## Rationale for rejection

**1** is what the dyad procedure does and it would have produced confidently
wrong gold values for direction one, with nothing to signal it. This is the
near miss the cairn records: the dyad's structure made a whole-script capture
look like the general method, and it is not.

**2** requires editing the reference, which destroys the property the whole
harness rests on -- that the executed lines are byte-identical to the
published ones.

**3** was attractive because MATLAB genuinely has this notion. Rejected because
this script's banner rulers all begin with `%%`, so MATLAB sees section breaks
at nearly every banner line and the cell structure does not correspond to the
seven blocks the prose describes. Selecting by declared line range says what is
meant.

**4** re-introduces the shared-author step: the transcription and the R
implementation come from one reading. It is the thing being tested, not the
instrument to test it with.

## Forward cost

Every further script must be inspected for name collisions before a profile is
declared, rather than assumed to be linear. The five ENSO scripts are next and
have not been checked; they are ~2,000 lines each and the same pattern would be
easy to miss. The discriminating `Requires` condition is the cheap general
guard: name one variable whose shape or value differs between the intended
selection and any wrong one.

The payoff was immediate. With the two directions verified apart, both
committed predator-prey fixtures were re-graded against the authors' code at
the published N = 12000: `predprey_reference_predator_to_prey.csv` agrees
**exactly**, and `predprey_reference_prey_to_predator.csv` to 1.5e-13, which is
reassociation round-off over 12,000 steps. Had the transcription conflated the
directions -- the hazard this script's structure creates -- the disagreement
would have been of order one, not of order 1e-13. The fixtures keep the
directions correctly apart, and that is now measured rather than assumed.

## References

- `tools/oracle/parity/manifest/knobs.dcf` -- `predprey_dir1`, `predprey_dir2`
- `tools/oracle/parity/tools/reduce.R` -- section composition and its assertions
- `tools/oracle/parity/reports/regrade_predprey_dir1.csv`, `..._dir2.csv`
- `matlab_reference/noisy_predator_prey_model.m` lines 90, 119, 225, 276, 321
