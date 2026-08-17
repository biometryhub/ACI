---
id: 2026-08-14-minimal-substitution-capture-gate
schema_version: 1.4
date: 2026-08-14
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: infrastructure
project: aciR
status: accepted
title: "Reference workspaces captured under a machine-verified minimal substitution, not an unmodified run"
tags: [oracle_scope, parity_harness, independent_oracle, provenance]
triggers:
  - adr_class_commitment
  - methodology_change_downstream
reversal_cost: low
decision_pressure: null
review_due: null
review_trigger: "if a machine with enough memory and patience makes a full unmodified run practical, the substitution should be dropped rather than defended"
supersedes: []
superseded_by: null
related:
  - 2026-08-14-extraction-by-byte-verified-hoist
  - 2026-08-13-grading-register-as-build-gate
---

## Status

Decided and implemented 2026-08-14 as part of the MATLAB/aciR parity harness
(`tools/oracle/parity/`). The gold side of gate G1 is a workspace captured from the
reference script after applying the problem-size substitutions declared in
`manifest/knobs.dcf`, and `tools/reduce.R` asserts the substitution is minimal
before the script is run.

The original plan said "run the original script unmodified". That is not
available, and the arithmetic is not close. At its published settings
`dyad_interaction_model.m` allocates `RE_metric` as 20601 x 20601 (3.40 GB) and
three online-smoother cell triangles at N = 30001 (about 10.8 GB), then scans
513 epsilon values against 20601 reporting times over rows up to 20601 long --
of order 1e11 element operations. On the 32 GB workstation this project runs
on, that is days, not a gate.

Three knobs are substituted for the dyad -- `N`, `time_start_plot`,
`time_end_plot` -- and the reduction checks its own work: each knob must match
exactly one bare literal assignment, the set of lines that differ must equal
the set of lines declared, every changed line must be a bare
`<name> = <literal>;` on BOTH sides, and the file length must be unchanged so
that the line ranges the extraction manifest depends on still address the same
code. Line terminators are preserved byte for byte, so `diff` against the
reference shows three lines and nothing else.

Two profiles are declared rather than one, at different sizes, and
deliberately on opposite branches of the `last_idx` conditional: `small`
(T = 3.0) takes the lookahead branch the published settings take, `tiny`
(T = 1.2) takes the truncating branch.

## Alternatives considered

**1. Run the script unmodified and wait.** The honest ideal, and the one the
plan promised.

**2. Capture the cheap prefix with a debugger breakpoint** (`dbstop` before the
quadratic sections), leaving the script bytes untouched.

**3. Reduce the size but do not verify minimality** -- edit the knobs by hand,
state that only sizes changed.

**4. Abandon workspace capture; grade the extracted functions against the
existing hand-transcribed harnesses instead.**

## Rationale for rejection

**(1)** Rejected on arithmetic, not on impatience: roughly 14 GB of peak
allocation and of order 1e11 element operations in interpreted MATLAB. A gate
nobody can run is not a gate, and one run per week would not have caught the
degenerate profile described below.

**(2)** Attractive and genuinely non-invasive, and worth revisiting. Rejected
for now because the debugger's behaviour under `-batch` is not something this
project has established, and a gate whose mechanism is itself unverified buys
nothing over one whose modification is visible in a three-line diff. It also
captures only a prefix: the CIR sections would still be ungraded.

**(3)** Rejected because it is exactly the class of failure this package
exists to avoid. "Only sizes changed" is a claim about a file, and a claim
about a file that no program checks will eventually be false. The minimality
assertions cost about forty lines.

**(4)** Rejected because it defeats the purpose. The transcribed harnesses
share an author with the R implementation; grading extracted reference code
against them would re-introduce the shared step the whole harness exists to
remove.

## Forward cost

Each new reference script needs a knobs record with a justification and a
`Requires` list, and the profiles must be re-checked whenever a script is
re-cloned upstream. The reduction is re-run from the manifest, so the cost is
per-script, not per-run.

The load-bearing follow-on cost is the `Requires` field, and it earned its
place on the first run. The initial `tiny` profile used a reporting window of
0.2 while the reference's `lookahead_tolerance` is 0.6; the reference indexes
its exact objective CIR as `subjective_CIR(:, 1:end-lookahead_tolerance/dt)`,
which MATLAB evaluates to an EMPTY range rather than raising, so
`defn_objective_CIR` came back 1x0. Any pairing graded against it would have
compared nothing to nothing and passed. Every profile now declares conditions
that are evaluated against the captured workspace before anything is graded,
and a degenerate profile is an error rather than a quiet success.

## References

- `tools/oracle/parity/manifest/knobs.dcf` -- the declared substitutions and their
  justifications
- `tools/oracle/parity/tools/reduce.R` -- the minimality assertions
- `tools/oracle/parity/matlab/check_profile.m` -- the degeneracy gate
- Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal inference.
  *Nature Communications*, 17, 1854. \doi{10.1038/s41467-026-68568-0}
