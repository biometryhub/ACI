# API stability

This document states what a user of aciR may rely on, and how changes to it
will be handled. It applies from version 0.1.0.

## Lifecycle stages

aciR uses three stages, following the conventions of the
[lifecycle](https://lifecycle.r-lib.org/articles/stages.html) package.

* **Stable** — the signature and semantics will not change in a breaking way
  without a deprecation cycle: at least one minor release in which the old form
  still works and warns, and a `NEWS.md` entry with before/after code.
* **Maturing** — the design is settled and in use, but the signature may still
  change in a minor release. Any change is announced in `NEWS.md` with
  before/after code; a change that would silently alter results rather than
  raise an error gets a deprecation cycle regardless of stage.
* **Experimental** — the design is still being worked out and may change
  without a deprecation cycle.

## Stage per export

| Export | Stage | Notes |
|---|---|---|
| `aci()` | Maturing | The canonical entry point. `time` is the newest argument and the most likely to gain siblings. |
| `aci_cgns_model()` | Maturing | The general constructor. Will gain arguments as the core widens (a time-varying `L_y` is roadmapped); existing arguments are settled. |
| `aci_dyad_model()` | Maturing | Parameterisation follows the reference implementation and is not expected to change. |
| `aci_simulate()` | Maturing | See the reproducibility contract below. |
| `aci_filter()` | Maturing | Expert surface; see the components contract below. |
| `aci_smoother()` | Maturing | Expert surface. |
| `aci_metric()` | Maturing | Expert surface. |
| `aci_online_smoother()` | Experimental | New in the development version. The `lag` and `tol` arguments are settled in meaning, but the returned object may gain fields as the causal influence range and the adaptive-lag variant develop. |
| `aci_cir()` | Experimental | New in the development version. The returned object's shape, and in particular how an unresolved time is reported, may change: the saturation margin is this package's own device and has no counterpart in the reference implementation. |
| `aci_dyad_components()` | Maturing | The worked example of the components schema. |
| `aci_components` | Maturing | The components schema itself; see below. |
| `print.aci_model()`, `print.aci()` | Maturing | Printed output is for humans and may be reformatted in any release. Do not parse it; use `summary()` or `as.data.frame()`. |
| `summary.aci()` | Maturing | The returned object's fields are part of the contract; the printed form is not. |
| `as.data.frame.aci()` | Maturing | Column names and order are part of the contract. New columns may be appended. |
| `plot.aci()` | Experimental | The panel layout and styling may change. |

Nothing is Stable at 0.1.0. This is a research preview: the package has not yet
had enough external use for its design to have been tested by anyone but its
author, and declaring stability before that would be a claim about evidence
that does not exist. The intention is to promote the model-layer exports to
Stable at 1.0.0 once the roadmap items below have either landed or been
abandoned, since both would otherwise force a signature change.

## The components contract

`aci_filter()`, `aci_smoother()` and `aci_metric()` are a supported extension
surface, not implementation details: they exist so that a conditional Gaussian
system for which aciR supplies no constructor can still be filtered, smoothed
and scored. The schema they consume is documented at `?aci_components` and is
versioned with the package.

Within a minor release series:

* no entry will be removed or renamed;
* no entry's meaning, units or orientation will change;
* new *optional* entries may be added, and a components list built for an
  earlier version will keep working.

The one foreseeable change is the addition of a per-step `L_y` when the core
gains a time-varying self-drift. That will be additive: a scalar `L_y` will
continue to mean a self-drift constant in time.

## The reproducibility contract

`aci_simulate(seed = )` is contained. The generator state is saved before the
draw and restored when the function exits, so a seeded call does not consume
the caller's random stream. An unseeded call consumes the global stream in the
ordinary way.

The mapping from a seed to a path is **not** stable across versions. It changed
in 0.1.0, and it will change again if the integration scheme or the order of
draws changes. Reproducing a specific path across versions requires pinning the
package version; this is why the package's own validation runs on committed
fixtures rather than on a simulated path.

## The validation contract

The oracle tests are a release gate, not a convenience. They run from the
installed package fixtures, never skip, and must pass before any release. A
fixture is never regenerated in place: the refresh policy is recorded in
`inst/extdata/oracle-manifest.yml`, and a fixture whose hash changes without a
deliberate maintainer step is a defect rather than an update.

The manifest also records, per fixture, what it grades and what it does not. An
oracle grades the paths its parameters exercise and no others, and that scope is
part of what the package claims.

## Errors as interface

Error messages are part of the public interface and are locked by snapshot
tests. They may be reworded for clarity in any release; the condition that
triggers them will not change without a `NEWS.md` entry.

## Roadmap items that will force change

* A core admitting a **time-varying self-drift** `L_y`, needed for the noisy
  predator-prey model. Additive to the components schema.
* The **causal-influence-range** (`aci_cir()`). A new export; no change to
  existing ones.

Both wait on their own independent-oracle fixtures.
