# Production compiled-CGNS benchmark

`scalar-dyad.R` benchmarks the current production architecture; the historical
name is retained so old run commands remain recognizable. It is development
infrastructure, not a timing unit test.

Run a full record from the package root:

```sh
Rscript --vanilla benchmarks/scalar-dyad.R \
  --output=benchmarks/results/production-compiled-cgns
```

For a syntax and execution smoke run, use:

```sh
Rscript --vanilla benchmarks/scalar-dyad.R \
  --output=/tmp/aci-production-smoke --quick --no-alloc --no-profile
```

The script needs the development-only packages `pkgload`, `microbenchmark`, and
`digest`. They are not installed-package dependencies. `--no-alloc` and
`--no-profile` are intended only for smoke/debug runs; omit them for evidence
used in a performance report.

The scenarios are fixed and source-hashed:

- `dyad_3001`: the canonical seed-1, 3,001-step scalar dyad, whose observation
  SHA-256 is checked before any measurements;
- `affine_matrix`: a deterministic 401-step `k=2, l=2` affine model exercising
  authenticated batch realization and the matrix kernels;
- `conditioned_generic`: a deterministic 401-step `k=2, l=1` closure model with
  masked-innovation (`inflate`) conditioning; and
- `bounded_range`: the first 201 dyad points, a forward table capped at 25 lags,
  and streaming forward CIR without a retained triangular table.

Each run writes raw timings, timing summaries, stage contracts, fresh-process
allocation logs and summaries, selected raw profiles and profile summaries,
parity results, input hashes, complete R-source hashes, Git state when available,
and R/system/library details. `parity.csv` must pass before timings are accepted.

Interpret rows by `layer`; they are not interchangeable:

- `compilation` realizes and validates coefficients;
- `warm_kernel` consumes an existing coefficient bundle;
- `bundle_complete` performs filter, smoother, metric, likelihood, and result
  construction over an existing bundle;
- `public_workflow` includes the documented API boundary and fresh compilation;
  and
- `warm_lag_cir` consumes a precompiled bundle for bounded lag/CIR algorithms.

`stages.csv` records whether each expression includes compilation, predictive
likelihood, validation, result construction, or retained lag rows. Reports
should quote the raw median and IQR (plus allocation evidence) for the relevant
contract. Do not collapse these rows into one headline speed ratio.

The current benchmark never loads aci 0.0.21 into the same R process: both source
trees define package `aci`, so doing so produces an invalid namespace comparison.
Use the preserved raw data under `benchmarks/results/scalar-dyad-20260825-final/`
for the historical 0.0.21 machine run. If a fresh 0.0.21 measurement is needed,
run `aci-prev` from a separate `Rscript --vanilla` process and retain its own
source hashes and environment record. Never reintroduce retired private kernels
into the current package as a comparator.

The historical result directories are immutable evidence from the earlier
shootout. Their stage labels describe the code that produced them and are not a
description of current routing.
