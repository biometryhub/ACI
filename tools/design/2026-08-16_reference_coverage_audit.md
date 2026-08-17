# Reference coverage audit — what aciR does not implement, and why

2026-08-16. Package at 0.2.1. Reference clone `matlab_reference/`
(`github.com/marandmath/ACI_code`, MIT, Copyright (c) 2025 Marios Andreou).

## Why this exists

The question "what features of the MATLAB package are not implemented" had no
single answer in the tree. `PARITY_REPORT.md` §6 answers a *different* question
-- which quantities have been **graded** against the authors' code -- and the
two axes had been conflated. This audit differences the reference's capability
surface against the package's, and records the reason for each gap or the fact
that no reason was recorded.

## Method, and how to re-run it

Section banners were extracted from the six model scripts, not recalled:

```sh
cd matlab_reference
grep -oE '%{4,} +[A-Z][A-Z0-9 ,\-()/]+ +%{4,}' dyad_interaction_model.m |
  sed 's/%//g'
```

Executable content of a section was counted by stripping comment and blank
lines between its banner and the next, e.g. for the data-loading section:

```sh
for f in ENSO_model_cond_ACI_*.m; do
  s=$(grep -nE '%{4,} +LOADING THE OBSERVATIONAL DATA +%{4,}' $f | cut -d: -f1)
  e=$(grep -nE '%{4,} +MODEL SETUP +%{4,}' $f | cut -d: -f1)
  awk -v a=$s -v b=$e 'NR>a && NR<b' $f | grep -vcE '^\s*%|^\s*$'
done
```

All five ENSO scripts were checked individually rather than one being taken as
representative, on the F7 precedent that these five differ in exactly this way.

## Section-by-section difference

| Reference section | aciR | Status |
|---|---|---|
| MODEL SETUP | `aci_cgns_model`, `aci_dyad_model`, `aci_predprey_model`, `aci_enso_model` | covered, wider surface |
| GENERATING THE TRUE SIGNALS | `aci_simulate` (Euler-Maruyama, Milstein) | covered |
| FILTERING | `aci_filter` | covered |
| SMOOTHING | `aci_smoother`, `aci_online_smoother` | covered |
| ACI ANALYSIS | `aci_metric`, `aci`, `aci_cir`, `aci_conditional` | covered |
| LOADING THE OBSERVATIONAL DATA | -- | **not implemented** |
| OBTAINING MONTHLY CLIMATOLOGY | -- | **not implemented** |
| PLOTTING ACI ANALYSIS RESULTS | `plot.aci`, `plot.aci_cir` | **not implemented** |

Every computational section is covered. What is absent is the case-study
apparatus around the method.

## The gaps

### 1. Observational-data ingestion — not implemented, no reason recorded

The reference ships `ENSO_DATA/` with eight files (ERSST v5 Niño3, Niño3.4 and
Niño4 indices; GODAS `u`/`h_W`/`T_C`/`T_E`; GODAS `I`; NCEP-NCAR wind stress;
the Niño3.4 SSTA regression coefficients; a longitude grid) and a section of
roughly 700 lines describing how to load and align them.

**Measured: that section contains zero executable lines in all five ENSO
scripts.** It is commented guidance. The scripts run on simulated data.

So there is no executable upstream behaviour to reimplement, which is the same
situation as conditional ACI (F7) -- a documented manual procedure rather than
code. aciR records no statement of this either way.

### 2. Monthly climatology — designed out, reason in a code comment only

The reference post-processes simulated series with `movmean` at the monthly
step, subtracts the mean, and rescales each variable by a fixed factor
(`1.5`, `150`, `7.5`, `7.5`, `5`).

aciR does not have the step because it does not need it: `R/enso-model.R:77`
enforces zero climatology in the anomaly model, so the package works in anomaly
space throughout. The reason is a comment in the source and appears in no
user-facing document.

### 3. Figures — partly reasoned, partly unrecorded

The reference produces whisker CIR overlays, Hovmöller diagrams, and a spatial
SST reconstruction obtained by combining the loaded regression coefficients
with the longitude grid, together with a 12-month-window ENSO event
classification. aciR provides `plot()` methods for its two result classes and
nothing else.

`README.md` §"Reproducing the published causal-influence-range panels" records
that the defaults "are chosen for a reusable package rather than for
transcription" and gives the three arguments (`horizon`, `margin`, `epsilon`)
that recover the paper's CIR panels. That covers the CIR figures. Nothing
records why the spatial reconstruction and the event classification are absent.

### 4. Singular observation-noise covariances — deliberate, well recorded

The only place aciR declines work the reference will do. The reference forms
`pinv(S_x * S_x')` "for stability concerns"; aciR uses a Cholesky inverse and
requires positive-definiteness at every slice, so it rejects a path whose
Walker circulation reaches either end of `[0, 4]`.

Recorded in `PARITY_REPORT.md` F8 and, user-facing, in the
`aci_enso_components()` documentation: quietly nudging a singular covariance
"would be inventing an observation the data does not support".

## Where aciR is the wider surface

Recorded so the difference is not read as one-directional.

- `R0` is an argument; the reference hardcodes the initial covariance at `0.1`
  inside the filter range.
- Noise Grammians are accepted directly; the reference only forms them from
  feedback columns, so not every admissible Grammian pair is reachable there.
- `aci_conditional()` generalises what the reference does by hand-editing one
  `S_xoS_x_inv` line.
- Vector states run through the online smoother and the causal influence range.
- The latent self-drift may vary in time.
- `objective_exact` and `horizon` were added at 0.2.0 and close the last two
  capability gaps that did exist (F1, F2).

## Correction produced by this audit

The claim that the reference has "only `simps.m` and `legendUnq.m` callable,
both third-party" is **false**. There are three callable files, and
`progress_bar.m` is the authors' own; all six model scripts call it. Corrected
2026-08-16 in `PARITY_REPORT.md` §1, in the hoist cairn (as a dated amendment),
and in §02 of the public Development Ledger. The claim had not reached any
shipped package surface, so no version bump was required, and the conclusion it
supported is unaffected.

Worth noting for the future: the package's `test-retired-claims.R` guard covers
`R/`, `man/`, the vignettes, `DESCRIPTION`, `README.Rmd`, `NEWS.md` and the
oracle manifest. It does **not** cover `tools/oracle/parity/reports/`, `tools/ledger/` or
`public/`, which is where this claim lived for two days.

## Open actions

1. Add a scope paragraph to the assumptions article stating that the
   case-study apparatus -- data ingestion, climatology preparation, event
   classification, figure production -- is out of scope, and why. Gaps 1 and 3
   currently have no stated reason anywhere a user would look.
2. Consider extending the retired-claims guard beyond the package directory to
   the project-level documents that make claims about the reference.
