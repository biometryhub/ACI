# aci: Assimilative causal inference in R
[![R-CMD-check](https://github.com/aidanmoller/aci/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aidanmoller/aci/actions/workflows/R-CMD-check.yaml)

Dynamical causal inference for partially observed stochastic systems. This package
implements assimilative causal inference (ACI) (Andreou et al., 2026; Andreou & Chen, 2026;
Jiang et al., 2026; Moser et al., 2026) and companion papers,
*does* a latent process drive the observed ones, *how strongly over time*
(the ACI metric: relative entropy of smoother vs filter posteriors), and
*for how long* (forward/backward causal influence ranges).  with a formula
front-end, grid-based CIR estimators, conditional ACI, ensemble engines, model
discovery, and a documented partial implementation of the Moser et al. (2026)
extreme-event workflow.

This is a working R prototype adapted from three MATLAB codebases, and four papers
(see `?aci_references`). The package is currently in development, and may be 
different the next time you try to use it. Future updates are expected to refine 
rather than extend current features. If you would like to contribute please feel
free!

To get started, follow the steps below or see `vignette("vignette-1-intro")` for
an overview.

## Installation

```r
# install.packages("remotes")
remotes::install_github("aidanmoller/aci")
```

## Condensed example

```r
library(aci)
sim <- simulate(model_dyad(), seed = 8, T = 5, dt = 5e-3, burn_in = 1)
dat <- data.frame(t = sim$obs$t, x1 = sim$obs$x[, 1], y = sim$hidden[, 1])

fit <- aci_fit(y ~ x1, data = dat)   # learn the model, assimilate, score
summary(fit, verdict = TRUE)          # equations + ACI + significance
plot(fit, "latent")                   # reconstruction vs time
cir(fit)                              # how long the influence persists
```

Drop the `y` column and `aci_fit()` reconstructs the latent cause from the
observations alone (experimental). For mechanistic control, build models
directly (`model_dyad()`, `model_enso6()`, `cgns_model()`, ...) and use
`aci()`, `lag_table()`, `forward_cir()`, `backward_cir()`.


## Benchmark models

`model_dyad()`, `model_predator_prey()` (both causal directions),
`model_tipping_triad()`, `model_pathways()`, `model_l84()`, `model_l96()`,
`model_enso6()` (Chen, Fang and Yu (2022) parameters or the reference-code variant), 
`model_topographic()`. Topographic systems and four-state multiscale models based 
on Andreou and Chen (2026) / FBCIR multiscale and the MATLAB Codebase(s) have separate 
constructors (`model_multiscale_fbcir()` & `model_topographic_layered_fbcir()`), 
and do not all have the same validation status. 


## Numerics cheat sheet

* **Time step**: keep $\Delta t\,\lVert L_y\rVert \lesssim 0.1$; the filter
  warns otherwise. Sub-step with `nsub` rather than thinning data.
* **Stepper**: `"explicit"` is the reference scheme and the lag table
  contract; `"implicit"` is unconditionally SPD for stiff/learned models.
* **Priors**: pass an explicit `init` (the reference uses mean = truth,
  cov = 0.1 in synthetic studies). The default diffuse prior warns once:
  discard a burn-in window when you use it.
* **Burn-in**: both the metric and the CIRs inherit transient prior
  influence; `nil_surrogate_test()` drops 10% by default.
* **Suppressing known warnings**: they are classed; e.g.
  `withCallingHandlers(..., aci_warn_diffuse_init = function(w)
  invokeRestart("muffleWarning"))`.

## Limitations and roadmap

* O1: Partial-observation learning remains experimental.
* O2: Some helpers are explicit package extensions rather than results specified in 
the four papers.
* O3: Simulator/pathwise parity of the golden ENSO drift block is still open (R 
uses Euler-Maruyama, the scripts mix Euler and Milstein).
* O4: MATLAB and R random-number streams are not interchangeable (numerical 
comparisons must use identical saved trajectories or increments).
* O4: the Moser et al. (2026) event layer remains incomplete: representative/MAP event paths, 
weighted uncertainty quantities in Eqs. 3.21--3.29.
* O5: what's implemented of the event workflow (Moser et al., 2026) are checked 
against their equations and R regression tests, but do not produce MATLAB-parity 
results.

## References

Source material and literature for implementations are cited throughout the help pages 
(see `?aci_references` for a full list).

Pre-processing and calibration scripts use published reference implementations, 
where available (`screen_outliers()`, `as_uniform_trajectory()` and `calibrate_cgns()` 
are documented as package extensions, not translations).
