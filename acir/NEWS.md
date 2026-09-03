# acir 0.1.0.9000

* The vignettes show the package's refusals as the classed condition and its
  message, printed as ordinary output, instead of a halted chunk. A genuine
  error in those chunks now fails the vignette build.
* Five condition messages no longer cite internal specification sections or
  development notes: the singular observation-noise Gram, the non-uniform
  time grid, the noise-free observation contract, and the inadmissible
  `method = "reduce"` reduction.

# acir 0.1.0

The parity release: every quantity the reference MATLAB implementation
computes is reproduced to the tolerance recorded in
`inst/evidence/register.csv`, and every stage of the specification's
performance table is inside its budget.

## Breaking changes

* **The simulation horizon is `t_end`, not `T`.** `simulate()` and
  `aci_simulate()` take `t_end`; the name `T` shadowed the alias of `TRUE`
  in R and is retired. A call that still passes `T` works, with a warning of class
  `aci_warning_deprecated`, until acir 0.2.0; passing both is an error.


## Bug fixes

* **The built package carries its vignette index.** The build-ignore rule
  `^build($|/)`, inherited from aci, removed the `build/` directory that
  `R CMD build` creates, so an installed package listed no vignettes under
  `vignette()` or `browseVignettes()` and the CRAN incoming check reported a
  missing vignette index. The rule is gone.

## Performance

* **A model supplied as closures is realised once per record, not once per
  verb.** The generic route evaluated every coefficient closure at every grid
  point on each call, so a filter, a smoother, an online smoother and an
  influence range on one model and one record realised the same arrays four
  times. The last four realisations are now kept, keyed on the model object
  and the grid, and a stored one is reused after its coefficients are
  re-evaluated at three grid points and found unchanged, which catches a
  parameter the closures read from their environment moving between calls.
  The arrays are identical either way, so no number moves. The option
  `aci.realiser_cache = FALSE` bypasses the cache (`R/aci-realiser-cache.R`,
  the specification's fingerprint, plan v0.3 PR-7).

* For a scalar hidden state the Theorem 3 auxiliaries (the update factor,
  the gain row and the one-lag increments of every interval) are formed as
  vector expressions over the record (`.online_aux_scalar()`,
  `R/aci-online-scalar.R`) instead of one interval at a time, and feed the
  Theorem 3 smoother, the forward primitives of the range, the online window
  route and the lag table. On the dyad record the smoother falls from 0.17 s
  to 0.004 s and the primitives from 0.15 s to 0.02 s, within one rounding
  of the per-interval kernels (bit-identical where the BLAS's triangular
  solve divides); the covariance policy is reached by the per-interval route
  whenever a variance needs it. Ported from aciR 0.2.3 (`.aci_online_aux()`
  in `R/aci-online-smoother.R`, tag `parents-final`). The matrix path is
  unchanged.

* The forward influence range and the lag table form each anchor's row of
  divergences as one vector expression when the hidden state is scalar
  (`.cir_scalar_row()`, `R/aci-cir-rows.R`), instead of advancing every
  anchor one cell at a time. On the authors' dyad record (N = 3,000, 1,001
  reporting anchors) the range falls from 44 s to under a second and the
  adaptive lag table from 137 s to one second, with every reported quantity
  within 1e-13 of the cell-by-cell recursion and the freeze indices, tail
  estimates and statuses unchanged. Ported from aciR 0.2.3 (`.aci_cir_row()`
  in `R/aci-cir.R`, the cumulative logarithms of `.aci_online_aux()` in
  `R/aci-online-smoother.R`; tag `parents-final`), with one change of
  summation: the logarithms of the update factors are summed from the anchor
  outward in blocks of 512 cells rather than differenced from a record-length
  cumulative sum, which errs by at least the rounding of that sum on any
  platform (6e-13 at 100,000 steps of contracting factors on x86) and by
  1.7e-12 at 20,000 steps and 2.8e-11 at 100,000 where `cumsum` accumulates
  in double (arm64); the blocked form stays below 3e-13. The
  relative entropy of each cell is evaluated in the operations of
  `.kl_fast()`, so a cell is bit-identical to the recursion's given the same
  posterior. The matrix path (`l > 1`) is unchanged.
