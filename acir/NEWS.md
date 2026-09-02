# acir 0.0.0.9000

## Performance

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
