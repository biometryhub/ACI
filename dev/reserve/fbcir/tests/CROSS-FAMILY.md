# Cross-family test blocks that also cover FBCIR

Two excised blocks assert FBCIR constructors alongside EnKBS and extremes
ones and could not be split without rewriting them. Both are filed once,
under EnKBS:

- `reserve/enkbs/tests/test-02-zoo-models.R` — Z6-lite; its
  `model_tipping_triad(0.1)` line is the FBCIR part.
- `reserve/enkbs/tests/test-11-other-families.R` — its
  "tipping-triad partitions reproduce both conditional FBCIR questions"
  block and the `variant = "fbcir"` half of the Lorenz-84 block are the
  FBCIR parts.
