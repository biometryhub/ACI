# Scalar numerical corrections and package documentation

Valid finite inputs can produce an infinite scalar KL signal because an
intermediate square or mean difference overflows even when the result fits in
double precision. The dispersion corrections address separate small-ratio and
ratio overflow failures. The batch also includes the plotting argument,
online-window count and documentation corrections.

Preserve the existing source-comparison tolerances and ordinary finite signal
arithmetic. Verify computational cost. No API, dependency, fixture, tolerance
or version change is intended. The slower optional near-equality series is
excluded.
