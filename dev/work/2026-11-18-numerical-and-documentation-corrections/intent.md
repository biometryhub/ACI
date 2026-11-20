# Scalar numerical corrections and reviewed package documentation

The author approved the 20 default items in the 18 November package review and
requested correction of the remaining scalar signal overflow. Valid finite
inputs can produce an infinite scalar KL signal because an intermediate square
or mean difference overflows even when the result fits in double precision.
The approved dispersion corrections address separate small-ratio and ratio
overflow failures. The batch also includes the reviewed plotting argument,
online-window count and documentation corrections.

Preserve the existing source-comparison tolerances and ordinary finite signal
arithmetic. Verify computational cost before staging. No API, dependency,
fixture, tolerance or version change is intended. The slower optional
near-equality series requires the supervisor's decision and is excluded.

The author permits staging but has not authorised a commit, push or message.
The paper, raw inputs and prior evidence remain unchanged. Review artifacts,
test outputs and the deferred alternative are stored outside this checkout in
the ACI-review workspace under
outputs/implementation-v2/20261118-package-audit/signal-fix-adoption-20261118/.
