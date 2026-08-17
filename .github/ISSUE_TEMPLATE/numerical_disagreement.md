---
name: Numerical disagreement
about: A quantity this package computes differs from what you expected, or from another implementation
title: ''
labels: numerical
---

These are the most useful reports this project receives, and they are treated
as findings rather than as support requests. Conventions differ between
implementations in ways that look like disagreements and are not, so naming
the settings is what makes a report actionable.

**The system and its parameters**

<!-- Which model, and any parameter you changed from the default. -->

**The quantity that disagrees**

<!-- For example: the filter covariance, the objective range, the metric at a
     particular time. -->

**What you got, and what you expected**

<!-- The size of the disagreement matters more than the values alone. If you
     compared against another implementation, say which and at what settings. -->

**A reproducible example**

```r
# aci_simulate() takes a `seed`, so a simulated case can be made reproducible
# in one line.
```

**Session information**

<details>

```r
# output of sessionInfo()
```

</details>

---

Before filing, two conventions account for most apparent disagreements, and
both are documented:

- `objective` and `objective_exact` from `aci_cir()` are different
  functionals, not two quadratures of one. They coincide only where the
  divergence decreases with lag, which the `monotone` column reports.
- A range reported near the end of a record is a lower bound. Check the
  `status` column for `"censored"` before comparing it against anything.
