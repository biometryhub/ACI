# Validation and the independent oracle

## The problem this addresses

A test suite written by the author of the code it tests can only prove
the code is self-consistent. If the author misread an equation, the test
encodes the same misreading and passes. For a package whose entire
purpose is to reproduce a published numerical method, self-consistency
is close to worthless. The question is not whether aciR agrees with
itself, but whether it agrees with the method.

Every numerical claim aciR makes is therefore graded against a reference
the package did not produce. This article records what those references
are, how they were generated, what each one covers, and (the part most
easily overlooked) what each one does *not* cover.

``` r

library(aciR)
```

### What an oracle is

This article uses “oracle” and “grounded” throughout, and neither is
everyday vocabulary.

A test oracle is the mechanism that decides, for a given input, whether
a program’s output is correct (Barr et al., 2015,
[doi:10.1109/TSE.2014.2372785](https://doi.org/10.1109/TSE.2014.2372785)).
Most software has one close to hand, whether a specification, a worked
example or an answer that is obviously right. A program written to
compute a quantity that is not otherwise known has none, and Weyuker
(1982,
[doi:10.1093/comjnl/25.4.465](https://doi.org/10.1093/comjnl/25.4.465))
called such a program non-testable. The usual recourse is a
pseudo-oracle, a second implementation of the same specification,
written independently, whose disagreement with the first exposes an
error in one of them (Davis and Weyuker, 1981,
[doi:10.1145/800175.809889](https://doi.org/10.1145/800175.809889)).

We use “grounded” for something narrower, an oracle whose expected
values come from a source that did not also produce the code under test.
The term is ours rather than the literature’s, and the distinction it
draws is the one the rest of this article turns on.

## The three oracles

aciR uses three, because no one of them is sufficient.

The quantity all of them grade is, at each step, the relative entropy of
the smoother posterior $`\mathcal{N}(m_s, R_s)`$ from the filter
posterior $`\mathcal{N}(m_f, R_f)`$, which for scalar Gaussians has the
closed form

``` math
\mathrm{KL}\!\left(\mathcal{N}(m_s, R_s) \,\|\, \mathcal{N}(m_f, R_f)\right)
= \frac{1}{2}\left[
    \frac{(m_s - m_f)^2}{R_f}
    + \frac{R_s}{R_f} - 1
    - \log\,\frac{R_s}{R_f}
  \right],
```

together with the four posterior series that feed it. The first term is
the signal part, driven by the shift between the two means, and the
remainder is the dispersion part, driven by the covariance ratio;
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)
evaluates the dispersion part in a cancellation-resistant form, which is
one of the things the oracles hold it to.

### 1. The dyad fixture, from the authors’ reference implementation

The flagship oracle. A MATLAB harness reproduces the deterministic core
of `dyad_interaction_model.m` from the reference implementation at
[github.com/marandmath/ACI_code](https://github.com/marandmath/ACI_code)
(MIT licence, commit `733c49f`), using that reference’s own seed and
parameters. It writes two files. The first holds the observed signal.
The second holds the expected filter mean, filter covariance, smoother
mean, smoother covariance and ACI metric at every hundredth step.

aciR then runs its *own* filter, smoother and metric on that signal and
must reproduce all five series. The comparison is at 301 sampled indices
of a 30,001-step trajectory, and the gate is a maximum absolute error
below `1e-6`.

The observed agreement is far tighter than the gate.

| Oracle      | Maximum absolute error | Gate   |
|-------------|------------------------|--------|
| Dyad        | `4.574119e-14`         | `1e-6` |
| Cross-noise | `1.376677e-14`         | `1e-6` |

Both sit at the floating-point noise floor of a 30,001-step recursion.
This is not agreement within tolerance; it is agreement to the last few
bits the arithmetic can carry.

### 2. The cross-noise fixture, and why it had to exist

The dyad oracle is strong, but its scope is narrower than it first
appears.

Every scalar model in the reference implementation, the dyad and the
noisy predator-prey model alike, sets the noise decomposition so that
the noise cross-covariance is zero. The filter’s update carries a term
`aux = S_yoS_x + R * L_x`, and the smoother carries `A_j` and a
transport term, all of which involve that cross-covariance. On every run
of the dyad fixture those terms *execute*, and are multiplied by zero.

Two measurements therefore say less than they appear to. Line coverage
sees the arithmetic run and reports it covered. The oracle sees the
outputs agree and reports them graded. Neither measurement is wrong, and
neither one reaches the cross-covariance terms, which
`aci_cgns_model(S_yoS_x = ...)` exposes to users as a supported path.
Coverage of the graded path was therefore standing in as evidence about
an ungraded path beside it.

aciR ships a second harness for that path, `aci_oracle_cross.m`, which
is the reference dyad with correlated noise switched on. It has no
upstream counterpart, because upstream has no scalar model that
exercises these terms. It is therefore a second independent
implementation of the published equations in a different language and
runtime. Agreement with it refutes a transcription error on the R side
over the full transient. It is *not* an authors’-reference grounding,
and the manifest says so.

### 3. The analytic identities, which depend on no implementation

The primary grounding for the cross-noise path is not either fixture. It
is algebra.

For a constant-coefficient system the stationary filtered covariance is
the positive root of an algebraic Riccati equation, obtained by setting
the covariance derivative to zero.

``` math
L_x^2 R^2 \;-\; 2\left(L_y S_{xx} - S_{yx} L_x\right) R
\;-\; \left(S_{yy} S_{xx} - S_{yx}^2\right) \;=\; 0 .
```

Started at that root with a non-zero cross-covariance, the recursion
must stay there to machine precision, because the fixed point of the
differential equation is exactly the fixed point of its Euler map. This
oracle is derived from the governing equations and is independent of
aciR *and* of MATLAB. The package’s tests also check the smoother’s
analogous fixed point, the terminal identity, and the zero-information
limit.

The zero-information test is constructed with some care. With no
coupling and no cross-covariance the observed process is independent of
the unobserved one, so the smoother must equal the filter and the metric
must be zero. That identity is exact in continuous time only. Under an
explicit Euler scheme the forward and backward sweeps do not invert one
another exactly, so an assertion of equality would test a property of
the mathematics against an approximation of it. The test instead asserts
that the residual vanishes at second order in the step. That is the
claim the discretisation supports, and it is the stronger of the two,
because it checks the *rate* rather than only the size.

## Provenance and the manifest

Trusting a committed fixture requires knowing where it came from and
being able to tell an intentional change from an accidental one.
`inst/extdata/oracle-manifest.yml` records, for each fixture, the
upstream repository and immutable commit, the licence and copyright
holder, the MATLAB release, the exact generating command, the parameters
and seed, the sampled indices, both an MD5 and a SHA-256 hash, the
measured error, and a `grades` and `does_not_grade` field.

``` r

manifest <- system.file("extdata", "oracle-manifest.yml", package = "aciR")
cat(paste(readLines(manifest)[1:12], collapse = "\n"))
#> # Provenance of the independent-oracle fixtures shipped with aciR.
#> #
#> # Each fixture pair is (signal, reference). aciR runs its OWN filter, smoother
#> # and metric on the signal and must reproduce the reference at the sampled
#> # indices to the stated tolerance. The tests that enforce this are the
#> # tests/testthat/test-oracle-*.R family, one per fixture group;
#> # tests/testthat/test-oracle-manifest.R checks the shipped bytes against the
#> # hashes below, so a changed fixture cannot pass unnoticed, and
#> # tests/testthat/test-grading-matrix.R fails the build when a capability is
#> # neither graded nor declared ungraded.
#> #
#> # Read `grades` on each entry before citing a fixture as validation. An oracle
```

The recorded MATLAB release is grounded by reproduction rather than by
recollection. The harness was rerun on that release into a scratch
directory, and both files reproduced byte-for-byte, with identical
hashes to the committed copies. A test checks the shipped bytes against
the recorded hashes on every run, so a fixture cannot drift silently.

## What is graded, and what is not

An oracle grades the paths its parameters exercise and no others. The
two lists below separate the paths by the source of their evidence
rather than by whether they work, because the distinction that matters
here is whether the evidence comes from the authors. A third section
names a boundary that belongs on neither list.

### Graded, with the scope named

- **A time-varying self-drift**, in both causal directions, against the
  noisy predator-prey model, the system whose latent damping is set by
  the observed population. The dyad and cross fixtures cannot reach this
  path, since both pin the self-drift at a constant.
- **The whole predator-prey path, filter through causal influence
  range**, in both causal directions, against byte-verified extracts of
  the reference script. Twenty-six graded quantities, worst absolute
  difference `1.50e-12` and tightest headroom 18.8 times the declared
  tolerance. Two further comparisons report aciR’s default horizon
  against the reference’s, which is a designed difference rather than a
  disagreement, and they are excluded from those two figures. The
  subjective ranges agree closely enough that the two implementations
  select the same index at all 616,113 threshold-time cells in each
  direction, so what remains is round-off in converting counts to time.
- **The noise cross-covariance of the online smoother and the causal
  influence range**, against the `cir_cross` fixture. This is a second
  independent implementation rather than the authors’ own numbers. Every
  scalar model upstream sets the cross-covariance to zero, so no
  reference fixture can reach these terms.
- **Vector filter, smoother and metric**, against the ENSO fixture and
  against block-diagonal algebraic Riccati fixed points solved in closed
  form, which depend on neither implementation.
- **The vector online smoother and causal influence range**, by collapse
  onto the graded scalar path at every lag and by an independent
  transcription of equations (3.5) to (3.7) in full generality. Not
  against the authors’ ENSO scripts.
- **The Euler-Maruyama integrator**, against the reference’s own
  3001-step dyad capture. Path equality at a matching seed is
  impossible, since R and MATLAB draw normal variates by different
  algorithms, so that is not what is graded. The scheme is invertible.
  Subtract the drift from a captured transition and divide by the noise
  coefficient, and the variate that produced it comes back, so the
  reference’s own draws are recoverable from its own path. Driven with
  them,
  [`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
  reproduces the observed component to `5.1e-15` and the unobserved to
  `1.8e-15`, against a round-off bound of `1.15e-10` derived from the
  arithmetic before the comparison was run. What this grades is the
  update rule, and only that. The drift the recovery subtracts is the
  drift the simulator adds back, so an error common to both would cancel
  here and the path would still be reproduced. That common error is what
  the dyad oracle catches, by grading the filter against the reference’s
  own output.

The matrix noise-cross-covariance terms need their scope named
separately. They have **no upstream counterpart at all**. Every scalar
model in the reference sets the cross-covariance to zero, and the ENSO
scripts state in as many words that the cross-interaction terms are
absent. The analytic identities are the primary grounding there, and a
second implementation refutes a transcription error without being an
authors’-reference grounding.

### Not yet graded

One item, and it is work not yet done rather than work that cannot be
done.

- **The aciR side of the four scalar-latent ENSO configurations,
  including
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md).**
  The reference side is complete. All four scripts are extracted and
  gated, at a maximum absolute difference of zero, and three of the five
  ship with conditional ACI enabled, so a reference to grade
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
  against does exist. Two things stand between that and a comparison.
  Firstly, the harness exports its bundles in a square configuration
  with feedback-matrix noise, and these are three observed against one
  latent, or two against one, with scalar latent noise. The package
  itself accepts that shape, so the work is in the export rather than in
  aciR. Secondly, the reference leaves the first filter step of a
  conditional run unmasked where aciR masks every step, which moves the
  filtered mean by 0.238 at step two. That convention has to be matched
  deliberately, or the comparison will report a difference that is not a
  disagreement.

### A boundary, not a gap

Model adequacy is not on that list, and it never will be. Every oracle
in this article asks one question. Does the code compute the method it
claims to compute? Whether that method describes your data is a separate
question, and it belongs to the science rather than to the software. No
fixture, no identity and no reference implementation can reach it, and
none here pretends to. The companion article *Assumptions and
interpretation* is where that question is treated.

## Reproducing the grade

The oracle tests run from the installed package fixtures. They never
skip, on any platform, in any check environment. A validation gate that
can skip is not a gate, and an earlier version of these tests could skip
when run outside the development tree.

``` r

testthat::test_local(system.file(package = "aciR"))
```

Regenerating the fixtures needs MATLAB and the harnesses in the
project’s `tools/oracle/` directory. The commands are in the manifest.
