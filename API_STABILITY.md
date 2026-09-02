# API stability

acir is at version 0.x. What that means for a user:

- The numerical core is fixed. Every graded quantity reproduces the
  method authors’ reference implementation to the tolerance recorded in
  `inst/evidence/register.csv`, and a change that moved a number beyond
  round-off would be a change of method, not a release.
- The public interface may still change before 1.0: argument names, the
  shape of returned objects, and the names of the model coefficient
  functions. Every such change is announced in `NEWS.md` under a
  “Breaking changes” heading, with the old form kept working for one
  minor version where that is possible.
- From 1.0, the interface follows semantic versioning: breaking changes
  only at a major version, with a deprecation period.

Bugs and numerical disagreements:
<https://github.com/biometryhub/ACI/issues>.
