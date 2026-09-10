# Release acir 0.2.0

The corrections merged in #25 are ready for a numbered release. The package
still identifies itself as 0.1.0.9000, while its external citation metadata
describes 0.1.0 and an older R requirement.

Set DESCRIPTION and the current NEWS heading to 0.2.0. Align CITATION.cff
and codemeta.json with that version and the existing R 4.1.0 requirement.
Correct codemeta's package documentation paths and remove obsolete generated
runtime/size values. Preserve authorship, references and historical versions.

Retain the deprecated simulation argument `T`, its warning class and its
tests. Replace the expired removal promise in its warning, help and test
comment; direct callers to `t_end`. No argument removal is part of this release.

No numerical algorithms, defaults, fixtures, tolerances, exports or workflow
behaviour change. The agricultural simulation remains a separate example.

Validation: build and install a fresh candidate; run R CMD check with tests,
examples and vignettes; check release metadata and installed citation;
compare numerical source and fixture hashes with the merged main snapshot.
The existing alias tests must retain identical draws and reject simultaneous
`T` and `t_end` arguments.

After review and passing CI, use `acir-v0.2.0` on the final main commit.
The repository already has a distinct `v0.2.0` tag; do not replace it.
