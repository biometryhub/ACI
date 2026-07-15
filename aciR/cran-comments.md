# cran-comments

## Release summary

aciR 0.1.0 is the first research preview. The package is currently distributed
from a private repository and is **not being submitted to CRAN at this
version**; this file records the release-check state so that the notes are
adjudicated when they arise rather than at submission time.

## Test environments

* local macOS 26.5.2 (Tahoe), R 4.5.2 — `R CMD check --as-cran`, vanilla
  session, run on the built tarball
* GitHub Actions (configured, runs on push): ubuntu-latest (devel, release,
  oldrel-1, and 4.1 — the declared floor), macOS-latest (release),
  windows-latest (release), with `error-on: warning`

## R CMD check results

0 errors | 0 warnings | 2 notes

### Note 1 — new submission and unreachable URLs

```
Maintainer: 'Max Moldovan <max.moldovan@gmail.com>'
New submission
Found the following (possibly) invalid URLs:
  URL: https://github.com/max578/aciR  ... Status: 404
```

Expected and accurate. The declared repository is private, so its URL returns
404 to an anonymous checker. The URL is correct and is the canonical home of
the package; it is not reachable without authentication. This note will persist
for as long as the repository is private, and would be resolved -- not
suppressed -- by making it public, which is a decision for a later release
rather than something to work around in metadata.

### Note 2 — HTML manual validation skipped

```
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
```

Environmental. The local HTML Tidy is older than R's requirement; the note
concerns the checking machine, not the package. The CI matrix checks the manual
on current images.

## Notes deliberately not present

* The 0.0.0.9000 development version carried a "large version component" note.
  Resolved by the 0.1.0 version.
* The earlier check ran without the PDF manual. The manual now builds; the
  remaining note above is only about HTML validation tooling.

## Test suite behaviour worth stating

`R CMD check` reports 7 skipped tests. All 7 are the `expect_snapshot()` blocks
in `test-errors.R`, which testthat skips on CRAN by default because snapshot
output can vary with R version and locale. Those snapshots lock the *wording* of
error messages; the *conditions* that raise them are tested by ordinary
`expect_error()` assertions in `test-validate.R`, which run everywhere. No
contract is unenforced on CRAN.

The independent-oracle tests never skip, on any platform or environment. This
was verified explicitly by running them with `NOT_CRAN=false`: 87 assertions,
0 skipped. A validation gate that can skip is not a gate, and an earlier
version of these tests could skip when run outside the development tree; that
test has been removed and replaced with one that runs from the installed
package fixtures.

## Downstream dependencies

None. This is a new package with no reverse dependencies.
