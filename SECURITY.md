# Security policy

## Supported versions

The most recent release is supported. Fixes are issued as a new patch release
rather than as patches to earlier versions.

| Version | Supported |
|---|---|
| 0.2.x | yes |
| 0.1.x | no |

## Reporting a vulnerability

Report privately through GitHub's advisory form at
<https://github.com/biometryhub/ACI/security/advisories/new>, which opens a
channel visible only to the maintainers. Please do not open a public issue for
a vulnerability.

A report is most useful with the affected version, the R version and platform,
the steps that reproduce the problem, and what an attacker could achieve.

You can expect an acknowledgement within seven days and an assessment within
thirty. If a report is accepted, the fix, the release carrying it, and the
credit you would like are agreed with you before anything is published.

## Scope

This package computes on numeric input supplied by the caller. It opens no
network connections, reads no files unless the caller names one, and executes
no code supplied as data. The realistic exposure is therefore narrow, and a
report that identifies something outside that expectation is valuable
precisely because it contradicts it.

Reports of numerical error are handled as defects rather than as
vulnerabilities. See [`aciR/SUPPORT.md`](aciR/SUPPORT.md), which explains what
to include and why those reports are the most useful ones this project
receives.
