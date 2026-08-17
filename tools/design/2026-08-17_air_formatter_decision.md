# The Air formatter: evaluated, and not adopted

Date: 2026-08-17
Status: closed

## What Air is

Air is an R source-code formatter published by Posit, written in Rust. It is
the R equivalent of `gofmt` for Go, `black` for Python, or `prettier` for
JavaScript: a program that reads an R file and rewrites its whitespace, line
breaks and indentation into one canonical shape.

The point of such a tool is not that its chosen shape is the best one. The
point is that it removes formatting from the set of things anyone has to think
about or argue over. Every contributor's editor writes the same layout on save,
so a diff never contains a reindentation, and a review never spends a comment
on where a line should wrap.

Air formats layout only. It does not rename anything, reorder anything, change
control flow, or alter what the code computes. Its configuration surface is
deliberately small: a line width, an indent style and width, and little else.
That narrowness is the design, and it is why the tool is fast and predictable.

## Why it was considered here

The AAGI-AUS repository guidelines name Air as their current direction, and
several packages in that organisation reference it in their contributing
instructions. This package is destined for `biometryhub` rather than for
`AAGI-AUS`, so the guideline does not bind, but it is worth measuring against.

## What the measurement showed

Running `air format --check R/` against the package reported that Air would
reformat **15 of the 16 R source files**.

That result is the decision. The reformatting would not fix anything: the
package already passes `lintr` with zero lints under a configuration that
enforces the same 80-character width and two-space indent that Air would apply.
The difference is not correctness or consistency. It is that Air makes
different choices from the ones this code makes deliberately, chiefly in how it
breaks long argument lists and how it aligns the continuation of a wrapped
expression.

Adopting it would therefore produce a single commit touching almost every
source file, with no behavioural change and no defect fixed. That commit would
sit across the history at exactly the point where a reader tracing the origin
of a line would most want `git blame` to be informative. The cost is paid once
and then paid again by every future reader of the history.

## The decision

Air is not adopted, and this is closed rather than deferred.

The condition that would reopen it is a change in who writes the code. Air
earns its cost when several people commit to the same files, because then it
replaces a recurring coordination problem with a one-off reformatting. With a
single author and a clean `lintr` gate already in place, there is no
coordination problem for it to solve.

If this package later acquires regular contributors, the sequence to adopt it
is: add `air.toml` at the package root, run `air format R/ tests/` as one
isolated commit that changes nothing else, add that commit's hash to
`.git-blame-ignore-revs` so `git blame` skips over it, and note the convention
in `CONTRIBUTING.md`. The `.git-blame-ignore-revs` step is what makes the cost
above recoverable, and it is the step most often forgotten.

## What replaces it

`lintr`, configured in `aciR/.lintr` and run in CI by
`.github/workflows/lint.yaml`, which enforces the line width, the indent, the
naming convention and the absence of the constructs this project avoids. It
currently reports zero lints across the package. The distinction is that
`lintr` reports on layout and style without rewriting the source, which suits a
codebase where the layout is hand-set and intentional.
