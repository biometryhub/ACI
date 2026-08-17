#!/bin/sh
# Render the walkthrough and install it into the site's asset tree.
#
# The rendered page is self-contained, so it can be opened from disk, but its
# two links to the ledgers are relative and only resolve once the page sits
# beside them at the site root. `aciR/pkgdown/assets/` is copied verbatim into
# the built site, which puts `walkthrough.html` at `/walkthrough.html` and the
# ledgers at `/ledgers/`, and the links resolve.
#
# Rendering through this script rather than by hand is what keeps the published
# copy from drifting away from the source, and what stops the destination path
# from depending on remembering it.
#
# Run from the repository root:
#   sh tools/walkthrough/build.sh

set -eu

here=$(dirname "$0")
root=$(cd "$here/../.." && pwd)
dest="$root/aciR/pkgdown/assets"

if [ ! -d "$dest" ]; then
  echo "asset directory not found: $dest" >&2
  exit 1
fi

cd "$root/tools/walkthrough"
quarto render walkthrough.qmd --to html

cp walkthrough.html "$dest/walkthrough.html"
echo "installed: aciR/pkgdown/assets/walkthrough.html"

# The page must stay inside the self-contained size budget: it is served to
# readers on unknown connections and is embedded rather than streamed.
bytes=$(wc -c < "$dest/walkthrough.html" | tr -d ' ')
mb=$(awk -v b="$bytes" 'BEGIN { printf "%.2f", b / 1048576 }')
echo "size: ${mb} MB"
if [ "$bytes" -gt 5242880 ]; then
  echo "FAIL: over the 5 MB budget." >&2
  exit 1
fi
