# Superseded revisions of the development ledger

**None of the files in this directory is current.** The published development
ledger is at <https://biometryhub.github.io/ACI/ledgers/development_ledger.html>,
and its source is `aciR/pkgdown/assets/ledgers/development_ledger.html`.

These three are earlier renders, kept because the document was rewritten
substantially rather than edited, and the intermediate states record what the
review rounds changed. They are not duplicates of what `git log` holds: each is
a distinct state that was never committed as the live ledger.

| File | What it was |
|---|---|
| `aciR_development_ledger_v1.html` | First render. Structure only, before the review rounds. |
| `aciR_development_ledger_v2.html` | After the first review. Later found to break several of the documented voice rules, which is what prompted v3. |
| `aciR_development_ledger_v3.html` | After the voice pass. Superseded by the current revision, which reworked the closing sections from a defensive register into open problems with routes. |

They carry the same `<title>` as the live document, because they were rendered
from the same template and have not been edited since. Read the filename, not
the title.

Two of the corrections that landed across these revisions are worth naming,
since a reader comparing them will notice the difference:

- The count of callable files in the reference implementation was wrong in the
  earliest text. There are three, not two, and the third is the authors' own
  rather than third-party. The correction is present from v3 onward.
- The claim that the recovered variates' normality evidenced a correct drift
  was withdrawn. It does not: at the integration step used, a constant drift
  error displaces each variate too little to show up in that test.

If a claim in one of these files disagrees with the current ledger, the current
ledger is the one that has been checked.
