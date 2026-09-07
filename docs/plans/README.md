# Plan process

The TODO at the repo root is an unordered backlog; most items name a
plan file here, a few record instead in a design doc or a
differently-named plan. This document covers the plan-file format, the
citation rule, the RNG-class gates a change needs, and where a landing
gets recorded.

## Plan format

One file per TODO item, `<item>.md`. A small decision-gated stub
targets roughly 80 lines; a landed or multi-slice item's file commonly
runs to several hundred lines (occasionally thousands) once Landing
notes, amendments, and session logs accumulate - there is no enforced
cap on the final file. Front block, should carry:
`agent:`, `rng:`, `window:` (if any), `budget:` (expected diff size).
Sections, should carry:

- Goal: two or three sentences; what is true after the item lands.
- Context: durable references (a code symbol name, a doc section title
  or item label) and design-doc pointers; no narrative a pointer can
  replace. Never a bare line number - see Cross-references.
- Decision (decision-gated items only): the question, a recommendation,
  and what evidence would change it. VD signs off before implementation.
- Constraints: gates, contract freezes, explicit out-of-scope list.
- Steps: numbered; each independently verifiable.
- Verification: exact commands and expected outcomes.

## Cross-references

Documentation cites code by SYMBOL. Every citation - of code or of
another document - is written in one machine-checked grammar: an
ordinary markdown link, so it renders and resolves on GitHub. The target
is the file, written relative to the citing document's own directory;
the text is what in that file is meant:

    [`rollTreeResidual`](../../src/bartcore/chain.hpp)  current state, by symbol
    [`GaussianResponse`](../../src/bartcore/model.hpp), [`ProbitResponse`](../../src/bartcore/model.hpp)  several symbols, one cite
    [`bart2()`](../../R/bart.R)  a trailing () is stripped
    [`dbartsSampler$setResponse`](../../man/dbartsSampler-class.Rd)  split on $ as well as ::
    ["hurdle.lognormal"](../../inst/tinytest/test-argument-surface.R)  verbatim fragment
    [Implementation record](../design/data-ownership.md#implementation-record)  doc to doc, by heading slug
    [src/bartcore/chain.hpp:2581](https://github.com/vdorie/dbarts/blob/d477a46be658d885b769dc186c2832991348187d/src/bartcore/chain.hpp#L2581)  history, a line at a commit
    retired: [`refuseHostMutation`](../../R/spec.R)  the construct is gone
    unresolved: [inst/include/dbarts/R_C_interface.hpp:40](https://github.com/vdorie/dbarts/blob/d477a46be658d885b769dc186c2832991348187d/inst/include/dbarts/R_C_interface.hpp#L40)  cannot be placed

- A target is a path relative to the citing document, and it must
  resolve to a file that exists: a citation into another repository is
  prose, not a cite. A link whose text is the target's own path or
  basename is a plain file link, not a cite; any OTHER link into a
  tracked file fails, because it reads as a cite and is checked by
  nothing.
- A symbol is checked component-wise - split on `::` and `$`, each part a
  whole-word token somewhere in that file. Adjacency is not required and
  no line number is involved, so a definition that moves needs no re-pin.
  Several symbols of one file are several links separated by `, `; that
  run is one cite, and a marker in front of it covers all of it.
- A `.md` target carrying a `#` fragment is a heading cite: the fragment
  must be GitHub's slug of one of the target document's headings -
  lowercased, everything but letters, digits, spaces, hyphens and
  underscores dropped, spaces to hyphens - and the link text is that
  heading, so a renamed section breaks the cite instead of drifting
  quietly. A heading of the citing document itself is `(#slug)`, with no
  path.
- A quoted fragment is checked as a literal substring. That is the form
  for a file with no symbols to name - a tinytest script; a tests/cpp
  citation names its test function instead. A bracket inside the fragment
  is backslash-escaped so the link still parses.
- A line number appears ONLY in the history form, a link to
  `https://github.com/vdorie/dbarts/blob/<sha>/<path>#L<a>-L<b>` whose
  sha is the full 40 hex digits and must be an ancestor of HEAD, and
  whose text repeats `path:a-b`, the claim the URL makes. The claim is
  read AT THAT COMMIT: the path must exist in that commit's tree and the
  cited line must be within the file's length there. A file deleted since
  therefore still cites cleanly, while a line that never existed at the
  named commit does not - which is what catches a bare `:NNN` whose file
  was guessed wrong. Landing notes and a plan's Landing section cite this
  way, pinned to the commit the note describes.
- `retired:` before a cite says the named CONSTRUCT is gone, so its
  content is not checked; the location must still be real - the target
  resolves, and a history cite's sha, path and line are checked as usual.
  The prose around it must say the thing is gone.
- `unresolved:` before a HISTORY cite says the LOCATION itself cannot be
  established: the target of a frozen record's line reference could not
  be placed at any candidate commit. The link must still parse and the
  sha must still be an ancestor; the path and the line are not checked.
  It is the honest marker for a record that cannot be repaired, never a
  way to quiet a cite that can be.
- A marker binds only when it sits immediately before the link, so
  `un-retired:` and any other word ending in `retired:` disarms nothing.
- A cite parked inside a fenced code block is not checked at all: a fence
  is a code sample or a transcript. Do not put a live cite in one.
- A symbol is found by token search over the whole file, so a name
  occurring only in a comment or a string satisfies its cite. Qualify a
  name the target file defines more than once - `ProbitResponse::computeLogLikelihood`,
  not `computeLogLikelihood` - so the token that answers is the one the
  sentence means.

Backticks belong inside a symbol cite's link text and nowhere else in
the grammar. `tools/check-doc-freshness.R`
checks every cite in everything under docs/, the top-level README.md,
man/*.Rd and vignettes/*.Rmd, and fails on any bare line reference left
outside the history form, so an unconverted citation cannot pass
unnoticed. State each fact in one home doc; elsewhere link to it by
title, do not restate it (a copied fact is a second thing to keep in
sync, and the one that rots).

## RNG classes and their gates

- neutral: draws unchanged.
  Gates: tests/cpp component tests; full tinytest suite.
- shifting: draws change, the posterior does not.
  Gates: the above, plus regenerate RNG-locked snapshots by replaying
  whole test files, re-record the equivalence baseline, and pass the
  statistical (z) mode against the previous baseline.
- posterior-changing: the stationary distribution or a default changes.
  Gates: all of the above, plus the exact-posterior gates
  (.github/workflows/exact-gates.yaml's per-family list) and a design
  note in docs/design/.

Hot-path changes of any class additionally need bench-sampler.R compare
on a quiet machine (maintainer-run; never concurrent with other load).

## Landing

Append the plan's `## Landing` (or `## Landing note`) note, and bump the matching
`docs/design/<x>.md` `Status:` line to `LANDED <date> (<commit>)`. The
design Status line is the record most often missed; check it
explicitly at every landing.
