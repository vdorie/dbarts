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
cap on the final file. Front block:
`agent:`, `rng:`, `window:` (if any), `budget:` (expected diff size).
Sections:

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
another document - is written in one machine-checked grammar, delimited
by double brackets so a single regex finds all of them, with `#`
separating the path from what it names so a C++ `::` qualification never
collides with the path:

    [[src/bartcore/chain.hpp#rollTreeResidual]]     current state, by symbol
    [[MOD#GaussianResponse]]                        alias path
    [[R/bart.R#bart2()]]                            a trailing () is stripped
    [[sampler.Rd#dbartsSampler$setResponse]]        split on $ as well as ::
    [[test-argument-surface.R#"hurdle.lognormal"]]  verbatim fragment
    [[docs/design/data-ownership.md#Implementation record]]  doc to doc, by heading
    [[src/bartcore/chain.hpp:2581@d477a46b]]        history, a line at a commit
    retired: [[R/spec.R#refuseHostMutation]]        skipped; prose says it is gone

- A path is a repo path, a basename unique in the tracked inventory, or
  one of the aliases the head of feature-matrix.md defines. It must
  resolve to a file that exists: a citation into another repository is
  prose, not a cite.
- A symbol is checked component-wise - split on `::` and `$`, each part a
  whole-word token somewhere in that file. Adjacency is not required and
  no line number is involved, so a definition that moves needs no re-pin.
  Several names in one cite are separated by commas.
- A `.md` target is always a heading cite: the text must appear in one of
  the target document's markdown headings.
- A quoted fragment is checked as a literal substring. That is the form
  for a file with no symbols to name - a tinytest script; a tests/cpp
  citation names its test function instead.
- A line number appears ONLY in the history form, `path:line@sha` or
  `path:a-b@sha`, whose sha is at least 8 hex digits and must be an
  ancestor of HEAD. Ancestry is all that is checked, because the line
  belongs to that commit rather than to the working tree. Landing notes
  and a plan's Landing section cite this way, pinned to the commit the
  note describes.
- `retired:` immediately before a cite skips it, and the prose around it
  must say the thing is gone.

Backticks around a cite are optional. `tools/check-doc-freshness.R`
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

Append the plan's `## Landing` note, and bump the matching
`docs/design/<x>.md` `Status:` line to `LANDED <date> (<commit>)`. The
design Status line is the record most often missed; check it
explicitly at every landing.
