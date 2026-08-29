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

Cite a durable landmark: name the symbol, section title, or item label
the change would move with - `rollTreeResidual (chain.hpp)`, not
`chain.hpp:2581`; `data-ownership.md "Implementation record" item 5`,
not `data-ownership.md:180`. Add a `file:NNN` line number alongside the
symbol as a locating convenience, never in its place.
`tools/check-doc-freshness.R` verifies docs/design anchors (both the
file:line and the cited symbol) and the docs/design and docs/plans
INDEX manifests; it does not check docs/plans citations themselves, so
name the symbol as well as the line - a moved anchor there is caught
only by re-reading, not by the checker. State each fact in one home
doc; elsewhere link to it by title, do not restate it (a copied fact is
a second thing to keep in sync, and the one that rots).

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
