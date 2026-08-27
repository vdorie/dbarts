# architecture-doc

agent: sonnet
rng: neutral (docs only)
budget: one new doc, ~300 lines

## Goal

A contributor who has never read the design docs' history can learn how
the engine works today from one document: docs/architecture.md.

## Context

- The design record (docs/design/) is proposal + landing notes -
  correct as history, poor as orientation; core-generalization.md's
  768 lines are chronological.
- readability-review step 3 adds a short current-architecture section
  to core-generalization.md; this doc is the full-length version and
  should absorb/replace that section if both land (coordinate).
- Source material: core-generalization.md (layers, dispatch tiers,
  concept decomposition, BartData), kernel-vocabulary.md,
  public-surface.md section 6, CLAUDE.local.md's layout paragraph.

## Constraints

- Present tense, current state only; history stays in the design docs
  (link, do not duplicate).
- Every claim checked against the code as written, not the docs (the
  docs describe some aspirations the code defers - per-column widths,
  arena allocation; the architecture doc must not repeat them as fact).
- Out of scope: user-facing documentation (vignette-refresh).

## Steps

1. Outline: package layering (R API -> R5 -> bridge -> facade ->
   engine -> kernel libs); dispatch tiers table; the concepts
   (ResponseModel, LeafModel, MoveStrategy, SplitSelector) with their
   shipped implementations; ColumnStore's actual layout; the mutation
   surface and its transaction semantics; tree storage (live, flat,
   state); RNG architecture; threading model; the gates
   (component/tinytest/equivalence/bench) and when each runs.
2. Write against the source, citing file paths not line numbers (lines
   rot); one diagram-as-ASCII for the layering.
3. Cross-link from CLAUDE.local.md's Layout section and the design
   docs' headers.

## Verification

- VD skims for wrongness; a fresh agent given only architecture.md can
  correctly answer five orientation questions (where do moves live,
  what is a facade, which layer owns cut points, what gates a
  sampling change, where does dispatch happen) - run that check.

## Status (2026-07-08)

Landed. docs/architecture.md (360 lines vs the ~300 budget; reviewed
line-for-line as load-bearing): layering diagram with the dbarts.h
side path, the four dispatch tiers, leaf concepts vs the runtime
ResponseFamily switch (incl. the GroupedResponse decorator), tree
moves + the rule/mask machinery, ColumnStore's real layout, the
transactional mutation surface, the three tree storage forms, RNG and
threading (incl. the callback restriction), and the gates table. The
implementer verified claims against code and correctly refused three
doc-borne aspirations (per-column code widths, arena allocation, the
MoveStrategy/SplitSelector names, which have no code symbols).
core-generalization.md's current-architecture section trimmed to a
pointer; kernel-vocabulary.md and public-surface.md gained pointer
lines. The five-orientation-questions check ran against a fresh
reader given only the doc: all five answered correctly. One review
fix: a "header-first" thinko about misc.a reworded to "stays a
compiled library". CLAUDE.local.md's cross-link is left to VD (his
private file). AWAITING: VD's skim for wrongness.
