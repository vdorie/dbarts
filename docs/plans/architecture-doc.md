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
