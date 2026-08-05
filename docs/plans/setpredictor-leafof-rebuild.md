# setpredictor-leafof-rebuild

status: OPEN (found 2026-08-04 by the constant-leaf-fits bench
        discharge, docs/plans/constant-leaf-fits.md compare record)
agent: opus
rng: neutral expected (leafOf is a derived cache; draws must not change)
budget: memo first, then ~100-200 engine lines

## Goal

setPredictor-accept returns to within the 5% gate of the pre-compaction
baseline. The constant-leaf-fits refactor made accepted predictor
mutations pay a full O(n * trees) leafOf rebuild (its own landing note
flags this), measured +22-28% on setPredictor-accept-n1000-t75 (x86,
reproducible, isolated to the refactor by an attribution run). This is
the embedded-Gibbs hot path - the package's distinguishing feature - so
the regression matters to per-sweep-mutation consumers (IRT-style
samplers, stan4bart).

## Context

- docs/plans/constant-leaf-fits.md: the leafOf/muByTree design, the
  per-tree leafOfStale lazy-rebuild mechanism (built for
  sampleTreesFromPrior), and the landing note conceding mutation keeps
  the full rebuild.
- The mutation transaction (docs/design/data-store.md): accepted
  setPredictor repartitions changed trees; a tree whose partitions did
  not change needs no leafOf rebuild, and one that changed locally may
  need only a subtree's rows remapped.

## Decision (memo first)

Measure where the +25% goes before building: (a) leafOf rebuild proper,
vs (b) rebuild running eagerly where the old code amortized, vs (c)
double work (rebuild plus a later stale-flagged rebuild). Then pick:
mark-stale-and-rebuild-lazily (reuse leafOfStale; cheapest, but the
next run() pays it), rebuild only trees whose partitions changed, or
per-subtree incremental remap. Evidence that would change the ranking:
if (a) dominates and most trees change on a typical accept, only the
incremental remap helps and the item should be weighed against its
complexity.

## Constraints

- Draws bitwise-unchanged (equivalence trio identical); leafOf equality
  against a freshly derived map stays asserted (the existing 25-sweep
  test).
- No public-surface change.
- Out of scope: the reject path (measured fine), sparse-column
  mutation internals.

## Steps

1. Memo: instrument the accept path on the bench box, attribute the
   +25%, pick the mechanism, record here.
2. Implement per the memo; extend the leafOf-consistency test with an
   accepted-mutation round.
3. Same-machine A/B on dbarts-bench vs 9047e05 and vs the pre-fix tip:
   setPredictor-accept within the gate, run metrics unchanged.

## Verification

- Gate battery per CLAUDE.local.md (install --preclean, tests/cpp from
  clean, full tinytest, equivalence trio bitwise).
- bench-sampler same-machine A/B as step 3; report the three-run table.
