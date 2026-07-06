# mutation-journal

agent: opus
rng: neutral
budget: ~250 lines

## Goal

Mutation rollback stops copying state it can rebuild or replace:
setPredictor builds the new codes into fresh storage and swaps on
success (no full-array copy-then-restore), the patch paths journal only
changed cells, and each entry point's cost contract is documented where
the design doc's O(changed cells) promise now holds or is explicitly
amended.

## Context

- Full-matrix setPredictor snapshots the entire n x p codes array for
  rollback (src/bartcore/sampler.hpp:614 `oldCodes = data_.codes;`) -
  the copy equals the operation's own cost, but build-new-and-swap
  gets rollback for free and halves peak traffic.
- Column form already snapshots only touched columns
  (sampler.hpp:656-701); the per-observation session
  (sampler.hpp:817-895) is the only path meeting the journal
  aspiration (docs/design/core-generalization.md:162-165).
- Tree re-routing and fit rebuilds dominate for small patches; the
  journal only helps the data-layer side. Do not oversell.

## Constraints

- Transaction semantics unchanged: validate-or-rollback, forceUpdate
  collapse, tri-state returns, all-or-none across chains.
- Bitwise-identical draws on the accept path (same codes, same order).
- Out of scope: data-ownership (separate item; rebase whichever lands
  second), per-column widths (hot-layer-u8).

## Steps

1. setPredictor: quantize into scratch, validate every chain against
   scratch, swap on accept, drop on reject; delete the oldCodes copy.
   Same for the cuts-refresh variant.
2. setCell/updatePredictor: journal (index, old value) pairs instead of
   column copies where the touched-cell count is small (threshold:
   copy when cells > column/4, measured not guessed).
3. Amend the design doc's mutation-contract paragraph with the actual
   per-entry-point costs.

## Verification

- Component tests: accept and reject paths bitwise-identical to
  current behavior on fixed inputs (existing mutation tests are the
  oracle; run before and after).
- Full tinytest; equivalence exact.
- bench-sampler setPredictor accept/reject metrics: no regression;
  expect improvement on the accept path.
