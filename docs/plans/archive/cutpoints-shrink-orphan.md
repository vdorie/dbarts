# cutpoints-shrink-orphan

agent: opus
rng: to be classified by the investigation (the fix may touch the
     setCutPoints remap, a live path)
budget: ~100 lines

## Goal

setCutPoints shrinking a column's grid never leaves a split indexing
past the new grid: flatten (state, getTrees) reads a real cut value
for every surviving rule.

## Context

- Found probing a strict state round-trip in the mutation fuzzer
  (fuzz-state-roundtrip.md landing note): after a grid-shrinking
  setCutPoints, a split can survive with an out-of-range splitIndex;
  flatten then reads past cutPoints[j], serializing a stale value.
  Potentially user-facing through getTrees, not just setState.
- mapCutPointsBelow/collapse are supposed to remap or collapse every
  affected split (tree.hpp); the fuzzer path that escapes them needs
  characterizing first - suspect saved trees flattened under the old
  grid, or the change/swap validity walk missing an interior rule.

## Constraints

- setState's reject-with-nothing-changed contract stays; do not relax
  buildFromFlat validation.
- Out of scope: mapOldCutPointsOntoNew's nearest-cut choice.

## Steps

1. Reproduce via the fuzzer (tighten OP_STATE locally) and reduce to
   a direct component test; characterize the escaping path.
2. Fix at the source (remap/collapse coverage, or refuse the shrink
   for affected saved trees); classify the rng impact honestly.
3. Tighten the fuzzer's OP_STATE round-trip assertion to strict now
   that all three edges are fixed.

## Verification

- Direct test + fuzzer >= 100 seeds with strict OP_STATE; full
  tinytest; equivalence vs current baseline (exact if the fix stays
  off live paths, else statistical + re-record).

## Landing note (2026-07-07, see merge commit)

Landed: collapseEmptyNodes also collapses ordinal splits whose index a
shrunken grid no longer addresses (ruleIsUnrepresentable); the
escaping path was missing values riding missingGoesRight keeping both
children occupied, so the empty-child test alone spared the orphan.
The fuzzer's OP_STATE round-trip assertion is now STRICT permanently.
Gates: direct component test, fuzzer 100 seeds (reviewer re-ran),
full suite, tinytest 2470 ok, equivalence exact 18/18 - the new
collapse is a no-op everywhere numCuts never shrinks. Source diff 51
lines. All three fuzzer-found edges are now closed.
