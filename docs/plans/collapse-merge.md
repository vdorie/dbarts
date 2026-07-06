# collapse-merge

agent: opus
rng: shifting (setData paths only)
budget: ~40 lines

## Goal

Both subtree-collapse sites merge leaf parameters the same way:
weighted by effective observation count.

## Context

- collapseEmptyNodesBelow: effective-observation-weighted mean
  (src/bartcore/tree.hpp:1091-1132).
- mapCutPointsBelow: plain unweighted mean over bottom nodes
  (src/bartcore/tree.hpp:1039-1053).
- The divergence is a verbatim classic port (main:src/dbarts/tree.cpp,
  collapseEmptyNodes vs mapCutPoints); no statistical rationale exists
  for the plain mean, and the rewrite's no-bit-parity license
  (docs/design/core-generalization.md) was available to unify them.

## Constraints

- The weighted form is the survivor; zero-total-weight subtrees (all
  descendants empty after remap) need the plain mean as the defined
  fallback - keep that branch explicit.
- Out of scope: any other setData remap behavior
  (mapOldCutPointsOntoNew's nearest-cut choice stays).

## Steps

1. mapCutPointsBelow accumulates (param * effectiveCount, effectiveCount)
   per bottom node; divide by the total, plain mean iff total == 0.
2. Component test in tests/cpp: a remap that starves an interval
   collapses to the weighted merge of its occupied descendants.
3. Regenerate any RNG-locked snapshots that exercise setData with
   collapsing subtrees (replay whole files); re-record equivalence
   baseline if the setData scenario shifts.

## Verification

- tests/cpp component tests pass, including the new one.
- Full tinytest suite after snapshot regeneration.
- Equivalence statistical mode passes against the prior baseline.
