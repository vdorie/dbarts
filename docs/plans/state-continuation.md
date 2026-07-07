# state-continuation

agent: opus
rng: shifting on the restore path only (fresh runs unchanged)
budget: ~300 lines, mostly deletions

## Goal

Saved state carries what the model needs (trees, parameters, sigma, k,
latents, DART state, RNG streams) and restore rebuilds the rest
canonically. The four fields that exist only to reproduce the last ulp
- per-tree observation index buffers, accumulated totalFits,
internal-scale sigma, sigmaPriorScale - are gone, along with the
bitwise-continuation contract that required them.

## Decision

Question: keep bitwise-exact restore? Recommendation: drop to semantic
restore. Classic never promised it (rebuildScratchFromState recomputed
totalFits and partitions on load); the contract is self-imposed
(src/bartcore/chain.hpp:175-206 documents the fields and why), costs
O(n * numTrees) state size, and stan4bart's splice-and-predict workflow
needs the same trees and sigma, not the same accumulation history.
Evidence that would change it: a consumer whose correctness (not
convenience) depends on run(saved) == run(never-saved) bitwise.
Signed off (VD, 2026-07-06): drop to semantic restore.

## Constraints

- Depends on constant-leaf-suffstat: with the order-insensitive
  suffstat landed, dropping the index buffers costs at most last-ulp
  differences; without it, restored chains diverge faster - land that
  first.
- setState validation (scratch-tree check, reject-with-nothing-changed)
  is behavior, not bitwise bookkeeping - keep.
- The state object stays opaque; version it (coordinate with
  flat-format-v2 if both are in flight - one version bump).
- Out of scope: keepTrees storage; getTrees.

## Steps

1. Remove indices/totalFits/sigmaPriorScale from ChainStateData and
   the get/setState paths; store sigma original-scale; restore
   recomputes totalFits from tree fits and partitions from tree
   structure + cut points (the classic recipe).
2. Component tests: the bitwise round-trip gates become (a) structural
   equality of trees/parameters after round trip, (b) statistical
   agreement of a continued run vs an uninterrupted one (calibrated
   bound), (c) setState rejection still mutates nothing.
3. Update chain.hpp's ChainStateData doc comment and the state
   paragraph in docs/design/core-generalization.md's landing notes
   (one-line amendment, not history rewriting).
4. test-capi.R state round-trip: predictions-identical assertion stays
   (predictions replay trees; unaffected).

## Verification

- Component tests as restructured; full tinytest.
- Equivalence exact on fresh runs (nothing outside restore changed).
- stan4bart consumer test passes.
