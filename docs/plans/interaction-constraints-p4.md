# interaction-constraints-p4

agent: opus (engine + bridge), sonnet (docs + R harness)
rng: default path POSTERIOR-NEUTRAL (no blocks(), no warm-start violation) ->
     equivalence trio BITWISE, NO re-record. An ACTIVE blocks() prior
     legitimately shifts the stream, like any constraint.
window: none (VD "also build variant A engine" 2026-07-21, after design+critique).
budget: F1 ~60-100 LOC; variant A ~150-220 engine + R + tests (~450-600 total).
        Design: docs/design/interaction-constraints.md + this session's memos.

## Goal (three deliverables from the P4 review)

1. F1 (bug fix, independent): a columnMaskStateFeasible warm-start/setState
   gate for the EXISTING columnMask (BCF moderators + variance forest).
   Reachable via $installTrees/$setState across configs; an out-of-mask tree
   otherwise mis-scores silently (for BCF, a corrupted causal decomp).
2. Variant A (engine): blocks(groups, trees.per.group) per-tree block-additive
   f = sum_G f_G via a static per-tree columnMask, per-forest (mu.blocks /
   tau.blocks); FIXED per-group capacity over the adaptive interactions(groups=).
3. Idiom docs: interactions(groups=<disjoint total partition>) already yields
   ADAPTIVE block-additivity; document it, cross-reference blocks().

DEFER (sharpened doors): soft path-dependent penalties, formal heredity.

## Decision (RESOLVED, VD 2026-07-21)

Build all three, overriding the critique's defer-the-consumerless-engine
call; the no-fixed-capacity-consumer caveat is noted.

## Context (read in code first)

- Static mask: setColumnMask / columnMask_ / columnAllowed (tree.hpp) already
  gates variableAvailable / collectAvailableVariables; installed per-tree in
  the single-forest ctor and buildBCFForest (chain.hpp).
- Warm-start containment: interactionStateFeasible / stateIsValid (chain.hpp)
  checked interaction_ but not columnAllowed; buildFromFlatBelow (tree.hpp)
  validates gauge, not mask.
- Forest-move safety: blockMasks / blockOfTree are POD vector members (heap
  data() survives a move), NOT unique_ptr - interaction_'s UAF (04ca425) was a
  borrowed OBJECT address, a vector buffer is not.

## Refinements from the critique (MUST fold in)

- F3 (BCF-tau collision): columnMask_ is ONE pointer per tree; tau already owns
  it for moderators. blocks() ANDs the block row with any existing columnMask
  (group INTERSECT moderators / variance cols) - never overwrite. Effective
  per-tree mask = block-row & base-columnMask.
- C1/C2: blocks() REQUIRES a TOTAL, disjoint partition, validated at fit time
  (factor dummies expanded). Unlike interactions(groups=)'s per-path allow-list
  (un-named = free), a static per-tree MASK makes an un-named predictor DEAD
  (masked out of every tree), so blocks() ERRORS on un-named predictors AND on
  overlap - the coherent contract for a fixed-capacity partition.
- C3: the block-additivity invariant test special-cases a never-split tree
  (empty variable set = subset of all blocks; the intercept).

## Steps (landed; see Landing)

1. F1 gate (Opus): columnMaskStateFeasible in installForests + stateIsValid.
2. Variant A engine (Opus): blocks()/resolveBlocks, Forest blockMasks /
   blockOfTree, chain.hpp install with the F3 intersection.
3. Docs + gates (Sonnet): Rd, NEWS, test-blocks.R, docs/design P4 landing.

## Verification (orchestrator re-runs independently before push)

Default-off equivalence trio BITWISE (no re-record); tests/cpp + ASAN clean;
full tinytest incl. the expanded test-blocks.R; air format --check; R CMD
check (codoc: blocks() usage).

## Landing

blocks(groups, trees.per.group) LANDED aadbbc8 (2026-07-21); F1
(columnMaskStateFeasible warm-start/setState gate) LANDED 103dbe2, follow-up
073d3db. Full narrative, the F3 tau-moderator intersection, the C1/C2
total-partition contract, and doors 2/3 as sharpened DEFER:
docs/design/interaction-constraints.md "P4 landing".

Gates: equivalence trio bitwise (no re-record); tests/cpp incl.
testBlockAdditiveConfinement, ASAN clean; full tinytest incl. test-blocks.R
(C3, fixed capacity, reproducibility, monotone, warm-start refusal); air
format clean; R CMD check 0/0/0 incl. check_pkgdown.
