# interaction-constraints-p4

agent: opus (engine + bridge), sonnet (docs + R harness)
rng: default path POSTERIOR-NEUTRAL (no blocks(), no warm-start violation) ->
     equivalence trio BITWISE, NO re-record. An ACTIVE blocks() prior
     legitimately shifts the stream, like any constraint.
window: none (VD "also build variant A engine" 2026-07-21, after design+critique).
budget: F1 ~60-100 LOC; variant A ~150-220 engine + R + tests (~450-600 total).
        Design: docs/design/interaction-constraints.md + this session's memos.

## Goal (three deliverables from the P4 review)

1. F1 (bug fix, independent): a columnMaskStateFeasible warm-start / setState
   gate for the EXISTING columnMask (BCF moderators + restricted variance
   forest). Reachable via $installTrees across BCF fits with different
   moderators, and cross-sampler $setState: installs a tree splitting on a
   masked column -> splitVariableLogProbability mis-scores over a menu that
   excludes it -> silent wrong posterior (for BCF, a corrupted causal decomp).
2. Variant A (engine): blocks(groups, trees.per.group) per-tree block-additive
   f = sum_G f_G via the static per-tree columnMask; per-forest (mu.blocks /
   tau.blocks). Adds FIXED per-group tree capacity + a "each tree in one block
   forever" guarantee over the already-adaptive interactions(groups=<partition>).
3. Idiom docs: interactions(groups=<disjoint total partition>) already yields
   ADAPTIVE block-additivity; document it, cross-reference blocks().

DEFER (sharpened doors): soft path-dependent penalties, formal heredity.

## Decision (RESOLVED, VD 2026-07-21)

Build all three. VD chose to build the engine over the critique's
defer-the-consumerless-engine call; the no-fixed-capacity-consumer caveat is
noted, overridden.

## Context (read in code first)

- Static mask: setColumnMask / columnMask_ / columnAllowed (tree.hpp);
  columnAllowed already gates variableAvailable AND collectAvailableVariables.
  Installed per-tree in the single-forest ctor and buildBCFForest (chain.hpp).
- Warm-start containment: interactionStateFeasible / stateIsValid /
  rebuildLiveForest (chain.hpp), installForests (sampler.hpp) check interaction_
  but NOT columnAllowed; buildFromFlatBelow (tree.hpp) validates gauge, not mask.
- Forest-move safety: blockMasks / blockOfTree are POD vector members (heap
  data() survives a move) - NOT unique_ptr; the interaction_ UAF (04ca425) was a
  borrowed OBJECT address, a vector buffer is not.
- R precedents: resolveInteractions, resolveColumns (dummy expansion via
  startsWith), the mu/tau split (R/bartcore.R), parsePriors, model attributes.

## Refinements from the critique (MUST fold in)

- F3 (BCF-tau collision): columnMask_ is ONE pointer per tree; tau already owns
  it for moderators. blocks() ANDs the block row with any existing columnMask
  (group INTERSECT moderators / variance cols) - never overwrite. Effective
  per-tree mask = block-row & base-columnMask.
- C1/C2: blocks() REQUIRES a disjoint partition (validated at fit time, factor
  dummies expanded); un-named predictors default to singleton blocks (additive
  main effect) so no variable goes dead.
- C3: the block-additivity invariant test special-cases a never-split tree
  (empty variable set = subset of all blocks; the intercept).

## Steps (SERIALIZED; orchestrator gates + commits between)

1. F1 gate (Opus): columnMaskStateFeasible (O(nodes) columnAllowed scan,
   mirror interactionStateFeasible) into installForests + stateIsValid, over the
   shipped columnMask. tests/cpp regression (mask-violating install -> refused)
   + R warm-start-refusal test (mirror interactionMismatch).
2. Variant A engine + bridge + R (Opus): blocks() + resolveBlocks (disjoint
   validation, dummy expansion, even default capacity, deterministic contiguous
   NO-RNG assignment); bridge into SamplerOptions / BCFForestSpec (NO dbarts.h
   change); Forest blockMasks (G*p) + blockOfTree; chain.hpp install at both
   sites with the F3 intersection; the F1 gate now covers block masks.
3. Docs + gates (Sonnet): blocks() Rd + NAMESPACE + _pkgdown.yml;
   interactions.Rd idiom cross-ref; NEWS; block-additivity invariant test
   (C3-aware, reuse test-interactions.R helpers) + reproducibility snapshot +
   monotone coexistence; update docs/design P4 + this Landing.

## Verification (orchestrator re-runs independently before push)

Default-off equivalence trio BITWISE (27/27 + bcf + multinomial, NO re-record);
tests/cpp cross-ISA from make clean; full tinytest; block-additivity invariant;
warm-start refusal (F1 + blocks); monotone coexistence; ASAN tests/cpp + watch
the CI sanitizer to green (F1 is R-reachable BCF, the 04ca425 class); air format
--check; full R CMD check (codoc: blocks() \usage).

## Landing

(to fill)
