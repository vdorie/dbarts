# interaction-constraints

agent: opus
rng: posterior-changing WHEN a constraint is set; the DEFAULT path (no
     interactions() prior) stays BITWISE-identical.
window: none (VD-approved build 2026-07-21; deliverable B).
budget: ~500-700 engine LOC + R surface + tests, staged into commit-sized
        steps. Full design: docs/design/interaction-constraints.md.

## Goal

Opt-in per-path interaction constraints: interactions(max.order = K,
groups = , forbid = ), per-forest for BCF (mu / tau). A constrained fit
forbids any root-to-leaf path from using more than K distinct variables
and / or from co-occurring a forbidden pair; an unconstrained fit is
byte-for-byte unchanged.

## Decision (RESOLVED, VD 2026-07-21)

Build deliverable B (per-path max-order + co-occurrence), not the cheaper
per-tree grouped variant A. Rationale and the A/B split: the design memo
"The open design call for VD".

## Context (read in code first)

- Availability primitive: collectAvailableVariables / variableAvailable /
  hasAnyAvailableVariable (tree.hpp), already ancestor-dependent via
  splitInterval (cut exhaustion). The predicate extends HERE.
- Selection + moves: CGMTreePrior::splitProbabilities /
  splitVariableLogProbability / drawSplitVariable (model.hpp);
  drawBirthableNode, the change and swap moves (moves.hpp); growTreeFromRoot
  (grow.hpp) is a FIFTH split consumer.
- Install precedent: columnMask_ / columnAllowed (tree.hpp) via setColumnMask,
  per-forest at the BCF sites (chain.hpp). The constraint installs the same.
- Subtree-validity precedent: categoricalSubtreeIsValid / findGoodOrdinalRules
  (moves.hpp) - but PER-VARIABLE; interaction needs a whole-subtree,
  all-variables walk (Step 3).

## Constraints

- DEFAULT path bitwise-identical: no interactions() prior => equivalence trio
  unchanged, NO re-record. Guard with the flag-off invariant.
- Exactness is the arbiter: the constrained posterior gates against an
  ENUMERABLE exact posterior, never a plausibility argument.
- No dbarts.h ABI change expected (interactions ride SamplerOptions + the
  bridge like DART / splitProbabilities). If a field must cross dbarts.h,
  STOP and flag (ABI checklist + stan4bart lockstep + MINOR bump).
- Out of scope (deferred in the memo): soft path-dependent penalties, formal
  heredity, the per-tree grouped variant A.

## Steps

1. Representation + install: a per-forest InteractionConstraint (max.order K;
   forbidden pairs as a p-bitset adjacency) in SamplerOptions, installed per
   forest mirroring setColumnMask; bridge wiring; NO behavior when unset.
   Gate: equivalence trio BITWISE (constraint off).
2. Predicate: extend collectAvailableVariables / variableAvailable /
   hasAnyAvailableVariable to read the node's ancestor variable-set (a bitset
   carried down the path) for max-order and co-occurrence. Birth / death /
   grow route through it (correct-by-construction). Component test:
   brute-force availability oracle.
3. Change / swap: interactionSubtreeIsValid - a whole-subtree walk carrying a
   running ancestor bitset, testing EVERY node's variable; a violation scores
   the -1.0 no-op (pi(T') = 0 auto-reject). DE-RISK FIRST: land the swap
   sibling-strand toy (forbid (x2, x3); root x1, children x2 / x3; swap root
   and left) in tests/cpp BEFORE the rest of the move work. Swap DB rests on
   proposal symmetry + treeLogProbability(parent) carrying the constrained
   prior - argue it in the commit.
4. Containment: an interaction-feasibility gate in buildFromFlat /
   stateIsValid; REFUSE a warm start whose constraint differs from the
   target; assert grow-from-root honors the predicate; assert monotone +
   interactions compose (orthogonal seams).
5. R surface: an interactions() prior object through parsePriors, bart2 +
   bcf(mu.interactions = , tau.interactions = ), Rd + NEWS; safe-over-fast
   validation (unknown names, empty groups, K < 1 error at fit time).

## Verification

- tests/cpp: brute-force availability oracle; feasibility invariant after
  every accepted move; the swap-strand toy scores 0; grow-from-root and
  monotone-coexistence.
- Exact-posterior gate: an enumerable constrained toy (small p, one forbidden
  pair / low K) - structure-visit frequencies + fitted means within Monte
  Carlo error.
- Equivalence trio BITWISE with constraints OFF (default unchanged, no
  re-record); full tinytest; ASAN tests/cpp for the new reachable engine code
  + watch the CI sanitizer green; air format --check on any R.
- R CMD check before push (codoc: interactions() in the Rd \usage).
