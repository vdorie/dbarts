# interaction-constraints

agent: opus
rng: posterior-changing WHEN a constraint is set; the DEFAULT path (no
     interactions() prior) stays BITWISE-identical.
window: none (VD-approved build 2026-07-21; deliverable B).
budget: ~500-700 engine LOC + R surface + tests, staged. Full design:
        docs/design/interaction-constraints.md.

## Goal

Opt-in per-path interaction constraints: interactions(max.order = K,
groups = , forbid = ), per-forest for BCF (mu / tau). A constrained fit
forbids any root-to-leaf path from using more than K distinct variables
and / or co-occurring a forbidden pair; an unconstrained fit is byte-for-byte
unchanged.

## Decision (RESOLVED, VD 2026-07-21)

Build deliverable B (per-path max-order + co-occurrence), not the cheaper
per-tree grouped variant A. Rationale: the design memo "open design call".

## Context (read in code first)

- Availability primitive: collectAvailableVariables / variableAvailable /
  hasAnyAvailableVariable (tree.hpp), already ancestor-dependent via
  splitInterval. The predicate extends HERE.
- Moves + selection: drawSplitVariable / splitVariableLogProbability /
  CGMTreePrior::splitProbabilities (model.hpp); drawBirthableNode, the change
  and swap moves (moves.hpp); growTreeFromRoot (grow.hpp) is a FIFTH consumer.
- Precedents: columnMask_ via setColumnMask installs per-forest at the BCF
  sites (chain.hpp) - same here; but the subtree-validity walks
  categoricalSubtreeIsValid / findGoodOrdinalRules (moves.hpp) are
  PER-VARIABLE, and interaction needs a whole-subtree all-variables walk.

## Constraints

- DEFAULT path bitwise-identical: no interactions() prior => equivalence trio
  unchanged, NO re-record. Guard with the flag-off invariant.
- Exactness is the arbiter: gate the constrained posterior against an
  ENUMERABLE exact posterior, never a plausibility argument.
- No dbarts.h ABI change expected (interactions ride SamplerOptions + the
  bridge like DART). If a field must cross dbarts.h, STOP and flag.
- Out of scope: soft path-dependent penalties, heredity, per-tree variant A.

## Steps

1. Representation + install: a per-forest InteractionConstraint (max.order K;
   forbidden pairs as a p-bitset adjacency) in SamplerOptions, installed per
   forest mirroring setColumnMask; bridge wiring; NO behavior when unset.
   Gate: equivalence trio BITWISE (off).
2. Predicate: extend the availability primitives to read the node's ancestor
   variable-set (a bitset carried down the path) for max-order and
   co-occurrence; birth / death / grow route through it. Test: brute-force
   availability oracle.
3. Change / swap: interactionSubtreeIsValid, a whole-subtree walk carrying a
   running ancestor bitset that tests EVERY node's variable; a violation
   scores the -1.0 no-op (pi(T') = 0). DE-RISK FIRST: land the swap
   sibling-strand toy (forbid (x2, x3); root x1, children x2 / x3; swap) in
   tests/cpp BEFORE the rest. Swap DB: proposal symmetry + treeLogProbability
   carrying the constrained prior - argue it in the commit.
4. Containment: an interaction-feasibility gate in buildFromFlat /
   stateIsValid; REFUSE a warm start whose constraint differs; assert
   grow-from-root honors the predicate and monotone + interactions compose.
5. R surface: interactions() through parsePriors, bart2 + bcf(mu.interactions
   = , tau.interactions = ), Rd + NEWS; fail at fit time on unknown names,
   empty groups, K < 1.

## Verification

- tests/cpp: availability oracle; feasibility invariant after every accepted
  move; the swap-strand toy scores 0; grow-from-root + monotone coexistence.
- Exact-posterior gate: an enumerable constrained toy (small p, one forbidden
  pair / low K) - structure-visit freq + fitted means within MC error.
- Equivalence trio BITWISE with constraints OFF; full tinytest; ASAN
  tests/cpp + watch CI sanitizer green; air format --check on any R; R CMD
  check before push (codoc: interactions() in the Rd \usage).

## Landing

LANDED f455d7c (2026-07-21), two phases. Phase 1 (0141054..9bb4432, engine):
the InteractionConstraint representation + per-forest install, the ancestor-
set availability predicate (max-order + co-occurrence), and the whole-subtree
change/swap interactionSubtreeIsValid walk (the per-variable precedent could
not see the swap sibling-strand break). Phase 2 (ac43ac9..f455d7c, surface):
the interactions(max.order, groups, forbid) R prior through parsePriors, the
.Call bridge, per-forest BCF (mu.interactions / tau.interactions), the
buildFromFlat / stateIsValid / warm-start containment gate, and the Rd + NEWS.
groups is a co-occurrence ALLOW-list lowered to the forbidden-pair primitive.
NO dbarts.h change (rides SamplerOptions + the bridge). ~1.1k LOC.

Gates (orchestrator re-ran independently): equivalence BITWISE with
constraints off 27/27 identical (no re-record); the exact-posterior
enumeration matches pi*m/Z with no illegal structure; the swap sibling-strand
and containment (9-donor) tests/cpp pass; grow-from-root + monotone
coexistence pass; ASAN clean (local, instrumented) and the Phase 1 sanitizer
CI is green; R feature tests 15/15; R CMD check 0/0/0. Phase 2 sanitizer CI
watched to green post-landing.

Deferred doors (design P4): soft path-dependent penalties, formal heredity,
the per-tree block-additive variant A. The blind critique that caught the two
must-fixes is docs/design/interaction-constraints.md "structure-move
exactness" and "Containment".
