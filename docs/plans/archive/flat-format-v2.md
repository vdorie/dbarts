# flat-format-v2

agent: opus
rng: neutral (storage representation only; draws untouched)
window: pre-release for the getTrees narrow-mask value column
budget: ~900 lines; the largest single item - split PRs by step if needed

## Goal

FlatNode stores its split payload in a tagged form (ordinal cut double
OR uint64 mask word OR pool offset) instead of bit-casting masks into
the double value field. The 53-category encoding boundary, the
54-63-band side channel, and the parallel mask channels exist only
because of the double encoding and are deleted with it. getTrees keeps
its reported vocabulary (value + directions).

## Context

- FlatNode {variable, value(double), flags}: [[src/bartcore/tree.hpp:129-135@6669dfe0]];
  mask bit-cast at [[src/bartcore/tree.hpp:1177@6669dfe0]], pool-offset-in-double at [[src/bartcore/tree.hpp:1166@6669dfe0]];
  maxValueEncodableCategories = 53 ([[src/bartcore/data.hpp:119-120@6669dfe0]]).
- Machinery that exists only for double-compatibility (from the
  pooled-masks audit): the flat two-tier split at 53 vs the engine's
  real 63/64 tier, the side channel for 54-63 where the engine mask is
  inline, savedTreeMasks_/treeMasks state slots, maskWordsForFlatTree
  restore-splitting, defaulted numCategories/maskWords parameters
  through the whole flat replay family.
- State objects are opaque and engine-specific
  (docs/design/public-surface.md section 2 DECIDED) - the format may
  change freely pre-release.
- gp/linear side channels (leaf params, function blocks) ride the same
  double-stream pattern; they are NOT in scope, but the tagged FlatNode
  must not preclude widening them later.

## Decision

getTrees currently reports a narrow categorical rule's raw mask in the
value column (documented). Recommendation: report value = NA +
directions for ALL categorical rules, unifying with the pooled path -
one vocabulary, no 53 special case. This is an R-surface change and
must land pre-release or not at all. Signed off (VD, 2026-07-06):
unify.

## Constraints

- Draws must be bitwise unchanged: this touches flatten/replay/state
  only, never proposal or acceptance code.
- dbarts_sampler_getTrees keeps returning the same data.frame shape.
- Out of scope: pooled storage for the engine's live rule words (the
  >63 pool is real and stays); leaf side channels; MIA flags semantics.

## Steps

1. FlatNode gains an explicit kind and a payload union; flags absorbs
   the kind tag. flattenBelow/buildFromFlat and the replay family
   (partitionFlatIndices, countFlatObservationsBelow,
   addFlat*PredictionsBelow, flatSubtreeIsWellFormed, printFlatSubtree)
   switch on the tag; delete the defaulted numCategories/maskWords
   threading.
2. Delete maxValueEncodableCategories and the 53-boundary branches;
   masks <= 63 categories are inline uint64, wider ones reference the
   per-tree pool directly.
3. State: serialize the tagged nodes (RAWSXP or widened REALSXP pair);
   version the state object; setState validates as today.
4. getTrees/decodeCategoricalSplits: implement the reporting decision;
   update man page and the C API header comment.
5. Regenerate only tests asserting raw mask values in getTrees output;
   all draw snapshots must NOT change (that is itself the gate).

## Verification

- Component tests: flatten -> replay reproduces code-based routing
  exactly (existing tests, plus pooled/inline boundary cases at 53, 54,
  63, 64 categories).
- Full tinytest: only getTrees-reporting tests change; any RNG-locked
  draw snapshot diff is a failure.
- Equivalence compare: exact match (bitwise-neutral claim is the gate).
- stan4bart consumer test (test-capi.R) passes unchanged.

## Landing note (2026-07-06, cb25c09)

Landed: FlatNode carries a FlatKind tag (flags bits 1-2) over a
value/mask/maskOffset union plus numMaskWords; the 53-category
boundary, its side channel for 54-63, and the numCategories/types
threading are gone (the real 63/64 pooled tier remains); state
serializes the payload as raw words and stateFormatVersion bumped to
2; getTrees reports value = NA + directions for every categorical
rule (the signed-off unify). Gates: component tests incl. new
53/54/63/64 boundary cases; tinytest 0/2465 failures with only the
mandated reporting and state-plumbing test edits; equivalence EXACT
18/18 identical draws vs 235bebc (reviewer re-verified); test-capi.R
untouched and passing; lint clean. Diff 704 lines vs 900 budget.
bench-sampler not owed (flatten/replay/state only; no hot-path
change). state-continuation shares format version 2.
