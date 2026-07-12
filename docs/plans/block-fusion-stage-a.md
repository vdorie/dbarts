# block-fusion-stage-a (b=1 exact refactor)

COMPLETE 2026-07-11: all seven commits landed, default flipped ON. Milestone gate held with NO re-record - equivalence 22/22 IDENTICAL, tinytest 2728/0 no snapshot regen, tests/cpp + full-vocabulary fuzzer clean, and bench-sampler NEUTRAL on the x86 box (interleaved A/B, every run-* metric within 4%, none over the 5% gate) after driving the b=1 overhead down (A cache + atomOf + S carry gated off at b=1, the three per-tree passes fused). See the Stage A landing note in docs/design/block-fusion.md section 9.

agent: opus
rng: NEUTRAL - bit-identical by construction. Stage A computes the SAME
  node sufficient statistics (sumWeights, sumWeightedResponse,
  sumWeightedResponseSq) by an atom-vocabulary path, feeds them to the
  UNCHANGED model math, in the SAME leaf/RNG order and the SAME
  summation order. No draw byte changes. Gate: equivalence 22/22
  IDENTICAL (exact-match mode) at EVERY commit; tinytest full pass with
  NO snapshot regeneration; tests/cpp bit-exact.
window: Stage A of docs/design/block-fusion.md (frontier 3.4, DESIGN
  PROPOSAL 2026-07-11, VD approved building + the one-time re-record).
  FOLLOWS data-ownership settling (design 5.2): branch from origin/bartcore (6086ccd). NOTE: the u8/width-templated
  partition (plan-1 steps 3-4) is PARKED, NOT landed; Stage A uses the
  current partitionChildren, which the atom split wraps (no u8 dependency). Stage B (b>1, the WIN, the re-record) builds on this map;
  Stage A must land bit-identical FIRST as the de-risk anchor.
budget: ~1400-2000 lines across ~7 files, seven staged commits. New
  header src/bartcore/atoms.hpp (~400-600) + new tests/cpp/test_atoms.cpp
  (~500-700); the rest are surgical hooks in chain.hpp / tree.hpp /
  moves.hpp behind a flag. Kill-switch cheap at every commit.
bench: VD grants the quiet machine. Fresh speed baseline is already
  current (bench-sampler-32fc7c8.csv). Stage A expects NEUTRAL (b=1 does
  the same O(n) work); a real regression = the atom-bookkeeping overhead
  is too high = kill signal for the WHOLE approach, caught cheaply
  before any FP change (this is the point of the stage). Confirmatory
  compare at commit (i), full compare at the milestone (vii).

---

## 0. What Stage A is, in one screen

Wire the atom-map apparatus - build, per-leaf aggregation, split (birth),
merge (death), snapshot/restore (change/swap), and cross-sweep persistence -
so the constant-leaf Gaussian sweep runs through it AT b=1, where an atom IS a
leaf and g IS treeY, and reproduces today's draws BITWISE. It proves the whole
machinery Stage B needs (design section 4.4, the primary correctness anchor)
with ZERO re-record before any floating-point regrouping happens in Stage B.

The insight that makes it clean and safe: the engine already funnels ALL
constant-leaf statistics through ONE seam - the node suffstat cache
`Node::{sumWeights, sumWeightedResponse, sumWeightedResponseSq}`
(tree.hpp:175-177). Every consumer reads ONLY that cache:
- move scoring: `logLikelihoodForBranch` -> `logIntegratedLikelihoodForNode`
  reads the node cache (moves.hpp:47-65, model.hpp:146-154);
- leaf draw: `drawFromPosteriorForNode` reads the node cache (model.hpp:156-162);
- k accumulation + fit scatter: `sampleParametersAndSetFits` reads the drawn
  param (chain.hpp:2105-2136).
Only FOUR sites WRITE that cache: `setNodeAverages`/`computeLeafStats`
(tree.hpp:492-521), `birth`'s two child `computeLeafStats` (tree.hpp:780-781),
`refreshSubtree`'s leaf `computeLeafStats` (tree.hpp:760), and
`orphanChildren`'s additive merge (tree.hpp:796-800).

Stage A introduces an atom path that PRODUCES those same cache values from an
atom SoA, behind a flag, and proves - per tree, bitwise - that it equals the
current writer, then flips the flag and lets the equivalence anchor confirm
end-to-end. Because the seam and every downstream reader are untouched, the
only thing that can differ is HOW the identical numbers were computed.

## 1. The b=1 wiring + bitwise-proof strategy (the crux)

### 1.1 Recommendation: option (a) - a compile-time flag with a tight per-tree oracle

Implement the atom path behind a compile-time flag (a `constexpr bool` /
build macro; see 1.5), default OFF, flipped ON at commit (vii). NOT option (b)
(replace outright, rely only on the end-to-end equivalence hash).

Justification:
1. LOCALIZATION. A single-bit divergence in the end-to-end equivalence hash
   tells you the sweep diverged but not WHERE. With the flag, tests/cpp runs
   BOTH writers on the same tree + same residual and asserts the node caches
   are bitwise-equal per leaf (`std::memcmp`/exact `==` on the three doubles),
   over a fuzzed move sequence. A divergence surfaces at the exact leaf, move,
   and op - debuggable, not a needle in a 2000-sweep haystack.
2. DETERMINISTIC ORACLE WITHOUT REPLAYING RNG. The differential is on the
   suffstat, which is a pure function of (residual, membership) - no RNG. Since
   `drawFromPosterior` is the UNCHANGED function of the suffstat, bitwise-equal
   suffstats => bitwise-equal draws given the same RNG stream. So the oracle
   never has to reset/replay the RNG to compare draws; it compares the
   deterministic caches, which is stronger and simpler.
3. CLEAN ABORT. If the milestone bench regresses (the kill signal) or a bug
   resists, flip the flag OFF: the engine is byte-for-byte today's, zero churn,
   ship nothing. Option (b) has no abort short of a revert.
4. NOT THROWAWAY. The flag is the natural seam Stage B needs anyway (b=1 vs
   b>1 vs off). Stage A's flag becomes Stage B's `b` dispatch. The dual path
   also keeps the fuzz differential (atom-vs-live) runnable as a permanent
   regression guard.

Cost: two coexisting paths for one stage, and a modest #ifdef/`if constexpr`
footprint. Acceptable and temporary.

### 1.2 The seam: the atom path WRITES the node suffstat cache

At b=1 the atom path does NOT replace the draw, the scoring, the scatter, the
roll, the sigma/latent stages, or `totalFits`. It replaces only the FOUR cache
writers with atom-vocabulary equivalents that produce bitwise-identical values:

| current writer (flag OFF)                    | atom writer (flag ON), b=1                                            |
|----------------------------------------------|-----------------------------------------------------------------------|
| `setNodeAverages` -> `computeLeafStats` all leaves | `AtomMap::aggregateTree`: for each leaf, aggregate its one atom's (A,G,Q) and store into the node cache |
| `birth` -> `computeLeafStats(left/right)`    | `AtomMap::splitAtom`: partition the parent atom, aggregate each child's (A,G,Q) into the two child node caches |
| `orphanChildren` additive merge              | `AtomMap::mergeAtoms`: A/G/Q of children add into the parent; the node-cache merge is ALREADY additive over children, so it stays as-is (see 1.4) |
| `refreshSubtree` -> `computeLeafStats(leaf)` | `AtomMap::refreshSubtree`: re-slice + re-aggregate the subtree's atoms into caches |

Because at b=1 an atom's members are exactly a leaf's members and its (A,G,Q)
are exactly (sumWeights, sumWeightedResponse, sumWeightedResponseSq) over the
residual g=treeY, "aggregate the atom" IS "call the same misc kernel over the
same member slice." Bitwise equality is by construction (see 1.3), not by luck.

### 1.3 The b=1 aggregation is the SAME kernel call - do NOT reimplement the reduction

The single most important bitwise rule: the atom aggregation must invoke the
EXISTING misc suffstat kernels over the EXISTING member order, NOT a hand-rolled
loop. `misc_computeIndexedSufficientStatisticsFast` (moments.c:334) and its
weighted/non-indexed siblings use a specific mod-5 remainder + stride-5
horizontal-sum accumulation tree (moments.c:315-332, 387-416). Any different
loop (even a mathematically identical sum) rounds differently and breaks the
anchor. So `AtomMap::aggregateAtom(leaf)` calls exactly the kernel
`computeLeafStats` would call, over the same slice.

Concretely, at b=1, for an atom c mapped to leaf L with member slice
`members[c.begin..c.end)`:
- A(c) = sumWeights, G(c) = sumWeightedResponse, Q(c) = sumWeightedResponseSq,
  obtained by the SAME `computeLeafStats` kernel dispatch (weighted vs
  unweighted, root vs non-root - see the root pin below) over the SAME member
  order. The three land in both the atom SoA AND the node cache.
- S(c) is initialized/updated by the draw (1.4); it is not part of the b=1
  suffstat and does not enter the cache.

Note Q(c) (= sum w g^2) is carried at b=1 because `logIntegratedLikelihood`
(model.hpp:117-118) needs `sumWeightedResponseSq`. Stage B DROPS Q (fact 1.2:
it cancels in every move ratio), but Stage A keeps it so the anchor is
bitwise. Keeping Q at b=1 is what makes the aggregation literally the
unchanged kernel.

### 1.4 The draw, roll, and merge stay byte-identical

- DRAW. The draw loop stays `sampleParametersAndSetFits` (chain.hpp:2105-2136)
  UNCHANGED: it iterates leaves via `fillBottom(0, bottoms)`, draws one normal
  per NON-EMPTY leaf from the node cache, forces 0.0 (no draw) for empty
  leaves, accumulates k over non-empty leaves, and scatters via
  `misc_setIndexedVectorToConstant`. Since the atom path already wrote the node
  cache, the draw reads identical inputs and consumes RNG identically. The
  ONLY atom-side action at draw time is bookkeeping: set S(c)=mu for c's atom
  (so Stage B's S carry is exercised and testable), which does NOT touch RNG or
  the cache.
- ROLL. At b=1 the residual roll is the existing per-tree fused O(n) pass
  (chain.hpp:731-744) and the block-exit `totalFits` rebuild (chain.hpp:779-787)
  - both UNCHANGED. g=treeY is used directly; NO separate g field is
  materialized at b=1 (design 4.4: "the block-entry g build is the current
  per-tree roll, unchanged"). Introducing a g scatter at b=1 would be a
  gratuitous FP-order risk; do not.
- MERGE (death). `orphanChildren` (tree.hpp:796-800) already forms the parent
  cache as left+right of the child caches, in that order. If the child caches
  hold atom-correct values, the parent cache is atom-correct with NO change to
  `orphanChildren`. The atom path's `mergeAtoms` only updates the atom SoA
  topology (child atoms -> parent atom); it must add A/G/Q in the SAME
  (left-then-right) order for the SoA to match, but the node cache is produced
  by the unchanged merge. Keep them consistent; the cache is the oracle.

### 1.5 Every place order matters, and how it is pinned (the RNG-neutrality traps)

Enumerate and pin all of these - each is a way to silently break the anchor:

1. LEAF / RNG DRAW ORDER. Both `setNodeAverages` and the draw iterate leaves
   via `fillBottom(0, out)` - a DFS, left child (`leftChild`) before right
   (`leftChild+1`), arena order (tree.hpp:282-286). The draw consumes one
   normal per non-empty leaf in this order. PIN: the atom path drives
   aggregation and draw off `fillBottom`, NEVER off an atom-id iteration order.
   At b=1 atom-id can equal leaf node id, but iterate the LEAF list.
2. WITHIN-LEAF SUMMATION ORDER. The kernel's mod-5 + stride-5 tree
   (moments.c). PIN: call the existing kernel over the member slice; never
   re-sum by hand. Verified: the four kernels (indexed/non-indexed x
   weighted/unweighted) are the only accumulators; `aggregateAtom` dispatches
   to the identical one `computeLeafStats` picks.
3. ROOT SPECIAL CASE. `computeLeafStats` uses the NON-indexed kernel
   (`misc_compute[Weighted]SufficientStatisticsFast`) when `node.parent ==
   invalidNode` (tree.hpp:494-508), i.e. a stump/root leaf, reading y[0..n)
   directly. For identity indices the non-indexed and indexed kernels are
   bitwise-equal (verified: identical mod-5+stride-5 structure, and a root leaf
   never had its indices permuted, so they are identity), but PIN the branch
   anyway: `aggregateAtom` replicates the `isRoot` dispatch. Zero-risk and
   future-proof against a kernel divergence.
4. BIRTH CHILD ORDER. `birth` calls `partitionChildren` then
   `computeLeafStats(left)` then `computeLeafStats(right)` (tree.hpp:779-781),
   children over their post-partition slices (always non-root -> indexed
   kernel). PIN: `splitAtom` aggregates left child then right child via the
   indexed kernel over the same slices.
5. DEATH MERGE ORDER. `orphanChildren`: parent = left + right (tree.hpp:796-800).
   PIN: SoA merge adds A/G/Q left-then-right; cache merge unchanged.
6. CHANGE / SWAP REFRESH ORDER. `refreshSubtree` does `partitionChildren` then
   recurses left, right, computing leaf stats in DFS order (tree.hpp:757-766).
   PIN: `AtomMap::refreshSubtree` mirrors the identical recursion order.
7. PARTITION PERMUTATION. `partitionChildren` (tree.hpp:672-739) is the
   width-templated two-pointer (u16 misc kernels; u8 scalar
   `partitionIndicesScalar`/`partitionRangeScalar`; categorical mask; sparse;
   MIA). Its permutation is integer-exact and ISA-stable. PIN: the atom split's
   member partition REUSES this exact primitive (design 5.2); it does NOT fork
   a second partitioner. (At b=1, if members alias tree.indices, the tree's own
   `partitionChildren` IS the atom split; see 2.4.)
8. EMPTY-LEAF HANDLING. Empty leaf => param forced 0.0, NO draw consumed, NOT
   counted toward k (chain.hpp:2109-2120). An empty leaf has zero atoms. PIN:
   the draw still visits it via `fillBottom` and forces 0.0 without drawing;
   the atom path must not skip the leaf list nor draw for empties.
9. k ACCUMULATION. `kSumSquaredParams += param*param` and `kNumLeaves += 1`
   over non-empty leaves in `fillBottom` order (chain.hpp:2117-2119). PIN:
   draw loop unchanged, so this is automatic; do not move it into the atom path.
10. SIGMA / SSR / LATENTS. `misc_computeSumOfSquaredResiduals` over totalFits
    (chain.hpp:908-910), `refreshLatents`, `drawSigma` - all once/sweep,
    OUTSIDE the block, over the totalFits rebuilt at block exit. PIN: untouched;
    the block-exit totalFits rebuild stays chain.hpp:779-787.
11. BCF / GROUPED-RE FOREST RESPONSE. `formForestResponse` builds each forest's
    residual (chain.hpp:700-704). PIN: at b=1 the block is per tree within the
    already-formed forest residual; nothing changes. (Blocks are per-forest in
    Stage B; irrelevant at b=1.)
12. GROW-FROM-ROOT WARM START (chain.hpp:938-1048) stays on the per-tree path
    (design 5.3): a fixed init few-sweeps, not hot. PIN: the flag engages only
    in the steady-state `run()` loop, never in the warm-start loop.

## 2. Where the atom map lives

### 2.1 New header: src/bartcore/atoms.hpp

A header-only `struct AtomMap` (or `BlockAtoms`), included by chain.hpp. It is
NOT templated on the leaf model (constant-leaf only, design 5.3); it is a plain
struct operating on doubles + a `Tree&` + the `ColumnStore&`. Keeping it in its
own header mirrors the layering and keeps chain.hpp's diff small.

### 2.2 Owned as a Forest member

One `AtomMap` per `Forest` (chain.hpp:257-307), added next to `treeFits`,
`treeY`, `scratch`. The Forest is the natural owner: it already holds the
per-tree index buffers (`indexBuffer`, sliced n per tree, chain.hpp:452-455),
`treeFits`, and `MoveScratch`. Sized once at forest init (the `initialize`
site, chain.hpp:452-459) and reused across sweeps (persistence, design 4.5).
For BCF/grouped, each Forest owns its own map (blocks are per-forest).

### 2.3 The struct (constant-leaf, b=1-capable, Stage-B-shaped)

```
struct AtomMap {
  // --- topology (persists across sweeps; patched by moves) ---
  // At b=1: atomId <-> leaf arena index, 1:1 for non-empty leaves.
  std::vector<uint32_t> atomOf;      // size n: obs -> atom id (Stage B move bookkeeping)
  std::vector<size_t>   atomBegin;   // per-atom: slice start into the member buffer
  std::vector<size_t>   atomEnd;     // per-atom: slice end
  std::vector<int32_t>  leafTuple;   // per-atom x b: leaf id per block tree (b=1 => 1 int per atom)
  // --- static-over-block (residual-independent; A cached, U'WU template) ---
  std::vector<double>   A;           // per-atom weight mass (count when unweighted)
  // --- dynamic-over-block (residual-dependent; rebuilt at block entry) ---
  std::vector<double>   G;           // per-atom residual mass = sum w*g
  std::vector<double>   Q;           // per-atom sum w*g^2 (b=1 anchor only; dropped Stage B)
  std::vector<double>   S;           // per-atom in-block fit (b=1: mu of the atom's leaf)
  // --- member buffer ---
  // b=1 recommendation (2.4): a VIEW onto the block-tree's tree.indices
  // (size_t* members + n), NOT a second owned n-buffer. Stage B introduces an
  // owned per-block buffer here.
  size_t* members = nullptr;         // aliases tree.indices at b=1
  // --- rollback workspace (mirrors Tree::SubtreeSnapshot) ---
  struct AtomSnapshot {
    std::vector<uint32_t> atomIds;
    std::vector<double>   A, G, Q, S;
    std::vector<size_t>   begin, end;
    std::vector<int32_t>  leafTuple;
    // member-segment snapshot only if members is owned (Stage B); at b=1 the
    // tree's own SubtreeSnapshot already restores tree.indices == members.
  } snapshot;
  // per-atom free list / recycling for split/merge id reuse.
};
```

Sizes: atomOf (n uint32) + A/G/Q/S (K doubles each, K = #leaves) + small
per-atom vectors. At b=1 with `members` aliasing tree.indices, the ONLY n-sized
add is `atomOf` (4n bytes) plus the K-sized SoA (tiny). No n*m addition. This is
the memory-neutral choice (2.4).

### 2.4 Member buffer: ALIAS tree.indices at b=1 (recommended), own it at Stage B

Two designs; both bit-identical. The choice is implementer-decidable (not a VD
call), decided by the bench/memory gate:

- DESIGN A (RECOMMENDED for Stage A): `members` aliases the block-tree's
  `tree.indices`. The atom's member partition at b=1 IS the tree's own
  `partitionChildren` (same buffer, same primitive), so the "atom split" wraps
  the tree's partition and then aggregates the child SoA; rollback of the member
  order is the tree's own `restoreSubtree`. The NEW, genuinely-tested code is
  the per-atom SoA lifecycle (aggregate/split/merge/refresh/snapshot of
  A/G/Q/S + topology). Pros: zero added n*m memory (critical at the n=1e6,
  m=200 bench point), zero added member-copy pass, so bench stays neutral. Cons:
  the owned-buffer two-pointer split is not exercised until Stage B.
- DESIGN B (alternative): `members` is an owned per-block buffer (n per block =>
  n*m at b=1, DOUBLING the index-buffer memory), maintained in lockstep with
  tree.indices by an explicit atom two-pointer split. Pros: exercises the
  owned-buffer split one stage earlier. Cons: n*m memory bump risks OOM/thrash
  at the largest bench point and an added O(n) build per block; the bench gate
  (the kill signal) is likely to flag it - defeating the "b=1 is neutral" anchor.

Recommend DESIGN A. Stage A's mandate (prove build / aggregation / split /
rollback / persistence) is fully met by DESIGN A because the SoA split, merge,
snapshot, and persistence are all real new code proven bitwise; only the member
BUFFER is deferred, and it is the one piece that pays for itself only at b>1.
Record this as a resolved implementer decision in the landing note.

### 2.5 Persistence via the U'WU-cache template (design 4.5)

Model the residual-INDEPENDENT persistence exactly on `LinearGaussianLeaf`'s
crossproduct cache (model.hpp:458-539): cache A(c) + the atom topology, keyed by
the atom's ordered member list, and re-validate on lookup by `std::memcmp`
against `tree.indices[begin..end)` (model.hpp:509-511). A structural move or a
rejected-move rollback that alters membership fails the compare and rebuilds A;
an unchanged leaf serves the cached A bitwise. This is "rollback-stable with no
per-move hook" (model.hpp:459-462) - the hardest part (4.5) inherits a landed,
proven pattern. G/Q/S are residual-DEPENDENT and always recomputed at block
entry from the current treeY (the block-entry aggregation), exactly as the U'WU
cache always rescans U'Wz / z'Wz. At b=1, A recompute is cheap, but exercising
the cache proves the Stage B mechanism.

## 3. Commit-by-commit breakdown

Each commit compiles, gates, and has a clean abort (flag OFF). The flag is
`AtomMap`-gated in chain.hpp; until commit (vii) the default build path is
byte-for-byte today's. Budgets are rough.

### (i) Atom-map struct + b=1 build + build oracle. ~250 lines. 1-2 days.

- FILES: new src/bartcore/atoms.hpp (struct + `buildForTree`); chain.hpp
  (Forest member + init sizing, no wiring yet); new tests/cpp/test_atoms.cpp;
  tests/cpp/Makefile + main.cpp (register the suite).
- WHAT: `AtomMap::buildForTree(tree, data)` builds the b=1 map from the tree's
  leaf partition: one atom per non-empty leaf, `atomBegin/atomEnd` = the leaf's
  `[begin,end)`, `leafTuple` = leaf id, `atomOf[i]` = the atom for obs i,
  `members` aliasing tree.indices, A left unfilled (commit ii).
- ORACLE (component, tests/cpp): build against a live tree grown by the real
  move code (reuse test_moves/test_tree helpers to grow a tree over synthetic
  data). Assert: every non-empty leaf has exactly one atom; the atom's member
  slice equals `tree.indices[leaf.begin..leaf.end)` (memcmp); `atomOf`
  round-trips (obs -> atom -> slice contains obs); atom count == non-empty leaf
  count. This is the snapshot-vs-live oracle from MEMORY sampler-internals; the
  live tree is the oracle, not a keepTrees snapshot.
- GATE: tests/cpp clean (new suite + existing); it compiles into the engine but
  is UNWIRED, so tinytest + equivalence are trivially unaffected. Confirmatory
  bench-sampler compare vs bench-sampler-32fc7c8.csv (should be identical -
  nothing on the hot path yet; this pins the harness/machine before real work).

### (ii) b=1 per-leaf aggregation -> (A, G, Q) bitwise vs the kernel. ~250 lines. 2-3 days.

- FILES: atoms.hpp (`aggregateAtom`, `aggregateTree`, the isRoot/weighted
  dispatch); test_atoms.cpp.
- WHAT: `aggregateAtom(tree, treeY, weights, leafId)` computes (A,G,Q) via the
  SAME `computeLeafStats` kernel dispatch (weighted vs unweighted, root
  non-indexed vs non-root indexed - pins 2,3 of 1.5) over the member slice, and
  fills the SoA. `aggregateTree` walks `fillBottom(0)` (pin 1). Still UNWIRED
  (writes only the SoA, not the node cache).
- ORACLE: for a fuzzed sequence of grown trees + random residuals, assert
  `aggregateAtom`'s (A,G,Q) == the node cache that `computeLeafStats` produced,
  BITWISE (exact `==` on the raw doubles) - INCLUDING the root/stump case and a
  weighted case. This is the tightest possible oracle for the aggregation.
- GATE: tests/cpp bit-exact; tinytest/equivalence unaffected (unwired). No
  bench (off hot path).

### (iii) Draw + residual roll via atoms wired at b=1 (S bookkeeping). ~250 lines. 3-4 days.

- FILES: chain.hpp (the flag + `if constexpr`/`if (useAtomPath)` in the sweep:
  `setNodeAverages` -> `aggregateTree` that ALSO writes the node caches; S(c)
  update at the draw); atoms.hpp (`writeNodeCaches`, `setInBlockFit`);
  test_atoms.cpp / test_sampler.cpp.
- WHAT: flag ON, route the SUFFSTAT source: `aggregateTree` writes the node
  cache (from the SoA) INSTEAD of `computeLeafStats`. The draw, scatter, roll,
  totalFits, sigma stay UNCHANGED (1.4). Add S(c)=mu bookkeeping after the draw
  (design 3.3 at b=1: S is the atom's single-tree fit) - carried but inert at
  b=1 except as a testable quantity. MOVES still use the live `computeLeafStats`
  in this commit (birth/change/swap not yet atom-routed), so only the
  per-sweep `setNodeAverages` is atom-sourced; this isolates the aggregation
  wiring from the move wiring.
- ORACLE: (a) differential - run one sweep with flag OFF and flag ON from the
  same seed/state; assert bitwise-equal treeFits, totalFits, and drawn params
  (test_sampler-style, small n). (b) S-consistency: after the draw+scatter,
  S(c) equals treeFits[i] for every member i of atom c (bitwise).
- GATE: tests/cpp bit-exact; with the flag ON in a dev build, run the FULL
  equivalence anchor (22/22 IDENTICAL, exact-match mode) and full tinytest with
  NO snapshot regeneration - since only the suffstat SOURCE changed and it is
  bitwise, draws are unchanged. This is the first end-to-end bit-identity
  checkpoint. (Ship default stays flag OFF until vii; but gate ON here.)

### (iv) Atom split (birth) + rollback + fuzz vs the live partition. ~350 lines. 4-5 days.

- FILES: atoms.hpp (`splitAtom`, `undoSplit`, `snapshotAtoms`/`restoreAtoms`);
  chain.hpp/moves.hpp hook so birth (flag ON) sources child suffstats from
  `splitAtom` and a rejected birth restores via the atom path; test_atoms.cpp.
- WHAT: on birth of leaf L, `splitAtom` partitions L's atom by the tree's
  `partitionChildren` (DESIGN A: same buffer/primitive, pin 7), aggregates left
  then right child SoA (pin 4), writes the two child node caches, updates
  topology (leafTuple, atomOf, begin/end, id recycling). Rejected birth:
  `undoBirth` already restores tree state; the atom path discards the child
  atoms and restores the parent SoA entry (A served from the U'WU-style cache;
  G/Q recomputed or snapshot-restored). Mirror `snapshotSubtree`/`restoreSubtree`
  (tree.hpp:827-842) for change/swap prep in commit (v).
- ORACLE: extend the mutation-fuzzer pattern (docs/plans/mutation-fuzzing.md,
  testMutationFuzzer): drive long seeded random move sequences; after EVERY
  move (accepted or rejected), assert the atom SoA matches a from-scratch
  `buildForTree`+`aggregateTree` rebuild (bitwise A/G/Q, matching topology, and
  `members`/atomOf consistent with the live tree.indices). A rejected move must
  leave the SoA bitwise-identical to before (the fingerprint pattern). Print
  seed + op trace on failure; sabotage-proof (break a restore, see it fail).
- GATE: tests/cpp bit-exact + fuzzer >= 20 seeds clean; equivalence 22/22
  IDENTICAL + full tinytest NO regeneration (flag ON dev build).

### (v) Death / change / swap in atom terms. ~300 lines. 3-4 days.

- FILES: atoms.hpp (`mergeAtoms`, `refreshSubtree`); moves.hpp/chain.hpp hooks
  for death (orphan/merge), change, swap under the flag; test_atoms.cpp.
- WHAT: death -> `mergeAtoms` (SoA topology merge, A/G/Q left+right per pin 5;
  node cache via unchanged `orphanChildren`, 1.4). change/swap ->
  `AtomMap::refreshSubtree` mirroring `refreshSubtree`'s DFS re-partition +
  re-aggregate (pin 6), with `snapshotAtoms`/`restoreAtoms` for rejection
  (mirror `snapshotSubtree`/`restoreSubtree`). At b=1 with DESIGN A, the member
  re-partition is the tree's own; the SoA topology + stats are the new work.
- ORACLE: the same fuzzer, now with the full move vocabulary weighted (birth,
  death, change, swap) as in testMutationFuzzer; invariant after every op as in
  (iv). Add the specific rejected-change/swap bitwise-stability assertion.
- GATE: tests/cpp bit-exact + fuzzer >= 20 seeds; equivalence 22/22 IDENTICAL +
  tinytest NO regeneration (flag ON dev build).

### (vi) Cross-sweep persistence + patch-on-accept (U'WU template). ~250 lines. 3-4 days.

- FILES: atoms.hpp (A-cache with member-list validation, `lookupA`/`storeA`
  mirroring `lookupCrossproduct`/`storeCrossproduct`); wire it into
  `aggregateTree`/`splitAtom`/`mergeAtoms`; test_atoms.cpp.
- WHAT: persist the atom topology + A across sweeps; on aggregation, serve A
  from the cache when the leaf's ordered member list matches
  tree.indices[begin..end) (memcmp), else rebuild - patch-on-accept and
  rollback-stability fall out of the compare, no explicit accept hook. G/Q/S
  recompute each block entry. Verify persistence saves work (no full rebuild
  per sweep) without changing any value.
- ORACLE: run N sweeps; assert the persisted-A path yields bitwise-identical
  node caches (hence draws) to a force-rebuild-every-sweep control; assert a
  membership change invalidates exactly the changed leaves' A. Fuzzer continues
  across sweeps (persistence stress).
- GATE: tests/cpp bit-exact; equivalence 22/22 IDENTICAL + tinytest NO
  regeneration (flag ON dev build).

### (vii) Full b=1 sweep through the atom path; flip default; MILESTONE. ~150 lines. 2-3 days.

- FILES: chain.hpp (flip the flag default ON for the constant-leaf path;
  ensure vector/function/GP leaf models and the warm-start loop stay on the
  legacy path - the flag is constant-leaf + steady-state-`run()` only);
  benchmarks (record the compare); design note (landing note); docs/plans note.
- WHAT: default the constant-leaf steady-state sweep through the atom path at
  b=1. Keep the flag so OFF remains a one-line abort and Stage B's b-dispatch
  seam is in place. Confirm the non-constant-leaf models are untouched.
- GATE (the Stage-A MILESTONE, ALL must hold, NO re-record):
  - equivalence: `Rscript benchmarks/R/equivalence.R compare
    benchmarks/baselines/equivalence-ac6ec2c.rds` => 22/22 EXACT match
    (identical RNG streams). Any non-exact scenario = defect, STOP.
  - tinytest: full `tinytest::test_package("dbarts")` from a `--preclean`
    install, NO snapshot regeneration (bitwise draws unchanged). All ok.
  - tests/cpp: `cd tests/cpp && make && ./test_bartcore` clean (delete stale
    binaries first - header edits, no dep tracking) + fuzzer >= 20 seeds.
  - bench-sampler: `Rscript benchmarks/R/bench-sampler.R compare
    benchmarks/baselines/bench-sampler-32fc7c8.csv` on the QUIET machine =>
    NEUTRAL (no metric > 5% slower). A real regression = the atom bookkeeping
    overhead is too high = KILL SIGNAL. This is the cheap early abort before any
    Stage B FP change.

## 4. Gates - per step and milestone

Standard gate procedure (CLAUDE.local.md): `R CMD INSTALL .` (`--preclean`
after any header/Makevars/facade change - atoms.hpp is a header, so preclean
every commit that touches it; ALWAYS preclean since chain.hpp includes it),
tinytest from the installed package, tests/cpp component tests (delete stale
`tests/cpp/*.o` + `test_bartcore` after header edits), equivalence, bench.

- EVERY commit (i-vii): tests/cpp clean (new + existing).
- Commits (iii-vii) [flag ON dev build]: equivalence 22/22 EXACT +
  full tinytest NO regeneration. Because the change is bit-identical by
  construction, these must pass at EVERY wired commit, not just the milestone -
  do not defer end-to-end checks.
- Commit (i): confirmatory bench compare (harness/machine pin).
- Milestone (vii): the full gate in 3(vii), bench on the quiet machine.
- rchk: not required (no PROT/SEXP machinery touched; pure engine + tests).

RNG-NEUTRALITY TRAPS to re-check at each wired commit (from 1.5): summation
order (call the kernel, never re-sum), leaf/draw order (`fillBottom`, not
atom-id order), the root non-indexed dispatch, birth left-then-right, death
left+right, the empty-leaf no-draw, and NO gratuitous g field at b=1. A single
misplaced accumulation is the whole failure mode; the per-tree differential
oracle (iii) and the fuzzer (iv-vi) are designed to catch it at the leaf level.

## 5. Risks and open items

### 5.1 Implementer-decidable (no VD call)

- Member buffer alias vs own (2.4): RECOMMEND alias (DESIGN A) for Stage A;
  decide by the bench/memory gate; record in the landing note.
- Flag mechanism: `constexpr bool` on the constant-leaf sweep vs a build macro
  vs a runtime Forest field. RECOMMEND a `constexpr`/`if constexpr` seam gated
  on the leaf-model trait so the legacy path is fully compiled out in the
  shipped default AND the atom path is compiled in dev/test builds for the
  differential oracle - i.e. a build knob (e.g. a config.hpp/Makevars define)
  that test builds set, defaulting to the atom path at (vii). Keep both
  compilable so tests/cpp can instantiate both writers side by side.
- Atom id recycling / free-list details, snapshot granularity, K-vector
  initial capacity: local, testable.

### 5.2 Risks

- BENCH REGRESSION AT b=1 (the designed kill signal). The atom SoA + its
  per-move maintenance is pure overhead at b=1 (no b>1 amortization). DESIGN A
  minimizes it (no member copy, aggregation REPLACES computeLeafStats). If it
  still regresses > 5%, that is the approach telling you the bookkeeping is too
  heavy - abort cheaply (flag OFF), before Stage B. Mitigation: keep the atom
  interior O(K); do not add any O(n) pass beyond the ones b=1 already pays.
- BITWISE FRAGILITY. The anchor lives or dies on summation/leaf/root order
  (1.5). Mitigation: the tight per-tree differential (iii) + the cross-sweep
  fuzzer (iv-vi) localize any divergence; the "call the kernel, never re-sum"
  rule removes the largest class.
- COMPLEXITY / BUG SURFACE. A second bookkeeping structure with
  split/merge/snapshot/rollback. Mitigation: DESIGN A defers the owned buffer;
  persistence reuses the landed U'WU-cache validate-by-memcmp (rollback-stable,
  no per-move hook); the fuzzer is sabotage-proofed (mutation-fuzzing precedent).
- SUBSTRATE DRIFT. Stage A follows data-ownership (plans 1-3 landed:
  width-templated `partitionChildren`, u8-on-scalar). Line numbers in this plan
  are from the post-data-ownership tree (worktree HEAD 64ef066); the design
  doc's refs are from 2e2b1c9 (pre). Mitigation: the implementer greps the named
  functions rather than trusting line numbers; the partition primitive the atom
  split reuses is the settled width-templated one (design 5.2).
- ENGINE FEATURES THAT MUST STAY ON THE LEGACY PATH at b=1: vector-param
  (linear) and function-param (GP) leaf models (design 5.3 - constant-leaf
  only), the grow-from-root warm start (design 5.3), and any path that is not
  the steady-state constant-leaf `run()` sweep. The flag must be scoped so these
  are provably untouched (they read the same node cache; do not route their
  suffstats through the constant-leaf atom SoA). Categorical / sparse / MIA /
  weighted / probit-logistic latents / BCF-grouped ARE covered at b=1 (design
  5.3) because the atom aggregation just calls the same kernel over the same
  members and g rebuilds from the working response each sweep.

### 5.3 Open questions for VD

NONE surface for Stage A. It is bit-identical, no re-record, no new public
surface, no ABI change; every choice above is implementer-decidable and gated
by bitwise equality. The VD-level questions in design 8.2 (default b, re-record
bundling, equivalence anchor under a b knob) all belong to Stage B, which
carries the one approved re-record. If the milestone bench regresses, that is a
KILL/scope decision, but it is a gate outcome, not a design question - surface
it to VD only if it fires.
