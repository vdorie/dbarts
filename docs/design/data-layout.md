# Structural contiguous-per-node data layout

Status: CLOSED - SHELVED after re-evaluation, 2026-08-04 (see the
post-mortem below; was DESIGN PROPOSAL 2026-07-11). Co-designed with the data-ownership program (docs/design/
data-ownership.md; plan-1 landed the owned predictor ColumnStore). This
note covers the RESPONSE-axis arrays (y/fits/weights/residual), a
different axis from plan-1's predictor CODES. It is the write-up of the
"store per-node y/fits in index-buffer (contiguous) order" lever flagged
as the highest-value structural direction by the two SIMD investigations
(docs/plans/x86-simd-plan.md, docs/plans/simd-survey.md).

POST-MORTEM (2026-07-21): the block-fused atom map this note points to as
the payoff (section 6; parallel-bart-frontier.md 3.4) is CLOSED, WONT-DO
(docs/design/block-fusion.md) - Stage B measured the fused path 4-6x
SLOWER than b=1, and the machinery was excised. So the "this layout is
the substrate for block-fusion" framing below no longer has a downstream
win to feed. What remains is the ~10% standalone single-tree reorder
(section 3): either it gets re-evaluated on its own merits (weighed
against its re-record and complexity costs, section 8), or this note is
shelved.

RE-EVALUATED 2026-08-04 (dbarts-bench re-profile at 06f73b0;
memory-wall-frontier.md section 10): SHELVED, CLOSED. The reorder's
touchable share is the residual-permutation-bound gathers, 46-49% of a
sweep - but the ~15% fit scatter this note planned to amortize the
reorder over is gone (the leafOf/muByTree compaction deleted it), the
21-22% partition/compare work is on predictor codes and does not
compound (section 5's own point), and the reorder still pays a gather
of the same magnitude as the one it removes. With one consumer fewer,
the realistic win falls below the already-marginal ~10% of section 3
against unchanged re-record and complexity costs. If the
histogram-fused suffstat (memory-wall-frontier.md sec 2a) ever ships,
its fused pass reads the residual in place and removes this note's
remaining premise entirely.

The headline of this note, stated up front because it corrects the
premise: the naive per-tree reorder is NOT a clean ~2.6x-on-32% win. The
"2.6x" is the per-element reduce ratio between contiguous and gather
access; it silently ignores that making the data contiguous costs a
gather (the reorder pass) of essentially the same magnitude as the
suffstat gather it removes. The residual's cross-tree coupling forbids
maintaining leaf order incrementally, so a per-tree O(n) permutation is
unavoidable every sweep. The realistic single-tree win is ~10% of sweep
time and ONLY materializes with the SIMD-lane reduce (which needs a
re-record). The genuinely large lever on this same primitive is the
block-fused atom map (parallel-bart-frontier.md 3.4), for which this
layout is the substrate.

## 1. The engine's current access pattern (verified in code)

Per forest (combiner.hpp:138-197):

- `indexBuffer` : `n * numTrees` size_t. Tree t owns the contiguous slice
  `indexBuffer + t*n`; `Tree::indices` points at it (tree.hpp:260,268).
  It is a PERSISTENT per-tree permutation P_t: initialized to identity
  (tree.hpp:224), partitioned in place by moves (tree.hpp:837
  partitionChildren), and carried across sweeps. Each node owns a
  contiguous RANGE `[node.begin, node.end)` of P_t (tree.hpp:221); the
  node's members are `indices[begin .. end)`, which are SCATTERED
  positions in 0..n-1.
- `treeFits` : `n * numTrees` doubles. Tree t's fitted contribution, in
  OBSERVATION order (treeFits[obs]). Persistent.
- `treeY` : ONE `n` double buffer, the running residual, OBSERVATION
  order, reused across all trees, rolled incrementally (chain.hpp:
  728-744).
- y, weights: obs order, owned by the response family (external to the
  forest).

The per-sweep, per-tree O(n) passes and their memory access shape:

1. Residual roll (chain.hpp:1457-1470): `treeY[i] += oldFits[i] -
   prevFits[i]` etc. Fully CONTIGUOUS (obs order); 3 streams; already
   near-optimal, bandwidth-bound (~15% inside Chain::run's 28.7% self).
2. setNodeAverages -> computeLeafStats (tree.hpp:725-760): per leaf,
   `misc_computeIndexedSufficientStatisticsFast(treeY, indices+begin,
   len, ...)`, i.e. read `treeY[indices[begin+k]]` for k in 0..len and
   accumulate (sumW, sumWZ, sumWZ2). This is a GATHER over the residual
   (32% of runtime, the top hotspot). It reads the RESIDUAL, which
   changes every tree, not the static y.
3. metropolisJumpForTree scoring (moves.hpp): a proposed birth partitions
   one node (partitionIndices GATHER, part of the 22.5% partition) and
   computes 2 child suffstats (more gathers); change/swap re-partition
   near the root (~n) and snapshot the affected index segment for
   rollback (tree.hpp:1016).
4. sampleParametersAndSetFits (chain.hpp:4896-4990): per leaf, draw a
   constant param, then `misc_setIndexedVectorToConstant(treeFits,
   indices+begin, len, param)` -> SCATTER of a per-leaf constant into
   obs-order treeFits (14.9%).

So of the four top hotspots, THREE (suffstat gather 32%, fit scatter 15%,
partitionIndices gather 8%) are permutation-bound on P_t, and one (roll
15%) is contiguous. Every one is memory-bound; SIMD width buys them
~0-1.15x (measured). The only structural escape is to change the LAYOUT
so the permutation-bound reads/writes become sequential.

## 2. Proposed layout

### 2.1 What moves into leaf-contiguous order

Define, for a tree t with permutation P_t and a leaf L owning range
`[begin, end)`, a "leaf-order" array A' with `A'[k] = A[P_t[k]]`. In
leaf order, leaf L's members are the CONTIGUOUS slice `A'[begin .. end)`.

Candidates and their character:

- residual (treeY): the hot gather target. CHANGES EVERY SWEEP. Cannot be
  held in leaf order incrementally (section 3). Must be produced in leaf
  order each sweep by a gather -> this is the unavoidable reorder pass.
- fits (treeFits): written per leaf as a constant. In leaf order the
  write is a sequential constant-fill (near free). BUT the roll needs
  fits in obs order (it sums across trees, section 3), so a leaf-order
  fit must be scattered back -> the unavoidable scatter pass.
- y, weights: STATIC within a sweep (change only on a latent refresh).
  Can be held persistently in leaf order, re-permuted only when the
  partition changes (O(changed), cheap). But the constant leaf's hot
  suffstat reads the RESIDUAL, not y; y/weights matter only for the
  weighted variant and the linear/GP leaves.

### 2.2 Two maintenance strategies

(1a) Physically permute alongside the partition. Keep leaf-order COPIES
of the arrays and, whenever partitionChildren swaps two index-buffer
slots, swap the corresponding double(s) too. y/weights (static) then live
permanently in leaf order at O(changed)-per-move maintenance. Costs:
+O(n*m) memory per array held this way (a leaf-order residual copy is
another n*m doubles == as large as treeFits itself: 160 MB at n=1e5,
m=200); partition swaps move 8-24 extra bytes each; change/swap rollback
must snapshot the double segments too (tree.hpp:1016 currently snapshots
only the size_t indexSegment). It does NOT solve the residual: the
residual VALUES change every sweep even when P_t is unchanged, so the
leaf-order residual must be re-gathered every sweep regardless -> 1a
helps only the static y/weights, which are not the hot path.

(1b) Gather-once-per-sweep into a reused n-scratch. Per tree, after the
roll, gather the residual into a single reused `rscratch[k] =
treeY[P_t[k]]`; run suffstat, move scoring, and the fit draw sequentially
over rscratch and a leaf-order fit scratch; scatter the fits back to
obs-order treeFits for the roll. Costs: +O(n) memory (negligible); one
gather + one scatter per tree. This is the tractable version and the one
recommended below. Its economics are analyzed in section 3.

Under either strategy the per-node suffstat becomes
`misc_computeSufficientStatisticsFast(rscratch + begin, len, ...)` (the
CONTIGUOUS kernel, no indices) and the fit write becomes
`misc_setVectorToConstant(fitScratch + begin, len, param)` (a sequential
fill). Both are sequential contiguous accesses.

## 3. THE ECONOMICS (why the naive framing overstates the win)

The residual is the crux. The roll computes tree t's residual as
`y - sum_{k != t} fits_k`, an inherently CROSS-TREE quantity. Consecutive
trees have DIFFERENT permutations (P_t != P_{t+1}), and the residual is
rolled between them, so there is no leaf order common to two trees and no
way to carry a leaf-order residual incrementally: any scheme that keeps a
per-tree leaf-order residual must, to roll from t to t+1, gather the
just-drawn fits_t (or the whole residual) into P_{t+1}'s order -> a
per-tree O(n) gather. The permutation is unavoidable every sweep.

Per-element costs (docs/plans/x86-simd-plan.md microbench, ns/elem, n=50000):

- suffstat, gather (current):        0.667
- suffstat, contiguous scalar:       0.258-0.284
- suffstat, contiguous SIMD (4x ILP):0.065-0.076
- reorder pass (gather-copy + write):~0.5 (gather overhead ~0.4 + write)
- fit scatter (current):             ~0.5 (setIndexedVectorToConstant)
- fit sequential fill:               ~0.1

Suffstat ALONE does not win. Reorder(0.5) + SIMD-contiguous-reduce(0.076)
= 0.576 vs current gather-reduce 0.667 -> only ~1.16x, because the reorder
pays the very gather the reduce avoided. With the SCALAR contiguous
reduce (order-preserving, section 4): 0.5 + 0.28 = 0.78 > 0.667 -> a LOSS.

The win requires AMORTIZING the single reorder across the three
permutation-bound consumers within one tree's turn (suffstat, move
scoring, fit write), and it requires the SIMD reduce:

  current (per elem-equiv, summing the permutation-bound passes):
    suffstat 0.667 + moves ~0.3 + fit 0.5            = ~1.47
  redesign (SIMD reduce):
    reorder 0.5 + scatter-back 0.5 + suffstat 0.076
    + moves 0.076 + fit 0.10                          = ~1.25
  -> ~1.18x on the ~57% of runtime these passes occupy = ~9% of sweep.

  redesign (order-preserving SCALAR reduce, no re-record):
    0.5 + 0.5 + 0.28 + 0.28 + 0.10                    = ~1.66  -> LOSS.

CONCLUSION: the realistic single-tree win is ~10% of sweep time, an order
of magnitude below the "2.6x on 32% = 20%" headline, and it exists ONLY
with the SIMD-lane reduce. The order-preserving (zero-re-record) variant
is a performance LOSS; it is a correctness-preserving fallback, not a
win.

### Break-even in n / p / tree-count

- n: the gather penalty (0.667 - 0.28 = ~0.39 ns/elem) is cache-miss
  driven and GROWS with n (measured 0.35 -> 0.61 -> 0.68 as N goes
  2e3 -> 5e4 -> 5e5). At small n (<~1e4, cache-resident) the penalty is
  small and the reorder overhead dominates -> NET NEGATIVE. Favorable
  only at n >= 1e5 (the confirmed-common regime).
- p: higher p means more candidate-split move scoring, more of which
  becomes sequential under the reorder -> MORE benefit. Favorable at
  high p.
- tree-count m: reorder passes and suffstat both scale with m ->
  ratio-neutral, but absolute memory (1a) scales n*m.

Verdict: NET POSITIVE only in the large-n, high-p, SIMD-reduce-enabled
corner, and even there ~10%. The scatter-back is not the dominant cost;
the reorder gather is. The honest lever here is not the single-tree
reorder but the block-fused atom map that reuses this layout (section 6).

## 4. THE RNG / BIT-IDENTITY DECISION (the crux)

dbarts's equivalence gate is machine-independent BY DESIGN: every
dispatched double kernel on the draw path is elementwise or a
permutation, and the hot REDUCTIONS (suffstat, SSR) are kept SCALAR and
undispatched so their fixed summation order is reproduced on every ISA
(simd-survey.md "key architectural invariant"). Physically reordering the
residual changes WHERE elements live, so we must ask precisely whether it
changes the ORDER in which they are summed.

### 4.1 The layout reorder is inherently ORDER-PRESERVING

The current gather kernel reads `treeY[indices[begin+k]]` for k = 0,1,2,...
in ascending k (with a fixed mod-5 unrolling) and accumulates in that
order. If the reorder fills `rscratch[k] = treeY[indices[begin+k]]` in
ascending k, then the contiguous kernel reads `rscratch[0,1,2,...]` in the
SAME order with the SAME unrolling. The reorder is a pure copy (it sums
nothing); the subsequent SCALAR contiguous reduce sums in identical
index-buffer order. Therefore:

  reorder + misc_computeSufficientStatisticsFast (scalar) is BITWISE
  IDENTICAL to the current gather kernel. Zero re-record.

The two-pointer partition is unstable, but that only sets what P_t IS;
whether we read via a gather or via a pre-gathered contiguous copy, the
summation visits the same elements in the same P_t order. Order is
preserved by construction.

### 4.2 ...but the order-preserving version does not win

From section 3, reorder + SCALAR reduce is a net LOSS (0.78 vs 0.667 on
suffstat; ~1.66 vs ~1.47 across the block). The 4x that makes the layout
worthwhile comes from the SIMD contiguous reduce, which sums in LANE
order (4 interleaved partial sums, combined by a fixed tree) != the
scalar sequential order. That reassociation changes the FP result ->
changes draws -> forces a one-time equivalence re-record, and (unlike the
partition) the result is not reproducible by a scalar fallback in the
same order.

### 4.3 Options and recommendation

(a) Accept a one-time re-record for the win. Ship reorder + the
    FIXED-LANE SIMD contiguous reduce. The lane structure (lane count,
    element->lane map, combine tree) is fixed across scalar/NEON/SSE2/
    AVX2 so the result stays machine-INDEPENDENT (same bytes on every
    ISA) even though it differs from today's scalar order. One re-record;
    machine independence PRESERVED. This is exactly the fast/reference
    TOGGLE VD already specified for the suffstat kernel (x86-simd-plan.md
    "bit-identity is a TOGGLE"): the scalar sequential kernel remains the
    reference for exact reproducibility; the SIMD path matches it only
    distributionally. The layout redesign SHOULD ride that same toggle
    and re-record.

(b) Order-preserving layout (scalar reduce, zero re-record). Correct and
    machine-independent-by-construction, but a PERFORMANCE LOSS
    (section 3). Not worth building for speed; only of interest if the
    layout were wanted for a non-speed reason (e.g. as the atom-map
    substrate, section 6, where the reduce is not the point).

(c) Hybrid. Land the layout plumbing behind the toggle: default install
    the SIMD reduce (path (a), re-record once); keep the scalar
    sequential reduce as the toggle's reference, which is bit-identical
    to TODAY's kernel (4.1) and therefore doubles as the migration
    anchor -- the re-record can be validated by flipping the toggle and
    confirming the reference reproduces the OLD baseline exactly.

RECOMMENDATION: (a)/(c). The re-record is REQUIRED for any win -- there is
no zero-re-record version that speeds this path up. Because the fixed-lane
kernel keeps machine independence, the re-record is a one-time snapshot
refresh (the "shifting" RNG class), NOT a surrender of cross-machine
reproducibility. That is a materially smaller ask than simd-survey.md's
worst case (ISA-dependent draws), and it is the ask VD already accepted in
principle via the toggle. SURFACE IT: this is a VD sign-off, not an
implementer call.

## 5. The scatter cost: verdict

The residual refresh forces a per-tree gather (reorder in) and the roll
forces a per-tree scatter (fits back to obs order). Neither is removable
while the roll stays cross-tree in obs order (section 3). Contrary to the
survey's framing, the SCATTER-back is NOT the dominant new cost -- the
reorder GATHER is (0.5 vs 0.5 per elem, but the gather also displaces the
suffstat's own gather, whereas the scatter is pure addition). Net across
the block: ~1.18x with SIMD, i.e. ~+10% of sweep, and only at large n/high
p. NET POSITIVE, but modestly, and conditional on the re-record. At small
n it is NET NEGATIVE (reorder overhead exceeds the shrunken gather
penalty). This is a real-but-minor lever as a standalone item; its value
is as the primitive the block-fused approach needs (section 6).

## 6. Co-sequencing with the data-ownership program

### 6.1 Scope boundary (clean, different axis)

Plan-1's ColumnStore owns predictor CODES (u8/u16), indexed by column and
consumed by the PARTITION. This redesign is about the RESPONSE-axis
per-observation arrays (treeY/treeFits/y/weights), owned by the Forest and
the response family, indexed by observation and permuted per-tree. The two
programs touch DISJOINT storage. Plan 2 (frame-direct ingestion, drop @x)
is entirely predictor-side; it has nothing to say about treeFits/treeY
layout, and this redesign has nothing to say about the frame.

### 6.2 The one shared touchpoint: partitionChildren

The single place the two axes meet is partitionChildren (tree.hpp:837),
which permutes the index buffer using the codes. Under strategy 1a it
would ALSO move the response-axis doubles alongside; and the data-
ownership program is ALREADY reworking this function (plan-1 step 3-4:
width-templated u8/u16 partition, the epi8 kernel family). If both land
independently, partitionChildren is reworked twice and the second rework
fights the first.

### 6.3 Recommendation: FOLLOW the program; do not fold into plan 2

1. This redesign FOLLOWS the data-ownership program (after plan 3/4),
   not folded into plan 2. Plan 2 is predictor ingestion; forcing the
   response-axis layout into it would conflate two orthogonal concerns
   and bloat a plan already at ~1200-1800 lines.
2. Plan 2's container need NOT anticipate the layout: the response-axis
   arrays are Forest-owned and the redesign is self-contained there. The
   ONE thing the program should preserve (it already does) is the clean
   index-buffer + node-range abstraction, and it should NOT bake in an
   assumption that response arrays are permanently obs-order in a way
   that is hard to relax.
3. Sequence the partition work so the width/code axis (plan-1 step 3-4)
   SETTLES FIRST, then this redesign extends the settled partition to
   optionally carry a response-axis payload (strategy 1a) or leaves
   partition untouched and gathers into scratch (strategy 1b -- which
   does NOT touch partitionChildren at all, a strong reason to prefer
   1b: it decouples entirely from the data-ownership partition rework).
4. Memory interaction: plan-1's u8 codes shrink predictor memory (~8x),
   but strategy 1a's leaf-order residual copy ADDS n*m doubles (160 MB
   at n=1e5,m=200), dwarfing that saving. Strategy 1b adds only O(n).
   Another vote for 1b.

Net: 1b (scratch gather) is the recommended strategy precisely because it
is ORTHOGONAL to the data-ownership partition rework -- it needs no change
to partitionChildren, no per-tree double copies, and can land any time
after plan-1 without waiting on plans 2-5.

## 7. Threading interaction

Contiguous-per-node layout HELPS within-chain threading (docs/plans/
within-chain-threading.md) and the block-fused frontier:

- The within-chain plan parallelizes the O(n) passes with a FIXED-BLOCK
  reduction for thread-count invariance. Sequential contiguous buffers
  partition cleanly into fixed blocks with predictable per-worker memory
  traffic and no false sharing; a gather/scatter over a shuffled buffer
  gives each worker scattered, unpredictable cache lines. So the reorder
  makes the suffstat/fit passes BETTER threading targets.
- The reorder gather and scatter-back are themselves O(n) memory-bound
  passes that would compete for DRAM bandwidth with the other threaded
  passes; they must be counted in the threading budget (do not double-
  count the win, per the within-chain plan's note on the u8 interaction).
- The fixed-lane SIMD reduce (section 4a) and the fixed-block threaded
  reduction share the SAME determinism discipline (fixed accumulator
  layout, combine in block order); building one eases the other, and
  both re-record together as one shifting-class change.
- Block fusion (frontier 3.4) COMPOSES with threading (fewer barriers AND
  cheaper passes) and with a GPU backend (resident contiguous data). The
  contiguous layout is the shared substrate for all three; a GPU seam
  (frontier section 4) wants exactly this resident-contiguous shape.

Verdict: contiguous layout is SYNERGISTIC with threading and the frontier,
not in tension. This raises its strategic value above its ~10% standalone
number.

## 8. Expected win, risks, open questions

### 8.1 Expected win (with uncertainty)

- Standalone single-tree reorder (strategy 1b + SIMD reduce + re-record):
  ~+10% of sweep time (range ~5-15%), at n >= 1e5 and moderate/high p.
  NET NEGATIVE below ~n=1e4. The survey's ~20% is optimistic; it counts
  the reduce ratio without the reorder cost.
- Order-preserving (no re-record): NET LOSS. Do not ship for speed.
- The compounding items the survey hoped for (fit scatter 15%, partition
  gather 8%) partially fold in: the fit write goes sequential (a real
  save) but the scatter-BACK reappears; the partitionIndices gather is on
  CODES, unaffected by residual layout, so it does NOT compound here.
- The LARGE win on this primitive is the block-fused atom map (frontier
  3.4, E1/E2 cleared 2026-07-08): ~6x DRAM reduction on the ~85% field,
  i.e. multi-x sweep speedups, by amortizing the reorder across a block
  of b~4-8 trees so the per-tree permutation is paid once per block, not
  once per tree. That is where the real perf lives; this note's layout is
  its substrate and its cheaper down-payment.

### 8.2 Risks

- Correctness of a physical-permutation partition (1a): partitionChildren
  and the snapshot/restore for rejected change/swap moves (tree.hpp:
  765-782) currently move only size_t indices; carrying doubles alongside
  doubles the state to snapshot and is a fertile bug source. Strategy 1b
  sidesteps this (rejected moves re-derive the affected rscratch range
  from treeY, or snapshot the small affected segment).
- The re-record (section 4): a one-time shifting-class snapshot refresh;
  RNG-locked tinytest snapshots regenerate, equivalence anchor re-records,
  z-mode passes. Machine independence preserved by the fixed lane layout.
- Memory blowup (1a): +n*m doubles per leaf-order array; avoid via 1b.
- Complexity: the elegant single shared-residual roll is wrapped in
  per-tree reorder machinery; the code gets harder to read for a ~10%
  win. Weigh against just building the atom map, whose complexity buys
  far more.
- Measurement risk: the ~10% is modeled from x86 microbenches, not an
  end-to-end sampler compare. A confirmatory bench-sampler run on the
  quiet machine is mandatory before committing (the plan-process hot-path
  rule).

### 8.3 Open questions for VD

1. STRATEGIC: is a ~10% standalone win worth the complexity + re-record,
   or should the effort skip straight to the block-fused atom map
   (frontier 3.4), which subsumes this layout and delivers the real
   multi-x win? Recommendation: treat the contiguous layout as the
   SUBSTRATE for 3.4 (and for threading and a GPU seam), and justify it
   by 3.4's payoff, not by the ~10% standalone. Do not ship the single-
   tree reorder as a standalone perf item unless 3.4 is deferred.
2. RE-RECORD: accept the one-time re-record (fixed-lane SIMD reduce,
   machine independence preserved)? There is no zero-re-record version
   that wins. This rides the fast/reference toggle VD already approved.
3. STRATEGY: 1b (scratch gather, O(n) memory, no partitionChildren
   change, orthogonal to data-ownership) vs 1a (physical permute, O(n*m)
   memory, entangled with the partition rework)? Recommendation: 1b.
4. SEQUENCING: confirm this follows the data-ownership program (after
   plans 3/4 settle the partition) rather than folding into plan 2, and
   that plan 2 need only preserve the index-buffer/node-range abstraction
   (no anticipatory container work).
5. SCOPE of leaf-order arrays: constant-leaf gaussian reads only the
   residual; linear/GP leaves read y/weights/covariates and have the U'WU
   cache. Is the layout scoped to the constant leaf first, or must it
   generalize across leaf models before it is worth the plumbing?
