# dbarts SIMD candidate survey (READ-ONLY)

Job 120e5d72. Engine at "step-3" (per-column u8/u16 code-width templating in
scan.hpp/tree.hpp). Baseline pre-step-3 = 2e2b1c9. Dev box arm64; x86 deferred.
No repo files modified.

Delivery vehicle: runtime ISA dispatch in src/misc/simd.c (function pointers
installed by misc_simd_init at load). Hand kernels live in src/misc/*.c and are
dispatched; bartcore/*.hpp headers get only CRAN-default -O2 (x86-64 baseline ==
SSE2, no AVX; arm64 baseline includes NEON).

## THE KEY ARCHITECTURAL INVARIANT (drives every verdict below)

The MCMC draw path is deliberately built so that **every dispatched double
kernel is elementwise or a permutation** - never a reduction:

- Dispatched today ([[simd.c:212-417@4a521760]]): partition (permutation, integer-exact),
  and the linalg family add/subtract/AXPY/addScalar/setConstant/transpose - all
  ELEMENTWISE (each output element = fixed arithmetic on the same operands,
  independent of vector width).
- Elementwise + permutation kernels are **bitwise-identical across ISAs** by
  construction. That is WHY the equivalence gate (equivalence-ac6ec2c.rds, 22
  scenarios, bitwise) is machine-independent: partition_*.c reproduce the scalar
  two-pointer permutation exactly ([[tree.hpp:626-628@4a521760]] "bitwise identical to the
  misc kernel"), and elementwise double ops don't reassociate.
- The hot REDUCTIONS that determine draws - leaf sufficient statistics
  (misc_compute*SufficientStatisticsFast, [[moments.c:311-416@4a521760]]) and the sigma-draw
  SSR (misc_computeSumOfSquaredResiduals, [[moments.c:2768@4a521760]]/2789) - are kept
  **SCALAR and UNDISPATCHED** on purpose. Comment [[moments.c:307-310@4a521760]] even
  pre-marks a SIMD suffstat as future work ("raw sums are order-insensitive").
- The dispatched moment reductions (Mean/Variance, moments_sse2.c) are NOT in
  the draw inner loop (see hotness table); they serve standalone/mt callers.

CONSEQUENCE: SIMD-izing + dispatching any draw-path reduction makes draws
**ISA-dependent** (kernel-vocabulary.md explicitly waives cross-variant
bitwise equality for the moment family). That is strictly bigger than a one-time
re-record: it ends machine-independent bitwise reproducibility unless the kernel
uses a FIXED cross-ISA accumulator layout (same lane count / combine tree on
scalar+NEON+SSE2+AVX2) plus a one-time anchor re-record. This is the central
tension: the clean elementwise wins are already taken; what's left hot is
reductions (bitwise-blocked) and scatters (hard to SIMD).

## Hotness map (default constant-leaf BART, per MCMC sweep)

Per tree, per sweep, the O(n) passes (n = numObservations) are:
1. Residual roll ([[chain.hpp:728-786@4a521760]]): 3 fused elementwise loops, HAND-WRITTEN in
   the header (not misc). `resid[i]=y-total+oldFits` / `resid[i]+=oldFits-prevFits`
   / `total[i]=y-resid+lastFits`.
2. setNodeAverages -> computeLeafStats over every leaf ([[tree.hpp:492-521@4a521760]]): a full
   GATHER+REDUCE pass over all n (y[indices[i]]), because backfitting changes the
   residual every sweep so cached leaf sums are stale. leafTracksNodeAverages =
   !hasVectorParams = TRUE for the constant leaf ([[chain.hpp:365@4a521760]], [[chain.hpp:747@4a521760]]).
3. Fit scatter (sampleParametersAndSetFits, [[chain.hpp:2076-2129@4a521760]];
   setTreeFitsFromParameters [[chain.hpp:2035@4a521760]]): per leaf, misc_setVectorToConstant (root,
   dispatched) or misc_setIndexedVectorToConstant (scatter of a constant, NOT
   dispatched).
Per-MOVE costs (metropolisJumpForTree) are O(n_node): birth partitions one node
(partition, dispatched u16 / scalar u8) + 2 child stats; death is incremental
(orphanChildren, [[tree.hpp:792-802@4a521760]]). So the O(n)-per-tree passes above dominate
for large n.
Once per sweep (not per tree): SSR for sigma ([[chain.hpp:908@4a521760]]), reduction, scalar.
Recorded sweeps only: test-fit descent (findBottomNodeForRow) - branchy tree
walk, not SIMD-able; accumulation already dispatched.

Cold in the default sampler: the cut-scan histogram (scan.hpp) - see #C.

---

## RANKED CANDIDATES

### #1  Fused residual-roll kernel family  ([[chain.hpp:728-736@4a521760]], 748-750, 780-786; grow twin [[chain.hpp:976-984@4a521760]], 1008-1009)
1. Loop/file: the three hand-written fused elementwise loops in run()'s per-tree
   backfit (and the duplicated pair in growForestFromRoot).
2. Frequency: **per tree, per sweep** - one of the two hottest O(n) passes.
3. Codegen: header-only -> CRAN-default -O2. On arm64 NEON is baseline so these
   ALREADY vectorize to NEON (no arm gap). On x86 they autovec to **SSE2 only**
   (2 doubles); AVX/AVX2 (4 doubles) is never emitted without -mavx.
4. Access: 2-3 contiguous loads + 1 contiguous store per element. No gather.
5. Why CRAN can't get it: not the compiler's fault on shape - it vectorizes
   fine; the gap is purely the missing -mavx width on x86. The existing misc
   vocabulary has no fused 3-input op (only y+=x, y-=x, y+=a*x), so the code was
   hand-fused to save a memory pass (comment [[chain.hpp:726@4a521760]] "One fused pass").
6. Benefit: order-of-magnitude MODEST and x86-only. Cache-resident regime
   (n<=~1e4, treeFits slab ~80KB fits L2): SSE2->AVX2 ~up to 2x on this pass.
   DRAM-bound (n>=1e5): bandwidth-limited, width barely helps. Unmeasurable on
   the arm64 dev box (already NEON); defer to x86.
7. ISA portability: maps cleanly to AVX/AVX2 (the only place it helps). NEON/SSE2
   give nothing over the existing autovec. So really an "AVX/AVX2 fill" kernel.
8. Draw-neutrality: **BITWISE-SAFE and ISA-independent** - pure elementwise, no
   reassociation. No re-record.
9. Effort: LOW. Add misc_fusedResidualRoll-style entries (dispatched AVX/AVX2 +
   scalar/SSE2/NEON = existing autovec), route the header loops through them.
10. u8 interaction: NONE (operates on doubles + no codes).
=> Cleanest safe win, but small and x86-gated. Recommend only after an x86
   microbench confirms a cache-resident delta.

### #2  Close the two known misc dispatch GAPS  (from docs/plans/archive/x86-simd.md)
Two "already-dispatched-with-a-gap" items, both bitwise-safe:
(a) misc_setIndexedVectorToConstant ([[linearAlgebra.c:168@4a521760]]) - the fit SCATTER of a
    constant to a leaf's member indices; used every sweep in the fit-scatter pass
    ([[chain.hpp:2048@4a521760]], [[chain.hpp:2125@4a521760]]). **Never entered the dispatch table** (x86-simd.md).
    - Access: scatter-store of one broadcast constant to data-dependent indices.
    - Why CRAN can't: data-dependent scatter; no SSE2/NEON scatter instruction.
    - Benefit: modest. Only AVX512 has scatter; on AVX2/NEON/SSE2 the best is
      scalar stores (indices within a leaf are NOT monotone, so no contiguous
      store fast-path in general). Bitwise-safe (pure stores).
    - Effort: LOW to just wire a scalar reference into the table (completeness);
      a real SIMD leaf is only worth it on AVX512, low priority.
(b) The three linalg arms whose SSE2/AVX bodies are scalar loops with intrinsics
    commented out (addVectorsInPlace/subtractVectorsInPlace/setVectorToConstant,
    linearAlgebra.c ~:173-215): the originals used non-temporal stores
    (_mm_stream_pd) that autovec won't emit - plausibly real in the DRAM-bound
    n>=1e5 regime. Restore-or-delete per x86 measurement. Bitwise-safe
    (elementwise). These are ON the hot residual/fit path (all dispatched).
Draw-neutrality: both bitwise-safe & ISA-independent. Effort: LOW-MED. x86
measurement required to choose restore-vs-delete; maintainer-run.

### #3  Leaf sufficient-statistic gather+reduce  ([[moments.c:311-416@4a521760]], via setNodeAverages/computeLeafStats)
1. misc_compute[Indexed][Weighted]SufficientStatisticsFast: scalar unrolled-by-5,
   3 accumulators (sumW, sumWx, sumWx^2).
2. Frequency: **per tree, per sweep** (setNodeAverages recomputes ALL leaves) -
   THE hottest reduction in the draw path, co-equal with #1 in raw O(n) work.
3. Codegen: scalar, undispatched (deliberately, per the invariant above).
4. Access: Indexed variant is a **GATHER** (x[indices[i]], y in obs-id order,
   indices a permutation) + reduce. Root/contiguous variant has no gather.
5. Why CRAN can't: (i) it's a reduction the compiler won't reassociate without
   -ffast-math; (ii) even hand-written, the GATHER dominates - AVX2 vgatherqpd
   is often no faster than scalar loads, and NEON/SSE2 have no gather at all, so
   the vectorizable multiply-add tail is a small fraction of the cost.
6. Benefit: LOW on most of the ISA family precisely because it's gather-bound;
   the contiguous (root) case would vectorize well but is rare/cheap. A real win
   needs an algorithmic change (maintain y permuted into index order to make the
   load contiguous) - out of scope for a drop-in kernel, and it trades the
   gather for a scatter on every residual refresh.
7. ISA portability: multiply-add reduction maps to all ISAs, but the gather does
   not; AVX2 gather is the only lane-parallel load and is microarch-dependent.
8. Draw-neutrality: **BITWISE-UNSAFE (reduction)** AND, if dispatched, makes
   draws ISA-DEPENDENT (breaks machine-independent equivalence). Acceptable only
   as a fixed cross-ISA accumulator-layout kernel (e.g. always 4 logical partial
   sums, identical combine tree on scalar+NEON+SSE2+AVX2) + one-time re-record.
9. Effort: MED-HIGH (fixed-layout kernel across 4 ISAs) + re-record + a VD
   decision on ISA-dependent draws.
10. u8 interaction: NONE (doubles + size_t indices).
=> Highest hotness, but blocked by draw-neutrality and gather-bound economics.
   NOT recommended without an explicit VD decision.

### #C  Cut-scan histogram  (unresolved: [[scan.hpp:61-144@4a521760]], histogramDenseCutScan + scanOrdinalCuts)
1. The scatter-reduce histogram: gather code=column[indices[i]] then SCATTER a
   (count,sumW,sumWz) double triple into bins[code]; then a small prefix scan.
2. Frequency: **COLD in the default sampler.** Only caller is growTreeFromRoot
   ([[grow.hpp:118@4a521760]]), i.e. the OPT-IN grow-from-root WARM START
   (bartcore_growFromRoot; a fixed few sweeps at init), NOT the per-iteration MH
   kernel. Default init is prior-grown (growSubtreeFromPrior). Becomes
   per-iteration hot ONLY if informed-proposal birth/death lands
   (parallel-bart-frontier.md 3.1, not yet implemented).
3. Codegen: scalar scatter-reduce in the header; -O2 will not vectorize a
   data-dependent scatter under SSE2 baseline (essentially scalar today).
4. Access: gather (codes) + **scatter-reduce** into bins keyed by code (a
   histogram), then contiguous prefix scan.
5. Why CRAN can't: data-dependent scatter index (code) with read-modify-write
   into doubles; no portable SIMD scatter; and it's a double reduction.
6. Benefit: this is the explicitly PLANNED kernel (kernel-vocabulary.md item 5,
   `misc_scanColumn` -> out[3*numCodes]). SIMD histograms (partial-histogram
   replication, or the cumulative-sum-over-sorted-codes / XGBoost trick) can win
   substantially WHEN hot - but it is not hot today.
7. ISA portability: partial-histogram replication maps to NEON/SSE2/AVX2;
   AVX2/AVX512 conflict-detection helps the scatter. Non-trivial per ISA.
8. Draw-neutrality: **BITWISE-UNSAFE** (double scatter-reduce; comment scan.hpp:
   123-126 depends on exact left-fold order, "bitwise total - left"). Any
   reordering changes grow-from-root draws -> re-record.
9. Effort: HIGH.
10. u8 interaction: u8 codes halve the code-load bandwidth, but the scan is
    bottlenecked by the double-triple SCATTER, not the code load, so u8 barely
    helps here.
=> Defer until informed proposals make it per-iteration hot; then it is the
   prime target. Not a current win.

### #D  SSR for the sigma draw  ([[moments.c:2768@4a521760]]/2789, [[chain.hpp:908@4a521760]], [[model.hpp:1763@4a521760]])
Reduction sum (y-yhat)^2, scalar, undispatched. **Once per sweep** (not per
tree) -> ~1/numTrees the hotness of #1/#3. Contiguous (no gather), so it WOULD
vectorize cleanly - but it is a reduction: bitwise-unsafe + ISA-dependence, same
regime as #3, for a fraction of the payoff. Low priority.

---

## moments ISA-gap verdict

moments_sse2.c ships **SSE2 only** for the Mean / VarianceForKnownMean family
(no NEON, no AVX2); dispatched via misc_stat_setSIMDInstructionSet
([[moments.c:1330@4a521760]]). BUT these dispatched Mean/Variance reductions are **NOT in the
sampler draw inner loop** - bartcore's hot reductions are the SEPARATE scalar
undispatched SufficientStatisticsFast family (#3) and SSR (#D). The Mean/Variance
kernels serve standalone/mt/utility callers. So the SSE2-only gap is REAL but
LOW PRIORITY: completing NEON/AVX2 for Mean/Variance does nothing for the
sampler hot path, and on arm64 they already fall back to a fine scalar unrolled
loop. Do NOT prioritize completing the moment ISA family for speed; only revisit
if a standalone/utility consumer profiles hot. (Also note: adding NEON/AVX2
Mean/Variance would introduce cross-ISA reduction-order divergence in those
outputs - fine for standalone callers under the waiver, but keep them out of any
draw path.)

## u8 disposition call

u8 code storage (step-3) is a **pure MEMORY lever, orthogonal to any SIMD-compute
win**, and it currently carries a partition REGRESSION:
- No hot SIMD kernel makes u8 a throughput win. The residual roll / suff-stat /
  fit scatter all operate on doubles + size_t indices (u8-invariant). The scan
  is scatter-bound (u8 barely helps) and cold anyway. Predict/partition is u8's
  only consumer, and phase-1 (hot-layer-u8.md) found NO arm64 partition
  win for u8.
- REGRESSION as shipped: misc has NO u8 SIMD partition kernel, so u8 dense
  ordinal columns take the SCALAR partition ([[tree.hpp:713-722@4a521760]]), while u16 columns
  get the SIMD partition (partition_*.c). Step-4 (u8 SIMD partition family) was
  deferred. So landing step-3 u8 as-is trades u16's 8-lane SIMD partition for a
  scalar loop on eligible columns - the memory saving must outweigh that.
- The one genuine SIMD story tied to u8 is the deferred u8 PARTITION port: 16
  lanes per 128-bit (double u16's 8) with SSE2 native unsigned u8 min/max - it is
  (a) NEEDED to avoid the partition regression and (b) MIGHT win on x86
  (unmeasured; the arm64 phase-1 NEON prototype was a coarser block-skip, not a
  full port, so its no-win is not conclusive).
Call: u8 should land ONLY as the container's structural memory optimization
(exactly VD's 2026-07-07 decision: fold into data-ownership, no standalone
phase-2), and MUST be paired with a proper u8 SIMD partition port (or accept the
partition regression on eligible columns, justified by DRAM-bound memory gains).
No scan/predict SIMD kernel flips u8 to a net compute win. u8's case rests on
end-to-end memory in the n>=1e5 regime, to be measured in the container.

## TOP RECOMMENDATIONS

1. **Fused residual-roll AVX/AVX2 kernel (x86)** - the hottest per-tree double
   pass currently capped at SSE2 on x86; bitwise-safe & ISA-independent
   (elementwise), zero re-record. Gate on an x86 cache-resident microbench;
   modest expected win, arm64 already NEON.
2. **Close the x86-simd.md dispatch gaps** - wire misc_setIndexedVectorToConstant
   into the table and restore-or-delete the three scalar-bodied SSE2/AVX linalg
   arms (non-temporal stores may matter at n>=1e5). Bitwise-safe; low effort;
   maintainer-run x86 measurement decides.
3. **Hold the suff-stat / cut-scan SIMD** pending a VD decision: both are
   reductions that would either break machine-independent bitwise draws
   (needs fixed cross-ISA layout + re-record) or, for the scan, are cold until
   informed proposals land. Do not implement speculatively.

## QUESTION FOR VD

The only large hotness left in the draw path is the leaf-suffstat reduction
(setNodeAverages, #3), and SIMD-izing it necessarily makes MCMC draws
ISA-dependent (or requires a fixed cross-ISA accumulator layout + a one-time
equivalence re-record). Recommendation: DO NOT trade away machine-independent
bitwise reproducibility for a gather-bound reduction whose SIMD ceiling is low on
NEON/SSE2 anyway. Keep suffstat scalar; take the safe elementwise x86 gaps
(#1/#2) instead. Confirm this stance before anyone hand-SIMDs a draw-path reduction.
