# Memory-wall frontier: attacking BART's per-sweep bottleneck (ideation)

Status: ideation panel, 2026-07-20. Three INDEPENDENT lenses (algorithms,
hardware, systems/memory), each pushed to go beyond the recorded designs, plus a
measured SSR-SIMD NO-GO. This is NOT a proposal to build - it is a ranked idea
map with falsifiers, for VD to pick the first experiment. Companion to
docs/design/parallel-bart-frontier.md, within-chain-threading.md, block-fusion.md
(excised), data-layout.md, gpu-bart.md.

## 0. The profile correction (load-bearing; the recorded profiles are STALE)

Two independent agents, reading current bartcore HEAD, converged on the same
finding: the leafOf / muByTree constant-leaf-fits refactor (which post-dates
data-layout.md and within-chain-threading.md) ALREADY eliminated the old ~15%
fit SCATTER, and made the residual roll and the totalFits rebuild STREAMING
passes (resid[i] uses mu[leafOf[i]], an L1-resident table lookup;
chain.hpp:2595-2596,2797-2828). The surviving latency-bound hotspot is now a
SINGLE shape: the random suffstat GATHER treeY[indices[k]] over each tree's
shuffled per-tree index buffer (misc_computeIndexedSufficientStatisticsFast,
moments.c:331; tree.hpp:507; plus the two move-phase child-suffstat gathers).
The docs' 32%-gather / 15%-scatter four-way split no longer describes reality -
the gather now dominates. CHEAPEST DECISIVE FIRST ACTION: re-profile HEAD on the
x86 quiet box (dbarts-bench) to re-baseline before committing to anything.

Both dead efforts routed AROUND the gather's access pattern: block-fusion
attacked pass VOLUME (ceiling ~1.3-2x, excised), within-chain-threading attacked
memory-level parallelism across CORES (collapsed on cross-CCX L3 traffic,
~0.9-1.67x). Neither touched the random ACCESS PATTERN. That is the unworked
axis, and every strong idea below attacks it - convert the random gather into a
streaming reduce, fewer elements, or sublinear-in-n work.

Exactness note used throughout: any change that regroups the suffstat sum is the
"shifting" RNG class - ungrouped fits stay bitwise, grouped/all draws shift once,
a one-time equivalence + snapshot re-record (VD has accepted this class before).
A SCALAR fixed-order regroup preserves the within-host cross-ISA bitwise contract
(no SIMD-lane reassociation); a precision change or a SIMD-lane reduce does not.

## 1. The convergent flagship: kill the random gather

### 1a. Histogram-fused suffstat (systems lens; LOWEST RISK, rides existing leafOf)
Fuse the suffstat into the residual roll as one obs-order pass: for i in 0..n,
roll resid[i] and then acc[leafOf[i]] += resid[i] into a tiny L1-resident
per-leaf accumulator; the node suffstat is then acc[nodeid] (unweighted
sumWeights = node count, already known). leafOf is ALREADY maintained obs-order
and current for the pre-move partition, so setNodeAverages collapses into the
roll and treeY is touched ONCE. Converts n random dependent loads (~64.n B random
traffic) into a sequential read (~12.n B, prefetched) + an L1 scatter-add - the
latency->bandwidth conversion neither prior effort attempted, and unlike
data-layout.md's reorder it REMOVES the gather rather than relocating it
(neutral-to-win even at small n). Expected ~2-4x on the now-dominant pass at
n>=1e5. Degenerate single-dominant-leaf tree -> one FP-add dependency chain;
mitigate with accumulator PRIVATIZATION (K fixed copies round-robined by i%K,
fixed combine order - also single-thread MLP for the scatter side). Scalar
fixed-order, so cross-ISA bitwise preserved. One re-record. Risk: low-moderate
(only the constant-leaf steady-state path; move-phase child gathers stay).
Falsifier: microbench roll+histogram vs roll+gather at n in {1e4,1e5,1e6}, K in
{1,4,8}, plus a worst-case single-dominant-leaf tree.

### 1b. Fixed predictor-space order + run-length leaves (algorithms lens; higher ceiling)
Sort all n observations ONCE into a global predictor-space order (Z-order /
Hilbert over the codes; the engine already holds per-column sorted orders) and
permute treeY/y/weights/totalFits/leafOf into it permanently. Represent each leaf
not as a shuffled index list but as a RUN-LIST: a few contiguous intervals into
the fixed order (leaves are axis-aligned boxes; the space-filling order clusters
neighbors). The suffstat becomes per-run contiguous streaming reduces, no
indices[]. Distinct from data-layout.md's rejected per-tree reorder: there is NO
per-tree reorder - one fixed order shared by all m trees, residual maintained in
it by the ordinary streaming roll, so the cost that doc feared does not exist.
Data-derived order is RNG-independent -> scalar fixed-order reduce stays cross-ISA
bitwise. One re-record. Risk: real - node rep becomes run-lists; the two-pointer
partition must split runs (fragmentation grows with depth; needs periodic O(n)
re-canonicalization; snapshot/restore snapshots run-lists). Falsifier (cheap,
decisive, RUN FIRST): on real post-burn-in forests, impose a Z-order and count
maximal runs per leaf; win if mean runs/leaf is small (< ~50 at n=1e5), dead if
fragmentation approaches n_leaf/2.

### 1c. Segment-tree / Fenwick over the fixed order (algorithms lens; SUBLINEAR, highest risk)
The aggressive sibling of 1b on the same fixed-order substrate. Hold the residual
in a lazy range-add / range-sum segment tree over the global order. A leaf redraw
is R range-adds; a leaf suffstat is R range-sum queries - O(#leaves.R.log n) per
sweep instead of O(m.n). The ONLY idea whose advantage GROWS with n (polylog vs
linear), so it targets the n=1e6+ hard case specifically. Same run-count census
gate as 1b + a standalone lazy-Fenwick throughput microbench. Highest constant
factors and implementation complexity; only escalate to it from 1b if the large-n
asymptotics justify the segment-tree constants.

## 2. Exactness-preserving accelerators (compose on top of section 1)

### 2a. Control-variate minibatch backfitting (SVRG-anchored)
Keep an exact per-leaf suffstat anchor from the last full sweep; on most sweeps
read only a fraction q of the residual and estimate the DELTA since the anchor
(low variance because the per-tree signal is O(1/sqrt m)). BART is unusually
subsample-friendly: shallow trees keep each leaf's members in the tens of
thousands even at large n, so subsample variance is smallest exactly where the
memory wall bites. ~1/q fewer gathered elements. Approximate as stated
(SBC-gated bias budget; anchor refresh every K bounds drift; restrict to the leaf
draw, keep move-scoring exact) - or made EXACT by wrapping in delayed acceptance
(2b). Falsifier: SBC coverage + posterior RMSE at matched wall-clock, q in
{0.05,0.1,0.25}, n in {1e5,1e6}; kill if move-decision bias distorts varcounts.

### 2b. Delayed acceptance with a compressed stage-1 surrogate
Twist on parallel-bart-frontier.md 3.2a (already falsifier-cleared: 0.97
frozen-vs-true, ~7-12% survivors). Score all moves' stage 1 against a COMPRESSED
frozen residual (subsample, block-average, count-sketch, or fp32); stage 2 pays
the exact full-resolution gather only on the ~10% survivors. Christen-Fox keeps
the posterior EXACT regardless of surrogate crudeness - this is what makes
subsampling/sketching exact rather than approximate. Smaller lever alone (the
move scan is a smaller band than the leaf-draw gather) but composes with 1/2 and
independently raises ESS/sweep. Falsifier: survivor-rate inflation vs surrogate
resolution; kill if a crude surrogate defeats the filter.

### 2c. Mixed-precision residual (minor; screening-only)
fp32/bf16 residual halves/quarters the BANDWIDTH-bound fraction (helps only at
n>=1e6, LLC-spilling; latency-bound moderate-n is untouched). Precision change ->
perturbs the posterior AND threatens cross-ISA bitwise, so use it ONLY as 2b's
stage-1 screen (exact fp64 stage-2 rescan), never on the exact path.

## 3. Cheap systems increments (mostly bit-identical, shippable independently)
- uint8/uint16 compacted leafOf (systems 2a): shallow trees rarely exceed 255
  leaves; 4x bandwidth cut on the leaf streams + 4x smaller n.m footprint (more
  LLC-resident at large n). Compose with 1a. Falsifier: live #leaves distribution.
- Software prefetch the surviving shuffled gathers (moments.c:344): __builtin_
  prefetch(&x[indices[i+D]]) D ahead. Bit-identical, no re-record. But the kernel
  already runs 5-wide MLP, so it only buys the latency the 5 in-flight loads miss
  (~1.15-1.4x if the LFBs are not already saturated). Falsifier: sweep D.
- Wider gather MLP (moments.c:338): 8-10-wide unroll + 2-4 accumulator banks to
  push outstanding misses toward the ~10-12 LFB limit; single-core MLP via ILP
  (dodges the cross-CCX collapse that no-go'd threading). One microbench falsifies
  both this and prefetch (if 10-wide ~ 5-wide, the LFBs are saturated).
- Huge pages / TLB on the big buffers (no doc mentions TLB): treeY at n=1e6 is
  2048 pages >> dTLB; 2 MB pages remove the page-walk from each gather miss.
  Bit-identical. Validate on the x86 Linux box (macOS masks it).
- Per-pass NT/temporal cache hints: mostly moot post-leafOf (the big fit-slab
  write is gone); a small once/sweep totalFits-rebuild NT store is the only
  survivor. Polish.

## 4. Hardware (the device answer, re-grounded)

The reframing that reorders the hardware question: for the flagship CONSTANT leaf,
the leaf-assignment matrix A is ONE-HOT, so suffstats = A^T r is a SEGMENTED SUM
(arithmetic intensity ~0), memory-bound, NOT a dense inner product. That kills the
tensor-core / NPU / AMX framing for the constant leaf on sight - those units
accelerate arithmetic intensity this op does not have. The lever is bandwidth +
MLP (the DRAM/HBM system), never the systolic array.

### 4a. Unified-memory APU (the real device answer; NEW premise)
gpu-bart.md 2f killed the CPU/GPU split on PCIe transfer granularity - a DISCRETE-
GPU assumption. On a unified-memory APU (Apple Silicon Metal/MPS, Grace-Hopper,
MI300A) the CPU and GPU share one physical pool with ZERO copy, so the residual/
fits/indices are visible to the host (branchy moves) and device (O(n) passes) with
no upload - the objection that closed 2f does not exist. Host keeps the pointer-
tree MH moves; device runs the O(n) passes at 5-50x the bandwidth a CPU core
sustains and 100-1000x the MLP (exactly what the latency-bound gather needs, which
8 cores structurally could not give). Uniquely PRESERVES dbarts's mutable-data
embedded-Gibbs contract (the reason dbartsSampler/dbarts.h exist) that a
from-scratch dense-heap GPU engine (bartz) abandons. CAVEAT: the hardware lens
built its launch-count math on the block-fusion atom map, which was EXCISED - the
physics (zero-copy, HBM bandwidth, MLP) stands, but the per-sweep launch structure
must be re-grounded to per-node or a fresh reduction blocking (compose with 1a's
histogram as the device kernel). Packaging: Metal/MPS are OS frameworks (no CUDA),
configure-probed optional backend behind the simd.c dispatch pattern; absent ->
CPU path unchanged; device determinism contract, never the CPU bitwise claim.
Falsifier: microbench MPS launch + strided-gather unified bandwidth on an M-series
part; if launches + readback exceed the DRAM time saved, same wall, higher ceiling.

### 4b. AMX / tensor for the ONE compute-bound corner: linear + GP leaves (NEW)
Linear leaves (U^T W U, p x p GEMM) and GP leaves (Cholesky up to gpMaxLeafSize
= 256, O(256^3)) are genuinely COMPUTE-bound - the batched-small-GEMM / batched-
Cholesky pattern Intel AMX, Apple AMX (Accelerate), and cuSOLVER-batched serve.
Best CRAN story of all: Intel/Apple AMX are CPU features reachable via the same
runtime dispatch as simd.c (no GPU, graceful fallback to scalar leaf math), and it
speeds up stan4bart-style and GP-leaf models where the O(n_leaf^3) leaf math
dominates. Risk: AMX's fast path is bf16/fp16 vs the fp64 draw/reproducibility
contract (use fp64/fp32 tiles, narrower win); pays off only if leaves near the 256
cap are common (shallow BART makes them rare). Falsifier: profile a GP run - what
fraction is Cholesky+Gram vs field maintenance.

### 4c. Reasoned NO's (each a substantive result)
- NPU / ANE / Hexagon: the WORST fit. NPUs want dense, low-precision, static-
  dataflow GEMM; the constant-leaf path is sparse fp64 + data-dependent segmented
  reduce + branchy control flow + an fp64 reproducibility contract - four-for-four
  against everything an NPU is good at. Structural no, not a tooling-maturity no.
- Tensor cores / AMX for the CONSTANT leaf: no (one-hot A -> intensity ~0). Only
  the p-wide dense-block linear/GP leaves of 4b qualify.
- PIM / near-memory (UPMEM, HBM-PIM): attacks the gather latency where data lives,
  but a leaf's members are shuffled across all banks -> all-to-all defeats bank
  locality unless observations are physically leaf/atom-sorted (composes with 1b);
  and no CRAN path (HPC-cluster hardware). Research/cluster only, n >> 1e6.
- CXL memory: a load/store fabric that ADDS ~100-300 ns/access - helps
  capacity-bound, not latency-bound; relevant only at n >> 1e8 dbarts does not
  target.
- Chains-as-batch on one device: off-target for the confirmed-common single-chain
  n>=1e5 flagship, and fights SIMT divergence.

## 5. SSR SIMD leaf: MEASURED NO-GO (2026-07-20)
Prototyped the fixed-lane SIMD SSR reduction and measured it (M1 Max, NEON 8-lane
vs current scalar, min-of-9): kernel 1.62x at n=1e3 (L1) -> 1.14x (L2) -> 1.05x
(n=1e6) -> 1.00x (n=1e7); bandwidth-bound not latency-bound (16-lane == 8-lane),
because the current scalar is already unrolled-by-5 and extracts the ILP. Fixed
8-lane / FMA-off layout is scalar==NEON bitwise (preserves cross-ISA) and differs
from current scalar by ~14 ULP (the re-record). SSR runs ONCE/sweep vs 75-200
trees' O(n) passes -> ~0.12-0.15% of a sweep -> implied whole-sweep speedup
<= 0.05% best case, ~0% at scale, under the +/-2% bench noise floor. Against a
one-time re-record + a permanent obligation to keep a fixed cross-ISA lane layout
+ no-FMA contract in sync across four ISAs (the first coupled draw-path
reduction). Clearly negative - keep SSR scalar.

## 6. Recommended experiment order (cheapest-decisive first)
1. RE-PROFILE HEAD on dbarts-bench (x86, quiet): confirm the scatter-gone /
   gather-dominant picture and re-baseline the % split. Everything else is scoped
   against this. (hours)
2. RUN-COUNT CENSUS under a Z-order on real post-burn-in forests: mean maximal
   runs per leaf. Gates 1b AND 1c together; the single highest-value, entirely-
   unrun experiment. (hours)
3. PROTOTYPE + microbench the histogram-fused suffstat (1a) - the lowest-risk
   flagship, rides existing leafOf, no new substrate. If it clears, it is the
   first thing to build.
4. SUBSAMPLE SBC/RMSE sweep (gates 2a and 2b's compressed stage-1). (day)

All of section 1-3 is single-core / draw-preserving-modulo-one-re-record and
CRAN-shippable. Section 4 is the device horizon (4a the real answer, 4b the
shippable CPU-AMX corner). None is committed; this doc is the record so the next
agent does not re-derive that (a) the profiles moved, (b) the gather's access
pattern is the unworked axis, and (c) the tensor/NPU framing is a mirage for the
constant leaf.

## 7. Post-panel corrections and measured results (2026-07-20; authoritative)

- BLOCK-FUSION RE-TREAD (VD's catch; corrects section 1). The "histogram-fused
  suffstat" flagship (1a) is NOT new: it is block-fusion Stage A, which fused the
  three per-tree passes at b=1 and measured BENCH-NEUTRAL (every run-* metric
  within 4% on the x86 box) before the whole project was excised WONT-DO. The
  "fixed-order run-list" (1b) re-derives the atom map, and its "run-count census"
  IS block-fusion's already-run E1 atom census (atoms far below n). So 1a and 1b
  are a re-derivation of an already-measured, bench-neutral, excised approach; the
  prior is discouraging. Caveat: Stage A predates the leafOf/muByTree refactor, so
  it is not a strict re-run - but the burden is on a re-profile to show the
  shifted profile changes the verdict. Do NOT re-open 1a/1b without it. The lesson
  the panel missed: atoms == binning, and block-fusion already spent it.

- METAL / UNIFIED-MEMORY GPU: MEASURED DEAD at BART sizes (corrects 4a's
  optimism). Probe on M1 Max (throwaway, discarded): (i) fp64 is UNSUPPORTED on
  the Apple GPU - `double` is a hard MSL compile error, not emulated; (ii)
  per-dispatch launch+readback floor ~140 us; (iii) the isolated gather-reduce is
  0.73x CPU at n=1e5 (LOSES - the 140 us launch dominates the 44 us of compute),
  4-8x at n>=1e6. Two killers at the common n=1e5: fp64-unsupported (which VD's
  fp32 tolerance dissolves), and - the structural one precision cannot fix -
  backfitting is SEQUENTIAL (tree t+1 needs tree t's fp64 draw CPU-side), so the
  140 us round trip is paid PER TREE PER SWEEP, un-batchable, and 140 us > the
  entire CPU gather (148 us) at n=1e5. GPU offload wins only at n>=1e6 AND only
  with fp32 suffstats. The ONLY escape is a restructure: batch the whole sweep's
  residual-independent stage-1 scoring into ONE launch (delayed acceptance, 2b),
  fp32 on GPU kept exact by the fp64 CPU survivor rescan, amortizing the 140 us
  over the sweep - speculative, heavy. Verdict: keep the reduce on the CPU; the
  GPU is an n>=1e6 / restructured-sweep horizon, not a near-term lever.

- PRECISION RE-ELEVATED to a TOP lever (VD 2026-07-20: fp64 is R's default, not a
  design commitment; reduced-precision STORAGE is acceptable for large gains).
  This is the one genuinely-NEW memory-wall attack, orthogonal to block-fusion's
  binning: store the residual/fits in fp32 (fp64-accumulated reductions) - halves
  the bytes on the bandwidth-bound streaming/gather passes and ~doubles the
  LLC-resident n (DRAM-bound -> cache-bound at the common n>=1e5), plausibly
  ~1.5-2x at large n, single-core, CRAN-shippable, one re-record, cross-ISA
  bitwise preserved (IEEE-deterministic). The gate is STATISTICAL not
  architectural: SBC / interval coverage vs the fp64 sampler, because MCMC
  stationarity can amplify a per-sweep bias a single-shot fp32 loss would not
  show. fp32-with-fp64-accumulation is the safe target (~1e-7, below the sampler
  noise floor); bf16 (~3 digits) is the aggressive version. This SUPERSEDES the
  histogram flagship (section 1) as the recommended first build. NEXT EXPERIMENT:
  confirm the fp32-vs-fp64 bandwidth/cache-crossover upside (microbench), then
  prototype the fp32 residual store and run it through SBC/coverage.

## 8. fp32 microbench: MEASURED (2026-07-20; refines section 7's prediction)

Ran a throwaway two-pass memory microbench (streaming reduce + shuffled-index
gather-reduce, fp64 store vs fp32 store with an fp64 accumulator, min-of-7,
~3e8 elements/measurement, n from 1e4 to 3e7) on BOTH the M1 Max (arm64,
~48 MB SLC, huge single-core BW) and the x86 bench box (Ryzen 7 3700X, 16 MB
L3/CCX, DDR4). Both index widths measured: size_t (8 B, the current engine) and
uint32 (4 B, a narrowing that is a SEPARATE lever). Raw numbers preserved off-
tree. The prediction ("top lever, ~1.5-2x at large n, ship it") is REFINED on
three points:

- THE WIN IS GATHER-ONLY AND CACHE-CROSSOVER-SHAPED, not asymptotic. fp32 speeds
  the shuffled gather by ~1.7x (M1) / ~2.0x (x86) in the band where fp32 stays
  LLC-resident but fp64 spills (n ~= 1e6 on both), then TAPERS to ~1.1x once both
  are deep in DRAM (n >= 1e7). It does not keep widening with n: the size_t index
  (8 B, common to both stores) caps the gather ratio near 16/12 = 1.33x, so at
  large n the win shrinks, not grows. Narrowing the index to uint32 lifts the
  asymptote toward 1.5x and is worth ~a modest extra gather gain (measured), but
  it is orthogonal engine work, not part of the fp32 store.

- THE STREAMING PASS IS MACHINE-DEPENDENT AND CAN LOSE. On the M1 the sequential
  reduce is COMPUTE-bound on the fp64-add accumulator chain (~0.19 ns/elem flat
  from L1 to 240 MB, ~41 GB/s single-core well under the BW ceiling) - so fp32
  buys ZERO on streaming there. On the x86 the streaming reduce IS bandwidth-bound
  at large n (fp64 GB/s collapses 57 -> 17.6 as it spills L3), and fp32 wins
  ~1.3-1.4x - BUT cache-resident (n <= 1e6) the x86 fp32 reduce is ~2x SLOWER than
  fp64, because the per-element fp32->fp64 conversion (cvtss2sd) is the bottleneck
  when the bytes are already free. So the fp32->fp64 conversion is NOT free, and
  an UNCONDITIONAL fp32 store REGRESSES the common small/medium-n cache-resident
  case. The store would have to be n-gated (opt-in or an auto-switch above a
  footprint threshold), which is real design complexity the section-7 note did
  not price in.

- NET: fp32 storage is a large-n (n >~ 1e6), gather-dominated, cache-residency
  lever worth ~1.3-2x ON THE BOUND PASSES IN THAT BAND, tapering at the very
  largest n and inverting to a loss for small n. Whether it yields a worthwhile
  END-TO-END sweep speedup now turns entirely on the gather's SHARE of a sweep
  (section 6 experiment 1, the re-profile) - the number that gates the engine
  prototype. [Share pending; verdict to follow.]
