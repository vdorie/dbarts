# Memory-wall frontier: attacking BART's per-sweep bottleneck (ideation)

Status: CLOSED as an idea map; its recommended lever LANDED. Ideation panel run
2026-07-20 (three lenses: algorithms, hardware, systems/memory), followed same
day by direct measurement of the strongest candidates. This was never a proposal
to build by itself - it is a ranked idea map with falsifiers, and the record of
which branches were tried, which closed, and which one shipped. Companion to
docs/design/parallel-bart-frontier.md, within-chain-threading.md, block-fusion.md
(excised), data-layout.md, gpu-bart.md, and docs/plans/blocked-jacobi-trees.md
(within-chain-parallelism verdicts).

Summary: the panel's own re-profile found that the leafOf/muByTree refactor had
already collapsed the old four-way 32%-gather/15%-scatter split into one
dominant hotspot - the random suffstat gather over each tree's shuffled index
buffer - and every strong idea below targets that access pattern. Of those
ideas, the "kill the gather" flagship (histogram-fused suffstat / fixed-order
run-lists) turned out to be a re-derivation of block-fusion, already measured
bench-neutral and excised; it was not reopened. A unified-memory GPU/APU was
measured dead at BART's common sizes (Apple GPUs do not support fp64, and
backfitting's per-tree sequentiality makes the ~140 us launch/readback floor
un-batchable at n=1e5). The lever that actually paid off was reduced-precision
storage - re-elevated mid-panel from a "minor, screening-only" idea to the
recommended first build, then measured, forked, and shipped in
reduced-precision-storage.md: a bitwise-preserving uint32 gather index (landed,
default-on) and an opt-in fp32 residual store (landed, ~1.1-1.3x depending on
thread/chain count). Post-landing re-profiling confirmed no other ruled-out
lever revives - within-chain threading and blocked-jacobi-trees remain NO-GO on
their own, independently measured grounds (see those docs); fp32 storage is the
memory-wall answer for this engine.

## 1. The bottleneck: profile correction

The profiles on record going into this panel were stale. Two independent agents,
reading current bartcore HEAD, converged on the same finding: the leafOf /
muByTree constant-leaf-fits refactor (which post-dates data-layout.md and
within-chain-threading.md) had ALREADY eliminated the old ~15% fit SCATTER, and
turned the residual roll and the totalFits rebuild into STREAMING passes
(resid[i] uses mu[leafOf[i]], an L1-resident table lookup; chain.hpp:2595-2596,
2797-2828). What survives as the latency-bound hotspot is now a SINGLE shape:
the random suffstat GATHER treeY[indices[k]] over each tree's shuffled per-tree
index buffer (misc_computeIndexedSufficientStatisticsFast, moments.c:331;
tree.hpp:507; plus the two move-phase child-suffstat gathers). The docs' old
32%-gather / 15%-scatter four-way split no longer describes reality - the
gather now dominates. This made the cheapest decisive first action obvious:
re-profile HEAD on the x86 quiet box (dbarts-bench) to re-baseline before
committing to anything (done later; see section 8b).

Both prior dead efforts routed AROUND the gather's access pattern rather than
through it: block-fusion attacked pass VOLUME (ceiling ~1.3-2x, excised),
within-chain-threading attacked memory-level parallelism across CORES (collapsed
on cross-CCX L3 traffic, ~0.9-1.67x). Neither touched the random ACCESS
PATTERN. That is the unworked axis, and every strong idea in this document
attacks it - convert the random gather into a streaming reduce, fewer elements,
or sublinear-in-n work.

Exactness note, used throughout: any change that regroups the suffstat sum is
the "shifting" RNG class - ungrouped fits stay bitwise, grouped/all draws shift
once, a one-time equivalence + snapshot re-record (VD has accepted this class
before). A SCALAR fixed-order regroup preserves the within-host cross-ISA
bitwise contract (no SIMD-lane reassociation); a precision change or a
SIMD-lane reduce does not.

## 2. Killing the gather: candidate designs

### 2a. Histogram-fused suffstat (systems lens; lowest risk, rides existing leafOf)

Fuse the suffstat into the residual roll as one obs-order pass: for i in 0..n,
roll resid[i] and then acc[leafOf[i]] += resid[i] into a tiny L1-resident
per-leaf accumulator; the node suffstat is then acc[nodeid] (unweighted
sumWeights = node count, already known). leafOf is ALREADY maintained obs-order
and current for the pre-move partition, so setNodeAverages collapses into the
roll and treeY is touched ONCE. Converts n random dependent loads (~64.n B
random traffic) into a sequential read (~12.n B, prefetched) + an L1
scatter-add - the latency->bandwidth conversion neither prior effort attempted,
and unlike data-layout.md's reorder it REMOVES the gather rather than
relocating it (neutral-to-win even at small n). Expected ~2-4x on the
now-dominant pass at n>=1e5. Degenerate single-dominant-leaf tree -> one FP-add
dependency chain; mitigate with accumulator PRIVATIZATION (K fixed copies
round-robined by i%K, fixed combine order - also single-thread MLP for the
scatter side). Scalar fixed-order, so cross-ISA bitwise preserved. One
re-record. Risk: low-moderate (only the constant-leaf steady-state path; the
move-phase child gathers stay). Falsifier: microbench roll+histogram vs
roll+gather at n in {1e4,1e5,1e6}, K in {1,4,8}, plus a worst-case
single-dominant-leaf tree.

### 2b. Fixed predictor-space order + run-length leaves (algorithms lens; higher ceiling)

Sort all n observations ONCE into a global predictor-space order (Z-order /
Hilbert over the codes; the engine already holds per-column sorted orders) and
permute treeY/y/weights/totalFits/leafOf into it permanently. Represent each
leaf not as a shuffled index list but as a RUN-LIST: a few contiguous intervals
into the fixed order (leaves are axis-aligned boxes; the space-filling order
clusters neighbors). The suffstat becomes per-run contiguous streaming reduces,
no indices[]. Distinct from data-layout.md's rejected per-tree reorder: there is
NO per-tree reorder - one fixed order shared by all m trees, residual
maintained in it by the ordinary streaming roll, so the cost that doc feared
does not exist. Data-derived order is RNG-independent -> scalar fixed-order
reduce stays cross-ISA bitwise. One re-record. Risk: real - node rep becomes
run-lists; the two-pointer partition must split runs (fragmentation grows with
depth; needs periodic O(n) re-canonicalization; snapshot/restore snapshots
run-lists). Falsifier (cheap, decisive, RUN FIRST): on real post-burn-in
forests, impose a Z-order and count maximal runs per leaf; win if mean
runs/leaf is small (< ~50 at n=1e5), dead if fragmentation approaches
n_leaf/2.

Neither of these is new. Block-fusion Stage A already fused the three per-tree
passes at b=1 - exactly the "histogram-fused suffstat" of 2a - and measured
BENCH-NEUTRAL (every run-* metric within 4% on the x86 box) before the whole
project was excised WONT-DO. The "fixed-order run-list" of 2b re-derives
block-fusion's atom map, and its "run-count census" IS block-fusion's
already-run E1 atom census (atoms far below n). So 2a and 2b are a
re-derivation of an already-measured, bench-neutral, excised approach, and the
prior is discouraging. Caveat: Stage A predates the leafOf/muByTree refactor,
so it is not a strict re-run - but the burden is on a re-profile to show the
shifted profile changes the verdict. Do NOT re-open 2a/2b without one. The
lesson the panel missed the first time around: atoms == binning, and
block-fusion already spent it.

### 2c. Segment-tree / Fenwick over the fixed order (algorithms lens; sublinear, highest risk)

The aggressive sibling of 2b on the same fixed-order substrate. Hold the
residual in a lazy range-add / range-sum segment tree over the global order. A
leaf redraw is R range-adds; a leaf suffstat is R range-sum queries -
O(#leaves.R.log n) per sweep instead of O(m.n). The ONLY idea whose advantage
GROWS with n (polylog vs linear), so it targets the n=1e6+ hard case
specifically. Gated by the same run-count census as 2b (2b and 2c share the
fragmentation gate) plus a standalone lazy-Fenwick throughput microbench.
Highest constant factors and implementation complexity; only escalate to it
from 2b if the large-n asymptotics justify the segment-tree constants. Given
the discouraging prior on 2b's underlying atom map (above), 2c inherits the
same caution.

## 3. Exactness-preserving accelerators (compose on top of section 2)

### 3a. Control-variate minibatch backfitting (SVRG-anchored)

Keep an exact per-leaf suffstat anchor from the last full sweep; on most sweeps
read only a fraction q of the residual and estimate the DELTA since the anchor
(low variance because the per-tree signal is O(1/sqrt m)). BART is unusually
subsample-friendly: shallow trees keep each leaf's members in the tens of
thousands even at large n, so subsample variance is smallest exactly where the
memory wall bites. ~1/q fewer gathered elements. Approximate as stated
(SBC-gated bias budget; anchor refresh every K bounds drift; restrict to the
leaf draw, keep move-scoring exact) - or made EXACT by wrapping in delayed
acceptance (3b). Falsifier: SBC coverage + posterior RMSE at matched
wall-clock, q in {0.05,0.1,0.25}, n in {1e5,1e6}; kill if move-decision bias
distorts varcounts.

### 3b. Delayed acceptance with a compressed stage-1 surrogate

Twist on the batched cross-tree stale-residual scoring of
parallel-bart-frontier.md 3.2a (already falsifier-cleared: 0.97
frozen-vs-true, ~7-12% survivors). Score all moves' stage 1 against a
COMPRESSED frozen residual (subsample, block-average, count-sketch, or fp32);
stage 2 pays the exact full-resolution gather only on the ~10% survivors.
Christen-Fox keeps the posterior EXACT regardless of surrogate crudeness - this
is what makes subsampling/sketching exact rather than approximate. Smaller
lever alone (the move scan is a smaller band than the leaf-draw gather) but
composes with 2a/2b and independently raises ESS/sweep. Falsifier:
survivor-rate inflation vs surrogate resolution; kill if a crude surrogate
defeats the filter.

### 3c. Mixed-precision residual (minor; screening-only)

fp32/bf16 residual halves/quarters the BANDWIDTH-bound fraction (helps only at
n>=1e6, LLC-spilling; latency-bound moderate-n is untouched). Precision change
-> perturbs the posterior AND threatens cross-ISA bitwise, so use it ONLY as
3b's stage-1 screen (exact fp64 stage-2 rescan), never on the exact path - or
so it seemed at panel time. Section 8 below is the story of how this
"minor, screening-only" idea got re-elevated to the lever that actually shipped.

## 4. Cheap systems increments (mostly bit-identical, shippable independently)

- uint8/uint16 compacted leafOf (systems lens): shallow trees rarely exceed 255
  leaves; 4x bandwidth cut on the leaf streams + 4x smaller n.m footprint (more
  LLC-resident at large n). Compose with 2a. Falsifier: live #leaves
  distribution. (The uint16 variant was later evaluated as a Track-1 follow-on
  in reduced-precision-storage.md and DECLINED by measurement - section 6 of
  that doc.)
- Software prefetch the surviving shuffled gathers (moments.c:344):
  __builtin_prefetch(&x[indices[i+D]]) D ahead. Bit-identical, no re-record.
  But the kernel already runs 5-wide MLP, so it only buys the latency the 5
  in-flight loads miss (~1.15-1.4x if the LFBs are not already saturated).
  Falsifier: sweep D.
- Wider gather MLP (moments.c:338): 8-10-wide unroll + 2-4 accumulator banks to
  push outstanding misses toward the ~10-12 LFB limit; single-core MLP via ILP
  (dodges the cross-CCX collapse that no-go'd threading). One microbench
  falsifies both this and prefetch (if 10-wide ~ 5-wide, the LFBs are
  saturated).
- Huge pages / TLB on the big buffers (no doc mentions TLB): treeY at n=1e6 is
  2048 pages >> dTLB; 2 MB pages remove the page-walk from each gather miss.
  Bit-identical. Validate on the x86 Linux box (macOS masks it).
- Per-pass NT/temporal cache hints: mostly moot post-leafOf (the big fit-slab
  write is gone); a small once/sweep totalFits-rebuild NT store is the only
  survivor. Polish.

## 5. Hardware (the device answer, re-grounded)

The reframing that reorders the hardware question: for the flagship CONSTANT
leaf, the leaf-assignment matrix A is ONE-HOT, so suffstats = A^T r is a
SEGMENTED SUM (arithmetic intensity ~0), memory-bound, NOT a dense inner
product. That kills the tensor-core / NPU / AMX framing for the constant leaf
on sight - those units accelerate arithmetic intensity this op does not have.
The lever is bandwidth + MLP (the DRAM/HBM system), never the systolic array.

### 5a. Unified-memory APU (the real device answer, on paper - then measured dead)

gpu-bart.md 2f killed the CPU/GPU split on PCIe transfer granularity - a
DISCRETE-GPU assumption. On a unified-memory APU (Apple Silicon Metal/MPS,
Grace-Hopper, MI300A) the CPU and GPU share one physical pool with ZERO copy,
so the residual/fits/indices are visible to the host (branchy moves) and
device (O(n) passes) with no upload - the objection that closed 2f does not
exist. Host keeps the pointer-tree MH moves; device runs the O(n) passes at
5-50x the bandwidth a CPU core sustains and 100-1000x the MLP (exactly what the
latency-bound gather needs, which 8 cores structurally could not give).
Uniquely PRESERVES dbarts's mutable-data embedded-Gibbs contract (the reason
dbartsSampler/dbarts.h exist) that a from-scratch dense-heap GPU engine (bartz)
abandons. The hardware lens built its launch-count math on the block-fusion
atom map, which was EXCISED - the physics (zero-copy, HBM bandwidth, MLP)
stood, but the per-sweep launch structure needed re-grounding to per-node or a
fresh reduction blocking (compose with 2a's histogram as the device kernel).
Packaging premise: Metal/MPS are OS frameworks (no CUDA), configure-probed
optional backend behind the simd.c dispatch pattern; absent -> CPU path
unchanged; device determinism contract, never the CPU bitwise claim.

This was measured, and the optimism did not survive contact. Probe on M1 Max
(throwaway, discarded): (i) fp64 is UNSUPPORTED on the Apple GPU - `double` is
a hard MSL compile error, not emulated; (ii) per-dispatch launch+readback floor
~140 us; (iii) the isolated gather-reduce is 0.73x CPU at n=1e5 (LOSES - the
140 us launch dominates the 44 us of compute), 4-8x at n>=1e6. Two killers at
the common n=1e5: fp64-unsupported (which VD's fp32 tolerance dissolves), and
- the structural one precision cannot fix - backfitting is SEQUENTIAL (tree t+1
needs tree t's fp64 draw CPU-side), so the 140 us round trip is paid PER TREE
PER SWEEP, un-batchable, and 140 us > the entire CPU gather (148 us) at n=1e5.
GPU offload wins only at n>=1e6 AND only with fp32 suffstats. The ONLY escape
is a restructure: batch the whole sweep's residual-independent stage-1 scoring
into ONE launch (delayed acceptance, 3b), fp32 on GPU kept exact by the fp64
CPU survivor rescan, amortizing the 140 us over the sweep - speculative,
heavy. Verdict: keep the reduce on the CPU; the GPU is an n>=1e6 /
restructured-sweep horizon, not a near-term lever.

### 5b. AMX / tensor for the ONE compute-bound corner: linear + GP leaves

Linear leaves (U^T W U, p x p GEMM) and GP leaves (Cholesky up to
gpMaxLeafSize = 256, O(256^3)) are genuinely COMPUTE-bound - the
batched-small-GEMM / batched-Cholesky pattern Intel AMX, Apple AMX
(Accelerate), and cuSOLVER-batched serve. Best CRAN story of all: Intel/Apple
AMX are CPU features reachable via the same runtime dispatch as simd.c (no
GPU, graceful fallback to scalar leaf math), and it speeds up stan4bart-style
and GP-leaf models where the O(n_leaf^3) leaf math dominates. Risk: AMX's fast
path is bf16/fp16 vs the fp64 draw/reproducibility contract (use fp64/fp32
tiles, narrower win), and it pays off only if leaves near the 256 cap are
common (shallow BART makes them rare). Falsifier: profile a GP run - what
fraction is Cholesky+Gram vs field maintenance.

### 5c. Reasoned NO's (each a substantive result)

- NPU / ANE / Hexagon: the WORST fit. NPUs want dense, low-precision,
  static-dataflow GEMM; the constant-leaf path is sparse fp64 + data-dependent
  segmented reduce + branchy control flow + an fp64 reproducibility contract -
  four-for-four against everything an NPU is good at. Structural no, not a
  tooling-maturity no.
- Tensor cores / AMX for the CONSTANT leaf: no (one-hot A -> intensity ~0).
  Only the p-wide dense-block linear/GP leaves of 5b qualify.
- PIM / near-memory (UPMEM, HBM-PIM): attacks the gather latency where data
  lives, but a leaf's members are shuffled across all banks -> all-to-all
  defeats bank locality unless observations are physically leaf/atom-sorted
  (composes with 2b); and no CRAN path (HPC-cluster hardware). Research/cluster
  only, n >> 1e6.
- CXL memory: a load/store fabric that ADDS ~100-300 ns/access - helps
  capacity-bound, not latency-bound; relevant only at n >> 1e8 dbarts does not
  target.
- Chains-as-batch on one device: off-target for the confirmed-common
  single-chain n>=1e5 flagship, and fights SIMT divergence.

## 6. SSR SIMD leaf: measured NO-GO (2026-07-20)

Prototyped the fixed-lane SIMD SSR reduction and measured it (M1 Max, NEON
8-lane vs current scalar, min-of-9): kernel 1.62x at n=1e3 (L1) -> 1.14x (L2)
-> 1.05x (n=1e6) -> 1.00x (n=1e7); bandwidth-bound not latency-bound (16-lane
== 8-lane), because the current scalar is already unrolled-by-5 and extracts
the ILP. Fixed 8-lane / FMA-off layout is scalar==NEON bitwise (preserves
cross-ISA) and differs from current scalar by ~14 ULP (the re-record). SSR
runs ONCE/sweep vs 75-200 trees' O(n) passes -> ~0.12-0.15% of a sweep ->
implied whole-sweep speedup <= 0.05% best case, ~0% at scale, under the
+/-2% bench noise floor. Against a one-time re-record + a permanent obligation
to keep a fixed cross-ISA lane layout + no-FMA contract in sync across four
ISAs (the first coupled draw-path reduction). Clearly negative - keep SSR
scalar.

## 7. Recommended experiment order (as set at panel time; cheapest-decisive first)

1. RE-PROFILE HEAD on dbarts-bench (x86, quiet): confirm the scatter-gone /
   gather-dominant picture and re-baseline the % split. Everything else is
   scoped against this. (hours)
2. RUN-COUNT CENSUS under a Z-order on real post-burn-in forests: mean maximal
   runs per leaf. Gates 2b and 2c together; the single highest-value,
   entirely-unrun experiment. (hours)
3. PROTOTYPE + microbench the histogram-fused suffstat (2a) - the lowest-risk
   flagship, rides existing leafOf, no new substrate. If it clears, it is the
   first thing to build.
4. SUBSAMPLE SBC/RMSE sweep (gates 3a and 3b's compressed stage-1). (day)

All of sections 2-4 is single-core / draw-preserving-modulo-one-re-record and
CRAN-shippable. Section 5 is the device horizon (5a the real answer, 5b the
shippable CPU-AMX corner). None of this was committed at panel time; this doc
is the record so the next agent does not re-derive that (a) the profiles
moved, (b) the gather's access pattern is the unworked axis, and (c) the
tensor/NPU framing is a mirage for the constant leaf.

What actually happened next diverged from this list: rather than working down
it in order, VD redirected the same day to a lever the panel had filed as
minor (section 3c) - reduced-precision storage. Item 1 above (the re-profile)
still got run, just folded into that investigation instead of standing alone.
Items 2-4 (run-count census, histogram-fusion prototype, subsample sweep) were
not run - the block-fusion prior recorded in section 2 above made 2a/2b
unattractive as a first build, and precision displaced them as the priority.

## 8. Precision: from a screening-only idea to the shipped lever

VD re-elevated mixed-precision storage from section 3c's "minor,
screening-only" framing to a top lever the same day the panel ran: fp64 is R's
default, not a design commitment, and reduced-precision STORAGE is acceptable
for large gains. The idea: store the residual/fits in fp32 (fp64-accumulated
reductions) - halving the bytes on the bandwidth-bound streaming/gather passes
and roughly doubling the LLC-resident n (DRAM-bound -> cache-bound at the
common n>=1e5). This is the one genuinely-NEW memory-wall attack in the whole
panel, orthogonal to block-fusion's binning. Plausible upside: ~1.5-2x at
large n, single-core, CRAN-shippable, one re-record, cross-ISA bitwise
preserved (IEEE-deterministic). The gate is STATISTICAL, not architectural: SBC
/ interval coverage vs the fp64 sampler, because MCMC stationarity can amplify
a per-sweep bias a single-shot fp32 loss would not show. fp32-with-fp64-
accumulation is the safe target (~1e-7, below the sampler noise floor); bf16
(~3 digits) is the aggressive version. This superseded the histogram flagship
(section 2) as the recommended first build. Two experiments followed
immediately: confirm the fp32-vs-fp64 bandwidth/cache-crossover upside
(microbench, 8a), then measure the gather's actual share of a sweep (8b) to
turn the isolated-kernel number into an end-to-end estimate (8c).

### 8a. fp32 microbench: measured

Ran a throwaway two-pass memory microbench (streaming reduce + shuffled-index
gather-reduce, fp64 store vs fp32 store with an fp64 accumulator, min-of-7,
~3e8 elements/measurement, n from 1e4 to 3e7) on BOTH the M1 Max (arm64,
~48 MB SLC, huge single-core BW) and the x86 bench box (Ryzen 7 3700X, 16 MB
L3/CCX, DDR4). Both index widths measured: size_t (8 B, the current engine)
and uint32 (4 B, a narrowing that is a SEPARATE lever). Raw numbers preserved
off-tree. The prediction above ("top lever, ~1.5-2x at large n, ship it") is
REFINED on three points:

- THE WIN IS GATHER-ONLY AND CACHE-CROSSOVER-SHAPED, not asymptotic. fp32
  speeds the shuffled gather by ~1.7x (M1) / ~2.0x (x86) in the band where
  fp32 stays LLC-resident but fp64 spills (n ~= 1e6 on both), then TAPERS to
  ~1.1x once both are deep in DRAM (n >= 1e7). It does not keep widening with
  n: the size_t index (8 B, common to both stores) caps the gather ratio near
  16/12 = 1.33x, so at large n the win shrinks, not grows. Narrowing the index
  to uint32 lifts the asymptote toward 1.5x and is worth a modest extra gather
  gain (measured), but it is orthogonal engine work, not part of the fp32
  store.
- THE STREAMING PASS IS MACHINE-DEPENDENT AND CAN LOSE. On the M1 the
  sequential reduce is COMPUTE-bound on the fp64-add accumulator chain (~0.19
  ns/elem flat from L1 to 240 MB, ~41 GB/s single-core well under the BW
  ceiling) - so fp32 buys ZERO on streaming there. On the x86 the streaming
  reduce IS bandwidth-bound at large n (fp64 GB/s collapses 57 -> 17.6 as it
  spills L3), and fp32 wins ~1.3-1.4x - BUT cache-resident (n <= 1e6) the x86
  fp32 reduce is ~2x SLOWER than fp64, because the per-element fp32->fp64
  conversion (cvtss2sd) is the bottleneck when the bytes are already free. So
  the fp32->fp64 conversion is NOT free, and an UNCONDITIONAL fp32 store
  REGRESSES the common small/medium-n cache-resident case. The store would
  have to be n-gated (opt-in or an auto-switch above a footprint threshold),
  which is real design complexity the earlier framing did not price in.
- NET: fp32 storage is a large-n (n >~ 1e6), gather-dominated, cache-residency
  lever worth ~1.3-2x ON THE BOUND PASSES IN THAT BAND, tapering at the very
  largest n and inverting to a loss for small n. Whether it yields a
  worthwhile END-TO-END sweep speedup now turns entirely on the gather's
  SHARE of a sweep - the re-profile in 8b.

### 8b. Gather sweep-share re-profile: measured (x86 bench box)

rdtsc-bracketed computeLeafStats (the constant-leaf node-average gather) vs a
full sample sweep, single-threaded, gaussian, 200 trees, 30 burn + 30 sample,
commit 2941808. rdtsc overhead <0.1% at these leaf sizes (tight, not an upper
bound):

    n=50k   gather share 35.4%   (reference)
    n=1e6   gather share 37.8%
    n=2e6   gather share 44.2%   (large AND growing with n)

The numerator is ONLY the node-average gather; the Metropolis proposal-time
candidate-suffstat scans also stream the data and would also benefit from an
fp32 store but are NOT counted - so the fp32-accelerable fraction is a LOWER
BOUND at ~38-44% for n>=1e6.

### 8c. Verdict: GREEN for n>=1e6, but n-conditional - a real fork

Applying end_to_end = 1/(1 - share*(1 - 1/gather_speedup)) with the measured
per-n gather speedups (M1 1.67x / x86 2.02x at n=1e6; ~1.7x at n=2e6):

    n=1e6   ~1.21x end-to-end   (lower bound; proposal scans lift it)
    n=2e6   ~1.22x end-to-end
    n=1e5   ~neutral            (gather still cache-resident; ~1.0-1.16x)

So the lever is REAL and above the +/-2% bench floor for the n>=1e6
single-chain workloads (the standing "n>=1e5 common" fact; the win
concentrates at the n>=1e6 end), and ~neutral-to-slightly-negative below. The
value proposition is MORE CONDITIONAL than the initial greenlight assumed:
~1.2-1.3x at large n (not the ~1.5-2x headline, which was the isolated-gather
kernel number, not end-to-end), n-gated to protect small n, one re-record, an
SBC/coverage gate, and permanent opt-in surface area. That is a genuine fork
for VD, not a default next step, because the next move is a multi-day
draw-changing engine arc.

The recommended prototype shape at this point was an OPT-IN control flag
(e.g. residual.precision = "single"), NOT an automatic n-switch (keeps the
repro/re-record story clean; small-n users simply do not opt in), starting
with FORK B (fp64 MASTER residual + fp32 SHADOW read by every tree's gather,
shadow refreshed by a once-per-sweep downcast): SBC-safe by construction
(master is exact, no accumulated roll bias), the gather reads the small
resident fp32 shadow, and the downcast is ~1 streaming pass amortized over
~200 tree-gathers (<1% overhead) - falling back to FORK A (sole fp32 store,
fp32 roll - max footprint win but accumulated-bias SBC risk) only if the
shadow eroded the speedup.

This Fork-B-first recommendation did not survive detailed design: the
residual is a per-tree incremental roll inside the backfit loop, so a
fp64-master + fp32-shadow scheme would have to refresh the shadow after every
tree's roll - an O(n) downcast on the same order as the gather it feeds, with
no amortization, strictly worse than rolling fp32 directly. Fork A shipped
instead. See reduced-precision-storage.md section 3b for the full argument.

## 9. Landed outcome: what shipped, and why no ruled-out lever revives

The precision lever from section 8 was BUILT and validated (full record:
reduced-precision-storage.md). Landed: a uint32 gather index (lossless,
default-on, ~400MB saved at n=5e5; commit 2980229) and an opt-in fp32
residual storage="single" (commit d384211). The re-profile section 7's
experiment list kept asking for is done (x86, n=1e6): the shipped default
gather share is 39% (unmoved from the historical 37.8% - the gather is
fp64-RESIDUAL-bandwidth-bound, so the index narrowing held its ratio); fp32
residual drops it to 33% and FLATTENS the sweep to near-even thirds (gather
33 / roll 32 / other 35). fp32 pays 1.10x single-thread -> 1.30x at 4 chains
on the box (parallel BART is DDR4-bandwidth-bound). Full detail:
reduced-precision-storage.md section 6.

Verdict on reviving ruled-out levers: NONE. The cuts made the CPU faster and
the data more cache-resident, RAISING the offload bar (GPU still dead) and
SHRINKING traffic-trick payoff (binning/atoms/block-fusion/NT-stores all
deader than before). fp32 loosens the memory bound on within-chain
parallelism (blocked-jacobi-trees) by ~30%, but that is a memory-side
mitigation only - its go/no-go rested on a separate statistical (ESS/sec)
experiment, unaffected by storage. That experiment has since run to
completion, independently of this doc: both within-chain threading and
blocked-jacobi-trees are NO-GO on the real engine (the binding wall there is
the serial fraction, not just bandwidth) - see within-chain-threading.md and
docs/plans/blocked-jacobi-trees.md for the full argument; it is not re-derived
here. fp32 storage IS the memory-wall answer for this engine.

## 10. Re-profile and census at 06f73b0 (2026-08-04, dbarts-bench)

The re-profile this doc gated its flagship on has run (x86 Ryzen 3700X;
perf blocked at paranoid=4, so a 1 kHz SIGPROF PC+stack sampler loaded
into R, symbolized offline - the 39.0% node-average share at n=1e6
reproduces section 2's rdtsc figure exactly, cross-validating the
method). Fresh split, gaussian Friedman p=10, 200 trees, single chain:

  pass                                       n=1e5   n=1e6
  random suffstat gather (total)             46.4%   48.8%
    node-average only (computeLeafStats)     40.2%   39.0%
    move/refresh child suffstats              6.1%    9.8%
  residual roll + totalFits (streaming)      27.9%   24.5%
  partition/compare on codes (avx2)          22.0%   20.8%
  bookkeeping / RNG / rest                    3.7%    5.9%

Three verdicts:

- HISTOGRAM-FUSED SUFFSTAT (sec 2a, the flagship): premise SURVIVES.
  The node-average gather is 39-40% of the sweep at both sizes and is
  exactly the pass 2a folds into the residual roll, which already loads
  leafOf[i]; ceiling ~1.67x on the sweep. Block-fusion Stage A's
  bench-neutral result stands as the counter-evidence, but its premise
  objection is removed: the gather grew, not shrank, since the fits
  compaction. Next artifact, behind a VD fork: a throwaway fused-kernel
  falsifier measured in situ (never a microbench alone;
  within-chain-threading.md sec 10a).
- Z-ORDER / LOCALITY REORDERING (secs 2b/2c): measured DEAD. Census on
  live post-burn trees (R-side replay validated against every node's
  reported n): best case (Z-order over the 5 relevant predictors) is
  2617 runs/leaf at n=1e5 against the recorded win threshold of ~50,
  and 15207 at n=1e6 - essentially the n_leaf/2 "dead" mark. With 3-5
  leaves per tree the identity gather already fetches 17-22 B/element
  against a contiguous floor of 8, capping ANY reordering at 2.1-2.8x
  on the residual stream alone before its own costs. Closed.
- DATA-LAYOUT REORDER (data-layout.md): re-evaluated and kept shelved;
  the dated close-out lives in that doc's post-mortem.

## 11. Fused-kernel falsifier: measured GREEN at K = 4 (2026-08-07, dbarts-bench)

The section-10 next-artifact ran (throwaway patch, chain.hpp only,
+147/-4, worktree wt/fused-falsifier): the roll and the node-average
suffstat fused into one obs-order pass, acc[leafOf[i]] += resid[i]
into a node-indexed accumulator, scalar and fixed-order, selected by
DBARTS_FUSED_SUFFSTAT with K banks combined in fixed order. The fused
path declines to stock for non-constant leaves, ResidT = float,
weights != nullptr, and a stale leafOf - the post-sampleTreesFromPrior
hazard the 2a sketch missed. In situ on Zen2 (200 trees, friedman
p=10, 1 chain/1 thread, min of 7 x 3 rounds, base and fused
back-to-back, 1-min loadavg 0.87-1.45 throughout); ratio = base/fused
ms per sample, per-round minima:

  K = 1: n=1e5 0.92-0.95 (a LOSS), n=1e6 1.15-1.20, no-signal worst
  case 0.67 - the single-accumulator FP dependency chain sec 2a
  predicted. NOT green against the pre-registered >= 1.25 line.
  K = 4: n=1e5 1.41-1.43, n=1e6 1.49-1.54, worst case a 1.11x WIN.
  GREEN in every round. 1.53x implies the fused pass costs ~7% of
  the sweep against the gather's 39%, independently corroborating
  the section-10 profile. Block-fusion Stage A's bench-neutral prior
  does NOT reproduce post-leafOf-refactor.

Privatization is the MECHANISM, not a mitigation, and K joins the
exactness contract: different K sum in different orders. Sanity vs
base (n=5000, 200 trees, 500+500): sigma posterior mean and rmse
agree to ~2e-15 relative, varcounts to 4 dp, node sums differ at
<= 1.4e-11 relative, and no MH accept flipped over 1000 sweeps - the
draws are NOT bitwise (snapshot re-record at landing, the sec-1
shifting class) but the measured shift is last-ULP, far milder than
the class name suggests. UNMEASURED: weighted, fp32 residual, the
move-phase child gathers, vector/function leaves, BCF/multinomial,
monotone, growForestFromRoot, K = 8, and any arm64 in-situ leg - the
arm A/B is the confirmation to demand before an arc lands
(within-chain-threading.md sec 10a). The real arc runs design +
refuting critique on that scope before implementation.

## 12. arm64 in-situ leg: direction confirmed, no-signal inverts (2026-08-07, M1 Max)

The section-11 arm A/B ran on an Apple M1 Max (8P+2E, macOS), same
throwaway patch (wt/fused-falsifier, K = 4) against pristine 63456a0,
mirroring the Zen2 protocol: friedman p=10, 200 trees, 1 chain/
1 thread, ms per sample after a 200-sweep burn-in, min of 7 reps,
base and fused back-to-back, 3 rounds. Per-round base/fused ratios:

  run-n1e5:   1.071 / 1.072 / 1.071  (base 43.4-43.8, fused 40.6-40.9)
  run-n1e6:   1.228 / 1.231 / 1.246  (base 456.7-463.7, fused 371.0-377.7)
  noise-n1e5: 0.951 / 0.948 / 0.949  (base 39.7-40.1, fused 41.9-42.2)

Direction confirmed at n >= 1e5 and the win still grows with n, but
far below Zen2 (1.07x/1.23-1.25x vs 1.41-1.43x/1.49-1.54x) -
consistent with the M-series' bandwidth headroom lowering the wall
the fusion attacks, though that attribution is unmeasured. The
no-signal worst case INVERTS: a consistent ~5% loss (0.948-0.951)
where Zen2 measured a 1.11x win. Mechanism unmeasured; the design
memo must weigh that regression explicitly (n-gate, arch-gate, or
accept) - this record does not decide it.

Protocol note: the rounds ran at 1-min loadavg 1.9-3.1 on a
GUI-loaded desktop - the armab loadavg < 2 gate never opened and was
relaxed to < 6 (macOS loadavg counts bursty runnable GUI threads; a
single-threaded R process on 8 P-cores sees little of them). Validity
rests on the alternating min-of-7 design plus cross-round agreement,
which is tight: n=1e5 ratio spread 1.071-1.072 across rounds taken at
different load. A fourth partial round (killed mid-run) corroborates
(1.065 at n=1e5, 1.260 at n=1e6). Gate future single-threaded A/Bs
on this box at loadavg < 6 and judge by cross-round spread, not
absolute quiet.

## 13. Landed: the fused pass ships (2026-08-07, c8f661a)

The production kernel landed as c8f661a: rollAndSetNodeAveragesFused
in chain.hpp, always on where eligible (constant leaf, ResidT double,
null weights, fresh leafOf map), the stock roll + setNodeAverages
pair everywhere else, K = 4 constexpr with the fixed positional-bank
contract, and NO root gate. The root gate the design first proposed
died on a census: stock noise-n1e5 carries only ~20% bare roots with
69% of trees single-split (friedman ~2%), so the gate could neither
capture the arm64 no-signal loss nor fire in signal fits. v1 accepts
the ~5% arm64 no-signal loss; the pre-registered arm re-run stands.
The roll expressions are bitwise identical to rollTreeResidual,
enforced by a Chain test hook that runs fused and stock per tree over
the same entering residual, so the suffstat summation association is
the only draw change anywhere. Probit, ordinal, AFT/survival, and
grouped random effects ride the fused path (null workingWeights);
BCF, multinomial, logistic, negbin, t, and the heteroscedastic mean
forest decline through their always-on weights.

Gates at landing, run independently of the implementer: tests/cpp
plain and ASAN/UBSAN clean from make clean; tinytest 3659/0;
equivalence vs 7903855 splits into exactly the predicted 11 bitwise
decline scenarios and 16 covered scenarios at max |z| = 0.00; BCF
and multinomial equivalence bitwise on every channel; air clean.
Baseline re-recorded as equivalence-c8f661a.rds (self-reproduces
27/27 under --strict-coverage); bcf-equivalence-99205ee and
multinomial-equivalence-ec2a3d0 deliberately not re-recorded - their
staying bitwise IS the scope gate. PENDING: the bench-sampler speed
baseline re-records on the next quiet-machine grant. The
wt/fused-falsifier throwaway probe is superseded by this landing and
removed.

The pre-registered arm64 no-signal re-run RAN 2026-08-08 (M1 Max,
production install of the landed kernel at tip d4e8c2a vs pristine
63456a0, sec-12 protocol: n=1e5, 200 trees, 1 chain/1 thread, min of
7, base/fused back-to-back, 3 rounds, 1-min loadavg 2.9-4.5).
Per-round base/fused ratios: noise-n1e5 0.964 / 0.970 / 0.950 - a
3-5% loss, matching the probe's 0.948-0.951 and far inside the
pre-registered acceptance (escalation line ~10%, memo sec 5.3);
run-n1e5 control 1.058 / 1.079 / 1.102, reproducing the sec-12 win
and corroborating build provenance (the fused install carries
rollAndSetNodeAveragesFused, the base install does not). The
no-signal loss is ACCEPTED and the re-run obligation is discharged;
reopening now requires a no-signal loss beyond ~10% or an inversion
in a real-signal scenario.
