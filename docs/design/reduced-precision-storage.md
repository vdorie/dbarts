# Optional reduced-precision / narrowed storage for the BART engine

Status: DESIGN (2026-07-20). Draw-changing where noted; design-first + an
independent blind critique precede any implementation. Motivation and the
measured upside live in memory-wall-frontier.md sections 7-8; this doc is the
concrete build design. Directive (VD 2026-07-20): "Build it, keep it optional,
think about other ways to optionally decrease storage sizes at the same time."

## 1. What and why

The per-sweep O(n) passes are memory-bound, dominated by the random suffstat
GATHER (x[indices[i]] over a shuffled index buffer). Measured (frontier sec 8):
an fp32 running residual buys ~1.2-1.3x end-to-end at n >= 1e6 (gather is
37.8-44.2% of a sweep and fp32 speeds it ~1.7-2.0x in the cache-crossover band),
tapering to neutral at n <= 1e5. The goal: make the hot per-observation storage
optionally narrower, shrinking memory traffic and doubling the LLC-resident n at
large scale, WITHOUT touching the default fp64 path.

## 2. Architecture: a compile-time storage axis behind the existing facade

INVARIANT (non-negotiable, the blind critique must verify it): the default
instantiation stays BITWISE-IDENTICAL and pays ZERO new runtime cost.

The facade already resolves runtime choices into monomorphic instantiations:
`createSampler` (facade.hpp:502) reads the family / IsCategorical /
ServesRawValues / leaf-covariate axes and dispatches ONCE to build a concrete
`SamplerFacade<L>` (facade.hpp:188), handing the R bridge a type-erased
`SamplerBase*`; every later call is "one virtual hop into fully typed code"
(facade.hpp:14). The sweep is fully monomorphic - no per-element dispatch.

Storage precision/width is ONE MORE axis on that factory. Thread a `Storage`
policy through `SamplerFacade<L, Storage>` and the Chain; the factory gains one
construction-time branch on the opt-in flag. Consequences (see the answer
recorded in frontier reasoning): the `Storage = default` instantiation is the
code the compiler emits today (byte-identical, no new branch/virtual on the hot
path); the cost is COMPILE-TIME + binary size (each storage type multiplies the
instantiations of whatever is templated on it) + one construction-time branch.
The user interface stays fp64 (R passes double; the downcast is internal at
ingestion), so dbarts.h / the ABI is UNCHANGED - no hash bump, no stan4bart
lockstep. Bound the compile-time blowup by instantiating narrowed storage only
for the branches where it pays (the continuous/gaussian residual path; latent
families - probit/logistic/ordinal/multinomial - keep fp64).

ANTI-PATTERN to avoid: a runtime precision TAG branched per element, or a
type-erased element accessor. That taxes the default path. Dispatch once,
monomorphize the core - exactly what the facade already does.

## 3. The levers, by draw-impact class

Grounded in the storage inventory (2026-07-20 sweep of combiner.hpp/chain.hpp/
tree.hpp/data.hpp/model.hpp). KEY INVENTORY FACTS: (a) the predictor store is
ALREADY uint16 cutpoint codes (xint_t, data.hpp:19/158) - splits run on integer
codes, raw doubles are quantized ONCE at setup, so the "fp32 predictor flips a
split" hazard does NOT exist on the hot path and the biggest n*p prize is
already spent; (b) the LARGEST hot array is the gather index (indexBuffer,
size_t, n*trees) - bigger than the predictor codes; (c) leafOf is already uint32,
categorical masks uint64 bitsets, several flags uint8 - no float anywhere today.

### 3a. PRESERVING (bitwise-identical draws; NO SBC, NO re-record; gated only by
the equivalence trio staying identical). DEFAULT-ON - a separate, lower-risk
TRACK 1 that lands first and independently of the fp32 work.
- **Index narrowing size_t -> uint32** (Forest::indexBuffer combiner.hpp:128, its
  Tree::indices alias tree.hpp:186, and VarianceForest::indexBuffer chain.hpp:311).
  THE TOP LEVER: footprint n*trees*8B (the single biggest hot array; ~1.6 GB at
  n=1e6/200 trees) -> halves to 4B, and it halves the streamed-index bytes inside
  the DOMINANT gather (microbench: index narrowing alone gives ~1.21x on the
  gather at n=1e6, g64->g32, and lifts the fp32 asymptote 1.33x->1.5x).

  PRESERVING - VERIFIED by the blind critique: the SIMD partition (partition_body.c)
  vectorizes the COLUMN VALUES fed to the code compare and never packs the index
  buffer into vectors (indices are only scalar subscripts), so its swap sequence /
  control flow is independent of index width -> a byte-identical permutation ->
  identical gather order -> bitwise-identical fp64 sums and draws. State save/
  restore does NOT serialize index width (ForestStateData/ChainStateData hold only
  FlatNodes + fp64; buffers are reconstructed by repartition). n < 2^31 always
  holds (predictor codes would be terabytes first); a guard throws above UINT32_MAX.

  BLAST RADIUS (critique correction - it is NOT a one-line typedef): misc_size_t
  is a single macro used for BOTH indices AND lengths, so there is no index typedef
  to flip. The change introduces a NEW misc_index_t DISTINCT from the length type
  and threads it through the whole misc.a index family: the two runtime-dispatched
  partition fn-pointers misc_partitionRange / misc_partitionIndices compiled across
  5 ISA TUs (avx2/sse2/sse4_1/neon + scalar body), misc_partitionIndicesSparse,
  the 16 misc_computeIndexed* suffstat/mean/variance kernels AND their misc_mt_
  variants, misc_sumIndexedVectorElements, misc_setIndexedVectorToConstant; plus
  Tree::indices (tree.hpp:186), the 8 scalar C++ partition kernels, and
  SubtreeSnapshot.indexSegment (std::vector<size_t> + a sizeof(size_t) memcpy,
  tree.hpp:836). THE TRAP: preservation holds ONLY if each SIMD kernel is RETYPED,
  not REWRITTEN - a fresh uint32 SIMD partition with a different swap sequence
  changes the permutation and silently turns a "preserving" lever into a draw-
  changing one. GUARD: static_assert(sizeof(index_t) == sizeof(misc_index_t))
  pinning the C++ buffer type to the C kernel param (mirroring the existing
  static_assert on misc_xint_t == uint16 at tree.hpp:618). ACCEPTANCE: the
  cross-ISA tests/cpp gate AND the bitwise equivalence trio must both stay
  identical - any drift means a kernel was rewritten, not retyped.
- **leafOf uint32 -> uint16** (combiner.hpp:143, n*trees). PRESERVING (node ids
  are exact) and a GENUINE gather input (mu[leafOf[i]] streamed every roll, so the
  reward is real), but CAPACITY-BOUNDED: a tree exceeding 65535 nodes overflows.
  Shallow BART trees (base=.95/power=2) astronomically never approach this, but the
  fallback must be a real TYPE spec, not a hand-wave: a uniform std::vector<uint16>
  cannot host a per-tree fallback, so the guard is a GLOBAL refuse (throw/cap at
  the first tree that would exceed 65535 nodes), not a per-tree type union. uint8
  is NOT safe (>255 nodes reachable). leafOf is internal scratch (not in
  ForestStateData) -> no state/ABI exposure. CONSIDER for Track 1 as a follow-on
  to the index change (same PRESERVING gate); defer if the global-refuse guard
  reads as a surprising failure mode.

### 3b. CHANGING (feeds precision-sensitive arithmetic; SBC/coverage gate + one
re-record; OPT-IN only). The primary target.

- fp32 running residual (treeY). DESIGN = FORK A (the master+shadow Fork B is
  DEAD - see below). Make treeY itself fp32: the per-tree roll writes fp32, the
  per-tree gather (setNodeAverages -> computeLeafStats) reads fp32, and the
  proposal-time suffstat scans read fp32. Both wins (half the roll-write bytes,
  half the gather bytes + cache-residency) come from treeY being fp32. Leaf draw
  stays fp64 (sumWeights / sumWeightedResponse are fp64 accumulators fed by the
  float-loading suffstat kernels); the four misc_compute*SufficientStatistics-
  Fast kernels get float-input variants (load float, accumulate double).

  WHY FORK B IS DEAD (chain.hpp:891-948): treeY is a running residual updated
  PER-TREE inside the backfit loop - `rollTreeResidual(forest, t, ...)` at :902
  rewrites it, then the gather at :906 reads tree t's freshly-rolled residual. A
  fp64-master + fp32-shadow scheme would have to refresh the shadow after every
  roll = a per-tree O(n) downcast, the same order as the gather it feeds -> no
  amortization, strictly worse than rolling fp32 directly. So the fp32 store
  cannot avoid the incremental fp32 roll.

  THE ROLL IS INCREMENTAL (rollTreeResidual chain.hpp:2581: t=0 recomputes
  resid = y - total + mu[leaf]; t>0 does resid += mu[leaf] - muPrev[leafPrev]),
  so fp32 rounding ACCUMULATES across trees (~200 updates/sweep/element) and,
  because it is never independently re-summed, across sweeps. This is the real
  and ONLY safety question for Track 2 (there is no "safe by construction").

  DRIFT CONTROL - corrected after the blind critique (the earlier "free
  finalizeTotalFits re-anchor" claim was WRONG):
  - finalizeTotalFits (chain.hpp:2613) computes total = y - resid + mu[leaf] -
    i.e. FROM the drifted resid - so it PROPAGATES the fp32 drift into totalFits,
    not resets it; next sweep's t=0 roll re-derives resid from that drifted
    total. Nothing in the existing per-sweep machinery re-anchors. A genuine
    re-anchor needs an INDEPENDENT total = sum_t mu_t[leafOf_t[i]] (a fresh
    per-sweep O(n*trees) re-sum). That re-sum gathers into the L1-resident mu
    tables (not the DRAM residual), so it is far cheaper per element than the
    treeY gather (~10-15% of it), i.e. affordable - but it is NOT free and NOT
    the existing pass. It bounds CROSS-sweep drift.
  - The draw-relevant error is WITHIN-sweep and IRREDUCIBLE below per-tree
    cadence: the gather at tree t (setNodeAverages) sees ~t accumulated fp32
    roundings since the last t=0 anchor. No re-anchor coarser than per-tree
    removes it. Magnitude ~sqrt(200)*1e-7 ~ 1e-6 relative, and it is ZERO-MEAN
    (round-to-nearest), and it is further suppressed by fp64 leaf-sum averaging
    (~1/sqrt(leaf size)); so the PRIOR is it sits below the MCMC noise floor -
    but that is a hypothesis SBC must confirm, not a guarantee.
  So Track 2's viability rests ENTIRELY on the STATISTICAL gate: prototype fp32
  treeY with NO re-anchor and measure (a) the drift max_i|treeY_fp32 - exact|
  growth over a long n=1e6 run, and (b) SBC / interval coverage vs fp64. If SBC
  passes un-anchored -> ship raw fp32 (cheapest). If drift/SBC fails -> add the
  independent per-sweep re-sum re-anchor (bounds cross-sweep; re-benchmark the
  net win, since it eats ~10-15% of the gather) and re-test; the within-sweep
  1e-6 is the floor either way. This falsification is the FIRST Track-2 action,
  BEFORE the full mixed-precision templating (which is large source churn).
- fp32 SCRATCH/FITS BUNDLE (same opt-in flag, same SBC gate as treeY, gaussian/
  continuous path only): Forest::totalFits (combiner.hpp:132, n), test-fit
  accumulators totalTestFits/currTestFits (nTest), the non-constant-leaf treeFits
  slab (n*trees, dense-leaf models only), the variance-forest arrays for HBART
  (factorByTree n*treesVar, combinedVariance/meanResidual/divisor/treeResidual n
  each, chain.hpp:314-319), and the gaussian working response yRescaled_
  (model.hpp, n). LATENT-family working buffers (probit/logistic/ordinal/
  multinomial) stay fp64 by design (little value, more instantiations). These
  extend the fp32 win to the streaming passes and the HBART weight channel; each
  is individually smaller than treeY but shares the flag and gate.

### 3c. CORRECTNESS-SENSITIVE (a discrete decision, not just rounding). EXCLUDED
from every track - flagged so no one bundles them as a transparent memory tier.
- Predictor codes uint16 -> uint8 (data.hpp:158): the largest n*p array but
  ALREADY uint16; narrowing further lowers maxNumCutsRepresentable (~254),
  capping cutpoints and silently changing which splits exist -> different tree
  structure, pathological on high-cardinality columns. Only ever an explicit
  low-resolution MODELING option with its own arc, never a transparent tier.
- fp32 cutpoints (data.hpp:221, tiny footprint) and fp32 ownedTestValues
  (data.hpp:251, cold re-quantize source): fp32 flips near-tie quantization ->
  correctness-sensitive; not worth the risk, and cutpoints are not n-scaled.

## 4. Opt-in surface (fork for VD)

- Option 1 (coarse): a single memory-tier flag on dbartsControl, e.g.
  `storage = c("double", "single")` (or a `low.memory` toggle) that flips a
  curated bundle. Matches how users think ("big-n, reduce memory"), not
  "narrow my index array". PRESERVING levers default-on regardless.
- Option 2 (fine): per-lever knobs. More control, more surface, more
  instantiations to compile and test.
- Lean: coarse flag for the CHANGING bundle (fp32 residual [+ fp32 scratch]);
  PRESERVING levers (uint32 index) default-on as a separate track.

## 5. Gates

- PRESERVING track: equivalence trio must stay BITWISE-IDENTICAL (no re-record);
  tests/cpp; sanitizer for new reachable code.
- CHANGING track (fp32 residual, opt-in): SBC / interval coverage vs the fp64
  sampler is the ARBITER (MCMC stationarity can amplify a per-sweep bias);
  same-machine A/B speed at n >= 1e6; the DEFAULT-path equivalence baselines
  stay identical (only the opt-in path shifts, so no default re-record); ASAN.
- The default fp64 instantiation must diff-clean against HEAD in generated
  assembly for the hot inner loop (spot-check) - the INVARIANT.

## 6. v1 scope: two tracks, sequenced (proposed; finalize after critique)

TRACK 1 (PRESERVING, default-on, no SBC, no re-record) - LANDED 2980229
(index_t = uint32; retype-not-rewrite across the misc.a index family; gates re-run
independently, equivalence trio BITWISE, ~400MB saved at n=5e5). leafOf uint16 is
a deferred smaller follow-on (needs a global-refuse >65535 guard).

TRACK 2 (CHANGING, opt-in `storage="single"`, gaussian constant-leaf path) -
falsification PASSED (6c) -> GO to implement, UN-ANCHORED:
- STEP 2.0 (falsification) DONE - GO. See 6c.
- STEP 2.1 v1 (implement): real fp32 treeY (std::vector<float>, FORK A,
  UN-ANCHORED - the re-anchor is deferred insurance, not built) + the float-input
  suffstat kernel variants (load float, accumulate double), behind the opt-in flag
  via the facade Storage axis (default=double byte-identical). Gaussian constant-
  leaf only; gate to it. NO re-anchor machinery.
- STEP 2.1 v2 (follow-on, if v1's measured win warrants): the fp32 scratch/fits
  bundle (totalFits, test fits, dense treeFits, variance-forest arrays, gaussian
  working response) on the same flag. The microbench says the GATHER (treeY)
  carries the win and streaming scratch gains little, so v1 (treeY) captures most
  of the value - measure before expanding.
- GATES for 2.1: full SBC/coverage (landing arbiter) + same-machine A/B at n>=1e6
  to confirm the REAL fp32 memory win (the falsification prototype used round32-on-
  double, so it did NOT realize the speedup - that is still to be measured) + ONE
  re-record of the opt-in path only (default path stays bitwise) + ASAN + tests/cpp
  + tinytest.

OUT: predictor codes uint8 and fp32 cutpoints/test-raw (correctness-sensitive,
sec 3c); bf16 (the aggressive tier; future); latent-family fp32 (low value).

## 6b. Blind-critique resolutions (2026-07-20; verdict SOUND-WITH-CAVEATS T1 /
NEEDS-REWORK-DONE T2)

- BLOCKER (T2 drift): the "finalizeTotalFits gives a free per-sweep re-anchor"
  claim was FALSE - it derives totalFits FROM the drifted residual, propagating
  the error; and the draw-relevant error is within-sweep, irreducible below
  per-tree cadence. FIXED in 3b: dropped the free-safety claim; Track 2 now rests
  on an explicit SBC falsification (step 2.0); a real re-anchor, if needed, is an
  independent O(n*trees) mu[leafOf] re-sum (cache-cheap but not free).
- SERIOUS (T1 scoping): index narrowing is genuinely bitwise-preserving (verified)
  but is NOT a one-line typedef - misc_size_t is shared with lengths; a new
  misc_index_t must thread through the whole misc.a index family across 5 ISA TUs,
  RETYPED-not-rewritten, with a sizeof static_assert pin and the cross-ISA +
  equivalence gates as acceptance. FIXED in 3a.
- CONFIRMED in the design's favor: the fp32 gather reduces in double (kernels
  already accumulate double); no hidden runtime dispatch (Storage is an if-constexpr
  compile-time axis, default codegen byte-identical); compile blowup bounded
  (+1 instantiation per offered L, not the cross-product); threading safe
  (chain-parallel, serial gather kernels); warm-start/state safe (fp64 wire);
  stan4bart/ABI unchanged; ownership unchanged.
- CAVEAT (T2 churn): "compile-time cost only" understates SOURCE churn - Storage
  threads through Forest<L>/Chain<L>/Sampler<L> and every hand-written mixed-
  precision loop (rollTreeResidual is where rounding enters). Default CODEGEN
  stays byte-identical (enforce via the hot-loop assembly spot-check), but the
  source around it is heavily rewritten - budget accordingly.

## 6c. Falsification result (2026-07-20; step 2.0) - GO, UN-ANCHORED

Throwaway fp32-treeY prototype (isolated worktree on the 2980229 tip; round32(x)
= (double)(float)x on every residual write, provably bit-identical to a true
std::vector<float>+float-load/double-accumulate since float->double promotion is
exact and every consumer already reduces in fp64). Flag-off byte-deterministic;
snapshot suite 186/186. Confirmed genuinely fp32-storing/fp64-reducing (drift
1e-15 off -> 1e-6 on, a 9-order gap).
- DRIFT: a slow RANDOM WALK (log-log slope ~0.5, technically unbounded because
  finalizeTotalFits propagates rather than re-anchors), but TINY - rms_rel ~
  7.9e-8*sqrt(sweeps); even 1e6 sweeps -> ~7.9e-5 rms / 7.5e-4 max relative to
  sd(y). Does NOT worsen with n (bigger leaves -> more fp64-averaging).
- SBC/COVERAGE (arbiter; fp32-vs-fp64 paired + MC-noise control, 30 seeds,
  n in {1e3,1e4}): fp32 INDISTINGUISHABLE from fp64 - no detectable bias in
  test RMSE, 90% CI coverage of true f, or the sigma posterior (all p >> 0.05);
  the fp32-vs-fp64 divergence is 3-6x SMALLER than an RNG reseed. (Absolute BART
  under-coverage/small-n overfit are identical across arithmetic modes - orthogonal
  to fp32.)
- RE-ANCHOR (built + measured as insurance): the independent per-sweep mu[leafOf]
  re-sum flattens the walk (slope -> ~0, rms ~1e-6) and makes fp32 near-bit-
  identical to fp64 - but costs ~+20% of a sweep (naive, unfused), a large fraction
  of the fp32 gather win. VERDICT: ship UN-ANCHORED (SBC passes, cheapest); keep
  the re-anchor as a documented, un-built fallback for >1e5-sweep chains or a
  stricter latent family.
- CAVEAT: gaussian constant-leaf only; SBC = paired-comparison proxy, not full
  rank-uniformity SBC (the full-implementation landing gate should run the fuller
  check). The prototype measured DRAWS faithfully but not the SPEED win (double
  storage + round op) - the real fp32 storage speedup is still to be measured.

## 6d. Track 2 v1 LANDED + validated (2026-07-20, d384211)

Opt-in `storage="single"` on dbartsControl/bart2 -> real fp32 treeY
(std::vector<ResidT>, ResidT=double default; +1 instantiation SamplerFacade<
ConstantGaussianLeaf,float> minted only at the gaussian-constant-leaf branch) +
four float-input suffstat kernels (load float, accumulate double). Un-anchored.
Errors clearly on probit/hetero/monotone/BCF/multinomial. Independently verified:
- INVARIANT (flag-off) BITWISE: equivalence 27/27 + bcf + multinomial -> default
  path provably byte-unchanged (the Storage axis defaulted to double).
- Genuine fp32 (not a silent double fallback): fp32-vs-fp64 posterior-mean
  divergence 1.9e-6, matching the ~1e-6 prediction.
- Statistically indistinguishable (real templated build, 10 seeds, n=1e4,
  200 trees): test-f RMSE Delta 5.7e-8 (p=.12), 90% coverage Delta 1e-5 (p=.34).
- codoc complete; ABI header untouched (no stan4bart lockstep).

SPEED A/B (n=1e6, single-thread, min-of-3, storage=single vs double):
  M1 Max (arm64):        1.17x (17% faster)  - matches the gather-share projection.
  x86 Ryzen 3700X (AVX2): 1.10x (10% faster)  - clean build, 0 warnings.
The x86 UNDERSHOOTS the ~1.24x projection because x86 pays a cvtss2sd conversion
tax on the CACHE-RESIDENT streaming roll (microbench sec 8: x86 fp32 streaming is
~2x slower cache-resident) that partly offsets the gather win; M1 (compute-bound
streaming) does not, so M1 matches projection. Both POSITIVE and above the
+/-2% floor. The win should GROW at n>1e6 on x86 (residual spills the 16MB L3 ->
streaming goes DRAM-bound, where fp32 helps streaming too - microbench sec 8).

## 6e. v2 (fp32 scratch/fits bundle): MEASURED NO-GO
The v1 A/B settles v2 against itself. The bundle (totalFits, test fits, dense
treeFits, variance-forest arrays, gaussian working response) is STREAMING, not
gathered - and the microbench + the x86 v1 undershoot show fp32 STREAMING gains
~nothing (M1 compute-bound) to NEGATIVE (x86 cache-resident conversion tax). So
fp32-ing the streaming scratch adds conversion overhead with no gather benefit ->
neutral-to-regressive. The gather (treeY, v1) carries the entire win. DO NOT build
v2 unless a re-profile at n>>1e6 shows a streaming pass that is both DRAM-bound AND
a large sweep share. leafOf uint16 (the deferred PRESERVING Track-1 follow-on)
remains the only other live lever, and it is small.

## 6f. Post-landing re-profile + multi-chain scaling (2026-07-20) - and why no ruled-out lever revives

x86 rdtsc re-profile at n=1e6 (single-thread, 200 trees), buckets = gather
(computeLeafStats) / roll (rollTreeResidual) / other (moves/scans/draws):
  historical (2941808, pre-cuts):  gather 37.8%
  storage=double (uint32 idx, fp64 resid):  39.2 / 29.9 / 30.9   sweep 42.4s
  storage=single (uint32 idx, fp32 resid):  33.4 / 32.0 / 34.7   sweep 37.1s
- The uint32 INDEX narrowing did NOT move the gather's share (39.2% vs 37.8%): the
  gather is fp64-RESIDUAL-bandwidth-bound, not index-bound (the 4B index is minor
  next to the 8B value), so it sped up the gather AND the other index-touching
  passes ~proportionally, holding the ratio. Track 1's value is the ~400MB saving +
  a proportional modest speedup, NOT a gather-concentrated win (corrects the
  earlier microbench framing).
- The fp32 RESIDUAL is what attacks the gather: -28% gather cycles, share 39->33%.
  Under fp32 the sweep FLATTENS to near-even thirds (33/32/35), "other" (moves/
  scans/draws - heterogeneous, no single lever) the largest slice.

MULTI-CHAIN A/B (n=1e6, fp32/fp64 speedup): the win GROWS with chains because
parallel BART shares memory bandwidth -
  x86 Ryzen (dual-channel DDR4):  1 thread 1.10x -> 4 chains/4 threads 1.30x
  M1 Max (huge BW):               1 thread 1.17x -> 4 chains/4 threads 1.19x
So fp32 pays MOST in the DEFAULT multi-chain config on commodity/bandwidth-limited
hardware (~30%), not the single-thread corner the first A/B measured. Document this
(NEWS/Rd): "storage='single' is most valuable for multi-chain runs and large n."

WHY NO RULED-OUT LEVER REVIVES: the reductions FLATTENED the profile (fp32) or held
it (index) rather than exposing a new fat target - no single pass to attack, and
the incumbent got faster + smaller/more-cache-resident, RAISING the bar for offload
(GPU: fp64-blocker dissolved but the sequential-backfit + launch killer stands and
the CPU baseline it must beat got faster) and SHRINKING the payoff of traffic-cutting
tricks (binning/atoms/block-fusion/NT-stores). The one shifted item - within-chain
parallelism (blocked-jacobi-trees) - has its memory-wall bound LOOSENED ~30% by
fp32, and the 4-chain A/B PROVES parallel BART is bandwidth-bound on commodity HW
(that IS why fp32 gives 30% there), so fp32 is a prerequisite MITIGATION for it -
but its go/no-go still rests on the statistical ESS/sec experiment, unchanged by
storage. Net: no revival; fp32 is itself the memory-wall answer, and it is
broader-value than first reported (the multi-chain default).

## 7. Open questions - status after the critique + falsification

- RESOLVED: index narrowing stays bitwise (SIMD partition is index-width-
  independent; state does not serialize width). See 3a/6b. LANDED 2980229.
- RESOLVED: Fork B dead (per-tree roll); Fork A is the design. See 3b.
- RESOLVED (step 2.0, sec 6c): un-anchored fp32 treeY drift stays far below the
  MCMC noise floor; SBC passes; re-anchor NOT needed (deferred insurance). GO.
- OPEN (step 2.1 landing): confirm the REAL fp32-storage speed win (same-machine
  A/B at n>=1e6) and pass a fuller SBC at the templated implementation.
- RESOLVED: compile blowup bounded (+1 instantiation per offered L); no hidden
  runtime dispatch. Source churn is the real cost (6b caveat).
- OPEN (minor): is n-gating needed beyond opt-in (x86 small-n conversion loss)?
  Opt-in likely covers it - large-n users opt in, small-n users do not; confirm
  the opt-in path is never auto-selected for small n.
