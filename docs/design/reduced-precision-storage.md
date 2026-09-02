# Optional reduced-precision / narrowed storage for the BART engine

Status: LANDED / COMPLETE (2026-07-20/21). Both tracks below shipped and were
independently validated; the storage arc is closed and no further lever is
queued. Motivation and the measured upside live in memory-wall-frontier.md
section 8 (the fp32-precision lever writeup); this doc is the concrete build design plus its full implementation
and measurement record. Directive (VD 2026-07-20): "Build it, keep it
optional, think about other ways to optionally decrease storage sizes at the
same time."

Summary: Track 1 (index narrowing, size_t -> uint32, bitwise-preserving,
default-on) landed at commit 2980229, saving ~400MB at n=5e5 with the
equivalence trio staying bitwise-identical - the largest hot array in the
engine, narrowed for free. Track 2 (opt-in fp32 residual storage, draw-
changing) landed at commit d384211 after a falsification gate (drift + SBC)
showed it safe to ship un-anchored; measured speedup is 1.10x-1.30x depending
on thread/chain count, and grows with n and with multi-chain bandwidth
pressure (the default multi-chain configuration on commodity hardware benefits
most). A follow-on fp32 scratch-and-fits bundle and a leafOf uint16 narrowing were
both designed and then declined by measurement - the residual gather carries
the entire win; streaming scratch does not. Bottom line: fp32 storage is the
memory-wall answer for this engine (see memory-wall-frontier.md section 9,
the landed-outcome umbrella verdict);
nothing further in the storage space remains open.

## 1. What and why

The per-sweep O(n) passes are memory-bound, dominated by the random suffstat
GATHER (x[indices[i]] over a shuffled index buffer). Measured (memory-wall-
frontier.md sec 8): an fp32 running residual buys ~1.2-1.3x end-to-end at
n >= 1e6 (gather is 37.8-44.2% of a sweep and fp32 speeds it ~1.7-2.0x in the
cache-crossover band), tapering to neutral at n <= 1e5. The goal: make the hot
per-observation storage optionally narrower, shrinking memory traffic and
doubling the LLC-resident n at large scale, WITHOUT touching the default fp64
path.

## 2. Architecture: a compile-time storage axis behind the existing facade

INVARIANT (non-negotiable, verified by the blind critique below): the default
instantiation stays BITWISE-IDENTICAL and pays ZERO new runtime cost.

The facade already resolves runtime choices into monomorphic instantiations:
`createSampler` ([[facade.hpp#createSampler]]) reads the family / IsCategorical /
ServesRawValues / leaf-covariate axes and dispatches ONCE to build a concrete
`SamplerFacade<L>` ([[facade.hpp#SamplerFacade]]), handing the R bridge a type-erased
`SamplerBase*`; every later call is "one virtual hop into fully typed code"
([[facade.hpp#SamplerBase]]). The sweep is fully monomorphic - no per-element dispatch.

Storage precision/width is ONE MORE axis on that factory. Thread a `Storage`
policy through `SamplerFacade<L, Storage>` and the Chain; the factory gains
one construction-time branch on the opt-in flag. The `Storage = default`
instantiation is the code the compiler emits today (byte-identical, no new
branch/virtual on the hot path); the direct cost is COMPILE-TIME + binary size
(each storage type multiplies the instantiations of whatever is templated on
it) + one construction-time branch. Independent verification later confirmed
this stays bounded to +1 instantiation per offered storage type, not a
cross-product, and that Storage is a compile-time if-constexpr axis with no
hidden runtime dispatch. The user interface stays fp64 (R passes double; the
downcast is internal at ingestion), so dbarts.h / the ABI is UNCHANGED - no
hash bump, no stan4bart lockstep. Bound the compile-time blowup by
instantiating narrowed storage only for the branches where it pays (the
continuous/gaussian residual path; latent families - probit/logistic/
ordinal/multinomial - keep fp64).

ANTI-PATTERN to avoid: a runtime precision TAG branched per element, or a
type-erased element accessor. That taxes the default path. Dispatch once,
monomorphize the core - exactly what the facade already does.

One caveat surfaced only once implementation started: "compile-time cost
only" understates the true cost. Storage threads through Forest<L>/Chain<L>/
Sampler<L> and every hand-written mixed-precision loop (rollTreeResidual is
where rounding enters, see 3b) - the default CODEGEN stays byte-identical
(enforced by a hot-loop assembly spot-check), but the SOURCE around it is
heavily rewritten. Budget for source churn, not just extra binary size.

## 3. The levers, by draw-impact class

Grounded in the storage inventory (2026-07-20 sweep of combiner.hpp/chain.hpp/
tree.hpp/data.hpp/model.hpp). KEY INVENTORY FACTS: (a) the predictor store is
ALREADY uint16 cutpoint codes ([[data.hpp#xint_t, codes]]) - splits run on
integer codes, raw doubles are quantized ONCE at setup, so the "fp32
predictor flips a split" hazard does NOT exist on the hot path and the
biggest n*p prize is already spent; (b) the LARGEST hot array is the gather
index (indexBuffer, size_t, n*trees) - bigger than the predictor codes; (c)
leafOf is already uint32, categorical masks uint64 bitsets, several flags
uint8 - no float anywhere today.

### 3a. PRESERVING (bitwise-identical draws; NO SBC, NO re-record; gated only
by the equivalence trio staying identical). DEFAULT-ON - a separate,
lower-risk TRACK 1 that lands first and independently of the fp32 work.

- **Index narrowing size_t -> uint32** (Forest::indexBuffer [[combiner.hpp#Forest::indexBuffer]],
  its Tree::indices alias [[tree.hpp#Tree::indices]], and VarianceForest::indexBuffer
  [[chain.hpp#VarianceForest::indexBuffer]]). THE TOP LEVER: footprint n*trees*8B (the single biggest hot
  array; ~1.6 GB at n=1e6/200 trees) -> halves to 4B, and it halves the
  streamed-index bytes inside the DOMINANT gather (microbench: index
  narrowing alone gives ~1.21x on the gather at n=1e6, g64->g32, and lifts the
  fp32 asymptote 1.33x->1.5x).

  PRESERVING - VERIFIED by the blind critique: the SIMD partition
  (partition_body.c) vectorizes the COLUMN VALUES fed to the code compare and
  never packs the index buffer into vectors (indices are only scalar
  subscripts), so its swap sequence / control flow is independent of index
  width -> a byte-identical permutation -> identical gather order ->
  bitwise-identical fp64 sums and draws. State save/restore does NOT
  serialize index width (ForestStateData/ChainStateData hold only FlatNodes +
  fp64; buffers are reconstructed by repartition). n < 2^31 always holds
  (predictor codes would be terabytes first); a guard throws above UINT32_MAX.

  BLAST RADIUS - the change is NOT a one-line typedef: misc_size_t is a single
  macro used for BOTH indices AND lengths, so there is no index typedef to
  flip. The change introduces a NEW misc_index_t DISTINCT from the length type
  and threads it through the whole misc.a index family: the two
  runtime-dispatched partition fn-pointers misc_partitionRange /
  misc_partitionIndices compiled across 5 ISA TUs (avx2/sse2/sse4_1/neon +
  scalar body), misc_partitionIndicesSparse, the 16 misc_computeIndexed*
  suffstat/mean/variance kernels AND their misc_mt_ variants,
  misc_sumIndexedVectorElements, misc_setIndexedVectorToConstant; plus
  Tree::indices ([[tree.hpp#Tree::indices]]), the 8 scalar C++ partition kernels, and
  SubtreeSnapshot.indexSegment ([[tree.hpp#SubtreeSnapshot::indexSegment]], std::vector<size_t> + a sizeof(size_t) memcpy).
  THE TRAP: preservation holds ONLY if each SIMD kernel is
  RETYPED, not REWRITTEN - a fresh uint32 SIMD partition with a different swap
  sequence changes the permutation and silently turns a "preserving" lever
  into a draw-changing one. GUARD: static_assert(sizeof(index_t) ==
  sizeof(misc_index_t)) pinning the C++ buffer type to the C kernel param
  (mirroring the existing [[tree.hpp#"static_assert(std::is_same_v<misc_xint_t, std::uint16_t>);"]]
  assert); that guard shipped with the lever and is live as
  [[tree.hpp#"static_assert(sizeof(index_t) == sizeof(misc_index_t));"]], six lines below the misc_xint_t assert it mirrors.
  ACCEPTANCE: the cross-ISA tests/cpp gate AND the bitwise
  equivalence trio must both stay identical - any drift means a kernel was
  rewritten, not retyped.

  LANDED (commit 2980229): retyped, not rewritten, across the whole misc.a
  index family described above; both gates re-run independently and the
  equivalence trio stayed BITWISE; ~400MB saved at n=5e5.

- **leafOf uint32 -> uint16** ([[combiner.hpp#Forest::leafOf]], n*trees). PRESERVING (node
  ids are exact) and a GENUINE gather input (mu[leafOf[i]] streamed every
  roll, so the reward is real), but CAPACITY-BOUNDED: a tree exceeding 65535
  nodes overflows. Shallow BART trees (base=.95/power=2) astronomically never
  approach this, but the fallback must be a real TYPE spec, not a hand-wave: a
  uniform std::vector<uint16> cannot host a per-tree fallback, so the guard is
  a GLOBAL refuse (throw/cap at the first tree that would exceed 65535 nodes),
  not a per-tree type union. uint8 is NOT safe (>255 nodes reachable). leafOf
  is internal scratch (not in ForestStateData) -> no state/ABI exposure.

  Considered as a Track 1 follow-on, gated the same PRESERVING way as the
  index change above - but ultimately DECLINED once the post-landing
  re-profile (section 6) showed the roll's cost is dominated by treeY, not
  leafOf, so the win would be negligible; see section 6 for the full
  disposition.

### 3b. CHANGING (feeds precision-sensitive arithmetic; SBC/coverage gate +
one re-record; OPT-IN only). The primary target.

- fp32 running residual (treeY). DESIGN = FORK A (a master+shadow Fork B was
  considered and is DEAD - see below). Make treeY itself fp32: the per-tree
  roll writes fp32, the per-tree gather (setNodeAverages -> computeLeafStats)
  reads fp32, and the proposal-time suffstat scans read fp32. Both wins (half
  the roll-write bytes, half the gather bytes + cache-residency) come from
  treeY being fp32. Leaf draw stays fp64 (sumWeights / sumWeightedResponse are
  fp64 accumulators fed by the float-loading suffstat kernels); the four
  misc_compute*SufficientStatisticsFast kernels get float-input variants
  (load float, accumulate double).

  WHY FORK B IS DEAD ([[chain.hpp#Chain::run]]): treeY is a running residual updated
  PER-TREE inside the backfit loop - `rollTreeResidual(forest, t, ...)`
  ([[chain.hpp#rollTreeResidual]]) rewrites it, then the gather
  ([[chain.hpp#Tree::setNodeAverages]]) reads tree t's freshly-rolled
  residual. A fp64-master + fp32-shadow scheme would have to refresh the
  shadow after every roll = a per-tree O(n) downcast, the same order as the
  gather it feeds -> no amortization, strictly worse than rolling fp32
  directly. So the fp32 store cannot avoid the incremental fp32 roll; this is
  why memory-wall-frontier.md's initial recommendation to prototype Fork B
  first did not survive detailed design (see that doc's section 8c).

  THE ROLL IS INCREMENTAL ([[chain.hpp#rollTreeResidual]]: t=0 recomputes
  resid = y - total + mu[leaf]; t>0 does resid += mu[leaf] - muPrev[leafPrev]),
  so fp32 rounding ACCUMULATES across trees (~200 updates/sweep/element) and,
  because it is never independently re-summed, across sweeps. This is the
  real and ONLY safety question for Track 2 (there is no "safe by
  construction").

  DRIFT CONTROL: an earlier claim that finalizeTotalFits gives a free
  per-sweep re-anchor turned out to be WRONG. finalizeTotalFits
  ([[chain.hpp#finalizeTotalFits]]) computes total = y - resid + mu[leaf] - i.e. FROM the drifted resid -
  so it PROPAGATES the fp32 drift into totalFits, not resets it; next sweep's
  t=0 roll re-derives resid from that drifted total. Nothing in the existing
  per-sweep machinery re-anchors. A genuine re-anchor needs an INDEPENDENT
  total = sum_t mu_t[leafOf_t[i]] (a fresh per-sweep O(n*trees) re-sum). That
  re-sum gathers into the L1-resident mu tables (not the DRAM residual), so it
  is far cheaper per element than the treeY gather (~10-15% of it), i.e.
  affordable - but it is NOT free and NOT the existing pass. It would bound
  CROSS-sweep drift, but the draw-relevant error is WITHIN-sweep and
  IRREDUCIBLE below per-tree cadence: the gather at tree t (setNodeAverages)
  sees ~t accumulated fp32 roundings since the last t=0 anchor, and no
  re-anchor coarser than per-tree removes it. Magnitude ~sqrt(200)*1e-7 ~ 1e-6
  relative, and it is ZERO-MEAN (round-to-nearest), and it is further
  suppressed by fp64 leaf-sum averaging (~1/sqrt(leaf size)); so the prior was
  that it sits below the MCMC noise floor - but that was a hypothesis for
  falsification to confirm, not a guarantee. So Track 2's viability rested
  ENTIRELY on a statistical gate, run BEFORE the full mixed-precision
  templating (large source churn): prototype fp32 treeY with NO re-anchor and
  measure (a) the drift max_i|treeY_fp32 - exact| growth over a long n=1e6
  run, and (b) SBC / interval coverage vs fp64. If SBC passes un-anchored,
  ship raw fp32 (cheapest); if drift/SBC fails, add the independent
  per-sweep re-sum re-anchor (bounds cross-sweep; re-benchmark the net win,
  since it eats ~10-15% of the gather) and re-test.

  This falsification was run (throwaway fp32-treeY prototype on the 2980229
  tip, round32(x) = (double)(float)x on every residual write - provably
  bit-identical to a true std::vector<float>+float-load/double-accumulate
  since float->double promotion is exact and every consumer already reduces
  in fp64). Flag-off stayed byte-deterministic; snapshot suite 186/186.
  Confirmed genuinely fp32-storing/fp64-reducing (drift 1e-15 off -> 1e-6 on,
  a 9-order gap). Results:
  - DRIFT: a slow RANDOM WALK (log-log slope ~0.5, technically unbounded
    because finalizeTotalFits propagates rather than re-anchors), but TINY -
    rms_rel ~ 7.9e-8*sqrt(sweeps); even 1e6 sweeps -> ~7.9e-5 rms / 7.5e-4 max
    relative to sd(y). Does NOT worsen with n (bigger leaves -> more
    fp64-averaging).
  - SBC/COVERAGE (the arbiter; fp32-vs-fp64 paired + MC-noise control, 30
    seeds, n in {1e3,1e4}): fp32 INDISTINGUISHABLE from fp64 - no detectable
    bias in test RMSE, 90% CI coverage of true f, or the sigma posterior (all
    p >> 0.05); the fp32-vs-fp64 divergence is 3-6x SMALLER than an RNG
    reseed. (Absolute BART under-coverage/small-n overfit are identical
    across arithmetic modes - orthogonal to fp32.)
  - RE-ANCHOR (built + measured as insurance, not shipped): the independent
    per-sweep mu[leafOf] re-sum flattens the walk (slope -> ~0, rms ~1e-6)
    and makes fp32 near-bit-identical to fp64 - but costs ~+20% of a sweep
    (naive, unfused), a large fraction of the fp32 gather win.
  - VERDICT: ship UN-ANCHORED (SBC passes, cheapest); keep the re-anchor as a
    documented, un-built fallback for >1e5-sweep chains or a stricter latent
    family.
  - CAVEAT: gaussian constant-leaf only; this SBC pass was a paired-comparison
    proxy, not full rank-uniformity SBC (the full-implementation landing gate
    ran the fuller check - see section 6). The prototype measured DRAWS
    faithfully but not the SPEED win (it used a double-storage round op, not
    real fp32 storage) - the real speedup is measured in section 6.

  Both wins (roll-write bytes, gather bytes/cache-residency) come from treeY
  being fp32; the actual implementation and its speed measurement are in
  section 6.

- fp32 SCRATCH/FITS BUNDLE (same opt-in flag, same SBC gate as treeY,
  gaussian/continuous path only): Forest::totalFits ([[combiner.hpp#Forest::totalFits]], n),
  test-fit accumulators totalTestFits/currTestFits (nTest), the
  non-constant-leaf treeFits slab (n*trees, dense-leaf models only), the
  variance-forest arrays for HBART (factorByTree n*treesVar,
  combinedVariance/meanResidual/divisor/treeResidual n each,
  [[chain.hpp#combinedVariance, meanResidual, divisor, treeResidual]]), and the gaussian working response yRescaled_ (model.hpp, n).
  LATENT-family working buffers (probit/logistic/ordinal/multinomial) stay
  fp64 by design (little value, more instantiations). This would extend the
  fp32 win to the streaming passes and the HBART weight channel; each piece
  is individually smaller than treeY but shares the flag and gate. This
  bundle was designed but never built - see section 6 for the measured NO-GO
  and why.

### 3c. CORRECTNESS-SENSITIVE (a discrete decision, not just rounding).
EXCLUDED from every track - flagged so no one bundles them as a transparent
memory tier.

- Predictor codes uint16 -> uint8 ([[data.hpp#codes, maxNumCutsRepresentable]]): the largest n*p array but
  ALREADY uint16; narrowing further lowers maxNumCutsRepresentable (~254),
  capping cutpoints and silently changing which splits exist -> different
  tree structure, pathological on high-cardinality columns. Only ever an
  explicit low-resolution MODELING option with its own arc, never a
  transparent tier.
- fp32 cutpoints ([[data.hpp#cutPoints]], tiny footprint) and fp32 ownedTestValues
  ([[data.hpp#ownedTestValues]], cold re-quantize source): fp32 flips near-tie quantization ->
  correctness-sensitive; not worth the risk, and cutpoints are not n-scaled.

Also out of scope, decided at design time rather than by later measurement:
bf16 residual storage (the aggressive tier beyond fp32; left for a future
arc) and fp32 for the latent-family working buffers (probit/logistic/
ordinal/multinomial - little value, more instantiations, so never pursued).

## 4. Opt-in surface (fork for VD; resolved)

- Option 1 (coarse): a single memory-tier flag on dbartsControl, e.g.
  `storage = c("double", "single")` (or a `low.memory` toggle) that flips a
  curated bundle. Matches how users think ("big-n, reduce memory"), not
  "narrow my index array". PRESERVING levers default-on regardless.
- Option 2 (fine): per-lever knobs. More control, more surface, more
  instantiations to compile and test.
- Chosen: the coarse flag for the CHANGING bundle (fp32 residual), with
  PRESERVING levers (uint32 index) default-on as a separate track. This is
  what shipped: opt-in `storage = "single"` on dbartsControl/bart2 (section 6).

## 5. Gates

- PRESERVING track: equivalence trio must stay BITWISE-IDENTICAL (no
  re-record); tests/cpp; sanitizer for new reachable code.
- CHANGING track (fp32 residual, opt-in): SBC / interval coverage vs the fp64
  sampler is the ARBITER (MCMC stationarity can amplify a per-sweep bias);
  same-machine A/B speed at n >= 1e6; the DEFAULT-path equivalence baselines
  stay identical (only the opt-in path shifts, so no default re-record); ASAN.
- The default fp64 instantiation must diff-clean against HEAD in generated
  assembly for the hot inner loop (spot-check) - the INVARIANT.

## 6. What landed

Before implementation, the design underwent an independent blind critique (per
the design-first process this doc committed to). It found two things worth
fixing - both already folded into section 3 above: Track 2's "free
re-anchor" claim was wrong (finalizeTotalFits propagates drift rather than
resetting it, so Track 2 rests on an explicit SBC falsification instead of an
architectural guarantee), and Track 1's index narrowing, while genuinely
bitwise-preserving, is not a one-line typedef (misc_size_t is shared with
lengths; a new misc_index_t has to thread through the whole misc.a index
family, retyped not rewritten, pinned by a sizeof static_assert). It also
CONFIRMED several properties held: the fp32 gather reduces in double (kernels
already accumulate double); no hidden runtime dispatch (Storage is a
compile-time if-constexpr axis; default codegen byte-identical); compile
blowup bounded to +1 instantiation per offered storage type, not a
cross-product; thread-safe (chain-parallel, serial gather kernels); warm-
start/state-safe (the state wire stays fp64); stan4bart/ABI unchanged;
ownership unchanged.

### Track 1: index narrowing - landed

Commit 2980229: index_t = uint32, retyped-not-rewritten across the misc.a
index family (section 3a). Gates re-run independently; equivalence trio
BITWISE; ~400MB saved at n=5e5.

### Track 2: real fp32 residual - landed and validated

Commit d384211: opt-in `storage="single"` on dbartsControl/bart2 -> real fp32
treeY (std::vector<ResidT>, ResidT=double default; +1 instantiation
SamplerFacade<ConstantGaussianLeaf,float> minted only at the
gaussian-constant-leaf branch) + four float-input suffstat kernels (load
float, accumulate double). Un-anchored, per the falsification result in 3b.
Errors clearly on probit/hetero/monotone/BCF/multinomial (not offered there).
Independently verified:

- INVARIANT (flag-off) BITWISE: equivalence 27/27 + bcf + multinomial ->
  default path provably byte-unchanged (the Storage axis defaulted to
  double).
- Genuine fp32 (not a silent double fallback): fp32-vs-fp64 posterior-mean
  divergence 1.9e-6, matching the ~1e-6 prediction from the falsification
  prototype.
- Statistically indistinguishable (real templated build, 10 seeds, n=1e4,
  200 trees): test-f RMSE Delta 5.7e-8 (p=.12), 90% coverage Delta 1e-5
  (p=.34) - a fuller pass than the paired-comparison proxy used at
  falsification time.
- codoc complete; ABI header untouched (no stan4bart lockstep).

SPEED A/B (n=1e6, single-thread, min-of-3, storage=single vs double):

    M1 Max (arm64):         1.17x (17% faster)  - matches the gather-share projection.
    x86 Ryzen 3700X (AVX2): 1.10x (10% faster)  - clean build, 0 warnings.

The x86 UNDERSHOOTS the ~1.24x projection because x86 pays a cvtss2sd
conversion tax on the CACHE-RESIDENT streaming roll (memory-wall-frontier.md
sec 8a: x86 fp32 streaming is ~2x slower cache-resident) that partly offsets
the gather win; M1 (compute-bound streaming) does not, so M1 matches
projection. Both results are POSITIVE and above the +/-2% floor. The win
should GROW at n>1e6 on x86 (residual spills the 16MB L3 -> streaming goes
DRAM-bound, where fp32 helps streaming too).

### v2 scope and leafOf: declined by measurement

The v1 A/B settles the fp32 scratch-and-fits bundle (section 3b) against itself.
The bundle (totalFits, test fits, dense treeFits, variance-forest arrays,
gaussian working response) is STREAMING, not gathered - and the microbench
plus the x86 v1 undershoot show fp32 STREAMING gains ~nothing (M1
compute-bound) to NEGATIVE (x86 cache-resident conversion tax). So fp32-ing
the streaming scratch adds conversion overhead with no gather benefit ->
neutral-to-regressive. The gather (treeY, v1) carries the entire win. DO NOT
build this bundle unless a re-profile at n>>1e6 shows a streaming pass that
is both DRAM-bound AND a large sweep share.

leafOf uint16 (the deferred PRESERVING Track-1 follow-on, section 3a) is
DECLINED: leafOf is only a minor input to the STREAMING roll (dominated by
the treeY read/write it sits beside, and the re-profile below puts roll at
~30% with treeY the bulk of it), so the speed gain is negligible; and it
cannot be done safely without a GLOBAL-REFUSE guard that makes the sampler
ERROR if any tree ever exceeds 65535 nodes (a uniform vector<uint16> cannot
host a per-tree fallback) - a surprising new failure mode for deep-tree
priors, not worth ~400MB. Reopen only if a concrete memory-constrained
large-n user needs the last halving AND accepts the node cap.

No other live storage lever remains; the arc is complete.

### Post-landing re-profile and multi-chain scaling

x86 rdtsc re-profile at n=1e6 (single-thread, 200 trees), buckets = gather
(computeLeafStats) / roll (rollTreeResidual) / other (moves/scans/draws):

    historical (2941808, pre-cuts):              gather 37.8%
    storage=double (uint32 idx, fp64 resid):      39.2 / 29.9 / 30.9   sweep 42.4s
    storage=single (uint32 idx, fp32 resid):      33.4 / 32.0 / 34.7   sweep 37.1s

The uint32 INDEX narrowing did NOT move the gather's share (39.2% vs 37.8%):
the gather is fp64-RESIDUAL-bandwidth-bound, not index-bound (the 4B index is
minor next to the 8B value), so it sped up the gather AND the other
index-touching passes ~proportionally, holding the ratio. Track 1's value is
the ~400MB saving + a proportional modest speedup, NOT a gather-concentrated
win (this corrects the earlier microbench-only framing, which suggested a
larger gather-specific effect). The fp32 RESIDUAL is what attacks the gather:
-28% gather cycles, share 39->33%. Under fp32 the sweep FLATTENS to
near-even thirds (33/32/35), with "other" (moves/scans/draws - heterogeneous,
no single lever) the largest slice.

MULTI-CHAIN A/B (n=1e6, fp32/fp64 speedup): the win GROWS with chains because
parallel BART shares memory bandwidth -

    x86 Ryzen (dual-channel DDR4):  1 thread 1.10x -> 4 chains/4 threads 1.30x
    M1 Max (huge BW):               1 thread 1.17x -> 4 chains/4 threads 1.19x

So fp32 pays MOST in the DEFAULT multi-chain configuration on
commodity/bandwidth-limited hardware (~30%), not the single-thread corner the
first A/B measured. Documented in NEWS/Rd: "storage='single' is most valuable
for multi-chain runs and large n."

WHY NO RULED-OUT LEVER REVIVES: the reductions FLATTENED the profile (fp32) or
held it (index) rather than exposing a new fat target - no single pass to
attack, and the incumbent got faster + smaller/more-cache-resident, RAISING
the bar for offload (GPU: the fp64-blocker dissolved, per memory-wall-
frontier.md sec 5a, but the sequential-backfit + launch killer stands and the
CPU baseline it must beat got faster) and SHRINKING the payoff of
traffic-cutting tricks (binning/atoms/block-fusion/NT-stores, memory-wall-
frontier.md sections 2/4). The one shifted item - within-chain parallelism
(blocked-jacobi-trees) - had its memory-wall bound LOOSENED ~30% by fp32, and
the 4-chain A/B PROVES parallel BART is bandwidth-bound on commodity HW (that
IS why fp32 gives 30% there), so fp32 was a prerequisite MITIGATION for it -
but its go/no-go rested on a separate statistical ESS/sec experiment,
unchanged by storage. That experiment has since concluded, independently of
this doc: both within-chain threading and blocked-jacobi-trees are NO-GO on
the real engine (see within-chain-threading.md and docs/plans/blocked-jacobi-
trees.md; not re-derived here). Net: no revival; fp32 is itself the
memory-wall answer, and it is broader-value than first reported (the
multi-chain default is where it pays most). See memory-wall-frontier.md
section 9 for the umbrella verdict.

## 7. Where the open questions landed

- RESOLVED: index narrowing stays bitwise (SIMD partition is index-width-
  independent; state does not serialize width). See 3a/6. LANDED 2980229.
- RESOLVED: Fork B dead (per-tree roll); Fork A is the design. See 3b.
- RESOLVED: un-anchored fp32 treeY drift stays far below the MCMC noise
  floor; SBC passes; re-anchor NOT needed (kept as a documented, un-built
  fallback). See 3b.
- RESOLVED: the real fp32-storage speed win was confirmed by the same-machine
  A/B at n>=1e6 (1.10x-1.30x depending on thread/chain count) and a fuller
  SBC pass on the templated implementation passed. See 6.
- RESOLVED: compile blowup bounded (+1 instantiation per offered storage
  type); no hidden runtime dispatch. Source churn was the real cost, budgeted
  for (section 2 caveat).
- RESOLVED: n-gating beyond opt-in is not needed - opt-in is manual
  (`storage="single"` on dbartsControl/bart2), never auto-selected, so
  small-n users who do not benefit simply do not opt in.

No open questions remain; the storage arc is complete.
