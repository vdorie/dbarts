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
  gather at n=1e6, g64->g32, and lifts the fp32 asymptote 1.33x->1.5x). PRESERVING:
  indices are exact integers and n < 2^31 always holds for real BART (the
  predictor codes alone would be terabytes first); a guard throws above UINT32_MAX.
  This is a plain `index_t = uint32_t` typedef swap across the partition/sort
  kernels (tree.hpp:530-706, the misc_size_t* two-pointer partitions) + the buffer
  decls - a straight DEFAULT change, NOT the opt-in template axis (uint32 is
  universally safe, so both precisions need not coexist). Validated by the
  equivalence trio staying BITWISE-IDENTICAL. Verify state save/restore does not
  serialize the index width (it should not - the index buffer is reconstructed
  scratch, not saved state).
- **leafOf uint32 -> uint16** (combiner.hpp:143, n*trees). PRESERVING (node ids
  are exact) but CAPACITY-BOUNDED: a tree exceeding 65535 nodes overflows. Shallow
  BART trees (base=.95/power=2 prior) essentially never approach this, but it is
  not structurally impossible, so it needs a guard (fall back / cap). uint8 is
  NOT safe (>255 nodes reachable). Halves another n*trees stream read every roll.
  CONSIDER for Track 1 if the guard is cheap; the critique weighs guard-vs-reward.

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

  THE ROLL IS INCREMENTAL (chain.hpp:894-896: "retires the previous tree's new
  fits and admits this tree's old ones"), so fp32 rounding ACCUMULATES across
  trees (~200 updates/sweep/element) and sweeps. This is the real SBC risk (not
  "safe by construction"). Drift control = a periodic fp64 RE-ANCHOR:
  recompute treeY exactly = forestY - sum_of_current_fits in fp64, store fp32,
  resetting drift; cadence (every sweep? every K?) trades accuracy vs cost. A
  per-sweep re-anchor is CHEAP because `finalizeTotalFits(forest, forestY)`
  (chain.hpp:965) ALREADY runs one exact fp64 total-fits pass per sweep for the
  latent/sigma updates - the re-anchor rides that existing pass rather than
  adding a new O(n*trees) one, so drift can be reset every sweep at ~no extra
  asymptotic cost. SBC decides whether the re-anchor is needed at all and how
  often. Within-sweep drift with a once-per-sweep anchor is
  ~sqrt(200)*1e-7 ~ 1e-6 relative, plausibly below the sampler noise floor;
  un-anchored multi-sweep drift is the thing to measure.
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

TRACK 1 (PRESERVING, default-on, no SBC, no re-record) - land FIRST and
independently; it de-risks and delivers the biggest footprint cut on its own:
- index_t = uint32 (the top lever; n*trees, halves the biggest hot array + helps
  the gather). Plain default typedef swap, gated by bitwise equivalence.
- leafOf uint16 CONSIDER (guarded; the critique decides if the guard earns it).

TRACK 2 (CHANGING, opt-in via a coarse `storage="single"` flag, SBC-gated, one
re-record of the opt-in path only) - land AFTER Track 1:
- fp32 treeY (FORK A, incremental fp32 roll + finalizeTotalFits-piggybacked
  re-anchor) + the four float-input suffstat kernels.
- fp32 scratch/fits bundle (totalFits, test fits, dense treeFits, variance-forest
  arrays, gaussian working response) behind the SAME flag.

OUT: predictor codes uint8 and fp32 cutpoints/test-raw (correctness-sensitive,
sec 3c); bf16 (the aggressive tier; future); latent-family fp32 (low value).

## 7. Open questions for the blind critique

- Does the uint32-index change actually stay bitwise on the equivalence trio, or
  does any code path serialize/hash the index width (state save/restore)?
- RESOLVED (Fork B is dead; treeY rolls per-tree - see 3b). The live question is
  the fp32-roll drift: does the incremental fp32 roll bias the posterior without
  a re-anchor, and if so what cadence clears SBC at least cost? Is a full fp64
  re-anchor per sweep too expensive (it is ~O(n*trees)), forcing a cheaper
  bound (e.g. re-anchor every K sweeps, or a fp64 running-sum-of-fits maintained
  cheaply)? This is the crux the prototype must settle.
- Compile-time / binary-size blowup: acceptable, or does it force restricting
  the storage axis further?
- Is n-gating needed beyond opt-in (the x86 small-n conversion loss), or does
  opt-in fully cover it?
