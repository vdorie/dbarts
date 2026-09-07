# Reading review: engine, bridge, support libs (b102e17c)

Read-only candidate list for the maintainer's judgement - nothing decided, nothing removed. Counts
re-derived in-tree; disagreements with review-lenses-memo.md F5 stated. Tags: `extension-point` (a
doc names what plugs in, or a consumer does), `compile-time-only` (no runtime indirection, no
reader cost), `sediment` (minted for a killed plan), `stale` (comment contradicts code), `diverge`
(two copies of one rule that answer differently).

## 0. Counts per class
- (1) abstractions with one user: 2, both type-erasure by necessity. Engine polymorphism is 5
  bases / 135 virtual declarations: ResponseModel 38 / 9 impls, ForestCombiner 30 / 2 + null
  default, ProgressSink 2 / 2, SamplerBase 61 / 1 source impl but 5 link-time instantiations,
  PredictorUpdateSession 4 / 1 (x5). Concepts 9, each satisfied by a shipped type and reaching a
  constraining template. Policies: Merge 2, Strategy 2, Columns 2. Template params: L (5
  arguments), ResidT (2).
- (2) duplication 12 sites, 5 DIVERGE, 3 of those correctness candidates. (3) dead or unreachable
  7. (4) bridge defense with no reachable trigger: 5 dead, 3 verified live. (6) runtime seams: 1
  violation of the architecture's own dispatch rule; the sweep hot path is otherwise virtual-free.
- (5) comments: 262 docs/ citations (R/ 109, bartcore 102, bridge 42, support libs 4, dbarts.h 2,
  Makevars 3) over 29 distinct paths - all 29 files exist, 38 anchors checked, 0 stale. 20
  narration sites, worst 15 below.

## 1. REFUTATIONS the memo's F5 rests on - read before deleting anything
- R1. `SamplerBase` is NOT "36 pure virtuals with one implementation". 61 virtual declarations
  ([src/bartcore/facade.hpp:139-401](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L139-L401)); the one source-level impl `SamplerFacade<L, ResidT>` ([src/bartcore/facade.hpp:405](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L405)) is a
  TEMPLATE instantiated five ways: ConstantGaussianLeaf/double, ConstantGaussianLeaf/float ([src/bartcore/facade.hpp:725](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L725)),
  MonotoneConstantGaussianLeaf ([src/bartcore/facade.hpp:793](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L793)), GPGaussianLeaf ([src/bartcore/facade.hpp:813](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L813), [src/bartcore/facade.hpp:846](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L846)), LinearGaussianLeaf ([src/bartcore/facade.hpp:817](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L817),
  [src/bartcore/facade.hpp:850](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L850)). The vtable is what lets the bridge hold `unique_ptr<SamplerBase>` without knowing L.
  `extension-point`, high.
- R2. `ResidT` has TWO instantiations; `float` is reachable from R (`control@storage == "single"`
  -> [src/R_interface_bartcore.cpp:398-404](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L398-L404) -> options.fp32Residual -> [src/bartcore/facade.hpp:725](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L725)).
  [src/bartcore/chain.hpp:4463-4480](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L4463-L4480) documents why both spellings exist - codegen, not taste.
  `compile-time-only`, high.
- R3. All nine leaf concepts are live. VectorLeafModel ([src/bartcore/model.hpp:60](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L60)) and FunctionLeafModel ([src/bartcore/model.hpp:85](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L85))
  head no template but are disjuncts of IntegrableLeafModel ([src/bartcore/model.hpp:98](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L98)), which heads
  chain/sampler/facade/combiner. ScaleLeafModel ([src/bartcore/model.hpp:111](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L111)) is the second disjunct of
  MoveScorableLeafModel ([src/bartcore/model.hpp:126](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L126)) and IS exercised: [src/bartcore/chain.hpp:4222](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L4222) instantiates moves.hpp with
  ConstantVarianceLeaf. `compile-time-only`, high.
- R4. The flat C API is fully exercised: of 62 symbols in inst/include/dbarts/dbarts.h every
  function is called by stan4bart (bartcore branch, 28) or by inst/tinytest/capi/consumer.c, and
  the 10 with neither are typedef/enum names and the `dbarts_stub_*` macro machinery. No dead
  exports.
- R5. Per-observation loops are monomorphic. The sweep's virtual boundaries are per-sweep
  ([src/bartcore/chain.hpp:1511-1521](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1511-L1521), [src/bartcore/chain.hpp:1358](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1358)) and per-sweep-per-forest ([src/bartcore/chain.hpp:1407-1408](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1407-L1408)); nothing virtual sits inside
  the per-tree loop ([src/bartcore/chain.hpp:1428-1503](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1428-L1503)) or below, exactly what [docs/design/core-generalization.md:70-76](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/core-generalization.md#L70-L76)'s dispatch
  table promises. The `misc_*` function-pointer kernels are per-node-op, sanctioned by
  [docs/design/kernel-vocabulary.md:12-16](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/kernel-vocabulary.md#L12-L16).
- R6. Comment rot is NOT in the docs/ citations: all 29 cited design docs exist and 38 sampled
  anchors resolve to what the citing comment claims. Only blemish: [src/R_interface_bartcore.cpp:6285](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6285)
  cites `public-surface.md 2`, whose TITLE is about factor ingestion though its body carries the
  registry rule past line 155.

## 2. DUPLICATION
DIVERGE first - a divergence outranks a large exact clone.

- D1. `alpha == 0.0` short-circuit exists only in the scalar kernels: [src/misc/linearAlgebra.c:67](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra.c#L67) and [src/misc/linearAlgebra.c:85](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra.c#L85)
  guard `if (length == 0 || alpha == 0.0) return;`; the same two functions guard on length alone
  at [src/misc/linearAlgebra_sse2.c:42](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra_sse2.c#L42), [src/misc/linearAlgebra_sse2.c:60](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra_sse2.c#L60), [src/misc/linearAlgebra_avx.c:53](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra_avx.c#L53), [src/misc/linearAlgebra_avx.c:75](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra_avx.c#L75), [src/misc/linearAlgebra_neon.c:150](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra_neon.c#L150), [src/misc/linearAlgebra_neon.c:220](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/linearAlgebra_neon.c#L220). The C fallback writes
  nothing where SIMD dispatch writes, in a library whose stated invariant is within-host bitwise
  identity across dispatch. Live callers ([src/bartcore/model.hpp:2833](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L2833), [src/bartcore/model.hpp:2912](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L2912), [src/bartcore/model.hpp:2968](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L2968), [src/bartcore/model.hpp:3033](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L3033)) pass `-min_`, i.e.
  `-0.0` when min_ is 0; `x += -0.0` is a round-to-nearest no-op for every finite x, so observable
  difference is confined to NaN payloads. `diverge`; high on the fact, low on impact.
- D2. Scale-leaf positivity applied on 3 of 4 state paths. [src/bartcore/chain.hpp:3301](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3301) (live variance trees)
  and [src/bartcore/chain.hpp:3334](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3334) (saved buffer) refuse a non-positive variance leaf, and [src/bartcore/sampler.hpp:919-923](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L919-L923) applies it
  to a warm start's SLOT-sourced trees with a comment naming the reason ("the buffer is
  hand-buildable and a rebuild scatters the leaf straight into a divisor"). The LIVE-sourced arm,
  [src/bartcore/sampler.hpp:891](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L891) (`dst.varianceTrees = src.varianceTrees;`), applies nothing and neither donor
  parser checks ([src/R_interface_bartcore.cpp:7233-7248](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L7233-L7248)); `rebuildVarianceForest`
  ([src/bartcore/chain.hpp:4356-4368](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L4356-L4368)) then scatters it into a divisor. Reachable from a `.Call` with a
  hand-built state, not from R's own state objects. High. CORRECTNESS CANDIDATE.
- D3. `applyNewData` ([src/bartcore/chain.hpp:2453](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2453)) and `recoverTreeParameters` ([src/bartcore/chain.hpp:2421](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2421)) take `forests_[0]` only
  where siblings loop `forests_` ([src/bartcore/chain.hpp:2329](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2329), [src/bartcore/chain.hpp:2401](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2401), [src/bartcore/chain.hpp:2526](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2526), [src/bartcore/chain.hpp:2540](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2540)). A whole-data replacement on a
  BCF/multinomial chain would leave forests 1..K-1 on the old grid; guarded two layers away by
  `refuseMultiForestMutation` ([src/R_interface_bartcore.cpp:4632](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L4632)), with no engine-level assert. High.
  CORRECTNESS CANDIDATE.
- D4. One enumerator, two opposite meanings, same file. `augmentationFamily`
  ([src/R_interface_bartcore.cpp:6151](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6151)) rejects "gaussian" and maps "student" onto
  `ResponseFamily::gaussian`; `drawAugmentationLaws` reads `case RF::gaussian:` as "the Student-t
  scale mixer" ([src/R_interface_bartcore.cpp:6225](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6225)), while `resolveFamily` ([src/R_interface_bartcore.cpp:1582](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L1582)) maps "gaussian" to the same enumerator
  meaning the Gaussian law and never accepts "student". `diverge`/`sediment`, high; a dedicated
  token costs nothing.
- D5. The shipped-BCF glue draw and the general amplitude draw disagree on the scale-mixture
  refresh: [src/bartcore/combiner.hpp:1171](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L1171) refreshes only `prior[0].variance`, [src/bartcore/combiner.hpp:1213](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L1213) refreshes every forest with
  `halfCauchyScale > 0`. The selector at [src/bartcore/combiner.hpp:986](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L986) routes on `shippedShape()` ([src/bartcore/combiner.hpp:1474](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L1474)), which tests
  basis widths and canonicality only, NOT the half-Cauchy flag - so a two-forest spec with
  canonical bases and a per-forest prior scale (admitted at [src/R_interface_bartcore.cpp:2178-2180](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2178-L2180))
  gets a fixed-variance prior where the general path would sample it. The duplication is
  deliberate and well documented ([src/R_interface_bartcore.cpp:940-984](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L940-L984), with a stated deletion trigger); the PREDICATE's
  coverage is the finding. Medium-high. CORRECTNESS CANDIDATE.

AGREE-but-restated (reading cost and drift risk, no defect today):

- D6. Six augmentation laws restated in the bridge: [src/R_interface_bartcore.cpp:6145-6270](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6145-L6270) (126 lines)
  re-implements Probit/Ordinal/AFT/Logistic/NB/T response draws ([src/bartcore/model.hpp:3085](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L3085), [src/bartcore/model.hpp:3220](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L3220), [src/bartcore/model.hpp:3548](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L3548) and
  the NB/T drawers). [src/R_interface_bartcore.cpp:6146-6149](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6146-L6149) declares it deliberate (a different generator, citing
  r-c-division.md); nothing enforces lockstep - no shared kernel, no cross-check test named at the
  site. High.
- D7. `linearAlgebra_sse2.c` and `_avx.c` contain ZERO intrinsics (`grep -c '_mm'` = 0 for both;
  `_neon` has 89) - they are the `_c` bodies with a different unroll factor (4 vs 8), relying on
  auto-vectorization. With `_neon` (303 lines, real intrinsics) and linearAlgebra.c (286) that is
  five routines written four times, while the sibling partition family shares one 350-line
  `partition_body.c` behind 28-line wrappers, an idiom used nowhere else in the directory.
  [src/misc/moments.c:251-344](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/moments.c#L251-L344) vs [src/misc/moments.c:354-447](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/moments.c#L354-L447) repeats it: four fp32 suffstat kernels are exact clones of the
  fp64 four, 94 lines kept in lockstep by hand ([src/misc/moments.c:349-353](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/moments.c#L349-L353) makes the summation order a correctness
  contract). High.
- D8. Five "build a scratch tree from flat, check containment" walks ([src/bartcore/chain.hpp:3199-3230](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3199-L3230),
  [src/bartcore/chain.hpp:3303-3324](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3303-L3324), [src/bartcore/chain.hpp:3365-3382](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3365-L3382), [src/bartcore/chain.hpp:3409-3428](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3409-L3428), [src/bartcore/chain.hpp:3441-3450](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3441-L3450)): same eight lines, but the first two `return
  false` on a failed build and the last three `continue`, and the mask selection is spelled two
  ways ([src/bartcore/chain.hpp:3206-3210](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3206-L3210) vs [src/bartcore/chain.hpp:3411-3415](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3411-L3415)). Each comment names the others as its mirror; no doc covers it.
  High.
- D9. Eight near-parallel saved-vs-live replay functions ([src/bartcore/chain.hpp:2799](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2799)/[src/bartcore/chain.hpp:2825](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2825), [src/bartcore/chain.hpp:2865](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2865)/[src/bartcore/chain.hpp:2898](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2898),
  [src/bartcore/chain.hpp:2938](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2938)/[src/bartcore/chain.hpp:2962](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2962), [src/bartcore/chain.hpp:904](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L904)/[src/bartcore/chain.hpp:932](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L932)). `predictVariance` ([src/bartcore/chain.hpp:934](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L934)) guards `if (!varianceForest_) return;`; its
  saved twin ([src/bartcore/chain.hpp:906](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L906)) dereferences immediately, unreachable today only because the bridge gates on
  `shape.hasVarianceForest` ([src/R_interface_bartcore.cpp:5732](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L5732)). High.
- D10. Leaf-shape flatten switch three times: [src/bartcore/chain.hpp:2683-2716](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2683-L2716), [src/bartcore/chain.hpp:2652-2677](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2652-L2677), [src/bartcore/chain.hpp:3004-3040](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L3004-L3040) - the
  first two a documented cache-vs-recompute pair, getState's a different format. High on
  duplication, medium on intent.

REFUTED - the memo's "the variance forest re-implements the mean forest's lifecycle as 12 parallel
members though ConstantVarianceLeaf already satisfies ScaleLeafModel". The parallelism is real and
larger (15 shadowed members, [src/bartcore/chain.hpp:407-484](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L407-L484) vs [src/bartcore/combiner.hpp:138-231](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L138-L231), plus ~15 paired Chain
methods) but the premise does not carry: [src/bartcore/model.hpp:322-324](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L322-L324) static_asserts
`!IntegrableLeafModel<ConstantVarianceLeaf>`, and Forest ([src/bartcore/combiner.hpp:137](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L137)) and Chain
([src/bartcore/chain.hpp:532](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L532)) are both constrained on IntegrableLeafModel. What CAN be shared is -
MoveScorableLeafModel admits it, `logLikelihoodForBranch` scores it unchanged. The blocker is the
combination law (Forest's treeFits/totalFits/treeY are additive; the variance forest is
multiplicative, [src/bartcore/chain.hpp:400-407](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L400-L407)), logged as open debt in [docs/design/forest-combiner.md:207-211](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/forest-combiner.md#L207-L211), [docs/design/forest-combiner.md:252](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/forest-combiner.md#L252) and
[docs/design/heteroscedastic.md:59-62](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/heteroscedastic.md#L59-L62), [docs/design/heteroscedastic.md:586-590](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/heteroscedastic.md#L586-L590). `extension-point`, high.

## 3. DEAD OR UNREACHABLE
- U1. `Tree::rightChildOf` ([src/bartcore/tree.hpp:366](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/tree.hpp#L366)) and `Sampler::setCurrentSampleNum` ([src/bartcore/sampler.hpp:485](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L485)):
  zero references in src/, tests/cpp, benchmarks/, R/, inst/. `sediment`, high.
- U2. [src/bartcore/chain.hpp:766](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L766) `default:` on a 3-arm family switch. `createAmplitudeSampler`
  ([src/bartcore/facade.hpp:869-873](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L869-L873)) refuses every family but gaussian/probit/logistic first and the comment at
  [src/bartcore/facade.hpp:750-755](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L750-L755) says so. Unreachable today; the cost is that `case gaussian: default:` suppresses
  -Wswitch, so a 7th enumerator would silently fit Gaussian here. `sediment` - defense that
  disables the compiler's own check.
- U3. Three further -Wswitch suppressions on `ResponseFamily`: [src/bartcore/chain.hpp:5033](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L5033) (folding
  gaussian/aft/ordinal/nbinom), [src/R_interface_bartcore.cpp:2298](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2298), [src/R_interface_bartcore.cpp:2842](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2842). Of 7 family ladders only 3
  are exhaustive. Deleting each `default:` and reading the compiler is the cheapest check on this
  list.
- U4. `refusedAmplitudeFamilyReason` ([src/R_interface_bartcore.cpp:2268-2284](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2268-L2284)): exhaustive switch,
  every arm returns, so `return "this response family";` ([src/R_interface_bartcore.cpp:2283](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2283)) is unreachable - but required to
  silence -Wreturn-type. `compile-time-only`; touch only alongside U3.
- U5. `misc_simd_getMaxSIMDInstructionSet` ([src/misc/simd.c:142-183](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/simd.c#L142-L183)) detects SSE, SSE3, SSSE3,
  SSE4_2, AVX512F, AVX512VL, AVX512BW; the dispatcher ([src/misc/simd.c:284-340](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/simd.c#L284-L340)) tests only AVX2, SSE4_1, SSE2,
  AVX, NEON. Seven levels select nothing, and the AVX512 block ([src/misc/simd.c:170-183](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/simd.c#L170-L183), with its own second
  `__cpuidex`) computes state no branch reads. One-time cost; inherited from pre-bartcore misc
  (main carries simd.c at 9b0ae65b) but 125 of its lines changed here. `extension-point` if an
  AVX512 kernel is planned, else `sediment`.
- U6. [src/misc/partition_body.c:106](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L106), [src/misc/partition_body.c:123](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L123), [src/misc/partition_body.c:149](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L149) - a commented-out PARTITION_RANGE if/else/endif
  triple wrapping two abandoned NEON load strategies ([src/misc/partition_body.c:107-123](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L107-L123), inside a block comment): ~18 dead
  lines, the same shape as `#if 0`, which the memo's sec 1.0 reported as zero occurrences.
- U7. `XINT_TYPE` width generality has no user and a wrong-answer failure mode.
  [docs/design/kernel-vocabulary.md:24-27](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/kernel-vocabulary.md#L24-L27) says the code type is "configure-selected via `--with-xint-size`";
  [configure.ac:21](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/configure.ac#L21) hard-wires `uint16_t` with no such option and no width-suffixed kernel exists.
  The `#ifndef XINT_TYPE` guard (src/include/misc/types.h.in:6) still lets a CPPFLAGS define
  override it, and every SIMD partition kernel hard-codes `epi16`/`u16` ([src/misc/partition_body.c:11](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L11), [src/misc/partition_body.c:73](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L73),
  [src/misc/partition_body.c:151](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/misc/partition_body.c#L151)) with no `static_assert(sizeof(misc_xint_t) == 2)` anywhere. `sediment`, high.

## 4. DEFENSIVE CODE WITH NO REACHABLE TRIGGER (bridge)
- V1. The state-format compatibility window is empty by construction: `stateFormatVersion` and
  `minReadableStateFormatVersion` are both 3 ([src/R_interface_bartcore.cpp:6312](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6312), [src/R_interface_bartcore.cpp:6321](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6321)), and the
  comment at [src/R_interface_bartcore.cpp:6310-6311](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6310-L6311) concedes "no release ever shipped format 3". The two floor checks ([src/R_interface_bartcore.cpp:6658](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6658),
  [src/R_interface_bartcore.cpp:7282](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L7282)) can only fire on a state with no version attribute (reads 0). Forward-looking; the
  registry rule at [src/R_interface_bartcore.cpp:6285-6298](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6285-L6298) is the durable part. `extension-point`, high.
- V2. A check byte-identical to its sole caller's: `createMultinomialCountsHolder` re-runs `if
  (!Rf_isInteger(countsExpr) || Rf_xlength(dimsExpr) != 2)` ([src/R_interface_bartcore.cpp:3379-3380](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3379-L3380)) which
  `createMultinomialDataHolder` already ran at [src/R_interface_bartcore.cpp:3457-3458](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3457-L3458), with only `parseMultinomialData` (which
  never touches countsExpr) between. The helper has C++ linkage only, so it is not a flat-C
  backstop. `dead-defense`, high. Its surrounding PROTECT is separately annotated ([src/R_interface_bartcore.cpp:3374-3377](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3374-L3377)) as
  deliberate analyzer bait - decide the two separately.
- V3. Author-acknowledged unreachable: `resolveCscCategoricalReferences`'s `if (source >=
  numSparseColumns) Rf_error("%s", boundMessage)` ([src/R_interface_bartcore.cpp:628-629](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L628-L629)). Its own doc comment at [src/R_interface_bartcore.cpp:609-610](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L609-L610) says
  so, and `mapColumnSources` ([src/R_interface_bartcore.cpp:570](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L570)) bounds every CSC entry at all three call sites.
  `dead-defense`, high; the odd part is that the comment names the redundancy and the check stays.
- V4. PROTECT convention applied inconsistently inside one function. Four pairs guard an attribute
  of the already-rooted `stateExpr`/`donorStateExpr` argument with nothing allocating in the
  window: [src/R_interface_bartcore.cpp:6649-6653](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6649-L6653) (formatVersion in setState), [src/R_interface_bartcore.cpp:7277-7281](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L7277-L7281) (its copy in installForests),
  [src/R_interface_bartcore.cpp:6731-6740](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6731-L6740) (sampleNum), [src/R_interface_bartcore.cpp:6742-6753](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6742-L6753) (recordedDraws). In the SAME functions, sibling reads of the
  identical pattern skip PROTECT ([src/R_interface_bartcore.cpp:6699-6700](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6699-L6700), [src/R_interface_bartcore.cpp:6711](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6711), [src/R_interface_bartcore.cpp:7079](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L7079), [src/R_interface_bartcore.cpp:7094-7095](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L7094-L7095)). None is among the four
  sites docs/plans/release-candidate-review.md's rchk note ([src/R_interface_bartcore.cpp:1138-1173](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L1138-L1173), [src/R_interface_bartcore.cpp:1300-1330](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L1300-L1330)) records, and
  all are distinct from the two `setState` rchk BAILOUTS the note calls already-balanced. The
  finding is the inconsistency, not the PROTECTs. Medium.
- V5. `rbart_getFitted`'s two PROTECTs ([src/R_interface_rbart.cpp:16-17](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_rbart.cpp#L16-L17)): both dims are reduced to raw
  `int*` at [src/R_interface_rbart.cpp:19-20](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_rbart.cpp#L19-L20) and last read at [src/R_interface_rbart.cpp:40-42](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_rbart.cpp#L40-L42), and the only allocation is `rc_newReal(n)` at [src/R_interface_rbart.cpp:44](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_rbart.cpp#L44).
  Untouched by the rchk commit. Medium - defensible only if `rc_getDims` can itself allocate.
- V6. NOT dead, keep. `setState`'s three "already-non-null" guards ([src/R_interface_bartcore.cpp:6799](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6799), [src/R_interface_bartcore.cpp:6813](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6813), [src/R_interface_bartcore.cpp:6830](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6830)) make the
  null branch in `readFunctionTreeParams` ([src/R_interface_bartcore.cpp:5490](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L5490)), `readTreeParams` ([src/R_interface_bartcore.cpp:5469](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L5469)) and `readTreeMasks`
  ([src/R_interface_bartcore.cpp:5438](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L5438)) unreachable on setState's path, but the same helpers are called unguarded from
  `readWarmStartState` ([src/R_interface_bartcore.cpp:7144-7164](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L7144-L7164)).
- V7. NOT dead, keep - this closes the "R already refused it" class. No C-side check on the .Call
  path is provably unreachable: every `dbarts_bartcore_*` symbol is an ordinary DL_FUNC in
  `R_callMethods[]` ([src/R_interface.cpp:180-259](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface.cpp#L180-L259)), reachable by any R code holding the internal `C_`
  name; two sites say so themselves ([src/R_interface_bartcore.cpp:3987-3993](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3987-L3993), [src/R_interface_bartcore.cpp:4878-4882](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L4878-L4882)). The S4 route confirms it:
  `dbartsControl`'s `setValidity` ([R/A_class.R:272-383](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/A_class.R#L272-L383)) is not re-run on plain `@<-` mutation and
  `setControl` ([R/dbarts.R:1155-1206](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/dbarts.R#L1155-L1206)) never calls `validObject()`. ~13 more sites carry the same
  explicit "backstop" annotation.
- V8. rc constraint API re-derived: 42 constrained calls (35 `rc_getInt`/`rc_getDouble`/
  `rc_getBool` carrying RC_LENGTH or RC_VALUE, 7 `rc_assert*Constraints`, 0 bare `rc_get*0`)
  against 111 unconstrained `rc_getListElement` fetches and 168 open-coded `Rf_is<Type>` checks in
  R_interface_bartcore.cpp (3 in C_interface.cpp, deliberately
  - [src/C_interface.cpp:3-5](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/C_interface.cpp#L3-L5)). The memo's "7" counted only `rc_assert*Constraints`. Of 50 `dbarts_bartcore_*` .Call
    entries (59 rows total), 20 call one of the 14 shared `refuse*`/`validate*` helpers directly;
    the rest go through `createHolder`/`parseData`/`parseModel`. 304 `Rf_error` sites in the
    bridge, 42 in the C API.

## 5. COMMENTS
- C1. THREE `stale` "Not yet exposed through the R surface" claims - all three features ARE
  reachable from R today: [src/bartcore/chain.hpp:88](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L88) (monotone) vs [R/model.R:99](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/model.R#L99), [R/model.R:526-531](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/model.R#L526-L531); [src/bartcore/chain.hpp:141](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L141)
  (interaction constraints) vs [R/spec.R:415-416](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/spec.R#L415-L416) and [src/R_interface_bartcore.cpp:1254](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L1254); [src/bartcore/chain.hpp:166](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L166)
  (variance forest) vs [R/spec.R:532](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/spec.R#L532), [R/dbarts.R:816](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/dbarts.R#L816), [src/R_interface_bartcore.cpp:2017](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2017). A reader who
  trusts these will not look for the R-side refusals that guard them. High. Best comment finding
  here.
- C2. [src/R_interface_bartcore.cpp:6272-6311](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6272-L6311), a 40-line block comment of which ~28 lines narrate three
  pre-release format iterations that never shipped. It opens "The shipped format (version 2)"
  while `stateFormatVersion = 3` sits 40 lines below - a contradiction inside one comment. The
  load-bearing part is the registry rule ([src/R_interface_bartcore.cpp:6285-6298](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6285-L6298), ~13 lines). `stale` + narration, high. Worst
  site by volume.
- C3. The other worst narration sites (criterion: the payload is a comparison to code that no
  longer exists and cannot be seen), most-useless first - [src/bartcore/chain.hpp:1041](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1041) ("the family conjunct
  that used to stand beside it was there because setResponse handed forest 0's bare totals as
  though combined"); [R/generics.R:1154](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/generics.R#L1154) and [R/generics.R:1287](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/generics.R#L1287) ("the same keepTrees gate a deleted `$bc` field
  used to"); [R/utility.R:119](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/utility.R#L119) ("gone now that the rename has landed"); [R/utility.R:49](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/utility.R#L49) ("silent
  before bart2 could forward a resid.prior object at all"); [src/R_interface_bartcore.cpp:2766](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2766) (quotes
  an old, replaced error message); [src/C_interface.cpp:624](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/C_interface.cpp#L624) ("which this entry used to drop on the
  floor"); [R/data.R:453](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/data.R#L453) (three appeals to an invisible prior implementation in one comment);
  [R/data.R:342](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/data.R#L342) ("this used to be a function evaluated in the caller's frame"); [R/spec.R:490](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/spec.R#L490),
  [R/model.R:1759](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/model.R#L1759) and [R/spec.R:652](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/spec.R#L652) (removed flat formals and a vanished literal); [R/bart.R:769](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/bart.R#L769) and
  [R/dbarts.R:627](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/R/dbarts.R#L627) (fixed-bug narration). 20 sites total, 14 in R/ and 6 in C++; src/misc,
  src/external, src/rc and src/include are narration-free.
- C4. Forward-looking PLAN narration (a plan is not a constraint): [src/bartcore/model.hpp:2589-2594](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L2589-L2594) ("v1 ships
  the exact integer envelope... A later real-shape mode routes a fractional b"), [src/bartcore/model.hpp:2774](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/model.hpp#L2774),
  [src/bartcore/facade.hpp:784](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L784) ("v1 keeps the mean leaf constant"), [src/R_interface_bartcore.cpp:3008](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3008), [src/R_interface_bartcore.cpp:3014](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3014), [src/R_interface_bartcore.cpp:3030](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L3030).
  Six sites; each also states a live refusal, so only the forward half is narration. Medium.
- C5. 262 docs/ citations no installed user can follow (`.Rbuildignore` has `^docs$`). The
  decision is not "are they rotten" - they are not - but whether a shipped comment should cite a
  stripped path.

## 6. SEAMS THAT COST AT RUNTIME
- S1. PER-OBSERVATION VIRTUAL DISPATCH, contradicting the architecture's own rule.
  `PredictorUpdateSession` ([src/bartcore/sampler.hpp:89-100](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L89-L100)) declares `observationWouldRemainValid(i)` and
  `commitObservation(i)`; `updatePredictorPerObservationJointly` ([src/bartcore/facade.hpp:694-703](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L694-L703)) calls both
  inside `for (i = 0; i < numObservations; ++i)`, once per sampler per observation.
  [docs/design/core-generalization.md:69-76](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/core-generalization.md#L69-L76) states "nothing dispatches per observation" and "Per obs | none:
  monomorphic loops/kernels". The erasure is NOT removable by templating - the joint sweep takes
  `SamplerBase* const*` over samplers of possibly different L, its whole purpose
  (R/updatePredictorPerObservationJointly.R; the bairrtt consumer) - though the single-sampler
  path ([src/bartcore/facade.hpp:217](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L217)) pays it without needing it. Frequency PER-OBSERVATION, 2 virtual calls x
  numSamplers x n; the work inside is a per-tree descent so the ratio is likely fine, but either
  the code or the doc's absolute rule should move. `extension-point`, high; not benchmarked.
- S2. `std::function` on the run path: pollInterrupt/shouldCancel ([src/bartcore/facade.hpp:147](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L147), [src/bartcore/sampler.hpp:274](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L274),
  [src/bartcore/sampler.hpp:362-370](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L362-L370), [src/bartcore/chain.hpp:1345](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1345)) and SweepCallback ([src/bartcore/chain.hpp:396](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L396)), both called once per sweep
  ([src/bartcore/chain.hpp:1358](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1358), [src/bartcore/chain.hpp:1362](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1362)). PER-SWEEP, negligible; listed to close the class.
- S3/S4. `combiner_->` virtuals: 44 sites in chain.hpp, the two inside the sweep's forest loop
  being drawForestGlue/formForestResponse ([src/bartcore/chain.hpp:1407-1408](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1407-L1408), prior-sampling copy at [src/bartcore/chain.hpp:1973-1974](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1973-L1974)) -
  PER-SWEEP-PER-FOREST, two impls plus a null default. `response_->` virtuals in `run`
  ([src/bartcore/chain.hpp:1511-1521](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L1511-L1521)) - PER-SWEEP, nine impls.

## 7. LOWER-VALUE, FOR COMPLETENESS
- L1. Fourteen `*ForTesting` accessors on the shipped engine that no production code calls (all
  reached from tests/cpp, three by exactly one assertion): they widen the engine's member surface
  for testability only. `extension-point` (test seam).
- L2. 113 named helpers appear exactly twice in the 30488-line shipped C++ corpus (definition plus
  one call site); 84 are never referenced from tests/cpp. Re-derived; the memo says 144/103 over
  the engine alone. Decomposition, not sediment - the 84 with no test reach belongs to lens F2.
- L3. core-generalization.md, the arbiter document, names two extension points no shipped code
  carries: `SplitSelector` ([src/bartcore/chain.hpp:129-133](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L129-L133)) and `MoveStrategy` ([src/bartcore/chain.hpp:118-126](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L118-L126), template rule at [src/bartcore/chain.hpp:86](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L86)), zero
  occurrences of either in src/. It cannot justify a seam by naming it alone; check the code.

## 8. THE TEN TO DECIDE FIRST
Ranked by maintainer time saved per decision, not by size.

1. C1 `stale` - three "Not yet exposed through the R surface" comments ([src/bartcore/chain.hpp:88](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L88), [src/bartcore/chain.hpp:141](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L141), [src/bartcore/chain.hpp:166](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L166))
   on features R reaches today. One-line fix, misleading now.
2. D2 `diverge` - variance-leaf positivity checked on 3 of 4 state paths ([src/bartcore/sampler.hpp:891](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L891) vs
   [src/bartcore/sampler.hpp:919-923](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/sampler.hpp#L919-L923)). Correctness candidate.
3. D3 `diverge` - `forests_[0]` hardcoded in applyNewData ([src/bartcore/chain.hpp:2453](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2453)) and
   recoverTreeParameters ([src/bartcore/chain.hpp:2421](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L2421)) where siblings loop. Correctness candidate.
4. D5 `diverge` - `shippedShape()` ([src/bartcore/combiner.hpp:1474](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/combiner.hpp#L1474)) routes on basis shape but not the
   half-Cauchy flag, so two amplitude specs get two different models.
5. S1 `extension-point` - per-observation virtual dispatch ([src/bartcore/facade.hpp:694-703](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/facade.hpp#L694-L703)) against
   [docs/design/core-generalization.md:69-76](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/docs/design/core-generalization.md#L69-L76)'s "nothing dispatches per observation"; the code or the doc's
   absolute rule moves.
6. D4 `diverge` - `ResponseFamily::gaussian` means "Student-t" in the augmentation surface
   ([src/R_interface_bartcore.cpp:6225](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6225)), "Gaussian" in resolveFamily ([src/R_interface_bartcore.cpp:1582](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L1582)).
7. C2 `stale` - the 40-line state-format comment ([src/R_interface_bartcore.cpp:6272-6311](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L6272-L6311)), ~28 lines
   of pre-release history opening "the shipped format (version 2)" 40 lines above
   `stateFormatVersion = 3`.
8. U2/U3 `sediment` - four `default:` arms suppressing -Wswitch on `ResponseFamily`
   ([src/bartcore/chain.hpp:766](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L766), [src/bartcore/chain.hpp:5033](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/bartcore/chain.hpp#L5033), [src/R_interface_bartcore.cpp:2298](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2298), [src/R_interface_bartcore.cpp:2842](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L2842)). Delete each and read the
   compiler; cheapest decision here.
9. V3 `dead-defense` - a refusal whose own comment says it is unreachable
   ([src/R_interface_bartcore.cpp:628-629](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L628-L629), annotated at [src/R_interface_bartcore.cpp:609-610](https://github.com/vdorie/dbarts/blob/b102e17cc7226338e7c7b26f73eedae617c8026f/src/R_interface_bartcore.cpp#L609-L610)).
10. D7 - `linearAlgebra_sse2.c` and `_avx.c` carry zero intrinsics: five routines written four
    times while the sibling partition family shares one body; D1's divergence rides on it.

Deliberately NOT here: SamplerBase, ResidT, the leaf concepts, the variance forest's parallel
lifecycle, the flat C API's 34 stan4bart-unused symbols - each checked, each with a second user or
a named plan (sec 1, sec 2's REFUTED paragraph).
