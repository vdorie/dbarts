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
  (facade.hpp:139-401); the one source-level impl `SamplerFacade<L, ResidT>` (facade.hpp:405) is a
  TEMPLATE instantiated five ways: ConstantGaussianLeaf/double, ConstantGaussianLeaf/float (:725),
  MonotoneConstantGaussianLeaf (:793), GPGaussianLeaf (:813, :846), LinearGaussianLeaf (:817,
  :850). The vtable is what lets the bridge hold `unique_ptr<SamplerBase>` without knowing L.
  `extension-point`, high.
- R2. `ResidT` has TWO instantiations; `float` is reachable from R (`control@storage == "single"`
  -> R_interface_bartcore.cpp:398-404 -> options.fp32Residual -> facade.hpp:725).
  chain.hpp:4463-4480 documents why both spellings exist - codegen, not taste.
  `compile-time-only`, high.
- R3. All nine leaf concepts are live. VectorLeafModel (model.hpp:60) and FunctionLeafModel (:85)
  head no template but are disjuncts of IntegrableLeafModel (:98), which heads
  chain/sampler/facade/combiner. ScaleLeafModel (:111) is the second disjunct of
  MoveScorableLeafModel (:126) and IS exercised: chain.hpp:4222 instantiates moves.hpp with
  ConstantVarianceLeaf. `compile-time-only`, high.
- R4. The flat C API is fully exercised: of 62 symbols in inst/include/dbarts/dbarts.h every
  function is called by stan4bart (bartcore branch, 28) or by inst/tinytest/capi/consumer.c, and
  the 10 with neither are typedef/enum names and the `dbarts_stub_*` macro machinery. No dead
  exports.
- R5. Per-observation loops are monomorphic. The sweep's virtual boundaries are per-sweep
  (chain.hpp:1511-1521, :1358) and per-sweep-per-forest (:1407-1408); nothing virtual sits inside
  the per-tree loop (:1428-1503) or below, exactly what core-generalization.md:70-76's dispatch
  table promises. The `misc_*` function-pointer kernels are per-node-op, sanctioned by
  kernel-vocabulary.md:12-16.
- R6. Comment rot is NOT in the docs/ citations: all 29 cited design docs exist and 38 sampled
  anchors resolve to what the citing comment claims. Only blemish: R_interface_bartcore.cpp:6285
  cites `public-surface.md 2`, whose TITLE is about factor ingestion though its body carries the
  registry rule past line 155.

## 2. DUPLICATION
DIVERGE first - a divergence outranks a large exact clone.

- D1. `alpha == 0.0` short-circuit exists only in the scalar kernels: linearAlgebra.c:67 and :85
  guard `if (length == 0 || alpha == 0.0) return;`; the same two functions guard on length alone
  at linearAlgebra_sse2.c:42, :60, _avx.c:53, :75, _neon.c:150, :220. The C fallback writes
  nothing where SIMD dispatch writes, in a library whose stated invariant is within-host bitwise
  identity across dispatch. Live callers (model.hpp:2833, :2912, :2968, :3033) pass `-min_`, i.e.
  `-0.0` when min_ is 0; `x += -0.0` is a round-to-nearest no-op for every finite x, so observable
  difference is confined to NaN payloads. `diverge`; high on the fact, low on impact.
- D2. Scale-leaf positivity applied on 3 of 4 state paths. chain.hpp:3301 (live variance trees)
  and :3334 (saved buffer) refuse a non-positive variance leaf, and sampler.hpp:919-923 applies it
  to a warm start's SLOT-sourced trees with a comment naming the reason ("the buffer is
  hand-buildable and a rebuild scatters the leaf straight into a divisor"). The LIVE-sourced arm,
  sampler.hpp:891 (`dst.varianceTrees = src.varianceTrees;`), applies nothing and neither donor
  parser checks (R_interface_bartcore.cpp:7233-7248); `rebuildVarianceForest`
  (chain.hpp:4356-4368) then scatters it into a divisor. Reachable from a `.Call` with a
  hand-built state, not from R's own state objects. High. CORRECTNESS CANDIDATE.
- D3. `applyNewData` (chain.hpp:2453) and `recoverTreeParameters` (:2421) take `forests_[0]` only
  where siblings loop `forests_` (:2329, :2401, :2526, :2540). A whole-data replacement on a
  BCF/multinomial chain would leave forests 1..K-1 on the old grid; guarded two layers away by
  `refuseMultiForestMutation` (R_interface_bartcore.cpp:4632), with no engine-level assert. High.
  CORRECTNESS CANDIDATE.
- D4. One enumerator, two opposite meanings, same file. `augmentationFamily`
  (R_interface_bartcore.cpp:6151) rejects "gaussian" and maps "student" onto
  `ResponseFamily::gaussian`; `drawAugmentationLaws` reads `case RF::gaussian:` as "the Student-t
  scale mixer" (:6225), while `resolveFamily` (:1582) maps "gaussian" to the same enumerator
  meaning the Gaussian law and never accepts "student". `diverge`/`sediment`, high; a dedicated
  token costs nothing.
- D5. The shipped-BCF glue draw and the general amplitude draw disagree on the scale-mixture
  refresh: combiner.hpp:1171 refreshes only `prior[0].variance`, :1213 refreshes every forest with
  `halfCauchyScale > 0`. The selector at :986 routes on `shippedShape()` (:1474), which tests
  basis widths and canonicality only, NOT the half-Cauchy flag - so a two-forest spec with
  canonical bases and a per-forest prior scale (admitted at R_interface_bartcore.cpp:2178-2180)
  gets a fixed-variance prior where the general path would sample it. The duplication is
  deliberate and well documented (:940-984, with a stated deletion trigger); the PREDICATE's
  coverage is the finding. Medium-high. CORRECTNESS CANDIDATE.

AGREE-but-restated (reading cost and drift risk, no defect today):

- D6. Six augmentation laws restated in the bridge: R_interface_bartcore.cpp:6145-6270 (126 lines)
  re-implements Probit/Ordinal/AFT/Logistic/NB/T response draws (model.hpp:3085, :3220, :3548 and
  the NB/T drawers). :6146-6149 declares it deliberate (a different generator, citing
  r-c-division.md); nothing enforces lockstep - no shared kernel, no cross-check test named at the
  site. High.
- D7. `linearAlgebra_sse2.c` and `_avx.c` contain ZERO intrinsics (`grep -c '_mm'` = 0 for both;
  `_neon` has 89) - they are the `_c` bodies with a different unroll factor (4 vs 8), relying on
  auto-vectorization. With `_neon` (303 lines, real intrinsics) and linearAlgebra.c (286) that is
  five routines written four times, while the sibling partition family shares one 350-line
  `partition_body.c` behind 28-line wrappers, an idiom used nowhere else in the directory.
  `moments.c:251-344` vs `:354-447` repeats it: four fp32 suffstat kernels are exact clones of the
  fp64 four, 94 lines kept in lockstep by hand (:349-353 makes the summation order a correctness
  contract). High.
- D8. Five "build a scratch tree from flat, check containment" walks (chain.hpp:3199-3230,
  :3303-3324, :3365-3382, :3409-3428, :3441-3450): same eight lines, but the first two `return
  false` on a failed build and the last three `continue`, and the mask selection is spelled two
  ways (:3206-3210 vs :3411-3415). Each comment names the others as its mirror; no doc covers it.
  High.
- D9. Eight near-parallel saved-vs-live replay functions (chain.hpp:2799/:2825, :2865/:2898,
  :2938/:2962, :904/:932). `predictVariance` (:934) guards `if (!varianceForest_) return;`; its
  saved twin (:906) dereferences immediately, unreachable today only because the bridge gates on
  `shape.hasVarianceForest` (R_interface_bartcore.cpp:5732). High.
- D10. Leaf-shape flatten switch three times: chain.hpp:2683-2716, :2652-2677, :3004-3040 - the
  first two a documented cache-vs-recompute pair, getState's a different format. High on
  duplication, medium on intent.

REFUTED - the memo's "the variance forest re-implements the mean forest's lifecycle as 12 parallel
members though ConstantVarianceLeaf already satisfies ScaleLeafModel". The parallelism is real and
larger (15 shadowed members, chain.hpp:407-484 vs combiner.hpp:138-231, plus ~15 paired Chain
methods) but the premise does not carry: model.hpp:322-324 static_asserts
`!IntegrableLeafModel<ConstantVarianceLeaf>`, and Forest (combiner.hpp:137) and Chain
(chain.hpp:532) are both constrained on IntegrableLeafModel. What CAN be shared is -
MoveScorableLeafModel admits it, `logLikelihoodForBranch` scores it unchanged. The blocker is the
combination law (Forest's treeFits/totalFits/treeY are additive; the variance forest is
multiplicative, chain.hpp:400-407), logged as open debt in forest-combiner.md:207-211, :252 and
heteroscedastic.md:59-62, :586-590. `extension-point`, high.

## 3. DEAD OR UNREACHABLE
- U1. `Tree::rightChildOf` (tree.hpp:366) and `Sampler::setCurrentSampleNum` (sampler.hpp:485):
  zero references in src/, tests/cpp, benchmarks/, R/, inst/. `sediment`, high.
- U2. chain.hpp:766 `default:` on a 3-arm family switch. `createAmplitudeSampler`
  (facade.hpp:869-873) refuses every family but gaussian/probit/logistic first and the comment at
  :750-755 says so. Unreachable today; the cost is that `case gaussian: default:` suppresses
  -Wswitch, so a 7th enumerator would silently fit Gaussian here. `sediment` - defense that
  disables the compiler's own check.
- U3. Three further -Wswitch suppressions on `ResponseFamily`: chain.hpp:5033 (folding
  gaussian/aft/ordinal/nbinom), R_interface_bartcore.cpp:2298, :2842. Of 7 family ladders only 3
  are exhaustive. Deleting each `default:` and reading the compiler is the cheapest check on this
  list.
- U4. `refusedAmplitudeFamilyReason` (R_interface_bartcore.cpp:2268-2284): exhaustive switch,
  every arm returns, so `return "this response family";` (:2283) is unreachable - but required to
  silence -Wreturn-type. `compile-time-only`; touch only alongside U3.
- U5. `misc_simd_getMaxSIMDInstructionSet` (src/misc/simd.c:142-183) detects SSE, SSE3, SSSE3,
  SSE4_2, AVX512F, AVX512VL, AVX512BW; the dispatcher (:284-340) tests only AVX2, SSE4_1, SSE2,
  AVX, NEON. Seven levels select nothing, and the AVX512 block (:170-183, with its own second
  `__cpuidex`) computes state no branch reads. One-time cost; inherited from pre-bartcore misc
  (main carries simd.c at 9b0ae65b) but 125 of its lines changed here. `extension-point` if an
  AVX512 kernel is planned, else `sediment`.
- U6. src/misc/partition_body.c:106, :123, :149 - a commented-out PARTITION_RANGE if/else/endif
  triple wrapping two abandoned NEON load strategies (:107-123, inside a block comment): ~18 dead
  lines, the same shape as `#if 0`, which the memo's sec 1.0 reported as zero occurrences.
- U7. `XINT_TYPE` width generality has no user and a wrong-answer failure mode.
  kernel-vocabulary.md:24-27 says the code type is "configure-selected via `--with-xint-size`";
  configure.ac:21 hard-wires `uint16_t` with no such option and no width-suffixed kernel exists.
  The `#ifndef XINT_TYPE` guard (src/include/misc/types.h.in:6) still lets a CPPFLAGS define
  override it, and every SIMD partition kernel hard-codes `epi16`/`u16` (partition_body.c:11, :73,
  :151) with no `static_assert(sizeof(misc_xint_t) == 2)` anywhere. `sediment`, high.

## 4. DEFENSIVE CODE WITH NO REACHABLE TRIGGER (bridge)
- V1. The state-format compatibility window is empty by construction: `stateFormatVersion` and
  `minReadableStateFormatVersion` are both 3 (R_interface_bartcore.cpp:6312, :6321), and the
  comment at :6310-6311 concedes "no release ever shipped format 3". The two floor checks (:6658,
  :7282) can only fire on a state with no version attribute (reads 0). Forward-looking; the
  registry rule at :6285-6298 is the durable part. `extension-point`, high.
- V2. A check byte-identical to its sole caller's: `createMultinomialCountsHolder` re-runs `if
  (!Rf_isInteger(countsExpr) || Rf_xlength(dimsExpr) != 2)` (:3379-3380) which
  `createMultinomialDataHolder` already ran at :3457-3458, with only `parseMultinomialData` (which
  never touches countsExpr) between. The helper has C++ linkage only, so it is not a flat-C
  backstop. `dead-defense`, high. Its surrounding PROTECT is separately annotated (:3374-3377) as
  deliberate analyzer bait - decide the two separately.
- V3. Author-acknowledged unreachable: `resolveCscCategoricalReferences`'s `if (source >=
  numSparseColumns) Rf_error("%s", boundMessage)` (:628-629). Its own doc comment at :609-610 says
  so, and `mapColumnSources` (:570) bounds every CSC entry at all three call sites.
  `dead-defense`, high; the odd part is that the comment names the redundancy and the check stays.
- V4. PROTECT convention applied inconsistently inside one function. Four pairs guard an attribute
  of the already-rooted `stateExpr`/`donorStateExpr` argument with nothing allocating in the
  window: :6649-6653 (formatVersion in setState), :7277-7281 (its copy in installForests),
  :6731-6740 (sampleNum), :6742-6753 (recordedDraws). In the SAME functions, sibling reads of the
  identical pattern skip PROTECT (:6699-6700, :6711, :7079, :7094-7095). None is among the four
  sites docs/plans/release-candidate-review.md's rchk note (:1138-1173, :1300-1330) records, and
  all are distinct from the two `setState` rchk BAILOUTS the note calls already-balanced. The
  finding is the inconsistency, not the PROTECTs. Medium.
- V5. `rbart_getFitted`'s two PROTECTs (R_interface_rbart.cpp:16-17): both dims are reduced to raw
  `int*` at :19-20 and last read at :40-42, and the only allocation is `rc_newReal(n)` at :44.
  Untouched by the rchk commit. Medium - defensible only if `rc_getDims` can itself allocate.
- V6. NOT dead, keep. `setState`'s three "already-non-null" guards (:6799, :6813, :6830) make the
  null branch in `readFunctionTreeParams` (:5490), `readTreeParams` (:5469) and `readTreeMasks`
  (:5438) unreachable on setState's path, but the same helpers are called unguarded from
  `readWarmStartState` (:7144-7164).
- V7. NOT dead, keep - this closes the "R already refused it" class. No C-side check on the .Call
  path is provably unreachable: every `dbarts_bartcore_*` symbol is an ordinary DL_FUNC in
  `R_callMethods[]` (R_interface.cpp:180-259), reachable by any R code holding the internal `C_`
  name; two sites say so themselves (:3987-3993, :4878-4882). The S4 route confirms it:
  `dbartsControl`'s `setValidity` (R/A_class.R:272-383) is not re-run on plain `@<-` mutation and
  `setControl` (R/dbarts.R:1155-1206) never calls `validObject()`. ~13 more sites carry the same
  explicit "backstop" annotation.
- V8. rc constraint API re-derived: 42 constrained calls (35 `rc_getInt`/`rc_getDouble`/
  `rc_getBool` carrying RC_LENGTH or RC_VALUE, 7 `rc_assert*Constraints`, 0 bare `rc_get*0`)
  against 111 unconstrained `rc_getListElement` fetches and 168 open-coded `Rf_is<Type>` checks in
  R_interface_bartcore.cpp (3 in C_interface.cpp, deliberately
  - :3-5). The memo's "7" counted only `rc_assert*Constraints`. Of 50 `dbarts_bartcore_*` .Call
    entries (59 rows total), 20 call one of the 14 shared `refuse*`/`validate*` helpers directly;
    the rest go through `createHolder`/`parseData`/`parseModel`. 304 `Rf_error` sites in the
    bridge, 42 in the C API.

## 5. COMMENTS
- C1. THREE `stale` "Not yet exposed through the R surface" claims - all three features ARE
  reachable from R today: chain.hpp:88 (monotone) vs R/model.R:99, :526-531; chain.hpp:141
  (interaction constraints) vs R/spec.R:415-416 and R_interface_bartcore.cpp:1254; chain.hpp:166
  (variance forest) vs R/spec.R:532, R/dbarts.R:816, R_interface_bartcore.cpp:2017. A reader who
  trusts these will not look for the R-side refusals that guard them. High. Best comment finding
  here.
- C2. R_interface_bartcore.cpp:6272-6311, a 40-line block comment of which ~28 lines narrate three
  pre-release format iterations that never shipped. It opens "The shipped format (version 2)"
  while `stateFormatVersion = 3` sits 40 lines below - a contradiction inside one comment. The
  load-bearing part is the registry rule (:6285-6298, ~13 lines). `stale` + narration, high. Worst
  site by volume.
- C3. The other worst narration sites (criterion: the payload is a comparison to code that no
  longer exists and cannot be seen), most-useless first - chain.hpp:1041 ("the family conjunct
  that used to stand beside it was there because setResponse handed forest 0's bare totals as
  though combined"); R/generics.R:1154 and :1287 ("the same keepTrees gate a deleted `$bc` field
  used to"); R/utility.R:119 ("gone now that the rename has landed"); R/utility.R:49 ("silent
  before bart2 could forward a resid.prior object at all"); R_interface_bartcore.cpp:2766 (quotes
  an old, replaced error message); C_interface.cpp:624 ("which this entry used to drop on the
  floor"); R/data.R:453 (three appeals to an invisible prior implementation in one comment);
  R/data.R:342 ("this used to be a function evaluated in the caller's frame"); R/spec.R:490,
  R/model.R:1759 and R/spec.R:652 (removed flat formals and a vanished literal); R/bart.R:769 and
  R/dbarts.R:627 (fixed-bug narration). 20 sites total, 14 in R/ and 6 in C++; src/misc,
  src/external, src/rc and src/include are narration-free.
- C4. Forward-looking PLAN narration (a plan is not a constraint): model.hpp:2589-2594 ("v1 ships
  the exact integer envelope... A later real-shape mode routes a fractional b"), model.hpp:2774,
  facade.hpp:784 ("v1 keeps the mean leaf constant"), R_interface_bartcore.cpp:3008, :3014, :3030.
  Six sites; each also states a live refusal, so only the forward half is narration. Medium.
- C5. 262 docs/ citations no installed user can follow (`.Rbuildignore` has `^docs$`). The
  decision is not "are they rotten" - they are not - but whether a shipped comment should cite a
  stripped path.

## 6. SEAMS THAT COST AT RUNTIME
- S1. PER-OBSERVATION VIRTUAL DISPATCH, contradicting the architecture's own rule.
  `PredictorUpdateSession` (sampler.hpp:89-100) declares `observationWouldRemainValid(i)` and
  `commitObservation(i)`; `updatePredictorPerObservationJointly` (facade.hpp:694-703) calls both
  inside `for (i = 0; i < numObservations; ++i)`, once per sampler per observation.
  core-generalization.md:69-76 states "nothing dispatches per observation" and "Per obs | none:
  monomorphic loops/kernels". The erasure is NOT removable by templating - the joint sweep takes
  `SamplerBase* const*` over samplers of possibly different L, its whole purpose
  (R/updatePredictorPerObservationJointly.R; the bairrtt consumer) - though the single-sampler
  path (facade.hpp:217) pays it without needing it. Frequency PER-OBSERVATION, 2 virtual calls x
  numSamplers x n; the work inside is a per-tree descent so the ratio is likely fine, but either
  the code or the doc's absolute rule should move. `extension-point`, high; not benchmarked.
- S2. `std::function` on the run path: pollInterrupt/shouldCancel (facade.hpp:147, sampler.hpp:274,
  :362-370, chain.hpp:1345) and SweepCallback (chain.hpp:396), both called once per sweep
  (chain.hpp:1358, :1362). PER-SWEEP, negligible; listed to close the class.
- S3/S4. `combiner_->` virtuals: 44 sites in chain.hpp, the two inside the sweep's forest loop
  being drawForestGlue/formForestResponse (:1407-1408, prior-sampling copy at :1973-1974) -
  PER-SWEEP-PER-FOREST, two impls plus a null default. `response_->` virtuals in `run`
  (chain.hpp:1511-1521) - PER-SWEEP, nine impls.

## 7. LOWER-VALUE, FOR COMPLETENESS
- L1. Fourteen `*ForTesting` accessors on the shipped engine that no production code calls (all
  reached from tests/cpp, three by exactly one assertion): they widen the engine's member surface
  for testability only. `extension-point` (test seam).
- L2. 113 named helpers appear exactly twice in the 30488-line shipped C++ corpus (definition plus
  one call site); 84 are never referenced from tests/cpp. Re-derived; the memo says 144/103 over
  the engine alone. Decomposition, not sediment - the 84 with no test reach belongs to lens F2.
- L3. core-generalization.md, the arbiter document, names two extension points no shipped code
  carries: `SplitSelector` (:129-133) and `MoveStrategy` (:118-126, template rule at :86), zero
  occurrences of either in src/. It cannot justify a seam by naming it alone; check the code.

## 8. THE TEN TO DECIDE FIRST
Ranked by maintainer time saved per decision, not by size.

1. C1 `stale` - three "Not yet exposed through the R surface" comments (chain.hpp:88, :141, :166)
   on features R reaches today. One-line fix, misleading now.
2. D2 `diverge` - variance-leaf positivity checked on 3 of 4 state paths (sampler.hpp:891 vs
   :919-923). Correctness candidate.
3. D3 `diverge` - `forests_[0]` hardcoded in applyNewData (chain.hpp:2453) and
   recoverTreeParameters (:2421) where siblings loop. Correctness candidate.
4. D5 `diverge` - `shippedShape()` (combiner.hpp:1474) routes on basis shape but not the
   half-Cauchy flag, so two amplitude specs get two different models.
5. S1 `extension-point` - per-observation virtual dispatch (facade.hpp:694-703) against
   core-generalization.md:69-76's "nothing dispatches per observation"; the code or the doc's
   absolute rule moves.
6. D4 `diverge` - `ResponseFamily::gaussian` means "Student-t" in the augmentation surface
   (R_interface_bartcore.cpp:6225), "Gaussian" in resolveFamily (:1582).
7. C2 `stale` - the 40-line state-format comment (R_interface_bartcore.cpp:6272-6311), ~28 lines
   of pre-release history opening "the shipped format (version 2)" 40 lines above
   `stateFormatVersion = 3`.
8. U2/U3 `sediment` - four `default:` arms suppressing -Wswitch on `ResponseFamily`
   (chain.hpp:766, :5033, R_interface_bartcore.cpp:2298, :2842). Delete each and read the
   compiler; cheapest decision here.
9. V3 `dead-defense` - a refusal whose own comment says it is unreachable
   (R_interface_bartcore.cpp:628-629, annotated at :609-610).
10. D7 - `linearAlgebra_sse2.c` and `_avx.c` carry zero intrinsics: five routines written four
    times while the sibling partition family shares one body; D1's divergence rides on it.

Deliberately NOT here: SamplerBase, ResidT, the leaf concepts, the variance forest's parallel
lifecycle, the flat C API's 34 stan4bart-unused symbols - each checked, each with a second user or
a named plan (sec 1, sec 2's REFUTED paragraph).
