# Slice spec: out-of-sample per-forest predict replay
# (priced 2026-08-20; AMENDED per blind critique, verdict EXECUTE-WITH-AMENDMENTS)

Door: docs/plans/archive/bcf-bartcause-relocation.md:1271-1276. Verdict MEDIUM-PLUS:
expect ~700-800 raw insertions, >= 13 files. Comparables: 88ffe12
(multinomial C2 predict, 347/62/9) and 77013d72 (one new read method, 656
insertions / 15 files - the closer calibration for the Rd/bridge surface).

## What exists / what is missing (anchors critique-verified at 0f5a2285)

- Chain::predictFromSavedSample (src/bartcore/chain.hpp:2759-2777) reads
  forests_[0] only. Live twin predictFromCurrentTrees (:2784-2803).
- Chain::predictFromSavedSampleMulti (:2823-2850) is the multinomial
  precedent (loops K forests, softmaxes; live twin :2856-2884). It reports
  softmax probability channels keyed off numReportedLocations, not
  per-forest raw totals.
- Sampler::predictColumns (src/bartcore/sampler.hpp:540-567) dispatches on
  numReportedLocations() (:1196); AmplitudeForestCombiner does not override
  it (base default 1, combiner.hpp:664), so amplitude fits take the
  forest-0 path. That path is unreachable from R today: K=1 amplitude
  construction is refused ("a multi-forest model needs at least two
  forests") and refuseUndefinedTestFits guards K>=2
  (R_interface_bartcore.cpp:2849-2857, :5732).
- Storage is K-wide BY CODE: initializeSavedTrees loops every forests_
  entry (chain.hpp:2547-2580); savedTreeCapacity() (:2586) reads forest
  0's, set equal for all. (Do not cite the getTrees probe - $getTrees
  hardcodes forest 0L, R/dbarts.R:1860-1861, so it cannot establish this.)
- Serialization is ALREADY per-forest: ForestStateData
  (combiner.hpp:35-46) carries per-forest savedTrees/Params/Masks;
  getState sizes state.forests by forests_.size() (chain.hpp:2892-2896).
  stateFormatVersion (= 2, R_interface_bartcore.cpp:6171) needs NO bump.
  Restore rebuilds totalFits by summation, documented "equivalently, not
  bitwise" (combiner.hpp:60-66).
- In-sample already ships: forestFits/glue/bases (R/bart.R:216-241,
  368-381), extract(type="forest") (R/generics.R:389-435, extractForest
  :395, resolveForestSelection :367), identity documented man/bart.Rd:296,
  pinned inst/tinytest/test-argument-surface.R:743-800 (tolerance bar at
  :742-745).
- Reusable gates: shape.forestReportingIsDefined (facade.hpp:116-120;
  :114 is testFitsAreDefined - do not confuse).

## Binding conventions (BLOCKER fixes - implement exactly this)

- SCALE: the engine entry emits RAW totalFits per forest - NO fitScale,
  NO fitShift (response.scale IS fitScale, chain.hpp:1174; the run's
  forestFits channel is the same raw memcpy, chain.hpp:5209-5219, and
  R/bart.R:220-222 already multiplies by response.scale). The C++ test
  compares replay-at-training-rows against forestTotalFits (internal
  scale). R multiplies by response.scale exactly as in-sample.
- TOLERANCE, NOT IDENTITY: both new tests use tolerance 1e-12 (the
  test-argument-surface.R:742-745 bar). Bitwise is defeated twice over:
  the amplitude ridge rescales totalFits as c*sum but saved leaves as
  sum(c*) (combiner.hpp:1367 vs :1378-1383), and a measured single-forest
  keepTrees gap of 1.44e-15 exists already. Never assert identical().

## Slice sketch

engine: chain.hpp predictPerForestFromSavedSample + live twin (fork of
  :2823-2884 minus softmax, RAW - no scale) ~55; sampler.hpp
  predictPerForest + predictPerForestColumns (dense/CSC pair, out =
  nTest x K x capacity x nChains) ~45; facade.hpp one new virtual +
  forwarder ~10 (--preclean MANDATORY; delete benchmarks/kernels
  binaries); combiner.hpp 0 (gate exists).
  While editing chain.hpp: sweep its stale "BCF is a constant-leaf
  model" static_assert string (chain.hpp:745) to amplitude vocabulary;
  combiner.hpp's twin (:746) is NOT touched by this slice and stays.
bridge: R_interface_bartcore.cpp new bartcore_predictPerForest ~120-125:
  reuses parseTestSource / validateTestContainerAgainstStore /
  refuseSparseLeafCovariate, refuses on !forestReportingIsDefined, but
  CANNOT reuse predictFromSource (:5605-5711) - the K margin changes the
  dim vector and a per-forest raw fit takes NO offset (the offset belongs
  to the combination, as the shift does) - so it forks the allocation/
  wrapper logic. Plus declaration in src/R_interface_bartcore.hpp
  (pattern: :69) and R_interface.cpp DEF_FUNC.
R: R/bartcore.R wrapper ~15; R/dbarts.R R5 $predictForests ~30;
  R/generics.R predict.bart: ADD a `forest` formal (extract.bart's
  :254-259 is the pattern), new type value through
  validateType(type, eval(formals(predict.bart)$type)) (:208), a
  type="forest" arm reusing resolveForestSelection, BYPASSING the
  ci.level / probabilityFromLatents / s-attribute paths (:228-250) ~85.
Rd (~100 total): man/bart.Rd - new type value, new formal in \usage,
  Value shape, \seealso tie to extract's type="forest";
  man/dbartsSampler-class.Rd - rewrite the "Five methods report a fitted
  quantity" sentence (:321) + table row (:323-330), \alias, \S4method
  usage entry, \arguments, Value paragraph. Wording: say
  "amplitude-coupled", never "multi-forest" (a multinomial fit's raw f_k
  are defined but stay refused by the forestReportingIsDefined gate -
  deliberate scope; do not imply otherwise). The getTrees doc must not
  imply per-forest tree access exists (that gap stays, unpriced).
docs: feature-matrix.md:960-962 prose ("remains a door - S11 did not
  open it") is falsified by this slice - update that line (and its row).
messages: refuseUndefinedTestFits's text stays SURFACE-NEUTRAL (five
  call sites incl. flat C_interface.cpp:728/:775 - C consumers have no
  predict(type="forest")); extract's sample="test" refusal
  (R/generics.R:402-407) gets a message-only reword so it no longer
  claims the capability is absent (the stored channel still is).
tests: tests/cpp/test_sampler.cpp - per-forest replay at training rows
  vs forestTotalFits, tolerance-based, twin of :4729-4800 ~70;
  inst/tinytest - out-of-sample twin of the test-argument-surface.R:743
  identity at 1e-12, PLUS a storeState/loadState round trip through the
  new read (precedent test-bcf-r5-surface.R:214) ~110.
inst/NEWS.Rd ~8.
NOT touching inst/include/dbarts/dbarts.h: no ABI/hash move,
test-capi.R:56-60 inert, stan4bart unaffected. No new Rd topic, so no
_pkgdown.yml entry. predict.bart is the only R target (dbarts ships no
bcf()/bartBCF class; the door text's consumer list is stale).

## Draws

Predict-only, EXPECTED no draw movement (replay paths take no ext_rng,
mutate no chain state; forest-reporting channel rng-free,
chain.hpp:5203-5219). But this is an expectation, not a derivation:
combiner.hpp:921-930 records ulp movement from expression arrangement
alone. RUN the trio gates and treat any non-bitwise result as a defect
to diagnose, not to re-record. Expect equivalence 42 / bcf 12 /
multinomial 11 bitwise. Vtable hazard: new SamplerBase virtual
bus-errors against stale objects (--preclean; delete benchmarks/kernels
binaries).

## Scope boundary (binding)

Blocked on the test-basis modelling decision
(bart2-argument-consolidation.md:883-901; door
bcf-bartcause-relocation.md:1274-1276) - DO NOT TOUCH: resident test
basis, setTestPredictors/setTestOffset on an amplitude sampler,
run()$yhat.test for such fits, extract(type="forest", sample="test")
BEHAVIOR (message reword only, above), bcf() missing-response handling,
testFitsAreDefined itself.
NOT blocked: per-forest raw f_k at new rows. Full recombination
(shift + sum_k (B_k %*% g_k) * (scale * f_k)) is pure R arithmetic over
shipped $glue once per-forest fits exist; document as a 3-line idiom in
man/bart.Rd. A bases= argument on predict is EXCLUDED (reopens "what
predict means on a BCF fit"; belongs with the modelling decision).
No basis re-expansion is possible (data@bases is a bare matrix list,
R/A_class.R:514-524; no terms/xlev/contrasts retained) - any
caller-supplied basis arrives pre-expanded; retaining basis terms is a
separate formula-machinery item, not in this slice.

## Mutation proofs (prove the gates discriminate; record which move)

M1: make the engine loop read forest 0 into every channel -> the
per-forest C++ test and the R out-of-sample twin must move; in-sample
extract cells stay green (they never touch the new entry).
M2: apply fitScale inside the engine entry -> both new tests move
(squared-scale on the R side, scale mismatch on the C++ side).
M3: add the offset in the new bridge path -> the R twin with a non-null
offset moves.
M4: drop the forestReportingIsDefined refusal -> the multinomial
refusal cell moves (in-bounds silent wrong value otherwise).

## Risks

- Per-forest numTrees may differ; loop must use each forest.numTrees
  (the Multi path already does, chain.hpp:2835).
- rchk: new bridge entry with a fresh allocation shape.
- refuseUndefinedTestFits message: surface-neutral update only (above).
