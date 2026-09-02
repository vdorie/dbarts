# multinomial-test-surface

agent: opus (a REPORTING/REPLAY arc, not a posterior arc - the softmax test
  blend and the out-of-sample prediction are both deterministic functions of
  draws the shipped multinomial sampler already produces, so no sampled
  quantity is added. The care is entirely in keeping the seam byte-identical on
  the single-forest and BCF paths and getting the K-forest test/predict shape,
  the level-shift invariance, and the surface threading right).
rng: ONE CLASS, throughout: RNG-NEUTRAL. Test fits are formed from each forest's
  totalTestFits (populated in the sweep from tree fits, no draw); predict
  replays saved trees (no draw); combinedTestFits and the softmax consume no
  rng; flipping testFitsAreDefined() and un-refusing keepTrees move no draw
  (storeSavedTreeRecord is pure serialization). Every commit stays bitwise on
  ALL THREE standing anchors: single-forest equivalence 22/22 IDENTICAL vs
  equivalence-ac6ec2c.rds, the BCF fixture 5 scenarios x 6 channels IDENTICAL vs
  bcf-equivalence-99205ee.rds, and the multinomial fixture's EXISTING train
  channels IDENTICAL vs multinomial-equivalence-bb8855e.rds (the new test
  channel is additive). If any anchor moves, the seam touched a draw path it
  must not.
window: NONE. The public surface already ships (bart2(family = "multinomial"),
  bb29d00); this arc removes by-name refusals on it (test=, predict, keepTrees)
  and adds the outputs behind them. No new family slot, no dbarts.h change, no
  pre-release window of its own.
budget: ~500-600 lines. C1 test-at-creation ~290 (combinedTestFits virtual +
  storeSample test seam + testFitsAreDefined accessor + bridge test-wiring +
  refusal relax + the internal and bart2 R surface + one exact-gate arm + the
  fixture test channel + tinytest). C2 predict-on-newdata ~220 (K-forest
  saved-tree replay + the predict bridge shape + un-refuse keepTrees +
  predict.bartMultinomial + gate/cpp). C3 docs ~50. Chiefly
  src/bartcore/combiner.hpp, chain.hpp, facade.hpp, sampler.hpp,
  src/R_interface_bartcore.cpp, R/bartcore.R, R/bart.R, R/generics.R,
  benchmarks/R/multinomial-exact.R + multinomial-equivalence.R + a re-recorded
  baseline, tests/cpp/test_sampler.cpp. Header edits -> --preclean; delete stale
  tests/cpp binaries.

## Goal

Give the shipped multinomial softmax model the out-of-sample capability it
deferred at bb29d00's close (docs/plans/archive/multinomial.md C7 landing note, the
scope-narrowings paragraph): (1) TEST-AT-CREATION - x.test supplied at fit
time reports per-draw test-set softmax probabilities in the run's test channel
(the mbart2 prob.test pattern), and (2) PREDICT-ON-NEWDATA - predict(fit,
newdata) replays the K saved forests to K-column probabilities. Both are
identified deliverables (probabilities, not the non-identified raw f_ik). The
whole arc is RNG-neutral: it reports and replays existing draws, so every
posterior path - single-forest, BCF, and the multinomial train channel - stays
bitwise. dbarts.h stays frozen (all internal + R surface).

## Why the blend is well-defined here (unlike BCF)

BCF leaves testFitsAreDefined() false for a real reason: a mu + b_z tau
off-sample needs a test treatment vector the API never carries, so only the
bare prognostic forest could be recorded and it would misreport
(refuseBCFTestSurface, [[R_interface_bartcore.cpp:1403-1409@6a48351b]];
BCFForestCombiner::testFitsAreDefined, [[combiner.hpp:490@6a48351b]]). The multinomial test
blend has no such gap: softmax over the K forests' totalTestFits is the SAME
map combinedFits already applies to the K forests' totalFits ([[combiner.hpp:597@6a48351b]]-
611), and every one of the K forests already accumulates its own totalTestFits
in the sweep (resizeTestStorage loops all forests, [[chain.hpp:491-498@6a48351b]]; the
per-forest test accumulation runs inside the `for f` sweep loop,
[[chain.hpp:585-641@6a48351b]]). What is missing is a combined-TEST-fits seam: storeSample's
test loop reads the reported forest's totalTestFits directly (single-forest
shaped, [[chain.hpp:2191-2194@6a48351b]]) instead of asking the combiner for a K-location
test blend.

### The level-shift invariance (why the test channel is correct without touching it)

MultinomialForestCombiner::afterCombine's grand-shift c is added to every f_ik
uniformly across all K forests ([[combiner.hpp:655-663@6a48351b]]) but is NOT applied to
totalTestFits (nor need it be: keepTrees/test were out of scope at C4). This is
fine and must be stated as a correctness fact, not a bug: the softmax is
invariant to a common shift of all K categories, and c is common, so
softmax({totalTestFits_k}) == softmax({totalTestFits_k + c}). At storeSample the
train channel reads totalFits WITH this sweep's c (afterCombine at
[[chain.hpp:673@6a48351b]] precedes storeSample at 692) and the test channel reads
totalTestFits WITHOUT it; both recover the softmax of the identified log-odds
and agree. Same argument makes the predict replay correct although saved
(flattened) tree leaves carry no c.

## Context (seams, read in code)

- The combiner base virtuals ([[combiner.hpp:230-292@6a48351b]]): formForestResponse (236),
  combinedFits (241), the reporting map reportedForest/testFitsAreDefined/
  logLikelihoodIsDefined (276-278), numReportedLocations (286),
  drawForestGlue/afterCombine (257-269). combinedTestFits is the ONE new virtual
  this arc adds; it slots beside combinedFits, base inert, MultinomialForest-
  Combiner overriding it as combinedFits' totalTestFits analog. BCF leaves
  testFitsAreDefined false so it is never called on BCF.
- MultinomialForestCombiner ([[combiner.hpp:540-688@6a48351b]]): combinedFits softmaxes over
  totalFits log-sum-exp-safe (597-611); testFitsAreDefined() returns false
  (559) - this is the flip. numReportedLocations() = K (558). The test analog is
  the same loop over forests[k].totalTestFits into an n_test x K
  location-major scratch.
- storeSample's test loop ([[chain.hpp:2179-2198@6a48351b]]): if a combiner has
  testFitsAreDefined() false it NaN-fills nTest * numLocations (2182-2189, BCF
  and today's multinomial); else it loops `loc` over numLocations writing the
  REPORTED forest's totalTestFits into every channel (2191-2197) - correct at
  L = 1 (single-forest), WRONG for K locations (it would copy forest 0's test
  fits into all K channels). The seam: introduce a `const double* testSrc` that
  is forest.totalTestFits.data() for a single-location test-defined path and
  combiner_->combinedTestFits(forests_) (n_test x K location-major) for a
  multi-location one; the existing `scale * testSrc[loc*nTest+i] + shift` loop
  is then byte-identical at L = 1 (loc 0 reads totalTestFits[i]).
- The reporting accessor gap: chain has numReportedLocations() ([[chain.hpp:438@6a48351b]])
  surfaced on SamplerBase ([[facade.hpp:145@6a48351b]]/339) but NO testFitsAreDefined()
  accessor - it is combiner-internal, read only in storeSample. To let
  refuseBCFTestSurface allow multinomial while still refusing BCF, surface a
  chain testFitsAreDefined() (combiner_ ? combiner_->testFitsAreDefined() :
  true) and a SamplerBase override beside numReportedLocations.
- The run allocation ALREADY handles L > 1 test (the C2 location seam, commit
  9030d93): allocFitsArray sizes nTest x numLocations x numSamples (x numChains)
  and the test slot picks it when numLocations > 1 (R_interface_bartcore.cpp:
  2099-2116, 2147-2157); results.testFits is wired whenever numTestObservations
  > 0 (2203) and results.numReportedLocations is set (2212). NOTHING in the run
  allocation changes - it is already K-test-shaped and dormant only because
  multinomial creation wires no test data and NaN-flags it.
- Creation refuses test only by omission: createMultinomialHolder
  ([[R_interface_bartcore.cpp:1667-1729@6a48351b]]) parses the data (parseData reads x_test/
  numTestObservations) but never calls setTestPredictors - unlike the
  single-forest createHolder, which wires it at 1552-1567. The multinomial path
  carries no offset (design note "The surface"), so test offset stays refused.
  refuseBCFTestSurface guards setTestPredictor (2456), setTestOffset (2501),
  setTestPredictorAndOffset (2527), and predict (3269) on numForests >= 2 -
  currently catching multinomial; it must gate on !testFitsAreDefined() instead.
- Predict addresses forest 0 only: predictFromSavedSample ([[chain.hpp:1389-1408@6a48351b]],
  `const Forest<L>& forest = forests_[0]`) and predictFromCurrentTrees (1413-)
  sum a single forest; sampler.predict loops chains/slots into an nTest x
  capacity (x nChains) buffer ([[sampler.hpp:491-504@6a48351b]]); the predict bridge
  allocates that single-location shape ([[R_interface_bartcore.cpp:3283-3304@6a48351b]]). The
  saved trees themselves ARE per-forest already: storeSavedTreeRecord runs
  inside the `for f` sweep loop ([[chain.hpp:631@6a48351b]]), so every kept sample stores all
  K forests' trees (Forest<L>::savedTrees, [[combiner.hpp:126-132@6a48351b]]). K-forest
  replay + softmax is the only new logic - the harder half.
- R surface refusals to lift:
  - R/bartcore.R: bartcoreMultinomialSampler (556-580) passes no x.test to
    createMultinomial; bartcoreSetTestPredictor (672-675)/bartcorePredict
    (780-786) exist and route to the same .Call entries BCF uses. Comment at 555
    "Out-of-sample prediction and test fits are refused this arc" is the doc to
    correct.
  - R/bart.R: the multinomial branch refuses test= (429-434), keepTrees
    (452-457); bart2Multinomial (666-714) builds host+engine and runs;
    packageMultinomialResults (731-761) reshapes only samples$train, names
    levels on the trailing K margin - the test channel is the same reshape on
    samples$test.
  - R/generics.R: extract.bartMultinomial refuses sample="test" (407-412);
    predict.bartMultinomial refuses outright (456-461); fitted.bartMultinomial
    (439-454) is train-only.
- The exact gate (benchmarks/R/multinomial-exact.R): three train-only arms
  (intercept K=3, K=2==logistic, covariate K=3), each matching the posterior
  mean of identified train probabilities to quadrature (header 1-56). A test arm
  is cheap and strong: pass test rows duplicating an arm's cells and assert the
  test-channel probabilities equal the train/quadrature target to MC error.
- The bitwise fixture (benchmarks/R/multinomial-equivalence.R): records train
  (K softmax channels), per-category forestFits, per-category varcount, over a
  K=3 covariate and a K=2 scenario (recordChannels 63-75). It byte-guards the
  K-channel OUTPUT seams the anchors miss (its own header 3-14). A test channel
  is additive: record result$test and per-category test fits; the existing three
  channels stay bitwise (adding x.test moves no train draw).

## Constraints

- dbarts.h FROZEN. Multinomial is internal (bartcore .Call) + the bart2 R
  surface; no public C-API test-fits path exists for it, so no ABI concern.
- The single-forest AND BCF byte layouts are UNTOUCHED. combinedTestFits is a
  new inert base virtual (BCF never calls it, testFitsAreDefined stays false);
  storeSample's L = 1 test write stays the bare totalTestFits path; the run
  allocation is unchanged (already K-test-capable). Both anchors IDENTICAL and
  the multinomial train channels IDENTICAL are the per-commit proof.
- RNG-NEUTRAL every commit. No sampled quantity is added; test fits and predict
  are pure functions of existing draws. This is the whole neutrality claim.
- Levels/dimnames threading on the new K-shaped test/predict outputs follows the
  C7 conventions (packageMultinomialResults names levels on the trailing K
  margin; every K-shaped output threads them for the object's lifetime).
- Single-trial scope unchanged (labels 0..K-1, n_i = 1). No offset, no weights,
  no count matrix - those stay the recorded multi-forest-models follow-ups.
- OUT of scope (stay refused, by name, as C7 left them): the latent escape hatch
  (type="bart"; raw f_ik non-identified and unrecorded), per-sample per-category
  varcount, the count-matrix likelihood, the formula interface.

## C1 - Test-at-creation (the mbart2 prob.test pattern), one gated commit

The engine seam is unreachable from R without the bridge and surface wiring, so
this lands end to end as one gated commit (the multinomial arc's C4 discipline:
nothing lands before the gate that exercises it). RNG-NEUTRAL.

Sub-parts:

(a) ENGINE. Add combinedTestFits(forests) to ForestCombiner (combiner.hpp,
   beside combinedFits): base inert (returns nullptr; never called because base
   testFitsAreDefined() is true only for the single-forest null-combiner path,
   which storeSample handles combiner-free). MultinomialForestCombiner overrides
   it - the combinedFits softmax over forests[k].totalTestFits into an
   n_test x K location-major combinedTest_ scratch (reads data_.numTest-
   Observations). Flip MultinomialForestCombiner::testFitsAreDefined() to true
   ([[combiner.hpp:559@6a48351b]]). storeSample ([[chain.hpp:2179-2198@6a48351b]]): in the else (test
   defined) branch, set `const double* testSrc` = (combiner_ &&
   combiner_->numReportedLocations() > 1) ? combiner_->combinedTestFits(forests_)
   : forest.totalTestFits.data(), then the existing loop writes
   scale*testSrc[loc*nTest+i]+shift - byte-identical at L = 1. Surface chain
   testFitsAreDefined() ([[chain.hpp:438@6a48351b]] neighbor) and a SamplerBase override
   ([[facade.hpp:145@6a48351b]]/339 neighbor).

(b) BRIDGE. createMultinomialHolder ([[R_interface_bartcore.cpp:1667-1729@6a48351b]]): after
   sampler creation, wire test predictors when data.numTestObservations > 0
   (mirror 1552-1565's dense/mixed setTestPredictors; refuse a leaf-covariate
   sparse test column is moot - multinomial is constant-leaf). Multinomial
   carries no offset, so do NOT wire testOffset; if data.testOffset != NULL,
   Rf_error. Change refuseBCFTestSurface (1403-1409) to gate on
   sampler.testFitsAreDefined() (refuse only numForests>=2 AND test-undefined) -
   BCF still refused, multinomial and single-forest allowed; keep the BCF-worded
   message but reached only by BCF. No run-allocation change (verified
   K-test-shaped already).

(c) INTERNAL R SURFACE (R/bartcore.R). bartcoreMultinomialSampler (556-580):
   the data object already carries x.test from dbarts()/the data spec, so
   createMultinomial receives it through sampler$data - confirm the K-forest
   engine picks it up, else pass through explicitly. bartcoreSetTestPredictor
   (672-675) and bartcorePredict (780-786) already route to the shared .Call
   entries; they now succeed on a multinomial bc (the refusal moved to
   test-defined). Correct the 555 comment.

(d) BART2 SURFACE (R/bart.R, R/generics.R). Remove the test= refusal
   ([[bart.R:429-434@6a48351b]]); bart2Multinomial (666-714) threads test through the host
   sampler creation (the host already accepts test=; the label-response host
   carries x.test) and bartcoreRun reports samples$test. packageMultinomial-
   Results (731-761): reshape samples$test the same way as samples$train, name
   levels on the trailing K margin, add yhat.test to the result. extract.bart-
   Multinomial: remove the sample="test" refusal ([[generics.R:407-412@6a48351b]]), return
   yhat.test; keep the type="bart" refusal. fitted stays train-only (binary
   convention). print gains a test-rows line if present.

(e) EXACT-GATE ARM (benchmarks/R/multinomial-exact.R). Add a test arm reusing
   arm 3's covariate K=3 setup with test rows duplicating the two cells: the
   test-channel posterior-mean probabilities must equal the train/quadrature
   target per cell to the same MC tolerance. This is the only gate that the
   combined-TEST-fits blend equals the (already-gated) train blend.

(f) FIXTURE TEST CHANNEL (benchmarks/R/multinomial-equivalence.R). Give each
   scenario an x.test (a held-out slice), record result$test (n_test x K x
   n.samples) in recordChannels (63-75). Re-record the baseline as
   multinomial-equivalence-<C1-hash>.rds; add a MANIFEST entry (mark
   bb8855e historical, superseded); the 5-channel identical() compare guards it.
   NEUTRALITY: the three pre-existing channels (train, forestFits, varcount)
   reproduce bb8855e bitwise (adding x.test consumes no rng, moves no train
   draw) - RECORD that as the per-scenario proof; only the test channel is new.

(g) TINYTEST (inst/tinytest/test-multinomial-surface.R): a bart2(family=
   "multinomial", test=x.test) fit produces a levels-named n x K test array;
   a public-fit test channel reproduces the internal bartcoreRun test channel
   bit-for-bit at a fixed seed (the C7 reproduction-gate pattern extended to
   test); extract(sample="test") returns it; the K=2 case.

(h) COMPONENT TEST (tests/cpp/test_sampler.cpp, beside testMultinomial:1843):
   combinedTestFits softmax positivity + sum-to-one + the K=2 reduction, and
   that it equals combinedFits when test rows equal train rows.

Files: combiner.hpp, chain.hpp, facade.hpp, sampler.hpp,
R_interface_bartcore.cpp, R/bartcore.R, R/bart.R, R/generics.R,
tests/cpp/test_sampler.cpp, benchmarks/R/multinomial-exact.R,
benchmarks/R/multinomial-equivalence.R + re-recorded baseline, MANIFEST,
inst/tinytest/test-multinomial-surface.R. Gate: all THREE anchors identical
(single-forest 22/22, BCF 5x6, multinomial's 3 existing channels) + the new
multinomial test channel identical to itself + the exact gate's four arms to MC
error + tests/cpp from make clean + full tinytest (grows; no regen - RNG-neutral
so no snapshot shifts) + air. Size: L. --preclean; delete tests/cpp binaries.
Abort: any anchor or any existing multinomial channel moves = the seam touched a
draw path (it must not).

### C1 landing (2026-07-16)

C1 = bcefa63. The combinedTestFits blend landed exactly as planned; both
bitwise anchors and the old multinomial baseline's three channels identical
(the neutrality claim held: supplying x.test moves no train draw), the new
four-channel baseline recorded as multinomial-equivalence-bcefa63.rds with
bb8855e demoted, the exact gate's fourth arm pins the test channel to the
train channel's quadrature target on duplicated cells (bit-for-cell equal),
tinytest 2916 (2905 + 11, zero regen), cpp clean, air clean. Also swept the
C7-era docs/plans references out of R/bart.R and R/generics.R comments.
Deviations recorded from the implementer: the fixture guards result$test
rather than per-forest test fits (no per-forest test accessor exists), and
the multinomial creation refuses test offsets and mixed/sparse test stores
before rng creation (dense-only makes the mixed branch unreachable).

## C2 - Predict-on-newdata (the harder half; separable, RECOMMENDED with a defer-abort)

RNG-NEUTRAL (replay + keepTrees serialization consume no draw). Lands after C1.

Sub-parts:

(a) ENGINE. Generalize the saved-tree replay to K forests. predictFromSaved-
   Sample ([[chain.hpp:1389-1408@6a48351b]]) and predictFromCurrentTrees (1413-) loop
   forests_[f] instead of forests_[0], accumulating each forest's per-row total
   into a K-location-major slab, then apply the combiner's softmax (reuse
   combinedTestFits' map on the replayed totals, or a small predict-side softmax
   over the K totals). sampler.predict ([[sampler.hpp:491-504@6a48351b]]) sizes the out slab
   nTest x K x capacity (x nChains) for a multi-location sampler; single-location
   samplers keep the exact current shape and code path. The level-shift is absent
   from saved leaves but harmless (softmax invariance, above).

(b) BRIDGE. bartcore_predict ([[R_interface_bartcore.cpp:3262-3314@6a48351b]]): allocate the
   K-location shape when numReportedLocations > 1 (an allocFitsArray-style
   insert of the K dimension), refuse per the testFitsAreDefined gate (BCF still
   refused). refuseBCFTestSurface already relaxed in C1.

(c) R SURFACE. Un-refuse keepTrees in the bart2 multinomial branch
   ([[bart.R:452-457@6a48351b]]) - storeSavedTreeRecord already saves all K forests, so
   keepTrees just needs enabling and the K-forest replay from (a).
   predict.bartMultinomial ([[generics.R:456-461@6a48351b]]): require keepTrees at fit
   (mirror predict.bart's guard, [[generics.R:153@6a48351b]]), call the K-forest predict, and
   return a levels-named n x K probability array (per-chain shape matching the
   binary predict convention). bartcorePredict ([[R/bartcore.R:780-786@6a48351b]]) needs no
   change beyond the K-shaped .Call return.

(d) GATE. A reproduction test: predict on the fit-time test rows equals the run
   test channel from C1 to the bit (both replay the same trees/softmax). A
   tests/cpp K-forest replay case. All three anchors + the C1-updated
   multinomial fixture identical (the fixture does not exercise predict/keepTrees
   so it is inert here - state that).

DEFER-ABORT (present the fork honestly): predict is genuinely more than C1. If
the K-forest replay entangles the recently-hardened getTrees / saved-tree-count
ABI (a73ca50) or keepTrees-for-K-forests beyond a bounded commit, DEFER predict
to its own follow-up: C1 (test-at-creation) is a complete, shippable deliverable
that matches the mbart2 convention (mbart2 ships prob.test, has NO predict
method). See Q1.

Files: chain.hpp, sampler.hpp, facade.hpp, R_interface_bartcore.cpp, R/bart.R,
R/generics.R, tests/cpp/test_sampler.cpp, inst/tinytest. Gate: predict==test
reproduction to the bit + cpp replay case + all three anchors identical + full
tinytest + air + rchk note (bartcore_predict shape change on a new path). Size:
L. Abort/defer per above.

### C2 landing (2026-07-16)

C2 = 88ffe12. The defer-abort was NOT taken: the K-forest saved-tree replay
stayed a bounded generalization rather than entangling the just-hardened
getTrees/saved-tree-count ABI (a73ca50) - predictFromSavedSample and
predictFromCurrentTrees loop all K forests into a location-major slab and
softmax it through the same map combinedTestFits uses, sampler.predict sizes
the K-location out slab, the bridge allocates it, and keepTrees/predict.bart-
Multinomial are un-refused on the R surface. Implementation surfaced a real
bug along the way, not just a generalization: the multinomial chain
constructor never called the saved-tree storage initialization every other
keepTrees-capable constructor calls, so keepTrees on a multinomial fit
silently allocated no storage - fixed in this commit. The strongest gate
evidence: predicting at the fit-time test matrix reproduces C1's run test
channel bitwise (both replay/report the same saved per-forest leaf fits
through the same softmax), which the tinytest reproduction check pins;
single-forest and BCF replay stayed byte-for-byte untouched. All three
anchors identical, tests/cpp gained the K-forest replay case, tinytest 2924
(2916 + 8, zero regen), air clean.

### Arc close

C1 (bcefa63) and C2 (88ffe12) deliver Q1's recommended staged BOTH: test-at-
creation and predict-on-newdata are both landed, RNG-neutral throughout,
every anchor identical across both commits. C3 (docs) closes the arc - no
further commits queued here.

## C3 - Docs

Update docs/design/multinomial.md "The surface", the "Test fits are DEFINED"
bullet: it now DEFINES test fits (combinedTestFits) and (if C2 lands)
out-of-sample predict; correct the "No probit path" section's neighbor claims
only if touched. Record the level-shift invariance as the correctness fact
licensing the untouched totalTestFits. Update docs/plans/archive/multinomial.md C7
landing follow-up list (predict/test now landed off it). Update
docs/plans/archive/multi-forest-models.md. This plan's own landing notes per commit.
Files: docs/design/multinomial.md, docs/plans/*. Gate: full tinytest; R CMD
check man unaffected (no man/ topic added - predict/extract/fitted are existing
topics). Size: S.

## Verification

- Every commit, all three standing anchors: equivalence 22/22 IDENTICAL vs
  equivalence-ac6ec2c.rds, bcf-equivalence 5x6 IDENTICAL vs
  bcf-equivalence-99205ee.rds, and the multinomial fixture's EXISTING train/
  forestFits/varcount channels IDENTICAL (vs bb8855e for C1's neutrality proof,
  then vs the C1 re-recording thereafter). A moved draw on any anchor means a
  reporting/replay path touched a sampled quantity - stop.
- C1 additionally: the multinomial fixture's NEW test channel identical to
  itself across the compare, and the exact gate's four arms (three train + the
  new test arm) to MC error.
- C2 additionally: predict-on-newdata reproduces C1's run test channel bit-for-
  bit on the same fit, and the tests/cpp K-forest replay case.
- No bench-sampler run needed: this arc adds no per-sweep work (combinedTest-
  Fits and predict run once per stored sample / once per predict call, off the
  hot path; the single-forest and BCF sweeps are untouched). Note it and skip
  the quiet-window bench unless a reviewer disputes the off-hot-path claim.
- dbarts.h unchanged -> no stan4bart lockstep; C1's relaxed setTestPredictor
  entry and C2's predict shape earn a "rchk on next scheduled run" note.

## Open questions for VD

- Q1 (arc scope: test-at-creation only, or + predict-on-newdata). RECOMMEND
  BOTH, staged (C1 then C2), with C2 defer-abortable. The ecosystem fork: mbart2
  ships prob.test but has NO predict method - test data at fit is the
  multinomial-BART convention, which makes C1 alone a complete, precedent-
  matching deliverable and argues C2 is optional. AGAINST that: dbarts' OWN
  bart2 ships predict() (keepTrees) for every other family, so refusing it for
  multinomial is a dbarts-internal asymmetry a user meets directly, not just an
  mbart2 gap; and the raw material (per-forest saved trees) already exists, so
  C2 is a bounded replay generalization, not new machinery. The cost that could
  flip it: C2 touches the just-hardened getTrees/saved-tree-count ABI (a73ca50)
  and keepTrees-for-K-forests, which if larger than a bounded commit is better
  as its own arc - hence the defer-abort rather than a hard commitment. C1 ships
  regardless.
- Q2 (fixture: re-record vs a second baseline). RECOMMEND re-record
  multinomial-equivalence with the test channel to one new baseline
  (multinomial-equivalence-<C1-hash>.rds), the three old channels proving
  bitwise inside the same compare; the alternative (a separate test-only
  fixture) duplicates the harness for no added guard. State the neutrality claim
  in the MANIFEST entry as bcf-equivalence does.
