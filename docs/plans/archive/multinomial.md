# multinomial

agent: opus (the first multi-forest MODEL. The softmax coupling is an
  INTERLEAVED one-vs-rest Polya-Gamma cycle - omega_k drawn against the CURRENT
  margins immediately before forest k's tree update - plus a likelihood-invariant
  level-centering move; the K-forest combiner, the response family, and the
  exact-posterior gate are all posterior-defining. The Phase-A seams are
  RNG-neutral relocations but touch the sweep, storeSample, and the run bridge,
  where a leak is a silent posterior change).
rng: TWO CLASSES, cleanly gated. Phase A (C1-C3): NEUTRAL - the combiner.hpp
  extraction, the combined-output/reporting widening to K locations, and the
  multi-forest mutation guard, every step bitwise on BOTH anchors (single-forest
  equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds + the BCF fixture's five
  scenarios x six channels IDENTICAL vs bcf-equivalence-99205ee.rds). Phase B
  (C4-C6): POSTERIOR-CHANGING (new model) - the model, its response family, its
  creation/reporting surface, its component tests, its exact benchmark, and its
  bitwise fixture land as ONE gated commit (C4; nothing model-bearing lands
  before its gate exists), with every pre-existing path still bitwise on both
  anchors.
window: ENGAGED by the resolved Q5 - the public surface ships (C7), so the
  family slot lands on dbartsModel per family-on-model.md, taking that
  pre-release decision. The engine work (C1-C6) is internal and takes no
  window of its own.
budget: ~1400-1800 lines. Phase A ~450 (combiner.hpp motion + co-movers + the
  location-count seam + the run-bridge allocation + the mutation guard). Phase B
  ~950-1350 (MultinomialForestCombiner with the interleaved PG hook and the
  centering move, the response family, the creation surface, reporting, the exact
  gate with three arms + component tests + the mandatory bitwise fixture).
  Chiefly src/bartcore/combiner.hpp (new), chain.hpp, model.hpp,
  R_interface_bartcore.cpp, src/C_interface.cpp, R/bartcore.R,
  benchmarks/R/multinomial-exact.R + baseline, a bitwise multinomial fixture +
  baseline, tests/cpp/test_sampler.cpp. Header edits -> --preclean; delete stale
  tests/cpp binaries.

## Status

LANDED. All seven steps closed: C1-C3 (Phase A, RNG-neutral) on 2026-07-14,
C4-C7 (Phase B, the model itself plus its public surface) on 2026-07-15, final
commit bb29d00. Two follow-up arcs grew directly out of scope narrowings
recorded at close: docs/plans/archive/multinomial-test-surface.md and
docs/plans/archive/multinomial-varcounts.md.

## Summary

This arc built dbarts' first K-forest sampler: a K-category multinomial
classifier where every category gets its own forest, the K forests are
coupled through a softmax likelihood via an interleaved one-vs-rest
Polya-Gamma augmentation, and a likelihood-invariant level-centering move pins
the resulting non-identified flat direction. It rides the forest combiner BCF
introduced, extracting Forest<L>/ForestCombiner<L> into their own header
(combiner.hpp) as the second combiner instance was always meant to
(docs/design/forest-combiner.md). The work split into RNG-neutral
seam-widening (Phase A, C1-C3) that touches every existing single-forest and
BCF code path but changes no draw, and the posterior-changing model itself
(Phase B, C4-C7), landed as one gated commit plus a public bart2 surface.
Every phase gates on two bitwise-equivalence anchors (single-forest, BCF) plus
an exact-posterior benchmark with three arms; design forks opened while
drafting the plan (Q1-Q5) were resolved by VD before Phase B began.

## Goal

An internal K-category multinomial sampler: K forests coupled through a softmax
(log-linear) likelihood with an interleaved one-vs-rest Polya-Gamma augmentation,
riding the forest combiner the BCF refactor introduced. The combiner is BCF's
second instance and the occasion to extract Forest<L>+ForestCombiner<L> to
combiner.hpp (the deferral recorded in docs/design/forest-combiner.md). Chain
gains a location-count on its combined-output/reporting seams (1 for every
existing model, K for multinomial) so the widening is bitwise-neutral on the
single-forest and BCF paths. Creation mirrors BCF's internal .Call surface;
dbarts.h stays frozen.

Five forks were opened while drafting this plan. VD resolved Q1 (augmentation),
Q2 (identification), Q3 (output widening), and Q5 (public surface) together on
2026-07-14 by an "expectations rule" - two independent ecosystem surveys (BART
packages; general R multinomial tooling), with expectations winning absent a
performance issue. Q4 (state wire format) was resolved separately that same
day, on pre-release-compatibility grounds rather than the survey. Each fork is
discussed where it arises below: Q1 and Q2 under Model choice, Q3 under Key
finding, Q4 at step C5, Q5 at step C7.

## Model choice (the settled literature; defer to it per VD's standing guidance)

- Q1 - AUGMENTATION: log-linear / softmax with Polya-Gamma, one forest per
  category. Given the K forest fits f_ik, P(y_i = k) = softmax(f_i)_k. Category
  k's ONE-VS-REST conditional is a binomial logistic with linear predictor
  eta_ik = f_ik - C_ik, where C_ik = log sum_{j != k} exp(f_ij) is the
  log-sum-exp margin; so omega_ik ~ PG(n_i, eta_ik) and forest k sees working
  response (y_ik - n_i/2)/omega_ik + C_ik under precision omega_ik - the EXACT
  logistic PG machinery already shipped, one binomial per category.
  CRITICAL - THE AUGMENTATION IS INTERLEAVED, NOT JOINT. The product of K
  one-vs-rest kernels integrates to a product of binomials, NOT the multinomial:
  each omega_k is a TEMPORARY latent valid only for category k's conditional. So
  omega_k must be drawn against the CURRENT margins C_ik immediately BEFORE forest
  k's tree update, cycling the categories (Held and Holmes 2006; Polson-Scott-
  Windle 2013 sec 4; the BayesLogit mlogit interleaved one-vs-rest cycle). Drawing
  all K omegas in a single post-loop draw (as an earlier draft proposed, in the
  combiner's drawGlue) is a Jacobi-style simultaneous update that does NOT target
  the posterior: forest k>0 would get a working-response MEAN from fresh C_ik but
  a PRECISION drawn against stale margins.

  The fork: PG log-linear softmax reuses the shipped PG(1, psi) sampler
  (ext_rng_simulatePolyaGamma) and the weighted-conjugate kernels verbatim - one
  binomial-logistic per category, drawn INTERLEAVED via the per-forest hook -
  costing no new numerical code beyond that hook; its price is the softmax
  non-identification (handled like BCF's a-sign, comparing probabilities) plus
  the interleaving discipline. Multinomial probit (Kindo et al. 2016) is
  identified via a reference category but needs truncated multivariate-normal
  latent utilities and a covariance sampler - a whole new latent + move
  machinery reusing none of the PG code, and a heavier exact gate.

  ATTRIBUTION: the algorithm is the Held-Holmes / PSW interleaved one-vs-rest
  conditional cycle; Murray 2021 ("Log-Linear Bayesian Additive Regression
  Trees", JASA) is RELATED WORK, not the source - its construction uses a
  per-observation normalizer latent, a different device.

  RESOLVED 2026-07-14 (VD): PG softmax, interleaved. The ecosystem does surface
  probit-vs-logit as a user choice (mbart's pbart/lbart), but multinomial
  probit needs truncated-MVN latents plus a covariance sampler - real new
  machinery - and the modern many-category precedent (mbart2) defaults to
  logit. The design note states plainly that no probit path exists: a one-way
  door versus the mbart convention, taken deliberately.

- Q2 - IDENTIFICATION: K SYMMETRIC forests (softmax is invariant to a
  per-observation constant shift of all f_ik, so the raw f_ik are
  non-identified - exactly BCF's a-sign situation). K gives every category its
  own forest and its own variable-count channel (symmetric reporting), and the
  non-identification is handled the way BCF's is - the gate and the reporting
  compare IDENTIFIED functionals (softmax probabilities), never the raw f_ik.
  K-1 with a reference category removes the invariance (and the need for a
  centering move) but makes the reference category's reporting implicit and
  the varcount asymmetric.

  THE FLAT DIRECTION NEEDS AN EXPLICIT MOVE: the likelihood is invariant to a
  per-observation common shift of all f_ik, a flat additive direction pinned
  ONLY by the prior, which mixes as a slow random walk under K symmetric
  forests. ADD a likelihood-invariant LEVEL-CENTERING Gibbs move (the
  BCF-ridge-interweave analog) in afterCombine: draw the per-observation
  common shift from its Gaussian full conditional under the forest priors and
  recenter f_ik.

  RESOLVED 2026-07-14 (VD): K symmetric forests, with the centering move.
  mbart2, the closest precedent, is fully symmetric with K-length
  varcount/treedraws; mpbart's user-facing base= reference argument is the
  recorded counterexample wart; BART has no coefficient table for a reference
  category to hide in, and the user-visible K-shaped objects (probabilities,
  per-category variable counts) are symmetric in every surveyed package. The
  reporting symmetry is worth the one added centering move, which the gate and
  R surface already sidestep by reporting probabilities.

  AS BUILT (see the C4 landing note below): the centering move ended up being
  a GLOBAL scalar shift, not the per-observation shift described above - an
  adversarial review caught that a per-observation shift is not representable
  by shared-leaf trees, so it would bias the identified probabilities.

- SCOPE (M5) - SINGLE-TRIAL this arc (n_i = 1 classification). The one-vs-rest
  binomial for category k needs the SUCCESS COUNT y_ik, not just a category
  label; borrowed integer labels support only n_i = 1. Scope this arc to
  single-trial classification and record the grouped-count generalization
  (store the n x K count matrix, n_i > 1) as a follow-up - PG(n_i, .) is the
  sum of n_i PG(1, .) draws (weighted-logistic.md), so grouped counts add no
  numerical code, only the count-matrix data model, and the exact gate below
  is single-trial anyway. Real-shape (non-integer) counts stay out of scope
  (the deferral weighted-binary/negative-binomial carry).

## Context (seams, read in code)

- The combiner landed in chain.hpp beside Forest<L>: ForestCombiner<L> base
  (chain.hpp:394-439) and BCFForestCombiner<L> (447-...). The virtuals a second
  combiner needs already exist: formForestResponse (per-forest (response,
  weights) over the whole forest vector, called per forest INSIDE the sweep loop
  at chain.hpp:963 with the partially-updated forests_), combinedFits (405),
  drawGlue (421) + afterCombine (423, called post-loop at 1067-1069), the
  reporting map reportedForest/testFitsAreDefined/logLikelihoodIsDefined
  (431-433), serializeGlue/restoreGlue (437-438), setTreatment/bcfGlue (410-415).
  The per-forest call site (963) is exactly where the interleaved omega_k draw
  belongs; extraction to combiner.hpp is the deferred decision
  (docs/design/forest-combiner.md "Header location").
- The single n-vector combined output, the seam the design note flags
  (docs/design/forest-combiner.md "What still re-carves"): Chain::combinedFits()
  returns const double* (chain.hpp:2498, single-forest arm 2499-2500 returns the
  bare forests_[0].totalFits.data() pointer); storeSample writes ONE trainingFits
  and ONE testFits channel (2518-2554); the run bridge allocates trainFits as
  n x numSamples (R_interface_bartcore.cpp:1987) and testFits as nTest x numSamples
  (1996). refreshLatents/drawSigma take one location (chain.hpp:1054-1065) - but
  see the finding below: multinomial does NOT widen those.
- Per-forest fits already reach K forests generically: getForestFits(forest) ->
  forestTotalFits(chain, forestIndex, out) validates forestIndex <
  numForests() (R_interface_bartcore.cpp:1924-1938; facade.hpp:146). The
  BCF-shaped reporting (NaN-flag the scalar train/test channels, expose per-forest
  fits) is the precedent multinomial reuses (chain.hpp:2537-2554,
  bcf-testfits-guard).
- The PG sampler: ext_rng_simulatePolyaGamma(rng, psi) draws PG(1, psi)
  (Devroye, src/external/random.c:581); LogisticResponse::refreshLatents sums
  lround(w_i) draws and forms working response (y-0.5)/omega - offset
  (model.hpp:2231-2244). MultinomialForestCombiner's per-category draw is this,
  with the log-sum-exp margin C_ik as the per-category "offset" - but drawn
  INTERLEAVED, one category at a time against the current margins immediately
  before that category's forest update (M1; Q1 above), never batched post-loop.
- Creation mirrors BCF: R/bartcore.R:482 bartcoreBCFSampler packs a params vector
  + moderators, calls C_dbarts_bartcore_createBCF (R_interface_bartcore.cpp:1879),
  which builds a BCFSpec (chain.hpp:348) and calls createBCFSampler
  (facade.hpp:499, single ConstantGaussianLeaf instantiation). The multinomial
  analog: bartcoreMultinomialSampler + C_dbarts_bartcore_createMultinomial +
  MultinomialSpec + createMultinomialSampler.
- ResponseModel is per-observation-single-location (model.hpp:1788); multinomial
  has no sigma (fixed, like the binary families) and K latent channels, so its
  response family is thin (Q3, discussed below): a MultinomialResponse whose
  single-location seams are vestigial (no-op refreshLatents/drawSigma) and whose
  log-likelihood channel is NaN-flagged - storeSample scores one forest's fits via
  computeLogLikelihood (chain.hpp:2586-2593) and cannot see the K-blend, so
  logLikelihoodIsDefined() = false is BCF's exact choice (chain.hpp:638). The
  combiner owns the K interleaved PG draws, the per-forest working responses, and
  the level-centering move.
- State: the by-name block reader (R_interface_bartcore.cpp:3317
  stateFormatVersion = 3, minReadable = 3) already serializes a per-forest tree
  list of any length (SLOT_FORESTS), so K forests serialize with no format work;
  only combiner-specific scalar/latent state would need a block (Q4, resolved
  at step C5 below).

## Key finding (contradicts the design note's re-carve framing; resolves Q3)

docs/design/forest-combiner.md lumps three seams as forced to widen for a
non-BCF combiner: combinedFits -> const double*, refreshLatents/drawSigma take
one location, results.trainingFits is one channel. Read in code, only the
REPORTING/combined-OUTPUT seam genuinely widens, and only if we want the engine
to report probabilities:

- refreshLatents/drawSigma do NOT widen. The K-category PG augmentation lives in
  the combiner's INTERLEAVED per-forest hook (M1; Q1 above) - each omega_k drawn
  against the current margins just before forest k's update - and multinomial
  has no sigma. So MultinomialResponse's refreshLatents/drawSigma are no-ops
  and the single-location `combined` they are handed is ignored - no signature
  change, and every existing family is untouched. GUARD (M6): the K x n
  combinedFits buffer is handed to the no-op refreshLatents as a single-location
  pointer; carry a guard note (assert / comment) so a future non-no-op
  MultinomialResponse::refreshLatents cannot silently misread a K-channel
  buffer as one channel.
- combinedFits + trainingFits/testFits widen to K ONLY under the
  engine-reports-probabilities choice. The minimal alternative keeps them at
  one location, NaN-flags the scalar channels exactly as BCF does, exposes the
  K forests via getForestFits (already K-capable), and lets the R surface
  compute softmax probabilities. That path needs NO combined-output widening
  at all - the n x K softmax object stays internal to
  MultinomialForestCombiner (it reads forests_ directly for the margins) - and
  under it Phase A would shrink to just C1 and C3.

This is Q3: where the K-location generalization stops. PRIMARY: widen the
combined-output seam to a location count (C2) so the engine reports the K
softmax PROBABILITIES directly in trainingFits/testFits - identified, and the
natural deliverable for a classification model, where the raw per-category
fits are non-identified nuisance a user should not have to softmax by hand;
the widening is also the reusable multi-location channel
heteroscedastic/hurdle will want later, and it stays fully bitwise on both
anchors because L = 1 everywhere today. MINIMAL: the alternative just above -
NaN-flag the scalar train/test channels, compute softmax R-side, and skip the
combined-output widening entirely.

RESOLVED 2026-07-14 (VD): PRIMARY - the engine reports the K softmax
probabilities directly. Unanimous across the survey: BART::mbart2 returns
prob.train/prob.test directly (upgraded from user-side softmax in its own
history), every surveyed general package returns a labeled n x K probability
object, and dbarts' binary generics already default to probability scale
(type = "ev"). No performance issue: the L = 1 widening is byte-neutral and
the softmax is O(nK) per stored sample.

## Constraints

- dbarts.h FROZEN. Multinomial is internal-only (bartcore .Call), like BCF; the
  bartcore::Results widening (Phase A) is the internal struct, NOT the public
  dbarts_results the LinkingTo consumers (stan4bart) call - no header change, no
  lockstep, no rchk.
- Single-forest AND BCF hot paths pay nothing new. The location-count is 1 for
  every existing model, combiner_ stays null off multi-forest; every Phase-A
  touchpoint keeps its `if (combiner_)` shape, and Chain::combinedFits()'s
  single-forest arm (chain.hpp:2498-2501) stays the bare forests_[0].totalFits
  pointer (no scratch write, no memcpy - the widening is the combiner VIRTUAL
  only, G8). Equivalence 22/22 IDENTICAL and the BCF fixture IDENTICAL are the
  per-step proof, both anchors, every commit.
- Per-forest-per-SWEEP virtual dispatch only (core-generalization.md). The
  interleaved omega_k draw is one O(n) hook per forest per sweep (K hooks/sweep),
  the margin/centering loops O(n*K) inside afterCombine once per sweep, never in
  the per-observation partition/suffstat kernels.
- fast over safe in C/C++; header-only C++20; Doxygen/LLVM/ASCII; the
  comment/code delta bounded (README brevity).
- OUT of scope: heteroscedastic and hurdle (their own notes); grouped counts
  (n_i > 1, the n x K count-matrix follow-up above) and real-shape counts; a
  per-observation sigma; ordinal outcomes (ordinal-outcomes.md, cumulative
  probit - a different model); grouped-x-multinomial composition (the combiner
  makes it expressible via GroupedResponse below the combiner, but this arc does
  not build it).

## Phase A - RNG-neutral seams (C1-C3)

Every Phase-A step is gated: equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c
+ the BCF fixture five-scenarios-x-six-channels IDENTICAL vs
bcf-equivalence-99205ee + full tinytest 2832 no regen + tests/cpp (delete stale
binaries after header edits). If any step cannot be made bitwise on both
anchors, the seam design is wrong - stop and rethink.

C1. Extract combiner.hpp (pure motion; NEUTRAL). Move Forest<L>, ForestResponse,
   ForestCombiner<L>, and BCFForestCombiner<L> from chain.hpp into a new
   src/bartcore/combiner.hpp; chain.hpp includes it. The move set is BIGGER than a
   naive read suggests (size it honestly): ForestCombiner::serializeGlue/
   restoreGlue take ChainStateData& (chain.hpp:247), BCFForestCombiner holds
   BCFState (361) and is built from BCFSpec (348), and ForestStateData (224) and
   BCFForestSpec (333) are the other co-movers - name them or resolve the ordering
   with a mid-file include so combiner.hpp sees the state/spec structs. Pure
   motion, caught at compile time. This is the second-consumer extraction
   docs/design/forest-combiner.md deferred - do it FIRST so C4 adds
   MultinomialForestCombiner beside BCFForestCombiner in the new home. Note:
   tests/cpp's Makefile HAS -MMD/-MP dependency tracking (tests/cpp/Makefile:8) -
   binary deletion is belt-and-suspenders. Files: src/bartcore/combiner.hpp (new),
   chain.hpp. Gate: NEUTRAL - both anchors identical + tinytest + tests/cpp. Size:
   M. --preclean; delete tests/cpp binaries.

C2. Location-count on the combined-output/reporting seam + the multi-forest
   mutation guard (NEUTRAL; 1 location everywhere today). Give the combiner a
   numReportedLocations() (default 1; BCF 1) and widen the combinedFits VIRTUAL to
   write numReportedLocations x n into caller-owned scratch (BCF keeps writing one
   channel) - but Chain::combinedFits()'s single-forest arm stays the bare
   forests_[0].totalFits.data() pointer, no scratch write, no memcpy (G8: the
   widening is the combiner virtual only). storeSample's trainingFits/testFits
   loop over the reported locations (one iteration off multinomial).
   bartcore::Results grows a location stride the run bridge reads (Q3-primary); the
   bridge allocates trainFits as n x L x numSamples and testFits as
   nTest x L x numSamples with L = numReportedLocations (L = 1 -> the exact current
   n x numSamples shape and byte layout), and sampler.hpp's chain-major stride
   collapses at L = 1. Add C_interface.cpp (dbarts_sampler_run, the Results
   consumer, the stan4bart path - covered by test-capi.R) to the file list.
   MUTATION GUARD (G5): add refuseMultiForestMutation to bartcore_setData
   (R_interface_bartcore.cpp:2218, which guards only refuseViewSampler at 2221),
   bartcore_setPredictor (2557/2560), setResponse, and setWeights (setTreatment
   stays allowed) - applyNewData/revalidateTrees/rebuildFitsFromParameters rebuild
   forests_[0] ONLY, so whole-data mutation on a numForests >= 2 sampler silently
   corrupts the other forests; this closes a PRE-EXISTING BCF gap multinomial would
   widen. ANCHOR-COVERAGE CAVEAT (state honestly): equivalence reads only
   test/varcount/sigma/k as summaries and never reads train, and its poolChains
   tolerates reshapes, so tinytest (NOT the anchors) is what guards the train/test
   OUTPUT shapes at L = 1. Files: combiner.hpp, chain.hpp (combinedFits virtual,
   storeSample, Results), facade.hpp, sampler.hpp, src/C_interface.cpp,
   R_interface_bartcore.cpp (bartcore_run allocation + the mutation guards). Gate:
   NEUTRAL - both anchors identical (the L=1 layout is unchanged) + tinytest +
   tests/cpp + test-capi.R. Size: L. Abort: any equivalence or BCF-fixture channel
   moves = the L=1 path is not byte-identical.

C3. Per-forest reporting query for K forests (NEUTRAL; BCF-inert; foldable into
   C1/C2 per G1). Add a per-forest variable-count query beside getForestFits
   (getForestVariableCounts(forest) -> forestVariableCounts(chain, forestIndex,
   out)), so a K-forest model can report each category's split usage without
   widening the single storeSample varcount channel (which keeps addressing
   reportedForest()). BCF exposes forest 0 as today; the new query is additive and
   unused off multinomial. Files: facade.hpp, sampler.hpp, chain.hpp,
   R_interface_bartcore.cpp, R/bartcore.R. Gate: NEUTRAL - both anchors identical +
   tinytest + tests/cpp. Size: S.

### Phase A landing (2026-07-14)

C1 = 300e6c4 (combiner.hpp extraction, pure motion + co-movers), C2 = 9030d93
(location-count seam, byte-identical at L = 1; multi-forest mutation guards
hazard-exact - blanket on setData/setResponse/setWeights and the
per-observation sessions, non-force-only on setPredictor/updatePredictor, the
force paths kept per VD's "support whatever was supported"), C3 = 7c83419
(per-forest variable-count query, storeSample counting factored byte-neutrally).
Every commit gated: both anchors identical, full tinytest (2865 at C3),
tests/cpp from make clean.

Phase-A-close bench (G8), quiet window 2026-07-14, at 7c83419 vs
bench-sampler-4008675.csv: run-n10000-p10-t75 flagged at 1.066 (re-run 1.073;
the re-run read uniformly ~1.5% warmer across every arm including the
untouched setPredictor-reject path - thermal). All other arms 0.974-1.050.
CONTROL: the pre-Phase-A tip f859b2a, rebuilt and benched in the same window
on the same machine, shows the SAME flag at 1.056 with the same profile.
Differenced against the control, Phase A's arm-by-arm deltas are <= ~1.2%,
mixed sign: Phase A is NEUTRAL by measurement, the G8 claim holds. The
absolute flag decomposes as the accepted forest-combiner run-arm cost
(recorded 1.018-1.043 at its close, zero flags then) plus machine drift since
the 4008675 recording; their sum now trips the 5% gate on the largest arm.
The 4008675 baseline predates the forest-combiner close and will keep
tripping n = 10000: re-record it in a COLD quiet window (not this one - the
warm-machine runs above would bake a generous reference that masks small
regressions) before it gates another arc. DISCHARGED 2026-07-15:
bench-sampler-60a13b6.csv recorded cold at the arc close; its values sit in
the band of the four same-code runs above, so the drift vs 4008675 was
persistent machine state, not thermal.

## Phase B - the multinomial model (C4-C6)

POSTERIOR-CHANGING. C4 is THE MODEL AS ONE GATED COMMIT: no model-bearing code
lands before the gate that proves it exists, so the engine combiner, the response
family, the creation/reporting surface, the cpp component tests, the exact
benchmark, and the bitwise fixture all land together. Equivalence 22/22 + the BCF
fixture stay IDENTICAL at every Phase-B step (the new model must not perturb any
existing path) + full tinytest.

C4. The multinomial model, one gated commit. Sub-parts:

   (a) ENGINE - MultinomialForestCombiner<L> in combiner.hpp (constant-leaf,
   static_assert as BCF): holds the borrowed single-trial category labels
   (n_i = 1 this arc; the count matrix is the M5 follow-up), the n x K PG omega
   scratch (cold-started at 1/4 per category, the logistic seed), and the
   K-channel combined scratch. THE INTERLEAVED PG DRAW (M1; Q1): omega_k is drawn
   against the CURRENT margins C_ik = log sum_{j != k} exp(f_ij) immediately BEFORE
   forest k's tree update - via a per-forest pre-update hook (either add ext_rng*
   to formForestResponse, called at forest k's turn chain.hpp:963 with the
   partially-updated forests_, OR a new virtual drawForestGlue(f, rng, forests)
   fired just before each forest's update; base no-op, BCF unchanged - BCF's glue
   draw stays post-loop in drawGlue, and BCF's a/aVariance one-sweep lag is the
   benign recorded exception, chain.hpp:610). formForestResponse(k, ...) then forms
   category k's PG working response (y_ik - 1/2)/omega_ik + C_ik under weight
   omega_ik, IGNORING the passed y (the combiner owns the labels). The hook
   addition is BITWISE-NEUTRAL (both anchors prove it - default no-op). drawGlue's
   post-loop role for multinomial shrinks or empties. THE LEVEL-CENTERING MOVE
   (M3; Q2) lives in afterCombine: draw the per-observation common shift from its
   Gaussian full conditional under the forest priors and recenter f_ik (the
   BCF-ridge-interweave analog). combinedFits writes the K softmax probabilities
   per observation (numReportedLocations() = K); reportedForest addresses per
   category via the C3 query.

   (b) RESPONSE - a thin MultinomialResponse (model.hpp): no sigma
   (sigmaIsFixed), no-op refreshLatents/drawSigma (its single-location seams are
   vestigial - the combiner owns the K PG draws), and logLikelihoodIsDefined() =
   false (G6, BCF's exact choice: storeSample scores via
   response_->computeLogLikelihood(forest.totalFits..., chain.hpp:2586-2593) - ONE
   forest's fits, which cannot see the K-blend; the earlier "carries the labels for
   computeLogLikelihood" claim is DROPPED as self-contradictory with the NaN-flag).
   Carry the M6 guard note on the K x n buffer.

   (c) SURFACE - MultinomialSpec (chain.hpp: K, the per-category forest spec, the
   leaf-scale calibration below), createMultinomialSampler (facade.hpp, single
   ConstantGaussianLeaf instantiation as BCF), createMultinomialHolder +
   C_dbarts_bartcore_createMultinomial (R_interface_bartcore.cpp: K forests from
   the model spec, labels validated as 0..K-1 integers, refuse a non-multinomial
   response), and bartcoreMultinomialSampler + getForestFits/getForestVariableCounts
   wrappers (R/bartcore.R). The run output reports the K softmax-probability
   channels (Q3-primary) or NaN-flags them and relies on getForestFits (Q3-minimal);
   refuse the single-forest test surface as refuseBCFTestSurface does.
   LEAF-SCALE CALIBRATION (M2): the logistic precedent is node.scale =
   pi*sqrt(3) ~ 5.44 (dbarts.R:401, xbart.R:262; pi/sqrt(3) ~ 1.81 is the WIDENING
   FACTOR over probit's 3.0, dbarts.R:394-395; negative-binomial.md confirms
   "the logistic pi*sqrt(3) precedent"). The IDENTIFIED quantity is a DIFFERENCE of
   forests (a pairwise log-odds f_ik - f_ij has SD s*sqrt(2) for per-forest SD s),
   so the K = 2 anchor is s = pi*sqrt(3)/sqrt(2) ~ 3.85, and C_ik's variance is
   K-DEPENDENT, so NO single per-forest scale matches binary calibration across K.
   Calibrate per-forest node.scale to the K = 2 anchor pi*sqrt(3)/sqrt(2); record
   the implied prior-on-probabilities curve for K = 2, 3, 5 as the calibration
   artifact in the design note; note the K-dependence honestly.

   (d) COMPONENT TESTS (tests/cpp/test_sampler.cpp). testMultinomial: a
   deterministic fixed-seed fit recovering a monotone K = 3 category signal; the
   interleaved PG margin/omega draw checked for positivity and the K = 2 reduction
   (a two-category multinomial fit must reproduce the logistic path's working-
   response/weights construction); PG-moment tests at the trial counts used
   (mean/variance of PG(n, psi), following weighted-logistic's component test); a
   state round-trip. MANDATORY growForestFromRoot case (G4): the multinomial
   combiner branch (formForestResponse chain.hpp:1222-1227, drawGlue/afterCombine
   1285-1287) is UNREACHABLE from R - the BCF lesson (testBCFGrowForestFromRoot)
   verbatim.

   (e) EXACT-POSTERIOR GATE (benchmarks/R/multinomial-exact.R; M4, G3), in
   bcf-exact.R / categorical-exact.R style:
   - INTERCEPT-ONLY K = 3 softmax (each forest a single leaf via a constant
     predictor -> no valid cutpoints -> root-only trees; the trick is sound but the
     implementer must confirm the engine tolerates a no-cutpoint column). The
     quadrature integrates the sampler's ACTUAL symmetric N(0, tau^2)^K prior,
     reducing dimension by marginalizing the LEVEL analytically - the induced prior
     on the identified differences (d_k = f_k - f_K, k = 1..K-1) is a CORRELATED
     Gaussian, covariance tau^2(I + 11'), NOT independent and NOT a flat simplex.
     Do NOT fix f_K = 0 (that zeroes Var(f_K), breaks exchangeability, and compares
     against the wrong posterior). Match the posterior mean of the IDENTIFIED
     softmax probabilities to the sampler's to MC error (raw f_k sign/shift-ill-posed
     by construction, compared never - BCF's a precedent).
   - K = 2 == LOGISTIC arm (REQUIRED): the multinomial fit matches the shipped
     logistic sampler DISTRIBUTIONALLY to MC error, NOT draw-for-draw (different
     rng consumption - state it so). The strongest internal-consistency check:
     K = 2 multinomial IS binary logistic.
   - COVARIATE-DEPENDENT K >= 3 arm (REQUIRED, G3): one categorical predictor,
     small cells, per-cell quadrature (categorical-exact.R style) - the ONLY exact
     gate on tree-growth-under-softmax for K >= 3 (K = 2 == logistic covers
     covariate K = 2; intercept-only covers the margin/coupling).
   Failure means the softmax coupling, the interleaved PG draw, the centering move,
   or the margin formation is wrong.

   (f) THE BITWISE MULTINOMIAL FIXTURE (MANDATORY, G2; not "if wanted"):
   bcf-equivalence-shaped, recording result$train (the K softmax-probability
   channels), per-category getForestFits, and per-category getForestVariableCounts;
   two scenarios - one K = 3 covariate-dependent + one K = 2 (logistic-equivalent).
   Rationale to record IN the fixture: equivalence byte-guards NEITHER train nor
   test shape (it reads only test/varcount/sigma/k as summaries and never reads
   train; its poolChains tolerates reshapes), and heteroscedastic/hurdle will
   re-touch exactly these seams - the fixture inherit-guards them, the
   bcf-equivalence lesson verbatim. Add to MANIFEST/equivalence.yaml and record the
   baseline.

   Files: combiner.hpp, model.hpp, chain.hpp, facade.hpp, sampler.hpp,
   R_interface_bartcore.cpp, R/bartcore.R, tests/cpp/test_sampler.cpp,
   benchmarks/R/multinomial-exact.R + baseline, the multinomial bitwise fixture +
   baseline. Gate: the exact benchmark passes to MC error (all three arms) + the
   cpp cases (incl. growForestFromRoot) + BOTH bitwise anchors still identical +
   the new multinomial fixture identical to itself + full tinytest + rchk note (new
   .Call entry points; "rchk on next scheduled run" per README). Size: XL (this is
   the model).

### C4 landing (2026-07-15)

Landed as bb8855e, one gated commit. Gates verified twice (implementer, then
an independent orchestrator re-run): equivalence 22/22 identical draws,
bcf-equivalence five scenarios x six channels identical, the exact gate's
three arms at max gaps 0.0000 / 0.0008 / 0.0012 against tolerances
0.008 / 0.008 / 0.015, tests/cpp from make clean, full tinytest 2865 all ok
with zero regenerated snapshots, air check clean. The bitwise fixture is
recorded as multinomial-equivalence-bb8855e.rds and registered in the
MANIFEST (which also gained the missing bcf-equivalence-99205ee and
bench-sampler-4008675 entries). rchk on the next scheduled run
(bartcore_createMultinomial is a new .Call entry point).

An adversarial Opus review confirmed the interleaved draw, the exact-gate
construction, the bridge lifetimes, and BOTH deviations from the plan text
as correctness-forced; C6's design note must record them:

- The LEVEL-CENTERING move (M3; Q2) is a GLOBAL scalar shift, not the
  per-observation shift the plan described. A per-observation shift is not
  representable by shared-leaf trees, so the next backfit projects the
  mismatch and biases the identified probabilities (~2-4 percent Jensen bias,
  measured against the exact gate). The global shift is the one flat
  direction a forest carries exactly (absorbed into totalFits plus tree 0's
  fit slab, keeping the residual roll consistent), moves only the
  non-identified level, and is drawn from its exact Gaussian conditional
  under the symmetric per-forest priors - valid Gibbs, mixing-only in the
  level direction.
- The K = 2 == logistic exact-gate arm is INTERCEPT-ONLY. A covariate K = 2
  multinomial log-odds is a difference of two m-tree ensembles, logistic a
  single m-tree ensemble: equal prior covariance at the sqrt(2) calibration
  but a different function-space prior, so equality is exact only at
  intercept. Covariate tree growth under softmax is gated by arm 3's
  per-cell quadrature.

C5. State round-trip.

The K forests already serialize through the existing per-forest SLOT_FORESTS
list (any length), so no format work is needed there. The remaining question
(Q4) is whether the combiner's own latents need a wire format at all.

PRIMARY: redraw the PG latents on restore, matching the "restore is
structural" contract (test-bcf.R, state-format-policy.md). serializeGlue/
restoreGlue stay no-ops for the multinomial combiner (it carries no
un-recoverable scalar glue - softmax has no a/b coefficients), omega is
cold-restarted and redrawn on the first restored sweep, and re-creation is
driven by the stored R-side creation call as BCF's is. THE INTERLEAVED DRAW
(M1; Q1) resolves the restore wrinkle for free: the first restored sweep seeds
omega_k against the RESTORED forests correctly, because omega_k is drawn just
before forest k's update against whatever margins the restored forests
present. ALTERNATIVE: serialize the n x K omega in one append-only by-name
block, bumping stateFormatVersion 3 -> 4 and holding
minReadableStateFormatVersion at 3 (additive, per c-api-growth.md's registry
rule), for a tighter restore. This CONTRADICTS the design note's blanket claim
that "a non-BCF combiner needs a format bump": softmax has nothing of its own
shape to write.

RESOLVED 2026-07-14 (VD): PRIMARY, no bump - and further, the BCF-era format
increment can be UNDONE outright: none of the bartcore work is released, so
there is no reader of any earlier version to protect. C5 therefore flattens
the state format numbering: stateFormatVersion and
minReadableStateFormatVersion collapse to one pre-release version, with the
hasBCF glue block absorbed into the base format (no versioned append history
until a release creates a compat obligation).

Files: chain.hpp (serializeGlue/restoreGlue - untouched, already no-ops),
R_interface_bartcore.cpp (no new SLOT needed). Gate: a state round-trip
component test (structural, NOT bitwise - the BCF fixture's rule) + both
anchors + tinytest. Size: S.

### C5 landing (2026-07-15)

Q4-primary held with nothing to serialize: the multinomial combiner carries
no wire glue (serializeGlue/restoreGlue untouched no-ops), omega is a
per-sweep latent the interleaved draw re-seeds against the restored forests,
and the structural K-forest round-trip landed with C4's component tests. The
remaining content was the VD-directed flatten: stateFormatVersion and
minReadableStateFormatVersion collapsed from 3 to 1 - the 1.0-0 encoding is
the first shipped format, so the pre-release development increments
(the forests-list/BCF restructure included) fold into it; a state with no
version attribute reads as 0 and still refuses at the floor. Pins updated in
test-sampler-state-format.R (current 1, future 2, below-floor 0);
test-bartcore.R's relative below-floor construction needed nothing. Gates:
both anchors and the multinomial fixture identical, tests/cpp from make
clean, tinytest 2865 all ok, air clean.

C6. Docs. Design note docs/design/multinomial.md (likelihood; the INTERLEAVED
   one-vs-rest PG cycle with the per-forest hook; the LEVEL-CENTERING move;
   identification; the prior/leaf-scale calibration curve for K = 2, 3, 5; the
   exact-gate construction incl. the correlated-difference prior tau^2(I + 11');
   the surface; the Q resolutions incl. the ecosystem-survey grounding; the
   PLAINLY-STATED no-probit-path door (mbart exposes pbart/lbart; this model is
   PG-softmax only, the deliberate performance carve-out); ATTRIBUTION -
   Held-Holmes 2006 / PSW 2013 sec 4 as the algorithm, Murray 2021 as related
   work, not the source). Update multi-forest-models.md (multinomial landed) and
   docs/design/forest-combiner.md's "Anticipated, not built" (built) + the
   combiner.hpp extraction record + the correction that the softmax combiner
   needs no wire state. Files: docs/design/multinomial.md, docs/plans/*. Gate:
   full tinytest; R CMD check man unaffected (no man/ touched). Size: M.

### C6 landing (2026-07-15)

docs/design/multinomial.md written against the landed code (both C4
deviations recorded as design facts; the calibration table computed at
K = 2, 3, 5); multi-forest-models.md's multinomial entry marked landed;
forest-combiner.md updated (second consumer built, combiner.hpp extraction,
the drawForestGlue hook on the surface, the re-carve list corrected -
refreshLatents/drawSigma did not widen and the softmax combiner serializes
no glue). Docs-only.

C7. Public surface.

Q5 - internal-only this arc, matching BCF, or a public surface now? BCF ships
internal-only (bartcoreBCFSampler, an unexported .Call path; no
family = "bcf"). Multinomial could do the same - internal creation, exercised
by the exact gate and by future consumers via the bartcore surface - or grow a
public bart2(family = "multinomial") / a multinomial front-end this arc.
Internal-only keeps the arc focused on the engine + the coupling correctness
and takes no pre-release window; a public surface adds argument design,
factor-response ingestion, a prediction/probability front-end, and Rd +
tinytest coverage, and (per family-on-model.md) the family slot lives on
dbartsModel. The original recommendation, matching BCF's own staging, was
internal-only this arc, filing the public surface as its own follow-up once
the coupling is proven.

RESOLVED 2026-07-14 (VD): the public surface ships, staged as this step (C7),
landing only after C4-C6 prove the engine. BART users expect a callable
top-level fitting function; internal-only departs from the ecosystem's notion
of "shipped," and only staging (not a performance issue) argued for
deferral - staging survives, just in a different form: C7 lands after the
engine is proven (C4), rather than the public surface being deferred to a
separate arc entirely. This ENGAGES family-on-model's pre-release window (the
family slot lands on dbartsModel).

Surface shape per the ecosystem surveys: bart2 gains family = "multinomial"
(the existing family seam). The response is a FACTOR, K inferred from
levels(y), no explicit K argument (mbart2's convention); levels(y) is captured
at fit time and threaded onto every K-shaped output for the object's lifetime
(probability array dimnames, per-category variable counts, class predictions);
fitted/predict/extract return PROBABILITY-scale K-column output under the
existing type = "ev" default with the latent scale as the escape hatch
(dbarts' binary two-tier convention); a class-prediction convenience (argmax
as a factor over the original levels); per-category variable counts through
the C3 query. Ingestion: factor-response validation, single-trial this arc
(the M5 scope). Files: R/dbarts.R, R/A_class.R (family slot), R/generics.R,
R/bartcore.R, man/*.Rd, inst/tinytest. NEW EXPORTED Rd TOPICS (if any beyond
existing-topic edits) need _pkgdown.yml entries + check_pkgdown. Gate: full
tinytest (grows; a public-surface fit reproduces the internal fixture's
probabilities identically) + R CMD check man + both bitwise anchors + the
multinomial fixture identical. Size: L. Abort: any public-path fit diverges
from the internal path on the same seed.

### C7 landing (2026-07-15) - ARC CLOSED

Landed as bb29d00. bart2(family = "multinomial") over the matrix interface:
factor y, K from levels(y), levels on every K-shaped output, fit class
bartMultinomial (deliberately not "bart" - no single-forest generic can
misread the K-widened array), fitted ev/class, extract ev/ppd, and the
reproduction gate pinning the public fit bit-for-bit to the internal
pattern at K = 2 and K = 3. Gates: tinytest 2900 (2865 + 35 new, zero
regen), both anchors and the multinomial fixture identical, air clean,
codoc clean; no new Rd topic, pkgdown untouched.

Scope narrowings vs this plan's optimistic phrasing, all because
Q3-primary records only the probability channels (each refused BY NAME in
the surface, recorded as follow-ups in the TODO's multi-forest-models
entry): the latent escape hatch (raw per-category fits are
non-identified and unrecorded), the count-matrix likelihood, and the
formula interface (matrix interface only this arc). The predict/test
surface named here (a combined-test-fits location channel) landed as
its own arc, docs/plans/archive/multinomial-test-surface.md (C1 = bcefa63, C2 =
88ffe12); per-sample per-category variable counts (the channel then
addressed category 0 only, so the fit omitted varcount rather than
mislabel it) landed as its own arc too, docs/plans/archive/multinomial-varcounts.md
(C1 = 5afb09a).

## Verification

- Phase A (C1-C3), every commit: equivalence 22/22 IDENTICAL vs
  equivalence-ac6ec2c.rds AND the bcf-equivalence fixture IDENTICAL vs
  bcf-equivalence-99205ee.rds (both anchors - the single-forest and BCF paths are
  the neutrality claim) + full tinytest 2832 no regen + tests/cpp (delete stale
  binaries after header edits) + test-capi.R for C2 (the C_interface.cpp touch). A
  moved draw on either anchor means the L=1 / null-combiner path was not preserved.
  NOTE the anchor-coverage caveat (C2): tinytest, not the anchors, guards the
  train/test output shapes.
- Phase B, every commit: the two bitwise anchors STILL identical (the new model
  must not perturb any existing path) + full tinytest. C4 is one gated commit; it
  lands only with multinomial-exact.R passing to MC error (all three arms - the
  primary correctness gate) + the new tests/cpp cases (incl. the mandatory
  growForestFromRoot case) + the new bitwise multinomial fixture identical to
  itself.
- Phase-A-close bench-sampler compare vs bench-sampler-4008675.csv is
  UNCONDITIONAL (G8): the location-count seam (C2) touches storeSample and the run
  allocation, so the single-forest hot path is CONFIRMED clean by measurement, not
  asserted. dbarts.h unchanged (internal model) -> no stan4bart lockstep; the
  bridge's new .Call entry points earn a "rchk on next scheduled run" note (README
  review step 2).
