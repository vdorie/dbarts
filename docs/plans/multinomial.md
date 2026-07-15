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
window: none for the engine (internal reorg + a new internal model). IF a
  public R surface ships this arc (open question Q5), family lives on
  dbartsModel and the slot is pre-release (family-on-model.md); the
  recommendation is internal-only like BCF, so no window is taken.
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

## Model choice (the settled literature; defer to it per VD's standing guidance)

- AUGMENTATION: log-linear / softmax with Polya-Gamma, one forest per category.
  Given the K forest fits f_ik, P(y_i = k) = softmax(f_i)_k. Category k's
  ONE-VS-REST conditional is a binomial logistic with linear predictor
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
  a PRECISION drawn against stale margins. RECOMMEND this over multinomial probit
  (Kindo et al. 2016): probit needs truncated-MVN latent utilities and a
  covariance sampler, a whole new latent-draw + move machinery, and reuses none of
  the PG code; the log-linear route reuses ext_rng_simulatePolyaGamma and the
  weighted-conjugate kernels untouched. ATTRIBUTION: the algorithm is the
  Held-Holmes / PSW interleaved one-vs-rest conditional cycle; Murray 2021
  ("Log-Linear Bayesian Additive Regression Trees", JASA) is RELATED WORK, not the
  source - its construction uses a per-observation normalizer latent, a different
  device. Q1 records the fork.
- IDENTIFICATION: K SYMMETRIC forests (softmax is invariant to a per-observation
  constant shift of all f_ik, so the raw f_ik are non-identified - exactly BCF's
  a-sign situation). RECOMMEND K over K-1-with-reference: K gives every category
  its own forest and its own variable-count channel (symmetric reporting), and the
  non-identification is handled the way BCF's is - the gate and the reporting
  compare IDENTIFIED functionals (softmax probabilities), never the raw f_ik. K-1
  removes the invariance but makes the reference category's reporting implicit and
  the varcount asymmetric. THE FLAT DIRECTION NEEDS AN EXPLICIT MOVE: the
  likelihood is invariant to a per-observation common shift of all f_ik, a flat
  additive direction pinned ONLY by the prior, which mixes as a slow random walk.
  ADD a likelihood-invariant LEVEL-CENTERING Gibbs move (the BCF-ridge-interweave
  analog) in afterCombine: draw the per-observation common shift from its Gaussian
  full conditional under the forest priors and recenter f_ik. Q2 records the fork.
- SINGLE-TRIAL this arc (n_i = 1 classification). The one-vs-rest binomial for
  category k needs the SUCCESS COUNT y_ik, not just a category label; borrowed
  integer labels support only n_i = 1. Scope this arc to single-trial
  classification and record the grouped-count generalization (store the n x K
  count matrix, n_i > 1) as a follow-up - PG(n_i, .) is the sum of n_i PG(1, .)
  draws (weighted-logistic.md), so grouped counts add no numerical code, only the
  count-matrix data model, and the exact gate below is single-trial anyway.
  Real-shape (non-integer) counts stay out of scope (the deferral
  weighted-binary/negative-binomial carry).

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
  before that category's forest update (M1 above), never batched post-loop.
- Creation mirrors BCF: R/bartcore.R:482 bartcoreBCFSampler packs a params vector
  + moderators, calls C_dbarts_bartcore_createBCF (R_interface_bartcore.cpp:1879),
  which builds a BCFSpec (chain.hpp:348) and calls createBCFSampler
  (facade.hpp:499, single ConstantGaussianLeaf instantiation). The multinomial
  analog: bartcoreMultinomialSampler + C_dbarts_bartcore_createMultinomial +
  MultinomialSpec + createMultinomialSampler.
- ResponseModel is per-observation-single-location (model.hpp:1788); multinomial
  has no sigma (fixed, like the binary families) and K latent channels, so its
  response family is thin (Q3): a MultinomialResponse whose single-location seams
  are vestigial (no-op refreshLatents/drawSigma) and whose log-likelihood channel
  is NaN-flagged - storeSample scores one forest's fits via
  computeLogLikelihood (chain.hpp:2586-2593) and cannot see the K-blend, so
  logLikelihoodIsDefined() = false is BCF's exact choice (chain.hpp:638). The
  combiner owns the K interleaved PG draws, the per-forest working responses, and
  the level-centering move.
- State: the by-name block reader (R_interface_bartcore.cpp:3317
  stateFormatVersion = 3, minReadable = 3) already serializes a per-forest tree
  list of any length (SLOT_FORESTS), so K forests serialize with no format work;
  only combiner-specific scalar/latent state would need a block (Q4).

## Key finding (contradicts the design note's re-carve framing; grounds Q3)

docs/design/forest-combiner.md lumps three seams as forced to widen for a non-BCF
combiner: combinedFits -> const double*, refreshLatents/drawSigma take one
location, results.trainingFits is one channel. Read in code, only the
REPORTING/combined-OUTPUT seam genuinely widens, and only if we want the engine
to report probabilities:

- refreshLatents/drawSigma do NOT widen. The K-category PG augmentation lives in
  the combiner's INTERLEAVED per-forest hook (M1) - each omega_k drawn against the
  current margins just before forest k's update - and multinomial has no sigma. So
  MultinomialResponse's refreshLatents/drawSigma are no-ops and the single-location
  `combined` they are handed is ignored - no signature change, and every existing
  family is untouched. GUARD (M6): the K x n combinedFits buffer is handed to the
  no-op refreshLatents as a single-location pointer; carry a guard note (assert /
  comment) so a future non-no-op MultinomialResponse::refreshLatents cannot
  silently misread a K-channel buffer as one channel.
- combinedFits + trainingFits/testFits widen to K ONLY under the
  engine-reports-probabilities choice. The minimal alternative (Q3) keeps them at
  one location, NaN-flags the scalar channels exactly as BCF does, exposes the K
  forests via getForestFits (already K-capable), and lets the R surface compute
  softmax probabilities. That path needs NO combined-output widening at all - the
  n x K softmax object stays internal to MultinomialForestCombiner (it reads
  forests_ directly for the margins).

The plan's PRIMARY line takes the widening (engine reports identified
probabilities - see Q3 rationale: for a classification model probabilities are
the deliverable and the raw per-category fits are non-identified nuisance). Q3
carries the minimal alternative with its full tradeoff.

## Constraints

- dbarts.h FROZEN. Multinomial is internal-only (bartcore .Call), like BCF; the
  bartcore::Results widening (Phase A) is the internal struct, NOT the public
  dbarts_results the LinkingTo consumers (stan4bart) call - no header change, no
  lockstep, no rchk.
- Single-forest AND BCF hot paths pay nothing new. The location-count is 1 for
  every existing model; combiner_ stays null off multi-forest; every Phase-A
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
   tests/cpp's Makefile HAS -MMD/-MP dependency tracking (tests/cpp/Makefile:8;
   the CLAUDE.md "no dependency tracking" note is outdated for tests/cpp) - binary
   deletion is belt-and-suspenders. Files: src/bartcore/combiner.hpp (new),
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
   K-channel combined scratch. THE INTERLEAVED PG DRAW (M1): omega_k is drawn
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
   (M3) lives in afterCombine: draw the per-observation common shift from its
   Gaussian full conditional under the forest priors and recenter f_ik (the
   BCF-ridge-interweave analog), pinning the flat additive direction the softmax
   likelihood is invariant to. combinedFits writes the K softmax probabilities per
   observation (numReportedLocations() = K); reportedForest addresses per category
   via the C3 query.

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
   FACTOR over probit's 3.0, dbarts.R:394-395; negative-binomial.md:39 confirms
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

C5. State round-trip (separable under Q4-primary). The K forests already serialize
   through SLOT_FORESTS (the per-forest list, any length). Under Q4-primary (redraw
   the PG latents on restore, matching the "restore is structural" contract
   test-bcf.R/state-format-policy.md), NO wire-format bump is needed: serializeGlue/
   restoreGlue stay no-ops for the multinomial combiner (it carries no
   un-recoverable scalar glue - softmax has no a/b coefficients), omega is
   cold-restarted and redrawn on the first restored sweep, and re-creation is
   driven by the stored R-side creation call as BCF's is. THE INTERLEAVED DRAW
   (M1) also resolves the restore wrinkle: the first restored sweep seeds omega_k
   against the RESTORED forests correctly, because omega_k is drawn just before
   forest k's update against whatever margins the restored forests present. If Q4
   resolves to serialize the n x K omega, add ONE append-only by-name block, bump
   stateFormatVersion to 4, hold minReadableStateFormatVersion at 3 (additive, per
   c-api-growth.md's registry rule). Files: chain.hpp (serializeGlue/restoreGlue if
   Q4-serialize), R_interface_bartcore.cpp (a new SLOT only if Q4-serialize). Gate:
   a state round-trip component test (structural, NOT bitwise - the BCF fixture's
   rule) + both anchors + tinytest. Size: S (Q4-primary) or M (Q4-serialize).

C6. Docs. Design note docs/design/multinomial.md (likelihood; the INTERLEAVED
   one-vs-rest PG cycle with the per-forest hook; the LEVEL-CENTERING move;
   identification; the prior/leaf-scale calibration curve for K = 2, 3, 5; the
   exact-gate construction incl. the correlated-difference prior tau^2(I + 11');
   the surface; the Q3/Q4 resolutions; ATTRIBUTION - Held-Holmes 2006 / PSW 2013
   sec 4 as the algorithm, Murray 2021 as related work, not the source). Update
   multi-forest-models.md (multinomial landed) and docs/design/forest-combiner.md's
   "Anticipated, not built" (built) + the combiner.hpp extraction record + the
   correction that the softmax combiner needs no wire state. Files:
   docs/design/multinomial.md, docs/plans/*. Gate: full tinytest; R CMD check man
   unaffected (no man/ touched). Size: M.

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

## Open questions for VD

- Q1 (augmentation: interleaved PG softmax vs multinomial probit). PG log-linear
  softmax reuses the shipped PG(1, psi) sampler and the weighted-conjugate kernels
  verbatim - one binomial-logistic per category, drawn INTERLEAVED (omega_k against
  the current margins just before forest k's update, via the per-forest hook) -
  and costs no new numerical code beyond that hook; its price is the softmax
  non-identification (handled like BCF's a-sign, comparing probabilities) plus the
  interleaving discipline. Multinomial probit (Kindo 2016) is identified via a
  reference category but needs truncated multivariate-normal latent utilities and a
  covariance sampler - a whole new latent + move machinery reusing none of the PG
  code, and a heavier exact gate. The algorithm is Held-Holmes 2006 / PSW 2013
  (Murray 2021 is related work). RECOMMEND PG softmax: it reuses everything and the
  non-identification is a solved problem here.
- Q2 (K vs K-1 forests). K symmetric forests give every category its own forest
  and variable-count channel, at the cost of a per-observation additive
  non-identification of the raw f_ik (the probabilities are identified) that
  requires the explicit level-centering Gibbs move in afterCombine (M3). K-1 with a
  reference category removes the invariance (and the centering move) but makes the
  reference category's effect implicit, its reporting asymmetric, and the
  combiner's per-category symmetry uneven. RECOMMEND K symmetric: the reporting
  symmetry is worth the one added centering move the gate and R surface already
  sidestep by reporting probabilities.
- Q3 (where the K-location generalization stops). PRIMARY: widen the combined-
  OUTPUT seam to a location count (C2) so the engine reports the K softmax
  PROBABILITIES directly in trainingFits/testFits - identified, and the natural
  deliverable for a classification model, where the raw per-category fits are
  non-identified nuisance a user should not have to softmax by hand. MINIMAL: keep
  every seam at one location, NaN-flag the scalar train/test channels exactly as
  BCF does, expose the K forests via getForestFits (already K-capable), and compute
  softmax R-side - the n x K object stays internal to the combiner and Phase A
  shrinks to C1 and C3. RECOMMEND PRIMARY: multinomial's deliverable is
  probabilities, unlike BCF where the per-forest fits ARE the estimands; the
  widening is the reusable multi-location channel heteroscedastic/hurdle will also
  want, and it stays fully bitwise on both anchors because L = 1 everywhere today.
- Q4 (state wire format). PRIMARY: NO bump. The K forests serialize through the
  existing per-forest SLOT_FORESTS list, the softmax combiner carries no
  un-recoverable scalar glue (no a/b), and the PG omega latents are redrawn on
  restore - the interleaved draw (M1) seeds them against the restored forests
  correctly on the first restored sweep - a structural (not bitwise) continuation,
  already the contract (test-bcf.R, state-format-policy.md). This CONTRADICTS the
  design note's blanket "a non-BCF combiner needs a format bump": softmax has
  nothing of its own shape to write. ALTERNATIVE: serialize the n x K omega in one
  append-only by-name block (stateFormatVersion 3 -> 4, minReadable held at 3) for
  a tighter restore. RECOMMEND PRIMARY: no bump, redraw omega; it is the cheapest,
  matches the structural-restore contract, and the finding that the softmax
  combiner needs no wire state should be recorded in the design note as a
  correction.
- Q5 (public R surface this arc, or internal like BCF). BCF ships internal-only
  (bartcoreBCFSampler, an unexported .Call path; no family = "bcf"). Multinomial
  could do the same - internal creation, exercised by the exact gate and by future
  consumers via the bartcore surface - or grow a public bart2(family =
  "multinomial") / a multinomial front-end this arc. Internal-only keeps the arc
  focused on the engine + the coupling correctness and takes no pre-release window;
  a public surface adds argument design, factor-response ingestion, a
  prediction/probability front-end, and Rd + tinytest coverage, and (per
  family-on-model.md) the family slot lives on dbartsModel. RECOMMEND
  INTERNAL-ONLY this arc, matching BCF: land the engine + the gate now, file the
  public surface as its own follow-up once the coupling is proven - the same
  staging forest-split-bcf/BCF used.
