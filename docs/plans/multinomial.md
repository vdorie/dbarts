# multinomial

agent: opus (the first multi-forest MODEL. The softmax coupling draw (a
  Polya-Gamma augmentation per category), the K-forest combiner, the response
  family, and the exact-posterior gate are all posterior-defining; the Phase-A
  seams are RNG-neutral relocations but touch the sweep, storeSample, and the
  run bridge, where a leak is a silent posterior change).
rng: TWO CLASSES, cleanly gated. Phase A (steps 1-3): NEUTRAL - the combiner
  extraction and the combined-output/reporting widening to K locations, every
  step bitwise on BOTH anchors (single-forest equivalence 22/22 IDENTICAL vs
  equivalence-ac6ec2c.rds + the BCF fixture's five scenarios x six channels
  IDENTICAL vs bcf-equivalence-99205ee.rds). Phase B (steps 4-8): POSTERIOR-
  CHANGING (new model) - gated by a NEW multinomial exact-posterior benchmark
  and component tests, with every pre-existing path still bitwise on both
  anchors.
window: none for the engine (internal reorg + a new internal model). IF a
  public R surface ships this arc (open question Q5), family lives on
  dbartsModel and the slot is pre-release (family-on-model.md); the
  recommendation is internal-only like BCF, so no window is taken.
budget: ~1200-1600 lines. Phase A ~400 (combiner.hpp motion + the location-
  count seam + the run-bridge allocation). Phase B ~800-1200
  (MultinomialForestCombiner, the response family, the creation surface,
  reporting, the exact gate + component tests). Chiefly src/bartcore/
  combiner.hpp (new), chain.hpp, model.hpp, R_interface_bartcore.cpp,
  R/bartcore.R, benchmarks/R/multinomial-exact.R + baseline,
  tests/cpp/test_sampler.cpp. Header edits -> --preclean; delete stale
  tests/cpp binaries.

## Goal

An internal K-category multinomial sampler: K forests coupled through a softmax
(log-linear) likelihood with Polya-Gamma augmentation, riding the forest
combiner the BCF refactor introduced. The combiner is BCF's second instance and
the occasion to extract Forest<L>+ForestCombiner<L> to combiner.hpp (the
deferral recorded in docs/design/forest-combiner.md). Chain gains a
location-count on its combined-output/reporting seams (1 for every existing
model, K for multinomial) so the widening is bitwise-neutral on the
single-forest and BCF paths. Creation mirrors BCF's internal .Call surface;
dbarts.h stays frozen.

## Model choice (the settled literature; defer to it per VD's standing guidance)

- AUGMENTATION: log-linear / softmax with Polya-Gamma, one forest per category
  (Murray 2021, "Log-Linear Bayesian Additive Regression Trees", JASA). Given
  the K forest fits f_ik, P(y_i = k) = softmax(f_i)_k. Conditional on the other
  categories, category k's likelihood is a binomial logistic with linear
  predictor eta_ik = f_ik - C_ik, where C_ik = log sum_{j != k} exp(f_ij) is
  the log-sum-exp margin. So omega_ik ~ PG(n_i, eta_ik) and forest k sees
  working response (y_ik - n_i/2)/omega_ik + C_ik under precision omega_ik -
  the EXACT logistic PG machinery already shipped, one binomial per category
  per sweep, cycling the categories. RECOMMEND this over multinomial probit
  (Kindo et al. 2016): probit needs truncated-MVN latent utilities and a
  covariance sampler, a whole new latent-draw + move machinery, and reuses none
  of the PG code; the log-linear route reuses ext_rng_simulatePolyaGamma and the
  weighted-conjugate kernels untouched. Q1 records the fork.
- IDENTIFICATION: K SYMMETRIC forests (softmax is invariant to a per-observation
  constant shift of all f_ik, so the raw f_ik are non-identified - exactly BCF's
  a-sign situation). RECOMMEND K over K-1-with-reference: K gives every category
  its own forest and its own variable-count channel (symmetric reporting), matches
  Murray's formulation, and the non-identification is handled the way BCF's is -
  the gate and the reporting compare IDENTIFIED functionals (softmax
  probabilities), never the raw f_ik. K-1 removes the invariance but makes the
  reference category's reporting implicit and the varcount asymmetric. Q2 records
  the fork.
- Integer counts only (n_i trials): PG(n_i, .) is the sum of n_i PG(1, .) draws
  (weighted-logistic.md), so single-trial classification (n_i = 1) and grouped
  counts both ride the shipped sampler with no new numerical code. Real-shape
  counts are out of scope (the same deferral weighted-binary/negative-binomial
  carry).

## Context (seams, read in code)

- The combiner landed in chain.hpp beside Forest<L>: ForestCombiner<L> base
  (chain.hpp:394-439) and BCFForestCombiner<L> (447-...). The virtuals a second
  combiner needs already exist: formForestResponse (per-forest (response,
  weights) over the whole forest vector, 400-405), combinedFits (405), drawGlue
  (421) + afterCombine (423), the reporting map reportedForest/testFitsAreDefined
  /logLikelihoodIsDefined (431-433), serializeGlue/restoreGlue (437-438),
  setTreatment/bcfGlue (410-415). Extraction to combiner.hpp is the deferred
  decision (docs/design/forest-combiner.md "Header location").
- The single n-vector combined output, the seam the design note flags
  (docs/design/forest-combiner.md "What still re-carves"): Chain::combinedFits()
  returns const double* (chain.hpp:2498); storeSample writes ONE trainingFits and
  ONE testFits channel (2518-2554); the run bridge allocates trainFits as n x
  numSamples (R_interface_bartcore.cpp:1987) and testFits as nTest x numSamples
  (1996). refreshLatents/drawSigma take one location (chain.hpp:1054-1065) - but see the
  finding below: multinomial does NOT widen those.
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
  with the log-sum-exp margin C_ik as the per-category "offset".
- Creation mirrors BCF: R/bartcore.R:482 bartcoreBCFSampler packs a params vector
  + moderators, calls C_dbarts_bartcore_createBCF (R_interface_bartcore.cpp:1879),
  which builds a BCFSpec (chain.hpp:348) and calls createBCFSampler
  (facade.hpp:499, single ConstantGaussianLeaf instantiation). The multinomial
  analog: bartcoreMultinomialSampler + C_dbarts_bartcore_createMultinomial +
  MultinomialSpec + createMultinomialSampler.
- ResponseModel is per-observation-single-location (model.hpp:1788); multinomial
  has no sigma (fixed, like the binary families) and K latent channels, so its
  response family is thin (Q3): a MultinomialResponse carrying the category
  labels for the log-likelihood channel and cold-starting nothing on the
  single-location seams, while the combiner owns the K-category PG draws and the
  per-forest working responses.
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

- refreshLatents/drawSigma do NOT widen. The recommended architecture puts the
  K-category PG augmentation in the combiner's drawGlue (it already receives the
  full forest vector), and multinomial has no sigma. So MultinomialResponse's
  refreshLatents/drawSigma are no-ops and the single-location `combined` they are
  handed is ignored - no signature change, and every existing family is
  untouched.
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
  touchpoint keeps its `if (combiner_)` shape. Equivalence 22/22 IDENTICAL and
  the BCF fixture IDENTICAL are the per-step proof, both anchors, every commit.
- Per-SWEEP virtual dispatch only (core-generalization.md). The combiner's PG
  loop and margin loop are O(n*K) inside one drawGlue call per sweep, never in
  the per-observation partition/suffstat kernels.
- fast over safe in C/C++; header-only C++20; Doxygen/LLVM/ASCII; the
  comment/code delta bounded (README brevity).
- OUT of scope: heteroscedastic and hurdle (their own notes); real-shape counts;
  a per-observation sigma; ordinal outcomes (ordinal-outcomes.md, cumulative
  probit - a different model); grouped-x-multinomial composition (the combiner
  makes it expressible via GroupedResponse below the combiner, but this arc does
  not build it).

## Phase A - RNG-neutral seams (steps 1-3)

Every Phase-A step is gated: equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c
+ the BCF fixture five-scenarios-x-six-channels IDENTICAL vs
bcf-equivalence-99205ee + full tinytest 2832 no regen + tests/cpp (delete stale
binaries after header edits). If any step cannot be made bitwise on both
anchors, the seam design is wrong - stop and rethink.

1. Extract combiner.hpp (pure motion; NEUTRAL). Move Forest<L>, ForestResponse,
   ForestCombiner<L>, and BCFForestCombiner<L> from chain.hpp into a new
   src/bartcore/combiner.hpp; chain.hpp includes it. This is the second-consumer
   extraction docs/design/forest-combiner.md deferred - do it FIRST so steps 4+
   add MultinomialForestCombiner beside BCFForestCombiner in the new home. No
   behavior change; the include order (combiner.hpp needs Forest<L>, so Forest<L>
   moves with it) is the whole content of the step. Files: src/bartcore/
   combiner.hpp (new), chain.hpp. Gate: NEUTRAL - both anchors identical +
   tinytest + tests/cpp. Size: M (motion is mechanical but wide). --preclean;
   delete tests/cpp binaries.

2. Location-count on the combined-output/reporting seam (NEUTRAL; 1 everywhere
   today). Give the combiner a numReportedLocations() (default 1; BCF 1) and
   widen combinedFits to write numReportedLocations x n into caller-owned scratch
   (BCF keeps writing one channel). storeSample's trainingFits/testFits loop over
   the reported locations (one iteration off multinomial). bartcore::Results
   grows a location stride the run bridge reads (Q3-primary); the bridge allocates
   trainFits as n x L x numSamples and testFits as nTest x L x numSamples with
   L = numReportedLocations (L = 1 -> the exact current n x numSamples shape and
   byte layout). The facade forwards numReportedLocations. Because L = 1 for every
   existing model the allocation, the storeSample writes, and the draw values are
   byte-identical - equivalence (which runs through bartcore_run) is the proof.
   Files: combiner.hpp, chain.hpp (combinedFits, storeSample, Results),
   facade.hpp, sampler.hpp, R_interface_bartcore.cpp (bartcore_run allocation).
   Gate: NEUTRAL - both anchors identical (the L=1 layout is unchanged) +
   tinytest + tests/cpp. Size: L. Abort: any equivalence or BCF-fixture channel
   moves = the L=1 path is not byte-identical.

3. Per-forest reporting query for K forests (NEUTRAL; BCF-inert). Add a per-forest
   variable-count query beside getForestFits (getForestVariableCounts(forest) ->
   forestVariableCounts(chain, forestIndex, out)), so a K-forest model can report
   each category's split usage without widening the single storeSample varcount
   channel (which keeps addressing reportedForest()). BCF exposes forest 0 as
   today; the new query is additive and unused off multinomial. Files:
   facade.hpp, sampler.hpp, chain.hpp, R_interface_bartcore.cpp, R/bartcore.R.
   Gate: NEUTRAL - both anchors identical + tinytest + tests/cpp. Size: S.
   (Foldable into step 6 if VD prefers; kept separate so the neutral query lands
   under the bitwise gate before any model reads it.)

## Phase B - the multinomial model (steps 4-8)

POSTERIOR-CHANGING. Gated by a NEW exact-posterior benchmark
(benchmarks/R/multinomial-exact.R) + new component tests, WITH equivalence 22/22
+ the BCF fixture still IDENTICAL at every step (the new model must not perturb
any existing path) + full tinytest.

4. MultinomialForestCombiner<L> + the response family (engine). Add
   MultinomialForestCombiner<L> to combiner.hpp (constant-leaf, static_assert as
   BCF): holds the borrowed integer category labels and trial counts, the n x K
   PG omega scratch (cold-started at n_i/4 per category, the logistic seed), and
   the K-channel combined scratch. drawGlue draws omega_ik ~ PG(n_i, eta_ik) via
   summed ext_rng_simulatePolyaGamma over the log-sum-exp margins C_ik computed
   from forests_ (its own internal draw order over the K categories, fixed and
   documented - the Q3-resolved "combiner owns its draw order"). formForestResponse
   (k, forests, y, w) forms category k's PG working response (y_ik - n_i/2)/
   omega_ik + C_ik under weight omega_ik, IGNORING the passed y (the combiner owns
   the labels). combinedFits writes the K softmax probabilities per observation
   (numReportedLocations() = K). reportedForest addresses per category via the
   step-3 query; testFitsAreDefined/logLikelihoodIsDefined per Q3. A thin
   MultinomialResponse (model.hpp): no sigma (sigmaIsFixed), no-op
   refreshLatents/drawSigma, carries the labels for computeLogLikelihood
   (categorical log P(y_i)); its single-location seams are vestigial. Files:
   combiner.hpp, model.hpp, chain.hpp (the multinomial Chain constructor building
   K forests + the combiner, mirroring the BCF constructor). Gate: exact benchmark
   (step 8) + component tests (step 7) + both bitwise anchors + tinytest. Size: L.

5. Creation + reporting surface (bridge + R). MultinomialSpec (chain.hpp: K, the
   per-category forest spec, the leaf-scale calibration on the PG latent scale -
   pi/sqrt(3) per category, following logistic's node.scale precedent and
   negative-binomial.md's latent-scale paragraph), createMultinomialSampler
   (facade.hpp, single ConstantGaussianLeaf instantiation as BCF),
   createMultinomialHolder + C_dbarts_bartcore_createMultinomial
   (R_interface_bartcore.cpp: K forests from the model spec, labels validated as
   0..K-1 integers, refuse a non-multinomial response), and
   bartcoreMultinomialSampler + getForestFits/getForestVariableCounts wrappers
   (R/bartcore.R). The run output reports the K softmax-probability channels
   (Q3-primary) or NaN-flags them and relies on getForestFits (Q3-minimal); refuse
   the single-forest test surface as refuseBCFTestSurface does. Files:
   chain.hpp, facade.hpp, sampler.hpp, R_interface_bartcore.cpp, R/bartcore.R.
   Gate: as step 4 + rchk note (the bridge is touched - new .Call entry points;
   "rchk on next scheduled run" per README). Size: L.

6. State (de)serialization. The K forests already serialize through SLOT_FORESTS
   (the per-forest list, any length). Under Q4-primary (redraw the PG latents on
   restore, matching the "restore is structural" contract test-bcf.R/
   state-format-policy.md), NO wire-format bump is needed: serializeGlue/
   restoreGlue stay no-ops for the multinomial combiner (it carries no
   un-recoverable scalar glue - softmax has no a/b coefficients), omega is
   cold-restarted and redrawn on the first restored sweep, and re-creation is
   driven by the stored R-side creation call as BCF's is. If Q4 resolves to
   serialize the n x K omega, add ONE append-only by-name block, bump
   stateFormatVersion to 4, hold minReadableStateFormatVersion at 3 (additive,
   per c-api-growth.md's registry rule). Files: chain.hpp (serializeGlue/
   restoreGlue if Q4-serialize), R_interface_bartcore.cpp (a new SLOT only if
   Q4-serialize). Gate: a state round-trip component test (structural, NOT
   bitwise - the BCF fixture's rule) + both anchors + tinytest. Size: S
   (Q4-primary) or M (Q4-serialize).

7. Component tests (tests/cpp). testMultinomial: a deterministic fixed-seed fit
   recovering a monotone K=3 category signal; the PG margin/omega draw checked
   for positivity and the K=2 reduction (a two-category multinomial fit must
   reproduce the logistic path's working response/weights construction); the
   combiner's growForestFromRoot branch if reachable; a state round-trip.
   PG-moment tests at the trial counts used (mean/variance of PG(n, psi),
   following weighted-logistic's component test). Files: tests/cpp/
   test_sampler.cpp. Gate: tests/cpp passes; both anchors still identical. Size:
   M.

8. Exact-posterior gate (benchmarks/R/multinomial-exact.R) + design note + docs.
   The gate, in bcf-exact.R / categorical-exact.R style: (a) INTERCEPT-ONLY K = 3
   softmax (each forest a single leaf via a constant predictor), y ~
   Multinomial(softmax(f)), f_k ~ N(0, tau^2) - the posterior mean of the
   IDENTIFIED softmax probabilities computed by (K-1)-dimensional grid quadrature
   over the counts, matched to the sampler's mean category probabilities to MC
   error (raw f_k are sign/shift-ill-posed by construction, compared never -
   BCF's a precedent); (b) a K = 2 arm asserting the multinomial fit matches the
   shipped logistic sampler (the strongest internal-consistency check - K = 2
   multinomial IS binary logistic); optionally (c) a small-cell (one categorical
   predictor, 2 cells) extension if the per-cell quadrature stays tractable.
   Failure means the softmax coupling, the PG draw, or the margin formation is
   wrong. Record the baseline; add it to MANIFEST/equivalence.yaml if a bitwise
   multinomial fixture is also wanted (Q3-primary makes the softmax-probability
   channel a natural bitwise fixture like bcf-equivalence). Design note
   docs/design/multinomial.md (likelihood, the K-forest softmax coupling, the PG
   augmentation, identification, prior/leaf-scale calibration, the surface, the
   Q3/Q4 resolutions); update multi-forest-models.md (multinomial landed) and
   docs/design/forest-combiner.md's "Anticipated, not built" (built) + the
   combiner.hpp extraction record. Files: benchmarks/R/multinomial-exact.R,
   benchmarks/baselines/multinomial-exact* or a bitwise fixture baseline,
   docs/design/multinomial.md, docs/plans/*. Gate: the exact benchmark passes to
   MC error; both bitwise anchors identical; full tinytest. Size: M.

## Verification

- Phase A, every commit: equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds
  AND the bcf-equivalence fixture IDENTICAL vs bcf-equivalence-99205ee.rds
  (both anchors - the single-forest and BCF paths are the neutrality claim) +
  full tinytest 2832 no regen + tests/cpp (delete stale binaries after header
  edits). A moved draw on either anchor means the L=1 / null-combiner path was
  not preserved.
- Phase B, every commit: the two bitwise anchors STILL identical (the new model
  must not perturb any existing path) + full tinytest. On landing:
  multinomial-exact.R passes to MC error (the primary correctness gate) + the new
  tests/cpp cases + the K=2==logistic arm.
- No bench-sampler owed unless a Phase-A step measurably touches the single-forest
  hot path (the location-count is compile-time-1 and the null-combiner
  short-circuit is unchanged, so none should); one confirmatory bench-sampler
  compare vs bench-sampler-4008675.csv at Phase-A close if step 2's storeSample
  edit is non-trivial. dbarts.h unchanged (internal model) -> no stan4bart
  lockstep; the bridge's new .Call entry points earn a "rchk on next scheduled
  run" note (README review step 2).

## Open questions for VD

- Q1 (augmentation: PG softmax vs multinomial probit). PG log-linear softmax
  (Murray 2021) reuses the shipped PG(1, psi) sampler and the weighted-conjugate
  kernels verbatim - one binomial-logistic per category per sweep - and costs no
  new numerical code; its price is the softmax non-identification (handled like
  BCF's a-sign, comparing probabilities). Multinomial probit (Kindo 2016) is
  identified via a reference category but needs truncated multivariate-normal
  latent utilities and a covariance sampler - a whole new latent + move machinery
  reusing none of the PG code, and a heavier exact gate. RECOMMEND PG softmax: it
  is the settled PG-based BART multinomial, reuses everything, and the
  non-identification is a solved problem here.
- Q2 (K vs K-1 forests). K symmetric forests give every category its own forest
  and variable-count channel and match Murray's log-linear form, at the cost of a
  per-observation additive non-identification of the raw f_ik (the probabilities
  are identified). K-1 with a reference category removes the invariance but makes
  the reference category's effect implicit, its reporting asymmetric, and the
  combiner's per-category symmetry uneven. RECOMMEND K symmetric: the reporting
  symmetry is worth more than removing an invariance the gate and the R surface
  already sidestep by reporting probabilities.
- Q3 (where the K-location generalization stops). PRIMARY: widen the combined-
  OUTPUT seam to a location count (Phase A step 2) so the engine reports the K
  softmax PROBABILITIES directly in trainingFits/testFits - identified, and the
  natural deliverable for a classification model, where the raw per-category fits
  are non-identified nuisance a user should not have to softmax by hand.
  MINIMAL: keep every seam at one location, NaN-flag the scalar train/test
  channels exactly as BCF does, expose the K forests via getForestFits (already
  K-capable), and compute softmax R-side - the n x K object stays internal to the
  combiner and Phase A shrinks to steps 1 and 3. The minimal path is smaller and
  strictly BCF-shaped; the primary path costs step 2 (a bitwise-neutral L=1
  widening of storeSample + the run allocation) and buys engine-computed
  identified probabilities. RECOMMEND PRIMARY: multinomial's deliverable is
  probabilities, unlike BCF where the per-forest fits ARE the estimands; the
  widening is the reusable multi-location channel heteroscedastic/hurdle will also
  want, and it stays fully bitwise on both anchors because L = 1 everywhere today.
- Q4 (state wire format). PRIMARY: NO bump. The K forests serialize through the
  existing per-forest SLOT_FORESTS list, the softmax combiner carries no
  un-recoverable scalar glue (no a/b), and the PG omega latents are redrawn on
  restore - a structural (not bitwise) continuation, which is already the
  contract (test-bcf.R, state-format-policy.md; the BCF fixture explicitly does
  NOT bitwise-compare a getState->setState round-trip). This CONTRADICTS the
  design note's blanket "a non-BCF combiner needs a format bump": softmax has
  nothing of its own shape to write. ALTERNATIVE: serialize the n x K omega in one
  append-only by-name block (stateFormatVersion 3 -> 4, minReadable held at 3 per
  c-api-growth's additive rule) for a tighter restore. RECOMMEND PRIMARY: no bump,
  redraw omega; it is the cheapest, matches the structural-restore contract, and
  the finding that the softmax combiner needs no wire state should be recorded in
  the design note as a correction.
- Q5 (public R surface this arc, or internal like BCF). BCF ships internal-only
  (bartcoreBCFSampler, an unexported .Call path; no family = "bcf"). Multinomial
  could do the same - internal creation, exercised by the exact gate and by
  future consumers via the bartcore surface - or grow a public bart2(family =
  "multinomial") / a multinomial front-end this arc. Internal-only keeps the arc
  focused on the engine + the coupling correctness and takes no pre-release
  window; a public surface adds argument design, factor-response ingestion, a
  prediction/probability front-end, and Rd + tinytest coverage, and (per
  family-on-model.md) the family slot lives on dbartsModel. RECOMMEND
  INTERNAL-ONLY this arc, matching BCF: land the engine + the gate now, file the
  public surface as its own follow-up once the coupling is proven - the same
  staging forest-split-bcf/BCF used.
