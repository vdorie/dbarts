# GP leaves and the non-conjugate move strategy (phase 6)

Proposal, 2026-07-04; since LANDED through stage 4 - see ## Status at
the end. Phase 6 of docs/design/core-generalization.md is
"Non-conjugate MoveStrategy: GP leaves, general likelihoods". This
document separates those two axes, argues they land in that order, and
stages the work; the proposal text is preserved as written.

## The two axes are separable

The conjugate engine samples tree structure with leaf parameters
marginalized out: metropolisJumpForTree (moves.hpp) scores branches via
logIntegratedLikelihoodForNode, and sampleParametersAndSetFits draws
leaf parameters given the accepted structure. That machinery requires a
closed-form marginal - the IntegrableLeafModel concept - not a cheap
one.

A GP leaf under any working-Gaussian response (gaussian, probit,
logistic via Polya-Gamma, grouped decorations) IS integrable: with
f_leaf ~ GP(0, tau^2 K) over the leaf's rows and working response z
with weights w,

    log p(z_leaf) = log N(z_leaf; 0, sigma^2 W^-1 + tau^2 K),

a closed form costing one O(n_leaf^3) Cholesky. So GP leaves ride the
EXISTING move strategy; what they need is a third leaf-parameter shape
(function-valued, variable length per leaf) and cost guardrails.

The non-conjugate MoveStrategy is for models with no working-Gaussian
form at all - general likelihoods where no marginal exists and
structure must be sampled jointly with parameters. No concrete consumer
is scheduled; GP leaves do not need it. Recommendation: land GP leaves
first (part 1), keep the non-conjugate strategy as a designed-for
follow-up (part 2), and let a real consumer pull it in.

## Part 1: GP leaves on the conjugate engine

### Model

node.prior = gp(columns, k, ...) designates q continuous columns, as
linear() does. Per bottom node, the leaf function over the standardized
designated covariates u is

    f ~ GP(0, (scale / k)^2 C_theta),    scale = nodeScale / sqrt(m),

with C_theta a squared-exponential kernel over u (per-column
lengthscales theta, plus a small nugget for conditioning). The prior
variance ties to k exactly as the constant and linear leaves do, so the
forest shrinkage story is unchanged; a constant-kernel limit recovers
the constant leaf.

The leaf "parameter" is the vector of function values at the leaf's
training rows. Two consequences that simplify the engine work:

- Training fits ARE the parameters: like the scalar leaf (whose
  recoverParametersFromFits reads fits), a GP leaf's drawn values live
  directly in treeFits_. No fixed-stride paramsByTree_ blocks apply.
- Test and saved-tree prediction is the GP conditional mean given the
  drawn training values: f(x*) = c(x*)' C^-1 f, which requires the
  leaf's training covariate rows and per-leaf weights alpha = C^-1 f at
  prediction time (see formats below).

### Concept extension

VectorLeafModel assumes a fixed numParams() stride; GP parameters are
per-observation. Add a third shape:

    concept FunctionLeafModel = LeafModelCore<L> && L::hasFunctionParams
      && requires { drawFromPosteriorForNode(rng, tree, y, weights, k,
                    sigmaSq, node, fitsOut);   // n_leaf values
                    predictForNode(...); ... };

Chain code branches if constexpr on the third flag next to the existing
hasVectorParams branches; scalar and vector paths compile untouched
(the linear-leaves precedent: equivalence must report identical draws
for existing models).

Suffstat shape is the roadmap's "universal fallback": the node's
observation-index span (tree.indices + node.begin, node.end). Nodes
already carry it; accumulate is O(1) and all cost lives in the leaf's
own math. Covariate access copies the linear-leaf pattern exactly: the
leaf gathers standardized designated columns at initialize/regather
(rawColumn serves raw matrices, views' gathered copies, and mixed
stores' dense-backed columns; standardization constants inherit through
buildFromParent), so views/xbart and mixed input compose for free and
the facade refuses designations the store cannot serve raw - all
machinery that exists today.

### Cost model and guardrails

Every move scores affected bottoms: a birth/death costs a Cholesky of
each affected child, O(n_leaf^3). BART leaves are usually small, but
nothing bounds them. Guardrails, in order of preference:

1. A hard cap on GP-leaf size (option, default a few hundred): the
   integrated likelihood of an over-cap branch returns the veto value
   (the empty-leaf precedent in logLikelihoodForBranch), so trees
   simply do not grow leaves the exact math cannot afford. Simple,
   honest, and matches how smooth-BART is used (small m, smooth
   local pieces). [SUPERSEDED at landing: the veto deadlocks from the
   root; over-cap leaves fall back to the constant leaf's exact math
   instead - see the stage 1 landing notes.]
2. Cache per-leaf Cholesky factors keyed on the node's span, reusing
   them across the four move types within a sweep (moves re-score the
   same unchanged siblings repeatedly). Optimization only; land after
   correctness.
3. Low-rank/inducing-point kernels: out of scope for v1; the concept
   accommodates them later behind the same interface.

### Hyperparameters

- Lengthscales theta: v1 fixes them at creation from the data scale
  (e.g. the median pairwise distance heuristic per designated column,
  computed once on standardized values), exposed as an argument.
  Slice-sampling theta per iteration is a mechanical follow-up (the
  grouped-effects sliceSampleOnce is a free function already) but
  changes the draw count per sweep, so it lands separately if wanted.
- k: fixed-k works exactly as the other leaves (prior variance scales
  by (scale/k)^2). The chi-k hyperprior needs a sum-of-squares
  contribution; the natural quantity is f' C^-1 f per leaf (the
  standardized magnitude of the drawn function). Whether chi-k over GP
  leaves is in v1 is an open decision; refusing updateK + gp() at
  first is the conservative default.
- Nugget: fixed small constant (1e-6 relative) for conditioning only;
  not exposed.

### Formats, state, and predict

The flat format's side channels extend from the linear-leaf precedent
(fixed numParams - 1 slopes per leaf) to variable-length blocks:

- flatten emits, per GP leaf, alpha = C^-1 f (n_leaf doubles) plus the
  n_leaf x q standardized covariate rows, with a per-leaf count vector
  alongside (lengths are recoverable from counts; validation mirrors
  the slopes split by (sizes + 1) / 2 leaves).
- State slots tree.params/saved.params carry the concatenation, as
  linear states do; restore validates counts against the tree records.
- addFlatPredictionsBelow computes c(x*)' alpha per row, standardizing
  on the fly with the stored constants - the same replay structure as
  the linear leaf's slope path.
- getTrees reports leaves with value = NA plus a gp marker column, or
  simply the leaf's mean value with the function riding only predict;
  open decision (reporting a whole function per row does not fit a
  data frame).

Storage note: keepTrees with GP leaves stores O(n_leaf (q + 1)) per
leaf per kept sample. That is the honest cost of a function-valued
leaf; the man page documents it next to the existing keepTrees memory
warning.

### Mutation surface

setPredictor/updatePredictor/per-observation sessions re-route
observations; the GP leaf regathers covariates under existing constants
(regatherTrainingCovariates precedent) and the transactional
snapshot/rollback of fits already covers function-valued leaves because
fits ARE the parameters. setData reinitializes constants like the
rebuilt cut grid. No new refusals expected beyond what linear leaves
impose today; stage 3 verifies each path with component tests rather
than assuming.

### Staging (each gated: tests/cpp, suite, equivalence identical draws
for existing models, speed, check)

1. Engine: FunctionLeafModel concept + GPGaussianLeaf (marginal, draw,
   fit paths; span suffstats; size cap; fixed theta/k). Component
   tests against R reference GP math (mvtnorm-style constants
   hardcoded, computed by an independent R implementation); recovery test on a
   smooth 1-d function where constant-leaf BART visibly steps.
2. Formats: flatten/state/keepTrees/predict with per-leaf counts;
   bitwise state round trip; saved-replay == recorded test fits.
3. R surface: node.prior = gp(columns, k, lengthscale) on dbarts();
   resolveLeafCovariates reuse; views/xbart composition; tinytest.
4. Follow-ups by demand: Cholesky caching, chi-k coupling, sampled
   lengthscales, low-rank kernels.

## Part 2: the non-conjugate MoveStrategy (designed-for, not scheduled)

For likelihoods with no working-Gaussian representation the chain's
(z, w) backbone and integrated branch scores both disappear. The
designed-for shape (core-generalization.md):

- Concept: NonIntegrableLeafModel exposes logLikelihoodForNode(params)
  and prior draws; no marginal.
- Moves: joint structure+parameter proposals - birth grows a rule from
  the prior and draws child parameters from the prior (or a Laplace
  approximation centered at the MLE), acceptance is likelihood ratio
  times prior ratio times proposal correction. Prior-draw proposals
  make the reversible-jump Jacobian trivial; Laplace centering is a
  mixing upgrade behind the same interface.
- The sweep: per-tree residuals no longer decompose; the likelihood is
  evaluated at total fits with the tree's contribution swapped, over
  the affected span. That is a second Chain::run variant, selected at
  compile time by the leaf concept - the largest single piece of this
  part.

Landing this without a consumer would be speculative engine work; the
recommendation is to write the second sweep only when a concrete model
(heteroscedastic BART, survival, count responses beyond Polya-Gamma)
is scheduled, using part 1's concept plumbing.

## Open decisions

1. Part 1 v1 scope: fixed lengthscales + fixed k + hard size cap
   (recommended), or any of those sampled/soft from the start.
2. Kernel family: squared-exponential only in v1 (recommended), or a
   Matern option up front.
3. getTrees reporting for function-valued leaves (NA value + marker
   recommended).
4. Whether part 2 waits for a named consumer (recommended) or lands as
   engine-only capability in this major-version window.
5. Timing relative to release: GP leaves are additive (new leaf model,
   existing paths untouched, equivalence must stay identical), but
   1.0-0 is otherwise release-ready; Vincent decides whether phase 6
   precedes or follows the CRAN submission.

## Stage 1 landing notes (2026-07-04)

Engine only; nothing is exposed through the facade or bridge, so the
shipped .so is semantically unchanged (equivalence: identical draws).
Deltas and facts vs the plan above:

- DELTA, over-cap guardrail: NOT the veto value. The veto deadlocks -
  every tree starts as a single root leaf holding all n observations,
  and a birth splitting one over-cap leaf into two over-cap children
  compares veto against 2 x veto, which no draw accepts, so trees could
  never grow at all. Instead leaves above gpMaxLeafSize delegate their
  marginal AND their draw exactly to ConstantGaussianLeaf. The scoring
  rule is a deterministic function of leaf membership, so MH
  comparisons stay coherent; oversized regions run constant-BART
  (data-informed splits throughout) and the GP refinement switches on
  once a leaf falls under the cap. The component test pins over-cap
  scores bitwise to the constant leaf's.
- Concepts (model.hpp): LeafModelCore additionally requires the
  hasFunctionParams trait (both existing models gained "= false");
  ScalarLeafModel adds !hasFunctionParams; new FunctionLeafModel
  requires beginTreeDraw(tree), drawFromPosteriorForNode(rng, tree, y,
  weights, k, sigmaSq, node, fits) -> FunctionLeafDrawStats
  (sumSquaredParams, numParams - the chi-k contribution),
  drawFromPriorForNode, and fitForTestObservationForNode;
  IntegrableLeafModel is the three-way disjunction.
- GPGaussianLeaf: squared-exponential correlation with fixed 1e-6
  nugget over standardized designated columns. The covariate gather
  duplicates the linear leaf's pattern deliberately (LinearGaussianLeaf
  stays textually untouched; factor a shared struct only if a third
  consumer appears). Lengthscales: supplied per column (standardized
  scale) or the median pairwise-distance heuristic - deterministic
  subsample of at most 512 rows (stride ceil(n/512)), R median
  convention, degenerate columns fall back to 1.
- Math as proposed: score = -0.5 log det V + 0.5 log det(sigma^2 W^-1)
  - 0.5 z' V^-1 z with V = (scale/k)^2 C + sigma^2 W^-1, one Cholesky
  (reduces to the constant leaf's formula for a constant kernel, the
  component test checks the limit numerically). Posterior draw by
  Matheron's rule f = f0 + s^2 C V^-1 (z - f0 - e0), consuming exactly
  2 n_leaf standard normals (f0's in row order first); empty leaves
  consume none. Prior draw s L_C eps, n_leaf normals. Prediction
  weights alpha = C^-1 f are cached per node (arena-indexed;
  beginTreeDraw resets per tree sweep) and test rows evaluate
  c(x*)' alpha; chi-k accumulates f' alpha over n_leaf coordinates
  (over-cap: value^2 over 1, the scalar accounting).
- SamplerOptions += gpLengthscales (borrowed; consumed) and
  gpMaxLeafSize (default 256). The designation shares
  leafCovariateColumns; the leaf model TYPE picks linear vs gp, so the
  facade needs a discriminator only at stage 3.
- Chain: the function-leaf branches next to the existing if constexpr
  pairs (scalar/vector text unchanged). Draws land in currFits_
  directly; test rows route per tree through the draw cache. Mutation
  flows: revalidate/recover leave params empty (fits ARE the
  parameters and stay in place; rebuildFitsFromParameters only
  regathers covariates), forceRefreshTrees collapses structure keeping
  per-observation fits, applyNewData cold-starts fits at zero through
  the scalar broadcast path (old fits are meaningless over replaced
  data; structures still carry through the split remap).
- Stage-1 refusals (graceful, the linear stage-2 precedent): run()
  skips flattening kept samples, getState returns an empty-tree state,
  stateIsValid/setState report false, flattenTree/printTree/
  printSavedTree come back empty or print nothing, and both predict
  paths zero their output. Sampler::numLeafCovariates and
  leafCovariateColumns serve function leaves.
- Tests (tests/cpp, constants computed in R by direct evaluation of the
  GP marginal and posterior formulas):
  testGPLeafMarginal (five marginals + median thetas vs R at 1e-9 /
  1e-12, constant-kernel limit, over-cap bitwise), testGPLeafDraw
  (posterior mean/variance vs R, E[f' C^-1 f], prior moments,
  duplicated-test-row identity, over-cap broadcast, empty leaf),
  testGPLeafEndToEnd (sin(4 pi x), n 250, m 10, cap 100: sse ratio
  0.014 beats the constant leaf's 0.016, sigma 0.211 vs truth 0.2,
  test-fit gap 0.0013, updateK sane, refusals, prior-draw
  continuation, and a SamplerFacade<GPGaussianLeaf> instantiation for
  full compile coverage of the surface).

## Stage 2 landing notes (2026-07-04)

Formats: flatten, keepTrees, state, and prediction now serve
function-valued leaves; still engine-only (no facade or bridge
exposure; equivalence identical draws). Facts vs the plan:

- Saved side channel (savedTreeParams, parallel to savedTrees): one
  variable-length block per leaf in pre-order - [count, constant] when
  count is zero (over-cap and empty leaves replay a constant),
  otherwise [count, alpha (count), plain standardized covariate rows
  (count x q, row-major, member order)]. Counts ride the stream as
  doubles (exact to 2^53); computeFunctionBlockOffsets (tree.hpp)
  walks and validates the channel and produces per-leaf offsets.
  Lengthscales are NOT baked into the rows - the replay divides - so
  the format survives a future sampled-theta follow-up unchanged.
- Replay: addFlatFunctionPredictionsBelow (tree.hpp) is the linear
  analogue with per-leaf block offsets threaded exactly like
  leafOffset. It standardizes test rows on the fly and repeats the
  live evaluation's arithmetic order (((x - mean)/sd - row)/theta,
  members in draw order), so saved replays BIT-MATCH the recorded test
  fits - the component test asserts equality. Its stack scratch bound
  maxFunctionLeafCovariates = 8 (GPGaussianLeaf::maxNumCovariates
  aliases it; the stage-3 factory validates).
- run() keepTrees flatten reads the DRAW CACHE
  (appendLeafBlockFromCache): the alpha weights recorded are the exact
  values the recorded test fits used. Between runs (no fresh cache),
  appendLeafBlock recomputes alpha = C^-1 f from the persisted fits
  against current covariates (flattenTree, predictFromCurrentTrees).
- FlatNode.value for a gp leaf = the leaf's mean fit
  (functionLeafValues, chain-side): reporting only - prediction reads
  the blocks, zero-count blocks carry their own constant - and the
  natural value column for stage 3 getTrees.
- Live-tree state channel: treeParams[t] is simply the tree's fits
  slab (n doubles, observation order; the fits ARE the parameters).
  getState copies it out, setState memcpys it back - bitwise by
  construction, no per-leaf bookkeeping. stateIsValid demands one
  n-slab per tree and walks every saved side channel against its
  tree's leaf count; restores of truncated channels or short slabs
  refuse without mutating (component test).
- initializeSavedTrees default function slot: {0, 0} (zero-constant
  block), the analogue of the scalar zero leaf. printTree/
  printSavedTree print the mean values through the scalar path.
- All stage-1 refusals are lifted; what remains engine-refused is
  nothing - stage 3 adds creation-time validation (factory) and the R
  surface, including the bridge state slots for the new channels.
- Tests: testGPLeafFormats (bitwise saved replay vs recorded test
  fits; bitwise state round trip into a fresh sampler with identical
  sigma/train/test continuation; truncated-channel and short-slab
  refusals; flatten emits walkable blocks); the end-to-end test's
  refusal checks became positive checks (state restores into its own
  sampler, flatten emits records, fresh live prediction sits at the
  response center).

## Stage 3 landing notes (2026-07-04)

The R surface: node.prior = gp(columns, k, lengthscale, max.leaf.size)
on dbarts() and xbart(). Facts vs the plan:

- Dispatch: SamplerOptions.gpLeaves selects SamplerFacade<GPGaussianLeaf>
  in createSampler/createSamplerOverStore behind the shared designation
  validation (both leaf models bound at 8 covariates); the third facade
  instantiation's size cost stayed within speed noise, as the linear
  landing's did. SamplerBase gained usesFunctionLeaves() (vtable
  change - --preclean), which the bridge keys its state and reporting
  layouts on.
- R: dbartsGPPrior (A_class.R; k + columns like linear, plus
  lengthscale and max.leaf.size); gp() constructor in dbartsPriors;
  resolveLeafCovariates is shared with linear (label-parameterized
  messages) and additionally recycles a scalar lengthscale per resolved
  column. Open decision 1 implemented per recommendation: v1 REFUSES a
  non-fixed k for gp priors at parsePriors ("gp node priors require a
  fixed 'k'"), which binary responses hit through their chi default -
  they must pass an explicit k. The engine's chi-k accumulation
  (f' C^-1 f) is implemented and component-tested, so lifting the
  refusal later is an R-level change only.
- Bridge: parseModel accepts either prior class into the shared
  designation, plus gp lengthscales and cap. GOTCHA hit: a NULL S4 slot
  arrives via Rf_getAttrib as the pseudo-NULL SYMBOL, not R_NilValue -
  test positively for the resolved REALSXP. setModel's designation
  refusal also compares the leaf-model kind. setState branches on
  usesFunctionLeaves(): live channels read as one n-slab per tree,
  saved channels split by walking each tree's per-leaf blocks
  (readFunctionTreeParams/readFunctionSavedParams); the engine's
  stateIsValid re-validates after.
- getTrees (open decision 3 per recommendation): function-leaf
  samplers emit NO coefficient columns and NA values on leaf rows (the
  function rides prediction only); rules keep their split values.
  plotTree labels gp leaves "gp".
- xbart: the node.prior mini-env gained gp; folds inherit
  standardization constants from the full data (the linear precedent)
  while HEURISTIC lengthscales recompute over each fold's rows -
  supply lengthscale explicitly for exact cross-fold calibration.
- rbart_vi/bart/bart2 build normal(k) internally, as with linear.
  Sparse and sparse-backed mixed columns refuse through the shared
  resolveLeafCovariates paths.
- Docs: dbartsPriors.Rd gained the gp entry (and linear's stale
  "unavailable in xbart" sentence was corrected), dbarts.Rd/xbart.Rd
  mention the specification, NEWS carries the feature bullet.
- Tests: tinytest test-gp-leaves.R (26 results): validation and
  fixed-k refusals, smooth-signal recovery, bitwise saved replay
  through predict, getTrees/plotTree semantics, bitwise state round
  trip via storeState/setState, predictor mutation on the designated
  column, probit composition, full-rows view bitwise vs raw path, fold
  prediction quality, xbart, sparse refusal.

## Stage 4 landing notes: chi-k coupling (2026-07-04)

The fixed-k restriction was R-side only - the engine has accumulated
f' C^-1 f over n_leaf (value^2 over 1 for over-cap leaves) into the
chi-k update since stage 1, component-tested in testGPLeafEndToEnd.
Lifting it removed the parsePriors refusal and nothing else; the
bridge's updateK/kHyperprior plumbing is leaf-agnostic. Binary
responses now get their usual chi(1.25, Inf) default under gp leaves,
and dbarts run results carry the sampled k as for the other leaf
models. xbart continues to reduce a k found inside the prior to a
per-cell fixed value only when it is numeric; a chi object passes
through to a sampled-k fit, matching normal(chi(...)) semantics.
test-gp-leaves.R's refusal checks became positive coverage
(continuous chi + binary default, sampled k varies and stays
positive).

## Stage 4 landing notes: kernel caching (2026-07-05)

The correlation kernel and its Cholesky factor depend only on the
leaf's member list (values and order) and the fixed lengthscales - not
on sigma, k, weights, or the response - so they are reusable across
sweeps for any tree whose structure survives its move. GPGaussianLeaf
now keeps a per-(tree, node) cache of both, and every lookup
re-validates by memcmp of the stored member list against
tree.indices - an O(p) check guarding O(p^3) work that makes a hit
bitwise identical to a fresh build and makes staleness structurally
impossible: there are no invalidation hooks in the move code to miss.
Specifics:

- The cached kernel is the exact buildKernel output (nugget included);
  a hit skips gather + buildKernel in the marginal and additionally
  chol(K) in the draws. chol(V) recomputes every draw - sigma moves
  every sweep - and is the remaining cubic term. Score-to-draw reuse
  of chol(V) within a sweep (sigma and k are fixed across one
  iteration's tree updates, and V does not depend on z) is measured
  second-order headroom, left undone.
- Trees are keyed by address (the chain's tree vector never
  reallocates after construction); nodes by arena index. Entries whose
  node died, stopped being a bottom node, or changed membership are
  evicted in beginTreeDraw - hygiene for the budget, not correctness,
  since lookups re-validate regardless.
- Covariate-value changes are the one thing member comparison cannot
  see, so reinitialize and regatherTrainingCovariates clear the cache;
  those are exactly the mutation flows that touch u_. The component
  test pins both surfaces bitwise against a cold state-restored clone:
  a WITHIN-BIN designated-column perturbation (members identical,
  kernels stale - only the clear saves you) and a cross-bin
  non-designated change (members re-route - only the comparison saves
  you). Note the test-design trap: mutations with updateCuts map the
  old grid, so a clone built fresh over mutated data has a different
  grid and legitimately diverges; the clone must be built over the
  ORIGINAL data, restored, and handed the identical mutation.
- Memory is hard-bounded: 256 MB total split evenly across chains at
  initialize, leaves under 32 observations never cached, and a spent
  budget degrades to recomputation. Entries charge members + kernel +
  factor at insert.
- Measured on n = 2000, one designated column, 50 trees, 200
  iterations: cap 256 5.85s -> 4.58s, cap 128 0.53s -> 0.41s, cap 64
  0.12s -> 0.12s; constant-leaf fits unchanged (0.06s). The draw
  phase's chol(V) dominates what remains. The stored speed baselines
  were unusable under an interactively loaded machine (same-binary
  reruns swung 1.04x-1.34x); the gate ran as a same-climate A/B
  against a fresh pre-cache recording and came in 0.963-1.038, no
  flags.

## Stage 4 addendum: sampled lengthscales are not mechanical
(2026-07-05, analysis only - nothing implemented)

The proposal calls slice-sampling theta "a mechanical follow-up".
Costing it while landing the kernel cache says otherwise; recording
here so the eventual decision is informed:

- Theta's conditional given trees and function values is
  p(theta | f) ~ prod over every gp leaf of N(f_leaf; 0, s^2 C_theta)
  times the prior: one O(p^3) Cholesky per under-cap leaf per
  evaluation, across ALL trees, times the 5-10 evaluations a slice
  step takes - several multiples of the whole draw phase, every
  iteration. A moving theta also empties the kernel cache each sweep,
  giving back the caching win.
- The saved format stores plain standardized rows with theta divided
  out at replay, so it survives - but replay needs the theta CURRENT
  AT EACH SAVED DRAW, i.e. a per-sample theta channel in the state
  and flatten formats, plus predict-side plumbing.
- Open decisions with no recommendation on record: the prior on theta
  (log-normal centered at the heuristic is the natural default), one
  theta shared per column vs per-column independence (current fixed
  code is per-column), and update cadence (every sweep vs thinned) as
  the cost lever.

## Status

Stages 1 (engine), 2 (formats), and 3 (R surface) LANDED, and stage
4's chi-k coupling and Cholesky caching are in; see the landing notes
above. Remaining follow-ups (sampled lengthscales - see the addendum,
they need decisions and are NOT cheap - low-rank kernels, and
score-to-draw chol(V) reuse) wait on demand, and the part-2
non-conjugate strategy remains designed-for; the open decisions were
implemented per their recommendations, with the over-cap fallback
replacing the veto in decision 1's scope.
