# GP leaves and the non-conjugate move strategy (phase 6)

Proposal, 2026-07-04. Phase 6 of docs/design/core-generalization.md is
"Non-conjugate MoveStrategy: GP leaves, general likelihoods". This
document separates those two axes, argues they land in that order, and
stages the work. Nothing here is implemented.

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
   local pieces).
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
   hardcoded, the linear_leaf_reference.R pattern); recovery test on a
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

## Status

PROPOSAL. Nothing implemented.
