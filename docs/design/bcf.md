# Bayesian Causal Forests: design

Status: LANDED. Proposed 2026-07-07; the two-forest sampler and the Forest
split it forced shipped 2026-07-07, with mixing refinements through 2026-07-10
(see "Landing" below). Model and defaults follow Hahn, Murray, and Carvalho
(2020), "Bayesian Regression Tree Models for Causal Inference" (Bayesian Analysis
15(3), 965-1056), and the `bcf` R package; deviations from either are flagged.

BCF fits two ensembles at once: a prognostic forest mu(x, pihat) and a
treatment forest tau(x), combined as

    y_i = mu(x_i, pihat_i) + b_{z_i} tau(x_i) + eps_i,   eps ~ N(0, sigma^2),

with z_i binary treatment. It is the concrete case the multi-forest provision
(core-generalization.md:138-144) was designed for, and bartCause is the intended
consumer, driving from R over the sampler mutation API. Scope: binary z and a
Gaussian response. Continuous treatment is out of scope (it needs a coefficient
function b(z), not a two-level scalar).

## The Forest member split

Chain today (src/bartcore/chain.hpp) holds the trees, their fits, the
backfitting residual, the leaf model, the split selector, and the response state
in one flat object - there has only ever been one forest. The split moves the
per-forest members into a `Forest` struct `Chain` holds one-or-more of
(`std::vector<Forest>`, size 1 for every current sampler).

Recommendation. A `Forest` owns: `trees_`/`indexBuffer_`; the leaf model and its
scale, plus the per-forest k / node-scale state (`k_`, the k accumulators and
hyperprior); the split selector (`treePrior_`, `fixedSplitProbabilities_` or the
DART state); backfitting state (`treeY_` residual, `treeFits_`/`totalFits_`
contribution, `paramByNode_`/`paramsByTree_`, the test-fit contribution); the
per-forest keepTrees channel (`savedTrees_` and companions); and the per-forest
options (tree count, base/power, move probabilities, node scale).

`Chain` keeps what is shared: the RNG (`rng_`), the response model (`response_`:
sigma, latents, the y-scaling, the working response/weights, the sigma draw),
`sigmaIsFixed_` / `family_`, and chain-level orchestration (chain count,
thinning, verbose). This is the design's division of labor
(core-generalization.md:107-116): a Forest carries its own working response and
backfitting state, the ResponseModel the global variance draw. BCF's glue
scalars (b0, b1, the optional prognostic scale a) live on the combining response
model, not any Forest.

The sweep becomes a loop over forests: each backfits against its residual
`treeY_` (for one forest, `y*` minus nothing, so bit-identical to today), then
the response model refreshes latents, draws sigma from the combined fit, and
draws the glue. The single-forest paths stay draw-neutral - the equivalence exact
mode gates it - and the refactor lands before any BCF math.

State serialization. `ChainStateData`'s tree channels (`trees`, `savedTrees`,
`treeParams`, `treeMasks`, per-forest k) gain a forest dimension, and BCF adds
b0/b1/a to the glue. One format-version bump, on the single version scheme
storeState already stamps (public-surface.md:142-150) - shared with
flat-format-v2 and state-continuation, so a reader checks one version and refuses
a mismatch by name. A single-forest state is a length-1 forest vector; no
cross-version migration.

## Propensity score

Recommendation. Propensity enters as an ordinary covariate of the prognostic
forest's view (`include_pi = "control"`, bcf's default and recommendation for
observational data): pihat is one more column mu(x, pihat) sees, an ordinary
ordinal predictor; the treatment forest does not see it by default. The caller
computes it - bartCause already fits richer propensity models than the engine
would offer, and caller-side keeps the embeddable-sampler design (a wrapper
resampling it between sweeps reinstalls the column through the mutation surface
below). The engine never estimates propensities.

Alternative, not default: pihat in both forests (`"both"`), for callers who want
the effect surface to bend on it too. Excluding it entirely is discouraged with
observational data - the covariate-dependent prior it induces is the paper's
mechanism against regularization-induced confounding.

## Per-forest defaults

BCF calibrates the forests differently; the note adopts its conventions, since
bartCause expects bcf-like behavior:

- Tree count: prognostic 200, treatment 50. The smaller treatment forest
  regularizes tau harder - the paper's point that homogeneous or mildly
  heterogeneous effects should shrink toward a constant. This prognostic
  default (200) is bcf's, above dbarts's own single-forest default of 75
  (prior-defaults.md); BCF carries bcf's counts, not dbarts's.
- Tree prior: prognostic base 0.95, power 2 (standard CGM decay); treatment
  base 0.25, power 3 - a much stronger pull to shallow trees, so tau defaults to
  near-constant effects and deepens only under real heterogeneity.
- Leaf scale: prognostic a half-Cauchy with median 2 sd(y), treatment a
  half-Normal with median sd(y) (`use_muscale` / `use_tauscale`) - adaptive scale
  hyperpriors, the role dbarts's `chi(1.25, Inf)` hyperprior on k plays for
  binary responses (prior-defaults.md).

Scale anchoring. The package range-anchors globally: y net of offset maps to
[-0.5, 0.5] by observed min/max, every constant calibrated against that range so
it transfers across the BART lineage (prior-defaults.md, "Response scaling").
`bcf` sd-anchors instead - it standardizes y by its weighted sd and states those
leaf scales in sd(y) units. Recommendation: keep range-anchoring globally and
express BCF's per-forest scales through the existing per-forest k / node-scale
and k hyperprior, converting the sd(y) medians to a range unit (approximate for
non-Gaussian y). VD scoped to this item whether the treatment forest should
instead sd-anchor exactly - a local deviation reproducing bcf's numbers
bit-for-bit; first open question below.

## Moderators: the treatment forest's column subset

The treatment forest is often given a subset of columns (effect moderators)
while the prognostic forest sees all of them plus pihat. This is the column-
subset view over a shared store from data-ownership.md's Sharing section: a view
is a column-index list (kernels consume one column at a time, so no contiguity is
needed), the forests attach through their own views, and a shared column mutated
once is visible to both under the single-writer rule.

Recommendation. mu attaches over (all x, pihat); tau over the declared moderator
subset (defaulting to all of x, no pihat). The store is built once and shared;
only a per-forest cut-grid override would allocate a diverging column, which BCF
does not need. BCF is the first multi-forest consumer of handle-plus-views.

## How the treatment enters

z is a 0/1 indicator; the treatment forest contributes only through the
treated-minus-control contrast. mu-versus-tau identification is the usual BART
centering (mean-zero leaf priors, the offset carrying the grand mean) plus a
scalar expansion `bcf` uses to fix the split of level between forests and mix
well.

Recommendation. Ship the bscale expansion: b_{z_i} = b0 when z_i = 0, b1 when
z_i = 1, with b0, b1 ~ N(0, 1/2), so the estimated effect is (b1 - b0) tau(x).
The two-coefficient form materially improves mixing on the treatment forest
(treated and control scales move independently in the Gibbs step) and is bcf's
shipped default. The prognostic scalar a (half-Cauchy) is a smaller win, and may
be deferred at a = 1.

Alternative for the first cut and the gate: fix b0 = 0, b1 = 1 (tau enters
directly as z * tau(x), no expansion) - the minimal coherent model, making the
gate a pure two-forest enumeration; the expansion is added afterward as the
default. Continuous z would replace the two-level b_z with a function b(z), out
of scope.

## Exact-posterior gate

The existing gates (benchmarks/R/logistic-reference.R, categorical-exact.R) take
a single-tree, one-predictor problem, enumerate the tree space under the CGM
prior, and quadrature the leaf marginals; the sampler must match the exact
posterior predictive to Monte Carlo error. BCF extends this to two forests: one
predictor, one prognostic and one treatment tree, binary z over a fixed design.

Recommendation. Enumerate each forest's single-tree structures independently (the
existing `enumerate` recursion, once per forest); the joint tree space is their
product. For each (mu-tree, tau-tree) pair the model is linear in the leaf
parameters,

    y_i = sum_l mu_l 1[i in mu-leaf l] + sum_m tau_m z_i 1[i in tau-leaf m] + eps_i,

so with Gaussian leaf priors and error the integrated likelihood over all leaf
parameters is a closed-form multivariate Gaussian block marginal - conjugacy
replaces the per-leaf quadrature the binary gates needed with an exact block. The
remaining quadrature is 1-D over sigma (low-dimensional over b0/b1 with the
expansion on); the first gate fixes b0 = 0, b1 = 1. The matched quantity is the
posterior predictive of (b1 - b0) tau(x) at the test cells; failure means the
two-forest backfit or the glue draw is wrong.

## Mutation surface

bartCause swaps response-side quantities between samples; the split decides
which entry points fan per forest.

- Chain-level (touch the shared response, fan to recompute each forest's
  residual): setResponse, setOffset, setWeights, setSigma, getLatents. y, sigma,
  and the offset are shared; each forest's `treeY_` is derived.
- Per-forest: setModel (each forest has its own base/power, k, node scale, split
  probabilities), setTreeStorage / getTrees / printTrees / predict (addressed by
  forest - a caller wants tau's trees or mu's), and the leaf-scale draw hook.
- Predictor mutation over shared views: setPredictor and family target a forest's
  view, but a mutation to a column both forests reference fans to both under the
  single-writer rule (collapsing the two-copy setPredictorJointly workaround,
  data-ownership.md Sharing).
- New for BCF: setTreatment(z). z is the one response-side quantity with no
  single-forest analog; updating it between sweeps re-forms b_{z_i} and both
  residuals. Installing pihat as a mutable prognostic column uses the ordinary
  predictor path.

C exposure stays internal first (bartcore helpers plus the bridge, as the data
handle does); bartCause drives from R, dbarts-level public exposure a later
decision (public-surface.md section 5's handle deferral).

## Open questions

1. sd- versus range-anchoring for the treatment forest. The recommendation keeps
   range-anchoring and maps bcf's sd(y) scales onto per-forest k hyperpriors -
   approximate. Sd-anchor the treatment forest exactly (a local deviation) to
   match bcf numerically, or accept the approximate map?
2. Prognostic tree count: adopt bcf's 200 for mu, above dbarts's own 75
   (recommended, since BCF carries its own calibration)?
3. The prognostic scalar a (half-Cauchy muscale): ship it with the b0/b1
   expansion, or defer at a = 1?
4. Should tau see pihat by default (`"both"`), or stay moderator-only?

## Resolutions (VD, 2026-07-07)

Questions 1 and 2: keep range-anchoring and map bcf's sd(y) scales
onto per-forest k hyperpriors (the approximate map); adopt bcf's
200/50 tree convention for the BCF entry point, leaving dbarts's
single-forest default at 75. Questions 3 (muscale expansion vs a = 1)
and 4 (tau seeing pihat) remain open; neither blocks the Forest
split, which is unblocked as of these resolutions.

Questions 3 and 4 (VD, 2026-07-07): ship the half-Cauchy prognostic
scalar a - parity with bcf is the goal, so both expansions
(use_muscale and use_tauscale equivalents) are defaults, as in bcf.
tau stays moderator-only by default (bcf's include_pi = "control");
pihat in tau's view is opt-in through the moderator subset. With both
expansions shipped, the adaptive per-forest scaling lives in the glue
(a for mu, b0/b1 for tau), which refines question 1's mechanism: the
per-forest k values are FIXED constants converted from bcf's range
units - no k hyperprior runs for BCF defaults; the approximate map is
the unit conversion only. The exact-posterior gate gains a low-
dimensional quadrature over the glue (its first mode fixes a = 1,
b0 = 0, b1 = 1; a second mode integrates the expansions). All BCF
steps are unblocked.

## Calibration (2026-07-07)

The map converts bcf's sd(y)-unit scales to the engine's range unit at
BCF creation. Let s be the sample sd of the range-scaled working
response (y net of offset mapped to [-0.5, 0.5]); range-anchoring makes
s the sd(y) magnitude in internal units.

- Prognostic. The mu forest's node scale is s, so the leaves of its
  trees sum to mu ~ N(0, s^2). The half-Cauchy scalar a has median
  |a| = aPriorScale = sd.control (default 2), putting the prognostic
  total a mu at a median sd.control sd(y) - bcf's use_muscale.
- Treatment. The tau forest's node scale is sd.moderate s / 0.674, so
  tau ~ N(0, (sd.moderate s / 0.674)^2). With b0, b1 ~ N(0, 1/2) the
  contrast b1 - b0 ~ N(0, 1) has half-normal median 0.674, so the
  effect (b1 - b0) tau sits at a median sd.moderate sd(y) - bcf's
  use_tauscale correction.
- Both forests fix k = 1; the map overrides the host model's node prior
  and k for mu, since the adaptive magnitude lives entirely in the glue.

The R wrapper exposes sd.control and sd.moderate (the two magnitudes in
sd(y) units) and converts internally. benchmarks/R/bcf-exact.R
reproduces the map and the range/sigma calibration to validate the
implementation end to end.

## Burn-in under strong prognostic signal (2026-07-10)

BCF fits on data whose prognostic amplitude is large relative to sigma
carry a long burn transient: the glue a reaches the right amplitude
within ~10 sweeps, but the mu forest's split STRUCTURE mixes slowly at
high SNR, and the shape misfit sits in sigma until it does. Measured
settle time scales with |a|/sigma and reaches ~72k sweeps in the
Cauchy(0, 2) prior's tail (|a| > 5; bias in E[sigma] up to ~1.2x at
burn = 18k). Amplitude-aware initialization and no-glue warm starts
were measured and do NOT help (raising early SNR freezes structure
further); burn length is the working lever. Normal use (|a| ~ O(1))
settles in ~2k sweeps; the strong-|a| regime arises mainly when a
response is swapped in against a stale build scale (the in-Gibbs
setResponse(updateScale = FALSE) path). The SBC harness pins its BCF
burn to absolute sweeps accordingly. The |a| >= 40 extreme tail (~3%
of the prior) is a structure-mixing limit no burn fixes; an engine
remedy (tempered early sweeps) would be its own item. Records:
docs/plans/bcf-sigma-residual.md.

## Landing (2026-07-07 to 2026-07-10)

The Forest member split and the two-forest sampler above shipped
2026-07-07: the sampler itself (d1ffb92), calibrated and gated
against the exact posterior, closing the core item (c9fd2fe;
docs/plans/forest-split-bcf.md); the state format gained a forest
dimension (4e6b206); and a sampler can warm-start from a donor fit's
forests (933eed8). Mixing refinements landed 2026-07-10: an
interweaving rescale move on the glue ridge (9617c94;
docs/plans/bcf-ridge-interweaving.md) and the sigma burn-in
calibration recorded above (docs/plans/bcf-sigma-residual.md).
bartCause is the intended consumer, driving from R over the mutation
surface described above.

## Status

LANDED. Two-forest sampler and Forest split 2026-07-07; mixing
refinements 2026-07-10 (see "Landing" above).
