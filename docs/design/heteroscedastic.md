# Heteroscedastic BART (HBART): design

Status: LANDED 2026-07-20 (arc 3775437..994ec7e; records 31c0204); predictor
mutation (setPredictor, transactional and per-observation) on a
heteroscedastic sampler LANDED 2026-08-12
(docs/plans/multiforest-predictor-mutation.md S3, a825263; see section 14).
Plan: docs/plans/multi-forest-models.md (the heteroscedastic half). Users
declare a second ensemble modeling the residual
variance as a function of predictors, so y_i = f(x_i) + s(x_i) eps_i with
eps_i ~ N(0, 1): the mean forest f is the ordinary constant-leaf ensemble, and
the variance forest s^2(x) = product of multiplicative-positive trees (Pratola,
Chipman, George, McCulloch, "Heteroscedastic BART via Multiplicative Regression
Trees", JCGS 2020, 29(2):405-417; arXiv:1709.07542). The variance forest is a
SECOND, differently-typed forest coupled to the mean forest through the
per-observation precision channel. rng class: posterior-changing when a variance
forest is declared, BYTE-neutral (no variance forest instantiated, mean and sigma
paths untouched) when it is not (section 8).

Landscape (surveyed 2026-07-18): rbart, the reference HBART R package, was
archived from CRAN 2024-01; stochtree ships maintained log-linear variance
forests for bart() and bcf() (Murray 2021, "Log-Linear BART"); bartMachine issue
#40 (an hBART merge request) has sat unresolved since 2021. No maintained R
package cleanly offers this today except stochtree's variance-forest feature -
concrete demand, partial coverage elsewhere.

## 0. The load-bearing finding: the variance forest is CONJUGATE

The TODO and the forest-combiner note (docs/design/forest-combiner.md) both
frame the variance forest's leaf as "multiplicative-positive, NON-integrable,"
requiring "the non-conjugate MoveStrategy the design anticipated but never built"
and specializing the monotone M-access seams (ParamScoringLeafModel,
TreeDrawLeafModel). Verified against the two live implementations, that premise
is WRONG, and correcting it collapses most of the anticipated work:

- HBART (arXiv:1709.07542 Section 3.3) puts a SCALED-INVERSE-CHI-SQUARED prior on
  each variance leaf, s^2_lk ~ chi^-2(nu', lambda'). Conditional on the mean fit
  and every OTHER variance tree, tree j's problem is a pure partition-of-variance
  with a conjugate prior: the leaf marginal is closed form (paper eq. 7, "depends
  on the data only via sum_i e_i^2") and the leaf full conditional is again
  scaled-inverse-chi-squared (paper eq. 6). This is EXACTLY the conditional
  conjugacy the mean forest already exploits (backfit one tree against the others),
  transposed from a location problem to a scale problem.
- stochtree (Murray 2021) writes the same model on the log scale, sigma_0^2
  exp{h(x)} with exp(leaf) ~ inverse-gamma - the identical conjugate family
  (inverse-chi-squared is inverse-gamma), differing only in calibration/reporting
  convention.

So the variance forest fits the EXISTING conjugate move machinery
(metropolisJumpForTree, moves.hpp:777; the per-leaf marginal sum in
logLikelihoodForBranch, moves.hpp:61-78) with a new INTEGRABLE leaf model. The
consequences, developed below:

- No non-conjugate MoveStrategy is built (fork 3). The monotone seams
  (ParamScoringLeafModel, TreeDrawLeafModel) are NOT specialized - they solve a
  DIFFERENT problem (dependent leaves), which the variance forest does not have.
  The "payoff of building them generically" does not land here; the honest reuse
  is the ForestCombiner and the weight channel, not the move seams.
- The genuinely new engine work is elsewhere: a fourth INTEGRABLE leaf kind with a
  scale suffstat and a scale marginal (section 3), the MULTIPLICATIVE within-forest
  residual roll (section 6), and - the real structural blocker - a Chain that holds
  TWO different leaf types (section 6), the "single-leaf-type forest vector" re-carve
  forest-combiner.md named as heteroscedastic's remaining work.

## 1. The HBART model and its reduction to our engine

**Model (paper eq. 2-3).** y_i = f(x_i) + s(x_i) eps_i, eps_i ~ N(0,1), so
y_i | f, s ~ N(f(x_i), s^2(x_i)) with

    f(x)   = sum_{t=1}^{m}  g(x; T_t, M_t)        (mean forest, constant leaves)
    s^2(x) = prod_{l=1}^{m'} h(x; T'_l, M'_l)     (variance forest, positive leaves)

Each variance tree contributes a positive multiplicative factor; the leaf value of
tree l on the region x falls into is s^2_lk > 0.

**Backfitting decouples both forests into conditional-conjugate pieces.** With user
case weights w_i, Var(y_i) = s^2(x_i) / w_i (the dbarts precision-weight convention;
GaussianResponse borrows w_i, model.hpp:2486,2498). The mean forest is untouched:
conditional on s(.), each observation has a known precision w_i / s^2(x_i), so the
mean update is the ordinary WEIGHTED Gaussian backfit the engine already runs
(weighted node suffstats, tree.hpp:510-518). The variance forest backfits
multiplicatively: to update tree j, define the mean residual e_i = y_i - f(x_i) and
the other-trees product s^2_{-j}(x_i) = prod_{l != j} h(x_i); then

    r_i := e_i sqrt(w_i / s^2_{-j}(x_i)) ~ N(0, s^2_jk)  for i in leaf k of tree j,

so leaf k's likelihood is chi^-2-conjugate through the single statistic
sum_{i in k} w_i e_i^2 / s^2_{-j}(x_i) (the w_i carries the case weight into the
scale suffstat, section 3). So both forests reduce to the two phases the
engine already runs per tree - a birth/death Metropolis move on structure, then a
leaf draw - and BOTH are conjugate. The mean tree's move/draw are byte-identical to
today (only the weights differ); the variance tree's move/draw use a new leaf model
(section 3) over the existing conjugate machinery.

**HBART's own sweep (paper Algorithm 1).** Per iteration: update all m mean trees
(conditioning on the current s), then all m' variance trees (conditioning on the
current f). Our run() sweep already loops forests in order (chain.hpp:668); the
variance forest is a second phase in that loop, guarded like BCF's combiner
(section 6).

## 2. Decision (fork 1) - variance-forest parameterization

The load-bearing fork. Three parameterizations:

- **A. HBART multiplicative, conjugate (recommend).** s^2(x) = prod of leaves,
  each leaf s^2_lk ~ chi^-2(nu', lambda') (paper Section 3.3). Conjugate: the
  structure move reuses the existing conjugate marginal machinery, the leaf draw is
  a scaled-inverse-chi-squared Gibbs step, and it targets the EXACT posterior. The
  prior calibrates to dbarts's existing sigma prior by mean-matching across the m'
  trees (paper Section 3.4): lambda' = lambda^(1/m'), nu' from matching the leaf-
  product's prior spread to the homoscedastic chi^-2(nu, lambda) - so a
  heteroscedastic fit with a constant s(x) reproduces the homoscedastic sigma prior
  exactly (this is the reduction gate, section 9). This is what rbart shipped and
  what the exact-posterior gate can close.

- **B. stochtree log-linear, conjugate.** log s^2(x) = sum of leaves, sigma_0^2
  exp{h(x)}, exp(leaf) ~ inverse-gamma (Murray 2021). Mathematically the same
  conjugate family as A up to the sigma_0^2 anchor and log-scale calibration.
  Reporting is on the log-variance scale and it carries an explicit global
  sigma_0^2. Equivalent in exactness and mixing; the only substantive difference is
  whether a global scale sigma_0 is drawn separately (B) or absorbed into the
  product's calibration (A) - the sub-fork in section 5.

- **C. Log-linear with a NORMAL leaf prior (non-conjugate).** log s^2(x) = sum of
  Gaussian leaves, matching the mean forest's leaf prior aesthetically (symmetric on
  the log scale). This is the ONLY parameterization that is non-integrable and would
  specialize the monotone M-access seam with a Laplace / pseudo-marginal move - i.e.
  the one the TODO's framing describes. Rejected: it is strictly harder (a
  non-conjugate move dbarts has never built, gp-leaves.md "Part 2: the non-conjugate MoveStrategy"), mixes worse
  (structure moves lose the Rao-Blackwellization the conjugate marginal buys), and
  neither live package uses it. There is no reason to pay for a non-conjugate move
  when the conjugate form is available, exact, and precedented.

**Recommendation: A.** Conjugate, exact, reuses the existing move, matches the
reference implementation the exact gate validates against, and mean-matches the
current sigma prior so the homoscedastic reduction is clean. B is an acceptable
alternative if VD prefers the log-variance reporting scale and an explicit global
sigma_0 (some users reason more naturally about a multiplicative log-variance
surface); the engine work is nearly identical, so this can stay a reporting/
calibration choice rather than an architecture one. C is rejected outright.

RESOLVED (VD, 2026-07-19): **A**, and proceed with the implementation now (chosen
over locking the 1.0-0 release first; heteroscedastic is the last big planned
feature and goes into 1.0-0). Fork 4 resolves to the weight channel, fork 6 to the
nullable distinctly-typed variance forest, all per the recommendations below. The
design was blind-critiqued (REJECT: the roll gate did not exercise the roll; the
variance suffstat dropped case weights), amended, and re-critiqued (ACCEPT, with the
conjugacy foundation verified against HBART eq. 6-7 / Section 3.4). Two implementation
riders from the re-critique: (i) the roll component test must also perturb the tree's
OWN leaf and assert the suffstat does not move (perturbing another tree scales both
sides equally and cannot alone distinguish s^2 from s^2_{-j}); (ii) the non-identified
a_c d_c ridge may mix slowly, so gate (b) may need long runs to hit MC tolerance - a
runtime/tolerance matter, budget for it.

## 3. Decision (fork 2) - the fourth leaf kind

RESOLVE as an INTEGRABLE scale leaf, `ConstantVarianceLeaf`, a new leaf model in
model.hpp beside ConstantGaussianLeaf (model.hpp:125). It is NOT a
ParamScoringLeafModel and NOT a TreeDrawLeafModel (model.hpp:106-122) - the monotone
seams solve dependent-leaf coupling the variance forest does not have. What it
declares:

- **Parameter.** One positive scalar per leaf, s^2_k (the multiplicative variance
  factor). Reported as such; predicted as a factor of the product.

- **Prior.** s^2_k ~ chi^-2(nu', lambda'), the calibrated per-tree scaled-inverse-
  chi-squared (section 2, paper Section 3.4). The existing chi^-2 draw primitive
  covers it: a scaled-inverse-chi-squared draw is scale / chisq with
  ext_rng_simulateChiSquared (random.h:137, = gamma(df/2, 2)). No new rng primitive.
  SCALE CONVENTION (pin, reconcile with dbarts): paper eq. 6 writes the posterior
  scale as nu' lambda'^2 + (data), so the paper's lambda' is SD-scale (lambda'^2 is
  variance-scale). dbarts's existing sigma prior ChiSquaredScalePrior
  (model.hpp:2287-2307) parameterizes on the VARIANCE scale: its `scale` field is a
  variance and its posterior scale is degreesOfFreedom*scale + SSR. The variance leaf
  reuses that variance-scale convention (its calibrated per-tree variance-scale is
  lambda'^2 = the Section-3.4 mean-matched analog of dbarts's `scale`), so at m'=1 it
  collapses to ChiSquaredScalePrior's own (nu, scale) exactly (section 9 gate a).

- **Sufficient statistic (the one genuinely new node stat).** The leaf marginal and
  draw depend on the node only through (n_k, sum_{i in k} w_i e_i^2 / s^2_{-j}(x_i)) -
  a COUNT and a WEIGHTED SUM OF SCALED SQUARED RESIDUALS (the w_i is the user case
  weight, dropped from an earlier draft; without it a weighted Gaussian fit targets
  the wrong variance scale - section 5). The mean leaf's cached node stat (sumWeights,
  sumWeightedResponse; Node, tree.hpp:172-175) does NOT carry a sum of squares - it
  was deliberately dropped because it cancels in the mean leaf's MH ratio
  (model.hpp:132-134). The variance leaf needs it, so it accumulates its own suffstat
  over the node's index span, the LinearGaussianLeaf precedent (per-call O(n_leaf)
  accumulation of its own blocks, model.hpp:704-707) - NOT a change to the shared Node
  struct. This keeps the mean path's node stat byte-identical (section 8). The
  per-observation divisor w_i / s^2_{-j}(x_i) is supplied through the sweep's weight
  scratch (section 6), so the leaf reads it as an ordinary per-obs weight.

- **Marginal (the structure-move score).** The closed-form scaled-inverse-chi-
  squared marginal of the leaf (paper eq. 7): an independent factor per leaf, so the
  branch score is the ordinary per-leaf SUM in logLikelihoodForBranch
  (moves.hpp:61-78) - the variance leaf provides logIntegratedLikelihoodForNode over
  its own (n_k, weighted sum of scaled squared residuals) suffstat and falls through
  the SAME conjugate path, no ParamScoringLeafModel override. The empty-leaf veto
  (moves.hpp:74) applies unchanged.

- **Draw.** drawFromPosteriorForNode returns a scaled-inverse-chi-squared draw
  chi^-2(nu' + n_k, [nu' lambda'^2 + sum_{i in k} w_i e_i^2 / s^2_{-j}(x_i)] /
  (nu' + n_k)) (paper eq. 6, lambda' SD-scale per the convention above) - the ordinary
  independent per-node draw the ScalarLeafModel path already dispatches
  (chain.hpp:4857-4874), NOT a TreeDrawLeafModel coupled sweep. drawFromPrior draws
  chi^-2(nu', lambda'^2).

So the leaf is close to `ScalarLeafModel` in SHAPE (independent per-node draw,
per-leaf marginal sum) but its math is a scale problem, not the (scale/k) location
problem the ScalarLeafModel signature bakes in (model.hpp:47-55: it passes k and
residualVariance, threads a Gaussian mean-explained marginal). The clean move is a
sibling concept, `ScaleLeafModel`, with the scale suffstat and the chi^-2
marginal/draw, dispatched at compile time like the location concepts - rather than
forcing the Gaussian signature to double as a scale leaf. This is a modest concept
addition, not a seam reuse; be honest in the plan that monotone's seams are not the
mechanism.

## 4. Decision (fork 3) - the MoveStrategy is CONJUGATE, not non-conjugate

RESOLVE: the variance forest uses the EXISTING conjugate move (metropolisJumpForTree
over an integrable leaf), NOT a new non-conjugate MoveStrategy. Section 0 and
section 3 establish the conjugacy; the move machinery requires only that the leaf
expose a closed-form marginal, which it does. No Laplace, no pseudo-marginal, no
prior-grown reversible jump.

The M-access seam (MoveContext::leafParams, moves.hpp:44-47;
logLikelihoodForBranchWithParams, model.hpp:571) is UNUSED by the variance forest.
That seam exists so a leaf whose leaves are mutually DEPENDENT (monotone's ordered
neighbors) can score a move against the frozen values of the other leaves. The
variance leaves are conditionally INDEPENDENT given the other trees - the coupling
to the other variance trees rides the multiplicative residual (the weight/divisor
the suffstat accumulates against), not a read of sibling leaf parameters. So the
score integrates every touched leaf out, per-leaf, exactly as the Gaussian leaf
does. This is a direct conflict with the TODO's framing, recorded in section 12.

## 5. Decision (fork 4) - per-observation residual variance channel

RESOLVE: fold s^2(x_i) into the WEIGHT channel for the mean forest, the route
forest-combiner.md explicitly reserved ("a future variance forest can route
into the WEIGHT channel rather than the additive location - the seam exists,
unused"). The mean forest sees effective weight

    w_i^mean = user_w_i / s^2(x_i),   with global residualVariance = 1 (sigma fixed),

so posteriorPrecision = sum_i w_i / s^2(x_i) - EXACTLY the heteroscedastic Gaussian
precision - with ZERO change to the constant leaf's math (it already forms
posteriorPrecision = sumWeights / residualVariance, model.hpp:141,157). The weight
channel is already per-sweep and already refreshed when it varies
(workingWeights + workingWeightsVaryPerSweep, model.hpp:2353-2358; the logistic
omega precedent, model.hpp:3537-3539). It does add TWO n-length per-sweep scratch
vectors (not "no new array", an earlier overstatement): the divided mean weights
w_i^mean, and the variance forest's own per-tree divisor w_i / s^2_{-j}(x_i) - the
BCF forestWeights/combined precedent (combiner.hpp:235,388-396), per-SWEEP scratch,
never per-observation dispatch inside a kernel.

**Requires Gaussian; refuses latent families (construction-time).** The weight
channel is not free real estate: the latent families ALREADY own
response_->workingWeights() - logistic/PG returns omega (model.hpp:3537-3539),
NBResponse and (future) ordinal likewise - and the run loop pulls
weights = response_->workingWeights() every sweep (chain.hpp:783). Routing s^2
through the same channel would COLLIDE with the family's own latent precisions. So a
variance forest requires ResponseFamily::gaussian and is REFUSED at construction for
probit, logistic, nbinom, and ordinal, checked at the factory against the family
argument (facade.hpp:765,846,868) before any Chain is built - the same
construction-time refusal shape monotone uses for categorical columns
(facade.hpp:743). Only the Gaussian family, whose workingWeights() returns the
borrowed user weights (model.hpp:2498), leaves the channel available to carry
w_i / s^2(x_i). (Heteroscedastic probit/logistic - modeling a variance surface under
a latent link - is a real model but needs a distinct latent-plus-scale plumbing out
of v1 scope, section 11.)

Tradeoffs vs a DEDICATED sigma-vector (a per-observation residualVariance threaded
through every leaf marginal/draw):

- Mean-forest draw: weight channel = no signature change, reuses every weighted
  kernel and the cached node suffstat. Sigma-vector = every logIntegratedLikelihood
  / drawFromPosterior takes a per-obs variance, breaking the cached-suffstat
  reduction (the node can no longer cache sumWeightedResponse under one scalar
  variance) and touching the hot path. Weight channel wins decisively.
- Sigma update: there is NO global sigma draw under heteroscedastic - the variance
  forest IS the variance (parameterization A, section 2), so residualVariance is
  fixed at 1 and drawSigma is skipped (sigmaIsFixed_ true for the mean side,
  chain.hpp:790). Under parameterization B a global sigma_0 IS drawn, from the
  scaled residuals, reusing drawSigma; this is the section-2 sub-fork. Recommend A
  (no global sigma draw): it avoids the sigma_0-vs-variance-forest-level
  identification that B must pin, and the calibration (Section 3.4) already anchors
  the product's overall level. If B is chosen, a level-centering move on the
  variance forest (the multinomial precedent, multinomial.md) pins sigma_0 against
  the product's geometric mean.
- Prediction: the mean fit f(x) is unperturbed; s(x) is the variance forest's own
  combined fit, reported on its own channel (section 6). Both replay from saved
  trees independently (the saved variance trees became part of the serialized
  state in section 16; before that they replayed only through a live pointer).

## 6. Decision - two leaf types in one Chain, and the coupling owner

This is the real structural work, and it is NOT solved by the ForestCombiner as it
stands. `Chain` is `template <IntegrableLeafModel L>` with `forests_` a
`std::vector<Forest<L>>` (chain.hpp:286-287) - EVERY forest shares ONE leaf type.
BCF and multinomial exploit this (all their forests are ConstantGaussianLeaf); the
combiner is `ForestCombiner<L>` over that single L. Heteroscedastic breaks it: the
mean forest is ConstantGaussianLeaf, the variance forest is ConstantVarianceLeaf.
Three ways to hold both:

- **(a) Second Chain template parameter** `Chain<LMean, LVar>`. Precise but the
  widest blast radius: every Chain<L> reference, the facade instantiations, and the
  combiner become two-parameter. Risk to byte-neutrality (section 8) is high unless
  LVar defaults to a `void`-like tag that compiles the whole variance path out.
- **(b) A distinct typed member** `std::optional<VarianceForest> varianceForest_`
  beside forests_, where VarianceForest is a concrete (non-template) object over
  ConstantVarianceLeaf holding its own trees, leaf, and sweep. Chain<L> stays
  monomorphic on the MEAN leaf; the variance object is a fixed type (the variance
  leaf is always constant-scale in v1). Null off heteroscedastic - the BCF
  combiner_ null-short-circuit shape (combiner.hpp:87-101), byte-neutral by
  construction.
- **(c) Type-erased** `std::unique_ptr<VarianceForestBase>`: (b) behind a small
  virtual surface (score-and-move a tree, draw leaves, form combined variance). Its
  O(n_leaf) work lives inside its own methods, so the per-tree/per-node virtual call
  is amortized within the dispatch tiers (core-generalization.md, "per node
  op: virtual or switch"). Cleanest for keeping the variance leaf's translation unit
  and move instantiation self-contained.

**Recommend (b) or (c)** - a nullable variance-forest member, not a second template
parameter. The variance forest is a fundamentally different object (a scale ensemble
with a multiplicative sweep), not another Forest<L>; forcing it into the
single-L vector is what (a) fights. (c) is preferred if the variance leaf's move
instantiation is best isolated in its own TU; (b) if the extra vtable is unwanted
and the compile-time coupling is acceptable.

**The coupling owner.** BCF/multinomial route coupling through ForestCombiner<L>,
but that hierarchy is templated on the SINGLE mean L and its vector holds Forest<L>
only - it cannot hold the variance forest. So heteroscedastic's coupling is
Chain-level, guarded like `if (combiner_)`:

- A new `if (varianceForest_)` phase in run() (chain.hpp:668 loop): after the mean
  forest sweep, run the variance forest's own tree loop, then refresh the combined
  variance s^2(x_i) and hand the mean forest w_i^mean = user_w_i / s^2(x_i) for the
  next sweep (a ResponseModel decorator or a direct weight-vector formation, the
  formForestResponse weight route of section 5).
- The variance forest's OWN sweep needs a MULTIPLICATIVE running residual: the
  engine's rollTreeResidual (chain.hpp:713) maintains an ADDITIVE residual (treeY =
  y minus other trees' fits). The variance tree instead needs, per tree j, the
  divisor s^2_{-j}(x_i) = s^2(x_i) / h_j(x_i) - maintained by dividing the combined
  variance by the current tree's factor, a multiplicative analog. This is the one
  new per-tree sweep operation, local to the variance forest object.

Whether this lives as a HeteroscedasticForestCombiner subclass (widening the
combiner to hold a differently-typed forest - a real generalization of the
hierarchy) or as a dedicated Chain-level integration (the BCF-era shape, a member
plus guarded phases) is a design sub-fork. Recommend the dedicated Chain-level
integration for v1: the combiner's L-templated vector is the wrong container for a
second leaf type, and forcing heteroscedastic through it would generalize the
hierarchy for one consumer before hurdle (the other two-leaf-type model) sharpens
the shape - the same "let the second consumer shape the split" discipline the
combiner extraction followed (forest-combiner.md).

## 7. Decision (fork 5) - the R surface

Declare the variance forest on dbarts()/bart2() with a dedicated argument group; no
collision with the mean-side sigma/tree args (sigest, sigdf, sigquant, k, n.trees,
power, base, R/dbarts.R).

- **Which predictors drive variance.** Default the SAME x as the mean forest (the
  common case and stochtree's default). Allow a SUBSET or a distinct set via a
  formula/spec, the BCF-moderator precedent (a second design matrix): recommend a
  `variance` argument taking a one-sided formula or a column selection,
  `variance = ~ x1 + x2`; absent, it defaults to all mean predictors. HBART allows
  a distinct set and it is cheap to support (a second column view over the shared
  store).

- **Number of variance trees.** `n.trees.variance` (default a small ensemble, e.g.
  40, HBART/stochtree use fewer variance than mean trees). Clear, parallels
  `n.trees`, no collision. stochtree's `num_trees` under `variance_forest_params`
  is the ecosystem analog; `n.trees.variance` is the dbarts lower.case.dotted form.

- **Variance prior hypers.** The (nu', lambda') calibration derives from a
  homoscedastic-style spec: recommend `power.variance` / `base.variance` for the
  variance tree prior and a `k.variance`-or-df spec for the leaf prior spread,
  calibrated through paper Section 3.4 from the same sigquant/sigdf vocabulary the
  mean sigma prior uses. Keep the names suffixed `.variance` so they never shadow
  the mean-side hypers.

- **Naming rationale (collision trap).** Do NOT overload `sigma`/`sigest` (scalar,
  mean-side) or add a bare `heteroscedastic = TRUE` flag that then needs every hyper
  anyway. A single `variance` (predictors) + `n.trees.variance` + `.variance`-
  suffixed hypers is the least surprising, and `variance` as a predictor selector
  reads naturally ("model the variance as a function of ..."). bart2 takes the same
  arguments; the fit object gains an `s.test`/`sigma.test` channel (below).

- **Reporting.** The fit carries the mean fit as today plus a per-observation s(x)
  (or s^2(x)) channel - train and test - and per-variance-tree variable counts. This
  is a NEW reporting channel, NOT the multinomial numReportedLocations widening
  (combiner.hpp:328): that seam widens the SAME-typed response location K-fold through
  combinedFits/refreshLatents (multinomial.md), whereas the variance surface is a
  SEPARATELY-typed forest routed through the weight channel, never through response_,
  so it is reported from the variance forest's own combined fit directly. predict
  returns both f and s.

## 8. Bitwise neutrality for homoscedastic fits (fork 6)

Binding requirement (rng class: neutral when no variance forest is declared). A
homoscedastic fit MUST be byte-identical to today. Construction-time mechanism,
mirroring monotone (monotone.md) and BCF (combiner.hpp:87):

- No `variance` / `n.trees.variance` -> the factory (facade.hpp:756 createSampler)
  builds the UNCHANGED single-forest ConstantGaussianLeaf Chain with
  varianceForest_ null (or the combiner absent). Not a variance-capable type with
  the forest switched off - the SAME object code, so the identical rng stream. The
  variance sweep phase, the multiplicative roll, and the weight-division are all
  behind `if (varianceForest_)`, never entered.
- The mean leaf's node suffstat is unchanged (the variance leaf accumulates its own,
  section 3), so no shared hot-path struct widens. The weight channel already exists;
  homoscedastic weights (user weights or null) flow exactly as today.
- Binary size grows by the variance-leaf instantiation and (option a) any second
  template parameter - bounded like the existing leaf matrix
  (core-generalization.md); option (b)/(c) add one concrete type, not a matrix
  row.

Gates that prove neutrality (section 9): the equivalence trio (equivalence.R vs
equivalence-<hash>.rds, 22/22 IDENTICAL) reproduces byte-identically with NO
re-record - heteroscedastic adds no draw to any existing family or leaf; and
bench-sampler compare shows no regression on the homoscedastic constant-leaf paths.

## 9. Gates (fork 7)

Heteroscedastic has no single cheap exact-posterior enumeration, but the CONJUGACY
(section 0) is what makes exact and reduction gates possible at all - both leaf
marginals are closed form. The gate that MATTERS is (b): it is the only one that
exercises the top correctness risk (the multiplicative divisor s^2_{-j}, section 13).
Established style: logistic-reference.R, categorical-exact.R, multinomial-exact.R.

**(a) Homoscedastic reduction (a sanity floor, NOT a roll gate).** A variance forest
with n.trees.variance = 1 over a CONSTANT variance predictor can only stay a root, so
s^2(x) = a single chi^-2(nu', lambda') leaf = the homoscedastic sigma^2 prior under
the Section-3.4 calibration. At m'=1 that calibration collapses ALGEBRAICALLY -
lambda' = lambda^(1/1) = lambda and nu' = 2/[1-(1-2/nu)^1] = nu - so the
heteroscedastic sampler's sigma^2(x) = const posterior must match the homoscedastic
sigma^2 posterior distributionally (confirm the algebra NUMERICALLY, not just on
paper). Honest limitation: at m'=1 the divisor s^2_{-j} is identically 1, so this
gate does NOT exercise the multiplicative roll at all - it validates only the
calibration, the scale suffstat, and the single-leaf draw. It is a floor, not the
closing gate.

**(b) The m'=2 CLOSING exact gate (the top-risk gate).** TWO variance trees over one
binary predictor (two cells), plus a mean forest, sigma-free (the variance forest IS
the variance). Now s^2(x) = h_1(x) * h_2(x), so scoring tree 1 must divide by tree
2's contribution and vice versa - the divisor s^2_{-j} is genuinely NON-trivial and
the roll is on the critical path. Enumerate the joint tree-structure space (each
variance tree a root or a two-cell split; the mean tree likewise); for each, integrate
the variance leaves (v1 for tree 1, v2 for tree 2) by per-cell quadrature over their
product, with the mean leaves - Gaussian-conjugate given the per-observation variance
= the covering product - integrated ANALYTICALLY (closed form), the
multinomial-exact.R nested-quadrature precedent (multinomial.md). Match the
sampler's posterior mean of f(x) AND s(x) to MC + quadrature error, never widened.
This is the repo's bar - a CLOSING exact gate for the top risk, the class of gate that
caught the stationary-bias in the monotone arc (monotone.md). A self-consistently
wrong divisor (using s^2 or s^2_{-j'} for the wrong j') fails here, where it passes
(a) and (c).

**(c) Single-tree coupling gate (m'=1, coupling and growth only).** ONE mean tree +
ONE variance tree, one binary predictor: enumerate the 2x2 joint structure and match
posterior f(x), s(x) by nested quadrature. Validates the mean<->variance coupling and
tree growth cheaply, but - like (a) - at m'=1 the divisor is 1, so it does not guard
the roll; (b) subsumes its coverage and adds the divisor.

**(d) Simulation recovery / SBC (validates the sampler against ITS OWN model).**
Simulate from a known heteroscedastic truth (s(x) a ramp or step over x); the fit
recovers f(x) and s(x) with calibrated intervals, and s(x) tracks the truth where a
homoscedastic fit reports a flat inflated sigma. SBC ranks on the variance-leaf draws
check the scaled-inverse-chi-squared update is self-consistent. HARD LIMITATION: SBC
validates the sampler against the model IT implements, so a divisor bug that is
self-consistent between the data-generating and inference code passes SBC silently -
this is why (b), an INDEPENDENT quadrature, is the load-bearing gate and SBC is only
family-level smoke.

**(e) Reference cross-check.** Against stochtree (maintained) and/or the archived
rbart on a shared dataset - the used-car / Boston-housing examples from the HBART
paper. Informational (calibration conventions differ), not a bitwise gate, but a
gross s(x) mismatch falsifies the parameterization or calibration.

**Component tests (rng-free where possible).** The scale marginal and draw match
direct quadrature on fixed (n_k, weighted sum of scaled squared residuals); the
scaled-inverse-chi-squared draw's moments match closed form (a deep-tail df exercises
the gamma primitive); and - the divisor guard - a test that asserts the marginal/draw
actually USE s^2_{-j}: hold tree j's data fixed, change ANOTHER variance tree's leaf
values, and confirm tree j's suffstat and posterior draw move by the predicted
1/s^2_{-j} factor (not merely that a full s^2 = h_1 h_2 reconstruction is numerically
right, which a divisor that cancels its own error would also satisfy).

**What would FALSIFY correctness:** (b) an m'=2 exact-gate f/s mismatch beyond
MC+quadrature -> the coupling, the scale marginal, or the divisor is wrong (the
primary falsifier); the component divisor test failing -> the draw ignores s^2_{-j};
(a) a constant variance forest that does NOT reduce to homoscedastic -> the
calibration or suffstat is wrong.

**Neutrality trail.** equivalence.R compare IDENTICAL for every existing scenario
(no re-record); bench-sampler compare unchanged on the homoscedastic paths; a NEW
heteroscedastic equivalence fixture records f, s, and the reported channels as the
posterior-changing baseline for this arc.

## 10. Commit decomposition (fork 8)

Four independently gate-able units. Line estimates are engine C++ before R unless
noted, and are honest (this arc is larger than monotone's - the structural
two-leaf-type change dominates):

- **C1 - the fourth leaf kind (the scale leaf, in isolation).** ScaleLeafModel
  concept + ConstantVarianceLeaf: the scale suffstat (own-span accumulation), the
  chi^-2 marginal, the chi^-2 draw, the Section-3.4 calibration. Gate: tests/cpp
  component tests vs quadrature (marginal, draw moments), NO Chain wiring yet.
  ~350-500 lines.

- **C2 - Chain integration + neutrality (RE-BUDGETED UP).** The nullable
  variance-forest member (section 6 option b/c), the multiplicative residual roll,
  the variance sweep phase, and the weight-channel coupling (w_i^mean =
  user_w_i / s^2) with its two n-scratch vectors and the gaussian-only refusal.
  HONEST re-budget: only metropolisJumpForTree/logLikelihoodForBranch reuse across
  leaf types (section 4); the variance forest reimplements its OWN per-sweep
  plumbing - the multiplicative roll (no rollTreeResidual reuse, that is additive),
  scale-leaf draw dispatch, saved-tree flattening, test-fit routing, and the
  mutation/rebuild paths (revalidateTrees/applyNewData are forest-0-only today,
  forest-combiner.md) - because that surrounding chain.hpp plumbing is templated
  on the mean L and does not transfer to a differently-typed forest. Homoscedastic
  byte-neutral. Gate: equivalence 22/22 IDENTICAL (no re-record) + a new
  heteroscedastic fixture identical to itself + tests/cpp (incl. the divisor guard).
  ~800-1300 lines, not 500-800.

- **C3 - reporting + the exact/reduction gates.** The s(x) train/test channels (a
  NEW separately-typed-forest reporting channel, NOT the numReportedLocations
  same-typed widening, section 7), predict for both forests, state (a NEW scale-leaf
  FlatNode branch encoding a positive variance value - NOT same-type per-forest-list
  reuse; the value semantics and validation differ from a Gaussian leaf value), and
  benchmarks/R/heteroscedastic-exact.R (gates a, b, c - the m'=2 closing gate is
  here). ~400-600 lines + gate scripts.

- **C4 - the R surface.** `variance`, `n.trees.variance`, the `.variance`-suffixed
  hypers on dbarts()/bart2(); the bridge SamplerOptions fields and the gaussian-only
  refusal message; the fit class and probability/interval generics; tinytest
  (including the recovery test, gate d). ~400-600 lines R + bridge.

Total ~1950-3000 lines (up from the earlier estimate once C2's own-plumbing reality
is counted). C1 is a clean standalone; C2 is the risk-bearing structural change and
the widest review; C3/C4 are additive on top.

## 11. Out of scope, and the doors

- **Non-constant variance leaves** (linear or GP variance leaves). Constant scale
  leaves only in v1; a smooth log-variance surface is a door, not v1 scope.
- **Parameterization B's global sigma_0 with level-centering** (section 2/5): if VD
  prefers the log-variance reporting scale, this is a bounded follow-up, not a
  rearchitecture.
- **The non-conjugate log-Gaussian variance leaf** (fork 1 option C): rejected;
  reopen only if a use case needs a Gaussian log-variance prior the conjugate form
  cannot express. It is the one path that would specialize the monotone move seam.
- **Heteroscedastic composed with BCF / grouped effects / multinomial**: the
  weight-channel routing composes in principle (a GroupedResponse decorates below the
  variance weighting), but building any cross-product is out of scope. CORRECTION
  (docs/design/hurdle.md, landed 2026-07-20): hurdle did NOT end up sharpening
  section 6's structural choice - it composed R-side as two independent ordinary
  fits (occupancy + lognormal positive part), no engine change, so section 6's
  Chain-level two-leaf-type integration remains a single-consumer design with no
  second consumer to generalize against.
- **dbarts.h exposure:** none in v1 (the ordinal/monotone precedent); heteroscedastic
  is reachable only through the R surface and internal bridge.

## 12. Plan-vs-code note (the TODO conflict)

The TODO's multi-forest-models entry and forest-combiner.md state heteroscedastic
needs "a fourth leaf kind (multiplicative-positive, non-integrable) plus the
non-conjugate MoveStrategy the design anticipated but never built," and that "the
generic move-scoring-reads-M seam ... now EXISTS ... the non-conjugate MoveStrategy
specializes it (Laplace/pseudo-marginal)." FINDING: for the two LIVE
parameterizations (HBART's multiplicative inverse-chi-squared and stochtree's
log-linear inverse-gamma) the variance leaf is CONJUGATE and INTEGRABLE (section 0),
so:

- The fourth leaf kind IS needed, but as an INTEGRABLE scale leaf, not a
  non-integrable one (section 3).
- The non-conjugate MoveStrategy is NOT needed and should not be built for
  heteroscedastic (section 4). The variance forest reuses the EXISTING conjugate move.
- The monotone M-access seam (ParamScoringLeafModel) and the tree-granularity draw
  seam (TreeDrawLeafModel) are NOT specialized - they solve dependent-leaf coupling
  the variance forest does not have. The "payoff of building them generically" does
  not land here; a non-conjugate move would only be needed for the rejected
  Gaussian-log-leaf variant (fork 1 option C).

The real, under-anticipated work is the Chain-level two-leaf-type change and the
multiplicative residual roll (section 6), which the TODO/forest-combiner note DID
flag as heteroscedastic's remaining re-carve (forest-combiner.md) but did not
connect to the (mistaken) non-conjugacy framing. Net: the arc is simpler in its
statistics (no new MoveStrategy) and harder in its plumbing (two leaf types, a
multiplicative sweep) than the TODO reads.

This note SUPERSEDES forest-combiner.md ("a non-integrable leaf and its own
MoveStrategy") and the TODO's matching clause on the non-conjugacy point; those are
corrected at landing (the orchestrator owns the edit to forest-combiner.md and the
TODO). Cite this note, not them, for the variance-forest's conjugacy.

## 13. Costs, risks, and confidence

- **Per-fit cost.** A second ensemble roughly doubles the tree work, plus the
  multiplicative roll and the weight refresh - honest, and paid only when a variance
  forest is declared; homoscedastic fits are untouched (section 8).
- **Biggest correctness risk: the multiplicative residual roll (section 6).** A wrong
  divisor (using s^2 instead of s^2_{-j}, or a stale combined variance) silently fits
  the wrong variance model, and unlike a crash it produces plausible-looking s(x).
  Neither m'=1 gate (a, c) exercises it (divisor 1), and SBC (d) passes a
  self-consistent divisor bug. The GUARDS are the m'=2 closing exact gate (section 9
  gate b, an INDEPENDENT quadrature) and the component divisor test that asserts the
  draw actually reads s^2_{-j}; both must run with m' >= 2. This gate must exist and
  close before C2 lands, matching the repo bar the monotone arc met (monotone.md).
- **Second risk: the two-leaf-type Chain change** perturbing the homoscedastic path.
  The equivalence 22/22 bitwise gate is the guard; option (b)/(c) over (a) keeps the
  blast radius off Chain<L>.
- **Implementation-time risk to MEASURE: variance-forest mixing under param A.**
  Param A carries no global sigma_0 and no level-centering move, so the variance
  forest's overall level is pinned only by the calibrated prior. Whether that mixes
  adequately at m' >= 2 without a centering move is ASSERTED, not confirmed here -
  the HBART paper's mixing evidence at its own m' was not verified for this note.
  Measure IACT/ESS of s(x) on the recovery data at landing; if the level random-walks,
  the multinomial-style level-centering move (multinomial.md) or param B's explicit
  sigma_0 draw is the reserved remedy. Cheap to add later; flagged now so C2/C3 leave
  room for it.
- **Confidence in the variance-forest DRAW as specified: HIGH.** The scaled-inverse-
  chi-squared leaf marginal and draw are textbook conjugate updates
  (paper eq. 6-7), the primitive exists (random.h:137), and both live packages agree
  on the family. There is no non-conjugate approximation to get wrong; the residual
  numeric care is the SD-vs-variance-scale convention (section 3), pinned against
  dbarts's own ChiSquaredScalePrior.
- **Confidence in the MoveStrategy as specified: HIGH.** It is the EXISTING conjugate
  move over a new integrable leaf - no new acceptance math, only a new per-leaf
  marginal in the same per-leaf-sum path (moves.hpp:61-78). The confidence rests on
  the conjugacy claim (section 0), which is verified against two independent
  implementations rather than recalled; the residual risk is engineering (the roll,
  the weighted suffstat, and the two-leaf plumbing), not algorithmic.

## 14. Post-landing: predictor mutation on a heteroscedastic sampler

Sections 0-13 cover the model's own landing (2026-07-20). Everything below is a
later arc, docs/plans/multiforest-predictor-mutation.md S3 (2026-08-12), which
did not touch the model - no new draw, no changed prior, no changed sweep -
only widened which store mutations the variance forest can survive.

Before S3, `$setPredictor` (transactional, whole-matrix or column-subset via
its `column` argument, and per-observation via `forceUpdate = "partial"`)
refused outright on any sampler carrying a variance forest; only the FORCED
whole-matrix swap, `setCutPoints`, and `setData` re-routed it (section 6's
Chain-level integration covers those paths, which are unchanged). S3 splits
`refreshVarianceForest` - collapse-and-remap, still exactly what
`forceRefreshTrees` and `setData`'s `applyNewData` use - into two halves for
the transactional paths only: `revalidateVarianceTrees` recovers each tree's
node-indexed factors through the LIVE partition, repartitions, and reports
occupancy without collapsing or scattering anything; `rebuildVarianceFactors`
then drops stale missing directions, scatters the recovered factors through the
new partition, and recomputes `combinedVariance`. `Chain::repartitionTrees`
(the transaction's rollback route) gained a variance arm to match, so a
rejected transaction restores the variance trees' partitions exactly rather
than leaving them re-routed by the declined proposal.

The acceptance criterion is the one this document's variance forest has always
needed and did not fully have: `Chain::stateIsValid`'s variance branch checked
tree well-formedness and strict leaf positivity but not occupancy, so
`$setState` could install a variance state with an unoccupied bottom - one that
carries no drawn scale - where its mean-forest sibling had refused the
analogous state all along. S3 closes that gap: the variance branch now builds
the same scratch-tree-and-repartition check the mean branch already ran, so
`$setState` refuses such a state too. This is a behavior change (a state that
used to install now does not) and ships its own NEWS bullet and tinytest.

Net effect: a heteroscedastic sampler now accepts `$setPredictor` - whole
matrix, column subset, and per-observation alike - and the per-observation
session across samplers (`updatePredictorPerObservationJointly`) on the same
terms as an ordinary single-forest sampler - a row installs only if it
empties no leaf in any mean tree AND no variance tree, rolling back (or
declining, per row) otherwise. See
docs/design/empty-leaf-veto.md for why that criterion is not new, and
docs/plans/multiforest-predictor-mutation.md for the full mechanism, the
falsifiers that gate it (F9, F10 cover the variance-specific rollback and
recovery-ordering claims), and the equivalence baselines.

## 15. Post-landing: Student-t residuals and grouped random effects refused (2026-08-17)

Two more gaussian-family compositions reach a variance forest's construction
path unadjudicated, since neither the scale-mixture reweighting `student()`
installs nor rbart_vi's per-group Gibbs blocks were checked against the
variance forest's weight-channel routing (section 5) when either landed:
`resid.dist = student()` residuals declared alongside `variance` constructed
and ran on the pre-refusal build and is now a validation error at the family
gate's site (R/spec.R); a grouped (rbart_vi) fit alongside `variance` was
already unreachable (rbart_vi declares no `variance` formal and its `...` is
rejection-only), so no new refusal was needed there. Both are validation
errors only, not a withdrawal of the `variance` formal from any surface -
support arrives by adjudicating whether the weight-channel divisor composes
with the scale-mixture weights or the grouped Gibbs blocks and then dropping
the corresponding refusal (or, for grouped, adding `variance` to rbart_vi's
formals); no new surface is needed for either, since the formal already
exists on dbarts()/bart2() and rbart_vi needs only the same formal plus the
collision check section 5 already applies to the latent families.

## 16. Post-landing: the saved variance trees ride the state (2026-08-18)

Section 5's "both replay from saved trees independently" was true of a LIVE
sampler only. `Chain::getState`/`setState` carried the variance forest's LIVE
trees and not its SAVED (keepTrees) buffer, so any state round trip -
`storeState()` + `saveRDS`/`readRDS`, a new session, or `$copy()` - left a
re-created sampler predicting from whatever `initializeSavedTrees` had put in
the slots. That fill was a leaf of 0.0, and `predictVarianceFromSavedSample`
forms s^2(x) as a PRODUCT over the variance trees, so every restored slot
reported s(x) == 0 with no error and no warning while yhat survived bitwise.

What now rides the state, in `ChainStateData` beside `varianceTrees` (the
variance forest sits outside `forests_`, so these are siblings rather than a
`ForestStateData` entry, and they carry no capacity scalar - the block size
against a construction-fixed tree count already implies it):

- `savedVarianceTrees`, slot-major capacity x numTrees, serialized as the
  append-only chain-level blocks `variance.saved.vars/.values/.sizes/.flags`.
- `varianceTreeMasks` / `savedVarianceTreeMasks`, the pooled categorical side
  channels for the live and saved trees, as `variance.masks` and
  `variance.saved.masks`. `Tree::flatten` requires a channel whenever the store
  holds a pooled column, and the variance forest was passing none: a variance
  tree splitting on a >63-level factor crashed inside `storeVarianceSavedTrees`
  under keepTrees and inside `Chain::getState` without it. The channel now runs
  end to end - flatten, rebuild, validation, warm-start install, and both
  replay paths.

Blocks were APPENDED, so neither `stateFormatVersion` nor
`minReadableStateFormatVersion` moved (the registry rule at their declaration).

Semantics:

- An ABSENT saved block against a live capacity > 0 is REFUSED. `getState`
  fills the block whenever keepTrees is on, so an empty one can only be a state
  written before the channel existed; accepting it would restore the
  destination's own identity fill and report a plausible constant s(x) where the
  defect it replaces at least reported an obvious zero. The same size gate
  refuses a capacity mismatch, as the mean side already does. `installForests`
  does not route through `stateIsValid`, so it carries its own gate on the same
  posture (below).
- The default fill is now the MULTIPLICATIVE identity 1.0. A 0.0 leaf is
  user-visible with no round trip at all - `predict()` over a slot no sweep has
  written yet multiplies a zero into the product - and, once validation applies
  the scale-leaf positivity law to saved trees, it would make the engine's own
  `storeState` emit a state its own `stateIsValid` rejects.
- Saved trees are held to form (`flatTreeIsWellFormed`, masks included) and leaf
  positivity, but NOT to the occupancy pass section 14 added for the live trees.
  A saved slot is a historical replay target routed over NEW rows, never over
  this sampler's partition; the mean side does not occupancy-check its saved
  trees either. A slot-sourced warm start is the one path that makes a saved
  slot LIVE, and it is occupancy-checked there, by `installVarianceForest` -
  so a donor that kept sweeps and then had its rows moved can hold a slot the
  destination refuses. Refusing is right: an unoccupied scale leaf reports a
  scale the data never supported.
- `setTreeStorage` re-runs `initializeSavedTrees` and resets the sample counter,
  so a capacity change followed by `$setState` refuses under the size gate -
  already the mean-side behavior.

Warm start pairs the two halves at ONE sample. `installForests` slices the
donor's saved variance buffer at the slot its mean forests come from, so
`samples = k` means sample k's mean forest AND sample k's scale surface; the two
buffers are index-aligned by construction (one slot base drives both, and both
index by the sample number), so the pairing is exact rather than nominal. A
live-sourced start (slot < 0, the pool a keepTrees = FALSE donor offers) takes
the donor's live pair, unchanged. This matters most where the argument is
`samples = NULL`: the chains spread across the donor pool for overdispersion,
which the old reassembly delivered on the mean half while handing every chain
the same final scale surface.

The buffer must cover the named sample, and its STRIDE must be the donor's own
variance tree count: a state whose live variance block was replaced by a shorter
one would otherwise slice across two sweeps' trees and still present the right
count downstream. Both are refused as `varianceSlotMismatch`, before any live
state is touched. The sliced trees are additionally held to the scale-leaf
positivity law state validation applies to a saved slot - the fix multiplies the
hand-buildable surface by capacity, and a rebuild scatters the leaf straight
into a divisor.
