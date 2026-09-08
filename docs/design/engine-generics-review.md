# What the engine's generic axes should be

Status: REVIEW (analysis only; nothing proposed for build).

A memo on the shape of the C++ core. It asks what must vary independently in a
BART engine that covers the models this package fits or has designed for, works
that out from the models rather than from the code, and only then reads the code
against it. Section 7 argues against the answer.

## 1. The models the engine is meant to cover

Shipped. Continuous responses with Gaussian errors, and the same with Student-t
errors by scale mixture. Binary responses through a probit latent and through a
logistic Polya-Gamma augmentation, the latter taking trial counts as weights.
Ordered categorical responses by cumulative probit with sampled cutpoints.
Counts by negative binomial with an integer dispersion. Log-normal accelerated
failure time with right censoring, and discrete-time hazard as a person-period
recoding of the binary families. Multinomial responses as one forest per
category coupled by a softmax. Semicontinuous two-part responses composed in R
from two ordinary fits. Alongside those: a two-forest causal specification and
its generalization, K forests each carrying a per-row basis whose amplitudes are
drawn jointly; a variance ensemble modelling the residual scale as a second
product of trees; leaf models that hold a constant, a monotone-constrained
constant, a small linear regression, or a Gaussian process over designated
covariates; a sparsity-inducing prior on split variables; per-forest interaction
and column restrictions.

Designed but unbuilt, or plausibly next. Real-valued negative-binomial
dispersion. Exact serially-correlated or spatially-correlated errors, which need
a banded rather than diagonal precision. Multivariate responses, which today are
composed outside the engine as one sampler per outcome exchanging offsets and
scales. Varying-coefficient models, which the multiplier family expresses except
that its forests must carry constant leaves. Shape constraints beyond
monotonicity. Leaf models that draw their own hyperparameters each sweep -
sampled kernel lengthscales, random feature bases, global-local shrinkage on
leaf scales - which is where three current peer packages differ from this one.
Subset splits on unordered factors carrying a network structure. Per-forest row
restriction, which a mixture-cure specification wants.

Two things are out on model grounds and should stay out. Soft trees replace the
partition with probabilistic routing, so every row contributes to every leaf and
per-leaf sums stop being the sufficient statistics; nothing above the tree
representation survives that. A non-diagonal error precision couples the leaves
of a tree, so the per-leaf independent draw and the per-leaf marginal both go.
Those are different engines, not missing axes.

## 2. What varies from one model to the next

Asked of every model above, six questions have distinct answers: what a leaf
holds and how it is drawn; what one row contributes to the likelihood; what
augmentation the sweep carries and redraws; how many forests there are and how
their outputs combine; what global quantities the chain draws each sweep and
must carry in its saved state; and what a structural proposal needs to score
itself.

Reading the answers down the list, one fact dominates. Every shipped model
reaches the tree-sampling machinery through exactly one interface: a per-row
working response and a per-row precision. The family is what turns the observed
data and the current fit into that pair - a latent draw for the binary and
censored families, a Polya-Gamma draw for the logistic and count families, the
identity for the Gaussian one. The tree machinery then reads that pair only
through sums over a node's members. Nine response families, four leaf models and
three couplings all ride that one interface, which is the engine's genuine
central abstraction and is worth stating plainly because everything else in this
memo is a footnote to it.

Three things do not fit through it, and they are precisely where the trouble is.
The monotone leaf's score reads its neighbours' current values, so its leaves are
not independent given the tree. The variance ensemble's residual is
multiplicative rather than additive, and it enters the mean model through the
precision rather than through the location. And the per-row multiplier family
scales each forest's output by a row-specific number before the sum. All three
are statements about how a forest's output reaches the likelihood, not about
what a leaf holds or what the family does with a location.

## 3. The axes the model space implies

Six, with the granularity at which each is decided.

**A. The leaf parameter and its marginal.** Decided per node. Varies in the
parameter's shape - a scalar, a fixed-length vector, one value per member row, a
positive scale - in which statistic of the node the marginal reads, and in
whether leaves within a tree are independent.

**B. The channel by which a forest's output reaches the likelihood.** Decided
per forest, applied per row. Three values occur in the shipped space: added to a
location; multiplied by a row-specific scalar and then added; multiplied into
the row's precision. A fourth is designed for: a per-forest row restriction.

**C. The observation model.** Decided per row, refreshed per sweep. Consumes the
combined location or locations and produces the working response and precision,
owning whatever augmentation that takes.

**D. The global blocks.** Decided per sweep, carried in state. The residual
scale, each forest's leaf-prior scale, the sparsity weights, the ordinal
cutpoints, the dispersion, the residual degrees of freedom, the amplitudes.
Every one of these is the same kind of thing: a small named block a family or a
coupling draws each sweep, serializes, and reports.

**E. The split-rule vocabulary and the availability law.** Per node, per column.

**F. The structural proposal.** Per node.

Independence. A, B and C are close to a full product in the model space and are
never treated as one. A varying-coefficient model is B crossed with A: per-forest
multipliers over linear leaves. A monotone treatment forest is B crossed with A
again. Heteroscedastic survival is B crossed with C. The refusals that block
these are, on inspection, refusals of bookkeeping rather than of models: the
reason a censored or ordinal response cannot take the multiplier family is that
nobody has threaded that family's global block through the coupling's draw
order, not that the conditional does not exist.

Cost tiers, which is what decides compile-time against run-time. Only one tier
is per-observation: the residual roll, the partition compare, the sufficient-
statistic accumulation and the fit scatter, which run for every row of every tree
of every sweep. The per-node tier - scoring a proposal, drawing a leaf - is a few
thousand calls a sweep, each amortizing a loop over the node's members. The
per-forest and per-sweep tiers are free at any dispatch. So the rule the engine
states for itself is right: compile-time only where per-observation work or the
sufficient-statistic type is at stake.

What that rule implies, though, is narrower than what was built. It says the
leaf's own accumulation and scatter loops must be monomorphic. It does not say
the whole sampler must be instantiated once per leaf kind. What forces that is
not the arithmetic but the storage: a constant leaf keeps per-tree leaf values
plus a row-to-leaf map, a vector leaf keeps a dense per-tree fit slab, a function
leaf keeps the fits themselves as the parameters. Three chain-level data
structures, so the chain is templated. Had the leaf owned its own storage behind
the same interface it owns its math behind, the leaf would be a per-node runtime
choice and the sampler would be instantiated once.

## 4. How the current shape lines up

What matches. The working-response-and-precision backbone is exactly the right
abstraction and is what lets nine families ride four leaf models with no
per-family tree code. Choosing the family at run time is right; it is decided
once per chain and touched a handful of times per sweep. The forest as the
composable unit - its own tree count, move probabilities, tree prior, leaf-prior
scale, column restriction and interaction constraint - is right, and is already
heterogeneous across the forests of one chain in every respect but one. The
coupling object is the right shape for axis B: it is asked, per forest per sweep,
for that forest's response and precisions and for the combined location. And the
interface a proposal needs from a leaf turns out to be two numbers - a veto rank
and a log marginal - which is as small as it could be.

What does not match, in descending order of consequence.

*The leaf kind is a property of the sampler rather than of the forest.* This is
the single misplacement everything else follows from. Because a chain's forests
must share one leaf type, the variance ensemble cannot be one of them: it is a
hand-written second forest kind carrying its own trees, index buffer, move
probabilities, tree prior, saved-tree buffer, categorical mask pool and scratch -
a near-duplicate of the ordinary forest, existing because the type system forbids
the obvious spelling. For the same reason every coupling is constant-leaf only,
so varying coefficients with linear leaves, a monotone treatment forest, a
smooth prognostic forest are all unreachable, none of them for a model reason.

*The channel is not an axis.* Additive is implicit in the forest; the multiplier
lives inside one coupling subclass; multiplicative-into-precision lives in a
bespoke chain member swept at a fixed point after the family's own draws. So
"two forests, one into the mean and one into the variance" is a special case
rather than an instance of anything, and the row-restriction door has nowhere to
land.

*The global blocks have no shape.* Each is a bespoke pair of accessors on the
response base class, plus a flag on the sampler's shape record, plus a named
state slot, plus a bridge channel, plus an R accessor. There are three of them
now and they are three transcriptions of one pattern; a fourth costs what the
third did.

*The type-erasure boundary is above everything.* Because the leaf kind is a
compile-time parameter, the one concrete class every host holds a pointer to must
project every capability of every family and every coupling onto itself: some
sixty flat virtual functions, a third of them meaningful for exactly one model.
That bloat is the direct price of putting the compile-time axis at the top of the
stack instead of the bottom, and it is also the cheapest thing to fix
independently, since shrinking that class to a common core plus a queried
capability object touches no sampling code at all.

*A leaf model cannot draw its own hyperparameters.* The leaf prior is one scale
per forest with a chi hyperprior on it. A leaf that wants a per-sweep draw of its
own - kernel lengthscales, a random feature map, a global-local shrinkage law
over leaf scales - has no place to make it. This is the axis on which several
recent peer packages differ, and the current design has no seam for it at all.

The operational measure of how far the leaf axis has leaked out of the leaf: the
chain reads leaf-shape traits at compile time in over a hundred places, and the
shared library carries four full copies of the sampler stack plus a fifth for the
reduced-precision residual variant.

## 5. The price of the sixth leaf model and the tenth family

A sixth leaf model today costs a concept and a struct - the honest part - and
then the chain's three storage shapes and their hundred-odd compile-time
branches; the flatten, saved-tree, state and prediction channels; a factory arm;
a fifth instantiation of a five-thousand-line chain, adding on the order of a
hundred exported symbols; an entry in an enumeration published on the C header;
and an R prior constructor. Under the derived shape it costs the struct, its own
storage and its own format serializer behind the same interface, and one registry
entry: no new instantiation, no new branch in the chain.

A tenth family today costs a response class, an enumerator, a construction arm,
and - for each global block it carries - the five-place pattern above, which for
the most recent such block reaches roughly nineteen files. Under the derived
shape it costs the response class plus a declaration of the blocks it carries,
which state, reporting and the R accessors read generically.

Both figures should be read with the counter-cost in view: the derived shape pays
a per-node virtual call, gives up the constant-leaf whole-sweep specialization
unless that is re-expressed as a specialization behind the interface, and is a
rewrite of a working engine that is gated bitwise against recorded baselines.

## 6. What the engine actually needs from its host

Less than one would guess, and it is worth recording because it bounds any
future portability question. The engine calls into R's mathematical library in
four places only: the normal density and distribution function and the t density,
used by the per-row log-likelihood channel and by the censoring and ordinal
calculations. The chi-squared quantile that calibrates the residual prior is
computed at the boundary, not in the engine. Progress output goes through a sink
the host installs; cancellation is a flag the host sets after polling; random
numbers come from the package's own generator, one per chain, which never touches
R's stream after seeding. The engine reports every refusal as a null or a false
return rather than an error, with exactly two exceptions - two rejection samplers
that give up after a bounded number of attempts, both reachable only from
main-thread entry points.

So a host-hook version is a vtable of four scalar functions, an error callback
that returns instead of jumping, and the two sinks that already exist. The
residual value of doing it is not portability for its own sake; it is that the
engine could then be exercised, and its component tests run, with no R at all.
The caveat is that the response transform and the prior calibration still speak
in terms the host derives, so this is a hook set, not a clean severance.

## 7. Three places a reasonable designer would land elsewhere

**The leaf may deserve to be compile-time after all.** The strongest argument for
the current shape is a measured one: the constant leaf's fused residual roll and
sufficient-statistic pass is a whole-sweep, per-observation transformation that
is only expressible when the leaf kind and its storage are known statically, and
it is worth roughly half again on the shipped hot path. If the plain Gaussian
constant-leaf fit is the case that matters - and by usage it plainly is - then
five instantiations to protect it is a defensible trade, and the derived shape
has to re-earn that path through a specialization hook whose adequacy this memo
asserts rather than demonstrates. A designer who weights the shipped hot path
above the unbuilt model space chooses what exists.

**The variance ensemble may not be a channel.** The derivation calls it "a forest
whose output multiplies the precision". But its running residual is
multiplicative, its leaf prior is a different conjugate family, its calibration
is a product across trees, and it is swept at a different point of the iteration
for a real reason - it conditions on a settled mean. A designer could hold that
"channel" is a false generalization over two objects sharing only the word
forest, and that the honest count of forest kinds is two, written out. The reply
is that the coupling machinery already exists in general form for the multiplier
family and only the shared-leaf-type restriction keeps the variance ensemble
outside it - but that cannot be shown without building it, which is exactly the
kind of claim this memo should not be trusted on.

**The product may not be worth reaching.** The derivation leans on combinatorics:
leaf times channel times family is nearly a full product in the model space, so
the axes should be independent. Demand is not combinatorial. Nobody has asked for
a Gaussian-process treatment forest, a monotone ordinal fit, or a
negative-binomial multiplier model; the cells with named consumers are a handful,
and each could be opened by hand more cheaply than the general mechanism costs. A
designer could reasonably keep the current axes, shrink the refusal matrix by
declaring the empty cells out of scope by decision rather than by accident, and
spend the saved effort on the two axes with named consumers: leaf-level
hyperparameter draws, and a leaf-prior law richer than one scale per forest.

## Appendix A. Where the claims can be checked

Leaf models and their concept hierarchy, including the parameter shapes of
section 3's axis A and the two optional coupling seams:
[`LeafModelCore`, `ScalarLeafModel`, `VectorLeafModel`, `FunctionLeafModel`, `ScaleLeafModel`, `MoveScorableLeafModel`, `TreeDrawLeafModel`, `ParamScoringLeafModel`, `ConstrainedLeafModel`](../../src/bartcore/model.hpp).
The four shipped location leaves and the scale leaf:
[`ConstantGaussianLeaf`, `MonotoneConstantGaussianLeaf`, `LinearGaussianLeaf`, `GPGaussianLeaf`, `ConstantVarianceLeaf`](../../src/bartcore/model.hpp).

The response family as a runtime interface, and the working-response-and-precision
backbone of section 2:
[`ResponseModel`, `ResponseFamily`, `workingResponse`, `workingWeights`, `refreshLatents`](../../src/bartcore/model.hpp).
The three global blocks of section 4 and their per-family accessors:
[`carriesOrdinalThresholds`, `carriesDispersion`, `carriesResidualDf`](../../src/bartcore/model.hpp).

The forest as the composable unit, and the coupling object that answers per-forest
response, precisions and combined location:
[`Forest`, `ForestCombiner`, `formForestResponse`, `combinedFits`](../../src/bartcore/combiner.hpp).
The two shipped couplings and the per-forest basis:
[`AmplitudeForestCombiner`, `MultinomialForestCombiner`, `ForestBasis`](../../src/bartcore/combiner.hpp).

The variance ensemble as a hand-written second forest kind, its multiplicative
roll, and its fixed sweep point:
[`VarianceForest`, `sweepVarianceForest`](../../src/bartcore/chain.hpp).
The chain templated on the leaf, and the constant-leaf whole-sweep specialization
of section 7:
[`Chain`, `rollAndSetNodeAveragesFused`](../../src/bartcore/chain.hpp).

The erasure boundary, its capability record, the per-leaf instantiation set, and
the composition refusals of section 3:
[`SamplerBase`, `SamplerFacade`, `SamplerShape`, `createSampler`, `createAmplitudeSampler`, `createMultinomialSampler`, `varianceForestIsRefused`](../../src/bartcore/facade.hpp).

What a proposal asks of a leaf - a veto rank and a log marginal:
[`metropolisJumpForTree`, `logLikelihoodForBranch`, `BranchScore`, `resolveVetoRank`](../../src/bartcore/moves.hpp),
[`leafVetoRank`, `Rule`](../../src/bartcore/tree.hpp),
[`ColumnStore`](../../src/bartcore/data.hpp).

The host dependencies of section 6: the four mathematical functions
([`Rf_pnorm5`, `Rf_dnorm4`, `Rf_dt`](../../src/bartcore/model.hpp)), the two
terminating rejection samplers ([`ext_throwError`](../../src/bartcore/chain.hpp)),
and the output sink ([`ext_printf`](../../src/bartcore/sampler.hpp),
[`ProgressSink`](../../src/bartcore/chain.hpp)).

The user-facing vocabulary the model space of section 1 is stated in:
[`dbartsPriors`](../../R/model.R), [`resolveSamplerSpec`](../../R/spec.R).
