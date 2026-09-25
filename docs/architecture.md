# dbarts architecture

A contributor orientation to the engine as it exists today: what runs where,
and which layer owns which decision. The proposals and landing records behind
it are in docs/design/, anchored by core-generalization.md.

## Layering

    R/                        bart (the front door; bart2 is its one-release
                               alias), bartBT (0.9-x's BayesTree-style door),
                               dbarts, xbart; the dbartsSampler reference
                               class (R/dbarts.R)
                               |
                               v  method calls
    R/bartcore.R               engine-specific method bodies the reference
                               class delegates to
                               |
                               v  .Call, C_-prefixed entry points
    src/R_interface_bartcore.cpp   the bridge: SEXP <-> engine types,
    src/R_interface_bartcore_common.hpp   argument validation, PROTECT
                               bookkeeping, error conversion
                               |
                               v  bartcore::createSampler(...) [one call, at creation]
    src/bartcore/facade.hpp    SamplerBase (virtual) / SamplerFacade<L>
                               the type-erased boundary
                               |
                               v  one vtable hop per call thereafter
    src/bartcore/{sampler,chain,combiner,model,moves,grow,scan,tree,data}.hpp
                               the header-only C++20 engine (bartcore.hpp
                               is the umbrella include)
                               |
                               v  calls into
    src/misc/*, src/external/*, src/rc/*
                               compiled support libraries: linear algebra,
                               partitioning and moment kernels with runtime
                               SIMD dispatch (misc.a); RNG and IO (external.a);
                               bounds-checked SEXP extraction (rc.a)

A second path down from compiled code: `inst/include/dbarts/dbarts.h` is the
only shipped header, a flat C API (`dbarts_sampler_*` entry points reached
through `R_GetCCallable`) that `LinkingTo: dbarts` consumers such as
stan4bart build against without seeing any C++ types. A consumer creates the
sampler in R and hands the C API the handle from its external pointer. The
implementation (`src/C_interface.cpp`) shares the bridge core with the
R-facing `.Call` entry points.

Above the reference class, sampler creation has one resolution step.
`bart` and `bartBT` build their samplers by calling `dbarts()`, and
[`resolveSamplerSpec`](../R/spec.R) turns `dbarts()`'s arguments - the
family object, which carries the residual prior
(`gaussian(sigma = chisq(3, 0.9))`), the tree and node prior objects, and the
control - into the `(control, model, data)` triple plus the resolved response
family that construction consumes. It is exported as `dbartsSpec()`
(docs/design/consumer-spec-surface.md). `xbart` builds its per-fold samplers
itself.

The sampler lives in C++ and the R object holds an external pointer to it, a
`bartcore::SamplerBase`. `src/bartcore/` is header-only: it compiles into
whichever translation unit includes it - the R bridge, the C API, `tests/cpp`
- and is not itself a library. `misc.a` is a library because its SIMD kernels
are compiled once per instruction set in separate translation units that
cannot move into a consumer's.

## Dispatch tiers

Four tiers, chosen by how much work each dispatch decision amortizes over:

| Granularity   | Mechanism                              | Example                              |
|---------------|-----------------------------------------|--------------------------------------|
| Per R call    | one virtual call through `SamplerBase`  | run, setPredictor, setOffset         |
| Per iteration | virtual calls on the response family and the forest combiner; plain struct calls on the tree and DART priors | latent refresh, sigma/k/DART draws, glue draws |
| Per node op   | template instantiation selected at compile time; switch over a closed set | move type (moves.hpp), kernel table lookup |
| Per observation | monomorphic loops and kernel calls, no dispatch | partition compare, suffstat accumulation |

The leaf model (`ConstantGaussianLeaf` and its constrained variant
`MonotoneConstantGaussianLeaf`, `LinearGaussianLeaf`, `GPGaussianLeaf`, and
the variance forest's `ConstantVarianceLeaf`, all in
`src/bartcore/model.hpp`) is a compile-time template parameter `L`, threaded
through `Sampler<L>` -> `Chain<L>` -> `Forest<L>` and the free functions in
`src/bartcore/moves.hpp` because it sits inside the per-observation and
per-node-op tiers: `accumulate` and `logIntegratedLikelihoodForNode` must
inline. The response family is chosen once per chain and called through a
handful of per-iteration virtual hops, so runtime polymorphism costs it
nothing.

## The facade

`src/bartcore/facade.hpp` defines `SamplerBase`, an abstract class exposing
every operation a host needs (run, the setters, tree and state serialization,
prediction, queries), and `SamplerFacade<L, ResidT = double>`, a `final`
class template that implements `SamplerBase` by forwarding each call to a
`Sampler<L, ResidT> impl_`. It exists because the leaf model is a
compile-time parameter, which forces type erasure: it is the one concrete
type the rest of the package holds a pointer to, so no caller of the engine
has to be templated on `L` itself.

`ResidT` is a second compile-time parameter, carrying the opt-in fp32
residual storage (`storage = "single"`;
docs/design/reduced-precision-storage.md). It mints exactly one extra
instantiation, the gaussian constant leaf with `ResidT = float`, alongside
the `double` path everything else uses.

Selection of `L` happens exactly once, in the free factory functions at the
bottom of facade.hpp (`createSampler` and its siblings), called from the
bridge at sampler creation. Every factory that can build a variance forest
asks one shared predicate, [`varianceForestIsRefused`](../src/bartcore/facade.hpp),
so they cannot drift apart: a variance forest is refused for any family other
than gaussian and aft, with designated leaf covariates, with Student-t
residuals, or with an active monotone constraint. Nothing re-dispatches on
`L` thereafter.

## Model concepts and their shipped implementations

`src/bartcore/model.hpp` defines two independent extension points, dispatched
by two different mechanisms.

**Leaf model** - compile-time, via C++20 concepts, each a named set of
requirements a type must satisfy to instantiate a template.
`LeafModelCore` is the base requirement, a closed-form
`logIntegratedLikelihoodForNode`. It refines into `ScalarLeafModel` (one
parameter per leaf: `ConstantGaussianLeaf`), `VectorLeafModel`
(`numParams()` doubles per leaf, fits evaluated per observation:
`LinearGaussianLeaf`), and `FunctionLeafModel` (one drawn value per member
observation: `GPGaussianLeaf`).
`IntegrableLeafModel` is the union of the three, and it is what the
templates in `chain.hpp` and `sampler.hpp` are constrained on; the free
functions in `moves.hpp` take the wider `MoveScorableLeafModel`, which also
admits `ScaleLeafModel` and so the variance forest's `ConstantVarianceLeaf`,
a leaf that scores moves without satisfying the mean-model concepts.

Every shipped leaf model requires a closed-form marginal likelihood, and one
kernel samples the tree structure: the conjugate Metropolis-Hastings moves in
`moves.hpp` (see "Tree moves"). `src/bartcore/grow.hpp` is a second tree
builder, XBART-style root-down construction, but it initializes a forest
rather than sampling from its posterior: stationarity belongs to the exact
moves that follow, so a grown forest need only be a legal chain state. It is
reached as `dbartsSampler$growFromRoot` (docs/design/grow-from-root.md).

**Response family** - runtime, via a virtual base class. `ResponseModel`
(model.hpp) owns the working response and weights the sweep reads, latent
refresh, sigma draws, and the response-side mutation entry points. The
concrete class is chosen with a `switch` on
`ResponseFamily` (`enum class ResponseFamily { gaussian, probit, logistic,
aft, ordinal, nbinom }`) inside `Chain`'s constructor
(`src/bartcore/chain.hpp`). Two paths sit off the enum: the gaussian arm
yields `TResponse` instead of `GaussianResponse` when Student-t residuals are
requested, and the multinomial model installs `MultinomialResponse` through
its own construction path. Every chain in a sampler shares the same family.

**Split-variable selection**: `CGMTreePrior` (model.hpp) owns the
depth-decaying growth probability and the split-variable log-probability -
uniform over available variables by default, or weighted by
`splitProbabilities` when DART is active, in which case `DartPrior`
Gibbs-updates those weights from each iteration's per-variable split counts.
Both are plain, non-virtual structs held by value inside `Forest<L>`: one
instance per forest, no dispatch, since a chain's split-selection policy
never changes after construction.

## Forests and combiners

A `Forest` (`src/bartcore/combiner.hpp`) is one ensemble of the backfitting
sampler: its trees, their fits, its running residual, its leaf model
instance, its split selector, and its own tree count, move probabilities and
`k` hyperprior. A chain holds a vector of forests - one for most models, two
for BCF, one per category for multinomial.

When a chain holds more than one forest it delegates their coupling to a
`ForestCombiner<L, ResidT>` (combiner.hpp), which answers three questions per
sweep: `formForestResponse` gives forest f the response and precisions its
own leaf draws see, the residual net of every other forest's scaled
contribution; `formForestVetoWeights` gives the precisions forest f's
empty-leaf veto reads; and `combinedFits` returns the per-observation
location all the forests together imply, which the response family's latent
and sigma draws consume. A single-forest chain carries no combiner.

`AmplitudeForestCombiner` is the multiplier family. Each forest carries a
`ForestBasis` - an n x q row-major matrix - and an amplitude block of q
coefficients; the forest's multiplier for observation i is the dot product of
the two, and the combined location is the sum over forests of multiplier
times forest fit. A forest whose multiplier is a plain scalar carries a dense
all-ones column, so there is exactly one multiplier path. Amplitudes are
drawn jointly per forest from a Gaussian full conditional (`drawGlue`),
optionally under a half-Cauchy scale mixture, and `setForestBasis` is the
only route that replaces a basis on a live sampler. BCF is the two-forest
instance: the prognostic forest takes a one-column basis and a single
amplitude, the treatment forest the (1 - z, z) indicator pair and the two
treatment scales. `MultinomialForestCombiner` instead couples K symmetric
category forests through a softmax likelihood, with a one-vs-rest
Polya-Gamma augmentation drawn against the current margins immediately
before each category's own forest updates (`drawForestGlue`).

A forest's cached fits (`totalFits`) are kept by difference updates, so they
may differ from the forest's leaves gathered in tree order by additive
rounding only. A transform that writes leaf values in bulk - the multinomial
level shift in `afterCombine`, the level-fibre shift - updates the cache in
place only when it is additive; a multiplicative one must re-derive the cache
from the leaves before the sweep ends, since a gap it multiplies compounds.
`Chain::run` checks the rule in debug builds.

Two per-observation channels ride alongside, only one of them per-forest.
`Chain::setForestWeights` installs a precision factor composed into forest f's
leaf conditionals alone (`composeForestWeights`), admitted only by a combiner
whose `supportsForestWeights` is true. `Chain::setActiveRows` is global: a 0/1
mask saying which rows are in the data set this sweep, admitted by the
response family (`supportsActiveRows`) rather than by the combiner. An
all-ones mask installs nothing and any element other than 0 or 1 is refused.
The mask is also where probit and ordinal fits put 0/1 case weights, since
neither family has a weight channel (`enforceWeightPolicy`,
docs/design/active-rows-mask.md).

Split availability is restricted per forest, not per chain: a forest may
carry a column mask (BCF moderators, a column-restricted variance forest) and
an `InteractionConstraint` (`src/bartcore/tree.hpp`), a maximum interaction
order plus a forbidden co-occurrence adjacency that every tree of that forest
borrows. Since a donor grown under a different restriction may hold a tree
the destination forbids, warm start refuses on containment before touching
live state (`Chain::interactionStateFeasible`, `columnMaskStateFeasible`).

Design notes: docs/design/forest-combiner.md for the hierarchy,
multiplier-combiner.md for the amplitude family, and
interaction-constraints.md for containment.

## Tree moves

`src/bartcore/moves.hpp` holds five structural moves: birth/death, change,
swap, perturb and rule_gibbs, the conjugate Metropolis-Hastings proposals and
their acceptance ratios, as free functions templated on
`MoveScorableLeafModel`. `metropolisJumpForTree` (a free function in
moves.hpp, called from `Chain`) is the per-iteration, per-tree entry: it
draws a step type (`StepType::birth/death/swap/change/perturb/ruleGibbs`) and
dispatches to the corresponding move function. Swap, perturb and rule_gibbs
all ship at probability 0. Swap - at production forest sizes it is nearly
all no-op - is the only move that rotates a child's rule up the tree, so a
single-tree fit wanting to cross between rootings sets it positive through
`dbartsControl(proposal.probs = )`. Perturb moves one interior node's
ordinal cut by a single grid position, keeping its split variable and the
tree's shape; it ships at zero pending a benefit measurement
(docs/design/perturb-move.md).
Rule_gibbs replaces the split rule at a nog node - an interior node whose two
children are both leaves - with an exact draw from that rule's own full
conditional over the available ordinal variables and their admissible cuts,
at acceptance one; it ships at zero pending a benefit measurement
(docs/design/nog-gibbs.md). All five structural probabilities exactly zero
freezes the structures: [`Chain::run`](../src/bartcore/chain.hpp)
reads [`structureIsFrozen`](../src/bartcore/moves.hpp) once per forest and
skips `metropolisJumpForTree` for every tree of the mean and variance
forests, so no move is proposed and none of its randomness is drawn, while
leaf values, sigma and the family's latents keep sampling - a fitted forest
re-sampled as a fixed basis.

Every candidate branch's empty-leaf veto is ranked
(`Tree::leafVetoRank`, [`resolveVetoRank`](../src/bartcore/moves.hpp)): rank 2 is
a leaf with no member at all, rank 1 a leaf whose members all carry zero
weight, rank 0 a leaf a likelihood term reaches. Comparing a (current,
proposal) pair, the worse-ranked branch takes `-HUGE_VAL` outright; when both
ranks are equal the comparison runs on the finite log-likelihoods as usual.
Only rank 2 is absolute - no move may install a leaf with no member, from any
state - so a chain sitting on a rank-1 branch still mixes under the prior and
transition kernel at constant likelihood rather than freezing. Member-empty
leaves stay out of the chain state entirely; a weight-emptied leaf is
penalized rather than forbidden (docs/design/empty-leaf-veto.md).

Rules themselves are typed by column: ordinal rules compare a code against a
threshold, categorical rules test a bit of a direction mask. Masks up to 63
categories live inline in the rule's 64-bit word; wider ones store an offset
into a per-tree mask pool (`Tree::maskPool`, compacted between moves once
garbage passes a high-water mark), machinery that turns on only when some
column has more than 63 categories (`ColumnStore::hasPooledCategorical`).
Both kinds share one word, so bit 63 doubles as the missing-value direction
flag for either.

## One sweep

`dbartsSampler$run` (R/dbarts.R) calls `bartcoreSamplerRun` (R/bartcore.R),
which `.Call`s `bartcore_run` (`src/R_interface_bartcore.cpp`). The bridge
calls `run` on the `SamplerBase` it holds; the facade forwards to
`Sampler::run` (sampler.hpp), which hands each chain to `Chain::run`
(chain.hpp). One iteration of `Chain::run` is:

1. Under a variance forest, form the mean weights `w_i / s^2(x_i)`
   (`formMeanWeights`); the global sigma stays fixed at 1.
2. When `levelGibbs` is `TRUE`, or `NA` (the default) and that forest's
   structural mixture is frozen, draw the constant-leaf forest's level shift
   (`drawLevelShift`): a constant added to every occupied leaf of a tree, the
   constants summing to zero across the forest's trees, leaving the fitted
   function unchanged. It runs ahead of every channel the rest of the sweep
   writes; the mixture is read once per forest here and reused by step 3's
   frozen skip, and a forest that skips consumes no generator draw.
3. For each forest in turn, and for each of its trees: roll the running
   residual so `treeY` holds the response net of every other tree's current
   fits, propose one move with `metropolisJumpForTree` and accept or reject
   it, then draw the tree's leaf values and write its fits
   (`sampleParametersAndSetFits`). A multi-forest chain asks the combiner for
   this forest's own response and precisions first.
4. Rebuild the forest's `totalFits` once the tree loop ends
   (`finalizeTotalFits`).
5. Refresh the response family's latents against the combined location
   (`ResponseModel::refreshLatents`) and, where sigma is a free parameter,
   draw it (`drawSigma`).
6. Draw the combiner's glue and its post-combine move; sweep the variance
   forest against the mean residual (`sweepVarianceForest`); draw each
   forest's `k` and, under DART, its split weights.
7. Record the sample if this iteration is a kept one.

Every draw there reads the chain's own generator. Two optional host hooks
run inside the loop, a per-sweep hook before each iteration and a per-draw
callback at each kept draw; [Callbacks and errors](#callbacks-and-errors)
covers both.

## ColumnStore

`src/bartcore/data.hpp` defines `ColumnStore`, the predictor container every
chain in a sampler shares; chains read it and never mutate it directly. Its
governing idea is that the engine owns quantized codes rather than raw
predictors: `std::vector<xint_t> codes` (`xint_t` is `std::uint16_t`) against
per-column cut points, with the store the sole owner of cut construction and
re-quantization. docs/design/data-store.md is the standing reference for the
cut grid, the code blocks, the source descriptors, the borrowed view's value
channels and the mutation transaction. Four facts bear on the rest of this
document:

- A column is `numeric`, `categorical` or `orderedFactor`, but rules, scans,
  masks and the flat replay branch on the derived predicate
  `splitsBySubset(j)` rather than on the kind.
- A factor column of either kind keeps no doubles on either side; what a leaf
  model reads instead is gathered into store-owned copies
  (`gatheredRawValues`).
- A CSC-built column at or below 20% nonzero density
  (`sparseDensityThreshold`) takes a rank-bitmap representation; a denser one
  densifies into `codes`. A dense build never takes rank-bitmap storage,
  however sparse its values.
- `hasMissing[j]` gates the extra missing-direction draw in rules; a reserved
  code marks a missing cell.

None of this is arena-allocated; every array is an ordinary `std::vector`
sized at build or mutation time.

## The mutation surface

The contract that lets BART sit inside a larger sampler - predictors,
response, offset and weights swappable between MCMC iterations - is
implemented on `bartcore::Sampler<L>` (`src/bartcore/sampler.hpp`) and
exposed through `SamplerBase`; docs/design/bart-as-a-component.md states which
mutations are legal when. Three kinds:

- **Response-side** (`setResponse`, `setOffset`, `setWeights`, `setSigma`):
  passed to the response family; no tree structure changes, and a length
  mismatch fails before any state does. `setResponse` and `setOffset` keep
  the response transform - the location and scale that map the response onto
  the tree prior's scale - as it was at creation unless `updateScale` is set.
  When it is set under a variance forest, the variance forest's prior and its
  current surface are restated in the new working units as well, so the chain
  matches one created on the new response.
- **Whole-data replacement** (`setData`): predictors, response and
  optionally weights, offset and test rows, with a possibly different number
  of observations. It first checks every factor cell of the replacement and
  of any test matrix against each column's fixed level table, and refuses the
  whole call, leaving the sampler untouched, if any cell is not a known
  level. Past that check it is not transactional: cut points are rebuilt,
  each existing split moves to the nearest new cut, and any subtree left
  empty or invalid collapses. It always re-derives the response transform and
  re-anchors a variance forest with it.
- **Predictor-side** (`setPredictor`, `updatePredictor`, and the
  per-observation session API): transactional, with three outcomes -
  `enum class PredictorUpdateResult { accepted, rolledBack,
  invalidCutPoints }`. With `forceUpdate` set, every chain force-refreshes
  its trees, collapsing any split that would empty a leaf into its parent
  with an effective-observation-weighted parameter merge, and the call always
  returns `accepted`. Without it, `revalidateAllChains()` validates every
  tree of every chain first and only rebuilds fits if all stay valid,
  restoring the snapshot and repartitioning every tree on failure. A
  quantile-mode cut refresh that would induce fewer cuts than an existing
  split needs is rejected before any mutation happens. The snapshot and
  rollback mechanics are [Predictor mutation transaction](design/data-store.md#predictor-mutation-transaction).
- **Per-observation updates** use a `PredictorUpdateSession`
  (`beginPredictorUpdate`/`updatePredictorPerObservation`): stage one
  observation's leaf moves against running per-leaf occupancy counts, test
  validity, then commit or skip, with one fits rebuild at the end rather
  than one per observation. `updatePredictorPerObservationJointly`
  (facade.hpp) sweeps several samplers sharing an index-aligned column in one
  randomized scan order, installing an observation everywhere or nowhere so
  their fits never diverge.

The predictor-side transactions are all-or-none across every tree of every
chain: validation runs over the whole sampler before any chain's fits are
rebuilt.

## Tree storage forms

Four representations, used for different purposes. `Tree` and `FlatNode` are
the two encodings; `SamplerStateData` and the wire format are aggregates
built on top of `FlatNode`.

- **Live**: `Tree` (`src/bartcore/tree.hpp`) is a flat arena -
  `std::vector<Node> nodes`, children allocated as adjacent pairs so
  `rightChild == leftChild + 1` always holds. A `Node` carries its `Rule`,
  its `[begin, end)` span into the tree's external observation-index buffer,
  and, for scalar and function leaves, the constant-leaf sufficient
  statistic. The moves and the sweep operate on this directly.
- **Flat**: `FlatNode` is one node of a pre-order-flattened tree - an
  ordinal cut point, a categorical direction mask (inline or pooled, tagged
  by `FlatKind`), or a leaf parameter - replayable against raw predictors
  without the `ColumnStore` that quantized them. It is the one format shared
  by saved-tree storage (`keepTrees`), external reporting (`getTrees`) and
  state serialization.
- **State**: `SamplerStateData` (sampler.hpp) is the whole sampler's
  in-process serializable state - one `ChainStateData` per chain, itself one
  `ForestStateData` per forest, plus the store's cut points and the
  saved-tree write cursors. `stateFormatVersion` is 1, as is
  `minReadableStateFormatVersion` (`src/R_interface_bartcore.cpp`); blocks
  are read by name and an absent optional block is defaulted, so adding one
  bumps neither number - an older reader ignores the name it does not know,
  and a newer reader defaults it when an older state omits it. Restore is
  semantic, not bitwise: a restored chain rebuilds partitions from tree
  structure and cut points and `totalFits` by summing tree fits, and
  restates the residual prior - a variance forest's included - on the
  restored response transform, so it continues equivalently but not
  bit-for-bit.
- **Wire**: what actually leaves the process. `storeState`
  (`src/R_interface_bartcore.cpp`) flattens `SamplerStateData` into a
  struct-of-arrays SEXP that `setState` reads back, one list per chain, whose
  node values ride a RAWSXP at 8 bytes per node, so an inline categorical
  mask's bit pattern survives verbatim.

## RNG architecture

Every chain owns its own Mersenne Twister (`ext_rng*`), created by
`createChainRngs` (`src/R_interface_bartcore.cpp`) at sampler construction,
a single-chain sampler included. Seeding is the only place a chain's stream
meets R's. With an explicit `control@seed`, a separate seed generator seeded
from that value hands each chain its seed, so a single-chain run with seed S
reproduces chain 0 of any multi-chain run with the same seed. Without one,
the chain seeds are drawn once from R's stream (`unif_rand()` inside a
`GetRNGstate()`/`PutRNGstate()` bracket), so a prior `set.seed()` determines
them.

After that, nothing a sampler does advances R's stream. Sweeps, prior draws,
grow-from-root, the probit latent redraw on `setResponse` and the
per-observation update's scan order all draw from the chains' generators (the
scan order from chain 0's). Several bridge entries still bracket those calls
with `GetRNGstate()`/`PutRNGstate()`; the brackets draw nothing. The one
bridge entry that does draw from R's stream is `dbartsDrawLatents`, which
draws a family's augmentation variables for a caller's fit rather than for a
sampler.

## Threading model

`Sampler::run` (sampler.hpp) runs chains on as many `std::thread` workers as
`numThreads` and the chain count allow. While running, a chain reads the
shared `ColumnStore` but never writes it, never calls R, and never touches
another chain, so results are bitwise identical at any thread count. Two execution paths across chains, and one pool
inside a chain:

- **Inline** (one worker or fewer): chains run sequentially on the calling
  thread; progress prints directly and the interrupt poll runs between
  sweeps, throttled to roughly every 100ms.
- **Worker threads** (more than one): each of W workers runs every W-th
  chain. Workers never call into R, so progress lines are queued
  (`QueuedProgressSink`) for the main thread to print, and cancellation is an
  atomic flag the main thread sets after polling for an interrupt. The main
  thread waits on a condition variable the last chain signals, so the call
  returns as soon as the chains finish; the wait wakes every 100ms only so
  the main thread can print progress and poll. On POSIX platforms `SIGINT` is blocked in worker threads, so a Ctrl-C
  reaches only the main thread, whose poll turns it into a cooperative cancel
  rather than an R longjmp across threads.
- **routeTestRows** (`chain.hpp`): the one pool that runs inside a chain
  during sampling rather than across chains - a `misc_mt` pool
  (`testFitPool_`) fanning a tree's test-row routing across this chain's
  share of the thread budget, serial below `testFitParallelCutoff` test rows.
  General within-chain parallelism is deliberately not used
  (docs/design/within-chain-threading.md); this pool is the one sanctioned
  exception, because it is RNG-free and output-disjoint, and so bitwise
  identical at any thread count.

**Prediction** (`Sampler::predictColumns`, sampler.hpp) mirrors `run`'s
worker-thread design: a per-call `n.threads` partitions the (chain, draw)
work list across a `std::thread` fan-out with no cross-thread reduction, so
replay is bitwise identical at every thread count, including the inline
below-cutoff path. `numThreads = 0` means the sampler's own thread count
(docs/design/threaded-predict.md).

## Callbacks and errors

**Per-sweep hook.** `SweepCallback`, passed to `Sampler::run` as `onSweep`,
fires on the calling thread before every sweep and lets the host change
conditioning state between sweeps without a round trip; which changes are
legal there is docs/design/bart-as-a-component.md's subject. It runs only
when chains run inline, and `Sampler::run` requires the caller not to set it
alongside worker-thread chains. The flat C API has no entry for it. The one
bridge entry that installs it, `bartcore_runWithCallback`, evaluates an R
closure through `R_tryEval` and refuses more than one chain. No function in
dbarts calls it; it is kept as the entry host packages build on when they
move work into per-sweep callbacks, for example to keep memory down.

**Per-draw callback.** A compiled function of the shipped
`dbarts_draw_callback` type, registered through
`dbarts_sampler_setDrawCallback` in the C API or passed from R as the
`callback` argument of `bart` and the sampler's `run`. It fires once per kept
draw, on the thread running that chain - a worker thread when chains run in
parallel - and sees everything recorded for that draw. A nonzero return
stops the run (docs/design/per-draw-callbacks.md).

**Errors.** The engine reports every error as a C++ exception. Each bridge
entry catches it, leaves the handler, and only then raises the R error, so
no R longjmp crosses an engine frame. A per-draw callback may itself raise an
R error only when its run is inline on the thread that entered it: there the
call is made under `R_UnwindProtect`, which turns the jump into a C++ unwind
through the engine's frames, and the entry resumes the jump once they are
gone. A C++ exception the callback throws is held until `R_UnwindProtect`'s
own frame has returned and then rethrown along the same path. On a worker
thread the callback must not raise; it returns nonzero instead. A run that
stops early - an error, an interrupt, a nonzero return - does not advance
the sample cursors, so the draws and saved trees it wrote are dropped, and
the sampler stays usable for a later run.

## Reproducibility contract

Within one host, the same seed gives bitwise-identical draws whichever SIMD
instruction set the runtime dispatch selects. That holds by construction:
every dispatched double-precision kernel is elementwise or a permutation
(partition, vector add and subtract, AXPY, transpose), and the draw-path
reductions - the per-node sufficient statistics of
[`misc_computeSufficientStatisticsFast`](../src/misc/moments.c) and the
residual sum of squares behind the sigma draw - are scalar with a fixed
order. The gate is
["C_dbarts_getMaxSIMDInstructionSet"](../inst/tinytest/test-simd.R),
which fits the same data at every forced dispatch level and compares at
tolerance zero. It compares dispatch levels of one binary, so a header
change that moves every level together is caught by the seed-locked
snapshot tests and the equivalence baselines instead. tests/cpp sets
dispatch once, at the host maximum.

A second build mode, `--enable-reference-build` at configure time, is
reported at runtime by [`buildInfo`](../R/buildInfo.R) (mode, compiled
instruction sets, dispatch level). It exists to keep a scalar, fixed-order
draw path should the shipped build ever vectorize those reductions; none
does, so today the two builds compile the same kernels. The seed-locked
snapshot tests and the recorded equivalence baselines are keyed to the
reference build all the same, and exit on the shipped one; [CI](plans/README.md#ci)
says where they run.

Across hosts the guarantee is never bitwise, for two reasons outside the
engine: the equivalence scenarios generate their data through the platform
libm (the Friedman function calls `sin`), so the inputs already differ in
the last bit between macOS and Linux, and compilers differ in how they
generate transcendental and floating-point code. The equivalence baselines
are therefore recorded and compared bitwise on one machine. On another
host the main corpus is compared statistically only. The BCF and
multinomial harnesses' `--cross-host` mode, run on the reference build,
requires the draws to match within a relative deviation of 1e-8
([`crossHostRtol`](../benchmarks/R/bcf-equivalence.R)); a scenario outside
that bound still passes if a weak statistical comparison cannot tell the
runs apart. A cross-host bitwise mismatch is not a regression.

For kernel work this means: vectorizing or FMA-contracting a draw-path
reduction reorders its sum, which breaks the within-host property and turns
the gates into a statistical comparison plus a re-record. Measured, it buys
almost nothing - vector sufficient-statistics kernels gained under one
percent and FMA nothing on the load-bound hot loops - so the performance
lever is data layout, not SIMD reductions.

## Further reading

- docs/design/INDEX.md - every design note, with its status.
- docs/design/core-generalization.md - the proposal behind the engine.
- docs/design/data-store.md - the standing reference for the predictor data
  layer and its mutation transaction.
- docs/design/forest-combiner.md, multiplier-combiner.md - the multi-forest
  coupling and the amplitude family.
- docs/design/bart-as-a-component.md - what a host may mutate between sweeps.
- docs/design/interaction-constraints.md - the containment predicates.
- docs/design/empty-leaf-veto.md - the ranked veto and its stationarity
  argument.
- docs/design/feature-matrix.md - what each response family and extension
  supports today.
- docs/design/kernel-vocabulary.md - the compiled-kernel contract (`misc.a`)
  the engine dispatches into.
- docs/design/public-surface.md - what the R and C surfaces expose;
  consumer-spec-surface.md - the exported `dbartsSpec()` resolution surface.
- docs/design/pooled-masks.md, sparse-columns.md, mia-missingness.md,
  linear-leaves.md, gp-leaves.md, bcf.md,
  multinomial.md, monotone.md, heteroscedastic.md, grow-from-root.md,
  reduced-precision-storage.md - design and landing notes for each extension
  mentioned above.
