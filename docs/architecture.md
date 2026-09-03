# dbarts architecture

A contributor orientation to the engine as it exists today: what runs where,
and which layer owns which decision. The proposals and landing records behind
it are in docs/design/, anchored by core-generalization.md.

## Layering

    R/                        bart, bart2, dbarts, xbart, rbart_vi;
                               the dbartsSampler reference class (R/dbarts.R)
                               |
                               v  method calls
    R/bartcore.R               engine-specific method bodies the reference
                               class delegates to
                               |
                               v  .Call, C_-prefixed entry points
    src/R_interface_bartcore.cpp   the bridge: SEXP <-> engine types,
    src/R_interface_bartcore_common.hpp   argument validation, PROTECT bookkeeping
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

A second, independent path down from the C ABI: `inst/include/dbarts/dbarts.h`
is the only shipped header, a flat C API (`dbarts_sampler_*` entry points
reached through `R_GetCCallable`) that `LinkingTo: dbarts` consumers such as
stan4bart build against without seeing any C++ types. Its implementation
(`src/C_interface.cpp`) shares the same bridge core as the R-facing `.Call`
entry points.

Above the reference class, sampler creation has one resolution step:
`resolveSamplerSpec` (`R/spec.R`) turns the user-facing arguments into the
`(control, model, data)` triple plus the resolved response family that
construction consumes. `dbarts()` delegates to it, and it is also exported as
`dbartsSpec()` (docs/design/consumer-spec-surface.md).

Two conventions are worth stating. The sampler is created and held C-side via
an external pointer on the R5 object, targeting a `bartcore::SamplerBase`.
And `src/bartcore/` is header-only, compiling into whichever translation unit
includes it - the R bridge, the C ABI, `tests/cpp` - so it is not itself a
compiled library. `misc.a` is, because its SIMD kernels are compiled per
instruction set in separate translation units that cannot move into a
consumer's.

## Dispatch tiers

Four tiers, chosen by how much work each dispatch decision amortizes over:

| Granularity   | Mechanism                              | Example                              |
|---------------|-----------------------------------------|--------------------------------------|
| Per R call    | one virtual call through `SamplerBase`  | run, setPredictor, setOffset         |
| Per iteration | virtual calls on the response family; plain struct calls on the tree and DART priors | latent refresh, sigma/k/DART draws |
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
asks one shared predicate, `varianceForestIsRefused`, so the two cannot
drift: the request is refused under a non-gaussian family, under designated
leaf covariates, under a finite `residualDf`, or under an active monotone
constraint. Nothing re-dispatches on `L` thereafter.

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
its own construction path. When `options.numGroups > 0` the chosen family is
wrapped in `GroupedResponse`, a decorator that Gibbs-samples per-group
intercepts into the offset between tree sweeps and forwards everything else
to the wrapped model. Every chain in a sampler shares the same family.

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

Two per-forest channels ride alongside: `Chain::setForestWeights` installs a
per-observation precision factor composed into forest f's leaf conditionals
alone (`composeForestWeights`), admitted only by a combiner whose
`supportsForestWeights` is true, and `Chain::setActiveRows` installs a global
0/1 mask saying which rows are in the data set this sweep, validated and
normalized in one pass, with an all-ones mask installing nothing.

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

`src/bartcore/moves.hpp` holds birth/death, change, and swap: the conjugate
Metropolis-Hastings proposals and their acceptance ratios, as free functions
templated on `MoveScorableLeafModel`. `metropolisJumpForTree` (a free
function in moves.hpp, called from `Chain`) is the per-iteration, per-tree
entry: it draws a step type (`StepType::birth/death/swap/change`) and
dispatches to the corresponding move function.

Every candidate branch's empty-leaf veto is ranked rather than a flat
`-HUGE_VAL` (`Tree::leafVetoRank`, [[moves.hpp#resolveVetoRank]]): rank 2 is
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
2. For each forest in turn, and for each of its trees: roll the running
   residual so `treeY` holds the response net of every other tree's current
   fits, propose one move with `metropolisJumpForTree` and accept or reject
   it, then draw the tree's leaf values and write its fits
   (`sampleParametersAndSetFits`). A multi-forest chain asks the combiner for
   this forest's own response and precisions first.
3. Rebuild the forest's `totalFits` once the tree loop ends
   (`finalizeTotalFits`).
4. Refresh the response family's latents against the combined location
   (`ResponseModel::refreshLatents`) and draw sigma (`drawSigma`).
5. Draw the combiner's glue and its post-combine move; sweep the variance
   forest against the mean residual (`sweepVarianceForest`); draw each
   forest's `k` and, under DART, its split weights.
6. Record the sample if this iteration is a kept one.

Every draw there reads the chain's own generator; nothing in the loop touches
R.

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

The embeddable-in-a-larger-sampler contract (predictors, response, offset and
weights swappable between MCMC iterations) is implemented on
`bartcore::Sampler<L>` (`src/bartcore/sampler.hpp`) and exposed through
`SamplerBase`; docs/design/bart-as-a-component.md states which mutations are
legal when. Two families:

- **Response-side** (`setResponse`, `setOffset`, `setWeights`, `setSigma`,
  `setData`): forwarded straight to the response family; no tree structure is
  touched, so these either apply or, for length mismatches, fail before
  touching state. `setData` also replaces the predictors, and so answers a
  status: it checks every factor cell of the replacement and of any test
  matrix against each column's fixed level table before anything moves, and
  refuses the whole call rather than ingesting one side and dropping the other.
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
  rollback mechanics are [[data-store.md#Predictor mutation transaction]].
- **Per-observation updates** use a `PredictorUpdateSession`
  (`beginPredictorUpdate`/`updatePredictorPerObservation`): stage one
  observation's leaf moves against running per-leaf occupancy counts, test
  validity, then commit or skip, with one fits rebuild at the end rather
  than one per observation. `updatePredictorPerObservationJointly`
  (facade.hpp) sweeps several samplers sharing an index-aligned column in one
  randomized scan order, installing an observation everywhere or nowhere so
  their fits never diverge.

Data transactions are all-or-none across every tree of every chain:
validation runs over the whole sampler before any chain's fits are rebuilt.

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
  saved-tree write cursors. `stateFormatVersion` is 3, as is
  `minReadableStateFormatVersion` (`src/R_interface_bartcore.cpp`); blocks
  are read by name and defaulted when absent, so an additive block loads into
  an older reader. Restore is semantic, not bitwise: a restored chain rebuilds
  partitions from tree structure and cut points, `totalFits` by summing tree
  fits, and the variance-prior anchor by re-running the same transform
  construction does, so it continues equivalently but not bit-for-bit.
- **Wire**: what actually leaves the process. `storeState`
  (`src/R_interface_bartcore.cpp`) flattens `SamplerStateData` into a
  struct-of-arrays SEXP that `setState` reads back, one list per chain, whose
  node values ride a RAWSXP at 8 bytes per node - so an inline categorical
  mask's bit pattern survives rather than being normalized by a REALSXP.

## RNG architecture

Every chain owns a dedicated Mersenne Twister (`ext_rng*`), created by
`createChainRngs` (`src/R_interface_bartcore.cpp`) at sampler construction -
including a single-chain sampler; there is no path where a chain samples
directly off R's generator. Seeding is the only place the two streams meet:
with an explicit `control@seed`, a dedicated seed-generator (itself a
fresh Mersenne Twister seeded from the control value) hands each chain a
seed, so a single-chain run with seed S reproduces chain 0 of any
multi-chain run with the same seed. Without an explicit seed, chain seeds
are drawn from R's stream once via `unif_rand()` inside a
`GetRNGstate()`/`PutRNGstate()` bracket, so a prior `set.seed()` in R
determines them.

After that one-time seeding, sampling itself never advances R's stream -
`Chain::run` and everything it calls draws exclusively from the chain's own
`ext_rng*` - and the bridge's remaining `GetRNGstate()`/`PutRNGstate()`
brackets cover only draws made outside a chain's sweep loop: a probit latent
redraw, the per-observation session's scan-order permutations, and the
sample-from-prior, grow-from-root and draw-latents entry points.

## Threading model

`Sampler::run` (sampler.hpp) runs chains on as many `std::thread` workers as
`numThreads` and the chain count allow; chains never touch the `ColumnStore`,
the R API, or each other while running, so results are bitwise identical
regardless of thread count. Two execution paths:

- **Inline** (one worker or fewer): chains run sequentially on the calling
  thread; progress prints directly and the interrupt poll runs between
  sweeps, throttled to roughly every 100ms.
- **Worker threads** (more than one): each worker claims chains in a
  round-robin. Workers must never call into R, so progress lines are queued
  (`QueuedProgressSink`) for the main thread to flush and cancellation is a
  relaxed atomic flag the main thread sets after polling. On POSIX platforms
  `SIGINT` is blocked in worker threads, so a Ctrl-C reaches only the main
  thread, whose poll turns it into a cooperative cancel rather than an R
  longjmp across threads.
- **routeTestRows** (`chain.hpp`): the one pool that runs inside a chain
  during sampling rather than across chains - a `misc_mt` pool
  (`testFitPool_`) fanning a tree's test-row routing across this chain's
  share of the thread budget, serial below `testFitParallelCutoff` test rows.
  General within-chain parallelism is deliberately not used
  (docs/design/within-chain-threading.md); this pool is the one sanctioned
  exception, because it is RNG-free and output-disjoint, and so bitwise
  identical at any thread count.

**Callback restriction**: the per-sweep conditioning callback (`SweepCallback`
in the engine, `dbarts_sampler_callback` on the C ABI) fires on the calling
thread before every sweep and lets the host mutate conditioning state between
sweeps without a round trip. It is usable only when chains run inline:
`Sampler::run` requires the caller not to set `onSweep` alongside
worker-thread chains, and `dbarts_sampler_setCallback` raises an error when
`numThreads > 1 && numChains > 1`. Which mutations are legal there is
docs/design/bart-as-a-component.md's subject.

**Prediction** (`Sampler::predictColumns`, sampler.hpp) mirrors `run`'s
worker-thread design: a per-call `n.threads` partitions the (chain, draw)
work list across a `std::thread` fan-out with no cross-thread reduction, so
replay is bitwise identical at every thread count, including the inline
below-cutoff path. `numThreads = 0` floors to the sampler's own thread count
(docs/design/threaded-predict.md).

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
  linear-leaves.md, gp-leaves.md, grouped-random-effects.md, bcf.md,
  multinomial.md, monotone.md, heteroscedastic.md, grow-from-root.md,
  reduced-precision-storage.md - design and landing notes for each extension
  mentioned above.
