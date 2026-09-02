# dbarts architecture

This is a contributor orientation to the engine as it exists today: what
runs where, which layer owns which decision, and what has to pass before a
change lands. It does not narrate how the engine got here; that history,
including abandoned alternatives and the landing record, is in docs/design/
(core-generalization.md is the anchor document).

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
reached through `R_GetCCallable`) that `LinkingTo: dbarts` consumers (e.g.
stan4bart) build against without seeing any C++ types. Its implementation
(`src/C_interface.cpp`) shares the same bridge core as the R-facing
`.Call` entry points (`R_interface_bartcore_common.hpp`).

Above the reference class, sampler creation has one resolution step:
`resolveSamplerSpec` (`R/spec.R`) turns the user-facing arguments into the
`(control, model, data)` triple plus resolved family token that
construction consumes. `dbarts()` delegates to it, and it is exported as
`dbartsSpec()` so a `LinkingTo` consumer holding its sampler C-side can
produce the same triple through supported surface
(docs/design/consumer-spec-surface.md).

Non-obvious conventions: R reaches the bridge entries as `C_dbarts_*`
(`C_rbart_fitted` is the one exception, `useDynLib(dbarts, .fixes = "C_")`
in NAMESPACE); `dbartsSampler`
is created and held C-side via an
external pointer stored on the R5 object, targeting a
`bartcore::SamplerBase`; `src/bartcore/` is header-only and compiles into
whichever translation unit includes it - the R bridge
(`src/R_interface_bartcore.cpp`), the C ABI (`src/C_interface.cpp`), and
`tests/cpp` - so it is not itself a compiled library - unlike `misc.a`,
which stays a compiled library because its SIMD kernels are compiled
per-ISA in separate translation units that cannot move into a consumer's
TU.

## Dispatch tiers

Four tiers, chosen by how much work each dispatch decision amortizes over:

| Granularity   | Mechanism                              | Example                              |
|---------------|-----------------------------------------|--------------------------------------|
| Per R call    | one virtual call through `SamplerBase`  | run, setPredictor, setOffset         |
| Per iteration | virtual calls on `ResponseModel`; plain struct calls on the tree/DART priors | latent refresh, sigma/k/DART draws |
| Per node op   | template instantiation selected at compile time; switch over a closed set | move type (moves.hpp), kernel table lookup |
| Per observation | monomorphic loops and kernel calls, no dispatch | partition compare, suffstat accumulation |

The leaf model (`ConstantGaussianLeaf` and its constrained variant
`MonotoneConstantGaussianLeaf`, `LinearGaussianLeaf`, `GPGaussianLeaf`,
plus the variance forest's `ConstantVarianceLeaf`, all in
`src/bartcore/model.hpp`) is a compile-time template
parameter `L` threaded through `Sampler<L>` -> `Chain<L>` -> `Forest<L>` and
the free functions in `src/bartcore/moves.hpp`, because it sits inside the
per-observation and per-node-op tiers (`accumulate`, `logIntegratedLikelihoodForNode`
must inline). The response model, by contrast, is chosen once per chain and
called through a handful of per-iteration virtual hops, so runtime
polymorphism costs nothing that matters; see "Model concepts" below.

## The facade

`src/bartcore/facade.hpp` defines `SamplerBase`, an abstract class exposing
every operation a host needs (run, the setters, tree/state serialization,
prediction, queries) and `SamplerFacade<L, ResidT = double>`, a `final`
class template that implements `SamplerBase` by forwarding each call to a
`Sampler<L, ResidT> impl_`. The leaf model is the parameter that matters for
dispatch; `ResidT` is a second compile-time parameter carrying the
opt-in fp32 residual storage (`options.fp32Residual`, `storage = "single"`),
which mints exactly one extra instantiation - the gaussian constant leaf
with `ResidT = float` - alongside the default `double` path every other
leaf model and family uses. This exists because the leaf model is a compile-time parameter: without a
type-erasure layer, every caller of the engine (the R bridge, the C ABI,
`tests/cpp`) would itself need to be templated on `L`. The facade collapses
that back down to one concrete type the rest of the package can hold a
pointer to.

Selection of `L` happens exactly once, in the free factory functions at the
bottom of facade.hpp - `createSampler`, `createSamplerOverStore`,
`createConstantLeafSampler`, `createAmplitudeSampler`,
`createMultinomialSampler` - called from the bridge at sampler creation
(`src/R_interface_bartcore.cpp`, `bartcore::createSampler(...)` and
siblings). Every factory that can build a variance forest asks one shared
predicate, `varianceForestIsRefused`, so the two cannot drift: a
variance-forest request (`options.numVarianceTrees > 0`) is refused under a
non-gaussian family, designated leaf covariates, a finite `residualDf` (the
Student-t augmentation shares the weight channel while reporting the
gaussian family, so it is its own term rather than a family test), or an
active monotone constraint. `createSampler` picks
`MonotoneConstantGaussianLeaf` when a monotone constraint is active without
leaf covariates, `ConstantGaussianLeaf` when no leaf covariates are
designated, `GPGaussianLeaf` when covariates are designated and
`options.gpLeaves` is set, otherwise `LinearGaussianLeaf`;
`createSamplerOverStore` (a store-backed variant used for row-subset views)
never selects the monotone leaf - it goes straight to `ConstantGaussianLeaf`
whenever no leaf covariates are designated. Every call
after construction goes through the chosen `SamplerFacade<L>` and pays one
virtual hop; nothing re-dispatches on `L` per call, per iteration, or per
observation.

## Model concepts and their shipped implementations

`src/bartcore/model.hpp` defines two independent extension points, dispatched
by two different mechanisms:

**LeafModel** - compile-time, via C++20 concepts. `LeafModelCore` is the
base requirement (a closed-form `logIntegratedLikelihoodForNode`); it
refines into `ScalarLeafModel` (one parameter per leaf: `ConstantGaussianLeaf`),
`VectorLeafModel` (`numParams()` doubles per leaf, fits evaluated per
observation: `LinearGaussianLeaf`), and `FunctionLeafModel` (one drawn value
per member observation, no per-leaf parameter storage: `GPGaussianLeaf`).
`IntegrableLeafModel` is the union of the three, and it is what the
templates in `chain.hpp` and `sampler.hpp` are constrained on; the free
functions in `moves.hpp` take the wider `MoveScorableLeafModel`
(`IntegrableLeafModel` or `ScaleLeafModel` - the latter admits the
heteroscedastic variance forest's `ConstantVarianceLeaf`, which scores
moves without satisfying the mean-model concepts).
Every shipped leaf model requires a closed-form marginal; there is
exactly one tree-structure sampling strategy today, the conjugate
Metropolis-Hastings moves in `moves.hpp` (see "Tree moves" below). A
non-conjugate strategy for leaf models without a working-Gaussian marginal
is designed for (docs/design/gp-leaves.md) but not implemented - GP leaves
did not need it, because a GP leaf under any working-Gaussian response is
itself integrable.

**ResponseModel** - runtime, via a virtual base class. `ResponseModel`
(model.hpp) owns the working response and weights the backfitting loop
reads, latent refresh, sigma draws, and the embedded-Gibbs mutation
entry points (setResponse, setOffset, setData, ...). The concrete class is
chosen with a `switch` on `ResponseFamily` (`enum class ResponseFamily
{ gaussian, probit, logistic, aft, ordinal, nbinom }`) inside `Chain`'s
constructor (`src/bartcore/chain.hpp`); the gaussian arm yields
`TResponse` instead of `GaussianResponse` when Student-t residuals are
requested (`resid.dist`), and the K-forest multinomial model installs
`MultinomialResponse` through its own construction path rather than the
enum. When `options.numGroups > 0` the chosen response is
wrapped in `GroupedResponse`, a decorator that Gibbs-samples per-group
intercepts into the offset between tree sweeps and forwards everything
else to the wrapped model. This choice is made once per chain at construction; every chain in a
sampler shares the same family.

**Split-variable selection**: `CGMTreePrior` (model.hpp) owns the
depth-decaying growth probability and the split-variable log-probability -
uniform over available variables by default, or weighted by
`splitProbabilities` when DART is active. `DartPrior` Gibbs-updates those
weights from each iteration's per-variable split counts (a normalized-gamma
Dirichlet draw), with an optional concentration (alpha) update sampled off
a fixed lambda grid with precomputed log-gamma terms. Both are plain,
non-virtual structs held by value inside `Forest<L>` - one instance per
forest, no dispatch at all, since a chain's split-selection policy never
changes after construction.

## Tree moves

`src/bartcore/moves.hpp` holds birth/death, change, and swap: the conjugate
Metropolis-Hastings proposals and their acceptance ratios, as free functions
templated on `MoveScorableLeafModel`. `metropolisJumpForTree` (a free
function in moves.hpp, called from `Chain`) is the per-iteration, per-tree
entry: it draws a step type
(`StepType::birth/death/swap/change`) and dispatches to the corresponding
move function. Every branch's empty-leaf veto is RANKED, not a flat
`-HUGE_VAL` (`Tree::leafVetoRank`, `moves.hpp#resolveVetoRank`,
docs/design/empty-leaf-veto.md): rank 2 is a leaf with no member at all,
rank 1 a leaf whose members all carry zero weight, rank 0 a leaf a
likelihood term reaches. Comparing a (current, proposal) branch pair, the
worse-ranked branch takes `-HUGE_VAL` outright; when both ranks are equal
the comparison runs on the finite log-likelihoods as usual. Only rank 2 is
absolute - no move may install a leaf with no member, from any state - so a
chain whose CURRENT branch is rank-1 vetoed (weights installed after the
tree was grown can strand a leaf this way) still mixes under the prior and
transition kernel at constant likelihood, rather than freezing. This keeps
member-empty leaves out of the chain state entirely, while a
weight-emptied leaf is penalized rather than forbidden.

Rules themselves are typed by column: ordinal rules compare a code against
a threshold, categorical rules test a bit of a direction mask. Masks up to
63 categories live inline in the rule's 64-bit word; wider ones store an
offset into a per-tree mask pool (`Tree::maskPool`, compacted between moves
once garbage passes a high-water mark). Both kinds share one word so a
single bit (position 63) doubles as the missing-value direction flag for
either rule kind.

## ColumnStore

`src/bartcore/data.hpp` defines `ColumnStore`, the predictor container
every chain in a sampler shares (chains never mutate it directly). Layout:

- **Codes**: `std::vector<xint_t> codes`, per-column quantized integer codes
  against per-column cut points (ordinal) or the raw category value
  (categorical). `xint_t` is `std::uint16_t` for every column.
- **Column kinds**: `std::vector<ColumnKind> types` marks each column
  `numeric`, `categorical`, or `orderedFactor` (`ColumnKind`, data.hpp). A
  categorical column's raw double IS its integer level index - `codeFor`
  casts it directly rather than binning it against cut points, so raw value
  and code coincide for a categorical column. An ordered factor carries the
  0..K-1 level codes and is binned against cut points over them, its order
  being meaningful. Rules, scans, masks and the flat replay branch on the
  DERIVED predicate `splitsBySubset(j)` rather than on the kind, so the
  mechanic (subset mask vs. threshold) is stated once.
- **Cut points**: `std::vector<std::vector<double>> cutPoints` plus
  `numCuts`/`categoryCounts`/`maxNumCuts` per column. `numCuts[j]` is the
  threshold count - 0 for a categorical column, whose `cutPoints[j]` stays
  empty, and K - 1 for an ordered factor; `categoryCounts[j]` is the fixed
  level count K of a factor column of EITHER kind, 0 for a numeric one.
  ColumnStore is the sole owner of cut construction and re-quantization, in
  three modes: uniform-over-range or quantile for a numeric column, selected
  by `useQuantiles`, and the K - 1 midpoints between consecutive declared
  level codes for an ordered factor, selected by the kind - so `n.cuts`
  bounds a numeric column's grid and applies to neither factor kind. Every
  other layer (moves, the tree prior) reads cuts through
  `ColumnStore::codeFor`/`column()`/`cutPoints`; nothing above the data
  layer computes or caches its own cut grid.
- **Ingestion channels**: the borrowed view a host builds from
  (`PredictorSource`) carries its dense values in two channels - doubles in
  `denseValues` and int32 level codes in `denseCodes`, with `denseChannels`
  naming which holds each column and where. A host whose factor columns are
  integers (which is what an R factor is) hands them over as integers, so a
  factor column's transient costs 4 bytes a cell rather than 8 and no cell
  is widened and narrowed again. The channel is a storage fact and the kind
  a reading one: both `build` arms and `buildTest` consume either, widening
  a coded column the store reads as numeric. The flat replay's reader
  (`PredictorSourceColumnReader`) reads either too and picks the channel
  once per node, so `predict` and the saved-tree replay route a factor test
  set off its own codes without widening a cell. Only the mutation entrances
  still need one dense block, column-major, and they materialize their own.
- **Raw values**: per-column storage and re-quantize source is one
  descriptor, `ColumnSource` (carried in `sources[j]` on each `CodeBlock`),
  discriminated by `ColumnSourceKind` - which names WHICH POOL the column's
  raw lives in, not who owns the codes: `denseCallSupplied` (the build-reset
  default; the raw arrives with the call - the call-time `x` on the train
  side, nothing on test), `denseResident` (re-quantizes from the
  descriptor's `residentRaw`, a slice of store-owned memory -
  `ownedDenseValues` on the train side of a mixed build, `ownedTestValues`
  on test - never a live host pointer), `denseCodesOnly` (no raw anywhere),
  and `cscRank`/`cscDensified` (both re-quantize from the descriptor's
  retained `slice`, a genuinely borrowed CSC values/rows pair that a
  mutation repoints at store-owned `ownedCscValues`/`ownedCscRows` (train
  side) on first write).

  **A column's stored representation is the one its values actually have.**
  A FACTOR column of either kind takes `denseCodesOnly` wherever its
  storage is dense: its cells are level codes, which the codes already
  carry, and its grid follows the level table rather than its values, so it
  never re-quantizes and there is nothing a double would serve. Both owned
  dense blocks therefore hold the REAL-VALUED columns alone, packed per
  predictor over the columns they serve. A dense build retains no raw at
  all past the build call. What a leaf model reads instead is gathered: the
  store owns a copy of each designated leaf-covariate column
  (`gatheredRawValues`) plus its standardization constants
  (`gatheredMeans`/`gatheredSds`), and `buildTest` gathers the test-side
  twin (`gatheredRawTestValues`, over its own column list) for the
  designated columns its block no longer holds. A mixed build gathers what
  its own block cannot serve on the same rule. The test side owns every raw
  it keeps whatever the view's shape: `ownedTestValues` (real-valued dense)
  and `ownedTestCscValues`/`ownedTestCscRows` (packed nonzeros of a
  mixed/CSC test build), so no borrowed pointer survives the call.
- **Sparse storage**: a CSC-built column at or below 20% nonzero density
  (`sparseDensityThreshold`) takes a rank-bitmap representation
  (`SparseColumnData`: a bitset, per-word popcount ranks, and packed
  nonzero codes); a denser CSC-built column densifies into the regular
  `codes` array. The tiering applies only to CSC-built columns - a dense
  build never takes rank-bitmap storage, however sparse its values.
- **Missingness**: `hasMissing` per column gates the extra missing-direction
  draw in rules; a reserved code (`naCode` for ordinal, a reserved category
  position for categorical) marks a missing cell.
- **Categorical tiering**: `hasPooledCategorical` is true when any column
  has more than 63 categories, which is the only condition that turns on
  the mask-pool machinery in `Tree`.

None of this is arena-allocated; every array is an ordinary `std::vector`
sized at build or mutation time.

## The mutation surface and its transaction semantics

The embeddable-in-a-larger-sampler contract (predictors/response/offset/
weights swappable between MCMC iterations) is implemented on
`bartcore::Sampler<L>` (`src/bartcore/sampler.hpp`) and exposed through
`SamplerBase`. Two families of mutation:

- **Response-side** (`setResponse`, `setOffset`, `setWeights`, `setSigma`,
  `setData`): forwarded straight to `ResponseModel`; no tree structure is
  touched, so these either apply or, for length mismatches, fail before
  touching state. `setData` also replaces the predictors, and so answers a
  status: it checks every factor cell of the replacement AND of any test
  matrix against each column's fixed level table before anything moves, and
  refuses the whole call rather than ingesting one side and dropping the
  other.
- **Predictor-side** (`setPredictor`, `updatePredictor`, and the
  per-observation session API): transactional. `Sampler` snapshots what it
  needs to restore (`oldCodes`, `oldHasMissing`, `oldCuts` when cut points
  would refresh), calls `ColumnStore::setPredictors`/`setColumns` to install
  the new values and re-quantize, then either force-updates or validates:
  with `forceUpdate` set, every chain force-refreshes its trees, collapsing
  any split that would empty a leaf into its parent with an
  effective-observation-weighted parameter merge, and always returns
  `accepted`; without it, `revalidateAllChains()` validates every tree of
  every chain first and only rebuilds fits if all stay valid, restoring the
  snapshot and repartitioning every tree on failure (`rolledBack`). A
  quantile-mode cut refresh that would induce fewer cuts than an existing
  split needs is rejected before any mutation happens (`invalidCutPoints`).
  Return type: `enum class PredictorUpdateResult { accepted, rolledBack,
  invalidCutPoints }`.
- **Per-observation updates** use a `PredictorUpdateSession`
  (`beginPredictorUpdate`/`updatePredictorPerObservation`): stage one
  observation's leaf moves against running per-leaf occupancy counts, test
  validity, then commit or skip, with one fits rebuild at the end rather
  than one per observation. `updatePredictorPerObservationJointly`
  (facade.hpp) sweeps several samplers that share an index-aligned column in
  the same randomized scan order, installing an observation everywhere or
  nowhere so their fits never diverge.

Data transactions are all-or-none **across every tree of every chain**:
validation runs over the whole sampler before any chain's fits are rebuilt.

## Tree storage forms

Four representations, used for different purposes. The last two are not
independent tree encodings - both are aggregates built on top of the flat
form below:

- **Live**: `Tree` (`src/bartcore/tree.hpp`) is a flat arena -
  `std::vector<Node> nodes`, children allocated as adjacent pairs so
  `rightChild == leftChild + 1` always holds. A `Node` carries its `Rule`,
  its `[begin, end)` span into the tree's external observation-index buffer,
  and (for scalar/function leaves) the constant-leaf sufficient statistic.
  This is what the moves and the backfitting loop operate on directly.
- **Flat**: `FlatNode` is one node of a pre-order-flattened tree - an
  ordinal cut point, a categorical direction mask (inline or pooled, tagged
  by `FlatKind`), or a leaf parameter, replayable against raw predictors
  without the `ColumnStore` that quantized them. This is the one format
  shared by saved-tree storage (`keepTrees`), external reporting
  (`getTrees`), and state serialization - `flattenTree`, `savedTree`, and
  `ChainStateData::forests[*].trees` all traffic in it.
- **State**: `SamplerStateData` (sampler.hpp) is the whole sampler's
  in-process serializable state - one `ChainStateData` per chain (itself
  one or more `ForestStateData`, since BCF chains carry two forests), the
  store's cut points, and two saved-tree write cursors. It is an aggregate,
  not a third tree encoding: `ForestStateData`'s tree fields, `trees` and
  `savedTrees`, are both `std::vector<std::vector<FlatNode>>` (combiner.hpp).
  A live forest's own saved-tree circular buffer (`Forest::savedTrees`,
  combiner.hpp) is flat too, which is why `getState` re-flattens every live
  tree into `ForestStateData::trees` but straight-copies the already-flat
  buffer into `ForestStateData::savedTrees`. The two write cursors are
  `currentSampleNum` (the next saved-tree slot, wrapping circularly) and
  `recordedDraws` (slots written since the last reset, capped at capacity -
  required rather than inferred, since an unwritten slot holds a zero-leaf
  tree that would read as a legitimate draw). `stateFormatVersion` is 3.
  Restoring a chain rebuilds everything else canonically (partitions from
  tree structure and cut points, `totalFits` by summing tree fits, the
  variance-prior anchor by re-running the same transform construction does)
  rather than replaying history, so a restored chain continues equivalently
  but not bit-for-bit - the last ulp of the original accumulation order is
  not reproduced.
- **Wire**: what actually leaves the process. `storeState`
  (`src/R_interface_bartcore.cpp`) flattens `SamplerStateData` into a
  struct-of-arrays SEXP that `setState` reads back, one list per chain. Each
  tree block's core is four parallel R vectors, built by `storeFlatTrees` -
  `tree.vars` (INTSXP), `tree.values`
  (RAWSXP, 8 bytes per node, so an inline categorical mask's bit pattern
  survives rather than being normalized by a REALSXP), `tree.sizes`
  (per-tree node counts), and `tree.flags` (the missing direction plus the
  `FlatKind` tag) - plus two optional side channels alongside them,
  `tree.params` (vector-leaf slopes) and `tree.masks` (the wide-categorical
  mask pool), present when the fit needs them.

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
`ext_rng*`. Nine further `GetRNGstate()`/`PutRNGstate()` pairs exist in the
bridge (`R_interface_bartcore.cpp`): two bracket `bartcore_run`'s own calls
into `sampler.run(...)`, even though nothing inside `sampler.run` touches
R's stream; the other seven bracket genuine R-level draws - a probit latent
redraw, two scan-order permutations for the per-observation update session,
and four entry points that can run outside a chain's own sweep loop
(`bartcore_sampleTreesFromPrior`, `bartcore_sampleNodeParametersFromPrior`,
`bartcore_growFromRoot`, `bartcore_drawLatents`).

## Threading model

`Sampler::run` (sampler.hpp) runs chains on
`min(options.numThreads, numChains)` `std::thread` workers; chains never
touch the `ColumnStore`, the R API, or each other while running, so results
are bitwise identical regardless of thread count (asserted by a `tests/cpp`
component test). Two execution paths:

- **Inline** (`numWorkers <= 1`): chains run sequentially on the calling
  thread. Progress prints directly; the interrupt poll runs between sweeps,
  throttled to roughly every 100ms.
- **Worker threads** (`numWorkers > 1`): each worker claims chains in a
  round-robin (`c = w; c < numChains; c += numWorkers`). Workers must never
  call into R: progress lines are queued (`QueuedProgressSink`) for the main
  thread to flush, cancellation is a relaxed atomic flag the main thread
  sets after polling, and on POSIX platforms `SIGINT` is blocked in worker
  threads (so a Ctrl-C is delivered only to the main thread, whose poll
  turns it into a cooperative cancel rather than an R longjmp across
  threads).
- **routeTestRows** (`chain.hpp`): the one pool that runs inside a chain
  during sampling rather than across chains - a `misc_mt` pool
  (`testFitPool_`) that fans a tree's test-row routing across this chain's
  share of the thread budget (`options.numThreads / numChains`). Routing
  draws no RNG and each row writes its own output slot, so results are
  bitwise identical at any thread count; below `testFitParallelCutoff`
  (65536 test rows) it runs serially, since dispatch would outweigh the
  routing it saves. General within-chain parallelism is deliberately not
  used (`docs/design/within-chain-threading.md`); this pool is the one
  sanctioned exception because it is RNG-free and output-disjoint.

**Callback restriction**: the per-sweep conditioning callback (`SweepCallback`
in the engine, `dbarts_sampler_callback` on the C ABI) fires on the calling
thread before every sweep and lets the host mutate conditioning state
(`setSigma`, `setOffset`, ...) between sweeps without a round trip. It is
only usable when chains run inline: `Sampler::run` requires the caller not
to set `onSweep` alongside worker-thread chains, and
`dbarts_sampler_setCallback` raises an error outright when
`numThreads > 1 && numChains > 1`. With several chains running inline, the
callback still only sees one chain progress at a time - chain c completes
every sweep before chain c + 1 starts.

**Prediction** (`Sampler::predictColumns`, sampler.hpp) mirrors `run`'s
worker-thread design: a per-call `n.threads` partitions the (chain, draw)
slab list across a `std::thread` fan-out with no cross-thread reduction, so
replay is bitwise identical at every thread count, including the inline
below-cutoff path. `numThreads = 0` floors to the sampler's own thread
count; the R5 `predict`/`predictForests` default to the fit's thread count,
and the formal is appended last on every generic (docs/design/
threaded-predict.md).

## Gates

Changes are classified by RNG effect and gated accordingly:

| Class | Meaning | Gates |
|---|---|---|
| neutral | draws unchanged | `tests/cpp` component tests; full tinytest suite |
| shifting | draws change, posterior does not | neutral's gates, plus regenerated RNG-locked tinytest snapshots, a re-recorded equivalence baseline, and the statistical (z) equivalence mode against the previous baseline |
| posterior-changing | stationary distribution or a default changes | shifting's gates, plus the exact-posterior gates (`.github/workflows/exact-gates.yaml`'s per-family list) and a design note in docs/design/ |

Any hot-path change, regardless of class, additionally needs
`benchmarks/R/bench-sampler.R compare` against a saved baseline, run on an
otherwise idle machine.

## Further reading

- docs/design/core-generalization.md - the generalization proposal and the
  full landing history behind everything above.
- docs/design/kernel-vocabulary.md - the compiled-kernel contract
  (`misc.a`) the generic engine dispatches into.
- docs/design/public-surface.md - what the R and C surfaces expose, and the
  cutover sequencing that retired the classic engine.
- docs/design/consumer-spec-surface.md - the exported `dbartsSpec()`
  resolution surface for embedding packages.
- docs/design/pooled-masks.md, sparse-columns.md, mia-missingness.md,
  linear-leaves.md, gp-leaves.md, grouped-random-effects.md, bcf.md,
  monotone.md, heteroscedastic.md -
  design and landing notes for each extension mentioned above.
