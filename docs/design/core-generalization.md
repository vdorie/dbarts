# Generalized BART core: design

Status: accepted plan, 2026-07-02. This document is the reference for the
core-generalization work; update it when decisions change.

## Goals and scope decisions

Generalize dbarts over data types (heterogeneous data.frame-like input,
categorical rules, sparse columns) and probability models (logistic responses,
non-constant leaf models, correlated observations, Dirichlet split-variable
priors), while preserving the package's distinguishing feature: the sampler is
embeddable in outer Gibbs/MH samplers, with data mutable between iterations.

Decisions:

- **Parallel new core.** The generalized engine is a new header-first C++
  library developed inside this repo. The existing R API and the
  `R_C_interface.hpp` ABI keep working throughout; the old engine is retired
  once the new one reaches parity.
- **C++20 floor.** Concepts define the model/data extension points.
  Requires R >= 4.3 toolchains.
- **Zero performance regression** on the classic Gaussian/probit path;
  cutover is gated on benchmarks.
- **No bit-parity requirement.** The port need not preserve the
  Node/index-array tree representation or the old engine's RNG stream.
  Validation is by exact RNG-free component tests, a statistical-equivalence
  harness, and regenerated test snapshots at cutover.
- **First shipped models** beyond Gaussian and probit: logistic via
  Polya-Gamma, and DART (Dirichlet split-variable prior). Linear leaves,
  correlated observations/random effects, and GP leaves are designed-for but
  implemented later.

## Current-state constraints (what the design starts from)

- Predictors are quantized once into per-column integer codes
  (`xint_t`, configure-selected width, default uint16); all tree-structure
  work is integer comparison. Hot loops are per-column: partitioning takes a
  column pointer and an integer cut (src/dbarts/node.cpp, partition kernels
  in src/misc/partition_body.c), and leaf sufficient statistics are vector
  reductions (src/misc/moments.c).
- SIMD kernels are compiled per-ISA in separate translation units and
  dispatched through function pointers set at load (src/misc/simd.c). This is
  why the core is header-*first*, not header-only: per-ISA compilation cannot
  move into consumers' translation units.
- The Gaussian conjugacy is behind a virtual interface
  (`EndNodePrior::computeLogIntegratedLikelihood`/`drawFromPosterior`,
  inst/include/dbarts/model.hpp) but its signature hardcodes the sufficient
  statistic as (mean, effective count, variance).
- Categorical rules exist (`Rule::categoryDirections` bitmask,
  src/dbarts/node.hpp; drawn in src/dbarts/treePrior.cpp) but are outside the
  active flow and capped at 32 categories.
- The response model is the least factored part: `responseIsBinary` is
  special-cased inline in the chain loop; y-rescaling, sigma draws, and probit
  latents are scattered rather than owned by a model object.
- `Data` borrows R memory directly (`REAL(slot)`); the quantized codes are a
  derived cache rebuilt from it.

## Architecture

Layers, bottom to top:

1. **Kernel library** (compiled; `misc.a` extended): a closed vocabulary of C
   functions compiled per (operation x code-width x ISA), dispatched through
   function-pointer tables. See kernel-vocabulary.md. Everything hot must be
   expressible as kernel calls.
2. **Generic core** (headers, C++20): the backfitting engine, tree storage,
   move machinery, written against concepts.
3. **Instantiations**: a curated, explicit list of `Sampler<Config>`
   specializations compiled in dedicated translation units.
4. **Facade**: a virtual `SamplerBase` exposing the mutation vocabulary
   (run, setPredictor, setOffset, ...). A factory maps the R-side model spec
   to an instantiation exactly once, at allocation.
5. **R package and C ABI**: `.Call` bridge and `R_C_interface.hpp` wrap the
   facade; LinkingTo consumers never see templates.

### Dispatch tiers

Dispatch is free when amortized over the work it gates; nothing dispatches
per observation.

| Granularity   | Mechanism                        | Gates          | Examples                          |
|---------------|----------------------------------|----------------|-----------------------------------|
| Per R call    | one vtable hop via `SamplerBase` | ms-scale work  | run, setPredictor, setOffset      |
| Per iteration | virtual hooks                    | O(n * trees)   | latent refresh, sigma/k/DART draws|
| Per node op   | virtual or switch on closed set  | O(n_leaf)      | move type, kernel table lookup    |
| Per obs       | none: monomorphic loops/kernels  | O(1)           | partition compare, suffstat accum |

The rule for what is a template parameter: **compile-time only if it touches
per-observation work or the `SufficientStat` type.**

- Template (header code): `LeafModel` (its `SufficientStat` is a value type
  living on the stack in branch scoring; `accumulate` must inline;
  `logMarginal` is called thousands of times per iteration), branch scoring,
  `MoveStrategy` internals (generic over the leaf model; selection among
  strategies is a per-tree switch over a closed set).
- Runtime virtual (shared across instantiations): `ResponseModel`,
  `SplitSelector`, variance/hyper priors, the facade. Concrete classes own
  their O(n) loops, so runtime polymorphism still yields monomorphic inner
  loops.

The instantiation matrix therefore collapses to roughly #leaf models
(constant, linear, index-span universal); each instantiation is one explicit
template instantiation plus a factory registry entry.

### Concept decomposition

- **ResponseModel**: owns the working weighted-Gaussian representation
  (z_i, w_i), latent refresh, y-scaling, and global variance draws. The
  backfitting engine is written once against (z, w) and is not templated on
  the response model. Gaussian: identity, hooks are no-ops. Probit:
  truncated-normal latents (current behavior). Logistic: Polya-Gamma omega_i
  become weights plus working response, riding the existing weighted kernels
  untouched. Grouped random effects: a decorator Gibbs-sampling per-group
  intercepts into the offset between tree sweeps (rbart_vi moved in-core,
  composing with any response model). Models with no working-Gaussian form
  use a non-conjugate MoveStrategy instead.
- **LeafModel** (concept; refined by `IntegrableLeaf`): declares
  `SufficientStat`, `accumulate`, `logMarginal` (if integrable), parameter
  type, tree-level parameter draw, `predict(param, x_row)`. Constant leaf:
  suffstat (sum w, sum wz, sum wz^2), delegating to existing moment kernels.
  Linear leaf: small (X'WX, X'Wz) blocks, still integrable. Universal
  fallback: the leaf's observation-index span (O(1) accumulate; cost lives in
  the model's own math). Leaf draws are at **tree granularity** (vector of
  leaf parameters given the tree), defaulting to independent per-leaf draws;
  this is the hook for monotone/shape-constrained BART. Leaf models consume
  raw column values, not codes: codes exist only for split structure.
- **MoveStrategy**: tree-structure sampling. Implementation 1 is the
  conjugate MH port (birth/death/change/swap on integrated branch scores;
  requires `IntegrableLeaf`, enforced at compile time). Designed-for later:
  joint structure+parameter proposals with prior-grown subtrees
  (reversible-jump style, optionally Laplace/pseudo-marginal scoring) for
  non-integrable leaves, and grow-from-root style samplers (enabled by the
  cut-scan kernel).
- **SplitSelector**: variable-selection prior; log-probabilities enter
  acceptance ratios; per-iteration hyper update hook. Uniform and DART
  (Dirichlet posterior over split counts) are the first two. Extensible to
  grouped/structured sparsity and interaction constraints.
- **Forest is the composable unit.** A Forest = trees + leaf model + split
  selector + its own working response and backfitting state. A `Sampler`
  orchestrates one *or more* forests plus glue parameters per iteration.
  v1 samplers have exactly one forest, but BCF (prognostic + treatment),
  multinomial (K forests), heteroscedastic (mean x variance forests), and
  hurdle models all reduce to several forests plus a response model that
  combines them.

### Data model: BartData

A purpose-built container (analogous to XGBoost's DMatrix); ingestion from R
types is a real step and the sampler never sees R types.

- **Cold layer**: canonical typed columns as provided (double, factor codes +
  level table, CSC sparse). Needed for cutpoint (re)computation, raw-value
  views for leaf models, and round-tripping to R. May borrow R memory while
  the R object lives; allowed to own.
- **Hot layer**: per-column bin-code arrays (per-column integer width chosen
  by cardinality: uint8 for <= 255 cuts, uint16 above), per-column cutpoint
  tables (shared with test containers so test rows bin identically), column
  metadata (kind, bin count, allowed rule vocabulary, NA policy, mutability),
  and both orientations where access patterns need them. Aligned,
  arena-allocated. A reserved NA code per column and a missing-direction bit
  on rules are designed in from the start (MIA-style missingness later).
- **Transactional mutation API**: column-granular replacement (with or
  without re-cutting), per-cell patches, and the per-observation
  propose/validate/commit-or-rollback cycle as an undo journal on the hot
  layer. Contract: common mutations cost O(changed cells), never O(n*p).
  Per-column-type semantics (factor level changes, sparse pattern changes)
  are defined here once.
- **Standalone R handle**: an external-pointer object constructible
  independently of any sampler; built once, shared read-only across chains
  and samplers (single-writer rule); row-subset views make CV folds a
  restricted initial index set over the same hot layer (xbart stops copying);
  serializable.

### Rules

Trees store compact tagged rules; the tag resolves per node operation
(one kernel-table lookup amortized over the leaf), so heterogeneous columns
cost nothing per observation. Vocabulary: ordinal (code > cut), categorical
membership (bitset over level codes; uncapped via inline bitset or pooled
storage past 64 levels), each with an optional missing-direction bit. Sparse
columns partition via a dedicated kernel (index set vs. nonzero rows);
order-preservation requirements to be settled by prototype.

## Extensions provisioned for

| Extension | Provision |
|---|---|
| BCF, multinomial, heteroscedastic, hurdle | multi-forest Sampler |
| Monotone / shape constraints | tree-level leaf draws |
| Grow-from-root / particle Gibbs, warm starts | cut-scan kernel; tree state import/export (FlattenedTrees) |
| Survival (AFT, discrete-time hazard), ordinal, quantile | ResponseModel latents; person-period expansion at ingest |
| Grouped/structured DART, interaction constraints | SplitSelector |
| Missingness (MIA) | reserved NA code + rule direction bit |
| Spatial/GP residuals | partially: grouped effects in-core; full GLS is a stretch goal |

Explicitly out of scope: soft BART (probabilistic routing replaces index-set
partitioning entirely; a different library sharing only the tree structure).

## Performance strategy

- Kernel vocabulary is the perf floor; the generic layer only decides which
  kernel on which column over which index set (O(1) per node op).
- Degradation ladder when no specialized kernel exists: (1) specialized SIMD
  kernel; (2) generic scalar kernel in the same table slot; (3) monomorphic
  header-template fallback. Never per-observation virtual dispatch.
- Verification: the classic instantiation calls the same kernel functions
  through the same pointers, so per-observation code is identical by
  construction; orchestration overhead is covered by per-operation
  microbenchmarks and the end-to-end zero-regression gate (benchmarks/).
  Diff the generated assembly of the classic tree-update loop once during
  phase 1.
- Honest costs, managed: binary size per instantiation (small matrix,
  measured; CRAN cares), template compile times (explicit-instantiation TUs),
  template debuggability (concepts at the extension points).

## Validation

1. Exact RNG-free component tests: branch scores, suffstat accumulation,
   acceptance ratios, posterior-draw parameters on fixed inputs, old vs. new.
2. Statistical-equivalence harness: matched synthetic scenarios, posterior
   summaries compared across seeds within Monte Carlo error
   (benchmarks/R/equivalence.R).
3. Regenerated tinytest snapshots at cutover (established convention for
   RNG-shifting changes).
4. The zero-regression benchmark gate (benchmarks/).

## Phases

0. **Groundwork** (DONE 2026-07-02): benchmark + equivalence harness
   (benchmarks/); kernel vocabulary documented as a stable interface.
1. **Core skeleton + classic parity** (first cut DONE 2026-07-02): C++20
   header library in src/bartcore/ (not yet wired into the package build);
   flat arena tree storage; conjugate MH moves generic over an
   IntegrableLeafModel concept; runtime-virtual ResponseModel with Gaussian
   and probit; ordinal rules. Gates passed: C++ component tests
   (tests/cpp/), statistical equivalence vs the old engine (20 seeds x 3
   scenarios, max |z| = 2.55), and throughput 12-16% FASTER than the old
   engine across n in {100..10000}, trees in {75, 200}. Caveat (VD): that
   is an end-to-end engine comparison, NOT a generics-cost measurement --
   the engines call different kernel variants (unrolled-vs-online moments
   before the fix, aligned-vs-unaligned axpy still). Isolating abstraction
   cost requires pinning both engines to identical kernels; do this before
   claiming the dispatch design is free.
   Not yet ported (phase 2+): multiple chains/threads, thinning, keepTrees,
   the k hyperprior (ChiHyperprior), quantile-based cut points, sampler
   mutation API (setPredictor and friends), callbacks, MATCH_BAYES_TREE.
   Availability is recomputed from ancestor walks instead of cached
   per-node flags; Forest is not yet split out of Sampler (do this when the
   facade lands in phase 2).
2. **R boundary + DART** (DART half DONE 2026-07-02): DART implemented in
   bartcore as the split-selector seam: CGMTreePrior carries optional
   split-variable weights (the classic engine's splitprobs semantics, ported
   and cross-engine equivalence-tested), and DartPrior Gibbs-updates them
   from forest split counts (normalized-gamma Dirichlet draw) with an
   optional concentration update sampled on a fixed lambda grid with
   precomputed lgamma terms (O(grid) multiply-adds per iteration).
   Validated: component tests of the Dirichlet update, alpha adaptation
   under sparsity, and a sparsity-recovery behavioral test (signal mass
   0.24 uniform -> 0.99 DART on p=25 with 2 signal variables).
   Boundary work (DONE 2026-07-02): SamplerBase facade + factory
   (bartcore/facade.hpp); embedded-Gibbs mutation API on the new engine
   (setOffset with/without rescale, setResponse incl. probit latent redraw,
   setSigma, latents, setTestPredictors); the package now builds as C++20
   (CXX_STD in Makevars.in/.win) with an internal .Call bridge
   (src/R_interface_bartcore.cpp, R/bartcore.R, not exported) constructed
   from a dbartsSampler's validated control/model/data; tinytest
   inst/tinytest/test-bartcore.R covers the surface. Full existing suite
   (1724 results) passes under C++20 -- the classic engine is unchanged.
   Chi hyperprior on k (DONE 2026-07-02): ported with the sum of squared
   leaf parameters accumulated during the parameter sweep; the default
   binary specification now runs on bartcore, run results include k
   samples, and a sampler-API equivalence scenario gates it cross-engine.
   Remaining for phase 2: predictor mutation (setPredictor and friends,
   incl. per-observation rollback) on the new engine, quantile-based cut
   points, multiple chains/threads, and the user-facing dbartsSampler
   opt-in flag once those reach parity.
3. **Logistic (Polya-Gamma)**: PG sampler in external.a; exercises latent
   hooks and the weighted path.
4. **Data generalization**: BartData container, data.frame ingestion,
   categorical rules enabled (> 32 levels), sparse columns, new kernels,
   per-column-type mutation semantics, standalone R data handle.
5. **Wave 2 models**: linear leaves; in-core grouped random effects
   (retiring the rbart_vi R loop).
6. **Non-conjugate MoveStrategy**: GP leaves, general likelihoods.
7. **Cutover**: new core default, R_C_interface mapped onto it, old engine
   deleted, headers published for LinkingTo.

## Risks

- Instantiation bloat / compile times / .so size: bounded by the small
  matrix; measure from phase 2.
- Sparse partitioning performance is open; prototype before committing to a
  representation (order-preserving partition likely needed).
- Polya-Gamma draw cost (one per observation per iteration) dominates
  logistic; choose sampler variant by benchmark.
- setPredictor/rollback semantics for factors and sparse patterns need
  careful per-column-type definitions (phase 4).
- MATCH_BAYES_TREE and similar compat flags: decide in phase 1 whether the
  new core carries them or they remain old-engine-only until cutover.
