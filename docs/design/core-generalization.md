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
   Re-measured at full R5 parity (d31c0f8, 2026-07-03) with
   bench-sampler.R engine=new, after multiple chains, the facade,
   categorical branches in the hot loops, keepTrees, and the progress
   plumbing all landed: the zero-regression cutover bar HOLDS. Sampling
   throughput 6-12% faster than classic (continuous n in {1000, 10000},
   trees in {75, 200}; binary 6-9% faster); embedded-Gibbs
   (setOffset + run(0, 1)) and both setPredictor paths at par (accept
   1.01, reject 0.98, interleaved A/B medians). The harness's original
   random-replacement setPredictor metric conflated the engines'
   accept/reject mixes (acceptance depends on the chain state - one
   single-observation leaf rejects most candidate columns); it now times
   the accept path (identity swap) and reject path (degenerate column)
   separately. Same caveat as above: this is the user-facing cutover
   gate, not a generics claim.
   Not yet ported (phase 2+): callbacks deliberately wait for the new
   public C surface (see phase 4 notes). MATCH_BAYES_TREE is not preserved
   (old-engine-only until cutover; see Risks). (Multiple chains/threads,
   thinning, quantile-based cut points, the k hyperprior, the sampler
   mutation API, and keepTrees/serialization were all ported later; see
   below.)
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
   Predictor mutation (DONE 2026-07-02): full setPredictor family ported
   with the classic transaction semantics. ColumnStore gains in-place
   mutation (setPredictors pointer swap, setColumns, setCell); the sampler
   snapshots x/codes/cuts, re-routes every tree, and rolls everything back
   if any leaf would empty (validate-or-rollback), or with forceUpdate
   collapses emptied splits into their parents with
   effective-observation-weighted parameter merges (collapseEmptyNodes
   port). Per-observation updates use a session object (occupancy counts +
   cached leaf assignments, randomized scan order, sequential commit) and a
   type-erased PredictorUpdateSession on the facade supports the joint
   all-or-none sweep across samplers sharing an index-aligned column
   (updatePredictorPerObservationJointly). Internal bridge + tinytests
   cover the whole surface.
   Port bug found by the mutation tests and fixed: the branch likelihood
   was missing the classic engine's empty-leaf veto (-1e7 for any branch
   containing an empty leaf), so bartcore chains could accept and carry
   empty leaves. With the veto restored, DART needed the established
   activation delay (the BART package starts DART halfway through burn-in);
   DartPrior.updateDelay holds s uniform for that many updates, restoring
   reliable sparsity recovery (10/10 seeds >= 0.997 signal mass, vs 3/10
   when updating from a cold forest).
   Cut-point generalization (DONE 2026-07-02): ColumnStore carries
   per-column cut counts and caps (heterogeneous n.cuts now accepted by the
   bridge) and a quantile mode (control.useQuantiles): cuts at sorted
   unique-value midpoints, thinned to the per-column cap with the reference
   engine's step/offset rule, counts fixed after construction. Predictor
   updates that refresh cuts pre-check quantile feasibility (a column
   inducing fewer cuts than existing splits require returns
   invalidCutPoints without mutating anything -- an improvement on the
   reference engine, which errors midway through installation); the
   transactional entry points now return
   accepted/rolledBack/invalidCutPoints. The explicit setCutPoints entry
   point installs externally chosen (strictly increasing) cuts per column,
   re-quantizes train and test codes, and force-refreshes trees so
   orphaned splits collapse.
   Multiple chains and threads (DONE 2026-07-02): the engine splits into
   Chain (bartcore/chain.hpp: trees, fits, response state, chain
   parameters, and its own rng over a shared read-only ColumnStore) and
   the Sampler coordinator (owns the store and the chains, fans mutation
   out, and runs chains on up to min(numThreads, numChains) std::thread
   workers -- chains never touch the store, the R API, or each other, so
   results are identical for any thread count, which a component test
   asserts bitwise). Data transactions are all-or-none across every tree
   of every chain like the classic engine's: validation runs over all
   chains before any chain's fits are rebuilt. Results arrays gain a
   trailing chain dimension (chain-major slabs), matching the classic
   run() shapes. The bridge honors control n.chains/n.threads; a single
   chain still draws through R's generator, while several chains each get
   a Mersenne twister seeded from R's stream. Within-chain (htm) kernel
   threading remains unported and is likely superseded by chain-level
   parallelism.
   Opt-in engine flag (DONE 2026-07-02): dbartsControl(engine =
   "bartcore") runs the standard dbartsSampler surface on the new engine.
   Supported methods delegate one line into R/bartcore.R helpers that
   mirror the classic semantics (run resolves control defaults and
   returns classic-named, chain-dimensioned results; setPredictor keeps
   the full/column/partial forms, forceUpdate defaults, data@x swap on
   success, and mask returns; setCutPoints, setTestPredictor(AndOffset),
   setResponse, setOffset with rescale, setSigma, getLatents, getSigmas).
   Unsupported methods (state serialization, keepTrees/tree extraction,
   predict, setWeights, test offsets, prior sampling) refuse
   loudly instead of passing a bartcore pointer to classic entry points.
   The engine gained classic-parity thinning (numThin: burn-in and
   samples count at the kept rate); the bridge honors keepTrainingFits
   and pins borrowed vectors in fixed slots of the external pointer's
   protection object so repeated setters do not accumulate. The
   equivalence harness drives both engines through the public sampler
   surface via the flag; flag-path results are RNG-identical to the
   internal bridge's.
   Whole-data replacement (DONE 2026-07-02): Sampler::setData swaps
   predictors, response, weights, offset, and test predictors at once,
   possibly resizing numObservations (the predictor count is fixed).
   Two phases: every chain recovers its leaf parameters against the old
   fits and partitions, then the store rebuilds cut points from scratch
   (quantile counts may shrink here, unlike the transactional updates)
   and each chain remaps split indices onto the value-nearest new cut
   within the ancestor-constrained interval
   (Tree::mapOldCutPointsOntoNew; interval-starved subtrees collapse
   with plain-mean parameter merges), re-routes, collapses anything left
   empty, and rebuilds fits. ResponseModel::setData re-borrows and
   resizes response state: Gaussian keeps sigma and the variance prior
   fixed on the original scale across the rescale, probit
   cold-initializes latents to 2y - 1, both as the reference engine
   does. The public sigma accessors were also brought onto the original
   scale, symmetric with setSigma (getSigmas previously leaked the
   internal scale). The R5 setData delegates with classic semantics
   (n.cuts/sigma slot preservation, data-slot rollback on error);
   phase 2 is complete.
3. **Logistic (Polya-Gamma)** (DONE 2026-07-02): ext_rng_simulatePolyaGamma
   in external.a draws PG(1, psi) exactly with Devroye's alternating-series
   method (component-tested against closed-form moments on both proposal
   branches). LogisticResponse rides the weighted conjugate path: given
   eta_i = f(x_i) + offset_i, omega_i ~ PG(1, eta_i) and the engine sees
   working response (y_i - 0.5) / omega_i - offset_i under per-iteration
   precision weights omega_i, with sigma fixed at 1; the cold start sets
   omega to its PG(1, 0) mean of 1/4. To support weights that change every
   iteration, the engine's precisions now come from
   ResponseModel::workingWeights() (user weights for gaussian, the omega
   draws for logistic) instead of a fixed chain member. The ctor bool
   became ResponseFamily {gaussian, probit, logistic}; the bridge takes a
   family argument, reachable via dbarts:::bartcoreSampler(sampler,
   family = "logistic") until the new public surface lands; latents()
   exposes the omega draws. Alongside, the classic engine's weighted
   probit was STRIPPED rather than ported (scaling the latent draws by
   1 / sqrt(w) is not a coherent weighted-likelihood model - Vincent);
   bartcore rejects weights with binary responses at creation and setData.
   No classic reference exists for a cross-engine gate; validated by exact
   PG moment tests, end-to-end log-odds recovery and calibration, and
   mutation smoke tests, plus an exact-posterior gate
   (benchmarks/R/logistic-reference.R): a single-tree, one-predictor,
   3-cut problem admits brute-force tree enumeration and quadrature leaf
   marginals, and both binary families match the exact posterior
   predictive to MC error (logistic within 0.001). The same script
   compares against the BART package under exactly matched priors, which
   surfaced a finding worth keeping: BART::lbart itself deviates from the
   exact posterior by ~0.03 in the anti-shrinkage direction on the tiny
   problem, and shows the same systematic ~0.02 probability gap
   (attenuation-signed, stable in chain length) on a 200-tree problem, so
   the lbart comparison is informational rather than a gate; pbart agrees
   with bartcore probit to ~0.005.
4. **Data generalization**: BartData container, data.frame ingestion,
   categorical rules enabled (> 32 levels), sparse columns, new kernels,
   per-column-type mutation semantics, standalone R data handle.
   Wide categorical masks (DONE 2026-07-03): direction masks are 64 bits,
   with the category cap raised from 32 to 53 - the widest mask whose every
   value the flattened format's double encoding represents exactly.
   Assignment patterns are now drawn bit by bit with the two all-same
   patterns rejected, exactly uniform for any width; the single range draw
   this replaced had only the generator's granularity, which for wide masks
   pinned low pattern bits to functions of the high ones (and shifted the
   categorical RNG stream in the process; ordinal paths are bit-identical,
   which the equivalence gate's exact reproduction confirms). The union
   with the ordinal split index zero-initializes at full width so rule
   equality can compare the wide member for both kinds. Uncapped categories
   need pooled mask storage and a wide-mask reporting format; deferred to
   the public-surface proposal (docs/design/public-surface.md).
   Data.frame ingestion (DONE 2026-07-03): `dbarts` and `dbartsData` take
   `factors = "categorical"`, coding unordered factors as single
   categorical columns (codes 0..K-1, K <= 53), ordered factors as
   ordinal, with per-column level tables retained on the model matrix so
   test data recode against the training levels (unseen levels error).
   The default `"indicators"` path is unchanged, and the classic engine
   refuses categorical columns at sampler creation (its categorical
   machinery is unreachable, untested code). The default flips at
   cutover per the public-surface plan.
   Priors as objects + DART exposure (DONE 2026-07-03, per the
   public-surface review): the cgm/normal/chisq/fixed/chi constructors
   are pure standard-evaluation functions returning validated S4 specs,
   exported only through the `dbartsPriors` container; parsePriors
   layers that vocabulary over the caller's environment (bare names keep
   working, immune to attach-order masking) and resolves the
   data/control-dependent parts at fit time - named split probabilities,
   the binary k default, and dart's update delay (half of n.burn, the
   startdart convention). `tree.prior = dart(...)` wires
   SamplerOptions::useDart through the bartcore bridge; the classic
   engine and setModel refuse DART. Per-sample varprobs landed
   2026-07-03: Results::splitProbabilities (filled only under DART),
   returned as "varprobs" from runs and packaged onto bart2 fits next
   to varcount; bart2 takes dart = TRUE or a full spec object.
   family exposure (DONE 2026-07-03): `family` on `dbarts()`/bart2
   resolves into `dbartsControl@family` ("auto" keeps the old dispatch;
   "gaussian" on 0/1 responses fits continuous; logistic is public on
   the bartcore engine with node.scale = pi * sqrt(3)); the slot drives
   the bridge's family argument including pointer re-creation after
   save/load, and setControl preserves it like `binary`. The wrappers'
   probability transforms are link-aware (2026-07-03): packaged bart and
   rbart results carry a `family` element, and
   predict/extract/fitted/plot transform latents through it
   (probabilityFromLatents in generics.R; missing element = probit, so
   old saved fits keep their meaning). rbart remains probit-only by
   construction - its ranef update is a normal-conjugate step.
   xbart backend (DONE 2026-07-03): the C++ crossvalidation monolith
   (crossvalidate.cpp, R_interface_crossvalidate.cpp) is deleted and
   xbart is an R-level driver over the dbartsSampler mutation API
   (setData for splits, setModel for cells), so it runs on either
   engine and inherits every engine feature; parallelism is a
   `parallel` cluster over replications. FINDING, fixed rather than
   replicated: the old driver warm-started chains ACROSS data splits.
   The held-out rows of a fold were training rows of the previous one,
   and slowly-mixing settings (a 2-tree ensemble's deep trees) never
   forgot them within the rep burn-in - measured 0.11 rmse of optimism
   at 400 burn-in iterations, enough to reverse a 2-versus-75-tree
   comparison. Warm starts now happen only across parameter cells over
   an unchanged split, which leaks nothing; every split starts a fresh
   chain with the full n.burn[1], and n.burn[3] is defunct. Also fixed:
   weighted rmse was sqrt(WSSR)/sum(w) (wrong units), binary losses on
   continuous responses silently became rmse, the multi-grid result
   array's dimension labels did not match its fill order, and per-fold
   offsets were dropped. The reproducibility contract is (seed,
   n.threads), with the caller's RNG stream restored when a seed is
   supplied; the driver draws splits through R's generator, so its
   reproducibility test pins sample.kind against leakage from other
   test files.
   THE FLIP (DONE 2026-07-03): dbartsControl defaults to engine =
   "bartcore" (class prototype included) and factors defaults to
   "categorical" on dbarts/dbartsData/bart2/xbart (xbart gained the
   argument and threads it to dbartsData); bart stays the frozen
   BayesTree shim - indicators and probit, explicitly requested - but
   runs the new engine. Running the whole wrapper surface on bartcore
   flushed out parity gaps the R5 harness never saw because it compares
   test fits, not training fits, and never exercised some paths:
   (1) recorded training fits omitted the offset while the classic
   engine includes it - rbart's ranef Gibbs reads train - ranef, so the
   ranef re-entered the residual and diverged geometrically (the
   "hang"); fixed engine-side (ResponseModel::offset(), added back in
   storeSample, symmetric with test fits). (2) control@rngSeed was
   ignored: now a single chain seeds R's generator with it (classic
   convention) and several chains draw their MT seeds from a dedicated
   generator, leaving R's stream untouched; seeded results are
   thread-count independent (classic's were not; its unseeded
   multichain clock-seeding is NOT replicated - unseeded chains seed
   from R's stream so set.seed suffices). rngKind/rngNormalKind are
   refused on bartcore. (3) test data could not be removed - bart2/bart
   null out test predictors for burn-in; buildTest(NULL, 0) transitions
   to the supported no-test state, offset cleared. (4) getLatents
   refused preallocated results, which rbart fills in place; the bridge
   now honors the classic storeLatents contract. (5) setResponse and
   setOffset skipped length validation (numeric(0) segfaulted).
   (6) copy() walked classic dbartsState slots; the bartcore branch
   installs the opaque state object in the duplicate (never mutated in
   place, so sharing is safe). (7) the state slot is now a delayedAssign
   like classic's, so forcing it before saveRDS captures the sampler.
   (8) fixed residual priors work on bartcore with the DOCUMENTED
   variance semantics (sigma = sqrt(value)) at create and setModel
   alike; classic ignored the value at create and installed it as an sd
   in setModel - not replicated. Mutation-error messages align with
   classic's wording. RNG-locked snapshots regenerated
   (binary/continuous regression, rbart, xbart); classic-mechanics
   tests (customMCMC's state cutPoints, rng kind semantics) pin
   engine = "classic" until removal; a knife-edge statistical bound
   (flat-hyperprior median k < 3; cross-seed medians 2.5-69 on BOTH
   engines) was replaced with a robust check. Suite 1975; equivalence
   reproduces all nine baseline z-stats exactly; exact-posterior gates
   pass.
   Classic removal, R surface (DONE 2026-07-03, cutover step 3b): the
   engine argument/slot, rngKind/rngNormalKind, and dbartsState are gone;
   every dbartsSampler method runs its bartcore body unconditionally
   (startThreads/stopThreads stay as no-ops for wrapper callers);
   binary+weights refuses at creation everywhere. Classic-pinned tests
   deleted (test-sampler-customMCMC) or reworked to engine-free forms
   (test-rng, test-bartcore's prior-KS and verbose byte-comparisons are
   now bartcore-only structural checks; binary-weights host creation is
   itself the refusal). Harnesses updated: equivalence's samplerApi path
   and bench-sampler run the installed engine; engine=new still selects
   the rshim path and old classic baselines remain the permanent
   cross-engine reference (re-recording classic is impossible). Suite
   1990 across 57 files; equivalence vs the saved baseline still
   statistically indistinguishable. The classic C++ is unreachable and
   deleted next (step 3c).
   Flat C API (DONE 2026-07-03, cutover step 3a): inst/include/dbarts/
   dbarts.h, the LinkingTo replacement for R_C_interface.hpp, sized by
   an audit of stan4bart 0.0-13 (the only compiled-boundary reverse
   dependency; findings and the resolved v1 surface in
   public-surface.md section 6). Opaque dbarts_sampler handle over the
   bartcore facade, dbarts_sampler_* CCallables (classic still exports
   the unprefixed dbarts_* symbols until removal), SEXPs only at
   creation (the R spec objects) and state/trees (R-serializable by
   purpose), run writes into a caller-owned dbarts_results with
   null-means-skip - which the engine already honored end to end.
   Implementation shares the bridge's cores, extracted into
   bartcore_bridge:: (createHolder, storeState, setState, getTrees,
   validateColumnValues; declared in R_interface_bartcore_common.hpp) -
   C_interface.cpp is a thin adapter. Tested by an actual consumer:
   test-capi.R compiles inst/tinytest/capi/consumer.c with R CMD SHLIB
   against the installed headers, resolves everything through
   R_GetCCallable, and drives the stan4bart workout (caller buffers,
   seeded bitwise reproducibility, offset/fixed-sigma conditioning,
   probit latents, predict == recorded test fits, state round trip with
   identical predictions); it self-gates on toolchain availability.
   Suite 2013 (was 1975); equivalence exact.
   Categorical splits (DONE 2026-07-02, engine-side): a finding first -
   the classic engine's categorical machinery is DEAD CODE from R (nothing
   ever assigns CATEGORICAL_VARIABLE; factors are dummy-expanded by
   makeModelMatrix, and the classic cut-map code never branches on type),
   so there is no cross-engine reference and no compatibility to keep.
   bartcore's design is clean-room: ColumnStore carries per-column
   ColumnType; categorical columns hold integer category codes 0..K-1
   directly (K <= 53, fixed at build), Rule holds a
   splitIndex/categoryDirections union, and observations route by bit
   test, one type dispatch per node operation. Rules live in a canonical
   gauge - direction bits confined to the categories reachable at the
   node (an ancestor-filtered bitmask), neither side empty - with the
   prior uniform over the 2^R - 2 such assignments, both orientations
   counted. Every kernel is closed on that space: birth/change draw in
   gauge, the swap validity walk checks gauge and nonemptiness for the
   swapped variables, and collapse/setData only enlarge reachability,
   which preserves gauge. The change move rejection-samples assignments
   whose descendant splits stay satisfiable (the good set depends only on
   the tree above and below the node, so proposals cancel and an
   exhausted draw aborts symmetrically, mirroring the ordinal
   findGoodOrdinalRules flow). Deviations from the classic design, all in
   unreachable code: the classic randomized unreachable categories'
   directions, pinned the first reachable category right, drew from half
   the assignment space while weighting rules by a positively-signed
   normalization in the tree prior, and never remapped rules below
   categorical nodes in setData. Gated by
   benchmarks/R/categorical-exact.R (tree-space enumeration + quadrature
   under the probit link; the sampler matches the exact posterior
   predictive to MC error) plus component tests over the mechanics,
   prior math, recovery, and mutation surface. Not yet exposed publicly:
   dbartsData still types every column ordinal, so tests flip varTypes by
   hand; factor ingestion without dummy expansion is public-surface work.

   Tree storage + serialization (DONE 2026-07-02; cutover-parity work for
   phase 7): one flattened representation serves everything - pre-order
   FlatNode records {variable, value} with splits value-encoded (an
   ordinal rule's cut point, a categorical rule's direction mask) and
   leaves holding their parameter. Cut values map back to split indices
   exactly (cuts are unique doubles), so the same format serves keepTrees
   storage, getTrees, predict, and state serialization; replay against
   raw predictors (x <= cut goes left, mask bit sends right) reproduces
   the engine's code-based routing exactly. keepTrees mirrors the classic
   circular buffer: capacity = n.samples at creation, trees flattened
   inside the tree loop while the freshly drawn parameters are live,
   currentSampleNum advancing per run. getTrees/predict are drop-in
   classic formats (pre-order n/var/value rows, var 1-based with -1
   leaves, internal-scale leaf values; predictions numTest x samples x
   chains with keepTrees, else per-chain live-tree fits; saved trees
   replay training data or newdata for the counts). State serialization
   (Sampler::get/setState through an opaque bartcoreState R object)
   restores BITWISE-exact continuation, which took three fields beyond
   the obvious ones: sigma on the internal scale (the original-scale
   round trip can drop a bit), the accumulated totalFits (its
   floating-point accumulation history differs from a fresh tree-order
   sum), and each tree's observation index buffer (within-leaf segment
   order feeds the suffstat sums; restored partitions are set by an
   order-preserving walk, Tree::setPartitionsFromOrderedIndices).
   setState validates every chain against the state's own cut points on
   scratch trees before touching anything, so a bad state rejects with
   nothing changed. The R5 surface honors updateState for bartcore, and
   getPointer transparently re-creates a save/load-ed sampler from its
   stored state, matching classic workflows; component tests gate the
   round trips (gaussian multi-chain + keepTrees, logistic latents +
   DART) bitwise. Also fixed in passing: bartcore_setData retained a
   leftover ordinal-only check that made its categorical validation
   unreachable.

   setWeights + test offsets (DONE 2026-07-02, cutover-parity cleanup):
   both mirror the classic engine exactly. setWeights is a bare pointer
   swap fanned to every chain's response model (gaussian only; the binary
   families reject, their reference weighting having been stripped) -
   nothing rescales, and the weighted residuals enter the next iteration's
   node statistics and sigma draw, so installing weights before any run is
   bitwise identical to creating with them. Test offsets live on the
   ColumnStore beside the test predictors and are added to recorded test
   fits in storeSample (predictions instead take an offset argument, as
   classic predict's NA sentinel never falls back to the stored offset).
   The R5 surface gets the full classic offset plumbing:
   testUsesRegularOffset sync in setOffset, setTestOffset breaking the
   sync, and setTestPredictorAndOffset for row-count changes; a lone
   setTestPredictor that would orphan the offset's length refuses (the
   classic engine silently keeps the stale pointer). Gated by exact
   component tests (offset shifts recorded fits by itself; setWeights
   equals creation with the new weights) and a ninth equivalence scenario
   (wtoffset) swapping weights and installing a test offset mid-chain.

   Prior sampling + printTrees (DONE 2026-07-02, cutover-parity cleanup):
   sampleTreesFromPrior replaces every tree's structure with the classic
   recursion (depth-decayed Bernoulli growth, rules from the tree prior,
   empty children keep growing until a final collapse; fits left stale,
   which run() tolerates since totalFits still sums the per-tree fits);
   sampleNodeParametersFromPrior draws every leaf from the node prior and
   rebuilds tree/total/test fits. printTrees reproduces both classic
   formats: Node::print for live trees (occupancy, TBN flags, per-variable
   availability recomputed from rule intervals, rules by index=value) and
   SavedNode::print for saved slots (values only), with the classic
   chain/sample headers and indent bookkeeping; categorical rules print
   their direction bits (the classic categorical print branch never
   terminates, being dead code). Gated by a component test (growth
   frequencies at depth 0/1 match base/(1+d)^power, leaf-parameter moments
   match N(0, (scale/k)^2), occupied well-formed trees, run continues from
   a prior-drawn state) and cross-engine tinytests (leaf-count KS and
   node-prior spread, classic vs bartcore).

   Callbacks (DECIDED 2026-07-02): not ported. Control::callback is
   reachable only through the C ABI - the R bridge never sets it - and the
   old engine serves that ABI until cutover, so a bartcore equivalent has
   no possible consumer today. Design one with the new public C surface
   instead.

   setControl + setModel + getSumsOfSquaredResiduals (DONE 2026-07-02,
   final R5 parity slice; every dbartsSampler method now runs on the
   flag). setControl covers what the bart/rbart/pdbart wrappers change
   mid-flight - keepTrainingFits, keepTrees flips (saved-tree storage
   reallocates and the write position resets; a no-op reconfiguration
   preserves stored samples), n.threads, n.thin, and the R-side-only
   verbose/updateState/n.burn/n.samples - and refuses what creation fixed
   (engine, chain and tree counts, generators, cut-grid settings), where
   the classic engine instead resizes chain state nobody resizes. setModel
   maps a dbartsModel onto a ModelParameters install fanned to every
   chain: tree prior base/power, move probabilities, node scale, k or its
   chi hyperprior, fixed split probabilities (DART samplers keep their
   Dirichlet machinery - a dbartsModel cannot express DART), and the
   variance prior re-anchored to the original-scale sigma estimate exactly
   as construction, so installing before any run is bitwise identical to
   creating with the model (component-gated, fixed-k and chi-k variants).
   getSumsOfSquaredResiduals descales by range squared for continuous
   responses and reports binary responses on the latent scale; the
   classic version multiplied the internal SSR by one factor of the range
   (and by the binary placeholder range of 2), a units slip with no
   consumers - fixed there too (Vincent, 2026-07-02), so the engines
   agree.

   Verbose output (DONE 2026-07-02): full classic parity. The creation
   summary (BARTFit::printInitialSummary) is reproduced in the bridge
   from the same Control/Model/Data structures and is byte-identical
   across engines over the same inputs (tinytest-gated); the quantile
   and internal-scale lines reduce at creation to expressions over the
   raw prior scale and the response range (a new fitScale facade
   accessor). Runs print the loop header, "iteration: k (of N)" every
   printEvery kept iterations ("[c] "-prefixed under multiple chains),
   and the terminal summary (accumulated loop seconds, leaf counts, and
   variable usage). Chains hand formatted lines to a ProgressSink:
   inline runs print directly; worker threads queue lines behind a mutex
   and the main thread flushes every 0.1 seconds, since workers must
   never call into R (the classic engine's htm-buffered scheme, without
   the thread manager). setControl toggles verbose/printEvery live.

   Remaining unsupported surface: weights + binary (by design) and
   setControl changes to creation-fixed settings.
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
- MATCH_BAYES_TREE and similar compat flags: DECIDED (Vincent, 2026-07-02)
  - not preserved; they remain old-engine-only and die at cutover.
