# sampler-core-and-family-records: one constant record per kind of response, and a sampler interface of a core plus a queried capability object

Status: PLANNED 2026-11-10, not built. Two TODO items planned together, as their rulings ask: family-fact-records
(dec-B767) and sampler-core-and-capabilities (dec-B763). Awaits a blind critique, then collapsed-move-rebuild's
landing.

agent: opus implementer, one per slice; opus reviewer running the mutants below, effort high.
rng: NEUTRAL in every slice. No draw, default, message or fit object moves; the equivalence trio is bitwise.
window: pre-release, before the host-neutral engine's S3 (the C creation entry and the facts struct), and best
before its S2a (below). dbarts.h does not change; its hash does not move.
budget: ~1,600 lines written, ~1,700 removed, ~37 agent-hours over three slices (Cost). The rulings estimated
~1,300 lines and 23 to 27 hours; the difference is the bridge's sites, counted below, and the forwarding test.

## Summary for the maintainer

Today a sampler's kind of response is tested by name 129 times in the engine and the R and C bridge, and nine of those
tests ("gaussian or aft") stand for three different facts. After this work each of the eight kinds (gaussian,
Student-t, probit, logistic, aft, ordinal, negative binomial, multinomial) states its fixed facts once, in a constant
record written beside its class, and a table gathers the eight so the build fails when one is missing. The engine,
the factories and the bridge read the record; 84 of the 129 tests become reads of a named fact, 36 remain in switches
that list every kind because their arms really differ (building the object, or R's wording), 8 parse what the user
asked for, and 1 is family behavior kept as it is.

Then the class the bridge holds shrinks. Its 85 virtual functions become a core of about 60 that every sampler
answers, plus ten small interfaces for what only some samplers have (a residual scale to set, a variance forest, a
basis, a count matrix, ...). The bridge asks a sampler for its capability object and gets each interface or nothing;
having the interface is the capability, so the bridge's check and the engine's refusal cannot disagree. The
forwarding class, which repeats every call and the default arguments of 13 virtual functions, is deleted.

Nothing a user or a package linking to dbarts.h sees changes, and no draw moves. Two calls are yours (Open calls):
what exactly the capability object is, and whether the fold may come before the shrink.

## Goal

After the three slices land:

- `ResponseKind` names eight kinds; each response class carries one `constexpr` record of its fixed facts beside its
  code, gathered into one table whose static checks fail the build when a kind lacks a row, a row sits at the wrong
  index, or a row contradicts its class (an optional step flagged without its interface, or the reverse).
- The engine, the factories and the bridge read the record wherever they now test a family's name for a fact; the
  remaining switches over the kind list every enumerator with no default arm. The engine's multinomial path no
  longer calls itself logistic, and no test of a Student-t sampler passes through "gaussian".
- `Sampler<L, ResidT>` is the final implementation of `SamplerBase`; `SamplerFacade` and `makeIngestingFacade` are
  gone; no virtual function carries a default argument.
- `SamplerBase` is a core plus one query, `capabilities()`, whose answer holds a typed pointer for each optional
  function set, null where the sampler has none. Presence is fixed for the sampler's life.
- A tests/cpp check proves each class's virtual answers match its row; a second proves each fixture's capability
  object matches its record and shape.

## Context

- Rulings: dec-B763 (option (b), the generics memo's shape, "Oh, wait, do (b)"); dec-B767 (option (a), "OK, let's do
  it then"); dec-B768 (no latent-kind field in the C facts struct at 1.0); dec-B615 and dec-B639 (the facts struct,
  built in host-neutral S3); dec-B92 (leaf compile-time, family run-time; the axes doubted, not rebuilt); dec-B85 (the
  host-neutral work precedes 1.0); dec-B764 (chain-split's acceptance bar). TODO items family-fact-records,
  sampler-core-and-capabilities, capi-sampler-info, chain-split.
- The memo: [4. How the current shape lines up](../design/engine-generics-review.md#4-how-the-current-shape-lines-up),
  "shrinking that class to a common core plus a queried capability object touches no sampling code at all".
- The boundary today: [`SamplerBase`](../../src/bartcore/facade.hpp), [`SamplerFacade`](../../src/bartcore/facade.hpp),
  [`SamplerShape`](../../src/bartcore/facade.hpp), [`makeIngestingFacade`](../../src/bartcore/facade.hpp); its
  conformance gate [`runFacadeTests`](../../tests/cpp/test_facade.cpp) and
  [`checkShapeMatchesImpl`](../../tests/cpp/test_shape.cpp). The bridge fetches `shape()` at 106 sites (83 in
  R_interface_bartcore.cpp, 23 in C_interface.cpp); the family and global-block fields SC1 replaces are read at 28 of
  them, the presence fields SC3 moves at 29.
- The facts today: [`ResponseFamily`](../../src/bartcore/model.hpp) has six values for eight `final` classes under
  [`ResponseModel`](../../src/bartcore/model.hpp). [`TResponse`](../../src/bartcore/model.hpp) is built from the
  gaussian family plus a finite residual df; [`MultinomialResponse`](../../src/bartcore/model.hpp) is built under
  `ResponseFamily::logistic` "and family() is not read on this path". Class facts are already virtual queries
  ([`supportsActiveRows`](../../src/bartcore/model.hpp), true in all eight classes and false only in the base;
  [`workingWeightsVaryPerSweep`](../../src/bartcore/model.hpp); [`carriesOrdinalThresholds`](../../src/bartcore/model.hpp);
  [`scaleLatents`](../../src/bartcore/model.hpp), whose refusing default no caller reaches), and four unit-scale
  classes repeat the same four overrides (`initialSigma`, `fitScale`, `fitShift`, `sigmaScale` returning 1, 1, 0, 1).
- Every name test, line by line on main 0580ff13a and on collapsed-move-rebuild 3b304fdff, classified:
  scratch/plan-core/name-tests.md. The by-symbol table is [Appendix: the name tests](#appendix-the-name-tests).
- Branches that bear: collapsed-move-rebuild (worktree kmp-step0, being built) adds a logistic name test in its new
  `collapsedRescaleApplies` and a NaN-default virtual
  (`logLikelihoodAtScaledFits`, "which no caller reaches"); chain.hpp still holds 45 references. This plan is written
  against main plus that branch. host-neutral-engine (docs/plans/host-neutral-engine.md on that branch) rewrites the
  bridge's creation and refusal code in S2a and S2b and builds the facts struct in S3. apx-size edits model.hpp,
  combiner.hpp and the bridge.
- 0.9-34, run (scratch/famtraits/probe-cpp-09.R against scratch/libs/old-ecosystem, rerun 2026-11-10): a probit
  sampler accepted `setControl` with `binary = FALSE`, still read binary afterward and kept sigma at 1, 1, 1. The
  family was fixed at creation in practice; a record of facts fixed at creation is no regression. 0.9-34 had no
  type-erased boundary and no C header of this shape.

## Constraints

- Draws identical, bit for bit, in every slice. The latent-scale anchor, the forest sd unit and the basis-scale
  expressions move verbatim: a record stores an enumerator, never a precomputed constant (the logistic anchor stays
  `std::numbers::pi / std::sqrt(3.0)` at its reader; the bridge's logistic leaf scale stays
  `std::numbers::pi * std::sqrt(3.0)`).
- Messages identical: every tinytest message pin unchanged, the order of the bridge's refusal cascades included.
- The record is internal. It never appears in dbarts.h (its rows would become ABI); the `latent` column is the
  engine's and is not published (dec-B768). dbarts.h is untouched and `tools/check-api-hash.sh` passes with
  `DBARTS_C_API_HASH` unchanged.
- The engine stays R-agnostic: the record holds no R wording; the bridge keeps its words and keys them by kind.
- Out of scope: the R layer's own "gaussian or aft" tests (R/generics.R, R/spec.R, R/bart.R), which read the record
  after 1.0 through an internal call if ever (the ruling names the engine, the bridge and the factories); the forest
  couplings' base class [`ForestCombiner`](../../src/bartcore/combiner.hpp) and its capability virtuals (chain-split's
  territory, below); the factory's half-built-sampler and shared-tail findings of the C++ review's slice 5 (not ruled);
  the header text on family growth, which is host-neutral S3's (TODO capi-sampler-info); the aft fallback question
  (Not verified), whose fix would move draws.

## Design

### The record

```cpp
enum class ResponseKind : std::uint8_t {
  gaussian, studentT, probit, logistic, aft, ordinal, nbinom, multinomial
};
struct ResponseFacts {
  ResponseKind kind;
  bool drawsResidualScale;       // a sigma the sampler draws or a caller sets
  bool precisionChannelFree;     // the per-row precision channel is free for a variance forest
  bool carriesCensoring;         // a per-row event status (aft)
  bool takesTrials;              // a binomial response's trials may be given (logistic)
  bool carriesResidualDf;        // the global blocks: Student-t's nu,
  bool carriesShape;             //   negative binomial's r,
  bool carriesOrdinalThresholds; //   ordinal's cutpoints
  bool admittedByAmplitudeCoupling;
  bool admitsSinglePrecision;    // the fp32 running residual
  bool rescalesLatents;          // optional step built: the probit rescaling step
  bool scoresScaledFits;         // optional step built: the collapsed scale move
  ResponseUnits units;           // none, range, countShift
  LatentScale latentScale;       // fromResponse, standardNormal, logisticDistribution, notStated
  LatentKind latent;             // none, location, precision
  WeightsMeaning weights;        // precisionOrFrequency, membership, counts, none
  ResponseSupport support;       // real, binary, categoryIndex, count, categoryCounts
  DefaultNodeScale defaultNodeScale; // responseUnits, unitLatent, logisticLatent, pairwiseLogistic
  LeafRowsRead leafRowsRead;     // none, responseLessOffset, response
};
```

The names are the implementer's to refine; the columns are the facts the three expert reviews named plus those the
sites below read. Each class declares its row as a `static constexpr ResponseFacts` member beside its code; a virtual
`facts()` returns it; `responseFacts(ResponseKind)` reads the table for code that has no object yet (the factories,
the bridge). The table is built from the list of eight classes, and its static checks are: one row per kind; row i
holds kind i; `rescalesLatents` and `scoresScaledFits` equal `std::is_base_of_v` of the step's interface (Optional
steps); and the model invariants `carriesResidualDf` implies `latent == precision`, `precisionChannelFree` implies
`latent != precision`, `drawsResidualScale` implies `units == range`, `admittedByAmplitudeCoupling` implies
`latentScale != notStated`.

The rows, each value checked against the site that reads it today:

| column | gaussian | studentT | probit | logistic | aft | ordinal | nbinom | multinomial |
|---|---|---|---|---|---|---|---|---|
| drawsResidualScale | yes | yes | | | yes | | | |
| precisionChannelFree | yes | | | | yes | | | |
| carriesCensoring | | | | | yes | | | |
| takesTrials | | | | yes | | | | |
| carriesResidualDf | | yes | | | | | | |
| carriesShape | | | | | | | yes | |
| carriesOrdinalThresholds | | | | | | yes | | |
| admittedByAmplitudeCoupling | yes | | yes | yes | | | | |
| admitsSinglePrecision | yes | yes | | | | | | |
| rescalesLatents | | | yes | | | | | |
| scoresScaledFits | | | | yes | | | | |
| units | range | range | none | none | range | none | countShift | none |
| latentScale | fromResponse | fromResponse | standardNormal | logisticDistribution | fromResponse | standardNormal | notStated | notStated |
| latent | none | precision | location | precision | location | location | precision | none |
| weights | precisionOrFrequency | precisionOrFrequency | membership | counts | none | membership | counts | none |
| support | real | real | binary | binary | real | categoryIndex | count | categoryCounts |
| defaultNodeScale | responseUnits | responseUnits | unitLatent | logisticLatent | responseUnits | unitLatent | unitLatent | pairwiseLogistic |
| leafRowsRead | responseLessOffset | responseLessOffset | none | none | none | none | response | none |

Notes on the rows:

- Student-t follows every "gaussian" site today, so each column keeps what that site does for it now (it draws sigma,
  carries units, admits single precision, reads its weights as gaussian's, reads the response for its leaf rows),
  except where the code already treats it apart: the variance-forest refusal's finite-df term becomes
  `precisionChannelFree` false, the df block becomes `carriesResidualDf`, and amplitude admission becomes false. That
  last one is a new engine refusal (the amplitude factory refuses gaussian with a finite residual df, which the
  amplitude chain silently ignores today); the bridge refuses the combination first with "Student-t residuals", so no
  message moves (Calls made, 6).
- `latentScale` is read only by the calibration map, which only admitted kinds reach. Ordinal's value is the model's
  (a unit-variance latent, as its default node scale of 3 already says), not today's unreachable empirical arm;
  nbinom and multinomial state none, and admitting either later fails the build until it does (Calls made, 7).
- `defaultNodeScale` mirrors R's `defaultLeafScale` in R/model.R, multinomial's pairwise value included, which no C++
  reader uses until R reads the record.
- aft's `leafRowsRead` is none, as [`withResponseFallback`](../../src/bartcore/sampler.hpp) has it, although
  [`carriesUnits`](../../src/bartcore/chain.hpp) counts aft as carrying units: the record keeps today's draws and
  the question goes to Not verified.
- `latent` is what `latents()` returns: multinomial's Polya-Gamma variates live in its coupling, so its row says
  none. `workingWeightsVaryPerSweep` is `latent == precision` (logistic, Student-t, nbinom, the three overriding it
  today).

### The kind and the request

`ResponseFamily` stays the factories' request type (six values; Student-t is gaussian with a finite residual df).
One exhaustive function, `responseKindFor(ResponseFamily, double residualDf)`, holds the rule "a finite residual df
selects Student-t", which today lives in the chain's construction switch and the factory's variance-forest refusal;
the multinomial factory names its kind directly. The 684 tests/cpp lines that build samplers by family keep working
(Calls made, 4). `Chain` drops its `family_` member and holds a pointer to its row, set once at construction;
`SamplerShape` gains `responseKind` and loses `family`, `carriesShape` and `carriesResidualDf`, which the row answers.
The C bridge maps the kind to `dbarts_family` in one exhaustive switch (Student-t to GAUSSIAN, as dec-B615 has it;
multinomial to MULTINOMIAL), replacing the counts-mutation fingerprint in [`dbarts_sampler_family`](../../src/C_interface.cpp).

### Optional steps

Each optional move becomes an interface its class implements: `LatentRescaling` (the probit rescaling step, now
[`scaleLatents`](../../src/bartcore/model.hpp) with a refusing default) and `ScaledFitScoring` (the collapsed scale
move's evaluator, now a NaN default on the rebuild branch). `ProbitResponse` and `LogisticResponse` inherit them; the
chain takes the typed pointer once, at construction, from the class it builds. Each move keeps one eligibility
predicate evaluated every sweep, unchanged in shape: [`forestRescaleApplies`](../../src/bartcore/chain.hpp) reads
"the step is built" (the pointer) where it reads the probit name, and `collapsedRescaleApplies` likewise for logistic,
each still conjoined with the scale-orbit conditions, the toggle and the finite k-prior scale it has now. No class
carries an evaluator no move calls (Calls made, 3).

### The core and the capability object

A function joins a facet when its whole refusal condition is a creation fact; the facet's presence is that fact.
A function whose answer off a shape is a true value (a count of 0, an empty tally, a null `latents()`) stays in the
core. A refusal narrower than its facet (an index naming no forest, a non-map forest, no test rows) stays a returned
refusal inside the facet. By that rule the expected partition, which the implementer confirms against each
function's documented refusal:

| facet | functions | present when |
|---|---|---|
| `ResidualScaleControl` | setSigma, setResidualPrior | row.drawsResidualScale and no variance forest |
| `CensoringControl` | setSurvivalStatus | row.carriesCensoring |
| `ResidualDfControl` | setResidualDf | row.carriesResidualDf |
| `ShapeControl` | setNegativeBinomialShape, shapeParameter | row.carriesShape |
| `BinomialCountsControl` | setBinomialCounts | built with trials (row.takesTrials) |
| `VarianceForestSurface` | predictVariance, currentVarianceFits, sampleVarianceForestFromPrior | a variance forest is built |
| `AmplitudeCoupling` | numForestAmplitudes, amplitudes, setForestBasis, setForestBasisRescaled, predictPerForest, setForestWeights, setForestMapSd, setForestCoefficientPrior | the amplitude coupling |
| `CountsCoupling` | setCounts, setCategoryOffset, setCategoryTestOffset | the counts (softmax) coupling |
| `AdditiveLocation` | fitsWithoutOffset | every coupling but the counts coupling |
| `LeafCovariateControl` | rederiveStandardization, setLeafLengthscales (GP's own refusal kept) | a linear or GP leaf |

That moves 24 virtuals off the core; `totalAmplitudes` goes (it duplicates `SamplerShape::numAmplitudes`) and
`capabilities()` comes, so the core holds about 60. `SamplerShape` loses `hasVarianceForest`, `holdsTrials`,
`supportsForestWeights` (unless its predicate proves narrower than the coupling, when it stays a value),
and `supportsCountsMutation`. (`supportsActiveRows`, true for every class, goes in SC1 with its virtual and the
bridge's unreachable refusal, Calls made, 5.)

The capability object, `SamplerCapabilities`, is a trivially destructible struct of facet pointers, filled once when
the sampler is built (every presence is a creation fact; no setter moves one) and returned by `capabilities()`. The
facets are interfaces with protected non-virtual destructors that `Sampler<L, ResidT>` implements; a const sampler's
query hands facets whose mutators it cannot call (two views of the struct, or one templated on constness: the
implementer's choice, the reviewer checks there is no const hole). The bridge asks:

```cpp
if (auto* control = sampler.capabilities().residualScale)
  control->setSigma(value);
else
  refusePinnedSigmaChange(sampler, "$setSigma");
```

How each question is asked afterward: a count, a value or a setting that a setter moves comes from `shape()`, fetched
once per entry; a fixed fact of the kind from `responseFacts(shape.responseKind)`; a function set from
`capabilities()`. The C header's capability statuses (dbarts_sampler_setSigma's UNSUPPORTED, getLatents' 0) read the
facet and the row; dbarts.h's "test it once at setup" holds by construction, since presence never changes.

### The fold

`Sampler<L, ResidT>` derives `final` from `SamplerBase`, each method that implements the interface marked
`override`; the shape fill moves from the forwarder into `Sampler::shape()`; `Sampler::shape(chainNum)`, which
collides with it, becomes `shapeParameter`, the interface's name; the factories return `std::make_unique<Sampler<L,
ResidT>>` (the ingestion check kept, `makeIngestingFacade` renamed). Virtuals take no default arguments; the short
spellings callers use become non-virtual overloads on the base, re-exposed in `Sampler` by using-declarations as the
dense convenience spellings already are (Calls made, 9). Calls inside `Sampler` devirtualize; every bridge call is
still one virtual hop, as today.

### Relation to chain-split (dec-B764)

Neither plan blocks the other; both cite by symbol, and whichever lands second rebases.

- SC1 touches `Chain` at its constructors, the 33 fact sites and the `family_` member; chain-split's per-concern
  headers move those definitions. The conflict surface is the three constructors, `setResidualPrior`,
  `carriesUnits`, the two rescale predicates, `buildVarianceForest`'s gate and the anchor trio
  ([`latentScaleAnchor`](../../src/bartcore/chain.hpp), [`derivedForestSdUnit`](../../src/bartcore/chain.hpp),
  [`rowsForBasisScale`](../../src/bartcore/chain.hpp)).
- The record serves chain-split's acceptance bar: the calibration-map class can take the row's `latentScale` and the
  response's scale accessors as inputs, instead of reading `Chain`'s family, so it reaches nothing back in `Chain`.
- The forest-structure base is untouched here. `VarianceForestSurface` is implemented by `Sampler` over the chain's
  variance-forest functions; if chain-split renames them, the facet's implementation follows.
- The couplings' base class keeps its capability virtuals here. SC3 reads them once, at construction, to fill the
  capability object, so a later change of that class (a coupling kind in place of its booleans) touches one reader.
- Expected order: SC1 to SC3 before chain-split, by the rulings' timing (before host-neutral S3); chain-split waits for
  the move-settings and replay slices anyway.

### Relation to the host-neutral engine

- SC1's bridge half best lands before S2a, so S2a's R-free validators (`enforceBinaryWeightPolicy`,
  `validateResponseSupport`, `requireResolvedSigmaEstimate`, `sigmaIsPinned`) read the record from their first
  version; if S2a goes first, SC1 converts the same functions where S2a leaves them.
- S2b's `buildSampler` calls the factories, whose signatures SC1 to SC3 keep (`ResponseFamily` in, `SamplerBase`
  out). S2b may later pass a `ResponseKind` straight from `dbarts_family` (STUDENT maps one to one).
- S3 reads the new shape: `usesVarianceForest` from `capabilities().varianceForest`, `usesStudentResidual` from the
  row's `carriesResidualDf`, `family` from the kind map. Its plan's facts-struct section names
  `SamplerShape::carriesResidualDf` and `SamplerShape::hasVarianceForest`; those two sentences change at its next
  amendment. The header sentence on family growth (TODO capi-sampler-info) is S3's.

## Slices

Each is its own engine feature branch, serial behind collapsed-move-rebuild, landed on its own after the reviewer's
pass and its gates. Tier: changes draws if wrong, so the engine tier of [Process by risk](README.md#process-by-risk)
in each, an opus review that runs the mutants.

### SC1: the records (engine, factories, bridge)

Steps:

1. Recount the name tests on the slice's base by the rule in scratch/plan-core/name-tests.md; any test a landing
   added since (k-move-pricing, approximate-draws) joins the slice.
2. model.hpp: `ResponseKind`, the column enumerations, `ResponseFacts`, a row in each of the eight classes, the
   table and its static checks, `responseFacts`, `responseKindFor`; the optional-step interfaces; `facts()`.
   The fact virtuals go (`supportsActiveRows`, `workingWeightsVaryPerSweep`, the three `carries*` queries become
   non-virtual reads of `facts()`), and the four unit-scale classes drop their four scale overrides for base
   defaults of 1, 1, 0, 1 (Calls made, 8). Remove the stale comment "only the binary families fix it".
3. chain.hpp, facade.hpp, sampler.hpp: the 43 F sites read the row; the two construction switches key on the kind
   (Student-t its own arm in the single-forest one); `family_` goes from `Chain`; the multinomial chain and sampler
   carry kind multinomial.
4. The bridge: the 41 F sites read `responseFacts(kind)`; `refusedAmplitudeFamilyReason`,
   `enforceBinaryWeightPolicy` and `supportFamily` key on the kind with their wording kept, and `supportFamily`
   maps Student-t to its own kind (its "documentary" gaussian arm goes); [`refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp)
   keeps its cascade order, the family reason read from the request's base kind so "Student-t residuals" still
   fires where it does; `dbarts_sampler_family` maps the kind. `SamplerShape` gains `responseKind`, loses `family`,
   `carriesShape`, `carriesResidualDf` and `supportsActiveRows` (with the bridge's unreachable refusal). Wherever a
   switch gains a Student-t arm, that arm does what the gaussian arm does today.
5. Before each bridge site is converted, its implementer names the tinytest pin that would fail if the site read
   the wrong column, in a table in the landing note; where none exists (the unresolved-sigma refusal of a flat
   creation and the single-precision refusal looked thin), a pin joins the slice.
6. tests/cpp, a new suite `test_response_facts.cpp` (added to the Makefile's sources and main.cpp): for each of
   the eight classes, built on a small fixture, `facts()` is its table row and its virtual answers match the row: `initialSigma`, `fitScale`, `fitShift`,
   `sigmaScale` and `getScale` against `units` and `drawsResidualScale`; `latents()` non-null exactly where `latent`
   is not none, after a refresh; working weights changing across `refreshLatents` exactly where `latent` is
   precision; the global blocks present exactly where their column says; `setBinomialCounts` accepted exactly on a
   trials-built logistic; `setSurvivalStatus` reaching the class exactly on aft. Plus a bitwise pin of the logistic
   anchor at `std::numbers::pi / std::sqrt(3.0)` through the test peer.

Size: ~300 engine lines written and ~120 removed, ~180 bridge lines touched, ~260 tests/cpp. Cost: agent ~8, review
~4, gates ~3.

### SC2: the fold

Steps: `Sampler` final over `SamplerBase` as in [The fold](#the-fold); `SamplerFacade` deleted; the factories build
`Sampler` directly; default arguments off every virtual, the short forms as non-virtual base overloads;
test_facade.cpp's forwarding rows go (each compared a call through the base with the same call on `Sampler`, which
`override` now makes one function); the rows that check engine behavior no other suite checks (the partly filled
and wrapped saved-tree store, the two-chain and two-forest confusables) move to test_sampler.cpp, listed in the
landing note; test_shape.cpp reads `Sampler::shape()`. A clang-tidy run with `google-default-arguments` over
facade.hpp and sampler.hpp, clean.

Size: ~250 lines written; ~400 engine and ~1,100 tests/cpp removed. Cost: agent ~5, review ~2.5, gates ~2.

### SC3: the core and the capability object

Steps: the facet interfaces and `SamplerCapabilities` as in [The core and the capability object](#the-core-and-the-capability-object);
`Sampler` implements the facets and fills the object at construction from the row, the coupling's existing
predicates, the leaf type and the variance forest; the 24 functions leave `SamplerBase`; `totalAmplitudes` and the
presence fields go; the bridge's 29 presence reads and about 35 facet calls ask `capabilities()`, the C entries'
capability statuses among them; tests/cpp `test_capabilities.cpp`, built on test_facade.cpp's fixtures (gaussian,
probit, logistic, binomial, Student-t, aft, ordinal, nbinom, BCF, multinomial, variance forest, linear, GP,
monotone): each fixture's object has exactly the facets its record and shape imply, and each facet's functions
answer through the facet as through `Sampler`.

Size: ~600 lines touched (facade.hpp ~150, sampler.hpp ~120, bridge ~180, tests/cpp ~250, net ~100 removed).
Cost: agent ~7, review ~3, gates ~2.

## Verification

Per slice, against the slice's own private library, every install `--preclean` (vtables and class layouts change):

    R CMD INSTALL --preclean -l <lib> .                       # buildInfo()$mode the shipped build
    cd tests/cpp && make && ./test_bartcore                   # plain
    make clean && make OPT="-O2 -g -fsanitize=address,undefined" && ASAN_OPTIONS=detect_container_overflow=0 ./test_bartcore
    R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts", lib.loc = "<lib>")'   # one process, default ncpu
    R CMD INSTALL --preclean --configure-args=--enable-reference-build -l <reflib> .
    R_LIBS=<reflib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/<MANIFEST current> --bitwise --strict-coverage
    R_LIBS=<reflib> Rscript benchmarks/R/bcf-equivalence.R compare benchmarks/baselines/<MANIFEST current> --bitwise
    R_LIBS=<reflib> Rscript benchmarks/R/multinomial-equivalence.R compare benchmarks/baselines/<MANIFEST current> --bitwise
    # the four test-reproducibility-*.R files on the reference build, run_test_file each
    Rscript -e 'lintr::lint_package()' && air format --check . && Rscript tools/check-rc-codoc.R . && \
      Rscript tools/check-win-drift.R . && Rscript tools/check-doc-freshness.R .
    Rscript benchmarks/R/mutation-battery.R verify-anchors
    sh tools/check-api-hash.sh                                # DBARTS_C_API_HASH unchanged
    R CMD check --as-cran <tarball from a clean copy staged outside the tree>

Expected: every tests/cpp suite passes plain and under ASan with no diagnostic; the suite passes with every message
pin unchanged; each trio compare prints its full count of "identical draws (same RNG stream)" lines and no "max |z|"
line (Gate hygiene); the four snapshots pass. SC1 and SC3, which edit the bridge's validators and capability checks,
add the R-loaded ASan run over the creation, mutation and test-capi.R suites; SC2, which changes ownership in the
factories, adds the x86 leg with LeakSanitizer. The orchestrator's own check, each slice: `bench-sampler.R record`
on the base and `compare` on the tip, same quiet machine, interleaved; expected neutral (the record is read a handful
of times per sweep per chain, never per row).

## Mutants

The reviewer runs each on the slice's tree, rebuilds `--preclean`, and records which gate failed; a mutant no gate
catches is a finding, answered with a pin. `touch` each reverted file.

SC1:
1. Student-t's `precisionChannelFree` true: test_response_facts fails, and a variance forest with a Student-t
   residual is no longer refused (test-heteroscedastic.R's "variance forest does not support" pin).
2. Student-t's `drawsResidualScale` false: sigma pinned on a Student-t fit; equivalence.R's Student-t scenario
   fails bitwise (no arm64 snapshot holds a Student-t fit).
3. `responseKindFor` ignores the residual df: every Student-t fit becomes gaussian; equivalence.R's Student-t
   scenario and the tinytest files that read a df channel fail.
4. The logistic anchor written as `std::numbers::pi * (1.0 / std::sqrt(3.0))`: the bitwise anchor pin fails; record
   whether bcf-equivalence's logistic scenario also catches it.
5. nbinom's `units` range: a constant count response's transform is accepted where it is refused; the state and
   anchor tests fail.
6. Two rows swapped in the table: the build fails on the index check.
7. `scoresScaledFits` true on probit, or `ScaledFitScoring` dropped from `LogisticResponse`: the build fails on
   the interface check.
8. The bridge's `refusedAmplitudeFamilyReason` read from the Student-t row: test-bcf-creation.R's
   "Student-t residuals" pin fails on message order.
9. `requireResolvedSigmaEstimate` read for aft as not drawing sigma: an aft flat creation with an NA sigma is
   accepted; the pin step 5 names or adds fails.
10. aft's `leafRowsRead` set to responseLessOffset (the suspected fix): expected to pass every gate, which measures
    the coverage gap of Not verified; recorded, not answered here.

SC2:
11. A default argument re-added to a virtual with a value different from the base overload's: clang-tidy flags it.
12. `Sampler::shape()` fills `numSavedDraws` from the store's capacity: test_shape fails.
13. `createSampler` builds the GP instantiation for a linear designation: the leaf tests fail.

SC3:
14. `varianceForest` present on a homoscedastic gaussian: test_capabilities fails; predictVariance reaches a sampler
    with no variance forest.
15. `residualScale` present under a variance forest: the setSigma capability arms of test-capi.R and the tinytest
    pinned-sigma refusal on a heteroscedastic sampler fail (the reviewer confirms both arms exist, else adds them).
16. `countsCoupling` absent on multinomial: `dbarts_sampler_family` and the R count setters' pins fail.
17. `binomialCounts` keyed on the row alone (present on a binary logistic): test-count-response-doors.R's refusal
    fails.

## Open calls

1. **What the capability object is.** (a) One object built at creation, holding a typed pointer per optional
   function set, null where absent, read by one query; presence leaves `SamplerShape`. Recommended: it is the
   memo's "queried capability object", the answer and the refusal are one fact, and the C header's "test once"
   holds by construction. (b) One query per facet on `SamplerBase`, const and non-const, about 20 virtuals in
   place of one. (c) One second interface holding every non-core function, each still refusing: the class is
   split but the union is not, so it does not meet the ruling. Cost of (a) as planned; (b) about the same.
2. **The fold before the shrink.** The ruling's words are "shrink ... then fold". Folding first (SC2 before SC3)
   reaches the same end state without rewriting the forwarding class and its test table into facets only to delete
   them, about 2 agent-hours saved. Recommended: fold first. The alternative is SC3 before SC2, the facets first
   implemented by `SamplerFacade`.

## Calls made

1. Three slices, the records first (engine and bridge together), as the expert synthesis advised ("Q3 first, so the
   capability object is filled from the record"). Alternative: four, the bridge's half its own slice (one more gate
   battery and review, about 4 hours).
2. Every family name test classified F, S, P or K, and only F replaced; exhaustive switches whose arms build objects
   or carry R's wording stay, keyed by the kind. Alternative: a name column in the record and the wording composed
   from it, which puts text in the engine's record.
3. Optional steps as interfaces tied to their flags by `std::is_base_of_v` (the MCMC review's sketch). Alternatives:
   a flag plus a refusing default virtual (the default is dead code); a required evaluator on every class (writes
   evaluators no move calls).
4. `ResponseFamily` stays the factories' request type, the kind derived by one function. Alternative: the factories
   take `ResponseKind`, which touches 684 tests/cpp lines and the bridge's creation routes S2b rewrites anyway.
5. `supportsActiveRows` deleted, every class answering true. Alternative: a record column that is true in every row.
6. Student-t not admitted by the amplitude coupling, the engine refusing what it now ignores. Alternative: admit it
   and keep building a gaussian response, dropping the residual df in silence.
7. `latentScale` states the model's value for ordinal and `notStated` for nbinom and multinomial, all unreachable.
   Alternative: today's empirical arm for every non-binary row, which states a scale no model gives ordinal.
8. The four scale accessors default to the unit scale in the base, the tests/cpp row check catching a class that
   should override. Alternative: keep them pure, the four unit-scale classes repeating them.
9. No default arguments on virtuals; short forms as non-virtual base overloads. Alternative: the non-virtual
   interface idiom, a public non-virtual function over a private virtual for every entry.
10. The bridge's "logistic restates its weights after setData" stays a kind test: it is family behavior, and the
    nearest fact (counts and precision latents) also covers nbinom, whose setData does not need it.
11. The forwarding test's rows deleted, not kept as smoke tests, the behavior rows moved.
12. The couplings' base class and its capability virtuals left to chain-split.

## Not verified

- aft and the leaf-row fallback: [`carriesUnits`](../../src/bartcore/chain.hpp) counts aft as carrying units while
  [`withResponseFallback`](../../src/bartcore/sampler.hpp) does not read its response. An aft fit with a linear or
  GP leaf, an active-row mask and one distinct response value among the rows in may read leaf rows other than its
  transform's. It needs a constructed fit to show, which this planner could not run. SC1 keeps today's draws; a
  probe belongs before SC1 starts, and a fix, if one is owed, is its own draw-changing slice.
- The facet partition: each function's refusal condition was read from its documentation, not traced; the
  implementer confirms it and the critique should check the table.
- Whether bcf-equivalence's scenarios reach the logistic anchor (mutant 4 says).

## Cost

| slice | lines | agent-hours |
|---|---|---|
| SC1 records | ~740 written, ~120 removed | ~15 |
| SC2 fold | ~250 written, ~1,500 removed | ~9.5 |
| SC3 core and capabilities | ~600 touched, ~100 net removed | ~12 |
| total | ~1,600 written, ~1,700 removed | ~37 |

Plans have run 1.5 to 2 times low.

## Appendix: the name tests

Counted by `grep -c 'ResponseFamily::\|RF::'` on main 0580ff13a, chain.hpp on collapsed-move-rebuild 3b304fdff:
chain.hpp 45, facade.hpp 7, sampler.hpp 3, combiner.hpp 1, R_interface_bartcore.cpp 67, C_interface.cpp 6; 129.

| where | references | becomes |
|---|---|---|
| [`Chain`](../../src/bartcore/chain.hpp) single-forest and amplitude constructors, the construction switches | 12 | switches over the kind |
| `sigmaIsFixed_` at both constructors, [`setResidualPrior`](../../src/bartcore/chain.hpp) | 6 | drawsResidualScale |
| [`buildVarianceForest`](../../src/bartcore/chain.hpp)'s gate, [`varianceForestIsRefused`](../../src/bartcore/facade.hpp) | 4 | precisionChannelFree |
| [`carriesUnits`](../../src/bartcore/chain.hpp) | 3 | units |
| [`latentScaleAnchor`](../../src/bartcore/chain.hpp), [`derivedForestSdUnit`](../../src/bartcore/chain.hpp), [`rowsForBasisScale`](../../src/bartcore/chain.hpp) | 18 | latentScale |
| [`forestRescaleApplies`](../../src/bartcore/chain.hpp), `collapsedRescaleApplies` (rebuild) | 2 | the optional-step interfaces |
| the multinomial chain and sampler's logistic, `Chain`'s `family_` member | 3 | the multinomial row, the row pointer |
| [`createConstantLeafSampler`](../../src/bartcore/facade.hpp) | 1 | admitsSinglePrecision |
| [`trialsAreRefused`](../../src/bartcore/facade.hpp) | 1 | takesTrials |
| [`createAmplitudeSampler`](../../src/bartcore/facade.hpp) | 3 | admittedByAmplitudeCoupling |
| [`withResponseFallback`](../../src/bartcore/sampler.hpp) | 2 | leafRowsRead |
| [`AmplitudeSpec`](../../src/bartcore/combiner.hpp)'s family default | 1 | kept (a request) |
| [`resolveFamily`](../../src/R_interface_bartcore.cpp), the resid.df request check | 7 | kept (parsing) |
| [`applySurvivalStatus`](../../src/R_interface_bartcore.cpp), the setResponse status and setData refusals | 3 | carriesCensoring |
| [`refusedAmplitudeFamilyReason`](../../src/R_interface_bartcore.cpp) | 6 | admittedByAmplitudeCoupling, wording by kind |
| [`defaultLeafScale`](../../src/R_interface_bartcore.cpp) | 6 | defaultNodeScale |
| the creation weights reading, the aft weights refusal, [`familyCarriesNoWeights`](../../src/R_interface_bartcore.cpp) | 4 | weights |
| [`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp) | 5 | a switch over the kind, wording kept |
| the nbinom shape and weighted-count caps at creation, setResponse, setData, setWeights | 5 | carriesShape |
| [`requireResolvedSigmaEstimate`](../../src/R_interface_bartcore.cpp), [`sigmaIsPinned`](../../src/R_interface_bartcore.cpp), [`fillResidualParameters`](../../src/R_interface_bartcore.cpp) | 6 | drawsResidualScale |
| [`familyHasResponseScale`](../../src/R_interface_bartcore.cpp), the count family's supplied and stored ranges | 8 | units |
| [`validateResponseSupport`](../../src/R_interface_bartcore.cpp), setData's counts on a binary logistic | 8 | support, takesTrials |
| the single-precision refusal | 1 | admitsSinglePrecision |
| setData's logistic weights restated | 1 | kept (family behavior) |
| [`supportFamily`](../../src/R_interface_bartcore.cpp) | 7 | a map to the kind |
| [`dbarts_sampler_family`](../../src/C_interface.cpp) | 6 | a map from the kind |

Replaced by record reads 84; kept in switches over the kind 36; parsing and request defaults 8; family behavior 1.
