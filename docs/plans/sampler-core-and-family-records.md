# sampler-core-and-family-records: one constant record per kind of response, and a sampler interface of a core plus a queried capability object

Status: PLANNED 2026-11-10, amended the same day from the blind critique (verdict AMEND FIRST,
scratch/plan-core-crit/critique.md, A1 to A7 taken); not built. Two TODO items planned together, as their rulings
ask: family-fact-records (dec-B767) and sampler-core-and-capabilities (dec-B763), plus one small fix the critique's
probe found (SC0). Awaits collapsed-move-rebuild's landing.

agent: opus implementer, one per slice; opus reviewer running the mutants below, effort high.
rng: NEUTRAL in SC1 to SC3: no draw, default, message or fit object moves; the equivalence trio is bitwise. SC0 is
POSTERIOR-CHANGING in one corner only (aft with a linear or gp leaf and a creation mask keeping two or more rows that
hold one distinct log-time): the leaf covariates' standardization, and so the slopes' prior, moves there, and the
creation warning is raised; with one row in only the warning's words change. The trio and the snapshots reach
neither.
window: pre-release, before the host-neutral engine's S3 (the C creation entry and the facts struct), and best
before its S2a (below). dbarts.h does not change; its hash does not move.
budget: ~1,680 lines written, ~1,700 removed, ~42.5 agent-hours over four slices (Cost). The rulings estimated
~1,300 lines and 23 to 27 hours; the difference is the bridge's sites, counted below, the forwarding test and SC0.

## Summary for the maintainer

Today a sampler's kind of response is tested by name 129 times in the engine and the R and C bridge, and nine of those
tests ("gaussian or aft") stand for three different facts. After this work each of the eight kinds (gaussian,
Student-t, probit, logistic, aft, ordinal, negative binomial, multinomial) states its fixed facts once, in a constant
record written beside its class, and a table gathers the eight so the build fails when one is missing. The engine,
the factories and the bridge read the record; 78 of the 129 tests become reads of a named fact, 42 remain in switches
that list every kind because their arms really differ (building the object, R's wording, or R's defaults), 8 parse
what the user asked for, and 1 is family behavior kept as it is.

Then the class the bridge holds shrinks. Its 85 virtual functions become a core of about 65 that every sampler
answers, plus ten small interfaces for what only some samplers have (a residual scale to set, a variance forest, a
basis, a count matrix, ...). The bridge asks a sampler for its capability object and gets each interface or nothing;
having the interface is the capability, so through that object the bridge's check and the engine's answer cannot
disagree. The forwarding class, which repeats every call and the default arguments of 13 virtual functions, is
deleted. The records land on their own; the fold and the shrink land together, once.

Nothing a user or a package linking to dbarts.h sees changes in SC1 to SC3, and no draw moves. One small fix comes
first: an aft sampler with a linear or gp leaf, created with an active-row mask whose rows hold one distinct
log-time, reads its leaf covariates over those rows where every other family with a response scale reads every row
and says so (run on main, below); SC0 makes aft do what gaussian does. No call is left for the maintainer; the two
that were open are made (Calls made, 13 and 14).

## Goal

After the four slices land:

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
  object matches its record and shape; a third pins, through `Sampler`, each returned refusal that a facet's
  absence now also covers.
- An aft sampler's leaf rows fall back to every row where its rows in hold fewer than two distinct log-times, as
  gaussian's do (SC0).

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
- Every name test, line by line on main 0580ff13a and on collapsed-move-rebuild as of 2026-11-10 (not yet on main; its hashes change at its landing), classified:
  scratch/plan-core/name-tests.md. The by-symbol table is [Appendix: the name tests](#appendix-the-name-tests).
- Branches that bear: collapsed-move-rebuild (worktree kmp-step0, being built) adds a logistic name test in its new
  `collapsedRescaleApplies` and a NaN-default virtual
  (`logLikelihoodAtScaledFits`, "which no caller reaches"); chain.hpp still holds 45 references. This plan is written
  against main plus that branch. host-neutral-engine (docs/plans/host-neutral-engine.md on that branch) rewrites the
  bridge's creation and refusal code in S2a and S2b and builds the facts struct in S3. apx-size edits model.hpp,
  combiner.hpp and the bridge. The multinomial `$setData` slice (dec-B640, TODO multinomial-setdata, due before
  host-neutral S2a) makes a multinomial sampler reach setData's family sites (the aft refusal, counts on a binary
  logistic, the nbinom cap, the supplied range, the logistic weights restatement), which today see "logistic" and
  after SC1 see the multinomial row: whichever of the two lands second checks each such site (SC1, step 5).
- 0.9-34, run (scratch/famtraits/probe-cpp-09.R against scratch/libs/old-ecosystem, rerun 2026-11-10): a probit
  sampler accepted `setControl` with `binary = FALSE`, still read binary afterward and kept sigma at 1, 1, 1. The
  family was fixed at creation in practice; a record of facts fixed at creation is no regression. 0.9-34 had no
  type-erased boundary and no C header of this shape, and neither aft nor linear leaves (run: its namespace holds
  neither), so SC0's corner has no 0.9-34 behavior.
- The aft leaf-row probe, run 2026-11-10 on a main build (scratch/libs/dobj-main, 0580ff13a, shipped mode; script
  and output in scratch/plan-core/aft-probe/): 200 rows, a linear leaf on x1, a creation mask keeping rows 1 to 100,
  which hold one response value (log-time 1.5, every row an event). Gaussian falls back: its leaf center and scale
  are those of all 200 rows (5.194, 5.169), `standardization.all.rows` is TRUE, and it warns "fewer than two distinct
  values of the response are in the likelihood (100 of 200 observations are in it), so the response's scale, the
  prior's center and the center and scale of "x1" are taken from all 200". aft does not: its leaf numbers are those
  of the 100 rows in (0.493, 0.280), with no warning. Both read the response range from every row, the creation mask
  being installed after the transform as the help of `new("dbartsSampler", ...)` says; with distinct log-times among
  the rows in, both read the leaf over the rows in. So aft departs from the rule
  [`LeafRows`](../../src/bartcore/model.hpp) documents ("fewer than two distinct response values are in ... every
  row is read in their place") and from the setData help's "or the response's scale is read from every row":
  [`withResponseFallback`](../../src/bartcore/sampler.hpp) tests gaussian and nbinom by name, while
  [`carriesUnits`](../../src/bartcore/chain.hpp) and `AFTResponse::lacksDistinctRowsIn` treat aft as a response
  with a scale. The only path reaching it is creation (setData refuses aft).

## Constraints

- Draws identical, bit for bit, in every slice. The latent-scale anchor, the forest sd unit and the basis-scale
  expressions move verbatim: a record stores an enumerator, never a precomputed constant (the logistic anchor stays
  `std::numbers::pi / std::sqrt(3.0)` at its reader; the bridge's logistic leaf scale stays
  `std::numbers::pi * std::sqrt(3.0)`).
- Messages identical in SC1 to SC3: every tinytest message pin unchanged, the order of the bridge's refusal cascades
  included. SC0 changes aft's creation warning in its corner to gaussian's words (SC0).
- The record is internal. It never appears in dbarts.h (its rows would become ABI); the `latent` column is the
  engine's and is not published (dec-B768). dbarts.h is untouched and `tools/check-api-hash.sh` passes with
  `DBARTS_C_API_HASH` unchanged.
- The engine stays R-agnostic: the record holds no R wording; the bridge keeps its words and keys them by kind.
- Out of scope, SC0's two R sites apart: the R layer's own "gaussian or aft" tests (R/generics.R, R/spec.R, R/bart.R), which read the record
  after 1.0 through an internal call if ever (the ruling names the engine, the bridge and the factories); the forest
  couplings' base class [`ForestCombiner`](../../src/bartcore/combiner.hpp) and its capability virtuals (chain-split's
  territory, below); the factory's half-built-sampler and shared-tail findings of the C++ review's slice 5 (not ruled);
  the header text on family growth, which is host-neutral S3's (TODO capi-sampler-info); R's per-family defaults
  (the default node scale), which stay with R and its bridge, as host-neutral's "What stays with R" has it.

## Design

### The record

```cpp
enum class ResponseKind : std::uint8_t {
  gaussian, studentT, probit, logistic, aft, ordinal, nbinom, multinomial
};
struct ResponseFacts {
  ResponseKind kind;
  bool hasResidualScale;         // a free residual scale (drawn or held is a setting)
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
};
```

The names are the implementer's to refine; the columns are the facts the three expert reviews named plus those the
sites below read. Each class declares its row as a `static constexpr ResponseFacts` member beside its code; a virtual
`facts()` returns it; `responseFacts(ResponseKind)` reads the table for code that has no object yet (the factories,
the bridge). The table is built from the list of eight classes, and its static checks are: one row per kind; row i
holds kind i; `rescalesLatents` and `scoresScaledFits` equal `std::is_base_of_v` of the step's interface (Optional
steps); and the model invariants `carriesResidualDf` implies `latent == precision`, `precisionChannelFree` implies
`latent != precision`, `hasResidualScale` implies `units == range`, `admittedByAmplitudeCoupling` implies
`latentScale != notStated`. Two values follow from others and are read, not stored: the leaf-row fallback's reading
(below), and `workingWeightsVaryPerSweep` (`latent == precision`).

The rows, each value checked against the site that reads it today:

| column | gaussian | studentT | probit | logistic | aft | ordinal | nbinom | multinomial |
|---|---|---|---|---|---|---|---|---|
| hasResidualScale | yes | yes | | | yes | | | |
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
| latentScale | fromResponse | fromResponse | standardNormal | logisticDistribution | fromResponse | notStated | notStated | notStated |
| latent | none | precision | location | precision | location | location | precision | none |
| weights | precisionOrFrequency | precisionOrFrequency | membership | counts | none | membership | counts | counts |
| support | real | real | binary | binary | real | categoryIndex | count | categoryCounts |

Notes on the rows:

- Student-t follows every "gaussian" site today, so each column keeps what that site does for it now (it draws sigma,
  carries units, admits single precision, reads its weights as gaussian's, reads the response for its leaf rows),
  except where the code already treats it apart: the variance-forest refusal's finite-df term becomes
  `precisionChannelFree` false, the df block becomes `carriesResidualDf`, and amplitude admission becomes false. That
  last one is a new engine refusal (the amplitude factory refuses gaussian with a finite residual df, which the
  amplitude chain silently ignores today); the bridge refuses the combination first with "Student-t residuals", so no
  message moves (Calls made, 6).
- `latentScale` is read only by the calibration map, which only admitted kinds reach. Ordinal, nbinom and
  multinomial state none (`notStated`): today's ordinal arm is unreachable, and admitting any of the three later fails
  the build until its row states a scale (Calls made, 7).
- Multinomial's `weights` is counts: a multinomial sampler takes count weights
  ([`enforceMultinomialWeightPolicy`](../../src/R_interface_bartcore.cpp) at creation and in setWeights), and today it
  reports logistic, so `familyCarriesNoWeights` answers false for it in `bartcore_setWeights` before the
  counts-mutation branch. A row of none would refuse it there with aft's message. The implementer confirms the
  policy matches logistic's reading of counts, else gives multinomial a value of its own.
- The leaf-row fallback reads `units`: range reads the response less its offset, countShift the response, none
  nothing. Seven rows agree with today; aft (range) does not, which SC0 fixes. If SC1 lands before SC0, aft's
  exception is written as one commented kind test in [`withResponseFallback`](../../src/bartcore/sampler.hpp) naming
  SC0, and SC0 deletes it.
- `latent` is what `latents()` returns: multinomial's Polya-Gamma variates live in its coupling, so its row says
  none. `workingWeightsVaryPerSweep` is `latent == precision` (logistic, Student-t, nbinom, the three overriding it
  today).

### The kind and the request

`ResponseFamily` stays the factories' request type (six values; Student-t is gaussian with a finite residual df).
One exhaustive function, `responseKindFor(ResponseFamily, double residualDf)`, holds the rule "a finite residual df
selects Student-t", which today lives in the chain's construction switch and the factory's variance-forest refusal;
the multinomial factory names its kind directly. A finite df is ignored off gaussian, as the chain's construction switch
ignores it today, and the bridge's "resid.df requires gaussian" refusal runs before the function is reached. The 684
tests/cpp lines that build samplers by family keep working
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
core; so does a function a caller sweeping every sampler calls unconditionally, where it is a no-op off a shape. A
refusal narrower than its facet (an index naming no forest, a non-map forest, no test rows) stays a returned refusal
inside the facet. By that rule the partition, checked function by function by the critique, which the implementer
confirms against each function's documented refusal:

| facet | functions | present when |
|---|---|---|
| `ResidualScaleControl` | setSigma | row.hasResidualScale and no variance forest |
| `CensoringControl` | setSurvivalStatus | row.carriesCensoring |
| `ResidualDfControl` | setResidualDf | row.carriesResidualDf |
| `ShapeControl` | setNegativeBinomialShape, shapeParameter | row.carriesShape |
| `BinomialCountsControl` | setBinomialCounts | built with trials (row.takesTrials) |
| `VarianceForestSurface` | predictVariance, currentVarianceFits | a variance forest is built |
| `AmplitudeCoupling` | amplitudes, setForestBasis, setForestBasisRescaled, predictPerForest, setForestWeights, setForestMapSd, setForestCoefficientPrior | the amplitude coupling |
| `CountsCoupling` | setCounts, setCategoryOffset, setCategoryTestOffset | the counts (softmax) coupling |
| `AdditiveLocation` | fitsWithoutOffset | every coupling but the counts coupling |
| `LengthscaleControl` | setLeafLengthscales | a GP leaf |

Kept in the core by the rule: `setResidualPrior` (`$setModel` calls it on every sampler; a no-op off a residual
scale), `rederiveStandardization` (a constant leaf's empty answer is the true one, and `rederiveAfterInstall` and
`$setPredictor` call it whatever the leaf), `sampleVarianceForestFromPrior` (documented as a no-op on a homoscedastic
sampler, which R's prior-draw sweep calls unconditionally) and `numForestAmplitudes` (0 off the coupling, as
`SamplerShape::numAmplitudes` answers).

That moves 20 virtuals off the core; `totalAmplitudes` goes (it duplicates `SamplerShape::numAmplitudes`) and
`capabilities()` comes, so the core holds about 65 (SC3 recounts). `SamplerShape` loses `hasVarianceForest`, `holdsTrials`,
`supportsForestWeights` and `forestReportingIsDefined` (both true exactly on the amplitude coupling, so its presence
answers them; either stays a value if the implementer finds a sampler where it is narrower), and
`supportsCountsMutation`. (`supportsActiveRows`, true for every class, goes in SC1 with its virtual and the bridge's
unreachable refusal, Calls made, 5.)

The capability object, `SamplerCapabilities`, is a trivially destructible struct of facet pointers, filled once when
the sampler is built (every presence is a creation fact; no setter moves one) and returned by `capabilities()`
(Calls made, 13). Whether `capabilities()` is a virtual computing the struct per call or a non-virtual base function
over a struct the derived constructor fills is the implementer's choice; the stored form makes "fixed for life"
structural and holds self-pointers, safe while `Sampler`'s copy is deleted (so no move is declared), and the class
says that an added move constructor would leave them dangling. The
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

The engine's own refusals stay. `Sampler<L, ResidT>` implements every facet in every instantiation, so presence is a
run-time fact, and code holding `Sampler` directly (tests/cpp, internal callers) can still call an absent facet's
function. Each returned refusal or guard that exists today stays and is pinned in tests/cpp through `Sampler`: among
them [`fitsWithoutOffset`](../../src/bartcore/chain.hpp)'s refusal off one location, `currentVarianceFits`' false
with no variance forest, and [`Chain::setSigma`](../../src/bartcore/chain.hpp)'s return under a variance forest
(Calls made, 15). `Chain::setSigma` writes sigma today on a kind with no residual scale (probit), which no bridge path
reaches; it gains the same guard read from the row, unless a tests/cpp caller relies on the write, which the
implementer then reports. Its engine callers were read: `Sampler::setSigma` (the bridge and C entries, which the facet
now gates); `Chain::setResidualPrior` (already only on a kind with a residual scale and no variance forest); and
`Chain::moveScale`, which writes back `sigma()`, that is 1 times a unit scale of 1 on a unit-scale kind, so the
guard skips a write of 1 and is neutral. The state install and the constructors write `sigma_` directly (the install
under `drawsSigma()`), not through `setSigma`, so the guard does not reach them; any other effect would show in the
trio. The claim that the bridge's
check and the engine's answer agree is a claim about calls through the capability object.

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

Three features, serial behind collapsed-move-rebuild. SC0 is its own small feature, before SC1 or beside it. SC1
is its own feature (dec-B767). SC2 and SC3 are one feature (dec-B763 is one TODO item): one branch, each slice
reviewed and gated on it, merged to main once, so main never holds `Sampler` final over the full 85-function
interface, the shape the ruling revised away from. Tier: changes draws if wrong, so the engine tier of
[Process by risk](README.md#process-by-risk) in each, an opus review that runs the mutants.

### SC0: aft's leaf-row fallback (posterior-changing in its corner; engine and R)

The probe (Context) shows aft departing from the documented rule. Steps:

1. Engine: [`withResponseFallback`](../../src/bartcore/sampler.hpp) treats aft as gaussian, reading the values aft's
   transform reads (the log-times less any offset, a censored row at its censoring time; the implementer confirms
   against `AFTResponse` at creation), so an aft creation whose rows in hold fewer than two distinct log-times reads
   its leaf covariates over every row and sets `standardization.all.rows`.
2. R, the creation warning: the `dbartsSampler` initializer in R/dbarts.R calls
   [`warnRowsFallBack`](../../R/bartcore.R) after `leafFromAllRows` with `range = model@family %in% c("gaussian",
   "nbinom")` and `byMask = model@family %in% c("logistic", "aft")`. aft joins `range` (it reads a response with a
   scale) and leaves `byMask`, which is read only where `range` is false. Without this the engine's new flag would
   raise the leaf-only lead ("... is kept by the active-row mask ...") in place of the response's.
3. R, [`scaleRowsFallBack`](../../R/spec.R): it excludes aft by name, its comment saying aft takes no weights and
   meets a mask only after creation. No path reaches it with an aft row out (dbartsSpec's creation path has no aft
   row out, setData refuses aft, xbart refuses aft); the implementer confirms by tracing and leaves it, or adds aft
   if a path is found.
4. The words aft gets, as gaussian's: (a) two or more rows in holding one distinct log-time (the probe's case, today
   silent): "fewer than two distinct values of the response are in the likelihood (100 of 200 observations are in
   it), so the response's scale, the prior's center and the center and scale of "x1" are taken from all 200". (b)
   fewer than two rows in (one row in): today the leaf-only "one observation is in the likelihood (1 of 200 is kept
   by the active-row mask), so the center and scale of "x1" are taken from all 200"; after SC0 gaussian's and
   Student-t's one-row text, "fewer than two distinct values of the response are in the likelihood (1 of 200
   observations is in it), so the response's scale, the prior's center and the center and scale of "x1" are taken
   from all 200", as rule 9 of the rows-in design asks ("one text with the response's where the response falls back
   too, as one row in always makes it under the families that read it"). test-leaf-standardization-rows.R's
   "one row in, aft" pin changes to that text.
5. Pins: the probe's setup beside "one row in, aft", expecting every row's numbers, the flag TRUE and exactly one
   warning, counted with `withCallingHandlers` (Gate hygiene: never `suppressWarnings`); the changed one-row pin.
6. Design note: one sentence in rule 9 of docs/design/response-scale-rows.md, saying aft follows rule 4 there as
   gaussian does, its log-times the response read, censored rows at their censoring times. No NEWS entry: 0.9-34
   had neither aft nor linear leaves.

Gates, the posterior-changing class's (README, RNG classes and their gates): the neutral gates; the trio compared
bitwise against the current baselines, which stand in for re-recording and the z mode because the corner is reached
by none of their scenarios (a full count of identical-draws lines shows nothing outside it moved); the four
snapshots; and exact-gates.yaml's aft entry, `Rscript benchmarks/R/aft-exact.R quick` on the reference build, which
must pass. aft-exact has no linear leaf and no mask, so it too shows only that aft outside the corner kept its
exact posterior; inside the corner the posterior moves by design, to the standardization the rule names, and the
pins state that standardization exactly.

An agent call at the engine tier: it restores a documented rule, and no reasonable maintainer keeps the departure
(the critique's reading, A3). If SC1 lands first, SC0 also deletes SC1's commented aft kind test.

Size: ~15 engine lines, ~10 R lines, ~70 tinytest, a help sentence where the help lists the families, one design
sentence. Cost: agent ~3, review ~1.5, gates ~2.

### SC1: the records (engine, factories, bridge)

Steps:

1. Recount the name tests on the slice's base by the rule in scratch/plan-core/name-tests.md; any test a landing
   added since (k-move-pricing, approximate-draws, multinomial-setdata) joins the slice. The recount also names the
   tests that stand for a kind without its name: the counts-mutation fingerprint for multinomial
   (`supportsCountsMutation`, about six bridge sites, which SC3 moves) and `model.responseIsBinary`.
2. model.hpp: `ResponseKind`, the column enumerations, `ResponseFacts`, a row in each of the eight classes, the
   table and its static checks, `responseFacts`, `responseKindFor`; the optional-step interfaces; `facts()`.
   The fact virtuals go (`supportsActiveRows`, `workingWeightsVaryPerSweep`, the three `carries*` queries become
   non-virtual reads of `facts()`), and the four unit-scale classes drop their four scale overrides for base
   defaults of 1, 1, 0, 1 (Calls made, 8). Remove the stale comment "only the binary families fix it".
3. chain.hpp, facade.hpp, sampler.hpp: the 43 F sites read the row, `withResponseFallback` reading `units` (with
   aft's commented kind test if SC0 has not landed); the two construction switches key on the kind (Student-t its
   own arm in the single-forest one); `family_` goes from `Chain`; the multinomial chain and sampler carry kind
   multinomial.
4. The bridge: the 35 F sites read `responseFacts(kind)`; `refusedAmplitudeFamilyReason`,
   `enforceBinaryWeightPolicy`, `defaultLeafScale` and `supportFamily` key on the kind with their wording or constant
   kept at the reader (`defaultLeafScale` mirrors R's default, a model setting the record does not hold), and `supportFamily`
   maps Student-t to its own kind (its "documentary" gaussian arm goes); [`refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp)
   keeps its cascade order, the family reason read from the request's base kind so "Student-t residuals" still
   fires where it does; `dbarts_sampler_family` maps the kind. `SamplerShape` gains `responseKind`, loses `family`,
   `carriesShape`, `carriesResidualDf` and `supportsActiveRows` (with the bridge's unreachable refusal). Wherever a
   switch gains a Student-t arm, that arm does what the gaussian arm does today.
5. Before each bridge site is converted, its implementer names the tinytest pin that would fail if the site read
   the wrong column, in a table in the landing note; where none exists (the unresolved-sigma refusal of a flat
   creation and the single-precision refusal looked thin), a pin joins the slice. Every site a multinomial sampler
   reaches today reading `logistic` gets the same answer from the multinomial row, or is shown unreachable; the
   landing note lists each, `bartcore_setWeights`' `familyCarriesNoWeights` among them, and a bridge-level or
   tests/cpp pin holds multinomial's count weights accepted there.
6. tests/cpp, a new suite `test_response_facts.cpp` (added to the Makefile's sources and main.cpp): for each of
   the eight classes, built on a small fixture, `facts()` is its table row and its virtual answers match the row: `initialSigma`, `fitScale`, `fitShift`,
   `sigmaScale` and `getScale` against `units` and `hasResidualScale`; `latents()` non-null exactly where `latent`
   is not none, after a refresh; working weights changing across `refreshLatents` exactly where `latent` is
   precision; the global blocks present exactly where their column says; `setBinomialCounts` accepted exactly on a
   trials-built logistic; `setSurvivalStatus` reaching the class exactly on aft. Plus a bitwise pin of the logistic
   anchor at `std::numbers::pi / std::sqrt(3.0)` through the test peer.

Size: ~300 engine lines written and ~120 removed, ~180 bridge lines touched, ~270 tests/cpp. Cost: agent ~8, review
~4, gates ~3.

### SC2: the fold

Steps: `Sampler` final over `SamplerBase` as in [The fold](#the-fold); `SamplerFacade` deleted; the factories build
`Sampler` directly; default arguments off every virtual, the short forms as non-virtual base overloads;
test_facade.cpp's forwarding rows go (each compared a call through the base with the same call on `Sampler`, which
`override` now makes one function); the rows that check engine behavior no other suite checks (the partly filled
and wrapped saved-tree store, the two-chain and two-forest confusables) move to test_sampler.cpp, listed in the
landing note; test_shape.cpp reads `Sampler::shape()`. A clang-tidy run with `google-default-arguments` over
facade.hpp and sampler.hpp, clean: clang-tidy is not on PATH here but at /opt/homebrew/opt/llvm/bin/clang-tidy, run
with tests/cpp's include flags after `--` (or a compile database). It does not see a wrong default in a non-virtual
overload, which mutant 13 covers. Reviewed and gated on the SC2-SC3 branch; not merged alone.

Size: ~250 lines written; ~400 engine and ~1,100 tests/cpp removed. Cost: agent ~5, review ~2.5, gates ~2.

### SC3: the core and the capability object

Steps: the facet interfaces and `SamplerCapabilities` as in [The core and the capability object](#the-core-and-the-capability-object);
`Sampler` implements the facets and fills the object at construction from the row, the coupling's existing
predicates, the leaf type and the variance forest; the 20 functions leave `SamplerBase`; `totalAmplitudes` and the
presence fields go; the bridge's 29 presence reads and about 30 facet calls ask `capabilities()`, the C entries'
capability statuses among them; tests/cpp `test_capabilities.cpp`, built on test_facade.cpp's fixtures (gaussian,
probit, logistic, binomial, Student-t, aft, ordinal, nbinom, BCF, multinomial, variance forest, linear, GP,
monotone): each fixture's object has exactly the facets its record and shape imply, each facet's functions answer
through the facet as through `Sampler`, and each returned refusal or guard an absent facet also covers is pinned
through `Sampler` ([The core and the capability object](#the-core-and-the-capability-object)). The SC2-SC3 branch
merges once, after SC3's review and the branch's full gates.

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
pin unchanged (SC0: its new pin passes, every other pin unchanged); each trio compare prints its full count of "identical draws (same RNG stream)" lines and no "max |z|"
line (Gate hygiene); the four snapshots pass. SC1 and SC3, which edit the bridge's validators and capability checks,
add the R-loaded ASan run over the creation, mutation and test-capi.R suites; SC2, which changes ownership in the
factories, adds the x86 leg with LeakSanitizer. The SC2-SC3 branch runs the full battery again at its tip before its
one merge. The orchestrator's own check, each slice: `bench-sampler.R record`
on the base and `compare` on the tip, same quiet machine, interleaved; expected neutral (the record is read a handful
of times per sweep per chain, never per row).

## Mutants

The reviewer runs each on the slice's tree, rebuilds `--preclean`, and records which gate failed; a mutant no gate
catches is a finding, answered with a pin. `touch` each reverted file.

SC1:
1. Student-t's `precisionChannelFree` true: test_response_facts fails, and a variance forest with a Student-t
   residual is no longer refused (test-heteroscedastic.R's "variance forest does not support" pin).
2. Student-t's `hasResidualScale` false: sigma pinned on a Student-t fit; equivalence.R's Student-t scenario
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
9. `requireResolvedSigmaEstimate` read for aft as having no residual scale: an aft flat creation with an NA sigma
   is accepted; the pin step 5 names or adds fails.
10. Student-t's `admitsSinglePrecision` false: a single-precision Student-t creation is refused. The row check does
    not see it (the mutant edits the row itself); the pin step 5 names or adds, a single-precision Student-t
    creation, fails.
11. Multinomial's `weights` none: `bartcore_setWeights` refuses a multinomial sampler's count weights with aft's
    message; the bridge-level or tests/cpp pin step 5 adds fails (R's `$setWeights` refuses multinomial first, so no
    existing tinytest sees it).

SC0:
12. aft left out of the fallback again: SC0's tinytest pin fails (every-row numbers, the flag and the warning).
12a. aft left in `warnRowsFallBack`'s `byMask` and out of its `range`: the pins' warning text fails.

SC2:
13. A default argument re-added to a virtual with a value different from the base overload's: clang-tidy flags it.
14. A non-virtual short form passing `SplitPlacement::byValue` where the default is `byPosition`: clang-tidy does not
    see it; the tests/cpp row that calls the short form and compares with the long form's `byPosition` fails (the
    reviewer confirms one exists, else adds it).
15. `Sampler::shape()` fills `numSavedDraws` from the store's capacity: test_shape fails.
16. `createSampler` builds the GP instantiation for a linear designation: the leaf tests fail.

SC3:
17. `varianceForest` present on a homoscedastic gaussian: test_capabilities fails; predictVariance reaches a sampler
    with no variance forest.
18. `residualScale` present under a variance forest: the setSigma capability arms of test-capi.R and the tinytest
    pinned-sigma refusal on a heteroscedastic sampler fail (the reviewer confirms both arms exist, else adds them).
19. `countsCoupling` absent on multinomial: `dbarts_sampler_family` and the R count setters' pins fail.
20. `binomialCounts` keyed on the row alone (present on a binary logistic): test-count-response-doors.R's refusal
    fails.
21. `fitsWithoutOffset`'s refusal off one location removed, now that the facet is absent there: the tests/cpp pin
    through `Sampler` fails.

## Open calls

None. Nothing here changes what a user sees beyond SC0's corner, which restores a documented rule; the two calls
that were open are made (Calls made, 13 and 14), for close-out.

## Calls made

1. The records first (engine and bridge together, one slice), as the expert synthesis advised ("Q3 first, so the
   capability object is filled from the record"). Alternative: the bridge's half its own slice (one more gate battery
   and review, about 4 hours).
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
7. `latentScale` is `notStated` for ordinal, nbinom and multinomial, all unreachable, so admitting any of them forces
   the question (the critique's advice; the first draft gave ordinal `standardNormal`). Alternatives: ordinal's
   model value, which changes an unreachable arm away from verbatim; today's empirical arm for every non-binary row.
8. The four scale accessors default to the unit scale in the base, the tests/cpp row check catching a class that
   should override. Alternative: keep them pure, the four unit-scale classes repeating them.
9. No default arguments on virtuals; short forms as non-virtual base overloads. Alternative: the non-virtual
   interface idiom, a public non-virtual function over a private virtual for every entry.
10. The bridge's "logistic restates its weights after setData" stays a kind test: it is family behavior, and the
    nearest fact (counts and precision latents) also covers nbinom, whose setData does not need it.
11. The forwarding test's rows deleted, not kept as smoke tests, the behavior rows moved.
12. The couplings' base class and its capability virtuals left to chain-split.
13. The capability object is one struct of typed facet pointers, null where absent, read by one query (option (a)),
    because it is the memo's "queried capability object" with the fewest virtuals, the answer and the refusal are one
    fact, and dbarts.h's "test once" holds by construction. Alternatives: (b) a query virtual per facet, the same
    design with about 20 virtuals in place of one; (c) a second all-in interface, each function still refusing,
    which does not meet the ruling.
14. SC2 (the fold) before SC3 (the shrink) on the one branch, because the end state is the same and the forwarding
    class and its test table are not rewritten into facets only to be deleted (about 2 agent-hours). With one merge
    main never sees the order. Alternative: SC3 first, the facets implemented by `SamplerFacade`.
15. Each returned refusal or guard an absent facet also covers stays in `Sampler` and is pinned in tests/cpp through
    it. Alternative: the facet overrides private in `Sampler`, reachable only through the facet pointer, and the
    whole-facet refusals turned into asserts, which churns test_sampler.cpp's direct calls.
16. The default node scale stays in the bridge, a switch over the kind with its constants at the reader: it is R's
    default, a model setting, not a fact of the kind. Alternative: a record column, which puts a front end's default
    in the engine and has no C++ reader for multinomial's value.
17. No record column for the leaf-row fallback: it derives from `units`, and aft's one departure is fixed in SC0
    rather than recorded as a fact of the kind. Alternative: a column that keeps the departure, which the record
    would present as part of the model.
18. SC0 is an agent call at the engine tier, planned here: it restores the rule `LeafRows` and the help state, in a
    corner the trio does not reach. Alternative: hold it for the maintainer, which keeps a documented rule broken
    meanwhile.

Implementation choices the critique leaves to the implementer, each reported in the landing note: the capability
object's const view (two structs or one templated on constness, the reviewer checking for a const hole); stored or
computed `capabilities()`; record names; `admittedByAmplitudeCoupling` as a column or as the coupling's predicate
over the row (the factory's stated reasons are `notStated`, censoring and a global block); the scale-accessor defaults
(call 8); refining `latent` to the augmentation (none, truncated normal, Polya-Gamma, gamma mixture), which
approximate-draws' method refusal will need.

## Not verified

- The facet partition was checked by reading (the critique traced each presence against its function's refusal),
  not by a build; the implementer confirms it.
- Whether multinomial's count-weight policy matches logistic's reading exactly (A1): if not, multinomial's
  `weights` gets a value of its own.
- What value aft's fallback reads for a censored row at creation (SC0 names it after reading `AFTResponse`).
- Settled since the first draft: the aft fallback (run, Context); bcf-equivalence reaches the logistic anchor (its
  logistic scenarios build from a fresh data object, so `latentScaleAnchor` sets every forest's node scale; the
  critique's reading).

## Cost

| slice | lines | agent-hours |
|---|---|---|
| SC0 aft fallback | ~95 written | ~6.5 |
| SC1 records | ~750 written, ~120 removed | ~15 |
| SC2 fold | ~250 written, ~1,500 removed | ~9.5 |
| SC3 core and capabilities | ~580 touched, ~100 net removed | ~11.5 |
| total | ~1,680 written, ~1,700 removed | ~42.5 |

Plans have run 1.5 to 2 times low: 64 to 85 hours.

## Appendix: the name tests

Counted by `grep -c 'ResponseFamily::\|RF::'` on main 0580ff13a, chain.hpp on collapsed-move-rebuild as of 2026-11-10 (not yet on main; its hashes change at its landing):
chain.hpp 45, facade.hpp 7, sampler.hpp 3, combiner.hpp 1, R_interface_bartcore.cpp 67, C_interface.cpp 6; 129.

| where | references | becomes |
|---|---|---|
| [`Chain`](../../src/bartcore/chain.hpp) single-forest and amplitude constructors, the construction switches | 12 | switches over the kind |
| `sigmaIsFixed_` at both constructors, [`setResidualPrior`](../../src/bartcore/chain.hpp) | 6 | hasResidualScale |
| [`buildVarianceForest`](../../src/bartcore/chain.hpp)'s gate, [`varianceForestIsRefused`](../../src/bartcore/facade.hpp) | 4 | precisionChannelFree |
| [`carriesUnits`](../../src/bartcore/chain.hpp) | 3 | units |
| [`latentScaleAnchor`](../../src/bartcore/chain.hpp), [`derivedForestSdUnit`](../../src/bartcore/chain.hpp), [`rowsForBasisScale`](../../src/bartcore/chain.hpp) | 18 | latentScale |
| [`forestRescaleApplies`](../../src/bartcore/chain.hpp), `collapsedRescaleApplies` (rebuild) | 2 | the optional-step interfaces |
| the multinomial chain and sampler's logistic, `Chain`'s `family_` member | 3 | the multinomial row, the row pointer |
| [`createConstantLeafSampler`](../../src/bartcore/facade.hpp) | 1 | admitsSinglePrecision |
| [`trialsAreRefused`](../../src/bartcore/facade.hpp) | 1 | takesTrials |
| [`createAmplitudeSampler`](../../src/bartcore/facade.hpp) | 3 | admittedByAmplitudeCoupling |
| [`withResponseFallback`](../../src/bartcore/sampler.hpp) | 2 | units (aft's departure fixed in SC0) |
| [`AmplitudeSpec`](../../src/bartcore/combiner.hpp)'s family default | 1 | kept (a request) |
| [`resolveFamily`](../../src/R_interface_bartcore.cpp), the resid.df request check | 7 | kept (parsing) |
| [`applySurvivalStatus`](../../src/R_interface_bartcore.cpp), the setResponse status and setData refusals | 3 | carriesCensoring |
| [`refusedAmplitudeFamilyReason`](../../src/R_interface_bartcore.cpp) | 6 | admittedByAmplitudeCoupling, wording by kind |
| [`defaultLeafScale`](../../src/R_interface_bartcore.cpp) | 6 | a switch over the kind, constants at the reader |
| the creation weights reading, the aft weights refusal, [`familyCarriesNoWeights`](../../src/R_interface_bartcore.cpp) | 4 | weights |
| [`enforceBinaryWeightPolicy`](../../src/R_interface_bartcore.cpp) | 5 | a switch over the kind, wording kept |
| the nbinom shape and weighted-count caps at creation, setResponse, setData, setWeights | 5 | carriesShape |
| [`requireResolvedSigmaEstimate`](../../src/R_interface_bartcore.cpp), [`sigmaIsPinned`](../../src/R_interface_bartcore.cpp), [`fillResidualParameters`](../../src/R_interface_bartcore.cpp) | 6 | hasResidualScale |
| [`familyHasResponseScale`](../../src/R_interface_bartcore.cpp), the count family's supplied and stored ranges | 8 | units |
| [`validateResponseSupport`](../../src/R_interface_bartcore.cpp), setData's counts on a binary logistic | 8 | support, takesTrials |
| the single-precision refusal | 1 | admitsSinglePrecision |
| setData's logistic weights restated | 1 | kept (family behavior) |
| [`supportFamily`](../../src/R_interface_bartcore.cpp) | 7 | a map to the kind |
| [`dbarts_sampler_family`](../../src/C_interface.cpp) | 6 | a map from the kind |

Replaced by record reads 78; kept in switches over the kind 42; parsing and request defaults 8; family behavior 1.
