# host-neutral-engine: a plain-C specification the engine is built from, build-time math and printing, a status from every entry point, and a creation entry in dbarts.h

Status: PLANNED 2026-11-03 on main (c3d268c8f), for dec-B85; amended the same day after a blind critique (BUILD AFTER
AMENDMENTS, then BUILD), for the maintainer's rulings on every open call, and rebased onto main at 8fce28568 after forest-door
landed. Not built. Rulings, quoted where they settle a call: dec-B85 (the staged route), dec-B84 (the header is pure C),
dec-B119 (a status from every entry point, the unwind to the entry point staying), dec-A74 (printing chosen at build time, the
run-time print pointers removed), dec-A73 (a plain-C sigma getter rides this header change), dec-B124 (the specification's
predictor part is `dbarts_predictor_source`, its dense channel stated per column), dec-A50 (the engine's R touchpoints), dec-A20
(creation's refusal channel), with dec-B86 (the header ships what a consumer calls), dec-B87 and dec-A33 (copy-on-set), dec-B112
(destroy stays; "A C consumer should have a path to invalidate a C object"), dec-B56 and dec-B57 (the list, the hash, per-symbol
lookup), dec-B60, dec-A139 and dec-B424 (the version pair stays 1/0 before the release, the hash moves), dec-B111 (the
hash-to-minor CI assertion), dec-B160 and dec-B163 (multinomial run and predict through the flat API, its creation a post-release
door), dec-B186 (the C consumer test compiles and runs on CRAN), dec-A36 (consumers check the boundary under sanitizers
themselves), dec-A32 (a returned integer is a value or a capability status), dec-C01 (full symbol names), dec-C04 (the support
library's print pointers), dec-B571, dec-B572, dec-B575, dec-B582, dec-B577, dec-B578, dec-B579 and dec-B581 (this plan's
open calls, all ruled, and the status words; see [Ruled](#ruled)), and the frequency-weights rulings dec-B563 to dec-B570 and
dec-B580 (the reading of a weight is fixed at creation and the header gains its field in this change; real weights where the
draw is exact; probit and ordinal membership as 0/1 weights; a zero weight alike on every family, its predictive draws NA; a
binomial response for logistic). The runs behind it: reads of the tree at
c3d268c8f, of stan4bart main (27ef163) and treatSens master (268b2db), and of frequency-weights (da210ce2d); the critic's probe of
R 4.6.2's nmath against Rmath (17,304 points, below). No prototype was built here.

agent: S1 is its own feature, one opus implementer (effort high). The main feature: S2a, S2b and S3 one opus implementer each
(effort high), S2b the riskiest; S4 one sonnet implementer per consumer, in that consumer's own checkout; S5 sonnet. One opus
review per engine slice (S1, S2a, S2b, S3) that runs the mutants below, the same reviewer checking the fix round unless the fixes
change engine behavior. Tier "Changes draws" for every engine slice ([Process by risk](README.md#process-by-risk)): nothing is
meant to move a draw, but every creation route and every flat entry is rewritten, and the gates are the draw-changing tier's.
rng: neutral, bit for bit, on the R build: the equivalence trio, the four snapshot files, tests/cpp and the full tinytest suite
do not move, and no R message changes. The standalone build is new; its math is bitwise equal to Rmath on the reference host and
within a relative floor elsewhere; see [RNG class](#rng-class).
window: before the 1.0-0 submission (dec-B85: "post merging into main, but before 1.0 since it might change a number of packages
and is an interface change"). Two features. S1 (build-time math and printing, branch host-neutral-math) touches no header and no
creation code and lands on its own whenever it is green, serial with other engine merges. The main feature (branch
host-neutral-engine: S2a, S2b, S3, then S5) starts after forest-door and coefficient-prior-live land and merges once (dec-B494),
the consumer ports (S4) merging right after it.
budget: ~5,650 lines written in this tree plus ~6,000 copied unedited from R's math library and ~400 of reduced R headers: S1
~350, S2a ~900, S2b ~1,600, S3 ~2,800; S4 ~80 in stan4bart and ~50 in treatSens. Agent time, estimated and not measured: S1 ~8,
S2a ~20, S2b ~29, S3 ~26, S4 ~6, S5 ~1, rebases onto main after the pending plans ~3: ~93. Plans here run 1.5 to 2 times low,
so 140 to 185. The callback statuses and the status words (dec-B579, dec-B581) are about 300 of S3's lines and 3 of its hours. Counted under dec-B575's middle course (51 fields), with the forward-compatibility work it asks for
([Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial)):
about 300 lines and 4 hours over the bare middle course, which itself took 400 lines and 10 hours off every one-forest model.

## Summary for the maintainer

You ruled in September that before 1.0-0 the engine stops depending on R below the layer that talks to R, so that a program
other than R, a Python package for instance, can build and run a sampler through the one C header dbarts ships. This plan builds
that as two pieces of work.

The first is the math and the printing, and it can land on its own as soon as it is ready. The engine and its support code call
six of R's numerical routines (the normal distribution's density, probability and quantile; the t distribution's density and
quantile; the chi-squared quantile) and R's console. Under R they keep calling exactly those, so every draw an R user makes stays
bit for bit the same. Any other build uses a copy of R's own routines, kept in the repository but left out of the package R users
install, and writes to standard output. Your ruling named four routines; the engine has since come to use two more, the t
quantile and a chi-squared median, both read only when a sampler is built or its prior is changed, and this plan copies them too.
A check of the copy found it gives exactly R's answers when compiled the way R compiles, and differs in the last few digits
otherwise. The two unused printing hooks go.

The second is the rest. Today the code that turns R's control, model and data objects into an engine reads those objects
directly. It will instead fill one plain C description of the sampler - its data, family, priors, constraints, how its weights
are read, and run settings, each already resolved - and the engine is built from that description alone. R still does everything
R does now: reading the formula, choosing the family from the response, estimating the starting sigma, turning the residual
prior's quantile into a scale, parsing the prior vocabulary, and naming columns and forests in its messages. Nothing an R user
sees changes. Then the header changes once, so the two compiled packages rebuild once: a program can create a sampler from that
description, every entry that can fail returns a status in place of raising an R error, the failure's message can be read back,
each chain's current sigma can be read, and the predictor description says column by column which block a column's values sit
in, as you ruled. stan4bart and treatSens are then changed to check those statuses; with the helper you chose, an R user of
either sees the same errors, interrupts and time limits as today.

You have ruled every choice this plan raised. Entries that cannot fail keep returning their value. A failure's message is read
from one call, per thread, and a second call gives fixed words for each status, the same text for every program, which the R
helper falls back on and never replaces with words of its own. A C program creates the middle of the model in 1.0, with the
description shaped so that models of several forests and the multinomial model can be added later without changing what any
1.0 field means; the plan says what that addition looks like and tests that a 1.0 description reads the same under it. The
description is one struct with a size field and an initializer. Creating from C requires a seed. And a callback's return now
says what happened: zero continues, a positive value stops the run on purpose, a negative value reports the callback's own
failure, and a C++ exception thrown from a callback is caught on whatever thread ran it and reported the same way.

What it costs: about 5,650 lines here plus the copied math code, and roughly 140 to 185 agent-hours, most of it in rewriting how
samplers are created. The forest door has landed; the main piece still waits for the live coefficient prior, which reworks the
same creation code. Frequency weights' field joins the new description whichever lands first.

## Goal

The engine, the support libraries and the flat C entry file compile and run with no R header and no R library: a tests/cpp
recipe builds them into a plain C consumer that creates, runs, predicts, reads sigma and refuses through statuses alone, and no R
symbol is in its binary. Under R every draw, every message and every printed line is what main produces at the base of each
branch. Every creation route, R's and the new flat one, fills one plain-C specification and passes one translate-and-validate
stage. Every flat entry that can fail returns a `dbarts_status`, its message readable on the calling thread. stan4bart and
treatSens build against the new header and pass their suites.

## Context

- The rulings above, read whole in docs/decisions.md; TODO's "host-neutral engine" item, whose list of touchpoints this plan
  completes (below) and whose "Rewords the Python-binding entry" S5 carries out.
- The engine's creation inputs are already plain C++ structs with borrowed pointers: [`SamplerOptions`](../../src/bartcore/chain.hpp),
  [`AmplitudeSpec`](../../src/bartcore/combiner.hpp), [`ForestSpec`](../../src/bartcore/combiner.hpp),
  [`MultinomialSpec`](../../src/bartcore/combiner.hpp), [`PredictorSource`](../../src/bartcore/data.hpp), consumed by the
  factories [`createSampler`](../../src/bartcore/facade.hpp), [`createSamplerOverStore`](../../src/bartcore/facade.hpp),
  [`createAmplitudeSampler`](../../src/bartcore/facade.hpp) and [`createMultinomialSampler`](../../src/bartcore/facade.hpp). What
  is R-shaped is the stage above them.
- The R parse: [`parseSamplerSpecification`](../../src/R_interface_bartcore.cpp) fills three bridge structs from the S4
  objects, [`optionsFromParsed`](../../src/R_interface_bartcore.cpp) turns them into options (reading the DART prior off the model
  object directly), and [`createHolder`](../../src/R_interface_bartcore.cpp) adds the survival status, the variance record, the
  forest records and bases, the move record, the chain generators ([`createChainRngs`](../../src/R_interface_bartcore.cpp)), the
  factory call through [`createSamplerOrRaise`](../../src/R_interface_bartcore.cpp), the test predictors and the verbose summary
  ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)). Two more routes build engines,
  [`createMultinomialCountsHolder`](../../src/R_interface_bartcore.cpp) and
  [`createMultinomialDataHolder`](../../src/R_interface_bartcore.cpp); [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp)
  builds a `ColumnStore`, not an engine, and [`bartcore_createFromHandle`](../../src/R_interface_bartcore.cpp) builds xbart's
  fold engines over it.
- The flat entries: [`translateSource`](../../src/C_interface.cpp) validates a caller's predictor source (on R's transient
  stack); every refusal is an `Rf_error`; engine exceptions become R errors through
  [`callConvertingExceptions`](../../src/R_interface_bartcore_common.hpp); a run's callback raise travels as an
  [`UnwindJump`](../../src/R_interface_bartcore_common.hpp) the entry resumes; the poll
  [`userInterrupted`](../../src/R_interface_bartcore.cpp) consumes either an interrupt or an R error raised while polling (a
  `setTimeLimit` expiry, pinned in test-monotone.R), and [`raiseInterrupt`](../../src/R_interface_bartcore.cpp) re-raises the
  one it consumed; a slow monotone count warns through [`raiseSlowCountWarning`](../../src/C_interface.cpp). The handle is the
  bridge's holder, [`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp), which carries no R type.
- R frames and C++ exceptions: [`unwindProtect`](../../src/R_interface_bartcore.cpp) protects a body against R's jump and catches
  no C++ exception, so an exception must never leave a body it wraps or cross any R frame. The shared validators in
  R_interface_bartcore_common.hpp (`validateResponseSupport`, `enforceBinaryWeightPolicy`, `refuseNonFinite` and thirteen more)
  raise `Rf_error` today and are called at about 60 sites in the R bridge (counted by name at c3d268c8f, definitions included),
  besides the flat file's; the bridge's other 379 `Rf_error` sites are R's own and stay.
- Process-wide state on the predict path: [`predictPartition`](../../src/bartcore/sampler.hpp), a plain struct holding a
  `std::vector`, is written by every predict and read by a test hook; it is safe today only because predict is main-thread-only.
  [`testFitPartition`](../../src/bartcore/chain.hpp) is relaxed atomics, last writer wins, a diagnostic no answer reads.
- Engine exception types: `std::invalid_argument` (ten sites), `std::length_error` (two, saved-tree storage past the addressable
  size), `std::bad_alloc` (one), and `std::runtime_error` at three: the two prior-draw exhaustions in chain.hpp
  (sampleTreesFromPrior) and predict's join, which rethrows any worker failure, an allocation failure included, as
  `runtime_error("predict worker w failed: ...")`.
- The R touchpoints below the bridge, as found at c3d268c8f (TODO lists the first group only):
  - in the engine: `Rf_pnorm5` (eleven calls), `Rf_dnorm4` (two) and `Rf_dt` (one) called by name in
    src/bartcore/model.hpp; `ext_quantileOfStudentT` in [`forestMapEntry`](../../src/bartcore/combiner.hpp) (a coefficient's
    t prior, at creation and at a live restatement); `ext_quantileOfChiSquared` in
    [`ChiKHyperprior::medianOfChi`](../../src/bartcore/model.hpp) (the chi prior's median, read by the calibration map's leaf
    prior size); `ext_printf` throughout sampler.hpp and tree.hpp.
  - in the support library: `Rf_pnorm5` and `Rf_qnorm5` declared by hand in src/external/random.c for the truncated normal
    draws ([`ext_rng_simulateTruncatedNormalScale1`](../../src/external/random.c)), and `ext_quantileOfNormal` in
    [`ext_rng_simulateStandardNormal`](../../src/external/randomNorm.c), which every normal draw goes through (inversion is the
    default algorithm), so qnorm is on the hot path; R's own generators, `.Random.seed` and `set.seed` in
    src/external/randomBase.c; `Rf_error`, `Rf_warning` and `Rprintf` in src/external/io.c; the unused
    [`misc_printf`](../../src/include/misc/io.h) and `misc_flushOutput` pointers in misc and their five installers.
  - in the bridge's own core: `R_alloc`, `R_finite`, `R_UnwindProtect` (the callback leaf, [`ShippedDrawHook`](../../src/R_interface_bartcore_common.hpp)),
    `R_CheckUserInterrupt` (the poll).
  - The residual prior's chi-squared quantile is already the bridge's (`sigmaRawScale`), as dec-B85 has it.
- Build: src/Makevars.in sets no `OBJECTS`, so every .c and .cpp at the top of src/ is compiled into dbarts.so. configure (which
  reads `R CMD config`) generates src/config.hpp, src/misc/config.h, src/external/config.h and src/include/misc/types.h from
  their `.in` files; Windows uses checked-in `.win` variants. R's compiler on the reference host is `clang -arch arm64` with
  `-falign-functions=64 -Wall -g -O2`, so contraction of a * b + c is clang's default, on.
- tests/cpp links libR and has no standalone target. The rshim the TODO's predecessors name was deleted in c0ddd93e8; nothing in
  the tree builds the engine without R.
- The shipped header: [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h), [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)
  checked against [`dbarts_apiToken`](../../src/C_interface.cpp), [`dbarts_predictor_source`](../../inst/include/dbarts/dbarts.h)
  (today a dense column's index is read in the block its STORE kind selects, the rule dec-B124 replaces),
  [`DBARTS_SWAP_STATUS_LIST`](../../inst/include/dbarts/dbarts.h). The shipped test consumer
  [`capi_hash`](../../inst/tinytest/capi/consumer.c) and its handshake arms.
- Consumers, read-only: stan4bart main calls run (2 sites), sampleTreesFromPrior, setOffset (2), setSigma (2), getLatents (2),
  predict, setTreeStorage (2), printTrees, setVerbose, numTrees, destroy, the version pair and six accessors, all in
  src/init.cpp, plus `dbarts_results` in src/bart_util.hpp and `dbarts_dense_predictor_source`. treatSens master calls run (3),
  setResponse (4), setOffset, setSigma, setVerbose (2), setNumThreads (2), destroy (2) and the version pair, in
  src/R_interface.cpp, src/bartTreatmentModel.cpp and src/sensitivityAnalysis.cpp. Both define `DBARTS_USE_STUBS`; neither
  defines `DBARTS_REQUIRE_EXACT_ABI`; neither creates through C or registers a draw callback. bartCause and bairrtt link nothing
  of dbarts.
- docs/plans/pure-c-header.md (the last header change, its slices and gates), docs/design/consumer-spec-surface.md (why one
  resolution home), docs/design/r-c-division.md, error-style.md (R7, a flat entry names itself; R8, the R message is canonical
  and a C backstop is worded apart), docs/design/frequency-weights.md on its branch, [Callbacks and errors](../architecture.md#callbacks-and-errors),
  [Reproducibility contract](../architecture.md#reproducibility-contract), docs/plans/python-bindings.md.

## Constraints

- Neutral on the R build: no draw, message, warning or printed line changes; the gates in [Tests and gates](#tests-and-gates).
- No C++ exception crosses an R frame or leaves an `unwindProtect` body; no R jump crosses a C++ frame that owns anything.
- The header stays pure C (dec-B84): `gcc -std=c99 -pedantic` and `g++ -std=c++11` on its prototype view with no R include
  path (the tests/cpp header-compiles target). The version pair stays 1/0 and the hash is re-baked (dec-B60, dec-B424).
- No run-time hook for math, printing or interrupts (dec-B85, dec-A74): each is chosen when the code is compiled.
- Nothing only a standalone build compiles sits where the R package build compiles it or ships in the source tarball.
- The R surface is untouched: no R function, argument, return shape or help text changes, except the help's C section
  (dbarts-embedding.Rd).
- Out of scope: a C creation path for multi-forest or multinomial models (dec-B575 and dec-B160: a later append, whose
  shape this plan fixes); an opaque C state buffer (dec-B84, conditional on a consumer); a shipped standalone library or build
  system beyond the tests/cpp recipe; any R wrapper for a C-created handle; the SamplerShape fields that size the bridge's arrays
  (dec-A50 names them; no ruling asks this arc to move them); frequency weights' and the binomial response's own work (below);
  compatibility for development consumer binaries built against today's header (below).

## What stays with R

The R parse keeps every step that needs R or that R alone resolves, and hands the specification values, never R objects:

- ingestion: the formula and matrix interfaces, `dbartsData`, factors and their levels, sparse and mixed containers, `Surv`,
  the hazard expansion, `cbind(successes, failures)` once dec-B567 lands, `na.action`;
- resolution ([`resolveSamplerSpec`](../../R/spec.R)): the family from the response shape (the spec takes a stated family,
  never `DBARTS_FAMILY_AUTO`), ordinal recoding, the weight policy and `weights.type = NULL`'s reading (dec-B564: the spec takes
  the resolved reading), the starting sigma estimate, offset normalization, the prior vocabulary and its evaluation, the
  per-family defaults;
- the residual prior's chi-squared quantile (dec-B85): the spec carries the resulting scale;
- seeds drawn from R's stream when the control states none: R hands the chain seeds (see [Seeding](#seeding));
- every refusal that names an R thing - a column or forest by its R name, a slot, `$setData`, `dbartsSpec()` - stays raised in
  the R parse before the spec is filled, worded as now; the shared stage repeats the predicate in its own words (R8);
- the verbose creation summary, which prints the residual prior's quantile the spec does not carry;
- the records a re-created sampler reads (cut grids and their mass, the columns seen missing, a leaf standardization, the
  calibration map's anchor and unit, a numeric basis column's scale), the multi-forest records, multinomial counts and the data
  handle's prebuilt store: these ride the internal extension (below), which only the R bridge fills.

## Design

### The specification

Two plain-C structs, one shipped and one internal:

- `dbarts_sampler_spec`, in dbarts.h: a caller-filled, `structSize`-led struct (the `dbarts_predictor_source` contract in the
  caller-to-library direction) holding the middle of the model dec-B575 names, every value resolved: 51 fields after
  `structSize` (the ruling's "about 45" was an estimate). Its fields, grouped (the order below is the struct's):
  - data: `x` (a `dbarts_predictor_source`, its typing channel READ at creation, which gives `columnTypes` and `categoryCounts`
    the reader dec-B124 owes them), `y`, `weights`, `weightsType` (a `dbarts_weights_type`, below, read only when `weights` is non-NULL), `offset`, `survivalStatus`
    (aft, double 0/1 per row as R holds it), `xTest` (a second source, 0 rows for none), `testOffset`, `numCuts` and
    `maxNumCuts` (one or per column), `useQuantiles`;
  - family: `family` (a `dbarts_family`, 1 to 7), `sigmaEstimate`, `sigmaIsFixed`, `sigmaPriorDf`, `sigmaPriorScale` (R's
    `qchisq(1 - q, df) / df`), `residualDf` (read under STUDENT: positive fixes it, 0 estimates it), `numCategories` (ordinal
    K), `shape` (nbinom: positive fixes r, non-positive estimates it);
  - leaf prior: `k`, `nodeScale`, `namedSd` (NaN none), `kIsDrawn`, `kPriorDf`, `kPriorScale`;
  - tree prior: `numTrees`, `base`, `power`, `splitProbabilities`, `useDart` and DART's six (`dartA`, `dartB`, `dartRho`,
    `dartAlpha`, `dartUpdateAlpha`, `dartUpdateDelay`);
  - run: `numChains`, `numThreads`, `numThin`, `keepTrees`, `numSamplesToStore`, `verbose`, `printEvery`,
    `singlePrecisionResidual`, `hasSeed`, `seed`, `chainSeeds`;
  - engine limits: `categoricalExhaustiveCap`, `testFitParallelCutoff`, `predictParallelCutoff`, `sparseDensityThreshold`.

  Appended in a later minor when a binding asks (dec-B575), each group taking effect only when its fields are non-zero: the leaf
  models (`leafModel`, giving `dbarts_leaf_model` its first reader, leaf covariates, the GP settings), monotone fits, the
  constraints (columns, interactions, blocks), the move settings and the variance forest; then several forests and multinomial
  ([Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial)). Until then the R parse hands these to the build stage
  through the internal extension, in the layout the append will ship.

  `DBARTS_SAMPLER_SPEC_INIT` and `dbarts_sampler_spec_init()` (dec-B582) set `structSize` and `dbarts()`'s defaults for a
  gaussian response wherever they are numbers, NULL for every vector, `weightsType` PRECISION, and no seed (`hasSeed` 0). The
  header lists the values R resolves per family (probit's node scale, for one), since the spec carries no "the family's default"
  marker: one resolution home, as dbartsSpec()'s design argued.
- The weights' reading: `DBARTS_WEIGHTS_TYPE_LIST`, PRECISION 1 and FREQUENCY 2, a field of the family part as dec-B564 has it
  ("the shipped header gains one field in dec-B85's change"). The spec holds the resolved reading, never "the family's own": R
  resolves `weights.type = NULL` before filling it. With `weights` NULL the field is not read, so the initializer's PRECISION
  refuses no spec without weights (a logistic one, for instance). With weights, the build checks the pair against dec-B563,
  dec-B564 and dec-B568 to dec-B570 (FREQUENCY on the gaussian-type families at real weights and the Polya-Gamma families at
  whole numbers; PRECISION only where there is a residual scale; 0/1 weights under either, membership on probit and ordinal;
  a 0 a row out of the likelihood on every family, logistic included) with the words frequency weights gives R, prefixed.
  Until frequency weights' engine work lands, FREQUENCY is taken only where it reads as the family's own (the Polya-Gamma
  families' counts, 0/1 weights). The spec says nothing of predictive draws: a 0 row is still drawn, so the stream is what it
  was, and its draws are set to NA where dec-B580 puts it, on every family, its fitted value still reported.
- `bartcore_bridge::SpecExtension`, internal (src/spec_internal.h, plain C structs and arrays), in two parts. The PUBLIC-SHAPED
  part holds what later appends will ship, already in their sketched layout: the one-forest groups 1.0 leaves out (leaf models,
  monotone, constraints, moves, variance forest), the several-forest records (`ForestRecord`, field for field the sketched
  `dbarts_forest_spec`) and the multinomial record (`MultinomialRecord`, the sketched multinomial fields). The R-ONLY part holds
  what no C caller will state: the re-creation records (cut grids and mass, columns seen missing, a leaf standardization, the
  calibration map's anchor and unit, a numeric basis column's scale as first computed, a sharer's share as first set), the data
  handle's prebuilt store, the creation mask and the unit count (`tooFewUnitsIn`). Empty, it is the flat entry's.

One stage, `buildSampler(const dbarts_sampler_spec&, const SpecExtension&)`, in a new R-free src/sampler_build.cpp: validates
the spec with neutral words (the 1.0 reading rules of
[Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial)
included), translates both sources through the one `translateSource` (dec-B124: the flat predict entry,
creation from C, and the R parse's views all pass it), fills the engine structs, makes the chain generators from the seeds,
selects and calls the factory, installs the move record, the test predictors and the test offset, and returns a holder or a
refusal. The engine creation routes (createHolder, the two multinomial holders, createFromHandle over the handle's store) and
`dbarts_sampler_create` call it; createDataHandle fills the spec's predictor part, passes `translateSource` and builds its store
from the result. What createHolder does today in that order moves there whole; only the SEXP reads stay behind.

### Forward compatibility with several forests and multinomial

dec-B575 asks the 1.0 shape to be ready for creating several-forest and multinomial models from C, so that adding them later is
an append under dec-B86's minor bump and no 1.0 field changes meaning. The later append, sketched (not shipped at 1.0):

```c
/* appended to dbarts_sampler_spec, after every 1.0 field and every one-forest group appended before it */
size_t numForests;                    /* 0: the one forest the 1.0 tree and leaf fields describe */
size_t forestStride;                  /* sizeof(dbarts_forest_spec) as the caller compiled it */
const dbarts_forest_spec* forests;    /* numForests elements, forestStride bytes apart */
const int* counts;                    /* MULTINOMIAL: n x K counts, category-major; y unread */
const int* trials;                    /* MULTINOMIAL: per-row trial count n_i */
const double* categoryOffset;         /* MULTINOMIAL: n x K offset, category-major; offset unread */

typedef struct {                      /* one forest of a several-forest sampler */
  size_t numTrees; double base, power;                    /* tree prior */
  double sd, sdMultiple; int updateAmplitude; double amplitudeDf;   /* size: the forest statement */
  const size_t* columns; size_t numColumns;               /* columns it may split on, NULL all */
  int kind;                           /* dbarts_forest_kind: NONE 0, LEVELS 1, NUMERIC 2 */
  const double* basis; size_t numBasisColumns;            /* the forest's coefficient basis */
  /* then, as the one-forest groups append them: interactions, blocks, moves, per forest */
} dbarts_forest_spec;
```

The element is the engine's [`ForestSpec`](../../src/bartcore/combiner.hpp) (a `ForestStatement`, a `ForestStructureSpec`, a
kind and a basis) in C. It carries no share: forest-door, landed, keeps `ForestSpec::share` as a sharer's part of its group's
prior variance "as a first creation set it", which "No caller states"; the build derives the group (`forestSharers`: the forests
with no basis that state no sd and no df), and the share rides the extension's R-only part as a re-creation record. The
multinomial part is the engine's
[`MultinomialSpec`](../../src/bartcore/combiner.hpp): K is the 1.0 field `numCategories`, the shared category forest is the 1.0
tree and leaf fields (`numTrees`, `base`, `power`, `k`, `nodeScale`), and counts, trials and the category offset are the three
appended pointers. The result side is already there: the flat run sizes `dbarts_results` by F
(`dbarts_sampler_numFittedValuesPerObservation`, K or the forest count) and V (`dbarts_sampler_numVariableCountForests`), and
multinomial runs and predicts through the flat API today (dec-B163), so the append adds nothing to `dbarts_results`.

Invariants every 1.0 description keeps, under the append and any later one:

1. Zero means absent. Every appended field's all-zero bytes (0, NULL) mean "not used": the build behaves as the build before the
   append did. `DBARTS_SAMPLER_SPEC_INIT` and `dbarts_sampler_spec_init` zero everything they do not set.
2. A 1.0-sized description to a later build: the build reads `structSize` bytes and treats every field past them as zero, so a
   1.0 description is one forest (`numForests` 0) of a 1.0 family, read exactly as at 1.0. This is the existing `structSize`
   contract of `dbarts_results` and `dbarts_predictor_source`, in the caller-to-library direction.
3. A later-sized description to an earlier build: the build accepts it only if every byte past the `structSize` it knows is
   zero, and otherwise returns UNSUPPORTED naming the first unknown nonzero byte's offset, the rule Linux's extensible structs
   (`openat2`, `clone3`) follow. A newer caller using only 1.0 features runs on a 1.0 library; one asking for several forests
   is refused, never silently given one.
4. The 1.0 tree prior fields (`numTrees`, `base`, `power`) describe THE forest when there is one: the single forest, or under
   MULTINOMIAL the forest every category forest copies, as `MultinomialSpec` shares one `MultinomialForestSpec`. With
   `numForests` > 0 they are not read; every forest states its trees in its element, and nothing is inherited from them. The
   1.0 leaf prior fields (`k`, `nodeScale`, `namedSd`, `kIsDrawn`, `kPriorDf`, `kPriorScale`) keep one meaning under every
   structure: the model's leaf prior, the single forest's, or under several forests the total the sharing forests divide, "the
   single-forest chain's leaf prior with one k for the group", as forest-door built it; a forest stating its own sd or df is
   sized by its element instead. A sampler-wide field (data, family part, sigma, DART, run settings, engine limits) keeps its
   meaning under every structure.
5. `family` keeps naming the likelihood. Several forests is a structure over GAUSSIAN, PROBIT or LOGISTIC, not a family value;
   MULTINOMIAL (8) is a family. Values are never renumbered: AUTO 0 is refused at C creation; 1 to 7 are created at 1.0;
   MULTINOMIAL and every value past 8, reserved for families appended later, return UNSUPPORTED at 1.0 ("this build creates no
   multinomial sampler from C" or "family %d is not one this build creates").
6. A field is read only under the families or structures that name it, and that set never shrinks: `numCategories` under ORDINAL
   at 1.0 and MULTINOMIAL later; `residualDf` under STUDENT; `shape` under NBINOM; `y` and `offset` under every 1.0 family, not
   under MULTINOMIAL, which states the appended counts and category offset instead.
7. Nested arrays version by a stride, not a size per element: the parent states `forestStride`, the caller's
   `sizeof(dbarts_forest_spec)`, once; the build indexes by it and reads each element under invariants 1 to 3 at that size. The
   elements of one array share a type, so a size in each would repeat the stride; Win32's `SendInput(nInputs, pInputs, cbSize)`
   versions an array of structs the same way. A stride of 0 with `numForests` > 0 is INVALID_ARGUMENT.

What 1.0's validation refuses now, so the append is never a change of meaning: a nonzero byte past the known `structSize`
(invariant 3); MULTINOMIAL and any family value past 8 (5); and nothing else, since every other 1.0 field already has its final
reading. Placeholder fields are not shipped (see Calls made).

The later work is mostly exposure. From S2b the R path already builds several-forest and multinomial samplers through
`buildSampler` from the extension's `ForestRecord` array and `MultinomialRecord`, laid out as the sketch, and the equivalence
trio's BCF and multinomial baselines hold those routes bitwise. The append moves the two records' public prefixes into dbarts.h
as `dbarts_forest_spec` and the multinomial fields, adds `dbarts_forest_spec_init(element, size)` for bindings, and has
`dbarts_sampler_create` copy the caller's array into the extension at its stride before the same `buildSampler` call; the
R-only part (a sharer's share and a basis scale as first set, the anchor and unit) stays R's.

The gate (S3): a tests/cpp target builds the library with `-DDBARTS_TEST_SPEC_APPEND`, under which dbarts.h appends the sketch
above to the spec, and the build stage honors it. A C99 translation unit compiled WITHOUT the define fills a 1.0-sized
description (structSize the 1.0 size) of a gaussian and a probit sampler; a second, compiled WITH it, fills the same values with
the append zeroed; both create through the appended build, and their draws equal the shipped build's at the same spec and seed,
bit for bit. The same target creates a two-forest sampler through the test append and compares it bitwise with the R-created
BCF fit of the bcf-equivalence scenario it copies, and hands the shipped build a description carrying a nonzero appended
field, which returns UNSUPPORTED. The two struct definitions never meet in one translation unit, and both consumers are C.

### Seeding

The spec takes either `hasSeed` 1 with `seed` (every 32-bit value valid, 0 included) and `chainSeeds` NULL - one generator
seeded from it hands each chain its seed, which is the `control@seed` path today and gives the same draws - or `chainSeeds`, one
per chain, used as they stand, which is how the R parse passes the seeds it draws from R's stream when the control states none.
R's draws are therefore unchanged on both paths. A spec with neither is refused (dec-B577), and a binding fills the seed from
its own environment's stream, as dbarts' R code fills it from R's.

### Refusals inside the shared code

A refusal is a value, not an exception, wherever it can meet an R frame: `bartcore_bridge::Refusal`, a status and a fixed
1,024-byte message, trivially destructible (a static_assert, so an R jump past one strands nothing), empty on success. The
shared validators return one, marked `[[nodiscard]]`; an R site raises it with `raiseRefusal(r)`, which calls
`Rf_error("%s", ...)` with exactly today's text, and a flat site returns it. Below the validators, `buildSampler` and the engine
throw, and one catch per creation route (the shape [`createSamplerOrRaise`](../../src/R_interface_bartcore.cpp) has) turns any
exception into a refusal value before any R call, so no exception crosses an R frame. The about 60 validator sites in the R
bridge change from a call to a wrapped call; the 379 R-only `Rf_error` sites do not change. The review checklist: every
validator call is consumed, every `buildSampler` call sits under its route's catch, and no `unwindProtect` body lets an
exception out.

### Statuses

`dbarts_status`, an ABI enum generated like the others (`DBARTS_STATUS_LIST`), replacing `DBARTS_SWAP_STATUS_LIST`. The header's
rule: zero is success, any positive value is a success with a note, any negative value a failure; codes may be added in a later
minor, so a caller tests the sign and treats an unknown negative as a failure.

| value | name | meaning |
|---|---|---|
| 0 | `DBARTS_OK` | done |
| 1 | `DBARTS_SCALE_KEPT` | setResponse or setOffset swapped and kept the transform in force (today's SCALE_KEPT) |
| 2 | `DBARTS_STOPPED` | a draw callback's positive return stopped a run on purpose (dec-B581) |
| -1 | `DBARTS_UNSUPPORTED` | the sampler cannot do this at all (today's capability 0) |
| -2 | `DBARTS_INVALID_ARGUMENT` | this call is wrong (today's `Rf_error` "this call is wrong") |
| -3 | `DBARTS_OUT_OF_MEMORY` | an allocation failed or a size exceeds what can be addressed |
| -4 | `DBARTS_INTERRUPTED` | the build's poll stopped a run: an interrupt, or under R an error raised while polling (a time limit) |
| -5 | `DBARTS_CALLBACK_FAILED` | a draw callback returned a negative value or threw a C++ exception (dec-B581) |
| -6 | `DBARTS_FAILED` | the engine could not complete the call; the message says why |

Engine exceptions map by type: `std::invalid_argument`, `std::domain_error` and `std::out_of_range` to INVALID_ARGUMENT,
`std::bad_alloc` and `std::length_error` to OUT_OF_MEMORY, the rest (the two prior-draw exhaustions among them) to FAILED.
Predict's join (S3) rethrows a failed worker's exception in its own category, an allocation failure as an allocation failure,
with the message text it has today, so R's words hold and a C caller sees OUT_OF_MEMORY.

Per entry (the infallible ones keep returning their value, dec-B571), with what a failure leaves:

| entry | returns | can fail with | after a failure |
|---|---|---|---|
| `dbarts_sampler_create` | status, handle through an out pointer | INVALID_ARGUMENT, UNSUPPORTED (multinomial), OUT_OF_MEMORY, FAILED | `*sampler` NULL, nothing allocated |
| `dbarts_sampler_run` | status | INVALID_ARGUMENT (structSize 0), INTERRUPTED, CALLBACK_FAILED, OUT_OF_MEMORY, FAILED; STOPPED | the chains stand where they stopped, a valid state to run on from; this call's results and saved trees discarded; after OUT_OF_MEMORY or FAILED, valid to read and destroy |
| `dbarts_sampler_sampleTreesFromPrior` | status | FAILED (every draw of a tree left an empty leaf), OUT_OF_MEMORY | trees already drawn stand; drawing again is valid |
| `dbarts_sampler_setResponse`, `setOffset` | status | UNSUPPORTED, INVALID_ARGUMENT (support, finiteness, updateScale refusals); SCALE_KEPT | unchanged |
| `dbarts_sampler_setSigma` | status | UNSUPPORTED (pinned), INVALID_ARGUMENT (not finite and positive) | unchanged |
| `dbarts_sampler_getSigmas` (new) | status | UNSUPPORTED (heteroscedastic, as R's getSigmas() returns NULL) | out untouched |
| `dbarts_sampler_getLatents` | status | UNSUPPORTED (augments nothing) | out untouched |
| `dbarts_sampler_predict` | status | UNSUPPORTED (blend undefined), INVALID_ARGUMENT (source, offset, empty store), OUT_OF_MEMORY | sampler unchanged, out unspecified |
| `dbarts_sampler_setTreeStorage` | status | OUT_OF_MEMORY | the store empty, as any change of storage leaves it |
| `dbarts_sampler_printTrees` | status | INVALID_ARGUMENT (an index out of range) | nothing printed |
| `dbarts_sampler_setVerbose` | status | INVALID_ARGUMENT (printEvery 0) | unchanged |
| `dbarts_sampler_numTrees` | status, count through an out pointer | INVALID_ARGUMENT (forest out of range) | out untouched |
| `destroy`, `setDrawCallback`, `setNumThreads`, the other accessors, the version triple, `dbarts_lastErrorMessage`, `dbarts_statusString`, `dbarts_sampler_spec_init` | value, void, or OK | nothing (spec_init: INVALID_ARGUMENT on a structSize below the first release's) | - |

The message: `const char* dbarts_lastErrorMessage(void)`, the message of the most recent call on this thread that returned a
negative status, a capability probe included, self-prefixed with the entry's name as today (`"dbarts_sampler_run:
results.structSize is 0 - ..."`), valid until the next failing call on the thread; a success leaves it (dec-B572). The record
is a fixed thread-local buffer, so writing it never allocates; every failing call overwrites it, with the empty string when no
message could be formed (memory ran out before one was), so an empty message always belongs to the latest failure. The record
also notes, for INTERRUPTED under R, whether the poll consumed an interrupt or an error, which the R raise reads.

The words: `const char* dbarts_statusString(int status)` (dec-B579), fixed words per status, as zlib's `zError`, SQLite's
`sqlite3_errstr` and libcurl's `curl_easy_strerror` give. It cannot fail, needs no handle and no thread record, and returns
static strings:

| status | words |
|---|---|
| `DBARTS_OK` | `success` |
| `DBARTS_SCALE_KEPT` | `success; the response transform was kept` |
| `DBARTS_STOPPED` | `stopped by the draw callback` |
| `DBARTS_UNSUPPORTED` | `not supported by this sampler` |
| `DBARTS_INVALID_ARGUMENT` | `invalid argument` |
| `DBARTS_OUT_OF_MEMORY` | `out of memory` |
| `DBARTS_INTERRUPTED` | `interrupted` |
| `DBARTS_CALLBACK_FAILED` | `the draw callback failed` |
| `DBARTS_FAILED` | `the sampler could not complete the call` |
| any other | `unknown status` |

A binding names its exceptions from it, and it is the text a caller falls back on when `dbarts_lastErrorMessage()` is empty.
The table is generated from `DBARTS_STATUS_LIST`, whose X-macro gains the words beside each value, so a status cannot be added
without its words, and the words are folded into the hash with the values.

Two things stay as they are, because a status cannot carry them: the stub view's handshake, consumer code in an R-only view,
still raises; and the slow-count warning goes out through the build's warning routing, last, as now (under R an R warning, which
`options(warn = 2)` can still turn into a jump; standalone, a line on standard output). A draw callback's R error on an inline run
under R stays R's jump too ([The draw callback](#the-draw-callback)).

### The draw callback

dec-B581 makes the callback's return sign decide, and catches what a callback throws on every thread and in every build. The
contract, which replaces the RETURN and RAISING paragraphs of `dbarts_draw_callback` today:

- RETURN. 0 continues. A POSITIVE return stops the run on purpose: every chain stops at its next sweep boundary and
  `dbarts_sampler_run` returns STOPPED, a positive status, with this call's results and saved trees discarded (the sample
  cursors have not advanced past draws already written into their slots, exactly as on the interrupt path). A NEGATIVE return
  reports the callback's own failure: the run stops the same way and returns CALLBACK_FAILED, the message
  `dbarts_sampler_run: the draw callback for chain %zu returned %d`. When chains return differently in one sweep, a failure
  wins over a stop, and the first failing chain by index is the one named, so the status does not depend on thread timing.
- EXCEPTIONS. A C++ exception escaping the callback is caught AT THE CALL, on the main thread and on every worker, in the R
  build and the standalone one, and the run stops and returns CALLBACK_FAILED with `dbarts_sampler_run: the draw callback for
  chain %zu threw: <what()>` (`an exception not derived from std::exception` where there is no `what()`). An R
  package sees it through `DBARTS_R_CHECK` as an R error with that message; today an inline exception reached R as an error
  carrying `what()` alone, and no consumer registers a callback. Worker threads called the callback unprotected until now, so a
  worker's exception ended the process; the catch is new there.
- R ERRORS. Under R, on an inline run (min(numThreads, numChains) <= 1), the callback MAY raise an R error: the call is made
  under `R_UnwindProtect`, the jump becomes a C++ unwind at the callback, every frame to the entry is destroyed, and the entry
  resumes R's jump, exactly as now; the run returns nothing, and the sampler is left as a stop leaves it. From a worker an R
  raise stays undefined, as now: the header says to return negative instead.
- A BINDING catches its own exception in its wrapper, keeps it in the context, returns -1, and re-raises it once the run returns
  CALLBACK_FAILED, as C callbacks from Python are written; libgit2's callbacks abort with a negative value the call returns.

The catch is in the build's callback unit (src/host_r.cpp beside the `R_UnwindProtect` leaf, standalone/host_standalone.cpp's
plain call), shared by the inline call and the worker trampoline, so one wrapper serves every thread. About 40 engine and bridge
lines, as the ruling estimated; the three-way stop record (continue, stop, failure with chain and value or text) replaces
today's one stop flag.

### Threads

A handle is used by one thread at a time; distinct handles may be used at once from distinct threads; under the R build, `run`,
`sampleTreesFromPrior` and `printTrees` stay main-R-thread-only (the poll, R's warning, and printing through `Rprintf`), and every other entry may be called from any
thread once `predict` no longer allocates through `R_alloc`. That is the contract C libraries with handles state (libcurl: a
handle in one thread at a time, different handles in different threads; SQLite's multi-thread mode: any thread, no connection in
two at once), and what a Python binding that releases the GIL around predict on a thread pool would assume. `predictPartition`
becomes `thread_local`: predict writes it and the test hook reads it on one thread, so the hook still works under R and two
handles predicting on two threads no longer race (see Calls made for keeping the main-thread rule instead).
`testFitPartition` stays as it is: relaxed atomics, no race, read by nothing that answers.

### Build-time math, printing and interrupts (S1, and S3 for the poll)

One preprocessor switch, `DBARTS_HOST_STANDALONE`, defined by the standalone recipe and never by the R package's Makevars, so
the R build compiles as it does today.

- Math: src/include/external/stats.h gains the six names with every flag spelled: `ext_pnorm(x, mu, sd, lowerTail, logP)`,
  `ext_qnorm`, `ext_dnorm(x, mu, sd, logP)`, `ext_dt(x, df, logP)`, `ext_qt(p, df, lowerTail, logP)`, `ext_qchisq`. Under R each
  is a macro for the Rmath name the code calls today (`Rf_pnorm5` and so on), so the object code is the same; standalone each is a
  macro for a renamed copy, `dbarts_nmath_pnorm5` and so on. The existing convenience names (`ext_quantileOfNormal`, ...) are
  defined over them. A slim src/include/external/nmath.h declares the six without pulling Rmath.h, for random.c, which today
  declares two by hand. model.hpp's fourteen direct calls become the ext_ names.
- The copy: standalone/nmath/, a new top-level directory listed in .Rbuildignore, so neither the package build nor the tarball
  sees it. It holds the closure of the six in R 4.6.x's nmath (4.6.2 probed: 25 files, 5,963 lines, bd0.c to toms708.c, with
  nmath.h and dpq.h), unedited except a force-included prefix header that renames every external symbol so a program that also
  links libRmath sees no clash; reduced copies of Rmath.h, Rconfig.h and the eight R_ext headers it includes (about 400 lines);
  and a README naming the R version, the file list, R Core's copyright and the GPL (>= 2), the package's own license. It is
  compiled with `MATHLIB_STANDALONE` and `-include errno.h`. Since it is not in the tarball, DESCRIPTION gains no author or
  Copyright field. The pin moves when the reference build's R moves minor, and the probe below is rerun then.
- Printing: src/include/external/io.h maps `ext_printf` and `ext_fflush_stdout` to `Rprintf` and `R_FlushConsole` under R as
  now; standalone, to `printf` and `fflush(stdout)`. standalone/io_standalone.c gives `ext_printMessage` and
  `ext_issueWarning` ("Warning: " first) on standard output, as dec-A74 has it ("standard output for a standalone build"), and
  `ext_throwError`, which prints there, flushes and calls `abort()` (its two callers are unreachable invariants, an unknown
  generator kind). The misc pointers and their five installers go (dec-A74, dec-C04).
- R's generators: the parts of randomBase.c that use `unif_rand`, `norm_rand`, `.Random.seed` and `set.seed` compile under R
  only; the Mersenne Twister and the rest compile everywhere.
- Interrupts and the callback leaf (S3): a unit per build, src/host_r.cpp (the R_ToplevelExec poll and its record, the
  R_UnwindProtect leaf, the slow-count warning, `dbarts_R_raiseStatus`, all moved from the bridge; compiled into dbarts.so) and
  standalone/host_standalone.cpp (no poll, a plain call to the callback with C++ exceptions caught). A binding supplies its own
  unit in its own build, as dec-A74 has it for printing.
- Configuration: the standalone recipe cannot run configure, which reads `R CMD config`. standalone/config/ holds checked-in
  config.hpp, misc/config.h, external/config.h and misc/types.h for a C++20 compiler on a 64-bit POSIX system, written from the
  `.in` files as the `.win` variants are and kept in step by tools/check-win-drift.R, extended to that set. Their SIMD settings
  are configure's choice on the reference host, and the recipe builds misc through misc's own Makefile with the flag variables
  configure would pass (its per-file ISA flags: partition_avx2.c at the AVX2 flag, partition_sse4_1.c at SSE4.1, and the
  rest), so the standalone kernels are the R build's and the bitwise standalone-versus-R gate holds on the kernels too.
- The recipe: a tests/cpp target (tests/cpp is already build-ignored) that compiles the engine, src/C_interface.cpp,
  src/sampler_build.cpp, the misc and external sources, and standalone/ with `-DDBARTS_HOST_STANDALONE`, standalone/config ahead
  of src/ on the include path and no R include path, into objects under tests/cpp/build-standalone/, never into src/, where a
  later `R CMD INSTALL` would link a stray object into dbarts.so; it links no libR. It is the recipe a binding copies; nothing
  else is shipped.

### The flat entry file without R

src/C_interface.cpp and the shared helpers it calls move off R: `R_alloc` scratch becomes owned vectors, `R_finite` becomes
`std::isfinite`, the shared validators return refusals ([Refusals inside the shared code](#refusals-inside-the-shared-code)), and
R_interface_bartcore_common.hpp splits into an R-free core header (holder, `Refusal`, validators, `buildSampler`) and the
SEXP-facing declarations that stay R-only.

### Handles

A C-created handle is owned by its caller: `dbarts_sampler_destroy` frees it whole, engine and holder, and the handle is then
gone, so a second destroy is the double free C's `free` is. An R-owned handle behaves as today: destroy releases the engine, the
R object keeps the holder, and a second destroy does nothing. The holder records which it is.

### The sigma getter

`dbarts_sampler_getSigmas(const dbarts_sampler*, double* out)` writes each chain's current residual sd on the original response
scale, numChains values, as `$getSigmas()` returns them ([`bartcore_getSigmas`](../../src/R_interface_bartcore.cpp)): the pinned
value on the fixed-unit-scale families, and UNSUPPORTED, writing nothing, on a heteroscedastic sampler, where R returns NULL.

### The header diff

The list after the change (prose and Doxygen are listed after it; `int` returns marked V are values, every other non-void return
is a `dbarts_status`):

```c
#define DBARTS_C_API_LIST(X) \
  X(int, dbarts_apiMajorVersion, (void), ())                                   /* V */ \
  X(int, dbarts_apiMinorVersion, (void), ())                                   /* V */ \
  X(uint64_t, dbarts_apiHash, (void), ())                                      /* V */ \
  X(const char*, dbarts_lastErrorMessage, (void), ())                          /* new */ \
  X(const char*, dbarts_statusString, (int status), (status))                  /* new, dec-B579 */ \
  X(dbarts_status, dbarts_sampler_spec_init, \
    (dbarts_sampler_spec* spec, size_t structSize), (spec, structSize))        /* new */ \
  X(dbarts_status, dbarts_sampler_create, \
    (const dbarts_sampler_spec* spec, dbarts_sampler** sampler), (spec, sampler)) /* new */ \
  X(void, dbarts_sampler_destroy, (dbarts_sampler* sampler), (sampler)) \
  X(dbarts_status, dbarts_sampler_run, \
    (dbarts_sampler* sampler, size_t numBurnIn, size_t numSamples, dbarts_results* results), \
    (sampler, numBurnIn, numSamples, results))                                 /* was void */ \
  X(dbarts_status, dbarts_sampler_sampleTreesFromPrior, (dbarts_sampler* sampler), (sampler)) /* was void */ \
  X(void, dbarts_sampler_setDrawCallback, \
    (dbarts_sampler* sampler, dbarts_draw_callback fn, void* context), (sampler, fn, context)) \
  X(dbarts_status, dbarts_sampler_setResponse, \
    (dbarts_sampler* sampler, const double* y, int updateScale), (sampler, y, updateScale)) /* was int */ \
  X(dbarts_status, dbarts_sampler_setOffset, \
    (dbarts_sampler* sampler, const double* offset, int updateScale), (sampler, offset, updateScale)) \
  X(dbarts_status, dbarts_sampler_setSigma, (dbarts_sampler* sampler, double sigma), (sampler, sigma)) \
  X(dbarts_status, dbarts_sampler_getSigmas, \
    (const dbarts_sampler* sampler, double* out), (sampler, out))              /* new */ \
  X(dbarts_status, dbarts_sampler_getLatents, (const dbarts_sampler* sampler, double* out), (sampler, out)) \
  X(dbarts_status, dbarts_sampler_predict, \
    (dbarts_sampler* sampler, const dbarts_predictor_source* xTest, const double* offsetTest, \
     size_t numThreads, double* out), (sampler, xTest, offsetTest, numThreads, out)) \
  X(dbarts_status, dbarts_sampler_setTreeStorage, \
    (dbarts_sampler* sampler, int keepTrees, size_t numSamplesToStore), \
    (sampler, keepTrees, numSamplesToStore))                                   /* was void */ \
  X(dbarts_status, dbarts_sampler_printTrees, (... unchanged parameters ...), (...)) /* was void */ \
  X(void, dbarts_sampler_setNumThreads, (dbarts_sampler* sampler, size_t numThreads), (sampler, numThreads)) \
  X(dbarts_status, dbarts_sampler_setVerbose, \
    (dbarts_sampler* sampler, int verbose, size_t printEvery), (sampler, verbose, printEvery)) /* was void */ \
  X(size_t, dbarts_sampler_numObservations, (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numPredictors, (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numTestObservations, (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numChains, (const dbarts_sampler* sampler), (sampler)) \
  X(dbarts_status, dbarts_sampler_numTrees, \
    (const dbarts_sampler* sampler, size_t forest, size_t* numTrees), (sampler, forest, numTrees)) /* was size_t */ \
  X(size_t, dbarts_sampler_numSavedSamples, (const dbarts_sampler* sampler), (sampler)) \
  X(int, dbarts_sampler_kIsSampled, (const dbarts_sampler* sampler), (sampler)) /* V */ \
  X(int, dbarts_sampler_usesDart, (const dbarts_sampler* sampler), (sampler))   /* V */ \
  X(int, dbarts_sampler_family, (const dbarts_sampler* sampler), (sampler))     /* V */ \
  X(size_t, dbarts_sampler_numFittedValuesPerObservation, (const dbarts_sampler* sampler), (sampler)) \
  X(size_t, dbarts_sampler_numVariableCountForests, (const dbarts_sampler* sampler), (sampler))

/* the status list carries each status's fixed words, read by dbarts_statusString */
#define DBARTS_STATUS_LIST(X) \
  X(DBARTS_OK, 0, "success") \
  X(DBARTS_SCALE_KEPT, 1, "success; the response transform was kept") \
  X(DBARTS_STOPPED, 2, "stopped by the draw callback") \
  X(DBARTS_UNSUPPORTED, -1, "not supported by this sampler") \
  X(DBARTS_INVALID_ARGUMENT, -2, "invalid argument") \
  X(DBARTS_OUT_OF_MEMORY, -3, "out of memory") \
  X(DBARTS_INTERRUPTED, -4, "interrupted") \
  X(DBARTS_CALLBACK_FAILED, -5, "the draw callback failed") \
  X(DBARTS_FAILED, -6, "the sampler could not complete the call")

/* R-hosted view only (DBARTS_USE_STUBS), dec-B578 */
#define DBARTS_C_API_R_LIST(X) \
  X(void, dbarts_R_raiseStatus, (dbarts_status status), (status))
#define DBARTS_R_CHECK(call) \
  do { dbarts_status dbarts_s_ = (call); if (dbarts_s_ < 0) dbarts_R_raiseStatus(dbarts_s_); } while (0)

/* dbarts_draw_callback keeps its type; its return now has three meanings (dec-B581) */
typedef int (*dbarts_draw_callback)(void* context, const dbarts_draw* draw);  /* 0 go, > 0 stop, < 0 failed */
```

Types: `DBARTS_STATUS_LIST` and `dbarts_status` added (the list takes three arguments, so it gets its own enumerator macro beside
`DBARTS_ENUMERATOR`, and the library's word table expands the third), `DBARTS_SWAP_STATUS_LIST` and `dbarts_swap_status` removed;
`DBARTS_WEIGHTS_TYPE_LIST` and `dbarts_weights_type` added; `dbarts_sampler_spec` and `DBARTS_SAMPLER_SPEC_INIT` added; in
`dbarts_predictor_source` the field `denseChannels` (per column, k >= 0 double column k of `denseValues`, k < 0 code column ~k
of `denseCodes`, NULL every dense column in `denseValues` at `columnSources[j]`) is added above the boundary beside `denseCodes`
and `numDenseCodeColumns`, which move above it too, and the kind-selected reading of `columnSources` for a dense column goes
(dec-B124): the struct then is the engine's [`PredictorSource`](../../src/bartcore/data.hpp) field for field, `columnTypes` and
`categoryCounts` read at creation and checked and ignored by predict. `DBARTS_PREDICTOR_SOURCE_INIT` and
`dbarts_dense_predictor_source` follow. `dbarts_results`, `dbarts_draw`, `dbarts_column_type`, `dbarts_leaf_model` and
`dbarts_family` keep their layout and values (the family and leaf enums gain a reader in the spec).

`dbarts_R_raiseStatus` raises, never returns, and writes no words of its own (dec-B578): it chooses R's condition and forwards
text it reads through the neutral calls. For INTERRUPTED it re-raises the condition the poll consumed: R's interrupt condition as
the flat run raises it today, or an R error whose text is `dbarts_lastErrorMessage()`, where the R build's poll recorded the R
condition's own message (a time limit's "reached elapsed time limit"). For any other negative status it raises
`Rf_error("%s", m)`, m being `dbarts_lastErrorMessage()`, or `dbarts_statusString(status)` when that is empty. Its list is folded
into the hash and stubbed only under `DBARTS_USE_STUBS`.

Prose: the contract block's error paragraph is rewritten from "An Rf_error means THIS call is wrong ... longjmps" to the status
rule and the table above; dec-A32's two meanings become "a value entry cannot fail; every other entry returns a status"; THE
HANDLE gains the C-created origin and its lifetime; the thread paragraph states [Threads](#threads); the RETURN and RAISING
paragraphs of `dbarts_draw_callback` are replaced by the four of [The draw callback](#the-draw-callback) (RETURN, EXCEPTIONS,
R ERRORS, A BINDING), and the WORKER-thread paragraph says to return negative; `dbarts_sampler_run`'s interrupt, callback-stop
and slow-count paragraphs follow the statuses, the stop paragraph naming STOPPED and CALLBACK_FAILED; `dbarts_statusString` and
`dbarts_lastErrorMessage` document the fallback; `dbarts_sampler_predict` drops "main-R-thread-only ...
R_alloc"; `dbarts_sampler_getLatents` gains frequency weights' sentence on a frequency row's latent (dec-B564's note); every
per-entry "raises" becomes the status it returns; the spec's block states invariants 1 to 7 of the forward-compatibility
section, the reserved family values and the zero-tail rule. About +600/-180 lines in the header.

The hash: [`dbarts_foldLayout`](../../src/C_interface.cpp) gains the spec's size and every field's name and offset in pointer
units, as `dbarts_results`, and `DBARTS_STATUS_LIST` (values and words) and `DBARTS_WEIGHTS_TYPE_LIST` join the enum fold;
`dbarts_statusString` enters the hash as a list entry like any other; both tokens are re-baked
([`dbarts_apiSignatureToken`](../../src/C_interface.cpp) and the literal). Several by-value doubles in the spec make its offsets
differ by platform in pointer units, as `dbarts_draw`'s do, so the spec folds as `dbarts_draw` does: field count, and doubles by
name, position and width. tools/check-api-hash.sh is unchanged and still skips (no release tag).

Development binaries: a stan4bart or treatSens binary built against today's header still loads against the new dbarts, since
the version pair stays 1/0 and neither defines `DBARTS_REQUIRE_EXACT_ABI`; it ignores every status, and its two-argument
`numTrees` call leaves the new entry writing through whatever the third argument register holds. Only shipped versions are owed
compatibility, so nothing guards this; S4's ports merge with the main feature, and any development install of either is rebuilt.

## Refusals

R's words do not change: every refusal an R user can reach is raised by the R parse before the spec is filled, or by the shared
stage with the text it has today, and the tinytest pins hold them. New words reach a C caller of `dbarts_sampler_create` or
`dbarts_sampler_spec_init` only, read back through `dbarts_lastErrorMessage`, each prefixed with the entry's name (R7):

- `spec is NULL`; `spec.structSize is 0 - set it to sizeof(dbarts_sampler_spec) (e.g. dbarts_sampler_spec s =
  DBARTS_SAMPLER_SPEC_INIT, or dbarts_sampler_spec_init(&s, sizeof s))`; `structSize %zu is smaller than any dbarts_sampler_spec`;
- `the response is NULL`; `the response must be finite`; then the support rule's existing words, family by family;
- `family 0 (AUTO) is not a family; state a dbarts_family from GAUSSIAN to STUDENT` (INVALID_ARGUMENT); `this build creates no
  multinomial sampler from C` and `family %d is not one this build creates` (UNSUPPORTED, for 8 and every value past it);
- `the description is larger than this build knows and byte %zu past its end is not zero` (UNSUPPORTED);
- `weightsType %d is not a dbarts_weights_type`; the frequency and precision refusals in frequency weights' words;
- `a gaussian, Student-t or aft sampler needs a finite positive sigmaEstimate`; `sigmaPriorDf and sigmaPriorScale must be
  finite and positive`;
- `numChains must be at least 1`; `printEvery must be at least 1`; `no seed: set spec.hasSeed and spec.seed, or
  spec.chainSeeds`;
- the source refusals of `translateSource` with this entry's name, the test source's under `xTest`;
- `block tree counts must sum to numTrees`; and the engine's own creation refusals (leaf covariates, variance forest, monotone
  with a proposal mix, storage single), worded as the engine words them, with the entry's prefix.

The run's callback failures (CALLBACK_FAILED), under `dbarts_sampler_run`'s prefix: `the draw callback for chain %zu returned
%d`; `the draw callback for chain %zu threw: %s`; `the draw callback for chain %zu threw an exception not derived from
std::exception`. A deliberate stop (STOPPED) is not a failure and writes no message. Every status's fixed words are the
`dbarts_statusString` table above, which is also the text of a failure that left no message.

## RNG class

Neutral on the R build, by construction and by gate:

- S1 changes which macro spells each math call, not which function runs: under R every `ext_` name is the Rmath symbol the code
  calls today. S2a, S2b and S3 move code without changing what is computed or in which order the generators are touched; the
  one place order could move, the chain seeds, keeps both of today's paths (one seed through a seed generator, or seeds drawn
  from R's stream) and passes them through the spec unchanged.
- Proven by: the equivalence trio bitwise with `--strict-coverage` against the baselines the MANIFEST names as current at the
  branch's base, counting the "identical draws (same RNG stream)" lines ([Gate hygiene](README.md#gate-hygiene)), on the
  reference build and again on the shipped build; the four seeded-drift snapshot files on the reference build; the full tinytest
  suite in one process, every message pin unchanged; tests/cpp; and the x86 leg, recording the base tip on the box and comparing
  --bitwise there (the arm64 baselines are cross-host only on x86). A new arm proves the two creation routes agree: under R, a
  sampler created by `dbarts_sampler_create` from the spec the R parse would fill, and the R-created one, give the same draws bit
  for bit (consumer.c).
- The standalone build. The critic's probe compiled R 4.6.2's nmath closure with clang at R's flags (-O2, contraction on) and
  compared it with this R's Rmath on 17,304 points over the six routines and their flag combinations: every value bitwise equal.
  With contraction off, qnorm differed by up to 6 ulp, qt 10 and qchisq 32. So the copy is compiled with R's own C flags on the
  reference host, and gated there bitwise against Rmath on that grid (a tests/cpp suite in the libR-linked target), and the
  standalone consumer's draws equal the R-created sampler's bit for bit at the same spec and seed. Off the reference host the
  grid is gated at a relative difference of 1e-13 and the consumer's draws by the equivalence harness's statistical z.

## Ruled

The maintainer ruled every open call of this plan on 2026-11-03, and added two rulings it now builds (dec-B579 and the words
clause of dec-B578); the numbers are the plan's open-call numbers, which the ledger's records cite:

1. dec-B571: in the shipped header only an entry that can fail returns a status; the counts and the other accessors return their
   value, and numTrees moves to a status since its forest index can be out of range. The maintainer: "For the host-neutral
   decision, use (a)."
2. dec-B572: a C caller reads a failure's message from `dbarts_lastErrorMessage()`, the last failing call's on the calling
   thread; a capability probe that fails overwrites an earlier message, which the header states. The maintainer: "For
   host-neutral decision 2, use (a)."
3. dec-B575: C creates the middle of the model in 1.0 - every family, constant leaves, the standard tree prior and DART, weights
   and their reading, offset, test predictors and run settings - and the 1.0 shape must already admit several forests and
   multinomial as an append. The maintainer: "You can do (b), but it would be a good time to make sure that whatever is
   necessary for (c) will be forwards compatible. Basically, make sure we're ready to do (c) when the time comes." The plan
   answers it in [Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial).
4. dec-B582: the description is one struct the caller fills, led by `structSize`, with `DBARTS_SAMPLER_SPEC_INIT` for C and
   `dbarts_sampler_spec_init(spec, sizeof *spec)` for bindings that cannot expand a macro; its field order is fixed at 1.0 and
   grows only at the end, the several-forest append coming as an array of element structs at a stated stride, with an element
   initializer. The maintainer: "I think a struct with a size field."
5. dec-B577: the creation entry refuses a description with no seed, saying so, and a binding fills the seed from its own
   environment's stream. The maintainer: "Use (a)."
6. dec-B578: `DBARTS_R_CHECK(call)`, backed by one R-only entry, raises a failure as an R error, an interrupt as R's interrupt
   and a time limit as its error; "the R helper only chooses R's condition and forwards text it reads through the neutral
   calls, and writes none of its own." The maintainer: "For the host-neutral engine call, you can use (a), but if there is any
   translation of error codes into strings, I assume that that might also want to be platform independent too, right?"
   Built in [The header diff](#the-header-diff).
   - dec-B579, from that question: `const char* dbarts_statusString(int status)`, fixed words per status, which cannot fail,
     returns fixed words for an unknown status, and is the fallback of a binding and of `DBARTS_R_CHECK` when a failure left
     no message. The maintainer: "Yes, makes sense. Do (a)." Built in [Statuses](#statuses).
7. dec-B581, revising the call: the callback's return sign decides - 0 continues, a positive return stops on purpose and the run
   returns STOPPED, a negative return returns CALLBACK_FAILED with the value in the message - and a C++ exception escaping a
   callback is caught at the call on every thread and in every build and returned as CALLBACK_FAILED with its `what()`; an R
   error raised in the callback on a single-threaded run under R stays R's jump. The maintainer: "Use your revised
   recommendation." Built in [The draw callback](#the-draw-callback).

## Calls made

1. The vendored set is the closure of all six routines the engine and the support library call, not the four TODO names: qt
   ([`forestMapEntry`](../../src/bartcore/combiner.hpp)) and qchisq at 0.5
   ([`ChiKHyperprior::medianOfChi`](../../src/bartcore/model.hpp)) arrived after the ruling and are read at creation and at live
   restatements (setLeafPrior, and coefficient-prior-live's `$setForest(amplitude = )`). Alternative: compute them on R's side as
   the residual prior's quantile is, which cuts the copy by about 4,000 lines but adds an R-resolved number to every live
   restatement entry and makes a pure function of the prior the caller's. dec-B85's words, "a copy of R's routines kept in the
   support library", fit the copy.
2. The switch is one preprocessor macro, `DBARTS_HOST_STANDALONE`, that only a standalone build defines. Alternative: a configure
   option, which the R build does not need and Windows' checked-in Makevars would have to mirror.
3. The copy keeps R's file layout unedited and renames symbols by a force-included prefix header. Alternative: edit the names in
   place, which makes the next refresh from R a diff to hand-merge.
4. Everything only a standalone build compiles (the copy, io_standalone.c, host_standalone.cpp, the config headers) lives in a
   build-ignored top-level standalone/ directory. Alternative: keep it in src/ and set `OBJECTS` in Makevars.in and
   Makevars.win, then name R Core in DESCRIPTION (`Authors@R` or a `Copyright:` field) for the copy the tarball would carry;
   more to maintain for code no R user runs. Nothing in the R build needs the copy.
5. The standalone build's math gate is bitwise at R's own flags on the reference host and a 1e-13 relative floor elsewhere, the
   copy pinned to R 4.6.x. Alternative: a fixed ulp bound everywhere, which the probe shows is either too loose to catch a wrong
   flag (2 ulp passes nothing at contraction off, so it fails; 32 passes a miscompile) or too tight for another compiler.
6. The shipped spec is one struct and the R-only parts an internal extension beside it, both read by one build stage. Alternative:
   one internal spec of which the shipped struct is a projection; two homes for the same fields, the drift dbartsSpec() was made
   to end.
7. Every route passes the one translate stage: the engine routes through `buildSampler`, the data handle through the spec's
   predictor part. Alternative: leave the handle and multinomial routes on their own parse, keeping two translate stages against
   dec-B124's one.
8. Refusals are values wherever they can meet an R frame, exceptions only below `buildSampler` and the engine, each caught once
   per route before any R call. Alternative (the first draft): a thrown `Refusal` everywhere, which needs a catch inside about 60
   R bridge call sites' frames and fails silently (an exception through `unwindProtect`) wherever one is missed.
9. Exception types map to statuses as in [Statuses](#statuses), with FAILED for what is neither an argument nor an allocation.
   Alternative: name -6 INTERNAL_ERROR, which mislabels the prior-draw exhaustion, a data and prior condition, as a bug.
10. Predict's join rethrows in the failing worker's own category with today's text. Alternative: leave it a `runtime_error`, so
    a C caller sees FAILED for an allocation failure.
11. Two handles may be used at once from two threads, and `predictPartition` becomes `thread_local`. Alternative: keep predict
    main-thread-only under every build; free, but "the main thread" means nothing to a C program, and a binding predicting with
    several samplers on a thread pool would race on one global without a word in the header.
12. When chains of one run end differently, a callback failure wins over a deliberate stop and the lowest failing chain index is
    named. Alternative: report whichever chain's callback returned first, which makes the status and its message depend on
    thread timing.
13. The poll's consumed error is INTERRUPTED with the error's message and a note on the thread's record, so the R raise restores
    it as an error. Alternative: a separate status for "an error raised while polling", a code that only the R build can return.
14. The slow-count warning goes through the build's warning routing, raised last as today. Alternative: a status note, which no
    ruling asks for and dec-B167 left open.
15. Standalone messages, warnings and the abort path print to standard output, as dec-A74 says. Alternative: warnings on standard
    error, conventional for C tools but against the ruling's words.
16. The return type is the enum `dbarts_status`, so a status entry's type says it is one. Alternative: `int`, as the swap status
    is today; the enum costs nothing on R's platforms, none of which builds with short enums.
17. `dbarts_sampler_getSigmas` is named after `$getSigmas()` and refuses where R returns NULL. Alternative: write NaN on a
    heteroscedastic sampler; a status is how this header says "no such channel" (getLatents).
18. A C-created handle is freed whole by `dbarts_sampler_destroy` (dec-B112). Alternative: a second entry, `dbarts_sampler_free`,
    which keeps "a second destroy is a no-op" for both origins at the cost of an entry whose only use is to differ.
19. The verbose creation summary stays in the R bridge (it prints the residual prior's quantile, which the spec does not carry); a
    C-created sampler prints only run progress. Alternative: move it into the build stage, which changes R's printed summary
    unless the quantile rides the spec as well.
20. The specification states a family; AUTO is refused at C creation. Alternative: dispatch on the response shape as R's empty
    family string does, which consumer-spec-surface.md records as silently wrong for aft and logistic.
21. Student-t is `family = DBARTS_FAMILY_STUDENT` with `residualDf`, mirroring the enumeration's existing meaning. Alternative:
    GAUSSIAN with a finite `residualDf`, as the engine's options hold it, which gives the STUDENT enumerator no use at creation.
22. The spec holds the resolved weights reading, two values. Alternative: a third value, "the family's own", which would put a
    second resolution of dec-B564's NULL in the engine beside R's.
23. A seed is stated by `hasSeed`, every 32-bit seed valid. Alternative: 0 as "no seed", which takes a valid seed from R users'
    `control@seed` path.
24. The standalone build is a tests/cpp recipe with checked-in config headers, not a shipped library or a CMake project.
    Alternative: a top-level build for a libdbarts, a maintained surface with no consumer yet; or a configure run, which needs R.
25. The two static_asserts that no shape or source field owns storage stay, and `Refusal` gains one: the R parse still raises
    around them.
26. S1 is its own feature; the rest is one branch merged once, the header change landing whole in S3 so stan4bart and treatSens
    rebuild once. Alternative: S1 on the main branch, which holds a byte-neutral change with no dependency behind forest-door.
27. Under several forests the 1.0 tree prior fields are not read; every forest states its trees in its element, while the 1.0
    leaf prior fields keep their one meaning, the model's leaf prior the sharing forests divide, as forest-door built it.
    Alternatives: (i) the 1.0 tree fields are forest 1's, the array holding forests 2 onward, which makes the array's first
    element differ from the rest and puts one forest's settings in two places; (ii) they are defaults each element inherits
    where it leaves a field unset, which needs an "unset" value for every element field and makes a forest's value depend on
    two places, the drift dbartsSpec() ended.
32. The status words live beside the values in `DBARTS_STATUS_LIST`, so the table is generated and a status cannot be added
    without its words. Alternative: a separate switch in the library, which can fall behind the list.
33. The thread's message record is overwritten by every failure, empty when none could be formed. Alternative: leave the last
    message in place when a failure cannot write one, which shows a stale, unrelated message beside the new status.
34. The callback catch sits in the build's callback unit, one wrapper for the inline call and the worker trampoline.
    Alternative: a catch in each of the engine's call sites, several places to keep alike, and an engine that would then know
    the callback's failure text.
28. A nested array versions by one stride on the parent. Alternative: a `structSize` in each element, which repeats one number
    per element and lets elements of one array disagree, a state with no meaning.
29. Forward compatibility rests on "zero means absent" and on refusing a nonzero byte past the known size, with no placeholder
    fields at 1.0. Alternative: ship `numForests`, `forestStride` and `forests` now, refused unless zero; it fixes their names,
    types and order before the element is designed, for a guarantee the zero-tail rule already gives.
30. K is the 1.0 field `numCategories`, read under ORDINAL at 1.0 and MULTINOMIAL later. Alternative: a separate appended K for
    multinomial, two fields for "the number of response categories".
31. The R path builds several-forest and multinomial samplers from extension records laid out as the sketched append from S2b,
    so the append exposes a tested route. Alternative: keep the R path's own records and translate at the append, a second
    layout and a second test of it later.

## Open calls

None. Every call this plan raised is ruled ([Ruled](#ruled)).

## Slices

### S1 (its own feature, branch host-neutral-math): build-time math and printing

Steps: the six ext_ names in stats.h and the slim nmath.h; model.hpp's fourteen calls and random.c's two declarations onto them;
standalone/nmath with the copy, its prefix header, the reduced R headers and README; .Rbuildignore's `^standalone$`; io.h's two
mappings and standalone/io_standalone.c; randomBase.c's R-only parts guarded; misc/io.c, misc/io.h and their five installers
(R_interface.cpp, tests/cpp's main.cpp, three benchmark kernels) removed; standalone/config and the drift check's extension; a
tests/cpp `nmath` suite comparing each copy with Rmath on the probe's grid (run by the default libR-linked target, bitwise on
the reference host); a `standalone-math` target that compiles external, misc and the copy with the switch into its own
directory and links no libR.

Gates: tests/cpp (both targets); full tinytest; the trio bitwise; the four snapshots on the reference build; R CMD check
--as-cran from a clean tarball, which must contain no standalone/ file; sanitizers on tests/cpp; mutation-battery
verify-anchors; the x86 leg. Diagnostic, not a gate: `objdump -d` of the installed dbarts.so before and after, expecting only
R_init_dbarts and the removed misc symbols to differ.

Budget: ~350 lines written, ~6,000 copied, ~400 reduced headers. Agent ~5, review ~1.5, gates ~1.5.

### S2a: refusals as values, one translate stage (bridge)

Steps: `Refusal` and `raiseRefusal`; the sixteen shared validators, and the flat file's `translateSource`, `validateTestSource`
and `refuseTestMissingness`, returning refusals, their about 60 R sites and the flat file's wrapped (createDataHandle and the R
test-predictor setters call the last three inside `unwindProtect` bodies, outside any route's catch); the translated source,
which owns storage, held as an init-capture wherever it lives in an R frame, by `unwindProtect`'s own rule; R_interface_bartcore_common.hpp split into the R-free core and the SEXP-facing half; `translateSource` moved
into the core, used by the R parse's views and the data handle, carrying the engine's per-column channels (the shipped field
waits for S3); the engine routes' catch turned into a refusal value.

Gates: as S1, plus the R-loaded ASan run over the creation and mutation suites; every message pin in tinytest unchanged (the
whole suite in one process).

Budget: ~900 lines. Agent ~12, review ~5, gates ~3.

### S2b: the plain-C specification (bridge)

Steps: `dbarts_sampler_spec`, its INIT and `DBARTS_WEIGHTS_TYPE_LIST` in dbarts.h (no entry yet); src/spec_internal.h with
`SpecExtension`; src/sampler_build.cpp with `buildSampler`, taking over the engine-building half of createHolder,
optionsFromParsed, the two multinomial routes and createFromHandle; the R parse filling spec and extension (the parse functions
keep reading SEXPs; their outputs become spec fields); seeds through `hasSeed`/`seed` or `chainSeeds`; the weights reading; a
tests/cpp suite that fills every spec field with a distinct value and checks the engine options built from it field by field,
so a field the build forgets fails by name; the extension's public-shaped part (`ForestRecord`, `MultinomialRecord` and the
one-forest groups 1.0 leaves out) laid out as the sketched append, the several-forest and multinomial R routes reading it; the
1.0 reading rules (zero tail, reserved families, fields read only under the families that name them).

Gates: as S2a.

Budget: ~1,600 lines (header ~200, extension ~250, build stage ~650 moved, parse rework ~400, tests/cpp ~250). Agent ~18,
review ~7, gates ~4.

### S3: statuses, the creation entry, the sigma getter, the header (bridge, flat entry file, tests)

Steps: `DBARTS_STATUS_LIST` with its words and `dbarts_statusString`; the list as in [The header diff](#the-header-diff);
`dbarts_sampler_create`, `spec_init`, `getSigmas`, `lastErrorMessage` with the overwrite-on-every-failure record; C_interface.cpp
off R (owned scratch, statuses, the thread's record); `predictPartition` thread_local; predict's join keeping its category;
src/host_r.cpp and standalone/host_standalone.cpp, with the callback catch for the inline call and the worker trampoline; the
engine's three-way stop record (continue, stop, failure with chain and value or text) in place of the one stop flag, and the
lowest-failing-chain rule; the holder's origin and destroy; `dbarts_R_raiseStatus` and `DBARTS_R_CHECK` (dec-B578), forwarding
only neutral text; `denseChannels` in the shipped source; the hash folder and both tokens; the header's prose, the callback
contract rewritten. Tests: consumer.c and test-capi.R rewritten for statuses (each entry's failure reached and its
status and message read; INIT plus data and seed against `dbarts()`'s draws bitwise; a C-created sampler against the R-created
one bitwise; `denseChannels` with packed blocks against the dense source; getSigmas on two chains and its UNSUPPORTED on a
heteroscedastic sampler; `DBARTS_R_CHECK` raising an error, through the interrupt test hook an interrupt, and under
`setTimeLimit` the time limit's error, as test-monotone.R pins for R's run; `DBARTS_R_CHECK` on a failure whose message is
empty raising `dbarts_statusString`'s words; two samplers predicting on two threads; a C-created handle destroyed; and the four
callback paths, each run on a gaussian sampler of two chains: a positive return on chain 1 (STOPPED, no message, results
discarded, a second run valid from where it stopped), a negative return (CALLBACK_FAILED, the message naming the chain and the
value), a `std::runtime_error` thrown on a worker (two threads, CALLBACK_FAILED with its `what()`, the process alive), and the
same thrown inline (one thread, the same status and message); and, under R on an inline run, an `Rf_error` from the callback
still reaching R as its own error), the five handshake arms kept; a tests/cpp unit for `dbarts_statusString`, every listed status
against its words and an unlisted one against "unknown status"; tests/cpp's `standalone-consumer` target: a C99 file including the prototype
view, creating gaussian, probit and Student-t samplers, running, predicting, reading sigma, refusing, destroying, built and
linked with no R, `nm` showing no undefined R symbol, and its draws compared with the R-created sampler's at the same spec and
seed; the `spec-append` target of [Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial)
(the test-only append, the 1.0-sized and the appended consumer, the two-forest arm and the nonzero-tail refusal); the
cpp-tests workflow runs both. Docs: dbarts-embedding.Rd's C section, the architecture's layering and callbacks sections,
NEWS.

Gates: as S2a, plus test-capi.R under ASan on the R-loaded path, tests/cpp's standalone consumer under ASan and ThreadSanitizer
for the two-thread and worker-exception arms, on the x86 leg LeakSanitizer (C-created handles are new owners); the header-compiles target;
`sh tools/check-api-hash.sh`.

Budget: ~2,800 lines (header +600/-180 with the test-only append ~80, C_interface.cpp ~530, build units ~290 with the callback
catch, engine stop record ~40, bridge ~120, consumer.c ~600, test-capi.R ~340, tests/cpp ~580, workflow ~20, docs ~130). Agent
~16, review ~6, gates ~4.

### S4: the consumer ports (stan4bart, treatSens)

Each in its own checkout, on a branch, against a private library of the S3 tip, its CI's dbarts reference set to the branch
while in review. stan4bart: the fallible calls in src/init.cpp (run twice, sampleTreesFromPrior, setOffset twice, setSigma twice,
getLatents twice, predict, setTreeStorage twice, printTrees, setVerbose) wrapped in `DBARTS_R_CHECK`, numTrees through its out
pointer; nothing else changes (the dense helper is source-compatible). treatSens: run three times, setResponse four times,
setOffset, setSigma and setVerbose twice, likewise. Each suite run (stan4bart's tinytest, recompute-exactness and
compare-posterior; treatSens' testthat), R CMD check, and revdep-smoke dispatched on the dbarts branch tip with both branches,
bartCause included, all green. They merge right after dbarts' merge, their dbarts reference back to main.

Budget: stan4bart ~80 lines, treatSens ~50. Agent ~3 each, review folded into S3's reviewer reading both diffs.

### S5: records

TODO's host-neutral entry and python-bindings entry (the engine is R-free below the bridge; what a binding still owns:
ingestion, resolution, its own build unit and the standalone recipe); docs/plans/python-bindings.md's Context; the feature
matrix's reach column for the flat header; INDEX rows; NEWS's 1.0-0 entry for compiled consumers. Doc-freshness only.

## Tests and gates

Per slice above; across both features, run independently of the implementer on each slice's own library, `--preclean` on every
install (the facade's virtuals do not change, but every creation route does):

```
R CMD INSTALL --preclean -l <lib> .                       # and --configure-args=--enable-reference-build for the snapshots
cd tests/cpp && make && ./test_bartcore                   # header-compiles, the libR suite, the nmath grid
make standalone-math standalone-consumer && ./build-standalone/standalone_consumer
nm build-standalone/standalone_consumer | grep -E ' U _?(Rf_|R_|Rprintf|unif_rand)'   # expect nothing
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts", lib.loc = "<lib>")'
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare --bitwise --strict-coverage benchmarks/baselines/<current>.rds
R_LIBS=<lib> Rscript benchmarks/R/bcf-equivalence.R compare --bitwise benchmarks/baselines/<current>.rds
R_LIBS=<lib> Rscript benchmarks/R/multinomial-equivalence.R compare --bitwise benchmarks/baselines/<current>.rds
  # every scenario "identical draws (same RNG stream)", none skipped, no "max |z|"
# the four test-reproducibility-*.R files on the reference build, run_test_file each
Rscript -e 'lintr::lint_package()' && air format --check . && Rscript tools/check-rc-codoc.R . \
  && Rscript tools/check-win-drift.R . && Rscript tools/check-doc-freshness.R . \
  && Rscript benchmarks/R/mutation-battery.R verify-anchors
R CMD build . && tar tzf dbarts_*.tar.gz | grep standalone   # expect nothing
R CMD check --as-cran <tarball from a clean copy>
sh tools/check-api-hash.sh                                 # "no release tag, skipped"
```

plus sanitizers on tests/cpp and on the R-loaded path ([Gate hygiene](README.md#gate-hygiene)), and the x86 leg before landing
each slice (tests/cpp plain and with ASan, UBSan and LeakSanitizer, tinytest, the trio recorded on the base tip there and
compared --bitwise, the nmath grid at the relative floor).

## Mutants

The reviewer of each slice runs these; each must fail the named gate, and is reverted with a `touch`:

1. S1: one model.hpp pnorm site's lower-tail flag flipped - the probit and aft equivalence scenarios.
2. S1: a constant in the copied qnorm changed - the nmath grid.
3. S1: the copy compiled with `-ffp-contract=off` - the nmath grid's bitwise gate on the reference host.
4. S1: `ext_qnorm` mapped to the copy in the R build - record which scenarios move; with the copy exact at R's flags, none
   should, so this proves the gate reads the R build's symbol by `nm` rather than by draws.
5. S2a: one validator's refusal dropped unconsumed at an R site - the `[[nodiscard]]` build error, and that refusal's tinytest pin.
6. S2a: a refusal thrown, not returned, from a validator called inside an `unwindProtect` body - on arm64 macOS the exception
   crosses R's frame silently (the critic's probe: no crash, the cleanup never runs, the heap closure leaks), so the gate is
   LeakSanitizer on the x86 leg, and a test-build counter of `unwindProtect` cleanups against bodies entered, read by a
   tinytest hook, catches it on every host.
7. S2a: a dense column read in the block its kind selects instead of through the channels - consumer.c's packed-blocks arm.
8. S2b: one spec field (perturbProbability) not copied into the options - the field coverage suite, and the perturb scenario.
9. S2b: `chainSeeds` ignored in favor of `seed` - every unseeded tinytest snapshot and the C-versus-R identity arm.
10. S2b: one R-worded refusal moved behind the spec fill, reworded by the shared stage - its tinytest pin.
11. S3: setSigma's refusal returning OK - consumer.c's status arm.
12. S3: the thread's message not written on one entry - consumer.c's message arm.
13. S3: `dbarts_R_raiseStatus` raising an error for an interrupt, or an interrupt for a time limit - test-capi.R's two arms.
14. S3: `predictPartition` left process-wide - the two-thread arm under ThreadSanitizer.
15. S3: predict's join rethrowing an allocation failure as FAILED - a tests/cpp unit whose worker throws `std::bad_alloc`.
16. S3: getSigmas writing chain 0 to every slot - its two-chain arm.
17. S3: two spec fields swapped in the header - the build's static_assert against the re-baked hash.
18. S3: destroy on a C-created handle releasing the engine only - LeakSanitizer on the x86 leg.
19. S2b: the build reading a field past a 1.0-sized `structSize` instead of treating it as zero - the `spec-append` target's
    1.0-sized consumer (its draws move or it is refused).
20. S2b: the zero-tail check dropped, a nonzero appended byte accepted - the `spec-append` target's refusal arm.
21. S2b: the 1.0 tree fields read under `numForests` > 0 - the two-forest arm against the BCF scenario.
22. S2b: MULTINOMIAL accepted at C creation - consumer.c's reserved-family arm.
23. S3: a negative callback return mapped to STOPPED - consumer.c's negative-return arm.
24. S3: a positive return mapped to CALLBACK_FAILED - the positive-stop arm (and `DBARTS_R_CHECK` raising for a deliberate stop).
25. S3: the worker trampoline's catch removed - the worker-exception arm (the process ends).
26. S3: the inline catch dropping `what()` - the inline-exception arm's message check.
27. S3: the first chain to return named instead of the lowest failing one - a two-chain arm where both fail, run 20 times.
28. S3: a failure leaving the previous message in place when it cannot form one - a tests/cpp unit that arms an allocation
    failure in the message path and reads an empty message, then `DBARTS_R_CHECK` raising the status words.
29. S3: `dbarts_statusString` missing one status's words - its tests/cpp unit (and the X-macro makes the omission a build error).
30. S3: `dbarts_R_raiseStatus` writing its own text for UNSUPPORTED - test-capi.R's arm comparing the R message with the neutral
    one.

## Interaction with the other pending plans, and order

- S1 depends on nothing; it lands when green, serial with coefficient-prior-live's merge (engine features merge serially),
  rebasing onto whichever lands first (two qt lines in combiner.hpp are its only overlap).
- forest-door (slices 1 to 3, with setmodel-parameters' slice 4) LANDED on main; this plan was rebased onto 8fce28568. What it
  moved: `ForestSpec` gained `share`, documented "No caller states it", the group derived by `forestSharers`, so the share left
  the sketched public element for the extension's R-only part; the leaf prior of a several-forest chain is the total the sharers
  divide, so the 1.0 leaf prior fields keep one meaning under several forests (invariant 4 and call 27 rewritten); the bridge
  gained the leaf-prior install on several forests and `setForestTreePrior`, which S2a's validator count and S2b's parse rework
  take as they stand. The math touchpoints, the three `runtime_error` sites and the two process-wide channels are unchanged.
- coefficient-prior-live lands before S2a starts. It reworks the live forest writes in the same bridge and engine files, and
  S1's qt and chi-median copies serve its live restatements.
- Frequency weights (dec-B563 to dec-B566) and the binomial response (dec-B567) are ruled; their design note and plan are not
  built. The spec carries the resolved weights reading from S2b whichever lands first (dec-B564: "the shipped header gains one
  field in dec-B85's change"); if frequency weights' engine work lands later, its plan widens what the build takes for
  FREQUENCY and adds the getLatents sentence if S3 has merged. The binomial response enters the spec in the form its design
  gives the engine (likely a per-row trial count beside y under logistic): in S2b if it lands first, else appended by its own
  work before 1.0-0, a hash move that rebuilds no one, since neither consumer creates through C.
- multinomial-doors' C creation path stays post-release (dec-B160, dec-B575); this plan ships the shape it will append to and
  the R route it will expose.
- Frequency weights' dec-B580 (a 0 row's predictive draws NA on every family) asks nothing of the spec: the row is drawn and its
  draws set to NA after, so the stream holds, and the spec states only the weights and their reading.
- The sister re-verification at the freeze follows this arc (TODO), with stan4bart's sanitizer job dispatched against the final
  header (dec-A36).
- Nothing else touching src/ merges while S2a, S2b or S3 is in review.

## Cost

About 5,650 lines written here, 6,000 copied and 400 reduced headers, none of them in the tarball; ~130 in the two consumers;
~93 agent-hours estimated, 140 to 185 by the record of plans here, S2a and S2b most of it. The forward-compatibility work is
about 300 of those lines and 4 of those hours; the callback statuses and the status words about 300 and 3.

## Not verified

- The probe's bitwise match held for one compiler, one R (4.6.2) and one host; the x86 leg's agreement is gated at the floor, not
  measured.
- That `objdump` of dbarts.so before and after S1 differs only where named; relative call offsets may make the diff noisy.
- What a failed run (OUT_OF_MEMORY, FAILED) and a failed sampleTreesFromPrior leave in the chains beyond "valid to read and
  destroy"; S3's tests pin what the table states, and the header promises no more than they show.
- That no other process-wide mutable state is on the predict path beyond the two channels named; a ThreadSanitizer arm in S3 is
  the check.
- How coefficient-prior-live leaves the forest records, which S2b's extension must mirror; forest-door's, landed, were read on
  main at 8fce28568 for `ForestSpec` and `forestSharers` only, not every bridge record.
- Where the engine's sweep loop can carry a chain index and a callback value to the stop record without a cost per sweep; the
  ruling's "about 40 lines" is its estimate, not measured.
- The binomial response's engine form, not yet designed.
- Whether the engine's several-forest chain takes DART, split probabilities and the k prior per forest; the sketch leaves them
  sampler-wide and the append refuses what the factory refuses.
- stan4bart's and treatSens' call sites were read on their main and master; their 1.0 compat branches may differ.
- That no R-hosted consumer other than the two links dbarts.h; the reverse LinkingTo list on CRAN was not queried.
