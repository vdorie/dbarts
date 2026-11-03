# host-neutral-engine: a plain-C specification the engine is built from, build-time math and printing, a status from every entry point, and a creation entry in dbarts.h

Status: PLANNED 2026-11-03 on main (c3d268c8f), for dec-B85. Not built. Rulings, quoted where they settle a call: dec-B85 (the
staged route), dec-B84 (the header is pure C), dec-B119 (a status from every entry point, the unwind to the entry point
staying), dec-A74 (printing chosen at build time, the run-time print pointers removed), dec-A73 (a plain-C sigma getter rides
this header change), dec-B124 (the specification's predictor part is `dbarts_predictor_source`, its dense channel stated per
column), dec-A50 (the engine's R touchpoints), dec-A20 (creation's refusal channel), the errors-unwind ruling of 2026-09-27
(dec-B119's record, "Errors unwind instead of jumping"), with dec-B86 (the header ships what a consumer calls), dec-B87 and
dec-A33 (copy-on-set), dec-B112 (destroy stays; "A C consumer should have a path to invalidate a C object"), dec-B56 and dec-B57
(the list, the hash, per-symbol lookup), dec-B60, dec-A139 and dec-B424 (the version pair stays 1/0 before the release, the hash
moves), dec-B111 (the hash-to-minor CI assertion), dec-B160 and dec-B163 (multinomial run and predict through the flat API, its
creation a post-release door), dec-B186 (the C consumer test compiles and runs on CRAN), dec-A36 (consumers check the boundary
under sanitizers themselves), dec-A32 (a returned integer is a value or a capability status), dec-C01 (full symbol names) and
dec-C04 (the support library's print pointers). The runs behind it: reads of the tree at c3d268c8f, of stan4bart main (27ef163)
and treatSens master (268b2db), and of frequency-weights (da210ce2d); no prototype was built.

agent: S1 one opus implementer (effort high); S2 one opus implementer (effort high), the riskiest slice; S3 the same or a fresh
opus implementer; S4 one sonnet implementer per consumer, in that consumer's own checkout; S5 sonnet. One opus review per engine
slice (S1, S2, S3) that runs the mutants below, the same reviewer checking the fix round unless the fixes change engine behavior.
Tier "Changes draws" for S1 to S3 ([Process by risk](README.md#process-by-risk)): nothing is meant to move a draw, but every
creation route and every flat entry is rewritten, and the gates are the draw-changing tier's.
rng: neutral, bit for bit, on the R build: the equivalence trio, the four snapshot files, tests/cpp and the full tinytest suite
do not move, and no R message changes. The standalone build is new and owes no bitwise match to R; see [RNG class](#rng-class).
window: before the 1.0-0 submission (dec-B85: "post merging into main, but before 1.0 since it might change a number of packages
and is an interface change"). After forest-door and coefficient-prior-live land on main; one branch, host-neutral-engine, five
slices, merged once (dec-B494), the consumer ports merging right after it.
budget: ~5,000 lines written in this tree plus ~6,000 copied unedited from R's math library: S1 ~350 (+6,000 copied), S2
~1,900, S3 ~2,300, S5 ~150; S4 ~80 in stan4bart and ~50 in treatSens. Agent time, estimated and not measured: S1 ~7.5, S2 ~17,
S3 ~17, S4 ~6, S5 ~1, rebases onto main after the pending plans ~2: ~50. Plans here run 1.5 to 2 times low, so 75 to 100.

## Summary for the maintainer

You ruled in September that before 1.0-0 the engine stops depending on R below the layer that talks to R, so that a program
other than R, a Python package for instance, can build and run a sampler through the one C header dbarts ships. This plan builds
that in five steps.

First, the math and the printing. The engine and its support code call six of R's numerical routines (the normal distribution's
density, probability and quantile; the t distribution's density and quantile; the chi-squared quantile) and R's console. Under R
they keep calling exactly those, so every draw an R user makes stays bit for bit the same. Every other build uses a copy of R's
own routines kept in the package and writes to standard output. Your ruling named four routines; the engine has since come to
use two more, the t quantile and a chi-squared median, both read only when a sampler is built or its prior is changed, and this
plan copies them too rather than asking every caller to compute them. The two unused printing hooks go.

Second, the specification. Today the code that turns R's control, model and data objects into an engine reads those objects
directly. It will instead fill one plain C description of the sampler - its data, family, priors, constraints and run settings,
each already resolved - and the engine is built from that description alone. R still does everything R does now: reading the
formula, choosing the family from the response, estimating the starting sigma, turning the residual prior's quantile into a
scale, parsing the prior vocabulary, and naming columns and forests in its messages. Nothing an R user sees changes.

Third, the header change, made once so the two compiled packages rebuild once. A program can create a sampler from that
description, every entry that can fail returns a status in place of raising an R error, the failure's message can be read
back, each chain's current sigma can be read, and the predictor description says column by column which block a column's values
sit in, as you ruled. A sampler created this way has no R object behind it and is freed by the destroy entry.

Fourth, stan4bart and treatSens are changed to check those statuses; with the helper recommended below, an R user of either
sees the same errors and the same interrupt behavior as today. Fifth, the records.

Six choices only you can make, each a shape the two compiled packages and any future binding would code against: whether the
accessors that cannot fail keep returning their value; where a failure's message is read; whether the C description covers
one-forest models only for now; whether the description is one struct or built up by calls; whether creating from C requires a
seed; and how an R package that calls dbarts' C code turns a status back into R's error or interrupt. The recommendations keep the
header small and keep every behavior an R user sees.

What it costs: about 5,000 lines here and a copied 6,000 of R's math code, and roughly 75 to 100 agent-hours. It has to wait for
the forest door and the live coefficient prior, which rework the same creation code, and it should go before frequency weights,
which adds one field to the new description.

## Goal

The engine, the support libraries and the flat C entry file compile and run with no R header and no R library: a tests/cpp
target builds them into a plain C consumer that creates, runs, predicts, reads sigma and refuses through statuses alone, and no R
symbol is in its binary. Under R every draw, every message and every printed line is what main produces at the base of this
branch. Every creation route, R's five and the new flat one, fills one plain-C specification and passes one translate-and-validate
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
  factory call, the test predictors and the verbose summary ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)). Four
  more routes build engines: [`createMultinomialCountsHolder`](../../src/R_interface_bartcore.cpp),
  [`createMultinomialDataHolder`](../../src/R_interface_bartcore.cpp), [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp)
  and [`bartcore_createFromHandle`](../../src/R_interface_bartcore.cpp) (xbart's folds).
- The flat entries: [`translateSource`](../../src/C_interface.cpp) validates a caller's predictor source (on R's transient
  stack); every refusal is an `Rf_error`; engine exceptions become R errors through
  [`callConvertingExceptions`](../../src/R_interface_bartcore_common.hpp); a run's callback raise travels as an
  [`UnwindJump`](../../src/R_interface_bartcore_common.hpp) the entry resumes; an interrupt is caught by
  [`userInterrupted`](../../src/R_interface_bartcore.cpp) and re-raised by [`raiseInterrupt`](../../src/R_interface_bartcore.cpp);
  a slow monotone count warns through [`raiseSlowCountWarning`](../../src/C_interface.cpp). The handle is the bridge's holder,
  [`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp), which carries no R type.
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
    [`misc_printf`](../../src/include/misc/io.h) and `misc_flushOutput` pointers in misc.
  - in the bridge's own core: `R_alloc`, `R_finite`, `R_UnwindProtect` (the callback leaf, [`ShippedDrawHook`](../../src/R_interface_bartcore_common.hpp)),
    `R_CheckUserInterrupt` (the poll).
  - The residual prior's chi-squared quantile is already the bridge's (`sigmaRawScale`), as dec-B85 has it.
- tests/cpp links libR ("external.a's rng object references the R API even for non-native algorithms") and has no standalone
  target. The rshim the TODO's predecessors name was deleted in c0ddd93e8; nothing in the tree builds the engine without R.
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
  defines `DBARTS_REQUIRE_EXACT_ABI`; neither creates through C. bartCause and bairrtt link nothing of dbarts.
- docs/plans/pure-c-header.md (the last header change, its slices and gates), docs/design/consumer-spec-surface.md (why one
  resolution home), docs/design/r-c-division.md, error-style.md (R7, a flat entry names itself; R8, the R message is canonical
  and a C backstop is worded apart), [Callbacks and errors](../architecture.md#callbacks-and-errors),
  [Reproducibility contract](../architecture.md#reproducibility-contract), docs/plans/python-bindings.md.

## Constraints

- Neutral on the R build: no draw, message, warning or printed line changes; the gates in [Tests and gates](#tests-and-gates).
- The header stays pure C (dec-B84): `gcc -std=c99 -pedantic` and `g++ -std=c++11` on its prototype view with no R include
  path (the tests/cpp header-compiles target). The version pair stays 1/0 and the hash is re-baked (dec-B60, dec-B424).
- No run-time hook for math, printing or interrupts (dec-B85, dec-A74): each is chosen when the code is compiled.
- The R surface is untouched: no R function, argument, return shape or help text changes, except the help's C section
  (dbarts-embedding.Rd).
- Out of scope: a C creation path for multi-forest or multinomial models (open call 3; dec-B160 keeps multinomial's a
  post-release door); an opaque C state buffer (dec-B84, conditional on a consumer); a shipped standalone library or build
  system beyond the tests/cpp recipe; any R wrapper for a C-created handle; the SamplerShape fields that size the bridge's arrays
  (dec-A50 names them; no ruling asks this arc to move them); frequency weights' own work (below).

## What stays with R

The R parse keeps every step that needs R or that R alone resolves, and hands the specification values, never R objects:

- ingestion: the formula and matrix interfaces, `dbartsData`, factors and their levels, sparse and mixed containers, `Surv`,
  the hazard expansion, `na.action`;
- resolution ([`resolveSamplerSpec`](../../R/spec.R)): the family from the response shape (the spec takes a stated family,
  never `DBARTS_FAMILY_AUTO`), ordinal recoding, the weight policy, the starting sigma estimate, offset normalization, the prior
  vocabulary and its evaluation, the per-family defaults;
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
  caller-to-library direction) holding what a one-forest sampler needs, every value resolved. Its fields, grouped (the order
  below is the struct's):
  - data: `x` (a `dbarts_predictor_source`, its typing channel READ at creation, which gives `columnTypes` and `categoryCounts`
    the reader dec-B124 owes them), `y`, `weights`, `offset`, `survivalStatus` (aft, double 0/1 per row as R holds it),
    `xTest` (a second source, 0 rows for none), `testOffset`, `numCuts` and `maxNumCuts` (one or per column), `useQuantiles`;
  - family: `family` (a `dbarts_family`, 1 to 7), `sigmaEstimate`, `sigmaIsFixed`, `sigmaPriorDf`, `sigmaPriorScale` (the
    host's `qchisq(1 - q, df) / df`), `residualDf` (read under STUDENT: positive fixes it, 0 estimates it), `numCategories`
    (ordinal K), `shape` (nbinom: positive fixes r, non-positive estimates it);
  - leaf prior: `leafModel` (a `dbarts_leaf_model`, giving the enumeration its first reader), `k`, `nodeScale`, `namedSd` (NaN
    none), `kIsDrawn`, `kPriorDf`, `kPriorScale`, `leafCovariates` and their count, `gpLengthscales`, `gpMaxLeafSize`,
    `monotoneDirections`, `monotonePrior`;
  - tree prior: `numTrees`, `base`, `power`, `splitProbabilities`, `useDart` and DART's six (`dartA`, `dartB`, `dartRho`,
    `dartAlpha`, `dartUpdateAlpha`, `dartUpdateDelay`);
  - constraints: `forestColumns`, `interactionMaxOrder`, `interactionForbiddenPairs`, `numBlocks`, `blockOfColumn`,
    `blockTreeCounts`;
  - moves: the six probabilities, `perturbWidth`, `levelGibbs` (-1 automatic, 0 off, 1 on);
  - variance forest: `numVarianceTrees` (0 none), `varianceBase`, `variancePower`, `varianceForestColumns`;
  - run: `numChains`, `numThreads`, `numThin`, `keepTrees`, `numSamplesToStore`, `verbose`, `printEvery`,
    `singlePrecisionResidual`, `seed`, `chainSeeds`;
  - engine limits: `categoricalExhaustiveCap`, `testFitParallelCutoff`, `predictParallelCutoff`, `sparseDensityThreshold`.

  `DBARTS_SAMPLER_SPEC_INIT` sets `structSize` and `dbarts()`'s defaults for a gaussian response wherever they are numbers,
  NULL for every vector, and no seed. The header lists the values R resolves per family (probit's node scale, for one), since
  the spec carries no "the family's default" marker: one resolution home, as dbartsSpec()'s design argued.
- `bartcore_bridge::SpecExtension`, internal (src/spec_internal.h, plain C structs and arrays): the re-creation records, the
  multi-forest records (each forest's tree count and settings, sd, multiple, held flag, df, share as forest-door leaves it,
  kind, basis, basis scale), the amplitude glue settings, the multinomial counts, trials and K, the data handle's prebuilt store,
  the creation mask and the host's unit count (`tooFewUnitsIn`). Empty, it is the flat entry's.

One stage, `buildSampler(const dbarts_sampler_spec&, const SpecExtension&)`, in a new R-free src/sampler_build.cpp: validates
the spec with host-neutral words, translates both sources through the one `translateSource` (dec-B124: the flat predict entry,
creation from C, and the R parse's views all pass it), fills the engine structs, makes the chain generators from the seeds,
selects and calls the factory, installs the move record, the test predictors and the test offset, and returns a holder or throws
a `Refusal`. Every one of R's five routes and `dbarts_sampler_create` call it. What createHolder does today in that order moves
there whole; only the SEXP reads stay behind.

### Seeding

The spec takes either `seed` with `chainSeeds` NULL - one generator seeded from it hands each chain its seed, which is the
`control@seed` path today and gives the same draws - or `chainSeeds`, one per chain, used as they stand, which is how the R
parse passes the seeds it draws from R's stream when the control states none. R's draws are therefore unchanged on both paths.
A spec with neither is refused (open call 5).

### Statuses

`dbarts_status`, an ABI enum generated like the others (`DBARTS_STATUS_LIST`), replacing `DBARTS_SWAP_STATUS_LIST`. Zero is
success, a positive value a success with a note, a negative value a failure that touched nothing the caller can rely on:

| value | name | meaning |
|---|---|---|
| 0 | `DBARTS_OK` | done |
| 1 | `DBARTS_SCALE_KEPT` | setResponse or setOffset swapped and kept the transform in force (today's SCALE_KEPT) |
| 2 | `DBARTS_STOPPED` | a run a callback's nonzero return stopped; results and saved trees of the call are discarded |
| -1 | `DBARTS_UNSUPPORTED` | the sampler cannot do this at all, nothing touched (today's capability 0) |
| -2 | `DBARTS_INVALID_ARGUMENT` | this call is wrong (today's `Rf_error` "this call is wrong") |
| -3 | `DBARTS_OUT_OF_MEMORY` | an allocation failed (a tree store, a build) |
| -4 | `DBARTS_INTERRUPTED` | the host's interrupt stopped a run; the sampler as after a callback's stop |
| -5 | `DBARTS_CALLBACK_ERROR` | a draw callback threw a C++ exception; its `what()` is the message |
| -6 | `DBARTS_INTERNAL_ERROR` | anything else the engine threw |

Inside the shared core every refusal is a C++ exception (`bartcore_bridge::Refusal`, a status and a message) thrown where it is
found, so the unwind to the entry point stays (dec-B119) and only the last step changes: the flat entry catches, records the
message on the thread and returns the status; an R `.Call` entry catches and raises `Rf_error` with the message it raises today.
Engine exceptions map by type: `std::invalid_argument`, `std::domain_error` and `std::out_of_range` to INVALID_ARGUMENT,
`std::bad_alloc` and `std::length_error` to OUT_OF_MEMORY, the rest to INTERNAL_ERROR.

Per entry (the infallible ones keep returning their value under open call 1's recommendation):

| entry | returns | can fail with |
|---|---|---|
| `dbarts_sampler_create` | status, handle through an out pointer | INVALID_ARGUMENT, UNSUPPORTED (multinomial, AUTO), OUT_OF_MEMORY |
| `dbarts_sampler_run` | status | INVALID_ARGUMENT (structSize 0), INTERRUPTED, CALLBACK_ERROR, OUT_OF_MEMORY, INTERNAL_ERROR; STOPPED |
| `dbarts_sampler_sampleTreesFromPrior` | status | INVALID_ARGUMENT (a tree whose every draw left an empty leaf), OUT_OF_MEMORY |
| `dbarts_sampler_setResponse`, `setOffset` | status | UNSUPPORTED, INVALID_ARGUMENT (support, finiteness, updateScale refusals); SCALE_KEPT |
| `dbarts_sampler_setSigma` | status | UNSUPPORTED (pinned), INVALID_ARGUMENT (not finite and positive) |
| `dbarts_sampler_getSigmas` (new) | status | UNSUPPORTED (heteroscedastic, as R's getSigmas() returns NULL) |
| `dbarts_sampler_getLatents` | status | UNSUPPORTED (augments nothing) |
| `dbarts_sampler_predict` | status | UNSUPPORTED (blend undefined), INVALID_ARGUMENT (source, offset, empty store) |
| `dbarts_sampler_setTreeStorage` | status | OUT_OF_MEMORY |
| `dbarts_sampler_printTrees` | status | INVALID_ARGUMENT (an index out of range) |
| `dbarts_sampler_setVerbose` | status | INVALID_ARGUMENT (printEvery 0) |
| `dbarts_sampler_numTrees` | status, count through an out pointer | INVALID_ARGUMENT (forest out of range) |
| `destroy`, `setDrawCallback`, `setNumThreads`, the other accessors, the version triple, `dbarts_lastErrorMessage` | value or void | nothing |

The message: `const char* dbarts_lastErrorMessage(void)`, the message of the most recent call on this thread that returned a
negative status, self-prefixed with the entry's name as today (`"dbarts_sampler_run: results.structSize is 0 - ..."`), valid
until the next failing call on the thread; a success leaves it (open call 2).

Three things stay as they are, because a status cannot carry them: a draw callback that RAISES an R error on an inline run is
R's own jump, which the entry resumes after the unwind exactly as now (abandoning it would drop the target of an enclosing
`tryCatch`); the stub view's handshake, consumer code in an R-only view, still raises; and the slow-count warning goes out through
the build's warning routing, last, as now (under R an R warning, which `options(warn = 2)` can still turn into a jump; standalone,
a line on standard error).

### Build-time math, printing and interrupts

One preprocessor switch, `DBARTS_HOST_STANDALONE`, defined by a standalone build and never by the R package's Makevars, so the R
build compiles as it does today.

- Math: src/include/external/stats.h gains the six names with every flag spelled: `ext_pnorm(x, mu, sd, lowerTail, logP)`,
  `ext_qnorm`, `ext_dnorm(x, mu, sd, logP)`, `ext_dt(x, df, logP)`, `ext_qt(p, df, lowerTail, logP)`, `ext_qchisq`. Under R each
  is a macro for the Rmath name the code calls today (`Rf_pnorm5` and so on), so the object code is the same; standalone each is a
  macro for a renamed copy, `dbarts_nmath_pnorm5` and so on. The existing convenience names (`ext_quantileOfNormal`, ...) are
  defined over them. A slim src/include/external/nmath.h declares the six without pulling Rmath.h, for random.c, which today
  declares two by hand. model.hpp's fourteen direct calls become the ext_ names.
- The copy: src/external/nmath/, the closure of the six in R's nmath at a pinned R release (4.5.x, the version named in a
  README there): about 25 files, toms708.c (pbeta, for qt) and qgamma.c and pgamma.c (for qchisq) the large ones. Unedited
  except a prefix header force-included to rename every external symbol, so a host that also links libRmath sees no clash;
  compiled with `MATHLIB_STANDALONE` by the standalone recipe only. inst/COPYRIGHTS names R Core's copyright and the GPL (>= 2),
  the package's own license.
- Printing: src/include/external/io.h maps `ext_printf` and `ext_fflush_stdout` to `Rprintf` and `R_FlushConsole` under R as
  now; standalone, to stdout. src/external/io_standalone.c gives `ext_printMessage` (stdout), `ext_issueWarning` (stderr,
  "Warning: " first) and `ext_throwError` (stderr, then `abort()`; its two callers are unreachable invariants, an unknown
  generator kind). The misc pointers and their four installers go (dec-A74, dec-C04).
- R's generators: the parts of randomBase.c that use `unif_rand`, `norm_rand`, `.Random.seed` and `set.seed` compile under R
  only; the Mersenne Twister and the rest compile everywhere.
- Interrupts and the callback leaf: a host unit per build, src/host_r.cpp (the R_ToplevelExec poll, the R_UnwindProtect leaf,
  the slow-count warning, all moved from the bridge) and src/host_standalone.cpp (no poll, a plain call to the callback with
  C++ exceptions caught). A binding supplies its own unit in its own build, as dec-A74 has it for printing.
- The standalone recipe: a tests/cpp target that compiles the engine, src/C_interface.cpp, src/sampler_build.cpp,
  src/host_standalone.cpp, external (Mersenne Twister, io_standalone, nmath) and misc with `-DDBARTS_HOST_STANDALONE` and no R
  include path, and links no libR. It is the recipe a binding copies; nothing else is shipped.

### The flat entry file without R

src/C_interface.cpp and the shared helpers it calls move off R: `R_alloc` scratch becomes owned vectors, `R_finite` becomes
`std::isfinite`, the shared validators (`validateResponseSupport`, `enforceBinaryWeightPolicy`, `refuseNonFinite`,
`refuseMultiForestResponseMutation` and the rest of the common header's raisers) throw `Refusal` instead of calling `Rf_error`,
and the R entries that share them catch and raise the same text. R_interface_bartcore_common.hpp splits into an R-free core
header (holder, `Refusal`, validators, `buildSampler`) and the SEXP-facing declarations that stay R-only. `predict` stops being
main-thread-only for R's sake (no `R_alloc`); the main-thread rule stays for `run`, whose poll is the host's.

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
```

Types: `DBARTS_STATUS_LIST` and `dbarts_status` added, `DBARTS_SWAP_STATUS_LIST` and `dbarts_swap_status` removed;
`dbarts_sampler_spec` and `DBARTS_SAMPLER_SPEC_INIT` added; in `dbarts_predictor_source` the field `denseChannels` (per column,
k >= 0 double column k of `denseValues`, k < 0 code column ~k of `denseCodes`, NULL every dense column in `denseValues` at
`columnSources[j]`) is added above the boundary beside `denseCodes` and `numDenseCodeColumns`, which move above it too, and
the kind-selected reading of `columnSources` for a dense column goes (dec-B124): the struct then is the engine's
[`PredictorSource`](../../src/bartcore/data.hpp) field for field, `columnTypes` and `categoryCounts` read at creation and checked
and ignored by predict. `DBARTS_PREDICTOR_SOURCE_INIT` and `dbarts_dense_predictor_source` follow. `dbarts_results`,
`dbarts_draw`, `dbarts_column_type`, `dbarts_leaf_model` and `dbarts_family` keep their layout and values (the family and leaf
enums gain a reader in the spec).

Stub view: `DBARTS_R_CHECK(call)`, a macro that evaluates the call and, on a negative status, hands it to
`dbarts_R_raiseStatus`, an R-only callable the R build registers (under open call 6's recommendation): it raises R's interrupt
condition for INTERRUPTED, as the flat run does today, and otherwise `Rf_error("%s", dbarts_lastErrorMessage())`. It sits in a
second list, `DBARTS_C_API_R_LIST`, folded into the hash and stubbed only under `DBARTS_USE_STUBS`.

Prose: the contract block's error paragraph is rewritten from "An Rf_error means THIS call is wrong ... longjmps" to the status
rule and the table above; dec-A32's two meanings become "a value entry cannot fail; every other entry returns a status"; THE
HANDLE gains the C-created origin and its lifetime; the RAISING paragraph of `dbarts_draw_callback` says a thrown exception
becomes CALLBACK_ERROR and an R raise is still R's jump; `dbarts_sampler_run`'s interrupt, callback-stop and slow-count
paragraphs follow the statuses; `dbarts_sampler_predict` drops "main-R-thread-only ... R_alloc"; every per-entry "raises" becomes
the status it returns. About +500/-150 lines in the header.

The hash: [`dbarts_foldLayout`](../../src/C_interface.cpp) gains the spec's size and every field's name and offset in pointer
units, as `dbarts_results`, and `DBARTS_STATUS_LIST` joins the enum fold; both tokens are re-baked
([`dbarts_apiSignatureToken`](../../src/C_interface.cpp) and the literal). Several by-value doubles in the spec make its offsets
differ by platform in pointer units, as `dbarts_draw`'s do, so the spec folds as `dbarts_draw` does: field count, and doubles by
name, position and width. tools/check-api-hash.sh is unchanged and still skips (no release tag).

## Refusals

R's words do not change: every refusal an R user can reach is raised by the R parse before the spec is filled, or by the shared
stage with the text it has today, and the tinytest pins hold them. New words reach a C caller of `dbarts_sampler_create` only,
read back through `dbarts_lastErrorMessage`, each prefixed `dbarts_sampler_create: ` (R7):

- `spec is NULL`; `spec.structSize is 0 - set it to sizeof(dbarts_sampler_spec) (e.g. dbarts_sampler_spec s =
  DBARTS_SAMPLER_SPEC_INIT)`;
- `the response is NULL`; `the response must be finite`; then the support rule's existing words, family by family;
- `family %d is not one this entry builds; state a dbarts_family from GAUSSIAN to STUDENT` (AUTO and out-of-range values);
  `a multinomial sampler is built through the R interface` (UNSUPPORTED);
- `a gaussian, Student-t or aft sampler needs a finite positive sigmaEstimate` (the backstop createHolder words in R's terms);
  `sigmaPriorDf and sigmaPriorScale must be finite and positive`;
- `numChains must be at least 1`; `printEvery must be at least 1`; `no seed: set spec.seed or spec.chainSeeds`;
- the source refusals of `translateSource` with this entry's name, the test source's under `xTest`;
- `block tree counts must sum to numTrees`; and the engine's own creation refusals (leaf covariates, variance forest, monotone
  with a proposal mix, storage single), worded as the engine words them, with the entry's prefix.

## RNG class

Neutral on the R build, by construction and by gate:

- S1 changes which macro spells each math call, not which function runs: under R every `ext_` name is the Rmath symbol the code
  calls today. S2 and S3 move code without changing what is computed or in which order the generators are touched; the one place
  order could move, the chain seeds, keeps both of today's paths (one seed through a seed generator, or seeds drawn from R's
  stream) and passes them through the spec unchanged.
- Proven by: the equivalence trio bitwise with `--strict-coverage` against the baselines the MANIFEST names as current at the
  branch's base, counting the "identical draws (same RNG stream)" lines ([Gate hygiene](README.md#gate-hygiene)), on the
  reference build and again on the shipped build; the four seeded-drift snapshot files on the reference build; the full tinytest
  suite in one process, every message pin unchanged; tests/cpp; and the x86 leg, recording the base tip on the box and comparing
  --bitwise there (the arm64 baselines are cross-host only on x86). A new arm proves the two creation routes agree: under R, a
  sampler created by `dbarts_sampler_create` from the spec the R parse would fill, and the R-created one, give the same draws bit
  for bit (consumer.c).
- The standalone build is new. Its math is R's code compiled by another compiler with other flags (clang contracts a * b + c by
  default on arm64, R's binary may not), so it owes no bitwise match to R. Its gates: tests/cpp compares each copied routine with
  Rmath on a grid of a few thousand points, at most 2 ulp apart, and reports how many agree exactly; the standalone consumer's
  draws are compared with the R-hosted consumer's at the same spec and seed, bitwise where every grid point agreed exactly and
  otherwise by the equivalence harness's statistical z.

## Calls made

1. The vendored set is the closure of all six routines the engine and the support library call, not the four TODO names: qt
   ([`forestMapEntry`](../../src/bartcore/combiner.hpp)) and qchisq at 0.5
   ([`ChiKHyperprior::medianOfChi`](../../src/bartcore/model.hpp)) arrived after the ruling and are read at creation and at live
   restatements (setLeafPrior, and coefficient-prior-live's `$setForest(amplitude = )`). Alternative: move them host-side as the
   residual prior's quantile is, which cuts the copy by about 4,000 lines but adds a host-resolved number to every live
   restatement entry and makes a pure function of the prior the host's. dec-B85's words, "a copy of R's routines kept in the
   support library", fit the copy.
2. The switch is one preprocessor macro, `DBARTS_HOST_STANDALONE`, that only a standalone build defines. Alternative: a configure
   option, which the R build does not need and Windows' checked-in Makevars would have to mirror.
3. The copy keeps R's file layout unedited and renames symbols by a force-included prefix header. Alternative: edit the names in
   place, which makes the next refresh from R a diff to hand-merge.
4. The shipped spec is one struct and the R-only parts an internal extension beside it, both read by one build stage. Alternative:
   one internal spec of which the shipped struct is a projection; two homes for the same fields, the drift dbartsSpec() was made
   to end.
5. Every creation route, the data handle's included, fills the spec. Alternative: leave the handle and multinomial routes on their
   own parse, keeping two translate stages against dec-B124's one.
6. Inside the core every refusal is a thrown `Refusal`, caught at the entry. Alternative: thread a status return through every
   helper, which rewrites every call path and replaces the unwind dec-B119 keeps.
7. Exception types map to statuses as in [Statuses](#statuses). Alternative: one generic failure code, which leaves a host unable
   to tell an allocation failure from a wrong argument.
8. A callback's R raise on an inline run is still R's jump, resumed at the entry. Alternative: swallow it into CALLBACK_ERROR,
   which drops the target of an enclosing tryCatch and leaves R's condition stack inconsistent.
9. The slow-count warning goes through the build's warning routing, raised last as today. Alternative: a status note, which no
   ruling asks for and dec-B167 left open.
10. The return type is the enum `dbarts_status`, so a status entry's type says it is one. Alternative: `int`, as the swap status
    is today; the enum costs nothing on R's platforms, none of which builds with short enums.
11. `dbarts_sampler_getSigmas` is named after `$getSigmas()` and refuses where R returns NULL. Alternative: write NaN on a
    heteroscedastic sampler; a status is how this header says "no such channel" (getLatents).
12. A C-created handle is freed whole by `dbarts_sampler_destroy` (dec-B112). Alternative: a second entry, `dbarts_sampler_free`,
    which keeps "a second destroy is a no-op" for both origins at the cost of an entry whose only use is to differ.
13. The verbose creation summary stays in the R bridge (it prints the residual prior's quantile, which the spec does not carry); a
    C-created sampler prints only run progress. Alternative: move it into the build stage reading the spec, which changes R's
    printed summary unless the quantile rides the spec as well.
14. The specification states a family; AUTO is refused at C creation. Alternative: dispatch on the response shape as R's empty
    family string does, which consumer-spec-surface.md records as silently wrong for aft and logistic.
15. Student-t is `family = DBARTS_FAMILY_STUDENT` with `residualDf`, mirroring the enumeration's existing meaning; `residualDf`
    is read under STUDENT only. Alternative: GAUSSIAN with a finite `residualDf`, as the engine's options hold it, which gives the
    STUDENT enumerator no use at creation.
16. The standalone build is a tests/cpp recipe, not a shipped library or a CMake project. Alternative: a top-level build for a
    libdbarts, a maintained surface with no consumer yet.
17. The two static_asserts that no shape or source field owns storage stay: the R parse still raises around both structs.
18. The arc is one branch merged once, the header change landing whole in S3 so stan4bart and treatSens rebuild once; S2 places
    the spec struct in dbarts.h on the branch already, with no entry, since nothing outside the branch sees it.

## Open calls

Each is something stan4bart, treatSens and any future binding code against, and none is settled by a ruling.

1. **Do the entries that cannot fail keep returning their value?** Background: dec-B119 has every entry point return a status
   instead of raising. Most accessors - the counts, kIsSampled, usesDart, family - cannot fail at all (a bad handle crashes, by
   the header's rule, rather than being reported). Options: (a) every entry that can fail returns a status and every entry that
   cannot returns its value, numTrees moving to a status since its forest index can be out of range; cost: a reader learns which
   is which from the return type, and stan4bart's twelve accessor calls do not change. (b) every entry returns a status and
   values come back through pointers; cost: one rule with no exceptions, every accessor call in both consumers rewritten, and
   `dbarts_sampler_numChains(s)` becomes two lines. Recommendation: (a).
2. **Where is a failure's message read?** Options: (a) one call, `dbarts_lastErrorMessage()`, giving the message of the last
   failing call on this thread, as `dlerror` and `errno` work; cost: a caller reads it before the next failing call on that
   thread. (b) an error struct passed to every entry; cost: every signature grows a parameter. (c) a message kept on each handle;
   cost: creation, which has no handle yet, needs a second channel. Recommendation: (a).
3. **What can a C program create in 1.0?** Background: R creates every model; dec-B86 holds the header to what a consumer or a
   named design calls, and dec-B160 left multinomial creation from C for after the release. Options: (a) every model of one
   forest: all six families and Student-t, every leaf model, monotone, DART, the constraints and the variance forest; multi-forest
   and multinomial models stay R's, and a later release can append them; cost: a Python user fits BCF through R or waits. (b) also
   the multi-forest and multinomial models; cost: about 40 more fields frozen at 1.0, with no consumer asking. (c) gaussian and
   the binary families only; cost: a smaller struct that grows field by field as families are asked for. Recommendation: (a).
4. **Is the description one struct or built up by calls?** Options: (a) one struct the caller fills, with a size field and an
   initializer holding dbarts()'s defaults, as the results and predictor structs already work; cost: its field order is fixed at
   1.0, appends only after. (b) an opaque description object with a creation call and one setter per setting; cost: about 40
   entries frozen at 1.0 and checked one call at a time, against dec-B86's small header. Recommendation: (a).
5. **Does creating from C require a seed?** Background: R draws chain seeds from its own stream when none is given; a C program
   has no such stream. Options: (a) required: a spec with no seed is refused; cost: one more field a caller sets. (b) a seed from
   the clock and process id when none is given; cost: a run that cannot be repeated unless the caller thought to seed it, and a
   host-specific source in a host-neutral entry. Recommendation: (a).
6. **How does an R package turn a status back into R's error or interrupt?** Background: today a failing flat call raises the R
   error itself, and an interrupted run raises R's interrupt condition. Under statuses stan4bart and treatSens must do it, and R's
   interrupt has already been consumed by the poll, so only dbarts can raise it again. Options: (a) the R-only part of the header
   gives `DBARTS_R_CHECK(call)`, backed by one R-only entry dbarts registers, raising the message as an R error and an
   interrupt as an interrupt; cost: one more symbol in the R-hosted view, and the two consumers wrap each fallible call. (b)
   each consumer writes `if (status < 0) Rf_error("%s", dbarts_lastErrorMessage())`; cost: Ctrl-C during a stan4bart fit then
   reaches the user as an error, which `tryCatch(interrupt = )` and the quiet interrupt behavior miss. Recommendation: (a).

## Slices

### S1: build-time math and printing (engine, support libraries, tests/cpp)

Steps: the six ext_ names in stats.h and the slim nmath.h; model.hpp's fourteen calls and random.c's two declarations onto them;
the copied nmath closure with its prefix header, README (R version, file list) and COPYRIGHTS entry; io.h's two mappings and
io_standalone.c; randomBase.c's R-only parts guarded; misc/io.c, misc/io.h, their installer in R_interface.cpp, tests/cpp's
main.cpp and the three benchmark kernels' installers removed; a tests/cpp `nmath` component suite comparing each copy with Rmath
on its grid (run by the default libR-linked target) and a `standalone-math` target that compiles external and misc with the
switch and links no libR.

Gates: tests/cpp (both targets); full tinytest; the trio bitwise; the four snapshots on the reference build; R CMD check
--as-cran from a clean tarball (the copied files must not be compiled by the package build and must not draw a NOTE); sanitizers
on tests/cpp; mutation-battery verify-anchors; the x86 leg. Diagnostic, not a gate: `objdump -d` of the installed dbarts.so
before and after, expecting only R_init_dbarts and the removed misc symbols to differ.

Budget: ~350 lines written, ~6,000 copied. Agent ~4.5, review ~1.5, gates ~1.5.

### S2: the plain-C specification (bridge)

Steps: `dbarts_sampler_spec` and its INIT in dbarts.h (no entry yet); src/spec_internal.h with `SpecExtension`;
src/sampler_build.cpp with `buildSampler`, taking over the engine-building half of createHolder, optionsFromParsed and the two
multinomial and two handle routes; `translateSource` moved into the core and used by the R parse's views, with `denseChannels`
added to the shipped source and the R parse's code channel passed through it (dec-B124); the R parse filling spec and extension
(the parse functions keep reading SEXPs; their outputs become spec fields); `Refusal` and the R entries catching it and raising
the same text; seeds through `seed` or `chainSeeds`. A tests/cpp suite fills every spec field with a distinct value and checks
the engine options built from it field by field, so a field the build forgets fails by name.

Gates: as S1, plus the R-loaded ASan run over the creation suites; every message pin in tinytest unchanged (run the whole suite
in one process).

Budget: ~1,900 lines (header ~250, extension ~150, build stage ~650 moved, parse rework ~450, tests/cpp ~300, the rest
`translateSource`). Agent ~10, review ~4, gates ~3.

### S3: statuses, the creation entry, the sigma getter, the header (bridge, flat entry file, tests)

Steps: `DBARTS_STATUS_LIST`; the list as in [The header diff](#the-header-diff); `dbarts_sampler_create`, `getSigmas`,
`lastErrorMessage`; C_interface.cpp off R (owned scratch, statuses, the thread's message); host_r.cpp and host_standalone.cpp;
the holder's origin and destroy; `dbarts_R_raiseStatus` and `DBARTS_R_CHECK` (open call 6); the hash folder and both tokens; the
header's prose; R_interface.cpp's registration (the list does it) and the R-only callable. Tests: consumer.c and test-capi.R
rewritten for statuses (each entry's failure reached and its status and message read; INIT plus data and seed against
`dbarts()`'s draws bitwise; a C-created sampler against the R-created one bitwise; `denseChannels` with packed blocks against
the dense source; getSigmas on two chains and its UNSUPPORTED on a heteroscedastic sampler; `DBARTS_R_CHECK` raising an error
and, through the interrupt test hook, an interrupt; a C-created handle destroyed), the five handshake arms kept; tests/cpp's
`standalone-consumer` target: a C99 file including the prototype view, creating gaussian, probit and Student-t samplers, running,
predicting, reading sigma, refusing, destroying, built and linked with no R, and `nm` showing no undefined R symbol; the
cpp-tests workflow runs it. Docs: dbarts-embedding.Rd's C section, the architecture's layering and callbacks sections, NEWS.

Gates: as S2, plus test-capi.R under ASan on the R-loaded path, tests/cpp's standalone consumer under ASan and, on the x86 leg,
LeakSanitizer (C-created handles are new owners); the header-compiles target; `sh tools/check-api-hash.sh`.

Budget: ~2,300 lines (header +500/-150, C_interface.cpp ~500, host units ~250, bridge ~100, consumer.c ~450, test-capi.R ~300,
tests/cpp ~350, workflow ~15, docs ~120). Agent ~10, review ~4, gates ~3.

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
ingestion, resolution, its own host unit); docs/plans/python-bindings.md's Context; the feature matrix's reach column for the
flat header; INDEX rows; NEWS's 1.0-0 entry for compiled consumers. Doc-freshness only.

## Tests and gates

Per slice above; across the arc, run independently of the implementer on each slice's own library, `--preclean` on every
install (the facade's virtuals do not change, but every creation route does):

```
R CMD INSTALL --preclean -l <lib> .                       # and --configure-args=--enable-reference-build for the snapshots
cd tests/cpp && make && ./test_bartcore                   # header-compiles, the libR suite
make standalone-math standalone-consumer && ./standalone_consumer
nm standalone_consumer | grep -E ' U _?(Rf_|R_|Rprintf|unif_rand)'   # expect nothing
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts", lib.loc = "<lib>")'
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare --bitwise --strict-coverage benchmarks/baselines/<current>.rds
R_LIBS=<lib> Rscript benchmarks/R/bcf-equivalence.R compare --bitwise benchmarks/baselines/<current>.rds
R_LIBS=<lib> Rscript benchmarks/R/multinomial-equivalence.R compare --bitwise benchmarks/baselines/<current>.rds
  # every scenario "identical draws (same RNG stream)", none skipped, no "max |z|"
# the four test-reproducibility-*.R files on the reference build, run_test_file each
Rscript -e 'lintr::lint_package()' && air format --check . && Rscript tools/check-rc-codoc.R . \
  && Rscript tools/check-win-drift.R . && Rscript tools/check-doc-freshness.R . \
  && Rscript benchmarks/R/mutation-battery.R verify-anchors
R CMD check --as-cran <tarball from a clean copy>
sh tools/check-api-hash.sh                                 # "no release tag, skipped"
```

plus sanitizers on tests/cpp and on the R-loaded path ([Gate hygiene](README.md#gate-hygiene)), and the x86 leg before landing
each slice (tests/cpp plain and with ASan, UBSan and LeakSanitizer, tinytest, the trio recorded on the base tip there and
compared --bitwise).

## Mutants

The reviewer of each slice runs these; each must fail the named gate, and is reverted with a `touch`:

1. S1: one model.hpp pnorm site's lower-tail flag flipped - the probit and aft equivalence scenarios.
2. S1: a constant in the copied qnorm changed - the nmath suite.
3. S1: `ext_qnorm` mapped to the copy in the R build - nothing should fail bitwise only if the copy is exact; record which
   scenarios move (it measures how near the copy is, and proves the R build reads Rmath).
4. S2: one spec field (perturbProbability) not copied into the options - the field coverage suite, and the perturb equivalence
   scenario.
5. S2: `chainSeeds` ignored in favor of `seed` - every unseeded tinytest snapshot and the C-versus-R identity arm.
6. S2: a dense column read in the block its kind selects instead of through `denseChannels` - consumer.c's packed-blocks arm.
7. S2: one R-worded refusal moved behind the spec fill, reworded by the shared stage - its tinytest pin.
8. S3: setSigma's refusal returning OK - consumer.c's status arm.
9. S3: the thread's message not written on one entry - consumer.c's message arm.
10. S3: `dbarts_R_raiseStatus` raising an error for INTERRUPTED - test-capi.R's interrupt arm.
11. S3: getSigmas writing chain 0 to every slot - its two-chain arm.
12. S3: two spec fields swapped in the header - the build's static_assert against the re-baked hash.
13. S3: destroy on a C-created handle releasing the engine only - LeakSanitizer on the x86 leg.

## Interaction with the other pending plans, and order

- forest-door (slices 1 to 3) and coefficient-prior-live land first, as their plans say. They rework the forest records, the
  calibration map and the live forest writes in the same bridge and engine files; S2's extension carries the forest records as
  they leave them (the share among them), and S1's qt and chi-median copies serve their live restatements. S1 could be built
  beside them and rebased, being nearly file-disjoint (two qt sites in combiner.hpp); S2 starts after both land, from a main
  whose engine CI is green.
- frequency-weights (docs/design/frequency-weights.md, PROPOSED with open calls) adds "one enumeration field on its family part
  in dec-B85's change" and a getLatents comment. Neither consumer creates through C, so the field moves the hash and rebuilds no
  one: if frequency weights' engine and R slices land first, S2's spec carries the field; otherwise its slice 4 appends it,
  before 1.0-0, with a re-bake. Recommended order: this arc first, since it waits on nothing but the two plans above and the
  note still has calls open; the getLatents comment rides whichever header edit comes second.
- multinomial-doors' C creation path stays post-release (dec-B160), consistent with open call 3 (a).
- The sister re-verification at the freeze follows this arc (TODO), with stan4bart's sanitizer job dispatched against the final
  header (dec-A36).
- Engine features merge serially: nothing else touching src/ merges while S2 or S3 is in review.

## Cost

About 5,000 lines written here and 6,000 copied, ~130 in the two consumers; ~50 agent-hours estimated, 75 to 100 by the record
of plans here. The copy adds about 250 KB to the source tarball, compiled by nothing in the R build.

## Not verified

- That the copied routines agree with the Rmath a given R binary carries bit for bit; the flags R's own build used are not known
  here, so the standalone gate is stated within 2 ulp and the bitwise count only reported.
- The exact file closure of qt and qchisq in R's nmath; about 25 files is a reading of the call graph, not a build.
- That `objdump` of dbarts.so before and after S1 differs only where named; relative call offsets may make the diff noisy.
- How forest-door slice 3 and coefficient-prior-live leave the forest records, which S2's extension must mirror; their
  branches were read for the files they touch, not field by field.
- stan4bart's and treatSens' call sites were read on their main and master; their 1.0 compat branches may differ.
- That no R-hosted consumer other than the two links dbarts.h; the reverse LinkingTo list on CRAN was not queried.
