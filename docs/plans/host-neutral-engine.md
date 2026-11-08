# host-neutral-engine: a plain-C specification the engine is built from, build-time math and printing, a status from every entry point, and a creation entry in dbarts.h

Status: PLANNED 2026-11-03 on main (c3d268c8f), for dec-B85; amended the same day after a blind critique (BUILD AFTER
AMENDMENTS, then BUILD), for the maintainer's rulings on every open call, and rebased onto main at 8fce28568 after forest-door
landed. S1 LANDED on main 2026-11-03 as its own feature (dec-A260). Rebased onto main at 7937d252f and amended 2026-11-05 for
the work landed or planned since ([Amendments of 2026-11-05](#amendments-of-2026-11-05)): S1 as built, coefficient-prior-live
(dec-A261), hurdle-completion (dec-A262), latents-first (dec-B586, dec-B589), frequency-weights as amended and
approximate-draws (dec-B583 to dec-B585). Rebased onto main at f79874ea8 and amended 2026-11-07
([Amendments of 2026-11-07](#amendments-of-2026-11-07)): dec-B615 (the sampler's facts in one struct), the landings since
7937d252f (coefficient-prior-live, hurdle-completion, latents-first, and the fourth review's slices S-D, S-F and S-A, with
docs-10's R share), the fourth review's S-C, S-B and S-G and frequency-weights' S1 to S3 on their branches, and one new open
call. Rebased onto main at 7e7edb526 and amended again 2026-11-07
([Amendments of 2026-11-07, second](#amendments-of-2026-11-07-second)): dec-B639 rules open call 8, the per-forest facts
reported through an array of per-forest structs at a stated element size; the fourth review's S-C (dec-A269), S-G (dec-A267) and S-B
(dec-A270) landed, and the work in flight placed (rchk-findings, frequency-weights' S7 and S8, state-basis-levels for
dec-B647, small-rulings-1107, approximate-draws, the dec-B640 slice). No call is open. The main feature is not built.
Rulings, quoted where they settle a call: dec-B85
(the staged route), dec-B84 (the header is pure C), dec-B119 (a status from every entry point, the unwind to the entry point staying), dec-A74 (printing chosen at build time, the
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
binomial response for logistic), dec-B615 (the header describes a sampler through one struct led by its size, filled by one
call, the one-fact entries removed; "Yes, let's do that."), dec-B639 (each forest reported through the caller's array of
per-forest fact structs at a stated element size; "Let's just do (d) now, and think ahead."), with dec-B601 (sigma NaN on a sampler with a variance forest, and
a fact saying whether it has one) and dec-B600 (a registered callback receives every channel the fit carries), dec-B593
(approximate-draws' setting is a method name and a separate tolerance). The runs behind it: reads of the tree at
c3d268c8f, of stan4bart main (27ef163) and treatSens master (268b2db), and of frequency-weights (da210ce2d); the critic's probe of
R 4.6.2's nmath against Rmath (17,304 points, below). No prototype was built here. The 2026-11-05 amendment read main 7937d252f,
coefficient-prior-live at 0c485ed44, hurdle-completion at d27be2cc7, and the plans latents-first (f3a818b2b), frequency-weights
(72914ff71) and approximate-draws (e0c62e566) on their branches, with S1's builder and review reports; it built nothing. The
2026-11-07 amendment read main f79874ea8 against 7937d252f, the branches review4-prior-predictive (59c3e59a5),
review4-sigest-bases (c170d0f9e), frequency-weights (831e0f434) and approximate-draws (94aed03c1), slice S-B's build
prompt with its plans and accepted critique (its branch, review4-draw-channels, not yet made), stan4bart main (8e23a1b) and
treatSens master (268b2db); it built nothing. The second 2026-11-07 amendment read main 7e7edb526 against f79874ea8 (the
records of S-C, S-G and S-B with their Landing notes and dec-A267 to dec-A270), the branches rchk-findings (843470f35),
frequency-weights (ab739157c), state-basis-levels (5b068cceb, its plan), small-rulings-1107 (2e170e9f1) and
approximate-draws (e50e1145a), dec-B628, dec-B632, dec-B639, dec-B640 and dec-B647, stan4bart main (98db852) and treatSens
master (268b2db); it recounted by grep and built nothing.

agent: S1 landed (one opus implementer, one opus review with its mutants and a fix round). The main feature: S2a, S2b and S3
one opus implementer each (effort high), S2b the riskiest; S4 one sonnet implementer per consumer, in that consumer's own
checkout; S5 sonnet. One opus review per engine slice (S2a, S2b, S3) that runs the mutants below, the same reviewer checking
the fix round unless the fixes change engine behavior. Tier "Changes draws" for every engine slice ([Process by risk](README.md#process-by-risk)): nothing is
meant to move a draw, but every creation route and every flat entry is rewritten, and the gates are the draw-changing tier's.
rng: neutral, bit for bit, on the R build: the equivalence trio, the four snapshot files, tests/cpp and the full tinytest suite
do not move against the baselines main names at the branch's base (S-G's re-recorded equivalence and bcf-equivalence
baselines on main now, and any frequency weights, dec-B640 or approximate-draws records after them), and no R message changes. The standalone build is new; its math is bitwise equal to Rmath on the reference host and within a relative
floor elsewhere; see [RNG class](#rng-class).
window: before the 1.0-0 submission (dec-B85: "post merging into main, but before 1.0 since it might change a number of packages
and is an interface change"). Two features. S1 (build-time math and printing) LANDED on main 2026-11-03. The main feature
(branch host-neutral-engine: S2a, S2b, S3, then S5) starts after frequency-weights, the dec-B640 slice and
approximate-draws merge (the fourth review's S-C, S-G and S-B landed on 2026-11-07), the release map's order
([Interaction with the other pending plans, and order](#interaction-with-the-other-pending-plans-and-order));
it merges once (dec-B494), the consumer ports (S4) merging right after it.
budget: S1 as built ~910 lines written, 6,329 copied unedited from R's math library and 395 of reduced R headers. The rest
~6,120 lines written: S2a ~1,000, S2b ~1,790, S3 ~3,330; S4 ~115 in stan4bart and ~50 in treatSens. Agent time for what
is left, estimated and not measured: S2a ~22.5, S2b ~33, S3 ~35, S4 ~7, S5 ~1, rebases onto main ~5: ~103.5. Plans here run
1.5 to 2 times low, so 155 to 205 (S1's builder reported about 3 hours against the ~5 planned, the simplest slice). Counted
with frequency weights, the dec-B640 slice and approximate-draws landed first (52 fields, `trials` among them). The
sampler-facts struct (dec-B615, with dec-B639's per-forest element) is about 420 of S3's lines and 5 of its hours, and 35
lines and an hour of S4's. The callback
statuses and the status words (dec-B579, dec-B581) are about 300 of S3's lines and 3 of its hours. Counted under dec-B575's
middle course, with the forward-compatibility work it asks for
([Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial)):
about 300 lines and 4 hours over the bare middle course, which itself took 400 lines and 10 hours off every one-forest model.

## Summary for the maintainer

You ruled in September that before 1.0-0 the engine stops depending on R below the layer that talks to R, so that a program
other than R, a Python package for instance, can build and run a sampler through the one C header dbarts ships. This plan builds
that as two pieces of work.

The first, the math and the printing, has landed. The engine and its support code call six of R's numerical routines (the
normal distribution's density, probability and quantile; the t distribution's density and quantile; the chi-squared quantile)
and R's console. Under R they still call exactly those, and every draw an R user makes was checked to be bit for bit the same.
Any other build uses a copy of R's own routines, kept in the repository but left out of the package R users install, and
writes to standard output. The copy gives exactly R's answers on the two machines it was checked on. The two unused printing
hooks are gone. Frequency weights adds a seventh routine the same way, R's log of a binomial coefficient, one more file of
R's copied; that is built on its branch.

The second is the rest. Today the code that turns R's control, model and data objects into an engine reads those objects
directly. It will instead fill one plain C description of the sampler - its data, family, priors, constraints, how its weights
are read, a binomial response's trial counts, and run settings, each already resolved - and the engine is built from that
description alone. R still does everything R does now: reading the formula, choosing the family from the response, estimating
the starting sigma, turning the residual prior's quantile into a scale, parsing the prior vocabulary, and naming columns and
forests in its messages. Nothing an R user sees changes. Then the header changes once, so the two compiled packages rebuild
once: a program can create a sampler from that description, every entry that can fail returns a status in place of raising an
R error, the failure's message can be read back, each chain's current sigma can be read, the predictor description says column
by column which block a column's values sit in, and a sampler's facts - its sizes, its family, whether k is drawn, whether it
uses DART, a variance forest or a Student-t residual - are read in one call into one struct, in place of the eleven separate
calls that report them today, as you ruled on 2026-11-06. stan4bart and treatSens are then changed to check those statuses,
and stan4bart reads its sizes from the struct; with the helper you chose, an R user of either sees the same errors,
interrupts and time limits as today.

You have ruled every choice this plan raised. Entries that cannot fail keep returning their value. A failure's message is read
from one call, per thread, and a second call gives fixed words for each status, the same text for every program, which the R
helper falls back on and never replaces with words of its own. A C program creates the middle of the model in 1.0, with the
description shaped so that models of several forests and the multinomial model can be added later without changing what any
1.0 field means; the plan says what that addition looks like and tests that a 1.0 description reads the same under it. The
description is one struct with a size field and an initializer. Creating from C requires a seed. A callback's return says what
happened: zero continues, a positive value stops the run on purpose, a negative value reports the callback's own failure, and a
C++ exception thrown from a callback is caught on whatever thread ran it and reported the same way. The facts struct follows
the description's pattern: a size field first, an initializer, fields added only at the end. On 2026-11-07 you chose how it
reports each forest: the caller hands it an array of small per-forest structs and says how large each is, and the call fills
each one up to the number the caller gave. In 1.0 each holds its forest's tree count. Later facts about a forest are added to
the end of that small struct, and it is shaped beside the per-forest part of the creation description sketched for later, so
what a forest reports matches what creating it states. stan4bart hands it an array of one.

What is left costs about 6,100 lines and roughly 155 to 205 agent-hours, most of it in rewriting how samplers are created. The
fourth review's fixes have all landed: the prior predictive draw, the callback receiving every channel the fit carries (sigma
reported as not a number where a variance forest gives each row its own scale), and the starting sigma that reads each
forest's basis. None changes what this plan builds; the callback's channels and the k flag that follows the model are carried
into the new header as they stand. Several pieces of work come first, in the release map's order: a fix to how the R layer
protects its objects from R's memory manager, which touches the code this plan rewrites and lands next; frequency weights;
the change that lets a multinomial sampler take new data, which you ruled on 2026-11-07; then the approximate draws. None
makes this plan write anything twice: the description and the header carry the final weights reading and the binomial trial
counts from the start, and the approximate draws' setting gets no field in the 1.0 description (a C program gets the default
an R user gets, and the field can be added later without changing the meaning of anything already there). Two smaller pieces
of R work may land any time: a stored state that matches a factor basis's levels by name, and a batch of small rulings.
Other work can still land between this plan's pieces; the plan rebases and rechecks after each.

## Goal

The engine, the support libraries and the flat C entry file compile and run with no R header and no R library: a tests/cpp
recipe builds them into a plain C consumer that creates, runs, predicts, reads sigma and the sampler's facts and refuses through
statuses alone, and no R symbol is in its binary. Under R every draw, every message and every printed line is what main produces at the base of each
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
  one it consumed, at five sites on main: the flat run, the R run, its burn-in route, the sweep-closure entry and S-C's
  [`bartcore_sampleHyperparametersFromPrior`](../../src/R_interface_bartcore.cpp), which polls once after a call that made
  joint draws; a slow monotone count warns through [`raiseSlowCountWarning`](../../src/C_interface.cpp). The handle is the
  bridge's holder, [`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp), which carries no R type.
- R frames and C++ exceptions: [`unwindProtect`](../../src/R_interface_bartcore.cpp) protects a body against R's jump and catches
  no C++ exception, so an exception must never leave a body it wraps or cross any R frame. The shared refusals declared in
  R_interface_bartcore_common.hpp, seventeen on main (`validateResponseSupport`, `enforceBinaryWeightPolicy`, `refuseNonFinite`,
  thirteen more validators, and the `[[noreturn]]` `refuseCountOverCap` the first two call), raise `Rf_error` today; their names
  stand on 74 lines of R_interface_bartcore.cpp, definitions included, and 12 of the flat file, at main 7937d252f and at
  f79874ea8 alike (the first count, "about 60" at c3d268c8f, missed some; comment lines are not counted). frequency-weights'
  S2 and S3, on its branch at 831e0f434, add three (`refuseWeightedCountsOverCap`, `refuseWeightedTrialsOverCap`,
  `refuseBinomialResponseSwap`), twenty in all on 90 bridge lines and 13 flat. The bridge's other `Rf_error` sites are R's own
  and stay: 421 lines of the file at f79874ea8, 379 at c3d268c8f, forest-door adding 32, coefficient-prior-live 5 and the
  fourth review's S-D and S-A 5; 428 on frequency-weights' tip and 423 on S-C's. At main 7e7edb526 (the second 2026-11-07
  amendment, a grep of the twenty names in the two files with comment lines dropped, a wider match than the counts above,
  which matched calls and definitions only): S-C, S-G and S-B moved no shared-refusal line (90 bridge and 20 flat lines at
  f79874ea8 and at 7e7edb526 alike, 108 and 21 on frequency-weights' tip ab739157c), and added three R-only `Rf_error` lines,
  S-C's two (a captured engine refusal and the attempt count) and S-B's setCutPoints refusal, about 424 by the earlier count;
  frequency-weights' tip adds eight over its base. S2a recounts on the base it starts from, by one rule.
- Process-wide state on the predict path: [`predictPartition`](../../src/bartcore/sampler.hpp), a plain struct holding a
  `std::vector`, is written by every predict and read by a test hook; it is safe today only because predict is main-thread-only.
  [`testFitPartition`](../../src/bartcore/chain.hpp) is relaxed atomics, last writer wins, a diagnostic no answer reads.
- Engine exception types: `std::invalid_argument` (thirteen sites at main, forest-door's one and S-C's two, the
  hyperparameter draw's refusals, included; S-B added none), `std::length_error` (two, saved-tree storage past the addressable size),
  `std::bad_alloc` (one), `std::logic_error` (one, S-A's: [`Sampler::run`](../../src/bartcore/sampler.hpp) started while a
  run of the same sampler is in progress, as from inside its sweep hook or a draw callback), and `std::runtime_error` at
  three: the two prior-draw exhaustions in [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) and predict's join, which rethrows any worker
  failure, an allocation failure included, as `runtime_error("predict worker w failed: ...")`. Since S-A an allocation
  failure inside a sweep, a prior tree draw or a rebuild leaves each chain as its state describes, its fits re-summed from its
  trees, a valid state to run on from (dec-A266; residuals in TODO allocation-failure-residuals).
- The R touchpoints below the bridge, as S1 left them on main (dec-A260):
  - math: the engine and the support library call six names spelled with every flag, `ext_pnorm`, `ext_qnorm`, `ext_dnorm`,
    `ext_dt`, `ext_qt` and `ext_qchisq` ([nmath.h](../../src/include/external/nmath.h)), R's own routines in the R package and a
    renamed copy of R 4.6.2's nmath under `DBARTS_HOST_STANDALONE` (standalone/nmath, build-ignored). The qt serves
    [`forestMapEntry`](../../src/bartcore/combiner.hpp) (a coefficient's t prior, at creation and at a live restatement,
    coefficient-prior-live's `Chain::setForestCoefficientPrior` among them, on its branch) and the qchisq
    [`ChiKHyperprior::medianOfChi`](../../src/bartcore/model.hpp); qnorm is on the hot path through
    [`ext_rng_simulateStandardNormal`](../../src/external/randomNorm.c), pnorm and qnorm in the truncated normal draws
    ([`ext_rng_simulateTruncatedNormalScale1`](../../src/external/random.c)).
  - printing: [`ext_printf`](../../src/include/external/io.h) and `ext_fflush_stdout` are R's console in the R package and
    standard output standalone, where standalone/io_standalone.c gives the message, warning and error routines; misc's print
    pointers and their five installers are gone.
  - R's own generators, `.Random.seed` and `set.seed` in src/external/randomBase.c compile in the R package only.
  - in the bridge's own core, still: `R_alloc`, `R_finite`, `R_UnwindProtect` (the callback leaf,
    [`ShippedDrawHook`](../../src/R_interface_bartcore_common.hpp)), `R_CheckUserInterrupt` (the poll).
  - The residual prior's chi-squared quantile is already the bridge's (`sigmaRawScale`), as dec-B85 has it.
- Build: src/Makevars.in sets no `OBJECTS`, so every .c and .cpp at the top of src/ is compiled into dbarts.so. configure (which
  reads `R CMD config`) generates src/config.hpp, src/misc/config.h, src/external/config.h and src/include/misc/types.h from
  their `.in` files; Windows uses checked-in `.win` variants, and a standalone build the checked-in standalone/config (S1).
  R's compiler on the reference host is `clang -arch arm64` with `-falign-functions=64 -Wall -g -O2`, so contraction of
  a * b + c is clang's default, on.
- tests/cpp (S1): the default target links libR and holds the nmath suite (the copy against libR on a 36,093-point grid);
  `standalone-math` builds misc (through misc's own Makefile), external and the copy with the switch into build-standalone and
  checks the binary references nothing of R; `r-math-symbols` checks dbarts.so calls R's six routines and carries no copy;
  `standalone-config` compares the checked-in configuration headers' macro values with configure's. No target yet compiles the
  engine, the bridge's core or the flat entry file without R: that is S3's `standalone-consumer`.
- The shipped header: [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h), [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)
  checked against [`dbarts_apiToken`](../../src/C_interface.cpp), [`dbarts_predictor_source`](../../inst/include/dbarts/dbarts.h)
  (today a dense column's index is read in the block its STORE kind selects, the rule dec-B124 replaces),
  [`DBARTS_SWAP_STATUS_LIST`](../../inst/include/dbarts/dbarts.h). The shipped test consumer
  [`capi_hash`](../../inst/tinytest/capi/consumer.c) and its handshake arms.
- Consumers, read-only: stan4bart main calls run (2 sites), sampleTreesFromPrior, setOffset (2), setSigma (2), getLatents (2),
  predict, setTreeStorage (2), printTrees, setVerbose, numTrees (forest 0), destroy, the version pair and six accessors
  (numObservations, numTestObservations, kIsSampled, numPredictors 3 times, numChains and numSavedSamples twice each; eleven
  one-fact calls with numTrees, forest 0's count read once to bound printTrees' indices; the same at 98db852 as at 8e23a1b
  and 27ef163), all in src/init.cpp, plus `dbarts_results` in src/bart_util.hpp and `dbarts_dense_predictor_source`. treatSens master calls run (3),
  setResponse (4), setOffset, setSigma, setVerbose (2), setNumThreads (2), destroy (2) and the version pair, in
  src/R_interface.cpp, src/bartTreatmentModel.cpp and src/sensitivityAnalysis.cpp. Both define `DBARTS_USE_STUBS`; neither
  defines `DBARTS_REQUIRE_EXACT_ABI`; neither creates through C or registers a draw callback. bartCause and bairrtt link nothing
  of dbarts.
- docs/plans/pure-c-header.md (the last header change, its slices and gates), docs/design/consumer-spec-surface.md (why one
  resolution home), docs/design/r-c-division.md, error-style.md (R7, a flat entry names itself; R8, the R message is canonical
  and a C backstop is worded apart), docs/design/frequency-weights.md on its branch, [Callbacks and errors](../architecture.md#callbacks-and-errors),
  [Reproducibility contract](../architecture.md#reproducibility-contract), docs/plans/python-bindings.md.
- The work landed since the fourth amendment: docs/plans/coefficient-prior-live.md, docs/plans/hurdle-completion.md and
  docs/plans/latents-first.md (Landing notes), the fourth review's records in docs/plans/review-2026-11-06/ (r-surface.md,
  r-surface-2.md, engine-robustness.md, docs.md), and since the fifth its prior-predictive.md (S-C, dec-A269),
  sigest-bases.md (S-G, dec-A267) and draw-channels.md (S-B, dec-A270), each with its Landing note. The work that lands
  around this plan, each on its branch until merged: TODO rchk-findings (branch rchk-findings);
  docs/plans/frequency-weights.md (Interaction with the other pending plans, and order); docs/plans/state-basis-levels.md
  (dec-B647); TODO multinomial-setdata (dec-B640, no plan yet); docs/plans/approximate-draws.md (What is built, Decision);
  TODO capi-sampler-info (dec-B615, dec-B639); S1's builder and review reports in scratch/hnm-s1 and scratch/hnm-rev.

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
  caller-to-library direction) holding the middle of the model dec-B575 names, every value resolved: 52 fields after
  `structSize` with frequency weights' `trials`, 51 without it (the ruling's "about 45" was an estimate). Its fields, grouped
  (the order below is the struct's):
  - data: `x` (a `dbarts_predictor_source`, its typing channel READ at creation, which gives `columnTypes` and `categoryCounts`
    the reader dec-B124 owes them), `y`, `trials` (frequency weights' binomial response, read under LOGISTIC only: a per-row
    trial count with `y` the successes, NULL a 0/1 response; a non-NULL `trials` under any other family is refused),
    `weights`, `weightsType` (a `dbarts_weights_type`, below, read only when `weights` is non-NULL), `offset`,
    `survivalStatus` (aft, double 0/1 per row as R holds it), `xTest` (a second source, 0 rows for none), `testOffset`, `numCuts` and
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
  constraints (columns, interactions, blocks), the move settings, approximate-draws' latent sampler setting (zero the default
  rule, as R's NULL; see Calls made, 35) and the variance forest; then several forests and multinomial
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
  a 0 a row out of the likelihood on every family, logistic included) with the words frequency weights gives R, prefixed:
  the rows of that plan's resolution table, as built, which a later engine change widens by the build taking what the engine
  then takes (approximate-draws' S3, real weights on logistic and the binomial under its route, among them), never by a new
  field. With frequency weights landed first, as recommended, S2b translates that table as it stands; were this plan first,
  FREQUENCY would be taken only where it reads as the family's own (the Polya-Gamma families' counts, 0/1 weights) until
  frequency weights widens it. The spec says nothing of predictive draws: a 0 row is still drawn, so the stream is what it
  was, and its draws are set to NA where dec-B580 puts it, on every family, its fitted value still reported.
- `bartcore_bridge::SpecExtension`, internal (src/spec_internal.h, plain C structs and arrays), in two parts. The PUBLIC-SHAPED
  part holds what later appends will ship, already in their sketched layout: the one-forest groups 1.0 leaves out (leaf models,
  monotone, constraints, moves, approximate-draws' latent sampler setting once that is built, variance forest), the
  several-forest records (`ForestRecord`, field for field the sketched
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
(`numFittedValuesPerObservation`, K or the forest count) and V (`numVariableCountForests`), both facts-struct fields from S3,
and multinomial runs and predicts through the flat API today (dec-B163), so the append adds nothing to `dbarts_results`. The
reported side is there from 1.0 too: the facts struct reports each forest through a caller's array of `dbarts_forest_info`
at a caller-stated stride (dec-B639; [The sampler's facts (dec-B615)](#the-samplers-facts-dec-b615)), forest f of it the
forest f of `forests` above, of `dbarts_draw` and of printTrees' index.

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
   at 1.0 and MULTINOMIAL later; `residualDf` under STUDENT; `shape` under NBINOM; `trials` under LOGISTIC; `y` and `offset`
   under every 1.0 family, not under MULTINOMIAL, which states the appended counts and category offset instead.
7. Nested arrays version by a stride, not a size per element: the parent states `forestStride`, the caller's
   `sizeof(dbarts_forest_spec)`, once; the build indexes by it and reads each element under invariants 1 to 3 at that size. The
   elements of one array share a type, so a size in each would repeat the stride; Win32's `SendInput(nInputs, pInputs, cbSize)`
   versions an array of structs the same way. A stride of 0 with `numForests` > 0 is INVALID_ARGUMENT.

What 1.0's validation refuses now, so the append is never a change of meaning: a nonzero byte past the known `structSize`
(invariant 3); MULTINOMIAL and any family value past 8 (5); a non-NULL `trials` off LOGISTIC (6), so the probit binomial that
dec-B567 leaves for later can be read under PROBIT without a field meaning something new; and nothing else, since every other
1.0 field already has its final reading. Placeholder fields are not shipped (see Calls made).

The later work is mostly exposure. From S2b the R path already builds several-forest and multinomial samplers through
`buildSampler` from the extension's `ForestRecord` array and `MultinomialRecord`, laid out as the sketch, and the equivalence
trio's BCF and multinomial baselines hold those routes bitwise. The append moves the two records' public prefixes into dbarts.h
as `dbarts_forest_spec` and the multinomial fields, adds `dbarts_forest_spec_init(element, size)` for bindings, and has
`dbarts_sampler_create` copy the caller's array into the extension at its stride before the same `buildSampler` call; the
R-only part (a sharer's share and a basis scale as first set, the anchor and unit) stays R's.

The two per-forest elements are designed side by side (dec-B639: "so a forest's reported facts mirror what creating it
states"). `dbarts_forest_info` ships at 1.0 with one field, `size_t numTrees`, the name, type, meaning and first position of
`dbarts_forest_spec`'s first field. Its rules are invariant 7's in the library-to-caller direction: the parent states the
array, its length and the caller's element size once (`forests`, `forestsLength`, `forestStride` in the facts struct, the
spec's `forests` and `forestStride` mirrored); the library writes element i at byte i times the stride, each field whose end
lies within the stride, and never a byte past it. A fact a later minor appends to the element takes `dbarts_forest_spec`'s
name, type and meaning where the spec states the same thing (a forest's `kind`, `numBasisColumns`, and the tree prior as it
now stands after a live restatement, `base` and `power`), and reads 0 where it is absent, so an older library's 0 is
"absent". The element reports what holds now; the spec element states what creation asks for; where they share a name the
values agree at creation.

The gate (S3): a tests/cpp target builds the library with `-DDBARTS_TEST_SPEC_APPEND`, under which dbarts.h appends the sketch
above to the spec, and the build stage honors it. A C99 translation unit compiled WITHOUT the define fills a 1.0-sized
description (structSize the 1.0 size) of a gaussian and a probit sampler; a second, compiled WITH it, fills the same values with
the append zeroed; both create through the appended build, and their draws equal the shipped build's at the same spec and seed,
bit for bit. The same target creates a two-forest sampler through the test append and compares it bitwise with the R-created
BCF fit of the bcf-equivalence scenario it copies, and hands the shipped build a description carrying a nonzero appended
field, which returns UNSUPPORTED. The two struct definitions never meet in one translation unit, and both consumers are C.
Under the same define the header appends a test field to `dbarts_forest_info` (`size_t testForestIndex`, the forest's
index plus one), which the appended library writes. The unit compiled WITHOUT the define reads the facts of the two-forest
sampler the appended unit created, through the appended library at the 1.0 stride: each element's `numTrees` is the forest's, and a sentinel the caller placed after the
last element survives, so the library wrote no field past the stride the caller gave. The unit compiled WITH it reads both
fields of each element at its own stride.

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
`std::bad_alloc` and `std::length_error` to OUT_OF_MEMORY, then any other `std::logic_error` to INVALID_ARGUMENT (S-A's run
started inside its own run, as from a draw callback that calls `dbarts_sampler_run` on its own sampler; caught after
`std::length_error`, which derives from it), the rest (the two prior-draw exhaustions among them) to FAILED.
Predict's join (S3) rethrows a failed worker's exception in its own category, an allocation failure as an allocation failure,
with the message text it has today, so R's words hold and a C caller sees OUT_OF_MEMORY.

Per entry (the infallible ones keep returning their value, dec-B571), with what a failure leaves:

| entry | returns | can fail with | after a failure |
|---|---|---|---|
| `dbarts_sampler_create` | status, handle through an out pointer | INVALID_ARGUMENT, UNSUPPORTED (multinomial), OUT_OF_MEMORY, FAILED | `*sampler` NULL, nothing allocated |
| `dbarts_sampler_run` | status | INVALID_ARGUMENT (structSize 0; a run started inside a run of the same sampler), INTERRUPTED, CALLBACK_FAILED, OUT_OF_MEMORY, FAILED; STOPPED | the chains stand where they stopped, a valid state to run on from; this call's results discarded and the saved trees it wrote dropped, the store serving the earlier runs' draws it did not overwrite (S-A); after OUT_OF_MEMORY each chain's fits re-summed from its trees, so it runs on from there (S-A); after FAILED, valid to read and destroy |
| `dbarts_sampler_sampleTreesFromPrior` | status | FAILED (every draw of a tree left an empty leaf), OUT_OF_MEMORY | trees already drawn stand, the tree in flight settled (S-A); drawing again is valid |
| `dbarts_sampler_setResponse`, `setOffset` | status | UNSUPPORTED, INVALID_ARGUMENT (support, finiteness, updateScale refusals); SCALE_KEPT | unchanged |
| `dbarts_sampler_setSigma` | status | UNSUPPORTED (pinned), INVALID_ARGUMENT (not finite and positive) | unchanged |
| `dbarts_sampler_getSigmas` (new) | status | UNSUPPORTED (heteroscedastic, as R's getSigmas() returns NULL) | out untouched |
| `dbarts_sampler_getLatents` | status | UNSUPPORTED (augments nothing) | out untouched |
| `dbarts_sampler_predict` | status | UNSUPPORTED (blend undefined), INVALID_ARGUMENT (source, offset, empty store), OUT_OF_MEMORY | sampler unchanged, out unspecified |
| `dbarts_sampler_setTreeStorage` | status | OUT_OF_MEMORY | the store empty, as any change of storage leaves it |
| `dbarts_sampler_printTrees` | status | INVALID_ARGUMENT (an index out of range) | nothing printed |
| `dbarts_sampler_setVerbose` | status | INVALID_ARGUMENT (printEvery 0) | unchanged |
| `dbarts_sampler_getInfo` (new, dec-B615, dec-B639) | status | INVALID_ARGUMENT (structSize 0 or below the first release's; a non-NULL `forests` with `forestsLength` above 0 and a `forestStride` below the first release's element size, 0 included) | info and the caller's array untouched |
| `destroy`, `setDrawCallback`, `setNumThreads`, the version triple, `dbarts_lastErrorMessage`, `dbarts_statusString`, `dbarts_sampler_spec_init`, `dbarts_sampler_info_init` | value, void, or OK | nothing (the two inits: INVALID_ARGUMENT on a structSize below the first release's) | - |

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
  `dbarts_sampler_run` returns STOPPED, a positive status, with this call's results discarded and the saved trees it wrote
  dropped, the store serving the earlier runs' draws it did not overwrite, exactly as on the interrupt path (S-A's guard in
  [`Sampler::run`](../../src/bartcore/sampler.hpp), landed, which the header already states for a nonzero return). A NEGATIVE return
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

### Build-time math, printing and interrupts (S1, landed; S3 for the poll)

One preprocessor switch, `DBARTS_HOST_STANDALONE`, defined by the standalone recipe and never by the R package's Makevars, so
the R build compiles as it does today.

As built (dec-A260; the bullets below are the design it followed): the six names live in
[nmath.h](../../src/include/external/nmath.h), which stats.h's convenience names sit over; four of R's headers are needed
under `MATHLIB_STANDALONE` and reduced (Rmath.h, Rconfig.h, R_ext/RS.h, R_ext/Print.h, 395 lines), not eight; a prefix header
stands in for R's config.h and renames every symbol the copy defines, the build failing if one escapes; the copy is 6,329
lines, byte-identical to CRAN's R-4.6.2 source; misc builds through its own Makefile from another directory (vpath); R's
unused chi-squared probability went. The review added the `standalone-config` check, which on the first x86 leg found
configure's SSE4.1 macro never defined (fixed on main, no draw moving), and flag checks of every ext_ name in the standalone
driver. Under R the shared library differs from main's before S1 only in two print-pointer stores in its initializer and one
data offset in five run instantiations.

A seventh name, as frequency weights' S3 adds `ext_lchoose` (S1's review gave the recipe): R 4.6.2's choose.c, 158 lines, whose
callees the copy already holds; three renames in prefix.h (`choose`, `lchoose`, `Rf_lfastchoose`); the reduced Rmath.h's
`choose` and `lchoose` declarations restored; `ext_lchoose` in both branches of nmath.h; `Rf_lchoose` in the `r-math-symbols`
list once the engine calls it; the clean recipe and the README. Its grid, about 3,000 points: integer n from 0 to about 1e7
with k at 0, 1, n/2, n - 1 and n; k > n and k < 0 giving -Inf; every k for small n; non-integer and negative n; NaN. No
non-integer k: R's rounding warning goes through `Rf_warning`, which ends the test binary. The nbinom research plan's pchisq
is the same size (pchisq.c, its pgamma already copied).

- Math: src/include/external/stats.h gains the six names with every flag spelled: `ext_pnorm(x, mu, sd, lowerTail, logP)`,
  `ext_qnorm`, `ext_dnorm(x, mu, sd, logP)`, `ext_dt(x, df, logP)`, `ext_qt(p, df, lowerTail, logP)`, `ext_qchisq`. Under R each
  is a macro for the Rmath name the code calls today (`Rf_pnorm5` and so on), so the object code is the same; standalone each is a
  macro for a renamed copy, `dbarts_nmath_pnorm5` and so on. The existing convenience names (`ext_quantileOfNormal`, ...) are
  defined over them. A slim src/include/external/nmath.h declares the six without pulling Rmath.h, for random.c, which today
  declares two by hand. model.hpp's fourteen direct calls become the ext_ names.
- The copy: standalone/nmath/, a new top-level directory listed in .Rbuildignore, so neither the package build nor the tarball
  sees it. It holds the closure of the six in R 4.6.x's nmath (4.6.2 probed: 25 files, 5,963 lines, bd0.c to toms708.c, with
  nmath.h and dpq.h), unedited except a force-included prefix header that renames every external symbol so a program that also
  links libRmath sees no clash; reduced copies of the R headers it needs (four as built, 395 lines);
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
The facts struct's `usesVarianceForest` predicts that refusal, as it predicts setSigma's, and the NaN in every draw's sigma
and in `results.sigma` there, which S-B landed (dec-B601, dec-A270).

### The sampler's facts (dec-B615)

One struct, filled by one call, in place of the eleven one-fact entries (the eight counts, `kIsSampled`, `usesDart` and
`family`), with dec-B601's variance forest and the Student-t residual as fields of their own, and each forest reported
through the caller's array of per-forest fact structs at a stated element size (dec-B639, which ruled open call 8 as (d)):

```c
typedef struct dbarts_forest_info_t {  /* one forest's facts; no size field: the parent states the stride */
  size_t numTrees;             /* the forest's tree count; dbarts_forest_spec's first field */
  /* 1.0-0 field boundary: appends go below, never above */
} dbarts_forest_info;
#define DBARTS_FOREST_INFO_INIT { 0 }

typedef struct dbarts_sampler_info_t {
  size_t structSize;           /* the caller sets it: sizeof(dbarts_sampler_info) as compiled */
  size_t forestsLength;        /* the caller sets it: elements its array holds; the first min(this, numForests) are written */
  size_t forestStride;         /* the caller sets it: sizeof(dbarts_forest_info) as compiled; the INIT macro sets it */
  dbarts_forest_info* forests; /* the caller's array, NULL none: forest f at byte f * forestStride */
  size_t numObservations;
  size_t numPredictors;
  size_t numTestObservations;
  size_t numChains;
  size_t numForests;           /* 1; 2 on BCF; K on multinomial; a variance forest is not counted */
  size_t numSavedSamples;
  size_t numFittedValuesPerObservation;  /* F */
  size_t numVariableCountForests;        /* V */
  int family;                  /* a dbarts_family, 1 to 6 or 8: never AUTO, never STUDENT */
  int kIsSampled;              /* a run writes results.k: k drawn under the model as it now stands */
  int usesDart;
  int usesVarianceForest;      /* dec-B601: sigma is NaN, setSigma and getSigmas UNSUPPORTED */
  int usesStudentResidual;     /* a Student-t residual: a run writes results.residualDf */
  /* 1.0-0 field boundary: appends go below, never above */
} dbarts_sampler_info;
```

- The entry: `dbarts_status dbarts_sampler_getInfo(const dbarts_sampler* sampler, dbarts_sampler_info* info)`, with
  `DBARTS_SAMPLER_INFO_INIT` (every member spelled, structSize and forestStride set to the caller's sizes, the rest zero)
  and `dbarts_sampler_info_init(info, structSize)` for bindings, as dec-B582 gives the description; the function sets
  structSize and zeroes the rest, `forestStride` included, since it cannot know the element size the binding compiled
  (Calls made, 47). The direction is `dbarts_results`': the caller sizes, the library writes each value field whose end lies
  within `structSize` and never writes the four caller fields. A `structSize` of 0 or below the first release's is
  INVALID_ARGUMENT with `info` untouched (`dbarts_sampler_getInfo: info.structSize is 0 - set it to
  sizeof(dbarts_sampler_info) (e.g. dbarts_sampler_info i = DBARTS_SAMPLER_INFO_INIT)`).
- The forests (dec-B639): with `forests` NULL or `forestsLength` 0 no element is written and the other fields are. Otherwise
  `forestStride` must be at least the first release's element size, `sizeof(size_t)`, else INVALID_ARGUMENT with nothing
  written (`dbarts_sampler_getInfo: info.forestStride is 0 - set it to sizeof(dbarts_forest_info)`, or `info.forestStride %zu
  is smaller than any dbarts_forest_info`); the checks come before any write. Then for each f below min(`forestsLength`,
  `numForests`) the library writes, at byte f times `forestStride`, each element field it knows whose end lies within
  `forestStride`, by `memcpy` so a stride a binding computed need not be aligned, and leaves every other byte of the array:
  the elements past `numForests`, and the tail of each element past the fields it knows. A caller that hands a short array
  reads `numForests` in the same call and calls again with a longer one. A single-forest caller hands an array of one.
- Zero means absent, in this direction: a caller built against a later header that initialized with the INIT reads 0 in a
  field the installed library does not know, so every field appended after 1.0 is defined so that 0 is what a library
  predating it implies (no such feature), invariant 1 of the description read the other way. The same holds per element: the
  library never writes past the fields it knows, so the caller zeroes each element (`DBARTS_FOREST_INFO_INIT`, or zeroed
  memory) and reads 0 in a field an older library left.
- It reads the values as they stand at the call (`SamplerBase::shape` and `numTreesInForest`, nothing cached), so a caller
  reads it again after a call that changes one: a run (`numSavedSamples`), setTreeStorage (the store emptied), R's data and
  test conduits (`numObservations`, `numTestObservations`; a binomial sampler's setData under dec-B628 and, once dec-B640's
  slice lands, a multinomial one's among them), and R's setLeafPrior, setModel and setForest (`kIsSampled`, which S-B made
  follow the model: `Sampler::kIsSampled` reads chain 0's forest calibrations, dec-A270). The header names those; every other
  field is fixed for the sampler's life as far as traced (setForest and setVarianceForest refuse `n.trees` by name, so an
  element's `numTrees` holds; see Not verified).
- `family` keeps the one-fact entry's meaning: a Student-t sampler reports GAUSSIAN with `usesStudentResidual` 1, as dec-B615
  has the residual a field of its own (its option (c), reporting STUDENT, was not taken), while a C caller creates one with
  `family = DBARTS_FAMILY_STUDENT` (Calls made, 21); the header says both. `usesStudentResidual` reads
  `SamplerShape::carriesResidualDf` and `usesVarianceForest` `SamplerShape::hasVarianceForest`.
- The hash: the struct folds as `dbarts_draw` does, its field count, the size_t and pointer fields by offset in pointer units
  and the five ints by name, declaration position and width, since an int's offset in pointer units is not one number on
  every platform; `dbarts_forest_info` folds the same way, its field count and `numTrees` by offset. `DBARTS_C_API_LIST`
  gains `dbarts_sampler_getInfo` and `dbarts_sampler_info_init` and loses the eleven; the element adds no entry.
- The prose that names a removed entry names the field instead: the contract's VALUE list (the version accessors remain),
  the forest paragraph (an index past the last forest is INVALID_ARGUMENT on printTrees, the counts read from the elements'
  `numTrees`),
  `dbarts_results`' F, V and the k, varprobs and residualDf conditions, `dbarts_draw`'s LAYOUT, predict's draw count and
  printTrees' bounds, setSigma's paragraph and S-B's sigma sentence (`usesVarianceForest`), getSigmas.
- What it replaces owes nothing: the eleven entries never shipped and 0.9-34 had no C header (compatibility is owed to shipped
  versions only, dec-B461). A stan4bart binary built against today's header fails at its first removed call with
  `R_GetCCallable`'s "not provided" error, its stubs resolving by name; S4's port merges with the main feature. treatSens
  calls none of them.

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
  X(dbarts_status, dbarts_sampler_info_init, \
    (dbarts_sampler_info* info, size_t structSize), (info, structSize))       /* new, dec-B615 */ \
  X(dbarts_status, dbarts_sampler_getInfo, \
    (const dbarts_sampler* sampler, dbarts_sampler_info* info), (sampler, info)) /* new, dec-B615 */
  /* removed (dec-B615): numObservations, numPredictors, numTestObservations, numChains, numTrees, numSavedSamples,
     kIsSampled, usesDart, family, numFittedValuesPerObservation, numVariableCountForests */

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

/* the R view only (DBARTS_USE_STUBS), dec-B578 */
#define DBARTS_C_API_R_LIST(X) \
  X(void, dbarts_R_raiseStatus, (dbarts_status status), (status))
#define DBARTS_R_CHECK(call) \
  do { dbarts_status dbarts_s_ = (call); if (dbarts_s_ < 0) dbarts_R_raiseStatus(dbarts_s_); } while (0)

/* dbarts_draw_callback keeps its type; its return now has three meanings (dec-B581) */
typedef int (*dbarts_draw_callback)(void* context, const dbarts_draw* draw);  /* 0 go, > 0 stop, < 0 failed */
```

Types: `DBARTS_STATUS_LIST` and `dbarts_status` added (the list takes three arguments, so it gets its own enumerator macro beside
`DBARTS_ENUMERATOR`, and the library's word table expands the third), `DBARTS_SWAP_STATUS_LIST` and `dbarts_swap_status` removed;
`DBARTS_WEIGHTS_TYPE_LIST` and `dbarts_weights_type` added; `dbarts_sampler_spec` and `DBARTS_SAMPLER_SPEC_INIT` added;
`dbarts_sampler_info` and `DBARTS_SAMPLER_INFO_INIT`, and before them `dbarts_forest_info` and `DBARTS_FOREST_INFO_INIT`
(dec-B639), added ([The sampler's facts (dec-B615)](#the-samplers-facts-dec-b615)); in
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
R_alloc"; `dbarts_sampler_getLatents` keeps latents-first's sentence (dec-B586: the values the last sweep's trees were drawn
from) and gains frequency weights' sentence on a frequency row's latent (dec-B564's note); every
per-entry "raises" becomes the status it returns; the spec's block states invariants 1 to 7 of the forward-compatibility
section, the reserved family values and the zero-tail rule; the eleven one-fact entries' docs go, the info struct's block
taking their sentences (S-B's kIsSampled sentence among them, and its "re-reads dbarts_sampler_kIsSampled before sizing it"
in `dbarts_results`' paragraph becoming a re-read of the info struct) and every cross reference re-pointed to its fields;
`dbarts_forest_info`'s block states the stride rule, the zeroed element and the mirror of the creation element
([The sampler's facts (dec-B615)](#the-samplers-facts-dec-b615)); `dbarts_draw_callback` says a run started from inside it on its own
sampler returns INVALID_ARGUMENT. The text S-A and latents-first landed (a stopped run drops the saved trees it wrote, an
allocation failure leaves the chains runnable, the draws that open a sweep) and S-B's (the channels a callback receives,
sigma NaN under a variance forest, no "host" for the caller) is carried into the status form as it stands. S-B's docs-10 share
left no "host" in dbarts.h (a grep at 7e7edb526 finds none), so the release map's host-wording item has nothing left for S3;
the new prose says "the caller" or "the code that runs the BART sampler inside a larger sampler", never "host". About
+700/-240 lines in the header.

The hash: [`dbarts_foldLayout`](../../src/C_interface.cpp) gains the spec's size and every field's name and offset in pointer
units, as `dbarts_results`, and `DBARTS_STATUS_LIST` (values and words) and `DBARTS_WEIGHTS_TYPE_LIST` join the enum fold;
`dbarts_statusString` enters the hash as a list entry like any other; both tokens are re-baked
([`dbarts_apiSignatureToken`](../../src/C_interface.cpp) and the literal). Several by-value doubles in the spec make its offsets
differ by platform in pointer units, as `dbarts_draw`'s do, so the spec folds as `dbarts_draw` does: field count, and doubles by
name, position and width. The info struct folds the same way, its ints by name, position and width, and `dbarts_forest_info`
by its field count and offsets; the INIT macros are not folded, as `DBARTS_RESULTS_INIT` is not. tools/check-api-hash.sh is
unchanged and still skips (no release tag).

Development binaries: a stan4bart or treatSens binary built against today's header still loads against the new dbarts, since
the version pair stays 1/0 and neither defines `DBARTS_REQUIRE_EXACT_ABI`; it ignores every status. A stan4bart binary then
stops at its first call of a removed accessor, which its stub resolves by name and `R_GetCCallable` no longer finds (dec-B615
removes them); a treatSens binary calls none and runs on. Only shipped versions are owed
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
- `trials is read under DBARTS_FAMILY_LOGISTIC only; family %d takes no trial counts`; the binomial count refusals (successes
  above the trials, a negative or fractional count) in frequency weights' words;
- `a gaussian, Student-t or aft sampler needs a finite positive sigmaEstimate`; `sigmaPriorDf and sigmaPriorScale must be
  finite and positive`;
- `numChains must be at least 1`; `printEvery must be at least 1`; `no seed: set spec.hasSeed and spec.seed, or
  spec.chainSeeds`;
- the source refusals of `translateSource` with this entry's name, the test source's under `xTest`;
- `block tree counts must sum to numTrees`; and the engine's own creation refusals (leaf covariates, variance forest, monotone
  with a proposal mix, storage single), worded as the engine words them, with the entry's prefix.

The build stage also repeats, in its own words, a refusal the R parse raises first in R's: S-D's variance calibration
(dec-A264's models-01: a variance forest of two or more trees under a chi-squared residual prior of df 2 or below), `a
variance forest of %zu trees needs sigmaPriorDf above 2`, which no 1.0 description reaches (the variance forest is an append)
and the extension's variance group does. `dbarts_sampler_getInfo` and `dbarts_sampler_info_init` refuse a `structSize` as
the spec's entries do (`info.structSize is 0 - ...`, `structSize %zu is smaller than any dbarts_sampler_info`), and getInfo a
`forestStride` below the first release's element size where it would write an element (`info.forestStride is 0 - set it to
sizeof(dbarts_forest_info)`, `info.forestStride %zu is smaller than any dbarts_forest_info`). A run
started inside a run of the same sampler returns INVALID_ARGUMENT with S-A's engine words under the entry's prefix.

The run's callback failures (CALLBACK_FAILED), under `dbarts_sampler_run`'s prefix: `the draw callback for chain %zu returned
%d`; `the draw callback for chain %zu threw: %s`; `the draw callback for chain %zu threw an exception not derived from
std::exception`. A deliberate stop (STOPPED) is not a failure and writes no message. Every status's fixed words are the
`dbarts_statusString` table above, which is also the text of a failure that left no message.

## RNG class

Neutral on the R build, by construction and by gate:

- S1, landed, changed which macro spells each math call, not which function runs: under R every `ext_` name is the Rmath
  symbol the code called before, which the trio bitwise and the `r-math-symbols` check confirmed. S2a, S2b and S3 move code
  without changing what is computed or in which order the generators are touched; the one place order could move, the chain seeds, keeps both of today's paths (one seed through a seed generator, or seeds drawn
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
  grid is gated at a relative difference of 1e-13 and the consumer's draws by the equivalence harness's statistical z. As
  built, S1's grid (36,093 points from its own generator, every flag pair) held bitwise per routine on the reference build and
  on the x86 leg alike; the floor, about 450 ulp, catches only gross errors there, so small ones are caught by the bitwise gate
  on arm64 macOS, which CI's cpp-tests runs on every push.

## Ruled

The maintainer ruled open calls 1 to 7 of this plan on 2026-11-03, and added two rulings it now builds (dec-B579 and the words
clause of dec-B578); the numbers are the plan's open-call numbers, which the ledger's records cite. Item 8 was ruled on
2026-11-06 from the fourth review, outside this plan's calls, and revises item 1's accessors; item 9 is this plan's open call
8, the one choice item 8 left to it that a caller sees, ruled on 2026-11-07:

1. dec-B571: in the shipped header only an entry that can fail returns a status; the counts and the other accessors return their
   value, and numTrees moves to a status since its forest index can be out of range. The maintainer: "For the host-neutral
   decision, use (a)." Revised by dec-B615 (item 8): the counts and the facts are no longer entries of their own.
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
8. dec-B615, from the critique of the fourth review's slice S-B, which found no entry saying whether a run writes
   results.residualDf: the header describes a sampler through one struct led by structSize, with an initializer macro and
   function as dec-B582's description has, filled by one entry that reads the current values; it carries every count and fact
   the one-fact entries report, dec-B601's variance forest and the Student-t residual as fields of their own; its field order
   is fixed at 1.0 and grows only at the end; the one-fact entries are removed, stan4bart's calls moving to it in the same
   landing. "The entry's and struct's names, and how the per-forest tree counts are carried, come with the plan." The
   maintainer: "Yes, let's do that." Built in [The sampler's facts (dec-B615)](#the-samplers-facts-dec-b615), S3 and S4; the names are
   Calls made, 38, the tree counts item 9.
9. dec-B639, open call 8 (how the facts struct carries each forest's tree count): the facts struct carries the forest count
   and a pointer to the caller's array of per-forest fact structs, with their count and element size, versioned as the
   sketched several-forest creation array is (invariant 7, a stride stated once by the parent); the call fills each element
   up to the count the caller gave, at the size the caller compiled, and leaves the rest; at 1.0 an element holds a forest's
   tree count, and later per-forest facts append to the element under the header's size rules; the element's shape is
   designed beside the sketched `dbarts_forest_spec`, so a forest's reported facts mirror what creating it states. The plan
   had recommended (a), a caller's array of tree counts with its length. The maintainer asked "Is (d) inevitable, based on
   planned work?", was shown that no plan or TODO schedules creating several-forest samplers from C and that a per-forest
   array appended later would report a forest's tree count twice, and ruled: "Let's just do (d) now, and think ahead." Built
   in [The sampler's facts (dec-B615)](#the-samplers-facts-dec-b615) and
   [Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial), S3 and S4;
   the names and the remaining details are Calls made, 45 to 49.

## Calls made

1. The vendored set is the closure of all six routines the engine and the support library call, not the four TODO names: qt
   ([`forestMapEntry`](../../src/bartcore/combiner.hpp)) and qchisq at 0.5
   ([`ChiKHyperprior::medianOfChi`](../../src/bartcore/model.hpp)) arrived after the ruling and are read at creation and at live
   restatements (setLeafPrior, and coefficient-prior-live's `$setForest(amplitude = )`, built). Built so in S1. Alternative:
   compute them on R's side as the residual prior's quantile is, which cuts the copy by about 4,000 lines but adds an
   R-resolved number to every live
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
26. S1 is its own feature (landed so); the rest is one branch merged once, the header change landing whole in S3 so stan4bart
    and treatSens rebuild once. Alternative: S1 on the main branch, which holds a byte-neutral change with no dependency behind
    forest-door.
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
35. approximate-draws' latent sampler setting gets no 1.0 field. Its shape is ruled (dec-B593): a method by name, NULL the
    automatic rule (dec-B590's per-row rule until the defaults study, the study's after, dec-B592), and a separate tolerance
    read only by an approximate method; the names are settled at that plan's S5. It is a setting of how the sampler draws,
    beside the move settings dec-B575 leaves for an append, and its zeros read as R's NULL and the default tolerance, so the
    append is a plain one under invariants 1 to 3; a C caller at 1.0 gets the automatic rule, the default an R user gets.
    approximate-draws lands first in the order below, so S2b carries the method and the tolerance into the extension's
    public-shaped part, laid out as the append will ship. Alternative: a 1.0 field now, freezing names the maintainer settles
    at that plan's S5.
36. `refuseCountOverCap`, the `[[noreturn]]` helper beneath `enforceBinaryWeightPolicy` and `validateResponseSupport`, returns
    a refusal value in S2a like the other validators (nineteen on frequency weights' tip), its R sites wrapped. Alternative: leave it raising, which keeps an
    `Rf_error` inside two validators the R-free core holds, so the flat entries and the standalone build could not compile them.
37. `trials` is a 1.0 field read under LOGISTIC only and refused elsewhere, in frequency weights' form (`y` the successes, NULL a
    0/1 response). Alternative: ignore it off LOGISTIC, as `weightsType` is ignored without weights; but a `trials` given to a
    probit sampler states a model the build does not fit, and refusing it is what lets a later probit binomial read it.
38. The names (dec-B615 leaves them to the plan): the struct `dbarts_sampler_info`, the entry `dbarts_sampler_getInfo`, as
    `getSigmas` and `getLatents` read a sampler, and `DBARTS_SAMPLER_INFO_INIT` and `dbarts_sampler_info_init`, as the
    description's. The two new facts are `usesVarianceForest`, dec-B601's own name, and `usesStudentResidual`, beside
    `usesDart`. Alternative: `dbarts_sampler_describe`, which reads as making a description to create from, the spec's job.
39. The struct carries `numForests`, which no one-fact entry reports (`dbarts_draw` does, inside a run), since the tree counts
    are read against it and a several-forest sampler built from R is otherwise not countable from C before a run.
    Alternative: leave it to `dbarts_draw`, so a caller of a BCF sampler sizes nothing per forest until its first draw.
40. `family` keeps the one-fact entry's meaning, never STUDENT, and `usesStudentResidual` says a Student-t residual is in
    force, while creation takes `DBARTS_FAMILY_STUDENT` (call 21); the header states the pair. Alternative: report STUDENT,
    dec-B615's option (c), not taken.
41. Field order: the caller's fields first, then the counts in today's entry order with `numForests` after `numChains`, then
    the int facts, so the size_t block folds by offset and the ints by name. Alternative: every fact as a size_t, folding as
    `dbarts_results` does, which spells a yes-or-no as a count.
42. Fields appended after 1.0 read 0 from a library predating them, through the INIT's zeroes, so each is defined so that 0
    means what that library implies. Alternative: the library writes its own size into `structSize`, as `dbarts_draw` does,
    which changes an input the caller passes again on its next read.
43. A plain `std::logic_error` (S-A's run started inside its own run) maps to INVALID_ARGUMENT, caught after
    `std::length_error`. Alternative: FAILED by the catch-all, which tells a caller who called the entry wrongly that the
    engine failed.
44. approximate-draws' close-call count, if the maintainer rules at its S2 that dbarts.h reports it (that plan's Decision),
    is a field S3 adds, to `dbarts_results` if it is per run and to the info struct if it is the sampler's, never an entry
    of that plan's S2 that S3 then moves, so the header changes once. Alternative: an entry at its S2, two ABI events and an
    entry dec-B615's shape would retire.
45. The names under dec-B639: the element `dbarts_forest_info`, as `dbarts_sampler_info` reports what `dbarts_sampler_spec`
    states and `dbarts_forest_info` what the sketched `dbarts_forest_spec` states; the parent's `forests` and `forestStride`,
    the sketch's own names; `forestsLength` for the elements the caller's array holds, since `numForests` is the sampler's
    count the library writes; `DBARTS_FOREST_INFO_INIT`. Its one 1.0 field is `numTrees`, `dbarts_forest_spec`'s first
    field's name and type. Alternatives: `dbarts_forest_facts`, which breaks the spec-info pairing; `numForestsGiven` for the
    length, a second "numForests" a reader must tell from the library's.
46. The caller's fields follow `structSize` in the sketch's order, `forestsLength`, `forestStride`, `forests`, then the
    library's fields as call 41 has them. Alternative: the pointer first, as the earlier tree-count buffer was, which orders
    the same three facts differently in the two structs.
47. `DBARTS_SAMPLER_INFO_INIT` sets `forestStride` to `sizeof(dbarts_forest_info)` as the caller compiles it;
    `dbarts_sampler_info_init` leaves it 0, since a library cannot know the element size a binding compiled, and a binding
    that forgets it is refused, never overrun. Alternatives: neither sets it, one more line for every C caller; or the
    function writes the library's element size, which on a library newer than the binding writes past each of the binding's
    elements.
48. No `dbarts_forest_info_init` function: every element default is 0 and a binding zeroes its memory, while the sketched
    `dbarts_forest_spec_init` exists for the creation element's nonzero defaults (base, power). Alternative: the function
    for symmetry, an entry and a hash event whose body is a memset.
49. The library writes each element field by `memcpy` at its offset and takes any stride at least the first release's
    element size; the stride is checked only where an element would be written (`forests` non-NULL and `forestsLength`
    above 0), so a caller that wants no forest leaves all three at 0. Alternatives: refuse a stride that is not a multiple
    of the element's alignment, a refusal no compiled caller meets; check the stride always, which refuses the function
    path's zeroed struct with no array.

## Open calls

None. Open call 8, the last, was ruled on 2026-11-07 as (d) (dec-B639; [Ruled](#ruled), 9). The second 2026-11-07
amendment raised no call a user sees that the rulings leave open: the per-forest element's names and details are the plan's
under dec-B615's "The entry's and struct's names ... come with the plan" and dec-B639's "designed beside the sketched
dbarts_forest_spec" (Calls made, 45 to 49), and the landings and work in flight it placed change no R surface and no
header meaning beyond what their own rulings state.

## Slices

### S1 (its own feature): build-time math and printing - LANDED

Landed on main 2026-11-03, cbed014fb to 789ad3ff1, recorded in d6ef1cc13 (dec-A260); built as
[Build-time math, printing and interrupts (S1, landed; S3 for the poll)](#build-time-math-printing-and-interrupts-s1-landed-s3-for-the-poll)
says, with the review's fix round (the `standalone-config` check, the standalone driver's flag checks, prefix.h's note on
cospi, the README's warnings and license lines) and the orchestrator's configure fix. Gates as dec-A260 records them: tests/cpp
426 ok, plain and under ASan and UBSan, on the reference build and on x86; the grid bitwise per routine on both; the
standalone target from a clean archive referencing nothing of R; tinytest; the trio bitwise on the reference build; the four
snapshot files; R CMD check --as-cran, no standalone/ file in the tarball; CI. Lines: about 910 written, 6,329 copied, 395 of
reduced headers, against ~350, ~6,000 and ~400 planned. The review ran mutants 1 to 4 ([Mutants](#mutants)).

What S2a starts from: the engine and the support library call R only through the six ext_ names and io.h's two print macros,
and randomBase.c's R generators compile under R only; standalone/ holds the copy, io_standalone.c and the checked-in
configuration headers; tests/cpp builds misc, external and the copy standalone but nothing above them. Still on R, for S2a to
S3: the shared refusals' `Rf_error`, the flat file's `R_alloc` and `R_finite`, the poll, the callback leaf and the slow-count
warning in the bridge's core, and `predictPartition`'s main-thread rule. No host_r.cpp or host_standalone.cpp exists yet (S3).

### S2a: refusals as values, one translate stage (bridge)

Steps: `Refusal` and `raiseRefusal`; the shared validators, nineteen with frequency weights' three
(`refuseWeightedCountsOverCap`, `refuseWeightedTrialsOverCap`, `refuseBinomialResponseSwap`), and `refuseCountOverCap` beneath
them (Calls made, 36), and the flat file's `translateSource`, `validateTestSource` and `refuseTestMissingness`, returning
refusals, their named lines in the bridge and the flat file wrapped (74 and 12 on main f79874ea8, 90 and 13 on frequency
weights' S3 tip by the first count; 90 and 20 on main 7e7edb526 and 108 and 21 on frequency weights' tip ab739157c by the
second amendment's wider grep, Context; the implementer recounts on the base S2a starts from, by one rule) (createDataHandle and the R test-predictor setters call the last
three inside `unwindProtect` bodies, outside any route's catch); the translated source, which owns storage, held as an
init-capture wherever it lives in an R frame, by `unwindProtect`'s own rule; R_interface_bartcore_common.hpp split into the
R-free core and the SEXP-facing half; `translateSource` moved into the core, used by the R parse's views and the data handle,
carrying the engine's per-column channels (the shipped field waits for S3); the engine routes' catch turned into a refusal
value. What it takes as it stands: forest-door's and coefficient-prior-live's live forest entries
(`bartcore_setForestSd`, `bartcore_setForestCoefficient`, `bartcore_setForestTreePrior`), which call no shared refusal and keep
their R-only `Rf_error`s; frequency weights' lifted and widened refusals (`enforceBinaryWeightPolicy` taking 0/1 on probit and
ordinal and a logistic 0, `refuseCountOverCap`'s product cap, `validateResponseSupport` taking successes up to the trials),
converted as that feature leaves them; latents-first's setState without its digest and zero-weight reconciliation, which
holds no shared refusal. The fourth review's bridge code, none of it a shared refusal, keeps its R-only `Rf_error`s and moves
to the SEXP-facing half as it stands: S-D's `refuseVarianceCalibration` and parseModel's infinite-df refusal (raised in the
R parse before the spec fill), S-A's sweep-closure checks (`checkSweepBuffer` and the shape stop in
`bartcore_runWithCallback`, an R entry, not a flat one) and the warm-start value refusal, S-C's
`bartcore_sampleHyperparametersFromPrior` (landed, dec-A269: its captured engine refusal and attempt-count refusal, raised
after every C++ object of its frame is gone, and its poll), and S-B's setCutPoints floor and repeated-cut refusal (landed,
dec-A270). The work landing before S2a, each taken as it stands: rchk-findings' protection of the model's slot reads
(PROTECT and UNPROTECT around `slotOrNull` in parseModel, parseDartSettings, the move-record readers, the variance and
forest record readers, `refuseVarianceCalibration` and `forestMessageName`), which S2a's wrapping must keep balanced where a
wrapped call sits between them; frequency weights' binomial setData (dec-B628) and nbinom's counts reading under every
creation (dec-B632, the bridge's creation-time switch removed), R-only entries and no shared refusal; state-basis-levels'
`bartcore_setState` argument and glue reorder (dec-B647), moved to the SEXP-facing half as it stands, wherever its builder
puts the level record (its plan has an R attribute; the engine's state, its alternative, holds no shared refusal either);
and the dec-B640 slice's multinomial setData (not yet planned), which must let a multinomial sampler past
`refuseMultiForestMutation` at setData while BCF stays refused, and replace the counts and trials the holder owns
(`ownedCounts`, `ownedTrials`), which the softmax combiner borrows; S2a converts that shared refusal and moves the holder
into the R-free core as the slice leaves them.

Gates: tests/cpp, every target; the full tinytest suite in one process, every message pin unchanged (frequency weights' and
hurdle-completion's new pins among them); the trio bitwise; the four snapshots on the reference build; R CMD check --as-cran
from a clean tarball, no standalone/ file in it; sanitizers on tests/cpp and the R-loaded ASan run over the creation and
mutation suites; mutation-battery verify-anchors; rchk dispatched on the branch tip, no finding (main's are zero once
rchk-findings lands); the x86 leg.

Budget: ~1,000 lines. Agent ~13.5, review ~5, gates ~4.

### S2b: the plain-C specification (bridge)

Steps: `dbarts_sampler_spec`, its INIT and `DBARTS_WEIGHTS_TYPE_LIST` in dbarts.h (no entry yet); src/spec_internal.h with
`SpecExtension`; src/sampler_build.cpp with `buildSampler`, taking over the engine-building half of createHolder,
optionsFromParsed, the two multinomial routes and createFromHandle; the R parse filling spec and extension (the parse functions
keep reading SEXPs; their outputs become spec fields); seeds through `hasSeed`/`seed` or `chainSeeds`; the weights reading,
frequency weights' resolved model slot becoming `weightsType` and its resolution table the build's check; frequency weights'
binomial counts slot under logistic becoming `y` (the successes) and `trials`; a tests/cpp suite that fills every spec field
with a distinct value and checks the engine options built from it field by field, so a field the build forgets fails by name
(`weightsType` and `trials` among them); the extension's public-shaped part (`ForestRecord`, `MultinomialRecord` and the
one-forest groups 1.0 leaves out, approximate-draws' method and tolerance with them, Calls made 35) laid out
as the sketched append, the several-forest and multinomial R routes reading it; the 1.0 reading rules (zero tail, reserved
families, fields read only under the families that name them, `trials` off LOGISTIC refused). The forest records are as on
main: coefficient-prior-live writes `held.amplitude`, `amplitude.df`, `sd` and `share` at its live door and adds no slot, so
`ForestRecord` mirrors forest-door's records. The build stage repeats S-D's variance calibration in its own words over the
extension's variance group (Refusals), R's words staying in the parse before the fill. What it takes as it stands, all
landed: S-G's starting sigma regressed on each forest's basis (dec-A267), which R resolves and the spec's `sigmaEstimate`
carries; S-B's cut-grid floor (dec-A270), which `cutGridIsValid` holds in the engine, so the R-only cut-grid record the
extension carries needs no check of its own; S-C's hyperparameter draw (dec-A269), no creation route and no creation field,
whose exact-prior tests in test-prior-predictive.R read k's prior, DART's six and the response hyperpriors off a created
sampler, so a prior field the build forgets moves them. And as they land before it: rchk-findings' protected parse, which
the rework keeps balanced as its outputs move into spec fields (its own gate, rchk on the branch tip, catches a slip);
nbinom's counts reading under every creation (dec-B632), so under NBINOM every pair the weights check accepts reads
`weights` as counts, 0/1 weights included, with no creation record for it; the dec-B640 slice's multinomial setData, which replaces the holder's
counts and trials after creation, so `MultinomialRecord` borrows them from the holder as the combiner does today, never a
copy taken at the fill.

Gates: as S2a, plus a binomial and a frequency-weighted sampler created through the spec against the same fits created on the
base, bitwise, and a fit under each of approximate-draws' methods likewise; a multinomial sampler created through the spec,
then given `$setData` at another row count (dec-B640), against the same calls on the base.

Budget: ~1,790 lines (header ~220, extension ~300, build stage ~720 moved, parse rework ~430, tests/cpp ~290). Agent ~20.5,
review ~7.5, gates ~5.

### S3: statuses, the creation entry, the sigma getter, the facts struct, the header (bridge, flat entry file, tests)

Steps: `DBARTS_STATUS_LIST` with its words and `dbarts_statusString`; the list as in [The header diff](#the-header-diff);
`dbarts_sampler_create`, `spec_init`, `getSigmas`, `lastErrorMessage` with the overwrite-on-every-failure record;
`dbarts_sampler_info`, its INIT, `info_init` and `getInfo` (dec-B615), with `dbarts_forest_info`, its INIT and the stride
rule (dec-B639), the eleven one-fact entries and their provider-side bindings removed, `std::logic_error` mapped (Calls
made, 43); C_interface.cpp
off R (owned scratch, statuses, the thread's record); `predictPartition` thread_local; predict's join keeping its category;
src/host_r.cpp and standalone/host_standalone.cpp, with the callback catch for the inline call and the worker trampoline; the
engine's three-way stop record (continue, stop, failure with chain and value or text) in place of the one stop flag, and the
lowest-failing-chain rule; the holder's origin and destroy; `dbarts_R_raiseStatus` and `DBARTS_R_CHECK` (dec-B578), forwarding
only neutral text; `denseChannels` in the shipped source; the hash folder and both tokens; the header's prose, the callback
contract rewritten, carrying latents-first's sentences for `dbarts_sampler_getLatents` (what it returns after a run and
between a setter and the next run) and for a stop inside the latent draws, which then open a sweep, into the status form,
and frequency weights' sentence on a frequency row's latent; the spec's block stating `weightsType` and `trials`; the info
struct's block ([The sampler's facts (dec-B615)](#the-samplers-facts-dec-b615)). The stop record sits where latents-first leaves the
sweep: the callback still fires at the record, the sweep's last step, and the interrupt's cancel point inside the latent draws
keeps its arm. In [`Sampler::run`](../../src/bartcore/sampler.hpp) it replaces the `stoppedByCallback` out-parameter and the
relay's `stopped` flag beside what S-A and S-B put there, both landed: S-A's stopped-run guard and `running_` flag, and S-B's
observer scratch (`provideObserverScratch`, one block per channel the fit carries, allocated in the run's frame only when a
callback is registered, from `Chain::carriedChannels`), with the per-draw strides `Results` gained for it (dec-B600,
dec-A270). The flat run keeps S-B's `results.k` rule (left NULL to the engine where `kIsSampled` is 0) under its status
return, and the poll moved to src/host_r.cpp serves its five call sites, S-C's among them (Context). Tests: consumer.c and test-capi.R rewritten for statuses
(each entry's failure reached and its status and message read; INIT plus data and seed against `dbarts()`'s draws bitwise; a
C-created sampler against the R-created one bitwise, a binomial one with `trials` and a gaussian one under FREQUENCY among
them; `denseChannels` with packed blocks against the dense source; getSigmas on two chains and its UNSUPPORTED on a
heteroscedastic sampler; `DBARTS_R_CHECK` raising an error, through the interrupt test hook an interrupt, and under
`setTimeLimit` the time limit's error, as test-monotone.R pins for R's run; `DBARTS_R_CHECK` on a failure whose message is
empty raising `dbarts_statusString`'s words; two samplers predicting on two threads; a C-created handle destroyed; and the four
callback paths, each run on a gaussian sampler of two chains: a positive return on chain 1 (STOPPED, no message, results
discarded, a second run valid from where it stopped), a negative return (CALLBACK_FAILED, the message naming the chain and the
value), a `std::runtime_error` thrown on a worker (two threads, CALLBACK_FAILED with its `what()`, the process alive), and the
same thrown inline (one thread, the same status and message); and, under R on an inline run, an `Rf_error` from the callback
still reaching R as its own error; a callback that calls `dbarts_sampler_run` on its own sampler getting INVALID_ARGUMENT
and the outer run completing), the five handshake arms kept, and S-B's recording callback with its channel arms over the
four routes kept under the status return; the facts struct, every accessor call in consumer.c (50 at 7e7edb526, S-B's arms
added six) read from it through one helper: each field on gaussian, probit, Student-t, BCF, multinomial and variance-forest
samplers against R's view of the same sampler; the forest elements (dec-B639) on BCF and multinomial against R's tree
counts, each array placed between sentinel words: a short array (`forestsLength` 1 on BCF: forest 0 written, the next
slot's sentinel intact, `numForests` 2), a long array (`forestsLength` 4 on BCF: two written, slots 2 and 3 intact), a
larger stride than the library knows (`sizeof(dbarts_forest_info)` plus 16: each element at its offset, its tail intact),
and an older element size (the `spec-append` target's test field, below); `forests` NULL with the other fields written;
`forestsLength` 0 with a stride of 0 accepted; a stride of 0 and a stride of 4 with an array refused with their messages,
the struct and the array untouched; `dbarts_sampler_info_init`'s struct with an array and no stride refused, the INIT's
accepted; `usesVarianceForest` 1 on gaussian and aft with a
variance forest and 0 on gaussian, probit, multinomial and BCF, unchanged after `$setVarianceForest(power = 3)` and on the
re-read handle after `s$copy()` (critique amendment 1's arms, built here); `usesStudentResidual` equal to whether a run
writes a preset `results.residualDf`; S-B's `kIsSampled` transitions through setLeafPrior and setModel read from the struct;
a re-read after a run, setTreeStorage, R's setTestPredictor and a setData of another row count (a gaussian, a binomial
under dec-B628 and a multinomial under dec-B640's slice); a struct whose `structSize`
ends before its last field leaving a sentinel there; `structSize` 0 and below the first release's refused with their
messages; stan4bart-shaped reads (one forest, an array of one); a tests/cpp unit for `dbarts_statusString`, every listed status
against its words and an unlisted one against "unknown status"; tests/cpp's `standalone-consumer` target: a C99 file including the prototype
view, creating gaussian, probit and Student-t samplers, running, predicting, reading sigma, refusing, destroying, built and
linked with no R, `nm` showing no undefined R symbol, and its draws compared with the R-created sampler's at the same spec and
seed; the `spec-append` target of [Forward compatibility with several forests and multinomial](#forward-compatibility-with-several-forests-and-multinomial)
(the test-only append, the 1.0-sized and the appended consumer, the two-forest arm, the nonzero-tail refusal, and the
forest element's test field read at the 1.0 stride and at the appended one); the
cpp-tests workflow runs both. Docs: dbarts-embedding.Rd's C section, the architecture's layering and callbacks sections,
NEWS.

Gates: as S2a, plus test-capi.R under ASan on the R-loaded path, tests/cpp's standalone consumer under ASan and ThreadSanitizer
for the two-thread and worker-exception arms, on the x86 leg LeakSanitizer (C-created handles are new owners); the header-compiles target;
`sh tools/check-api-hash.sh`.

Budget: ~3,330 lines (header +700/-240 with the test-only append ~90, the weights and trials prose and the facts struct with
its forest element, C_interface.cpp ~585 with getInfo, the stride writes and the eleven accessors gone, build units ~290 with
the callback catch, engine stop record ~40, bridge ~120, consumer.c ~790, test-capi.R ~440, tests/cpp ~640, workflow ~20,
docs ~140). Agent ~21, review ~8, gates ~6. The facts struct is ~420 of those lines and ~5 of those hours: the ~270 and 3.5
of the fifth amendment, and dec-B639's element, ~150 lines and 1.5 hours (header ~25, library ~25, consumer.c and
test-capi.R ~60, the spec-append arm ~40), above the ruling's "about 70 lines" by its test arms; with S-B landed, the flat
run's scratch and kIsSampled are read, not estimated, and consumer.c's six more accessor calls are counted.

### S4: the consumer ports (stan4bart, treatSens)

Each in its own checkout, on a branch, against a private library of the S3 tip, its CI's dbarts reference set to the branch
while in review. stan4bart: the fallible calls in src/init.cpp (run twice, sampleTreesFromPrior, setOffset twice, setSigma twice,
getLatents twice, predict, setTreeStorage twice, printTrees, setVerbose) wrapped in `DBARTS_R_CHECK`, and its eleven one-fact
calls (numObservations, numTestObservations, kIsSampled, numPredictors three times, numChains and numSavedSamples twice each,
numTrees of forest 0) read from one `dbarts_sampler_getInfo` per function, printTrees' tree count through an array of one
`dbarts_forest_info` (dec-B639): `dbarts_forest_info forest = DBARTS_FOREST_INFO_INIT;`, the INIT'd info's `forests = &forest`
and `forestsLength = 1`, `DBARTS_R_CHECK(dbarts_sampler_getInfo(fit, &info))`, then `forest.numTrees`, where it calls
`dbarts_sampler_numTrees(fit, 0)` today (src/init.cpp, printTrees, read at 98db852); nothing else changes (the dense helper
is source-compatible). Its R code's empty-design check, which the fourth review's S-F
moves ahead of the call (TODO stan4bart-empty-bart-guard), is R code beside this port, not in its files. treatSens: run three times, setResponse four times,
setOffset, setSigma and setVerbose twice, likewise. Each suite run (stan4bart's tinytest, recompute-exactness and
compare-posterior; treatSens' testthat), R CMD check, and revdep-smoke dispatched on the dbarts branch tip with both branches,
bartCause included, all green. They merge right after dbarts' merge, their dbarts reference back to main.

Budget: stan4bart ~115 lines, treatSens ~50. Agent ~4 and ~3, review folded into S3's reviewer reading both diffs.

### S5: records

TODO's host-neutral entry and python-bindings entry (the engine is R-free below the bridge; what a binding still owns:
ingestion, resolution, its own build unit and the standalone recipe); TODO capi-sampler-info closed (dec-B615 and dec-B639,
built in S3); TODO host-wording-sweep closed if still open (S-B's docs-10 share left no "host" in dbarts.h, and the help, the
vignettes and the R strings were done before it); docs/plans/python-bindings.md's Context; the feature
matrix's reach column for the flat header; INDEX rows; NEWS's 1.0-0 entry for compiled consumers; the Landing note and
dec-A entry for the main feature's calls. Doc-freshness only.

## Tests and gates

Per slice above; across both features, run independently of the implementer on each slice's own library, `--preclean` on every
install (the facade's virtuals do not change, but every creation route does):

```
R CMD INSTALL --preclean -l <lib> .                       # and --configure-args=--enable-reference-build for the snapshots
cd tests/cpp && make && ./test_bartcore                   # S1's all: header-compiles, the libR suite and nmath grid,
                                                          # r-math-symbols, standalone-config, standalone-math
make standalone-consumer && ./build-standalone/standalone_consumer   # S3
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
gh workflow run rchk --ref host-neutral-engine             # no finding; the bridge's PROTECT balance after each slice
```

plus sanitizers on tests/cpp and on the R-loaded path ([Gate hygiene](README.md#gate-hygiene)), and the x86 leg before landing
each slice (tests/cpp plain and with ASan, UBSan and LeakSanitizer, tinytest, the trio recorded on the base tip there and
compared --bitwise, the nmath grid at the relative floor).

## Mutants

The reviewer of each slice runs these; each must fail the named gate, and is reverted with a `touch`. S1's four were run by
its review and are recorded as run:

1. S1: one model.hpp pnorm site's lower-tail flag flipped. The plan named the probit and aft equivalence scenarios, a wrong
   premise: in probit and aft the engine's pnorm feeds only the log-likelihood channel, which the equivalence harnesses do not
   record, and their draws come from random.c's truncated normal samplers. As run: (a) the probit log-likelihood's flag
   swapped left the equivalence scenarios identical and was killed by tests/cpp and test-capi.R's log-likelihood arms; (b) the
   near-tail flag in `logStandardNormalMass`, which the monotone leaves read, moved monotoneleaf and monotonejoint (|z| 3.28
   and 1.83) and failed six tests/cpp checks. A later slice that moves a pnorm site gates it by those, not by the probit or
   aft scenarios.
2. S1: a constant in the copied qnorm changed - the nmath grid (as run: qnorm 9 points at 1 ulp, qt 4 at 3 ulp).
3. S1: the copy compiled with `-ffp-contract=off` - the nmath grid's bitwise gate on the reference host (as run: 8,673 of
   36,093 points off; the 1e-13 floor passes it).
4. S1: `ext_qnorm` mapped to the copy in the R build - no scenario moves (as run: 55 of 55 identical), so the gate is
   `r-math-symbols` reading the R build's symbols by `nm`, which fails.
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
31. S2b: `trials` not handed to the engine, a binomial read as 0/1 - the field coverage suite, and S2b's binomial creation twin.
32. S2b: a non-NULL `trials` accepted off LOGISTIC - a tests/cpp unit of the build stage (S2b) and consumer.c's trials refusal
    arm (S3).
33. S2a: `refuseCountOverCap`'s value dropped at the logistic weights site - the `[[nodiscard]]` build error, and the cap's
    tinytest pin.
34. S3: `usesVarianceForest` answering `usesDart` - test-capi.R's variance-forest arm (critique amendment 1's mutant, moved
    to the field).
35. S3: `usesStudentResidual` reading `carriesShape` - the Student-t arm and the nbinom arm, against the preset
    `results.residualDf`.
36. S3: getInfo writing a field past the caller's `structSize` - the short-struct arm's sentinel.
37. S3: a forest element written for every forest whatever `forestsLength` says - the short-array arm's sentinel (and
    R-loaded ASan).
38. S3: getInfo reading `numSavedSamples` or `kIsSampled` from a copy taken at creation - the re-read arms after a run and
    after setLeafPrior.
39. S3: two info fields swapped in the header - the build's static_assert against the re-baked hash, as mutant 17.
40. S3: a plain `std::logic_error` mapped to FAILED - the re-entering callback's arm, which reads INVALID_ARGUMENT.
41. S2b: the build stage's variance calibration check dropped - a tests/cpp unit of the build stage over the extension's
    variance group at df 2 and two trees; R's pin from S-D holds either way, the parse raising first.
42. S3: the elements indexed by the library's `sizeof(dbarts_forest_info)` instead of the caller's `forestStride` - the
    larger-stride arm (element 1 found at the wrong offset, a tail sentinel overwritten).
43. S3: every element field the library knows written whatever `forestStride` says - the `spec-append` target's 1.0-stride
    unit (the next element's `numTrees` and the closing sentinel overwritten), and R-loaded ASan on the larger-stride arm's
    neighbor.
44. S3: a stride of 0 accepted with an array - the stride-0 refusal arm (every forest written onto element 0, the status OK).
45. S3: the stride checked when no element would be written - the `dbarts_sampler_info_init` arm with no array, refused
    where it must succeed.
46. S3: elements written past `numForests` up to `forestsLength` (a zeroed count in the long array's slots 2 and 3) - the
    long-array arm's sentinels.
47. S2a: one wrapped validator site in a protected parse function returning early on success before its UNPROTECT - rchk on
    the branch tip, and R's "stack imbalance" warning counted in the suite run that reaches it.

## Interaction with the other pending plans, and order

- S1 LANDED on main 2026-11-03 (dec-A260), before coefficient-prior-live, which rebased onto it with no conflict.
- forest-door (slices 1 to 3, with setmodel-parameters' slice 4) LANDED on main; this plan was rebased onto 8fce28568. What it
  moved: `ForestSpec` gained `share`, documented "No caller states it", the group derived by `forestSharers`, so the share left
  the sketched public element for the extension's R-only part; the leaf prior of a several-forest chain is the total the sharers
  divide, so the 1.0 leaf prior fields keep one meaning under several forests (invariant 4 and call 27 rewritten); the bridge
  gained the leaf-prior install on several forests and `setForestTreePrior`, which S2a's validator count and S2b's parse rework
  take as they stand. The math touchpoints, the three `runtime_error` sites and the two process-wide channels are unchanged.
- coefficient-prior-live LANDED on main 2026-11-05 (dec-A261), on this branch's base: the bridge entry `bartcore_setForestCoefficient`, registered in R_interface.cpp, with
  five R-only `Rf_error` lines, one of them `forestSdUnitMessage`, which creation raises too; `Chain::setForestCoefficientPrior`
  and its check `takesForestCoefficientPrior`; the facade and sampler virtual `setForestCoefficientPrior` and the combiner's
  `setForestAmplitudePrior`; in R, setForest's amplitude (R/forestDoor.R), `$setModel`'s coefficient refusal
  (R/fittingLevel.R), and R/model.R's `leafPriorGivesNoSize` and `refuseSizelessCoefficient`, which
  [`resolveSamplerSpec`](../../R/spec.R) now calls at creation with unchanged words. S2a: it adds no shared refusal call (74 and 12 lines, as on main); its entry is neither a creation route nor
  a flat entry and keeps raising, its lines joining the R-only set. S2b: no record slot and no creation read change ("No draw
  moves for any call taken before"), so `ForestRecord` mirrors forest-door's records, and the sizeless refusal stays R's,
  raised before the spec fill. S3: nothing, no flat entry writing a live forest. S1's qt serves its live restatement through
  `forestMapEntry`.
- hurdle-completion LANDED on main (dec-A262; R only): it changes no
  file under src/, inst/include, tests/cpp or standalone/ and nothing any slice here rewrites. Its new tinytest files and pins
  join the full suite S2a to S3 must hold unchanged; inst/NEWS.Rd and TODO are textual neighbors only. Its 0/1 weights reach the
  probit part as R's mask, through no bridge refusal S2a converts.
- latents-first LANDED on main 2026-11-05, slices 1 to 4 (dec-B586, dec-B589, dec-A263), as its plan read. Its overlaps,
  slice by slice:
  - S2a: R_interface_bartcore.cpp's state region. Its 2b retires [`storeState`](../../src/R_interface_bartcore.cpp)'s
    `weights.digest`, `weights.zero` and `survival.digest` writes, the digest encoders, and
    [`setState`](../../src/R_interface_bartcore.cpp)'s reads and reconciliation, about 140 lines. None calls a shared refusal and
    setState is no creation route, so S2a inherits the shorter setState and moves its declaration to the SEXP-facing half as is.
  - S2b: none in code. Creation keeps its deterministic latent starts and draws nothing, and the retired state attributes are
    no creation record the extension holds. It inherits the baselines latents-first records again (one snapshot file, the
    trio's three).
  - S3: facade.hpp and sampler.hpp, from which its 2b removes `weightsDigest`, `survivalDigest`, `reapplyWeights`,
    `zeroWeightRows` and `reapplySurvivalStatus` with their virtuals; S3 adds no virtual, and its predict join,
    `predictPartition` and the stop record sit beside them. [`Chain::runSweeps`](../../src/bartcore/chain.hpp), which its 2a
    reorders, the latent draws first with their cancel point kept: the draw callback still fires at the record, the sweep's
    last step, where S3's stop record and catch go. dbarts.h: its slice 3 writes the getLatents sentence and the run's
    interrupt paragraph (a stop inside the latent draws, which now open a sweep, leaves no tree of that sweep drawn), which S3
    carries into the status form. tests/cpp's test_facade.cpp, test_sampler.cpp and test_state.cpp, where it retires spies and
    twins and S3 adds its units in files of their own. It landed first, so S3 rewrites settled text.
- The fourth review's slice S-D, R surface (dec-A264), LANDED: the bridge's `refuseVarianceCalibration` (a variance forest's
  calibration refused at small df and few trees) and parseModel's infinite-df refusal, both R-worded in the creation parse.
  S2a: no shared refusal; both stay in the SEXP-facing half, raised before the fill, their lines among the 421. S2b: the build
  stage repeats the variance calibration check in neutral words over the extension's variance group (Refusals, mutant 41),
  so a C creation meets it. S3: nothing.
- The fourth review's slice S-F, second R surface (dec-A265), LANDED: R only, nothing under src/, inst/include or tests/cpp.
  Its pins join the suite; its TODO stan4bart-empty-bart-guard is S4's note.
- The fourth review's slice S-A, engine robustness (dec-A266), LANDED; the hash did not move. What it put where this plan
  writes:
  - S2a: the state floors (the engine's `drawnScalarsAreValid`, its k checks among them), which the bridge's setState and
    warm start report in R's words with the warm-start value refusal (`WarmStartResult::valueMismatch`), and the
    sweep-closure entry's `checkSweepBuffer` (an R entry): none a shared refusal, the R words moved to the SEXP-facing half
    as they stand. The generator's state check (src/external/random.c) is compiled as is by the standalone build.
  - S2b: none; creation's parse is untouched.
  - S3: [`Sampler::run`](../../src/bartcore/sampler.hpp)'s re-entry guard and running flag, a `std::logic_error`, which S3
    maps to INVALID_ARGUMENT (Calls made, 43; mutant 40), its stop record sitting beside them; the moves taking back their
    proposal and the rebuild settle paths in chain.hpp and moves.hpp, which make the statuses table's run and
    sampleTreesFromPrior rows state what S-A pins (a stopped run drops the saved trees it wrote; after OUT_OF_MEMORY every
    forest's totals equal its trees and the chain runs on); its dbarts.h text (the stopped run, the state and warm-start value
    refusals, the sweep-closure buffer rules), which S3 carries into the status form; tests/cpp's test_alloc.cpp, whose
    harness S3's OUT_OF_MEMORY arm for the flat run reuses.
- docs-10 LANDED its R share (strings and code comments name the code that runs the BART sampler inside a larger sampler, or
  the caller); its dbarts.h share is S-B's. This plan's new prose follows it. The build-environment names the plan holds
  (`DBARTS_HOST_STANDALONE`, src/host_r.cpp, standalone/host_standalone.cpp) are a different sense, which build, not who
  calls, and stay.
- The fourth review's slice S-C, prior predictive, LANDED on main 2026-11-07 (9a144199b, records 980bc50f9, dec-A269): the
  engine entry `sampleHyperparametersFromPrior` (chain, sampler and facade virtual; `DartPrior` drawing alpha and the split
  weights from its prior), with two `std::invalid_argument` sites; the bridge entry `bartcore_sampleHyperparametersFromPrior`
  (an R entry, registered in R_interface.cpp; two R-only `Rf_error` lines, one forwarding the engine's captured refusal
  after its frame's C++ objects are gone, and a poll after a call that made joint draws); tests/cpp's test_alloc.cpp DART
  joint-draw arm and test_moves.cpp units. S2a: its lines join the R-only set, moved to the SEXP-facing half as is; no
  shared refusal. S2b: none in code; it is no creation route and adds no field, and test-prior-predictive.R's exact-prior
  arms check the prior fields a spec-built sampler carries. S3: no flat entry; its virtual sits in facade.hpp beside the
  predict join, its `invalid_argument` sites are reached from no flat entry, and its poll is one of the five the moved poll
  serves.
- The fourth review's slice S-G, sigest bases, LANDED on main 2026-11-07 (records 61f610bae, dec-A267): R only (R/utility.R,
  R/spec.R), nothing under src/, inst/include or tests/cpp; the equivalence and bcf-equivalence baselines re-recorded and
  renamed after the draw-changing commit as rebased, which the MANIFEST names. S2b takes the starting sigma R resolves, now
  regressed on each forest's basis, as the spec's `sigmaEstimate`, unchanged; every slice's gates run against those
  baselines or later ones. S2a and S3: nothing.
- The fourth review's slice S-B, draw channels, LANDED on main 2026-11-07 (fd78973ee, records cdb3fd26f, dec-A270); not an
  ABI event, dbarts.h changing in comments only, the hash unmoved. What it built: `Results` carries a per-draw stride for
  every draw-indexed channel (`logLikelihoodStride` among five new); `Chain::carriedChannels` states the channels the fit
  carries; `Sampler::provideObserverScratch`, called by `Sampler::run` only when a callback is set, points each carried
  channel the caller left NULL at one draw of per-chain scratch owned by the run's frame (dec-B600); `Sampler::kIsSampled`
  reads chain 0's forest calibrations, so it follows setLeafPrior and setModel, and `bartcore_run` and the flat run read it
  (dec-A270's bridge-03); a draw's sigma and `results.sigma` are NaN under a variance forest (dec-B601); the uniform grid
  stays finite over a span wider than `DBL_MAX`; `cutGridIsValid` requires finite ends and `bartcore_setCutPoints` refuses a
  non-finite cut in R's words (dec-B604); dbarts.h's PRESENCE AND ABSENCE paragraph, its results paragraph's k and
  log-likelihood sentences, `dbarts_draw::sigma`'s and `dbarts_sampler_kIsSampled`'s comments, and "caller" for every
  "host"; consumer.c's recording callback and test-capi.R's channel arms over four routes. S-B is done: nothing of it is
  left to land. What it hands this plan, all S3's but one:
  - the facts: the header's sigma sentence names no accessor, and S3 adds that `usesVarianceForest` predicts the NaN; S3's
    getInfo reads `kIsSampled` from the current computation; the one-fact `dbarts_sampler_kIsSampled` and its new comment
    go with the other ten, the comment's sentence moving to the info block and `dbarts_results`' "re-reads
    dbarts_sampler_kIsSampled" becoming a re-read of the struct; TODO capi-sampler-info is this plan's, closed at S5;
  - the flat run: S3 rewrites `dbarts_sampler_run` with the scratch in the engine, so the flat file holds none, and keeps its
    `results.k` rule under the status return; the stop record replaces the `stoppedByCallback` out-parameter and the relay's
    flag beside S-A's guard and the scratch in `Sampler::run`; the callback catch on the worker trampoline sees the scratch
    pointers S-B hands it, which tests/cpp's `testDrawCallbackCarriesEveryChannel` (each chain's blocks distinct and stable)
    holds under the catch;
  - consumer.c (1,369 lines, 50 one-fact calls) and test-capi.R (2,258 lines), which S3 rewrites for statuses with S-B's
    recording arms kept;
  - the caller's name: no "host" left in dbarts.h, so the release map's host-wording item folded into S3 is done, and S5 closes
    TODO host-wording-sweep if it is still open;
  - S2a (the one not S3's): the setCutPoints refusal is R-worded and no shared refusal, moved to the SEXP-facing half as is;
    the state's and donor's infinite-grid refusals reuse existing messages through the engine's floor.
- frequency-weights (dec-B563 to dec-B570, dec-B580; S1 to S3 built and S4 building, read at 831e0f434). On its tip the
  shared refusals are twenty, on 90 bridge lines and 13 flat (`refuseWeightedCountsOverCap`, `refuseWeightedTrialsOverCap`
  and `refuseBinomialResponseSwap` new, the last called from C_interface.cpp's setResponse), and the R-only `Rf_error` lines
  428; `ext_lchoose` and choose.c are built. Its call 37 refuses a binomial response swap as R's `$setResponse` does, so S3's
  support row gives UNSUPPORTED for it and INVALID_ARGUMENT for a value out of support; its call 36's setData refusal is
  R-only. Its overlaps:
  - S2a: the shared refusals it lifts or widens, about 200 bridge lines in all: `enforceBinaryWeightPolicy` (its S4: 0/1 on
    probit and ordinal, a logistic 0), `refuseCountOverCap` (its S2: the product cap at w > 1, its Refusal 8),
    `validateResponseSupport` (its S3: successes up to the trials), `refuseBinaryWeightChange` as its S4 leaves it, and the
    R-worded creation refusals of nbinom, ordinal and multinomial weights, which its S2 and S4 remove.
  - S2b: its S1's reading of the resolved `weights.type` model slot becomes the spec's `weightsType`; its S3's counts slot read
    under logistic becomes `y` and `trials`; the build's weights check is its resolution table; the field coverage suite takes
    both.
  - S3: the header's `weightsType` and `trials` and the getLatents frequency sentence; consumer.c's binomial and frequency arms.
    Its S3 has the flat setResponse refuse a binomial sampler's response swap through `refuseBinomialResponseSwap`, as R's
    `$setResponse` refuses it (its Refusal 15, call 37), so S3 gives that UNSUPPORTED and a value out of support
    INVALID_ARGUMENT.
  - S1's copy: its S3 added `ext_lchoose` and choose.c by the recipe in S1's section; S3's standalone consumer then links seven
    names.
  - No overlap: its engine work in model.hpp and the multinomial loop in combiner.hpp; predict's weights and predict's doors
    (R only, and the flat predict takes no weights); its reworded two-column refusals (R's data and xbart code); its R slices S5
    to S8.
  - Since 831e0f434 (read again at ab739157c, S1 to S7 built, S8's records written; S7 in its fix round, then the x86 leg,
    the SBC dispatch and its landing): no dbarts.h change and no new shared-refusal declaration; its bridge gains a binomial
    sampler's setData with counts (dec-B628) and loses the creation-time switch for nbinom's weights reading (dec-B632);
    S7 is gate-only (frequency-exact.R, the SBC arm), S8 records, and dec-B649's fit-time refusal of a factor count column
    is R's. By the wider grep its tip has 108 bridge and 21 flat shared-refusal lines and eight more R-only `Rf_error`
    lines than its base. S2a converts the refusals as it leaves them; S2b reads dec-B632 as the build's NBINOM rule; S3
    reads a binomial setData's row count in the facts re-read arm.
- rchk-findings (TODO rchk-findings; branch at 843470f35, two commits on c4ba5427a, reviewed LAND, CI and rchk running,
  landing next): the bridge's model parser protects its slot reads (`slotOrNull` results PROTECTed across allocating calls in
  parseModel, parseDartSettings, the move-record readers, the variance and forest record readers, `refuseVarianceCalibration`
  and `forestMessageName`, about 120 lines of R_interface_bartcore.cpp), and a debug-only count in chain.hpp moves into its
  check block. No draw, message or header change. S2a: those functions sit in the creation parse that stays R-side; where
  S2a wraps a validator call inside one, the protection stays balanced (mutant 47), and rchk joins every slice's gates.
  S2b: its parse rework reads exactly these functions and moves their outputs into spec fields, keeping each PROTECT pair;
  about half an hour more than on the unprotected form, counted. S3: chain.hpp's moved count is beside no S3 code.
  Rebase cost: a docs rebase now, the branch holding the plan only.
- state-basis-levels (dec-B647, revising dec-A210's call 6; plan at 5b068cceb, building): a stored state records each factor
  forest's level names and `$setState`, `copy()` and a reload put each coefficient on its own level by name. Its plan's
  call 1 puts the record in an R attribute, R matching the names and handing `bartcore_setState` one integer order per
  forest, the bridge reordering each chain's glue block; its alternative, the engine's state, is not taken there (the engine
  keeps no names). S2a: `setState` and `readAmplitudeGlue` move to the SEXP-facing half with the added argument as they
  stand; no shared refusal. S2b: none; a state is no creation input. S3: none; the flat API carries no state. Were the record
  put in the engine's state instead, it would still meet no slice here but S2a's move of `storeState`, at no added cost.
  Rebase cost: a docs rebase.
- small-rulings-1107 (dec-B642, dec-B643, dec-B646; branch at 2e170e9f1, in review): R only, nothing under src/,
  inst/include or tests/cpp. Its pins join the suite every slice must hold. Rebase cost: a docs rebase.
- The dec-B640 slice (TODO multinomial-setdata; no plan or branch yet; after frequency-weights lands): a multinomial sampler
  takes `$setData`, the new n x K counts giving the response and the number of rows free; about 200 engine lines (the
  whole-data replacement, `Chain::recoverTreeParameters` and `Chain::applyNewData`, over the K category forests and their
  saved draws, and `MultinomialForestCombiner`'s per-row state resized and re-pointed), 40 bridge, 20 R and 150 test lines.
  S2a: it changes where `refuseMultiForestMutation` (a shared refusal) refuses and replaces the holder's owned counts and
  trials, both of which S2a converts and moves; it should land before S2a starts, about half an hour of S2a's recount and
  conversion if so, a between-slice rebase and re-gate (about an hour) if it lands during S2a to S3. S2b: the multinomial
  route's `MultinomialRecord` borrows the counts from the holder, which setData replaces after creation, and a
  spec-created multinomial sampler's setData joins S2b's gates. S3: no flat entry; the facts re-read arm reads its row
  count. Recommended: land it after frequency-weights and before approximate-draws' S2, or at the latest before S2a.
- approximate-draws (dec-B583 to dec-B585, dec-B593; S1 built and accepted at e50e1145a, its tests/cpp density units only;
  S2 waits for frequency weights and an independent proof of its call 28's bound): dec-B593 rules its setting a method
  name and a separate tolerance. Its S2 adds their engine side, read from the control in the parse S2b moves; it lands
  first, so S2b carries both into the extension, with no 1.0 field (Calls made, 35). Its S2's close-call reader becomes a
  field S3 appends to the facts struct, never an entry (Calls made, 44). Its S3's real weights on logistic and the binomial
  widen what the build takes under FREQUENCY, with no field. Its S1 (the density, src/bartcore/pgdensity.hpp and
  ddouble.hpp, header-only and R-free) and S4 (the series) touch nothing here; the standalone build compiles them as they
  are.
- Landings between this plan's slices: the DART split rule (dec-B605, TODO dart-split-rule) and k-move-pricing if its move
  wins (dec-B530) are engine work that may merge while S2a to S3 build; each costs this plan one rebase and re-gate between
  slices, below.
- multinomial-doors' C creation path stays post-release (dec-B160, dec-B575); this plan ships the shape it will append to and
  the R route it will expose.
- Frequency weights' dec-B580 (a 0 row's predictive draws NA on every family) asks nothing of the spec: the row is drawn and its
  draws set to NA after, so the stream holds, and the spec states only the weights and their reading.
- The sister re-verification at the freeze follows this arc (TODO), with stan4bart's sanitizer job dispatched against the final
  header (dec-A36).

Order, the release map's critical path with the work since placed: S-A, S-C, S-G and S-B (all landed by 2026-11-07); then
rchk-findings (landing next); state-basis-levels and small-rulings-1107 whenever each is ready (R and bridge-state only,
neither on the critical path); frequency weights S1 to S8 (S7's fix round, then its landing); the dec-B640 slice (after
frequency weights, ruled "its own slice after frequency-weights lands"; placed here, before approximate-draws' S2, so it
lands before S2a); approximate-draws S2 to S5 (its S1 accepted; S2 after frequency weights and the proof of its call 28);
then S2a, S2b and S3, with S4 right after the merge; then the approximate draws' defaults study, the freeze (the hash
re-baked, the sisters re-verified, stan4bart's sanitizers) and the release battery. That is the 2026-11-05 amendment's
recommended order (frequency weights before S2a); the other order's cost, 6 to 10 hours of re-targeting in frequency
weights and a second header change, is no longer in play. The header and its hash change once, with the final weights
reading, the approximate draws' settings in the extension, S-B's channels and the facts struct with its forest element.

What each rebase costs this plan:

- Before S2a, onto rchk-findings, state-basis-levels, small-rulings-1107, frequency weights, the dec-B640 slice and
  approximate-draws as each lands: the branch holds this plan only, so each is a docs rebase, minutes, any conflict in
  docs/plans/INDEX.md resolved by keeping both rows. No re-gate: nothing is built. The rebase onto 7e7edb526 (S-C, S-G and
  S-B landed) was such a rebase, with no conflict.
- At S2a's start, one recount on the base it then has: shared refusal lines, the R-only `Rf_error` lines, the
  `invalid_argument` sites, and the bridge regions frequency weights, rchk-findings, state-basis-levels and the dec-B640
  slice moved, about half an hour, counted in S2a.
- The features' own code, counted in the slices: frequency weights' twenty refusals and ~200 bridge lines, about an hour
  each in S2a and S2b; approximate-draws' two settings, about half an hour in S2b; rchk-findings' protected parse, about
  half an hour in S2b and a quarter in S2a's gates; the dec-B640 slice's refusal and holder counts, about half an hour in
  S2a and half in S2b; S-B's run scratch and kIsSampled, landed and read, about half an hour in S3.
- If the dec-B640 slice or state-basis-levels slips into S2a to S3's build, it becomes a between-slice landing: the rebase
  and re-gate below, about an hour, and for dec-B640 the S2a conversion of its refusal done on the rebased tree.
- Between slices, each landing that touches src/ (the DART split rule, k-move-pricing, any fix round): one rebase and the
  re-gate below, about half an hour to an hour each; the ~5 hours of rebases count four such and the docs rebases.



The plan said "Nothing else touching src/ merges while S2a, S2b or S3 is in review". Restated so other features can land
between slices:

- Other features merge to main whenever they are ready, engine features serially as before; this branch, merged once,
  absorbs them.
- The branch rebases onto main only between slices: never while a slice is in review (from the implementer's report through
  the reviewer's verdict on the fix round), and always before the next slice's implementer starts.
- After each rebase the slices already built are re-gated on the rebased tree, against the baselines the MANIFEST then names:
  the trio bitwise, the four snapshot files, tests/cpp and the full suite in one process; the R-loaded ASan run and rchk when
  the rebase touched the bridge; and, once S3 is built, the hash re-bake, header-compiles and test-capi.R when it touched dbarts.h, with
  S4's consumer runs repeated before the merge.
- A feature that rewrites what a slice in build rewrites (the creation parse, the shared refusals, the flat entry file,
  dbarts.h) merges before that slice starts or after its review ends; if it merges during the build, the slice rebases before
  its report.

This works because every gate here is neutral: a rebase is checked against the base it lands on, not a fixed one.

## Cost

S1, landed: about 910 lines written, 6,329 copied and 395 of reduced headers, none of them in the tarball. What is left: about
6,120 lines written here and ~165 in the two consumers; ~103.5 agent-hours estimated (S2a ~22.5, S2b ~33, S3 ~35, S4 ~7,
S5 ~1, rebases ~5), 155 to 205 by the record of plans here, S2b and S3 most of it, in the release map's order. The rise from
~99.5 (the fifth amendment): dec-B639's per-forest element (about 150 lines and 1.5 hours in S3, 5 lines in S4), S3's
consumer.c carrying S-B's six more accessor calls and its recording arms (a quarter hour), rchk-findings' protected parse
(about 40 lines and three quarters of an hour across S2a's gates and S2b's rework, rchk added to every slice's gates), the
dec-B640 slice's shared refusal and holder counts (about 20 lines and an hour across S2a and S2b) and S2b's multinomial
setData twin (a quarter hour). S-C, S-G and S-B landing as read cost nothing more. The forward-compatibility work is about
300 of those lines and 4 of those hours; the callback statuses and the status words about 300 and 3.

| slice | state | lines written | agent-hours (was) |
|---|---|---|---|
| S1 | LANDED 2026-11-03 (dec-A260) | ~910 (+6,329 copied, 395 reduced) | builder ~3 (8 planned in all) |
| S2a | planned, after rchk-findings, frequency weights, the dec-B640 slice and approximate-draws | ~1,000 | ~22.5 (21.5) |
| S2b | planned | ~1,790 | ~33 (32) |
| S3 | planned, with the facts struct and its forest element | ~3,330 | ~35 (33) |
| S4 | planned, in stan4bart and treatSens | ~165 | ~7 (7) |
| S5 | planned, records | - | ~1 (1) |
| rebases | docs rebases before S2a, the recount, and between slices | - | ~5 (5) |

## Not verified

- The copy's bitwise match is measured on two hosts (arm64 macOS with R 4.6.2 and the x86 Linux box with R 4.6.1, dec-A260);
  other compilers and C libraries are gated only at the floor, and on a C library with C23's cospi R's qt can differ from the
  copy's at df within 1e-12 of 1 for p other than a quartile (prefix.h's note), which the engine's qt at 0.75 does not reach.
- What a failed run (OUT_OF_MEMORY, FAILED) and a failed sampleTreesFromPrior leave in the chains beyond "valid to read and
  destroy"; S3's tests pin what the table states, and the header promises no more than they show.
- That no other process-wide mutable state is on the predict path beyond the two channels named; a ThreadSanitizer arm in S3 is
  the check.
- Every bridge record S2b's extension must mirror: coefficient-prior-live adds no slot (its diff against main read), and
  forest-door's were read for `ForestSpec` and `forestSharers` only.
- Where the engine's sweep loop can carry a chain index and a callback value to the stop record without a cost per sweep; the
  ruling's "about 40 lines" is its estimate, not measured, read against latents-first's loop as landed, with S-A's guard and
  S-B's scratch (landed) beside it.
- The binomial response's engine form is frequency weights' plan (`trials` beside `y` under logistic); not built, and its flat
  setResponse on a binomial sampler is that plan's to settle.
- frequency weights' bridge changes (about 200 lines) are its plan's estimate, counted on its tip with S4 still building;
  S2a's implementer recounts the shared refusals' lines on the base it starts from.
- S-B, S-C and S-G are landed and their code read on main 7e7edb526; what S3 does with S-B's scratch under the worker
  catch is planned, not prototyped.
- The dec-B640 slice has no plan: its overlap with S2a (`refuseMultiForestMutation`, the holder's owned counts) and S2b
  (`MultinomialRecord` borrowing them) is read from the ruling's sizing, not from its code, and its own estimate goes back
  to the maintainer if it sizes well above it.
- state-basis-levels is building: where its builder puts the level record (its plan: an R attribute) was read from the
  plan; either place meets no slice here beyond S2a's move of the state code.
- rchk-findings' protection was read on its branch before landing; rchk's own CI result on it was not read here.
- The forest element's stride writes by `memcpy`, a stride that is not a multiple of the element's alignment, and a
  32-bit platform's `sizeof(size_t)` of 4 (the first release's element size there) are designed, not probed.
- approximate-draws' setting names are its S5's to settle and the close-call reader's shape its S2's (dec-B593 rules only a
  method name and a separate tolerance); Calls made 35 and 44 hold for any names whose zero means the default.
- Which facts setters change beyond those the facts section names (kIsSampled through setLeafPrior and setModel, the saved
  sample count after a run): whether setModel can change `usesDart`, or a setter a tree count, was not traced through every
  setter (setForest and setVarianceForest refuse `n.trees` by name, read in their help); S3's re-read arms pin the ones
  named, and getInfo reads current values, each element's `numTrees` included, whichever change.
- S-A's settle is now pinned by its tests for the in-sweep allocation failure; what a FAILED run leaves beyond that is still
  "valid to read and destroy", as above.
- Whether the engine's several-forest chain takes DART, split probabilities and the k prior per forest; the sketch leaves them
  sampler-wide and the append refuses what the factory refuses.
- stan4bart's and treatSens' call sites were read on their main (98db852) and master (268b2db); their 1.0 compat branches
  may differ.
- That no consumer under R other than the two links dbarts.h; the reverse LinkingTo list on CRAN was not queried.

## Amendments of 2026-11-05

Brought up to date for S2a onward, on main 7937d252f; nothing was built or run beyond reads, recounts and the doc check.

| # | change | where it went |
|---|---|---|
| 1 | S1 landed on main (dec-A260): its as-built differences (four reduced headers, the 36,093-point grid, misc's relocatable Makefile, `r-math-symbols` and `standalone-config`, the configure fix), its gates, and what S2a starts from | Status; front block; Summary; Context (touchpoints, build, tests/cpp); Build-time math, As built; RNG class; Slices, S1; Cost; Not verified |
| 2 | Mutant 1's premise corrected: the engine's pnorm feeds only probit's and aft's log-likelihood channel, not their draws; S1's four mutants recorded as run with their killers | Mutants, 1 to 4 |
| 3 | The seventh math name for frequency weights' `lchoose`: one file, choose.c, with S1's review's recipe and grid | Build-time math, A seventh name; Interaction, frequency-weights |
| 4 | coefficient-prior-live (dec-A261): its bridge entry, chain and facade methods and R door placed; no shared refusal call and no record slot added; the forest-records gap closed | Context; Interaction; S2a; S2b; Calls made, 1; Not verified |
| 5 | `refuseCountOverCap` counted as a seventeenth shared refusal, returning a value in S2a; the shared refusals recounted at 74 bridge lines and 12 flat (the old "about 60" low), R-only `Rf_error` lines at 416, `invalid_argument` sites at eleven | Context; S2a; Calls made, 36; Mutants, 33 |
| 6 | hurdle-completion (dec-A262): R only, touches nothing here; its pins join the suite | Interaction; S2a gates |
| 7 | latents-first: the state region S2a inherits trimmed, nothing for S2b but the re-recorded baselines, and for S3 the sweep loop, facade and sampler removals, tests/cpp files, and the header's getLatents and interrupt sentences carried into the status form | Interaction; S2a; S2b; S3; The header diff, prose; front block, rng; Not verified |
| 8 | frequency-weights: `trials` joins the spec (52 fields) read under LOGISTIC only and refused elsewhere; the weights check is its resolution table, the interim FREQUENCY clause only for the other order; its lifted refusals converted as it leaves them; consumer.c's binomial and frequency arms; the flat setResponse on a binomial sampler mapped to the support row | The specification; Forward compatibility, invariant 6 and the 1.0 refusals; Refusals; S2a; S2b; S3; Calls made, 37; Mutants, 31 and 32; Interaction |
| 9 | approximate-draws: its latent sampler setting gets no 1.0 field, its zero the default rule, carried in the extension until a binding asks; its real weights widen the build's check, not the spec | The specification; Calls made, 35; S2b; Interaction; Not verified |
| 10 | Order: frequency weights before S2a to S3 recommended, with the costs of each order and a stacked variant; "nothing else touching src/ merges while a slice is in review" restated as rebases between slices with re-gates | Interaction, Order and the restated constraint; front block, window; Summary |
| 11 | Re-estimated: ~5,550 lines and ~91 agent-hours left (135 to 180), was ~93 in all with S1 | front block, budget; Slices; Cost; Summary |
| 12 | The retired `misc_printf` cite (S1 deleted its header) dropped | Context |

## Amendments of 2026-11-07

Rebased onto main f79874ea8 and brought up to date for dec-B615 and the work landed or in flight since 7937d252f; nothing was
built or run beyond reads, recounts and the doc check.

| # | change | where it went |
|---|---|---|
| 1 | dec-B615: `dbarts_sampler_info`, led by its size, with `DBARTS_SAMPLER_INFO_INIT`, `dbarts_sampler_info_init` and `dbarts_sampler_getInfo` returning a status; its fields and their order, `usesVarianceForest` (dec-B601) and `usesStudentResidual` as fields of their own, growth only at the end, current values on every read | front block; Summary; Goal; Statuses; The sampler's facts; The header diff; Ruled, 1 and 8; Calls made, 38 to 42 |
| 2 | The eleven one-fact entries removed in S3, no compatibility owed (0.9-34 had no such header); the hash re-baked with the info fold; stale stan4bart binaries fail at symbol lookup | The sampler's facts; The header diff, list, hash and development binaries; S3 |
| 3 | Open call 8: how the struct carries each forest's tree count, (a) recommended, a caller array and its length | Summary; Open calls, 8; The sampler's facts |
| 4 | The facts' tests (every field, tree counts, re-reads, size refusals, the short-struct sentinel) and mutants 34 to 39 | S3; Mutants, 34 to 39 |
| 5 | Consumers: stan4bart's eleven calls become one getInfo per function, treatSens unchanged in this | Context; S4 |
| 6 | S-B's handover: no facts entry, the sigma sentence naming no accessor, kIsSampled's current value, the run's callback scratch and stride fields, the cut-grid refusals R-only | Interaction, S-B; The sigma getter; S3 |
| 7 | S-A (dec-A266): run re-entry's `std::logic_error` mapped to INVALID_ARGUMENT, the stopped and failed run rows as S-A pins them, the settled tree in sampleTreesFromPrior, its state and sweep-closure refusals moved as is | Context, exception types; Statuses; The draw callback; Calls made, 43; Mutants, 40; Interaction; S2a; S3 |
| 8 | S-D (dec-A264): the variance calibration refusal repeated in the build stage's own words; S-F (dec-A265) and docs-10 R only; the plan's new prose names the caller | Refusals; S2b; Mutants, 41; Interaction |
| 9 | coefficient-prior-live, hurdle-completion and latents-first recorded as landed | Interaction; Summary |
| 10 | S-C, S-G and frequency weights' built slices placed; shared refusals recounted (17 on main, 20 on frequency weights' tip), R-only lines 421 on main; the binomial swap's UNSUPPORTED row settled by its call 37 | Context; Interaction; S2a; S3; Calls made, 36; Not verified |
| 11 | approximate-draws under dec-B593: method and tolerance carried in S2b's extension, no 1.0 field; its close-call count a later facts field | Calls made, 35 and 44; S2b; Interaction; Not verified |
| 12 | Order restated as the release map's, the 2026-11-05 options retired; each rebase's cost named | Interaction, Order; front block, window; Summary |
| 13 | Re-estimated: ~5,920 lines and ~99.5 agent-hours left (150 to 200), was ~91 | front block, budget; Slices; Cost; Summary |
| 14 | Not verified updated: S-B unbuilt, setters that change a fact, approximate-draws' names, S-A's settle | Not verified |
| 15 | S5 closes TODO capi-sampler-info | S5 |

## Amendments of 2026-11-07, second

Rebased onto main 7e7edb526 (from f79874ea8, no conflict) and brought up to date for dec-B639, the landings of S-C, S-G and
S-B, and the work in flight; nothing was built or run beyond reads, recounts by grep and the doc check.

| # | change | where it went |
|---|---|---|
| 1 | dec-B639 rules open call 8 as (d), not the (a) recommended: Open call 8 moves to Ruled with the maintainer's words; no call is open | Ruled, 9; Open calls; Status; Summary |
| 2 | The per-forest element `dbarts_forest_info` (one field, `numTrees`, `dbarts_forest_spec`'s first), `DBARTS_FOREST_INFO_INIT`, and the parent's `forestsLength`, `forestStride` and `forests` in place of `treeCounts` and `treeCountsLength`; the stride rule (checks before writes, `memcpy` per field within the stride, nothing past it, elements past `numForests` left), the stride refusals, the INIT setting the stride and the function not | The sampler's facts; Statuses, getInfo row; Refusals; The header diff, types and prose; S3 |
| 3 | The element designed beside the sketched creation element: shared names, types and meanings, appends mirroring it, zero meaning absent per element | Forward compatibility, the reported side and the element paragraph; The sampler's facts |
| 4 | The hash folds the element by field count and offsets, the INIT macros unfolded; the token re-baked once as before | The sampler's facts, hash; The header diff, hash |
| 5 | Calls made 45 to 49: the names, the caller fields' order, the INIT and the init function on the stride, no element init function, `memcpy` writes and the stride checked only where an element is written | Calls made |
| 6 | Tests: a short array, a long array, a larger stride, an older element size through the `spec-append` target's test field, NULL and zero-length arrays, the stride refusals, the init function's zeroed stride | S3 tests; Forward compatibility, the gate |
| 7 | Mutants 37 rewritten and 42 to 46 added for the element; 47 for rchk's balance under S2a | Mutants |
| 8 | stan4bart's port: printTrees' `dbarts_sampler_numTrees(fit, 0)` becomes an array of one element read through getInfo; its call sites re-read at 98db852 | S4; Context, consumers; Not verified |
| 9 | S-C (dec-A269) recorded as landed: engine and facade virtual, the bridge entry's two R-only lines and its poll, the poll's five sites, thirteen `invalid_argument` sites | Context; Interaction, S-C; S2a; S2b; S3 |
| 10 | S-G (dec-A267) recorded as landed: R only, the re-recorded baselines | Interaction, S-G; front block, rng; S2b |
| 11 | S-B (dec-A270) recorded as landed and done: what it built, what it hands S3 (the facts, the flat run's scratch and k rule, the stop record's place, consumer.c's 50 calls and recording arms, no "host" left) and S2a (the setCutPoints refusal) | Interaction, S-B; The sigma getter; The sampler's facts; The header diff, prose; S2a; S3; S5 |
| 12 | Recount at 7e7edb526 by a wider grep: shared-refusal lines 90 and 20, unmoved since f79874ea8, 108 and 21 on frequency weights' tip; three R-only `Rf_error` lines added | Context; S2a |
| 13 | rchk-findings placed: its protected parse functions inherited by S2a and S2b, rchk added to every slice's gates and to the rebase re-gate | Interaction; S2a; S2b; Tests and gates; Order |
| 14 | frequency-weights since 831e0f434: dec-B628's binomial setData, dec-B632's nbinom counts reading under every creation (the build's NBINOM rule, no creation record), S7 gate-only and S8 records | Interaction, frequency-weights; S2a; S2b; S3 |
| 15 | state-basis-levels (dec-B647), small-rulings-1107 and the dec-B640 slice placed, each slice's overlap and each rebase's cost named; dec-B640 recommended before S2a | Interaction; S2a; S2b; S3; Order |
| 16 | approximate-draws: S1 accepted, header-only and R-free, S2 waiting on frequency weights and its call 28's proof | Interaction, approximate-draws |
| 17 | Order restated: rchk-findings, then frequency weights, the dec-B640 slice, approximate-draws, then S2a to S3, the R-only work any time | Interaction, Order; front block, window; Summary |
| 18 | Re-estimated: ~6,120 lines and ~103.5 agent-hours left (155 to 205), was ~99.5 | front block, budget; Slices; Cost; Summary |
| 19 | S5 also closes TODO host-wording-sweep if still open | S5 |
| 20 | Not verified updated: S-B, S-C and S-G landed, the dec-B640 slice unplanned, state-basis-levels' record, rchk's CI, the element's stride writes | Not verified |
