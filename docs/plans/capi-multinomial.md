# capi-multinomial

Status: LANDED 2026-10-01 (5233f0ae, 58bface9, ae1f46cd; landing note at
EOF); questions ruled (dec-B160, dec-B163 to dec-B169 in
[decisions.md](../decisions.md)).

agent: opus (commits 1 and 2: the header, the C entry file, the bridge
  poll and the test consumer; step 4's stan4bart review); sonnet (commit 3,
  docs; step 4's consumer runs; commit 5, records). Serial: 1 to 5.
rng: NEUTRAL. No engine file changes and no draw moves. What changes is the
  shipped contract: two accessors, one renamed callback field, the varcount
  width the flat caller declares, the predict offset rule, and the flat
  run's interrupt poll and slow-count warning. A gaussian or probit run with
  no interrupt and no slow count is byte for byte what it is today.
window: nothing else edits `inst/include/dbarts/dbarts.h` while this is
  open. Lands before the 1.0-0 tag, or bumps the minor version (Constraints).
budget: header ~+110 -30; C entry file ~+110 -15; bridge and common header
  ~+50 -25; test consumer ~+280 -10; test files ~+340; tests/cpp, vignette
  and test helper renames ~+15 -15; docs ~+60 -25; records ~+50. About 1050
  lines; plan on 1600-2100.

## Goal

A multinomial sampler built in R runs and predicts through the flat C API,
and a caller who sizes its buffers by the documented layouts gets bitwise
what R's `$run` and `$predict` return on every channel the flat struct
carries, on every model, the split counts of every multi-forest model
included.
- The header gains `dbarts_sampler_numFittedValuesPerObservation` (K on
  multinomial, 1 elsewhere) and `dbarts_sampler_numVariableCountForests` (K on
  multinomial, the mean-forest count on an amplitude-coupled model - 2 on
  BCF, 3 on a three-forest `forests = list(...)` sampler - and 1 elsewhere).
- The [`dbarts_draw`](../../inst/include/dbarts/dbarts.h) field `numReportedLocations` is renamed
  `numFittedValuesPerObservation`.
- Flat predict reads a multinomial offset as a rows x K matrix where R does.
- [`dbarts_sampler_run`](../../inst/include/dbarts/dbarts.h) honours R's interrupt poll and raises
  R's slow-count warning once per sampler, closing TODO
  monotone-count-host-interrupt.

Out of scope: a flat creation path, and any flat count or category-offset
setter. Both stay with the R object (TODO multinomial-doors).

## Context

The defect.
- Every handle is an R sampler's pointer
  ([Landing note, S2 (2026-09-08)](pure-c-header.md#landing-note-s2-2026-09-08) removed flat creation),
  and `dbarts()` builds multinomial samplers.
- The header treats them as in contract:
  [`dbarts_sampler_family`](../../inst/include/dbarts/dbarts.h) reports `DBARTS_FAMILY_MULTINOMIAL`,
  [`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) documents its multinomial answer,
  and the [`dbarts_results`](../../inst/include/dbarts/dbarts.h) comment describes multinomial
  log-likelihood. But its buffer sizes assume one fitted value per
  observation.
- Run overflows. [`dbarts_sampler_run`](../../src/C_interface.cpp) passes the sampler's
  location count to the engine, and [`Chain::storeSample`](../../src/bartcore/chain.hpp) writes n x K
  per draw whatever the caller sized. Train and test therefore come back
  n x K x S x C into buffers documented as n x S x C, a heap overflow.
  The verification probe (n = 30, K = 3, S = 4, C = 2) wrote 720 of a
  documented 240 train entries and crashed later in unrelated R code.
- Predict overflows the same way. [`dbarts_sampler_predict`](../../src/C_interface.cpp) writes
  nTest x K x S x C. It adds a non-null offset slab by slab at stride nTest,
  an additive shift on probability channels at the wrong entries. It also
  ignores a category offset the sampler carries from R and returns the
  offset-free surface, where R refuses.
- varcount is safe today. The run leaves `numVariableCountForests` at 1, so
  the engine writes one slab: BCF's prognostic forest, or a multinomial's
  first category.
- No current entry reports K before a run.
  [`dbarts_sampler_numTrees`](../../inst/include/dbarts/dbarts.h) raises past the last forest, and
  the callback field arrives only inside a run.

What R returns (probed against the installed build: n = 30, K = 3, nTest =
10, four draws, two chains). [`bartcore_run`](../../src/R_interface_bartcore.cpp) gives:

| channel | R shape | content |
|---|---|---|
| sigma | S x C | the pinned 1 |
| train | n x K x S x C | softmax probabilities; each row of K sums to 1 |
| test | nTest x K x S x C | the same, under any category test offset set from R |
| varcount | p x K x S x C | slab k is category k's forest (BCF: p x 2 x S x C) |
| k, varprobs | NULL | a k hyperprior and DART are refused on multinomial |

Predict in R:
- [`predictFromSource`](../../src/R_interface_bartcore.cpp) returns nTest x K x S x C with saved
  trees and nTest x K x 1 x C without.
- Its offset must be a finite nTest x K per-category matrix entering
  before the softmax.
- A sampler carrying a train or test category offset set from R (installed
  by [`bartcore_setCategoryOffset`](../../src/R_interface_bartcore.cpp) or
  [`bartcore_setCategoryTestOffset`](../../src/R_interface_bartcore.cpp)) refuses a predict that names
  no offset, in [`bartcore_predict`](../../src/R_interface_bartcore.cpp).
- Category k is column k of the sampler's count matrix.

The other entries on a multinomial handle already behave.
- setResponse and setOffset return 0
  ([`responseConduitIsFixed`](../../src/R_interface_bartcore.cpp)), as do setSigma and getLatents.
- printTrees and numTrees take forests 0..K-1.
- The callback sees K fitted values per observation and, on the flat
  route, `numVariableCountForests` = 1.

The monotone gaps (TODO monotone-count-host-interrupt):
- The flat run passes the engine an empty interrupt poll and reads no
  slow-count tally.
- R's run passes `bartcore_userInterrupted`, a static in
  [R_interface_bartcore.cpp](../../src/R_interface_bartcore.cpp). That function runs
  `R_CheckUserInterrupt` under `R_ToplevelExec`, so the sampler joins its
  workers before the interrupt becomes the error "sampler run interrupted".
  The engine relays the poll into a leaf-order count through its count
  cancel function, which tests/cpp exercises in `testMonotoneCountInterrupt`.
- After every run, R's run attaches the tally
  ([`attachSlowCountTally`](../../src/R_interface_bartcore.cpp)), and the R function
  [`warnOnSlowCount`](../../R/bartcore.R) raises `dbartsSlowCountWarning`.
- The slow threshold has a test hook,
  [`bartcore_setMonotoneCountHooks`](../../src/R_interface_bartcore.cpp), driven from
  ["countHooks"](../../inst/tinytest/test-monotone.R). No R-level test of an interrupt exists on
  either route.

Consumers today.
- stan4bart (branch bartcore at 9a3be93) calls run, predict, setOffset,
  setSigma, getLatents, setTreeStorage, sampleTreesFromPrior, printTrees,
  destroy, the size queries and the version pair. It builds only gaussian and
  probit samplers.
  - Its `run` .Call (src/init.cpp) loops over iterations calling
    `dbarts_sampler_run(..., 0, 1, ...)` once per sweep. It holds a raw
    `new`ed `IterableBartResults` and PROTECTed callback vectors across the
    loop.
  - Its init path holds a `std::unique_ptr<Sampler>` and a `std::vector`
    across its first run.
  - It never polls for interrupts itself.
- treatSens (dbarts-1.0 at 7cc6a0f) is also a C consumer. It calls run,
  setResponse, setOffset, setSigma, setNumThreads, setVerbose and destroy,
  on gaussian and probit samplers only.
- Both use `DBARTS_USE_STUBS`. Neither defines
  [`DBARTS_REQUIRE_EXACT_ABI`](../../inst/include/dbarts/dbarts.h) (dropped under dec-B111), neither
  pins the hash, and neither reads `numReportedLocations` or `dbarts_draw`
  at all (`git -C <repo> grep`, both branches).
- bartCause (dbarts-1.0) has no compiled code against dbarts. bairrtt's
  LinkingTo names Rcpp and RcppEigen but not dbarts.

## Design

### The accessors and the rename

```
/// The fitted values a draw carries per observation: K, the category count,
/// on a multinomial sampler, and 1 on every other. A VALUE, never 0.
/// dbarts_sampler_run's train and test, and dbarts_sampler_predict's out,
/// are observations x this x draws x chains.
size_t dbarts_sampler_numFittedValuesPerObservation(
  const dbarts_sampler* sampler);
/// The sets of split counts a run writes per draw, one per forest that keeps
/// them: K on multinomial, the mean-forest count on an amplitude-coupled
/// model (2 on BCF), 1 elsewhere (a variance forest keeps none). A VALUE, never 0. dbarts_results' varcount is numPredictors x this
/// x draws x chains.
size_t dbarts_sampler_numVariableCountForests(const dbarts_sampler* sampler);
```

- Both are appended to [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h) after
  `dbarts_sampler_family`, with readable prototypes in the non-stub branch.
- Registration ([`DBARTS_API_REGISTER`](../../src/R_interface.cpp)), the stubs and the binding
  asserts expand from the list.
- The bodies return `shape().numReportedLocations` and
  `shape().numVariableCountForests`. The engine's own field names do not
  change; only shipped names do (dec-B163, dec-B165).
- The [`dbarts_draw`](../../inst/include/dbarts/dbarts.h) field `numReportedLocations` becomes
  `numFittedValuesPerObservation` in the same position.
  [`fillShippedDraw`](../../src/R_interface_bartcore_common.hpp) assigns the new name. The compiled
  readers move with it: the test consumer, `tests/cpp/test_capi.cpp` and the
  callback recipe in
  [vignettes/dbarts-as-a-component.Rmd](../../vignettes/dbarts-as-a-component.Rmd).
- No existing entry could carry these counts.
  - `dbarts_sampler_family` already says multinomial but not K.
  - `dbarts_sampler_numTrees` raises past the last forest.
  - A shape-struct query would add a third size-first struct to the layout
    fold and duplicate the seven size queries.
  - The callback carries both counts only during a run.

### Layouts (header text, per field)

F = `dbarts_sampler_numFittedValuesPerObservation` and V =
`dbarts_sampler_numVariableCountForests`. Observation (or predictor) varies
fastest, then the F (or V) axis, then draw, then chain. That is R's array
order, so `as.vector` of R's result equals the flat buffer.

- `train`: numObservations x F x numSamples x numChains. On multinomial each
  draw's F columns are category probabilities, rows summing to 1: the
  softmax of the K forests plus any category offset set from R.
- `test`: numTestObservations x F x numSamples x numChains. A category test
  offset set from R applies here, never to predict.
- `varcount`: numPredictors x V x numSamples x numChains, slab j forest j's
  counts. On BCF slab 0 is the prognostic forest (dec-B164). The run sets
  `numVariableCountForests` to V.
- `sigma`: numSamples x numChains, unchanged. On multinomial it is the pinned 1.
- `logLikelihood`: numObservations x numSamples x numChains, unchanged. On
  multinomial it is filled with the engine's quiet NaN.
- `k`, `varprobs`, `dispersion`, `residualDf`: unchanged, and untouched on
  multinomial, since no multinomial sampler carries them.
- predict `out`: xTest->numRows x F x S' x numChains, where S' is
  `numSavedSamples` with tree storage and 1 without. Probability scale on
  multinomial.
- predict `offsetTest` when F > 1 (dec-B166): xTest->numRows x F,
  column-major, entering before the softmax as R's matrix does. Every
  entry must be finite; a non-finite one raises with R's sentence. A null
  offset raises with R's sentence when the sampler carries a train or test
  category offset set from R. F = 1 keeps today's rule: numRows values,
  added after the fit.

Header edits beyond the fields:
- The run and predict paragraphs state the layouts above.
- The common-contracts bullet "Result and prediction layouts put samples and
  then chains in trailing dimensions" gains "after the fitted-value (or
  split-count) axis".
- The forest-index bullet's "states its count nowhere here" gains an
  exception: on multinomial the K forests are the K fitted values, and
  `dbarts_sampler_numVariableCountForests` states the count on both
  multi-forest models.
- The `dbarts_results` varcount paragraph ("this struct declares no forest
  count ... whatever the sampler's forest count is") is rewritten to V.
- `dbarts_draw`'s layout paragraph is restated. The F line names the
  accessor. The varcount sentence ("one slab for a single-forest model, K
  for a multinomial or multi-forest one"), false on the flat route today,
  becomes "V slabs, as `dbarts_sampler_numVariableCountForests` reports, on
  either route".

In [C_interface.cpp](../../src/C_interface.cpp):
- The stale stride comment in [`dbarts_sampler_run`](../../src/C_interface.cpp) ("1 for every
  dbarts.h-created sampler") is replaced.
- The [`dbarts_sampler_family`](../../src/C_interface.cpp) comment ("no entry here builds one
  yet") is restated for R-built handles.
- Predict's offset checks read only the shape, the holder's two owned
  category offsets and the caller's buffer. They run at the top of the
  entry, before [`callConvertingExceptions`](../../src/R_interface_bartcore_common.hpp) opens, and
  never as an `Rf_error` inside the captured body.
- For F > 1 the offset reaches the engine's predict as its category offset,
  and the slab add stays for F = 1. The sentences are shared with the R
  route's [`validateCategoryOffset`](../../src/R_interface_bartcore.cpp) and the
  [`bartcore_predict`](../../src/R_interface_bartcore.cpp) refusal, lifted into the bridge namespace
  if needed so the two routes cannot word them differently.

### The interrupt poll and the slow-count warning (dec-B167 to dec-B169)

Where the poll lives.
- `bartcore_checkInterrupt` and `bartcore_userInterrupted` move out of
  their file-static scope into the `bartcore_bridge` namespace. They are
  declared in [R_interface_bartcore_common.hpp](../../src/R_interface_bartcore_common.hpp) and stay
  defined in R_interface_bartcore.cpp.
- [`bartcore_run`](../../src/R_interface_bartcore.cpp) and [`dbarts_sampler_run`](../../src/C_interface.cpp) pass
  the same function. The flat run replaces its `{}` with it.
- The flat run passes `&stoppedByCallback` to the engine as
  [`bartcore_run`](../../src/R_interface_bartcore.cpp) does. A run the engine reports cancelled
  raises only when the stop was NOT the callback's. A callback's nonzero
  return keeps its documented contract (the entry returns normally), which
  the existing `capi_draw_reset(2L)` stop arms in test-capi.R pin.
- After the jump handling, a real cancel raises "dbarts_sampler_run:
  sampler run interrupted". That is the R route's sentence with the entry's
  prefix, raised where nothing the library owns is live.
- The sampler is then in the state a nonzero callback return leaves: cursors
  not advanced past written slots, so results and saved trees are
  discarded. The handle stays valid.
- The poll also fires at sweep boundaries, not only inside a count, exactly
  as on the R route. Every flat run becomes interruptible.

The interrupt test hook.
- A process-wide `std::atomic<int>` in R_interface_bartcore.cpp, armed
  through a third argument to
  [`bartcore_setMonotoneCountHooks`](../../src/R_interface_bartcore.cpp) (`interruptAfterPolls`; its
  registration arity in [R_interface.cpp](../../src/R_interface.cpp) moves 2 -> 3).
- When armed, the shared poll counts down and reports an interrupt on the
  armed poll without touching R's signal state. It disarms itself.
- Every test that arms it resets it to 0 in its cleanup, whatever the
  arm's outcome: it is process-wide, and a leftover count would interrupt
  an unrelated later run.
- Both routes get their first R-level interrupt test from it.
- That test asserts the wiring, the "interrupted" error and a usable sampler
  afterwards. It does not assert where the interrupt landed. The engine
  throttles the poll to one per ~100 ms after an immediate first call, so
  an armed count of 1 fires at a sweep boundary. Landing one inside a
  leaf-order count from R would need a count lasting over 100 ms, which is
  timing-dependent, and the count hooks carry no poll interval. The
  in-count relay is engine mechanics, already pinned without timing by
  tests/cpp `testMonotoneCountInterrupt` through its own poll-interval
  knob, so it stays there.

Where the once-per-sampler flag lives.
- A `bool slowCountWarned = false` member appended at the END of the
  holder, `dbarts_sampler_t` in
  [R_interface_bartcore_common.hpp](../../src/R_interface_bartcore_common.hpp). The creation sites build
  the holder by positional aggregate initialization, so a member placed
  anywhere else would shift every initializer after it. The struct is opaque
  in the shipped header, so this is no ABI change.
- After a run that returns normally, the flat entry reads
  `sampler.slowCountTally()` (the last run's, summed over chains). If it
  counted a slow count and the flag is false, the entry sets the flag and
  raises the warning.
- A holder lives as long as the engine handle. An R-side re-creation from a
  stored state builds a new holder, so a re-created sampler may warn once
  more; the header says so.
- The R route keeps warning per `$run` call, unchanged, and neither reads
  nor sets the flag.

How the warning is raised.
- The entry builds the same carrier [`bartcore_run`](../../src/R_interface_bartcore.cpp) builds.
  [`attachSlowCountTally`](../../src/R_interface_bartcore.cpp) is lifted into the bridge namespace
  beside the poll.
- It evaluates [`warnOnSlowCount`](../../R/bartcore.R) from the dbarts namespace on that
  carrier, so the class (`dbartsSlowCountWarning`), the sentence and the
  `tally` field are R's own.
- This happens last in the entry, after every buffer is released and every
  engine frame unwound.
- Any handler that exits on the warning turns it into a jump out of the
  entry: `options(warn = 2)`, `tryCatch(warning = )`, `expect_warning`, or
  a calling handler that stops or invokes a restart. The run's results are
  already complete and the sampler consistent when that happens. The header
  states this.

What stan4bart sees.
- A Ctrl-C during a fit used to be ignored until its `run` .Call returned,
  since stan4bart polls nowhere. Now the next poll inside a
  `dbarts_sampler_run` call converts it into an R error. That error
  longjmps through stan4bart's sweep loop and its .Call frame, so the
  newed `IterableBartResults` leaks and the stan4bart sampler stops
  mid-iteration.
- The same jump already follows any engine error from that call, so the
  path is not new; it becomes reachable by the user.
- The init path's `std::unique_ptr` and `std::vector` are skipped the same
  way.
- The warning arrives inside the .Call as an ordinary deferred R warning,
  shown at top level once per fit. Any exiting handler around stan4bart's
  fit (`warn = 2`, `tryCatch(warning = )`, `expect_warning` in its tests, a
  calling handler that stops) turns it into the interrupt's jump.
- Step 4 reviews these. The recommended fix is stan4bart-side and small:
  hold the run's results in an R-allocated, PROTECTed buffer in place of
  the raw `new`, so a jump leaks nothing. The alternative is a loop body
  that catches every C++ exception before any R jump can cross it.
  `R_UnwindProtect` around the whole loop is not recommended: the loop runs
  WALNUTS code that can throw C++ exceptions, and those must not cross the
  unwind-protect boundary. The fix lands as its own stan4bart commit on
  bartcore.
- A stan4bart test that drives `dbarts:::` hooks (the interrupt hook) skips
  when the installed dbarts hook's arity differs from the one it was written
  against, so stan4bart's suite does not break on an internal dbarts
  change.
- treatSens runs its fits through the same entry and gets the same review.

### Wrong sizes

No output size is declared anywhere: `dbarts_results` carries bare pointers,
and predict's `out` and `offsetTest` are bare pointers. A short buffer
remains the caller's crash under the header's "Validation is deliberately
partial" rule. The declared sizes that exist, a predictor source's numRows
and numColumns, are checked as today. Considered and not taken: a
caller-declared width field in `dbarts_results`. It covers run only, and
predict would need a signature change that forces consumer source edits.

### The ABI event

- Re-bake. The two appended entries move the signature token. The renamed
  field moves the layout fold, which folds field names. Both literals are
  re-baked in commit 1: `DBARTS_C_API_HASH` in the header and the
  [`dbarts_apiSignatureToken`](../../src/C_interface.cpp) assert, by the probe procedure in
  [5. Hash re-bake](dbarts-h-freeze.md#5-hash-re-bake). Commit 2 changes header comments only and
  moves neither.
- Version pair. It is held at 1/0: no version has shipped, and the header
  says the constants do not move before the first release.
  [tools/check-api-hash.sh](../../tools/check-api-hash.sh) prints its skip line until a tag exists.
- If the 1.0-0 tag exists before this lands, the additive parts (the two
  entries) alone would bump `DBARTS_C_API_MINOR`. The field rename would
  not fit a minor bump: renaming a field breaks a consumer's source. It
  must then either bump `DBARTS_C_API_MAJOR`, which every stub consumer's
  handshake refuses until rebuilt, or be dropped, keeping
  `numReportedLocations` as the field name and adding the new name only as a
  documented alias macro. Neither is wanted, which is why this lands before
  the tag. If the tag wins the race, the choice goes back to the
  maintainer.
- Pin sites: the header, C_interface.cpp, and test-capi.R's
  ["expect_identical(hashes$text"](../../inst/tinytest/test-capi.R) line. The outgoing
  `0x6380bf095d5cae3f` joins the file's stale-token block with a one-line
  reason. Dated mentions elsewhere stay.
- Existing entries whose writes or answers change:
  - Run's varcount writes widen from one slab to V on multinomial (K) and
    BCF (2). A stale binary that sized one slab and passed no train or test
    was safe on those handles and now overflows. No known consumer drives
    either through the flat API.
  - The callback's `numVariableCountForests` on the flat route moves from
    1 to V on the same handles.
  - Predict on a multinomial handle reads its offset K wide, and it raises
    on a null or non-finite one where R does.
  - Every run can now raise on an interrupt and emit one warning per
    sampler.
- What a stale consumer binary sees.
  - With stubs and no exact-ABI flag (stan4bart, treatSens): the version
    handshake passes, and every entry resolves by name with an unchanged
    signature. Neither reads `dbarts_draw`, so the field rename does not
    reach them. On gaussian and probit samplers they run bitwise as
    before, plus the interrupt and warning behaviour, which arrives with the
    library and not with the rebuild. Step 4 runs exactly this case first.
  - With `DBARTS_REQUIRE_EXACT_ABI`: the first stub call raises "dbarts C ABI
    mismatch ... rebuild".
  - A binary that reads the renamed field still reads the same offset, so a
    stale one keeps working. A rebuild of its source fails to compile until
    the name is updated.
  - A consumer built against the new header with an older dbarts installed
    fails only when it first calls a new entry ("not provided by package
    'dbarts'").

## Constraints

- Gates (neutral class):
  - Full tinytest on a `--preclean` private-library install, and tests/cpp
    (test_capi.cpp carries the rename).
  - The equivalence trio, labelled a formality (no harness drives a flat
    entry), IDENTICAL against the baselines in benchmarks/baselines/MANIFEST.
  - The R-loaded AddressSanitizer run of test-capi.R per the plans README.
  - A C99 `-pedantic -fsyntax-only` compile of the header's prototype view.
  - `R CMD check --as-cran` for commit 3 (man/ and the vignette touched).
- Frozen: no existing entry signature or enumerator changes. The one struct
  change is the field rename. No engine file is touched.
- Out of scope: flat creation; flat setters for counts and category offsets;
  a flat multinomial log-likelihood; changing the R route's per-call
  slow-count warning.

## Steps

1. Commit 1, the ABI event. Code, header and tests go together, because
   the hash pin makes them inseparable.
   - Header: the two accessors, the field rename, the layout and offset text,
     the bullet edits.
   - [C_interface.cpp](../../src/C_interface.cpp):
     - The accessor bodies.
     - `engineResults.numVariableCountForests` set to V.
     - Predict's pre-capture offset checks and the F > 1 offset pass-through.
     - The two comments.
   - [`fillShippedDraw`](../../src/R_interface_bartcore_common.hpp) and the shared offset sentences.
   - The re-bake of both literals.
   - [consumer.c](../../inst/tinytest/capi/consumer.c):
     - `capi_run_canaried(ptr, burn, samples, F, V, tailFactor)` sizes every
       buffer from the F and V it is handed, never from the accessors. The
       test computes them in R: F as `ncol` of the counts matrix the
       sampler was built from (1 for gaussian), V as the sampler's own
       forest count from R.
     - All nine `dbarts_results` pointers are passed. Each buffer has a
       tail of `tailFactor` x body (the test uses K) filled with a quiet NaN
       with payload 0x7FF8DEADBEEF0001, or 0xDEADBEEF for varcount.
     - It returns per channel the body, "tail intact" (memcmp), "body
       untouched" and "body fully written".
     - `tailFactor` 0 `malloc`s each body at exactly its size, for the
       AddressSanitizer case.
     - `capi_predict_canaried(ptr, x, offset, F, tailFactor)` does the same
       for predict and also returns the status.
     - [`capi_dims`](../../inst/tinytest/capi/consumer.c) gains both accessors.
     - [`capi_draw_report`](../../inst/tinytest/capi/consumer.c) renames its field and gains
       `numVariableCountForests`.
     - The mean callback reads the new field name.
   - `tests/cpp/test_capi.cpp` and the vignette recipe take the new field
     name.
   - [test-capi.R](../../inst/tinytest/test-capi.R):
     - Hash pins.
     - The accessors checked separately against the R-side counts: F is K, 1,
       1, 1 and V is K, 1, 1, 2 on multinomial, gaussian, probit and BCF.
     - Gaussian arm: two identically seeded two-chain samplers with test
       rows. Flat run on one and `$run` on the other.
       - `expect_identical` on `as.vector` of sigma, train, test and
         varcount.
       - logLikelihood is finite.
       - k, varprobs, dispersion and residualDf are body-untouched.
       - Every tail is intact.
     - Multinomial arm (K = 3, two chains), the same comparison:
       - sigma is exactly 1.
       - Train and test match R, with rows summing to 1.
       - varcount matches R's p x K x S x C.
       - Every logLikelihood word is bitwise
         `std::numeric_limits<double>::quiet_NaN()`, distinct from the
         canary.
       - k, varprobs, dispersion and residualDf are body-untouched.
       - Tails are intact.
     - BCF arm: flat varcount equals R's p x 2 x S x C.
     - Three-forest arm: a `forests = list(...)` sampler with three mean
       forests. The accessor answers 3, and flat varcount equals R's
       p x 3 x S x C.
     - Predict on gaussian and multinomial, against `$predict` on the same
       sampler, with and without tree storage.
     - Offset arms:
       - A gaussian offset still adds.
       - A multinomial nTest x K matrix equals `$predict(x, offset = m)`
         bitwise.
       - A non-finite entry raises.
       - A null offset raises on a sampler given `$setCategoryOffset`, and
         separately on one given only `$setCategoryTestOffset`.
       - An all-zero matrix on those samplers equals R's.
     - setResponse and setOffset return 0 on the multinomial handle.
     - The callback reports F and V equal to the accessors on multinomial,
       BCF and the three-forest sampler.
     - The existing `capi_draw_reset(2L)` callback-stop arms still return
       normally from the flat run, with no "interrupted" error.
   - Run the gaussian parity check first. If flat and R runs differ there,
     stop: that is a separate finding.
   - Mutation proofs, each reverted and the file `touch`ed:
     - An accessor returning 1 fails its check.
     - Dropping the V line fails varcount parity.
     - Dropping each offset check fails its arm.
   - AddressSanitizer: on the pre-slice build, a multinomial run through
     `tailFactor` 0, sized F = V = 1, reports a heap-buffer-overflow. On the
     slice there is none.
   - Gates: Constraints, plus `air format --check .` and lintr on the touched
     R files.
2. Commit 2, the flat run's monotone gaps (no hash move).
   - The bridge poll and `attachSlowCountTally` lifted into `bartcore_bridge`.
   - The `interruptAfterPolls` hook and the arity change.
   - The holder flag.
   - The flat run's poll, `&stoppedByCallback`, the raise on a real cancel
     only, and the warning tail.
   - Header comments on [`dbarts_sampler_run`](../../inst/include/dbarts/dbarts.h): it can raise on an
     interrupt, leaving the callback-abort state; it raises
     `dbartsSlowCountWarning` at most once per sampler (holder), and once
     more after an R-side re-creation; any exiting warning handler turns
     it into a jump after complete results.
   - consumer.c gains `capi_run_plain(ptr, burn, samples)` (no result
     buffers), for the warning and interrupt arms.
   - Tests, in test-monotone.R beside its slow-count block, using
     ["countHooks"](../../inst/tinytest/test-monotone.R) and that file's `slowSampler`. A
     `compileCapiConsumer` call there skips as test-capi.R does.
     - Interrupt, R route: `$run` with `interruptAfterPolls` armed raises
       "sampler run interrupted". A second `$run` afterwards succeeds.
     - Interrupt, flat route: the same through `capi_run_plain`, on a
       "leaf" sampler. It asserts only the "interrupted" error and a usable
       sampler afterwards. Where the interrupt lands is left to
       `testMonotoneCountInterrupt` (see "The interrupt test hook").
     - Each interrupt arm resets `interruptAfterPolls` to 0 in its cleanup
       (`on.exit` or a final `countHooks` call reached on every path).
     - A callback-stop flat run with the hook unarmed returns normally.
     - Slow count, flat route: with the threshold at -1, the first
       `capi_run_plain` on a "leaf" sampler yields exactly one warning,
       counted with `withCallingHandlers`, inheriting
       `dbartsSlowCountWarning` with the R route's sentence. Three more
       calls on the same sampler yield none. A fresh sampler warns once
       again. A "joint" sampler never warns.
     - The R route still warns on every `$run` with a slow count, pinned so
       the flag cannot leak into it.
   - Mutation proofs: dropping the flag fails once-per-sampler, and passing
     `{}` again fails the flat interrupt arm.
   - Gates: as commit 1, with tests/cpp's `testMonotoneCountInterrupt`
     named in the run.
3. Commit 3, docs.
   - [2. Reach](../design/feature-matrix.md#2-reach): the multinom flat cell
     goes M -> S, citing the accessor.
   - [Footnotes](../design/feature-matrix.md#footnotes): [f4] is rewritten.
   - [Gaps](../design/feature-matrix.md#gaps): the "Flat C reach for the K-forest
     softmax family" row is deleted.
   - docs/plans/bartcore-landing/changes.md:
     - chg-C03 goes 25 -> 27 entry points.
     - chg-C25 becomes six entries neither consumer calls, adding both
       accessors.
     - chg-C30's field list carries the renamed field.
     - A new chg-C34 records the accessors, the rename, the K-wide layouts,
       the per-forest varcount, the predict offset rule, the re-bake with
       the pair at 1/0, and dec-B160, dec-B163 to dec-B166.
     - A new chg-C35 records the flat run's interrupt and once-per-sampler
       warning (dec-B167 to dec-B169).
   - docs/design/per-draw-callbacks.md: its struct listing takes the new
     field name.
   - TODO:
     - monotone-count-host-interrupt is removed, closed by this slice.
     - The multinomial-doors dbarts.h clause says run and predict are open
       (dec-B160) and creation stays the door.
   - man/dbarts-embedding.Rd: one sentence on the two accessors, the
     interrupt and the warning in the "What the header reaches" paragraph.
   - No NEWS text. The 1.0-0 C API bullet introduces the header as new;
     monotone and multinomial are new in 1.0-0 too.
   - Gates: `tools/check-doc-freshness.R`, `tools/check-rc-codoc.R`, and
     `R CMD check --as-cran` from a clean tarball, each on its own exit
     status.
4. Consumer runs and the stan4bart review (no dbarts commit).
   - stan4bart bartcore, first with the stale binary:
     - Install the post-slice dbarts with `--preclean` into a library that
       still holds the PRE-slice stan4bart binary.
     - Run its tinytest suite at_home and expect a pass.
   - Then rebuild stan4bart with `--preclean`, run the suite again and
     expect a pass with no source change.
   - The stan4bart review (opus), against "What stan4bart sees":
     - Interrupt a long fit once by hand. Confirm the error, and confirm
       that a fresh fit afterwards works.
     - List what leaks, and land the stan4bart-side cleanup as its own
       commit on bartcore, with a test using dbarts' interrupt hook.
     - Confirm one slow-count warning per fit on a monotone BART component
       if stan4bart exposes one; otherwise record that it cannot reach a
       monotone leaf.
   - treatSens dbarts-1.0: the same stale and rebuilt runs, and the same
     review of its run call sites.
   - bartCause and bairrtt: nothing (neither links dbarts).
   - A failed stale run, or any needed consumer edit beyond the review's
     cleanup, stops the slice.
5. Commit 5, records.
   - This plan's Landing note, pinned to the landed shas, with the consumer
     runs and the stan4bart commit.
   - Status line; INDEX row status; dec-B160 and dec-B163 to dec-B169 marked
     landed as the register's practice is.
   - Plain commit messages throughout.

## Verification

- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e
  'tinytest::run_test_file("inst/tinytest/test-capi.R")'` and the same for
  test-monotone.R, test-callback-example.R: zero failures. Then
  `tinytest::test_package("dbarts")`: zero failures.
- `cd tests/cpp && make && ./test_bartcore`: all pass.
- The equivalence trio compares print the full count of "identical draws"
  lines and no "max |z|".
- The AddressSanitizer driver over test-capi.R and test-monotone.R reports
  zero diagnostics on the slice, and the overflow on the pre-slice case.
- The C99 syntax check of `inst/include/dbarts/dbarts.h` is clean.
- `sh tools/check-api-hash.sh` prints its skip line, or passes with the
  minor bump if a tag exists.
- Each consumer suite passes stale and rebuilt; the stan4bart interrupt
  test passes on its branch.

## Rulings

All ruled 2026-10-01 in [decisions.md](../decisions.md); this plan carries
them out.
- dec-B160: flat run and predict support multinomial samplers built in R; an
  ABI event every LinkingTo consumer rebuilds against; no flat creation.
- dec-B163: the K accessor is `dbarts_sampler_numFittedValuesPerObservation`,
  and the callback field `numReportedLocations` is renamed to match in the
  same event.
- dec-B164: varcount holds one set per forest that keeps split counts (K on
  multinomial, the mean-forest count on an amplitude model, 2 on BCF),
  matching R's run on every model.
- dec-B165: that width's accessor is `dbarts_sampler_numVariableCountForests`.
- dec-B166: flat predict reads a multinomial offset as rows x K, requires one
  where R does, and refuses non-finite entries as R does.
- dec-B167: the flat run's interrupt and slow-count gaps ride this change.
- dec-B168: the flat run raises R's `dbartsSlowCountWarning` itself; no
  results field.
- dec-B169: that warning fires once per sampler.

No question remains open.

## Landing note (2026-10-01)

Landed on top of the review-3 bridge fixes (39b76a6b), as four dbarts
commits and one stan4bart commit.

- 5233f0ae, step 1, the ABI event. The two accessors, the field rename,
  the layout and offset text, the predict offset checks ahead of the
  captured body, V in the flat run, and the re-bake: `DBARTS_C_API_HASH`
  0x6380bf095d5cae3f -> 0xa7415a6f1bcc93c3, signature token
  0xb6f41cfcbd996897 -> 0xfbf29fc67c22558b, the pair held at 1/0. The
  gaussian parity check passed first, every channel bitwise. One change
  the plan did not list: the flat run wrote k on a sampler with no k
  hyperprior (the R run passes null there, and the header already said
  "left untouched otherwise"), so the gaussian and multinomial arms'
  k-untouched checks failed; the flat run now passes null k off
  `kIsSampled`. Neither consumer passes k without a hyperprior.
  Mutation proofs: an accessor returning 1 fails 3 checks, dropping the V
  line 6, each offset refusal (train, test, non-finite) 1, and dropping
  the offset pass-through 2. AddressSanitizer: the slice-1 build, a
  multinomial run sized F = V = 1 through `tailFactor` 0, reports a
  heap-buffer-overflow; the slice build, sized by the accessors, reports
  none.
- 58bface9, step 2. Mutation proofs: dropping the holder flag fails the
  once-per-sampler arm; passing `{}` again fails the flat interrupt arm.
- ae1f46cd, step 3, docs. TODO lands with these records (the brief put
  it here rather than in step 3).
- fdddeff2, a lint fix to step 1's test helper.

Gates at fdddeff2, `--preclean` private library: tinytest 204 files,
11043 tests, 0 failures; tests/cpp all pass, the C99 and C++ header
compiles included, and again under `-fsanitize=address,undefined` with
no diagnostic; the R-loaded AddressSanitizer run of test-capi.R,
test-monotone.R and test-callback-example.R, zero diagnostics; the
equivalence trio bitwise (55/55 under `--strict-coverage`, BCF 15/15,
multinomial 11/11); every exact gate in exact-gates.yaml in quick mode,
the monotone successive-conditional and enumeration gates under both
priors included, PASS; lintr, air, check-rc-codoc, check-win-drift and
check-doc-freshness clean; check-api-hash skips (no tag); `R CMD check
--as-cran` one NOTE (Date), run with `--no-tests` beside the full
tinytest run above.

Step 4, consumers, each installed from `git archive` into the private
library.
- Stale binaries: stan4bart bartcore (9a3be93) and treatSens dbarts-1.0
  (7cc6a0f), built against the slice-1 dbarts, then the post-slice dbarts
  installed `--preclean` beside them. stan4bart tinytest 566/0, treatSens
  testthat 194/0.
- Rebuilt `--preclean` with no source change: the same counts, 0
  failures.
- stan4bart review. An interrupt injected through the count hooks (the
  fifth poll, inside warmup) ends the fit with "dbarts_sampler_run:
  sampler run interrupted", and a fresh fit afterwards matches one made
  before it. stan4bart reaches a monotone component through `bart_args`;
  under prior "leaf" with the threshold lowered, one fit warns exactly
  once. The cleanup lands on stan4bart as a2a94d4 on branch wt/capi
  (worktree .claude/worktrees/capi, off bartcore 9a3be93, not pushed): the
  run's draws land in the PROTECTed R list it returns, in place of the raw
  `new`, and creation's first draw in R's transient storage; test-24
  drives the hook and skips when the installed hook's arity is not 3.
  stan4bart's suite on that commit: 570/0. Left as is: creation's
  `std::unique_ptr<Sampler>` still leaks if an interrupt lands in its one
  sweep, and an error from the R callback's `Rf_eval` longjmps past the
  run as before.
- treatSens review: its sensitivity loop holds raw `new[]` buffers (the
  grid cells, the train and test stores, the per-cell state arrays)
  across `dbarts_sampler_run`, so an interrupt there now leaks them, the
  same class as stan4bart's; no slow-count warning is reachable (no
  monotone component). Not fixed here.
