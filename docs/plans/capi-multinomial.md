# capi-multinomial

Status: PROPOSED 2026-10-01. Ruling: dec-B160 in [decisions.md](../decisions.md).

agent: opus (commit 1, the header, the C entry file and the test consumer);
  sonnet (commit 2, docs; step 3, the consumer runs; commit 4, records).
  Serial: 1, 2, 3, 4.
rng: NEUTRAL. No engine file changes. The run entry already hands the
  engine the K-wide location stride. What moves is the contract, one or two
  accessors (Q2), the predict offset rule (Q3) and the varcount width the
  flat caller declares (Q2). Every sampler with one location and one varcount
  slab runs byte for byte what it runs today.
window: nothing else edits `inst/include/dbarts/dbarts.h` while this is
  open. Lands before the 1.0-0 tag, or bumps the minor version (Constraints).
budget: header ~+70 -20; C entry file ~+45 -12; test consumer ~+230; test
  file ~+200; docs ~+40 -15; records ~+40. About 650 lines; plan on 950-1250.

## Goal

A multinomial sampler built in R runs and predicts through the flat C API.
The header gains `dbarts_sampler_numReportedLocations`, which reports K
before any run, and it documents the K-wide layouts of
[`dbarts_sampler_run`](../../inst/include/dbarts/dbarts.h) and
[`dbarts_sampler_predict`](../../inst/include/dbarts/dbarts.h). A caller who sizes its buffers by
the documented layouts gets bitwise what R's `$run` and `$predict` return on
every channel the flat struct carries. One exception remains unless Q2
takes (c): varcount on a multi-forest amplitude sampler (BCF), where R
returns p x 2 x S x C and flat returns the prognostic p x S x C. Under Q3(a), a
multinomial sampler carrying a category offset set from R predicts through R
only. Out of scope: a flat creation path, and any flat count or
category-offset setter. Both stay with the R object (TODO multinomial-doors).

## Context

The defect. Every handle is an R sampler's pointer
([Landing note, S2 (2026-09-08)](pure-c-header.md#landing-note-s2-2026-09-08) removed flat creation), and
`dbarts()` builds multinomial samplers. The header treats them as in
contract: [`dbarts_sampler_family`](../../inst/include/dbarts/dbarts.h) reports
`DBARTS_FAMILY_MULTINOMIAL`, [`dbarts_sampler_getLatents`](../../inst/include/dbarts/dbarts.h) documents
its multinomial answer, and the [`dbarts_results`](../../inst/include/dbarts/dbarts.h) comment
describes multinomial log-likelihood. But its buffer sizes assume one
location. [`dbarts_sampler_run`](../../src/C_interface.cpp) passes the sampler's location
count to the engine, and [`Chain::storeSample`](../../src/bartcore/chain.hpp) writes n x K per draw
whatever the caller sized. So train and test come back n x K x S x C into
buffers documented as n x S x C, a heap overflow. The verification probe
(n = 30, K = 3, S = 4, C = 2) wrote 720 of a documented 240 train entries and
crashed later in unrelated R code. [`dbarts_sampler_predict`](../../src/C_interface.cpp) likewise
writes nTest x K x S x C. A non-null offset there is added slab by slab at
stride nTest: an additive shift on probability channels, at the wrong
entries. Flat predict also ignores a category offset the sampler carries
from R (see below) and returns the offset-free surface where R refuses.
varcount is safe today: the run leaves `numVariableCountForests` at 1, so the
engine writes only the first category's slab. No current entry reports K
before a run. [`dbarts_sampler_numTrees`](../../inst/include/dbarts/dbarts.h) raises past the last
forest, and the [`dbarts_draw`](../../inst/include/dbarts/dbarts.h) field `numReportedLocations`
arrives only inside a run, too late to size buffers.

What R returns (probed against the installed build: n = 30, K = 3, nTest =
10, four draws, two chains). [`bartcore_run`](../../src/R_interface_bartcore.cpp) gives:

| channel | R shape | content |
|---|---|---|
| sigma | S x C | the pinned 1 |
| train | n x K x S x C | softmax probabilities; each row of K sums to 1 |
| test | nTest x K x S x C | the same, under any category test offset set from R |
| varcount | p x K x S x C | slab k is category k's forest (BCF: p x 2 x S x C) |
| k, varprobs | NULL | a k hyperprior and DART are refused on multinomial |

[`predictFromSource`](../../src/R_interface_bartcore.cpp) returns nTest x K x S x C with saved
trees and nTest x K x 1 x C without. Its offset must be an nTest x K
per-category matrix entering before the softmax. A flat vector is refused,
because after the blend it moves values off the simplex and before it a
common shift is the softmax's null direction. A sampler carrying a train or
test category offset set from R (the holder's owned category offsets,
installed by [`bartcore_setCategoryOffset`](../../src/R_interface_bartcore.cpp) and
[`bartcore_setCategoryTestOffset`](../../src/R_interface_bartcore.cpp)) refuses a predict that names no
offset, in [`bartcore_predict`](../../src/R_interface_bartcore.cpp). The predicted rows are not its
rows, and an all-zero matrix is how a caller asks for the offset-free
surface. Category k is column k of the sampler's count matrix.

The other entries on a multinomial handle already behave, and nothing here
changes them. setResponse and setOffset return 0
([`responseConduitIsFixed`](../../src/R_interface_bartcore.cpp)), as do setSigma and getLatents.
printTrees and numTrees take forests 0..K-1. The per-draw callback sees
`numReportedLocations` = K with n x K train, and `numVariableCountForests` = 1
on the flat route.

Consumers today. stan4bart (branch bartcore at 9a3be93) calls run, predict,
setOffset, setSigma, getLatents, setTreeStorage, sampleTreesFromPrior,
printTrees, destroy, the size queries and the version pair. It builds only
gaussian and probit samplers, sizes run's train as n, and passes a non-null
offsetTest to predict. treatSens is also a C consumer, contrary to the
brief. On dbarts-1.0 at 7cc6a0f it carries LinkingTo dbarts and calls run,
setResponse, setOffset, setSigma, setNumThreads, setVerbose and destroy, on
gaussian and probit samplers only. Both use `DBARTS_USE_STUBS`; neither
defines [`DBARTS_REQUIRE_EXACT_ABI`](../../inst/include/dbarts/dbarts.h) (dropped under dec-B111), and
neither pins the hash. bartCause (dbarts-1.0) has no compiled code against
dbarts. bairrtt's LinkingTo names Rcpp and RcppEigen but not dbarts. No
shipped consumer reaches a multinomial or BCF handle through the flat API.

## Design

### The accessor

```
/// The per-observation locations each draw carries, L: K, the category
/// count, on a multinomial sampler, and 1 on every other. A VALUE, never 0.
/// dbarts_sampler_run's train and test, and dbarts_sampler_predict's out,
/// are sized by it.
size_t dbarts_sampler_numReportedLocations(const dbarts_sampler* sampler);
```

It is appended to [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h) after
`dbarts_sampler_family`, with a readable prototype in the non-stub branch.
Registration ([`DBARTS_API_REGISTER`](../../src/R_interface.cpp)), the stubs and the binding
asserts all expand from the list, so nothing else is hand-added. The body
returns `shape().numReportedLocations`. The name matches the [`dbarts_draw`](../../inst/include/dbarts/dbarts.h)
field that already carries L (Q1). It is a location count, not a forest
count: BCF has two forests and answers 1. Under Q2(c) a second VALUE entry,
`dbarts_sampler_numVariableCountForests`, mirrors the other draw field and
returns `shape().numVariableCountForests`: K on multinomial, the forest count
on BCF, 1 otherwise.

An existing entry cannot carry K:
- `dbarts_sampler_family` already says multinomial but not K. Packing K into
  its int would make a VALUE entry mean two things.
- `dbarts_sampler_numTrees` raises past the last forest, so probing for K
  means catching a longjmp. A forest count is also not L (BCF).
- A shape struct filled by one query would subsume the seven size queries,
  but it adds a third size-first struct to the layout fold for one number and
  duplicates entries that stay.
- The draw callback carries L only during a run.

### Layouts (header text, per field)

L = `dbarts_sampler_numReportedLocations`. V = the varcount slab count (Q2).
Observation is fastest, then location, then draw, then chain. That is R's
array order, so `as.vector` of R's result equals the flat buffer.

- `train`: numObservations x L x numSamples x numChains. On multinomial each
  draw's L columns are category probabilities, rows summing to 1: the
  softmax of the K forests plus any category offset set from R.
- `test`: numTestObservations x L x numSamples x numChains. A category test
  offset set from R applies here, never to predict.
- `varcount`: numPredictors x V x numSamples x numChains, forest-major within
  a draw. Q2(a): V = L, so K on multinomial and 1 on BCF. Q2(c): V = the
  sampler's varcount forest count, so K on multinomial and 2 on BCF. Q2(b):
  V = 1 everywhere, today's rule. The run sets `numVariableCountForests` to V.
- `sigma`: numSamples x numChains, unchanged. On multinomial it is the pinned 1.
- `logLikelihood`: numObservations x numSamples x numChains, unchanged. On
  multinomial it is filled with the engine's quiet NaN (its existing
  sentence stays).
- `k`, `varprobs`, `dispersion`, `residualDf`: unchanged, and left untouched
  on multinomial. No multinomial sampler carries any of them.
- predict `out`: xTest->numRows x L x S' x numChains, where S' is
  `numSavedSamples` with tree storage and 1 without. Probability scale on
  multinomial; otherwise as today.
- predict `offsetTest` when L > 1 (Q3):
  - Q3(a): a non-null offset raises, an argument error, since null works. A
    sampler carrying a category offset set from R returns capability 0
    with `out` untouched, since no argument the flat entry can take would
    work.
  - Q3(b): the offset is xTest->numRows x L, entering before the softmax
    exactly as R's matrix does. A null offset on a sampler carrying an
    R-set category offset raises with R's sentence, an argument error,
    since an all-zero buffer works.
  - Either way the check reads only the shape and the holder. It is
    raised (or returned) before [`callConvertingExceptions`](../../src/R_interface_bartcore_common.hpp)
    opens, never by an `Rf_error` inside the captured body.

Header edits beyond the fields:
- The run and predict paragraphs state the layouts above.
- The common-contracts bullet "Result and prediction layouts put samples and
  then chains in trailing dimensions" gains "after L locations".
- The forest-index bullet's "states its count nowhere here" gains an
  exception: a multinomial sampler's K forests are its L locations.
- The `dbarts_results` varcount paragraph ("this struct declares no forest
  count ... whatever the sampler's forest count is") is rewritten to the
  chosen V rule.
- `dbarts_draw`'s L line points at the accessor. Its varcount sentence
  ("one slab for a single-forest model, K for a multinomial or multi-forest
  one") is false on the flat route today, which declares 1. It becomes "the
  slabs the run declared": V on the flat run and one per forest on the R route.

In [C_interface.cpp](../../src/C_interface.cpp):
- The stale stride comment in [`dbarts_sampler_run`](../../src/C_interface.cpp) ("1 for every
  dbarts.h-created sampler") is replaced.
- The comment in [`dbarts_sampler_family`](../../src/C_interface.cpp) ("no entry here builds one
  yet") is restated for R-built handles.

### Wrong sizes

No output size is declared anywhere: `dbarts_results` carries bare pointers
and predict's `out` and `offsetTest` are bare pointers. So there is nothing
to check, and a short buffer remains the caller's crash under the header's
"Validation is deliberately partial" rule. The declared sizes that exist, a
predictor source's numRows and numColumns, are checked as today. Considered
and not taken: a caller-declared L field in `dbarts_results` that the run
would refuse on mismatch. It covers run only. Predict would need a new
parameter, a signature change that forces consumer source edits, where the
ruling asks only for a rebuild.

### The ABI event

- Appending to the list moves the signature token and the full token. Both
  literals are re-baked in the same commit: `DBARTS_C_API_HASH` in the
  header and the [`dbarts_apiSignatureToken`](../../src/C_interface.cpp) assert, by the probe
  procedure in [5. Hash re-bake](dbarts-h-freeze.md#5-hash-re-bake). No struct, enumerator or
  layout moves, so the fold in [`dbarts_apiToken`](../../src/C_interface.cpp) is otherwise
  unchanged.
- Version pair: held at 1/0. No version has shipped, and the header says the
  constants do not move before the first release. If the 1.0-0 tag exists
  when this lands, the same commit bumps `DBARTS_C_API_MINOR` to 1 instead,
  an additive entry. [tools/check-api-hash.sh](../../tools/check-api-hash.sh) prints its skip line until
  a tag exists and enforces exactly that afterwards.
- Pin sites: the header, C_interface.cpp, and test-capi.R's
  ["expect_identical(hashes$text"](../../inst/tinytest/test-capi.R) line. The outgoing
  `0x6380bf095d5cae3f` joins the file's stale-token block with a one-line
  reason. Other mentions of it (TODO, per-draw-callbacks.md,
  bartcore-landing) are dated history and stay.
- Existing entries whose writes or answers change (not only the new entry):
  - Q2(a) or (c) widens run's varcount writes on a multinomial handle from
    one slab to K. A stale binary that sized varcount p x S x C, and passed
    no train or test, was safe on that handle and now overflows. Q2(c) does
    the same on a BCF handle (1 to 2 slabs).
  - The callback's `numVariableCountForests` on the flat route moves from 1
    to V on the same handles.
  - Predict on a multinomial handle changes under either Q3 option: a
    non-null offset raises or is read K wide, and an R-set category offset
    makes it return 0 or raise.
  - Gaussian, probit and every other one-location, one-slab handle is
    unchanged.
- What a stale consumer binary sees:
  - With stubs and no exact-ABI flag (stan4bart, treatSens): the version
    handshake passes, and every entry it calls resolves by name with an
    unchanged signature and struct layout. On its gaussian and probit
    samplers it runs bitwise as before. Step 3 runs exactly this case before
    rebuilding.
  - With `DBARTS_REQUIRE_EXACT_ABI`: the first stub call raises "dbarts C
    ABI mismatch ... rebuild".
  - A consumer built against the new header with an older dbarts installed
    fails only when it first calls a new entry ("not provided by package
    'dbarts'"). The pre-1.0 version pair cannot catch that; after 1.0-0 the
    minor bump does.
  - So the "rebuild once" in the ruling is lockstep hygiene. No consumer
    source changes, and a skipped rebuild breaks nothing a known consumer
    does.

## Constraints

- Gates (neutral class):
  - Full tinytest on a `--preclean` private-library install, and tests/cpp.
  - The equivalence trio, labelled a formality (no harness drives a flat
    entry), comes back IDENTICAL against the baselines named in
    benchmarks/baselines/MANIFEST.
  - The R-loaded AddressSanitizer run of test-capi.R per the plans README:
    the defect is a heap overflow, and the new code is reachable only
    through `.Call`.
  - A C99 `-pedantic -fsyntax-only` compile of the header's prototype view.
- Frozen: no existing signature, struct layout or enumerator changes. No
  engine file is touched.
- Out of scope: flat creation; flat setters for counts, the category offset
  or the category test offset; a flat log-likelihood for multinomial; the
  monotone interrupt and slow-count work (Q4).

## Steps

1. Commit 1, the ABI event. Code, header and tests go together, because
   the hash pin makes them inseparable.
   - Header: the accessor(s), the layout text above, the bullet edits.
   - [C_interface.cpp](../../src/C_interface.cpp):
     - The accessor bodies.
     - `engineResults.numVariableCountForests` set to V in run.
     - The Q3 checks, made at the top of predict from the shape and the
       holder before the captured body. The empty-store refusal stays where
       it is.
     - Under Q3(b), the offset is passed to the engine's predict as its
       category offset when L > 1, and the slab add is kept for L = 1.
     - The two comments.
   - The re-bake of both literals.
   - [consumer.c](../../inst/tinytest/capi/consumer.c):
     - `capi_run_canaried(ptr, burn, samples, K, tailFactor)` sizes every
       buffer from the K it is handed, never from the accessor. That K is
       computed in R as `ncol` of the counts matrix the sampler was built
       from, 1 for gaussian. V is handed the same way.
     - All nine `dbarts_results` pointers are passed. Each buffer has a
       tail of `tailFactor` x body (the test uses K, so an overrun up to K
       times lands in the canary) filled with a quiet NaN carrying a nonzero
       payload (0x7FF8DEADBEEF0001), and 0xDEADBEEF for varcount.
     - It returns per channel: the body, "tail intact" (memcmp), "body
       untouched" (every word still canary) and "body fully written" (no
       canary word left).
     - `tailFactor` 0 allocates each body with `malloc` at exactly its size
       and no tail. That is the AddressSanitizer case.
     - `capi_predict_canaried(ptr, x, offset, K, tailFactor)` does the same
       for predict and also returns the status.
     - [`capi_dims`](../../inst/tinytest/capi/consumer.c) gains the accessor(s).
     - [`capi_draw_report`](../../inst/tinytest/capi/consumer.c) gains `numVariableCountForests`.
   - [test-capi.R](../../inst/tinytest/test-capi.R):
     - Hash pins.
     - The accessor is checked separately against the R-side K: K on
       multinomial and 1 on gaussian, probit and BCF (the two-forest case).
       Under Q2(c), the second accessor gives K, 2 and 1.
     - Gaussian arm: two identically seeded two-chain samplers with test
       rows. Flat run on one and `$run` on the other.
       - `expect_identical` on `as.vector` of sigma, train, test and
         varcount.
       - logLikelihood is finite.
       - k, varprobs, dispersion and residualDf come back body-untouched.
       - Every tail is intact.
     - Multinomial arm (K = 3, two chains), the same comparison:
       - sigma is exactly 1.
       - Train and test match R, with rows summing to 1.
       - varcount matches R's p x K x S x C under Q2(a) or (c), and is the
         first slab under (b).
       - Every logLikelihood word is bitwise
         `std::numeric_limits<double>::quiet_NaN()`, which differs from the
         canary payload.
       - k, varprobs, dispersion and residualDf are body-untouched.
       - Every tail is intact.
     - Predict on each, against `$predict` on the same sampler after the
       run, with and without tree storage (nTest x K x 1 x C without).
     - Offset arms:
       - Gaussian with an offset still adds.
       - Q3(a): multinomial with a non-null offset raises. A multinomial
         sampler given a category offset through `$setCategoryOffset`, and
         separately one given only `$setCategoryTestOffset`, returns 0 with
         `out` body-untouched.
       - Q3(b): an nTest x K matrix matches `$predict(x, offset = m)`
         bitwise. A null offset on those two samplers raises.
     - setResponse and setOffset return 0 on the multinomial handle (the
       arms capi-shape.md recorded as unreachable are now reached).
     - The callback's `num.reported.locations` equals K on multinomial.
       Its `numVariableCountForests` equals V on multinomial, and on BCF
       under (c), and is 1 on gaussian.
   - Run the gaussian parity check first. If flat and R runs differ there,
     stop: that is a separate finding, not this slice's.
   - Mutation proofs, each reverted and the file `touch`ed:
     - The accessor returning 1 fails the accessor check, and the K-sized
       arm with it.
     - Dropping the varcount line fails varcount parity.
     - Dropping each Q3 check fails its arm.
   - AddressSanitizer discrimination: on the pre-slice build, a multinomial
     run through `tailFactor` 0, sized by the old contract (K = 1), must
     report a heap-buffer-overflow. On the slice it must report none.
   - Gates: Constraints, plus `air format --check .` and lintr on
     test-capi.R.
2. Commit 2, docs.
   - [2. Reach](../design/feature-matrix.md#2-reach): the multinom flat cell
     goes M -> S, citing the accessor.
   - [Footnotes](../design/feature-matrix.md#footnotes): [f4] is rewritten. The
     sampler is built by `dbarts()` and driven flat by run and predict, with
     K from the accessor. Counts and category offsets go through the R
     object, and an offset-carrying sampler's flat predict follows Q3.
   - [Gaps](../design/feature-matrix.md#gaps): the "Flat C reach for the
     K-forest softmax family" row is deleted.
   - docs/plans/bartcore-landing/changes.md:
     - chg-C03's entry count goes from 25 to 26, or 27 under Q2(c).
     - chg-C25 becomes five (six) entries neither consumer calls, adding
       the accessor(s).
     - A new chg-C34 row records the accessor, the K-wide layouts, the
       varcount and predict-offset behaviour on multinomial handles, the
       re-bake with the pair at 1/0, and dec-B160.
   - TODO multinomial-doors: the dbarts.h clause says run and predict are
     open (dec-B160) and creation stays the door.
   - man/dbarts-embedding.Rd: one sentence on L and the accessor in the
     "What the header reaches" paragraph.
   - No NEWS text. The 1.0-0 C API bullet introduces the header as new and
     names no family or query. Per-family layout is header and Rd detail,
     and NEWS lists what a user would act on.
   - Gates: `tools/check-doc-freshness.R`, `tools/check-rc-codoc.R`, and
     `R CMD check --as-cran` from a clean tarball (man/ touched), each on
     its own exit status.
3. Consumer runs (no dbarts commit).
   - stan4bart bartcore: install the post-slice dbarts with `--preclean`
     into a library that still holds the PRE-slice stan4bart binary. Run
     stan4bart's tinytest suite at_home with that stale binary. Expect a
     pass, which is the stale-binary claim above. Then rebuild stan4bart
     with `--preclean` and run the suite again. Expect a pass with no source
     change.
   - treatSens dbarts-1.0: the same two runs with its suite.
   - bartCause and bairrtt: nothing, since neither links dbarts (confirmed
     by `git -C <repo> grep`).
   - Any needed edit, or a failure of the stale run, is a finding that
     stops the slice.
4. Commit 4, records.
   - This plan's Landing note, pinned to the landed shas, with both consumer
     runs.
   - Status line; INDEX row status.
   - Plain commit messages throughout.

## Verification

- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e
  'tinytest::run_test_file("inst/tinytest/test-capi.R")'`: zero failures,
  new arms included. Then `tinytest::test_package("dbarts")`: zero failures.
- `cd tests/cpp && make && ./test_bartcore`: all pass.
- The equivalence trio compares print the full count of "identical draws"
  lines and no "max |z|".
- The AddressSanitizer driver over test-capi.R reports zero diagnostics on
  the slice and the overflow on the pre-slice discriminating case.
- The C99 syntax check of `inst/include/dbarts/dbarts.h` is clean.
- `sh tools/check-api-hash.sh` prints its skip line, or passes with the
  minor bump if a tag exists.

## Open questions (VD)

Q1. The accessor's name. Background: the header already has one word for
the per-observation channels a draw carries, "reported locations", used by
the per-draw callback's struct. It is 1 for every model but multinomial,
where it is K.
- (a) Name it after that word, `..._numReportedLocations`. It matches the
  callback struct, stays right if another multi-location model ever ships,
  and pairs with `..._numVariableCountForests` under Q2(c). But "reported
  locations" is jargon to a multinomial user.
- (b) Name it `..._numCategories`. A multinomial user reads it at once. But
  it answers 1 on models that have no categories, differs from the
  callback's word, and would mislead if a non-categorical multi-location
  model arrived.

Recommendation: (a), with the doc comment saying "K, the category count, on
multinomial".

Q2. The split-count channel. Background: the flat run reports each draw's
per-predictor split counts, today from one forest only. On BCF that is the
prognostic forest; on multinomial it is the first category's. R's run
reports every forest: K slabs on multinomial and 2 on BCF.
- (a) Make the width L. Multinomial gets all K, equal to R's. BCF keeps one
  slab and stays unequal to R. No second accessor is needed. But a stale
  binary that sized one slab on a multinomial handle now overflows.
- (b) Keep one slab everywhere and document it. Nothing changes, but
  multinomial reports a quantity no one wants: no category is privileged,
  unlike BCF's prognostic forest.
- (c) Make the width the sampler's own varcount forest count, read through a
  second accessor that mirrors the callback struct's
  `numVariableCountForests`. Flat then equals R on every channel, BCF
  included, under one rule. Its costs: a second new entry, a second
  forest-related name, and the same stale-binary overflow on BCF handles
  as well (no known consumer drives one).

Recommendation: (c). Locations and split-count slabs are different
quantities that agree only on multinomial, and (c) gives one parity rule
instead of a BCF exception.

Q3. A prediction offset on multinomial. Background: R's predict on a
multinomial sampler takes an offset only as a rows x K matrix added before
the softmax. If the sampler carries a category offset set from R, R's
predict also refuses to run without one, since the new rows' offset cannot
be inferred. The flat entry's offset is a bare pointer.
- (a) Refuse. A non-null offset raises, and a sampler carrying an R-set
  category offset returns 0. Such samplers then have NO flat predict at
  all, and their host must call R's predict. The upside is that it can be
  opened later without an ABI change.
- (b) Read the offset as a rows x K matrix, as R does, and raise on a null
  offset where R does. That gives full parity with R's predict, including
  offset-carrying samplers, at about fifteen lines plus tests. But the
  matrix's finiteness is unchecked, where R refuses non-finite entries
  (they give NaN probabilities). Its shape is unchecked too, as `out`'s is.
  It is also the only category-offset channel the flat API has.

Recommendation: (b). Under (a) the samplers most likely to sit inside a
host's Gibbs step, those with an offset, cannot predict from C.

Q4. Bundle the monotone run gaps (TODO monotone-count-host-interrupt)?
Background: the flat run neither polls for a user interrupt nor reports the
slow leaf-order count tally, so a host's monotone fit cannot be interrupted
mid-count and never warns. R's run does both.
- (a) Bundle now.
  - The interrupt poll is behaviour only, with no ABI change: the flat run
    passes R's poll and raises "sampler run interrupted" after the engine
    unwinds, as R's run does.
  - The tally needs a channel. One option is an R warning raised by the
    entry itself, with no ABI change. The other is an appended
    `dbarts_results` field (four values per run, plus `DBARTS_RESULTS_INIT`),
    which moves the layout fold and shares this re-bake.
  - Cost: about 60 lines plus tests, and a review on stan4bart's side. An
    interrupt then unwinds through stan4bart's sweep loop, and a warning can
    fire once per call in its Gibbs loop.
- (b) Leave it to its own item. This slice stays multinomial-only. Only a
  field-shaped tally would share this re-bake, and before 1.0-0 a re-bake
  costs one test literal and one consumer rebuild.

Recommendation: (b). Most of it needs no ABI event at all, and the part
that matters is a stan4bart behaviour change that deserves its own review.
