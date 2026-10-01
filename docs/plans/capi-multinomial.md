# capi-multinomial

Status: PROPOSED 2026-10-01. Ruling: dec-B160 in [decisions.md](../decisions.md).

agent: opus (commit 1, the header, the C entry file and the test consumer);
  sonnet (commit 2, docs; step 3, the consumer rebuilds; commit 4, records).
  Serial: 1, 2, 3, 4.
rng: NEUTRAL. No engine file changes. The run entry already hands the
  engine the K-wide location stride; what moves is the contract, one
  accessor, one refusal and (Q2) the varcount width the flat caller declares.
  Every sampler with one location runs byte for byte what it runs today.
window: nothing else edits `inst/include/dbarts/dbarts.h` while this is
  open. Lands before the 1.0-0 tag, or bumps the minor version (Constraints).
budget: header ~+50 -15; C entry file ~+25 -12; test consumer ~+170; test
  file ~+130; docs ~+25 -10; records ~+40. About 450 lines; plan on 700-900.

## Goal

A multinomial sampler built in R runs and predicts through the flat C API.
The header gains `dbarts_sampler_numReportedLocations`, which reports K
before any run, and it documents the K-wide layouts of
[`dbarts_sampler_run`](../../inst/include/dbarts/dbarts.h) and
[`dbarts_sampler_predict`](../../inst/include/dbarts/dbarts.h). A caller who sizes its buffers by
the documented layouts gets bitwise what R's `$run` and `$predict` return.
Out of scope: a flat creation path, and any flat count or category-offset
setter. Both stay with the R object (TODO multinomial-doors).

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
whatever the caller sized, so train and test come back n x K x S x C into
buffers documented as n x S x C. That is a heap overflow. The verification probe
(n = 30, K = 3, S = 4, C = 2) wrote 720 of a documented 240 train entries and
crashed later in unrelated R code. [`dbarts_sampler_predict`](../../src/C_interface.cpp) likewise
writes nTest x K x S x C. A non-null offset there is added slab by slab at
stride nTest: an additive shift on probability channels, landing on the wrong
entries. varcount is safe: the run leaves `numVariableCountForests` at 1, so
the engine writes only the first category's slab. No current entry reports K
before a run. [`dbarts_sampler_numTrees`](../../inst/include/dbarts/dbarts.h) raises past the last
forest, and the [`dbarts_draw`](../../inst/include/dbarts/dbarts.h) field `numReportedLocations` arrives
only inside a run, too late to size buffers.

What R returns (probed against the installed build: n = 30, K = 3, nTest =
10, four draws, two chains). [`bartcore_run`](../../src/R_interface_bartcore.cpp) gives:

| channel | R shape | content |
|---|---|---|
| sigma | S x C | the pinned 1 |
| train | n x K x S x C | softmax probabilities; each row of K sums to 1 |
| test | nTest x K x S x C | the same, under any category test offset installed from R |
| varcount | p x K x S x C | slab k is category k's forest |
| k, varprobs | NULL | a k hyperprior and DART are refused on multinomial |

[`predictFromSource`](../../src/R_interface_bartcore.cpp) returns nTest x K x S x C with saved
trees and nTest x K x 1 x C without. Its offset must be an nTest x K
per-category matrix entering before the softmax. A flat vector is refused,
because after the blend it moves values off the simplex and before it a
common shift is the softmax's null direction. Category k is column k of the
sampler's count matrix.

The other entries on a multinomial handle already behave, and nothing here
changes them. setResponse and setOffset return 0
([`responseConduitIsFixed`](../../src/R_interface_bartcore.cpp)), as do setSigma and getLatents.
printTrees and numTrees take forests 0..K-1. The per-draw callback sees
`numReportedLocations` = K, with n x K train.

Consumers today. stan4bart (branch bartcore at 9a3be93) calls run, predict,
setOffset, setSigma, getLatents, setTreeStorage, sampleTreesFromPrior,
printTrees, destroy, the size queries and the version pair. It builds only
gaussian and probit samplers, sizes run's train as n, and passes a
non-null offsetTest to predict. treatSens is also a C consumer, contrary to
the brief: on dbarts-1.0 at 7cc6a0f it carries LinkingTo dbarts and calls run,
setResponse, setOffset, setSigma, setNumThreads, setVerbose and destroy, on
gaussian and probit samplers only. Both use `DBARTS_USE_STUBS`, and neither
defines [`DBARTS_REQUIRE_EXACT_ABI`](../../inst/include/dbarts/dbarts.h) (dropped under dec-B111).
Neither pins the hash. bartCause (dbarts-1.0) and bairrtt have no LinkingTo
and no flat calls. No shipped consumer reaches a multinomial handle.

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
field that already carries L. It is a location count, not a forest count:
BCF has two forests and answers 1 (Q1 weighs the multinomial-only name).

An existing entry cannot carry K:
- `dbarts_sampler_family` already says multinomial but not K. Packing K into
  its int would make a VALUE entry mean two things.
- `dbarts_sampler_numTrees` raises past the last forest, so probing for K means
  catching a longjmp. A forest count is also not L (BCF).
- A shape struct filled by one query would subsume the seven size queries, but
  it adds a third size-first struct to the layout fold for one number. It
  would also duplicate entries that stay.
- The draw callback carries L only during a run.

### Layouts (header text, per field)

L = `dbarts_sampler_numReportedLocations`. Observation is fastest, then
location, then draw, then chain. That is R's array order, so `as.vector` of
R's result equals the flat buffer.

- `train`: numObservations x L x numSamples x numChains. On multinomial each
  draw's L columns are category probabilities, rows summing to 1. That is the
  softmax of the K forests plus any category offset installed from R.
- `test`: numTestObservations x L x numSamples x numChains. A category test
  offset installed from R applies here, never to predict.
- `varcount`: numPredictors x L x numSamples x numChains, slab k category k's
  forest (Q2). The run sets `numVariableCountForests` to L. That is 1 for
  every other sampler, BCF included, so their single prognostic slab is
  unchanged.
- `sigma`: numSamples x numChains, unchanged. On multinomial it is the pinned 1.
- `logLikelihood`: numObservations x numSamples x numChains, unchanged and
  NaN-filled on multinomial (its existing sentence stays).
- `k`, `varprobs`, `dispersion`, `residualDf`: unchanged and left untouched on
  multinomial. No multinomial sampler carries any of them.
- predict `out`: xTest->numRows x L x S' x numChains, where S' is
  `numSavedSamples` with tree storage and 1 without. The scale is
  probability on multinomial and otherwise as today.
- predict `offsetTest`: must be null when L > 1 (Q3). A non-null one raises
  before anything is written. It is an argument error, not a capability 0,
  since a null offset works.

Header edits beyond the fields:
- The run and predict paragraphs state the layouts above.
- The common-contracts bullet "Result and prediction layouts put samples and
  then chains in trailing dimensions" gains "after L locations".
- The forest-index bullet's "states its count nowhere here" gains an exception:
  a multinomial sampler's K forests are its L locations.
- The `dbarts_results` paragraph about varcount and "this struct declares no
  forest count" is rewritten to the L rule.
- `dbarts_draw`'s L line points at the accessor.

In [C_interface.cpp](../../src/C_interface.cpp):
- The stale stride comment in [`dbarts_sampler_run`](../../src/C_interface.cpp) ("1 for every
  dbarts.h-created sampler") is replaced.
- The comment in [`dbarts_sampler_family`](../../src/C_interface.cpp) ("no entry here builds one yet")
  is restated for R-built handles.

### Wrong sizes

No output size is declared anywhere: `dbarts_results` carries bare pointers
and predict's `out` is a bare pointer. So there is nothing to check, and a
short buffer remains the caller's crash under the header's "Validation is
deliberately partial" rule. The declared sizes that exist, a predictor
source's numRows and numColumns, are checked as today. Considered and not
taken: a caller-declared L field in `dbarts_results` that the run would
refuse on mismatch. That covers run only. Predict would need a new parameter,
a signature change that forces consumer source edits, where the ruling asks
only for a rebuild.

### The ABI event

- Appending to the list moves the signature token and the full token. Both
  literals are re-baked in the same commit: `DBARTS_C_API_HASH` in the header
  and the [`dbarts_apiSignatureToken`](../../src/C_interface.cpp) assert. Use the probe procedure in
  [5. Hash re-bake](dbarts-h-freeze.md#5-hash-re-bake). No struct, enumerator or layout moves, so
  the fold in [`dbarts_apiToken`](../../src/C_interface.cpp) is otherwise unchanged.
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
- What a stale consumer binary sees. With stubs and no exact-ABI flag
  (stan4bart, treatSens): the version handshake passes, and every entry it
  calls resolves by name with an unchanged signature and struct layout. It
  loads and runs as before, bitwise on its gaussian and probit samplers. With
  `DBARTS_REQUIRE_EXACT_ABI`: the first stub call raises "dbarts C ABI
  mismatch ... rebuild". A consumer built against the new header with an
  older dbarts installed fails only when it first calls the new entry ("not
  provided by package 'dbarts'"). The pre-1.0 version pair cannot catch
  that; after 1.0-0 the minor bump does. The "rebuild once" in the ruling
  is therefore lockstep hygiene: no consumer source changes, and a skipped
  rebuild breaks nothing a consumer does today.

## Constraints

- Gates (neutral class): full tinytest on a `--preclean` private-library
  install, and tests/cpp. The equivalence trio, labelled a formality (no
  harness drives a flat entry), must come back IDENTICAL against the
  baselines named in benchmarks/baselines/MANIFEST. Also the R-loaded
  AddressSanitizer run of test-capi.R per the plans README (the defect is a
  heap overflow and the new code is reachable only through `.Call`). Also
  a C99 `-pedantic -fsyntax-only` compile of the header's prototype view.
- Frozen: no existing signature, struct layout or enumerator changes. No
  engine file is touched.
- Out of scope: flat creation; flat setters for counts, the category offset
  or the category test offset; a category offset on flat predict (Q3); a
  flat log-likelihood for multinomial.

## Steps

1. Commit 1, the ABI event (code, header, tests together; the hash pin makes
   them inseparable).
   - Header: the accessor, the layout text above, the bullet edits.
   - [C_interface.cpp](../../src/C_interface.cpp): the accessor body.
     `engineResults.numVariableCountForests = shape.numReportedLocations` in
     run (Q2). The offsetTest refusal in predict, raised beside the
     empty-store refusal and before the replay. The two comments.
   - The re-bake of both literals.
   - [consumer.c](../../inst/tinytest/capi/consumer.c):
     - `capi_run_canaried(ptr, burn, samples)` sizes every buffer from the
       accessors, adds a tail of the same length filled with a fixed NaN
       payload (0xDEADBEEF for varcount), runs, and returns the bodies, L,
       a per-channel "tail intact" flag (memcmp) and a per-channel "body
       fully written" flag (no canary word left).
     - `capi_predict_canaried(ptr, x, offset)` is the same for predict.
     - [`capi_dims`](../../inst/tinytest/capi/consumer.c) gains L.
   - [test-capi.R](../../inst/tinytest/test-capi.R):
     - Hash pins.
     - A gaussian arm: two identically seeded two-chain samplers with test
       rows. Flat run on one and `$run` on the other, `expect_identical` on
       `as.vector` of every channel; canaries intact and bodies written;
       L = 1.
     - The same for a two-chain multinomial (K = 3), plus row sums of 1 and
       logLikelihood all NaN.
     - Predict on each, against `$predict` on the same sampler after the
       run, with and without tree storage (nTest x K x 1 x C).
     - Flat predict with a non-null offset on multinomial raises. On
       gaussian it still adds.
     - L on BCF is 1 with two forests.
     - setResponse and setOffset return 0 on the multinomial handle (the
       arms capi-shape.md recorded as unreachable are now reached).
     - The callback's `num.reported.locations` equals the accessor on
       multinomial.
   - Run the gaussian parity check first. If flat and R runs differ there,
     stop: that is a separate finding, not this slice's.
   - Mutation proofs, each reverted and the file `touch`ed:
     - The accessor returning 1 fails the multinomial arm.
     - Dropping the varcount line fails varcount parity.
     - Dropping the refusal fails the offset expectation.
   - Gates: Constraints, plus `air format --check .` and lintr on
     test-capi.R.
2. Commit 2, docs.
   - [2. Reach](../design/feature-matrix.md#2-reach): the multinom flat cell
     goes M -> S, citing the accessor.
   - [Footnotes](../design/feature-matrix.md#footnotes): [f4] is rewritten (built
     by `dbarts()`, driven flat by run and predict, K from the accessor, counts
     and category offsets through the R object).
   - [Gaps](../design/feature-matrix.md#gaps): the "Flat C reach for the K-forest
     softmax family" row is deleted.
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
3. Consumer rebuilds (no dbarts commit).
   - stan4bart bartcore: `R CMD INSTALL --preclean` dbarts at the slice tip
     into its library, rebuild stan4bart with `--preclean`, and run its
     tinytest suite at_home. Expect a pass with no source change.
   - treatSens dbarts-1.0: the same with its suite.
   - bartCause and bairrtt: nothing (no LinkingTo, confirmed by
     `git -C <repo> grep`).
   - Any needed edit is a finding that stops the slice.
4. Commit 4, records.
   - This plan's Landing note, pinned to the landed shas, with the consumer
     results.
   - Status line; INDEX row status.
   - Plain commit messages throughout.

## Verification

- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e
  'tinytest::run_test_file("inst/tinytest/test-capi.R")'`: zero failures, new
  arms included. Then `tinytest::test_package("dbarts")`: zero failures.
- `cd tests/cpp && make && ./test_bartcore`: all pass.
- The equivalence trio compares print the full count of "identical draws"
  lines and no "max |z|".
- The AddressSanitizer driver over test-capi.R reports zero diagnostics.
  Against the pre-slice header sizes it must trip; run that once as the
  discriminating case.
- The C99 syntax check of `inst/include/dbarts/dbarts.h` is clean.
- `sh tools/check-api-hash.sh` prints its skip line (or passes with the minor
  bump if a tag exists).

## Open questions (VD)

Q1. The accessor's name. Background: the header already has one word for
the per-observation channels a draw carries, "reported locations", used by
the per-draw callback's struct. It is 1 for every model but multinomial,
where it is K.
- (a) Name it after that word, as `..._numReportedLocations`. It matches the
  callback struct and stays right if another multi-location model ever
  ships. "Reported locations" is jargon to a multinomial user.
- (b) Name it `..._numCategories`. A multinomial user reads it at once. But
  it answers 1 on models that have no categories, it differs from the
  callback's word, and it would mislead if a non-categorical multi-location
  model arrived.
Recommendation: (a), with the doc comment saying "K, the category count, on
multinomial".

Q2. The split-count channel on multinomial. Background: the flat run reports
each draw's per-predictor split counts. Today it reports one forest's: on
BCF the prognostic forest, which is meaningful alone, and on multinomial
the first category's.
- (a) Report all K, one per category, laid out like train. Flat equals R's
  run on every channel the struct carries. The width follows the same L the
  caller already reads, so no further accessor is needed. BCF is unchanged.
- (b) Keep the single first-category slab and document it. That changes
  nothing, but it reports a quantity no one wants: no category is
  privileged, unlike BCF's prognostic forest.
Recommendation: (a).

Q3. A prediction offset on multinomial. Background: R's predict on a
multinomial sampler takes an offset only as a rows x K matrix added before
the softmax, and refuses a plain vector. The flat entry's offset is a bare
pointer with no shape, so it cannot tell the two apart. Separately, the flat
API has no way to set the training-side category offset; that stays in R.
- (a) Refuse a non-null offset on multinomial with an error. This is safe and
  can be opened later without any ABI change. A C host wanting a category
  offset calls R's predict.
- (b) Read it as a rows x K category matrix, passed through to the engine as
  R does. That gives full parity with R's predict at about ten lines. But
  one pointer would then mean an additive shift on every other model and a
  pre-softmax per-category shift here, decided by the sampler rather than
  by anything the caller wrote. It would also be the only category-offset
  channel the flat API has.
Recommendation: (a).
