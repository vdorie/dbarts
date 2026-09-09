# pure-c-header

agent: opus (S1 header, bridge and C entry file; S2 the two consumer
  ports, worked from their own checkouts); sonnet (S3 the CI assertion,
  NEWS, Rd, records). Serialized: S1 lands before S2 starts; S3 may run
  beside S2.
rng: neutral on every slice. No engine file changes; the bridge copies
  the same values into the same buffers. The three bitwise equivalence
  baselines (gaussian, BCF, multinomial) expect IDENTICAL and any
  deviation is a leak.
window: first of the pre-release arcs. stan4bart and treatSens build
  against the result, so nothing else touches `dbarts.h` while this is
  open; the front-door arc (docs/plans/front-door.md) may run in parallel
  because it is R-only.
budget: ~-500 header, ~-700 C entry file, ~+120 bridge, ~+60 tools/CI,
  ~-900 +150 test consumer, ~+40 R (tinytest), ~60 docs; stan4bart ~70
  lines over four files, treatSens ~40 lines over three.

Decisions: dec-B84 (pure C, the four R-object entries removed), dec-B86
(trim to called entries plus sizing queries), dec-B87 (copy-on-set),
dec-B111 (exact-ABI flag off, hash-to-minor-bump CI assertion), all in
[docs/decisions.md](../decisions.md). dec-A20's null-pointer refusal
channel goes with the creation entry.

## Goal

The shipped header's prototype view compiles as C with no R header
included (the stub view still needs R's error and dynload headers, since
it resolves symbols through R). A compiled
consumer creates the sampler through dbarts's R interface, reads the
handle out of the sampler object's external pointer, and drives it with
the entries below. Every setter copies into a buffer the sampler owns
and allocated at creation. A change to the header's hash without a
version bump fails CI once 1.0-0 is tagged. stan4bart and treatSens
build and pass their tests against the new header on their compat
branches.

## Context

- The header today: 48 entries in
  [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h), the
  compile-time token
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) checked by
  `static_assert` in [`dbarts_apiToken`](../../src/C_interface.cpp), the
  stub view under [`DBARTS_USE_STUBS`](../../inst/include/dbarts/dbarts.h)
  and its opt-in
  [`DBARTS_REQUIRE_EXACT_ABI`](../../inst/include/dbarts/dbarts.h). The
  header includes Rinternals.h only for the four SEXP entries.
- The C handle wraps the bridge's holder:
  [`dbarts_sampler_t`](../../src/C_interface.cpp) holds a
  [`BartcoreHolder`](../../src/R_interface_bartcore_common.hpp) pointer,
  a preserved data SEXP and the callback pair. The R sampler holds the
  same holder behind
  [`dbartsSampler$getPointer`](../../man/dbartsSampler-class.Rd), created
  by [`createExternalHolder`](../../src/R_interface_bartcore.cpp) with
  [`holderFinalizer`](../../src/R_interface_bartcore.cpp) and read back by
  [`holderFromExpression`](../../src/R_interface_bartcore.cpp). Both
  creation routes reach
  [`createHolder`](../../src/R_interface_bartcore.cpp).
- Ownership today: the engine borrows. Its response models keep
  `const double*` members
  ([`ResponseModel`](../../src/bartcore/model.hpp),
  [`Chain::setWeights`](../../src/bartcore/chain.hpp),
  [`ColumnStore`](../../src/bartcore/data.hpp) for the test offset).
  The R bridge keeps the engine's pointers alive by pinning the R
  vectors (retired: [`PROT_RESPONSE`](../../src/R_interface_bartcore.cpp) and its
  siblings) on the main path and by moving a fresh vector into the
  holder's `ownedResponse`, `ownedOffset`, `ownedWeights`,
  `ownedTestOffset` ([`BartcoreHolder`](../../src/R_interface_bartcore_common.hpp))
  on the data-handle path; the C setters hand the caller's pointer
  straight through. Predictors already copy (data-ownership design,
  [The adopted design](../design/data-ownership.md#the-adopted-design)).
- Consumers, from their compat branches (stan4bart `bartcore`, treatSens
  `dbarts-1.0`): stan4bart calls create at two sites plus storeState,
  setState, getTrees, printTrees, getLatents, predict, run,
  sampleTreesFromPrior, setOffset, setSigma, setVerbose, setTreeStorage,
  destroy and seven sizing queries; treatSens calls create at two sites
  as well,
  plus run, setResponse, setOffset, setSigma, setNumThreads, setVerbose,
  destroy and the version pair. Both define `DBARTS_USE_STUBS` and
  `DBARTS_REQUIRE_EXACT_ABI` (stan4bart in its Makevars, treatSens per
  file). The shipped test consumer
  ["consumer.c"](../../inst/tinytest/test-capi.R) exercises 40
  entries including the whole multi-forest block, and five handshake
  arms (a consumer compiled with a wrong hash, with the wrong hash under
  the exact-ABI flag, with a wrong major, with a wrong minor, and with
  the exact-ABI flag alone).
- State format: [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)
  and the floor `minReadableStateFormatVersion` are read only by the
  bridge; a state with no `formatVersion` attribute reads as 0 and is
  refused there. dec-B109's message rides the front-door arc, not this
  one.
- Prior records this supersedes: public-surface design section
  [6. C API and callbacks](../design/public-surface.md#6-c-api-and-callbacks)
  ("SEXPs at exactly two boundaries"), the entry count and the
  getTrees reasoning in [capi-shape.md](capi-shape.md), the reader-name set
  in [dbarts-h-freeze.md](dbarts-h-freeze.md), and the "do not reopen" list
  in [dbarts-h-reshape.md](archive/dbarts-h-reshape.md). Each gets one
  Status-line amendment pointing here.

## Decision

No fork remains open. VD ruled on destroy (2026-09-08, recorded as
dec-B112): R-hosted consumers let the collector free the sampler, and a
C consumer keeps a path to invalidate its object. So
`dbarts_sampler_destroy` stays: it releases the engine sampler behind
the handle and marks the R object's pointer dead; the R object's
existing dead-pointer path then re-creates from a stored state or
refuses; a second destroy is a no-op.

Fixed by the decisions, not open: the handle is the address stored in
the sampler object's external pointer, read by the consumer with
`R_ExternalPtrAddr` in its own code, so the header names no R type and
no new entry is needed (the holder becomes the handle; the callback pair
and the preserved data SEXP move into it, a bridge-private change, and
the retained entries never touch the data SEXP); the trimmed list is
exactly the entries a known consumer calls plus the queries a host needs
to size the results struct (below); the exact-ABI flag comes off both
consumers; the hash assertion is CI, not a runtime check.

## Constraints

- Gates: neutral class. tests/cpp component tests; full tinytest suite
  including a rewritten test-capi.R; equivalence compare identical on
  all 50 gaussian, 12 BCF and 11 multinomial scenarios; sanitizers on
  the C entry file (test-capi.R has a sanitizers budget); R CMD check
  --as-cran; revdep-smoke by dispatch against both consumer branches.
- The header compiles under `gcc -std=c99 -pedantic -x c` and
  `g++ -std=c++11 -x c++` with no R include path; add that as a
  tests/cpp target.
- The hash and the version constants: pre-1.0-0 the hash moves and the
  constants stay at 1.0, per the header's own rule; the CI assertion is
  armed only when a release tag exists (step 8).
- Out of scope: any engine change; an opaque C state blob (dec-B84 adds
  it only when a consumer needs C-side store and restore, and none
  does); the plain-C specification stage (dec-B85, post-release); the
  memory audit (docs/plans/memory-footprint-audit.md) beyond stating
  what this arc allocates.

## Steps

S1, header and bridge:

1. Remove `dbarts_sampler_create`, `dbarts_sampler_getTrees`,
   `dbarts_sampler_storeState` and `dbarts_sampler_setState` from the
   list, the entry file
   and the registration table
   ([`DBARTS_API_REGISTER`](../../src/R_interface.cpp)). Delete the
   Rinternals.h include and its R_NO_REMAP dance from the prototype
   view; the stub view keeps `R_ext/Rdynload.h` and gains
   `R_ext/Error.h` for the `Rf_error` its handshake raises.
2. Trim the list to: the version triple; `run`, `sampleTreesFromPrior`;
   `setResponse`, `setOffset`, `setSigma`; `getLatents`; `predict` with
   `dbarts_predictor_source` and `dbarts_dense_predictor_source`;
   `setTreeStorage`, `printTrees`; `setNumThreads`, `setVerbose`; the
   nine sizing queries (`numObservations`, `numPredictors`,
   `numTestObservations`, `numChains`, `numTrees`, `numSavedSamples`,
   `kIsSampled`, `usesDart`, `family`; the last three say which
   `dbarts_results` fields the run will fill, which is how a host sizes
   that struct); `destroy` with the semantics above. Twenty-four
   entries. Delete
   with their entries: `dbarts_forest_calibration` and its init macro,
   `dbarts_sampler_callback` and `DBARTS_SAMPLER_CALLBACK_PARAMS`,
   `dbarts_drawLatents`, `dbarts_workingResponse`. Keep
   `dbarts_results` whole, `dbarts_column_type`, `dbarts_leaf_model` and
   `dbarts_family`.
3. Make the holder the handle: define `dbarts_sampler_t` as the holder
   type in the common bridge header so `R_ExternalPtrAddr` on the sampler
   object's pointer field is the handle. Write the contract as one
   header paragraph: obtained from the R sampler object, valid until
   that object is garbage collected or re-creates its pointer (a restore
   from a stored state), so a consumer re-reads it after any R-side
   restore and keeps the R object reachable while it holds the handle;
   `destroy` releases the engine sampler early and leaves the R object
   in its dead-pointer state.
4. Copy-on-set: `dbarts_sampler_setResponse` and `setOffset` copy into
   the holder's owned vectors and install those. Allocate `ownedResponse`,
   `ownedOffset`, `ownedWeights` and `ownedTestOffset` to length n (test
   offset to n.test) inside [`createHolder`](../../src/R_interface_bartcore.cpp)
   and make both R-path routes (the pin and the move) copies into those
   buffers, so no R vector is retained. The no-allocation-after-creation
   guarantee is the C handle's: the R conduits that change the
   observation count ([`bartcore_setData`](../../src/R_interface_bartcore.cpp))
   or the test count (`bartcore_setTestPredictor`) resize the buffers,
   and the header says so. Cost: three n-length doubles and one of
   test length per sampler, which the memory audit prices against the
   n-times-trees index buffers.
   One contract sentence in the header, above the setters.
5. Re-bake both tokens (the header literal and the signature token in
   [`dbarts_apiSignatureToken`](../../src/C_interface.cpp)); the build
   fails until they agree.
6. Tests: rewrite the shipped consumer to create through R (an R helper
   in the test builds a `dbarts` sampler and passes the object to the
   consumer's `.Call`), keep the five handshake arms, add a copy-on-set
   probe (set a response from a buffer, overwrite the buffer, run, and
   check the draws match a run from the unmodified copy), and a
   handle-after-restore probe (storeState, drop the pointer, getPointer,
   re-read the handle, run), and a destroy probe (destroy, then the R
   object's method re-creates from state or refuses, and a second
   destroy is a no-op). Delete the tests for removed entries. The
   R features those entries reached stay covered by their R tests.
   Add the two-compiler header-compiles target to tests/cpp.
7. Docs with the code: the header's own prose; `dbarts-embedding.Rd`
   (the C section names the handle route); NEWS 1.0-0 UPGRADING entry;
   Status-line amendments to the four prior records named in Context;
   the feature matrix's [2. Reach](../design/feature-matrix.md#2-reach)
   column for the flat C header, whose bcf row cites the creation entry.

S3, CI assertion (dec-B111):

8. `tools/check-api-hash.sh`: reads `DBARTS_C_API_HASH`, `_MAJOR` and
   `_MINOR` from the working tree and from the newest tag matching
   `v1.*`; fails when the hash differs and the pair does not; prints
   "no release tag, skipped" and exits 0 when no such tag exists. One
   step in check-standard.yaml (it already checks out with history).
   Prove it discriminates by running it against a scratch tree with the
   literal edited.
9. Record the tag-naming convention the check depends on in
   [docs/plans/README.md](README.md) under CI.

S2, consumer ports (each from its own checkout, own private library
built from the S1 tip):

10. stan4bart bartcore branch: the two creation sites build the
    sampler in R (stan4bart's R side already constructs the control,
    model and data objects; it calls `dbarts::dbarts` on them and hands
    the sampler object down) and the C++ reads the handle; the restore
    path uses the R object's `setState`; `getTrees` becomes the R
    method; keep its destroy calls; drop `DBARTS_REQUIRE_EXACT_ABI`
    from both Makevars. Run stan4bart's tinytest suite and its
    equivalence check against dbarts if it carries one.
11. treatSens dbarts-1.0 branch: the two creation sites likewise; drop
    the per-file exact-ABI define in three files. Run its tests.
12. revdep-smoke dispatch on the S1 tip, all three packages green
    (bartCause has no compiled boundary and needs no port).

## Verification

```
R CMD INSTALL --preclean -l <lib> .
cd tests/cpp && make && ./test_bartcore            # includes header-compiles
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-fbff1989.rds
  # 50 "identical draws (same RNG stream)", no "max |z|" line; same for bcf (12) and multinomial (11)
gcc -std=c99 -pedantic -Wall -x c -fsyntax-only inst/include/dbarts/dbarts.h   # no R include path
sh tools/check-api-hash.sh                          # "no release tag, skipped" today
gh workflow run revdep-smoke.yaml --ref <tip>       # three green
```

Expected: header size roughly halved; `grep -c SEXP inst/include/dbarts/dbarts.h`
is 0; both consumers build with no `DBARTS_REQUIRE_EXACT_ABI` and pass.
