# dbarts-h-reshape

agent: S0 opus (engine printer widening + bridge helper promotion, no dbarts.h
  change). S1 opus (the one dbarts.h edit, the C_interface bodies, the hash
  re-bake). S2 sonnet (consumer rebuilds, docs, records). Serialized: one
  implementer, each slice lands before the next starts.
rng: NEUTRAL, every slice. Nothing here touches a proposal, a suffstat, a cut
  builder or an RNG consumer. S0 widens two diagnostic printers and moves three
  refusal helpers between namespaces; S1 is dbarts.h + C_interface.cpp only.
  The trio is a FORMALITY and must be LABELLED one in the landing note: no
  equivalence harness drives a single dbarts.h entry point (measured -
  `benchmarks/` has zero `dbarts_sampler_` hits). `inst/tinytest/test-capi.R` +
  `inst/tinytest/capi/consumer.c` are the load-bearing neutrality gate (F1).
window: the LAST breaking change before the 1.0-0 freeze (TODO, "the
  typed-ingestion dbarts.h reshape (VD 2026-08-08) lands BEFORE the freeze").
  Runs AFTER bcf-public-surface (whose S3 re-bakes FIRST and whose five entries
  this arc adopts verbatim) and AFTER multiforest-predictor-mutation S0-S4
  (which adds no X-list entry and moves no hash - stated at four places in that
  plan). This arc's re-bake is the SECOND AND LAST of the window. NO version
  constant moves in either re-bake (VD 2026-08-10, binding decision 8).
  Everything is breakable; the sister packages migrate in lockstep at the
  freeze, once, against the final header (VD 2026-08-10).
budget: S0 ~15 engine + ~10 facade + ~65 bridge/common-header + ~45 tests/cpp.
  S1 ~150 header + ~230 C_interface + ~200 consumer.c + ~150 test-capi.R.
  S2 ~100 docs/records + 5 out-of-repo call-site edits (all stan4bart), 2
  optional handshake edits, and 4 rebuilds. Total in-repo ~965.
tip: every anchor below re-read against the live tree at 2e50cf1 (branch
  bartcore, local == origin). bcf-public-surface S0 has landed (994c161); its
  S1 and S3 have NOT. Design artifacts (memo, blind refuting critique, two
  probes) sit at `.claude/dbarts-h-reshape-design/`, which is GIT-IGNORED - this
  plan is the only tracked record and carries the load-bearing facts rather than
  pointing at them.

## Goal

The flat C surface stops being dense-only and forest-blind. The four predictor
entries take the same borrowed, self-describing source view the R ingestion
surface takes, so a C consumer mutates and predicts from CSC storage without
densifying and without a container class; the tree-facing queries take a forest
index, so they are unambiguous on the BCF sampler bcf-public-surface S1 makes
flat-creatable. Every predictor argument declares its own width and its own CSC
column count, which closes a live out-of-bounds read. The version handshake
gains an exact signature token and loses the packed integer that cannot express
one. ONE hash re-bake, NO version bump, ONE consumer migration.

## Binding decisions inherited (do not reopen)

1. **Reshape-by-replacement.** VD 2026-08-08 (TODO, typed-ingestion door
   status): the mutation entries are RE-SIGNED in place, not shadowed by
   parallel names. The `dbarts_results` precedent (c-api-growth.md, Design
   Part 1) is the sanctioned mode - a parallel run2/results2 duality was
   proposed and superseded the same day, once it was clear nothing on CRAN
   links the unreleased header.
2. **Consumer compatibility is a cost to ENUMERATE, never a design
   constraint** (VD 2026-08-10). Enumerated in "Migration costs" below.
3. **bcf-public-surface S3's entries are adopted VERBATIM, with ONE carve-out
   (amended 2026-08-10 by multiforest-extension-surface, which adjudicated a
   conflict between this item and bcf-public-surface's AMENDED Open decision
   block).** The widened `dbarts_sampler_setResponse(sampler, y,
   updateScale)`, `numForests` and `forestFits` are adopted verbatim and
   re-signed never; no guard S3 relaxed is re-relaxed. CARVE-OUT: the two
   BCF-SPECIFIC names `setTreatment` and `bcfGlue`, and the creation contract
   documented at `dbarts.h:348-357` plus the ownership sentence at `:43`, MAY
   be re-signed to engine vocabulary in this arc's S1, because
   bcf-public-surface's later VD amendment authorizes exactly that ("The flat
   C names S3 shipped (`setTreatment`, `bcfGlue`) are likewise renameable at
   the queued dbarts.h reshape re-bake if the settled surface uses engine
   vocabulary") and because this arc's S1 is the window's LAST re-bake, so
   declining costs a post-release MINOR bump plus a lockstep bump of
   stan4bart's floor - the exact cost resolved question 2 refused to pay for
   `setForestWeights`. The re-sign is CONDITIONAL: it happens only if
   multiforest-extension-surface's naming fork is answered before this slice
   starts, and its content is that plan's M3. Otherwise the shipped names go
   to release unchanged.

   **Fork-3 RESOLVED (VD 2026-08-11, multiforest-extension-surface).** The
   condition above is met before this slice starts. Rename content:
   `setTreatment` -> `setForestBasis(sampler, forest, basis, numColumns)`;
   `bcfGlue` -> `numForestAmplitudes` plus `forestAmplitudes` (S1 item 5b
   carries the mechanics, unchanged by this note). `setForestWeights` KEEPS
   its name - the precision channel and the mean channel stay two entries
   forever, and "weights" already means precisions in dbarts. The re-sign in
   this item is therefore no longer conditional; it is decided, and S1 item
   5b is unconditional content rather than a contingency.
4. **The source semantics are already signed off** - cheap-uniformity.md
   "Source-shaped mutation surface", rules 1-5. This plan spells them in C and
   adds nothing to them.
5. **The house return convention is 1 = accepted, 0 = refused.** Measured, not
   asserted: `dbarts_sampler_setPredictor` ends `return result ==
   bartcore::PredictorUpdateResult::accepted ? 1 : 0` (C_interface.cpp:285).
   bcf-public-surface S3 item 3 states the same ("1 on success, 0 on refusal -
   the shipped convention, not the inverse"). The `size_t` probes carry no error
   channel, so an out-of-range argument to one must `Rf_error`.
6. **Occupancy and veto semantics are untouched.** This arc changes no
   acceptance criterion; multiforest-predictor-mutation owns those.
7. **The flat BCF test-surface gap is NOT this arc's.** Orchestrator sequencing
   call, recorded: it belongs to bcf-public-surface S3, whose amendment text
   this plan carries verbatim under "TODO edits at landing".
8. **NO version increment, this arc or bcf-public-surface S3.** VD 2026-08-10,
   verbatim: *"No need to increment versions - no bartcore version has been
   released."* Application: no bartcore release has ever shipped, so
   `DBARTS_C_API_MAJOR` / `DBARTS_C_API_MINOR` have never been a compatibility
   contract with anyone, and whatever constants the header carries at the FIRST
   release simply ARE the initial contract. `DBARTS_C_API_MAJOR` stays 1 and
   `DBARTS_C_API_MINOR` stays 0 through both re-bakes of this window. Do not
   carry a "the minor floor moved to 1" claim forward from
   bcf-public-surface. **Consequence, load-bearing for this plan's design:** the
   `apiVersion`/`apiHash` machinery's job in the pre-release window is LOCKSTEP
   REBUILD DETECTION for the in-repo sister packages, not release
   compatibility - the hash is the only runtime signal that can move, so it is
   what must be designed for the job. The post-release rule is unchanged and
   still stands (TODO, "on any post-release dbarts.h change: bump
   DBARTS_C_API_MINOR (or MAJOR if incompatible) alongside the hash re-bake").

## Adjudication: corrected facts this plan is built on

The design memo was attacked by a blind refuting critique (verdict STANDS WITH
AMENDMENTS, seven blocking findings). Every finding was re-verified here against
2e50cf1; three were corrected further, in both directions. What survives:

- **A1 UPHELD, and independently reproduced. stan4bart uses 23 flat entries, not
  22**, across 38 call expressions, cross-checked by 23 emitted
  `dbarts_stub_fn` block-scope statics in `src/stan4bart.so` (unused
  `static inline` stubs are never emitted, so the emitted set IS the used set).
  The "22" collapsed `dbarts_apiMajorVersion` and `dbarts_apiMinorVersion` into
  one. The five-call-site conclusion is unaffected.
- **B1 UPHELD. `referenceCodes` needs a per-column sentinel.** `uint32_t` has
  none: 0 is a legal reference code, so `NULL = none` is an all-or-nothing
  absence bit, and a mixed source that legitimately declares a reference for its
  categorical CSC columns necessarily carries some value on its ordinal CSC
  columns. Refusing on non-NULL falsely refuses every mixed source; refusing
  only on nonzero passes "declared reference 0 against an ordinal column"
  silently - the exact hole cheap-uniformity rule 2 exists to close, at the one
  value `PredictorSource::referenceCodeOf` already returns for an absent array
  (data.hpp:220-222). The engine's own view carries no declared/absent bit
  either: the R layer keeps it one level up as a SEPARATE field,
  `ParsedMutationSource::referenceMeta` (R_interface_bartcore.cpp:204), "R
  INTEGER, NA_INTEGER == undeclared", which is exactly what the promoted helper
  keys on (`refuseCscReferenceAgainstStore`, :613-620). **Resolution:** the POD
  field becomes `const int32_t*`, `< 0` = column declared none, `>= 0` = the
  code, refused above `bartcore::maxCategories` (0xFFFF, data.hpp:95). This
  keeps the memo's decoupling goal (it still does not spell `xint_t`) and
  mirrors `columnSources`' own negative-means-something-else idiom in the same
  struct. F3 becomes constructible with a non-vacuous negative half.
- **NEW, found here, in neither document: the CSC triple is not
  self-describing.** The memo's POD carries `cscColumnPointers` ("length CSC
  columns + 1") with NO field giving the CSC column count. The library would
  have to INFER it as `max_j(~columnSources[j]) + 1` and could then never
  bounds-check the caller's pointer array or the `~v` decode. That is the same
  defect class as the missing width check the arc exists to close (F8).
  **Resolution:** add `size_t numCscColumns`. It also supplies
  `refuseCscReferenceAgainstStore`'s `numSparseColumns` argument directly, so
  the promoted helper is reused with its body UNCHANGED.
- **B2 UPHELD as a defect, remedy OVERTURNED.** `SamplerBase::setTestData(const
  PredictorSource&)` returns `bool` (facade.hpp:112-113): "Returns false without
  touching the test store when a designated leaf covariate column would be
  CSC-backed." Today's `void` flat entry is safe only because its argument is
  unconditionally dense - the header comment at facade.hpp:114-116 says so. The
  moment the entry takes a source the refusal is reachable, and a `void` entry
  discarding the bool leaves the store holding its PREVIOUS rows, so the next
  `dbarts_sampler_run` fills `results.test` for the old test set. Silent wrong
  answer. The critique preferred re-signing to `int`. **This plan does NOT.**
  Three reasons: (i) on the critique's own distinguishing criterion - "a
  legitimate caller condition rather than a malformed argument" - a CSC-backed
  leaf covariate is a FIXABLE MALFORMED ARGUMENT, which is why the R funnel
  raises it as an error with a repair instruction ("supply it as a dense test
  column", R_interface_bartcore.cpp:3672-3673); (ii) the entry already
  `Rf_error`s for its other malformed arguments (`validateColumnValues`), so an
  `int` for one failure mode and a longjmp for the rest is a worse contract than
  one uniform failure mode; (iii) `dbarts_sampler_predict` is a read-only replay
  that builds no store and so has NO bool to forward - it must use the explicit
  helper regardless, and one refusal implementation with one message across both
  flat test entries is what cheap-uniformity rule 5 asks for. **Resolution:**
  `setTestPredictors` stays `void`; it runs the promoted
  `refuseSparseLeafCovariate` BEFORE calling `setTestData`, and translates a
  `false` return to `Rf_error` anyway as defense in depth - the landed pattern
  from zero-weight-exactness S2 item 3 ("Translate a false return to Rf_error
  (defense in depth: the probe should already have refused)"). `predict` runs
  the same helper. F9 pins it with a mandatory negative half.
- **B3 UPHELD. The "zero consumer calls" justification for removing
  `dbarts_apiVersion` is FALSE for the fifth consumer.** True for stan4bart and
  treatSens (0 hits each, reproduced twice). False for
  `inst/tinytest/capi/consumer.c:16-27`, which resolves it BY HAND through
  `R_GetCCallable` with a written-out cast, labelled "Deliberate canary ... the
  un-stubbed per-symbol path a consumer that declines DBARTS_USE_STUBS (or a
  diagnostic tool) still relies on. Everything else goes through the stubs, so
  this one raw path guards that plain R_RegisterCCallable registration keeps
  working on its own." Removing it deletes the ONLY in-repo coverage of one of
  the two consumer paths the shipped header documents. `test-capi.R:43-50` pins
  `1000L` and the packed identity. `man/dbarts-package.Rd:43` documents the flat
  API as "versioned by `\code{DBARTS_C_API_VERSION}`"; `inst/NEWS.Rd:362` names
  it too. **Resolution:** the removal stands; the canary is RE-POINTED at
  `dbarts_apiHash` in the same commit (named in S1 item 6, gated by F6), and
  `man/dbarts-package.Rd` joins S2's doc list. Under binding decision 8 the
  re-point is doubly right: with no version constant moving, the hash is the
  only runtime lockstep signal there is, so the raw path must guard the symbol
  that carries it.
- **B4 UPHELD, with the probe, and now MORE load-bearing. `dbarts_apiHash()`
  does NOT catch "ANY signature change".** The token is FNV-1a over
  `DBARTS_C_API_DECLS`, which stringizes only return types, names and parameter
  lists (`#define DBARTS_API_STRINGIZE(ret, name, params, args) #ret " " #name
  #params ";"`, dbarts.h:277). It never sees a struct definition. Measured:
  three headers differing only in `dbarts_results`' layout - including
  `uint32_t* varcount` retyped to `uint64_t* varcount`, a hard ABI break - all
  hash to `0xf760898d116cb3a3`, identical, `declsBytes` 3023 in all three; the
  token text does not contain "varcount" at all. The same run independently
  re-derived the baked literal at dbarts.h:82. This arc makes the exposure worse
  in two ways at once: it adds a SECOND ABI struct, and binding decision 8
  promotes the token from one of two lockstep signals to the ONLY one.
  **Resolution:** the header says what the token actually covers - the
  entry-point SIGNATURES, not the layout of the structs they name - and points
  at `structSize` plus the exact-offset locks (C_interface.cpp:64-76) as the
  layout contract. **Do not ship "catches ANY signature change."** A struct
  layout change in this window is therefore NOT self-detecting and must be
  announced by hand to the sister packages; say so in the header and in
  c-api-growth.md.
- **B5 UPHELD, verbatim. The setForestWeights reservation states the INVERSE
  polarity, in two landed plans.** `docs/plans/c-api-growth.md:502-504`:
  "Nonzero return means refusal, matching `dbarts_sampler_setPredictor` and its
  siblings". `docs/plans/zero-weight-exactness.md:452`: "nonzero on refusal".
  Both are wrong on the fact they cite (binding decision 5). No falsifier
  covered the entry. **Resolution:** the entry is built with 1 = accepted /
  0 = refused; both plans get an ERRATUM edit at landing (text below, verbatim);
  F10 pins the polarity with an inverted negative half.
- **B6 UPHELD in substance, OVERTURNED in scope.** S0's "Bodies unchanged" is
  false as written, but the transitive closure is ONE type, not six symbols.
  (i) `refuseSparseLeafCovariate` (:761-772) takes `(const SamplerShape&, const
  PredictorSource&)` - both nameable in the common header. Body unchanged.
  (ii) `refuseCscReferenceAgainstStore` (:609-623) takes plain scalars and
  pointers, all header-nameable, and keys on `NA_INTEGER` per CSC column. Its
  body ALSO stays unchanged: `C_interface` already allocates O(p) scratch to
  convert `columnTypes` to `ColumnType` and `referenceCodes` to `xint_t`, so it
  fills one more `std::vector<int>` of length `numCscColumns`, `NA_INTEGER`
  where the POD says undeclared, and calls the helper as-is. One rule, one
  implementation, one message; the translation is explicit and F3 drives it.
  (iii) Only `validateTestContainerAgainstStore` (:1742-1760) is a real problem,
  and only because its parameter type `ParsedTestContainer` (:182-189) is an
  anonymous-namespace R-parsing struct that cannot be named in a header. Its
  body touches nothing but `parsed.view` (verified line by line), so it
  re-signs to `(const ColumnStore&, const PredictorSource&)` - a body change
  plus four call sites (:3666, :3741, :4655, :4720). **The critique's claim that
  `rawViewColumn`, `ParsedCscCodes`, `parsedCscCodes`, `refuseInvalidCategory
  Codes` and `categoricalTestMessage` must also be promoted is WRONG.** They are
  internal-linkage CALLEES defined earlier in the same translation unit
  (:1602-1643, anonymous namespace :37-2196); an external-linkage function
  defined lower in the same TU calls them legally. Nothing about them moves.
- **NEW, found here: S0's printer widening is under-scoped.** The memo names
  only `Chain::printTree` (chain.hpp:1916, `Forest<L, ResidT>& forest =
  forests_[0]`). `Chain::printSavedTree` (chain.hpp:1936, `const Forest<L,
  ResidT>& forest = forests_[0]`) has the same hardcode and is the branch
  `Sampler::printTrees` takes when `options_.keepTrees` (sampler.hpp:945-946).
  BOTH printers take the forest index, or `printTrees(forest = 1)` prints tau's
  trees without `keepTrees` and mu's with it.
- **B7 UPHELD. The defect handoff is not two lines, and the reachable guard is
  the wrong predicate.** `refuseBCFTestSurface` is defined at
  R_interface_bartcore.cpp:2097-2105 INSIDE the anonymous namespace (:37-2196):
  not declared in `R_interface_bartcore_common.hpp`, not reachable from
  `C_interface.cpp`. The one guard C_interface can see,
  `refuseMultiForestMutation` (common.hpp:115), is deliberately a DIFFERENT
  predicate: `refuseBCFTestSurface` fires on `shape.numForests >= 2 &&
  !shape.testFitsAreDefined`, gated on `testFitsAreDefined` "so a multi-forest
  model whose test blend IS defined (multinomial softmax over the K forests'
  totalTestFits) is allowed through". Guarding the flat entries with
  `refuseMultiForestMutation` would over-refuse a future flat-creatable
  multinomial - the exact case the engine comment says must pass.
  **Resolution:** the fix is the fourth namespace promotion of the same kind, it
  belongs to bcf-public-surface S3 by the orchestrator's sequencing call, and
  this plan hands that arc the exact mechanics under "TODO edits at landing".
- **A2 ADOPTED, narrowed.** "Under stubs a re-signed entry is a COMPILE ERROR at
  the consumer's call site" is TRUE for every re-sign this arc makes - all four
  predictor entries change arity or change a parameter from `const double*` to
  `const dbarts_predictor_source*` (no implicit conversion), and the three
  forest-indexed entries change arity; arity is reliable because `params` is
  also used as the pointer type `ret (*) params`, which makes C++ default
  arguments structurally impossible in the X-list. It is NOT a property of the
  stub mechanism: `R_GetCCallable` matches by name only and the cast is
  unchecked, so scalar retypes, return-type changes and same-type reorders all
  convert silently. Say it about these re-signs, never as a general rule.
- **A3 ADOPTED. multinomial-counts-mutation item 3 is SUPERSEDED, not adopted.**
  It asks that a tagged struct preserve the `refuseMultiForestMutation` guards
  on the flat RESPONSE-side entries; this arc declines the tagged struct (so
  there is nothing for it to bind) and bcf-public-surface S3 item 2 RELAXES
  exactly those guards. The memo's "adopted verbatim" pointed at a different
  guard (`refuseMultiForestTransactionalUpdate`) on different entries.
- **A6 ADOPTED. State the dense short-circuit.** The bridge does NOT materialize
  a plain dense argument: `if (Rf_isReal(xExpr)) { ...; values = REAL(xExpr); }
  else { ... materializeMutationSource(...); }`
  (R_interface_bartcore.cpp:3960-3975). The flat mutation entries take
  `source->isDenseBlock() ? source->denseValues : materialize(...)`, or the
  dense path acquires an unconditional O(n x p) copy that contradicts this
  plan's own byte-identical claim - and F1 would catch it.
- **A7 ADOPTED. Two validation obligations were unspecified**: a `columnTypes`
  value outside `{0, 1}` (the POD's `int32_t` admits 2^32 values; `ColumnType`
  is `enum class : uint8_t {ordinal, categorical}`, data.hpp:93) and a
  `categoryCounts` entry above `maxCategories`. Both are unchecked funnels into
  a quantizer. Folded into F8.
- **A8 ADOPTED with a different mechanism.** `DBARTS_PREDICTOR_SOURCE_DENSE(x,
  n, p)` expanding to `{ sizeof(...), (n), (p), (x) }` is reader-hostile and
  breaks first if anyone inserts rather than appends. The critique proposed
  designated initializers; this plan uses a `static inline` constructor instead,
  because a designated initializer is a warned extension below C++20 and a
  LinkingTo consumer building at CXX17 would take a `-Wc99-designator` hit that
  CRAN reports. `static inline` is already the header's own idiom (the stubs).
- **A9 ADOPTED.** S1's preconditions (a landed bcf-public-surface S1 and S3) are
  gated explicitly, not asserted in prose.
- **A5 ADOPTED as an OBSERVATION, not a gate.** The 34.5x / 9.96x sparse-predict
  payoff (cheap-uniformity S4 landing, 9929ede: n_test 1e5, p 1e4, density 0.01;
  0.075-0.082 s / ~1.0 GB against 2.59-3.47 s / 10.0-12.4 GB) is an R-surface
  measurement. F2 gates bitwise parity only. One order-of-magnitude C-side
  wall-clock and peak-RSS observation goes in the landing note; it is NOT a
  gate, needs no quiet machine, and must not be quoted as a benchmark.
- **The steelman of fork 1B is recorded and still loses.** Appending sparse
  siblings and leaving the dense entries alone would need no re-sign anywhere
  and no stale-binary window at all. It loses because the dense entries' absent
  width parameter is a LIVE memory-safety defect, measured (below), and B would
  preserve it forever behind a compatibility argument while needing a third
  entry or a break anyway to fix it. Fork 1A gets the width check for free
  because `numColumns` becomes explicit.

## The defect the arc closes, MEASURED

Reproduced live against installed dbarts 1.0.0 by driving
`dbarts_sampler_predict` on a gaussian sampler with p = 5, handing it the
caller's intended `nTest x 4` argument laid out inside an `nTest x 5`
allocation whose fifth column the caller never meant to pass (so nothing reads
outside an allocation):

    both calls returned normally: no refusal (no width parameter exists)
    max |predict(tail=0.10) - predict(tail=0.90)| = 1.0035
    identical answers: FALSE
    VERDICT: the library CONSUMED the tail -> NO width check
    == setTestPredictors, same narrow argument ==
    result: ACCEPTED without complaint, numTestObservations now 30

A one-unit swing on the response scale, silently, from bytes outside the
caller's intended matrix. Today both entries infer p from the sampler and index
`x + j * numTestObservations` for `j < shape.numPredictors` with no width check
(C_interface.cpp:288-300 and :310-328). This is the arc's strongest single
justification and F8's premise; it is not cosmetics.

## Verified seams (read at 2e50cf1)

- **The X-list is the single source; registration is fully mechanical.**
  `DBARTS_C_API_LIST` (dbarts.h:171-272) expands into the consumer stubs, the
  readable prototypes' bind-asserts (`DBARTS_BIND_ASSERT`, C_interface.cpp:441),
  the CCallable table (`DBARTS_API_REGISTER`, R_interface.cpp:256-265) and the
  FNV-1a token (`static_assert(dbarts_fnv1a(DBARTS_C_API_DECLS) ==
  DBARTS_C_API_HASH)`, C_interface.cpp:93). Adding or removing an entry costs
  ZERO registration lines. Baked literal `0x1a911c00bb26dcd7ULL` (dbarts.h:83);
  re-baking is a one-literal edit dbarts's own compile forces.
- **Both compiled consumers use `DBARTS_USE_STUBS`** (stan4bart
  `src/Makevars.in:1`, `src/Makevars.win:1`; treatSens `R_interface.cpp:27`,
  `bartTreatmentModel.cpp:10`, `sensitivityAnalysis.cpp:32`). See A2 for what
  that does and does not buy.
- **Both run a load-time handshake, both header-driven, and under binding
  decision 8 both are INERT for this window.** stan4bart
  `checkDbartsAPIVersion` (init.cpp:957-966, predicate `major !=
  DBARTS_C_API_MAJOR || minor < DBARTS_C_API_MINOR` at :961, invoked from
  `R_init_stan4bart` at :1058); treatSens `R_init_treatSens`
  (R_interface.cpp:446-448, same predicate). Both compare against
  `DBARTS_C_API_MAJOR` / `DBARTS_C_API_MINOR` taken from the header, so a plain
  rebuild moves the baked expectation with ZERO source edits - and, since no
  constant moves, a STALE consumer binary also passes. The hash is the only
  signal that can distinguish them; see S2's optional handshake edit. Neither
  calls `dbarts_apiVersion`; neither uses `DBARTS_C_API_VERSION`.
- **No consumer carries an upper version bound.** All four declare
  `dbarts (>= 1.0-0)` and dbarts's own DESCRIPTION is `1.0-0`, so the R-level
  constraint cannot distinguish pre- from post-reshape dbarts either. Within
  this window the only runtime discriminator available is
  `dbarts_apiHash() != DBARTS_C_API_HASH`.
- **`PredictorSource`** (data.hpp:184-239) is the ten-field borrowed view:
  numRows, numColumns, denseValues, cscColumnPointers, cscRowIndices, cscValues,
  columnSources, columnTypes, categoryCounts, referenceCodes.
  `columnSources == nullptr` is the identity map; a negative entry names CSC
  column `~v` (`sourceOf`, :208). `isDenseBlock()` (:231) is the predicate the
  mutation kernels gate on. `xint_t` is `uint16_t` (data.hpp:20) and
  `ColumnType` is `enum class : uint8_t` (:93) - BOTH are storage decisions the
  reduced-precision work reserves the right to change (public-surface.md sec 7,
  per-column u8 code widths), so neither may be spelled into the shipped ABI.
  `maxCategories` is `0xFFFFu` (:95). The view is trivially destructible by
  static_assert (:243), because `Rf_error` longjmps past destructors.
- **The engine's mutation virtuals REFUSE a non-dense view.**
  `SamplerBase::setPredictor(const PredictorSource&, ...)` (facade.hpp:129-133)
  is documented "only a dense block is consumable"; sampler.hpp:1003 and :1017
  return `PredictorUpdateResult::unsupportedSource` otherwise. The R bridge
  therefore MATERIALIZES first (`materializeMutationSource`,
  R_interface_bartcore.cpp:837-856, over `bartcore::materializePredictorSource`,
  data.hpp:305). **This is the most important fact in the plan**: a
  source-shaped C mutation entry buys uniformity, an explicit shape and
  validation - NOT resident sparse mutation, and NOT a lower whole-matrix peak.
  cheap-uniformity S1 item 5 states the same honest limit for R.
- **The predict and test-ingestion virtuals DO consume a sparse view resident.**
  `SamplerBase::predict(const PredictorSource&, size_t, double*)` is itself a
  virtual (facade.hpp:191-192) with a dense convenience spelling over
  `densePredictorSource` (:200-206); `Sampler::predict` (sampler.hpp:502-510)
  branches on `isDenseBlock()` and otherwise routes rows through
  `PredictorSourceColumns` rank bitmaps; `setTestData(const PredictorSource&)`
  (facade.hpp:112) is the resident test store. So predict and test ingestion are
  where the reshape delivers a MEASURED capability; mutation delivers uniformity
  only. Say both, in the header.
- **`bartcore_bridge::getTrees` ALREADY takes a forest index**
  (R_interface_bartcore_common.hpp:96-103; `bartcore_getTrees` reads forestExpr
  and range-checks against `shape.numForests`, R_interface_bartcore.cpp:4688,
  message "forest index out of range" at :4689). `dbarts_sampler_getTrees`
  hardcodes 0 (C_interface.cpp:354). `SamplerBase::numTreesInForest(f)` exists
  (facade.hpp:260). `savedTree`, `savedTreeSlopes`, `savedTreeMasks` and
  `flattenTree` all take `forestIndex = 0` (facade.hpp:167-188).
- **printTrees is forest-0-only in the ENGINE, on BOTH branches.**
  `Sampler::printTrees` (sampler.hpp:924) calls `printTree` without keepTrees
  (:936) and `printSavedTree` with it (:946); both read `forests_[0]`
  (chain.hpp:1917 and :1937). A forest-indexed flat entry needs the engine
  widening on both.
- **`setTreeStorage` is per SAMPLER by construction.** `Sampler::setTreeStorage`
  (sampler.hpp:851) calls `chain->initializeSavedTrees(capacity)`, which
  allocates saved slots for the variance forest AND every entry of `forests_`
  (chain.hpp:1797-1806). There is no per-forest storage to select.
- **`shape.numTrees` is forest 0's count on every multi-forest sampler** - the
  BCF constructor sets `options_.numTrees = spec.mu.numTrees` and the
  multinomial one `spec.forest.numTrees`, each with the comment "single-forest
  queries (numTrees, savedTree, printTrees) address forest 0" (sampler.hpp:166,
  :195). So `dbarts_sampler_numTrees` is silently forest-0 on a BCF the moment
  BCF is flat-creatable.
- **Three helpers the C surface needs live in the ANONYMOUS namespace of
  R_interface_bartcore.cpp** (`namespace {` :37, closed :2196):
  `refuseCscReferenceAgainstStore` (:609), `refuseSparseLeafCovariate` (:761),
  `validateTestContainerAgainstStore` (:1742). `validateColumnValues` (:2290)
  and the four refusal helpers C_interface already uses sit in
  `bartcore_bridge` (:2198-2758) and are declared in
  `R_interface_bartcore_common.hpp`. **Commit 7299b8b is the exact promotion
  precedent** and its five mechanics are reused verbatim in S0 and quoted for
  bcf S3 below.
- **`setForestWeights` exists at the engine.** `SamplerBase::setForestWeights
  (size_t forestIndex, const double* weights)` returns `bool`, true on install
  (facade.hpp:246-249), with the `SamplerFacade` fan-out at :451-453 and the
  shape predicate at :74. Bridge `bartcore_setForestWeights` and R
  `bartcoreSetForestWeights` landed at 153d1dd. Only the flat entry is missing.

## S0. Plumbing. No dbarts.h change, no hash move.

1. Promote three helpers from the anonymous namespace of
   `R_interface_bartcore.cpp` into `bartcore_bridge`, by the 7299b8b mechanics:
   move the definition DOWN below `} // namespace` (:2196) into the
   `bartcore_bridge` block; copy its doc comment into
   `R_interface_bartcore_common.hpp` beside `validateColumnValues`; append one
   sentence to the moved comment naming the flat C API as the second consumer;
   add a `using bartcore_bridge::<name>;` at the top of
   `R_interface_bartcore.cpp` so the existing unqualified call sites still
   resolve. The three:
   - `refuseSparseLeafCovariate` - body UNCHANGED.
   - `refuseCscReferenceAgainstStore` - body UNCHANGED (S1 adapts C-side; see
     the B6 adjudication).
   - `validateTestContainerAgainstStore` - RE-SIGNED to
     `(const bartcore::ColumnStore&, const bartcore::PredictorSource&)`, since
     `ParsedTestContainer` is anonymous-namespace and unnameable in a header,
     and the body reads nothing but `parsed.view`. Four call sites become
     `validateTestContainerAgainstStore(<store>, parsed.view)`: :3666, :3741,
     :4655, :4720. Its callees stay put - they are internal-linkage functions
     called from the same translation unit and nothing about them moves.
2. Widen the printers with a forest index, defaulted nowhere - explicit at every
   caller: `Chain::printTree(size_t t, int indentation, size_t forestIndex)`
   (chain.hpp:1916) and `Chain::printSavedTree(size_t slot, size_t t, int
   indentation, size_t forestIndex)` (:1936) read `forests_[forestIndex]`;
   `Sampler::printTrees` (sampler.hpp:924) takes and forwards it; the
   `SamplerBase` virtual (facade.hpp:230) and the `SamplerFacade` override
   (:433) gain the parameter. Existing callers pass 0. Range-check in the
   BRIDGE, not the engine (house convention; the engine is fast over safe).
3. `tests/cpp`: `printTrees(forest = 1)` on a BCF configuration asserting it
   prints tau's trees, on BOTH branches (keepTrees off and on); and forest 0
   printing today's output byte-for-byte on both.

rng: NEUTRAL. Gates: `R CMD INSTALL --preclean` into
`.claude/dbarts-h-reshape-design/privlib` with `R_LIBS` set (a facade virtual
moves - a stale object bus-errors; never the user library); delete the
`benchmarks/kernels` binaries (no header dependency tracking); `tests/cpp` from
clean, plain AND ASAN (`ASAN_OPTIONS=detect_container_overflow=0`); full
tinytest with NO snapshot regenerated; the trio bitwise, LABELLED a formality;
`air format --check .`.
ABORT: any trio divergence; any change to forest-0 printed output on either
branch.

## S1. The one dbarts.h edit and the hash re-bake.

PRECONDITION, gated not assumed: bcf-public-surface S1 and S3 have LANDED
(`dbarts_sampler_numForests` must exist for the range checks, a flat BCF for
F5's second leg, and a flat BCF for F10's coupling refusal). Verify by grepping
the live X-list for `numForests` before starting; if it is absent, STOP and
report - do not build a private substitute.

1. **The POD**, size-first, in the ABI-types block beside `dbarts_results`:

       typedef enum { DBARTS_COLUMN_ORDINAL = 0, DBARTS_COLUMN_CATEGORICAL = 1 }
         dbarts_column_type;

       typedef struct dbarts_predictor_source_t {
         size_t structSize;              /* set to sizeof(dbarts_predictor_source) */
         size_t numRows;
         size_t numColumns;
         const double* denseValues;      /* column-major, numRows x its own columns */
         size_t numCscColumns;           /* CSC columns; 0 when there is no CSC part */
         const int* cscColumnPointers;   /* length numCscColumns + 1 */
         const int* cscRowIndices;
         const double* cscValues;
         const int32_t* columnSources;   /* NULL = identity; >= 0 dense col; < 0 CSC col ~v */
         const int32_t* columnTypes;     /* dbarts_column_type per column; NULL = all ordinal */
         const uint32_t* categoryCounts; /* declared K per column; NULL/0 = infer */
         const int32_t* referenceCodes;  /* per column; < 0 = declared none */
         /* 1.0-0 field boundary: appends go below, never above. */
       } dbarts_predictor_source;

       #define DBARTS_PREDICTOR_SOURCE_INIT { sizeof(dbarts_predictor_source) }

       static inline dbarts_predictor_source
       dbarts_dense_predictor_source(const double* values, size_t numRows,
                                     size_t numColumns);

   `DBARTS_PREDICTOR_SOURCE_INIT` initializes only the first member, so C
   zero-initializes and C++ value-initializes the rest; both give NULL/0 and
   neither depends on field order. The dense constructor replaces the memo's
   argument-reordering macro (A8). Add exact-`offsetof` and exact-`sizeof`
   static_asserts in `C_interface.cpp` beside the `dbarts_results` ones
   (:64-76). Generalize `DBARTS_RESULTS_HAS` (dbarts.h:137) into
   `DBARTS_HAS_FIELD(type, ptr, field)` and re-express the old spelling over it;
   zero consumer uses either (verified, all four repos). A zero `structSize`
   errors, as `dbarts_sampler_run` already does (C_interface.cpp:135; the
   comment explaining why is at :131-134).
2. **Re-sign the four predictor entries.**

       X(int,  dbarts_sampler_setPredictor,
         (dbarts_sampler*, const dbarts_predictor_source*, int forceUpdate,
          int updateCutPoints), ...)
       X(int,  dbarts_sampler_updatePredictor,
         (dbarts_sampler*, const dbarts_predictor_source*, const size_t* columns,
          size_t numColumns, int forceUpdate, int updateCutPoints), ...)
       X(void, dbarts_sampler_setTestPredictors,
         (dbarts_sampler*, const dbarts_predictor_source*), ...)  /* NULL removes */
       X(void, dbarts_sampler_predict,
         (dbarts_sampler*, const dbarts_predictor_source*,
          const double* offsetTest, double* out), ...)

   Bodies, in one shared C_interface validation helper so the four cannot
   diverge: refuse `structSize == 0`; refuse `numColumns` disagreeing with
   `shape.numPredictors` (`updatePredictor` compares against its own
   `numColumns` argument, whose columns are in ARGUMENT order exactly as
   `bartcore_updatePredictor` treats them, R_interface_bartcore.cpp:4035-4042);
   refuse a `columnTypes` entry outside `{0, 1}`, a `categoryCounts` entry above
   `maxCategories`, a `referenceCodes` entry above `maxCategories`, and a
   `columnSources` entry naming a CSC column `>= numCscColumns` or a dense
   column `>= numColumns`; build the O(p) `ColumnType`, `uint32_t` and `xint_t`
   scratch plus the `NA_INTEGER`-encoded per-CSC-column adapter, inside the
   unwind-protected frame; run the promoted `refuseCscReferenceAgainstStore`
   (cheap-uniformity rule 5 - "the malformed-reference refusal applies at EVERY
   funnel", and the flat surface has no `validateXTest`); run the existing
   per-column `validateColumnValues` loop.
   Then: **mutation MATERIALIZES**, `isDenseBlock() ? source->denseValues :
   materializePredictorSource(...)` (A6), and NEVER hands a non-dense view to
   the engine - the engine would return `unsupportedSource`, which the
   `accepted ? 1 : 0` mapping reports to the caller as an ordinary rollback, a
   silent lie (F7). **predict and setTestPredictors pass the view straight
   through**, after `refuseSparseLeafCovariate`; `setTestPredictors` also
   translates `setTestData`'s `false` to `Rf_error` as defense in depth (F9).
   One header sentence states the mutation/predict asymmetry so no reader
   over-reads the sparse capability.
3. **Forest-index the three tree queries**, trailing `size_t forest`, 0 meaning
   exactly today's behavior:

       X(size_t, dbarts_sampler_numTrees, (const dbarts_sampler*, size_t forest), ...)
       X(SEXP,   dbarts_sampler_getTrees, (..., int useLiveTrees, size_t forest), ...)
       X(void,   dbarts_sampler_printTrees, (..., size_t numTreeIndices, size_t forest), ...)

   `forest >= dbarts_sampler_numForests(sampler)` raises `Rf_error` naming the
   entry, matching the bridge's "forest index out of range". `numTrees` returns
   `size_t` with no error channel (binding decision 5), so its out-of-range case
   must error rather than return 0 - a 0 tree count is indistinguishable from a
   legitimate answer to a caller that does not check.
   **Do NOT build a forest-indexed `setTreeStorage`**: storage is per sampler
   (chain.hpp:1797-1806), so its only legal value would be "all forests". Record
   CLOSED BY FACT. **Do NOT build a forest-indexed `predict`**: it needs
   per-forest saved-tree replay, which bcf-public-surface itself holds as a door
   with named consumers. Shipping the symbol without the engine behind it is
   worse than not shipping it.
4. **Versioning: the token moves, the constants do NOT.** Append
   `X(uint64_t, dbarts_apiHash, (void), ())` returning `DBARTS_C_API_HASH`.
   Remove `dbarts_apiVersion` from the X-list and `DBARTS_C_API_VERSION` from
   the header. **Leave `DBARTS_C_API_MAJOR` at 1 and `DBARTS_C_API_MINOR` at 0**
   (binding decision 8; VD 2026-08-10, "No need to increment versions - no
   bartcore version has been released"). Re-bake `DBARTS_C_API_HASH` - that
   re-bake is the whole acknowledgment. Header text, three sentences: (i) the
   supported post-release handshake remains major-equality plus a minor floor,
   and the constants become a contract only at the first release; (ii) within
   the pre-release window the constants do not move, so
   `dbarts_apiHash() == DBARTS_C_API_HASH` is the lockstep-rebuild check and the
   only runtime signal a stale consumer binary trips; (iii) the token covers the
   entry-point SIGNATURES, not the layout of the structs they name - `structSize`
   plus the exact-offset locks carry the layout, a struct layout change is NOT
   self-detecting, and it must be announced to the sister packages by hand
   (B4, measured).
5. **Append `dbarts_sampler_setForestWeights`**:

       X(int, dbarts_sampler_setForestWeights,
         (dbarts_sampler*, size_t forest, const double* weights),
         (sampler, forest, weights))

   at the END of the X-list. **1 = accepted, 0 = refused** (binding decision 5;
   the two recorded reservations say the inverse and are corrected at landing).
   `weights == NULL` clears. Body: capability probe on
   `shape.supportsForestWeights` FIRST (never a forest count, so a multinomial
   cannot slip through), then `forest < numForests`, then every element finite
   and `>= 0`; forward `SamplerBase::setForestWeights`' bool. Ownership: the
   flat entry BORROWS - it retains nothing and copies nothing, so the caller
   owns the buffer for the sampler's life. State that in the Doxygen; it is the
   one place the flat entry's contract differs from the bridge's, which copies
   into a holder-owned buffer. The precision channel and the mean channel stay
   TWO entries forever: `s_{f,i}` scales forest f's own leaf conditionals and
   never enters `combinedFits` (measured, 1.8e-15), while a basis scales
   forest f's contribution to the mean. Widening this entry with a basis is
   refused on the rule zero-weight-exactness applied to
   `dbarts_sampler_setWeights`.
5b. **The mean channel, CONDITIONAL on multiforest-extension-surface's fork 3
   (see binding decision 3's carve-out).** Re-sign
   `dbarts_sampler_setTreatment` as
   `X(int, dbarts_sampler_setForestBasis, (dbarts_sampler*, size_t forest,
   const double* basis, size_t numColumns), (sampler, forest, basis,
   numColumns))`, `basis` column-major n x numColumns; replace
   `dbarts_sampler_bcfGlue` with
   `X(size_t, dbarts_sampler_numForestAmplitudes, (const dbarts_sampler*,
   size_t forest), (sampler, forest))` and
   `X(int, dbarts_sampler_forestAmplitudes, (const dbarts_sampler*, size_t
   forest, double* out), (sampler, forest, out))` filling
   numForestAmplitudes(forest) x numChains. 1 = accepted, 0 = refused. The
   body accepts only what today's engine honours - forest 1, a two-column
   complementary 0/1 basis, Gaussian family - and refuses the rest naming the
   capability, so the family relaxes guard bodies later and moves no header.
   Ownership: the entry COPIES, matching `setTreatment` (`dbarts.h:43`); state
   it in the Doxygen, because a continuous basis cannot be coerced-and-copied
   incidentally the way a 0/1 z can. Re-word the creation Doxygen
   (`dbarts.h:348-357`) in engine vocabulary. Budget ~40 header + ~80
   C_interface + ~70 consumer.c + ~60 test-capi.R.
6. **`consumer.c` + `test-capi.R`**: every entry above; the refusal matrix; the
   `structSize` canary; the sparse-predict parity oracle; and the RE-POINTED raw
   canary - `p_apiHash_raw` resolved by hand through `R_GetCCallable` with a
   written-out cast, replacing `p_apiVersion_raw`, keeping the comment's
   explanation of why one raw path must survive (B3). Drop the `1000L` and
   packed-identity assertions; assert instead that `apiMajorVersion()` and
   `apiMinorVersion()` are STILL 1 and 0, that the raw and stubbed hash agree,
   and that the hash literal differs from the literal recorded at
   slice start (F6).
7. **Nameable-calibration's mid-chain footprint** (docs/plans/nameable-calibration.md
   sec 8 and its S0, "Signature freeze. No code."). **AMENDABLE until this
   slice starts**: if S2's implementation falsifies a signature choice, the
   plan is corrected rather than frozen. Creation half: NONE - the model
   crosses as SEXP (`dbarts_sampler_create`, `dbarts.h:175-177`) and the
   `prior.scale` -> `node.scale` conversion is engine-side, so a flat-C
   consumer reaches it with no header change. Mid-chain half: one output
   POD, one enum, two X-list entries appended at the END of
   `DBARTS_C_API_LIST`:

       typedef enum {
         DBARTS_LEAF_CONSTANT = 0,
         DBARTS_LEAF_MONOTONE = 1,
         DBARTS_LEAF_LINEAR   = 2,
         DBARTS_LEAF_GP       = 3
       } dbarts_leaf_model;

       /// Caller-owned output buffers for dbarts_sampler_forestCalibration, the
       /// dbarts_results contract: set structSize; EVERY member is a pointer and
       /// is filled only when both present-by-size and non-null, each over
       /// numChains; a zero structSize errors. Fields append below the marked
       /// boundary and never reorder. All quantities are in RESPONSE units (the
       /// family's latent units where the response is not rescaled).
       typedef struct dbarts_forest_calibration_t {
         size_t  structSize;      ///< caller sets to sizeof(dbarts_forest_calibration)
         double* priorScale;      ///< numChains; forest-total prior sd at k = 1
         double* priorSd;         ///< numChains; priorScale / k at the current k
         double* priorMean;       ///< numChains; prior mean of the forest total
         double* k;               ///< numChains
         double* responseScale;   ///< numChains; internal-to-response multiplier
         double* responseShift;   ///< numChains; internal-to-response offset
         int*    kHasHyperprior;  ///< numChains; THIS FOREST's own k law (not the
                                  ///< sampler-wide dbarts_sampler_kIsSampled,
                                  ///< which reads the sampler option and
                                  ///< disagrees on BCF and multinomial)
         int*    leafModel;       ///< numChains; dbarts_leaf_model, qualifying
                                  ///< priorSd and priorMean (see below)
         /* 1.0-0 field boundary: appends go below, never above. */
       } dbarts_forest_calibration;

       #define DBARTS_FOREST_CALIBRATION_INIT { sizeof(dbarts_forest_calibration) }

       X(int, dbarts_sampler_forestCalibration, \
         (const dbarts_sampler* sampler, size_t forest, \
          dbarts_forest_calibration* out), \
         (sampler, forest, out)) \
       X(int, dbarts_sampler_setForestPriorScale, \
         (dbarts_sampler* sampler, size_t forest, double priorScale), \
         (sampler, forest, priorScale))

   **Why the two trailing members are `int*` and not `int`** (compiled and
   measured, `cc -std=c99`, arm64 LP64): with two trailing `int`s, tail
   padding makes `sizeof` identical with and without the last field - 56
   through `responseShift`, 64 with one `int`, 64 with two,
   `offsetof(leafModel) = 60` - so a caller omitting `leafModel` sets
   `structSize = 64` and `DBARTS_HAS_FIELD(..., leafModel)` returns TRUE for
   a field it does not carry, while the exact-`sizeof` static assert cannot
   see a future sub-word append at all because it lands in existing padding.
   That would be the first POD to break the header's stated invariant
   (`dbarts_results`' own Doxygen: the library fills only fields whose end
   offset falls within `structSize`). With pointers the arithmetic is
   honest: 64 with one, 72 with two, and the omitting caller's `HAS` is
   FALSE. The uniform pointer shape also gives the two members the "null
   member skips" spelling the POD's Doxygen claims for every member, and the
   per-chain shape the rest of the struct already has.

   Prototype-view Doxygen (the `#else` branch) states, in this order: what
   the getter fills and that `priorScale` is the quantity the setter
   writes; that `priorSd` is `priorScale / k` per chain and moves every
   sweep under `kHasHyperprior` while `priorScale` does not; the
   LEAF-PARAMETER sentence (`prior.scale`/`prior.sd` describe the
   leaf-parameter scale of the forest total, equal to the prior sd of
   `f(x)` at every x for the constant leaf only) with equality for
   `DBARTS_LEAF_CONSTANT` only; then per tag - `DBARTS_LEAF_LINEAR` a LOWER
   bound attained at the standardized covariate origin, larger by
   `sqrt(1 + ||z(x)||^2)`; `DBARTS_LEAF_GP` an UPPER bound attained at rows
   reproducing a leaf member and on over-cap leaves, elsewhere
   `priorSd^2 c(x)' C^-1 c(x)` decaying to 0 as x leaves the leaf's data
   cloud, where every draw equals `priorMean`; `DBARTS_LEAF_MONOTONE` a
   LOWER bound in the interior (realized sd a few per cent to ~20% above
   it) whose `priorMean` is NOT the prior mean of `f(x)` under an active
   constraint, that marginal being skew with an x-dependent mean spanning
   several `priorSd` along the constrained axis; and that `priorMean` is
   exact for the constant, linear and gp leaves. Returns 1, or 0 without
   touching `out` when `forest` names no forest; errors on a zero
   `structSize`.

   The setter's Doxygen states: it restates forest `forest`'s leaf prior on
   EVERY chain so the forest total's prior sd at `k = 1` is `priorScale`, in
   response units; `k`, the response transform, sigma and the tree prior
   are untouched; it takes effect on the NEXT sweep and never reinterprets
   leaf values already drawn; a write reproducing the current internal
   scale bitwise is skipped, so a read-then-write is inert; to move the
   prior MEAN, shift the reported fit with `dbarts_sampler_setOffset`; the
   leaf model qualifies the write exactly as it qualifies the read. TWO
   ERROR CHANNELS, because the header's global contract raises on invalid
   arguments: a CAPABILITY answer is a RETURN VALUE - 0, touching nothing,
   when `forest` names no forest or a combiner owns this forest's
   calibration (a two-forest or multinomial sampler) - while a MALFORMED
   VALUE RAISES, namely a non-finite or non-positive `priorScale`.
   1 = accepted. The flat surface deliberately carries no `prior.sd`
   spelling, so it has no sampled-k refusal; that sugar and its refusal are
   R-side only.

   Forward compatibility: `forest` is a parameter, so the general basis
   family (multiforest-extension-surface M4) relaxes the refusal in its own
   guard body and moves no header - record that in the reshape landing
   note.

   Signature assumption the design carries, NOT a reshape obligation:
   `Chain::setModel` RE-DERIVES `leaf.scale = model.priorScale /
   response_->fitScale()` against the CURRENT transform whenever
   `priorScale` is finite - the same conversion creation runs - so that a
   no-op `$setModel(sampler$model)` does not silently REVERT a named
   calibration (MEASURED 1.5 -> 12.0, 8x, on a range-24 gaussian, absent
   this rule). The signature above assumes that rule holds; S0 records it
   as an assumption the signature carries, not as work item S1 owes.

   PRECONDITION, verified live in this worktree (4f0aeab8):
   `dbarts_sampler_numForests` IS in the X-list, at `dbarts.h:264` exactly
   as the calibration plan cites - reshape S1's own start condition is met.
   The reshape plan's zero-`structSize` anchor sentence above (item 1) reads
   "(C_interface.cpp:135; the comment explaining why is at :131-134)",
   confirmed correct against the live file: the comment sits at :131-134 and
   the check at :135. No errata to carry forward.
8. **Latent-subset-mask's flat entry** (docs/plans/latent-subset-mask.md,
   "The dbarts.h footprint (carried by dbarts-h-reshape S1)"). **AMENDABLE
   until this slice starts**, same as item 7. ONE entry, appended at the END
   of the X-list:

       X(int, dbarts_sampler_setActiveRows,
         (dbarts_sampler* sampler, const double* active), (sampler, active))

   Contract (Doxygen): `active` is length `dbarts_sampler_numObservations`,
   each element exactly `0.0` or `1.0`; `NULL` clears (every row active); an
   all-ones vector is accepted and installs nothing. Returns **1 = accepted,
   0 = refused** - the shipped convention (`dbarts_sampler_setPredictor`
   ends `accepted ? 1 : 0`; the polarity erratum in
   `docs/plans/c-api-growth.md`). No version constant moves (binding
   decision 8). **Ownership: the entry RETAINS NOTHING.** The values are
   consumed (copied) into the sampler's own buffer during the call and the
   caller's array is free immediately after it returns. This is NEITHER
   `setForestWeights`'s borrow-and-retain (which is why THAT entry obliges
   the caller to keep the array alive) NOR a copy into a holder; it is the
   predictor setters' "retain no pointer" clause. No clause joins dbarts.h's
   keep-alive list. Reachable for gaussian, Student-t, probit, logistic,
   aft, ordinal and nbinom - the families `dbarts_sampler_create` builds by
   name; multinomial and BCF have no flat creation path, so their masking
   stays `dbarts:::`-only, as `bartcore_setCounts` and
   `bartcore_setForestWeights` already are.

   Body decision (subset-mask V4): the S1 body ACCEPTS gaussian (and
   Student-t, the same `shape.family`) and refuses every other family by
   name - the pattern this plan's own item 5b proposes ("The body accepts
   only what today's engine honours ... and refuses the rest naming the
   capability, so the family relaxes guard bodies later and moves no
   header"). Gaussian masking needs no new engine work at that point: it is
   `setWeights(w * a)` composed at the entry. Validation runs the
   capability probe on `shape.supportsActiveRows` FIRST, never a family
   switch; the length is implicit (no length argument - the entry reads
   `numObservations` from the shape); the exact-`{0,1}` scan is the
   ENGINE's, inherited by the flat entry, so the r-c-division defect-4 hole
   does not reopen.

   `test-capi.R` positive-arm obligation (this plan's own S1 item 6 requires
   coverage of every entry plus the refusal matrix, calling `test-capi.R`
   "the load-bearing gate"): one positive arm (gaussian mask changes the
   fit, all-ones does not), one refusal arm per refused family reachable
   from `dbarts_sampler_create` (probit at S1 time), one fractional-value
   refusal, one `NULL`-clears arm. Because the body accepts gaussian, no
   assertion has to invert when this arc's S1 lands - only refusal arms are
   relaxed later, which is what item 5b was designed for. A body that always
   returns 0 is WITHDRAWN: it would ship a symbol that lies by omission and
   plant an assertion in another arc's gate file that must later invert.

   Two entries considered and NOT proposed, so a later reader does not
   reopen them: a reader (`dbarts_sampler_getActiveRows`) and a count
   (`dbarts_sampler_numActiveObservations`). The channel does not ride the
   state block, so the writer is the only source of the value and a reader
   can only echo it. If V6 is ever reversed and the mask DOES ride the
   state, a reader becomes necessary and must be added in the same re-bake.

rng: NEUTRAL. Gates: preclean install into the private library; `tests/cpp`
plain and ASAN from clean; `test-capi.R` (the load-bearing gate); full tinytest;
the trio bitwise, LABELLED a formality; `air format --check .`; dbarts.h
ASCII-clean and C99-clean (`consumer.c` compiles as C - the existing self-gating
check), and clean at CXX17 as well as CXX20 for the `static inline` constructor.
ABORT: a re-signed X-list that leaves `DBARTS_C_API_HASH` unchanged from
the literal recorded at slice start (the token failed to see the change);
any movement in `DBARTS_C_API_MAJOR` or `DBARTS_C_API_MINOR` (binding
decision 8); any sparse-vs-dense divergence at anything but bitwise; a missing
`dbarts_sampler_numForests` in the live X-list at slice start.

## S2. Consumers, docs, records.

1. stan4bart and treatSens rebuilt against the reshaped header - five stan4bart
   call sites, enumerated below; treatSens rebuild only.
2. **Optional, recommended, one line each: give the two handshakes a signal that
   can actually fire.** Under binding decision 8 no version constant moves, so
   `checkDbartsAPIVersion` (stan4bart init.cpp:961) and `R_init_treatSens`
   (treatSens R_interface.cpp:446) cannot distinguish a stale binary from a
   fresh one. Adding `|| dbarts_apiHash() != DBARTS_C_API_HASH` to each
   predicate restores lockstep detection for the rest of the window. **It must
   be removed or downgraded at the freeze**: post-release, a legitimate minor
   append moves the hash, and a hard equality check would fail every consumer
   until rebuild. Record that condition beside the edit, or do not make the
   edit.
3. Docs and records: `docs/design/public-surface.md` sec 6;
   `docs/plans/c-api-growth.md` (reservations closed, opened and corrected);
   `man/dbarts-package.Rd:43`; `inst/NEWS.Rd` (the new bullets and the
   historical :362 line); `TODO`; `docs/plans/INDEX.md` registration; the
   landing note in this file, carrying the F2 capability observation.

rng: NEUTRAL. Gates: both consumer packages install and pass their own suites
against the private library; full tinytest; `air format --check .`.

## Falsifiers (pre-registered)

- **F1 (S1), the neutrality gate.** Record `capi_run`'s sigma/train/varcount at
  the S0 tip; after S1 the same seed through the re-signed DENSE path must be
  BITWISE identical. The trio cannot see this - it drives no flat entry - so
  this is the only real neutrality evidence, and it is what catches an
  accidental unconditional materialization (A6).
- **F2 (S1), the sparse-predict oracle, BOTH halves.** From `consumer.c`: a CSC
  source through `dbarts_sampler_predict` is bitwise identical to a dense source
  holding the materialized same values, on {ordinal CSC, categorical CSC with a
  nonzero reference, mixed dense+CSC, all-implicit column}. NEGATIVE HALF,
  mandatory: force the implicit value to 0 for a categorical column and it must
  go red. Plus the A5 OBSERVATION (not a gate): one wall-clock and peak-RSS
  reading for a large sparse test set, recorded to an order of magnitude.
- **F3 (S1), rule 5 at the flat funnel, THREE legs.** A source declaring a
  reference code against a store-ORDINAL CSC column is REFUSED at each of
  setPredictor, updatePredictor, setTestPredictors and predict, from
  `consumer.c`, with the promoted helper's exact message. Leg 2, the one the
  `uint32_t` spelling made impossible: `referenceCodes[j] == 0` on that ordinal
  column must ALSO be refused - 0 is a legal code, not an absence. Leg 3, the
  negative half: `referenceCodes[j] < 0` (undeclared) on the same ordinal column
  in an otherwise identical mixed source must be ACCEPTED, and its result
  bitwise equal to the same source with `referenceCodes == NULL`.
- **F4 (S1), the input-struct guard.** Allocate a region larger than
  `dbarts_predictor_source`, set `structSize` to a boundary below `columnTypes`,
  and point the omitted fields at an unmapped page so a stray read segfaults;
  assert the library never reads them. NEGATIVE HALF: remove the guard and it
  must crash. The `capi_run_guard` idiom inverted for a read.
- **F5 (S1), forest addressing, BOTH halves.** `forest = 0` reproduces today's
  `getTrees` data.frame and `printTrees` output byte-for-byte on a single-forest
  sampler, on BOTH print branches (keepTrees off -> `printTree`, on ->
  `printSavedTree`); `forest >= numForests` errors on all three entries; on a
  flat-created BCF, `numTrees(0) != numTrees(1)` and `getTrees(1)` returns tau's
  trees. NEGATIVE HALF: wire the forest argument to a constant 0 and the BCF leg
  must go red on both print branches.
- **F6 (S1), the token moved and the raw path survives.** `dbarts_apiHash() ==
  DBARTS_C_API_HASH`; the baked literal differs from the literal recorded at
  slice start (a re-signed X-list that leaves it unchanged means the token is
  blind to the change - abort, do not proceed); `dbarts_apiMajorVersion() == 1`
  and `dbarts_apiMinorVersion() == 0`, i.e. the constants did NOT move (binding
  decision 8); the hand-resolved `R_GetCCallable("dbarts", "dbarts_apiHash")`
  pointer returns the same value as the stubbed call (the B3 canary,
  re-pointed); and a deliberately perturbed baked literal must still fail
  dbarts's own compile (the existing `static_assert`, with a one-line
  perturbation recipe in a comment so a later reader can re-check).
- **F7 (S1), the honesty falsifier.** A sparse mutation source must never return
  `unsupportedSource` to a caller: assert from `consumer.c` that a CSC
  `setPredictor` returns the SAME accept/reject decision as its dense
  equivalent, never a spurious 0. NEGATIVE HALF: bypass the materialization and
  it must go red.
- **F8 (S1), the self-description matrix.** Every entry taking a source refuses
  an argument that does not describe itself, one assertion per leg: (i)
  `numColumns` disagreeing with `dbarts_sampler_numPredictors`; (ii) a
  `columnTypes` entry outside `{0, 1}`; (iii) a `categoryCounts` or
  `referenceCodes` entry above `maxCategories`; (iv) a `columnSources` entry
  naming a CSC column `>= numCscColumns`. Leg (i) is written FAILING-FIRST
  against the pre-fix build, which consumed an out-of-matrix column with 1.0035
  of silent drift (measured above), so the commit shows the inversion. Leg (iv)
  runs under ASAN with a short `cscColumnPointers` array: without the bound the
  decode reads past it.
- **F9 (S1), the leaf-covariate refusal at the flat test surface, BOTH halves.**
  On a linear-leaf sampler with a designated leaf covariate, a source making
  that column CSC-backed must be REFUSED - by `dbarts_sampler_setTestPredictors`
  and by `dbarts_sampler_predict` - with the promoted helper's exact message,
  and the previously installed test store must be INTACT afterwards (assert by
  running and comparing `results.test` against the pre-call rows). NEGATIVE
  HALF, mandatory: remove both the C-side pre-refusal and the defensive
  translation of `setTestData`'s false, and the run must be shown reporting the
  OLD test rows - the silent stale store this falsifier exists to forbid.
- **F10 (S1), setForestWeights polarity, BOTH halves.** From `consumer.c` on a
  flat-created BCF: an accepted `setForestWeights(1, w)` returns **1**; a
  refused one returns **0** - drive the refusal twice, once with `forest >=
  numForests` and once on a single-forest sampler (no forest-weight coupling).
  A refused call must leave the sampler's subsequent draws BITWISE unchanged.
  NEGATIVE HALF: wire the return to the inverse polarity and both legs must go
  red. This falsifier exists because two landed plans specify the inverse.

## Migration costs per consumer (enumerated, never constraining)

Verified twice, independently, at stan4bart bartcore@6ce0440, treatSens
dbarts-1.0@bb9a121 (the worktree, not the `master` primary at d1da1dd, which
still binds the removed C++ ABI), bartCause dbarts-1.0@695c603, bairrtt
main@6167423 - whole-tree ripgrep plus a binary symbol census.

| consumer | entries used | this arc's call-site edits | bcf S3's edits | optional handshake | total |
|---|---|---|---|---|---|
| stan4bart | 23 (38 call expressions) | 5 | 0 | 1 | 5-6 + rebuild |
| treatSens | 10 (21 call expressions) | 0 | 4 | 1 | 4-5 + rebuild |
| bartCause | 0 (no `src/`) | 0 | 0 | - | rebuild only |
| bairrtt | 0 (LinkingTo Rcpp, RcppEigen) | 0 | 0 | - | rebuild only |
| in-repo `capi/consumer.c` | the gate | ~200 lines, in S1 | - | - | - |

- **stan4bart** - all five edits in `src/init.cpp`, all COMPILE ERRORS under the
  stubs: `dbarts_sampler_predict` x1 (:340) becomes a two-line dense-source
  construction plus the call; `dbarts_sampler_getTrees` x1 (:490),
  `dbarts_sampler_printTrees` x1 (:428) and `dbarts_sampler_numTrees` x2 (:390,
  :452) each gain a trailing `, 0`. UNCHANGED: `setTreeStorage` x2 (:201, :563 -
  no forest parameter), `dbarts_results` and the whole run loop, `setOffset`,
  `setSigma`, `getLatents`, `storeState`/`setState`, every `num*` probe. It
  calls none of `setPredictor`, `updatePredictor`, `setTestPredictors`,
  `setTestOffset`, `apiVersion`, `DBARTS_C_API_VERSION`, `DBARTS_RESULTS_HAS`.
  **Correction to bcf-public-surface's migration line:** stan4bart makes ZERO
  `dbarts_sampler_setResponse` calls (public-surface.md:347 recorded this for
  0.0-13; still true at 6ce0440), so bcf S3's setResponse widening costs it
  nothing but the rebuild.
- **treatSens** - ZERO call-site edits for THIS arc. Its window cost belongs to
  bcf S3: four `dbarts_sampler_setResponse` sites
  (`sensitivityAnalysis.cpp:173, :184, :278, :291`) gain an explicit
  `updateScale`, which must be **1**, not 0, to preserve today's behavior -
  `dbarts_sampler_setResponse` hardcodes `setResponse(y, true)`
  (C_interface.cpp:202). bcf S3's note ("0 for a BCF") is right for a BCF and
  wrong for treatSens's gaussian grid. Flagged for that arc, not decided here:
  three of the four sites (:173, :184, :291) are POST-BURN-IN, where
  dbarts.h:363-365 advises `false`, so preservation and desirability may diverge
  there.
- **bartCause / bairrtt** - no `dbarts_` or `DBARTS_` symbol anywhere, no
  `LinkingTo: dbarts`; they drive the R surface only.

VD's recorded estimate - "a stan4bart predict sweep only" (TODO, typed-ingestion
door status) - is CORRECT for the mutation half and INCOMPLETE for the arc as
scoped: the forest-indexed family adds stan4bart's four one-token edits.

**The number of re-bakes is NOT the number of migrations.** Consumers migrate
ONCE, at the freeze, against the final header (TODO, "re-verify stan4bart
bartcore against the reshaped header in the same window"; bcf-public-surface's
release-checklist edit). A re-bake is a one-literal edit only dbarts's own build
notices. Optimize the count for review clarity, not for migration cost, and do
not let it carry weight it does not have.

## Removals

Removal is argued on REDUNDANCY or AMBIGUITY only. Under VD's enabling-value
gate "no consumer calls it" is never the reason, so `setWeights`,
`setTestOffset`, `sampleNodeParametersFromPrior`, `setNumThin` and `usesDart` -
all with zero consumer calls - STAY.

1. **`dbarts_apiVersion()` and `DBARTS_C_API_VERSION`** - redundant with
   major/minor and unable to express the strict equality they advertise: the
   packed integer is `major * 1000 + minor`, so strict equality against it
   rejects every legitimate additive minor bump, and it does not move at all
   when a signature is REPLACED at the same major/minor - which, under binding
   decision 8, is exactly what every change in this window is. Replaced by
   `dbarts_apiHash()`. **CANARY REPLACEMENT, mandatory in the same commit:**
   `consumer.c`'s hand-resolved `R_GetCCallable` path re-points at
   `dbarts_apiHash`, keeping its comment; `test-capi.R`'s three version
   assertions are rewritten per S1 item 6; `man/dbarts-package.Rd:43` and
   `inst/NEWS.Rd:362` are edited. Without the re-point the removal deletes the
   only in-repo coverage of the un-stubbed per-symbol path the header documents.
2. **The standalone `numTestObservations` parameters** of `setTestPredictors`
   and `predict` - absorbed into `dbarts_predictor_source::numRows`. Replacement
   rather than deletion, and it is where "two sources of truth for the row
   count" dies. `numColumns` becomes explicit for the first time, which is what
   closes the measured out-of-bounds read.
3. **`DBARTS_RESULTS_HAS` as a distinct spelling** - folded into the generic
   `DBARTS_HAS_FIELD(type, ptr, field)`. Zero consumer uses (all four repos).
4. **NOT removed, recorded closed by fact:** the reserved forest-indexed
   `setTreeStorage`.

## Resolved questions (recorded, not open)

1. **The version question is MOOT. No constant moves. RESOLVED by VD ruling,
   not by argument.** VD 2026-08-10, verbatim: *"No need to increment versions -
   no bartcore version has been released."* No bartcore version has ever
   shipped, so `DBARTS_C_API_MAJOR` / `DBARTS_C_API_MINOR` have never been a
   compatibility contract with anyone; whatever the header carries at the first
   release simply IS the initial contract. The memo recommended MAJOR 1 -> 2 and
   the critique endorsed it; both are SUPERSEDED. `MAJOR` stays 1, `MINOR` stays
   0, here and in bcf-public-surface S3. What survives from that analysis, and
   why it now matters more: no consumer declares an upper version bound, both
   handshakes read their expectation from the header, and four sister packages
   sit at four branch tips with installed builds on this machine and on
   dbarts-bench - so with the constants pinned, `dbarts_apiHash()` is the ONLY
   runtime signal that can distinguish a stale consumer binary from a fresh one.
   Design the machinery for that job (S1 item 4, S2 item 2, F6), and keep B4's
   caveat in the header: the token sees signatures, not struct layout.
   The post-release rule is untouched and still binds (TODO, "on any
   post-release dbarts.h change: bump DBARTS_C_API_MINOR (or MAJOR if
   incompatible) alongside the hash re-bake, and bump stan4bart's
   Depends/LinkingTo dbarts floor in the same lockstep release").
2. **Build `dbarts_sampler_setForestWeights` in S1. RESOLVED on the evidence.**
   The reservation's stated ground was "the flat API has no BCF creation entry
   point to reach it from" (c-api-growth.md:484-490) - verifiably its whole
   rationale - and it EXPIRES at bcf-public-surface S1, which this arc runs
   after. Under VD's enabling-value gate the absence of a consumer today is
   never the gating fact, and named classes exist (stan4bart's multilevel BCF;
   the latent-treatment sensitivity class; zero-inflated log-linear BART, which
   keeps all n rows with augmented zero weights). The engine, bridge and R
   channels all landed at 153d1dd; only the flat entry is missing. Cost ~20
   C_interface + ~30 consumer.c lines inside a re-bake that is happening anyway.
   Deferring costs a POST-RELEASE append, which under the still-binding
   post-release rule means a minor bump PLUS a lockstep bump of stan4bart's
   Depends/LinkingTo floor in the same release - a recorded, non-trivial cost
   the in-window build avoids entirely. B5's polarity erratum is an argument for
   fixing the record, not for deferring the entry.
3. **The tagged response/offset reservation is DECLINED; new names are reserved
   instead. RESOLVED.** multinomial-counts-mutation.md's "What the reshape must
   reserve" item 1 asks for a tagged response and offset source expressing
   `{double* vector}` and `{int* counts, size_t numCategories}`. Declined,
   because: the header's own evolution rule is "function additions arrive under
   new names (a minor bump)" (dbarts.h:22-24); a counts response is a different
   data type, not a differently shaped y; `dbarts_sampler_create` routes only to
   `createHolder` (C_interface.cpp:104-110), so the flat surface cannot produce
   a multinomial sampler at all and every tag but the vector one is unreachable
   dead surface; building it would RE-SIGN `setResponse` a second time in the
   same window, moving treatSens's four sites twice; and the reservation's
   stated goal - no ABI break later - is met exactly as well by a new name.
   A plan may decline another plan's reservation when the substitute meets the
   reservation's own goal, the departure is recorded where the reserving author
   will see it, and the reservation is not yet in the registry. All three hold:
   multinomial-counts-mutation has not landed and c-api-growth.md carries no
   multinomial reservation today. Reserved in its place, both APPENDS:
   `dbarts_sampler_setCounts(dbarts_sampler*, const int* counts, size_t
   numCategories)` and `dbarts_sampler_setOffsetMatrix(dbarts_sampler*, const
   double* offset, size_t numCategories)`. Its item 2 (a size-first spec struct
   for a future flat creation) is ADOPTED. Its item 3 is SUPERSEDED, not
   adopted - see the A3 adjudication.
4. **The flat BCF test-surface gap belongs to bcf-public-surface S3.** Settled
   by the orchestrator's sequencing call; the amendment text is below.

## Open decisions (VD)

**None.** All four of the memo's questions are resolved above - three on
recorded principles and measured evidence, and the fourth (versioning) by VD's
2026-08-10 ruling, which made it moot.

## Doors held open (recorded, not scheduled)

- **Forest-indexed `predict`** - needs per-forest saved-tree replay. Named
  consumers: `predict.bartBCF` on new rows, bairrtt's MH filter.
- **Flat per-observation predictor update and the joint session** - the R5
  surface has both; the flat surface has neither. Additive later.
- **Flat `setCutPoints` and `setData`** - additive later; `setData` inherits the
  multi-forest refusal survey in `runsbcbcf-repair.md`.
- **Flat `predictVariance`** - the virtual exists (facade.hpp:195-197) with a
  dense convenience spelling; no flat entry. Additive later, and it takes the
  same source POD when it arrives.
- **Flat multinomial / counts creation, `setCounts`, `setOffsetMatrix`** -
  reserved names (resolved question 3).
- **A tagged response/offset source** - DECLINED in favour of new names; the
  reservation stays recorded so it is not re-proposed as fresh.
- **Per-column u8 code widths** (public-surface.md sec 7) - the reason the ABI
  must not spell `xint_t` or `ColumnType`.
- **A flat `setTestPredictorAndOffset`** - the bridge has one (the row count may
  change with the offset); the flat surface cannot change the test row count
  while an offset is installed. Additive later.
- **Struct-layout detection.** The API token is blind to it (B4, measured), so
  a layout change in a future window is announced by hand. A second token over
  the ABI structs' offsets is the obvious fix; unscheduled, and this arc
  prejudges nothing about it beyond recording that `structSize` plus the
  exact-offset locks are today's whole contract.
- **Resident sparse MUTATION** - not delivered here and must not be claimed. The
  engine's mutation virtuals consume only a dense block, so the flat entries
  materialize exactly as the bridge does.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

The 1.0-0 entry describes the flat C API as it SHIPS, not as it changed - no
prior version exists to have changed from (binding decision 8).

- S1: the flat C API's predictor entries take a single self-describing source
  struct, so a C consumer hands the sampler compressed-column storage for
  prediction and test data without densifying it and without an R container
  class, and every predictor argument declares its own width and CSC column
  count.
- S1: `dbarts_sampler_numTrees`, `getTrees` and `printTrees` take a forest
  index, so they are unambiguous on a BCF sampler; `dbarts_apiHash()` gives a
  consumer an exact signature-lockstep check alongside the major/minor
  handshake; `dbarts_sampler_setForestWeights` reaches the per-forest weight
  channel from C.

## TODO edits at landing

1. **The `bcf-public-surface` S3 amendment, to be applied VERBATIM by the
   orchestrator before that slice starts** (resolved question 4's mechanism, per
   the landed SL precedent at commit 7299b8b). Insert as a new S3 item, after
   item 2:

   > 2b. **Guard `dbarts_sampler_predict` and `dbarts_sampler_setTestPredictors`
   > with `refuseBCFTestSurface`.** S1 makes BCF flat-creatable and thereby
   > makes a SILENT WRONG ANSWER reachable: `Sampler::predict` routes to
   > `predictColumns`, whose `Chain::predictFromSavedSample` /
   > `predictFromCurrentTrees` both open `const Forest& forest = forests_[0]`
   > and loop `forests_[0].numTrees`, so a flat BCF consumer receives mu(x)
   > labelled as the fit. The R bridge already guards all four of its siblings
   > (`bartcore_predict` :4646, `bartcore_setTestPredictor` :3656,
   > `bartcore_setTestOffset` :3702, `bartcore_setTestPredictorAndOffset`
   > :3730); of the flat siblings only `setTestOffset` is guarded, and only
   > incidentally, by `refuseMultiForestMutation` (C_interface.cpp:306).
   > This is NOT two lines. `refuseBCFTestSurface` is defined at
   > R_interface_bartcore.cpp:2097-2105 INSIDE the anonymous namespace
   > (:37-2196): not declared in `R_interface_bartcore_common.hpp` and not
   > reachable from `C_interface.cpp`. The one guard C_interface can already
   > see, `refuseMultiForestMutation`, is the WRONG predicate - it fires on
   > `numForests >= 2`, while `refuseBCFTestSurface` fires on `numForests >= 2
   > && !shape.testFitsAreDefined`, deliberately gated on `testFitsAreDefined`
   > so a future flat-creatable multinomial (whose test blend IS defined) passes
   > through. Using it here would over-refuse exactly the case the engine
   > comment says must be allowed.
   > Mechanics, the five steps 7299b8b took for
   > `refuseMultiForestTransactionalUpdate`: (i) move the definition down out of
   > the anonymous namespace into the `bartcore_bridge` block (:2198-2758),
   > beside `refuseMultiForestMutation`; (ii) append one sentence to its comment
   > - "External linkage: the flat C API reuses this guard on its own predict
   > and test-predictor entries."; (iii) copy the doc comment as a declaration
   > into `R_interface_bartcore_common.hpp` inside `namespace bartcore_bridge`;
   > (iv) add `using bartcore_bridge::refuseBCFTestSurface;` at the top of BOTH
   > `R_interface_bartcore.cpp` (so its four existing unqualified call sites
   > still resolve) and `C_interface.cpp`; (v) call it at C_interface.cpp:288
   > (`dbarts_sampler_setTestPredictors`) and :310 (`dbarts_sampler_predict`),
   > naming each entry. Budget: ~35 bridge/common-header + ~25 consumer.c.
   > Extend F10 with two legs: on a flat-created BCF, `predict` and
   > `setTestPredictors` each refuse with the BCF message. NEGATIVE HALF: swap
   > the guard for `refuseMultiForestMutation` and a flat multinomial (when one
   > exists) must go red - until then, assert the predicate directly in
   > `tests/cpp` on a shape with `numForests >= 2 && testFitsAreDefined`.

   Also amend that plan in three further places:
   - S3 item 4 currently bumps `DBARTS_C_API_MINOR` 0 -> 1. Under VD 2026-08-10
     ("No need to increment versions - no bartcore version has been released")
     it does NOT: re-bake the hash and leave both constants alone.
   - S3's ABORT clause "a hash re-bake that is not accompanied by the minor
     bump" is void; replace with "a re-signed or appended X-list that leaves
     `DBARTS_C_API_HASH` unchanged, or any movement in `DBARTS_C_API_MAJOR` /
     `DBARTS_C_API_MINOR`".
   - Its "TODO edits at landing" line "note that the minor floor moved to 1" is
     void; drop it, and say instead that the constants stay 1.0 until the first
     release, when whatever they read becomes the initial contract.
   - S3 item 1's migration note: "0 for a BCF" is right for a BCF and wrong for
     treatSens's existing gaussian grid, whose four sites need **1** to preserve
     behavior; write both, and flag that three of the four (:173, :184, :291)
     are post-burn-in, where dbarts.h:363-365 advises `false`.
2. **`docs/plans/c-api-growth.md`, the setForestWeights reservation
   (:496-512).** ERRATUM, not an adoption. Replace constraint 1's first clause
   with: "**Erratum (dbarts-h-reshape, 2026-08-10):** this reservation
   originally read 'Nonzero return means refusal, matching
   `dbarts_sampler_setPredictor`'. That is the INVERSE of the shipped
   convention - `dbarts_sampler_setPredictor` ends `return result ==
   bartcore::PredictorUpdateResult::accepted ? 1 : 0` (C_interface.cpp), i.e.
   **1 = accepted, 0 = refused**, which is also what bcf-public-surface S3 item
   3 records. The entry is built with 1 = accepted / 0 = refused." In the same
   edit, strike "`DBARTS_C_API_MINOR` bumped when it lands" - it landed in this
   window with no version movement (VD 2026-08-10) - and mark the reservation
   BUILT, naming this arc's S1.
3. **`docs/plans/zero-weight-exactness.md`, S3 item 2.** Same erratum, two
   sentences: the phrase "nonzero on refusal" is wrong and reads "1 on
   acceptance, 0 on refusal", and the "`DBARTS_C_API_MINOR` bump" clause is void
   under VD 2026-08-10. The flat entry landed in dbarts-h-reshape S1. Do not
   restate the derivation - point at c-api-growth.md.
4. **`docs/plans/c-api-growth.md`, reservations closed and opened:**
   forest-indexed `setTreeStorage` CLOSED BY FACT (storage is per sampler,
   chain.hpp:1797-1806); forest-indexed `predict` DOOR, blocker = per-forest
   saved-tree replay; `dbarts_sampler_setCounts` / `dbarts_sampler_setOffsetMatrix`
   reserved as appends, SUPERSEDING multinomial-counts-mutation's
   tagged-response reservation item 1; the pre-release rule (VD 2026-08-10) that
   a re-bake carries no version movement, and the post-release rule that it must
   (REPLACE => MAJOR, APPEND => MINOR, plus the lockstep consumer floor); and
   the measured fact that `dbarts_apiHash()` is blind to struct layout, so
   `structSize` plus the exact-offset locks are the layout contract and a layout
   change is announced by hand.
5. **`docs/plans/multinomial-counts-mutation.md`, "What the reshape must
   reserve".** Item 1 DECLINED with the reasoning in resolved question 3, and
   the substitute names recorded; item 2 ADOPTED; item 3 SUPERSEDED - it asks
   that the flat response-side `refuseMultiForestMutation` guards be preserved,
   and bcf-public-surface S3 item 2 relaxes exactly those guards by decision.
6. **`TODO`.** Replace the typed-ingestion door-status reshape sentence with the
   landing record: plan path, artifacts under
   `.claude/dbarts-h-reshape-design/`, the slices landed, and the fact that no
   version constant moved. Correct the "consumer cost: a stan4bart predict sweep
   only" estimate to "five stan4bart call sites (predict, getTrees, printTrees,
   numTrees x2); treatSens, bartCause and bairrtt rebuild only". Add to the
   release checklist: the C ABI constants read 1.0 and become the initial
   contract at that release, and all four sister packages are re-verified
   against the final re-baked header once, at the freeze. Leave the existing
   post-release bump rule untouched - it still binds from the first release on.
7. **`man/dbarts-package.Rd:43`** - "versioned by `\code{DBARTS_C_API_VERSION}`"
   becomes the major/minor handshake, naming `dbarts_apiHash()` as the exact
   signature check. **`inst/NEWS.Rd:362`** - the historical 1.0-0 line names a
   macro that will not ship; edit deliberately.
8. **`docs/design/public-surface.md` sec 6** - the flat surface's shape, the
   source POD, the forest-indexed family, and the honest mutation-vs-predict
   asymmetry. Bump its `Status:` line per the plans README landing rule.
9. **`docs/plans/INDEX.md`** - one row for this plan.

## Departures from the memo and the critique (record)

1. **No version increment, against BOTH documents.** VD 2026-08-10, verbatim:
   *"No need to increment versions - no bartcore version has been released."*
   The memo's fork 4(a) recommended MAJOR 1 -> 2 / MINOR -> 0 and the critique
   independently endorsed it ("Recommendation A is correct"); both are
   SUPERSEDED. Neither constant moves, here or in bcf-public-surface S3. The
   consequence is not cosmetic: it makes `dbarts_apiHash()` the only runtime
   lockstep signal in the window, which is why S1 item 4, S2 item 2 and F6 are
   written the way they are, and why B4's "the token is blind to struct layout"
   caveat must ship in the header rather than being a footnote.
2. **B2's remedy is OVERTURNED.** The blocking finding is adopted in full - a
   `void` `setTestPredictors` taking a source would discard a reachable engine
   refusal and leave a stale test store. The critique's preferred fix (re-sign
   to `int`) is not taken. The entry stays `void` and pre-refuses C-side through
   the promoted `refuseSparseLeafCovariate`, with a defensive `Rf_error`
   translation of the bool. Reasons in the adjudication: the condition is a
   fixable malformed argument by the critique's own criterion, the entry already
   longjmps for its other malformed arguments, and `predict` must use the
   explicit helper regardless since it has no bool to forward.
3. **B6's transitive closure is OVERTURNED in scope.** Only the
   anonymous-namespace parameter TYPE `ParsedTestContainer` blocks a header
   declaration. The five named callees stay in the anonymous namespace: an
   external-linkage function defined lower in the same translation unit calls
   internal-linkage functions legally. The re-price is therefore ~65 bridge
   lines, not the larger, less mechanical figure the critique implied.
4. **`refuseCscReferenceAgainstStore`'s body does NOT change**, against the
   critique. The C entry builds an `NA_INTEGER`-encoded per-CSC-column adapter
   in the O(p) scratch it already allocates, so one rule keeps one
   implementation and one message, and F3 drives the translation.
5. **Two defects neither document found are folded in.** (i) The POD had no
   `numCscColumns`, so the CSC triple was not self-describing and no bound could
   be checked on the `~v` decode - the same defect class the arc exists to
   close; the field is added and F8 leg (iv) gates it. (ii) S0's printer
   widening missed `Chain::printSavedTree` (chain.hpp:1936), the branch
   `Sampler::printTrees` takes under `keepTrees`, which also hardcodes
   `forests_[0]`; both printers are widened and F5 gates both branches.
6. **A8 is adopted with a different mechanism**: a `static inline` constructor
   rather than designated initializers, which are a warned extension below
   C++20 and would cost a CXX17 LinkingTo consumer a CRAN-visible warning.
7. **The memo's four VD questions are RESOLVED, not asked.** Questions 2 and 3
   on recorded principles and measured evidence, 4 by the orchestrator's
   sequencing call, and 1 by VD's ruling above. The "Open decisions (VD)"
   section is deliberately empty.
8. **A5 is adopted as an observation, not a gate.** The 34.5x / 9.96x figure is
   an R-surface measurement and stays labelled one; the C-side reading is
   order-of-magnitude, recorded in the landing note, and is never a benchmark.
9. **The memo's anchor drift is drift, not error.** It read at 2ca0453; the tip
   advanced to 2e50cf1 with TODO and docs edits only, `src/` and `inst/include/`
   byte-identical. Every anchor in THIS plan was re-read at 2e50cf1 and several
   were corrected by a line or two in the process (`R_interface.cpp` register
   table :256-265; `sampler.hpp:1017`; `C_interface.cpp:354` for the getTrees
   hardcode; `treatSens R_interface.cpp:446-448`).
