# capi-dispatch-table

status: DECIDED (VD 2026-07-16): Decision 0 two-component version
encoding, Decision 1(i) X-macro single-source stubs, Decision 2(i)
in-header constexpr hash + baked static_assert; DESCRIPTION floor
lockstep adopted. The ASAN cross-repo contract job was landed and then
DROPPED by VD the same day (see the C0 landing note). Additive
escalations recorded, not built. Implementation plan at the end of this
file; landing notes appended per commit.
(Design pass 2026-07-15; adversarial review 2026-07-16; three blind
second opinions 2026-07-16, each designing from the raw problem before
reading this memo - all three rejected a get-api table independently and
reproduced the skew analysis below.)

Standing context for every claim here: NOTHING in dbarts.h is published.
The whole surface, DBARTS_C_API_VERSION included, first becomes a contract
at 1.0-0 submission. Every mechanism below is therefore a green-field
design choice, free today and permanent after the freeze; the incident
history is dev-time evidence about failure modes, not a shipped-contract
constraint.

## The question

VD's ask, verbatim intent: is it worth breaking with R's mechanism for
passing C function pointers by symbol (R_RegisterCCallable /
R_GetCCallable per function) and instead registering a single function
that returns a function table, indexed by a version string?

Motivating incident (dev-time, both trees unreleased): the getTrees ABI
break (2e2b1c9, restored a73ca50). dbarts_sampler_getTrees gained a
trainingData parameter under an unchanged symbol name and an unchanged
version constant. LinkingTo consumers call through their OWN declaration
of the signature, so stan4bart passed eight arguments to a nine-argument
function; the garbage ninth became the replay-predictor pointer and
corrupted extracted node counts stack-dependently (failed ~half of R CMD
check runs for a week). Not a link error, not a load error - call-time
stack garbage. This is exactly the class 1.0-0 must be designed against,
because after the freeze the same slip would ship.

## What the ground truth already is

Registration ([[src/R_interface.cpp:298-322@126fb2cd]]): ~33 entry points in a
C_callMethods table, each R_RegisterCCallable'd under its own name in
R_init_dbarts. Append-only-by-name discipline; a signature change is an ABI
event (docs/plans/README.md, the a73ca50 review-checklist step).
The current version constant is a single integer ([[dbarts.h:54@11888173]];
dbarts_apiVersion, [[C_interface.cpp:76@11888173]], returns one int) - an unpublished
dev convention with zero compatibility weight, replaceable at will until
submission.

Consumer (stan4bart, the SOLE reverse-LinkingTo package - CRAN confirms;
every other revdep is R-level, and stan4bart ships lockstep, same
maintainer):
- It already collapses the surface into ONE table struct: BARTFunctionTable
  (stan4bart's `src/init.cpp` lines 56-90), 22 of the ~33 symbols, populated by 22
  R_GetCCallable calls in lookupBARTFunctions (stan4bart's `src/init.cpp` lines 1001-1029), with 39
  invocation sites routed through bartFunctions.*. So "one lookup + a
  table" is a one-time ~21-line saving in the one consumer, not a per-site
  win.
- Its dev builds already perform a load-time version handshake:
  apiVersion() != DBARTS_C_API_VERSION errors before any pointer is used
  (stan4bart's `src/init.cpp` lines 1003-1006). Dev-time context, but it shows the consumer-side check is
  a habit this pair actually keeps.
- Each signature is stated THREE times on the consumer side even though it
  #include <dbarts/dbarts.h> (stan4bart's `src/init.cpp` line 23): the header prototype, the struct field
  type, and the bit_cast target type at the lookup. Neither of the latter
  two is derived from the first; a rebuild recompiles them verbatim. The
  getTrees skew lived in exactly that gap - the key mechanical fact for
  the analysis below.

Growth precedent already in this tree: dbarts_results is a
structSize-growable output struct ([[dbarts.h:83-102@11888173]], DBARTS_RESULTS_HAS) -
old callers are never written past, appends never reorder. The state format
reads its blocks BY NAME behind an encoding floor ([[dbarts.h:237-243@11888173]],
c-api-growth.md part 2). Both are self-describing-DATA designs; growth of
the FUNCTION surface is already handled by append-only names (a future
run2-style addition is just a new name, load-detectable when absent). A
dispatch table buys nothing for growth that these do not already provide.

## The decisive axis: mechanical vs disciplinary compatibility token

An early draft claimed "the table relocates the hazard, it does not remove
it." False as literally written. The honest statement:

- A HAND-VERSIONED table (a human bumps the version, get_api vets it)
  cannot catch change-without-bump - the token did not move, so no
  handshake anywhere can see the change. Same exposure as per-symbol.
- An AUTO-TOKEN table CAN catch it: both sides carry a token derived
  mechanically from the declaration set itself; the handshake compares
  tokens, so a change the author forgot to acknowledge is still visible.

So the decisive axis is not dispatch shape (per-symbol vs table) but
whether the compatibility token is DISCIPLINARY (a human bumps an integer)
or MECHANICAL (derived from the declarations). The auto-token table loses
anyway, on two grounds: timing - the same mechanical token enforced at
dbarts's OWN BUILD (Decision 2 below) catches the class before any header
ships, rather than at load in the field - and convention - it diverges
from the documented R mechanism for zero residual safety once the
build-time token exists. Dominated, not wrong.

## The rebuilt-consumer dimension

The incident was a REBUILD scenario: stan4bart was recompiled against the
nine-arg header repeatedly during that week (R CMD check builds it). The
operative staleness was not an old binary - it was the HAND-ROLLED
declaration (stan4bart's `init.cpp` lines 73-75, line 1017), which no rebuild re-derives. Had the
call at stan4bart's `init.cpp` line 531 gone through a type stated IN THE HEADER - a table
struct's field, or an inline stub - the rebuild would have produced a
COMPILE error (eight args passed to a nine-arg type). So for the motivating
incident, table-struct-in-header and header stubs are safety-EQUIVALENT:
both catch it at compile time; only hand-rolled declarations yield garbage.
The table must be rejected on the honest grounds - R-convention divergence
plus a registration rewrite for zero safety gain over stubs - not on a
false safety claim.

## Skew matrix

Designs: (a) status quo - per-symbol, hand-rolled consumer declarations,
hand-bumped single-integer version; (b) get_api table, struct in dbarts.h,
hand-versioned, structSize-growable, NULL on unknown version; (c)
per-symbol + header-derived consumer declarations (stubs) + mechanical
token. Outcomes: COMPILE (building the consumer), LOAD (before any BART
call), CALL (stack garbage at call time), BUILD (dbarts's own compile).

  scenario                        (a) status quo  (b) table       (c) stubs+token
  -------------------------------------------------------------------------------
  NOT-REBUILT consumer (binary built against header vN, library now vM):
  needs entry old lib lacks       LOAD            LOAD (missing   LOAD (stub's
  (vN > vM)                       (GetCCallable   get_api or NULL GetCCallable
                                  errors)         for unknown vN) errors)
  ADD only, append-only kept      SAFE            SAFE (past old  SAFE
                                                  structSize)
  CHANGE, not acknowledged        CALL            CALL (same      BUILD: dbarts
                                                  slot, stale     itself fails to
                                                  cast)           compile (D2)
  CHANGE, acknowledged (bumped)   LOAD            LOAD (welded:   LOAD
                                  (handshake)     NULL is forced) (handshake)
  -------------------------------------------------------------------------------
  REBUILT consumer (recompiled against current header, own decls as-is):
  CHANGE, bumped or not           CALL <-- the    COMPILE (header COMPILE (stub
                                  incident; hand- struct field is signature is
                                  rolled decl not new type; call  new type)
                                  re-derived      site mismatches)

Reading. The not-rebuilt CHANGE-without-acknowledgment row is CALL garbage
under BOTH (a) and (b) - only a mechanical token touches it, and enforcing
it at dbarts's own build stops the bad header from ever shipping. The
rebuilt row is where (b) and (c) both win and (a) alone loses - and it is
the row that actually happened. (b)'s only edge over (a) that (c) lacks is
welding the consumer-side NULL check to table access; the sole consumer
already performs that check voluntarily.

One trap the matrix exposes: a token alone does not close the rebuilt row.
Trace: the token forces the bump; a not-rebuilt consumer trips the
handshake at load - good. The maintainer rebuilds; versions realign; the
handshake now PASSES; but if the hand-rolled declaration was not also
widened (the original mistake - a header bump does not touch init.cpp),
the call is garbage again, this time WITH a passing handshake. The token
protects the mismatch window only. Header-derived declarations remove the
root cause: they are re-derived by every rebuild. Hence Decision 1 is
load-bearing, not optional polish.

## The handshake predicate (unanimous across the blind reviews)

All three blind designs found independently that a ONE-integer version
cannot carry a safe relaxed handshake: a bare installed >= built lets a
BREAKING bump through (a v1 consumer accepts a v3 library that changed
signatures), and strict equality forces a consumer rebuild on every
additive change. The safe relaxed predicate needs TWO components:

  major_installed == major_built && minor_installed >= minor_built

major = incompatible change, minor = additive. The encoding is a free
choice today - the current constant is unpublished and carries no weight -
but the freeze makes whatever ships at 1.0-0 permanent: adding a second
component afterward is itself an ABI event. So this is not a repair, it is
choosing the right green-field encoding while the choice is free
(Decision 0). Until a non-lockstep consumer exists, stan4bart may keep
strict equality on top; the two-component encoding is what makes the
relaxed predicate POSSIBLE later without another ABI event.

## Sub-choices, if a table were chosen anyway

- version string vs numeric: numeric major.minor (as Decision 0). Strings
  invite bespoke parsing.
- frozen table per major vs structSize-growable: growable, matching
  dbarts_results; frozen-per-major cannot absorb an additive minor without
  reverting to per-symbol registration for every in-major addition.
- unknown version: NULL, never best-effort - a table shaped for the wrong
  version re-opens call-time garbage.
- the struct lives in dbarts.h, which is also its limit: a NOT-rebuilt
  consumer still calls through a stale copy of it (the matrix row above).

## Declaration source: stubs, and the single-source (X-macro) form

Hand-written Matrix-style stubs cost ~33 inline wrappers, ~100-130 shipped
header lines, and each signature stated TWICE in the header (prototype +
stub) - a residual dual-statement the first draft's ledger accepted. Two of
the three blind designers converged independently on eliminating it:
declare the API ONCE as an X-macro list in dbarts.h,

  #define DBARTS_API_LIST(X) \
    X(sampler_create, dbarts_sampler*, (SEXP, SEXP, SEXP, const char*)) \
    X(sampler_run, void, (dbarts_sampler*, size_t, size_t, dbarts_results*)) \
    ...

and let the preprocessor expand it into the function-pointer typedefs, the
consumer stubs (or a bind-all helper), and dbarts's OWN registration table
- true single source, in-header, no external codegen, CRAN-safe. Doc
locality (the tooling designer's variant): keep the existing
Doxygen-documented prototypes as the human-readable surface and BIND the
X-list to them per entry (a static_assert via __builtin_types_compatible_p,
or assigning &dbarts_##name to the generated typedef), so prose stays next
to signatures while any drift between prototype and list fails dbarts's own
compile. This also hands any token mechanism a machine-readable declaration
set for free. numpy's single-source codegen (below) is the proof the
ingredient matters; the X-macro is its zero-toolchain form.

## Token mechanism, conditional on the declaration source

- If the X-macro is adopted: a header-embedded constexpr FNV hash over the
  stringized X-list (function signatures AND dbarts_results field
  spellings), with a baked #define that dbarts's own build static_asserts
  against. Change a signature and dbarts ITSELF fails to compile until the
  bake is updated - the update IS the mechanical acknowledgment. No parser,
  no committed golden, no CI dependency; it fires on every local R CMD
  INSTALL (VD installs locally constantly; CI-only guards miss hotfix and
  bypass paths). Layering is essential: the hash is the PROVIDER-side
  build gate (cannot forget); the major/minor handshake stays the
  CONSUMER-side load gate. Hash inequality alone must NOT be the load
  predicate - purely additive changes move the hash and must not break
  stale binaries. Provenance, honestly: two blind designers evaluated a
  declaration hash only as paired with the table or as release-side
  codegen; the in-header constexpr form coupled to the X-macro appeared
  only in the third design - but it dominates the alternatives on every
  axis this memo's constraints prioritize (no parser work, no CI coupling,
  catches local installs, one-person-maintainer cheap).
- If VD declines the X-macro: the fallback token is a CI guard comparing a
  committed golden against the header's normalized declaration set. Its
  cost must be stated honestly: NOT a raw text diff - a73ca50's header
  diff was ~3 declaration lines amid ~8 doc-comment lines, prototypes span
  lines ([[dbarts.h:214-220@11888173]]), and dbarts_results is field-order-sensitive
  (a mid-struct insertion must flag as an ORDER change, never wave through
  as an addition) - so it is a small normalizing C-declaration parser, real
  work under a one-person-maintainer constraint, and it only fires on CI.
  Both token forms share a correctness dependency: every ABI-crossing type
  stays defined inside dbarts.h (currently true - dbarts_results and the
  callback typedef are inline; control/model/data cross as SEXP); if an
  ABI type ever moves to another header, token and table alike go blind.
- Status quo token: the hand-bumped integer plus the a73ca50 checklist.
  Human discipline alone on the row that already failed once.

## Layers beyond the header (process; adopt without a fork)

- Pre-submission cross-repo contract job: build stan4bart's PINNED DEV
  SOURCE against the candidate dbarts and run its suite UNDER ASAN (plus
  stack-protector) - ASAN turns the incident's ~half-of-runs corruption
  deterministic. The existing revdep-smoke.yaml is structurally unable to
  catch this class: monthly, not per-PR ([[.github/workflows/revdep-smoke.yaml:10@126fb2cd]]); pulls stan4bart from CRAN,
  not dev ([[.github/workflows/revdep-smoke.yaml:49-55@126fb2cd]]); unsanitized R CMD check. This job is also the ONLY
  layer that catches SEMANTIC signature drift - same types, changed
  meaning or units - to which stubs, hashes, and diffs are all blind. Note
  that CRAN's own single-shot revdep check passes a 50%-flaky corruption
  half the time; it is not a backstop.
- DESCRIPTION floor asymmetry: stan4bart's LinkingTo/Depends dbarts (>= X)
  bumped in lockstep closes the CONSUMER-AHEAD direction at install/attach
  time; the load-time token closes LIBRARY-AHEAD. They must ship together;
  neither substitutes for the other.
- CRAN enforcement synergy for stubs: CRAN rebuilds reverse-LinkingTo
  packages on a dbarts submission, so once consumer calls go through
  header stubs, an arity change the consumer did not follow makes
  stan4bart FAIL TO COMPILE during CRAN's own rebuild - auto-blocking the
  dbarts update. Free strengthening of Decision 1 that hand-rolled
  declarations defeat.
- Recorded ADDITIVE escalations (adopt-later, no freeze cost): (i)
  signature-tagged registration names - register each function under
  name#sigtag derived from the X-list, keeping the plain alias for
  diagnostics - which converts stale-third-party-binary skew into a
  load-time lookup failure even when every other guard is bypassed; (ii)
  expose the baked content hash to consumers for load-time comparison if a
  non-lockstep consumer ever appears.

## Precedents

- Writing R Extensions documents this exact failure for R_GetCCallable:
  "this mechanism is fragile, as changes to the interface provided by
  packA have to be recognised by packB. The consequences of not doing so
  have included serious corruption to the memory pool of the R session" -
  and assigns the remedy to the CONSUMER's declaration. That is the
  canonical justification for moving declaration authorship into dbarts.h.
- Matrix - the largest LinkingTo C API in R - is per-function
  R_GetCCallable through an inline stub layer (Matrix_stubs.c,
  R_MATRIX_INLINE, e.g. M_cholmod_start), NOT a get-api table, with a
  numeric ABI version (inst/include/Matrix/version.h) bumped 1 -> 2 at
  1.7-0 and a documented rebuild-reverse-LinkingTo-on-bump rule.
- xts ships the exact cached-pointer static-inline stub idiom in
  inst/include/xtsAPI.h - a decade-proven, single-witness form cleaner
  than Matrix's macro layer.
- Counter-precedent, conceded: numpy's PyArray_API IS a function-pointer
  table - import_array() fetches it through a capsule and performs a
  runtime ABI check. Tables are proven outside CRAN. But numpy's safety
  comes from append-only slots + a handshake + SINGLE-SOURCE CODEGEN
  (generate_numpy_api.py emits both sides from one spec, so dual
  declarations cannot exist) - the same ingredients recommended here in
  per-symbol form. Ingredients, not shape; the X-macro is the
  zero-toolchain equivalent of numpy's generator.
- In-tree: dbarts_results structSize and the state format's by-name blocks
  are both self-describing-DATA precedents; neither argues for a
  self-describing FUNCTION table, since function growth is already
  append-only-by-name and load-detectable.

## Recommendation (the fork is VD's)

Do NOT adopt the get_api table. A hand-versioned table adds no safety on
either dangerous row; a mechanically-tokened table is dominated by the same
token enforced at dbarts's own build, without diverging from the R
convention; and on the rebuilt row - the incident that motivated the ask -
a table is exactly equivalent to header stubs, which do not disturb the
dispatch mechanism. All three blind designs reached this verdict
independently.

Keep per-symbol R_GetCCallable and decide, all green-field until 1.0-0:

- DECISION 0, version encoding (decide before the freeze; permanent
  after). Ship the 1.0-0 version constant as two components, major/minor -
  two macros + accessors, or a packed integer with documented
  decomposition - with the consumer handshake predicate
  major_installed == major_built && minor_installed >= minor_built.
  major = incompatible, minor = additive. The encoding is unconstrained
  today (the current single integer is unpublished); the freeze is what
  makes the choice permanent - a second component added later is itself
  an ABI event, so the richer form is free now and unreachable after.
  Ship cost: a handful of header lines. Cost of shipping one integer
  instead: every future additive release either forces a consumer rebuild
  (strict equality) or is indistinguishable from a breaking one (bare >=).
- DECISION 1, declaration source (the root-cause fix). Options:
  (i) RECOMMENDED - X-macro single-source: the API declared once as an
  X-list in dbarts.h, expanding into dbarts's own registration table, the
  consumer stubs/bind helper, and the typedefs; existing Doxygen
  prototypes kept and compile-time-bound to the list. Adopt: the macro
  machinery in the shipped header (~ the same 100-130 lines, but stated
  ONCE), a less conventional header idiom to teach; list discipline
  replaces prototype discipline.
  (ii) Hand-written Matrix/xts-style stub pairs: conventional and
  precedented, but each signature stated twice in the header, with
  stub-vs-prototype drift caught only by review.
  (iii) Status quo hand-rolled consumer declarations: zero dbarts cost;
  rebuilds never re-derive the consumer's declarations, so the incident
  class stays reachable even with Decisions 0 and 2 (the
  realigned-handshake trace above), and the CRAN-rebuild synergy is
  forfeited. Under (i) or (ii), stan4bart deletes BARTFunctionTable, its
  22 typedefs, and 22 bit_cast lookups (~90 lines; three statements per
  signature drop to zero) and the getTrees class becomes a compile error
  on every rebuild - structurally protecting any future second consumer
  the 1.0-0 freeze will bind.
- DECISION 2, token mechanism (conditional on Decision 1). Options:
  (i) RECOMMENDED IF X-MACRO - in-header constexpr hash + baked
  static_assert: dbarts's own compile fails on any unacknowledged
  declaration change, locally and on CI alike; no parser, no golden.
  (ii) Fallback if no X-macro - CI declaration-diff guard + committed
  golden: honest cost is a normalizing C-declaration parser (comment/
  whitespace stripping, multi-line prototypes, field-order sensitivity),
  and it fires only on CI.
  (iii) Status quo - hand-bumped integer + checklist: human discipline
  alone on the row that already failed once.
  Either mechanical token stays PROVIDER-side; the consumer load gate
  remains the Decision 0 handshake (hash inequality must not gate loads -
  additive changes move the hash and must not break stale binaries).

Adopt as process, no fork needed: the pre-submission ASAN cross-repo
contract job (the only semantic-drift catcher); the DESCRIPTION floor
bumped in lockstep (closes consumer-ahead; the token closes
library-ahead). Recorded for later, additive: signature-tagged
registration names; exposing the baked hash to consumers.

Timing: dbarts 1.0-0 freezes the dbarts.h compatibility contract at
submission; until then every choice above is free. Decision 0's richer
encoding is unreachable after the freeze; Decisions 1 and 2 are cheap now
and awkward after; the process layers are freeze-independent but cheapest
to stand up while stan4bart's port is already in flight.

## Implementation plan (VD sign-off 2026-07-16)

C0 - contract CI: LANDED e98d94f, then DROPPED by VD decision the same
day. Rationale: per-push CI building a pinned dev branch of a downstream
package is not an R-ecosystem practice (the convention is consumer-side
testing plus pre-submission revdep checks, which the release procedure
already has); the pin is a silent-rot hazard in the wrong repo; and the
goal it served is met by C1/C2 themselves - after the rework, an
out-of-date stan4bart fails NOISILY on its own (rebuild against a
changed header = stub compile error; stale binary = handshake load
error), which is all that was required. The DESCRIPTION-floor lockstep
practice stays on the release procedure.

C1 - dbarts header + bridge: the X-macro list in dbarts.h (name, return,
parameter list per entry); Doxygen prototypes KEPT and gated - consumer
TUs define DBARTS_USE_STUBS to get same-name static inline stubs
(cached-pointer R_GetCCallable per xts/Matrix) in place of the extern
prototypes; dbarts's own TUs get the prototypes and a provider-side
binding assert (decltype comparison in C_interface.cpp) so list-vs-
prototype drift fails dbarts's compile. Registration table in
R_interface.cpp generated from the same list. Version encoding ships as
major/minor (fresh at 1.0-0: major 1, minor 0) with accessors; the
constexpr FNV hash over the stringized list + dbarts_results field
spellings is static_assert'd against the baked DBARTS_C_API_HASH in
C_interface.cpp. inst/tinytest/capi/consumer.c consumes via the stubs
(one raw-GetCCallable canary retained). Gates: full tinytest, tests/cpp
from make clean, equivalence suites byte-identical (pure plumbing - no
sampler code moves), ABI checklist step satisfied by this doc.

C2 - stan4bart adoption (lockstep, after C1 lands): define
DBARTS_USE_STUBS, delete BARTFunctionTable + its typedefs + the 22
lookups, route the call sites through the stub names, handshake becomes
major-equality + minor-floor. Gates: full testthat, R CMD check tarball.

## Landing notes

C1 LANDED (the commit carrying this note): the 35-entry
DBARTS_C_API_LIST (33 existing + the major/minor accessors), four-field
entries (return, name, named parameter list, forwarding argument list);
Doxygen prototypes kept, mutually exclusive per TU with the same-name
DBARTS_USE_STUBS cached-pointer stubs (strict-C99 void dispatch, all
helper macros #undef'd after expansion); provider-side binding
static_asserts inside extern "C"; registration generated from the list;
dbarts_results offsets pinned per field; FNV-1a token over the
stringized list static_assert'd against the baked DBARTS_C_API_HASH;
version ships major 1 / minor 0 with DBARTS_C_API_VERSION packed as
major*1000+minor (value 1000) and dbarts_apiVersion retained so the
pinned stan4bart tip builds and loads unchanged until C2.
consumer.c consumes through the stubs with one deliberate raw
R_GetCCallable canary. Gates: R CMD INSTALL --preclean clean; tinytest
2926/0 (2924 + 2 version assertions); tests/cpp from clean binaries;
orchestrator re-ran all three equivalence suites - 22/22 identical
draws, BCF 5x6 bitwise, multinomial 2x4 bitwise; air format clean.
Anti-drift proofs fired live: stale hash bake, signature drift, and
mid-struct insertion each fail dbarts's own compile.

C2 LANDED (stan4bart d05a748): DBARTS_USE_STUBS on PKG_CPPFLAGS
(Makevars.in + Makevars.win; the generated Makevars untouched);
BARTFunctionTable, its 22 typedefs, and the 22 R_GetCCallable/bit_cast
lookups deleted (init.cpp 1127 -> 1077 lines); all 39 call sites
through the same-name stubs; handshake now major-equality + minor-floor
with an informative two-pair message, run before any other entry point.
Gates: R CMD INSTALL --preclean clean, testthat 314/0, R CMD check
Status OK with zero NOTEs. End-to-end proof: a wrong-arity stub call
against the installed header fails to compile with a diagnostic tracing
to DBARTS_C_API_LIST - the getTrees class is closed. DESCRIPTION
already carries the dbarts (>= 1.0-0) floor in Depends and LinkingTo.
ARC COMPLETE; the recorded additive escalations (signature-tagged
names, consumer-visible hash) remain future options.

Pre-release clarification (2026-07-22, 876a339): a hash re-bake before
the 1.0-0 submission DEFINES DBARTS_C_API_HASH; it does not consume a
MINOR bump. MAJOR/MINOR/VERSION stay 1/0/1000 across such a re-bake.
The MINOR-bump discipline (a signature-affecting change moves the hash
AND bumps MINOR, DESCRIPTION floor moved in lockstep) applies
POST-release, once some MINOR value has actually shipped to CRAN and a
consumer could be pinned against it. Before that there is nothing to
protect: stan4bart is the sole LinkingTo consumer and already
recompiles against dbarts in lockstep, every commit. Worked instance:
876a339 renamed setTestPredictors/setTestOffset/predict's x_test/
offset_test parameters to xTest/offsetTest - not an ABI change
(parameter names are not part of the call contract) but they feed the
stringized-signature FNV-1a, so DBARTS_C_API_HASH re-baked
0xc82cf27acefa5b81 -> 0xf760898d116cb3a3 with MAJOR/MINOR/VERSION held
at 1/0/1000.
