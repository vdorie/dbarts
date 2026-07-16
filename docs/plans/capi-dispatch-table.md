# capi-dispatch-table

status: memo revised after adversarial review, awaiting VD sign-off
(design pass 2026-07-15, adversarial second opinion 2026-07-16; no
implementation, no ABI change proposed to land without VD's fork decision
below).

## The question

VD's ask, verbatim intent: is it worth breaking with R's mechanism for
passing C function pointers by symbol (R_RegisterCCallable /
R_GetCCallable per function) and instead registering a single function
that returns a function table, indexed by a version string?

Motivating incident: the getTrees ABI break (2e2b1c9, restored a73ca50).
dbarts_sampler_getTrees gained a trainingData parameter under an unchanged
symbol name and an unchanged DBARTS_C_API_VERSION. LinkingTo consumers call
through their OWN declaration of the signature, so stan4bart passed eight
arguments to a nine-argument function; the garbage ninth became the
replay-predictor pointer and corrupted extracted node counts
stack-dependently (failed ~half of R CMD check runs for a week). Not a link
error, not a load error - call-time stack garbage.

## What the ground truth already is

Registration (src/R_interface.cpp:299-361): ~33 entry points in a
C_callMethods table, each R_RegisterCCallable'd under its own name in
R_init_dbarts. Append-only-by-name discipline; a signature change is an ABI
event (docs/plans/README.md:80-86, the a73ca50 review-checklist step).

Consumer (stan4bart, the SOLE reverse-LinkingTo package - CRAN confirms;
every other revdep is R-level, and stan4bart ships lockstep, same
maintainer):
- It already collapses the surface into ONE table struct: BARTFunctionTable
  (src/init.cpp:56-90), 22 of the ~33 symbols, populated by 22
  R_GetCCallable calls in lookupBARTFunctions (:1001-1029), with 39
  invocation sites routed through bartFunctions.*. So "one lookup + a
  table" is a one-time ~21-line saving in the one consumer, not a per-site
  win.
- It ALREADY does a load-time version handshake: apiVersion() !=
  DBARTS_C_API_VERSION errors before any pointer is used (:1003-1006).
- Each signature is stated THREE times on the consumer side even though it
  #include <dbarts/dbarts.h> (:23): the header prototype, the struct field
  type, and the bit_cast target type at the lookup. Neither of the latter
  two is derived from the first; a rebuild recompiles them verbatim. The
  getTrees skew lived in exactly that gap - the key mechanical fact for
  the analysis below.

Growth precedent already shipped in this tree: dbarts_results is a
structSize-growable output struct (dbarts.h:83-102, DBARTS_RESULTS_HAS) -
old callers are never written past, appends never reorder. The state format
reads its blocks BY NAME behind an encoding floor (dbarts.h:237-243,
c-api-growth.md part 2). Both are self-describing-DATA designs; growth of
the FUNCTION surface is already handled by append-only names (a future
run2-style addition is just a new name, load-detectable when absent). A
dispatch table buys nothing for growth that these do not already provide.

## The decisive axis: mechanical vs disciplinary compatibility token

The first draft claimed "the table relocates the hazard, it does not
remove it." False as literally written. The honest statement:

- A HAND-VERSIONED table (a human bumps DBARTS_C_API_VERSION, get_api vets
  it) cannot catch change-without-bump - the token did not move, so no
  handshake anywhere can see the change. Same exposure as per-symbol.
- An AUTO-TOKEN table CAN catch it: both sides generate, at their own
  build time, a hash of the normalized declaration set (signatures +
  struct layouts); get_api compares the consumer's baked-in hash at load
  and refuses a mismatch. No human bump required - change-without-bump
  becomes load-detectable in the field, mechanically.

So the decisive axis is not dispatch shape (per-symbol vs table) but
whether the compatibility token is DISCIPLINARY (a human bumps an integer)
or MECHANICAL (derived from the declarations themselves). The auto-token
table loses anyway, on two grounds: timing - a CI guard computing the same
normalized-declaration comparison catches the same class PRE-MERGE, before
any consumer can be built against the bad header, rather than at load in
the field - and convention - it diverges from the documented R mechanism
for zero residual safety once the CI guard exists. Dominated, not wrong.

## The rebuilt-consumer dimension (where the first draft's matrix erred)

The incident was a REBUILD scenario: stan4bart was recompiled against the
nine-arg header repeatedly during that week (R CMD check builds it). The
operative staleness was not an old binary - it was the HAND-ROLLED
declaration (init.cpp:73-75, :1017), which no rebuild re-derives. Had the
call at init.cpp:531 gone through a type stated IN THE HEADER - a table
struct's field, or an inline stub - the rebuild would have produced a
COMPILE error (eight args passed to a nine-arg type). So for the motivating
incident, table-struct-in-header and header stubs are safety-EQUIVALENT:
both catch it at compile time; only hand-rolled declarations yield garbage.
The table must be rejected on the honest grounds - R-convention divergence
plus a registration rewrite for zero safety gain over stubs - not on a
false safety claim.

## Skew matrix

Designs: (a) status quo - per-symbol, hand-rolled consumer declarations,
hand-bumped version; (b) get_api table, struct in dbarts.h, hand-versioned,
structSize-growable, NULL on unknown version; (c) per-symbol + inline stubs
in dbarts.h + mechanical CI guard. The auto-token table is covered in prose
above (it converts the not-rebuilt no-bump row to LOAD). Outcomes: COMPILE
(building the consumer), LOAD (before any BART call), CALL (stack garbage
at call time).

  scenario                        (a) status quo  (b) table       (c) stubs+guard
  -------------------------------------------------------------------------------
  NOT-REBUILT consumer (binary built against header vN, library now vM):
  needs entry old lib lacks       LOAD            LOAD (missing   LOAD (stub's
  (vN > vM)                       (GetCCallable   get_api or NULL GetCCallable
                                  errors)         for unknown vN) errors)
  ADD only, append-only kept      SAFE            SAFE (past old  SAFE
                                                  structSize)
  CHANGE, version NOT bumped      CALL            CALL (same      row PREVENTED:
                                                  slot, stale     CI fails the
                                                  cast)           merge un-bumped
  CHANGE, version bumped          LOAD            LOAD (welded:   LOAD
                                  (handshake)     NULL is forced) (handshake)
  -------------------------------------------------------------------------------
  REBUILT consumer (recompiled against current header, own decls as-is):
  CHANGE, bumped or not           CALL <-- the    COMPILE (header COMPILE (stub
                                  incident; hand- struct field is signature is
                                  rolled decl not new type; call  new type)
                                  re-derived      site mismatches)

Reading. The not-rebuilt CHANGE-without-bump row is CALL garbage under BOTH
(a) and (b) - only a mechanical token (auto-token at load, CI guard
pre-merge) touches it. The rebuilt row is where (b) and (c) both win and
(a) alone loses - and it is the row that actually happened. (b)'s only
edge over (a) that (c) lacks is welding the consumer-side NULL check to
table access; the sole consumer already performs that check voluntarily.

One trap the matrix exposes: the CI guard ALONE does not close the rebuilt
row. Trace: the guard forces the bump; a not-rebuilt consumer trips the
handshake at load - good. The maintainer rebuilds; versions realign; the
handshake now PASSES; but if the hand-rolled declaration was not also
widened (the original mistake - a header bump does not touch init.cpp),
the call is garbage again, this time WITH a passing handshake. The guard
protects the mismatch window only. Stubs remove the root cause: consumer
declarations become header-derived and are re-derived by every rebuild.
Hence stubs are load-bearing in the recommendation, not optional polish.

## Sub-choices, if a table were chosen anyway

- version string vs numeric: numeric major.minor. Strings invite bespoke
  parsing; integer major (incompatible) / minor (additive) maps onto the
  append-only rule and Matrix's integer ABI version.
- frozen table per major vs structSize-growable: growable, matching
  dbarts_results; frozen-per-major cannot absorb an additive minor without
  reverting to per-symbol registration for every in-major addition.
- unknown version: NULL, never best-effort - a table shaped for the wrong
  version re-opens call-time garbage.
- the struct lives in dbarts.h, which is also its limit: a NOT-rebuilt
  consumer still calls through a stale copy of it (the matrix row above).

## Costs

- (a) status quo: zero code. The CHANGE-without-bump row and the rebuilt
  row are both guarded only by the a73ca50 human checklist.
- (b) table: registration rewrite in dbarts, consumer rewrite in stan4bart,
  an idiom R does not document and no major LinkingTo package uses,
  permanent convention divergence. Zero safety gain over stubs on the
  rebuilt row; gains over status quo only the welded NULL check.
- (c-i) CI guard: a comparison script + committed golden + checklist line.
  NOT a raw text diff - a73ca50's header diff was ~3 declaration lines
  amid ~8 doc-comment lines, so naive diffs false-positive on comment
  edits; it needs comment/whitespace-stripped, symbol-name-keyed
  declaration comparison covering function signatures AND dbarts_results
  layout, with a mid-struct insertion flagged as a field-ORDER change,
  never waved through as an addition. Its correctness depends on every
  ABI-crossing type staying defined inside dbarts.h - currently true
  (dbarts_results and the callback typedef are inline; control/model/data
  cross as SEXP) - and it goes blind if an ABI type ever moves to another
  header (as would a table).
- (c-ii) stubs, both ledger columns. dbarts pays: ~33 inline wrappers,
  roughly 100-130 lines of shipped header, each signature stated twice in
  the header (prototype + stub), a stub naming choice (a prefix or macro
  gate, since LinkingTo consumers cannot link the real symbols), and a
  second documented consumption mechanism to teach alongside raw
  R_GetCCallable. The consumer gains: BARTFunctionTable, its 22 typedefs,
  and the 22 bit_cast lookups deleted (~90 lines); three signature
  statements per entry point drop to zero consumer-side; the incident
  class becomes a compile error on every rebuild.

## Precedents

- Matrix - the largest LinkingTo C API in R - is per-function
  R_GetCCallable through an inline stub layer (Matrix_stubs.c,
  R_MATRIX_INLINE, e.g. M_cholmod_start), NOT a get-api table, with a
  numeric ABI version (inst/include/Matrix/version.h) bumped 1 -> 2 at
  1.7-0 and a documented rebuild-reverse-LinkingTo-on-bump rule. This is
  recommendation (c) almost verbatim.
- Writing R Extensions documents only per-symbol R_GetCCallable; no major
  CRAN package was found shipping a get-api table.
- Counter-precedent, conceded: numpy's PyArray_API IS a function-pointer
  table - import_array() fetches it through a capsule and performs a
  runtime ABI check. Tables are proven outside CRAN. But numpy's safety
  comes from append-only slots + a handshake + SINGLE-SOURCE CODEGEN
  (generate_numpy_api.py emits both sides from one spec, so dual
  declarations cannot exist) - the same ingredients recommended here in
  per-symbol form. It reinforces that the ingredients matter, not the
  shape, and names a fourth option - generating consumer declarations
  from dbarts.h - which stubs approximate at a fraction of the machinery.
- In-tree: dbarts_results structSize and the state format's by-name blocks
  are both self-describing-DATA precedents; neither argues for a
  self-describing FUNCTION table, since function growth is already
  append-only-by-name and load-detectable.

## Recommendation (the fork is VD's)

Do NOT adopt the get_api table. A hand-versioned table adds no safety on
either dangerous row; an auto-token table does add field-level detection
but is dominated by a CI guard that mechanizes the same comparison
pre-merge without diverging from the R convention; and on the rebuilt row -
the incident that motivated the ask - a table is exactly equivalent to
stubs, which do not disturb the dispatch mechanism.

Keep per-symbol R_GetCCallable and decide two ORTHOGONAL questions:

- DECISION 1, compatibility token (recommend: adopt regardless of Decision
  2). Keep hand-bumped DBARTS_C_API_VERSION but make it mechanical: a CI
  guard that fails the build whenever the normalized declaration set of
  dbarts.h (function signatures and dbarts_results layout, field-order
  sensitive) changes without a same-diff version bump. Document the
  semantics - integer major = incompatible, minor = additive - and move
  stan4bart's handshake from != to >= major before the freeze.
  Adopt: a script, a committed golden, a checklist line. Skip: the
  CHANGE-without-bump row stays guarded by human discipline alone.
- DECISION 2, consumer declaration source (recommend: adopt - the
  root-cause fix). Ship Matrix-style inline stubs in dbarts.h and delete
  stan4bart's BARTFunctionTable; consumer declarations become
  header-derived, re-derived on every rebuild, making the getTrees class
  a compile error and structurally protecting any future second consumer
  the 1.0-0 freeze will bind.
  Adopt: the header cost above (~33 wrappers, ~100-130 lines, dual
  statement, second documented mechanism). Skip: rebuilds never re-derive
  hand-rolled declarations, so the incident class stays reachable even
  with Decision 1 (the realigned-handshake trace above).

Timing: dbarts 1.0-0 freezes the dbarts.h compatibility contract at
submission. Both decisions are cheap now and awkward-to-impossible after;
the version semantics (integer major = incompatible) should be written
down before release whichever way the fork goes, because that is what
gives any handshake teeth.
