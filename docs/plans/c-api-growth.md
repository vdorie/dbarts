# c-api-growth

agent: opus
rng: neutral - this is pure output-plumbing (a new run entry point) and a
  relaxation of the state reader's version gate. No draw, no accumulation,
  no RNG-consuming path is touched; every existing run and every existing
  state load produces bit-identical results. Equivalence is byte-identical
  (21/21), tinytest snapshots unchanged.
budget: dbarts.h ~35 lines (struct + macro + one prototype + version bump);
  C_interface.cpp ~60 lines (run2 core + run shim + static_asserts);
  R_interface_bartcore.cpp ~40 lines (setState/installForests gate swap +
  block-name refusal helper); R_interface.cpp 2 registration lines;
  inst/tinytest/capi/consumer.c ~70 lines (capi_run2 + guard canary);
  inst/tinytest/test-capi.R ~40 lines; one new state tinytest ~50 lines.

## Status: LANDED (bartcore, 2026-07-10)

**Summary.** This plan makes dbarts's two 1.0-0-bound compatibility surfaces
extensible without future ABI breaks.

Output channels: `dbarts_results` (inst/include/dbarts/dbarts.h) gained a
caller-set leading `structSize` field, and `dbarts_sampler_run` now fills
each output pointer only when it is both present-by-size (its offset+size
falls within the caller's declared `structSize`) and non-null - so fields
can append monotonically across releases and the library never writes past
an older caller's stack allocation. `DBARTS_C_API_VERSION` stays 1: the
struct is born size-first before its first release, so there is no prior
layout to version against. (A same-day initial proposal instead kept
`dbarts_results` byte-frozen forever and added parallel
`dbarts_sampler_run2`/`dbarts_results2` symbols; VD superseded it once it
was clear nothing on CRAN yet links the unreleased header, so the frozen
struct could simply be redesigned directly instead of shimmed - full record
and rationale in Design, Part 1.)

State format: `setState`/`installForests` were already reading per-chain
blocks by name; this plan demotes the `formatVersion` stamp from a strict-
equality gate to a `minReadableStateFormatVersion` encoding floor plus
provenance, and makes REQUIRED-block refusals name the missing or malformed
block. Additive block additions no longer orphan every state saved by a
prior release.

Both halves landed on bartcore 2026-07-10 as two independent commits. A
deferred first real append to the output channel - `logLikelihood`, feeding
the pointwise-loglik consumer - landed the same day as a follow-up, once
its engine-side scope was cut to fit budget (the initial landing shipped
mechanism-only). Full chronology, gate results, and the exact deviation
rationale are in the Landing log at the end of this document.

## Goal

Make the two frozen dbarts.h surfaces extensible without ABI breaks, per
the DECIDED design (VD, 2026-07-10, TODO c-api-growth, window pre-release):

1. Output channels: `dbarts_results` becomes size-first - a caller-set
   `structSize` as its leading member - so fields can append monotonically
   across releases while the library never writes past the caller's
   declared size. `dbarts_sampler_run` absorbs the guarded write core
   directly; there is no second entry point or second struct. (Originally
   proposed as a separate `dbarts_sampler_run2`/`dbarts_results2` pair that
   left `dbarts_results`/`dbarts_sampler_run` untouched behind a shim;
   changed the same day once it was clear the unreleased header has no
   installed base to protect - see Design, Part 1, for the full original
   design and the rationale for dropping it.)
2. State format: additive by-name blocks. Stop the every-encoding-change-
   orphans-all-states behavior by having the reader look up blocks by name,
   default absent OPTIONAL blocks, and refuse - naming the block - only when
   a REQUIRED one is missing. The `formatVersion` stamp is retained but
   demoted from an equality gatekeeper to an encoding-compatibility floor
   plus provenance. Absorbs the flat-format-v2 TODO item.

Land before 1.0-0 submission: `dbarts.h` and the state format both freeze
their compatibility contract at that release, so the extension mechanism
must be in place first (a `dbarts_results` grown after 1.0-0 breaks
stan4bart silently; a state gate left at strict-equality orphans every
saved fit on the next feature).

## Context

### The frozen surface, precisely

`dbarts_results` (inst/include/dbarts/dbarts.h:63-72) is eight caller-owned
output pointers, NULL-means-skip. It is an OUTPUT struct only. There is no
size or version field. `dbarts_sampler_run` (dbarts.h:89-90) takes a
`dbarts_results*`.

The freeze is real and demonstrable: the shipped C consumer stack-allocates
`dbarts_results results = {0};` and fills fields
(inst/tinytest/capi/consumer.c:204-210, :282-290, :335-343). `sizeof(dbarts_results)`
is baked into stan4bart's compiled object; if the library appended a ninth
field and wrote through it on a stan4bart-sized struct, it would scribble
past the caller's stack allocation. This is exactly why tau/groupEffects
(added for grouped ranef) forced stan4bart to re-value-initialize and
reinstall (TODO release note L447-449) - and that was the LAST field the
frozen struct can absorb for free, because those consumers were rebuilt in
lockstep. After 1.0-0 there is no lockstep rebuild for a CRAN release.

The mapping run does today (src/C_interface.cpp:60-92): it copies the eight
`dbarts_results` pointers into a stack `bartcore::Results` (the engine struct,
src/bartcore/chain.hpp:169-181, whose members are `sigma`, `trainingFits`,
`testFits`, `variableCounts`, `k`, `splitProbabilities`, `tau`,
`groupEffects` - same eight, engine names), brackets `GetRNGstate`/`run`/
`PutRNGstate`, and wires the optional per-sweep callback. The engine
`Results` already grows freely (it is internal, recompiled per build); the
ABI wall is only the public `dbarts_results`.

Entry points register two ways in src/R_interface.cpp: the `.Call` table
(R-facing `bartcore_*`) and the CCallable table (dbarts.h symbols, resolved
by consumers through `R_GetCCallable`). `dbarts_sampler_run` is at
R_interface.cpp:301; the CCallable block runs :295-345, registered in
`R_init_dbarts` :349-357. `DBARTS_C_API_VERSION` is 1 (dbarts.h:48), returned
by `dbarts_apiVersion()` (C_interface.cpp:43); the header instructs consumers
to check it before using the rest (dbarts.h:9-13).

### The state format, precisely

`storeState` (R_interface_bartcore.cpp:2787-2948) marshals a
`bartcore::SamplerStateData` (src/bartcore/sampler.hpp:63,
chain.hpp:196-235) into an R list, one element per chain, plus top-level
attributes. `setState` (R_interface_bartcore.cpp:2950-3181) reverses it.

Top-level attributes on the state list:
- `cutPoints` (REQUIRED; :2936, read :2981-2998, refuses if absent/wrong-length)
- `currentSampleNum` (REQUIRED; :2937, read :3000-3008)
- `formatVersion` (:2939, integer, currently == 3)
- `packageVersion` (:2941, provenance string, PACKAGE_VERSION "1.0-0",
  src/config.hpp:97)
- `class` = "bartcoreState" (:2943)

The KEY fact for this job: `setState` already reads every per-chain slot BY
NAME through `getListElement` (R_interface_bartcore.cpp:83-90, a linear
names scan), NOT by positional index, and already defaults absent OPTIONAL
slots. The write side names slots via a fixed enum only for its own
convenience (:2810-2821); the read side is name-driven and enum-free.

Per-chain slots, classified as the reader treats them today:
- `forests` (REQUIRED; :3015-3019 refuses if null/non-list). Each forest:
  - `tree.vars`/`tree.values`/`tree.sizes`/`tree.flags` (REQUIRED via
    readFlatTrees, :3029-3034)
  - `tree.params` (CONDITIONALLY REQUIRED: iff `usesFunctionLeaves()` or
    `numLeafCovariates() > 0`, :3046-3066)
  - `tree.masks` (CONDITIONALLY REQUIRED: iff `data().hasPooledCategorical`,
    :3069-3078)
  - `saved.*` (OPTIONAL: `saved.sizes` presence gates the block, :3035-3042;
    present iff tree storage was on)
  - `k` (REQUIRED, :3080-3085)
- `sigma` (REQUIRED, :3089-3094)
- `fit.scale` (REQUIRED, length 2, :3096-3102)
- `latents` (OPTIONAL: `if (!Rf_isNull(...))`, :3104-3112; present iff binary)
- `ranef` + `tau` (OPTIONAL pair: gated on `ranef` non-null, :3114-3125)
- `dart.probabilities` + `dart.alpha` + `dart.updates.skipped` (OPTIONAL
  triple, :3127-3145)
- `rng.state` (OPTIONAL, :3147-3157; absence forfeits bitwise continuation
  only - dbarts.h:194-196 already documents cross-kind unrestored streams)
- `bcf` (OPTIONAL, :3159-3170; present iff BCF)

So the additive-by-name reader is ~90% already built. The ONE thing that
orphans states across releases is the strict-equality gate at the top of
`setState` (:2959-2969) and `installForests` (:3337-3340):
`if (formatVersion != stateFormatVersion) Rf_error(...)`. Every past feature
bumped `stateFormatVersion` (now 3; history at :2773-2785:
v2 = flat-node tagging + dropped-slots ENCODING change, v3 = trees-into-
forests + bcf slot STRUCTURAL change) and thereby invalidated every state a
prior release wrote. public-surface.md records this gate as the "one
version scheme flat-format-v2, state-continuation, and forest-split-bcf all
bump" - c-api-growth is the item that ends the bump-orphans-everything half
of it while keeping the never-silently-misinterpret half.

### Why pre-1.0 states are not a compat target

The 1.0-0 cutover reset the compat contract (chi-k-runaway.md,
public-surface.md, MEMORY no-backwards-compat-constraint). Every state
in the wild after release is format 3 (the 1.0-0 encoding). v1/v2 states
cannot even structurally reach the v3 reader (they lack the `forests` slot,
so the REQUIRED-block check refuses them anyway). This is what makes it safe
to relax the gate: the only encodings the relaxed reader must never misread
are hypothetical FUTURE ones, and the registry rule below prevents those by
construction.

### Near-term output-channel consumers (survey)

TODO + docs/design, filtered to "named consumer or decided item AND the
engine can produce it now":
- pointwise-loglik (LANDED 2026-07-10, TODO L318-328): computes per-obs
  log-lik R-side from stored draws; its note "The engine-side Results field
  version rides c-api-growth" is the one DECIDED item that explicitly names
  this job for its channel. Producible now (engine has fits + sigma at
  storeSample). This is the single defensible demonstrator channel.
- negative-binomial (TODO L305-307): needs an `r` shape trace. Model NOT
  built - no engine to fill it. Future append; used as the worked example.
- ordinal-outcomes (TODO L308-309): needs a cutpoint block. Model NOT built.
  Future append.
- per-forest fits (forest-combiner/multi-forest, TODO L60-71, L122-130):
  BCF per-forest fits already ship via SEPARATE R-level entry points
  `bartcore_getForestFits`/`bartcore_getBCFGlue` (R_interface.cpp:229-230),
  not through `dbarts_results`. forest-combiner landed 2026-07-14
  (docs/design/forest-combiner.md), so folding them into `dbarts_results` is
  unblocked; do it when a consumer pulls. Future append.

JUDGMENT CALL (named): whether the growable `dbarts_results` ships with the
loglik channel at introduction or is mechanism-only (size + the frozen
eight, first real append lands with its model). Recommendation:
mechanism-only + loglik as the FIRST append IN THE SAME PR, so the
write-guard is exercised end-to-end by a real channel in the consumer test
rather than only by a canary. Rationale: adding a channel the engine cannot
yet fill (NB r, cutpoints) would be dead surface (violates "only what has a
consumer"); loglik is the one channel with both a decided item and a live
producer, so it doubles as the guard's live test. VD to confirm loglik-now
vs mechanism-only.

## Design

### Part 1: dbarts_results and dbarts_sampler_run (the growable output channel)

**Adopted design.** Header additions (inst/include/dbarts/dbarts.h): give
`dbarts_results` a leading caller-set `structSize`, ahead of the existing
eight output pointers, with a boundary comment marking where the 1.0-0
layout ends and future appends go. The library fills a field only when it
is both present-by-size and non-null:

```c
/// Caller-owned, growable results. structSize MUST be set by the caller to
/// sizeof(dbarts_results) as the caller compiled it; the library fills only
/// fields whose end offset is within structSize, so a caller built against
/// an older header (smaller struct) is never written past. Fields append
/// monotonically and never reorder across releases; a field is filled only
/// when both present-by-size and non-null. Value-initialize and set the
/// size: dbarts_results r = {0}; r.structSize = sizeof r;
typedef struct dbarts_results_t {
  size_t structSize;    ///< caller sets = sizeof(dbarts_results)
  double* sigma;        ///< numSamples x numChains
  double* train;        ///< numObservations x numSamples x numChains
  double* test;         ///< numTestObservations x numSamples x numChains
  uint32_t* varcount;   ///< numPredictors x numSamples x numChains
  double* k;            ///< numSamples x numChains
  double* varprobs;     ///< numPredictors x numSamples x numChains
  double* tau;          ///< numSamples x numChains
  double* groupEffects; ///< numGroups x numSamples x numChains
  /* --- 1.0-0 field boundary; appends after this release bump
     DBARTS_C_API_VERSION --- */
  double* logLikelihood; ///< numObservations x numSamples x numChains (opt.;
                         ///< landed as the first real append - see the
                         ///< Landing log)
} dbarts_results;

/// True when the caller's struct (per structSize) actually carries `field`.
/// sizeof is unevaluated, so this never dereferences past the caller buffer.
#define DBARTS_RESULTS_HAS(r, field) \
  ((r)->structSize >= offsetof(dbarts_results, field) + sizeof((r)->field))
```

Because nothing on CRAN links the unreleased 1.0-0 header yet,
`dbarts_results` had no installed base to protect: it could be redesigned in
place rather than frozen-and-shimmed. `DBARTS_C_API_VERSION` stays 1 - the
struct is BORN size-first, so there is no prior layout to distinguish;
future appends (after the 1.0-0 release ships) bump it from here. Adding a
new CCallable name would not by itself require the bump, but the bump is
the sanctioned "these fields/entry points exist" signal. `dbarts_sampler_run`
absorbs the guarded write core directly: there is no second entry point and
no second struct.

The write-guard idiom (library side, C_interface.cpp). A field is filled iff
present-by-size AND non-null:

```c
#define FILL(field, engineMember) \
  engineResults.engineMember = \
    (DBARTS_RESULTS_HAS(results, field) ? results->field : NULL)
```

`offsetof` is against the LIBRARY's (newest) layout; because fields only
append and never reorder, that offset equals the offset the caller would
compute for the same field, and `structSize` bounds the caller's buffer. A
smaller (older) caller struct makes `DBARTS_RESULTS_HAS` false for every
field past its end, so the library never reads the (absent) pointer slot nor
writes through it. A larger (newer-than-lib) caller struct is also safe: the
library only references fields it has code for (all within its own sizeof,
hence <= structSize), and value-init leaves the caller's newer fields NULL/0.

`dbarts_sampler_run` (C_interface.cpp) now runs the eight/nine-pointer map +
callback wiring + RNG bracket through this guarded fill for every field, old
and new alike - a caller that still declares only the frozen eight
(structSize ending at `logLikelihood`) gets exactly the historical behavior,
because every field past its declared size guards off.

Compile-time layout locks (C_interface.cpp, C++ side, NOT the shipped header
- the header must stay C99-clean for C consumers like consumer.c):
```cpp
static_assert(offsetof(dbarts_results, structSize) == 0);
static_assert(offsetof(dbarts_results, sigma)  == sizeof(size_t) /*+pad*/);
static_assert(offsetof(dbarts_results, sigma) < offsetof(dbarts_results, train));
// ... one monotonic-order assert per field pair ...
static_assert(offsetof(dbarts_results, groupEffects) <
              offsetof(dbarts_results, logLikelihood));
static_assert(sizeof(dbarts_results) == EXPECTED); // bump deliberately on append
```
The per-field-offset and final-sizeof asserts are the layout lock: any
reorder or mid-struct insertion fails the build, and appending forces the
author to update EXPECTED (an explicit acknowledgement that the ABI grew).

JUDGMENT CALL (named): pin exact literal offsets vs only-monotonic asserts.
Recommendation: monotonic-order asserts + one exact `sizeof` assert - exact
per-field offsets are padding/ABI-dependent and add churn without extra
safety once order + total size are pinned.

Thread-safety / ABI notes: main R thread only (RNG bracket internal),
callback refused while chains run on worker threads (C_interface.cpp:74-87).
No new global state. The struct is caller-owned and single-threaded per
call.

Registration: `dbarts_sampler_run`'s existing CCallable entry
(R_interface.cpp near :301) needs no new registration line - only its
signature (the now size-first `dbarts_results*`) changed.

**Originally proposed, and superseded: a parallel dbarts_sampler_run2 /
dbarts_results2.** Drafted 2026-07-10 as the initial plan, and landed as
Commit A's initial content (see Landing log), before VD's same-day design
change superseded it. The idea: leave `dbarts_results` and
`dbarts_sampler_run` byte-frozen forever, and add new-named siblings
`dbarts_results2` (`structSize` as its first member, the same eight
pointers, one appended field) and `dbarts_sampler_run2` (taking
`dbarts_results2*`). `dbarts_sampler_run` would become a thin shim: both
`run2` and the shim would call a shared static
`runCore(sampler, burnIn, numSamples, dbarts_results2*)` holding the
eight-pointer map, callback wiring, and RNG bracket. `run2` would call
`runCore` directly; `run` would build a local
`dbarts_results2 view = {0}; view.structSize = offsetof(dbarts_results2,
logLikelihood);` (pinned to the v1 boundary - the size an all-v1 caller
would declare), copy the eight frozen pointers field-by-field (not
`reinterpret_cast`, to keep `dbarts_results` and `dbarts_results2`
layout-DECOUPLED so v1 stays genuinely frozen and independent), and call
`runCore` - reproducing `run`'s exact behavior because every post-v1 field
then guards off. `run2` would share `run`'s contract verbatim: main R
thread only, callback refused on worker threads, no new global state.
`DBARTS_C_API_VERSION` would bump 1 -> 2 so a consumer could gate
`dbarts_apiVersion() >= 2` before resolving `run2` (the header already
mandates the apiVersion check, dbarts.h:9-13); `run2` would need its own
CCallable registration line (`DEF_FUNC("dbarts_sampler_run2",
dbarts_sampler_run2)`, R_interface.cpp near :301).

Changed because: nothing on CRAN links the unreleased 1.0-0 header, so
`dbarts_results` had no installed base to protect with a shim - it could
simply become size-first itself, and the one stan4bart lockstep rebuild
already planned for the 1.0-0 cutover absorbs the layout change for free.
Carrying `run2`/`results2` forward as a permanent parallel surface would
have meant a shim and a second registered symbol forever, for a
pre-vs-post-growable distinction that only matters before the struct ships
at all. VD made the call 2026-07-10 - the same day as the initial landing.
Commit B (state blocks, independent of this surface) was unaffected and
cherry-picked onto bartcore verbatim; Commit A was reworked in place to the
adopted design above before landing. The write-guard idiom, the
static_asserts (structSize at offset 0, monotonic order, one exact sizeof),
the boundary comment, and the mechanism-only adjudication (Landing log)
all carried over unchanged into the reworked `dbarts_sampler_run` - only
the duality itself (separate name, separate struct, version bump, second
registration) was dropped.

### Part 2: additive by-name state blocks

The reader is already by-name; the deliverables are (a) demote the version
gate, (b) name the missing block on refusal, (c) write down the registry
rule that keeps the relaxation safe.

(a) Version gate -> encoding floor + provenance. Introduce alongside
`stateFormatVersion` (R_interface_bartcore.cpp:2785):
```cpp
// The oldest ENCODING this reader still understands. Additive block
// additions do NOT raise it (they are read by name and defaulted when
// absent); only a non-additive change to an existing block's on-disk
// encoding raises it. Currently 3 (the 1.0-0 encoding); pre-1.0 states are
// not a compat target.
static const int minReadableStateFormatVersion = 3;
```
Replace the two `!= stateFormatVersion` gates (setState :2959, installForests
:3337) with `< minReadableStateFormatVersion`, keeping the packageVersion in
the message:
```cpp
if (formatVersion < minReadableStateFormatVersion)
  Rf_error("state encoding version %d (written by dbarts %s) predates the "
           "oldest this dbarts (%d) can read; re-fit with this version",
           formatVersion, packageVersion, minReadableStateFormatVersion);
```
`storeState` keeps stamping `formatVersion` (provenance + the floor input).
The stamp is bumped ONLY on a non-additive encoding change (which also bumps
the floor); additive block additions leave it at 3. A future release reading
a state written by an even-newer additive release passes the floor (>= 3) and
simply never looks up the unknown block names (getListElement is queried, not
enumerated), so forward-compat holds without code. This is the concrete
meaning of "the version stamp stays for provenance, not gatekeeping": it no
longer gates additive evolution, only genuinely-incompatible old encodings.

JUDGMENT CALL (named): floor (recommended) vs dropping the version check
entirely. Pure-provenance (no floor) is what the DECIDED text literally says;
the floor is defense-in-depth against a pre-1.0 or future-incompatible
encoding reaching the by-name reader and being misread. Recommendation: keep
the floor - it costs one comparison, changes no supported behavior (all real
states are >= 3), and is the only thing standing between a mis-versioned blob
and a silent misread. Frame as: the floor IS the "refuse a genuinely
incompatible encoding" half; the by-name reader IS the "additive evolution
is free" half.

(b) Name the missing REQUIRED block on refusal. Today required-absence gives
generic messages ("malformed forests in bartcore state", :3017; "malformed
parameters", :3082/:3091). Add a two-message convention distinguishing absent
from malformed, and use the block name:
```
"bartcore state is missing required block '%s'"        // getListElement null
"bartcore state block '%s' is malformed"               // present, wrong type
```
Apply to the REQUIRED and CONDITIONALLY-REQUIRED slots (`forests`, the tree.*
channels, `k`, `sigma`, `fit.scale`, and - conditioned on sampler config -
`tree.params`/`tree.masks`, and `latents` for binary if you promote it; see
worked example). OPTIONAL slots keep their silent-default behavior unchanged.

(c) Registry convention (document in the storeState comment block at :2773
and in a new short section of docs/design/public-surface.md section 2 or a
state-format design note):
- Block names are APPEND-ONLY. Once a release ships a block under a name,
  that name's on-disk ENCODING is FROZEN.
- A NEW capability adds a NEW optional block name; the old reader ignores it
  (forward-compat), the new reader defaults it when absent (backward-compat).
- Changing an existing block's encoding is FORBIDDEN as an in-place edit;
  it must introduce a new name and (if old states must still load) keep the
  old-name reader. Only a change that cannot be expressed that way bumps the
  floor.
- A block REQUIRED for a given sampler configuration (latents for binary,
  tree.params for vector leaves, tree.masks for pooled categoricals, and by
  extension a future `nb` for negative-binomial) is checked with the
  existing conditional pattern and refuses naming the block when the sampler
  needs it but the state omits it.

### Worked example: a future negative-binomial block (both surfaces)

Output channel (`dbarts_results`). NB adds an `r` shape trace. It appends
BELOW logLikelihood:
```c
  double* logLikelihood; // (already shipped)
  double* rTrace;        // numSamples x numChains; filled only under NB
```
Writer (C_interface.cpp): `FILL(rTrace, rTrace);` maps to a new
`bartcore::Results::rTrace` the NB response model fills at storeSample.
Reader behavior across versions:
- Old consumer (structSize ends at logLikelihood): `DBARTS_RESULTS_HAS(r,
  rTrace)` is false -> library never touches it. Old consumer keeps working
  against the NB-aware library, just without the trace.
- New consumer, non-NB sampler: `rTrace` non-null but the engine's NB
  producer is inactive -> the channel is simply never written (existing
  null-skip-by-model convention, dbarts.h:59-62). Consumer value-init leaves
  it 0; it must gate on `dbarts_sampler_...` model queries as with k/varprobs.
- New consumer, NB sampler: filled. `sizeof(dbarts_results)` EXPECTED assert
  bumped; `DBARTS_C_API_VERSION` bumped to 2 (the first bump after the
  1.0-0 release ships - logLikelihood's pre-release append did not bump it).

State block (save/load of the NB shape, which is chain state). Writer adds an
optional per-chain slot:
```cpp
if (chainState.hasNB)
  SET_VECTOR_ELT(chainExpr, SLOT_NB, Rf_ScalarReal(chainState.r));
```
named `"nb"` (append the enum + name at :2810-2821; no floor bump - it is
additive). Reader:
```cpp
SEXP nbExpr = getListElement(chainExpr, "nb");
if (sampler.usesNegativeBinomial()) {
  if (Rf_isNull(nbExpr) || !Rf_isReal(nbExpr) || Rf_xlength(nbExpr) != 1) {
    errorMessage = /* "missing required block 'nb'" or "block 'nb' malformed" */;
    break;
  }
  chainState.r = REAL(nbExpr)[0]; chainState.hasNB = true;
} // else: not an NB sampler, ignore the slot entirely
```
Cross-version behavior:
- NB state loaded by an OLD (pre-NB) reader: `nb` never queried, silently
  dropped. Harmless unless the destination is itself NB (an old reader has no
  NB sampler, so it never is) - the conditional-required check only fires on
  a sampler that needs it.
- Non-NB state loaded by a NEW reader on a non-NB sampler: `nb` absent,
  `usesNegativeBinomial()` false -> defaulted, no refusal.
- NB state loaded by a NEW reader on an NB sampler: required, present,
  restored bitwise.
This is exactly the latents/tree.params/tree.masks conditional-required
pattern already in setState (:3046-3078, :3104-3112), so NB adds no new
machinery - it instantiates the registry rule.

### Reserved: a future dbarts_sampler_setForestWeights (2026-08-10, not built)

docs/plans/zero-weight-exactness.md (S2) adds a caller-settable per-forest,
per-observation weight to the engine (`Chain::setForestWeights`) and the
bartcore bridge (`bartcore_setForestWeights`), `dbarts:::`-only - no
`dbarts.h` symbol, because the flat API has no BCF creation entry point to
reach it from (`dbarts_sampler_create` only ever routes to `createHolder`,
never a BCF holder). Reserving the name and signature here so a later BCF
creation entry does not collide with it:

    X(int, dbarts_sampler_setForestWeights,
      (dbarts_sampler*, size_t forest, const double* weights),
      (sampler, forest, weights))

appended at the END of the X-list (dbarts.h); no version constant moves when
it lands (VD 2026-08-10: no increments pre-release - the constants stay 1.0
until the first release). Three constraints recorded now so a later
implementer does not relitigate them:

1. **Erratum (dbarts-h-reshape, 2026-08-10):** this reservation originally
   read "Nonzero return means refusal, matching `dbarts_sampler_setPredictor`".
   That is the INVERSE of the shipped convention - `dbarts_sampler_setPredictor`
   ends `return result == bartcore::PredictorUpdateResult::accepted ? 1 : 0`
   (C_interface.cpp), i.e. **1 = accepted, 0 = refused**, which is also what
   bcf-public-surface S3 item 3 records. The entry is built with
   1 = accepted / 0 = refused; `weights == NULL` clears the installed weight
   rather than refusing.
2. Per-forest weights must NOT be folded into a future BCF creation struct -
   membership is typically resampled per sweep, so this is a mutation entry
   point, not a creation-time parameter.
3. `dbarts_sampler_setWeights` must NOT be widened with a forest index to
   carry this instead: the two channels differ on the sigma degrees of
   freedom (the case weight counts toward it; the per-forest weight must
   not), so collapsing them would need a call-site flag anyway and buys
   nothing.

### Landed: the BCF flat surface (bcf-public-surface S3, 1622eb9, 2026-08-10)

S3 appended four entries at the END of the X-list - `dbarts_sampler_numForests`,
`setTreatment`, `forestFits`, `bcfGlue` (inst/include/dbarts/dbarts.h:264-271)
- and widened `dbarts_sampler_setResponse` to take `updateScale`, re-baking
`DBARTS_C_API_HASH` from `0xf760898d116cb3a3ULL` to `0x1a911c00bb26dcd7ULL`
(dbarts.h:83); both version constants stayed at 1 and 0 (no bartcore release
has shipped). The three `int` entries return 1 = accepted, 0 = refused,
matching `dbarts_sampler_setPredictor` and this reservation's own convention.

**Reserved for the dbarts.h reshape, not built there:** forest-indexed
`getTrees`, `printTrees`, `predict` and `numTrees` - today's tree-facing
queries are forest-blind, ambiguous on a two-forest sampler
(bcf-public-surface.md S3, "Coordination with the queued dbarts.h reshape";
that section also named `setTreeStorage`, whose forest-indexed form the
reshape plan since CLOSED BY FACT - storage is per sampler).
`docs/plans/dbarts-h-reshape.md` S1 item 3 builds the three forest-indexed
tree queries (`numTrees`, `getTrees`, `printTrees`); forest-indexed `predict`
stays a recorded door there.
That plan's S1 item 5b is a CONDITIONAL carve-out re-signing `setTreatment`/
`bcfGlue` themselves to `setForestBasis`/`numForestAmplitudes`/
`forestAmplitudes`, engine vocabulary, if
`docs/plans/multiforest-extension-surface.md`'s fork 3 is answered before
that slice starts (dbarts-h-reshape.md binding decision 3's carve-out);
otherwise the BCF names S3 shipped go to release unchanged.

### No X-list change: multiforest-predictor-mutation (2026-08-10 to 08-12)

`docs/plans/multiforest-predictor-mutation.md` (SL through S4) widened
`setPredictor`, `updatePredictor`, and the per-observation session to accept a
BCF, multinomial, and heteroscedastic sampler, and retired
`refuseMultiForestTransactionalUpdate` and
`refuseVarianceForestPredictorMutation` from the R bridge. None of that
touched `dbarts.h`: the four transactional entries already existed and only
their guard bodies changed. No X-list entry was appended and
`DBARTS_C_API_HASH` did not move, so the queued dbarts.h reshape is
unaffected by this arc.

### No X-list change: multinomial-counts-mutation (2026-08-12), reservations for the reshape

`docs/plans/multinomial-counts-mutation.md` (S1-S5) added a counts mutation
channel (`bartcore_setCounts`) and an n x K / nTest x K category offset
(train, test, and a per-call predict offset) to the multinomial sampler, all
`dbarts:::bartcore*`-only - `dbarts_sampler_create` has no multinomial branch
to reach a flat entry from, so no `dbarts.h` symbol was added and
`DBARTS_C_API_HASH` did not move. S5 lifted `bart2(family = "multinomial")`'s
own `offset` refusal for an n x K matrix (a one-shot creation-time argument,
not a mutation surface); still no flat surface. Three reservations for
whenever the dbarts.h reshape scopes a flat multinomial creation entry:

1. **Source-shaped response and offset parameters.** The reshape already
   re-signs predictor entries onto a borrowed `PredictorSource` view; the
   response side needs the same treatment before a flat multinomial entry can
   exist - a tagged source expressing at minimum `{ double* vector }` and
   `{ int* counts, size_t numCategories }` for the response, and
   `{ double* vector }` / `{ double* matrix, size_t numCategories }` for the
   offset. Built this way, a later flat multinomial creation entry needs a new
   tag, not an ABI break.
2. **A size-first spec struct**, per the `dbarts_results` `structSize`
   precedent (Part 1, above), for any future flat multinomial creation entry.
3. **Whatever tagged struct the reshape adopts must PRESERVE the
   `refuseMultiForestMutation` / `refuseMultiForestResponseMutation` guards**
   already on the flat response-side entries (the D4 precedent,
   `docs/plans/multiforest-mutation-gaps.md`) - the reshape is exactly where a
   guard like this gets dropped by accident, widening a flat entry to silently
   accept a multi-forest sampler it was never exercised against.

## Verification

Gates run from a worktree against a private library (per the repo's install
gotcha: `R CMD INSTALL --preclean -l <lib>` after touching dbarts.h and the
bridge; delete tests/cpp + capi binaries so no stale-header link).

1. tests/cpp: rebuild clean, `./test_bartcore` green. No new engine case is
   strictly required (the growable-results mechanism is bridge-level), but
   if the loglik channel ships, add nothing here - the engine `Results`
   growth is covered by the C consumer round trip. (If VD wants an
   engine-level assert, a test_state.cpp case that the state round-trips an
   added optional block belongs there.)

2. C-consumer round trip (the load-bearing gate, mirrors the real ABI).
   Extend inst/tinytest/capi/consumer.c + test-capi.R:
   - a legacy-sized call: stack-allocate `dbarts_results r = {0};`, set
     `r.structSize` to the pre-append boundary (the offset of
     `logLikelihood` - what an all-v1 caller would have declared), fill
     sigma/train/varcount, run, and assert against a full-sized call
     (`structSize = sizeof r`, same fields plus logLikelihood if shipped)
     that the shared fields are IDENTICAL - proving the size guard is
     transparent to an older-sized caller.
   - The write-guard canary: allocate a struct region larger than
     `dbarts_results`, set `structSize` to a boundary below `test`/`varcount`
     (simulate an OLD, undersized caller), write a sentinel into the bytes
     past that offset, run, assert the sentinel is intact AND the guarded
     field was never dereferenced (leave it pointing at a poisoned address so
     a stray write would segfault). This exercises the ABI-safety claim from
     a genuinely-compiled consumer - the one test that would actually catch a
     guard regression.
   - `dbarts_apiVersion()` stays 1 (the struct is born size-first
     pre-release, so there is nothing to bump here - unlike the
     originally-proposed run2 duality, which would have moved this to 2).
   - Self-gates on toolchain availability exactly as today (test-capi.R:13-38).

3. State additive-load tinytest (new file, e.g.
   inst/tinytest/test-sampler-state-format.R): fit, `storeState`, then
   - hand-set `attr(state, "formatVersion")` to a FUTURE additive value (e.g.
     4) and confirm `setState`/predict still loads (floor >= 3 passes,
     by-name reader unaffected) - the anti-orphan proof;
   - hand-set it BELOW the floor (e.g. 2) and confirm refusal names both
     versions - the floor proof;
   - delete a REQUIRED per-chain slot (e.g. `sigma`) from a chain and confirm
     refusal names the block ("required block 'sigma'") - the naming proof;
   - delete an OPTIONAL slot on a fit that has it off-path and confirm it
     still loads - the default proof.
   States are opaque R lists with attributes, so this is pure R attribute
   surgery (no C needed); test-sampler-saveLoad.R is the sibling for style.

4. tinytest suite: `tinytest::test_package("dbarts")` - baseline + the new
   checks, 0 failures. No existing snapshot moves (no RNG path touched).

5. Equivalence: `benchmarks/R/equivalence.R compare
   benchmarks/baselines/equivalence-<current>.rds` - all scenarios IDENTICAL
   (same RNG stream). The growable `dbarts_sampler_run` adds no draw; a
   legacy-sized caller gets byte-identical output via the size guard; the
   state gate change only affects load-time refusal, not sampling.

6. air format + lintr on any touched R/tinytest files; dbarts.h ASCII-clean
   (grep for non-ASCII); the header still compiles as C99 (consumer.c is C -
   the CI compile in test-capi.R is the check; keep static_asserts OUT of the
   header).

### stan4bart lockstep

Adopted design: stan4bart's port adds `r.structSize = sizeof r` at its
`dbarts_results` allocations and gates on `dbarts_apiVersion() >= 1` as
before - the one lockstep rebuild for 1.0-0 that stan4bart already needs (it
is not yet released against the current header) absorbs this size-first
layout change for free. Record this in the release note (TODO L445-449
area).

Under the originally-proposed run2/results2 duality (Design, Part 1), this
would instead have been fully optional: `dbarts_sampler_run` and
`dbarts_results` would stay byte-frozen (via the shim), so stan4bart
0.0-13's compiled `dbarts_results results = {0}` + `dbarts_sampler_run`
calls would be unaffected, and a stan4bart that never touched `run2` would
need no rebuild at all. When stan4bart eventually wanted the loglik/NB
channels it would switch to `run2` with `structSize = sizeof(
dbarts_results2)` and gate on `dbarts_apiVersion() >= 2` - additive, at its
own pace, in a later lockstep window. That optionality was the tradeoff
given up when the duality was dropped: because stan4bart's 1.0-0 rebuild
was happening anyway, requiring it to also set `structSize` cost nothing
extra in practice.

The state gate relaxation (Part 2) is backward-compatible independent of
either design: states stan4bart writes today (format 3) still load (>=
floor), and this half needs no stan4bart change regardless of which
output-channel design shipped.

## Risks and sequencing

- Smallest-diff split (recommended): two independent landings under this one
  plan, since the two surfaces share nothing.
  - Landing A (output channels): dbarts.h struct/macro additions to
    `dbarts_results` (no version bump - the struct is born size-first
    pre-release); C_interface.cpp folds the write guard directly into
    `dbarts_sampler_run` plus the static_asserts; no new R_interface.cpp
    registration (`dbarts_sampler_run`'s existing CCallable entry carries the
    new signature); consumer.c + test-capi.R updated to size-set the struct
    at each call site. Self-contained.
  - Landing B (state blocks): R_interface_bartcore.cpp gate swap + block-name
    refusal helper + registry comment; docs note; state-format tinytest.
    Self-contained.
  They can land in either order; neither depends on the other. If VD prefers
  one commit, they still review as two logical halves.
- Risk: a consumer that forgets to set `structSize` (leaves it 0) gets an
  all-skip no-op run - safe but silently empty. Mitigation: loud header doc
  (the `= {0}; r.structSize = sizeof r;` line) and the consumer test models
  the correct idiom. Considered and rejected: a "0 means full" sentinel -
  it re-introduces exactly the size-blind write the design exists to prevent.
- Risk: appending a field later without bumping the `sizeof` static_assert or
  `DBARTS_C_API_VERSION`. Mitigation: the `sizeof(dbarts_results) ==
  EXPECTED` assert fails the build until EXPECTED is updated, forcing the
  author to notice the ABI grew (and, by convention, bump the version).
- Risk: someone bumps `stateFormatVersion` for a purely additive block out of
  habit, re-orphaning states. Mitigation: the registry comment states
  additive changes DO NOT bump; the floor stays 3; a review checklist line.
  The relaxed reader tolerates a spurious bump anyway (floor is >=, not ==),
  so the failure mode is degraded provenance, not orphaned states - a strict
  improvement over today.
- Risk: the loglik channel ships but no consumer fills-and-reads it yet
  (dead surface). Mitigation: this is the mechanism-only-vs-loglik-now
  judgment call above; if mechanism-only wins, Landing A ships size + the
  eight v1 fields and the canary alone proves the guard, with logLikelihood
  deferred to pointwise-loglik's engine-side follow-up.

## Landing log

- 2026-07-10: plan drafted (read-only). Not yet implemented. Open judgment
  calls for the implementer/VD, all named inline: (1) ship the loglik channel
  at run introduction vs mechanism-only; (2) monotonic-order static_asserts
  vs exact per-field offsets; (3) encoding floor vs pure-provenance (no gate).
  Recommendations recorded above: (1) loglik as first append in the same PR,
  (2) monotonic + one sizeof assert, (3) keep the floor.

- 2026-07-10: judgment calls adjudicated by the orchestrator (VD
  delegated sequencing; the design itself was decided 2026-07-10):
  (1) loglik ships as the first append in the same landing - it is
  the one channel with a decided item and a live producer, and it
  exercises the write guard end-to-end; if the engine-side per-
  family loglik producer balloons past budget, fall back to
  mechanism-only and report. (2) monotonic-order static_asserts
  plus the one exact sizeof assert. (3) keep the encoding floor.
  Two-landing split (A: output channels; B: state blocks) accepted;
  land as separate commits under this plan.

- 2026-07-10: LANDED as two commits on wt/c-api-growth.
  Commit A, "Add a growable results struct and run2 entry point" - its
  initial content shipped the then-current run2/results2 design (see
  Design, Part 1, "Originally proposed" for what this became a few hours
  later): dbarts.h (dbarts_results2 + DBARTS_RESULTS2_HAS + run2 prototype +
  DBARTS_C_API_VERSION 2); C_interface.cpp (runCore factoring, run
  as a v1-boundary shim, run2, and the monotonic-order + exact-sizeof
  static_asserts); R_interface.cpp registration; consumer.c capi_run2
  + capi_run2_guard canary; test-capi.R (apiVersion 1L->2L, run2==run
  bit-identity, guard).
  Commit B, "Read sampler state blocks by name behind an encoding floor" -
  landed as designed and unaffected by the later design change:
  minReadableStateFormatVersion floor replacing the two equality gates
  in setState/installForests; missing-vs-malformed block-name refusals
  on forests/tree.vars/tree.params/tree.masks/k/sigma/fit.scale; the
  registry-rule comment at storeState; the public-surface.md addition;
  new test-sampler-state-format.R; and test-bartcore.R's future-version
  refusal test retargeted to a below-floor refusal (the equality gate it
  asserted no longer exists).
  DEVIATION (judgment call 1, adjudicated fallback taken): shipped
  MECHANISM-ONLY, NOT loglik-now. The engine-side per-family loglik
  producer balloons past budget - it needs a new virtual across all four
  response classes (Gaussian/Probit/Logistic + the GroupedResponse
  decorator) with distinct family-specific math that must exactly match
  the just-landed R-side pointwiseLogLikelihood, plus per-observation
  group-index plumbing into storeSample for rbart (trainingFits omit the
  intercepts) and BCF NaN handling, with no gate verifying the numeric
  values. Per the adjudication, the results struct ships size + the eight
  v1 fields; logLikelihood is deferred to pointwise-loglik's engine-side
  follow-up (it landed the same day - see below). The write guard is still
  proven end-to-end: capi_run2_guard pins structSize below `test`/`varcount`
  (fields the gaussian sampler DOES produce) and poisons those slots, so a
  size-blind write would crash - a stronger canary than a null field.
  Judgment calls 2 (monotonic + one sizeof) and 3 (keep the floor) landed as
  recommended.
  Gates: preclean install clean; tests/cpp all green; test-capi 58/58
  (run2==run sigma/train/varcount, guard TRUE); tinytest 2602/0;
  equivalence 21/21 identical draws; air + lintr clean; dbarts.h ASCII
  and C99-clean (consumer compiles+runs as C).

- 2026-07-10: DESIGN CHANGE (VD), superseding Commit A's run2/results2
  duality (full adopted design and rationale in Design, Part 1). Nothing
  on CRAN links against the unreleased 1.0-0 header, so dbarts_results is
  still editable: it becomes size-first ITSELF (a leading size_t
  structSize, shifting the eight pointers), and the one lockstep stan4bart
  rebuild already planned for 1.0-0 absorbs the layout change. dbarts_results2,
  dbarts_sampler_run2, and the run shim are dropped; the macro is
  DBARTS_RESULTS_HAS; DBARTS_C_API_VERSION stays 1 (the struct is BORN
  size-first at 1.0-0 - there is no prior layout to distinguish; appends bump
  it from here). The write-guard idiom, the static_asserts (structSize at 0,
  monotonic order, one exact sizeof), the boundary comment, and the
  mechanism-only adjudication (above) all carry over verbatim into
  dbarts_sampler_run, which absorbs runCore's body directly. A zero
  structSize is an all-skip no-op (documented in the header).
  stan4bart lockstep action: its port adds r.structSize = sizeof r at its
  dbarts_results allocations and gates on dbarts_apiVersion() >= 1 as
  before.
  Commit B was approved unchanged and landed on bartcore as 2d09768
  (cherry-picked verbatim; gates re-run clean). Commit A reworked in place
  per the above; the guard canary retargets dbarts_sampler_run itself
  (structSize pinned to offsetof(test) with poisoned slots past it) and
  every consumer.c run call site sets structSize - the in-repo lockstep
  rebuild. Gates re-run after the rework; results recorded below.

- 2026-07-10: DEFERRED CHANNEL LANDED (the logLikelihood follow-up this
  plan named). dbarts_results gains a ninth pointer, logLikelihood
  (numObservations x numSamples x numChains), appended INSIDE the born v1
  layout - 1.0-0 is unreleased, so it ships below groupEffects with the
  "1.0-0 field boundary" comment moved beneath it and DBARTS_C_API_VERSION
  left at 1 (the same in-repo lockstep the size-first rework already took).
  C_interface.cpp maps it through the FILL guard and the exact-sizeof
  static_assert grew to 9 pointers (its message reworded: appends AFTER
  1.0-0 bump the version; this pre-release append does not). Engine:
  bartcore::Results gains logLikelihood (default null = skip, zero work);
  storeSample fills it per non-BCF observation via a new
  ResponseModel::computeLogLikelihood virtual. Family coverage: gaussian
  (dnorm(y, f(x)+offset, sigma/sqrt(w)) via Rf_dnorm4), probit
  (log dbinom(y, 1, Phi(eta)) via the stable Rf_pnorm5 log-tail), and
  logistic (w * log dbinom(y, 1, plogis(eta)) via a stable log1pExp, integer
  count weights). GROUPED samplers are FILLED, not NaN'd: the GroupedResponse
  decorator adds the per-observation intercept (the shiftFits pattern, ~6
  lines) then defers to the base family, matching rbart_vi's extract "ev" +
  ranef location - the <= 30-line clean threading the scope allowed. NaN
  POLICY: BCF NaN-fills the channel (the two-forest blended location is
  invisible to the response model, exactly the testFits precedent); the base
  ResponseModel default NaN-fills so any future family reports "unavailable"
  rather than a wrong value. The R bridge never sets Results::logLikelihood
  (default-initialized null), so R-side memory and every draw are unchanged;
  equivalence stays byte-identical. Verification: test-capi.R asserts the C
  channel == dnorm(y, train, sigma, log) recomputed in R (gaussian, exact to
  1e-12, same Rf_dnorm4) and == the stable/dbinom probit forms; tests/cpp
  testLogLikelihood asserts NULL costs nothing (bitwise sigma/train with and
  without the channel, same seed) and the per-family math on a fixed input;
  the write-guard canary still pins structSize below the new field and
  passes. Gates recorded below.
