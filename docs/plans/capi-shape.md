# dbarts.h shape freeze: the last breaking-if-deferred items

Status: PROPOSED. Read at bartcore 9d0ee10f; revised after an independent blind
critique (two blocking defects, several citation corrections, all folded). One
slice, one hash re-bake, one lockstep consumer rebuild - the shape of
docs/plans/dbarts-h-freeze.md (D3/D4), which this continues and does not reopen.

rng: NEUTRAL. Nothing here touches a proposal, a suffstat, a cut builder or an
RNG consumer. Ten of the 48 entry points change shape; every body below them is
the one that runs today. The equivalence trio is a FORMALITY and must be
labelled one in the landing note: `benchmarks/` has zero hits for
`dbarts_sampler_`, `DBARTS_USE_STUBS` or `dbarts/dbarts.h` (measured), so no
harness drives a flat entry point. The load-bearing gate is
`inst/tinytest/test-capi.R` + `inst/tinytest/capi/consumer.c`.

## 0. What a return value means, per entry

The header must not claim a universal reading of `int`. A full sweep of the 48
entries puts every non-void return in exactly one of three classes, and the
classes do not share a rule:

- **VALUE** (12): `apiMajorVersion`, `apiMinorVersion`, `apiHash`, the six
  `size_t` counts, `kIsSampled`, `usesDart`, `dbarts_sampler_family`. The number
  IS the answer; there is no refusal in it. `kIsSampled`/`usesDart`/`family`
  matter here because they are `int`-returning and NOT statuses.
- **TRANSACTION** (2): `setPredictor`, `updatePredictor`. 0 means the mutation
  was attempted and ROLLED BACK - every tree could not keep non-empty leaves -
  and a different argument would have worked (src/C_interface.cpp:753, :786;
  CAPI:869-870). This is the class the earlier draft's universal sentence would
  have contradicted.
- **CAPABILITY STATUS** (9 today, 16 after this slice): `getLatents`,
  `getDispersion`, `setForestBasis`, `getForestFits`, `getForestAmplitudes`,
  `setForestWeights`, `getForestCalibration`, `setForestPriorScale`,
  `setActiveRows`, plus the seven of item 1. Here, and only here:

> 0 means the SAMPLER cannot do this at all - no argument would have worked -
> and nothing was touched. An `Rf_error` means THIS call is wrong and the caller
> can fix it and call again.

Objects and voids carry no status: `dbarts_sampler_create` (pointer),
`getTrees`/`storeState` (SEXP), `printTrees` and the remaining voids refuse only
by raising.

Checked against every capability-status entry as it stands: `getLatents` 0
(family augments nothing), `getDispersion` 0 (no dispersion), `setForestBasis` 0
(no amplitudes / forest names no forest) vs raise (null basis, zero width,
non-finite value), `setForestWeights` 0 vs raise (non-finite, negative),
`setForestPriorScale` 0 vs raise (non-finite, non-positive),
`getForestCalibration` 0 vs raise (zero structSize), `getForestFits` 0 (no such
forest). All conform. `setActiveRows` is the single violation inside the
capability class - item 3.

The header states this as three named classes plus the capability rule, each
entry's own doc saying which class its return is in. It never says "0 means X"
of `int` as such.

## 1. Capability-refusal split (void -> int)

The five `void` setters refuse on a MIX of both classes today, so "unify" has to
say which refusal moves to which channel, not just which type the entry returns.
Taken from src/C_interface.cpp and the shared guards at
src/R_interface_bartcore.cpp:2613-2925:

| entry | capability refusals (-> 0) | call-is-wrong refusals (stay raise) |
| --- | --- | --- |
| setSigma | `refusePinnedSigmaChange`: variance forest owns the scale; family pins sigma. LIVE | none today |
| setWeights | `refuseBinaryWeightChange` (probit/ordinal/aft/nbinom), LIVE; multi-forest conduit fixed at creation, flat-unreachable | logistic counts not positive integers |
| setResponse | multi-forest conduit fixed at creation. FLAT-UNREACHABLE (see below) | y outside the family's support; `updateScale != 0` against a pinned transform (multi-forest, heteroscedastic, grouped) |
| setOffset | multi-forest conduit fixed at creation. FLAT-UNREACHABLE | `updateScale != 0` against a pinned transform |
| setTestPredictors | `refuseUndefinedTestFits`. LIVE (BCF) | every source refusal (`translateSource`, `validateTestSource`), the offset-length pairing, the sparse leaf covariate |
| setTestOffset | `refuseMultiForestMutation`. LIVE (BCF) | none today |
| predict | `refuseUndefinedTestFits`. LIVE (BCF) | `refuseEmptyTreeStore`, every source refusal |

FLAT-UNREACHABLE means exactly this: the capability arm of
`refuseMultiForestResponseMutation` needs `numForests >= 2 &&
!supportsResponseMutation` (src/R_interface_bartcore.cpp:2641-2642), and the only
combiner answering false is `MultinomialForestCombiner` (default
src/bartcore/combiner.hpp:693, overridden true for `AmplitudeForestCombiner` at
:1059), which has no flat creation path (src/C_interface.cpp:304-330). The BCF
legs confirm it empirically: `LEG_RESPONSE_PINNED`, `LEG_OFFSET_PINNED` and
`LEG_WEIGHTS` all ACCEPT (inst/tinytest/capi/consumer.c:1231, :1237, :1243). So
`setResponse` and `setOffset` gain a 0 arm no sampler this header can build will
fire, and `setWeights`' multi-forest arm is likewise dead while its family arm is
live. Kept anyway, on the precedent `DBARTS_FAMILY_MULTINOMIAL` set at the
freeze: the contract is stated before the path that reaches it, so multinomial
creation opening later is a minor bump and not a shape change. Recorded in
section 13.

The `updateScale` arms stay raises deliberately, and the consumer sweep shows
why: stan4bart drives `dbarts_sampler_setOffset(..., isWarmup && iter %
update_scale_mod == 0)` once per iteration (src/init.cpp:647) and discards the
return. Its own TODO `bart-args-forests-guard` records that a two-forest sampler
"dies at the first setOffset with a loud dbarts-side refusal" - that loudness is
the guard. `AmplitudeForestCombiner::supportsResponseMutation()` is true
(src/bartcore/combiner.hpp:1059), so a BCF sampler reaches only the `updateScale`
arm, and stan4bart's first call passes `updateScale = 1`. The guard survives this
slice unchanged. Verified, not assumed.

Alternatives.

- **A. Unify downward** (the capability ints raise instead). Kills probing
  outright: a host asking "does this sampler carry latents?" would have to
  longjmp to find out. Costs ~9 bridge sites + every consumer gaining an
  R_tryCatch. Reject.
- **B. Unify upward, five entries** (the endorsed remedy). Costs: 5 list entries,
  5 prototypes, 5 bodies, ~12 consumer.c lines, one hash re-bake.
- **C. Unify upward, seven entries** - B plus `setTestOffset`
  (`refuseMultiForestMutation`, a fixed property) and `predict`
  (`refuseUndefinedTestFits`, the same predicate `setTestPredictors` uses).
  RECOMMENDED. B leaves those two as the new outliers on the exact axis being
  fixed: `setTestPredictors` would answer 0 where `predict` raises for the same
  reason on the same sampler. Marginal cost over B is 2 entries, and both new 0
  arms are LIVE.
- **D. Additive `dbarts_sampler_supports*` predicates, signatures unchanged.**
  Non-breaking, so by the freeze rule it does not have to land now - which is
  also the argument against it as the ANSWER: it leaves the split frozen forever
  and states each rule twice (predicate + guard), the duplication the shared
  helpers exist to prevent. Reject as a substitute; it remains available as a
  post-1.0 convenience.

RECOMMENDATION: C, including the two entries whose 0 arm is unreachable today.
The contract is right where the code cannot yet exercise it, and a shape fixed
now cannot be fixed later.

Cost that must be stated, not hidden: a discarded 0 leaves the sampler
unconditioned and the run continues on whatever it held - a quieter failure than
today's longjmp. Three things carry it: the header says so in the doctrine
paragraph; the consumer migration below checks the answer ONCE at setup (a
capability answer is invariant over a sampler's life, so no per-sweep check is
needed); and `DBARTS_NODISCARD` on the non-void stubs is recorded as a door, not
taken.

## 2. `forest` argument position

Eight entries take it as argument 2 (`numTrees`, `setForestBasis`,
`getForestFits`, `numForestAmplitudes`, `getForestAmplitudes`,
`setForestWeights`, `getForestCalibration`, `setForestPriorScale`); `getTrees`
and `printTrees` take it last, after `useLiveTrees`, where the multi-forest work
appended it.

- **A. State the shape rule in a comment and leave the order.** Free, and freezes
  a rationalization: there is no rule that puts it last on exactly those two.
- **B. Move it to argument 2 on both.** RECOMMENDED. It is not only symmetry: the
  forest qualifies arguments that follow it, since `treeIndices` are read against
  THAT forest's tree count - `printTrees` checks it inline
  (src/C_interface.cpp:925-928) and `getTrees` delegates
  (src/C_interface.cpp:885-889) to the check at
  src/R_interface_bartcore.cpp:7734-7735. The R5 twins (`$getTrees`,
  `$printTrees`, R/dbarts.R:1968, :1999) carry no forest selector at all, so
  there is no R argument order to mirror and nothing pulls the other way.

Compile-visibility caveat: position 2 goes from `const size_t*` to `size_t`,
which no C++ call site converts, so both stan4bart sites fail to compile until
edited. dbarts.h is a C header, though, and a C consumer's pointer-to-integer
mismatch is a constraint violation some toolchains still diagnose as a warning
rather than an error. The REAL backstop for a C consumer is the hash re-bake:
`DBARTS_REQUIRE_EXACT_ABI` raises at first resolution, and a consumer without it
still fails the major/minor handshake path unchanged, so a mis-ordered binary
cannot silently run. Do B, and state the rule in the header so the next entry
that names a forest lands in the right place.

## 3. `setActiveRows`' two zeros

`Chain::setActiveRows` (src/bartcore/chain.hpp:1666) returns false both when the
family implements no mask and when an element is not exactly 0 or 1; the flat
entry maps both to 0. The R bridge already splits them: a capability raise first
(src/R_interface_bartcore.cpp:4038), then a post-hoc test of the engine's bool
with its own message (:4049-4052). By section 0's capability rule the value
refusal is recoverable and must raise.

Align. Write a shared `refuseNonBinaryMask` helper carrying the exact-{0,1} scan
and the R bridge's message text - NEW code, not a promotion, since the R bridge
has no scan of its own - and call it from both surfaces, leaving the engine's
bool contract alone. The R bridge's post-hoc test then becomes defense in depth.

## 4. `setForestBasis` and the column-major rule

Verified: the engine really is row-major (`SamplerBase::setForestBasis`,
src/bartcore/sampler.hpp:1528-38) and the R bridge transposes an R matrix on the
way in (src/R_interface_bartcore.cpp:3954-3968); the flat entry documents and
reads row-major (src/C_interface.cpp:1036-1042). The header is correct and the
global rule at CAPI:81 has an unnamed exception.

- **A. Transpose the flat contract to column-major.** Restores one rule, and the
  copy would absorb the transpose. Costs: a SILENT semantic break (same
  `const double*`, no diagnostic) for any caller; a per-call transpose on an
  entry a Gibbs host swaps every sweep; and it makes the flat surface state a
  layout the engine does not use, so both surfaces would transpose and neither
  would show the engine's own contract. No consumer calls the entry (verified:
  zero hits in stan4bart and in treatSens' dbarts-1.0 worktree), so the silent
  break is theoretical - but the other two costs are not.
- **B. Rename the parameter `basis` -> `basisRowMajor`, and name the exception in
  the global rule.** RECOMMENDED. The signature then documents itself at every
  call site, and the parameter name is inside the stringized signature the token
  folds, so the rename moves the hash and every `DBARTS_REQUIRE_EXACT_ABI`
  consumer gets a hard stop that sends it back to the entry's doc.

Consumer impact of B: none - no consumer calls it. In-repo: the live C test legs
already lay their bases row-major and correctly
(inst/tinytest/capi/consumer.c:1265-70, :1310-14, :1325-30); the R-facing wrapper
`capi_set_forest_basis` (:948-957) is the file's only orphan (no `.R` caller;
test-capi.R:57 resolves by name), its comment says "column-major", and it hands
`REAL(basisExpr)` through untransposed. RULING: FIX it, do not delete it -
consumer.c is a complete-coverage exemplar at 48/48 entries, `setForestBasis` is
otherwise covered only through the BCF legs, and fixing restores an R-facing leg.
Correct the comment to row-major, transpose properly, and rename its local to
match `basisRowMajor`.

## 5. `columnTypes` / `categoryCounts` above the boundary

VD ruling: KEEP, rationale in the header. Verified inert on every path a flat
entry takes: `translateSource` publishes `categoryCounts` into the engine view
(src/C_interface.cpp:217), but the only reader is
`ColumnSource::declaredCategoryCount` (src/bartcore/data.hpp:1186), set in the
TRAINING build, which no flat entry reaches; `buildTest` does not read it.
`columnTypes` is never published at all - the store's own types are gathered
instead. Both are checked for well-formedness (:144-152). The struct comment at
CAPI:290-297 stands; add, after it:

```
/// Both stay ABOVE the boundary rather than being dropped. They are the only
/// fields a source needs to describe a typing the SAMPLER does not already fix,
/// which is the one thing a predictor view cannot say about itself once they
/// are gone, and this is the last release in which they can be placed rather
/// than appended. They cost a filling caller nothing - null is "nothing
/// declared" - and what a caller does fill is checked, so a malformed
/// declaration is refused here rather than the first time an entry reads one.
```

and change the two field comments from "declared, not read" to "declared and
checked; no entry reads it" (same length, no offset moves).

## 6. `USE_FC_LEN_T`: remove it, do not balance it

The earlier draft proposed a symmetric push/pop. That remedy was built on a
mechanism that does not exist, and the real leak runs the other way.

Corrected facts. R's headers never read `USE_FC_LEN_T`; `R_ext/BLAS.h:38` gates
on `FC_LEN_T`, which `Rconfig.h:20` defines (R >= 4.3) unless
`DONT_USE_FC_LEN_T`. Measured on R 4.6.1: preprocessing `#include
<dbarts/dbarts.h>` then `<R_ext/BLAS.h>` yields `FC_LEN_T` DEFINED despite
CAPI:107, so the `#undef` strips nothing. In the only window where the macro
mattered (R 3.6.2-4.2), `Rconfig.h` is reached through CAPI:102's own
`<Rinternals.h>` and is include-guarded, so `FC_LEN_T` is frozen before :107
runs. The `#undef` closes nothing in either window.

What the DEFINE does is real and unwanted: when this header's `<Rinternals.h>` is
the translation unit's first, CAPI:96 turns `FC_LEN_T` ON for the consumer's
whole TU through `Rconfig.h`, and no restore of `USE_FC_LEN_T` undoes it -
`FC_LEN_T` is what `BLAS.h` reads and is not undoable. dbarts.h declares no
Fortran routine, so it has no business setting it.

REMEDY: delete CAPI:95-97 (the `<Rversion.h>` include and the `R_VERSION`-gated
`#define`) and CAPI:107 (the `#undef`) outright. The `R_NO_REMAP` dance at
:98-106 and its comment stay untouched. Hash-invisible.

In-repo verification that nothing relies on dbarts.h supplying the define: only
two TUs under `src/` include `<dbarts/dbarts.h>` (src/R_interface.cpp:17,
src/C_interface.cpp:7), and `src/` contains ZERO occurrences of `R_ext/BLAS.h`,
`R_ext/Lapack.h`, `FCONE`, `F77_CALL`, `dgemm`, `dpotrf` or `LAPACK` - no file
under `src/` declares or calls a Fortran routine at all. Both TUs additionally
obtain the macro from their own `external/*.h` headers, which carry the same
define/undef dance (src/include/external/R.h:9/:36,
src/include/external/Rinternals.h:9/:41). A future dbarts TU that does call
Fortran defines the macro itself, as every existing `external/*.h` already does.
Same for consumers: treatSens already defines its own
(src/R_interface.cpp:16-23), stan4bart never needs it. One migration note to
carry: a consumer that was relying on dbarts.h to pull `<Rversion.h>` in for it
now includes it itself.

## 7. Signature diff (inst/include/dbarts/dbarts.h)

Entry count stays 48. `DBARTS_C_API_MAJOR`/`MINOR` DO NOT MOVE (CAPI:114-118: no
version constant increments before 1.0-0). Ten `DBARTS_C_API_LIST` entries:

```
  X(int, dbarts_sampler_setResponse,
    (dbarts_sampler* sampler, const double* y, int updateScale),
    (sampler, y, updateScale))
  X(int, dbarts_sampler_setOffset,
    (dbarts_sampler* sampler, const double* offset, int updateScale),
    (sampler, offset, updateScale))
  X(int, dbarts_sampler_setWeights,
    (dbarts_sampler* sampler, const double* weights), (sampler, weights))
  X(int, dbarts_sampler_setSigma,
    (dbarts_sampler* sampler, double sigma), (sampler, sigma))
  X(int, dbarts_sampler_setTestPredictors,
    (dbarts_sampler* sampler, const dbarts_predictor_source* xTest),
    (sampler, xTest))
  X(int, dbarts_sampler_setTestOffset,
    (dbarts_sampler* sampler, const double* offsetTest),
    (sampler, offsetTest))
  X(int, dbarts_sampler_predict,
    (dbarts_sampler* sampler, const dbarts_predictor_source* xTest,
     const double* offsetTest, size_t numThreads, double* out),
    (sampler, xTest, offsetTest, numThreads, out))
  X(SEXP, dbarts_sampler_getTrees,
    (dbarts_sampler* sampler, size_t forest, const size_t* chainIndices,
     size_t numChainIndices, const size_t* sampleIndices,
     size_t numSampleIndices, const size_t* treeIndices,
     size_t numTreeIndices, int useLiveTrees),
    (sampler, forest, chainIndices, numChainIndices, sampleIndices,
     numSampleIndices, treeIndices, numTreeIndices, useLiveTrees))
  X(void, dbarts_sampler_printTrees,
    (dbarts_sampler* sampler, size_t forest, const size_t* chainIndices,
     size_t numChainIndices, const size_t* sampleIndices,
     size_t numSampleIndices, const size_t* treeIndices,
     size_t numTreeIndices, int useLiveTrees),
    (sampler, forest, chainIndices, numChainIndices, sampleIndices,
     numSampleIndices, treeIndices, numTreeIndices, useLiveTrees))
  X(int, dbarts_sampler_setForestBasis,
    (dbarts_sampler* sampler, size_t forest, const double* basisRowMajor,
     size_t numColumns),
    (sampler, forest, basisRowMajor, numColumns))
```

The readable prototypes at :807, :815, :826, :839, :905, :910, :933, :964, :974,
:1075 follow, each with its own doc gaining one sentence naming its return class
and, for the capability seven, what its 0 means. `printTrees` stays `void`: it
refuses only on out-of-range indices and the empty store, both recoverable.

Contracts block (CAPI:44-82), three edits:

- Replace the last sentence of the first bullet ("Where an entry does refuse...")
  with the three-class paragraph of section 0: a non-void return is a VALUE
  (`kIsSampled`, `usesDart`, `dbarts_sampler_family` and the counts - the number
  IS the answer, and an `int` here is not a status), a TRANSACTION result
  (`setPredictor`, `updatePredictor` - 0 is a completed rollback, and a different
  argument would have worked), or a CAPABILITY STATUS, whose rule is: 0 means the
  SAMPLER cannot do this at all and nothing was touched, so a host driving a
  sampler it did not build probes the channel instead of unwinding through its
  own frames; an `Rf_error` means this call is wrong - a value outside a family's
  support, a source whose declared shape disagrees, a malformed struct, a scale
  update against a transform the sampler pins - and longjmps through the caller's
  frames, so call from contexts that are safe to unwind. Each entry's own doc
  says which class its return is in. A DISCARDED capability 0 leaves the sampler
  unchanged and the run conditioned on what it held before; that answer is a
  fixed property of the sampler, so test it once at setup rather than every
  sweep.
- New bullet: where an entry names a forest, `forest` is the argument after the
  sampler, a 0-based index in [0, `dbarts_sampler_numForests`), and it qualifies
  every argument after it - a tree index list is read against THAT forest's tree
  count. An index past the last forest is a 0 only on the capability-status
  entries. It RAISES everywhere else: on the `size_t` probes (`numTrees`,
  `numForestAmplitudes`), whose value carries no refusal a caller could tell from
  a legitimate answer, and on `dbarts_sampler_getTrees` and
  `dbarts_sampler_printTrees`, where a caller assembling a data.frame could not
  tell an `R_NilValue` refusal from an empty answer.
- The column-major bullet names its one exception: `dbarts_sampler_setForestBasis`
  takes `basisRowMajor`, row i at `basisRowMajor + i * numColumns`, that
  contraction being the engine's only read of a basis.

Include block: delete CAPI:95-97 and CAPI:107 (section 6). Nothing replaces them.

## 8. Bridge changes (src/C_interface.cpp, and the shared guards)

The refusal texts live in `bartcore_bridge` so the two surfaces cannot state
different rules; a flat entry must therefore PREDICT a refusal, never restate it.
Declarations in src/R_interface_bartcore_common.hpp, bodies at
src/R_interface_bartcore.cpp:2613-2925.

FOUR whole-function guards take the mechanical form - a predicate stating the
condition, and the raiser rewritten as `if (!predicate(...)) return;` above its
existing message:

- `bool sigmaIsPinned(const SamplerBase&)` - `refusePinnedSigmaChange`'s
  disjunction (variance forest, or family not gaussian/aft).
- `bool familyCarriesNoWeights(const SamplerBase&)` - `refuseBinaryWeightChange`.
- `bool testFitsAreUndefined(const SamplerBase&)` - `refuseUndefinedTestFits`.
- `bool isMultiForest(const SamplerBase&)` - `refuseMultiForestMutation`, for
  `setTestOffset`.

`refuseMultiForestResponseMutation` does NOT take that form: it carries two arms,
and predicating the whole function would short-circuit exactly the BCF samplers
whose `updateScale` arm must fire, silently deleting the decalibration guard
stan4bart depends on. Only the capability block at
src/R_interface_bartcore.cpp:2642-2675 may be predicated. Exact restructure:

```
bool responseConduitIsFixed(const bartcore::SamplerShape& shape) {
  return shape.numForests >= 2 && !shape.supportsResponseMutation;
}

void refuseMultiForestResponseMutation(...) {
  bartcore::SamplerShape shape = sampler.shape();
  if (responseConduitIsFixed(shape)) {
    ... :2643-2674 unchanged ...
  }
  if (shape.numForests >= 2 && conduit != ResponseConduit::weights &&
      updateScale != FALSE)
    ... :2677-2680 unchanged ...
}
```

The `numForests >= 2` test moves from the deleted early return at :2641 onto the
`updateScale` arm, which is where it was doing its work; the R bridge's behaviour
is bit-for-bit what it is today.

Entry bodies: each of the seven gains its predicate test returning 0 ahead of the
guards it keeps, and `return 1;` at the end - including
`dbarts_sampler_setTestPredictors`' EARLY return on the null-`xTest` removal path
(src/C_interface.cpp:803), which becomes `return 1;`, and `predict`'s return
after the offset add. `setForestBasis` renames its parameter;
`getTrees`/`printTrees` reorder theirs, bodies otherwise unchanged.
`setActiveRows` calls the new `refuseNonBinaryMask` after the
`supportsActiveRows` probe, keeping `? 1 : 0` as defense in depth. Cost:
`setActiveRows` now scans the mask twice, entry and engine - one extra O(n) pass
per install.

Note for `setWeights`: `enforceBinaryWeightPolicy`'s probit arm becomes
unreachable on this path (probit already answered 0), which is correct - it stays
for creation.

## 9. Hash re-bake

Items 1, 2 and 4 all move the stringized signature list; items 3, 5 and 6 are
invisible to it. Two literals move together in the code commit:
`DBARTS_C_API_HASH` (CAPI:166, currently `0x0939c0224353505bULL`) and
`dbarts_apiSignatureToken` (src/C_interface.cpp:523, currently
`0x32e5b15aa6c88c69ULL`). No struct layout, enumerator or callback parameter
moves, so the layout fold is unchanged (src/C_interface.cpp:527 uses the macro,
not a literal). Derive both with the temporary `DbartsShowToken` probe procedure
recorded at dbarts-h-freeze.md section 5, including its `printf '0x%016xULL\n'`
gotcha.

Live-literal pin sites, complete: CAPI:166; src/C_interface.cpp:523;
inst/tinytest/test-capi.R:87; docs/design/threaded-predict.md:111 (signature
token) and :112 (ABI hash). Add the outgoing `"0x0939c0224353505b"` to
test-capi.R's stale-token block at :73-84 - FIVE `expect_false` lines, :73, :74,
:78, :81, :84 - and set :87 to the new literal.
`expect_equal(versions, c(1L, 0L))` is unchanged.

docs/plans/dbarts-h-freeze.md:188 records the transition `0x66d33f1613892406 ->
0x0939c0224353505b` and STAYS AS WRITTEN: it is that slice's history, true of the
commit it describes, not a claim about the live header. Decided explicitly, not
by omission.

## 10. Tests

`inst/tinytest/capi/consumer.c`:

- Return the status instead of `R_NilValue`: `capi_set_response` (:509),
  `capi_set_weights` (:515), `capi_set_test_offset` (:520), `capi_set_offset`
  (:560), `capi_set_sigma` (:567).
- `setTestPredictors` has THREE call sites in two wrappers - :764 and :767
  (null-clear and install) and :794 - all of which must propagate the status.
- `predict` has FOUR call sites - :807, :836, :860-862, :879-880 - each of which
  returns `R_NilValue` on a 0, the shape `capi_get_latents` (:571-579) already
  has.
- `sweepCallback` calls `setSigma` at :355 INSIDE a `dbarts_sampler_callback`:
  the status must not be discarded there, and it must not `Rf_error` out of a
  callback either - test it and return 0 from the callback to stop the run.
- `capi_print_trees` (:529) and `capi_get_trees` (:888, call at :917-919) pass
  `forest` second.
- `capi_set_forest_basis` (:948-957): correct the comment to row-major, transpose
  the R matrix, rename the local to match `basisRowMajor`, and give it a
  test-capi.R caller so the restored leg is exercised.
- BCF legs (:1222-): `LEG_TEST_OFFSET`, `LEG_TEST_PREDICTORS` and `LEG_PREDICT`
  flip from raising to returning 0 - each sets `legs->accepted = (entry(...) ==
  0)` and its `legRefusals[]` entry becomes NULL. `LEG_RESPONSE_RESCALED` and
  `LEG_OFFSET_RESCALED` keep their refusal strings: that pair is the
  discriminating half of the split.

`inst/tinytest/test-capi.R`, expectation flips (each currently `expect_error`):
:467 setSigma on probit -> `0L`; :480 setWeights on probit -> `0L`; :517
setWeights on nbinom -> `0L`. New: setSigma on the heteroscedastic sampler built
at :222-228 -> `0L`. Unchanged, and named in the file as the half that proves the
gate discriminates: :512 setResponse over the count bound still raises; :1703 and
:1707 logistic fractional/zero counts still raise; the `updateScale` refusals
still raise. :1409 `capi_set_active_rows(ptrMaskA, rep(0.5, n))` changes from
`expect_equal(..., 0L)` to `expect_error(..., "exactly 0 or 1")`; :1407 and :1410
are unchanged. Plus the hash pins of section 9.

`tests/cpp`: no facade virtual changes shape, so `test_facade.cpp` needs no edit
and the "always `--preclean` after changing virtuals in facade.hpp" bus-error
rule is NOT triggered (`--preclean` is still mandatory - a shipped header
changes). `test_facade.cpp:1035-1053` (`SPY_RET(bool, setActiveRows, ...)`)
already pins the facade's non-binary refusal and stays as written. One case is
worth adding, and it is the invariant the item-3 split rests on:
`test_sampler.cpp` gains a case asserting that a chain reporting
`supportsActiveRows()` still returns false from `setActiveRows` on a non-binary
mask, i.e. the engine's two false reasons stay separable by the probe.

## 11. Consumer migration (lockstep, after dbarts installs clean)

Both flat-API consumers carry `DBARTS_REQUIRE_EXACT_ABI`, so the re-bake alone
forces a rebuild of each; without the edits below neither would load. Neither has
a function-pointer table, a hand-rolled prototype, or a return value in a
type-sensitive context, so the void -> int change is source-compatible at every
call site.

**stan4bart**, /Users/vdorie/Repositories/stan4bart branch `bartcore`. No
vendored header (LinkingTo only); `DBARTS_USE_STUBS -DDBARTS_REQUIRE_EXACT_ABI`
in src/Makevars.in:3 and src/Makevars.win:3; includes at src/init.cpp:24 and
src/bart_util.hpp:10. Two MANDATORY edits:

- src/init.cpp:439 `dbarts_sampler_printTrees(fit, 0, chainIndices,
  numChainIndices, sampleIndices, numSampleIndices, treeIndices, numTreeIndices,
  0);`
- src/init.cpp:502-504 `dbarts_sampler_getTrees(fit, 0, chainIndices,
  numChainIndices, sampleIndices, numSampleIndices, treeIndices, numTreeIndices,
  useLiveTrees ? 1 : 0)`

Recommended, one check each at SETUP (the answer cannot change over the sampler's
life, so the per-sweep call sites stay bare): src/init.cpp:240 `setOffset` covers
:647 (:646 computes `update_scale_mod`; :648 is a commented-out duplicate); :242
`setSigma` covers :629; :268 `getLatents` covers :676; :349 `predict` now returns
int. No `setWeights`, `setResponse`, `setTestPredictors`, `setTestOffset`,
`setForestBasis` or `setActiveRows` call site exists. Its hand-rolled
major/minor check (:969-978) is unaffected. Its CI reinstall step
(feb8c29/f9bca65) is what keeps the runner from using a cached pre-re-bake
dbarts; keep it.

**treatSens**, the treatSens dbarts-1.0 branch
branch `dbarts-1.0` at 1db3d89. `DBARTS_REQUIRE_EXACT_ABI` + `DBARTS_USE_STUBS`
per translation unit (src/R_interface.cpp:29-30, src/bartTreatmentModel.cpp:10-11,
src/sensitivityAnalysis.cpp:32-33). NO mandatory edit - it calls no reordered
entry - so a rebuild alone compiles. Recommended setup checks:
src/sensitivityAnalysis.cpp:278 `setSigma`, :279 `setResponse` (covering :174,
:185, :292), src/bartTreatmentModel.cpp:85 `setOffset`. The
sensitivityAnalysis.cpp sampler is gaussian; the bartTreatmentModel.cpp one is
DBARTS_FAMILY_PROBIT (created at :63), and its `setOffset` still answers 1 -
`setOffset` has no family-keyed capability arm. Carry the hazard explicitly: a
probit sampler answers 0 to `setSigma`, so if that file ever gains one, an
unchecked call is a silent no-op. Its own `USE_FC_LEN_T` dance
(src/R_interface.cpp:16-23) is unaffected by item 6.

**R-API-only, confirmed by grep, not asserted**: bartCause (no `src/` at all;
DESCRIPTION has no `LinkingTo` field; zero `dbarts_` symbols repo-wide; uses
`bart2`/`xbart`/`extract`/`guessNumCores` via NAMESPACE plus `dbarts::bart`,
`dbartsData`, `dbartsControl`, `rbart_vi`, `forest`, `dbarts`) and bairrtt (has a
`src/`, but Rcpp/WALNUTS only - zero dbarts hits under `src/`; `LinkingTo: Rcpp,
RcppEigen`; uses `dbartsControl`, `dbarts`, `updatePredictorPerObservationJointly`,
`extract`).

CORRECTION to the task premise: treatSens is NOT R-API-only. Its MAIN worktree
still carries `LinkingTo: dbarts (>= 0.9-21)` against the deleted C++ header ABI
(`dbarts/bartFit.hpp` plus 14 `R_GetCCallable` lookups at
src/sensitivityAnalysis.cpp:734-747) - which is why the compat branch above
exists and why treatSens migrates in lockstep with stan4bart, not with bartCause.

## 12. Expected gate outcomes

`R CMD INSTALL . --preclean` MANDATORY (a shipped header changes; no facade
virtual moves, so the stale-object bus error is not in play). Deleting the
`benchmarks/kernels` binaries is NOT needed: dbarts.h is outside their include
graph. `tinytest::test_package("dbarts")` green. `cd tests/cpp && make &&
./test_bartcore` green. Equivalence trio BITWISE, 43/12/11 identical draws,
against equivalence-736bfb05.rds, bcf-equivalence-6e3b9fb8.rds and
multinomial-equivalence-4d9a3337.rds - a formality, since no harness drives a
flat entry point. `bench-sampler.R compare bench-sampler-ab1dc52.csv` on a quiet
machine: noise-level, no kernel or sweep code moves. `tools/check-rc-codoc.R`,
`air format --check .`, `R CMD check --as-cran` at the standing 0E/0W/1N. The
test-capi.R header compile matrix (default flags, stale hash alone, stale hash +
`DBARTS_REQUIRE_EXACT_ABI`, wrong MAJOR, wrong MINOR) all as before. Consumers:
`R CMD INSTALL` each, then `tinytest::test_package("stan4bart")` and treatSens
`testthat::test_local()`, both green before this slice lands. NEWS: one `\item`
in the 1.0-0 `\subsection{C API}` covering the seven return-type changes, the two
reordered entries, `basisRowMajor`, `setActiveRows`' value refusal becoming an
error, and the `USE_FC_LEN_T` removal.

Anchor budget, `tools/check-doc-freshness.R .`, re-aligned in the SAME commit
with each file's line count invariant. THREE files shift, not two:

- inst/include/dbarts/dbarts.h: 22 bare `dbarts.h:N` anchors plus 8 `CAPI:N`
  alias anchors (the alias is declared at tools/check-doc-freshness.R:147) = 30
  in docs/design, plus docs/plans/multiforest-extension-surface.md:3289.
- src/C_interface.cpp: 17 anchors.
- src/R_interface_bartcore.cpp: the predicates of section 8 insert into
  :2613-2925, shifting 90 `RIB:` anchors plus 26 bare
  `R_interface_bartcore.cpp:` anchors that sit at or below :2613. This is the
  largest of the three and the earlier draft omitted it.

Cross-slice coordination, since this slice is not file-disjoint from a
concurrently landing one: `inst/NEWS.Rd` takes an append under the same 1.0-0
`\subsection{C API}` (:1006) that any other pre-RC slice also appends to -
conflicts are expected at integration stacking and are resolved THERE, by the
practiced append-point resolution, not by pre-negotiating the text.
`src/R_interface_bartcore.cpp` is likewise shared: this slice touches the guard
block :2613-2925 and `bartcore_setActiveRows` :4030-4054, neither of which is
R-surface-disjoint. The `RIB:` anchor realignment therefore runs LAST, ONCE, over
the stacked tree - never concurrently with another slice re-anchoring the same
file, which would produce two line maps neither of which is true.

## 13. Residue (recorded, not fixed)

- `setResponse`' and `setOffset`'s new 0 arms are FLAT-UNREACHABLE: the only
  combiner with `supportsResponseMutation() == false` is the multinomial one, and
  multinomial has no flat creation path. So is `setWeights`' multi-forest arm
  (its family arm is live). Kept as stated contract, on the
  `DBARTS_FAMILY_MULTINOMIAL` precedent.
- `setActiveRows` can only ever return 1 after this slice. Every buildable family
  implements the mask (`ResponseModel::supportsActiveRows` default false at
  src/bartcore/model.hpp:2678, overridden true by all nine concrete families;
  GroupedResponse and AFTResponse forward), and every concrete `setActiveRows`
  returns true unconditionally, so src/bartcore/chain.hpp:1667 and :1677 are both
  dead and only the non-binary scan at :1672 can refuse - which this slice moves
  to a raise. The 0 arm stands as the declared capability channel.
- A heteroscedastic sampler answers `DBARTS_FAMILY_GAUSSIAN` while `setSigma`
  refuses it; the accessor still does not predict the refusal, but the probe now
  does without unwinding.
- `dbarts_sampler_kIsSampled`, `usesDart` and `dbarts_sampler_family` return
  `int` values that are NOT capability answers. The two-channel doctrine is a
  per-entry fact and never a property of the return type.
- `translateSource` publishes `categoryCounts` into the engine view although no
  flat path reads it (`declaredCategoryCount` is train-build only). If
  `buildTest` ever reads `categoryCountOf`, the header's "no entry reads it"
  sentence silently becomes false.
- `dbarts_sampler_setWeights` still does not check gaussian weights for
  finiteness or non-negativity at the flat entrance; R does.
- The R5 twins of `getTrees`/`printTrees` carry no forest selector at all, so a
  forest other than 0 is reachable only from the flat API.
- The flat `printEvery` entry still does not guard 0 (carried from the freeze).
- `setActiveRows` scans its mask twice after this slice, entry and engine.
- Pre-existing stale anchor found in passing: docs/design/feature-matrix.md:437
  cites test-capi.R:1286-1319 for the active-rows block, which lives at
  :1392-1433. The fix rides this slice's anchor pass.

## 14. Doors

- `DBARTS_NODISCARD` on the non-void stubs (the `DBARTS_IS_VOID` machinery in the
  stub branch already selects on return type, and the attribute stays out of the
  stringized list, so it is hash-invisible): turns a discarded capability 0 into
  a compile diagnostic. It cannot be taken as spelled, though - selecting on
  "non-void" would fire on all 23 int entries including the five VALUE ones
  (`kIsSampled`, `usesDart`, `family`, the version accessors), where discarding
  is legitimate, so it needs a per-entry annotation the four-column list does not
  carry. Additive and non-breaking either way, so post-1.0.
- Additive `dbarts_sampler_supports*` predicates remain available post-1.0 for a
  host that wants to probe without attempting the call.
- Drop `DBARTS_REQUIRE_EXACT_ABI` from both consumers at the coordinated 1.0
  merge (carried from the freeze).

## 15. Open decisions (VD)

1. Item 1 scope: seven entries (recommended) or the five named in the brief.
2. Item 4: rename (recommended) or transpose the contract.
3. Item 2: move `forest` (recommended) or state the rule and leave the order.
