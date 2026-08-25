# dbarts.h freeze fixes (D4) and the stub version check (D3)

Status: DESIGNED 2026-08-25, not started.

Spec: docs/plans/prerc-surface-freeze.md D3, D4, Sequencing. TODO `dbarts-h-freeze-fixes` and `stub-version-check` are ONE slice - one hash re-bake, one lockstep consumer
rebuild. Evidence: review-2026-08-24/memos/prerc-lens1-surface.md A1, A2, A6, A8, A9, A10, re-anchored live. Revised 2026-08-25 after an independent critique and two orchestrator rulings under the standing grant
(sections 2, 6). No sampling code moves; zero baseline re-records expected.

## 1. The family enum

Taken from code. Create (dbarts.h:417-419, :702-703 -> `resolveFamily`, src/R_interface_bartcore.cpp:1580-1618): `""` (dispatch on the response shape), `"gaussian"`,
`"probit"`, `"logistic"`, `"ordinal"`, `"nbinom"`, `"aft"`. `dbarts_drawLatents` (:535-541) and `dbarts_workingResponse` (:542-546) -> `augmentationFamily`
(src/R_interface_bartcore.cpp:6219-6228): `"probit"`, `"logistic"`, `"ordinal"`, `"aft"`, `"nbinom"`, `"student"`; `""` and `"gaussian"` are REFUSED there, and
`"student"` maps to `ResponseFamily::gaussian` internally while selecting a distinct law, so it needs its own enumerator. Engine `ResponseFamily` has six values
(src/bartcore/model.hpp:2580) and `SamplerShape::family` (src/bartcore/facade.hpp:76) is not total: multinomial pins `logistic` (src/bartcore/sampler.hpp:221), Student-t
reports `gaussian` (facade.hpp:104-107), heteroscedastic reports `gaussian` + `hasVarianceForest`. R's other tokens are not sampler families - `hazard*` remaps
pre-creation (R/dbarts.R:538), `twopart` aliases `hurdle.lognormal` (:409), which composes TWO samplers (:442-452), heteroscedastic is a variance forest on gaussian
(R/A_class.R:453-462) - so `multinomial` is the only engine family flat create cannot reach. Add between the leaf-model list and `#undef DBARTS_ENUMERATOR`
(dbarts.h:326-332); that `#undef` MUST move below the new list.

    #define DBARTS_FAMILY_LIST(X) \
      X(DBARTS_FAMILY_AUTO, 0) \
      X(DBARTS_FAMILY_GAUSSIAN, 1) \
      X(DBARTS_FAMILY_PROBIT, 2) \
      X(DBARTS_FAMILY_LOGISTIC, 3) \
      X(DBARTS_FAMILY_AFT, 4) \
      X(DBARTS_FAMILY_ORDINAL, 5) \
      X(DBARTS_FAMILY_NBINOM, 6) \
      X(DBARTS_FAMILY_STUDENT, 7) \
      X(DBARTS_FAMILY_MULTINOMIAL, 8)
    typedef enum { DBARTS_FAMILY_LIST(DBARTS_ENUMERATOR) } dbarts_family;

Prefix follows `DBARTS_COLUMN_*`/`DBARTS_LEAF_*`. Order: 0 = AUTO (what `""` means today; both consumers pass it), 1-6 mirror `ResponseFamily`'s declaration order so the
mapping switch reads down it, 7-8 are the two it lacks; values explicit, never renumbered, since they fold. Carried as `int`, the rule that no enum is ever a field or
parameter type holding: `dbarts_sampler_create(SEXP, SEXP, SEXP, int family)`, `dbarts_drawLatents(int family, ...)`, `dbarts_workingResponse(int family, ...)`, plus a
new entry appended LAST in `DBARTS_C_API_LIST` (48 after, 47 today) - `X(int, dbarts_sampler_family, (const dbarts_sampler* sampler), (sampler))`, prototype beside
`kIsSampled`/`usesDart` (:976-977), returned by value and so carrying no `get` (section 2). Admission is narrower per entry, every refusal naming entry and family: create
takes AUTO, GAUSSIAN, PROBIT, LOGISTIC, AFT, ORDINAL, NBINOM and refuses STUDENT, MULTINOMIAL and anything outside [0, 8]; drawLatents/workingResponse take PROBIT,
LOGISTIC, ORDINAL, AFT, NBINOM, STUDENT and refuse AUTO, GAUSSIAN, MULTINOMIAL (today's behaviour). Implement as two int -> family mappings in src/C_interface.cpp beside
`dbarts_sampler_create` (:478-484) and `augmentationArguments` (:1086-1126); `resolveFamily`/`augmentationFamily` keep their string forms for the R bridge
(R/dbarts.R:926, :1832 call `C_dbarts_bartcore_create`), so nothing below C_interface.cpp changes. `dbarts_sampler_family` is total over what the engine builds:
`shape().supportsCountsMutation` -> MULTINOMIAL (facade.hpp:110-114; dead until multinomial creation opens, the point being the ENUMERATOR existing pre-freeze), else
`shape().family` one-to-one. Never AUTO (creation resolved it), never STUDENT (a Student-t sampler's family IS gaussian - facade.hpp:107, dbarts.h:786-789). Hash:
enumerator names AND values fold (src/C_interface.cpp:442, :455-457), so add `hash = dbarts_fnv1a(hash, DBARTS_FAMILY_LIST(DBARTS_ENUMERATOR_TEXT));` to
`dbarts_apiToken()` (:453-459) after :456's leaf-model fold, without which the enum is invisible to the token. No R twin: `model` is a public R5 field (R/dbarts.R:894)
holding a `dbartsModel` whose `@family` slot carries the token (R/A_class.R:400, :414), so `sampler$model@family` is the read; residue, recorded not fixed, a hand-built
model can read `"auto"` where the engine resolved a family. Doc rewrites, all in enum terms: dbarts.h:654-663 (create's family paragraph); :124 "the two ABI enums'" and
:245 "Both ABI enums" become three; :1123-1124 and :1152-1156 name enumerators where they now spell quoted family strings. dbarts.h:774-775's "test the family before
calling" clause becomes "- read dbarts_sampler_family before calling on a sampler whose family the caller did not choose. A heteroscedastic sampler is flat-creatable (its
variance forest rides a control attribute), answers DBARTS_FAMILY_GAUSSIAN, and is still refused here, so the accessor does not by itself predict this refusal."
man/dbartsSpec.Rd:43, which documents `family` as the token for create's fourth argument and describes the empty string, names the `dbarts_family` enumerator instead,
`DBARTS_FAMILY_AUTO` for the empty string, keeping its caveat that shape dispatch is wrong for aft and logistic.

## 2. Reader-name form

Naming rule, stated once and applied everywhere here: a flat entry MIRRORS its R5 twin's name, and a reader delivering through an out pointer or an allocated object
carries `get`, while a scalar property returned by value (`numTrees`, `numForests`, `kIsSampled`, `usesDart`, `dbarts_sampler_family`) does not. That is D4's criterion,
and the R5 surface already carries `get` on every reader (R/dbarts.R:1663, 1681, 1691, 1707, 1730, 1932), so `get` is ADDED, not dropped. Four renames, old -> new
(dbarts.h list / prototype; C_interface.cpp definition; R5 twin):

- `dbarts_sampler_forestFits` -> `dbarts_sampler_getForestFits` (:512-514 / :1011; :954; `$getForestFits`, R/dbarts.R:1691)
- `dbarts_sampler_forestAmplitudes` -> `dbarts_sampler_getForestAmplitudes` (:517-519 / :1030; :975; `$getForestAmplitudes`, :1707)
- `dbarts_sampler_dispersion` -> `dbarts_sampler_getDispersion` (:533-534 / :800; :655; `$getDispersion`, :1681)
- `dbarts_sampler_forestCalibration` -> `dbarts_sampler_getForestCalibration` (:524-527 / :1083; :1015 plus the entry-naming error string at :1019; `$getCalibration`,
  :1730). The `Forest` stem stays: the struct is `dbarts_forest_calibration`.

`getLatents` (:443-444 / :790) and `getTrees` (:467-473 / :901) already comply. `dbarts_sampler_storeState` STAYS: its twin is `$storeState` (R/dbarts.R:1882), the
store/set pair's own verb, not a `get` reader. Other sites: registration and the binding asserts are X-macro expansions needing NO edit (src/R_interface.cpp:272-278,
src/C_interface.cpp:1169-1172); no Rd site, man/ never naming a flat entry; inst/NEWS.Rd:971, :973, :998, :1003; inst/tinytest/capi/consumer.c:612, :909, :993, :1036,
:1264, :1274-1276, :1287-1290, :1321, :1329, :1341 (the leg LABEL strings at :1092-1093, :1096 name channels, not entries - leave them). No consumer repo calls any of the
four.

## 3. Signature and type fixes

`printTrees` gains `useLiveTrees` in getTrees' exact position, after `numTreeIndices` and before `forest`: dbarts.h:474-480 (list), :909-917 (prototype + doc),
src/C_interface.cpp:830-861. Plumbing mirrors `bartcore_bridge::getTrees` (src/R_interface_bartcore.cpp:7700-7714): `useSaved = shape.savedTreeCapacity > 0 &&
!useLiveTrees`, with `refuseEmptyTreeStore` (C_interface.cpp:844) and the `sampleIndices` range check (:850-853) running ONLY when `useSaved`. The engine picks live vs
saved on `options_.keepTrees` (src/bartcore/sampler.hpp:1301-1334), so add a trailing `bool useLiveTrees` there and on the facade virtual and override
(src/bartcore/facade.hpp:321-327, :604-613) and change sampler.hpp:1312 to `if (useLiveTrees || !options_.keepTrees)`. The R bridge (:6095) passes `false`, so
`$printTrees` is unchanged. A facade virtual moves: `--preclean` MANDATORY. `const int*` -> `const int32_t*`, one spelling per role - index and enum arrays `int32_t`,
boolean flags `int`: dbarts.h:295 `cscColumnPointers` and :296 `cscRowIndices` (siblings at :298-302 are already `int32_t*`), with locals src/C_interface.cpp:134-135
following, and :353 `leafModel` `int*` -> `int32_t*` (the ABI's other enum-carrying array, `columnTypes`, is `int32_t*`) while `kHasHyperprior` (:349) is a FLAG and stays
`int*`; consumer cast consumer.c:972 `(int*)` -> `(int32_t*)`. All are hash-INVISIBLE (same width, offset, name - dbarts.h:132-135), which is why they land pre-freeze.
`bartcore::PredictorSource` (src/bartcore/data.hpp:188-189) stays `const int*` (`int32_t` is `int` on every supported platform), so C_interface.cpp:212-213 compiles
unchanged. `printEvery` `uint32_t` -> `size_t`: dbarts.h:489 (list), :959-960 (prototype), src/C_interface.cpp:885-888, widening the engine with it -
src/bartcore/facade.hpp:313, :591; sampler.hpp:1211; chain.hpp:212 (the `SamplerOptions` field), :1717. R side unchanged: `ParsedControl`
(src/R_interface_bartcore.cpp:178, :435) stays `uint32_t` and widens on assignment (:1846, :4920). This DOES move the hash (parameter types stringize -
src/C_interface.cpp:450-451, dbarts.h:551). Out of scope: the flat entry does not guard `printEvery == 0` (chain.hpp:1393 divides by it); the `>= 1` floor is the C++
bridge's slot read (src/R_interface_bartcore.cpp:431-433), and R/bart.R:811, R/rbart.R:100 compute `printEvery %/% n.thin`, which can reach 0 there.

## 4. D3: the stub check

Every stub calls `dbarts_stub_checkApiHash()` on its resolution branch (dbarts.h:610-621) and that check is hash EQUALITY (:572-587), so an additive MINOR would
hard-error every stubs consumer - the inverse of :29-33 and :144-148. The version accessors ARE list entries (:415-416) and are NOT stub-exempt, so the new check MUST
resolve them raw through `R_GetCCallable` as the current one resolves `dbarts_apiHash` (:577-578); calling the stubs would re-enter it. Replace :572-587 with:

    static inline void dbarts_stub_checkApi(void) {
      static int dbarts_stub_apiChecked = 0;
      if (dbarts_stub_apiChecked) return;
      int (*majorFn)(void) = (int (*)(void)) (void (*)(void)) R_GetCCallable("dbarts", "dbarts_apiMajorVersion");
      int (*minorFn)(void) = (int (*)(void)) (void (*)(void)) R_GetCCallable("dbarts", "dbarts_apiMinorVersion");
      int major = majorFn(), minor = minorFn();
      if (major != DBARTS_C_API_MAJOR || minor < DBARTS_C_API_MINOR)
        Rf_error("dbarts C API version mismatch: this package was built against %d.%d, the installed "
                 "dbarts provides %d.%d; rebuild this package against the installed dbarts",
                 DBARTS_C_API_MAJOR, DBARTS_C_API_MINOR, major, minor);
    #ifdef DBARTS_REQUIRE_EXACT_ABI
      /* today's hash-equality check and its :581-585 message, verbatim, with dbarts_apiHash
         resolved raw through R_GetCCallable exactly as majorFn/minorFn are above */
    #endif
      dbarts_stub_apiChecked = 1;
    }

Rename the call at :614. Wrap :114-115 in `#ifndef DBARTS_C_API_MAJOR`/`#ifndef DBARTS_C_API_MINOR` on the hash's own pattern (:152-154), carrying the sentence at
:150-151, so the handshake is testable. Header doc rewrites: :29-33 (the stubs enforce the two-component handshake, the hash is opt-in), :137-151, :563-571, :643-651 ("A
DBARTS_USE_STUBS consumer does not call it" -> "...unless it defines DBARTS_REQUIRE_EXACT_ABI"). Consumer docs: man/dbarts-package.Rd:43 and inst/NEWS.Rd:566-567 ("as an
exact signature check") -> "as an opt-in exact-ABI check, which moves on additive releases as well"; inst/NEWS.Rd:994-995 ("an exact signature-lockstep check") -> "an
opt-in exact-ABI lockstep check". Leave inst/NEWS.Rd:1020 (accurate as written).

## 5. Hash re-bake

Two literals move together, in the header edit's commit: `DBARTS_C_API_HASH` (dbarts.h:153) and the `dbarts_apiSignatureToken` assert (src/C_interface.cpp:461); the
combined assert is :465-470. `DBARTS_C_API_MAJOR`/`MINOR` DO NOT MOVE - no version constant increments before 1.0-0 (dbarts.h:111-113). Neither assert prints its value,
so derive both by temporarily adding, after `dbarts_apiToken()` (C_interface.cpp:459), `template <std::uint64_t> struct DbartsShowToken;` and
`DbartsShowToken<dbarts_apiSignatureToken> a; DbartsShowToken<dbarts_apiToken()> b;`; the compile diagnostic names both in decimal, convert with `printf '0x%016xULL\n'
<decimal>` (NOT `%016llx`: zsh's builtin and /usr/bin/printf both reject it), bake, remove the probe. Pin: add the outgoing `"0x66d33f1613892406"` to test-capi.R's
`expect_false` block (:73-81), a comment naming what moved, and set :84 to the new literal; :90's `expect_equal(versions, c(1L, 0L))` is unchanged. The outgoing literal
also appears at docs/design/threaded-predict.md:112, whose sentence says it is what dbarts.h:153 bakes LIVE - update that one. :274 records what the threaded-predict arc
itself re-signed to and stays, as do docs/design/multinomial-mutation-arc.md:334, :350, in its sections 1-4, declared at :6-8 to be the pre-arc proposal and "not
statements about the live code"; the docs/plans mentions are not anchor-checked, leave them.

## 6. Consumer migration (lockstep, after dbarts installs clean)

Both consumers keep lockstep by opting in: define `DBARTS_REQUIRE_EXACT_ABI` beside `DBARTS_USE_STUBS`, with the one-line comment "pre-release lockstep guard: MAJOR/MINOR
do not move before 1.0-0, so only the token catches a stale binary; dropped at the coordinated 1.0 merge". Both then DROP their hand-rolled `dbarts_apiHash() !=
DBARTS_C_API_HASH` clauses. Add to the dbarts RELEASE procedure entry in TODO: "drop DBARTS_REQUIRE_EXACT_ABI from stan4bart and treatSens Makevars". stan4bart,
/Users/vdorie/Repositories/stan4bart branch `bartcore`, 24 entries, 4 edits. (1) src/init.cpp:146-153 `getBARTFamily` returns `int`, not `const char*`: `""`/`"auto"` ->
`DBARTS_FAMILY_AUTO`, `"gaussian"` -> `DBARTS_FAMILY_GAUSSIAN`, probit/logistic/aft/ordinal/nbinom to theirs. :151-152 returns whatever string sits on the attribute, and
`dbartsModel@family` admits 8 tokens including `"multinomial"` (R/A_class.R:453-462), so the default branch MUST `Rf_error` naming the token, never fall through to AUTO.
Its call sites (:197-199, :366-367) need no text change. (2) src/init.cpp:432 `dbarts_sampler_printTrees`: insert `0` for `useLiveTrees` before the trailing `0` forest;
:495-497 already passes `useLiveTrees` to `getTrees`. (3) src/init.cpp:971-972: drop the hash clause and the stale comment at :966-970 scheduling its removal. (4)
`-DDBARTS_REQUIRE_EXACT_ABI` appended to `PKG_CPPFLAGS` in BOTH src/Makevars.in:1 and src/Makevars.win:1. No `printEvery` edit - the `100` literals at :200, :368 widen;
no `dbarts_results` edit (src/bart_util.hpp:27, :47); no renamed reader called. treatSens, /Users/vdorie/Repositories/treatSens/.claude/worktrees/dbarts-1.0 branch
`dbarts-1.0`, 11 entries, 4 edits: src/bartTreatmentModel.cpp:62 `"probit"` -> `DBARTS_FAMILY_PROBIT`; src/sensitivityAnalysis.cpp:259 `""` -> `DBARTS_FAMILY_AUTO`;
src/R_interface.cpp:452-453 drop the hash clause; and `#define DBARTS_REQUIRE_EXACT_ABI` immediately above each `#define DBARTS_USE_STUBS`, since this repo defines the
macro per translation unit rather than in Makevars - src/R_interface.cpp:27, src/bartTreatmentModel.cpp:10, src/sensitivityAnalysis.cpp:32, three sites with the comment
on the first. No `printEvery` edit (:64, :261); no printTrees, getTrees or renamed-reader call sites. Gate: `R CMD INSTALL . --preclean` in dbarts then each consumer,
then `tinytest::test_package("stan4bart")` and treatSens `testthat::test_local()` (runner tests/testthat.R), both green before this slice is landed.

## 7. Tests, NEWS and re-anchoring inside dbarts

consumer.c: `capi_create` (:130-139) takes the enum through a static string -> enum table, so test-capi.R's 46 `capi_create` call sites keep their string arguments; add
`capi_create_raw_family(control, model, data, familyInt)` passing an unmapped `int` for the refusal probes; the same table serves `capi_draw_latents` (:623-636) and
`capi_working_response` (:638-648), whose 3 call sites keep their strings; add `capi_family_constants()` (named integer vector of all nine) and
`capi_sampler_family(ptr)`; `capi_print_trees` (:471-479) gains `useLiveTrees`; apply the four renames and the `(int32_t*)` cast at :972. tests/cpp: the facade
`printTrees` virtual gains a parameter, so test_facade.cpp:184-188 (its `SPY_VOID` declaration) takes it while test_facade.cpp:496 and test_sampler.cpp:6345 (direct
calls) pass `false`; separately test_facade.cpp:179's `SPY_VOID(setVerbose, (bool v, std::uint32_t e), (v, e))` must widen `e` to `std::size_t` or the spy stays abstract.
test-capi.R additions. Enum round trip: `capi_family_constants()` equals 0:8 in header order, and `capi_sampler_family()` gives the matching enumerator on the gaussian,
probit, logistic, ordinal, nbinom and aft samplers the file already builds, with `resid.dist = student()` gaussian and the heteroscedastic sampler at :222-228 both
GAUSSIAN. Refusals: `capi_create_raw_family(..., 999L)` and `..., -1L` on the range, `..., 7L` and `..., 8L` naming the family, `capi_draw_latents` at AUTO and at
GAUSSIAN by name. `useLiveTrees`: on a `keepTrees` sampler TRUE and FALSE both print and differ; with storage on and nothing recorded TRUE prints while FALSE raises the
empty-store refusal (:598-620 is the model). D3 probes are built by REWRITING the existing stale-token block at inst/tinytest/test-capi.R:1710-1759, NOT from the :29-56
harness: that block compiles the consumer with `-DDBARTS_C_API_HASH=0x0123456789abcdefULL` alone (:1721) and expects "dbarts C ABI mismatch", which D3 makes pass
silently. Rewrite it into four arms, each with its own temp dir and object name (the loader caches by path) and the block's CI guard (`stop` under CI, else skip): (a)
that flag alone, now expecting NO error - the half proving the gate discriminates; (b) it plus `-DDBARTS_REQUIRE_EXACT_ABI` at :1721, expecting "dbarts C ABI mismatch";
(c) `-DDBARTS_C_API_MAJOR=99` and (d) `-DDBARTS_C_API_MINOR=99`, both expecting "dbarts C API version mismatch" (the second exercises the minor floor). The default-flags
consumer at :29-56 covers the clean path, so both macro modes compile. Plus section 5's hash pin edits. NEWS: one new `\item` in the 1.0-0 `\subsection{C API}`
(inst/NEWS.Rd:956-1023, before the `\itemize` close at :1022) covering `dbarts_family`, `dbarts_sampler_family`, the four renames, `printEvery` as `size_t`,
`printTrees`'s `useLiveTrees`, and the stub check now gating on major/minor with `DBARTS_REQUIRE_EXACT_ABI` as the opt-in - in addition to the three rewordings in
sections 2 and 4. Re-anchoring, in the SAME commit as the code edits: 22 docs/design anchors point into dbarts.h and 17 into src/C_interface.cpp, and the family list
alone pushes docs/design/public-surface.md:411 (`dbarts_sampler_setCallback` at dbarts.h:439, STRICT) onto `getLatents`. Run `Rscript tools/check-doc-freshness.R .` after
the edits and re-align every strict miss from the `git diff -U0` line map, editing the docs/design anchors in place so each file's line count is invariant; the
landed-hash stamps in docs/design/feature-matrix.md and threaded-predict.md belong to the records commit, not this one. Gate battery (CLAUDE.local.md): `R CMD INSTALL .
--preclean`; `tinytest::test_package("dbarts")`; `cd tests/cpp && make && ./test_bartcore`; `benchmarks/R/equivalence.R compare` against equivalence-736bfb05.rds,
bcf-equivalence-6e3b9fb8.rds and multinomial-equivalence-4d9a3337.rds, all bitwise since no sampling moves; `benchmarks/R/bench-sampler.R compare
bench-sampler-ab1dc52.csv` on a quiet machine; `tools/check-doc-freshness.R`; `air format --check .`; `R CMD check`. Delete the benchmarks/kernels binaries first.

## 8. Open sub-choices

None. Both are settled by VD, 2026-08-25: the naming rule and the fourth rename in section 2, the consumer lockstep macro and its release-time removal in section 6.
