One of four adversarial pre-review reports read ahead of the manual review; VD's rulings on its findings are recorded in docs/plans/pre-review-cleanup.md.

# Review-chain audit at bartcore 3080a9c5

READ-ONLY. Targets: docs/plans/bartcore-review-tour.md (975 lines, stamped
ae5b91d8) and docs/design/feature-matrix.md (1,037 lines, stamped 52c10e02).
Every disposition below was checked against the live tree at 3080a9c5, not
against a landing note.

Scope note the audit turned up first: the tour is stamped at ae5b91d8, but
**46 commits** have landed since, not four. Besides today's four slices it
missed the dbarts.h freeze (6446ddce), the predict-surface slice (7b3ac6bf,
71cc7133, befc8f45, ed43deef), composition-refusals (fe0b3292, de18ef2b) and
the fit-entry refusals (fcbbc478). Several "Known open" items the task expected
to be closed today were already closed before today.

## Citations of the tour by line

    [[docs/plans/review-2026-08-24/memos/prerc-lens2-backlog.md:73@fcbbc478]]   -> tour:290
    [[docs/plans/review-2026-08-24/memos/prerc-lens3-external.md:49@fcbbc478]]  -> tour:463-467
    [[docs/plans/INDEX.md:254@fcbbc478]]                                        -> no line
    [[docs/plans/release-candidate-review.md:594@fcbbc478]], [[docs/plans/release-candidate-review.md:1094@fcbbc478]]              -> no line

Only **:290** and **:463-467** are line-cited. [[docs/plans/release-candidate-review.md:290@fcbbc478]] still holds its cited
claim (the per-person-period LOO trap, "stated nowhere"). [[docs/plans/release-candidate-review.md:463-467@fcbbc478]] still holds
the bartCause paragraph; faf1d167 already amended it in place without moving a
line. Any tour edit outside those two spans is free to change line counts;
edits inside them should stay in-line, or the two memos re-pointed.

feature-matrix.md is line-cited from docs/design and its anchors are certified,
so every replacement below is in-line and character-count-compatible.

---

## 1. Summary counts

| disposition | tour | matrix | total |
|---|---|---|---|
| STALE-FALSE      | 18 | 6  | 24 |
| STALE-INCOMPLETE | 12 | 4  | 16 |
| OK (checked, holds) | 9 | 5 | 14 |
| unverified | 3 | 0 | 3 |

Broad finding, matrix: the 2026-08-26 re-align (bfc4571d) shifted anchors by
the diff line map rather than re-deriving them by content, so its certification
stamp at [[docs/plans/release-candidate-review.md:37@bfc4571d]] overstates what was done. Three anchor families do not hold - the
flat-C column of table 1, the pointwise-loglik column of table 3, and the
`$setCalibration` half of the calibration column - plus four scattered cites.
Every engine anchor I sampled (chain.hpp, model.hpp, combiner.hpp, sampler.hpp)
does hold; facade.hpp is uniformly +2.

Broad finding, tour: anchor rot is pervasive (~40 of ~120 cited line numbers
are wrong), and the entry maps for Stop 6 and Stop 9 are wrong wholesale.

---

## 2. Per-passage table

### 2a. docs/plans/bartcore-review-tour.md

| line | claim (quoted/abridged) | disposition | replacement |
|---|---|---|---|
| 3-4 | "anchored at branch tip **ae5b91d8**" | STALE-FALSE | `**3080a9c5**` |
| 7-9 | "cb290550 (main tip; bartcore is 1273 commits ahead) / 680 files, +229,948 / -32,374 / docs: +80,194 / -0 (205 files, ~35%...)" | STALE-FALSE | `merge base: cb290550 (main tip; bartcore is 1319 commits ahead)` / `scale:      828 files, +263,944 / -32,412` / `docs:       +108,721 / -0  (352 files, ~41% of the diff, zero deletions)` |
| 12-13 | "The file count fell since the previous refresh (733 -> 680)" | STALE-FALSE | "The file count rose again since the last refresh (680 -> 828), all of it docs." |
| 23 | "Every line number below was grep-verified against ae5b91d8 at write time." | STALE-FALSE | "Every line number below was grep-verified at 3080a9c5. Prefer the symbol name where one is given: the numbers rot within days." |
| 25-29 | "**One stamp to distrust.** `[[docs/design/feature-matrix.md:37@3080a9c5]]` says its anchors are verified by content 'against the tree at 54dec2ab'. That sha is not an ancestor..." | STALE-FALSE | The stamp now reads 52c10e02, a true ancestor, so the 54dec2ab paragraph is dead. Replace with: "**One stamp to read narrowly.** `[[docs/design/feature-matrix.md:37@23b9cde7]]` certifies its anchors by content at 52c10e02, but the 2026-08-26 re-align shifted them by the diff line map rather than re-deriving them. Its engine anchors hold; its flat-C column (table 1), its pointwise-loglik column (table 3) and its `$setCalibration` cites do not. Cell VALUES are separately adjudicated and were not re-checked." |
| 33-36 | "Two are reshaped this refresh: Stop 2 ... Stop 9 ..." | STALE-INCOMPLETE | Refresh framing is two refreshes old; drop the sentence or restate against the regeneration (section 5). |
| 38-49 | stop diff sizes | STALE-FALSE | 1: `+21,087 / -8,718`; 3: `3,027`; 6: `+8,339 / -385`; 7: `+6,140 / -1,295`; 8: `+23,432 / -3,251`; 10: `~89,948`. (2, 4, 5 unchanged: 4,923 / 1,857 / 2,085.) |
| 65 | "`docs/design/feature-matrix.md` (1,032 lines)" | STALE-FALSE | `(1,037 lines)` |
| 68 | "a 'Gaps' section ([[docs/design/feature-matrix.md:924@52c10e02]])" | STALE-FALSE | `([[docs/design/feature-matrix.md:929@52c10e02]])` |
| 100-103 | "Two doc-vs-code divergences still open" (ColumnStore 16-bit, not arena-allocated) | unverified | see section 4 |
| 105-113 | "**One document that contradicts itself where a reviewer will land.** `[[docs/plans/variance-forest-mutation-routing.md:415-426@52c10e02]]` opens ... with 'CLOSED: built at c95a5e83' and then ... continues with the original present-tense description of the hole" | STALE-FALSE | That door was rewritten: [[docs/plans/variance-forest-mutation-routing.md:416-426@c95a5e83]] now closes cleanly ("Nothing remains open under this bullet ... the description of the hole that preceded it is retired with the door"). Cut the whole paragraph. The `[[feature-matrix.md:1016-1022@c95a5e83]]` cite in it is also stale (now [[feature-matrix.md:1017-1027@c95a5e83]]). |
| 125-141 | Stop 1 read-order anchors | STALE-FALSE | `facade.hpp (908)` - `SamplerBase` [[feature-matrix.md:138-413@c95a5e83]] (the `sed -n '138,413p' ... grep -c '= 0;'` count is still 59), `SamplerFacade` [[feature-matrix.md:415@c95a5e83]], `createConstantLeafSampler` [[feature-matrix.md:732@c95a5e83]], `createSampler` [[feature-matrix.md:789@c95a5e83]], `createSamplerOverStore` [[feature-matrix.md:842@c95a5e83]], `createAmplitudeSampler` [[feature-matrix.md:870@c95a5e83]], `createMultinomialSampler` [[feature-matrix.md:896@c95a5e83]]. sampler.hpp and chain.hpp anchors all hold as written. |
| 169 | "`stateFormatVersion` is 3 (`[[src/R_interface_bartcore.cpp:6380@c95a5e83]]`)" | STALE-FALSE | `([[src/R_interface_bartcore.cpp:6446@c95a5e83]])` |
| 228 | "`setWeights`/`setSigma` refused by name (`[[src/R_interface_bartcore.cpp:2773@c95a5e83]]`, [[src/R_interface_bartcore.cpp:2898@c95a5e83]])" | STALE-FALSE | `[[src/R_interface_bartcore.cpp:2791@c95a5e83]]`, `[[src/R_interface_bartcore.cpp:2928@c95a5e83]]` |
| 272 | "DART is refused by name (`[[R/bart.R:901@c95a5e83]]`)" | STALE-FALSE | `[[R/bart.R:902@c95a5e83]]` |
| 298 | "`hurdleLogLik` (`[[R/generics.R:2045@c95a5e83]]`)" | STALE-FALSE | `[[R/generics.R:2344@c95a5e83]]` |
| 327-331 | "`heteroscedasticScale` (`[[R/generics.R:32@c95a5e83]]`) ... in `pointwiseLogLikelihood` ([[R/generics.R:95@c95a5e83]]) and, via `ppdNoiseScale` ([[R/generics.R:2754@c95a5e83]]), in `sampleFromPPD` ([[R/generics.R:2785@c95a5e83]])" | STALE-FALSE | `[[R/generics.R:32@c95a5e83]]` and `[[R/generics.R:95@c95a5e83]]` hold; `ppdNoiseScale` is `[[R/generics.R:3115@c95a5e83]]`, `sampleFromPPD` `[[R/generics.R:3146@c95a5e83]]` |
| 338 | "(`[[src/R_interface_bartcore.cpp:7118@c95a5e83]]`, [[src/R_interface_bartcore.cpp:7433@c95a5e83]])" | STALE-FALSE | `[[src/R_interface_bartcore.cpp:7184@c95a5e83]]`, `[[src/R_interface_bartcore.cpp:7499@c95a5e83]]` ([[feature-matrix.md:1027@c95a5e83]] already carries these) |
| 379 | Stop 3 Known open: "`monotone-leaf-branch-fill` is a live TODO entry." | STALE-FALSE | "`monotone-leaf-branch-fill` is CLOSED, not adopted: benched 2026-08-26 with no measurable effect (`fillBottom` 0.18% of a monotone run, within noise), bitwise; the arm is archived at `archive/monotone-branch-fill` (98067fb3) and the memo is `docs/plans/review-2026-08-24/memos/monotone-branch-fill-bench.md`. Door left: `coneProbability`'s quadrature ([[TODO:126-129@23b9cde7]])." |
| 451 | "`test-predict-blend.R` (390 lines)" | STALE-FALSE | `(417 lines)` |
| 466 | bartCause fix "(fixed 765a596 on its dbarts-1.0 branch)" | STALE-INCOMPLETE | Correct as far as it goes; add the verdict. In-line (this span is cited): `... - a bartCause-side follow-up (fixed 765a596 on its dbarts-1.0 branch; suite 692/0), which is why a grep here` |
| 478-490 | Stop 6 entry map ("`R_interface_bartcore.cpp` is 7,759 lines", `_create` [[TODO:3489@23b9cde7]] ... guards [[TODO:2265@23b9cde7]] ...) | STALE-FALSE | 7,825 lines. `bartcore_create` [[TODO:3541@23b9cde7]], `_createDataHandle` [[TODO:3565@23b9cde7]], `_createFromHandle` [[TODO:3619@23b9cde7]], `_createBCF` [[TODO:3791@23b9cde7]]; `_run` [[TODO:4307@23b9cde7]], `_runWithCallback` [[TODO:4551@23b9cde7]], `_growFromRoot` [[TODO:4632@23b9cde7]]; `_setOffset` [[TODO:4649@23b9cde7]], `_setResponse` [[TODO:4669@23b9cde7]], `_setSigma` [[TODO:4693@23b9cde7]], `_setData` [[TODO:4700@23b9cde7]], `_setWeights` [[TODO:4935@23b9cde7]]; `_setPredictor` [[TODO:5119@23b9cde7]], `_updatePredictor` [[TODO:5167@23b9cde7]], `_setActiveRows` [[TODO:4074@23b9cde7]]; `_setForestBasis` [[TODO:3981@23b9cde7]], `_getForestAmplitudes` [[TODO:4118@23b9cde7]], `_getForestVariableCounts` [[TODO:4280@23b9cde7]]; `_getCalibration` [[TODO:4205@23b9cde7]], `_setCalibration` [[TODO:4260@23b9cde7]]; `_storeState` [[TODO:5684@23b9cde7]], `_setState` [[TODO:5689@23b9cde7]], `_installForests` [[TODO:5702@23b9cde7]]; `_predict` [[TODO:5837@23b9cde7]], `_predictPerForest` [[TODO:5938@23b9cde7]], `_getTrees` [[TODO:6022@23b9cde7]], `_getLatents` [[TODO:6171@23b9cde7]]. Guards: `refusedAmplitudeFamilyReason` [[TODO:2266@23b9cde7]], `refuseMultiForestMutation` [[TODO:2618@23b9cde7]], `refuseMultiForestResponseMutation` [[TODO:2647@23b9cde7]], `refuseBinaryWeightChange` [[TODO:2791@23b9cde7]], `refuseUndefinedTestFits` [[TODO:2890@23b9cde7]], `refusePinnedSigmaChange` [[TODO:2928@23b9cde7]], and NEW `refuseNonBinaryMask` [[TODO:2946@23b9cde7]] (shared with the flat entry since 9df0cb50). |
| 498 | "`refusePinnedSigmaChange`'s own comment ([[TODO:2885-2897@9df0cb50]])" | STALE-FALSE | `([[TODO:2915-2927@9df0cb50]])` |
| 511-514 | "`dbarts.h` is 1,158 lines ...; `C_interface.cpp` 1,175. Constants: `DBARTS_C_API_MAJOR 1` ([[TODO:103@9df0cb50]]), `DBARTS_C_API_MINOR 0` ([[TODO:104@9df0cb50]]), `DBARTS_C_API_HASH 0x66d33f1613892406ULL` ([[TODO:142@9df0cb50]])." | STALE-FALSE | "`dbarts.h` is 1,307 lines (was 477); `C_interface.cpp` 1,315. Constants: `DBARTS_C_API_MAJOR 1` ([[TODO:143@9df0cb50]]), `DBARTS_C_API_MINOR 0` ([[TODO:146@9df0cb50]]), `DBARTS_C_API_HASH 0x5a32aa4cd3872d55ULL` ([[TODO:189@9df0cb50]])." Two re-bakes since the tour was written (freeze -> 0x0939c022...; capi-shape -> 0x5a32aa4c...). |
| 516-534 | Stop 7 body: threaded predict as the whole story | STALE-INCOMPLETE | Threaded predict is now the *older* of two ABI events. Must add the capi-shape slice (9df0cb50): (i) **a three-class return doctrine** stated at [[CAPI:53-69@9df0cb50]] - VALUE (12 entries), TRANSACTION (`setPredictor`/`updatePredictor`), CAPABILITY STATUS (16); (ii) **seven setters void -> int** (`setResponse`, `setOffset`, `setWeights`, `setSigma`, `setTestPredictors`, `setTestOffset`, `predict`), 0 meaning "this sampler cannot, nothing was touched", a raise meaning "this call is wrong"; (iii) **`forest` moved to argument 2** on `getTrees` ([[CAPI:1023@9df0cb50]]) and `printTrees` ([[CAPI:1034@9df0cb50]]), the last two that took it after `useLiveTrees` - a pointer-to-integer change at position 2 that fails to compile in C++ and may only warn in C, with the hash re-bake as the real backstop; (iv) `setForestBasis`'s `basis` renamed **`basisRowMajor`** ([[CAPI:1137@9df0cb50]]), naming the one documented exception to the column-major rule ([[CAPI:106-110@9df0cb50]]); (v) `setActiveRows` now **raises** on a non-binary mask instead of answering 0 ([[src/C_interface.cpp:1174-1189@9df0cb50]], shared `refuseNonBinaryMask`), so its 0 is purely a capability answer; (vi) dbarts.h no longer defines `USE_FC_LEN_T` (nor includes `<Rversion.h>`) - a consumer that relied on that pull-in now includes it itself. Entry count stays 48; MAJOR/MINOR do not move. Landing note: `docs/plans/capi-shape.md`. |
| 521-523 | "`dbarts_sampler_predict` gained `size_t numThreads` as its new last parameter ([[src/C_interface.cpp:861@9df0cb50]])" | STALE-FALSE | It is the fourth of five parameters (`out` is last), at [[src/C_interface.cpp:990@9df0cb50]], and the entry now RETURNS `int`: a capability status, 1 with `out` written or 0 leaving `out` untouched where the sampler's blend is undefined. |
| 526-527 | "Consumer migration: stan4bart's one call site gains a `0`; bartCause needs zero lines" | STALE-INCOMPLETE | Still the threaded-predict half only. capi-shape adds: both stan4bart `getTrees`/`printTrees` sites fail to compile until the forest argument moves to position 2; every consumer with `DBARTS_REQUIRE_EXACT_ABI` hard-stops on the new hash; a consumer relying on dbarts.h for `USE_FC_LEN_T`/`<Rversion.h>` must include them itself. treatSens is **not** R-API-only (its main worktree has `src/`); bartCause is (no `src/` at all). See capi-shape.md section 11. |
| 528 | "`test-capi.R` grew 379 -> 1,759 lines" | STALE-FALSE | `379 -> 2,012 lines` |
| 536-553 | Stop 7 "Known open" | STALE-FALSE (one item) + STALE-INCOMPLETE | Delete "a fractional `n.threads` truncates silently (`validatePredictThreads` `[[R/generics.R:234@9df0cb50]]` - `is.numeric(2.7)` passes every predicate and `as.integer` truncates)": it is fixed. `validatePredictThreads` ([[R/generics.R:234-248@9df0cb50]]) now refuses on `n.threads != round(n.threads)` and echoes the value, and `coerceOrError` ([[R/utility.R:173-180@9df0cb50]]) refuses a fractional integer at every other count. The rest of the list (D1, D3, unjoined workers on a `std::thread` throw, `Rf_error` after the join, `benchmarks/R/bench-predict.R` absent - re-confirmed absent, the tuned cutoff, and the four undocumented control attributes) all still hold. Add from capi-shape: a DISCARDED capability 0 now leaves the sampler unconditioned and the run continuing - a quieter failure than the old longjmp; `DBARTS_NODISCARD` on the non-void stubs is a recorded door, not taken. |
| 561-563 | R file sizes "R/generics.R 2,985, R/bart.R 2,896, R/dbarts.R 2,093, R/model.R 1,812, R/data.R 1,729, R/rbart.R 1,309, R/bartcore.R 1,107" | STALE-FALSE | `R/generics.R 3,346, R/bart.R 2,975, R/dbarts.R 2,178, R/data.R 1,828, R/model.R 1,812, R/rbart.R 1,309, R/bartcore.R 1,112` |
| 576-584 | "(`bartRedirectedFamilies` `[[R/bart.R:2621@9df0cb50]]`, `refuseBartRedirectedFamily` [[R/bart.R:2634@9df0cb50]])" and "`checkFamilyUnsupportedArgs` (`[[R/bart.R:620@9df0cb50]]`)" | STALE-FALSE | `[[R/bart.R:2699@9df0cb50]]`, `[[R/bart.R:2712@9df0cb50]]`; `[[R/bart.R:618@9df0cb50]]` |
| 586 | "`$n.chains` is set unconditionally on every fit that keeps a sampler." | STALE-FALSE | "`$n.chains` is on every fit, kept sampler or not: `packageRbartResults` sets it unconditionally ([[R/rbart.R:1295@9df0cb50]]) since 52c10e02, closing the gap where an `rbart_vi()` fit that kept its samplers carried none and `fitNChains` had to infer it." |
| 590-595 | "`tools/check-rc-codoc.R` (259 lines, wired into `lint.yaml`)" | OK | 259 lines, confirmed |
| 597-601 | Stop 8 Known open: "**`predict(bases =)`, `predict(sample =)` and `fitted(combineChains =)` on the own-class fits are still silent no-ops** - explicitly outside the judgement table ... do not assume the generics work closed every swallow." | STALE-FALSE | All three are refused by name since d48aef8a. Replace with: "**Foreign arguments are now refused by name across the surface.** `refuseUnusedGenericArgs`/`foreignArgsFor` ([[R/generics.R:2028@d48aef8a]], [[R/generics.R:2050@d48aef8a]]) compose four derived reason tables - `predictForeignReasons` [[R/generics.R:2109@d48aef8a]], `extractForeignReasons` [[R/generics.R:2118@d48aef8a]], `fittedForeignReasons` [[R/generics.R:2128@d48aef8a]], `residualsForeignReasons` [[R/generics.R:2141@d48aef8a]], plus `survivalProbabilitiesForeignReasons` [[R/generics.R:2154@d48aef8a]] - so `predict(bases =)`, `predict(sample =)` and `fitted(combineChains =)` on the own-class fits now name the channel that serves the caller instead of vanishing. One inert argument survives on purpose: `as_draws_array`/`as_draws_df.bartMultinomial`'s `vars` ([[R/diagnostics.R:226@d48aef8a]], [[R/diagnostics.R:229@d48aef8a]]), documented rather than refused because the `posterior` generics fix the signature." |
| 601-604 | "`setForestBasis(k, ~var)` still evaluates the formula in `environment(basis)` ..." | OK | Still true: `$setForestBasis` calls `evaluateForestBasis(basis)` with `data = NULL` ([[R/dbarts.R:1476@d48aef8a]] -> [[R/model.R:808-820@d48aef8a]]). |
| 604-606 | "`plotTree`'s dead padding branch stays out"; "`dbartsData(counts =)` with a formula still reaches its own later refusal" | OK | The counts/formula refusal is live at [[R/data.R:1287@d48aef8a]]. plotTree padding not re-checked in detail. |
| 606-608 | "**`par(mfrow)` restoration is now partial**: ... but `plot.bart` and `plot.rbart` set `mfrow` through the shared `plotSigmaTrace` (`[[R/plot.R:9-12@d48aef8a]]`) and still leak it." | STALE-FALSE | Fixed at fcbbc478. `plot.bart` ([[R/plot.R:67-68@fcbbc478]]) and `plot.rbart` ([[R/plot.R:133-134@fcbbc478]]) both `par(no.readonly = TRUE)` + `on.exit(par(oldpar))`; `plotSigmaTrace` gained `setLayout =` so a caller with its own layout is not reset. Replace with: "**`par(mfrow)` is restored by every plot method.** All eight save and restore (`R/plot.R` [[R/plot.R:67@fcbbc478]], [[R/plot.R:133@fcbbc478]], [[R/plot.R:232@fcbbc478]], [[R/plot.R:314@fcbbc478]], [[R/plot.R:374@fcbbc478]], [[R/plot.R:423@fcbbc478]], [[R/plot.R:521@fcbbc478]] and `plotTree`); `plotSigmaTrace`'s `setLayout = FALSE` (:9) is the one caller-owned-layout escape." |
| 597-608 | Stop 8 Known open, as a list | STALE-INCOMPLETE | Missing every surface-refusals fact a reviewer needs: `fitted()`'s positional slot 3 is now `ci.level` on all six methods, and `fitted.bartHurdle` lost `sample` ([[R/generics.R:910@fcbbc478]], 1185, 1515, 1804, 2363, 2889); `plot.bartMultinomial` reordered to `(x, plquants, cols, ...)` ([[R/plot.R:226@fcbbc478]]) matching its five siblings, with `plot.pdbart` deliberately a fourth shape; `type = "class"` with `ci.level`, `ci.level` on `type = "forest"`, forest selection outside the forest arm and `residuals(ci.level =)` all refused; `$getSigmas`/`$getSumsOfSquaredResiduals` refuse a vestigial `result` ([[R/dbarts.R:1675@fcbbc478]], [[R/dbarts.R:1695@fcbbc478]]); `coerceOrError` refuses a fractional integer ([[R/utility.R:173-180@fcbbc478]]). Landing note: `docs/plans/surface-refusals.md`. |
| 624 | "the ladder is `resolveTreePrior` (`[[R/bart.R:484@fcbbc478]]`)" | STALE-FALSE | The symbol does not exist. It is `resolveDartShorthand` (`[[R/bart.R:467@fcbbc478]]`), with the collision check at `[[R/bart.R:894-905@fcbbc478]]`. |
| 629 | "`refuseUnsupportedAmplitudeComposition`, `[[src/R_interface_bartcore.cpp:2323@fcbbc478]]`" | STALE-FALSE | `[[src/R_interface_bartcore.cpp:2316@fcbbc478]]` |
| 636-644 | setSigma: "R5 `[[R/dbarts.R:1504@fcbbc478]]`, bridge [[R/dbarts.R:4627@fcbbc478]] ... `[[src/R_interface_bartcore.cpp:2898@fcbbc478]]`" | STALE-FALSE | R5 `[[src/R_interface_bartcore.cpp:1504@fcbbc478]]` holds; bridge is `[[src/R_interface_bartcore.cpp:4693@fcbbc478]]`; the guard is `[[src/R_interface_bartcore.cpp:2928@fcbbc478]]`. Add that the flat entry now returns a capability status rather than raising. |
| 649-657 | "`estimateSigmaFromLinearModel` (`[[R/utility.R:542@fcbbc478]]`) ... floored at ... ([[R/utility.R:617@fcbbc478]]) ... re-checked *after* the fallback ([[R/utility.R:618@fcbbc478]]) ... the sparse branch returns at [[R/utility.R:562@fcbbc478]]" | STALE-FALSE | `estimateSigmaFromLinearModel` is `[[R/utility.R:551@fcbbc478]]`; the three interior cites shift with it - re-derive rather than offset. |
| 675-676 | "`[[R/dbarts.R:1398@fcbbc478]]` -> bridge [[R/dbarts.R:4022@fcbbc478]] -> `Chain::setActiveRows` [[chain.hpp:1666@fcbbc478]] -> [[facade.hpp:368@fcbbc478]]" | STALE-FALSE | bridge is `[[facade.hpp:4074@fcbbc478]]` (`[[facade.hpp:4022@fcbbc478]]` is `bartcore_setForestWeights`); facade is `[[facade.hpp:370@fcbbc478]]`. [[chain.hpp:1666@fcbbc478]] holds. |
| 679-684 | active-rows "Open:" list | STALE-INCOMPLETE | Add: the flat entry's value refusal now RAISES (shared `refuseNonBinaryMask`, `[[src/R_interface_bartcore.cpp:2946@fcbbc478]]`), so a flat-C 0 means only "this family carries no mask". |
| 688 | "R5 `installTrees(donor, samples)` (`[[R/dbarts.R:1888@fcbbc478]]`)" | STALE-FALSE | `[[R/dbarts.R:1973@fcbbc478]]` (`[[R/dbarts.R:1888@fcbbc478]]` is inside `reapplyForestWeights`) |
| 702-705 | "`extract.bart` (`[[R/generics.R:425@fcbbc478]]` ...), `.bartMultinomial` [[R/generics.R:966@fcbbc478]], `.bartOrdinal` [[R/generics.R:1297@fcbbc478]], `.bartNegbin` [[R/generics.R:1576@fcbbc478]], `.bartHurdle` [[R/generics.R:2001@fcbbc478]], `.rbart` [[R/generics.R:2379@fcbbc478]]" | STALE-FALSE | `extract.bart` [[R/generics.R:471@fcbbc478]], `.bartMultinomial` [[R/generics.R:1062@fcbbc478]], `.bartOrdinal` [[R/generics.R:1417@fcbbc478]], `.bartNegbin` [[R/generics.R:1715@fcbbc478]], `.bartHurdle` [[R/generics.R:2292@fcbbc478]], `.rbart` [[R/generics.R:2700@fcbbc478]]; shared `pointwiseLogLikelihood` [[R/generics.R:56@fcbbc478]] |
| 712 | "`$getCalibration`'s docstring (`[[R/dbarts.R:1731@fcbbc478]]`)" | STALE-FALSE | `[[R/dbarts.R:1794@fcbbc478]]` (method at [[R/dbarts.R:1793@fcbbc478]]) |
| 726-729 | "13 `test_*.cpp`, 26,424 lines ... `test_sampler.cpp` 6,841" | STALE-FALSE | `26,430 lines`; `test_sampler.cpp` 6,847 |
| 731-738 | "169 files, 43,098 lines ... `test-capi.R` 1,759 ... `test-predict-blend.R` 390" | STALE-FALSE | `169 files, 44,963 lines`; `test-capi.R` 2,012; `test-predict-blend.R` 417. Others (923 / 572 / 490 / 373 / 287 / 284 / 235) hold. |
| 755-768 | "**Five are `schedule` + `workflow_dispatch` only - `equivalence`, `rchk`, `sbc`, `valgrind`, `revdep-smoke` - and have never run**" | STALE-INCOMPLETE | Still literally true of the five workflow files. Add: since 3f532af2 `exact-gates.yaml` (push-triggered) carries the bcf and multinomial equivalence compares under `--cross-host` ([[R/dbarts.R:141-163@3f532af2]]), so two of `equivalence.yaml`'s three legs now run on every push under a different workflow. Only the gaussian `equivalence.R` leg, `rchk`, `sbc`, `valgrind` and `revdep-smoke` are dormant in substance. |
| 771-776 | "**exact-gates** runs 20 scripts." | STALE-INCOMPLETE | Still 20 exact-posterior scripts, plus a second step running `bcf-equivalence.R --cross-host` and `multinomial-equivalence.R --cross-host` on every push. |
| 805-815 | Stop 10 Known open: "five harnesses ... exactly one recorded (manual) run each - ... `composition-matrix.R`. All baselines are single-host (arm64 macOS); the x86 bench box is unavailable; and the bcf equivalence baseline predates the `summaries` field, so the most expensive dormant gate will degrade loudly on its first live cross-host run." | STALE-FALSE | Replace the last sentence entirely: "All baselines are recorded on arm64 macOS, but cross-host is no longer untested. `bcf-equivalence.R` and `multinomial-equivalence.R` gained a `--cross-host` compare mode (3f532af2) with a two-tier verdict - tier 1 LOCKED is the gate (`rtol = 1e-8`, `atol = rtol * max|a|` on continuous channels, `identical()` on combinatorial ones), tier 2 is a labelled-weak Welch z with an ESS-adjusted denominator, which the design shows cannot gate on its own. The absent `summaries` field turned out to be a pure function of what is stored, so the RC-tip re-record is a refresh, not a prerequisite. Validated on the x86 box before it retired: bcf 12/12 tier 1 PASS (worst ratio 3.7e-06), multinomial 11/11 (2.1e-06), gaussian 43 compared / 0 skipped all |z| = 0.00, tinytest 7478 ok; four discrimination probes fired as designed. `exact-gates.yaml` on ubuntu-latest is now the standing cross-host gate. Reproducibility within a host stays bitwise across SIMD dispatch." Also amend the harness list: `composition-matrix.R` was re-run and amended at 52c10e02 and its two recorded disagreements are closed (0). |
| 839-923 | "Appendix: recent arcs ... Everything landed since 2026-08-20, newest first" ending at NEWS coherence (ae5b91d8) | STALE-INCOMPLETE | Eleven arcs landed since and none appear: capi-shape (9df0cb50), surface-refusals (d48aef8a), rd-records (52c10e02), bcf-cross-host (3f532af2), composition-refusals (fe0b3292, de18ef2b), predict-surface (7b3ac6bf, 71cc7133, befc8f45, ed43deef), dbarts.h freeze (6446ddce), fit-entry refusals (fcbbc478), the monotone branch-fill bench (closed, not adopted), the bartCause `alignForestBasisToSubset` fix, and formal heredity ruled the first post-1.0 arc (3080a9c5, [[TODO:121-125@3080a9c5]]). Four of these now have tracked landing notes: `docs/plans/{capi-shape,surface-refusals,rd-records,bcf-cross-host}.md`. |
| 925-941 | Appendix: gate commands | STALE-INCOMPLETE | Add the two cross-host compares (they are the pre-landing check for anything touching the bcf/multinomial draw paths): `Rscript benchmarks/R/bcf-equivalence.R compare benchmarks/baselines/bcf-equivalence-6e3b9fb8.rds --cross-host` and the multinomial twin against `multinomial-equivalence-4d9a3337.rds`. The gaussian baseline (`equivalence-736bfb05.rds`) and `bench-sampler-ab1dc52.csv` are current. |
| 947-950 | "feature-matrix.md - ... a Gaps section ([[TODO:924@23b9cde7]])" | STALE-FALSE | `([[TODO:929@23b9cde7]])` |
| 955 | "`docs/plans/release-candidate-review.md` (3,831 lines)" | STALE-FALSE | `(3,894 lines)` |
| 960-966 | "docs/plans/review-2026-08-24/ - the second review's artifacts" | STALE-INCOMPLETE | Add `memos/monotone-branch-fill-bench.md` to the memo list; it is where a merge reader lands on the monotone door. |

### 2b. docs/design/feature-matrix.md (all replacements in-line)

| line | claim | disposition | replacement |
|---|---|---|---|
| 37 | "Anchors are verified BY CONTENT against the tree at 52c10e02" | STALE-FALSE (substance) | The 2026-08-26 re-align shifted by the diff line map, not by content; the eight cells below prove it. Either re-derive them and keep the stamp, or restate it. Minimal in-line honest form: `Anchors are line-map aligned to the tree at 52c10e02 - the cited line` |
| 79 | gaussian flat C `S [[CAPI:678@52c10e02]]` | STALE-FALSE | `S [[CAPI:748@52c10e02]]` (the `dbarts_sampler_create` family paragraph; [[CAPI:678@52c10e02]] is inside the `DBARTS_USE_STUBS` macro block, which the dbarts.h freeze inserted ahead of the prototypes) |
| 81 | probit flat C `S [[CAPI:676@52c10e02]]` | STALE-FALSE | `S [[CAPI:747@52c10e02]]` |
| 86 | aft flat C `S [[CAPI:679@52c10e02]]` | STALE-FALSE | `S [[CAPI:749@52c10e02]]` |
| 89 | bcf flat C `S [[CAPI:699-721@52c10e02]]` | STALE-FALSE | `S [[CAPI:770-792@52c10e02]]` - which is exactly what this file's own [f3] ([[CAPI:255@52c10e02]]) already cites for the K-forest paragraph, so the file contradicts itself |
| 161 | ordinal pointwise loglik `S [[generics.R:1425@52c10e02]]` | STALE-FALSE | `S [[generics.R:1437@52c10e02]]` (the `if (type == "loglik")` branch of `extract.bartOrdinal`; `ordinalLogLik` is [[generics.R:1488@52c10e02]], which the Gaps section at [[generics.R:956@52c10e02]] already cites) |
| 162 | nbinom `S [[generics.R:1741@52c10e02]]` | STALE-FALSE | `S [[generics.R:1750@52c10e02]]` (`negbinLogLik` [[generics.R:1778@52c10e02]], per Gaps [[generics.R:963@52c10e02]]) |
| 163 | multinom `S [[generics.R:1096@52c10e02]]` | STALE-FALSE | `S [[generics.R:1105@52c10e02]]` (`multinomialLogLik` [[generics.R:1131@52c10e02]], per Gaps [[generics.R:974@52c10e02]]) |
| 166 | hurdle `S [[generics.R:2300@52c10e02]] [f25]` | STALE-FALSE | `S [[generics.R:2317@52c10e02]] [f25]` (`hurdleLogLik` [[generics.R:2344@52c10e02]], which [f25] itself cites at [[generics.R:620@52c10e02]]) |
| 168 | grouped `S [[generics.R:2802@52c10e02]]` | STALE-FALSE | `S [[generics.R:2812@52c10e02]]` |
| 157-165, 169 | nameable calibration column, nine cells reading `[[dbarts.R:1793@52c10e02]], 1753` | STALE-FALSE | `[[dbarts.R:1793@52c10e02]], 1820`. `[[dbarts.R:1793@52c10e02]]` is `getCalibration`; `setCalibration` is `[[dbarts.R:1820@52c10e02]]`, not `[[dbarts.R:1753@52c10e02]]` (`ptr,`). [f16]'s own later text at [[dbarts.R:497@52c10e02]] already says `[[dbarts.R:1820-1850@52c10e02]]`, and [f23] at [[dbarts.R:604@52c10e02]] says `[[dbarts.R:1847@52c10e02]]`. |
| 175 | gaussian variance forest `S [[FAC:795@52c10e02]]` | STALE-FALSE | `S [[FAC:798@52c10e02]]` (the `if (options.numVarianceTrees > 0 && ...)` gate in `createSampler`) |
| 184 | hurdle DART `S [[bart.R:2346@52c10e02]], 2355 [f35]` | STALE-FALSE | `S [[bart.R:2359@52c10e02]], 2368 [f35]` |
| 185 | bcf variance forest `R [[FAC:874@52c10e02]] [f48]` | STALE-FALSE | `R [[FAC:879@52c10e02]] [f48]` |
| 115, 840 | "`createAmplitudeSampler` [[FAC:868-889@52c10e02]]" (twice) | STALE-FALSE | `[[FAC:870-891@52c10e02]]` |
| 211-215 | "inst/tinytest/test-capi.R drives only the `\"\"`/`\"probit\"` tokens plus grouped ([[FAC:697@52c10e02]]) and heteroscedastic ([[FAC:222@52c10e02]]) by control attribute, and BCF ([[FAC:1356@52c10e02]])" | STALE-FALSE | `grouped ([[FAC:785@52c10e02]]) and heteroscedastic ([[FAC:803@52c10e02]]) ... and BCF ([[FAC:1475@52c10e02]])`. The claim itself (logistic, aft, ordinal, nbinom reachable and untested through dbarts.h) still holds. |
| 342 | [f12] "`bart2Hurdle` ([[bart.R:2304@52c10e02]]) composes ... an occupancy probit ([[bart.R:2346@52c10e02]]) and a lognormal positive part ([[bart.R:2355@52c10e02]])" | STALE-FALSE | `([[bart.R:2317@52c10e02]])`, `([[bart.R:2359@52c10e02]])`, `([[bart.R:2368@52c10e02]])` |
| 355 | [f13] "The flat C API guards through the same call ([[C_interface.cpp:614@52c10e02]], 636)." | STALE-FALSE | `([[C_interface.cpp:658@52c10e02]], 682)` - the two `refuseGroupedScaleUpdate` calls in `dbarts_sampler_setResponse`/`setOffset`; [[C_interface.cpp:614@52c10e02]] and [[C_interface.cpp:636@52c10e02]] are now inside `run`/`setCallback` |
| 364-371 | [f15] "the facade's pure virtual [[FAC:368@52c10e02]] and its shape probe [[FAC:108@52c10e02]]" | STALE-FALSE | `[[FAC:370@52c10e02]]` ([[FAC:108@52c10e02]] holds) |
| 431-439 | [f15] "The flat-C entry, `dbarts_sampler_setActiveRows`, LANDED ... body `[[C_interface.cpp:1174-1187@52c10e02]]` ... (`[[inst/tinytest/test-capi.R:1401-1446@52c10e02]]` pins a genuine mask, an all-ones no-op, a fractional refusal, a NULL clear, and a probit mask moving draws)" | STALE-INCOMPLETE | Anchors hold ([[inst/tinytest/test-capi.R:1174-1189@52c10e02]] now; the test block still starts at [[inst/tinytest/test-capi.R:1401@52c10e02]]). But the semantics changed at 9df0cb50 and the footnote should say so: the fractional case now RAISES through the shared `refuseNonBinaryMask` ([[R_interface_bartcore.cpp:2946@9df0cb50]]) rather than answering 0, so the entry's `int` is a pure capability status under the header's three-class doctrine ([[CAPI:53-69@9df0cb50]]). |
| 464-465 | [f16] "facade [[FAC:354@9df0cb50]], 361" and "R5 ... [[dbarts.R:1793@9df0cb50]], 1753" | STALE-FALSE | `[[FAC:356@9df0cb50]], 363`; `[[dbarts.R:1793@9df0cb50]], 1820` |
| 756, 765 | [f34]/[f35] "[[bart.R:2346@9df0cb50]], 2355" | STALE-FALSE | `[[bart.R:2359@9df0cb50]], 2368` |
| 250-255 | [f3] header/attribute claims ([[CAPI:744-754@9df0cb50]], NOT [[CAPI:510-529@9df0cb50]], K-forest [[CAPI:770-792@9df0cb50]], four attributes undocumented) | OK | All four verified at 3080a9c5; the four selecting attributes are still absent from the header (only `"bartcore.survival"` appears, [[CAPI:753@3080a9c5]]) |
| 723-736 | [f30] "Split verdicts, measured by benchmarks/R/composition-matrix.R" | OK | Verdicts unchanged; the harness itself was amended at 52c10e02 (a `multinomBase` probe and a logistic integer-count `setWeights`) and now reports 0 disagreements. No cell value moves. |
| 787-790 | [f39] baselines | OK | `equivalence-736bfb05.rds` (43), `bcf-equivalence-6e3b9fb8.rds` (12), `multinomial-equivalence-4d9a3337.rds` (11), 11 files + MANIFEST - all confirmed |
| 1017-1027 | Gaps, heteroscedastic: "([[RIB:7184@23b9cde7]], [[RIB:7499@23b9cde7]])" and the `columnMaskStateFeasible` closure | OK | Anchors hold; this is the passage the tour's [[RIB:105-113@23b9cde7]] says is correct, and it still is |
| 929-1037 | Gaps section as a whole | STALE-INCOMPLETE | No gap changes value from today's four slices (they are surface/ABI/Rd shape, not capability). Worth one added line under **Cross-cutting**: the flat C API's capability-status doctrine landed at capi-shape (9df0cb50), so a `dbarts.h` consumer probes a family's capability instead of longjmping - the last flat-C shape item before 1.0-0. |

---

## 3. Reading order

The tour's twelve stops are ordered by blast radius through the *diff*. That is
the right order for someone auditing the change; it is the wrong order for
someone deciding to merge (section 5). Within the existing structure, the four
new landing notes slot as follows:

- **`docs/plans/capi-shape.md`** -> Stop 7, as the stop's primary document,
  ahead of `threaded-predict.md`. It is the newer and larger ABI event, it
  states the return doctrine the whole header now rests on, and its section 11
  is the consumer migration list a merge needs. `threaded-predict.md` section
  11 becomes the second read.
- **`docs/plans/surface-refusals.md`** -> Stop 8, replacing the "Known open"
  no-op paragraph. Read its section 4 (the refusal lists, method by method)
  beside `[[R/generics.R:2109-2166@9df0cb50]]`; its section 17 is the residue list that
  should become Stop 8's Known open verbatim.
- **`docs/plans/rd-records.md`** -> Stop 0, item 4, beside `inst/NEWS.Rd`. It
  is the only place the *documented* claims were checked against the code
  (print/xbart/`$predict` value shapes, `rbart` `$n.chains`), which is what a
  reader of NEWS is implicitly trusting. Its section 7 also belongs to Stop 10,
  as the record that `composition-matrix.R` was re-run.
- **`docs/plans/bcf-cross-host.md`** -> Stop 10, immediately after the CI
  dormancy paragraph, because it is the answer to that paragraph's own worst
  finding. Its sections 2.1-2.2 (why one |z| bar cannot be a gate) are the most
  transferable methodology in the branch and deserve a pointer from the "P17
  rule" sentence.

One reordering inside the existing scheme is worth making regardless: **Stop 7
should precede Stop 6**. The shipped C ABI is the only thing in the branch that
a third party can be broken by, it is frozen at the release, and it is 1,307
lines a reviewer can actually finish. The bridge (Stop 6) is 7,825 lines of
implementation behind it.

---

## 4. Claims I could not verify

1. **tour:100-103** - "`ColumnStore` ([[data.hpp:501@9df0cb50]]) uses uniform 16-bit codes,
   no per-column width ... and it is still not arena-allocated -
   `docs/architecture.md` says both in its own voice." I confirmed
   `ColumnStore` is at [[data.hpp:501@9df0cb50]] and that `hot-layer-u8.md` is a measured
   NO-GO, but did not read the allocation path or architecture.md's two
   sentences closely enough to certify "says both in its own voice".
2. **tour:604-605** - "`plotTree`'s dead padding branch stays out ('no
   judgement names it')." I did not locate the branch; `man/plotTree.Rd` took a
   +4 edit at d48aef8a, so the passage is at least adjacent to a change.
3. **matrix cell VALUES** (as distinct from anchors). The file's own preamble
   says a cell's value is adjudicated separately from its anchor. I verified
   anchors and the specific values named in the brief; I did not re-adjudicate
   the 26 x 13 grid. Two `?` cells (student loglik formula, grouped + variance
   forest) and the bcf pointwise-loglik `M`/`?` question at [[data.hpp:1003-1009@d48aef8a]] are
   unchanged by today's landings but were never adjudicated in the first place.

---

## 5. The tour as a merge document

VD's framing: *"I'm the consumer, and it needs to be worth my time to read
before we merge bartcore into main."* Judged against that, not against its own
stated purpose.

### (a) Length and reading time

7,351 words / 975 lines. Dense technical prose with ~120 file:line citations
reads at 150-200 wpm, so **the tour alone is 40-50 minutes**. But it does not
stand alone: Stop 0 budgets itself "~40 minutes" and directs the reader through
`docs/architecture.md` (2,842 w), `docs/design/feature-matrix.md` (9,743 w),
`inst/NEWS.Rd` (2,551 lines) and `TODO` (329 lines) *before* Stop 1. The entry
chain as written is a **4-6 hour commitment before the first line of diff**, and
the diff it then guides is 828 files / +263,944.

That is not a merge document. It is an audit syllabus.

### (b) What it must front-load, and whether it does

| a merge decision needs | where the tour puts it | verdict |
|---|---|---|
| Public R surface delta vs main | Stop 8, at line 557 of 975; and it is a *file-size* paragraph plus three prose bullets, not an inventory of what a user's script must change | LATE and INCOMPLETE - there is no list anywhere of "what breaks in a script that worked on main" |
| Shipped C ABI delta | Stop 7, line 506; correct in shape, two hash re-bakes stale, and the capi-shape return doctrine is absent entirely | LATE and STALE |
| Consumer migration cost | one sentence at [[data.hpp:526-527@d48aef8a]], threaded-predict only | ABSENT in substance |
| What the merge risks | scattered: Stop 10's dormant-workflow paragraph, eleven "Known open" lists, `TODO` | ABSENT as a single view |
| Doors needing VD's own ruling | scattered across eleven "Known open" lists and never distinguished from residue | ABSENT |
| Gate verdicts, and which are CI-enforced vs manual-only | Stop 10, line 718; genuinely good, and the single most merge-relevant section in the document | LATE - this is the section a merge reader should hit second |
| The one structural idea (leaf model = compile-time `L`) | Stop 0, [[data.hpp:81-85@d48aef8a]] | CORRECTLY front-loaded; 5 lines that pay for themselves |
| Three renames that make greps lie | Stop 0, [[data.hpp:91-98@d48aef8a]] | CORRECTLY front-loaded; prevents a false conclusion |

The tour front-loads *orientation for reading code* and back-loads *evidence
for a decision*. For a merge reader that is exactly inverted.

### (c) What to drop or hand off

- **The recent-arcs appendix ([[data.hpp:839-923@d48aef8a]], 579 w).** Eleven arcs of commit
  archaeology. Four now have tracked landing notes; the rest are in
  `release-candidate-review.md`. A merge reader never needs "which commit
  landed the xbart k-grid sort". CUT, replaced by a pointer.
- **Stop 2's thirteen-row narration (1,449 w, 20% of the document).** Every
  claim in it is a restatement of a `feature-matrix.md` row plus its footnote,
  and the tour says so itself (":the matrix cell and its footnote carry the
  rest", "it is the authority where they disagree"). Keeping both means one of
  them is always the stale one - which is exactly what this audit found. CUT to
  a half-page of what the matrix *cannot* say (the enum's six tokens do not
  enumerate the thirteen shipped rows; the four rows that arrive by attribute,
  ingestion rewrite, R composition or combiner).
- **Stop 11 (89 w) and Stop 4 (172 w).** Fold into one "sediment and the data
  layer" paragraph. Both already say "skim, do not read" / "prefer
  architecture.md's prose".
- **Process narration**: the citation convention ([[data.hpp:15-23@d48aef8a]]), the anchor
  re-derivation history ([[data.hpp:855-860@d48aef8a]]), the "two doc-vs-code divergences" and "one
  document that contradicts itself" paragraphs ([[data.hpp:100-113@d48aef8a]]). These are
  orchestration hygiene, not merge inputs. Hand to `docs/plans/README.md`.
- **Every count a reader would have to re-derive to trust** - file sizes, line
  totals, per-stop diff numbers. This audit found ~15 of them stale. Either
  print the command and not the number, or accept that the document must be
  regenerated at every tip. Prefer the former.
- **The eleven scattered "Known open" lists.** Merge them into one table with
  a column for *who decides*. Residue that needs no ruling goes to `TODO`.

### (d) Proposed structure and length

Target: **~2,400 words / ~320 lines / 25 minutes**, ending in a decision. The
matrix and the landing notes are cited as the places to go deeper, never
inlined.

| section | ~w | what goes in it |
|---|---|---|
| 0. The merge in one page | 200 | What bartcore replaces (the classic engine and `R_C_interface.hpp` are deleted; `inst/include/dbarts/dbarts.h` is the only shipped header), the one structural idea (leaf model is a compile-time `L`, response family a runtime virtual), and the three renames that make greps lie. Verbatim from the current [[data.hpp:81-98@d48aef8a]]. |
| 1. What a user's script sees | 350 | The R surface delta as an inventory, not prose: family vocabulary refused by name, one `offset` spelling, `...` removed from `bart2`/`rbart_vi`, foreign args refused across predict/extract/fitted/residuals, `fitted()`'s slot 3, the own-class generics' new methods. Every row cites `inst/NEWS.Rd`'s UPGRADING block, which is the user-facing contract. Deeper: `surface-refusals.md`, `rd-records.md`. |
| 2. What a linked consumer sees | 300 | `dbarts.h` 477 -> 1,307 lines; the three-class return doctrine; the seven void -> int setters; forest to argument 2 on getTrees/printTrees; `basisRowMajor`; `numThreads` on predict; `USE_FC_LEN_T` gone; MAJOR/MINOR frozen, hash re-baked twice. Then the migration table: stan4bart (N sites), treatSens (has `src/`), bartCause (R-API only, 692/0 after 765a596). Deeper: `capi-shape.md` s11, `threaded-predict.md` s11. |
| 3. Gates, and what each one actually proves | 450 | The current Stop 10, promoted and tightened into a table: gate / trigger / last verdict / date / CI-enforced or manual. Keeps the two headline honest findings verbatim - five workflows have never fired, and every clean rchk/valgrind/SBC claim rests on a manual local run - and adds the cross-host resolution (`bcf-cross-host.md`). This section is why the merge is safe or is not. |
| 4. Risks a merge carries | 300 | The short list of things that could be wrong and would not be caught: PROTECT balance rests on manual runs; no equivalence scenario or SBC arm reaches a latent K-forest; aft and heteroscedastic are uncovered at ensemble scale; warm start and grow-from-root are unrefused and untested at two forests; five harnesses call themselves permanent gates and have one manual run each. One line each, each with its evidence. |
| 5. Doors that need VD's ruling | 250 | One table, filtered from the eleven Known-open lists by "does this need a decision, or is it just unbuilt?": the K-forest batched front door's spelling; `updateScale` refused under every family on a multi-forest sampler against its own plan; real-`r` nbinom; weighted binary (parked); formal heredity as the first post-1.0 arc; the RC declaration itself. Everything else -> `TODO`. |
| 6. If you read code: the walk | 400 | The current Stops 1, 3, 5, 6, 7 compressed to entry anchors and the one thing to judge in each, in the order 7 (ABI) -> 1 (engine) -> 5 (multi-forest) -> 6 (bridge) -> 3/4 (moves, data). Symbol names, not line numbers. |
| 7. Where the rest lives | 150 | `architecture.md` (current state), `feature-matrix.md` (present-state grid, the one deep read), `docs/design/INDEX.md` + `docs/plans/INDEX.md`, the four landing notes, `review-2026-08-24/`, `TODO`. Plus the gate commands. |

The regenerated tour should also carry a one-line regeneration recipe rather
than baked counts, and a stamp of the form "current at <sha>, checked by
`tools/check-doc-freshness.R`" - which today checks anchor identity for
`feature-matrix.md` alone.

### (e) Fate of the existing passages

**Survive verbatim** (5 passages, ~700 words):
- [[data.hpp:81-85@d48aef8a]] the one structural idea -> new s0
- [[data.hpp:91-98@d48aef8a]] three renames that make greps lie -> new s0
- [[data.hpp:755-768@d48aef8a]] the CI dormancy finding, with the cross-host amendment -> new s3
- [[data.hpp:771-776@d48aef8a]] exact-gates, plus the two cross-host compares -> new s3
- [[data.hpp:796-803@d48aef8a]] the calibration lane's 11 families / 75-of-83 verdict -> new s3

**Rewritten** (7):
- Stop 7 -> new s2, restructured around capi-shape's doctrine rather than
  threaded predict, with the migration table added
- Stop 8 -> new s1, from prose into an inventory keyed to NEWS UPGRADING
- Stop 10's harness/baseline paragraphs -> new s3's table
- Stop 10's Known open -> split between new s3 (gate coverage) and s4 (risk)
- the eleven per-stop "Known open" lists -> merged into new s4 + s5
- Stops 1/5/6 -> compressed into new s6, symbols not line numbers
- the "where decisions live" appendix -> new s7, plus the four landing notes

**Cut** (5):
- Stop 2's thirteen-row narration (the matrix is the authority; keep three
  sentences on why thirteen rows exist behind a six-token enum)
- the recent-arcs appendix
- the citation-convention and anchor-history paragraphs
- the "one stamp to distrust" and "one document that contradicts itself"
  paragraphs (both now false, and both orchestration hygiene)
- Stop 11, folded into one sentence of new s6

Net: 7,351 -> ~2,400 words, and the first 850 words carry the whole merge
decision.
