# Composition refusals (D6) and the predict-time NA refusal (D8)

Status: LANDED 2026-08-25 at 936825d7 (design record 93f3155c; code tip de18ef2b; anchors in the doc body are 78f334c1's code).

Spec: docs/plans/prerc-surface-freeze.md D6, D8 and its Sequencing line ("then D1, D9, D5, D2, D6, D8, D7"). TODO
`composition-refusals` ([[TODO:63-65@78f334c1]]) and `predict-na-refusal` ([[TODO:181-183@78f334c1]]) - two commits, one slice; they share only the
gate battery and could land in either order. Evidence: review-2026-08-24/memos/prerc-lens2-backlog.md P5 and N1,
re-anchored live (the memos were written at 7a8c7286). Revised 2026-08-25 after an independent critique and orchestrator
rulings under the standing grant (sections 1, 3, 5, 6, 8, 9, 10, 12, 13). No sampling code moves, no RNG consumption
changes on any path an existing test walks; the equivalence trio is expected bitwise identical and zero baselines
re-record. Budget: ~65 R + ~20 C++ + 5 Rd/NEWS + ~95 test lines.

## 1. D6 census: how the two compositions are spelled

The two names in D6 - grouped + `variance =`, and heteroscedastic + `group.by` - are ONE model reached from its two ends:
a chain whose response model is wrapped by the grouped decorator AND which carries a variance forest. They are separate
cells in docs/design/feature-matrix.md's table 4 ([[TODO:186@7a8c7286]] `grouped` row x `variance forest` column, [[TODO:187@7a8c7286]] `hetero` row x
`grouped ranef` column, both `?` with note [f30]:727-731) because the matrix is read along both axes, and the note says
what the code says: they "still CONSTRUCT ([[CH:641@7a8c7286]] decorates before [[CH:742@7a8c7286]] builds the variance forest)". Both anchors hold
at this tip - [[src/bartcore/chain.hpp:641-646@7a8c7286]] wraps `response_` in `GroupedResponse`, [[chain.hpp:738-742@7a8c7286]] builds the variance
forest afterwards.

Where each half can be declared, taken from formals:

- Variance forest: `dbarts(variance = )` ([[R/dbarts.R:369@7a8c7286]], forwarded [[R/dbarts.R:726@7a8c7286]]), `bart2(variance = )` ([[R/bart.R:687@7a8c7286]], which
  builds its own control and calls `dbarts()`), `dbartsSpec(variance = )` ([[R/spec.R:783@7a8c7286]], forwarded [[R/spec.R:878@7a8c7286]]). All three funnel
  into `resolveSamplerSpec` ([[R/spec.R:102@7a8c7286]]), which resolves the selector at [[R/spec.R:507@7a8c7286]] and attaches `bartcore.variance` to the
  control at [[R/spec.R:532@7a8c7286]]. A LinkingTo consumer can attach that attribute itself; the header states the route at
  [[inst/include/dbarts/dbarts.h:836@7a8c7286]] ("is flat-creatable (its variance forest rides a control attribute)").
- Grouped random intercepts: NO public formal outside `rbart_vi(group.by = )` ([[R/rbart.R:17@7a8c7286]]), which writes
  `bartcore.groups` at [[R/rbart.R:385@7a8c7286]] and is the only shipped writer ([f8], [[feature-matrix.md:307-309@7a8c7286]]). The attribute IS
  the interface: `resolveSamplerSpec` only READS it ([[R/spec.R:449@7a8c7286]], [[R/spec.R:664@7a8c7286]]) and the bridge reads it at creation
  (`applyGroupAttribute`, [[src/R_interface_bartcore.cpp:1905@7a8c7286]], called once at [[src/R_interface_bartcore.cpp:3014@7a8c7286]]). The flat API documents the same route
  in `dbarts_sampler_setResponse` ([[src/C_interface.cpp:630-632@7a8c7286]], "this entry's control carries whatever bartcore.groups
  attribute the consumer put on it, so a flat-API sampler can be grouped") rather than at the create entry.

FIVE entrances construct the composition, not the three an earlier census claimed:

1. `dbarts(x, y, control = <control carrying bartcore.groups>, variance = ...)` - what
   benchmarks/R/composition-matrix.R's two `?` probes build (`extra:variance` on the `grouped` fixture at [[src/C_interface.cpp:664-669@7a8c7286]],
   `extra:groups` on the `hetero` fixture at [[src/C_interface.cpp:670-682@7a8c7286]], both through `do.call(dbarts, args)`). Resolution:
   `resolveSamplerSpec`.
2. `dbartsSpec(data, control = <same>, variance = ...)` ([[R/spec.R:772@7a8c7286]]) - the exported spec surface, same
   `resolveSamplerSpec` call, and the route a LinkingTo consumer is told to use.
3. `new("dbartsSampler", control, model, data)` - the class is exported (NAMESPACE:24 `exportClasses(dbartsSampler)`) and
   its `initialize` ([[R/dbarts.R:905-923@7a8c7286]]) calls `C_dbarts_bartcore_create` on the control it is handed, touching
   `resolveSamplerSpec` at no point. This is not a back door: it is the route `dbarts()` itself takes at [[R/dbarts.R:750@7a8c7286]].
   A caller holding a control with both attributes - one from `dbartsSpec()`, one added by hand - reaches the engine here.
4. `dbarts_sampler_create(control, model, data, family)` ([[src/C_interface.cpp:519-524@7a8c7286]]), which calls
   `bartcore_bridge::createHolder` ([[src/R_interface_bartcore.cpp:2991@7a8c7286]]) directly.
5. Re-creation from a stored control: `getPointer()` ([[R/dbarts.R:1895@7a8c7286]]) and `setState()` ([[R/dbarts.R:1928@7a8c7286]]) both call
   `C_dbarts_bartcore_create` with the R-side `control`, `model` and `data` after a session boundary invalidates the
   pointer. And the stored control can ACQUIRE an attribute: `$setControl`'s carry-forward loop ([[R/dbarts.R:1180-1196@7a8c7286]])
   iterates the names on the OLD control and copies those onto `newControl`, so a bartcore.* attribute the old control
   did NOT have survives on the replacement, is assigned to the field, and is read at the next re-creation. (The bridge's
   `setControl` itself ignores the attributes, [[src/R_interface_bartcore.cpp:1901-1903@7a8c7286]] - which is exactly why the
   acquisition is silent until the sampler is rebuilt.)

Still NOT constructible: through `rbart_vi`, which declares no `variance` formal and no `...` at all (its formals run
[[R/rbart.R:9-52@7a8c7286]]), so R's own wall raises "unused argument" - pinned at [[inst/tinytest/test-heteroscedastic.R:258-262@7a8c7286]];
through `bart2`, which declares no `control` formal; or through `xbart` (neither formal).

## 2. D6: what happens today on each composition

Derived by reading. Both cells CONSTRUCT and RUN, and the composed sampler is SILENTLY WRONG - not merely unadjudicated.
The mechanism, end to end:

- `buildVarianceForest` ([[chain.hpp:4176@7a8c7286]]) ends with `sigmaIsFixed_ = true; sigma_ = 1.0;` ([[chain.hpp:4203-4204@7a8c7286]]): under a variance
  forest sigma is not a parameter, and the residual scale lives in `s^2(x)` row by row - what `setSigma`'s own comment
  states, [[chain.hpp:1705-1712@7a8c7286]].
- Each sweep calls `response_->refreshLatents(rng_, combined, sigma_)` ([[chain.hpp:1535@7a8c7286]], and [[chain.hpp:2029@7a8c7286]] on the callback loop),
  so the grouped decorator is handed sigma = 1.0 forever.
- `GroupedResponse::refreshLatents` ([[model.hpp:4746-4756@7a8c7286]]) passes that sigma, and `base_->workingWeights()` - the USER
  weights, not the heteroscedastic precisions - to `drawGroupEffects` ([[model.hpp:4668-4691@7a8c7286]]), whose conditional is
  `prec_j = (sum_{i in j} w_i) / sigma^2 + 1/tau^2`, `mean_j = (sum_i w_i (z_i - F_i)) / (sigma^2 prec_j)`.
- The mean forest meanwhile backfits against `meanWeights_[i] = w_i / s^2(x_i)` ([[chain.hpp:1409-1416@7a8c7286]], formed at [[chain.hpp:4212@7a8c7286]]),
  and the variance forest backfits against `y - meanFits` where `y` is the grouped working response `z - b`
  ([[chain.hpp:1536@7a8c7286]], [[chain.hpp:1555-1556@7a8c7286]]).

So `b_j` is drawn from the conditional of a model whose residual variance is 1 at every row, while the fitted model's
residual variance is `s^2(x_i)`. Both moments are wrong, and the direction depends on the surface: illustrating at a
constant `s^2`, with `S_j` the group's weighted residual sum and `W_j` its weight sum, the model's own conditional is
mean `S_j / (W_j + s^2/tau^2)` and variance `s^2 / (W_j + s^2/tau^2)`, while the draw uses mean `S_j / (W_j + 1/tau^2)`
and variance `1 / (W_j + 1/tau^2)` - over-shrunk and too WIDE where `s^2 < 1`, under-shrunk and too NARROW where
`s^2 > 1`. At a genuinely varying `s^2(x)` the weights inside both sums are wrong as well. The variance-forest half of the
composition composes; the group block does not. Separately, `tau`'s prior scale is converted once at construction from
`rel.scale = sd(y)` ([[R/rbart.R:389@7a8c7286]], [[model.hpp:4712-4716@7a8c7286]]), a homoscedastic calibration a variance surface leaves undefined.

Nothing errors, nothing warns, and every reported channel has the right shape - which is why "acceptance IS the contract"
(P5) and why the refusal is owed before the release rather than an adjudication after it.

Confirming run, which this design does not depend on: `Rscript benchmarks/R/composition-matrix.R` prints the two `?`
cells' verdicts under "? resolutions" (the runner at [[model.hpp:706-733@7a8c7286]]); expect `grouped varianceForest -> CONSTRUCTS` and
`hetero groupedRanef -> CONSTRUCTS` before the change.

## 3. D6: the refusal sites, exactly

TWO sites, the two-layer pattern the sibling grouped-ordinal and grouped-nbinom refusals already use (R surface at
`rbart_vi`'s family token list, bridge backstop at [[src/R_interface_bartcore.cpp:3016-3025@7a8c7286]]). Five entrances, two sites: the
R check is the friendly early error on the documented path, and the C++ check is the CLOSING COVER for all five, because
every one of them ends in `createHolder`.

(a) R, inside `resolveSamplerSpec`'s variance block ([[R/spec.R:508-527@7a8c7286]]), after the Student-t refusal at [[R/spec.R:519-524@7a8c7286]] and
before the monotone one at [[R/spec.R:525-527@7a8c7286]] - the block's order is family, then residual law, then leaf constraint, and this is a
response-model composition, so it sits with the residual law. It reads the control the caller supplied, so it fires for
`dbarts()`, `bart2()` and `dbartsSpec()` (entrances 1 and 2):

    # rbart_vi's group block draws b_j at the SCALAR residual scale the chain
    # carries, which a variance forest pins at 1 while s^2(x) holds the residual
    # scale row by row; the group effects would condition on a residual variance
    # the fitted model does not have. Refuse rather than fit a composition whose
    # Gibbs blocks disagree.
    if (!is.null(attr(control, "bartcore.groups"))) {
      stop(
        "a variance forest does not support grouped random effects: the group ",
        "effects draw at a scalar residual scale, which the variance forest ",
        "replaces row by row"
      )
    }

The message names both parts, in the wording of its two siblings ("a variance forest does not support Student-t
residuals: the two are not yet shown to compose"), and it is reached from either spelling because there is only one
construction. Placement after the `family != "gaussian"` check at [[R/spec.R:509-515@7a8c7286]] means a grouped probit + `variance` still gets
the family refusal first, which is the more specific fact.

(b) C++, in `createHolder`'s creation lambda, immediately after `applyVarianceAttributes`
([[src/R_interface_bartcore.cpp:3031-3032@7a8c7286]]), in the shape of the two group refusals ten lines above:

    // grouped random effects and a variance forest is an unadjudicated
    // composition: the group block draws b at the scalar sigma a variance
    // forest pins at 1, so the effects condition on a residual variance the
    // model does not have. The backstop the R surface (spec.R) mirrors.
    if (options.numGroups > 0 && options.numVarianceTrees > 0)
      Rf_error("grouped random effects are not supported with a "
               "heteroscedastic variance forest");

Five lines, and they are the reason section 13's ruling 2 is (a) rather than optional: `createHolder` is the ONE creation
function that `new("dbartsSampler", ...)`, `dbarts_sampler_create`, `getPointer()` and `setState()` all reach
([[R/dbarts.R:922@7a8c7286]], [[R/dbarts.R:1895@7a8c7286]], [[R/dbarts.R:1928@7a8c7286]] and [[src/C_interface.cpp:521@7a8c7286]] all call `C_dbarts_bartcore_create` / `createHolder`), so the
three entrances the R check cannot see are closed here, including the acquire-then-rebuild path of entrance 5. No third
site in `$initialize`: it would be redundant with this one and would state the same rule twice.

Formals untouched, all of them: `dbarts(variance = )`, `bart2(variance = )`, `dbartsSpec(variance = )`,
`rbart_vi(group.by = , group.by.test = )`, `varianceForest()`'s own arguments, and the `bartcore.groups` /
`bartcore.variance` attribute contracts. Nothing is withdrawn from any surface; two compositions become validation errors.

## 4. D6 door memo (drafted here; its home at landing)

There is no doors ledger file. The practiced location is the owning docs/design doc's own post-landing section - e.g.
docs/design/heteroscedastic.md "15. Post-landing: Student-t residuals and grouped random effects refused (2026-08-17)"
([[src/C_interface.cpp:680-697@7a8c7286]]), which is where the sibling refusal's memo lives. So this memo lands as a NEW section 17 of
docs/design/heteroscedastic.md, and section 15 is CORRECTED in the same commit: its claim that "a grouped (rbart_vi) fit
alongside `variance` was already unreachable" is true of the rbart_vi spelling only and false of the four other
entrances, and its parenthetical reason ("its `...` is rejection-only") is wrong twice over - `rbart_vi` has no `...` at
all ([[R/rbart.R:9-52@7a8c7286]]), which is why R's own unused-argument wall is what raises. feature-matrix [f30] has recorded the
live behaviour correctly all along. docs/design/grouped-random-effects.md gets one cross-reference by section title, not
a restatement (docs/plans/README.md "Cross-references").

Drafted text:

> ## 17. Post-landing: grouped random effects refused with a variance forest (D6)
>
> Section 15 read the composition off `rbart_vi`'s formals alone and concluded it was unreachable. It is reachable
> wherever the two halves can be declared on one control, which is four entrances beside `rbart_vi` - `dbarts()`,
> `dbartsSpec()`, `new("dbartsSampler", ...)` and `dbarts_sampler_create`, plus re-creation from a control that acquired
> the attribute through `$setControl` - and it constructed and ran until D6. It is now a validation error at
> `resolveSamplerSpec`'s variance block with a backstop at `createHolder`, the one creation function all of them reach;
> both formals stay.
>
> What an adjudication would need, all three:
>
> - Engine. `drawGroupEffects` conditions on a scalar sigma and the base model's weights; under a variance forest the
>   per-row precision is `w_i / s^2(x_i)`, which the chain forms as `meanWeights_` and never shows the response model. A
>   correct b block needs that vector at the response-model interface - a change to `ResponseModel::refreshLatents`'s
>   contract, which every family implements, and the reason this is more than dropping a refusal.
> - Prior interaction and identifiability. `tau` and `s^2(x)` both explain dispersion, and they are not separated by the
>   data when the grouping is x-measurable (a group indicator among the predictors, a group that IS a region of x): the
>   variance forest can absorb between-group spread as a wider `s^2` on that region. The tau prior is also calibrated
>   once at construction from `rel.scale = sd(y)`, a homoscedastic quantity; under a variance surface it must be restated
>   against something the model still names.
> - Test oracle. The decisive cheap one is a conditional-exactness gate in the shape of
>   `benchmarks/R/heteroscedastic-exact.R` and `benchmarks/R/grouped-mixing.R`: hold `s^2(x)` fixed and check the b draw
>   against a brute-force posterior. It falsifies the pre-D6 code immediately and verifies a fix. A joint SBC arm is the
>   stronger form and needs the heteroscedastic SBC arm first, which is itself deferred ([f47]).
>
> A surface decision rides with it for the `group.by` spelling: the composition has no public spelling at all, so support
> means adding `variance` to `rbart_vi`'s formals (section 15's own suggestion) or `group.by` to `dbarts()`/`bart2()`
> (deferred post-1.0, VD 2026-08-24). Either is additive.

Records edited in the same commit, so no doc reads as if the cells were still open: [[feature-matrix.md:186@7a8c7286]] and [[feature-matrix.md:187@7a8c7286]] become
`R spec.R:NNN` (the two `?` cells; NNN is the new stop's line), [f30]:727-731 loses its "still CONSTRUCT" sentence and
gains one naming the refusal and pointing at heteroscedastic.md section 17 by title, and [[feature-matrix.md:1010-1011@7a8c7286]]'s grouped narrative
("Composition with a variance forest constructs unrefused and untested ([f30])") is rewritten. benchmarks/R/
composition-matrix.R needs no edit: it probes only `S` and `?` cells ([[feature-matrix.md:706@7a8c7286]], [[feature-matrix.md:715-731@7a8c7286]]), so the two drop out of the probe
set when the matrix says `R`, and the run must then report zero disagreements.

## 5. D8 census: newdata validation, and the training-NA record

Every predict route, and where its `newdata` is validated:

| method | R/generics.R | newdata path |
|---|---|---|
| `predict.bart` | [[R/generics.R:296@ed43deef]] | `object$fit$predict(newdata, offset, n.threads)` [[R/generics.R:386@ed43deef]] (and `predictForest` [[R/generics.R:668@ed43deef]], `predictBlend` [[R/generics.R:808@ed43deef]], through the same R5 methods) |
| `predict.bartMultinomial` | [[R/generics.R:1216@ed43deef]] | `validateXTest(newdata, object$fit$data@x)` [[R/generics.R:1240@ed43deef]], then `$predict` [[R/generics.R:1267@ed43deef]] (which validates again) |
| `predict.bartOrdinal` | [[R/generics.R:1502@ed43deef]] | `validateXTest` [[R/generics.R:1524@ed43deef]], then `bartcorePredict` [[R/generics.R:1534@ed43deef]] |
| `predict.bartNegbin` | [[R/generics.R:1754@ed43deef]] | `validateXTest` [[R/generics.R:1777@ed43deef]], then `bartcorePredict` [[R/generics.R:1786@ed43deef]] |
| `predict.bartHurdle` | [[R/generics.R:2147@ed43deef]] | `hurdleParts` [[R/generics.R:1990@ed43deef]] -> two `predict.bart` calls |
| `predict.rbart` | [[R/generics.R:2194@ed43deef]] | `object$fit[[i]]$predict(newdata, offset, n.threads)` [[R/generics.R:2249@ed43deef]], [[R/generics.R:2259@ed43deef]], [[R/generics.R:2269@ed43deef]] |

Plus `survivalProbabilities.bart` / `.rbart` ([[R/bart.R:2500@7a8c7286]], [[R/bart.R:2551@7a8c7286]]), which build an expanded newdata and forward to
`predict`, and the R5 surface itself: `$predict` ([[R/dbarts.R:1084@7a8c7286]]) and `$predictForests` ([[R/dbarts.R:1147@7a8c7286]]) both open with
`validateXTest(x.test, data@x)` ([[R/dbarts.R:1088@7a8c7286]], [[R/dbarts.R:1155@7a8c7286]]), as do `$setTestPredictorAndOffset` ([[R/dbarts.R:1586@7a8c7286]]) and `$getTrees(newdata)`
([[R/dbarts.R:2057@7a8c7286]]). `bartcorePredict` ([[R/bartcore.R:1085@7a8c7286]]) does NOT validate - it coerces to a double matrix and calls the bridge -
which is why the two families that use it validate in the generic instead.

So `validateXTest` ([[R/data.R:50-283@7a8c7286]]) is the single R-side chokepoint for out-of-sample predictors: every S3 predict,
every R5 predict/test-surface method, `dbartsData(test = )` ([[R/data.R:1492@7a8c7286]]), and the EXPORTED `makeTestModelMatrix`
([[R/data.R:46@7a8c7286]], `export(makeind, makeModelMatrixFromDataFrame, makeTestModelMatrix)` at NAMESPACE:16, documented
[[man/makeind.Rd:12@7a8c7286]]) all reach it. It returns with the test columns already reordered and renamed onto the training columns
([[R/data.R:265-279@7a8c7286]]), so column j of its value is column j of `x.train`.

Three routes do not pass through it, and the design states them as the contract's edge rather than chasing them: the
per-column `$setTestPredictor(x, column)` form ([[R/bartcore.R:571-600@7a8c7286]] takes the `column` branch), the internal store-view
build (section 8's residue), and DIRECT SLOT ASSIGNMENT - `data@x.test <- ...` on a `dbartsData` object, which bartCause
does at fit time (bartCause's `R/responseFit.R` lines 180-189). The contract this slice ships is therefore: the validated surfaces
plus the flat backstop; a raw slot poke is outside it, as it is for every other check `validateXTest` performs (column
count, category codes, factor levels).

The training-NA record EXISTS, on both sides, and neither is serialized:

- C++: `ColumnStore::hasMissing`, one `uint8_t` per column ([[src/bartcore/data.hpp:535@7a8c7286]], its invariant comment at [[src/bartcore/data.hpp:533-534@7a8c7286]]),
  assigned by every full quantize (`quantizeDenseInto` [[src/bartcore/data.hpp:954@7a8c7286]], `quantizeCscColumnInto` [[src/bartcore/data.hpp:999@7a8c7286]], both
  `hasMissingOut[j] = anyMissing`) and by a view build ([[src/bartcore/data.hpp:1465@7a8c7286]]). This is the record the ENGINE routes by, and it gates the
  extra missing-direction draw, so "a column flagged 0 consumes no missing-direction draw from the rng" is what makes the
  whole feature free on NA-free data.
- R: `data@x` keeps its NAs (pinned, [[inst/tinytest/test-data-missing.R:20@7a8c7286]]), and `validateXTest` already holds it as
  `x.train`.

Saved-fit consequence: NONE. The record is derived at build time from the stored data object, never written to state, so
`stateFormatVersion` stays 3 and `minReadableStateFormatVersion` stays 3 ([[src/R_interface_bartcore.cpp:6390@7a8c7286]], [[src/R_interface_bartcore.cpp:6399@7a8c7286]]); a
fit saved before this slice re-creates its sampler through `getPointer()` from the same `data@x` and gets the same
`hasMissing` it always had. Nothing to migrate, no version floor to move. This is the main reason NOT to add a new R-side
slot recording the NA columns: a new field would be absent on every pre-slice saved fit and would need a fallback that
reads `data@x` anyway.

One known divergence between the two records, stated because the R-side check inherits it: `setCell` marks
`hasMissing[j]` when a value becomes NA and never clears it ([[data.hpp:1755-1764@7a8c7286]], "the flag only clears on a full column
re-quantize"), so after a per-cell mutation removes a column's LAST NA the engine still calls the column missing-bearing
while `data@x` says it is complete, and the R check then refuses a newdata NA the engine might still route. The corner is
reachable: `updatePredictorPerObservationJointly` writes the mutated column back to `data@x` ([[R/updatePredictorPerObservationJointly.R:89@7a8c7286]])
while the engine's commit sets the flag through `setCell` ([[src/bartcore/sampler.hpp:1938@7a8c7286]]). Recorded as residue; a
whole-column `setPredictor` re-quantizes and the two records re-agree.

## 6. D8: what happens today

Two routing paths reach a test NA, and both end LEFT for the same reason.

- PREDICT (both `predictFromSavedSample`, [[chain.hpp:2832@7a8c7286]], and `predictFromCurrentTrees`, [[chain.hpp:2857@7a8c7286]], which flattens the live
  trees on the fly) replays through `addFlatPredictions` -> `partitionFlatIndices` ([[tree.hpp:1727@7a8c7286]]), which reads RAW
  values and branches on `isNA(value)` against `missingGoesLeft = (flat.flags & flatMissingGoesRight) == 0` ([[tree.hpp:1732@7a8c7286]]). The
  flag is set when a tree is flattened only if the rule carries the bit ([[tree.hpp:1465@7a8c7286]]), and `buildFromFlat` REFUSES to restore
  it on an NA-free column ([[tree.hpp:1558-1560@7a8c7286]]). `codeFor` is not on this path at all.
- The RUN's own test-fit channel (a test set installed on the sampler by `dbartsData(test = )` or
  `$setTestPredictor`) quantizes instead: `ColumnStore::codeFor` ([[data.hpp:695-706@7a8c7286]]) codes any NA to `naCode` for an
  ordinal column and to `missingCategoryCode(numCuts[j])` for a categorical one, and `Rule::sendsRight` ([[tree.hpp:155-162@7a8c7286]])
  returns `missingGoesRight()` for `naCode` while the categorical mask never sets the missing bit on an NA-free column
  ([[tree.hpp:479@7a8c7286]]).

Either way the bit is 0, because it is drawn only where `data.hasMissing[variableIndex]` ([[tree.hpp:450@7a8c7286]], [[tree.hpp:609@7a8c7286]], [[tree.hpp:875@7a8c7286]],
:889): on a training-complete column the row goes LEFT at every split on that column, in every tree, in every draw.
`validateColumnValues` ([[src/R_interface_bartcore.cpp:2911-2917@7a8c7286]]) passes the value through - `categoricalValueIsValid`
returns true for NA by construction ([[data.hpp:707-712@7a8c7286]]). So both surfaces return numbers, silently, from a route the model
never learned. [[man/dbarts.Rd:117@7a8c7286]] promises the opposite ("in training and test data alike") and [[inst/NEWS.Rd:720-725@7a8c7286]]
repeats the promise. A fit built with `missing = "error"` is the sharpest case: it has no NA-bearing column at all, so
today it answers every NA test row from the left branch of a model that was explicitly asked to reject NAs.

## 7. D8: what is refused, and what stays accepted

Refused: an NA in a test/newdata column whose TRAINING values carried none. Accepted, unchanged: an NA in a column that
carried training NAs (the learned route - the case [[test-data-missing.R:43-57@7a8c7286]] and [[test-data-missing.R:97-100@7a8c7286]] pin), and every complete
newdata. The rule is per COLUMN, not per fit, so a partially missing design keeps working exactly where it worked.

Family by family there is nothing to say, and the design says so once: every family shares one `ColumnStore` and one
routing rule, so the fact and the refusal are family-independent - gaussian, probit, logistic, aft, hazard, ordinal,
nbinom, multinomial (one store, K forests) and the grouped and heteroscedastic decorations all inherit it. Two
compositions have a visible consequence:

- Hurdle composes TWO samplers over different row sets ([[R/bart.R:1130@7a8c7286]], "occupancy probit on 1{y > 0} (all n) and a
  gaussian"; `hurdleParts`, [[R/generics.R:1990@7a8c7286]], forwards `newdata` to each). The positive part trains on `{i : y_i > 0}`
  only, so a column whose training NAs all sit in the zero rows is complete FOR THAT PART: the same newdata passes the
  occupancy predict and is refused by the positive one. That is correct - the positive part learned no route - and the
  message names the column, so the caller can see which half is speaking.
- rbart's per-chain samplers ([[R/generics.R:2249-2269@7a8c7286]]) share one data object, so all chains agree; the refusal fires on
  the first.

Other predict arguments are unaffected and stay complete-only where they already were: `offset` / `offset.test` are
refused with NAs at ingestion ([[R/data.R:1713-1716@7a8c7286]]) and a length-1 `NA` on the R5 flat path still reads as "no offset"
([[R/dbarts.R:1136-1141@7a8c7286]]); `weights` feeds a posterior-predictive draw and never a route; `bases` is a caller-supplied basis
for the amplitude blend, not a predictor, and goes through `validateForestBases`. None of them acquires a missingness
rule here.

## 8. D8: where the refusal fires

R primary, C++ backstop, on routes that do not overlap except where the R surface itself validates twice.

(a) R: one helper called at the tail of `validateXTest` (R/data.R, immediately before its `x.test` return at [[R/dbarts.R:283@7a8c7286]]). That
is one site for every route in section 5's table, including the R5 `$predict`, the fit-time `dbartsData(test = )` and the
exported `makeTestModelMatrix`, and it is the only place with column NAMES, which is what "refused by name" is worth.
Beside it, a per-column probe that reads either storage without densifying, on `predictorSourceColumn`
([[R/mixedMatrix.R:413-441@7a8c7286]]), which already answers for a plain matrix, a `dbartsMixedMatrix` and a bare `dgCMatrix`:

    ## A split rule learns a route for NA only on a column whose TRAINING values
    ## carried one: the missing direction is drawn only there and cannot be
    ## restored onto an NA-free column, so on a training-complete column every
    ## rule sends NA down one fixed branch. Refuse rather than answer from a
    ## route the model never learned.
    sourceColumnHasNA <- function(source, j, numColumns, numObservations) {
      column <- predictorSourceColumn(source, j, numColumns, numObservations)
      if (is.list(column)) {
        return(anyNA(column$x) || is.na(column$implicit))
      }
      anyNA(column)
    }

    sourceHasNA <- function(source) {
      if (inherits(source, "dbartsMixedMatrix")) {
        return(
          any(vapply(source$dense, anyNA, logical(1L))) ||
            (!is.null(source$sparse) && anyNA(source$sparse@x))
        )
      }
      if (inherits(source, "dgCMatrix")) {
        return(anyNA(source@x))
      }
      anyNA(source)
    }

    refuseTestMissingness <- function(x.test, x.train) {
      # the whole-object probe short-circuits, so complete test data - the usual
      # case - pays one scan and never touches the training side
      if (!sourceHasNA(x.test)) {
        return(invisible(NULL))
      }
      # the two sides carry the same p by now ([[R/mixedMatrix.R:220@7a8c7286]] refuses otherwise), so one
      # column count serves both
      numColumns <- NCOL(x.test)
      numTest <- NROW(x.test)
      numTrain <- NROW(x.train)
      offending <- integer(0L)
      for (j in seq_len(numColumns)) {
        if (!sourceColumnHasNA(x.test, j, numColumns, numTest)) next
        if (sourceColumnHasNA(x.train, j, numColumns, numTrain)) next
        offending <- c(offending, j)
      }
      # every NA sits in a column that carried training NAs: every one has a
      # learned route, and nothing is refused
      if (length(offending) == 0L) {
        return(invisible(NULL))
      }
      predictorNames <- colnames(x.train)
      labels <- if (is.null(predictorNames)) {
        paste0("column ", offending)
      } else {
        paste0("'", predictorNames[offending], "'")
      }
      shown <- head(labels, 5L)
      stop(
        "test predictors have missing values in ",
        toString(shown),
        if (length(labels) > 5L) {
          paste0(" and ", length(labels) - 5L, " more column(s)")
        },
        ", which carried none in training: a split rule learns a route for NA ",
        "only on a column that had missing values when the trees were grown, ",
        "so these rows have no route to take"
      )
    }

Every offending column is reported, not just the first, so a caller fixing an imputation gap sees the whole gap in one
message; the list is bounded at five plus a count so a wide design cannot print a wall. Column positions are 1-based, as
R indexes and as `colnames(x.train)` is subscripted. Cost: complete newdata pays one short-circuiting `anyNA` over the
test object and never touches the training side; NA-bearing newdata pays one training-column scan per offending test
column. `predict.bartMultinomial` validates twice - once in the generic ([[R/generics.R:1240@7a8c7286]]) and once inside `$predict`
([[R/dbarts.R:1088@7a8c7286]]) - so its scan runs twice; that is a pre-existing double coding this slice inherits rather than
introduces, and an `anyNA` short-circuit is cheap enough to leave alone rather than restructure a predict path for.

`validateXTest`'s own messages say `'test'` at every existing site ([[R/dbarts.R:190@7a8c7286]], [[R/dbarts.R:220@7a8c7286]], [[R/dbarts.R:258@7a8c7286]]), including when the caller is
`predict`, so the new one keeps that vocabulary rather than growing an `argument` parameter and eight call-site edits.

(b) C++: `validateTestSource` ([[src/C_interface.cpp:253-257@7a8c7286]]) is called at exactly two sites, [[src/C_interface.cpp:788@7a8c7286]]
(`dbarts_sampler_setTestPredictors`) and [[src/C_interface.cpp:831@7a8c7286]] (`dbarts_sampler_predict`), and by NO R-bridge entry - the R bridge
validates through `validatePredictorMatrix` ([[src/R_interface_bartcore.cpp:130-142@7a8c7286]]) and
`validateTestContainerAgainstStore` ([[src/R_interface_bartcore.cpp:2970@7a8c7286]]) instead. So a check placed inside it covers both flat test entrances and
never runs for an R caller. About ten lines, using the reader the replay itself builds (`PredictorSourceColumns`,
[[src/bartcore/data.hpp:401-440@7a8c7286]]), naming the column by INDEX since the flat caller supplies no names:

    // A test NA takes a rule's learned missing direction, and a rule learns one
    // only where the training column had NAs (ColumnStore::hasMissing gates the
    // draw), so on a complete column it would take one fixed branch at every
    // split. The R surface refuses first and names the column; this is the flat
    // entrances' backstop, which has only the index.
    void refuseTestMissingness(const bartcore::ColumnStore& store,
                               const bartcore::PredictorSource& source,
                               const char* caller) {
      bartcore::PredictorSourceColumns columns(source, store.types.data());
      for (size_t j = 0; j < source.numColumns; ++j) {
        if (store.hasMissing[j]) continue;
        bartcore::PredictorSourceColumnReader column = columns.column(j);
        for (size_t i = 0; i < source.numRows; ++i)
          if (bartcore::isNA(column.at(i)))
            Rf_error("%s: test column %zu has missing values but the training "
                     "column had none, so no rule routes them", caller, j + 1);
      }
    }

Placed beside `validateTestSource` in src/C_interface.cpp (not promoted to the bridge block: no R-side caller wants it),
called from it after `validateTestContainerAgainstStore`, so `caller` is the flat entry's own label already passed to
`translateSource`. Alternatives rejected: putting the check in `validateColumnValues` (nine call sites, most of them
TRAINING mutations - it would refuse `setPredictor` on an NA-free column, a different and recorded case,
docs/design/mia-missingness.md "Mutation surface"); putting it only in C++ (loses the column name, and misses
`dbartsData(test = )`, which never reaches a predict entry); adding a bridge reader for `hasMissing` so R could ask the
engine (a new entry, and R already holds `data@x`).

Residue, recorded not fixed, the same fact in places the validated surfaces cannot see:

- The internal store view (`ColumnStore::buildFromParent`, [[data.hpp:1392@7a8c7286]], built at [[src/R_interface_bartcore.cpp:3678@7a8c7286]]) is
  how xbart's folds are cut. A view recomputes `hasMissing` from ITS gathered training rows ([[src/R_interface_bartcore.cpp:1459-1465@7a8c7286]]) while the test
  side "gates no draws, so it tracks no missingness" ([[src/R_interface_bartcore.cpp:1050@7a8c7286]]), so a held-out row can be NA on a column complete within its
  fold and route left. Refusing it would mean a data-dependent per-fold error at ingestion; out of scope.
- `$setTestPredictor(x, column)`'s per-column form ([[R/bartcore.R:571-600@7a8c7286]]) and direct `data@x.test <-` assignment
  (section 5) are outside the validated surface by the same contract that governs every other `validateXTest` check.
- Mid-run mutation that ADDS NAs to a training-complete column keeps its documented behavior (left until moves revisit
  the rules) - unchanged by this slice.
- The sticky `hasMissing` after per-cell mutation (section 5), reachable through
  `updatePredictorPerObservationJointly`, can make the R check conservative by one corner.

## 9. Rd and NEWS

D6: no Rd changes. No formal moves, and neither `man/dbarts.Rd`'s `variance` item nor `man/rbart.Rd`'s `group.by` item
claims the composition works; adding a sentence to `\item{variance}` naming a composition that has no public spelling
would document a refusal a reader cannot trigger from the documented surface. The record lives in the design doc
(section 4).

D8: the rule goes where the promise it corrects lives.

- [[man/dbarts.Rd:116-118@7a8c7286]] `\item{missing}`, the one D8 names. Replace the clause "in training and test data alike" and add
  the rule: "... follows that rule's chosen branch (\dQuote{Missingness Incorporated in Attributes}, Twala et al. 2008).
  A route is learned per COLUMN and only where that column's TRAINING values were missing, so test predictors and
  \code{newdata} may carry \code{NA} in those columns and are refused, naming the column, in any column that was complete
  in training - there is no route for such a value to take. \code{"error"} rejects predictors containing \code{NA} ..."
- [[man/bart2.Rd:311-313@7a8c7286]] `\item{missing}` says "exactly as in \code{\link{dbarts}}" and needs one added sentence: the
  refusal on a training-complete column, since bart2 is where most users meet it.
- man/makeind.Rd, whose `makeTestModelMatrix(data, newdata)` entry ([[man/bart2.Rd:12@7a8c7286]]) is an exported caller of `validateXTest` and so
  acquires the refusal: one sentence in its description saying that `newdata` is refused, naming the column, where it
  carries \code{NA} in a column the fit's training data did not.
- [[man/xbart.Rd:99-101@7a8c7286]] defers wholly to `?dbarts` and needs no edit.
- [[inst/NEWS.Rd:720-725@7a8c7286]], the 1.0-0 NEW FEATURES entry that repeats "training and test data alike": same correction. It is
  unreleased text, so this is a fix, not a rewrite of history.
- Two new `\item`s in `\subsection{UPGRADING}` ([[inst/NEWS.Rd:5@7a8c7286]]), one per commit: D6 (the composition is a validation
  error, both formals stay) and D8 (NA in test predictors on a training-complete column is now an error). The file holds
  339 `\item` entries at this tip and should hold 341 after; that is a landing-note observation, not a gate.

## 10. Test plan

House battery per commit: `R CMD INSTALL .` (`--preclean` on commit 1, which touches no header but does touch
src/R_interface_bartcore.cpp; it costs nothing), `tinytest::test_package("dbarts")`, `tools/check-rc-codoc.R`,
`Rscript tools/check-doc-freshness.R .`, `air format --check .`, `lintr::lint_package()`, `R CMD check --as-cran`.
Neither commit changes an S4 generator's `methods = list(...)` or an S3 signature, so codoc has nothing new to compare;
both run anyway.

Existing tests, verified by reading:

- D6: nothing breaks. No in-repo test constructs both halves - the four files that mention `bartcore.groups` and
  `variance` keep them in separate cells ([[test-fits-without-offset.R:179-207@7a8c7286]], [[test-calibration-prior-draws.R:219-271@7a8c7286]],
  test-bcf-creation.R, test-capi.R) - and [[test-heteroscedastic.R:222-262@7a8c7286]] already pins the rbart_vi wall.
- D8: ONE test must be rewritten. [[inst/tinytest/test-sparse-factor.R:798-806@7a8c7286]], under the "NO FALSE REFUSALS" banner at
  [[inst/tinytest/test-sparse-factor.R:783@7a8c7286]], builds `test.dense.na <- cbind(rnorm(3L), c(0, NA, 2))` over columns `x1`/`f` and asserts that
  `dbarts(dbartsData(train.bound, y.bound, test = test.dense.na), control = boundControl)` constructs. `train.bound`
  ([[inst/tinytest/test-sparse-factor.R:616@7a8c7286]], `sparseFrame(labels.bound, levels.small, "a")`) carries no NA in `f`, so the new rule refuses it. Its intent -
  an NA code must not be read as an out-of-range CATEGORY - survives if the training column carries an NA too, so the
  fixture gains one: build a local NA-bearing training design and its own sampler for this assertion rather than
  mutating `train.bound`, which `sampler.bound` and every code-bound assertion above it share, and place the new block
  at the END of the file, since the file is seeded at [[inst/tinytest/test-sparse-factor.R:610@7a8c7286]] and a new `rnorm` call would otherwise shift every draw after
  it. The assertion itself keeps its shape (`expect_inherits(..., "dbartsSampler")`).
  Everything else holds: the missingness suite's test set carries NAs only on columns that carried them in training
  ([[test-data-missing.R:9-11@7a8c7286]] puts NAs in `x1` and `g`; `test.df <- df[1:20, c("x1","x2","g")]` at [[test-data-missing.R:45@7a8c7286]] inherits exactly
  those), so `dbarts(..., test = test.df)` at [[test-data-missing.R:52@7a8c7286]] and `sampler.keep$predict(test.df)` at [[test-data-missing.R:97@7a8c7286]] still pass. The implementer
  must still run the full suite: this is a refusal on data, and the census is a reading.

New tests.

D6, extending inst/tinytest/test-heteroscedastic.R's existing section "Student-t residuals and a grouped fit:
unadjudicated with 'variance'" ([[inst/tinytest/test-heteroscedastic.R:222@7a8c7286]]) rather than opening a file:
- a control carrying `bartcore.groups` (the four sibling files' idiom) plus `variance = varianceForest(n.trees = 3L)`
  through `dbarts()` raises, matching `"does not support grouped random effects"`;
- the same two halves through `dbartsSpec()` raise the same message;
- the same two halves through `new("dbartsSampler", control, model, data)`, built from a `dbartsSpec()` result with the
  group attribute added afterwards, raise the BRIDGE message - the entrance the R check does not see;
- each half ALONE still constructs (the assertion that this is a composition check, not a regression on either feature);
- the refusal fires whichever order the halves are declared in.
D6, in inst/tinytest/test-capi.R beside the existing `capi_create` refusal probes ([[test-data-missing.R:118-150@7a8c7286]]), reusing the grouped control
the file already builds at [[inst/tinytest/test-capi.R:776@fe0b3292]]: that control plus a `bartcore.variance` attribute through `capi_create` raises
`"not supported with a heteroscedastic variance forest"` - the flat entrance's only coverage.

D8, in inst/tinytest/test-data-missing.R, appended at the END of the file (several blocks there hardcode values that
depend on the file's full execution history, so a new block goes last or carries its own `set.seed`):
- `sampler.keep$predict(test.df)` with an NA written into `x2` (training-complete) raises, matching `"'x2'"` and
  `"carried none in training"`;
- the same NA in `x1` (training-incomplete) still predicts, shape unchanged - the acceptance half;
- NAs in BOTH `x2` and a second complete column: the message names both, which is the all-columns rule;
- `dbarts(y ~ x1 + x2 + g, df, test = <test.df with NA in x2>)` raises the same message, the fit-time twin;
- `sampler.keep$setTestPredictor(<matrix with NA in x2>)` raises it;
- `makeTestModelMatrix(data, <newdata with NA in x2>)` raises it, the exported caller;
- a `missing = "error"` fit refuses any NA newdata, naming the column - today it answers silently;
- `bart2(..., keepTrees = TRUE)` then `predict(fit, newdata)` raises through the S3 surface, so the refusal is proven to
  reach a fit and not only a sampler.
D8, in inst/tinytest/test-data-missing-matrix.R: the unnamed x/y matrix interface names the column by INDEX
(`"column 2"`, 1-based).
D8, in inst/tinytest/test-data-mixed.R or test-predict-sparse.R: a `dbartsMixedMatrix` newdata with an NA in a dense
column and one in a sparse column, against a container training design - the container probe's only coverage. A sparse
column stores missing as NaN ([[R/mixedMatrix.R:50-56@7a8c7286]]), so both storages must be exercised.
D8, in inst/tinytest/test-hurdle.R: a column whose training NAs lie only in the zero rows is complete for the positive
part, so `predict(hurdleFit, newdata)` with an NA there raises - the one composition with a surprising answer, pinned so
it is not later read as a bug.

Warning handling: none of these emits a warning, so no `expect_warning` count changes; every new assertion is
`expect_error(..., pattern = )` with a pattern short enough to survive a rewording (`"carried none in training"`,
`"does not support grouped random effects"`, `"heteroscedastic variance forest"`).

## 11. Consumer sweep (read-only, `git -C <repo> grep`)

- stan4bart, /Users/vdorie/Repositories/stan4bart branch `bartcore` (f9bca65). D6: CANNOT hit it - it resolves through
  `dbartsSpec` (stan4bart's `R/stan4bart_fit.R` lines 554-559) and explicitly refuses `variance` in `bart_args` (stan4bart's `R/stan4bart_fit.R` lines 544-548, "the parametric
  component draws the residual standard deviation"), and it never writes `bartcore.groups` (its only control attribute is
  `n.cuts`, stan4bart's `R/stan4bart_fit.R` line 524). D8: it calls the FLAT `dbarts_sampler_predict` (stan4bart's `src/init.cpp` line 349) with a dense source, so it is
  exactly the consumer the C++ backstop serves - but it forbids `na.action = na.pass` outright (stan4bart's `R/stan4bart.R` lines 41-42) and
  defaults to `na.omit` (:6), so no NA reaches a model frame it hands to dbarts. No migration cost, no expected hit.
- bartCause, /Users/vdorie/Repositories/bartCause branch `dbarts-1.0` (d825cfc). D6: reaches the grouped surface only
  through `rbart_vi` (bartCause's `R/responseFit.R` lines 121 and 210; its `R/treatmentFit.R` line 107), which has no `variance` formal - unreachable.
  D8: its predict calls (bartCause's `R/generics.R` lines 143-176) pass `x.new` built from the user's confounders; it validates the treatment
  column for NA (bartCause's `R/generics.R` line 209) but not the rest, so a bartCause user with NA confounders on a training-complete column is now
  refused instead of silently answered - the intended change, and BART with NA predictors is new in 1.0, so nothing there
  can be relying on it. Its fit-time `data@x.test <-` assignment (bartCause's `R/responseFit.R` lines 180-189) is outside the validated
  surface by section 5's contract, so that path is unchanged either way. No code change expected; run its suite as the
  gate.
- treatSens, branch `dbarts-1.0` (1db3d89). D6: creates through
  `dbarts_sampler_create` (treatSens's `src/bartTreatmentModel.cpp` line 63 and `src/sensitivityAnalysis.cpp` line 260) with neither attribute -
  unreachable. D8: it IS a test-surface user - `dbartsData(x, y, x.test)` (treatSens's `R/cibart.R` line 56) and `sampler$setTestPredictor`
  (treatSens's `R/treatSensBART.R` lines 228, 245 and 253) - so the R check runs on its data; its designs are simulated numeric confounders
  with no NAs, and the per-column `setTestPredictor(newColumn, i)` form is not on the checked route anyway. No expected
  hit.
- bairrtt, /Users/vdorie/Repositories/bairrtt branch `main` (6167423). D6: no dbarts composition surface at all. D8: it
  predicts through the R5 `$predict(frame)` six times (bairrtt's `R/irt_causal_bart.R` lines 568-713), which is on the checked route, but
  it refuses NA covariates at ingestion itself (`as_covariate_frame`, bairrtt's `R/engine.R` lines 303-305, "'x' must not contain NA") and
  its documented missingness is in the RESPONSE matrix, which is not a predictor. No hit.

Total expected lockstep migration cost: zero lines in any consumer. Rebuild stan4bart and treatSens against commit 1
(the bridge change) as the branch habit; bartCause and bairrtt need only their suites.

## 12. Commit plan and gates

Two commits, D6 then D8 (prerc-surface-freeze's own order; they are independent).

1. `composition-refusals`. R/spec.R (one refusal in the variance block); src/R_interface_bartcore.cpp (one refusal after
   `applyVarianceAttributes`); inst/tinytest/test-heteroscedastic.R and test-capi.R (new assertions);
   docs/design/heteroscedastic.md (section 15 corrected, section 17 added); docs/design/grouped-random-effects.md (one
   cross-reference); [[docs/design/feature-matrix.md:186@7a8c7286]], [[docs/design/feature-matrix.md:187@7a8c7286]], [[docs/design/feature-matrix.md:727-731@7a8c7286]], [[docs/design/feature-matrix.md:1010-1011@7a8c7286]]; one NEWS UPGRADING item.
   `R CMD INSTALL .`; `cd tests/cpp && make && ./test_bartcore` must stay at its current 268 ok (nothing under
   src/bartcore/ moves, so it cannot change).
2. `predict-na-refusal`. R/data.R (two helpers plus one call at `validateXTest`'s tail); src/C_interface.cpp (one helper
   plus one call in `validateTestSource`); [[man/dbarts.Rd:116-118@7a8c7286]]; [[man/bart2.Rd:311-313@7a8c7286]]; man/makeind.Rd;
   [[inst/NEWS.Rd:720-725@7a8c7286]] plus one UPGRADING item; inst/tinytest/test-data-missing.R, test-data-missing-matrix.R,
   test-data-mixed.R (or test-predict-sparse.R), test-hurdle.R, test-sparse-factor.R (the rewritten fixture).
   `R CMD INSTALL .`; tests/cpp unchanged.

Both commits: `Rscript tools/check-doc-freshness.R .` after the edits, re-aligning every strict miss from the
`git diff -U0` line map by editing the docs/design anchors in place so each file's line count is invariant. The scope is
larger than a first pass suggests, because docs/design cites this bridge by its `RIB:` abbreviation as often as by name:
70 anchors sit past [[src/R_interface_bartcore.cpp:3032@7a8c7286]] across 11 docs/design files (feature-matrix.md alone holds 53),
14 past [[R/spec.R:527@7a8c7286]] across 4, and 35 past [[R/data.R:283@7a8c7286]] and [[src/C_interface.cpp:257@7a8c7286]] across 11. The anchor pass runs LAST
in each commit, as practiced. Baseline at this tip: 0 FAIL / 68 WARN.

Equivalence, both commits: `benchmarks/R/equivalence.R compare` against `benchmarks/baselines/equivalence-736bfb05.rds`,
plus the sibling scripts `benchmarks/R/bcf-equivalence.R` against `bcf-equivalence-6e3b9fb8.rds` and
`benchmarks/R/multinomial-equivalence.R` against `multinomial-equivalence-4d9a3337.rds` (43 / 12 / 11 channels). ALL
BITWISE IDENTICAL, zero re-records. The expectation is verified, not assumed, for D8's C++ touch specifically: the new
`refuseTestMissingness` is called only from `validateTestSource`, reached only from `dbarts_sampler_predict` and
`dbarts_sampler_setTestPredictors` - flat entries no equivalence scenario calls - and it draws nothing, allocates only
the reader the replay builds anyway, and returns before any engine call. D6's C++ touch is a branch on two `size_t`
creation options, taken before the chain exists. Neither commit edits src/bartcore/, the RNG stream, or any recorded
channel; if any equivalence channel moves, something in this slice reached the sampler and the commit is wrong.

`bench-sampler` is not required: no sampling path is touched. The one measurable addition is D8's per-predict
`anyNA(newdata)` scan, which is on no benchmark's inner loop; run
`benchmarks/R/bench-sampler.R compare benchmarks/baselines/bench-sampler-ab1dc52.csv` once at the end of the slice on a
quiet machine if the RC tip wants a clean sheet.

One extra gate, D6 only: `Rscript benchmarks/R/composition-matrix.R` before AND after. Before, it reports the two `?`
cells resolving to CONSTRUCTS (section 2's reading). After, with the cells written as `R`, it probes neither and reports
zero disagreements - which is also the check that the matrix edit and the code agree.

## 13. Settled sub-choices

Three, settled by orchestrator ruling under the standing grant, 2026-08-25, and written into the sections above rather
than left open.

1. D8 covers `dbartsData(test = )` and `$setTestPredictorAndOffset`, not `predict` alone, by placing the check in
   `validateXTest` (section 8a). The fact is about the data and the trees, not about which function the rows arrived
   through, and `?dbarts`'s `missing` item speaks of "test data", not of `predict`. Two consequences ride with it and are
   in the plan: [[test-sparse-factor.R:798-806@8d80ab01]]'s fixture is rewritten (section 10) and man/makeind.Rd gains a line
   (section 9), since `makeTestModelMatrix` is an exported caller. Direct slot assignment (`data@x.test <- `) stays
   outside the contract, as it is for every other check `validateXTest` performs.
2. D6 ships the bridge backstop (section 3b), not the R refusal alone. Three of the five entrances - `new("dbartsSampler",
   ...)`, `dbarts_sampler_create`, and re-creation from a control that acquired the attribute through `$setControl` -
   never reach `resolveSamplerSpec`, and all five reach `createHolder`, so five lines there are what make the refusal
   true of the surface rather than of one path. No third site in `$initialize`: redundant with the backstop.
3. D8 keeps the flat backstop (section 8b), and the R-side message names ALL offending columns, 1-based, bounded at five
   plus a count, with the scan returning silently when every NA sits in a column that carried training NAs.

## Draft docs/plans/INDEX.md row

## Landing note (2026-08-25)

Landed as three commits, pushed together as 936825d7: fe0b3292 (D6: the resolveSamplerSpec refusal, the createHolder
backstop, the heteroscedastic.md section 15 in-line correction plus its new section 17 door memo, feature-matrix [f30]
rewrite, tests at the R and flat entrances), de18ef2b (D8: refuseTestMissingness at validateXTest's tail with the flat
backstop in validateTestSource, the hurdle construction-time pre-check, Rd wording in dbarts/bart2/makeind, tests),
936825d7 (docs/design anchor re-alignment by the 93f3155c..de18ef2b diff line map; stamps to de18ef2b). Implementer
gates per phase; independent battery on a git-archive snapshot of 936825d7: tinytest 7375/0, tests/cpp 268+69 ok,
equivalence 43/12/11 identical, check-rc-codoc 42, freshness 0 FAIL / 69 WARN (the +1 over the 68 baseline is an
advisory intrinsic to [f30]'s own prose), NEWS 341, air + lint_package clean, as-cran 1 NOTE, census zero, both
refusal messages probed live.

Deviations from the design, both ruled at implementation under the standing grant: (1) the D8 refusal reaches hurdle
CONSTRUCTION, not just predict - bart2Hurdle forces the positive part's test channel to the full design, so a
covariate NA only on zero-outcome rows now refuses at fit with its own message (refuseHurdlePositiveMissingness,
before either component builds; the generic replay refusal stands as backstop; a column NA on both row sets stays
constructible). The exemption branch was rejected: those replayed rows feed the fit's own combined channels, the
defect class D8 refuses. (2) head() was avoided in the two message builders (no utils import), spelled
labels[seq_len(min(5L, length(labels)))].

Residue for the record: sparseFactor() refuses NA outright, so test-sparse-factor.R's D8 fixture uses a dense factor
column; benchmarks/R/composition-matrix.R reports two pre-existing DISAGREEMENTS unrelated to this slice (multinom
dbarts5 "no base fixture recipe"; logistic setWeights integer-weights - fixed at HEAD, both harness bugs, not model gaps); multinomial-mutation-arc.md sections 5-8 and
model-space-survey.md hold drifted anchors inside frozen text ([[C_interface.cpp:460@936825d7]]/462, data.R ranges,
[[model-space-survey.md:368-369@936825d7]]), left by rule; direct data@x.test slot assignment bypasses the R refusal (validated
surfaces + flat backstop are the contract) - bartCause's `responseFit.R` lines 180-189 uses that route at fit time.
Consumers: zero migration lines (stan4bart refuses variance in bart_args and forbids na.pass; treatSens and bairrtt
feed complete designs).
