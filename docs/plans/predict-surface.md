# predict() surface unification (D1), surface smalls (D9), shim removal (D5), saved-tree refusal (D2)

Status: LANDED 2026-08-25 at 78f334c1 (design record 716acd1f; code tip ed43deef).

Spec: docs/plans/prerc-surface-freeze.md D1, D9, D5, D2 and its Sequencing line ("then D1, D9, D5, D2"). TODO
`predict-signature-unification`, `surface-smalls`, `deprecation-shim-removal`, `predict-refusal-names-cure` - four commits, one
slice, since D9's fifteen new usage entries and D5's and D2's edits all land inside functions D1 reshapes. Evidence:
review-2026-08-24/memos/prerc-lens1-surface.md B1, B2, B3, B5, B6, C1 and prerc-lens2-backlog.md P2, re-anchored live (the memos
were written at 7a8c7286; every line number below is b46add06's). Revised 2026-08-25 after an independent critique and
orchestrator rulings under the standing grant (sections 3, 5, 6, 7, 9, 11, 13). No sampling code moves, no RNG consumption
changes on any path an existing test walks; zero baseline re-records expected.

## 1. Scope, and what stays put

In scope: the six S3 `predict` methods, the six S3 `fitted` methods, the seven S3 `extract` methods, the two S3
`survivalProbabilities` methods that carry `group.by`, and the four R5 forest readers. Out of scope, stated so it is not
rediscovered: the reference class's own `$predict(x.test, offset.test, n.threads)` ([[R/dbarts.R:1084@b46add06]]) and
`$predictForests(x.test, offset.test, n.threads)` ([[R/dbarts.R:1147-1151@b46add06]]) keep `offset.test`. Those are the engine's own terms - the same
vocabulary as `$setTestOffset(offset.test)` ([[R/dbarts.R:1629@b46add06]]), `$setTestPredictorAndOffset` ([[R/dbarts.R:1573@b46add06]]) and `dbartsData(offset.test =)` - and
lens 1 B10 declined that rename with 201 in-repo uses behind it. The fit-time `bart2(offset.test =)`/`rbart_vi(offset.test =)`
formals and `dbartsData`'s slot keep the name for the same reason. D1 retires `offset.test` as a spelling on the S3 `predict`
surface ONLY, where it names a per-call argument, not a stored channel.

## 2. D1 census: every method as it stands

Taken from code, formals in declaration order. Registration is NAMESPACE:31-98 (67 S3 methods; `predict` on six classes,
`fitted` on six, `extract` on seven).

predict, R/generics.R:
- `predict.bart` [[R/generics.R:257-268@b46add06]] `(object, newdata, offset, weights, type = c("ev","ppd","bart","forest"), combineChains = TRUE,
  ci.level = NULL, forest = NULL, bases = NULL, n.threads = object$fit$control@n.threads, ...)`. `offset`/`weights` carry no
  default; [[R/generics.R:270-275@b46add06]] fills NULL through `missing()`.
- `predict.rbart` [[R/generics.R:2150-2161@b46add06]] `(object, newdata, group.by, offset, weights, type = c("ev","ppd","bart","ranef"), combineChains =
  TRUE, ci.level = NULL, n.threads = object$fit[[1L]]$control@n.threads, ...)`. `group.by` is positional THREE and has no
  default; `as.factor(group.by)` at [[R/generics.R:2199@b46add06]] is what raises when it is absent.
- `predict.bartMultinomial` [[R/generics.R:1181-1190@b46add06]] `(object, newdata, type = c("ev","ppd","bart","forest"), offset = NULL, combineChains =
  TRUE, ci.level = NULL, n.threads = object$fit$control@n.threads, ...)`.
- `predict.bartOrdinal` [[R/generics.R:1473-1481@b46add06]] `(object, newdata, type = c("ev","ppd","bart"), combineChains = TRUE, ci.level = NULL,
  n.threads = object$fit$control@n.threads, ...)` - no offset formal at all.
- `predict.bartNegbin` [[R/generics.R:1711-1720@b46add06]] `(object, newdata, type = c("ev","ppd","bart"), offset.test = NULL, combineChains = TRUE,
  ci.level = NULL, n.threads = object$fit$control@n.threads, ...)` - the one `offset.test` spelling; consumed at [[R/generics.R:1746-1751@b46add06]] as
  `bartcorePredict(list(ptr = ...), newdata, offset.test, n.threads)`.
- `predict.bartHurdle` [[R/generics.R:2107-2115@b46add06]] `(object, newdata, type = c("ev","ppd","prob","bart"), combineChains = TRUE, ci.level = NULL,
  n.threads = object$occupancy$fit$control@n.threads, ...)` - no offset formal; composes two `predict.bart` calls through
  `hurdleParts` ([[R/generics.R:1950-1986@b46add06]]), which forwards only `newdata`, `type` and `n.threads`, all by name ([[R/generics.R:1965-1979@b46add06]]), and no offset or
  weights at all.

fitted, R/generics.R: `.bart` [[R/generics.R:859-864@b46add06]], `.rbart` [[R/generics.R:2551-2556@b46add06]], `.bartMultinomial` [[R/generics.R:1092-1097@b46add06]], `.bartOrdinal` [[R/generics.R:1401-1406@b46add06]],
`.bartNegbin` [[R/generics.R:1665-1670@b46add06]], `.bartHurdle` [[R/generics.R:2073-2078@b46add06]]. None takes an offset or a `newdata`; all are `(object, type, [sample,]
[ci.level,] ...)`. `fitted.bartHurdle`'s `sample = "train"` is a bare string where every sibling uses a choice vector (lens 1 B8,
not in this slice).

extract, R/generics.R: `.bart` [[R/generics.R:434-442@b46add06]] is `(object, type, sample, combineChains, forest = NULL, contribution = FALSE, ...)` -
the only one with the per-forest pair. `.rbart` [[R/generics.R:2389-2394@b46add06]], `.bartMultinomial` [[R/generics.R:975-980@b46add06]], `.bartOrdinal` [[R/generics.R:1306-1311@b46add06]],
`.bartNegbin` [[R/generics.R:1585-1590@b46add06]] and `.bartHurdle` [[R/generics.R:2010-2015@b46add06]] are all `(object, type, sample, combineChains, ...)`.
`.dbartsSampler` [[R/generics.R:2543@b46add06]] is `(object, type = "predictors", ...)` - no `sample`, no `combineChains`, no `refuseUnusedGenericArgs`,
its own one-token validation at [[R/generics.R:2544-2546@b46add06]]; it documents the SAMPLER-class read (man/extract.dbartsSampler.Rd), not a fit's
channels, and is OUT of D1's scope. D1 changes no `extract` signature at all.

Offset spellings, package-wide: `offset` (predict.bart, predict.rbart, predict.bartMultinomial; the multinomial one is an
`nrow(newdata)` x K matrix, documented [[man/bart2.Rd:126-129@b46add06]]), `offset.test` (predict.bartNegbin; and the two R5 methods and every
fit-time channel, all out of scope), `binaryOffset` (bart()'s creation formal, out of scope). Two spellings on the S3 predict
surface, one after D1.

`group.by` appears in exactly one predict signature (predict.rbart [[R/generics.R:2153@b46add06]]) and CAN be matched positionally today - and is, at 21
in-repo test call sites and one consumer site (section 12). It appears once more on the S3 surface, positional FOURTH on
`survivalProbabilities.rbart` ([[R/bart.R:2551-2557@b46add06]]), which forwards it to predict.rbart at [[R/bart.R:2581-2587@b46add06]]; that method takes the same
named-only treatment (section 3), since it is the same argument on the same class.

Unknown-name refusal already exists and is partial: `refuseUnusedGenericArgs` ([[R/generics.R:1868-1883@b46add06]]) intersects
`names(reasons)` with `names(dots)` and stops on `supplied[1L]`, so the FIRST name in the reasons list wins when a call supplies
several - which is why section 3 gives ordinal and hurdle their own `offset.test` wording rather than composing another class's.
predict.bart calls it at [[R/generics.R:285@b46add06]] and predict.rbart at [[R/generics.R:2180@b46add06]], both with `predictOffsetUnusedArgs` ([[R/generics.R:253-255@b46add06]]), which holds
`offset.test` alone; the four own-class predicts call it at [[R/generics.R:1193-1198@b46add06]], [[R/generics.R:1483-1488@b46add06]], [[R/generics.R:1722-1727@b46add06]], [[R/generics.R:2117-2122@b46add06]] with lists holding
only `forest`/`contribution` ([[R/generics.R:964@b46add06]], [[R/generics.R:1301@b46add06]], [[R/generics.R:1580@b46add06]], [[R/generics.R:2005@b46add06]]). So today `predict(negbinFit, nd, offset = o)` and
`predict(ordinalFit, nd, offset = o)` still vanish into `...`. Memo B1's claim that "predict.bart/predict.rbart don't call
refuseUnusedGenericArgs at all" is STALE: [[R/generics.R:253-255@b46add06]] and its two call sites landed after the memo. What has not landed is the other
half - the offset and weights names on the four own-class lists.

## 3. D1 after: the signatures

The guarantee is the PREFIX `(object, newdata, type)`, identical on all six, plus one canonical relative order for the tail:
`offset, weights, combineChains, ci.level, forest, bases, n.threads`, a class simply omitting the names it has no channel for.
`n.threads` stays the last positional formal on every method. `group.by` moves AFTER `...`, which is R's own named-only
mechanism: a formal following `...` can be matched by its full name only, never positionally and never by partial match. It keeps
no default, so a missing one still raises - but by name, see below.

    predict.bart(object, newdata, type = c("ev", "ppd", "bart", "forest"),
                 offset = NULL, weights = NULL, combineChains = TRUE, ci.level = NULL,
                 forest = NULL, bases = NULL, n.threads = object$fit$control@n.threads, ...)

    predict.rbart(object, newdata, type = c("ev", "ppd", "bart", "ranef"),
                  offset = NULL, weights = NULL, combineChains = TRUE, ci.level = NULL,
                  n.threads = object$fit[[1L]]$control@n.threads, ..., group.by)

    predict.bartMultinomial(object, newdata, type = c("ev", "ppd", "bart", "forest", "class"),
                            offset = NULL, combineChains = TRUE, ci.level = NULL,
                            n.threads = object$fit$control@n.threads, ...)

    predict.bartOrdinal(object, newdata, type = c("ev", "ppd", "bart", "class"),
                        combineChains = TRUE, ci.level = NULL,
                        n.threads = object$fit$control@n.threads, ...)

    predict.bartNegbin(object, newdata, type = c("ev", "ppd", "bart"),
                       offset = NULL, combineChains = TRUE, ci.level = NULL,
                       n.threads = object$fit$control@n.threads, ...)

    predict.bartHurdle(object, newdata, type = c("ev", "ppd", "prob", "bart"),
                       combineChains = TRUE, ci.level = NULL,
                       n.threads = object$occupancy$fit$control@n.threads, ...)

    survivalProbabilities.rbart(object, times, newdata = NULL, combineChains = TRUE, ..., group.by)

The two `"class"` tokens are D9b's (section 6); they ride in commit 2 with the rest of that vocabulary work, NOT here, so each
commit's `\usage` matches its own formals. `survivalProbabilities.bart` ([[R/bart.R:2500-2506@b46add06]]) is already
`(object, times, newdata, combineChains, ...)`, so the rbart move makes that pair's positional prefix uniform too.

Body edits that follow:
- predict.bart: delete [[R/bart.R:270-275@b46add06]], the `missing(offset)`/`missing(weights)` block - the two defaults replace it exactly, since both
  helpers that receive them (`predictForest` [[R/bart.R:644-651@b46add06]], `predictBlend` [[R/bart.R:784-794@b46add06]]) are already given NULL on the missing path and
  are called positionally with the same locals ([[R/bart.R:319-326@b46add06]], [[R/bart.R:333-343@b46add06]]), so neither helper's own signature moves.
- predict.rbart: delete [[R/bart.R:2182-2187@b46add06]], the same block. Add, immediately after the saved-tree refusal and BEFORE `validateType`:

        if (missing(group.by)) {
          stop("'group.by' must be given by name: predict on an rbart fit needs the ",
               "test rows' grouping factor, and it is no longer the third positional ",
               "argument")
        }

  Ordering matters: `predict(fit, x, g)` now binds `g` to `type`, and without this check the caller would get `validateType`'s
  "type must be in 'ev', 'ppd', 'bart', 'ranef'" ([[R/bart.R:1845@b46add06]]) with nothing pointing at the real cause. With it, the old positional
  call gets a message that names the fix. This is the migration signpost for the 21 test sites and the one consumer site.
- predict.bartNegbin: rename `offset.test` -> `offset` at [[R/bart.R:1715@b46add06]] and at its single use [[R/bart.R:1749@b46add06]].
- survivalProbabilities.rbart: move `group.by` after `...` ([[R/bart.R:2551-2557@b46add06]]). Its own `missing(group.by)` refusal at
  [[R/bart.R:2576-2578@b46add06]] stays where it is, its text gaining the naming rule: `"'group.by' must be given by name when 'newdata' is
  given"`. Its forward at [[R/bart.R:2581-2587@b46add06]] already passes `group.by = as.factor(group.by)`, an exact full-name match, so it survives
  predict.rbart's own move untouched.
- [[R/generics.R:227-233@b46add06]], the comment above `validatePredictThreads`, justifies `n.threads`-last by "consumers call these methods
  positionally, so an earlier insertion would rebind their arguments" - the practice D1 retires. Rewrite the last sentence to the
  surviving constraint: "Every predict method takes it as its LAST positional formal, so the argument a caller is most likely to
  supply by position - `type` - stays third on every one of them."

Refusal lists, completing lens 1 B1's second half. `predictOffsetUnusedArgs` ([[R/generics.R:253-255@b46add06]]) keeps its name and one entry, with the
comment rewritten (it currently justifies itself by naming predict.bartNegbin's formal, which no longer exists):

    # One offset spelling is live across every predict method - 'offset'. The
    # fit-time channels keep 'offset.test' (dbartsData, bart2, rbart_vi, and the
    # sampler's own $predict), so a caller carrying that name here would otherwise
    # vanish into '...' with the offset silently dropped instead of applied.
    predictOffsetUnusedArgs <- list(
      offset.test = "this fit's out-of-sample offset argument is named 'offset'"
    )

Two more lists beside it. The no-offset one carries BOTH spellings with the same wording, rather than composing the list above:
`refuseUnusedGenericArgs` reports the first name in `names(reasons)`, so a composed list would answer
`predict(ordinalFit, nd, offset.test = o)` with "named 'offset'", pointing at an argument this class does not have.

    # A per-observation weight scales the posterior-predictive DRAW of a fit whose
    # noise the caller can rescale - gaussian sigma, a logistic trial count. A
    # count, category or two-part draw comes from its own law with no such factor,
    # so the argument has nothing to act on here.
    predictWeightsUnusedArgs <- list(
      weights = "this family's posterior-predictive draw takes no per-observation weight"
    )
    # An offset shifts the latent at rows the sampler never saw, and these two
    # families replay their trees with no offset channel at all, so either spelling
    # would be dropped rather than applied.
    noPredictOffsetReason <- paste0(
      "this fit has no out-of-sample offset channel; predict replays the ",
      "offset-free surface"
    )
    predictNoOffsetUnusedArgs <- list(
      offset = noPredictOffsetReason,
      offset.test = noPredictOffsetReason
    )

Wiring, at the six call sites: bart [[R/generics.R:285@b46add06]] and rbart [[R/generics.R:2180@b46add06]] unchanged (`predictOffsetUnusedArgs`; `weights` is a live formal on
both). multinomial [[R/generics.R:1193-1198@b46add06]] -> `c(multinomialUnusedArgs, predictOffsetUnusedArgs, predictWeightsUnusedArgs)`; negbin
:1722-1727 -> `c(negbinUnusedArgs, predictOffsetUnusedArgs, predictWeightsUnusedArgs)`; ordinal [[R/generics.R:1483-1488@b46add06]] ->
`c(ordinalUnusedArgs, predictNoOffsetUnusedArgs, predictWeightsUnusedArgs)`; hurdle [[R/generics.R:2117-2122@b46add06]] ->
`c(hurdleUnusedArgs, predictNoOffsetUnusedArgs, predictWeightsUnusedArgs)`. The predict-only names are combined AT the call site,
never folded into the shared `*UnusedArgs` lists, because those same lists serve extract/fitted/residuals ([[R/generics.R:988-993@b46add06]], [[R/generics.R:1100-1105@b46add06]],
:1140-1145, [[R/generics.R:1315-1320@b46add06]], [[R/generics.R:1408-1413@b46add06]], [[R/generics.R:1450-1455@b46add06]], [[R/generics.R:1594@b46add06]], [[R/generics.R:1672-1677@b46add06]], [[R/generics.R:1692-1697@b46add06]], [[R/generics.R:2019@b46add06]], [[R/generics.R:2081@b46add06]], [[R/generics.R:2094-2099@b46add06]]), which have no offset
or weights formal to protect.

Landed 74e2e050 ("the last predict slot-4 collision"), not post-1.0 residue - the integration tip this doc's
landing was branched from: `predict.bartOrdinal`/`predict.bartHurdle` gained an `offset = NULL` formal ahead of
`combineChains`, refused non-NULL by `refusePredictOffsetChannel`, reversing section 13's settled call below.
Superseded record only: Ordinal and hurdle do NOT gain a dummy refusing `offset` formal to close it (settled, section 13): the by-name
refusal above delivers the same message at the same moment, and a formal would buy positional uniformity for an argument D1 is
simultaneously telling callers to pass by name. The split that remains is "has an out-of-sample offset channel", not the
arbitrary one B2 found; every name means one thing everywhere, and every wrong guess is refused by name.

## 4. D1: the call sites that move

Inside R/: NONE. Every internal predict call passes everything past `newdata` by name - `hurdleParts` [[R/generics.R:1965-1979@74e2e050]],
`survivalProbabilities.bart` [[R/bart.R:2460@74e2e050]] and [[R/bart.R:2532@74e2e050]], `survivalProbabilities.rbart` [[R/bart.R:2581-2587@74e2e050]]. The three positional
`extract(object, type, sample, ...)` calls ([[R/generics.R:869@74e2e050]], [[R/generics.R:2570@74e2e050]], [[R/generics.R:2596@74e2e050]]) are extract's own order, which D1 does not touch.

Verified by parsing every .R file under R/, inst/tinytest/ and tests/ and reporting each `predict`/`fitted`/`extract`/`residuals`
call with three or more unnamed arguments; the complete result is section 11's list plus those three.

## 5. D9a: one `forest` argument on the four readers

The four are R5 methods on `dbartsSampler`, not the flat-C entries D4 renamed (`dbarts_sampler_getForestFits` and siblings) -
those already agree and are not touched here. Live names and conventions, R/dbarts.R: `getForestFits(forest)` [[R/generics.R:1691-1695@74e2e050]], no
default; `getForestVariableCounts(forest)` [[R/generics.R:1716-1729@74e2e050]], no default; `getForestAmplitudes(forest = NULL)` [[R/generics.R:1707-1715@74e2e050]], NULL = all
stacked; `getCalibration(forest = 1L)` [[R/generics.R:1730-1734@74e2e050]]. `setCalibration(..., forest = 1L, ...)` [[R/generics.R:1735-1804@74e2e050]] is a fifth site and a
WRITER: it is refused on every multi-forest sampler ([[src/R_interface_bartcore.cpp:4205-4210@74e2e050]], because a calibration map owns those
forests' scales), so `1L` is the only value that can succeed and NULL would name nothing writable. It keeps `forest = 1L`, and
the Rd item that pairs it with `getCalibration` ([[man/dbartsSampler-class.Rd:201@74e2e050]]) is rewritten to say why.

After: all four default to `forest = NULL` meaning every forest. The NULL path is NEW on three of them, so no existing shape
moves and there is no compatibility claim to make; what must hold is that a SINGLE-forest sampler's `forest = NULL` read is
bitwise what `forest = 1` returns, which is what lets the default change silently. The stacked shapes are stated one by one in
the Rd rather than folded into a single sentence, because the four readers return four different quantities:
- `getForestAmplitudes(NULL)`: unchanged - sum(q) x numChains, forest-major within the row margin, ragged by construction.
- `getForestFits(NULL)`: n x numChains at one forest; n x numForests x numChains above it, the forest margin between the rows
  and the chains, where `getForestAmplitudes` already puts it and where `run()`'s own multi-forest widening puts it
  (`forestFits` n.obs x n.forests x n.samples x n.chains, [[man/dbartsSampler-class.Rd:408@74e2e050]]).
- `getForestVariableCounts(NULL)`: numPredictors x numChains, then numPredictors x numForests x numChains, same placement; the
  predictor rownames stay on margin 1.
- `getCalibration(NULL)`: numChains x 12, then numChains x 12 x numForests - the forest margin appended LAST, because this
  reader's row margin IS the chain axis and there is no row margin to sit behind. The column dimnames and the `leaf.model`
  attribute ride unchanged.

The single-forest identity is what protects the bare `getCalibration()` reads already in the tree: 22 calls across
inst/tinytest and benchmarks (23 grep lines, one of them the comment at [[inst/tinytest/test-bcf-family.R:12@74e2e050]]) -
[[test-augmentation.R:217@74e2e050]], [[test-augmentation.R:245@74e2e050]], [[test-augmentation.R:246@74e2e050]]; [[test-calibration-midchain.R:50@74e2e050]], [[test-calibration-midchain.R:187@74e2e050]], [[test-calibration-midchain.R:262@74e2e050]], [[test-calibration-midchain.R:327@74e2e050]], [[test-calibration-midchain.R:330@74e2e050]], [[test-calibration-midchain.R:336@74e2e050]], [[test-calibration-midchain.R:350@74e2e050]], [[test-calibration-midchain.R:498@74e2e050]], [[test-calibration-midchain.R:547@74e2e050]], [[test-calibration-midchain.R:549@74e2e050]];
[[test-embedding-recipes.R:68@74e2e050]], [[test-embedding-recipes.R:217@74e2e050]], [[test-embedding-recipes.R:242@74e2e050]], [[test-embedding-recipes.R:266@74e2e050]], [[test-embedding-recipes.R:267@74e2e050]]; [[benchmarks/R/backfit-exact.R:143@74e2e050]]; [[benchmarks/R/geweke-mc.R:528@74e2e050]], [[benchmarks/R/geweke-mc.R:531@74e2e050]], [[benchmarks/R/geweke-mc.R:562@74e2e050]] -
plus two in vignettes R CMD check builds ([[vignettes/dbarts-as-a-component.Rmd:189@74e2e050]],
[[vignettes/gibbs_sampler_mixture_model.Rmd:247@74e2e050]]). Every one is on a single-forest sampler and indexes the result as a matrix. A
heteroscedastic sampler counts as single-forest here - its variance forest is a separate member ([[src/bartcore/chain.hpp:898@74e2e050]]
`forests_.size()`, [[src/bartcore/chain.hpp:906-912@74e2e050]] the variance forest beside it), so `numForests` is 1 - which is what makes the default safe for that
family too.

R needs one fact it does not have: the forest count. `forestIndexFrom` ([[src/R_interface_bartcore.cpp:3920-3927@74e2e050]]) is the only
forest decoder and it rejects NULL (`Rf_asInteger(R_NilValue)` is NA_INTEGER, which casts to an out-of-range `size_t`); the flat
API has `dbarts_sampler_numForests` ([[src/C_interface.cpp:997-999@74e2e050]]) but no bridge twin; and `data@bases`/`dataCounts` are
documented CAPABILITY probes, "deliberately not a forest count" ([[R/bartcore.R:16-23@74e2e050]], [[R/bartcore.R:45-52@74e2e050]]), so deriving one R-side would
misfire on exactly the samplers this serves. Add ONE bridge entry and do the stacking in R - the three per-forest readers keep
their current bodies untouched, which is the smaller and more testable change:

    // The sampler's forest count, the R twin of dbarts_sampler_numForests. The R5
    // readers stack their per-forest reads at forest = NULL and need the bound; a
    // capability probe cannot supply it, since a plain single-forest sampler
    // answers no to every one of them.
    SEXP bartcore_numForests(SEXP ptrExpr) {
      BartcoreHolder& holder(holderFromExpression(ptrExpr));
      return Rf_ScalarInteger(
        static_cast<int>(holder.sampler->shape().numForests));
    }

Placed beside `bartcore_getForestAmplitudes` ([[src/R_interface_bartcore.cpp:4052-4077@b46add06]]), declared in src/R_interface_bartcore.hpp beside [[src/R_interface_bartcore.hpp:26@b46add06]]
(`SEXP bartcore_numForests(SEXP ptr);`) and registered in src/R_interface.cpp beside [[src/R_interface.cpp:192-193@b46add06]] as
`DEF_FUNC("dbarts_bartcore_numForests", bartcore_numForests, 1)`. It stays an INTERNAL read - a package helper in R/bartcore.R,
not a new R5 method, so D9 adds no public surface the freeze would then lock:

    # The sampler's forest count. A COUNT, not a capability probe:
    # samplerCarriesAmplitudes and samplerCarriesCounts each answer only for their
    # own model, and neither sees a plain single-forest sampler.
    bartcoreNumForests <- function(ptr) .Call(C_dbarts_bartcore_numForests, ptr)

The R side of each reader, on `getForestAmplitudes`'s existing idiom `if (is.null(forest)) NULL else resolveForestIndex(forest)`
([[R/dbarts.R:1710-1714@74e2e050]]). `getForestFits`:

    getForestFits = function(forest = NULL) {
      "<docstring, gaining the forest = NULL shape>"
      ptr <- getPointer()
      if (!is.null(forest)) {
        return(.Call(C_dbarts_bartcore_getForestFits, ptr, resolveForestIndex(forest)))
      }
      # the bridge counts forests from 0, as resolveForestIndex converts to
      numForests <- bartcoreNumForests(ptr)
      blocks <- lapply(
        seq_len(numForests),
        function(f) .Call(C_dbarts_bartcore_getForestFits, ptr, f - 1L)
      )
      if (numForests == 1L) {
        return(blocks[[1L]])
      }
      result <- array(0.0, c(nrow(blocks[[1L]]), numForests, ncol(blocks[[1L]])))
      for (f in seq_len(numForests)) {
        result[, f, ] <- blocks[[f]]
      }
      result
    }

`getForestVariableCounts` is the same with `array(0L, ...)` (the per-forest reads are INTSXP and the assignment preserves the
type), keeping its `predictorNames <- colnames(data@x)` block ([[R/dbarts.R:1724-1727@74e2e050]]) applied to the stacked result unchanged and
UNGUARDED: `rownames<-` works on a 3-d array, filling `dimnames[[1]]` and leaving the other two NULL, so the predictor names land
on margin 1 in both shapes. `getCalibration` stacks into `array(0.0, c(nrow(first), ncol(first), numForests))` with
`dimnames(result) <- list(NULL, colnames(first), NULL)` and then
`attr(result, "leaf.model") <- attr(first, "leaf.model")` - the tag is a property of the sampler, identical on every forest, so
the first forest's is the sampler's.

Nothing below the bridge moves: `forestTotalFits`, `forestVariableCounts` and `forestCalibration` are read per (chain, forest)
exactly as now, and the three existing entries are not edited at all. `resolveForestIndex` ([[R/bartcore.R:1051-1057@74e2e050]]) keeps its
message and is called only on the non-NULL branch.

## 6. D9b: one `type` vocabulary per class

Current, taken from the choice vectors (predict / fitted / extract):
- bart: ev ppd bart forest / ev ppd bart / ev ppd bart loglik trees forest
- rbart: ev ppd bart ranef / ev ppd bart ranef / ev ppd bart loglik ranef trees
- bartMultinomial: ev ppd bart* forest* / ev class bart* / ev ppd bart* forest* loglik
- bartOrdinal: ev ppd bart / ev class bart / ev ppd bart loglik
- bartNegbin: ev ppd bart / ev bart / ev ppd bart loglik
- bartHurdle: ev ppd prob bart / ev prob bart / ev ppd prob bart loglik
(* named only to be refused by `refuseMultinomialLatentType`, [[R/generics.R:932-949@74e2e050]].)

After, and the rule that produces it: `predict` and `fitted` carry one vocabulary per class; `extract` carries the DRAW channels,
which is that vocabulary plus the draws-only ones and minus the reductions. Four exceptions, each stated in the Rd rather than
left to be inferred:
- `"trees"` and `"loglik"` are extract-only: neither has a per-observation posterior mean.
- `"forest"` is predict/extract-only: its value carries a forest margin rather than one per-observation channel, and
  `fitted.bart` [[R/generics.R:877-881@74e2e050]] reduces the LAST margin, which for `extractForest`'s value is the forest and not the observation. So
  `fitted` does not gain it (settled, section 13).
- `"class"` is predict/fitted-only: `extract` returns draw channels and a class is a reduction OVER draws (the argmax of the
  posterior-mean probability matrix), so there is no draws-shaped value for it to return.
- `"ppd"` is not on `fitted` for the two categorical families: their posterior-predictive draw is a category CODE
  (`multinomialPpdFromProbs`, [[R/generics.R:1067-1074@74e2e050]], reached at [[R/generics.R:1242-1244@74e2e050]] and [[R/generics.R:1543-1545@74e2e050]]), whose mean is not a quantity.

Four edits:
- `predict.bartMultinomial` [[R/generics.R:1184@74e2e050]] gains `"class"` LAST in its vector; `predict.bartOrdinal` [[R/generics.R:1476@74e2e050]] gains `"class"` last.
- `fitted.bartNegbin` [[R/generics.R:1667@74e2e050]] becomes `c("ev", "ppd", "bart")` and `fitted.bartHurdle` [[R/generics.R:2075@74e2e050]] becomes
  `c("ev", "ppd", "prob", "bart")` - `"ppd"` SECOND in both, since codoc compares the `\usage` default against the formal's and
  the two must be written the same way; second is where every sibling that has it puts it.

`fitted.bartHurdle` needs no body change - [[R/generics.R:2082@74e2e050]] already routes every type through `extract`, which has the arm.
`fitted.bartNegbin` [[R/generics.R:1678@74e2e050]] becomes a three-way switch:

    channel <- switch(
      type,
      bart = object$latent.train,
      ev = object$yhat.train,
      # the ppd arm is a draw, not a stored channel; extract pairs each mu with
      # its own draw's dispersion, and the mean over the observation margin below
      # is invariant to the chain layout it returns
      ppd = extract.bartNegbin(object, type = "ppd", sample = "train")
    )

extract.bartNegbin's ppd arm ([[R/generics.R:1623-1632@74e2e050]]) returns `dim(mu)`, so the observation margin stays last and [[R/generics.R:1682@74e2e050]]'s
`apply(channel, length(dim(channel)), mean)` is unchanged. Like `fitted.bart(type = "ppd")` before it, this consumes RNG - only
when a caller asks for it, so no existing test's stream moves.

`predict`'s `"class"` arm mirrors `fitted`'s, and the mirror is enforced by extracting the shared reduction rather than copying
it a third and fourth time. `fitted.bartMultinomial` [[R/generics.R:1110-1121@74e2e050]] and `fitted.bartOrdinal` [[R/generics.R:1425-1436@74e2e050]] are today the same seven
lines twice; replace both, and serve both new predict arms, with one helper pair placed beside `refuseMultinomialLatentType`:

    # The posterior-mean n x K probability matrix of a K-widened draws array
    # (observation margin next-to-last, category margin last in every chain
    # layout), and its argmax as a factor over the fit's own levels - the class
    # prediction fitted() and predict() share, so the two cannot drift.
    meanCategoryProbabilities <- function(probs, levels) {
      d <- length(dim(probs))
      meanProbs <- apply(probs, c(d - 1L, d), mean)
      dimnames(meanProbs) <- list(NULL, levels)
      meanProbs
    }
    categoryFromMeanProbabilities <- function(meanProbs, levels, ordered = FALSE) {
      factor(levels[max.col(meanProbs, ties.method = "first")],
             levels = levels, ordered = ordered)
    }

Placement in `predict`, and it is the one place the arm can go: AFTER the existing `ci.level` block and before the plain `probs`
return - multinomial after [[R/generics.R:1245-1251@74e2e050]] and before [[R/generics.R:1252@74e2e050]], ordinal after [[R/generics.R:1546-1549@74e2e050]] and before [[R/generics.R:1550@74e2e050]]. That reproduces exactly what
`fitted` does at [[R/generics.R:1107-1108@74e2e050]], where the `ci.level` early return precedes the mean/class reduction, so `type = "class"` WITH a
`ci.level` returns the band on the full probability draws rather than a factor, on the stated principle that the band is "taken
on the full probability draws before the class reduction so it is meaningful regardless of 'type'" ([[R/generics.R:1088-1091@74e2e050]]). Inserting the
arm before the block instead would return early and leave the widened `trailing` selection dead. With the placement right, the
two `trailing` selections at [[R/generics.R:1249@74e2e050]] and [[R/generics.R:1547@74e2e050]] do widen to `if (type %in% c("ev", "class")) 2L else 1L` and are live.
`combineChains` does not reach the class arm - the reduction is over every draw - and the Rd says so.

## 7. D9c: the alias-without-usage entries

Verified by matching NAMESPACE's `S3method` lines against every `\alias` and every `\method`/`\S3method` usage entry in man/.
67 registered S3 methods. TWO have no `\alias` anywhere (`print.bart`, `print.rbart` - lens 1 C3, not this slice). EIGHTEEN have
an `\alias` and no `\usage` entry, not fifteen. The plan's 15 are exactly C1's set, all in man/bart2.Rd: `extract`, `fitted`,
`predict`, `print`, `residuals` for each of `bartOrdinal` (aliases [[R/generics.R:10-14@74e2e050]]), `bartNegbin` ([[R/generics.R:17-21@74e2e050]]) and `bartHurdle` ([[R/generics.R:24-28@74e2e050]]) -
those get usage entries here. The other three are outside C1 and outside D9: `print.dbartsCompositionValidation`
(man/dbartsValidateComposition.Rd), `print.dbartsVarianceForest` and `format.dbartsVarianceForest` (man/varianceForest.Rd, lens 1
C8, deferred). Report the true count, take the 15.

Which of the 15 D1 makes moot: none of the entries, all of the divergences they would have exposed. `predict.bartNegbin`'s usage
line prints `offset` rather than `offset.test`; `predict.bartOrdinal`'s and `predict.bartHurdle`'s print the same
`(object, newdata, type, ...)` prefix as their siblings; `fitted.bartNegbin`'s prints the vocabulary D9b equalized. Writing these
entries at the OLD signatures would publish exactly the split the slice removes, which is why D1 is committed first.

## 8. D5: the two rbart shims

[[R/generics.R:2167-2178@74e2e050]], verbatim what is deleted:

    dotsList <- list(...)
    if (!is.null(dotsList[["value"]])) {
      warning("argument 'value' has been deprecated; use 'type' instead")
      type <- dotsList[["value"]]
      dotsList[["value"]] <- NULL
    }

    type <- foldTypeAliases(type)
    if (is.character(type) && length(type) > 0L && type[1L] == "post-mean") {
      warning("type of 'post-mean' for predict deprecated; use 'ev' instead")
      type[1L] <- "ev"
    }

:2179-2180 then become `type <- validateType(type, eval(formals(predict.rbart)$type))` and a
`refuseUnusedGenericArgs(list(...), "predict", "rbart", ...)` call - `dotsList` has no other reader, and `validateType` folds the
response/link aliases itself ([[R/generics.R:1843@74e2e050]]), so the standalone `foldTypeAliases` call goes with the block. One comment moves: [[R/generics.R:1826@74e2e050]]
("some also reject length-0 input; predict.rbart interposes a post-mean alias") loses its second clause.

The deletion is TOTAL but not silent (settled, section 13): `value` joins predict.rbart's reasons list, so the old spelling
errors by name instead of vanishing into `...`, which is the defect lens 1 B1 exists to close. A refusal is not a shim - it does
not accept the old spelling, it names it:

    # 'value' was predict.rbart's pre-1.0 name for 'type'. It is not accepted, only
    # refused by name, since a supplied one would otherwise choose the default
    # channel silently.
    rbartPredictValueUnusedArgs <- list(
      value = "predict's channel argument is named 'type'"
    )

so [[R/generics.R:2180@74e2e050]] reads `refuseUnusedGenericArgs(list(...), "predict", "rbart", c(predictOffsetUnusedArgs,
rbartPredictValueUnusedArgs))`. `"post-mean"` needs nothing: `validateType` reports "type must be in 'ev', 'ppd', 'bart',
'ranef'", which names the replacement.

Nothing else in the package references either shim: `git grep post-mean` finds only [[R/generics.R:1826@74e2e050]], [[R/generics.R:2175-2176@74e2e050]], this plan and the
memo, plus [[inst/NEWS.Rd:2394@74e2e050]], which is the 0.9-x-era history entry that introduced `"ev"` and stays (the TODO hit this
sentence once cited is gone; TODO carries no "post-mean" text now). No test and no Rd exercises
`value =` or `"post-mean"`.

One thing the memo does not say and the implementer must know: the shims are in the RELEASED 0.9.34 (checked against the
installed build - `predict.rbart`'s body there carries both, and its formals are `(object, newdata, group.by, offset, type,
combineChains, ...)`). "Nothing released to deprecate from" is true of 1.0-0 and false of 0.9-x, so the deletion is a 0.9-x-visible
removal and belongs in NEWS's UPGRADING subsection with D1's own positional change, not silently in NEW FEATURES.

## 9. D2: the saved-tree refusal

Every site that refuses a fit without stored trees, R/generics.R unless noted:
- [[inst/NEWS.Rd:277-283@74e2e050]] predict.bart, two arms: `"predict requires bart2 to be called with 'keepTrees' == TRUE"` and
  `"predict requires bart to be called with 'keeptrees' == TRUE"`, selected by `callName(object$call) == "bart2"`
  (`callName`, [[R/utility.R:187-189@74e2e050]]).
- [[R/generics.R:296-306@b46add06]] predict.bart's amplitude-coupled arm: `"predict on an amplitude-coupled fit requires 'keeptrees'/'keepTrees' == TRUE:
  ..."`.
- [[R/generics.R:1199-1204@b46add06]] predict.bartMultinomial, [[R/generics.R:1489-1494@b46add06]] predict.bartOrdinal, [[R/generics.R:1728-1733@b46add06]] predict.bartNegbin, [[R/generics.R:2123-2128@b46add06]]
  predict.bartHurdle: `"predict requires bart2(family = \"...\") to be called with 'keepTrees' == TRUE"`.
- [[R/generics.R:2162-2164@b46add06]] predict.rbart: `"predict requires rbart to be called with 'keepTrees' == TRUE"` - which also names the wrong
  function, since the entry point is `rbart_vi`.
On the same path and rewritten with them, because each is the same fact with the same cure:
- [[R/generics.R:445-455@b46add06]] extract.bart type = "trees", both arms; [[R/generics.R:2405-2409@b46add06]] extract.rbart type = "trees".
- [[R/generics.R:2648-2652@b46add06]] plotTree.bart (`"plotTree requires the trees to be kept: fit with keeptrees/keepTrees = TRUE"`) and [[R/generics.R:2672-2676@b46add06]]
  plotTree.rbart (`"... fit rbart_vi with keepTrees = TRUE"`). These two already name the cure in `= TRUE` form but on a third
  stem; folding them in is what leaves ONE message form for the fact across the whole surface, and plotTree.bart stops offering
  both spellings at once in favour of the one its own fit used.
Deliberately NOT rewritten: the `keepTrainingFits` refusals ([[R/generics.R:490@b46add06]], [[R/generics.R:2484@b46add06]], [[R/generics.R:2563@b46add06]]; [[R/plot.R:61@74e2e050]], [[R/plot.R:130@74e2e050]]) name a different argument and
a different fact, and the general error-message pass is lens 1 B9 / slice L's business. `[[R/bart.R:2415-2420@74e2e050]]`
(`hazardSurvivalProbabilities`) already says "requires the trees; refit with keepTrees = TRUE" and is the wording model; it stays
as written. `[[R/partialDependence.R:67@74e2e050]]`, [[R/partialDependence.R:73@74e2e050]] refuse on `keepSampler`, a different argument again.
`[[src/R_interface_bartcore.cpp:4127-4131@74e2e050]]`'s multinomial `getFitsWithoutOffset` refusal mentions `keepTrees` only in a trailing
caveat about what predict reports; it is not a refusal ON the tree store and stays.

Two helpers beside `predictOffsetUnusedArgs`:

    # predict, extract(type = "trees") and plotTree all read the fit's SAVED trees,
    # so a fit kept without them has nothing to read. The message names the one
    # argument that keeps them rather than restating the condition.
    refuseWithoutTrees <- function(what, keepTrees = "keepTrees") {
      stop(what, " requires the fit's saved trees; refit with ", keepTrees, " = TRUE")
    }

    # bart spells it 'keeptrees', bart2 and rbart_vi 'keepTrees'. A fit kept with
    # keepCall = FALSE stores call("NULL") and names neither, so it takes bart's
    # spelling, which is the surface such a fit most likely came from.
    bartKeepTreesArgument <- function(object) {
      if (callName(object[["call"]]) == "bart2") "keepTrees" else "keeptrees"
    }

Resulting text, exactly: `predict requires the fit's saved trees; refit with keepTrees = TRUE` (and the `keeptrees` twin);
`extract(type = "trees") requires the fit's saved trees; refit with keeptrees = TRUE`; `plotTree requires the fit's saved trees;
refit with keeptrees = TRUE`. The amplitude arm keeps its own reason, on the new stem:

    stop("predict requires the fit's saved trees; refit with ",
         bartKeepTreesArgument(object), " = TRUE: an amplitude-coupled fit pairs ",
         "each saved draw's forests with that draw's own amplitudes, and without ",
         "the tree store only the current trees replay, one set for every draw")

The four own-class arms drop the `bart2(family = "...")` naming: every one of those classes is reachable only through `bart2`, so
the spelling is fixed, and the family is not in doubt to whoever holds the fit. `predict.rbart`'s stops naming `rbart`.

A test that pinned the OLD `'keepTrees' == TRUE` text would break; in-repo none does (every pin is the bare word, section 11).
Downstream, bartCause's `tests/testthat/test-08-predict.R` line 208 pins `"keepTrees == TRUE"`, but against bartCause's OWN message
(bartCause's `R/generics.R` line 110), raised before the call ever reaches dbarts - unaffected.

## 10. Rd plan

man/bart.Rd. [[R/generics.R:45-53@74e2e050]], predict.bart's usage, after:

    \method{predict}{bart}(
        object, newdata,
        type = c("ev", "ppd", "bart", "forest"),
        offset = NULL, weights = NULL,
        combineChains = TRUE,
        ci.level = NULL,
        forest = NULL,
        bases = NULL,
        n.threads,
        \dots)

:187-189 `\item{offset}` gains the sentence that it is the shift at the PREDICTED rows and that the fit-time channel is
`offset.test`, so the two names are visibly different things rather than a typo. [[R/generics.R:202-204@74e2e050]] `\item{type}` gains the `"class"`
sentence for the two categorical families (it is the shared type item across the bart family pages) and states the four
vocabulary exceptions of section 6. [[R/generics.R:208-210@74e2e050]] `\item{forest}` is untouched.

man/rbart.Rd. [[R/generics.R:56-62@74e2e050]], after (the file uses `\S3method`, keep it):

    \S3method{predict}{rbart}(
        object, newdata,
        type = c("ev", "ppd", "bart", "ranef"),
        offset = NULL, weights = NULL,
        combineChains = TRUE,
        ci.level = NULL,
        n.threads,
        \dots, group.by)

:67-69 `\item{group.by}` gains: "For \code{predict} and \code{\link{survivalProbabilities}}, supplied by name only - it follows
\code{\dots} in the signature, so it is never matched positionally; a missing one is refused, naming itself." The shared
catch-all item at [[R/generics.R:85-87@74e2e050]] already covers `offset`, `offset.test` and `n.threads` and needs no edit.

man/survivalProbabilities.Rd. [[R/generics.R:35-42@74e2e050]], after:

    \method{survivalProbabilities}{rbart}(
      object,
      times,
      newdata = NULL,
      combineChains = TRUE,
      \dots,
      group.by
    )

with [[R/generics.R:71@74e2e050]]'s `\item{group.by}` gaining the same named-only sentence. [[R/generics.R:27-33@74e2e050]]'s `bart` entry and the four own-class entries are
unchanged.

man/bart2.Rd. [[R/generics.R:82-86@74e2e050]] predict.bartMultinomial gains `"class"` in its type vector (commit 2, with the formal). Fifteen new
`\method` entries (the file's own spelling) inserted so each class's five sit together, ordered as the aliases are:

    \method{extract}{bartOrdinal}(
        object, type = c("ev", "ppd", "bart", "loglik"),
        sample = c("train", "test"),
        combineChains = TRUE, \dots)

    \method{fitted}{bartOrdinal}(
        object, type = c("ev", "class", "bart"),
        ci.level = NULL, \dots)

    \method{predict}{bartOrdinal}(
        object, newdata,
        type = c("ev", "ppd", "bart", "class"),
        combineChains = TRUE, ci.level = NULL, n.threads, \dots)

    \method{print}{bartOrdinal}(x, \dots)

    \method{residuals}{bartOrdinal}(object, \dots)

    \method{extract}{bartNegbin}(
        object, type = c("ev", "ppd", "bart", "loglik"),
        sample = c("train", "test"),
        combineChains = TRUE, \dots)

    \method{fitted}{bartNegbin}(
        object, type = c("ev", "ppd", "bart"),
        ci.level = NULL, \dots)

    \method{predict}{bartNegbin}(
        object, newdata,
        type = c("ev", "ppd", "bart"),
        offset = NULL,
        combineChains = TRUE, ci.level = NULL, n.threads, \dots)

    \method{print}{bartNegbin}(x, \dots)

    \method{residuals}{bartNegbin}(object, \dots)

    \method{extract}{bartHurdle}(
        object, type = c("ev", "ppd", "prob", "bart", "loglik"),
        sample = c("train", "test"),
        combineChains = TRUE, \dots)

    \method{fitted}{bartHurdle}(
        object, type = c("ev", "ppd", "prob", "bart"),
        sample = "train", ci.level = NULL, \dots)

    \method{predict}{bartHurdle}(
        object, newdata,
        type = c("ev", "ppd", "prob", "bart"),
        combineChains = TRUE, ci.level = NULL, n.threads, \dots)

    \method{print}{bartHurdle}(x, \dots)

    \method{residuals}{bartHurdle}(object, type = "ev", \dots)

Every argument these introduce already has an `\item`: object [[R/generics.R:278@74e2e050]], newdata [[R/generics.R:281@74e2e050]], type [[R/generics.R:284@74e2e050]], sample [[R/generics.R:287@74e2e050]], ci.level [[R/generics.R:302@74e2e050]],
combineChains [[R/generics.R:178@74e2e050]], n.threads [[R/generics.R:173@74e2e050]], offset [[R/generics.R:126@74e2e050]], x [[R/generics.R:293@74e2e050]], `\dots` [[R/generics.R:275@74e2e050]] - so `checkDocFiles`'s undocumented-argument test is
satisfied without a new item. Prose edits in the same file: [[R/generics.R:338@74e2e050]], which says predict.bartNegbin takes "an optional log-exposure
\code{offset.test}", becomes `offset` (commit 1) and then gains fitted's `"ppd"` (commit 2); [[R/generics.R:330@74e2e050]] and [[R/generics.R:334@74e2e050]] gain the
`type = "class"` sentence for predict (commit 2). [[R/generics.R:126@74e2e050]] `\item{offset}` gains a third paragraph for the negbin log-exposure shape
(it already carries the multinomial matrix shape in its second). [[R/generics.R:131@74e2e050]] `\item{offset.test}` stays - it documents bart2's own
fit-time argument.

man/dbartsSampler-class.Rd. Usage [[R/generics.R:94@74e2e050]], [[R/generics.R:96@74e2e050]], [[R/generics.R:97@74e2e050]] become `(forest = NULL)`; [[R/generics.R:95@74e2e050]] and [[R/generics.R:98-100@74e2e050]] unchanged. `\item{forest}` [[R/generics.R:201@74e2e050]] is
rewritten: all four readers default to `NULL`, every forest, with the four stacked shapes named one by one (section 5) and the
statement that a single-forest sampler's `NULL` read is exactly its forest-1 read; `setForestWeights`/`setForestBasis` keep no
default (they are writers naming one target) and `setCalibration` keeps `1L` with the map refusal as its reason. `\value` [[R/generics.R:428@74e2e050]]
restates getForestFits's and getForestVariableCounts's stacked shapes beside the amplitudes shape it already carries; [[R/generics.R:432@74e2e050]]
restates getCalibration's. [[R/generics.R:210@74e2e050]]'s remedy text `setOffset(rep_len(-getCalibration()[1, "prior.mean"], n))` stays correct as
written - single forest, matrix result.

[[man/plotTree.Rd:48@74e2e050]] ("trees kept (\code{keeptrees}/\code{keepTrees} equal to \code{TRUE})") is prose about the requirement, not a
quotation of the message, and stays.

`tools/check-rc-codoc.R` parses the generator's `methods = list(...)` against these `\S4method` usage entries and compares names,
order and defaults, so [[man/plotTree.Rd:94@74e2e050]]/[[man/plotTree.Rd:96@74e2e050]]/[[man/plotTree.Rd:97@74e2e050]] must move in the SAME commit as [[R/dbarts.R:1691@74e2e050]], [[R/dbarts.R:1716@74e2e050]], [[R/dbarts.R:1730@74e2e050]]. `R CMD check`'s codoc covers
the S3 side the same way, which is why the `"class"`/`"ppd"` vocabulary edits and their `\usage` lines are one commit and D1's
reorder another. Keep each usage entry's defaults spelled exactly as the formals spell them. `n.threads` is shown with no default
in the three existing predict usage entries ([[man/bart.Rd:52@74e2e050]], [[man/bart2.Rd:86@74e2e050]], [[man/rbart.Rd:61@74e2e050]]) and passes today; the three new
predict entries follow that spelling for consistency within the file rather than introducing a fourth style. Spelling the real
defaults at all six sites instead is a free rider that also closes lens 1 C8 for these entries - implementer's call, either way
the six must agree.

## 11. Test plan

Files that MUST change, with the reason:
- [[inst/tinytest/test-rbart-generics.R:162@74e2e050]], [[inst/tinytest/test-rbart-generics.R:180@74e2e050]], [[inst/tinytest/test-rbart-generics.R:199@74e2e050]], [[inst/tinytest/test-rbart-generics.R:201@74e2e050]] (two calls), [[inst/tinytest/test-rbart-generics.R:203@74e2e050]], [[inst/tinytest/test-rbart-generics.R:204@74e2e050]] - positional `group.by`, 7 sites.
- [[inst/tinytest/test-rbart-groupby.R:183@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:187@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:195@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:231@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:235@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:241@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:277@74e2e050]], [[inst/tinytest/test-rbart-groupby.R:282@74e2e050]] - 8 sites ([[inst/tinytest/test-rbart-groupby.R:241@74e2e050]] and [[inst/tinytest/test-rbart-groupby.R:282@74e2e050]] are the
  `suppressWarnings(predict(` blocks, whose third argument sits on its own line).
- [[inst/tinytest/test-rbart-bartcore.R:68@74e2e050]], [[inst/tinytest/test-rbart-bartcore.R:73@74e2e050]] - 2 sites.
- [[inst/tinytest/test-generics-multithreaded.R:261@74e2e050]], [[inst/tinytest/test-generics-multithreaded.R:262@74e2e050]], [[inst/tinytest/test-generics-multithreaded.R:271@74e2e050]] - 3 sites ([[inst/tinytest/test-generics-multithreaded.R:272@74e2e050]] already passes `group.by = g` and is the
  positional/named equality check, which is exactly the test that must keep passing).
- [[inst/tinytest/test-generics-posteriorPredictiveDistribution.R:131@74e2e050]] - 1 site.
  All 21 become `group.by = g`. This is the complete positional list: parsing every R file under R/, inst/tinytest/ and tests/
  for `predict`/`fitted`/`extract`/`residuals` calls with three or more unnamed arguments returns these 21, plus
  `extract(fit, "trees", "train")` ([[test-sampler-trees.R:77@74e2e050]], extract's own unchanged order), plus
  `predict(fitted, x, bases = fitted$bases, ...)` ([[test-predict-blend.R:83@74e2e050]], a helper forwarding `...`), plus the three internal
  `extract(object, type, sample, ...)` calls in R/generics.R. No `survivalProbabilities` call passes `group.by` positionally:
  [[test-rbart-aft.R:100@74e2e050]] and [[test-rbart-aft.R:112@74e2e050]] already name it.
- [[inst/tinytest/test-nbinom.R:110@74e2e050]] - `offset.test = rep(log(2), 10L)` becomes `offset =`.
- [[inst/tinytest/test-predict-blend.R:376@74e2e050]], [[inst/tinytest/test-predict-blend.R:380@74e2e050]] - `pattern = "requires 'keeptrees'/'keepTrees'"` no longer matches; use
  `pattern = "saved trees"`.
- [[inst/tinytest/test-plot-generics.R:110-113@74e2e050]] - `pattern = "requires the trees to be kept"` no longer matches; use
  `pattern = "saved trees"`.
Files that need NO change though they look like they might, every one pinning the bare word rather than the sentence:
[[test-generics-errors.R:21@74e2e050]], [[test-generics-errors.R:23@74e2e050]] (`"keeptrees"`), [[test-generics-errors.R:39@74e2e050]], [[test-generics-errors.R:41@74e2e050]] (`"keepTrees"`); [[test-nbinom.R:125@74e2e050]], [[test-ordinal.R:121@74e2e050]],
[[test-hurdle.R:188@74e2e050]], [[test-multinomial-surface.R:961@74e2e050]] (`"keepTrees"`); [[test-hazard.R:234@74e2e050]], whose target message
([[R/bart.R:2415-2420@74e2e050]]) is not rewritten; [[test-fits-without-offset.R:241@74e2e050]], whose `"keepTrees"` is in the multinomial
`getFitsWithoutOffset` refusal, a different message that stays. Type-vocabulary pins are likewise safe: [[test-hurdle.R:241-242@74e2e050]],
[[test-multinomial-generics.R:265-268@74e2e050]], [[test-nbinom.R:410@74e2e050]], [[test-nbinom.R:440-441@74e2e050]], [[test-ordinal.R:420@74e2e050]], [[test-ordinal.R:444-445@74e2e050]] and [[test-pointwise-loglik.R:38-39@74e2e050]]
all probe types that stay refused.

New tests:
- test-generics-errors.R, one block per class: the `offset.test` spelling refused by name on all six (it already checks bart at
  [[test-pointwise-loglik.R:62-64@74e2e050]] and rbart at [[test-pointwise-loglik.R:80-88@74e2e050]] - extend to the four own-class fits, and check that ordinal's and hurdle's message says
  "no out-of-sample offset channel" rather than "named 'offset'"); `weights` refused on the four; `offset` refused on ordinal and
  hurdle; `value` refused on rbart, naming `type`; and, on rbart, `predict(fit, x, g)` raising the "given by name only" message
  and `predict(fit, x, type = "ev")` with no `group.by` raising it too.
- test-generics-errors.R: each of the six saved-tree refusals matches `"refit with keepTrees = TRUE"` (or `keeptrees` on the
  `bart` arm), and `extract(fit, type = "trees")` and `plotTree(fit)` match the same stem.
- test-rbart-generics.R: `predict(fit, x, group.by = g)` equals the old positional answer bitwise, so the reorder is proven inert.
- test-rbart-aft.R: `survivalProbabilities(fit, times, newdata = x.new, g)` now raises the named-only message rather than
  binding `g` to `combineChains`.
- test-calibration-midchain.R: on a single-forest sampler `getCalibration()` is `identical` to `getCalibration(1L)`; on the BCF
  sampler it already builds ([[test-pointwise-loglik.R:362-380@74e2e050]]) `getCalibration()` is a numChains x 12 x 2 array whose `[, , 1]` and `[, , 2]` are the two
  indexed reads, with the `leaf.model` attribute and the column dimnames preserved.
- test-bcf-family.R or test-fits-without-offset.R: the same three assertions for `getForestFits()` and
  `getForestVariableCounts()` - single-forest NULL identical to `1L`, BCF NULL an n x 2 x numChains (resp. p x 2 x numChains)
  array whose `[, f, ]` slices are the indexed reads, predictor rownames intact on margin 1 in both shapes.
- test-multinomial-r5-surface.R: `getForestFits()` on the K-forest sampler is n x K x numChains.
- test-multinomial-generics.R and test-ordinal.R: `predict(fit, newdata, type = "class")` is a factor over the fit's levels,
  `nrow(newdata)` long, ordered for ordinal, and equals `fitted(fit, type = "class")` when `newdata` is the training x; and
  `predict(fit, newdata, type = "class", ci.level = 0.9)` returns the SAME band `type = "ev"` does, which is the assertion that
  fails if the class arm is placed above the ci.level block.
- test-nbinom.R and test-hurdle.R: `fitted(fit, type = "ppd")` runs, is n long, and sits near `fitted(fit, type = "ev")` at a
  loose tolerance (it is a Monte Carlo estimate of the same mean).

RNG: no test's stream moves. Nothing in this slice touches the sampler, the predict replay, or the order in which any existing
arm draws. The two new draw paths - `fitted.bartNegbin(type = "ppd")` and `fitted.bartHurdle(type = "ppd")` - are reachable only
by asking for them, and the new tests that do must sit at the END of their files or carry their own `set.seed`, since several
regression tests in those files hardcode values that depend on the file's full execution history.

## 12. Consumer sweep (read-only, `git -C <repo> grep`)

- stan4bart, /Users/vdorie/Repositories/stan4bart branch `bartcore` (33b4aa8): ZERO hits. It calls no dbarts S3 predict/fitted/
  extract method and no forest reader; its own `predict.stan4bartFit` (stan4bart's `R/generics.R` lines 906-914) already has the
  `(object, newdata, type, ...)` shape and its comment at `R/generics.R` lines 907-908 records why offset must stay behind the defaulted
  arguments - independent confirmation of D1's ordering. Its `offset.test` hits (stan4bart's `R/mvbart.R` lines 38, 135, 150 and 195;
  `src/init.cpp` lines 327 and 329) are all the dbartsData/creation channel, which keeps the name. No migration cost.
- bartCause, /Users/vdorie/Repositories/bartCause branch `dbarts-1.0` (7ae6e83): ONE breaking site.
  bartCause's `R/generics.R` line 162 `p.score <- predict(object$fit.trt, x.new, group.by, combineChains = FALSE, ...)` on an rbart treatment fit
  becomes `group.by = group.by`. Its other dbarts reads are safe: `sampler$getCalibration(1L)` (bartCause's `R/bcf.R` line 244) and
  `fit$fit$getForestVariableCounts(1L|2L)` (bartCause's `tests/testthat/test-14-bcf.R` lines 162-164) all pass the index explicitly, so the new NULL
  default never fires; every `predict(fit, x, type = ...)` in bartCause's `R/generics.R` lines 143-176 and `tests/testthat/test-08-predict.R` passes
  `type` and `group.by` by name; and its `"keepTrees == TRUE"` pin at `tests/testthat/test-08-predict.R` line 208 is against its own
  message (bartCause's `R/generics.R` line 110), not dbarts's. Separately - not a cost, an observation - bartCause's own `predict.bartcFit`
  (bartCause's `R/generics.R` lines 87-91) has the same `(object, newdata, group.by, type, ...)` shape D1 removes from dbarts; whether it follows is
  its maintainer's call and not this slice's.
- treatSens, branch `dbarts-1.0` (1db3d89): ZERO hits. No
  `offset.test`, no `predict(`, no forest reader, no keepTrees text - it is a flat-C consumer only.
- bairrtt, /Users/vdorie/Repositories/bairrtt branch `main` (6167423): ZERO breaking hits. Its six predict calls
  (bairrtt's `R/irt_causal_bart.R` lines 568, 571, 636, 637, 709 and 713) are all `model$predict(frame)` on the reference class, whose signature
  this slice does not touch.
Total lockstep migration cost: one line, in bartCause.

## 13. Settled sub-choices

Five, all settled by orchestrator ruling under the standing grant, 2026-08-25, and written into the sections above rather than
left open. (1) Ordinal and hurdle refuse `offset`/`offset.test` by name and gain NO dummy formal, so the positional prefix
guarantee is the first three arguments and no further (section 3). (2) Those two classes carry their OWN `offset.test` wording
rather than composing `predictOffsetUnusedArgs`, since `refuseUnusedGenericArgs` reports the first name in `names(reasons)` and
the composed list would point at an argument they do not have (section 3). (3) D5 deletes both shims AND refuses `value` by name;
a refusal is not a shim (section 8). (4) `fitted` does not gain `"forest"`: a forest-margined value has no per-observation
posterior mean, and `fitted.bart` [[R/generics.R:877-881@b46add06]] reduces the last margin, which for that value is the forest (section 6). (5) D9's
`forest = NULL` is implemented by ONE new bridge entry (`bartcore_numForests`) plus R-side stacking, not by teaching three
existing bridge entries a NULL branch; the three per-forest readers are not edited at all (section 5).

## 14. Commit plan and gates

Four commits, in the plan's order. D1 first because D9's fifteen new usage entries would otherwise have to be written twice, and
D5 and D2 both land inside functions D1 has already reshaped. Every commit is codoc-clean on its own: a formal and its `\usage`
line always move together, which is why the `"class"` token is commit 2's and not commit 1's even though it lands in a signature
commit 1 also edits.

1. D1, `predict-signature-unification`. R/generics.R (three predict signatures actually move - bart, rbart and negbin's rename;
   multinomial, ordinal and hurdle already carry the target prefix and change only in the body - plus two `missing()` blocks
   deleted, the `group.by` by-name refusal, the negbin rename's one use, the [[R/generics.R:227-233@b46add06]] comment, `predictOffsetUnusedArgs`'s
   comment, two new reason lists, four call-site list compositions); [[R/bart.R:2551-2557@74e2e050]] and [[R/bart.R:2576-2578@74e2e050]]
   (survivalProbabilities.rbart); [[man/bart.Rd:45-53@74e2e050]], [[man/bart.Rd:187-189@74e2e050]]; [[man/rbart.Rd:56-62@74e2e050]], [[man/rbart.Rd:67-69@74e2e050]];
   [[man/survivalProbabilities.Rd:35-42@74e2e050]], [[man/survivalProbabilities.Rd:71@74e2e050]]; [[man/bart2.Rd:338@74e2e050]]; the 22 test call sites in six files; one NEWS UPGRADING item.
2. D9, `surface-smalls`. src/R_interface_bartcore.cpp (one new entry), src/R_interface_bartcore.hpp, src/R_interface.cpp
   (declaration and registration); R/bartcore.R (`bartcoreNumForests`); [[R/dbarts.R:1691@74e2e050]], [[R/dbarts.R:1716@74e2e050]], [[R/dbarts.R:1730@74e2e050]] and their docstrings;
   R/generics.R (four choice vectors, the negbin fitted switch, the two shared category helpers, two predict class arms, two
   `trailing` selections); [[man/dbartsSampler-class.Rd:94@74e2e050]], [[man/dbartsSampler-class.Rd:96@74e2e050]], [[man/dbartsSampler-class.Rd:97@74e2e050]], [[man/dbartsSampler-class.Rd:201@74e2e050]], [[man/dbartsSampler-class.Rd:428@74e2e050]], [[man/dbartsSampler-class.Rd:432@74e2e050]]; [[man/bart2.Rd:82-86@74e2e050]], [[man/bart2.Rd:330@74e2e050]], [[man/bart2.Rd:334@74e2e050]], [[man/bart2.Rd:338@74e2e050]] and
   the fifteen new usage entries; the new tests. This commit is the only one that compiles: `R CMD INSTALL .` suffices - no
   header, no facade virtual moves, so `--preclean` is not required, though it costs nothing.
3. D5, `deprecation-shim-removal`. [[R/generics.R:1826@74e2e050]], [[R/generics.R:2167-2180@74e2e050]] plus the new `value` reasons list; one NEWS UPGRADING item.
4. D2, `predict-refusal-names-cure`. R/generics.R (two helpers; [[R/generics.R:277-283@74e2e050]], [[R/generics.R:296-306@74e2e050]], [[R/generics.R:445-455@74e2e050]], [[R/generics.R:1199-1204@74e2e050]], [[R/generics.R:1489-1494@74e2e050]], [[R/generics.R:1728-1733@74e2e050]],
   [[R/generics.R:2123-2128@74e2e050]], [[R/generics.R:2162-2164@74e2e050]], [[R/generics.R:2405-2409@74e2e050]], [[R/generics.R:2648-2652@74e2e050]], [[R/generics.R:2672-2676@74e2e050]]); [[inst/tinytest/test-predict-blend.R:376@74e2e050]], [[inst/tinytest/test-predict-blend.R:380@74e2e050]];
   [[inst/tinytest/test-plot-generics.R:110-113@74e2e050]]; the new refusal tests; one NEWS item.

Re-anchoring, in the SAME commit as its own code edits, exactly as the dbarts.h freeze slice did: docs/design carries 21 anchors
into R/generics.R, 41 into R/dbarts.R, 27 into src/R_interface_bartcore.cpp and 7 into the man pages this slice moves. Run
`Rscript tools/check-doc-freshness.R .` after each commit's edits and re-align every strict miss from the `git diff -U0` line map,
editing the docs/design anchors in place so each file's line count is invariant. Baseline at this tip is 0 FAIL / 71 WARN.

Gate battery per commit (CLAUDE.local.md): `R CMD INSTALL .`; `tinytest::test_package("dbarts")`; `tools/check-rc-codoc.R`
(commit 2 especially); `tools/check-doc-freshness.R`; `air format --check .` and lintr; `R CMD check --as-cran`. Commits 1, 3 and
4 are R and Rd only, so `cd tests/cpp && make && ./test_bartcore` need only run on commit 2, where it must stay at its current
268 ok. Equivalence: `benchmarks/R/equivalence.R compare` against equivalence-736bfb05.rds, bcf-equivalence-6e3b9fb8.rds and
multinomial-equivalence-4d9a3337.rds, ALL BITWISE IDENTICAL on every commit - the harness calls predict only with named
arguments ([[equivalence.R:1229-1232@578ee440]]) and `getCalibration()` only on single-forest samplers ([[backfit-exact.R:143@578ee440]], [[geweke-mc.R:528@578ee440]],
:531, [[geweke-mc.R:562@578ee440]]), so it needs no edit. `bench-sampler` is not required: no sampling path is touched and the slice adds one scalar
bridge read that no benchmark calls in a loop; run it once at the end of the slice on a quiet machine against
bench-sampler-ab1dc52.csv if the RC tip wants a clean sheet. Zero baseline re-records expected; if any equivalence channel moves,
something in this slice reached the sampler and the commit is wrong.

Consumer gate, after commit 1: rebuild bartCause on `dbarts-1.0` with the one-line fix and run `testthat::test_local()`;
stan4bart, treatSens and bairrtt need no rebuild for this slice (no hits), though the branch's habit of rebuilding both flat-C
consumers costs little.

## Landing note (2026-08-25)

Landed as five commits, pushed together as 78f334c1: 7b3ac6bf (D1, commit 1), 71cc7133 (D9, commit 2), befc8f45 (D5, commit 3),
ed43deef (D2, commit 4), 78f334c1 (docs/design anchor re-alignment by the 716acd1f..ed43deef diff line map; stamps in
feature-matrix.md and threaded-predict.md moved to ed43deef). Implementer gates ran per commit; the independent battery re-ran on a
git-archive snapshot of 78f334c1: tinytest 7352/0, tests/cpp + sampler pass, equivalence 43 identical (equivalence.R
--strict-coverage), bcf 12 and multinomial 11 identical (their OWN sibling scripts bcf-equivalence.R/multinomial-equivalence.R -
section 12's invocation naming equivalence.R for all three is wrong), check-rc-codoc 42 methods, freshness 0 FAIL / 68 WARN,
NEWS.Rd 339 entries (three new UPGRADING items: D1, D5, D2), air + lintr::lint_package clean, R CMD check --as-cran 1 NOTE
(days-since-update), census greps zero. All baselines bitwise; zero re-records.

Deviations from the design, all content-neutral: the predict.rbart formals-order loop test in test-generics-multithreaded.R split
so rbart checks (n.threads, ..., group.by); man/bart2.Rd's class/ppd prose landed in the "Generics for a X fit" paragraphs the
design's line numbers actually named; commit 1's bart2.Rd edit was prose-only (the three own-class predicts already carried the
target prefix).

Consumers: bartCause dbarts-1.0 takes TWO named group.by fixes in R/generics.R (d825cfc), not the one line both sweeps counted:
the direct predict call at [[geweke-mc.R:162@78f334c1]] and a do.call args list at [[geweke-mc.R:193@78f334c1]] that neither parse-walk keyed on (it builds list(fit, x, group.by,
...) rather than calling predict directly) - caught by running bartCause's suite against the landed build; stan4bart, treatSens,
bairrtt unaffected (zero hits, re-verified by the critique). Residue: [[ordinal.md:364@78f334c1]] and [[survival.md:615@78f334c1]] cite
[[R/generics.R:269-275@78f334c1]], whose tail fell in a deletion - re-anchoring them needs intent, not arithmetic, so both stand as advisories;
[[multinomial-mutation-arc.md:835@78f334c1]]'s [[bart2.Rd:291@78f334c1]] cite drifted to 348 but sits in frozen section 5, left by rule;
[[model-space-survey.md:429@78f334c1]] was stale before the slice (frozen exempt). The two settled-sub-choice doors that remain open:
extract gains no "class" (a reduction over draws, not a channel) and fitted gains no "forest" (the last-margin reduction would
average over forests, not observations).
