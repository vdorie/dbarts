# Matrix review 2: fit classes x generics, and the sampler class

Branch `bartcore`, tip `b102e17c`. Read-only. Private lib built from `git archive HEAD`
(`R CMD INSTALL --preclean`); every cell below ran against that lib, never the checkout.
Scripts and raw grids: `generics/` beside this file.

## 0. Reproduction of the builder's grid

- `matrix.R` + `matrix-sampler-census.R` re-run against the fresh lib reproduce
  `matrix-grid.csv` BYTE-IDENTICALLY (877 rows, `diff` empty); F1-F9 stand on
  re-execution, and F3 and F6 are confirmed and sharpened below.

## 1. Cells added by this pass

- `generics/fill-grid.csv`: 4110 executed cells + 316 census/build rows over 40 fit
  variants (13 families x {base, keepTrees+keepSampler FALSE, n.chains 2 combined,
  n.chains 2 uncombined}) x 11 generics x every `type=`/`sample=` the method declares
  x `combineChains` TRUE/FALSE, run twice: live, and after `saveRDS`/`readRDS` in a
  FRESH process. Live 949 accepted / 1308 refused / 201 raw; rds 631 / 883 / 138.
- `generics/sampler-census4.csv`: 192 cells - 48 RC methods on each of a multinomial,
  an amplitude (two-forest), a heteroscedastic and a weighted-logistic sampler, reads
  before mutators so the builder's tree-store ordering artifact does not recur.
  143 accepted / 48 refused / 1 raw. Plus ~270 targeted probes
  (`generics/probe*.R`, `census4b.R`) for the Rd claims, the partial-matching class,
  and the argument-swallow class.

## 2. FINDINGS

### BLOCKER

**G1. `survivalProbabilities()` is dead on every discrete-time hazard fit whose
training predictor matrix carries column names.**
- `[[R/bart.R:2421-2424@b102e17c]]` (training branch) and `[[R/bart.R:2432-2435@b102e17c]]` (matrix-`newdata` branch)
  build the re-expanded design as `cbind(<named covariates>, rep(seq_len(K), each = n))`,
  leaving the period column UNNAMED. `predict` then matches by name against the fit's
  `x1,x2,period` and refuses. The data-frame branch at `[[R/bart.R:2427-2428@b102e17c]]` is correct
  (`bigX[["period"]] <- ...`), which is why only two of three branches are broken.
- Repro (`generics/probe15.R`):
  ```r
  x <- matrix(runif(120), 60, 2); colnames(x) <- c("x1","x2")
  fit <- bart2(x, Surv(t, s), family = "hazard", keepTrees = TRUE, ...)
  survivalProbabilities(fit)
  # Error: column names of 'test' do not match those of 'x': 'period' present in
  #        'x' but not in 'test' (whose columns are ...)
  ```
  Identical for `times=`, `newdata=<named matrix>` and `combineChains = FALSE`. Drop
  the `colnames` and all of it works (`6x60x60`, `6x60x4`); a `data.frame` `newdata`
  also works.
- Contradicts [[survivalProbabilities.Rd:11-13@b102e17c]], [[survivalProbabilities.Rd:45-47@b102e17c]], [[survivalProbabilities.Rd:87-96@b102e17c]] and [[bart.Rd:328@b102e17c]]. Hazard
  fits are matrix-interface only (the formula route is refused by name) and [[bart.Rd:75@b102e17c]]
  tells users column names are honoured, so the named matrix is the normal spelling.
- Why no test caught it: `[[inst/tinytest/test-hazard.R:12@b102e17c]]` builds
  `x <- matrix(runif(n * p), n, p)` with NO column names, so the file exercises only
  the working branch.
- Fix class: agent-fix (name the column `period` in both branches).

**G2. `extract(fit, type = "trees", <any own formal other than object/type>)` corrupts
its own call.** Confirms and broadens matrix-results F3.
- `[[R/generics.R:375-380@b102e17c]]` (`extract.bart`) rewrites `match.call()` onto
  `object$fit$getTrees` and strips only `object` and `type`, so the four surviving
  formals - `sample`, `combineChains`, `forest`, `contribution` - are forwarded into
  `getTrees(treeNums, chainNums, sampleNums, current, newdata)` ([[R/dbarts.R:1917@b102e17c]]).
  `extract.rbart` ([[R/dbarts.R:1857-1866@b102e17c]]) repeats it for `sample` and `combineChains`.
- Repro (`generics/probe3.R`), all on a `keepTrees = TRUE` fit:
  `extract(f, type = "trees", sample = "train")` and the positional
  `extract(f, "trees", "train")` -> `missing value where TRUE/FALSE needed`
  (`sample` partial-matches `sampleNums`, `as.numeric("train")` is NA, a
  `sampleNums <= 0` guard dies on it); `combineChains = FALSE`, `forest = 1L` and
  `contribution = TRUE` each -> `unused argument (...)`. 33 grid cells in
  matrix-grid.csv, plus 102 more from the `combineChains` axis this pass added.
- Contradicts [[bart.Rd:203-205@b102e17c]] (`sample` documented with no type qualification) and
  [[bart.Rd:219@b102e17c]] (only `\dots` is documented as reaching `getTrees`). Correction to
  matrix-results F3's last paragraph: with `sample` OMITTED nothing is forwarded and
  nothing is silently ignored - the bug is only the supplied-argument one.
- Fix class: agent-fix (strip the method's own formals from `treesCall`, or refuse
  them by name under `type = "trees"`).

### MAJOR

**G3. `extract(..., combineChains = FALSE)` is silently ignored on `bartMultinomial`,
`bartOrdinal` and `bartNegbin`, while `predict` on the same fits honours it.**
- `extract.bartMultinomial` ([[R/generics.R:862@b102e17c]]), `.bartOrdinal` ([[R/generics.R:1071@b102e17c]]) and
  `.bartNegbin` ([[R/generics.R:1237@b102e17c]]) have formals `(object, type, sample, ...)`, so it lands in
  `...` and is dropped; `extract.bart`, `.rbart` and `.bartHurdle` all carry it.
- Repro (`generics/probe11.R`, 2-chain fits): `extract(fO, combineChains = FALSE)`
  returns `12x60x3`, `identical()` to the combined result, where the matching
  `predict(fO, xTest, combineChains = FALSE)` returns `2x6x8x3`.
- Fix class: VD-judgement (group A below).

**G4. Same swallow for `sample=` and `ci.level=` on `fitted`, `ci.level=` on
`predict`, `type=` on `residuals`, `forest=` on `extract` and `vars=` on `summary`,
on those three classes.** Each verified `identical()` to the call without the argument
(`generics/probe11.R`): `fitted(fM, sample = "test")`, `fitted(fO, ci.level = 0.9)`,
`predict(fN, xTest, ci.level = 0.9)`, `residuals(fM, type = "bart")`,
`extract(fM, forest = 1L)`, `summary(fM, vars = "sigma")`. bart/rbart/bartHurdle
honour every one. Fix class: VD-judgement (A).

**G5. `$n.chains` is absent from every fit that keeps a sampler, i.e. from every
`keepTrees = TRUE` fit.** Confirms the lenses memo.
- `[[R/bart.R:390-394@b102e17c]]`: `if (keepSampler) result$fit <- fit else result$n.chains <- n.chains`.
  `keepTrees/keepSampler` in {TT, TF, FT} all give `is.null(fit$n.chains)`; only FF
  sets it. [[bart.Rd:315-317@b102e17c]] documents it unconditionally, and its stated purpose
  ("information that can be lost if `combinechains` is TRUE") is exactly the missing
  case.
- Asymmetric: rbart, bartMultinomial, bartOrdinal and bartNegbin DO carry `n.chains`
  under `keepTrees = TRUE`; bart and bartHurdle do not. Fix class: agent-fix (set it
  unconditionally, as the siblings do).

**G6. `extract(fit, type = "trees")` on a `keepTrees = FALSE, keepSampler = TRUE` fit
silently returns the CURRENT working trees, in a different column set, instead of
refusing.**
- The guard at [[R/generics.R:363-373@b102e17c]] tests `is.null(object$fit)`, i.e. the SAMPLER, not
  `keepTrees`; `getTrees`'s `useSaved <- control@keepTrees && !current` then falls back.
  `bart2(..., keepTrees = FALSE, keepSampler = TRUE)` then
  `extract(f, type = "trees")` -> an `11x4` frame with columns `tree,n,var,value` - no
  `sample` column, no warning - where [[bart.Rd:257@b102e17c]] scopes the feature to
  `keepTrees = TRUE` and [[bart.Rd:260-266@b102e17c]] documents `chain, sample, tree, n, var, value`.
- Fix class: VD-judgement (group F).

**G7. No `plot` method for `bartOrdinal`, `bartNegbin`, `bartHurdle`** (NAMESPACE:
31-34,56). Confirms matrix-results F6 and names the default: `plot(ordinalFit)` reaches
`plot.default` and raises `'x' is a list, but does not have components 'x' and 'y'`.
bart, rbart and bartMultinomial all have one; all six have summary and print.
Fix class: VD-judgement (group B).

**G8. `extract(type = "loglik")` exists only on `bart` and `rbart`.** All four
own-class families refuse by name (`type must be in 'ev', 'ppd', 'bart'` etc.); this is
the channel loo/WAIC consume, and [[bart.Rd:201@b102e17c]] documents `"loglik"` on the `extract`
generic without scoping it. Fix class: VD-judgement (group C).

**G9. A saved `bartHurdle` fit cannot `predict` after reload, and no Rd names the
recipe that fixes it.**
- The hurdle fit has no `$fit` of its own (`names()` = call, family, y, occupancy,
  positive), so [[bart.Rd:248@b102e17c]]'s `bartFit$fit$storeState()` does not apply, and [[bart.Rd:248@b102e17c]]
  names only multinomial/ordinal/nbinom as the extended-family cases.
- Repro (`generics/probe6a.R`/`6b.R`): `saveRDS` then a fresh session ->
  `predict(fH, xTest)` errors "samplers cannot be re-created without a stored state".
  `fH$occupancy$fit$storeState(); fH$positive$fit$storeState()` first makes it work
  (`generics/probe7.R`). Every other class survives its documented recipe - bart,
  rbart, multinomial and the amplitude fit all replay
  `predict`/`extract(type="forest")`/`plotTree` after reload.
- Fix class: agent-fix (extend [[bart.Rd:248@b102e17c]]).

**G10. [[dbartsSampler-class.Rd:328@b102e17c]] "Reading it forces the sampler's *current* state" is
false, and the stale value is rejected by the sampler's own `setState`.**
- `run()` ASSIGNS `state` (via `storeState`), replacing the lazy promise, so a later
  mutation leaves `state` stale. Repro (`generics/probe8.R`): two-forest sampler,
  `run()`, then `setForestBasis(2L, <1-column basis>)` narrowing the block from 2 to
  1; `s$state[[1]]$glue` has length 8 (pre-change) where a fresh `storeState()` gives
  7, and `s$setState(s$state)` errors `state is not consistent with this sampler`.
  With `s$storeState()` first it succeeds; a width-PRESERVING change succeeds while
  silently reinstalling the stale trees.
- The behavior follows the documented opt-in convention at [[dbartsSampler-class.Rd:116@b102e17c]]; the sentence at [[dbartsSampler-class.Rd:328@b102e17c]]
  is the wrong one. Fix class: agent-fix (Rd).

### MINOR

- M1. `sample` validation inconsistent: bart/rbart raise `sample must be in 'train',
  'test'` ([[R/generics.R:385-391@b102e17c]]); all four own-class `extract` methods use bare
  `match.arg` -> `'arg' should be one of "train", "test"`.
- M2. `residuals(fit, sample = "train")` -> raw `formal argument "sample" matched by
  multiple actual arguments` ([[R/generics.R:814@b102e17c]], [[R/generics.R:2078@b102e17c]], [[R/generics.R:1564@b102e17c]] all pin it).
- M3. `plotTree(fit, treeNum = 1L, sample = 2L)` silently partial-matches `sampleNum`
  via `do.call` ([[R/generics.R:2096-2103@b102e17c]]). Benign today; same mechanism as G2.
- M4. `setForestBasis` on an amplitude-free sampler says `forest index out of range`
  ([[R/dbarts.R:1456-1458@b102e17c]]) where `setForestWeights` says `requires a sampler that
  carries forest amplitudes`; [[dbartsSampler-class.Rd:189@b102e17c]] promises the latter.
- M5. `extract(type = "trees")` omits the `chain` column on a single-chain fit;
  [[bart.Rd:262@b102e17c]] lists `chain, sample, tree` unconditionally.
- M6. [[bart.Rd:248@b102e17c]] says only `predict` stops on a fit saved without `storeState()`;
  `extract(type = "trees")` and `plotTree` stop identically.
- M7. [[dbartsSampler-class.Rd:308@b102e17c]] offers "the return of `storeState`" to `setState`,
  but it returns NULL invisibly ([[dbartsSampler-class.Rd:423@b102e17c]]), so `s$setState(s$storeState())` refuses.
- M8. `setForestBasis(k, ~var)` evaluates the formula in `environment(basis)`, never
  against the sampler's data ([[R/model.R:810-821@b102e17c]]), so `~x1` naming a predictor column
  raises raw `object 'x1' not found`; Rd:186 silent on the scope. M9. An amplitude
  sampler's `setTestPredictor` after a refused `setTestPredictorAndOffset` raises raw
  `argument is of length zero`.
- M10. `names(fit)` lists NULL-valued `yhat.test`, `yhat.test.mean`, `s.train`,
  `s.test` entries even when absent, so `nm %in% names(fit)` misleads.

## 3. Empty cells: what the default does today

- `plot` on bartOrdinal/bartNegbin/bartHurdle -> `plot.default`. Rd SILENT. (G7)
- `plotTree` and `survivalProbabilities` on all four own-class fits -> `no applicable
  method`; each carries a `$fit` sampler under `keepTrees`, so plotTree would have
  something to draw. [[plotTree.Rd:26-29@b102e17c]] names only bart/rbart/dbartsSampler. SILENT.
- `as_draws_array`/`as_draws_df` on the four -> posterior's default, `All list
  elements must be lists themselves.` summary.bart.Rd:5-8 scopes them. SILENT.
- `xbart()` returns a BARE array (dimnames named `rep`,`n.trees`,`k`), no class:
  `predict`/`extract`/`plotTree` -> `no applicable method`; `fitted`/`residuals` ->
  raw `$ operator is invalid for atomic vectors`; `summary` -> `summaryDefault` of the
  loss numbers; `print`/`plot` -> the numeric defaults. SILENT.
- `dbartsSampler`: `predict`/`survivalProbabilities` -> `no applicable method`;
  `fitted` -> `'fitted' is not a valid field or method name for reference class
  "dbartsSampler"`; `residuals` -> the same about `'na.action'`; `plot` -> `cannot
  coerce type 'S4' to vector of type 'double'`; `summary` -> `summaryDefault`;
  `print` -> `$show()`, correct. [[dbartsSampler-class.Rd:355@b102e17c]] states the division, so
  the absences AGREE; only the raw messages are ugly.
- `pdbart`/`pd2bart`: `fitted`/`residuals` return NULL SILENTLY; `summary` gives a
  `summaryDefault` of the list; `predict`/`extract` -> `no applicable method`. SILENT.

## 4. Partial-matching and argument-swallow audit (every generic's formals read)

Hazard exists wherever a method forwards its own formals, or `...`, into a callee
with different names. Full formals table in `generics/probe11.R`'s first block.

- `extract.bart` / `extract.rbart` -> `sampler$getTrees`: 6 hazards, all live.
  `sample` -> `sampleNums` (partial match, then NA); `combineChains`, `forest`,
  `contribution` -> `unused argument`. THE F3 CLASS. (G2)
- `plotTree.{bart,rbart,dbartsSampler}` -> `sampler$plotTree(treeNum, chainNum,
  sampleNum, treePlotPars, ...)` via `do.call`: a user's `sample=`/`chain=`
  partial-matches `sampleNum`/`chainNum`. Live but benign. (M3)
- `residuals.{bart,rbart,bartHurdle}` -> `fitted.*(..., sample = "train", ...)`: a
  user-supplied `sample=` collides. (M2)
- `fitted.{bart,rbart}` -> `extract(object, type, sample, ...)` POSITIONALLY into
  `(object, type, sample, combineChains, forest, contribution, ...)`: safe today, but
  a positional coupling to another method's formal ORDER - inserting a formal into
  `extract.bart` before `combineChains` would silently misroute `fitted`'s dots.
- The four own-class extract/fitted/residuals/summary methods: no forwarding, but a
  bare `...` that ABSORBS the bart-family vocabulary without error. (G3, G4)
- `predict.bart`'s `\dots` is documented "Not used" ([[bart.Rd:219@b102e17c]]) and indeed ignores
  `sample=`: AGREES, but a live trap because extract/fitted take `sample`.
  `predict.rbart` reads `dotsList[["value"]]` by exact name only ([[bart.Rd:1619-1624@b102e17c]]), so a
  partial `val=` is dropped rather than warned about. No two formals within any one
  method are in a prefix relation, so no user abbreviation is ambiguous.

## 5. Rd sweep: testable claims, by topic

62 testable claims extracted; 50 AGREE, 7 CONTRADICT, 5 SILENT-where-behaviour-exists.
Only the non-AGREE rows and the load-bearing AGREEs are listed.

- [[bart.Rd:201@b102e17c]] type vocabulary, `"response"`/`"link"` synonyms on all four generics,
  loglik combined = samples-by-observations and uncombined = chains-first
  (`12x40` / `2x6x40`), `"forest"` refused on a multinomial fit - AGREE.
- [[bart.Rd:203-205@b102e17c]] `sample` - CONTRADICTS under `type = "trees"` (G2).
- [[bart.Rd:215-216@b102e17c]] `ci.level` -> `est, ci.lower, ci.upper` - AGREES on
  bart/rbart/bartHurdle; ignored on the other three (G4).
- [[bart.Rd:219@b102e17c]] `\dots` reaches `getTrees` (`chainNums`, `sampleNums`, `treeNums`,
  `newdata`, all four verified); "Not used in predict" - AGREE.
- [[bart.Rd:247-248@b102e17c]] Saving - AGREES for bart/rbart/multinomial/ordinal/nbinom and the
  amplitude fit; SILENT on hurdle (G9), understates the affected generics (M6).
- [[bart.Rd:257-258@b102e17c]] - CONTRADICTS (G6). [[bart.Rd:260-266@b102e17c]] columns - CONTRADICTS on `chain` (M5)
  and on the keepSampler-only shape (G6). [[bart.Rd:263@b102e17c]] `newdata` changes the `n` column -
  AGREES (960 train vs 192 for 8 test rows). [[bart.Rd:315-317@b102e17c]] `n.chains` - CONTRADICTS (G5).
- [[bart.Rd:277@b102e17c]] chain-major collapse with all collapsed fields sharing the row order -
  AGREES (combined `sigma` == the two chains concatenated, chain 1 first).
  [[bart.Rd:288-289@b102e17c]], [[bart.Rd:294-295@b102e17c]] (varcount `12x3`, variable-named columns), [[bart.Rd:300-301@b102e17c]]
  (`forestFits` `6x40x2` named `forest1,forest2`, `glue`, `bases`, `n.forests`),
  [[bart.Rd:303@b102e17c]], [[bart.Rd:306-307@b102e17c]], [[bart.Rd:312-313@b102e17c]], [[bart.Rd:318-319@b102e17c]], [[bart.Rd:324-326@b102e17c]], [[bart.Rd:327-328@b102e17c]], [[bart.Rd:331@b102e17c]], [[bart.Rd:335@b102e17c]] - AGREE.
- summary.bart.Rd:11-17 (mean.s on a heteroscedastic fit), [[bart.Rd:19-25@b102e17c]], [[bart.Rd:42-51@b102e17c]]
  (`field[column]` naming), [[bart.Rd:64-77@b102e17c]], [[bart.Rd:79-82@b102e17c]]; [[plotTree.Rd:18-22@b102e17c]], [[plotTree.Rd:26-29@b102e17c]], [[plotTree.Rd:39-42@b102e17c]], [[plotTree.Rd:51@b102e17c]]
  (including the current-trees fallback); extract.dbartsSampler.Rd:8, [[plotTree.Rd:15@b102e17c]], [[plotTree.Rd:27@b102e17c]]
  (`type must be 'predictors'`; `40x3` with the predictor column names) - all AGREE.
  The first two are SILENT on the own-class fits.
- [[survivalProbabilities.Rd:40-41@b102e17c]], [[survivalProbabilities.Rd:49-52@b102e17c]], [[survivalProbabilities.Rd:98-105@b102e17c]] (`12x2x60` / `2x6x2x60`) and the
  rbart branch - AGREE. [[survivalProbabilities.Rd:45-47@b102e17c]] and [[survivalProbabilities.Rd:87-96@b102e17c]] - CONTRADICTS on a named-predictor hazard
  fit (G1).
- [[dbartsSampler-class.Rd:351@b102e17c]] - AGREES on all ten documented multinomial refusals
  (`setResponse`, `setOffset`, `setWeights`, `setSigma`, `setData`, `setModel`,
  `setCalibration`, `setForestWeights`, `setForestBasis`, `getFitsWithoutOffset`),
  each naming the capability, and on every open channel (`setCounts`,
  `setCategoryOffset`/`setCategoryTestOffset` with an n x K matrix, the predictor
  family, `setActiveRows`, `$predict` -> `10x3x5`). The builder's untested paragraph
  is now tested and clean. [[dbartsSampler-class.Rd:189@b102e17c]] - CONTRADICTS for `setForestBasis` (M4).
  [[dbartsSampler-class.Rd:308@b102e17c]] - CONTRADICTS (M7). [[dbartsSampler-class.Rd:328@b102e17c]] - CONTRADICTS (G10). [[dbartsSampler-class.Rd:358-359@b102e17c]] - AGREES.
- [[bart2.Rd:310@b102e17c]], [[bart2.Rd:312@b102e17c]], [[bart2.Rd:314@b102e17c]] (`cutpoint[1]`, `cutpoint[2]`), [[bart2.Rd:318@b102e17c]] (`dispersion`), [[bart2.Rd:322@b102e17c]]
  (hurdle `prob`/`link`/`log`/`ppd`, `ci.level`, two-component summary) - AGREE. [[bart2.Rd:38-49@b102e17c]]
  gives `combineChains` to `predict.bartMultinomial` and not to
  `extract.bartMultinomial`, saying nothing about what extract does with it: SILENT
  (G3). No `\usage` for the aliased extract/fitted/predict/print/residuals of
  bartOrdinal/bartNegbin/bartHurdle: SILENT.

## 6. VD-judgement groups (the deciding principle, not the cells)

- A. "Do the own-class generics carry the bart-family argument vocabulary, or refuse
  it by name?" - `combineChains`/`forest` on extract, `sample`/`ci.level` on fitted,
  `ci.level` on predict, `type` on residuals, `vars` on summary. All silently
  swallowed today; honour or refuse, silence is the wrong option. (G3, G4)
- B. "Does every own-class fit get `plot()`?" - bartOrdinal/bartNegbin/bartHurdle are
  the three without one; bartMultinomial shows what one looks like. (G7)
- C. "Is `extract(type = 'loglik')` a bart/rbart channel or every family's?" (G8)
- D. "Which generics does an `xbart` result support?" - a bare array today, so
  `fitted`/`residuals` misfire on stats' defaults. Class it, or leave it. (sec 3)
- E. "Does `dbartsSampler` ever grow the bart-family S3 surface?" - dbartsSampler-
  [[class.Rd:355@b102e17c]] already answers no; the open part is only whether the six defaults
  should be refused by name instead of leaking RC-field errors. (sec 3)
- F. "Should `extract(type = 'trees')` read the CURRENT trees on a keepSampler-only
  fit, or refuse?" - [[plotTree.Rd:39-42@b102e17c]] already documents that fallback for plotTree,
  so the two surfaces should agree either way. (G6)
- G. "Do `pdbart`/`pd2bart` get the fit generics, or refuse?" (sec 3)

## 7. What this pass did NOT check

- Numeric correctness: cells record class, dims, dimnames and attributes only, so the
  heteroscedastic loglik/ppd/summary defect the 08-24 value scan reported is neither
  confirmed nor refuted here (shapes are consistent everywhere). Likewise
  `contribution = TRUE` arithmetic and the `bases =` identity at [[bart.Rd:301@b102e17c]].
- `plot` CONTENT: an accepted `plot` cell means only that it drew without error.
- `pdbart`'s own argument surface, and `run()`'s shapes at [[dbartsSampler-class.Rd:387@b102e17c]].
