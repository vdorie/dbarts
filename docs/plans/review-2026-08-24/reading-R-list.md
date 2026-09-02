# Reading review, R layer - candidate list

Read-only pass over R/*.R (25 files, 20623 lines) and NAMESPACE at b102e17c. Nothing here is a decision: every entry
carries a tag, the evidence, what removal would cost, and a confidence. Counts re-derived in-tree; method stated per
class. Base for classes 2/3: getParseData over all 25 files. 357 top-level function definitions - 27 exported, 55
registered S3, 275 internal; of the 275, 33 have ZERO call site in R/ and 79 exactly one. 596 stop() and 26 warning()
sites. 3573 comment lines (17.3% of R/).

## Top ten, by maintainer time saved per minute spent deciding

1. [test-surface] [[R/bartcore.R:1066-1535@b102e17c]], a second mutation API (31 functions, ~470 lines) parallel to the R5 methods,
    zero callers in R/ - 6.1.
2. [duplication] R/bart.R's four family drivers share ~55 lines of scaffold each; two are near-verbatim twins - 1.1,
    1.2.
3. [duplication] 19 of 55 C entry points reached both from an R5 inline .Call and a bartcore.R wrapper, with
    independent coercion - 1.3.
4. [option-no-caller] [[R/rbart.R:47@b102e17c]] rbart_vi(callback=) is wired, documented, and switches the sampler onto the R Gibbs
    loop; nothing anywhere exercises it - 3.1.
5. [comments] 102 of 129 process-narrating comment blocks cite a docs/ path .Rbuildignore strips - cl. 5.
6. [sediment] [[R/data.R:11-21@b102e17c]] dataSlotOrNULL guards 3 slots against a class version that also lacked @bases and
    @response.type, both read bare everywhere - 4.3.
7. [duplication] Mirror-then-.Call-then-rollback (10 copies, ~120 lines) - 1.4; keepTrees refusal (12 sites, 8
    wordings, 4 gates for one condition) - 1.5.
8. [dead formal] [[R/bart.R:2833@b102e17c]] makeind's `all` does nothing, and a test asserts as much - 3.2.
9. [sediment] [[R/bartcore.R:646@b102e17c]] bartcoreSampler says it is test-only; two production bart2 paths call it and build a
    second engine - 6.2. [[R/bartcore.R:1486@b102e17c]] bartcorePredictPerForest has no caller - 2.1.
10. [defensive] [[R/dbarts.R:295@b102e17c]], [[R/dbarts.R:321@b102e17c]], [[R/dbarts.R:344@b102e17c]], three is.null() branches already excluded upstream - 4.1.

## 1. Duplication across the five surfaces

Method: clustered all 596 stop() messages by literal text (41 recur verbatim), then read each cluster. Listed only
where the CODE repeats, not just the message.
- 1.1 [duplication] [[R/bart.R:1436-1524@b102e17c]] bart2Multinomial vs [[R/bart.R:1525-1602@b102e17c]] bart2MultinomialCounts - ~40 of their 77/69
  lines identical: same buildSamplerPriors, same buildHostSamplerCall(family="multinomial", sigest), same
  `samplerCall$data <- y`, same setCategoryOffset / samplerOnly / bartcoreRun / packageMultinomialResults / keepTrees
  tail. They differ only in how (y, levels, K) are derived, and the comment at [[R/bart.R:1573@b102e17c]] says so ("Mirrors
  bart2Multinomial's direct construction exactly"); docs/design/multinomial.md names no divergence. LOST: nothing.
  Confidence high.
- 1.2 [duplication + sediment] [[R/bart.R:1757-1894@b102e17c]] bart2Ordinal vs [[R/bart.R:2008-2136@b102e17c]] bart2Negbin - identical 12-argument
  signature, buildSamplerPriors call, buildHostSamplerCall/eval/samplerOnly block, n.chains/n.obs/n.test/n.samples
  block, bartcoreSampler+adoptPointer pair, varWidth computation, "drop the trailing singleton chain margin" block,
  keepTrees/keepSampler tail: ~55 shared lines each. ONE REAL DIVERGENCE - ordinal runs bartcoreRun(bc, n.burn,
  n.samples) once ([[R/bart.R:1824@b102e17c]]), nbinom loops per sample ([[R/bart.R:2073@b102e17c]]) for a per-sweep dispersion read; but
  [[src/R_interface_bartcore.cpp:4373-4376@b102e17c]] allocates dispersion as a full per-draw channel ("one scalar per draw, so the
  dispersion channel takes sigma's own shape") - bart2Ordinal's own argument for its single run - and
  docs/design/bart-as-a-component.md sec 1 makes the forms bitwise identical, so a fold is testable. LOST: n.samples-1
  .Call round trips. Loop claim med-high.
- 1.3 [duplication] Two R routes to one bridge entry. Of 55 distinct C_dbarts_* symbols, 19 are reached both from an
  inline .Call in a dbartsSampler R5 method (R/dbarts.R) and from a wrapper in R/bartcore.R: create, getCalibration,
  getFitsWithoutOffset, getForestAmplitudes, getForestFits, getForestVariableCounts, getLatents, getTrees, predict,
  predictPerForest, setActiveRows, setCalibration, setForestBasis, setForestWeights, setModel, setState,
  setTestOffset, setWeights, storeState. They validate differently: R5 setActiveRows ([[R/dbarts.R:1383-1403@b102e17c]]) checks
  length, NA and 0/1-ness; bartcoreSetActiveRows ([[R/bartcore.R:1317-1322@b102e17c]]) coerces and defers to the bridge. LOST: the
  wrapper route is how tests reach bridge refusals without R's validation - a real seam (6.1) - but it needs no
  independent coercion per entry. High (mechanical `comm` over both files' C_ symbol sets).
- 1.4 [duplication] Mirror-slot / .Call / restore-on-error, ten ~12-line copies: [[R/dbarts.R:1244@b102e17c]], [[R/dbarts.R:1340@b102e17c]], [[R/dbarts.R:1429@b102e17c]], [[R/dbarts.R:1473@b102e17c]],
  [[R/dbarts.R:1595@b102e17c]], [[R/dbarts.R:1635@b102e17c]] and [[R/bartcore.R:314@b102e17c]], [[R/bartcore.R:489@b102e17c]], [[R/bartcore.R:523@b102e17c]], [[R/bartcore.R:617@b102e17c]]. Each saves the old slot, writes the new, calls, restores in the
  handler, then `if (inherits(tryResult, "error")) stop(tryResult)`. LOST: the setModel copy alone rewrites e$call
  ([[R/bartcore.R:1248-1254@b102e17c]]), so a withMirrorRollback() needs one optional hook. Confidence high.
- 1.5 [duplication] keepTrees refusal: 12 sites ([[R/generics.R:229@b102e17c]], [[R/generics.R:231@b102e17c]], [[R/generics.R:247@b102e17c]], [[R/generics.R:367@b102e17c]], [[R/generics.R:371@b102e17c]], [[R/generics.R:984@b102e17c]], [[R/generics.R:1166@b102e17c]], [[R/generics.R:1300@b102e17c]], [[R/generics.R:1583@b102e17c]],
  [[R/generics.R:1616@b102e17c]], [[R/generics.R:1854@b102e17c]], [[R/generics.R:2093@b102e17c]]), 8 wordings, FOUR gates for one condition - `is.null(object[["fit"]])` ([[R/generics.R:227@b102e17c]]), `is.null(fit) ||
  !keepTrees` ([[R/generics.R:981@b102e17c]]), `is.null(object[["cutpoints.raw"]])` ([[R/generics.R:1163@b102e17c]]), `is.null(object[["dispersion.raw"]])` ([[R/generics.R:1297@b102e17c]]). Two
  sites duplicate the same bart-vs-bart2 callName() branch ([[R/generics.R:228-232@b102e17c]], [[R/generics.R:366-372@b102e17c]]), and predict.bartOrdinal ([[R/generics.R:1160-1174@b102e17c]])
  / predict.bartNegbin ([[R/generics.R:1294-1308@b102e17c]]) share a 12-line prologue differing only in family string and field name. LOST: the
  wording, which requireKeptTrees(object, family, field) carries.
- 1.6 [duplication] Numeric 0/1-binary detection, three near-copies: [[R/spec.R:197-206@b102e17c]], [[R/xbart.R:129-138@b102e17c]],
  [[R/rbart.R:353-361@b102e17c]]. The first two are line-for-line the same but for aft being exempted only in spec.R; rbart.R's
  copy computes gatedFamily and DROPS the refusal (the per-chain dbarts() call raises it later, stated at [[R/rbart.R:351-354@b102e17c]]).
  The 4-line core folds to numericResponseIsBinary(y). LOST: nothing - the exemption lists are family-set differences,
  not rule differences. Confidence high.
- 1.7 [duplication; ONE ARM INTENTIONAL] Thin-then-check, four surfaces, four behaviours: [[R/bart.R:799-808@b102e17c]] (bart2,
  isTRUE(n.samples <= 0L)), [[R/bart.R:2709-2726@b102e17c]] (bart, isTRUE(ndpost <= 0L)), [[R/rbart.R:101-106@b102e17c]] (rbart_vi, bare `==
  0L`, message "no posterior draws will be taken after thinning"), [[R/xbart.R:51-53@b102e17c]] + [[R/xbart.R:94-95@b102e17c]] (validates n.thin but
  NEVER divides - its n.samples is a kept-draw count, not a sweep budget). The xbart split is DOCUMENTED AND
  INTENTIONAL ([[man/bart2.Rd:152@b102e17c]] spells out both conventions; [[man/xbart.Rd:114@b102e17c]] says "n.samples are still returned
  regardless"). Residue: three copies of the division, four of the message, and rbart_vi's bare `== 0L` erroring on NA
  where the other two fall through.
- 1.8 [duplication] n.threads has four validators: [[R/A_class.R:295-296@b102e17c]] + [[R/A_class.R:349-350@b102e17c]] (S4 validity), [[R/rbart.R:82-84@b102e17c]]
  (coerceOrError(...)[1L], silently takes element 1, tests `< 1L`), [[R/xbart.R:398-400@b102e17c]] (coerceOrError with NO [1L],
  tests `<= 0L`), [[R/generics.R:294@b102e17c]] (as.integer(...)[1L], no validation). xbart's is the odd one: is.na(c(2L,3L)) has
  length 2, so xbart(n.threads = c(2,3)) raises R's "the condition has length > 1" rather than a dbarts message; same
  shape at [[R/xbart.R:44-53@b102e17c]] for n.cuts/useQuantiles/n.thin. High on the divergence; the length-2 failure is reasoned
  from coerceOrError ([[R/utility.R:182-206@b102e17c]]), not executed.
- 1.9 [duplication] Verbatim 9-line block twice: [[R/bart.R:1256-1264@b102e17c]] (bart2) and [[R/bart.R:2820-2828@b102e17c]] (bart), attaching
  $weights/$weights.test after packageBartResults; the second site's comment says "mirrors bart2's
  packageBartResults". packageBartResults already receives the sampler. LOST: nothing. Confidence high. (Neither field
  appears in man/bart.Rd or man/bart2.Rd.)
- 1.10 [duplication] The 13-element family token vector is written out verbatim twice ([[R/dbarts.R:376-389@b102e17c]],
  [[R/bart.R:679-693@b102e17c]]). dbartsSpec's is an 8-element subset ([[R/spec.R:787-795@b102e17c]]), xbart's 4 ([[R/xbart.R:26@b102e17c]]), rbart_vi's 3
  ([[R/rbart.R:49@b102e17c]]), bart's 3 ([[R/bart.R:2643@b102e17c]]); only the 13-element pair is literal duplication. LOST: nothing; a shared
  constant serves match.arg.
- 1.11 [duplication, small] Three twin pairs, nothing lost by folding any: [[R/generics.R:1400-1406@b102e17c]] validateType vs
  [[R/generics.R:1411-1422@b102e17c]] resolveHurdleType (the second adds a "log" -> "bart" fold and rewrites the refusal; an `aliases =`
  argument covers it); [[R/bart.R:1387-1400@b102e17c]] detectAutoMultinomial vs [[R/bart.R:1405-1418@b102e17c]] detectAutoOrdinal (13-line twins over
  one detectAutoResponse result, differing only in the type predicate, both single-caller at [[R/bart.R:735@b102e17c]]/[[R/bart.R:749@b102e17c]]);
  [[R/generics.R:907-926@b102e17c]] fitted.bartMultinomial vs [[R/generics.R:1113-1131@b102e17c]] fitted.bartOrdinal (shared meanProbs + max.col + factor()
  block; [[R/generics.R:936-946@b102e17c]] vs [[R/generics.R:1137-1146@b102e17c]] share the indicator construction - ordinal's `ordered = TRUE` and multinomial's count
  branch parameterise).
- 1.12 [duplication] [[R/A_class.R:272-390@b102e17c]] - dbartsControl's validity is a hand-unrolled table: 13 `length(x) != 1L`
  checks, 4 is.na checks, 7 range checks, each with its own string. ~118 lines a per-field spec renders in ~25. LOST:
  nothing; messages are mechanical from the field name. Pure reading cost, low risk. Confidence high.
- 1.13 [INTENTIONAL, stated] [[R/data.R:986-1022@b102e17c]] re-checks offset and weights lengths with the same
  messages validateXYOffset/validateXYWeights use later, so the formula path does not surface model.frame's "variable
  lengths differ"; the reason is at [[R/data.R:986-990@b102e17c]]. Listed so it is not re-flagged. Recommend keep.
- 1.14 [duplication, small] Same message, independent code, 2-4 sites each: "predictor columns cannot be entirely
  missing" [[R/data.R:1519@b102e17c]]/[[R/data.R:1525@b102e17c]]/[[R/data.R:1541@b102e17c]]; "newdata cannot be NULL" [[R/generics.R:989@b102e17c]]/[[R/generics.R:1171@b102e17c]]/[[R/generics.R:1305@b102e17c]]; "'times' must be finite
  and positive" [[R/bart.R:2410@b102e17c]]/[[R/bart.R:2500@b102e17c]]/[[R/bart.R:2544@b102e17c]]; "'n.samples' must be a positive integer" [[R/bart.R:807@b102e17c]], [[R/dbarts.R:792@b102e17c]],
  [[R/xbart.R:95@b102e17c]]; "'weights' must have the same length as 'y'" [[R/data.R:678@b102e17c]]/[[R/data.R:1020@b102e17c]], [[R/dbarts.R:1326@b102e17c]]/[[R/dbarts.R:1415@b102e17c]]; "multinomial
  counts must be a numeric matrix..." [[R/bartcore.R:979@b102e17c]]/[[R/bartcore.R:1069@b102e17c]].
## 2. Single-caller helpers

Method: parsed every top-level `name <- function`, counted SYMBOL/SYMBOL_FUNCTION_CALL tokens outside the definition's
own srcref, excluded the 27 exports and 55 S3 registrations. All 79 single-callers were read; the large majority are
`factoring` - a long entry point split for readability, called from the one place that needs it (bart2Ordinal,
ingestFormulaTerms, walkFormulaTerms, rbart_vi_fit, sliceSample, predictBlend, resolveForests, resolveMonotone,
packageOrdinalResults, ...). Those are fine and not listed. The residue:
- 2.1 [dead] [[R/bartcore.R:1486-1496@b102e17c]] bartcorePredictPerForest - zero callers in R/, inst/, man/, vignettes/,
  benchmarks/, tests/. The R5 predictForests ([[R/dbarts.R:1142-1154@b102e17c]]) inlines `.Call(C_dbarts_bartcore_predictPerForest,
  ...)` instead; the wrapper carries a 10-line comment duplicating the R5 docstring. LOST: nothing. Confidence high.
- 2.2 [test-surface] The other 30 zero-caller functions are all in [[R/bartcore.R:1066-1535@b102e17c]] and all ARE called - from
  inst/tinytest and benchmarks/R via `:::`. See 6.1; a seam decision, not 30 decisions.
- 2.3 [factoring, near-inline] [[R/bartcore.R:24-26@b102e17c]] samplerCarriesAmplitudes - a one-expression body with a 10-line
  comment, called once ([[R/bartcore.R:37@b102e17c]]); the same probe is open-coded at [[R/bart.R:182@b102e17c]], [[R/spec.R:286@b102e17c]]/[[R/spec.R:568@b102e17c]], [[R/dbarts.R:700@b102e17c]]. LOST:
  the comment's argument, worth keeping wherever it lands. 2.4 [sediment] [[R/data.R:11-13@b102e17c]] dataSlotOrNULL, one caller
  (dataCounts, [[R/data.R:19-21@b102e17c]]); see 4.3. ([[R/utility.R:317@b102e17c]] setDefaultsFromFormals is also single-caller but pre-existing.)
- 2.6 [extension-point] [[R/dbarts.R:945-952@b102e17c]] adoptPointer and [[R/dbarts.R:1790-1799@b102e17c]] reapplyForestWeights are the only 2 of 46 R5
  methods with no mention in any man/*.Rd and no call from inst/tinytest, benchmarks/ or vignettes/. Both carry full
  R5 docstrings, so ?dbartsSampler-class shows them. adoptPointer is a named seam
  ([[docs/design/multinomial-mutation-arc.md:1037@b102e17c]] calls it "a documented dbartsSampler$adoptPointer method");
  reapplyForestWeights is named at [[docs/design/bart-as-a-component.md:119@b102e17c]]. Internal callers: [[R/bart.R:1807@b102e17c]], [[R/bart.R:2056@b102e17c]];
  [[R/dbarts.R:1066@b102e17c]], [[R/dbarts.R:1827@b102e17c]], [[R/dbarts.R:1859@b102e17c]]. [[test-host-shell-pins.R:18-19@b102e17c]] pins both names; [[test-forest-weights-r5.R:108@b102e17c]] admits its
  test "stays green even with reapplyForestWeights deleted". LOST IF DROPPED FROM THE PUBLIC METHOD SET: adoptPointer
  is what bart2Ordinal/bart2Negbin rest on and what a host loop would want, and it is the most dangerous method on the
  class - it installs a raw externalptr whose provenance it does not check, as its docstring says. A surface decision,
  not a deletion.
## 3. Options with no caller

Method: formals of the five entries and their helpers read against inst/tinytest (167 files), man \examples
(brace-matched; 27 of 29 Rd carry one), vignettes/, benchmarks/R and the four sister repos (bartCause dbarts-1.0,
stan4bart bartcore, treatSens master, bairrtt main - branches confirmed, no checkout), plus the scan in 3.8. Every
zero-exercise claim was re-derived by reading its would-be call site, since same-name formals collide across entry
points.
- 3.1 [option with no caller; HIGHEST value here] [[R/rbart.R:47@b102e17c]] `rbart_vi(callback = NULL)` - wired (validated
  [[R/rbart.R:327-332@b102e17c]], threaded [[R/rbart.R:370@b102e17c]]), documented ([[man/rbart.Rd:82@b102e17c]], [[man/rbart.Rd:149@b102e17c]]), and it SWITCHES THE SAMPLER: [[R/rbart.R:384@b102e17c]]'s `builtinTauPrior
  && is.null(callback)` sends a callback-carrying fit down the R Gibbs loop rather than the engine loop. Zero uses in
  inst/tinytest, \examples, vignettes/, benchmarks/R or the four sister repos - test-rbart-loop-callback.R and
  test-capi.R exercise a DIFFERENT internal per-sweep engine callback. [[R/rbart.R:43@b102e17c]] `keepCall` is likewise never
  supplied to rbart_vi: the single "keepCall" hit, [[inst/tinytest/test-plot-generics.R:114@b102e17c]], is a bart2() call (read,
  not grepped). LOST: nothing - the coverage gap IS the finding. Confidence high.
- 3.2 [dead formal, exported, documented as dead] [[R/bart.R:2833@b102e17c]] `makeind(x, all = TRUE)`: the body is `ignored <- all
  ## for R check` then `makeModelMatrixFromDataFrame(x, TRUE)` - a hardcoded TRUE, never `all`; [[man/makeind.Rd:25-26@b102e17c]]
  says "Not currently implemented"; and [[inst/tinytest/test-makeModelMatrix.R:299@b102e17c]] asserts `makeind(df, all = FALSE)` is
  identical to `makeind(df)`, which can only pass BECAUSE the formal does nothing. LOST: BayesTree signature
  compatibility, the reason makeind exists. Confidence high.
- 3.3 [surface gap] `bartOwnClassFamilies` ([[R/bart.R:2587-2592@b102e17c]]) lists multinomial/ordinal/nbinom/ hurdle.lognormal and
  is checked at [[R/bart.R:2647-2653@b102e17c]] BEFORE match.arg so bart() can redirect to bart2() by name. "twopart" - bart2()'s and
  dbarts()'s own alias for hurdle.lognormal (resolved at [[R/bart.R:712@b102e17c]], [[R/dbarts.R:405@b102e17c]]) - is missing from that list, so
  `bart(family = "twopart")` falls through to R's generic match.arg message, naming neither the alias nor bart2(). The
  token is exercised with bart2()/dbarts() at [[inst/tinytest/test-hurdle-surface.R:35@b102e17c]], [[inst/tinytest/test-hurdle-surface.R:71@b102e17c]] and
  [[test-formula-terms.R:256@b102e17c]], never with bart(). High.
- 3.4 [accepted but undocumented tokens] [[R/model.R:554@b102e17c]], [[R/model.R:557@b102e17c]] - parseMonotoneSign accepts "inc" and "dec"; nothing in
  the repo or the sisters uses them, and its own refusal message ([[R/model.R:559-562@b102e17c]]) does not list them, nor does any Rd. The
  character "0" branch ([[R/model.R:558@b102e17c]]) is likewise unused - numeric 0 is well covered ([[inst/tinytest/test-monotone.R:38@b102e17c]], [[inst/tinytest/test-monotone.R:43@b102e17c]],
  [[inst/tinytest/test-monotone.R:61@b102e17c]]). LOST: nothing.
- 3.5 [zero exercise, low] [[R/bart.R:2624@b102e17c]] `printevery` and [[R/bart.R:2638@b102e17c]] `proposalprobs` (bart() formals, appearing only in
  man/bart.Rd \usage and \arguments); [[R/spec.R:782@b102e17c]] `blocks` and [[R/spec.R:797@b102e17c]] `dispersion` (dbartsSpec formals, never passed by
  keyword anywhere, stan4bart's dynamic bart_args forwarding included); [[R/spec.R:615@b102e17c]] and [[R/spec.R:623@b102e17c]] - the ordinal and
  multinomial arms of the forests=-plus-family refusal switch, whose aft and nbinom siblings ARE hit
  ([[inst/tinytest/test-bcf-family.R:457@b102e17c]], [[inst/tinytest/test-bcf-family.R:467@b102e17c]]).
- 3.6 [dead formal, pre-existing] [[R/utility.R:878@b102e17c]] `quoteInNamespace(name, character.only = FALSE)` - nine call sites,
  none passes it. Present on main, so not branch residue.
- 3.7 [refuse-by-name; recommend keep] [[R/generics.R:1523@b102e17c]] `extract.bartHurdle(sample = c("train","test"))` - `sample =
  "test"` always errors ([[R/generics.R:1529-1534@b102e17c]]); the formal exists so the refusal names the argument rather than R's "unused
  argument".
- 3.8 [NEGATIVE results, so they cost no time] (a) A mechanical dead-formal scan (formals minus every
  SYMBOL/SYMBOL_SUB/STR_CONST token in the body's srcref span) flags 13; all but 3.2/3.6 are false positives, because
  the five entries carry arguments through match.call()/redirectCall and a live formal need never appear as a symbol
  (pdbart/pd2bart's y.train reaches bart() via pdbart.prologue:43; bart2's 16 and dbarts's 6 likewise). rbart_vi_fit's
  `.chain.num.ignored` ([[R/rbart.R:872@b102e17c]]) is genuinely unread and says so in its name. (b) Eight S4 slots have no R-side
  `@` read outside their own validity (dart's b/rho/alpha/update.alpha, fixed's value, dbartsModel@p.birth,
  dbartsControl@useQuantiles/@printCutoffs) but ALL are read by the C++ bridge through the slot-as-attribute path
  ([[src/R_interface_bartcore.cpp:1209@b102e17c]], [[src/R_interface_bartcore.cpp:1832@b102e17c]], the REPROTECT_SLOT calls) - not findings. (c) dbartsSpec's
  tree.prior/resid.prior/proposal.probs/ parentEnv look unexercised inside dbarts, but stan4bart passes all four
  (stan4bart's `R/stan4bart_fit.R` lines 554-558) - a dbarts-only grep would have called them dead.
## 4. Defensive code with no reachable trigger
- 4.1 [unreachable] R/dbarts.R validateArgumentsInEnvironment, three sites. [[R/dbarts.R:295-297@b102e17c]] `if (is.null(n.samples))` runs
  after [[R/dbarts.R:285@b102e17c]]'s as.integer(n.samples) and [[R/dbarts.R:292@b102e17c]]'s `length != 1L` refusal; as.integer(NULL) is integer(0), so the length
  check always fires first. [[R/dbarts.R:321@b102e17c]] `is.null(sigma) || sigma <= 0.0` sits inside `!missing(sigma) && !is.na(sigma)`
  ([[R/dbarts.R:309@b102e17c]]); `TRUE && logical(0)` is NA and `if (NA)` raises "missing value where TRUE/FALSE needed" (executed, R 4.6.1),
  so a NULL sigma gets R's error, never this one; [[R/dbarts.R:344@b102e17c]] is the same shape for sigest. LOST: nothing - the NULL cases
  want an EARLIER check.
- 4.2 [partly unreachable] [[R/augmentation.R:44-52@b102e17c]] augRestrict's second arm (`if (!supplied && family %in% families)
  stop("family \"%s\" requires '%s'")`) is dead for two of its five call sites: [[R/augmentation.R:77@b102e17c]] and [[R/augmentation.R:80@b102e17c]] hardcode `supplied =
  TRUE`. So dbartsDrawLatents(family = "aft", ...) with sigma omitted silently uses the default 1 instead of raising
  "family \"aft\" requires 'sigma'", and family = "logistic" with weights omitted likewise. The comment at [[R/augmentation.R:38-41@b102e17c]]
  states a rule the code does not enforce there ("its own family REQUIRES it, the draw law carrying no default to fall
  back on") - sigma HAS a default ([[R/augmentation.R:66@b102e17c]]). LOST: nothing if the arm stays for the three reachable sites; the comment
  should drop the other two. High.
- 4.3 [sediment] [[R/data.R:4-21@b102e17c]]. dataSlotOrNULL exists so "a saved bart/bart2 fit or dbartsSampler [that] holds a
  dbartsData built under the class definition in force when it was written" reads @counts / @offset.category /
  @offset.category.test as NULL rather than erroring. Three problems. (a) The stated invariant is violated in-tree:
  [[R/generics.R:992@b102e17c]] reads `object$fit$data@offset.category` bare though the comment says "Every internal read of
  'counts', 'offset.category' and 'offset.category.test' goes through here". (b) Slots added in the same window are
  read bare everywhere - @bases at [[R/bart.R:182@b102e17c]], [[R/spec.R:286@b102e17c]]/[[R/spec.R:562@b102e17c]]/[[R/spec.R:568@b102e17c]]/[[R/spec.R:574@b102e17c]]/[[R/spec.R:700@b102e17c]], [[R/bartcore.R:25@b102e17c]], [[R/dbarts.R:700@b102e17c]];
  @response.type at [[R/data.R:502@b102e17c]]/[[R/data.R:615@b102e17c]], [[R/spec.R:196@b102e17c]], [[R/xbart.R:128@b102e17c]] - so any object old enough to lack @counts errors
  on those first and the guard cannot deliver. (c) `git log -S` dates the additions to this unreleased branch: @counts
  and dataSlotOrNULL at 5a3bc276 (2026-08-24), @bases at 983d7f0a (2026-08-14), @response.type at ee7bf84f
  (2026-07-17); none is on main, so under no-backwards-compat the only objects protected are ones serialized by an
  intermediate commit of this branch, and [[R/A_class.R:627-631@b102e17c]] concedes the half-measure ("stays READABLE ... but it is
  not revalidatable"). LOST IF REMOVED: nothing the branch owes anyone - but a real 0.9-x-fit migration story starts
  here and must then cover @bases and @response.type too. High.
- 4.4 [divergent, unconfirmed] [[R/rbart.R:104@b102e17c]] `if (control@n.samples == 0L)` where bart2 ([[R/bart.R:806@b102e17c]]) and bart
  ([[R/bart.R:2725@b102e17c]]) wrap the same test in isTRUE(... <= 0L). Reachable only if control@n.samples can be NA there; [[A_class.R:385@b102e17c]]
  permits NA_integer_ explicitly, but no concrete call was traced. Confidence low.
- 4.5 [stale precondition] [[R/bartcore.R:636-645@b102e17c]] - bartcoreSampler's comment says "Internal interface used by the tests
  and the equivalence harness". It is production code: [[R/bart.R:1802@b102e17c]] and [[R/bart.R:2051@b102e17c]] call it on every bart2 ordinal and
  nbinom fit. See 6.2.

## 5. Comments and docs that narrate process rather than state a constraint

Method: extracted every consecutive `#`/`##` block in all 25 files and judged each against the repo's own rule
(constraints, not provenance; literature citations are fine).
- COUNT 129 blocks. Per file, descending: bart.R 34, spec.R 16, dbarts.R 14, bartcore.R 11, data.R 11, generics.R 11,
  model.R 11, A_class.R 5, rbart.R 3, mixedMatrix.R 2, partialDependence.R 2, utility.R 2, xbart.R 2, augmentation.R
  1, diagnostics.R 1, formulaTerms.R 1, plot.R 1, sliceSample.R 1; nine at 0.
- 102 of the 129 cite a docs/ path. VERIFIED: 109 raw lines match `docs/` in R/ (152 more in src/, out of scope),
  naming 15 distinct .md paths, and ALL 15 EXIST. `.Rbuildignore` carries `^docs$`, so none resolves for an installed
  user. No `plans/` reference appears anywhere in R/.
- SAME DEFECT, NOT PREVIOUSLY COUNTED: 6 lines cite benchmarks/R/*.R ([[R/spec.R:476@b102e17c]], [[R/bart.R:817@b102e17c]], [[R/bart.R:1428@b102e17c]], [[R/bart.R:1523@b102e17c]], [[R/bart.R:2250@b102e17c]],
  [[R/bartcore.R:956@b102e17c]]); `.Rbuildignore` also carries `^benchmarks$`.
- STALE PATH: [[R/diagnostics.R:4@b102e17c]] says the posterior methods are registered "in zzz.R". There is no R/zzz.R; the
  registration is [[R/hooks.R:6-21@b102e17c]].
- Comment lines as a share of the file: spec.R 0.289, bartcore.R 0.272, diagnostics.R 0.247, bart.R 0.210, model.R
  0.202, mixedMatrix.R 0.202; tree-wide 0.173. spec.R is both densest and second-highest hit count; diagnostics.R is
  dense and nearly clean, so density is not itself the problem.

Worst 15 sites (quotes verified; three marked PRE-EXISTING are VD's own, not branch residue):
- [[R/data.R:341-342@b102e17c]] "this used to be a function evaluated in the caller's frame, but that causes warnings in R check so
  now it is just a block of code" - PRE-EXISTING (main:[[R/data.R:106@b102e17c]]).
- [[R/xbart.R:57@b102e17c]] "control = is no longer a formal: xbart builds its own control..."
- [[R/xbart.R:191-192@b102e17c]] "a custom control could once supply an n.trees a caller left unnamed; control = is gone".
- [[R/bart.R:769@b102e17c]] "used to silently take dbarts()'s own default rather than the token/value this signature advertises" -
  a fixed defect; [[R/rbart.R:66@b102e17c]] ships the same sentence.
- [[R/dbarts.R:627@b102e17c]] "this combination used to fit an ordinary single-forest model with the declaration silently
  discarded. Refuse it by name instead."
- [[R/spec.R:17@b102e17c]] "xbart once carried the weaker `family != \"gaussian\"`."
- [[R/spec.R:652@b102e17c]] "the gaussian-only literal 0.5 this used to compare against always meant" - unreadable without the
  removed comparison.
- [[R/data.R:1294-1296@b102e17c]] "(the bug: incorporation was only reachable through the formula path, since this branch always
  complete-cases-filtered first)".
- [[R/data.R:1031-1032@b102e17c]] "completeness is validated below (previous versions silently na.omit-dropped them)".
- [[R/generics.R:1154@b102e17c]] and [[R/generics.R:1287@b102e17c]] "rides the same keepTrees gate a deleted $bc field used to" - twice.
- [[R/utility.R:118-119@b102e17c]] "the passthrough that let it keep flowing through '...' is gone now that the rename has landed";
  [[R/bartcore.R:1018@b102e17c]] "the same creation-time checks the retired dedicated entries used".
- [[R/model.R:1759-1760@b102e17c]] and [[R/spec.R:489-490@b102e17c]] both narrate "the removed flat n.trees.variance / power.variance /
  base.variance formals" - shipped twice.
- [[R/sliceSample.R:304@b102e17c]] "leaves j == maxIter and used to raise the exhaustion error over a good sample".
  ([[R/rbart.R:483-484@b102e17c]] "To be polite ... we set the seed back" and [[R/data.R:1218@b102e17c]] "backwards compatibility of
  bart(x.train, y.train, x.test)" are both PRE-EXISTING.)

## 6. Over-engineering for one use
- 6.1 [test-surface / seam] [[R/bartcore.R:1066-1535@b102e17c]] - a low-level handle API (a bcSampler env carrying $ptr and $x)
  with 31 functions no other R/ file calls: bartcoreBCFSampler, bartcoreMultinomialSampler,
  bartcoreMultinomialCountSampler, bartcoreSet{Counts, CategoryOffset, CategoryTestOffset, ForestBasis, ForestWeights,
  ForestPriorScale, Offset, Response, ActiveRows, Weights, TestOffset, Data, TestPredictor, Predictor, CutPoints,
  State}, bartcoreUpdatePredictor{,PerObservation,PerObservationJointly}, bartcoreForest{Amplitudes, Fits,
  Calibration, VariableCounts}, bartcoreFitsWithoutOffset, bartcoreGet{Latents, Trees}, bartcoreStoreState,
  bartcorePredictPerForest. All but the last are reached from inst/tinytest and benchmarks/R via `:::`
  (bartcoreBCFSampler alone from 12 test files and 5 benchmarks). SIMPLER FORM: fold each into the R5 method that
  already reaches the same C entry (1.3) and let tests drive the R5 surface; the few bridge-refusal tests then need a
  deliberate escape hatch rather than a parallel API. LOST: the tests' ability to hit the bridge WITHOUT the R5
  layer's validation - which is what several assert - plus the equivalence harness's independence from R5. A real
  seam, gestured at by docs/design/bart-as-a-component.md's "graduation debt", but NO design doc names this API as it.
  RECOMMEND scoping the seam, not a deletion. High on the census; the call is VD's.
- 6.2 [sediment, recorded and priced] Two engine constructions per bart2 ordinal/nbinom fit. bart2Ordinal
  ([[R/bart.R:1802@b102e17c]]) and bart2Negbin ([[R/bart.R:2051@b102e17c]]) build a host through dbarts(), then call bartcoreSampler(sampler, family =
  ...) which issues a SECOND C_dbarts_bartcore_create ([[R/bartcore.R:647-654@b102e17c]]), then sampler$adoptPointer(bc$ptr) so
  $fit wraps the second and abandons the first. The comments ([[R/bart.R:1804-1807@b102e17c]], [[R/bart.R:2053-2056@b102e17c]]) give the reason: both creates in
  the same order keep the draw stream, i.e. avoid re-recording baselines. docs/design/multinomial-mutation-arc.md
  records the split - multinomial moved to direct construction (one create), ordinal/nbinom "adopt the engine that ran
  by pointer (bitwise, no re-record)". The simpler form already exists in-file at [[bart.R:1470-1491@b102e17c]]. LOST: a baseline
  re-record, the entire cost. High; a scheduling question, not a discovery.
- 6.3 [tidiness] [[R/A_class.R:272-390@b102e17c]], the hand-unrolled validity table; see 1.12.
- 6.4 [checked and CLEARED - need no maintainer time] R/validateComposition.R; R/diagnostics.R (every summary method
  delegates to summary.bart); R/partialDependence.R (pdbart.* helpers genuinely shared by pdbart and pd2bart);
  R/mixedMatrix.R (every helper has 2+ callers); R/formulaTerms.R (547 lines from two entry sites, one coherent
  module); R/model.R's prior class hierarchy (no virtual class has a single implementation); 4 setMethod calls total,
  none a single-method generic of dbarts' own making; R/utility.R's metaprogramming (ifelse_3, evalx, redirectCall,
  addCallArgument, subTermInLanguage, quoteInNamespace) is PRE-EXISTING on main; [[R/rbart.R:1-6@b102e17c]] rbart.priors is a
  2-entry registry with both exercised (prior = gamma, [[inst/tinytest/test-rbart-bartcore.R:95@b102e17c]],
  [[benchmarks/R/equivalence.R:1427@b102e17c]]).

## 7. Two notes outside the six classes
- [pre-existing, low] dcauchy and dgamma ([[R/rbart.R:2@b102e17c]], :4) and the `predict` generic ([[R/generics.R:1489-1490@b102e17c]],
  [[R/bart.R:2440@b102e17c]], [[R/bart.R:2512@b102e17c]], [[R/bart.R:2561@b102e17c]]) are used but appear in no `importFrom(stats, ...)`; both hold on main and 0.9-33 ships
  that way, so a standing condition, not a branch regression (derived mechanically).
- [for the cross-surface pass] Family tokens per entry: dbarts 13, bart2 13, dbartsSpec 8, xbart 4, rbart_vi 3, bart
  3; [[man/dbartsSpec.Rd:40@b102e17c]]'s "Both entry points call the same resolution" is true of the resolver, not the accepted
  set. dbarts() lists hurdle.lognormal then stop()s on it ([[R/dbarts.R:438-448@b102e17c]]) - accept-then-refuse, unlike bart()'s
  pre-match.arg redirect (3.3).
