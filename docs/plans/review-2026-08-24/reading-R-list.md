# Reading review, R layer - candidate list

Read-only pass over R/*.R (25 files, 20623 lines) and NAMESPACE at b102e17c. Nothing here is a decision: every entry
carries a tag, the evidence, what removal would cost, and a confidence. Counts re-derived in-tree; method stated per
class. Base for classes 2/3: getParseData over all 25 files. 357 top-level function definitions - 27 exported, 55
registered S3, 275 internal; of the 275, 33 have ZERO call site in R/ and 79 exactly one. 596 stop() and 26 warning()
sites. 3573 comment lines (17.3% of R/).

## Top ten, by maintainer time saved per minute spent deciding

1. [test-surface] R/bartcore.R:1066-1535, a second mutation API (31 functions, ~470 lines) parallel to the R5 methods,
    zero callers in R/ - 6.1.
2. [duplication] R/bart.R's four family drivers share ~55 lines of scaffold each; two are near-verbatim twins - 1.1,
    1.2.
3. [duplication] 19 of 55 C entry points reached both from an R5 inline .Call and a bartcore.R wrapper, with
    independent coercion - 1.3.
4. [option-no-caller] R/rbart.R:47 rbart_vi(callback=) is wired, documented, and switches the sampler onto the R Gibbs
    loop; nothing anywhere exercises it - 3.1.
5. [comments] 102 of 129 process-narrating comment blocks cite a docs/ path .Rbuildignore strips - cl. 5.
6. [sediment] R/data.R:11-21 dataSlotOrNULL guards 3 slots against a class version that also lacked @bases and
    @response.type, both read bare everywhere - 4.3.
7. [duplication] Mirror-then-.Call-then-rollback (10 copies, ~120 lines) - 1.4; keepTrees refusal (12 sites, 8
    wordings, 4 gates for one condition) - 1.5.
8. [dead formal] R/bart.R:2833 makeind's `all` does nothing, and a test asserts as much - 3.2.
9. [sediment] R/bartcore.R:646 bartcoreSampler says it is test-only; two production bart2 paths call it and build a
    second engine - 6.2. R/bartcore.R:1486 bartcorePredictPerForest has no caller - 2.1.
10. [defensive] R/dbarts.R:295, :321, :344, three is.null() branches already excluded upstream - 4.1.

## 1. Duplication across the five surfaces

Method: clustered all 596 stop() messages by literal text (41 recur verbatim), then read each cluster. Listed only
where the CODE repeats, not just the message.
- 1.1 [duplication] R/bart.R:1436-1524 bart2Multinomial vs :1525-1602 bart2MultinomialCounts - ~40 of their 77/69
  lines identical: same buildSamplerPriors, same buildHostSamplerCall(family="multinomial", sigest), same
  `samplerCall$data <- y`, same setCategoryOffset / samplerOnly / bartcoreRun / packageMultinomialResults / keepTrees
  tail. They differ only in how (y, levels, K) are derived, and the comment at :1573 says so ("Mirrors
  bart2Multinomial's direct construction exactly"); docs/design/multinomial.md names no divergence. LOST: nothing.
  Confidence high.
- 1.2 [duplication + sediment] R/bart.R:1757-1894 bart2Ordinal vs :2008-2136 bart2Negbin - identical 12-argument
  signature, buildSamplerPriors call, buildHostSamplerCall/eval/samplerOnly block, n.chains/n.obs/n.test/n.samples
  block, bartcoreSampler+adoptPointer pair, varWidth computation, "drop the trailing singleton chain margin" block,
  keepTrees/keepSampler tail: ~55 shared lines each. ONE REAL DIVERGENCE - ordinal runs bartcoreRun(bc, n.burn,
  n.samples) once (:1824), nbinom loops per sample (:2073) for a per-sweep dispersion read; but
  src/R_interface_bartcore.cpp:4373-4376 allocates dispersion as a full per-draw channel ("one scalar per draw, so the
  dispersion channel takes sigma's own shape") - bart2Ordinal's own argument for its single run - and
  docs/design/bart-as-a-component.md sec 1 makes the forms bitwise identical, so a fold is testable. LOST: n.samples-1
  .Call round trips. Loop claim med-high.
- 1.3 [duplication] Two R routes to one bridge entry. Of 55 distinct C_dbarts_* symbols, 19 are reached both from an
  inline .Call in a dbartsSampler R5 method (R/dbarts.R) and from a wrapper in R/bartcore.R: create, getCalibration,
  getFitsWithoutOffset, getForestAmplitudes, getForestFits, getForestVariableCounts, getLatents, getTrees, predict,
  predictPerForest, setActiveRows, setCalibration, setForestBasis, setForestWeights, setModel, setState,
  setTestOffset, setWeights, storeState. They validate differently: R5 setActiveRows (R/dbarts.R:1383-1403) checks
  length, NA and 0/1-ness; bartcoreSetActiveRows (R/bartcore.R:1317-1322) coerces and defers to the bridge. LOST: the
  wrapper route is how tests reach bridge refusals without R's validation - a real seam (6.1) - but it needs no
  independent coercion per entry. High (mechanical `comm` over both files' C_ symbol sets).
- 1.4 [duplication] Mirror-slot / .Call / restore-on-error, ten ~12-line copies: R/dbarts.R:1244, :1340, :1429, :1473,
  :1595, :1635 and R/bartcore.R:314, :489, :523, :617. Each saves the old slot, writes the new, calls, restores in the
  handler, then `if (inherits(tryResult, "error")) stop(tryResult)`. LOST: the setModel copy alone rewrites e$call
  (:1248-1254), so a withMirrorRollback() needs one optional hook. Confidence high.
- 1.5 [duplication] keepTrees refusal: 12 sites (R/generics.R:229, :231, :247, :367, :371, :984, :1166, :1300, :1583,
  :1616, :1854, :2093), 8 wordings, FOUR gates for one condition - `is.null(object[["fit"]])` (:227), `is.null(fit) ||
  !keepTrees` (:981), `is.null(object[["cutpoints.raw"]])` (:1163), `is.null(object[["dispersion.raw"]])` (:1297). Two
  sites duplicate the same bart-vs-bart2 callName() branch (:228-232, :366-372), and predict.bartOrdinal (:1160-1174)
  / predict.bartNegbin (:1294-1308) share a 12-line prologue differing only in family string and field name. LOST: the
  wording, which requireKeptTrees(object, family, field) carries.
- 1.6 [duplication] Numeric 0/1-binary detection, three near-copies: R/spec.R:197-206, R/xbart.R:129-138,
  R/rbart.R:353-361. The first two are line-for-line the same but for aft being exempted only in spec.R; rbart.R's
  copy computes gatedFamily and DROPS the refusal (the per-chain dbarts() call raises it later, stated at :351-354).
  The 4-line core folds to numericResponseIsBinary(y). LOST: nothing - the exemption lists are family-set differences,
  not rule differences. Confidence high.
- 1.7 [duplication; ONE ARM INTENTIONAL] Thin-then-check, four surfaces, four behaviours: R/bart.R:799-808 (bart2,
  isTRUE(n.samples <= 0L)), R/bart.R:2709-2726 (bart, isTRUE(ndpost <= 0L)), R/rbart.R:101-106 (rbart_vi, bare `==
  0L`, message "no posterior draws will be taken after thinning"), R/xbart.R:51-53 + :94-95 (validates n.thin but
  NEVER divides - its n.samples is a kept-draw count, not a sweep budget). The xbart split is DOCUMENTED AND
  INTENTIONAL (man/bart2.Rd:152 spells out both conventions; man/xbart.Rd:114 says "n.samples are still returned
  regardless"). Residue: three copies of the division, four of the message, and rbart_vi's bare `== 0L` erroring on NA
  where the other two fall through.
- 1.8 [duplication] n.threads has four validators: R/A_class.R:295-296 + :349-350 (S4 validity), R/rbart.R:82-84
  (coerceOrError(...)[1L], silently takes element 1, tests `< 1L`), R/xbart.R:398-400 (coerceOrError with NO [1L],
  tests `<= 0L`), R/generics.R:294 (as.integer(...)[1L], no validation). xbart's is the odd one: is.na(c(2L,3L)) has
  length 2, so xbart(n.threads = c(2,3)) raises R's "the condition has length > 1" rather than a dbarts message; same
  shape at R/xbart.R:44-53 for n.cuts/useQuantiles/n.thin. High on the divergence; the length-2 failure is reasoned
  from coerceOrError (R/utility.R:182-206), not executed.
- 1.9 [duplication] Verbatim 9-line block twice: R/bart.R:1256-1264 (bart2) and :2820-2828 (bart), attaching
  $weights/$weights.test after packageBartResults; the second site's comment says "mirrors bart2's
  packageBartResults". packageBartResults already receives the sampler. LOST: nothing. Confidence high. (Neither field
  appears in man/bart.Rd or man/bart2.Rd.)
- 1.10 [duplication] The 13-element family token vector is written out verbatim twice (R/dbarts.R:376-389,
  R/bart.R:679-693). dbartsSpec's is an 8-element subset (R/spec.R:787-795), xbart's 4 (R/xbart.R:26), rbart_vi's 3
  (R/rbart.R:49), bart's 3 (R/bart.R:2643); only the 13-element pair is literal duplication. LOST: nothing; a shared
  constant serves match.arg.
- 1.11 [duplication, small] Three twin pairs, nothing lost by folding any: R/generics.R:1400-1406 validateType vs
  :1411-1422 resolveHurdleType (the second adds a "log" -> "bart" fold and rewrites the refusal; an `aliases =`
  argument covers it); R/bart.R:1387-1400 detectAutoMultinomial vs :1405-1418 detectAutoOrdinal (13-line twins over
  one detectAutoResponse result, differing only in the type predicate, both single-caller at :735/:749);
  R/generics.R:907-926 fitted.bartMultinomial vs :1113-1131 fitted.bartOrdinal (shared meanProbs + max.col + factor()
  block; :936-946 vs :1137-1146 share the indicator construction - ordinal's `ordered = TRUE` and multinomial's count
  branch parameterise).
- 1.12 [duplication] R/A_class.R:272-390 - dbartsControl's validity is a hand-unrolled table: 13 `length(x) != 1L`
  checks, 4 is.na checks, 7 range checks, each with its own string. ~118 lines a per-field spec renders in ~25. LOST:
  nothing; messages are mechanical from the field name. Pure reading cost, low risk. Confidence high.
- 1.13 [INTENTIONAL, stated] R/data.R:986-1022 re-checks offset and weights lengths with the same
  messages validateXYOffset/validateXYWeights use later, so the formula path does not surface model.frame's "variable
  lengths differ"; the reason is at :986-990. Listed so it is not re-flagged. Recommend keep.
- 1.14 [duplication, small] Same message, independent code, 2-4 sites each: "predictor columns cannot be entirely
  missing" R/data.R:1519/:1525/:1541; "newdata cannot be NULL" R/generics.R:989/:1171/:1305; "'times' must be finite
  and positive" R/bart.R:2410/:2500/:2544; "'n.samples' must be a positive integer" R/bart.R:807, R/dbarts.R:792,
  R/xbart.R:95; "'weights' must have the same length as 'y'" R/data.R:678/:1020, R/dbarts.R:1326/:1415; "multinomial
  counts must be a numeric matrix..." R/bartcore.R:979/:1069.
## 2. Single-caller helpers

Method: parsed every top-level `name <- function`, counted SYMBOL/SYMBOL_FUNCTION_CALL tokens outside the definition's
own srcref, excluded the 27 exports and 55 S3 registrations. All 79 single-callers were read; the large majority are
`factoring` - a long entry point split for readability, called from the one place that needs it (bart2Ordinal,
ingestFormulaTerms, walkFormulaTerms, rbart_vi_fit, sliceSample, predictBlend, resolveForests, resolveMonotone,
packageOrdinalResults, ...). Those are fine and not listed. The residue:
- 2.1 [dead] R/bartcore.R:1486-1496 bartcorePredictPerForest - zero callers in R/, inst/, man/, vignettes/,
  benchmarks/, tests/. The R5 predictForests (R/dbarts.R:1142-1154) inlines `.Call(C_dbarts_bartcore_predictPerForest,
  ...)` instead; the wrapper carries a 10-line comment duplicating the R5 docstring. LOST: nothing. Confidence high.
- 2.2 [test-surface] The other 30 zero-caller functions are all in R/bartcore.R:1066-1535 and all ARE called - from
  inst/tinytest and benchmarks/R via `:::`. See 6.1; a seam decision, not 30 decisions.
- 2.3 [factoring, near-inline] R/bartcore.R:24-26 samplerCarriesAmplitudes - a one-expression body with a 10-line
  comment, called once (:37); the same probe is open-coded at R/bart.R:182, R/spec.R:286/:568, R/dbarts.R:700. LOST:
  the comment's argument, worth keeping wherever it lands. 2.4 [sediment] R/data.R:11-13 dataSlotOrNULL, one caller
  (dataCounts, :19-21); see 4.3. (R/utility.R:317 setDefaultsFromFormals is also single-caller but pre-existing.)
- 2.6 [extension-point] R/dbarts.R:945-952 adoptPointer and :1790-1799 reapplyForestWeights are the only 2 of 46 R5
  methods with no mention in any man/*.Rd and no call from inst/tinytest, benchmarks/ or vignettes/. Both carry full
  R5 docstrings, so ?dbartsSampler-class shows them. adoptPointer is a named seam
  (docs/design/multinomial-mutation-arc.md:1037 calls it "a documented dbartsSampler$adoptPointer method");
  reapplyForestWeights is named at docs/design/bart-as-a-component.md:119. Internal callers: R/bart.R:1807, :2056;
  R/dbarts.R:1066, :1827, :1859. test-host-shell-pins.R:18-19 pins both names; test-forest-weights-r5.R:108 admits its
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
- 3.1 [option with no caller; HIGHEST value here] R/rbart.R:47 `rbart_vi(callback = NULL)` - wired (validated
  :327-332, threaded :370), documented (man/rbart.Rd:82, :149), and it SWITCHES THE SAMPLER: :384's `builtinTauPrior
  && is.null(callback)` sends a callback-carrying fit down the R Gibbs loop rather than the engine loop. Zero uses in
  inst/tinytest, \examples, vignettes/, benchmarks/R or the four sister repos - test-rbart-loop-callback.R and
  test-capi.R exercise a DIFFERENT internal per-sweep engine callback. R/rbart.R:43 `keepCall` is likewise never
  supplied to rbart_vi: the single "keepCall" hit, inst/tinytest/test-plot-generics.R:114, is a bart2() call (read,
  not grepped). LOST: nothing - the coverage gap IS the finding. Confidence high.
- 3.2 [dead formal, exported, documented as dead] R/bart.R:2833 `makeind(x, all = TRUE)`: the body is `ignored <- all
  ## for R check` then `makeModelMatrixFromDataFrame(x, TRUE)` - a hardcoded TRUE, never `all`; man/makeind.Rd:25-26
  says "Not currently implemented"; and inst/tinytest/test-makeModelMatrix.R:299 asserts `makeind(df, all = FALSE)` is
  identical to `makeind(df)`, which can only pass BECAUSE the formal does nothing. LOST: BayesTree signature
  compatibility, the reason makeind exists. Confidence high.
- 3.3 [surface gap] `bartOwnClassFamilies` (R/bart.R:2587-2592) lists multinomial/ordinal/nbinom/ hurdle.lognormal and
  is checked at :2647-2653 BEFORE match.arg so bart() can redirect to bart2() by name. "twopart" - bart2()'s and
  dbarts()'s own alias for hurdle.lognormal (resolved at R/bart.R:712, R/dbarts.R:405) - is missing from that list, so
  `bart(family = "twopart")` falls through to R's generic match.arg message, naming neither the alias nor bart2(). The
  token is exercised with bart2()/dbarts() at inst/tinytest/test-hurdle-surface.R:35, :71 and
  test-formula-terms.R:256, never with bart(). High.
- 3.4 [accepted but undocumented tokens] R/model.R:554, :557 - parseMonotoneSign accepts "inc" and "dec"; nothing in
  the repo or the sisters uses them, and its own refusal message (:559-562) does not list them, nor does any Rd. The
  character "0" branch (:558) is likewise unused - numeric 0 is well covered (inst/tinytest/test-monotone.R:38, :43,
  :61). LOST: nothing.
- 3.5 [zero exercise, low] R/bart.R:2624 `printevery` and :2638 `proposalprobs` (bart() formals, appearing only in
  man/bart.Rd \usage and \arguments); R/spec.R:782 `blocks` and :797 `dispersion` (dbartsSpec formals, never passed by
  keyword anywhere, stan4bart's dynamic bart_args forwarding included); R/spec.R:615 and :623 - the ordinal and
  multinomial arms of the forests=-plus-family refusal switch, whose aft and nbinom siblings ARE hit
  (inst/tinytest/test-bcf-family.R:457, :467).
- 3.6 [dead formal, pre-existing] R/utility.R:878 `quoteInNamespace(name, character.only = FALSE)` - nine call sites,
  none passes it. Present on main, so not branch residue.
- 3.7 [refuse-by-name; recommend keep] R/generics.R:1523 `extract.bartHurdle(sample = c("train","test"))` - `sample =
  "test"` always errors (:1529-1534); the formal exists so the refusal names the argument rather than R's "unused
  argument".
- 3.8 [NEGATIVE results, so they cost no time] (a) A mechanical dead-formal scan (formals minus every
  SYMBOL/SYMBOL_SUB/STR_CONST token in the body's srcref span) flags 13; all but 3.2/3.6 are false positives, because
  the five entries carry arguments through match.call()/redirectCall and a live formal need never appear as a symbol
  (pdbart/pd2bart's y.train reaches bart() via pdbart.prologue:43; bart2's 16 and dbarts's 6 likewise). rbart_vi_fit's
  `.chain.num.ignored` (R/rbart.R:872) is genuinely unread and says so in its name. (b) Eight S4 slots have no R-side
  `@` read outside their own validity (dart's b/rho/alpha/update.alpha, fixed's value, dbartsModel@p.birth,
  dbartsControl@useQuantiles/@printCutoffs) but ALL are read by the C++ bridge through the slot-as-attribute path
  (src/R_interface_bartcore.cpp:1209, :1832, the REPROTECT_SLOT calls) - not findings. (c) dbartsSpec's
  tree.prior/resid.prior/proposal.probs/ parentEnv look unexercised inside dbarts, but stan4bart passes all four
  (stan4bart/R/stan4bart_fit.R:554-558) - a dbarts-only grep would have called them dead.
## 4. Defensive code with no reachable trigger
- 4.1 [unreachable] R/dbarts.R validateArgumentsInEnvironment, three sites. :295-297 `if (is.null(n.samples))` runs
  after :285's as.integer(n.samples) and :292's `length != 1L` refusal; as.integer(NULL) is integer(0), so the length
  check always fires first. :321 `is.null(sigma) || sigma <= 0.0` sits inside `!missing(sigma) && !is.na(sigma)`
  (:309); `TRUE && logical(0)` is NA and `if (NA)` raises "missing value where TRUE/FALSE needed" (executed, R 4.6.1),
  so a NULL sigma gets R's error, never this one; :344 is the same shape for sigest. LOST: nothing - the NULL cases
  want an EARLIER check.
- 4.2 [partly unreachable] R/augmentation.R:44-52 augRestrict's second arm (`if (!supplied && family %in% families)
  stop("family \"%s\" requires '%s'")`) is dead for two of its five call sites: :77 and :80 hardcode `supplied =
  TRUE`. So dbartsDrawLatents(family = "aft", ...) with sigma omitted silently uses the default 1 instead of raising
  "family \"aft\" requires 'sigma'", and family = "logistic" with weights omitted likewise. The comment at :38-41
  states a rule the code does not enforce there ("its own family REQUIRES it, the draw law carrying no default to fall
  back on") - sigma HAS a default (:66). LOST: nothing if the arm stays for the three reachable sites; the comment
  should drop the other two. High.
- 4.3 [sediment] R/data.R:4-21. dataSlotOrNULL exists so "a saved bart/bart2 fit or dbartsSampler [that] holds a
  dbartsData built under the class definition in force when it was written" reads @counts / @offset.category /
  @offset.category.test as NULL rather than erroring. Three problems. (a) The stated invariant is violated in-tree:
  R/generics.R:992 reads `object$fit$data@offset.category` bare though the comment says "Every internal read of
  'counts', 'offset.category' and 'offset.category.test' goes through here". (b) Slots added in the same window are
  read bare everywhere - @bases at R/bart.R:182, R/spec.R:286/:562/:568/:574/:700, R/bartcore.R:25, R/dbarts.R:700;
  @response.type at R/data.R:502/:615, R/spec.R:196, R/xbart.R:128 - so any object old enough to lack @counts errors
  on those first and the guard cannot deliver. (c) `git log -S` dates the additions to this unreleased branch: @counts
  and dataSlotOrNULL at 5a3bc276 (2026-08-24), @bases at 983d7f0a (2026-08-14), @response.type at ee7bf84f
  (2026-07-17); none is on main, so under no-backwards-compat the only objects protected are ones serialized by an
  intermediate commit of this branch, and R/A_class.R:627-631 concedes the half-measure ("stays READABLE ... but it is
  not revalidatable"). LOST IF REMOVED: nothing the branch owes anyone - but a real 0.9-x-fit migration story starts
  here and must then cover @bases and @response.type too. High.
- 4.4 [divergent, unconfirmed] R/rbart.R:104 `if (control@n.samples == 0L)` where bart2 (R/bart.R:806) and bart
  (:2725) wrap the same test in isTRUE(... <= 0L). Reachable only if control@n.samples can be NA there; A_class.R:385
  permits NA_integer_ explicitly, but no concrete call was traced. Confidence low.
- 4.5 [stale precondition] R/bartcore.R:636-645 - bartcoreSampler's comment says "Internal interface used by the tests
  and the equivalence harness". It is production code: R/bart.R:1802 and :2051 call it on every bart2 ordinal and
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
- SAME DEFECT, NOT PREVIOUSLY COUNTED: 6 lines cite benchmarks/R/*.R (R/spec.R:476, R/bart.R:817, :1428, :1523, :2250,
  R/bartcore.R:956); `.Rbuildignore` also carries `^benchmarks$`.
- STALE PATH: R/diagnostics.R:4 says the posterior methods are registered "in zzz.R". There is no R/zzz.R; the
  registration is R/hooks.R:6-21.
- Comment lines as a share of the file: spec.R 0.289, bartcore.R 0.272, diagnostics.R 0.247, bart.R 0.210, model.R
  0.202, mixedMatrix.R 0.202; tree-wide 0.173. spec.R is both densest and second-highest hit count; diagnostics.R is
  dense and nearly clean, so density is not itself the problem.

Worst 15 sites (quotes verified; three marked PRE-EXISTING are VD's own, not branch residue):
- R/data.R:341-342 "this used to be a function evaluated in the caller's frame, but that causes warnings in R check so
  now it is just a block of code" - PRE-EXISTING (main:R/data.R:106).
- R/xbart.R:57 "control = is no longer a formal: xbart builds its own control..."
- R/xbart.R:191-192 "a custom control could once supply an n.trees a caller left unnamed; control = is gone".
- R/bart.R:769 "used to silently take dbarts()'s own default rather than the token/value this signature advertises" -
  a fixed defect; R/rbart.R:66 ships the same sentence.
- R/dbarts.R:627 "this combination used to fit an ordinary single-forest model with the declaration silently
  discarded. Refuse it by name instead."
- R/spec.R:17 "xbart once carried the weaker `family != \"gaussian\"`."
- R/spec.R:652 "the gaussian-only literal 0.5 this used to compare against always meant" - unreadable without the
  removed comparison.
- R/data.R:1294-1296 "(the bug: incorporation was only reachable through the formula path, since this branch always
  complete-cases-filtered first)".
- R/data.R:1031-1032 "completeness is validated below (previous versions silently na.omit-dropped them)".
- R/generics.R:1154 and :1287 "rides the same keepTrees gate a deleted $bc field used to" - twice.
- R/utility.R:118-119 "the passthrough that let it keep flowing through '...' is gone now that the rename has landed";
  R/bartcore.R:1018 "the same creation-time checks the retired dedicated entries used".
- R/model.R:1759-1760 and R/spec.R:489-490 both narrate "the removed flat n.trees.variance / power.variance /
  base.variance formals" - shipped twice.
- R/sliceSample.R:304 "leaves j == maxIter and used to raise the exhaustion error over a good sample".
  (R/rbart.R:483-484 "To be polite ... we set the seed back" and R/data.R:1218 "backwards compatibility of
  bart(x.train, y.train, x.test)" are both PRE-EXISTING.)

## 6. Over-engineering for one use
- 6.1 [test-surface / seam] R/bartcore.R:1066-1535 - a low-level handle API (a bcSampler env carrying $ptr and $x)
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
  (R/bart.R:1802) and bart2Negbin (:2051) build a host through dbarts(), then call bartcoreSampler(sampler, family =
  ...) which issues a SECOND C_dbarts_bartcore_create (R/bartcore.R:647-654), then sampler$adoptPointer(bc$ptr) so
  $fit wraps the second and abandons the first. The comments (:1804-1807, :2053-2056) give the reason: both creates in
  the same order keep the draw stream, i.e. avoid re-recording baselines. docs/design/multinomial-mutation-arc.md
  records the split - multinomial moved to direct construction (one create), ordinal/nbinom "adopt the engine that ran
  by pointer (bitwise, no re-record)". The simpler form already exists in-file at bart.R:1470-1491. LOST: a baseline
  re-record, the entire cost. High; a scheduling question, not a discovery.
- 6.3 [tidiness] R/A_class.R:272-390, the hand-unrolled validity table; see 1.12.
- 6.4 [checked and CLEARED - need no maintainer time] R/validateComposition.R; R/diagnostics.R (every summary method
  delegates to summary.bart); R/partialDependence.R (pdbart.* helpers genuinely shared by pdbart and pd2bart);
  R/mixedMatrix.R (every helper has 2+ callers); R/formulaTerms.R (547 lines from two entry sites, one coherent
  module); R/model.R's prior class hierarchy (no virtual class has a single implementation); 4 setMethod calls total,
  none a single-method generic of dbarts' own making; R/utility.R's metaprogramming (ifelse_3, evalx, redirectCall,
  addCallArgument, subTermInLanguage, quoteInNamespace) is PRE-EXISTING on main; R/rbart.R:1-6 rbart.priors is a
  2-entry registry with both exercised (prior = gamma, inst/tinytest/test-rbart-bartcore.R:95,
  benchmarks/R/equivalence.R:1427).

## 7. Two notes outside the six classes
- [pre-existing, low] dcauchy and dgamma (R/rbart.R:2, :4) and the `predict` generic (R/generics.R:1489-1490,
  R/bart.R:2440, :2512, :2561) are used but appear in no `importFrom(stats, ...)`; both hold on main and 0.9-33 ships
  that way, so a standing condition, not a branch regression (derived mechanically).
- [for the cross-surface pass] Family tokens per entry: dbarts 13, bart2 13, dbartsSpec 8, xbart 4, rbart_vi 3, bart
  3; man/dbartsSpec.Rd:40's "Both entry points call the same resolution" is true of the resolver, not the accepted
  set. dbarts() lists hurdle.lognormal then stop()s on it (R/dbarts.R:438-448) - accept-then-refuse, unlike bart()'s
  pre-match.arg redirect (3.3).
