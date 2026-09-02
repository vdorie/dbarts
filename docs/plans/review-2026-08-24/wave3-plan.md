# Fix wave 3 - implementer spec (bartcore @ 0045507c)

Source: consolidated-report.md sections (b)(c)(d)(e)(g), the
per-cell evidence in matrix-results.md / matrix-review-entries.md / matrix-review-generics.md,
and VD's judgements J1-J8. Waves 1-2 are b102e17c..07ad73e4 (landing note at
[[docs/plans/release-candidate-review.md:543@07ad73e4]]). J7 (engine default: arms, in a
private worktree) and J9 (predict n.threads) are OUT of this wave.

Every count and line reference below was re-verified in this pass; section E lists the
command behind each. All line numbers are against 0045507c.

## A. Landed-vs-remaining ledger (R surface)

LANDED (do not re-plan)
- 8042cc2c: B1 extract(type="trees") refuses sample/combineChains/forest/contribution by
  name ([[R/generics.R:360-372@8042cc2c]], wired at [[R/generics.R:398@8042cc2c]] and [[R/generics.R:1901@8042cc2c]]); B2 survivalProbabilities names the
  hazard period column on both branches ([[R/bart.R:2421@8042cc2c]], [[R/bart.R:2434@8042cc2c]]); minor [M2] residuals(sample=)
  refused by name (refuseResidualsSample, [[R/bart.R:842-851@8042cc2c]], on bart/rbart/bartHurdle).
- b657e8ae: M3 ([[dbartsSpec.Rd:40@b657e8ae]],[[dbartsSpec.Rd:48@b657e8ae]]), M4 ([[bart.Rd:153@b657e8ae]]), M5 ([[bart.Rd:165@b657e8ae]]), M7 ([[bart.Rd:248@b657e8ae]]
  hurdle storeState recipe), M8 ([[dbartsSampler-class.Rd:328@b657e8ae]]), M12 (three chain.hpp comments),
  U7 ([[kernel-vocabulary.md:26@b657e8ae]]), minors [M4] setForestBasis wording, [M5] chain column, [M7]
  storeState's NULL return - all Rd/comment.
- 7318b266 + 66ac05b3: M13 in full (rchk gate, sanitizers [main, master], lint/pkgdown PR
  trigger, three unrunnable harnesses deleted, logistic-reference gating + BART in CI).
- e35c8797: B8 rbart fitted/residuals SEGFAULT, R-side refusal plus the C-side bound.
- 52d3b5ff / fe505ae3 / 07ad73e4: reach gaps A4-A12, B3-B9, M14-M21, and (g.5)'s C9
  (prior.sd*k), C10 (zero-weight loglik NaN) and C11 (n.burn %/% n.thin) - all three pinned.

REMAINING (this wave)
- M1 dbartsData's false "'x' must have the same number of observations as 'y'" on every
  multi-column response -> slice r-surface (J5 decides the shape).
- M2 bart()'s family redirect covers 4 of 10 tokens; "twopart" missing -> r-surface (J3).
- M6 $n.chains absent from every keepSampler fit -> r-surface.
- M9 dbartsDrawLatents refuses its own formal default -> r-surface.
- M10 xbart answers family-then-response and n.threads length-2 with raw R errors ->
  r-surface.
- E19 defaultNodeScale has no default arm (silent NULL) -> r-surface.
- minor [M6-gen] [[bart.Rd:248@07ad73e4]] understates which generics stop on a stateless fit ->
  r-surface (Rd only).
- minor [M10-gen] names(fit) lists NULL-valued yhat.test/s.train/s.test -> r-surface.
- minor [M1-gen] own-class extract methods give bare "'arg' should be one of" where
  bart/rbart say "sample must be in 'train', 'test'" -> generics.
- minor [M3-gen] plotTree(sample=) partial-matches sampleNum -> generics.
- minor [M8-gen] setForestBasis(k, ~var) evaluates the formula in environment(basis) ->
  NOT in any J; see D5.
- J1 dots removal, J3 family redirects/typed-token echo, J5 one offset name, J8 monotone
  vocabulary + n.samples = 0 message + makeind note -> r-surface.
- J2 plot/loglik on the four own-class families, extract(type="trees") keepSampler
  fallback, pdbart/pd2bart/xbart Rd sentences -> generics.
- J4 own-class generics honour combineChains/ci.level and refuse the rest by name ->
  generics.
- J6 handle-API move, dead engine members, adoptPointer/reapplyForestWeights docs,
  core-generalization.md amendment -> handle-api. dataSlotOrNULL deletion moves to
  r-surface (it guards the slots J5 reshapes).

## B. Slice specs

Shared conventions for every refusal written in this wave:
- Echo the token/argument the CALLER TYPED. Where an entry folds an alias
  ("twopart" -> "hurdle.lognormal", [[R/bart.R:712@07ad73e4]], [[R/dbarts.R:405@07ad73e4]]) or remaps a family
  ("hazard"/"hazard.logistic" -> "probit"/"logistic", [[R/dbarts.R:534@07ad73e4]]), capture
  `requestedFamily <- family` BEFORE the fold and interpolate that in every downstream
  stop(); never the resolved token.
- Name the KIND received, not just the expected shape ("a length-n vector was
  supplied", "an n x 2 matrix was supplied", "a Surv response was supplied").
- Single quotes around argument names, escaped double quotes around family tokens -
  the package's existing style ([[R/bart.R:849@07ad73e4]], [[R/data.R:628@07ad73e4]]).

### (i) Slice "r-surface" - J1 + J3 + J5 + J8 + remaining S3 + dataSlotOrNULL

R1. J1 - delete the rejection-only dots channel.
  Current: bart2 ([[R/bart.R:704@07ad73e4]] `...`) and rbart_vi ([[R/rbart.R:53@07ad73e4]] `...`) declare a dots
  formal used only to produce a nicer error; [[R/bart.R:764-765@07ad73e4]] and [[R/rbart.R:61-62@07ad73e4]] call
  rejectUnknownDotsArgs(argNames, fn) ([[R/utility.R:122-142@07ad73e4]]), which consults
  retiredDotsNames ([[R/utility.R:120@07ad73e4]], one row: rngSeed -> seed) and otherwise agreps the
  nearest formal.
  Required: delete the `...` formal from both signatures, delete the two
  argNames/rejectUnknownDotsArgs blocks ([[R/bart.R:762-765@07ad73e4]], [[R/rbart.R:59-62@07ad73e4]]), delete
  rejectUnknownDotsArgs and retiredDotsNames and the 12-line comment above them
  ([[R/utility.R:108-142@07ad73e4]]). R's own "unused argument" wall then serves all six entries
  (bart, bart2, dbarts, dbartsSpec, xbart, rbart_vi) identically - which also closes
  E15's asymmetry without writing anything. Partial matching is unaffected: both dots
  were the LAST formal, so no formal moves behind a `...`.
  Rd: delete `\dots)` from the bart2 usage block ([[man/bart2.Rd:68@07ad73e4]]) and the rbart_vi
  usage block ([[man/rbart.Rd:38@07ad73e4]]); delete the \item{\dots} paragraph at [[man/bart2.Rd:260-262@07ad73e4]];
  drop `\dots` from the shared \item list at [[man/rbart.Rd:85@07ad73e4]]. Leave every method usage's
  \dots alone (predict/fitted/extract/plot keep theirs).

R2. J3 - bart()'s by-name family redirects, and the typed-token rule.
  Current: [[R/bart.R:2648-2655@07ad73e4]] refuses only bartOwnClassFamilies ([[R/bart.R:2589-2594@07ad73e4]]:
  multinomial, ordinal, nbinom, hurdle.lognormal) through refuseBartOwnClassFamily
  ([[R/bart.R:2596-2607@07ad73e4]]); the other six of the ten tokens [[man/bart.Rd:174@07ad73e4]] names fall through to
  match.arg's "'arg' should be one of "auto", "logistic", "aft"".
  Required: extend the pre-match.arg by-name branch to all ten tokens with three
  reasons, exactly as [[bart.Rd:174@07ad73e4]] already states them:
    - multinomial, ordinal, nbinom, hurdle.lognormal, twopart ->
      'bart() does not fit family = "<typed>"; use bart2(x.train, y.train, family =
      "<typed>")'   (existing text, extended to carry "twopart" unfolded)
    - gaussian, probit, hazard.probit ->
      'bart() does not fit family = "<typed>" as a token; it is what family = "auto"
      already fits for this response - drop the argument, or use bart2(x.train,
      y.train, family = "<typed>")'
    - hazard, hazard.logistic ->
      'bart() does not fit family = "<typed>": the discrete-time expansion needs
      "breaks" and "max.rows", which bart() does not have - use bart2(x.train,
      y.train, family = "<typed>")'
  Implement as one named vector bartRedirectedFamilies (token -> reason class) beside
  bartOwnClassFamilies so the Rd and the code list the same ten.
  Typed-token echo, the other half of J3 (E16/E17): in R/dbarts.R capture the typed
  token before [[R/dbarts.R:405@07ad73e4]]'s twopart fold and before [[R/dbarts.R:534@07ad73e4]]'s hazard remap and use it in the
  downstream refusals that currently name the resolved one - the sites reached are
  R/model.R's resid.dist refusal ("family \"probit\" has its own fixed error scale")
  and the variance-forest refusal. Do the same in [[R/bart.R:712@07ad73e4]] for bart2. Rule to
  write in the comment: the resolved token is an implementation detail; the caller can
  only act on what they typed.
  rbart_vi / xbart / dbartsSpec keep their narrow vocabularies ([[R/rbart.R:52@07ad73e4]]
  c("auto","gaussian","aft"), [[R/xbart.R:26@07ad73e4]] c("auto","gaussian","probit","logistic"),
  R/spec.R's eight). Their Rd sentences are corrected to say the vocabulary is
  narrower BY DESIGN and that the wider set lives on bart2 - one sentence each in
  man/rbart.Rd, man/xbart.Rd, man/dbartsSpec.Rd. E18 (bart(keepevery = -1) refused as
  'n.thin') is the same rule one layer down: bart() already coerces by the typed name
  ([[R/bart.R:2670-2677@07ad73e4]]) but the positivity refusal comes from dbartsControl's validity,
  which names n.thin - add a bart()-side positivity check for keepevery/ndpost/nskip
  naming the typed spelling before the control is built.

R3. J5 - one offset name, shape follows family.
  Current spellings (E8, all re-verified): bart2(offset = <n x K>) accepted;
  bart2(offset.category = ) -> unknown argument; dbarts(offset.category = ) -> raw
  unused argument; dbarts(offset = <n x K>) -> "'offset' must have the same length as
  'y'" ([[R/data.R:765@07ad73e4]], a length() test, so an n x 1 matrix passes and an n x K does
  not); dbartsData(counts = , offset.category = ) accepted ([[R/data.R:899-900@07ad73e4]]).
  Required:
  a. dbartsData: DELETE the offset.category and offset.category.test formals
     ([[R/data.R:899-900@07ad73e4]]) and the countsIsMissing clauses that read them ([[R/data.R:907-909@07ad73e4]]).
     The n x K train shift arrives as `offset`, the nTest x K test shift as
     `offset.test`. The dbartsData SLOTS keep their names ([[R/A_class.R:536-538@07ad73e4]]) and so
     do the R5 methods $setCategoryOffset/$setCategoryTestOffset (method names, not
     argument names - see D2).
  b. Routing: in dbartsData and in the dbarts()/bart2() matrix paths, a matrix-valued
     `offset` on a counts-carrying data object installs the category offset
     (validateDataCategoryOffset, [[R/data.R:1430-1443@07ad73e4]], unchanged); a matrix-valued
     offset anywhere else is refused; a vector-valued offset on a multinomial fit is
     refused.
  c. Messages (one helper each, used by dbartsData, dbarts, bart2, predict):
     flat-on-multinomial:
       'family = "multinomial" requires an n x K matrix "offset", one column per
       category; a length-<n> vector was supplied, and a common per-observation shift
       is the softmax's own null direction - it cancels'
     matrix-on-anything-else:
       ''offset' must be a numeric vector of length n or a single number; a <r> x <c>
       matrix was supplied, which only family = "multinomial" accepts'
     wrong-K:
       ''offset' must have one column per category (K = <K>); a <r> x <c> matrix was
       supplied'  (validateCategoryOffset's existing wording, kept)
  d. predict.bartMultinomial: rename the formal offset.category.test -> offset
     ([[R/generics.R:1017@07ad73e4]], body reads at [[R/generics.R:1033-1056@07ad73e4]]). SPECIFIED HERE, EXECUTED BY THE
     GENERICS SLICE (section C: J4 rewrites the same method; one owner per function).
     predict.bart already spells the
     new-row shift `offset` ([[R/generics.R:210@07ad73e4]]), so this is the same name for the same
     thing on both classes. Update [[man/bart2.Rd:80@07ad73e4]] (usage) and [[man/bart2.Rd:272-273@07ad73e4]] (the \item),
     and the two prose paragraphs that name the old spelling ([[man/bart2.Rd:229@07ad73e4]], [[man/bart2.Rd:312@07ad73e4]]).
  e. M1's false message: [[R/data.R:1178-1180@07ad73e4]] (sparse branch) and [[R/data.R:1237-1239@07ad73e4]] (numeric
     branch) compare NROW(formula) against NROW(codeResponse(data)$y) after
     codeResponse ([[R/data.R:458-473@07ad73e4]]) has flattened an n x 2 response to length 2n.
     Insert refuseMultiColumnResponse(data) BEFORE codeResponse in both branches:
       Surv:      ''y' is a survival response (Surv); dbartsData() takes a
                  single-column response - fit through dbarts()/bart2() with family =
                  "aft" or "hazard", which extract time and status first'
       n x 2:     ''y' is an n x 2 matrix; dbartsData() takes a single-column
                  response - a (time, status) pair goes to dbarts()/bart2() with
                  family = "aft"/"hazard", per-category counts to
                  dbartsData(counts = )'
       n x K:     ''y' is an n x <K> matrix; dbartsData() takes a single-column
                  response - pass per-category counts as dbartsData(counts = ) and fit
                  with family = "multinomial"'
     Fixing this one helper fixes all four inheriting surfaces (dbartsData positional,
     rbart_vi's matrix route, xbart, dbarts) - E4's own repro list.
  f. [[R/bart.R:867-874@07ad73e4]]'s multinomial offset.test refusal names
     dbarts:::bartcoreSetCategoryTestOffset as the internal channel, and
     [[man/bart2.Rd:229@07ad73e4]] repeats it. Under the handle-api slice that symbol leaves the
     namespace; see C for who owns the edit.

R4. J8 - option vocabularies.
  a. monotone: [[R/model.R:548-570@07ad73e4]] parseMonotoneSign switches on tolower(value) and
     accepts "inc", "dec", "0" alongside the documented set. DELETE those three switch
     arms; keep tolower() and DOCUMENT the case fold in [[man/dbarts.Rd:72@07ad73e4]] and
     [[man/bart2.Rd:203@07ad73e4]] ("matching is case-insensitive, so "Increasing" is accepted").
     The numeric 0 arm is untouched - [[man/dbarts.Rd:72@07ad73e4]] documents 0 for unconstrained.
     No existing pin asserts "inc"/"dec"/"0" (verified: 0 hits in inst/tinytest).
  b. makeind(all = ): NO CODE CHANGE. [[R/bart.R:2835-2838@07ad73e4]] binds `ignored <- all`;
     [[man/makeind.Rd:26@07ad73e4]] already says "Not currently implemented". Add the reason to
     that item - "retained for signature compatibility with BayesTree::makeind" - and
     stop there.
  c. n.samples = 0: dbartsControl/dbarts keep accepting it (the host-loop shape).
     bart2 ([[R/bart.R:806-809@07ad73e4]]), xbart ([[R/xbart.R:94-96@07ad73e4]]) and rbart_vi ([[R/rbart.R:104-106@07ad73e4]],
     currently "no posterior draws will be taken after thinning") share ONE message
     from one helper in R/utility.R:
       ''n.samples' must leave at least one draw after thinning (n.samples %/% n.thin
       = 0); dbarts() and dbartsControl() accept a zero-draw run - a sampler driven by
       a host loop - but <caller>() returns posterior draws'
     with <caller> the entry point's own name. Add one sentence to
     [[man/dbartsControl.Rd:32-33@07ad73e4]] recording the split (it already contrasts the
     per-run() return count with bart2's sweep budget - this appends the zero case).
     The multinomial branch keeps its own family-named check ([[R/bart.R:923@07ad73e4]]).

R5. Remaining S3 items.
  - M6: [[R/bart.R:389-393@07ad73e4]] sets result$n.chains only in the else branch. Set it
    UNCONDITIONALLY (2 lines), as bartMultinomial/bartOrdinal/bartNegbin already do;
    bartHurdle's assembly needs the same. [[bart.Rd:315-317@07ad73e4]] then becomes true.
  - M9: [[R/augmentation.R:65@07ad73e4]] declares sigma = 1 and [[R/augmentation.R:79-81@07ad73e4]] guards on !missing(sigma).
    Change the default to NULL and guard on !is.null(sigma); the aft/student arms that
    consume it already treat NULL as absent ([[R/augmentation.R:52@07ad73e4]]'s augRestrict).
  - M10: (i) move `family <- match.arg(family)` ([[R/xbart.R:123@07ad73e4]]) ABOVE the data build
    so a bad family is named before the response is ingested, matching the other four
    entries; (ii) give n.threads the length/positivity check the others have -
    [[R/xbart.R:398-401@07ad73e4]] currently coerces then tests is.na() on a length-2 value, which
    raises R's raw "'length = 2' in coercion to 'logical(1)'". Add
    'n.threads' must be of length 1' ([[A_class.R:296@07ad73e4]]'s wording) before the positivity
    test.
  - E19: [[R/model.R:400-415@07ad73e4]] defaultNodeScale's switch() has no default arm, so
    ("hazard") and ("student") return NULL silently. Add the sibling's stop() text
    ([[R/model.R:451@07ad73e4]]): 'no node scale is defined for family "<family>"'.
  - minor [M6-gen]: [[man/bart.Rd:248@07ad73e4]]'s Saving subsection says only predict stops on a
    stateless fit; extract(type = "trees") and plotTree stop identically. One clause.
  - minor [M10-gen]: R/bart.R's result assembly writes yhat.test/s.train/s.test
    entries as NULL, so names(fit) lists absent components. Assign only when non-NULL
    (result[["yhat.test"]] <- ... inside the existing if), or drop the NULL entries at
    the end with result[!vapply(result, is.null, NA)] - the second is one line and
    order-preserving.
  - dataSlotOrNULL deletion (from J6, homed here): [[R/data.R:11-13@07ad73e4]] with its 8-line
    comment (:4-10). Its comment claims every internal read of counts/offset.category/
    offset.category.test goes through it; FALSE - it has exactly ONE in-package use
    ([[R/data.R:20@07ad73e4]], inside dataCounts), while @counts is read bare at [[R/A_class.R:634@07ad73e4]] and
    @offset.category bare at [[R/generics.R:1034@07ad73e4]]. Inline the slot read into dataCounts
    and delete the function. It protects only objects serialized by an intermediate
    commit of THIS branch, which nothing outside the branch holds.
    Test consequence: [[inst/tinytest/test-multinomial-r5-surface.R:469-479@07ad73e4]] exercises it
    by hand-stripping attributes - delete that block.

Files touched: R/bart.R, R/rbart.R, R/xbart.R, R/data.R, R/dbarts.R, R/model.R,
R/augmentation.R, R/utility.R, R/spec.R (family Rd cross-reference only if the
narrow-vocabulary sentence needs a code comment); man/bart.Rd, man/bart2.Rd,
man/rbart.Rd, man/xbart.Rd, man/dbarts.Rd, man/dbartsData.Rd, man/dbartsControl.Rd,
man/dbartsSpec.Rd, man/makeind.Rd, man/dbartsAugmentation.Rd; inst/NEWS.Rd.

Test pins to add or rewrite (inst/tinytest):
- test-argument-surface.R - REWRITE the 8 dots pins ([[inst/tinytest/test-argument-surface.R:387@07ad73e4]], [[inst/tinytest/test-argument-surface.R:389@07ad73e4]], [[inst/tinytest/test-argument-surface.R:392@07ad73e4]], [[inst/tinytest/test-argument-surface.R:399@07ad73e4]], [[inst/tinytest/test-argument-surface.R:403@07ad73e4]],
  [[inst/tinytest/test-argument-surface.R:549@07ad73e4]], [[inst/tinytest/test-argument-surface.R:553@07ad73e4]], [[inst/tinytest/test-argument-surface.R:557@07ad73e4]]) to expect R's "unused argument"; the rngSeed pair loses its
  retirement text entirely (that is the J1 trade). [[inst/tinytest/test-argument-surface.R:563@07ad73e4]]'s dbarts() "unused argument"
  pin is unchanged and becomes the shared expectation. [[inst/tinytest/test-argument-surface.R:407-412@07ad73e4]]'s partial-match pin
  must stay green (it will).
- [[test-heteroscedastic.R:259-262@07ad73e4]] - the rbart_vi(variance = ) pin moves from
  "unknown argument" to "unused argument".
- [[test-bart-bart2.R:108-122@07ad73e4]] - ADD six by-name refusal pins (gaussian, probit,
  hazard.probit, hazard, hazard.logistic, twopart), each asserting the typed token
  appears in the message; keep the four existing ones.
- [[test-error-quality.R:38-47@07ad73e4]] - the n.samples = 0 pin text changes to the shared
  message; ADD the rbart_vi and xbart arms so all three entries are pinned to ONE
  string. [[test-xbart-error.R:24@07ad73e4]] - same rewrite.
- test-data-compatibility.R or a new test-response-shape.R - ADD the three by-kind
  ingest refusals (Surv, n x 2, n x K), each asserting the message names the kind;
  plus the four inheriting surfaces (dbartsData positional, rbart_vi matrix,
  xbart, dbarts) reaching the same text.
- test-multinomial-r5-surface.R - REWRITE [[inst/tinytest/test-multinomial-r5-surface.R:132@07ad73e4]], [[inst/tinytest/test-multinomial-r5-surface.R:136@07ad73e4]], [[inst/tinytest/test-multinomial-r5-surface.R:140@07ad73e4]], [[inst/tinytest/test-multinomial-r5-surface.R:144@07ad73e4]], [[inst/tinytest/test-multinomial-r5-surface.R:434@07ad73e4]], [[inst/tinytest/test-multinomial-r5-surface.R:449@07ad73e4]] from
  offset.category=/offset.category.test= to offset=/offset.test=; DELETE [[inst/tinytest/test-multinomial-r5-surface.R:469-479@07ad73e4]]
  (dataSlotOrNULL).
- test-multinomial-generics.R - REWRITE the 5 offset.category.test= sites ([[inst/tinytest/test-multinomial-generics.R:288@07ad73e4]], [[inst/tinytest/test-multinomial-generics.R:295@07ad73e4]],
  [[inst/tinytest/test-multinomial-generics.R:302@07ad73e4]], [[inst/tinytest/test-multinomial-generics.R:310@07ad73e4]], [[inst/tinytest/test-multinomial-generics.R:321@07ad73e4]]) to offset=. Ships WITH the generics slice, per R3(d).
- test-monotone.R - ADD an "inc" refusal pin and a case-fold acceptance pin
  ("Increasing" equals "increasing").
- test-argument-surface.R or test-error-quality.R - ADD the matrix-offset and
  flat-offset refusal pins, and an xbart(n.threads = c(1, 2)) pin naming the argument,
  and an xbart(<Surv>, family = "zzz") pin showing the family is named first.
- test-plot-generics.R / test-nbinom.R - ADD a defaultNodeScale("hazard") refusal pin
  (dbarts:::defaultNodeScale, which now stops).
- test-augmentation.R - ADD dbartsDrawLatents(..., sigma = 1) on a probit fit returning
  the same vector as the omitted call (M9's regression).
- ADD an n.chains pin on a keepTrees fit (test-bart-bart2.R or
  test-control-valuesAreUsed.R): expect_equal(fit$n.chains, 2L) under keepSampler.

Dense-line budget: code ~150 lines net (dots removal is net negative: -35 in
R/utility.R, -6 at the two call sites; family redirects +25; offset routing +45;
by-kind ingest +20; J8 +20; S3 items +25; dataSlotOrNULL -15), Rd ~40 lines, tests
~180 lines (30 rewritten, 150 added).

Gates: full battery. Trio 43/12/11 must stay bitwise WITHOUT re-record - none of these
paths moves a draw. Call sites whose SPELLING changes and that the trio or the tinytest
suite executes: the 6 dbartsData(offset.category=) sites in
test-multinomial-r5-surface.R and the three n.samples = 0 pins (the 5
predict(offset.category.test=) sites ride the generics slice). benchmarks/ has ZERO
offset.category hits, so no harness spelling changes. J1 changes no call that any harness makes (benchmarks pass only real formals).
Discrimination proof: every new refusal shown failing on the 0045507c build.
Rd topics: bart, bart2, rbart, xbart, dbarts, dbartsData, dbartsControl, dbartsSpec,
makeind, dbartsAugmentation.
Anchor drift: branch off 0045507c. Only inst/tinytest/test-xbart-*.R and inst/NEWS.Rd
can move under it - the xbart-oracle worktree holds test-xbart-reproducibility.R (stay
out of that file) and any further wave-3 slice appends at [[NEWS.Rd:1854@0045507c]].

### (ii) Slice "generics" - J2 + J4 + trees fallback + the three Rd sentences

G1. J2 - plot on the four own-class families. Current: only plot.bartMultinomial
  exists ([[R/plot.R:188-210@0045507c]], registered NAMESPACE S3method(plot, bartMultinomial));
  plot(bartOrdinalFit) reaches plot.default and raises "'x' is a list, but does not
  have components 'x' and 'y'". Required: plot.bartOrdinal, plot.bartNegbin,
  plot.bartHurdle in R/plot.R + three NAMESPACE registrations.
  >>> SURVEY SLOT A: the per-family plot semantics (what is traced, and against what)
  come from the generics-survey agent; plot.bartMultinomial's "trace each category's
  training-mean predicted probability over the kept draws" ([[man/bart2.Rd:312@0045507c]]) is the
  shape to match. Implement whatever the survey returns; the file, registration and
  Rd obligations below are fixed regardless.

G2. J2 - extract(type = "loglik") on the four own-class families. Current: all four
  refuse by name through validateType ([[R/generics.R:1442-1448@0045507c]]) because "loglik" is
  absent from their type vocabularies ([[R/generics.R:906@0045507c]] multinomial, [[R/generics.R:1115@0045507c]] ordinal,
  [[R/generics.R:1281@0045507c]] negbin, [[R/generics.R:1564@0045507c]] hurdle); [[bart.Rd:201@0045507c]] documents "loglik" on the extract generic
  unscoped, and it is the channel loo/WAIC consume. Required: add "loglik" to each of
  the four vocabularies and give pointwiseLogLikelihood ([[R/generics.R:56-169@0045507c]]) an arm
  per family; it already switches on object[["family"]] and already handles the
  weights and chain-shape bookkeeping.
  >>> SURVEY SLOT B: the per-family density (multinomial: the count/one-hot
  multinomial log-pmf at the reported probabilities; ordinal: the observed category's
  cumulative-probit difference - the engine's own per-observation channel, model.hpp:
  3314, is the oracle; nbinom: dnbinom at the per-draw dispersion; hurdle: the
  two-part mixture, the zero spike plus the lognormal branch) is the survey's to fix.
  Whatever it returns must satisfy: shape identical to extract(type = "ev") minus any
  K margin, zero-weight rows flagged NaN as the gaussian arm does ([[R/generics.R:122@0045507c]]).

G3. J2 - extract(type = "trees") on a keepSampler-only fit. Current: the guard at
  [[R/generics.R:385-396@0045507c]] tests is.null(object$fit), i.e. the SAMPLER, so a
  keepTrees = FALSE / keepSampler = TRUE fit falls through and getTrees returns the
  CURRENT working trees as an 11 x 4 frame (tree, n, var, value - no chain, no sample
  column, no warning). Required: KEEP that behavior (it is plotTree's documented
  fallback, [[man/plotTree.Rd:39-42@0045507c]]) and DISCLOSE it in man/bart.Rd's "Extracting Trees"
  subsection ([[man/bart.Rd:257-266@0045507c]]): state that without keepTrees but with a kept sampler the
  frame holds the sampler's CURRENT trees and omits the chain/sample index columns,
  in the same words [[plotTree.Rd:39-42@0045507c]] uses. No message is emitted (plotTree emits
  none). Pin both shapes.

G4. J4 - the own-class argument vocabulary. Current formals are
  (object, type, sample, ...) on extract.bartMultinomial ([[R/generics.R:904@0045507c]]), .bartOrdinal ([[R/generics.R:1113@0045507c]]),
  .bartNegbin ([[R/generics.R:1279@0045507c]]), and the bare "..." silently swallows the bart-family vocabulary:
  ten verified cells, each identical() to the call without the argument -
  extract(combineChains=), extract(forest=), fitted(sample=), fitted(ci.level=),
  predict(ci.level=), residuals(type=), summary(vars=).
  Required:
  a. HONOUR combineChains on extract for the three classes: add the formal, reshape
     through combineOrUncombineChains(x, fitNChains(object), combineChains)
     ([[R/generics.R:196-205@0045507c]]), the helper predict on the same fits already uses.
  b. HONOUR ci.level on fitted and predict for the three classes, through
     posteriorInterval ([[R/generics.R:170-195@0045507c]]), returning est/ci.lower/ci.upper exactly
     as bart/rbart/bartHurdle do ([[bart.Rd:215-216@0045507c]]).
  c. REFUSE by name, one helper (refuseUnusedGenericArgs(dots, generic, class)) reading
     names(list(...)) and stopping on the first hit:
       forest:       "'forest' is not used by <generic> on a <class> fit: <single-forest
                     reason | a multinomial fit's K category forests are not identified
                     individually - the identified content is the reported probabilities>"
       contribution: "'contribution' is not used by <generic> on a <class> fit"
       sample:       "'sample' is not used by fitted on a <class> fit: fitted values are
                     always the training rows; use extract(object, sample = \"test\")"
       type:         "'type' is not used by residuals on a <class> fit: the residual is
                     <the per-category observed proportion minus the fitted probability
                     | ...>"
       vars:         "'vars' is not used by summary on a bartMultinomial fit: it pools
                     the per-category mean-probability channel, which selects nothing"
     summary.bartOrdinal/.bartNegbin already carry a real vars ([[R/diagnostics.R:229-247@0045507c]])
     and are untouched; only summary.bartMultinomial (R/diagnostics.R, the (object, ...)
     one) refuses.
  d. Same refusal helper covers the plotTree/survivalProbabilities absences: register
     plotTree.bartMultinomial/.bartOrdinal/.bartNegbin/.bartHurdle and
     survivalProbabilities.bartMultinomial/... (8 tiny methods + 8 NAMESPACE lines)
     that stop by name rather than leaving R's "no applicable method":
       "plotTree is defined for bart, rbart_vi and dbartsSampler fits; a <class> fit's
       trees live on its sampler - call plotTree(object$fit, ...) <or, for hurdle,
       on object$occupancy$fit / object$positive$fit>"
       "survivalProbabilities applies to a discrete-time hazard fit (bart2(family =
       \"hazard\")); a <class> fit has no hazard channel"

G5. minor [M1-gen] - the four own-class extract methods validate sample with bare
  match.arg ([[R/generics.R:911@0045507c]], [[R/generics.R:1120@0045507c]], [[R/generics.R:1286@0045507c]], [[R/generics.R:1570@0045507c]]) and give "'arg' should be one of
  ...", where bart/rbart say "sample must be in 'train', 'test'" ([[R/generics.R:415@0045507c]],
  [[R/generics.R:815@0045507c]], [[R/generics.R:1947@0045507c]], [[R/generics.R:2065@0045507c]]). Add validateSample() beside validateType and use it at all
  eight sites, so one wording serves every class.

G6. minor [M3-gen] - plotTree.bart ([[R/generics.R:2134-2149@0045507c]]) builds args via
  list(treeNum = treeNum, ...) and do.call, so a caller's sample=/chain= partial-
  matches sampleNum/chainNum and silently draws. Refuse sample/chain by name in the
  same helper: "'sample' is not used by plotTree; the saved sample is 'sampleNum'".
  plotTree.rbart ([[R/generics.R:2151@0045507c]]) takes the same treatment.

G7. J2 - the three unclassed results, one Rd sentence each.
  man/pdbart.Rd \value ([[man/pdbart.Rd:82@0045507c]]): "The result carries class \"pdbart\"/\"pd2bart\" for its
  plot method only; it is not a fit, so predict, extract, fitted and residuals are not
  defined for it (fitted and residuals fall through to stats' defaults and return
  NULL)."
  man/xbart.Rd \value ([[man/xbart.Rd:129-133@0045507c]]): "The result is a bare array with no class, so the fit
  generics - predict, extract, fitted, residuals - do not apply to it; it is a table of
  losses, not a fit."
  man/summary.bart.Rd (:5-8, where as_draws_* is scoped): one clause recording that the
  four own-class fits are not draws objects.

Files touched: R/generics.R, R/plot.R, R/diagnostics.R, NAMESPACE (3 plot + 4 plotTree
+ 4 survivalProbabilities registrations = 11 new S3method lines), man/bart.Rd,
man/bart2.Rd, man/pdbart.Rd, man/xbart.Rd, man/plotTree.Rd, man/survivalProbabilities.Rd,
man/summary.bart.Rd, inst/NEWS.Rd.

Test pins:
- test-multinomial-generics.R, test-ordinal.R, test-nbinom.R, test-hurdle.R - ADD, per
  class: a plot() smoke pin (expect_silent on a null device, plus a panel/content
  assertion the survey names - NOT expect_silent alone, which is the M14/print.bart
  shape that cannot fail); an extract(type = "loglik") pin against an independently
  coded R oracle (max abs diff 0 or < 1e-12), the test-heteroscedastic-channels.R
  idiom; a combineChains = FALSE shape pin that DIFFERS from the combined shape; a
  ci.level pin asserting names(est, ci.lower, ci.upper); and one refusal pin per
  refused argument naming the argument.
- test-sampler-trees.R (or test-plot-generics.R) - ADD the keepSampler-only
  extract(type = "trees") pin: column set is (tree, n, var, value) and the keepTrees
  fit's is (sample, tree, n, var, value); both non-empty.
- test-plot-generics.R - ADD the plotTree sample=/chain= refusals.
- test-multinomial-generics.R etc. - the sample-message change (G5) breaks any pin
  asserting "should be one of" for a bad sample: verified NONE exists (the 12 "should
  be one of" pins are family/factors/missing/storage/type pins, listed in E).

Dense-line budget: code ~230 lines (plot methods ~60 pending the survey, loglik arms
~70 pending the survey, J4 honour/refuse ~70, refusing method stubs ~30), Rd ~25,
tests ~220. This is the largest slice and the only one with an external dependency.

Gates: full battery. Draw-moving risk: extract(type = "ppd") and the new loglik arms
must not touch the RNG on the default path (the multinomial ppd already does, and is
unchanged). Trio 43/12/11 bitwise, no re-record - the trio constructs no own-class
generic call (verified: bcf/multinomial-equivalence call the handle API, not the S3
generics). Rd topics: bart, bart2, pdbart, xbart, plotTree, survivalProbabilities,
summary.bart.
Anchor drift: branch off 0045507c; rebase onto the r-surface slice before landing (they
share inst/NEWS.Rd and man/bart2.Rd). R/generics.R and R/plot.R are otherwise this
slice's alone.

### (iii) Slice "handle-api" - the rest of J6

H1. Move the 31 handle wrappers out of the namespace.
  The set is exactly the contiguous block [[R/bartcore.R:1066-1535@0045507c]], which holds 32
  definitions: bartcoreSetCounts, bartcoreSetCategoryOffset,
  bartcoreSetCategoryTestOffset, bartcoreSetForestBasis, bartcoreSetForestWeights,
  bartcoreForestAmplitudes, bartcoreForestFits, bartcoreFitsWithoutOffset,
  bartcoreForestCalibration, bartcoreSetForestPriorScale, bartcoreForestVariableCounts,
  bartcoreSetModel, bartcoreRun, bartcoreSetOffset, bartcoreSetResponse,
  bartcoreSetActiveRows, bartcoreSetWeights, bartcoreSetTestOffset, bartcoreSetData,
  bartcoreSetTestPredictor, bartcoreSetPredictor, bartcoreUpdatePredictor,
  bartcoreUpdatePredictorPerObservation,
  bartcoreUpdatePredictorPerObservationJointly, bartcoreSetCutPoints,
  bartcoreGetLatents, bartcorePredict, bartcorePredictPerForest, bartcoreGetTrees,
  bartcoreStoreState, bartcoreSetState (31 bartcore* wrappers) plus resolveForestIndex,
  which is NOT one of them and STAYS (7 uses in R/dbarts.R).
  Three of the 31 have live in-package callers and so keep their R/bartcore.R
  definitions: bartcoreRun ([[R/bart.R:1488@0045507c]], [[R/bart.R:1574@0045507c]], [[R/bart.R:1826@0045507c]], [[R/bart.R:2076@0045507c]] and [[R/xbart.R:695@0045507c]]),
  bartcorePredict ([[R/generics.R:1221@0045507c]], [[R/generics.R:1354@0045507c]]), bartcoreSetModel ([[R/xbart.R:691@0045507c]]).
  The other 28 are deleted from R/bartcore.R.
  New file inst/common/bartcoreHandle.R defines all 31 names so every test/benchmark
  call site is uniform: 28 verbatim bodies (comments and validation carried over -
  the R-side checks are the tests' own now), and 3 one-line aliases
  `bartcoreRun <- dbarts:::bartcoreRun` (likewise Predict, SetModel) so the two copies
  cannot drift. The moved bodies call exactly three package internals, which stay and
  are reached with ::: - asCountMatrix (1 site), validateCategoryOffset (3 sites),
  resolveForestIndex (1 site). Every .Call target is spelled dbarts:::C_dbarts_*.
  Call sites: 577 `dbarts:::bartcoreX` in inst/tinytest, 137 in benchmarks/R, 18
  `getFromNamespace("bartcoreX", "dbarts")` in benchmarks/R/sbc.R, 2 in
  inst/tinytest/test-fits-without-offset.R. The rewrite is mechanical - strip
  `dbarts:::` (and unwrap getFromNamespace) for the 31 moved names ONLY; leave the
  creators (bartcoreSampler, bartcoreDataHandle, bartcoreSamplerFromHandle,
  bartcoreBCFSampler, bartcoreMultinomial*Sampler) qualified, since they stay in R/.
  48 files gain one source line:
    source(system.file("common", "bartcoreHandle.R", package = "dbarts"))
  - 34 tinytest files (7 of which already source another inst/common helper: the
    directory holds 8 helpers today and no file is auto-sourced, so each file states
    its own) and 14 benchmark scripts (13 with dbarts::: sites plus sbc.R).
  Benchmarks run against an installed library, so system.file resolves there too;
  benchmarks/R/bartcore-shim.R is unrelated (it loads the C++ rshim, not these).
  NAMESPACE is UNCHANGED: none of the 31 was ever exported.

H2. Delete the two dead engine members: [[src/bartcore/tree.hpp:366@0045507c]]
  `Tree::rightChildOf` and [[src/bartcore/sampler.hpp:485@0045507c]]
  `Sampler::setCurrentSampleNum` (one definition each, zero callers in src/, tests/ or
  inst/). Header edits, so R CMD INSTALL --preclean and tests/cpp from make clean.

H3. adoptPointer / reapplyForestWeights: keep BOTH and document both. Evidence for
  reapplyForestWeights (the one J6 left open): it is NOT dead - three live call sites
  ([[R/dbarts.R:1066@0045507c]] in copy(), [[R/dbarts.R:1827@0045507c]] and [[R/dbarts.R:1859@0045507c]] in the getPointer/setState re-creation
  paths) - and it IS now pinned: [[inst/tinytest/test-forest-weights-r5.R:103-140@0045507c]] forces
  forestWeights to list() immediately before each of the three re-creation sites, which
  its own comment at [[inst/tinytest/test-forest-weights-r5.R:108@0045507c]] records as the only oracle that can see the call. What is
  missing is documentation: 0 man/ hits for either name. Add an "Infrastructure
  methods" paragraph to man/dbartsSampler-class.Rd covering adoptPointer,
  reapplyForestWeights and getPointer, lifted from their R5 docstrings ([[R/dbarts.R:946@0045507c]],
  [[R/dbarts.R:1791@0045507c]]). [[inst/tinytest/test-host-shell-pins.R:16-21@0045507c]] already classifies both as
  infrastructure and its census assertions (46 own / 5 infrastructure / 41 substantive,
  [[inst/tinytest/test-host-shell-pins.R:43-45@0045507c]]) are unchanged by this slice - do not perturb them.

H4. Amend [[docs/design/core-generalization.md:69-76@0045507c]]. The sentence "Dispatch is free when
  amortized over the work it gates; nothing dispatches per observation" and the table's
  "Per obs | none: monomorphic loops/kernels" row are both falsified by
  [[src/bartcore/facade.hpp:694-703@0045507c]], the joint per-observation predictor sweep, which
  calls two virtuals per observation through the per-sampler session objects
  (observationWouldRemainValid(j), commitObservation(j)). Amend both to name the joint
  sweep as the deliberate exception and say why: the sweep exists to let N samplers of
  DIFFERENT instantiations vote on one row, so the type erasure is its purpose and the
  per-row vtable hop is the price. docs/** fires no CI (paths-ignore), so this rides
  with the rest of the slice.

Files touched: R/bartcore.R, inst/common/bartcoreHandle.R (new),
src/bartcore/tree.hpp, src/bartcore/sampler.hpp, man/dbartsSampler-class.Rd,
docs/design/core-generalization.md, 34 inst/tinytest files, 14 benchmarks/R files.
No NAMESPACE change, no inst/NEWS.Rd entry (nothing user-visible moves: these were
never exported, and the C++ members were unreachable).

Dense-line budget: code ~+500 in inst/common, ~-470 in R/bartcore.R, -2 in src/;
docs/man ~25; test/benchmark edits ~750 touched lines but all mechanical (714 prefix
strips + 48 source lines + 18 getFromNamespace unwraps).

Gates: full battery, plus --preclean and tests/cpp from clean for H2. The trio is a
REAL gate here, not a formality: benchmarks/R/bcf-equivalence.R (50 bartcore refs) and
benchmarks/R/multinomial-equivalence.R (48) are themselves rewritten by this slice, so
43/12/11 bitwise is what proves the move inert. Also run the full tinytest suite twice
- once to confirm 0 failures, once with the helper deliberately unsourced in one file
to confirm the source lines are load-bearing (a file that silently found the symbol
elsewhere would hide a missed edit).
Anchor drift: land LAST, off whatever bartcore is by then, and re-run the prefix sed
against the tip - every wave-3 test edit and any further SBC work (benchmarks/R/sbc.R
moved 257 lines at 0045507c) lands in files this slice rewrites.

## C. Cross-slice overlap and execution shape

Files two slices both touch, and where:
- inst/NEWS.Rd - r-surface and generics both append to the 1.0-0 BUG FIXES itemize,
  which ends at [[inst/NEWS.Rd:1854@0045507c]] (waves 1 and 2 both appended there: 8042cc2c +34, e35c8797 +8).
  GUARANTEED textual conflict if both are written in parallel. handle-api adds nothing.
- man/bart2.Rd - r-surface edits the usage \dots ([[man/bart2.Rd:68@e35c8797]]), the \item{\dots} ([[man/bart2.Rd:260-262@e35c8797]]),
  the multinomial offset prose ([[man/bart2.Rd:229@e35c8797]]) and the offset.category.test \item ([[man/bart2.Rd:272-273@e35c8797]]);
  generics edits the own-class generic paragraphs ([[man/bart2.Rd:308-320@e35c8797]]) and the summary/vars
  items ([[man/bart2.Rd:279@e35c8797]]). Different regions of one file; git will usually merge, but the file is
  long-line Rd and a hand merge is likely.
- R/generics.R - owned by GENERICS ALONE, provided one item moves: predict.bartMultinomial
  ([[R/generics.R:1013-1070@e35c8797]]) is touched by J5 (rename offset.category.test -> offset) AND by J4 (add
  ci.level, refuse by name). MOVE the rename into the generics slice, together with the
  five predict(offset.category.test=) pins in inst/tinytest/test-multinomial-generics.R.
  r-surface then touches R/generics.R not at all (verified: dataSlotOrNULL has no
  generics.R use; the @offset.category read at [[R/generics.R:1034@e35c8797]] is already a bare slot read and is
  inside predict.bartMultinomial anyway).
- [[R/bart.R:867-874@e35c8797]] - the multinomial offset.test refusal names
  dbarts:::bartcoreSetCategoryTestOffset, which the handle-api slice moves out of the
  namespace; [[man/bart2.Rd:229@e35c8797]] repeats the same claim. RESOLUTION: r-surface rewrites
  both to name the R5 method $setCategoryTestOffset (which stays in the package),
  removing the coupling; handle-api then touches neither file. This is the one item
  that MUST move between slices to keep them disjoint.
- inst/tinytest/test-multinomial-r5-surface.R - r-surface (dbartsData spelling at
  [[inst/tinytest/test-multinomial-r5-surface.R:132@e35c8797]]/[[inst/tinytest/test-multinomial-r5-surface.R:136@e35c8797]]/[[inst/tinytest/test-multinomial-r5-surface.R:140@e35c8797]]/[[inst/tinytest/test-multinomial-r5-surface.R:144@e35c8797]]/[[inst/tinytest/test-multinomial-r5-surface.R:434@e35c8797]]/[[inst/tinytest/test-multinomial-r5-surface.R:449@e35c8797]], dataSlotOrNULL block [[inst/tinytest/test-multinomial-r5-surface.R:469-479@e35c8797]]) and handle-api (it is
  one of the 34 bartcore-calling files). Unavoidable; handle-api lands last and re-runs
  its sed over whatever the earlier slices left.
- inst/tinytest/test-xbart-*.R - r-surface (M10 pins) vs the in-flight xbart-oracle
  worktree, which has test-xbart-reproducibility.R modified and test-xbart-fold-oracle.R
  new. Disjoint files today; r-surface should put its xbart pins in
  test-xbart-error.R / test-error-quality.R and stay out of -reproducibility.R.
- benchmarks/R/sbc.R - handle-api (18 getFromNamespace lines at [[benchmarks/R/sbc.R:1055-1059@0045507c]], [[benchmarks/R/sbc.R:1261-1265@0045507c]],
  [[benchmarks/R/sbc.R:1596@0045507c]]) vs the SBC/leaf-scale arm, which landed 257 lines there at 0045507c and may
  come back. Land handle-api after any further SBC work, or coordinate that file.
- src/ - handle-api touches tree.hpp and sampler.hpp; the engine-default worktree has
  R_interface_bartcore.cpp, chain.hpp, model.hpp modified. Disjoint. E19's R-side fix
  (r-surface, R/model.R defaultNodeScale) is the twin of the engine-default worktree's
  [[R_interface_bartcore.cpp:2298@0045507c]] default arm - no shared file, but they should land in
  either order and then agree; note the pairing in both landing notes.

Recommended execution:
1. r-surface and generics in parallel worktrees off 0045507c. They share only
   inst/NEWS.Rd and man/bart2.Rd after the two reassignments above.
2. Land r-surface FIRST (it is the one with an outward-facing message contract other
   slices quote), then rebase generics onto it and land - the NEWS conflict is then a
   single append-point resolution, as in the residue burn-down's five hand-merged
   append points.
3. handle-api LAST and ALONE: its sed must run after every other test edit exists, and
   its trio run is the gate for the whole move. It is also the only slice that touches
   src/ and so the only one needing --preclean plus tests/cpp from clean.
4. The generics slice cannot start its plot/loglik bodies until the survey agent
   returns SURVEY SLOTS A and B; everything else in that slice (J4, G3, G5, G6, G7)
   is independent and can be written first.

## D. Open questions, with a default so nobody blocks

D1. Do the 714 qualified call sites really "stay untouched"? They cannot: `dbarts:::`
    stops resolving the moment a name leaves the namespace. DEFAULT: a scripted prefix
    strip over the 31 moved names only, reviewed as a pure prefix diff; argument lists,
    ordering and file structure stay byte-identical, which is what the instruction
    protects. If VD meant the stronger thing, the alternative is to keep the 31
    exported-internal and drop the move.
D2. Does "one offset name everywhere" reach the R5 methods $setCategoryOffset /
    $setCategoryTestOffset and their bartcore wrappers? DEFAULT: NO - method names are
    not argument names, their own arguments are already `offset`/`offset.test`, and the
    judgement names only dbartsData's formals and predict's. Merging them into
    $setOffset would be a separate, larger arc.
D3. Are the four creators part of "the 31"? DEFAULT: NO - the 31 is the contiguous
    handle block (E7). bartcoreBCFSampler, bartcoreMultinomialSampler,
    bartcoreMultinomialCountSampler and bartcoreMultinomialDataSampler have no
    in-package callers but are 200+ lines of model construction, not .Call wrappers;
    they stay in R/ and tests keep reaching them with :::.
D4. (g.5) C8, [[R/plotTree.R:9-35@0045507c]]'s dead padding branch, was folded into VD-H by the
    report but is not named in J6. DEFAULT: leave it, record it as residue.
D5. minor [M8-gen], setForestBasis(k, ~var) evaluating the formula in
    environment(basis) instead of against the sampler's data, is named by no judgement.
    DEFAULT: out of wave 3.
D6. as_draws_array/as_draws_df on the four own-class fits: refusing "by name" would
    mean registering methods for a Suggests-package generic. DEFAULT: Rd sentence in
    man/summary.bart.Rd, no registration.
D7. Should the keepSampler-only extract(type = "trees") fallback emit a message as well
    as an Rd sentence? DEFAULT: no - plotTree's identical fallback emits none, and J2
    says "follows plotTree's documented fallback".
D8. J8's shared zero-sample message names bart2/xbart/rbart_vi. bart()'s own ndpost = 0
    behaviour was not probed in this pass. DEFAULT: leave bart() as it is; if the
    implementer's probe shows it faults rather than refusing, add it to the same helper.
D9. J1 deletes the retiredDotsNames sunset mechanism along with the dots. If a future
    rename wants that channel back it must be rebuilt. DEFAULT: accepted, per the rule
    "delete when the only function is a different error message" - record the loss in
    the landing note so it is a decision, not an accident.

## E. Verification of every number above (all run this pass, at 0045507c)

E1. Landed commits and their file lists: `git log --oneline b102e17c..0045507c` (10
    commits) and `git show --stat <sha>` for each; the R-surface diffs read in full for
    8042cc2c, b657e8ae, e35c8797, 07ad73e4. S3 is unlanded: no commit in the range
    touches R/data.R, R/augmentation.R, R/xbart.R or R/model.R
    (`git log --oneline b102e17c..0045507c -- R/data.R R/augmentation.R R/xbart.R R/model.R`
    is empty).
E2. 8 dots pins in inst/tinytest/test-argument-surface.R ([[inst/tinytest/test-argument-surface.R:387@0045507c]], [[inst/tinytest/test-argument-surface.R:389@0045507c]], [[inst/tinytest/test-argument-surface.R:392@0045507c]], [[inst/tinytest/test-argument-surface.R:399@0045507c]], [[inst/tinytest/test-argument-surface.R:403@0045507c]],
    [[inst/tinytest/test-argument-surface.R:549@0045507c]], [[inst/tinytest/test-argument-surface.R:553@0045507c]], [[inst/tinytest/test-argument-surface.R:557@0045507c]]) + 1 in test-heteroscedastic.R ([[inst/tinytest/test-heteroscedastic.R:260@0045507c]]) = 9 expect_error sites on the
    dots channel: `grep -rn "unknown argument" inst/tinytest/*.R` (7 message lines, two
    of them comments) cross-read against the surrounding expect_error blocks.
E3. dots formals: `grep -n "^  \.\.\.$" R/bart.R R/rbart.R` -> [[R/bart.R:704@0045507c]] (bart2),
    [[R/rbart.R:53@0045507c]] (rbart_vi); the other two hits ([[R/bart.R:2487@0045507c]], [[R/bart.R:2539@0045507c]]) are method
    signatures, not entries. Helper: `grep -rn "rejectUnknownDotsArgs|retiredDotsNames"
    R/ man/ inst/tinytest/ docs/` -> [[R/utility.R:120@0045507c]],[[R/utility.R:122@0045507c]],[[R/utility.R:129@0045507c]],[[R/utility.R:130@0045507c]], call sites
    [[R/bart.R:765@0045507c]] and [[R/rbart.R:62@0045507c]], zero man/ hits, one test comment.
E4. bart()'s family vocabulary c("auto","logistic","aft") at [[R/bart.R:2645@0045507c]]; the
    four-token own-class list at [[R/bart.R:2589-2594@0045507c]]; the ten tokens [[bart.Rd:174@0045507c]] names,
    read in full. Existing pins at [[inst/tinytest/test-bart-bart2.R:108-122@0045507c]] (4 tokens).
E5. offset spellings: `grep -rn "offset.category" R/*.R man/*.Rd inst/tinytest/*.R
    benchmarks/R/*.R` -> 17 hits in R/data.R (formals at [[R/data.R:899-900@0045507c]], missing() clauses at
    [[R/data.R:907-909@0045507c]], validation at [[R/data.R:1430-1443@0045507c]]), R/A_class.R slots at [[R/A_class.R:536-538@0045507c]],
    [[R/generics.R:1017@0045507c]]/[[R/generics.R:1033-1056@0045507c]] (predict.bartMultinomial), [[man/bart2.Rd:80@0045507c]]/[[man/bart2.Rd:229@0045507c]]/
    [[man/bart2.Rd:272-273@0045507c]]/[[man/bart2.Rd:312@0045507c]], [[man/dbartsData.Rd:14@0045507c]]/[[man/dbartsData.Rd:26-27@0045507c]], [[man/dbarts.Rd:103@0045507c]],
    [[man/dbartsSampler-class.Rd:142@0045507c]]/[[man/dbartsSampler-class.Rd:150@0045507c]]; tests: test-multinomial-r5-surface.R 21 hits
    (6 of them dbartsData calls) and test-multinomial-generics.R 5 (all
    predict(offset.category.test=)); benchmarks/: ZERO
    (`grep -rc "offset.category" benchmarks/R/*.R` all 0).
E6. M1's two sites: [[R/data.R:1178-1180@0045507c]] and [[R/data.R:1237-1239@0045507c]], both `stop("'x' must have the
    same number of observations as 'y'")` after codeResponse ([[R/data.R:458-473@0045507c]]).
E7. The 31 handle wrappers: `awk 'NR>=1066 && NR<=1535' R/bartcore.R | grep -cE
    "^[a-zA-Z][A-Za-z0-9_.]* <- function"` -> 32 definitions, of which resolveForestIndex
    ([[R/bartcore.R:1252@0045507c]]) is a helper -> 31 bartcore* wrappers. Per-function caller census run over
    R/*.R, inst/tinytest/*.R and benchmarks/: bartcoreRun 5 in-package call sites
    ([[R/bart.R:1488@0045507c]],[[R/bart.R:1574@0045507c]],[[R/bart.R:1826@0045507c]],[[R/bart.R:2076@0045507c]]; [[R/xbart.R:695@0045507c]]), bartcorePredict 2
    ([[R/generics.R:1221@0045507c]],[[R/generics.R:1354@0045507c]]), bartcoreSetModel 1 ([[R/xbart.R:691@0045507c]]); the other 28 have
    zero in-package callers (the only other R/ hits are comments, listed by
    `grep -rn "\bbartcoreX\b" R/*.R | grep -v "^R/bartcore.R"`). bartcorePredictPerForest
    has zero callers anywhere.
E8. Call-site counts: `grep -rhoE "dbarts:::bartcore[A-Za-z0-9_]+" inst/tinytest/*.R |
    wc -l` -> 577; the same over benchmarks/R/*.R -> 137;
    `grep -rhoE "getFromNamespace\(\"bartcore[A-Za-z0-9_]+\"" benchmarks/R/*.R | wc -l`
    -> 18, all in sbc.R; 2 more in inst/tinytest/test-fits-without-offset.R. Files:
    `grep -rlE "\bbartcore[A-Za-z0-9_]+" inst/tinytest/*.R | wc -l` -> 34;
    `grep -rlE "dbarts:::bartcore" benchmarks/R/*.R | wc -l` -> 13 (+ sbc.R = 14).
    Of the 34 test files, 7 already source an inst/common helper
    (`comm -12` of the two file lists); inst/common holds 8 helpers today (`ls inst/common`).
E9. Trio composition: equivalence 43 scenarios, bcf-equivalence 12,
    multinomial-equivalence 11 ([[docs/plans/review-2026-08-24/gate-ledger.md:111@658869ac]]).
    bcf-equivalence.R holds 50 bartcore refs, multinomial-equivalence.R 48,
    equivalence.R 0 (`grep -ohE "dbarts:::bartcore[A-Za-z0-9_]+" <file> | sort | uniq -c`).
E10. Dead engine members: `grep -rn "rightChildOf|setCurrentSampleNum" src/ tests/ inst/
    docs/` -> [[src/bartcore/tree.hpp:366@0045507c]], [[src/bartcore/sampler.hpp:485@0045507c]], plus one
    docs/plans note. Zero call sites.
E11. reapplyForestWeights: 3 call sites ([[R/dbarts.R:1066@0045507c]], [[R/dbarts.R:1827@0045507c]], [[R/dbarts.R:1859@0045507c]]), definition at
    [[R/dbarts.R:1790@0045507c]], pin at [[inst/tinytest/test-forest-weights-r5.R:103-140@0045507c]] with the oracle comment
    at [[inst/tinytest/test-forest-weights-r5.R:108@0045507c]]; adoptPointer: definition [[R/dbarts.R:945@0045507c]], 2 call sites ([[R/bart.R:1807@0045507c]],
    [[R/bart.R:2056@0045507c]]); both have 0 man/ hits (`grep -rn "adoptPointer|reapplyForestWeights" man/`
    is empty). [[test-host-shell-pins.R:16-21@0045507c]] lists both as infrastructure; its counts are
    at [[test-host-shell-pins.R:43-45@0045507c]] (46 / 5 / 41).
E12. core-generalization.md's claim is at [[docs/design/core-generalization.md:69-71@0045507c]] plus the "Per obs" table row at [[docs/design/core-generalization.md:76@0045507c]];
    the falsifying loop is [[src/bartcore/facade.hpp:694-703@0045507c]] (two virtuals per
    observation).
E13. n.samples = 0: three refusal sites [[R/bart.R:806-809@0045507c]], [[R/xbart.R:94-96@0045507c]],
    [[R/rbart.R:104-106@0045507c]]; two pins ([[inst/tinytest/test-error-quality.R:46@0045507c]],
    [[test-xbart-error.R:24@0045507c]]); zero pins on rbart_vi's wording
    (`grep -rn "no posterior draws" inst/tinytest/ benchmarks/` is empty);
    [[man/dbartsControl.Rd:32-33@0045507c]] is the sentence to extend.
E14. monotone: parseMonotoneSign at [[R/model.R:548-570@0045507c]] with the "inc"/"dec"/"0" arms at
    [[R/model.R:553@0045507c]],[[R/model.R:557@0045507c]],[[R/model.R:559@0045507c]]; zero tests pin those spellings
    (`grep -rn '"inc"|"dec"' inst/tinytest/*.R` is empty); [[man/dbarts.Rd:72@0045507c]] and
    [[man/bart2.Rd:203@0045507c]] are the Rd items.
E15. The 12 "should be one of" pins that must survive G5:
    `grep -rn "should be one of" inst/tinytest/*.R` -> [[test-augmentation.R:376@0045507c]],
    [[test-argument-surface.R:259@0045507c]],[[test-argument-surface.R:260@0045507c]], [[test-bart-bart2.R:109@0045507c]] (a comment),
    [[test-hazard.R:339@0045507c]],[[test-hazard.R:346@0045507c]], [[test-hurdle-surface.R:77@0045507c]],[[test-hurdle-surface.R:85@0045507c]], [[test-nbinom.R:232@0045507c]],[[test-nbinom.R:236@0045507c]],
    [[test-prior-predictive.R:127@0045507c]], [[test-xbart-model.R:206@0045507c]] - none of them a `sample`
    argument, so the wording change is safe.
E16. NEWS append point: `grep -n "subsection" inst/NEWS.Rd` -> the 1.0-0 BUG FIXES
    subsection opens at [[inst/NEWS.Rd:995@0045507c]] and its itemize closes at [[inst/NEWS.Rd:1854@0045507c]]; both landed waves appended
    immediately above that line.
E17. In-flight worktrees: `git worktree list` -> engine-default and xbart-oracle, both
    at 0045507c; `git status --porcelain` in each -> engine-default has
    src/R_interface_bartcore.cpp, src/bartcore/chain.hpp, src/bartcore/model.hpp,
    tests/cpp/Makefile modified; xbart-oracle has inst/tinytest/test-xbart-
    reproducibility.R modified and test-xbart-fold-oracle.R new.
E18. Suite size: 167 files in inst/tinytest (`ls inst/tinytest/*.R | wc -l`); wave 2's
    recorded assertion count is 7040 (landing note [[docs/plans/release-candidate-review.md:619@0045507c]]).
