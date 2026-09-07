# Fix wave 3 - implementer spec (bartcore @ 0045507c)

Source: consolidated-report.md sections (b)(c)(d)(e)(g), the
per-cell evidence in matrix-results.md / matrix-review-entries.md / matrix-review-generics.md,
and VD's judgements J1-J8. Waves 1-2 are b102e17c..07ad73e4 (landing note at
[docs/plans/release-candidate-review.md:543](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/docs/plans/release-candidate-review.md#L543)). J7 (engine default: arms, in a
private worktree) and J9 (predict n.threads) are OUT of this wave.

Every count and line reference below was re-verified in this pass; section E lists the
command behind each. All line numbers are against 0045507c.

## A. Landed-vs-remaining ledger (R surface)

LANDED (do not re-plan)
- 8042cc2c: B1 extract(type="trees") refuses sample/combineChains/forest/contribution by
  name ([R/generics.R:360-372](https://github.com/vdorie/dbarts/blob/8042cc2c8a1d5afd86fbc59c6cd38b101bcae5b8/R/generics.R#L360-L372), wired at [R/generics.R:398](https://github.com/vdorie/dbarts/blob/8042cc2c8a1d5afd86fbc59c6cd38b101bcae5b8/R/generics.R#L398) and [R/generics.R:1901](https://github.com/vdorie/dbarts/blob/8042cc2c8a1d5afd86fbc59c6cd38b101bcae5b8/R/generics.R#L1901)); B2 survivalProbabilities names the
  hazard period column on both branches ([R/bart.R:2421](https://github.com/vdorie/dbarts/blob/8042cc2c8a1d5afd86fbc59c6cd38b101bcae5b8/R/bart.R#L2421), [R/bart.R:2434](https://github.com/vdorie/dbarts/blob/8042cc2c8a1d5afd86fbc59c6cd38b101bcae5b8/R/bart.R#L2434)); minor [M2] residuals(sample=)
  refused by name (refuseResidualsSample, [R/bart.R:842-851](https://github.com/vdorie/dbarts/blob/8042cc2c8a1d5afd86fbc59c6cd38b101bcae5b8/R/bart.R#L842-L851), on bart/rbart/bartHurdle).
- b657e8ae: M3 ([man/dbartsSpec.Rd:40](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/man/dbartsSpec.Rd#L40),[man/dbartsSpec.Rd:48](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/man/dbartsSpec.Rd#L48)), M4 ([man/bart.Rd:153](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/man/bart.Rd#L153)), M5 ([man/bart.Rd:165](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/man/bart.Rd#L165)), M7 ([man/bart.Rd:248](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/man/bart.Rd#L248)
  hurdle storeState recipe), M8 ([man/dbartsSampler-class.Rd:328](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/man/dbartsSampler-class.Rd#L328)), M12 (three chain.hpp comments),
  U7 ([docs/design/kernel-vocabulary.md:26](https://github.com/vdorie/dbarts/blob/b657e8ae41c308977f9550e80b9e109fae3aa0d1/docs/design/kernel-vocabulary.md#L26)), minors [M4] setForestBasis wording, [M5] chain column, [M7]
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
- minor [M6-gen] [man/bart.Rd:248](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart.Rd#L248) understates which generics stop on a stateless fit ->
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
  ("twopart" -> "hurdle.lognormal", [R/bart.R:712](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L712), [R/dbarts.R:405](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/dbarts.R#L405)) or remaps a family
  ("hazard"/"hazard.logistic" -> "probit"/"logistic", [R/dbarts.R:534](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/dbarts.R#L534)), capture
  `requestedFamily <- family` BEFORE the fold and interpolate that in every downstream
  stop(); never the resolved token.
- Name the KIND received, not just the expected shape ("a length-n vector was
  supplied", "an n x 2 matrix was supplied", "a Surv response was supplied").
- Single quotes around argument names, escaped double quotes around family tokens -
  the package's existing style ([R/bart.R:849](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L849), [R/data.R:628](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L628)).

### (i) Slice "r-surface" - J1 + J3 + J5 + J8 + remaining S3 + dataSlotOrNULL

R1. J1 - delete the rejection-only dots channel.
  Current: bart2 ([R/bart.R:704](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L704) `...`) and rbart_vi ([R/rbart.R:53](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/rbart.R#L53) `...`) declare a dots
  formal used only to produce a nicer error; [R/bart.R:764-765](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L764-L765) and [R/rbart.R:61-62](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/rbart.R#L61-L62) call
  rejectUnknownDotsArgs(argNames, fn) ([R/utility.R:122-142](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/utility.R#L122-L142)), which consults
  retiredDotsNames ([R/utility.R:120](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/utility.R#L120), one row: rngSeed -> seed) and otherwise agreps the
  nearest formal.
  Required: delete the `...` formal from both signatures, delete the two
  argNames/rejectUnknownDotsArgs blocks ([R/bart.R:762-765](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L762-L765), [R/rbart.R:59-62](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/rbart.R#L59-L62)), delete
  rejectUnknownDotsArgs and retiredDotsNames and the 12-line comment above them
  ([R/utility.R:108-142](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/utility.R#L108-L142)). R's own "unused argument" wall then serves all six entries
  (bart, bart2, dbarts, dbartsSpec, xbart, rbart_vi) identically - which also closes
  E15's asymmetry without writing anything. Partial matching is unaffected: both dots
  were the LAST formal, so no formal moves behind a `...`.
  Rd: delete `\dots)` from the bart2 usage block ([man/bart2.Rd:68](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L68)) and the rbart_vi
  usage block ([man/rbart.Rd:38](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/rbart.Rd#L38)); delete the \item{\dots} paragraph at [man/bart2.Rd:260-262](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L260-L262);
  drop `\dots` from the shared \item list at [man/rbart.Rd:85](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/rbart.Rd#L85). Leave every method usage's
  \dots alone (predict/fitted/extract/plot keep theirs).

R2. J3 - bart()'s by-name family redirects, and the typed-token rule.
  Current: [R/bart.R:2648-2655](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L2648-L2655) refuses only bartOwnClassFamilies ([R/bart.R:2589-2594](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L2589-L2594):
  multinomial, ordinal, nbinom, hurdle.lognormal) through refuseBartOwnClassFamily
  ([R/bart.R:2596-2607](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L2596-L2607)); the other six of the ten tokens [man/bart.Rd:174](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart.Rd#L174) names fall through to
  match.arg's "'arg' should be one of "auto", "logistic", "aft"".
  Required: extend the pre-match.arg by-name branch to all ten tokens with three
  reasons, exactly as [man/bart.Rd:174](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart.Rd#L174) already states them:
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
  token before [R/dbarts.R:405](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/dbarts.R#L405)'s twopart fold and before [R/dbarts.R:534](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/dbarts.R#L534)'s hazard remap and use it in the
  downstream refusals that currently name the resolved one - the sites reached are
  R/model.R's resid.dist refusal ("family \"probit\" has its own fixed error scale")
  and the variance-forest refusal. Do the same in [R/bart.R:712](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L712) for bart2. Rule to
  write in the comment: the resolved token is an implementation detail; the caller can
  only act on what they typed.
  rbart_vi / xbart / dbartsSpec keep their narrow vocabularies ([R/rbart.R:52](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/rbart.R#L52)
  c("auto","gaussian","aft"), [R/xbart.R:26](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/xbart.R#L26) c("auto","gaussian","probit","logistic"),
  R/spec.R's eight). Their Rd sentences are corrected to say the vocabulary is
  narrower BY DESIGN and that the wider set lives on bart2 - one sentence each in
  man/rbart.Rd, man/xbart.Rd, man/dbartsSpec.Rd. E18 (bart(keepevery = -1) refused as
  'n.thin') is the same rule one layer down: bart() already coerces by the typed name
  ([R/bart.R:2670-2677](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L2670-L2677)) but the positivity refusal comes from dbartsControl's validity,
  which names n.thin - add a bart()-side positivity check for keepevery/ndpost/nskip
  naming the typed spelling before the control is built.

R3. J5 - one offset name, shape follows family.
  Current spellings (E8, all re-verified): bart2(offset = <n x K>) accepted;
  bart2(offset.category = ) -> unknown argument; dbarts(offset.category = ) -> raw
  unused argument; dbarts(offset = <n x K>) -> "'offset' must have the same length as
  'y'" ([R/data.R:765](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L765), a length() test, so an n x 1 matrix passes and an n x K does
  not); dbartsData(counts = , offset.category = ) accepted ([R/data.R:899-900](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L899-L900)).
  Required:
  a. dbartsData: DELETE the offset.category and offset.category.test formals
     ([R/data.R:899-900](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L899-L900)) and the countsIsMissing clauses that read them ([R/data.R:907-909](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L907-L909)).
     The n x K train shift arrives as `offset`, the nTest x K test shift as
     `offset.test`. The dbartsData SLOTS keep their names ([R/A_class.R:536-538](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/A_class.R#L536-L538)) and so
     do the R5 methods $setCategoryOffset/$setCategoryTestOffset (method names, not
     argument names - see D2).
  b. Routing: in dbartsData and in the dbarts()/bart2() matrix paths, a matrix-valued
     `offset` on a counts-carrying data object installs the category offset
     (validateDataCategoryOffset, [R/data.R:1430-1443](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L1430-L1443), unchanged); a matrix-valued
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
     ([R/generics.R:1017](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/generics.R#L1017), body reads at [R/generics.R:1033-1056](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/generics.R#L1033-L1056)). SPECIFIED HERE, EXECUTED BY THE
     GENERICS SLICE (section C: J4 rewrites the same method; one owner per function).
     predict.bart already spells the
     new-row shift `offset` ([R/generics.R:210](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/generics.R#L210)), so this is the same name for the same
     thing on both classes. Update [man/bart2.Rd:80](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L80) (usage) and [man/bart2.Rd:272-273](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L272-L273) (the \item),
     and the two prose paragraphs that name the old spelling ([man/bart2.Rd:229](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L229), [man/bart2.Rd:312](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L312)).
  e. M1's false message: [R/data.R:1178-1180](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L1178-L1180) (sparse branch) and [R/data.R:1237-1239](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L1237-L1239) (numeric
     branch) compare NROW(formula) against NROW(codeResponse(data)$y) after
     codeResponse ([R/data.R:458-473](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L458-L473)) has flattened an n x 2 response to length 2n.
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
  f. [R/bart.R:867-874](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L867-L874)'s multinomial offset.test refusal names
     dbarts:::bartcoreSetCategoryTestOffset as the internal channel, and
     [man/bart2.Rd:229](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L229) repeats it. Under the handle-api slice that symbol leaves the
     namespace; see C for who owns the edit.

R4. J8 - option vocabularies.
  a. monotone: [R/model.R:548-570](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/model.R#L548-L570) parseMonotoneSign switches on tolower(value) and
     accepts "inc", "dec", "0" alongside the documented set. DELETE those three switch
     arms; keep tolower() and DOCUMENT the case fold in [man/dbarts.Rd:72](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/dbarts.Rd#L72) and
     [man/bart2.Rd:203](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart2.Rd#L203) ("matching is case-insensitive, so "Increasing" is accepted").
     The numeric 0 arm is untouched - [man/dbarts.Rd:72](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/dbarts.Rd#L72) documents 0 for unconstrained.
     No existing pin asserts "inc"/"dec"/"0" (verified: 0 hits in inst/tinytest).
  b. makeind(all = ): NO CODE CHANGE. [R/bart.R:2835-2838](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L2835-L2838) binds `ignored <- all`;
     [man/makeind.Rd:26](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/makeind.Rd#L26) already says "Not currently implemented". Add the reason to
     that item - "retained for signature compatibility with BayesTree::makeind" - and
     stop there.
  c. n.samples = 0: dbartsControl/dbarts keep accepting it (the host-loop shape).
     bart2 ([R/bart.R:806-809](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L806-L809)), xbart ([R/xbart.R:94-96](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/xbart.R#L94-L96)) and rbart_vi ([R/rbart.R:104-106](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/rbart.R#L104-L106),
     currently "no posterior draws will be taken after thinning") share ONE message
     from one helper in R/utility.R:
       ''n.samples' must leave at least one draw after thinning (n.samples %/% n.thin
       = 0); dbarts() and dbartsControl() accept a zero-draw run - a sampler driven by
       a host loop - but <caller>() returns posterior draws'
     with <caller> the entry point's own name. Add one sentence to
     [man/dbartsControl.Rd:32-33](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/dbartsControl.Rd#L32-L33) recording the split (it already contrasts the
     per-run() return count with bart2's sweep budget - this appends the zero case).
     The multinomial branch keeps its own family-named check ([R/bart.R:923](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L923)).

R5. Remaining S3 items.
  - M6: [R/bart.R:389-393](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/bart.R#L389-L393) sets result$n.chains only in the else branch. Set it
    UNCONDITIONALLY (2 lines), as bartMultinomial/bartOrdinal/bartNegbin already do;
    bartHurdle's assembly needs the same. [man/bart.Rd:315-317](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart.Rd#L315-L317) then becomes true.
  - M9: [R/augmentation.R:65](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/augmentation.R#L65) declares sigma = 1 and [R/augmentation.R:79-81](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/augmentation.R#L79-L81) guards on !missing(sigma).
    Change the default to NULL and guard on !is.null(sigma); the aft/student arms that
    consume it already treat NULL as absent ([R/augmentation.R:52](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/augmentation.R#L52)'s augRestrict).
  - M10: (i) move `family <- match.arg(family)` ([R/xbart.R:123](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/xbart.R#L123)) ABOVE the data build
    so a bad family is named before the response is ingested, matching the other four
    entries; (ii) give n.threads the length/positivity check the others have -
    [R/xbart.R:398-401](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/xbart.R#L398-L401) currently coerces then tests is.na() on a length-2 value, which
    raises R's raw "'length = 2' in coercion to 'logical(1)'". Add
    'n.threads' must be of length 1' ([R/A_class.R:296](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/A_class.R#L296)'s wording) before the positivity
    test.
  - E19: [R/model.R:400-415](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/model.R#L400-L415) defaultNodeScale's switch() has no default arm, so
    ("hazard") and ("student") return NULL silently. Add the sibling's stop() text
    ([R/model.R:451](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/model.R#L451)): 'no node scale is defined for family "<family>"'.
  - minor [M6-gen]: [man/bart.Rd:248](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/man/bart.Rd#L248)'s Saving subsection says only predict stops on a
    stateless fit; extract(type = "trees") and plotTree stop identically. One clause.
  - minor [M10-gen]: R/bart.R's result assembly writes yhat.test/s.train/s.test
    entries as NULL, so names(fit) lists absent components. Assign only when non-NULL
    (result[["yhat.test"]] <- ... inside the existing if), or drop the NULL entries at
    the end with result[!vapply(result, is.null, NA)] - the second is one line and
    order-preserving.
  - dataSlotOrNULL deletion (from J6, homed here): [R/data.R:11-13](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L11-L13) with its 8-line
    comment (:4-10). Its comment claims every internal read of counts/offset.category/
    offset.category.test goes through it; FALSE - it has exactly ONE in-package use
    ([R/data.R:20](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/data.R#L20), inside dataCounts), while @counts is read bare at [R/A_class.R:634](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/A_class.R#L634) and
    @offset.category bare at [R/generics.R:1034](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/R/generics.R#L1034). Inline the slot read into dataCounts
    and delete the function. It protects only objects serialized by an intermediate
    commit of THIS branch, which nothing outside the branch holds.
    Test consequence: [inst/tinytest/test-multinomial-r5-surface.R:469-479](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L469-L479) exercises it
    by hand-stripping attributes - delete that block.

Files touched: R/bart.R, R/rbart.R, R/xbart.R, R/data.R, R/dbarts.R, R/model.R,
R/augmentation.R, R/utility.R, R/spec.R (family Rd cross-reference only if the
narrow-vocabulary sentence needs a code comment); man/bart.Rd, man/bart2.Rd,
man/rbart.Rd, man/xbart.Rd, man/dbarts.Rd, man/dbartsData.Rd, man/dbartsControl.Rd,
man/dbartsSpec.Rd, man/makeind.Rd, man/dbartsAugmentation.Rd; inst/NEWS.Rd.

Test pins to add or rewrite (inst/tinytest):
- test-argument-surface.R - REWRITE the 8 dots pins ([inst/tinytest/test-argument-surface.R:387](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L387), [inst/tinytest/test-argument-surface.R:389](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L389), [inst/tinytest/test-argument-surface.R:392](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L392), [inst/tinytest/test-argument-surface.R:399](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L399), [inst/tinytest/test-argument-surface.R:403](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L403),
  [inst/tinytest/test-argument-surface.R:549](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L549), [inst/tinytest/test-argument-surface.R:553](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L553), [inst/tinytest/test-argument-surface.R:557](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L557)) to expect R's "unused argument"; the rngSeed pair loses its
  retirement text entirely (that is the J1 trade). [inst/tinytest/test-argument-surface.R:563](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L563)'s dbarts() "unused argument"
  pin is unchanged and becomes the shared expectation. [inst/tinytest/test-argument-surface.R:407-412](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-argument-surface.R#L407-L412)'s partial-match pin
  must stay green (it will).
- [inst/tinytest/test-heteroscedastic.R:259-262](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-heteroscedastic.R#L259-L262) - the rbart_vi(variance = ) pin moves from
  "unknown argument" to "unused argument".
- [inst/tinytest/test-bart-bart2.R:108-122](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-bart-bart2.R#L108-L122) - ADD six by-name refusal pins (gaussian, probit,
  hazard.probit, hazard, hazard.logistic, twopart), each asserting the typed token
  appears in the message; keep the four existing ones.
- [inst/tinytest/test-error-quality.R:38-47](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-error-quality.R#L38-L47) - the n.samples = 0 pin text changes to the shared
  message; ADD the rbart_vi and xbart arms so all three entries are pinned to ONE
  string. [inst/tinytest/test-xbart-error.R:24](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-xbart-error.R#L24) - same rewrite.
- test-data-compatibility.R or a new test-response-shape.R - ADD the three by-kind
  ingest refusals (Surv, n x 2, n x K), each asserting the message names the kind;
  plus the four inheriting surfaces (dbartsData positional, rbart_vi matrix,
  xbart, dbarts) reaching the same text.
- test-multinomial-r5-surface.R - REWRITE [inst/tinytest/test-multinomial-r5-surface.R:132](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L132), [inst/tinytest/test-multinomial-r5-surface.R:136](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L136), [inst/tinytest/test-multinomial-r5-surface.R:140](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L140), [inst/tinytest/test-multinomial-r5-surface.R:144](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L144), [inst/tinytest/test-multinomial-r5-surface.R:434](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L434), [inst/tinytest/test-multinomial-r5-surface.R:449](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L449) from
  offset.category=/offset.category.test= to offset=/offset.test=; DELETE [inst/tinytest/test-multinomial-r5-surface.R:469-479](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-r5-surface.R#L469-L479)
  (dataSlotOrNULL).
- test-multinomial-generics.R - REWRITE the 5 offset.category.test= sites ([inst/tinytest/test-multinomial-generics.R:288](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-generics.R#L288), [inst/tinytest/test-multinomial-generics.R:295](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-generics.R#L295),
  [inst/tinytest/test-multinomial-generics.R:302](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-generics.R#L302), [inst/tinytest/test-multinomial-generics.R:310](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-generics.R#L310), [inst/tinytest/test-multinomial-generics.R:321](https://github.com/vdorie/dbarts/blob/07ad73e403d77630c801a6705cf3afdd63f9bdd1/inst/tinytest/test-multinomial-generics.R#L321)) to offset=. Ships WITH the generics slice, per R3(d).
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
out of that file) and any further wave-3 slice appends at [inst/NEWS.Rd:1854](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L1854).

### (ii) Slice "generics" - J2 + J4 + trees fallback + the three Rd sentences

G1. J2 - plot on the four own-class families. Current: only plot.bartMultinomial
  exists ([R/plot.R:188-210](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/plot.R#L188-L210), registered NAMESPACE S3method(plot, bartMultinomial));
  plot(bartOrdinalFit) reaches plot.default and raises "'x' is a list, but does not
  have components 'x' and 'y'". Required: plot.bartOrdinal, plot.bartNegbin,
  plot.bartHurdle in R/plot.R + three NAMESPACE registrations.
  >>> SURVEY SLOT A: the per-family plot semantics (what is traced, and against what)
  come from the generics-survey agent; plot.bartMultinomial's "trace each category's
  training-mean predicted probability over the kept draws" ([man/bart2.Rd:312](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart2.Rd#L312)) is the
  shape to match. Implement whatever the survey returns; the file, registration and
  Rd obligations below are fixed regardless.

G2. J2 - extract(type = "loglik") on the four own-class families. Current: all four
  refuse by name through validateType ([R/generics.R:1442-1448](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1442-L1448)) because "loglik" is
  absent from their type vocabularies ([R/generics.R:906](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L906) multinomial, [R/generics.R:1115](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1115) ordinal,
  [R/generics.R:1281](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1281) negbin, [R/generics.R:1564](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1564) hurdle); [man/bart.Rd:201](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart.Rd#L201) documents "loglik" on the extract generic
  unscoped, and it is the channel loo/WAIC consume. Required: add "loglik" to each of
  the four vocabularies and give pointwiseLogLikelihood ([R/generics.R:56-169](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L56-L169)) an arm
  per family; it already switches on object[["family"]] and already handles the
  weights and chain-shape bookkeeping.
  >>> SURVEY SLOT B: the per-family density (multinomial: the count/one-hot
  multinomial log-pmf at the reported probabilities; ordinal: the observed category's
  cumulative-probit difference - the engine's own per-observation channel, model.hpp:
  3314, is the oracle; nbinom: dnbinom at the per-draw dispersion; hurdle: the
  two-part mixture, the zero spike plus the lognormal branch) is the survey's to fix.
  Whatever it returns must satisfy: shape identical to extract(type = "ev") minus any
  K margin, zero-weight rows flagged NaN as the gaussian arm does ([R/generics.R:122](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L122)).

G3. J2 - extract(type = "trees") on a keepSampler-only fit. Current: the guard at
  [R/generics.R:385-396](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L385-L396) tests is.null(object$fit), i.e. the SAMPLER, so a
  keepTrees = FALSE / keepSampler = TRUE fit falls through and getTrees returns the
  CURRENT working trees as an 11 x 4 frame (tree, n, var, value - no chain, no sample
  column, no warning). Required: KEEP that behavior (it is plotTree's documented
  fallback, [man/plotTree.Rd:39-42](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/plotTree.Rd#L39-L42)) and DISCLOSE it in man/bart.Rd's "Extracting Trees"
  subsection ([man/bart.Rd:257-266](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart.Rd#L257-L266)): state that without keepTrees but with a kept sampler the
  frame holds the sampler's CURRENT trees and omits the chain/sample index columns,
  in the same words [man/plotTree.Rd:39-42](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/plotTree.Rd#L39-L42) uses. No message is emitted (plotTree emits
  none). Pin both shapes.

G4. J4 - the own-class argument vocabulary. Current formals are
  (object, type, sample, ...) on extract.bartMultinomial ([R/generics.R:904](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L904)), .bartOrdinal ([R/generics.R:1113](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1113)),
  .bartNegbin ([R/generics.R:1279](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1279)), and the bare "..." silently swallows the bart-family vocabulary:
  ten verified cells, each identical() to the call without the argument -
  extract(combineChains=), extract(forest=), fitted(sample=), fitted(ci.level=),
  predict(ci.level=), residuals(type=), summary(vars=).
  Required:
  a. HONOUR combineChains on extract for the three classes: add the formal, reshape
     through combineOrUncombineChains(x, fitNChains(object), combineChains)
     ([R/generics.R:196-205](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L196-L205)), the helper predict on the same fits already uses.
  b. HONOUR ci.level on fitted and predict for the three classes, through
     posteriorInterval ([R/generics.R:170-195](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L170-L195)), returning est/ci.lower/ci.upper exactly
     as bart/rbart/bartHurdle do ([man/bart.Rd:215-216](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart.Rd#L215-L216)).
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
     summary.bartOrdinal/.bartNegbin already carry a real vars ([R/diagnostics.R:229-247](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/diagnostics.R#L229-L247))
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
  match.arg ([R/generics.R:911](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L911), [R/generics.R:1120](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1120), [R/generics.R:1286](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1286), [R/generics.R:1570](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1570)) and give "'arg' should be one of
  ...", where bart/rbart say "sample must be in 'train', 'test'" ([R/generics.R:415](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L415),
  [R/generics.R:815](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L815), [R/generics.R:1947](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1947), [R/generics.R:2065](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L2065)). Add validateSample() beside validateType and use it at all
  eight sites, so one wording serves every class.

G6. minor [M3-gen] - plotTree.bart ([R/generics.R:2134-2149](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L2134-L2149)) builds args via
  list(treeNum = treeNum, ...) and do.call, so a caller's sample=/chain= partial-
  matches sampleNum/chainNum and silently draws. Refuse sample/chain by name in the
  same helper: "'sample' is not used by plotTree; the saved sample is 'sampleNum'".
  plotTree.rbart ([R/generics.R:2151](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L2151)) takes the same treatment.

G7. J2 - the three unclassed results, one Rd sentence each.
  man/pdbart.Rd \value ([man/pdbart.Rd:82](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/pdbart.Rd#L82)): "The result carries class \"pdbart\"/\"pd2bart\" for its
  plot method only; it is not a fit, so predict, extract, fitted and residuals are not
  defined for it (fitted and residuals fall through to stats' defaults and return
  NULL)."
  man/xbart.Rd \value ([man/xbart.Rd:129-133](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/xbart.Rd#L129-L133)): "The result is a bare array with no class, so the fit
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
  The set is exactly the contiguous block [R/bartcore.R:1066-1535](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bartcore.R#L1066-L1535), which holds 32
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
  definitions: bartcoreRun ([R/bart.R:1488](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1488), [R/bart.R:1574](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1574), [R/bart.R:1826](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1826), [R/bart.R:2076](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2076) and [R/xbart.R:695](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/xbart.R#L695)),
  bartcorePredict ([R/generics.R:1221](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1221), [R/generics.R:1354](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1354)), bartcoreSetModel ([R/xbart.R:691](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/xbart.R#L691)).
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

H2. Delete the two dead engine members: [src/bartcore/tree.hpp:366](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/tree.hpp#L366)
  `Tree::rightChildOf` and [src/bartcore/sampler.hpp:485](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/sampler.hpp#L485)
  `Sampler::setCurrentSampleNum` (one definition each, zero callers in src/, tests/ or
  inst/). Header edits, so R CMD INSTALL --preclean and tests/cpp from make clean.

H3. adoptPointer / reapplyForestWeights: keep BOTH and document both. Evidence for
  reapplyForestWeights (the one J6 left open): it is NOT dead - three live call sites
  ([R/dbarts.R:1066](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1066) in copy(), [R/dbarts.R:1827](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1827) and [R/dbarts.R:1859](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1859) in the getPointer/setState re-creation
  paths) - and it IS now pinned: [inst/tinytest/test-forest-weights-r5.R:103-140](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-forest-weights-r5.R#L103-L140) forces
  forestWeights to list() immediately before each of the three re-creation sites, which
  its own comment at [inst/tinytest/test-forest-weights-r5.R:108](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-forest-weights-r5.R#L108) records as the only oracle that can see the call. What is
  missing is documentation: 0 man/ hits for either name. Add an "Infrastructure
  methods" paragraph to man/dbartsSampler-class.Rd covering adoptPointer,
  reapplyForestWeights and getPointer, lifted from their R5 docstrings ([R/dbarts.R:946](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L946),
  [R/dbarts.R:1791](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1791)). [inst/tinytest/test-host-shell-pins.R:16-21](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-host-shell-pins.R#L16-L21) already classifies both as
  infrastructure and its census assertions (46 own / 5 infrastructure / 41 substantive,
  [inst/tinytest/test-host-shell-pins.R:43-45](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-host-shell-pins.R#L43-L45)) are unchanged by this slice - do not perturb them.

H4. Amend [docs/design/core-generalization.md:69-76](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/core-generalization.md#L69-L76). The sentence "Dispatch is free when
  amortized over the work it gates; nothing dispatches per observation" and the table's
  "Per obs | none: monomorphic loops/kernels" row are both falsified by
  [src/bartcore/facade.hpp:694-703](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/facade.hpp#L694-L703), the joint per-observation predictor sweep, which
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
  which ends at [inst/NEWS.Rd:1854](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L1854) (waves 1 and 2 both appended there: 8042cc2c +34, e35c8797 +8).
  GUARANTEED textual conflict if both are written in parallel. handle-api adds nothing.
- man/bart2.Rd - r-surface edits the usage \dots ([man/bart2.Rd:68](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L68)), the \item{\dots} ([man/bart2.Rd:260-262](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L260-L262)),
  the multinomial offset prose ([man/bart2.Rd:229](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L229)) and the offset.category.test \item ([man/bart2.Rd:272-273](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L272-L273));
  generics edits the own-class generic paragraphs ([man/bart2.Rd:308-320](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L308-L320)) and the summary/vars
  items ([man/bart2.Rd:279](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L279)). Different regions of one file; git will usually merge, but the file is
  long-line Rd and a hand merge is likely.
- R/generics.R - owned by GENERICS ALONE, provided one item moves: predict.bartMultinomial
  ([R/generics.R:1013-1070](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/R/generics.R#L1013-L1070)) is touched by J5 (rename offset.category.test -> offset) AND by J4 (add
  ci.level, refuse by name). MOVE the rename into the generics slice, together with the
  five predict(offset.category.test=) pins in inst/tinytest/test-multinomial-generics.R.
  r-surface then touches R/generics.R not at all (verified: dataSlotOrNULL has no
  generics.R use; the @offset.category read at [R/generics.R:1034](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/R/generics.R#L1034) is already a bare slot read and is
  inside predict.bartMultinomial anyway).
- [R/bart.R:867-874](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/R/bart.R#L867-L874) - the multinomial offset.test refusal names
  dbarts:::bartcoreSetCategoryTestOffset, which the handle-api slice moves out of the
  namespace; [man/bart2.Rd:229](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/man/bart2.Rd#L229) repeats the same claim. RESOLUTION: r-surface rewrites
  both to name the R5 method $setCategoryTestOffset (which stays in the package),
  removing the coupling; handle-api then touches neither file. This is the one item
  that MUST move between slices to keep them disjoint.
- inst/tinytest/test-multinomial-r5-surface.R - r-surface (dbartsData spelling at
  [inst/tinytest/test-multinomial-r5-surface.R:132](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L132)/[inst/tinytest/test-multinomial-r5-surface.R:136](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L136)/[inst/tinytest/test-multinomial-r5-surface.R:140](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L140)/[inst/tinytest/test-multinomial-r5-surface.R:144](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L144)/[inst/tinytest/test-multinomial-r5-surface.R:434](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L434)/[inst/tinytest/test-multinomial-r5-surface.R:449](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L449), dataSlotOrNULL block [inst/tinytest/test-multinomial-r5-surface.R:469-479](https://github.com/vdorie/dbarts/blob/e35c8797ec5426de29c12dbc0761848b2d052453/inst/tinytest/test-multinomial-r5-surface.R#L469-L479)) and handle-api (it is
  one of the 34 bartcore-calling files). Unavoidable; handle-api lands last and re-runs
  its sed over whatever the earlier slices left.
- inst/tinytest/test-xbart-*.R - r-surface (M10 pins) vs the in-flight xbart-oracle
  worktree, which has test-xbart-reproducibility.R modified and test-xbart-fold-oracle.R
  new. Disjoint files today; r-surface should put its xbart pins in
  test-xbart-error.R / test-error-quality.R and stay out of -reproducibility.R.
- benchmarks/R/sbc.R - handle-api (18 getFromNamespace lines at [benchmarks/R/sbc.R:1055-1059](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/benchmarks/R/sbc.R#L1055-L1059), [benchmarks/R/sbc.R:1261-1265](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/benchmarks/R/sbc.R#L1261-L1265),
  [benchmarks/R/sbc.R:1596](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/benchmarks/R/sbc.R#L1596)) vs the SBC/leaf-scale arm, which landed 257 lines there at 0045507c and may
  come back. Land handle-api after any further SBC work, or coordinate that file.
- src/ - handle-api touches tree.hpp and sampler.hpp; the engine-default worktree has
  R_interface_bartcore.cpp, chain.hpp, model.hpp modified. Disjoint. E19's R-side fix
  (r-surface, R/model.R defaultNodeScale) is the twin of the engine-default worktree's
  [src/R_interface_bartcore.cpp:2298](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L2298) default arm - no shared file, but they should land in
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
D4. (g.5) C8, [R/plotTree.R:9-35](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/plotTree.R#L9-L35)'s dead padding branch, was folded into VD-H by the
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
E2. 8 dots pins in inst/tinytest/test-argument-surface.R ([inst/tinytest/test-argument-surface.R:387](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L387), [inst/tinytest/test-argument-surface.R:389](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L389), [inst/tinytest/test-argument-surface.R:392](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L392), [inst/tinytest/test-argument-surface.R:399](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L399), [inst/tinytest/test-argument-surface.R:403](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L403),
    [inst/tinytest/test-argument-surface.R:549](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L549), [inst/tinytest/test-argument-surface.R:553](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L553), [inst/tinytest/test-argument-surface.R:557](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L557)) + 1 in test-heteroscedastic.R ([inst/tinytest/test-heteroscedastic.R:260](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-heteroscedastic.R#L260)) = 9 expect_error sites on the
    dots channel: `grep -rn "unknown argument" inst/tinytest/*.R` (7 message lines, two
    of them comments) cross-read against the surrounding expect_error blocks.
E3. dots formals: `grep -n "^  \.\.\.$" R/bart.R R/rbart.R` -> [R/bart.R:704](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L704) (bart2),
    [R/rbart.R:53](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/rbart.R#L53) (rbart_vi); the other two hits ([R/bart.R:2487](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2487), [R/bart.R:2539](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2539)) are method
    signatures, not entries. Helper: `grep -rn "rejectUnknownDotsArgs|retiredDotsNames"
    R/ man/ inst/tinytest/ docs/` -> [R/utility.R:120](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/utility.R#L120),[R/utility.R:122](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/utility.R#L122),[R/utility.R:129](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/utility.R#L129),[R/utility.R:130](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/utility.R#L130), call sites
    [R/bart.R:765](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L765) and [R/rbart.R:62](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/rbart.R#L62), zero man/ hits, one test comment.
E4. bart()'s family vocabulary c("auto","logistic","aft") at [R/bart.R:2645](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2645); the
    four-token own-class list at [R/bart.R:2589-2594](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2589-L2594); the ten tokens [man/bart.Rd:174](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart.Rd#L174) names,
    read in full. Existing pins at [inst/tinytest/test-bart-bart2.R:108-122](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-bart-bart2.R#L108-L122) (4 tokens).
E5. offset spellings: `grep -rn "offset.category" R/*.R man/*.Rd inst/tinytest/*.R
    benchmarks/R/*.R` -> 17 hits in R/data.R (formals at [R/data.R:899-900](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/data.R#L899-L900), missing() clauses at
    [R/data.R:907-909](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/data.R#L907-L909), validation at [R/data.R:1430-1443](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/data.R#L1430-L1443)), R/A_class.R slots at [R/A_class.R:536-538](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/A_class.R#L536-L538),
    [R/generics.R:1017](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1017)/[R/generics.R:1033-1056](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1033-L1056) (predict.bartMultinomial), [man/bart2.Rd:80](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart2.Rd#L80)/[man/bart2.Rd:229](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart2.Rd#L229)/
    [man/bart2.Rd:272-273](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart2.Rd#L272-L273)/[man/bart2.Rd:312](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart2.Rd#L312), [man/dbartsData.Rd:14](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbartsData.Rd#L14)/[man/dbartsData.Rd:26-27](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbartsData.Rd#L26-L27), [man/dbarts.Rd:103](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbarts.Rd#L103),
    [man/dbartsSampler-class.Rd:142](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbartsSampler-class.Rd#L142)/[man/dbartsSampler-class.Rd:150](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbartsSampler-class.Rd#L150); tests: test-multinomial-r5-surface.R 21 hits
    (6 of them dbartsData calls) and test-multinomial-generics.R 5 (all
    predict(offset.category.test=)); benchmarks/: ZERO
    (`grep -rc "offset.category" benchmarks/R/*.R` all 0).
E6. M1's two sites: [R/data.R:1178-1180](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/data.R#L1178-L1180) and [R/data.R:1237-1239](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/data.R#L1237-L1239), both `stop("'x' must have the
    same number of observations as 'y'")` after codeResponse ([R/data.R:458-473](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/data.R#L458-L473)).
E7. The 31 handle wrappers: `awk 'NR>=1066 && NR<=1535' R/bartcore.R | grep -cE
    "^[a-zA-Z][A-Za-z0-9_.]* <- function"` -> 32 definitions, of which resolveForestIndex
    ([R/bartcore.R:1252](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bartcore.R#L1252)) is a helper -> 31 bartcore* wrappers. Per-function caller census run over
    R/*.R, inst/tinytest/*.R and benchmarks/: bartcoreRun 5 in-package call sites
    ([R/bart.R:1488](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1488),[R/bart.R:1574](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1574),[R/bart.R:1826](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1826),[R/bart.R:2076](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2076); [R/xbart.R:695](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/xbart.R#L695)), bartcorePredict 2
    ([R/generics.R:1221](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1221),[R/generics.R:1354](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1354)), bartcoreSetModel 1 ([R/xbart.R:691](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/xbart.R#L691)); the other 28 have
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
    multinomial-equivalence 11 ([docs/plans/review-2026-08-24/gate-ledger.md:111](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/review-2026-08-24/gate-ledger.md#L111)).
    bcf-equivalence.R holds 50 bartcore refs, multinomial-equivalence.R 48,
    equivalence.R 0 (`grep -ohE "dbarts:::bartcore[A-Za-z0-9_]+" <file> | sort | uniq -c`).
E10. Dead engine members: `grep -rn "rightChildOf|setCurrentSampleNum" src/ tests/ inst/
    docs/` -> [src/bartcore/tree.hpp:366](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/tree.hpp#L366), [src/bartcore/sampler.hpp:485](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/sampler.hpp#L485), plus one
    docs/plans note. Zero call sites.
E11. reapplyForestWeights: 3 call sites ([R/dbarts.R:1066](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1066), [R/dbarts.R:1827](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1827), [R/dbarts.R:1859](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1859)), definition at
    [R/dbarts.R:1790](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L1790), pin at [inst/tinytest/test-forest-weights-r5.R:103-140](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-forest-weights-r5.R#L103-L140) with the oracle comment
    at [inst/tinytest/test-forest-weights-r5.R:108](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-forest-weights-r5.R#L108); adoptPointer: definition [R/dbarts.R:945](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/dbarts.R#L945), 2 call sites ([R/bart.R:1807](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L1807),
    [R/bart.R:2056](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L2056)); both have 0 man/ hits (`grep -rn "adoptPointer|reapplyForestWeights" man/`
    is empty). [inst/tinytest/test-host-shell-pins.R:16-21](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-host-shell-pins.R#L16-L21) lists both as infrastructure; its counts are
    at [inst/tinytest/test-host-shell-pins.R:43-45](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-host-shell-pins.R#L43-L45) (46 / 5 / 41).
E12. core-generalization.md's claim is at [docs/design/core-generalization.md:69-71](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/core-generalization.md#L69-L71) plus the "Per obs" table row at [docs/design/core-generalization.md:76](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/design/core-generalization.md#L76);
    the falsifying loop is [src/bartcore/facade.hpp:694-703](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/bartcore/facade.hpp#L694-L703) (two virtuals per
    observation).
E13. n.samples = 0: three refusal sites [R/bart.R:806-809](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bart.R#L806-L809), [R/xbart.R:94-96](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/xbart.R#L94-L96),
    [R/rbart.R:104-106](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/rbart.R#L104-L106); two pins ([inst/tinytest/test-error-quality.R:46](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-error-quality.R#L46),
    [inst/tinytest/test-xbart-error.R:24](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-xbart-error.R#L24)); zero pins on rbart_vi's wording
    (`grep -rn "no posterior draws" inst/tinytest/ benchmarks/` is empty);
    [man/dbartsControl.Rd:32-33](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbartsControl.Rd#L32-L33) is the sentence to extend.
E14. monotone: parseMonotoneSign at [R/model.R:548-570](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/model.R#L548-L570) with the "inc"/"dec"/"0" arms at
    [R/model.R:553](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/model.R#L553),[R/model.R:557](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/model.R#L557),[R/model.R:559](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/model.R#L559); zero tests pin those spellings
    (`grep -rn '"inc"|"dec"' inst/tinytest/*.R` is empty); [man/dbarts.Rd:72](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/dbarts.Rd#L72) and
    [man/bart2.Rd:203](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/man/bart2.Rd#L203) are the Rd items.
E15. The 12 "should be one of" pins that must survive G5:
    `grep -rn "should be one of" inst/tinytest/*.R` -> [inst/tinytest/test-augmentation.R:376](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-augmentation.R#L376),
    [inst/tinytest/test-argument-surface.R:259](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L259),[inst/tinytest/test-argument-surface.R:260](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-argument-surface.R#L260), [inst/tinytest/test-bart-bart2.R:109](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-bart-bart2.R#L109) (a comment),
    [inst/tinytest/test-hazard.R:339](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-hazard.R#L339),[inst/tinytest/test-hazard.R:346](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-hazard.R#L346), [inst/tinytest/test-hurdle-surface.R:77](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-hurdle-surface.R#L77),[inst/tinytest/test-hurdle-surface.R:85](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-hurdle-surface.R#L85), [inst/tinytest/test-nbinom.R:232](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-nbinom.R#L232),[inst/tinytest/test-nbinom.R:236](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-nbinom.R#L236),
    [inst/tinytest/test-prior-predictive.R:127](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-prior-predictive.R#L127), [inst/tinytest/test-xbart-model.R:206](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/tinytest/test-xbart-model.R#L206) - none of them a `sample`
    argument, so the wording change is safe.
E16. NEWS append point: `grep -n "subsection" inst/NEWS.Rd` -> the 1.0-0 BUG FIXES
    subsection opens at [inst/NEWS.Rd:995](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L995) and its itemize closes at [inst/NEWS.Rd:1854](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L1854); both landed waves appended
    immediately above that line.
E17. In-flight worktrees: `git worktree list` -> engine-default and xbart-oracle, both
    at 0045507c; `git status --porcelain` in each -> engine-default has
    src/R_interface_bartcore.cpp, src/bartcore/chain.hpp, src/bartcore/model.hpp,
    tests/cpp/Makefile modified; xbart-oracle has inst/tinytest/test-xbart-
    reproducibility.R modified and test-xbart-fold-oracle.R new.
E18. Suite size: 167 files in inst/tinytest (`ls inst/tinytest/*.R | wc -l`); wave 2's
    recorded assertion count is 7040 (landing note [docs/plans/release-candidate-review.md:619](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/docs/plans/release-candidate-review.md#L619)).
