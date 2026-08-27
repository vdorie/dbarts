# Consolidated report - second whole-branch review, bartcore @ b102e17c

Consolidator pass, 2026-08-24. Read-only against the checkout. Every behavioral claim was re-run by THIS
pass against its own private library (`git archive HEAD` -> `R CMD INSTALL --preclean -l <d>/lib`, dbarts
1.0.0) and, for C++ items, its own `tests/cpp` build in a staged copy. Baselines first: the full 167-file
tinytest suite gives 0 failing files (89 s), `tests/cpp` gives "all tests passed". Probes in `consol/`;
mutation diffs and logs under the scratch dir named in section (e).

## (a) Summary

- CONFIRMED: 8 BLOCKERS, 44 MAJORS, 35 MINORS. REFUTED: 5 leg claims (2 of them by this addendum's
  replants). QUALIFIED: 6. NARROWED by a whole-suite replant: 6. FIXED BEFORE THE REVIEW BEGAN: 1.
- BINS: AGENT-FIX 44 canonical entries (8 BLOCKERS + 36 MAJORS) plus 15 minor one-liners; 9 VD-JUDGEMENT
  groups, deciding the other 8 MAJORS, ~260 executed grid cells and 8 deletion/seam items; 11 DEFER. Only
  M11's three engine divergences move draws; they land last.
- Three BLOCKERS are live user-visible defects (two wrong answers and one SEGFAULT); five are REACH GAPS - a
  planted defect the whole gate set (167 tinytest files plus 1355 C++ check() sites over 257 cases) leaves
  green. This pass re-planted 63 mutations against a full suite - 32 in sections (b) and the reach-gap
  table, 31 more for legs C and D in (g) - so every verdict here is whole-suite, and six leg findings were
  narrowed or refuted because of it.
- The five most consequential confirmed defects, in plain language: (1) `extract(fit, type = "trees", sample
  = 1)` silently returns a tree table filtered to a sample nobody asked for, and `sample = "train"` dies on
  an NA - three legs found it. (2) `survivalProbabilities()`, the point of a discrete-time hazard fit, fails
  on EVERY such fit whose training predictors carry column names, which is the normal spelling. (3)
  `fitted()` and `residuals()` on an `rbart` fit SEGFAULT the R session when the random-effect array has no
  dimnames - one assignment on an ordinary fit reaches it. (4) `dbartsData()` tells the user "'x' must have
  the same number of observations as 'y'" for every two-column response; the statement is false and four
  exported surfaces inherit it. (5) The facade - the one dispatch layer between the shipped flat C API and
  the engine - has no conformance test: 5 of 7 forwarding defects planted in it pass the entire C++ suite.

## (b) BLOCKERS and MAJORS

### BLOCKERS

- **B1. `extract(fit, type = "trees", ...)` corrupts its own call.** Asking for the trees of one posterior
  sample returns a table filtered to a DIFFERENT sample with no warning; asking by the documented
  `"train"`/`"test"` token dies on an internal NA; three other formals the method declares are refused with
  R's raw `unused argument`. EVIDENCE `R/generics.R:363-382` rewrites `match.call()` onto
  `object$fit$getTrees` and strips only `object` (:379) and `type` (:380), so `sample` partial-matches
  `getTrees`'s `sampleNums` (`R/dbarts.R:1917`); `extract.rbart` (:1851-1866) repeats it. REPRO (3 trees, 5
  samples, 1 chain, keepTrees): plain -> data.frame 75x5; `sample=1L` -> 17x5 SILENTLY FILTERED;
  `sample="train"` and the positional `extract(f,"trees","train")` -> `missing value where TRUE/FALSE
  needed`; `combineChains=FALSE`/`forest=1L`/`contribution=TRUE` -> `unused argument (...)`. CONTRADICTS
  `man/bart.Rd:203-205` and `:219`. LEGS F3, E1, G2 (33 grid cells + 102 on the combineChains axis).
  CONFIRMED. AGENT-FIX - strip the method's own formals from `treesCall`, or refuse them by name under
  `type="trees"`; Sonnet, ~15 lines. No draws.
- **B2. `survivalProbabilities()` is dead on every hazard fit whose training predictors carry column
  names.** Hazard fits are matrix-interface only and `man/bart.Rd:75` says column names are honoured, so the
  named design is the normal spelling - and the broken one. EVIDENCE `R/bart.R:2421-2424` (training) and
  `:2432-2435` (matrix `newdata`) build `cbind(<named covariates>, rep(seq_len(K), each=n))`, leaving the
  period column UNNAMED. REPRO (n=60, `Surv`, keepTrees): named training matrix -> `column names of 'test'
  do not match those of 'x': 'period' present in 'x' but not in 'test' (whose columns are 'x1, x2, x3, ')`;
  same for `newdata=<named matrix>` and for `times=`. UNNAMED training matrix -> array 5x60x60;
  `newdata=<data.frame>` -> 5x60x8. BROADER THAN THE LEG: a data.frame TRAINING input also fails (bart2
  converts it to a named matrix first), so the rule is "any named training design". CONTRADICTS
  `survivalProbabilities.Rd:45-47,:87-96`, `bart.Rd:328`. WHY NO TEST: `inst/tinytest/test-hazard.R:12`
  builds `matrix(runif(n*p),n,p)` with no column names (verified). LEGS G1 only. CONFIRMED and BROADENED.
  AGENT-FIX - name the period column in both branches, plus a named-column fixture; Sonnet, ~4 lines. No
  draws.
- **B3-B7. Five reach gaps.** Each was re-planted by this pass, alone, in a staged copy; R-side mutations
  installed into a private lib and scored over all 167 tinytest files, C++-side rebuilt in a staged
  `tests/cpp` and run whole. Pristine controls: 0 failing files / "all tests passed". Verdicts in the
  reach-gap table below. All five BIN AGENT-FIX (write the missing assertion); none moves draws.
  - B3 `R/dbarts.R:1025` `$growFromRoot(n.sweeps)` can discard `n.sweeps` [mutation-A 1].
  - B4 `combiner.hpp:906-913` the BCF near-zero multiplier snap can go [mutation-A 2, and mutation-B B5, the
    same site seen from the C++ side]. B5 `combiner.hpp:1818-1825` the multinomial veto's active-row mask
    can go [mutation-A 3]. Both feed `formForestVetoWeights`, which is what stops a prior tree draw from
    seating a leaf holding only weightless rows.
  - B6 `chain.hpp:1518` the latent refresh's vector-leaf cache invalidation can go - R GREEN and C++
    SURVIVED (planted twice); every Linear/GP-leaf fixture in tests/cpp is gaussian so the branch is never
    entered, while `facade.hpp:816` has no family gate and no tinytest pairs a vector leaf with a latent
    family [mutation-B 1].
  - B7 `facade.hpp` 5 of 7 forwarding defects survive (`:429` `filledSavedDraws` -> `savedTreeCapacity`;
    `:455` `setResponse(y,!updateScale)`; `:515` `savedSlotForDraw` -> identity; `:621` `setForestWeights`
    on forest 0 regardless; `:520` `savedTree` on forest 0 regardless). Every other test drives `Sampler<L>`
    directly, so the 61-virtual facade is exercised only through `shape()` [mutation-B 2].

### MAJORS - live behavior (all re-run by this pass)

- **M1. `dbartsData()` mis-reports EVERY multi-column response; four surfaces inherit it.**
  `R/data.R:1179-1181,:1238-1241` compare `NROW(formula)` against `NROW(codeResponse(data)$y)` and
  `codeResponse` flattens an n x 2 response to length 2n. REPRO with `NROW(x)==NROW(sv)==60`,
  `dim(sv)==60x2`: `dbartsData(x, Surv(t,s))`, `dbartsData(x, cbind(y,y))`, `dbartsData(x, <counts>)`
  positionally, `rbart_vi(x, Surv, family="aft")` (at n.samples 40, so thinning is not the cause) and
  `xbart(x, Surv)` ALL give `'x' must have the same number of observations as 'y'`; the rbart_vi FORMULA
  route succeeds, and `dbarts(x, sv, family="aft")` succeeds because it extracts log-time+status BEFORE
  calling `dbartsData()`. LEGS F2 (read as rbart_vi/aft-specific), E4 (generalized). CONFIRMED and
  GENERALIZED. AGENT-FIX for the message (it asserts something false); the refusal SHAPE is VD-C. Sonnet, ~8
  lines.
- **M2. `bart()`'s by-name family redirect covers 4 of 10 tokens; `man/bart.Rd:174` asserts all ten and
  `docs/design/feature-matrix.md` [f1] says the opposite.** REPRO: named redirect for `multinomial`,
  `ordinal`, `nbinom`, `hurdle.lognormal`; bare `'arg' should be one of "auto", "logistic", "aft"` for
  `gaussian`, `probit`, `twopart`, `hazard`, `hazard.probit`, `hazard.logistic`. `bartOwnClassFamilies`
  (`R/bart.R:2587-2592`) is exactly those four and is checked against the RAW token before `bart2`'s alias
  fold, so `"hurdle.lognormal"` redirects and its documented alias `"twopart"` (`dbarts.Rd:111`,
  `bart2.Rd:237`) does not. LEGS F9, E2, E3, reading-R 3.3. CONFIRMED. AGENT-FIX for the `twopart` half (1
  line); the other five are VD-A.
- **M3. `dbartsSpec()` and `dbarts()` resolve `family="multinomial"` differently for the same input; two
  sentences of one Rd paragraph are wrong.** REPRO: `dbartsSpec(dbartsData(x, <3-level factor>), ctl,
  "multinomial")` -> `family "multinomial" cannot fit a 3-level factor response` while `dbarts(x, <same
  factor>, family="multinomial")` is ACCEPTED (`model@family=="multinomial"`); the counts-built route is
  accepted on both. Cause: `dbarts`/`bart2` auto-convert a factor to counts before building the data.
  `dbartsSpec.Rd:40` ("a family can never resolve two ways") is falsified; `:48` calls `"multinomial"`
  unavailable as describing "more than one sampler" - false, and the token IS in its own vocabulary
  (verified: auto, gaussian, probit, logistic, aft, multinomial, ordinal, nbinom). `:25` is right. LEGS F1,
  E7. CONFIRMED. AGENT-FIX (Rd).
- **M4. `man/bart.Rd:153` states the `combinechains=TRUE` shape backwards.** REPRO (2 chains x 6 draws x 60
  obs): TRUE -> `yhat.train` 12x60 (collapsed matrix), FALSE -> 2x6x60; identical on `bart()`. Its own
  `\value` at `:277` and `bart2.Rd:164` are correct. LEG E5. CONFIRMED. AGENT-FIX (Rd).
- **M5. `man/bart.Rd:165` transposes the `proposalprobs` defaults** ("birth_death, change, and swap ...
  Defaults are 0.5, 0.1, 0.4, and 0.5 respectively"). REPRO: both `eval(formals(dbarts)$proposal.probs)` and
  `eval(formals(bart2)$proposal.probs)` are `c(birth_death=0.5, swap=0.1, change=0.4, birth=0.5)`, which
  `bart2.Rd:200` writes correctly. LEG E6. CONFIRMED. AGENT-FIX (Rd).
- **M6. `$n.chains` is absent from every fit that keeps a sampler, which `bart.Rd:315-317` and
  `bart2.Rd:164` promise unconditionally.** `R/bart.R:390-394`: `if (keepSampler) result$fit <- fit else
  result$n.chains <- n.chains`. REPRO (2 chains): keepTrees x keepSampler in {TT,TF,FT} all give
  `is.null(fit$n.chains)`; only FF sets it. The Rd's stated purpose - "information that can be lost if
  combinechains is TRUE" - is exactly the missing case. bartMultinomial/ bartOrdinal/bartNegbin DO carry it
  under keepTrees; bart and bartHurdle do not. LEG G5. CONFIRMED. AGENT-FIX (set it unconditionally, as the
  siblings do; 2 lines).
- **M7. A saved `bartHurdle` fit cannot `predict` after reload and no Rd names the recipe.** REPRO:
  `names(fH)` is `call, family, y, occupancy, positive` - no `$fit` - so `bart.Rd:248`'s
  `bartFit$fit$storeState()` does not apply, and `:248` names only multinomial/ordinal/nbinom. After
  `saveRDS` and a FRESH process, `predict(fH, xTest)` -> `samplers cannot be re-created without a stored
  state`, while bart (12x8), bartMultinomial (12x8x3), bartOrdinal (12x8x3) and bartNegbin (12x8) all
  replay. LEG G9. CONFIRMED. AGENT-FIX (extend `bart.Rd:248`).
- **M8. `dbartsSampler-class.Rd:328` "Reading it forces the sampler's *current* state" is false, and the
  stale value is then rejected by the sampler's own `setState`.** REPRO (two-forest sampler, `run()`, then
  `setForestBasis(2L, <1-column basis>)`): `s$state[[1]]$glue` reads length 8 AFTER the mutation where a
  fresh `storeState()` gives 7; `setState(<state read before>)` -> `state is not consistent with this
  sampler`; with `storeState()` first it is ACCEPTED. The behavior follows the documented opt-in at `:116`;
  `:328` is the wrong sentence. LEG G10. CONFIRMED. AGENT-FIX (Rd).
- **M9. `dbartsDrawLatents()` refuses its own formal default when written out.** REPRO:
  `formals(dbartsDrawLatents)$sigma` is `1`; omitting `sigma` returns numeric length 60, while `sigma = 1`
  -> `'sigma' applies only to family "aft" and "student", not "probit"`. `R/augmentation.R:65` declares
  `sigma = 1`, `:79-81` guards on `!missing(sigma)`. LEG E13. CONFIRMED. AGENT-FIX (default NULL, or guard
  on the value; 2 lines).
- **M10. `xbart()` answers two cross-entry inputs with raw R errors where four siblings name the argument.**
  REPRO: bad family AND a `Surv` response -> `dbarts`/`bart2`/`bart`/`rbart_vi` all report `'arg' should be
  one of ...` (family first), `xbart` reports `'x' must have the same number of observations as 'y'`;
  `n.threads=c(1,2)` -> `'n.threads' must be of length 1` on all four, `'length = 2' in coercion to
  'logical(1)'` on `xbart`. LEGS E9, E10, reading-R 1.8. CONFIRMED. AGENT-FIX (resolve family first; route
  n.threads through the shared coercion).
- **M11. Three engine divergences where one rule has two implementations that answer differently.** Verified
  by reading here, exactly as the engine leg states. `sampler.hpp:891` installs a LIVE-sourced donor's
  variance trees with NO positivity check where the slot-sourced arm at `:919-923` applies one and says why,
  and `rebuildVarianceForest` (`chain.hpp:4356-4368`) then scatters a non-positive leaf into a divisor
  (reachable from a `.Call` with a hand-built state). `chain.hpp:2421` (`recoverTreeParameters`) and `:2453`
  (`applyNewData`) hard-code `forests_[0]` where siblings at `:2329,:2401,:2526,:2540` loop `forests_`, so a
  whole-data replacement on a BCF/multinomial chain leaves forests 1..K-1 on the old grid - guarded two
  layers away by `refuseMultiForestMutation`, no engine assert. `combiner.hpp:986` routes `drawGlue` on
  `shippedShape()` (`:1474`), which tests basis widths and canonicality but NOT the half-Cauchy flag; the
  shipped path refreshes only `prior[0].variance` (`:1171`) where the general path loops every forest with
  `halfCauchyScale > 0` (`:1213`), so two admissible specs get two different MODELS. LEGS reading-engine D2,
  D3, D5. CONFIRMED BY READING, not executed (each needs a hand-built state or a spec the R surface
  refuses). AGENT-FIX; ALL THREE MOVE DRAWS on the paths they reach - adjudication plus a baseline
  re-record, land LAST.
- **M12. Three shipped comments say a feature is "Not yet exposed through the R surface" when R reaches all
  three today.** VERIFIED: `chain.hpp:88` (monotone) vs `R/model.R:99`; `chain.hpp:141-142` (interaction
  constraints) vs `R/spec.R:111`; `chain.hpp:166` (variance forest) vs `R/spec.R:113`. LEG reading-engine C1
  (its top item). CONFIRMED. AGENT-FIX, 3 lines.
- **M13. Two CI gates are mis-wired so a green check means less than it looks, and five declared gates have
  never run.** VERIFIED: 11 workflow files on bartcore, 1 on main (`git ls-tree main`); `sanitizers.yaml`
  push branches are `[bartcore]` alone, so the per-push memory-safety gate AND the only tinytest-count floor
  stop at the merge to main; `rchk.yaml:57-59` fails only on a NON-EMPTY `grep -E '[PB]|[UP]' ... || true`,
  so an OOM-killed analyzer leaving no `.bcheck` is a GREEN job - the cure is already recorded at
  `docs/plans/release-candidate-review.md:1170`; `lint.yaml` and `pkgdown.yaml` declare a BARE
  `pull_request:`, so a PR into bartcore draws style checks and no correctness checks. Five workflows
  (equivalence, sbc, rchk, valgrind, revdep-smoke) have zero runs ever, and three benchmark harnesses cannot
  execute (`DBARTS_CHANGE_LOG`, `DBARTS_CHANGE_STATS`, `BC_FALSIFIER` all 0 hits in `src/`). LEGS
  gate-ledger 1/7/9/12, gate-ledger-read B/C/D/F, remediations 1/2/6/7. CONFIRMED. AGENT-FIX for the rchk
  hardening and the sanitizers branch list (~12 lines); landing the five on main is VD's alone.

### MAJORS - reach gaps (a planted defect the gates leave green)

Each mutation re-planted ALONE in a staged copy; controls pristine; diffs at `mutlog/<id>.diff`. All
AGENT-FIX ("write the missing assertion"); none moves draws. R = the 167 tinytest files.

- A4 `model.hpp:3415` ordinal `updateCutpoints` freezes the last free cutpoint (at K=3, every free cutpoint)
  -> R GREEN, C++ KILLED (`ordinal gamma_2 posterior variance (actual 0, expected 0.1334)`).
  `test-ordinal.R:356`'s only gamma_2 assertion is a mean-within-0.35 band the frozen cold start sits
  inside. FIX `sd(cutpoints[,2]) > 0.02`.
- A5 `grow.hpp:186` `growth <= 0.0` -> `< 0.0` kills the depth/availability veto -> R GREEN, C++ KILLED (`a
  growth-vetoed node draws nothing`). `test-grow-from-root.R` is 11 of 18 smoke.
- A6-A9 four R-side refusals/constants with no behavioral cover, all R GREEN: `R/data.R:1495` precision
  threshold `1e-10` -> `1e-30` (both warn-cases have `diff(range(y))` EXACTLY 0, so the constant the file
  exists to defend is never exercised); `R/data.R:627` the `family="ordinal"` minimum-category refusal
  deleted; `R/xbart.R:608` the mcr loss ignoring case weights (the rmse arm IS covered); `R/rbart.R:270` the
  "survival (aft) fits do not support 'weights'" refusal deleted.
- A10-A12 three more, all R GREEN: `R/dbarts.R:1929` `getTrees(current=TRUE)` ignored, reading the saved
  store, though `test-empty-leaf-veto-weights.R`'s entire structural oracle is `getTrees(current=TRUE)`;
  `R/generics.R:2364` `print.bart` always printing family "gaussian", where `test-plot-generics.R:103-104`
  are `expect_true(is.character(capture.output(print(x))))`, which CANNOT FAIL; and (MINOR) `R/data.R:735`
  predict-path basis finiteness deleted.
- B3 `chain.hpp:4207` the variance forest binds `y + meanFits` for `y - meanFits` -> C++ SURVIVED (the only
  value-level variance assertion uses data with NO mean function). MAJOR not BLOCKER: tinytest backstops it.
- B4 `sampler.hpp:575` predict writes at the CAPACITY stride, reads at the draw count -> C++ SURVIVED;
  regressed that is an out-of-bounds heap WRITE for any chain >= 1, and every saved-tree test is
  single-chain so the sanitizer cannot help either.
- B5 `moves.hpp:635` and `:659`, the descendant-interval and categorical-gauge rule-validity predicates ->
  C++ SURVIVED (both); neither has a direct caller in tests/cpp.
- B6 `grow.hpp:222` `-log(numPredictors)` for `-log(numAvailable)` - identical at a root, wrong below it ->
  C++ SURVIVED; all three grow-from-root law tests measure ROOT rules only.
- B7 `model.hpp:3314` the ordinal log-likelihood loses its lower tail -> C++ SURVIVED. This is the EXPORTED
  per-observation loglik channel of every ordinal fit.
- B8 `model.hpp:2185` the `+ log 2` a missing-bearing ordinal column owes -> C++ SURVIVED; it cancels on a
  same-variable redraw, which is why the grow-side twin IS caught.
- B9 `model.hpp:3603` `LogisticResponse::setWeights`' cold start of INACTIVE rows deleted -> C++ SURVIVED -
  the line that IS the logistic weight channel's landing claim.

## (c) VD-JUDGEMENT groups, ordered by items decided

Nine groups. VD-I is stated in full at (g.3) because its evidence is this addendum's; the other eight follow
here.

- **VD-B. Does a fitting entry with no `...` owe a named unknown-argument diagnostic?** DECIDES 137 raw
  `unused argument` cells (leg-counted: dbartsSpec 42, xbart 32, dbarts 22, bart 12, plus resolver probes) -
  E15; F4 named only two entries. Mechanism verified here: `dbarts(..., offset.category=oc)` -> `unused
  argument (...)` where `bart2(...)` gives `unknown argument 'offset.category'`. ALTERNATIVES (a) give the
  four a `...` plus the shared `rejectUnknownDotsArgs` `bart2`/`rbart_vi` already run - ~5 lines per entry,
  and because the helper REFUSES rather than drops, no typo goes silent; (b) keep R's wall and say so once
  per Rd (4 Rd lines), asymmetry permanent. RECOMMEND (a): the window locks the public face and the helper
  exists.
- **VD-E. Which fit classes carry which generics, and what does an absent one do?** DECIDES `plot` on
  bartOrdinal/bartNegbin/bartHurdle (verified: all three fall to `plot.default` -> `'x' is a list, but does
  not have components 'x' and 'y'`; bart, rbart, bartMultinomial have one) [G7, F6];
  `extract(type="loglik")` on the four own-class families (verified refused by name on all four; it is the
  loo/WAIC channel `bart.Rd:201` documents unscoped) [G8]; `extract(type="trees")` on a keepSampler-only fit
  (verified: an 11x4 frame of the CURRENT trees, no `sample` column, no warning, where `bart.Rd:257` scopes
  the feature to keepTrees while `plotTree.Rd:39-42` documents exactly that fallback for plotTree) [G6]; and
  49 Rd-silent cells (plotTree/survivalProbabilities/as_draws_* on the own-class fits, six pdbart/pd2bart
  gaps, and whether an `xbart` result gets a class at all - a bare array today, so `fitted`/`residuals`
  misfire on stats' defaults). RECOMMEND add `plot` and `loglik` for the own-class families (the ordinal
  likelihood is already implemented and exported - reach gap B7); make `extract(type="trees")` agree with
  `plotTree`; leave pdbart/pd2bart and `xbart` unclassed with one Rd sentence each.
- **VD-A. May every entry SPELL every family its siblings spell, and does a refusal echo the token TYPED or
  the one it RESOLVED to?** DECIDES 30 bare-`match.arg` cells, leg-counted (bart 6, rbart_vi 10, xbart 9,
  dbartsSpec 5) plus E16/E17/E18 (`twopart` refused as `hurdle.lognormal`; `hazard` refused as `probit`;
  `bart(keepevery=-1)` refused as `n.thin`, an argument `bart()` lacks). `man/bart.Rd:174` asserts the
  redirects as done; `feature-matrix.md` [f1] says the opposite and the code follows the design doc.
  ALTERNATIVES (a) write the redirects on all four entries, ~40 R lines, and bart.Rd:174 becomes true; (b)
  rewrite bart.Rd:174 to match the design doc, 3 Rd lines, narrow vocabularies kept. RECOMMEND (a) for
  `bart()` - its Rd promises it and the `twopart` half is a plain bug (M2) - and (b) for the other three;
  separately make every refusal echo the token TYPED, the only one a caller can act on.
- **VD-D. Do the own-class generics carry the bart-family argument vocabulary, or refuse it by name?**
  DECIDES 10 verified silent swallows on bartMultinomial/bartOrdinal/bartNegbin - `combineChains` and
  `forest` on extract, `sample` and `ci.level` on fitted, `ci.level` on predict, `type` on residuals, `vars`
  on summary - each `identical()` to the call without it, because those methods' formals are `(object, type,
  sample, ...)` against `extract.bart`'s `(object, type, sample, combineChains, forest, contribution, ...)`.
  `predict` on the SAME fits honours `combineChains` (2x6x8x3 against a combined 12x60x3). [G3, G4]
  RECOMMEND honour `combineChains` and `ci.level`, refuse the rest BY NAME; silence is the wrong option, and
  `bart2.Rd:38-49` is silent, so nothing constrains the call.
- **VD-C. May every entry spell every `dbartsData` channel, and does the ingest name a survival/matrix
  response?** DECIDES E8's three offset.category spellings - verified: `bart2(offset=<n x K>)` ACCEPTED;
  `bart2(offset.category=)` -> unknown argument; `dbarts(offset.category=)` -> raw unused argument;
  `dbarts(offset=<n x K>)` -> `'offset' must have the same length as 'y'`; only the pre-built
  `dbartsData(counts=, offset.category=)` route works through `dbarts` - and M1's refusal shape. RECOMMEND
  give `dbarts`/`bart2` a real `offset.category=` formal, and have the ingest name a two-column response by
  kind.
- **VD-H. Deletions and seams.** DECIDES `R/bartcore.R`'s low-level handle API (31 functions, leg-counted;
  spot-verified: `bartcoreSetCounts`, `bartcoreGetTrees` and `bartcoreUpdatePredictorPerObservationJointly`
  have 0 callers in `R/` and 3/6/2 in tests+benchmarks, `bartcoreBCFSampler`'s only `R/` hit outside its own
  file is a COMMENT at `R/model.R:434`, and `bartcorePredictPerForest` has none anywhere);
  `adoptPointer`/`reapplyForestWeights` (0 `man/` hits each); `dataSlotOrNULL` (`R/data.R:11-21`), whose own
  comment claims every internal read of `counts`/`offset.category`/`offset.category.test` goes through it -
  verified FALSE (5 bare `@counts` reads elsewhere in `R/`, a bare `@offset.category` at `R/generics.R:992`) -
  and which protects only objects serialized by an intermediate commit of THIS branch; the two dead engine
  members `Tree::rightChildOf` and `Sampler::setCurrentSampleNum` (one reference each, the definition); and
  S1, per-observation virtual dispatch at `facade.hpp:694-703` against `core-generalization.md:69-76`'s
  absolute "nothing dispatches per observation". RECOMMEND fold the handle API into the R5 methods behind
  one deliberate escape hatch; keep `adoptPointer`, document or delete `reapplyForestWeights`; delete
  `dataSlotOrNULL` and the two dead members; amend the DOC for S1 - the erasure is the joint sweep's whole
  purpose.
- **VD-G. Should the four C++ `default:` arms on a `ResponseFamily` switch go, so a 7th enumerator is a
  compile error rather than a silent gaussian?** DECIDES `chain.hpp:766` (K-forest ctor ->
  GaussianResponse), `chain.hpp:5033` (latentScaleAnchor -> scaledResponseSd),
  `R_interface_bartcore.cpp:2298` (defaultNodeScale -> 0.5), `:2842` (validateResponseSupport -> NO
  validation), plus the open-coded chain at `:6254` - all five verified present. The mutation leg's
  7th-enumerator probe found exactly three sites warn today; `:2842`'s own comment enumerates the harm it
  prevents (a negative nbinom count underflowing into a ~1.8e19 allocation, "an uncatchable crash").
  RECOMMEND delete all four - the cheapest gate on this list, and one family already travels under another's
  name (`chain.hpp:862`).
- **VD-F. Undocumented option vocabularies: keep-and-document, or delete?** DECIDES `monotone`'s `"inc"`,
  `"dec"`, `"0"` and its case fold (verified accepted; `"up"` refused, and the refusal at
  `R/model.R:559-562` enumerates only the documented set) [E12]; `makeind(x, all=)` (verified inert -
  `body(makeind)` is `ignored <- all; makeModelMatrixFromDataFrame(x, TRUE)` - but `man/makeind.Rd:26` DOES
  document it as "Not currently implemented", so the options audit's framing is partly refuted) [E14]; and
  `n.samples = 0`, three answers (verified: `dbartsControl`/`dbarts` accept, `bart2`/`xbart` -> `'n.samples'
  must be a positive integer`, `rbart_vi` -> `no posterior draws will be taken after thinning`) [E11].
  RECOMMEND delete `"inc"/"dec"/"0"`, keep and document the case fold; keep `makeind(all=)` as BayesTree
  signature compatibility; keep the zero-run split (a zero-sample sampler is the host-loop shape) but give
  the three entries ONE message and say so in `dbartsControl.Rd`.

## (d) MINORS (V = spot-verified by this pass)

- V `extract(type="trees")` omits the `chain` column on a single-chain fit where `bart.Rd:262` lists it
  unconditionally. [M5] V bart/rbart raise `sample must be in 'train','test'` where all four own-class
  `extract` methods give bare `'arg' should be one of ...`. [M1]
- V `residuals(fit, sample="train")` -> raw `formal argument "sample" matched by multiple actual arguments`
  [M2]; V `plotTree(fit, treeNum=1L, sample=2L)` silently partial-matches `sampleNum` and draws - the B1
  mechanism, benign today [M3].
- V `setForestBasis` on an amplitude-free sampler says `forest index out of range` where
  `dbartsSampler-class.Rd:189` promises the `setForestWeights` wording. [M4]
- V `dbartsSampler-class.Rd:308` offers "the return of `storeState`" to `setState`, but it returns NULL
  invisibly, so `s$setState(s$storeState())` -> `'state' must inherit from bartcoreState`. [M7]
- V `setForestBasis(k, ~var)` evaluates the formula in `environment(basis)`, never against the sampler's
  data, so `~x3` naming a predictor raises raw `object 'x3' not found`; the same bites `forest(basis = ~x3)`
  at fit time. [M8]
- V `names(fit)` lists NULL-valued `yhat.test`, `s.train`, `s.test`, so `nm %in% names(fit)` misleads. [M10]
  `bart.Rd:248` understates which generics stop on a stateless fit. [M6]
- V `defaultNodeScale` has no default arm, so `dbarts:::defaultNodeScale("hazard")` and `("student")` return
  NULL SILENTLY while the sibling `defaultAmplitudePriorScale` 42 lines away `stop()`s by name. AGENT-FIX.
  [E19]
- V `resolveFamily` and `augmentationFamily` really do have disjoint vocabularies, but `R/augmentation.R:7`
  gates first, so the gap is unreachable from R. DEFER. [E21, E22]
- V `rbart_vi(callback=)` is wired, documented, and SWITCHES the sampler onto the R Gibbs loop with zero
  uses in inst/tinytest, benchmarks/R or vignettes - a coverage gap, not a defect. V C2 the 40-line
  state-format comment opens "The shipped format (version 2)" 40 lines above `stateFormatVersion = 3`; V1
  the compatibility window is empty (both versions 3); V3 a refusal whose own comment says it is
  unreachable.
- Reading-cost classes, all DEFER, fully enumerated in their leg files: reading-R's 14 duplication clusters;
  102 of 129 process-narrating comment blocks citing a stripped `docs/` path (262 tree-wide) plus
  `R/diagnostics.R:4`'s non-existent `R/zzz.R`; and the support-lib residue D1, D7, U5, U6, U7 - only U7's
  doc line (`--with-xint-size`, at `kernel-vocabulary.md:26`, absent from `configure.ac:21`) is AGENT-FIX.
- mutation-A 13-17 / mutation-B 10-11: an unpinned error string, an equivalent mutant on a dead accessor, a
  2-assertion unfalsifiable test file, a logistic PG mutation that HANGS rather than fails (CI's only signal
  is a timeout), seven reach gaps whose sole killer is an untouched file, and 16 of 63 C++ catches that
  crash naming nothing (fix: `setvbuf(stdout, NULL, _IOLBF, 0)` in `tests/cpp/main.cpp`, one line).

## (e) Fix-wave plan

Slices by file area, in landing order. Artifacts: private libs, staged sources and mutation logs under
`/private/tmp/claude-501/.../scratch/consol-78583/` (`mutlog/<id>.diff`, `.tinytest.log`, `.run.log`);
re-runnable probes in `consol/`.

- S1 `R/generics.R` + `R/bart.R` - the two live BLOCKERS: B1 (strip the trees-call formals in
  `extract.bart`/`extract.rbart`) and B2 (name the period column in both `survivalProbabilities` branches,
  plus a named-column hazard fixture). Sonnet, ~20 code + ~25 test lines. No answer needed, no draws.
- S2 Rd- and comment-only: M3 (`dbartsSpec.Rd:40,:48`), M4 (`bart.Rd:153`), M5 (`bart.Rd:165`), M7
  (`bart.Rd:248`), M8 (`dbartsSampler-class.Rd:328`), minors M4-gen/M5-gen/M7-gen, U7's
  `kernel-vocabulary.md:26` line, and M12's three "Not yet exposed" comments. Sonnet. No draws.
- S3 `R/bart.R` + `R/data.R` + `R/augmentation.R` + `R/xbart.R` + `R/model.R`: M1's message, M2's `twopart`
  line, M6, M9, M10 (both halves), E19. Sonnet, ~25 lines. M1's refusal SHAPE and M2's other five tokens
  WAIT ON VD-C and VD-A.
- S4 `.github/workflows/` + `tools/` - M13: harden `rchk.yaml` per its own recorded lesson, add `[main,
  master]` to `sanitizers.yaml`, fix the bare `pull_request:` on lint/pkgdown, delete the three unrunnable
  harnesses, wire `check-win-drift.R` and `check-doc-freshness.R`. Sonnet, ~40 lines. Landing the five
  schedule-bound workflows on main is VD's alone.
- S5 tests, R side - B3, the R halves of B4/B5, reach gaps A4-A12: one assertion each, in the file the leg
  names. Sonnet under an Opus read, ~80 lines. No draws.
- S6 tests, C++ side - B6, B7, reach gaps B3-B9: a facade conformance suite (the largest single item),
  direct unit tests for the four never-called `*IsValid` predicates, a depth >= 1 conditional-law test, a
  multi-chain partial-fill predict test, an ordinal-likelihood pin, a non-constant-mean variance fixture,
  the one-line `setvbuf`. Opus, ~400 lines. No draws.
- S7 C2's 28 narrative lines and the `docs/`-citation policy. Sonnet, AFTER VD answers the citation
  question. No draws.
- S9 `R/generics.R` + `src/R_interface_rbart.cpp` - BLOCKER B8: refuse an NA group match before
  `.Call(C_rbart_fitted, ...)` and bounds-check the group index C-side, plus one tinytest that strips the
  dimnames and expects an error rather than a crash. Sonnet R-side, Opus for the C bound; ~10 lines. No
  maintainer answer needed, no draws.
- S5 GROWS with legs C and D: M14-M21's assertions, the reach gaps the (g) replant table proves whole-suite,
  and the (g.6) one-liners. Same slice, same shape - one assertion each in the file that already claims the
  property. Sonnet under an Opus read; ~200 lines total now. No draws.
- S6 GROWS by one line, M16: `test-sampler-state-emptyLeafVeto.R` either rescales its fixture until it
  separates or drops its regression claim and points at `tests/cpp/test_moves.cpp:1212-1218`, which this
  pass proved does own the invariant (the C++ half is already proven).
- S10 CALIBRATION ENABLERS, from (h)'s not-covered pair. (i) A creation-time leaf-scale pin plus
  rebuild-per-replication - the FIX-A shape `runsbcbcf-repair.md` priced at ~30-60 engine lines and did not
  take - unblocks an aft SBC arm AND retires the grouped-gaussian sigma artifact, which is the same
  re-anchoring; no new R surface. (ii) A heteroscedastic arm lifted R-side through `setState` (build valid
  variance-forest structures in the flat state format, draw leaves from the calibrated inverse-chi-square,
  install, read the variance channel back), ~100-150 harness lines, no engine change, roughly a day.
  RECOMMEND (i) PRE-RC - it is cheap, engine-only, closes the sharpest calibration gap and fixes a standing
  harness artifact; DEFER (ii) unless (i) lands early, since heteroscedastic is gated at m = 1 under a
  CONSTANT predictor and so has the weakest ensemble-scale evidence of anything shipped. The aft status
  setter (~100-150 lines over three layers) is NEW PUBLIC SURFACE: under the pre-release window it lands
  before 1.0-0 or not at all, and (i) makes it unnecessary.
- S4 GROWS by one item, from (h)'s exact-gate census: `logistic-reference` part 2 - the only true
  ensemble-scale comparison against an independent implementation on the push path - does not gate, because
  `anyFailure` is fixed before it runs, and is skipped unless BART is installed. Make it gate, or say in the
  script why it does not.
- VD-I's branch gates nothing else: under the recommended (b) it is a `bart.Rd` item plus a refusing
  coercion (Sonnet, ~6 lines); under (a) it changes `inst/include/dbarts/dbarts.h` and must precede any
  consumer-facing freeze.
- S8 (LAST, DRAW-MOVING) `src/bartcore/` - M11's three divergences. Each needs adjudication of the intended
  model, then a re-record of the affected baselines, then the equivalence trio. Opus. DEPENDS ON a
  maintainer ruling, D5 above all: two admissible specs currently get two different models and which is
  intended is a modelling call.
- Deletion slice, GATED ON VD-H.

## (f) Refuted or qualified leg claims, with the probe

- REFUTED (lenses memo) "rbart_vi tests == 0L not <= 0L, so a NEGATIVE thinned count passes". PROBE
  `bart(keepevery=-1)` and `n.thin=-1` on the others: all refuse with `'n.thin' must be a positive integer`.
  The entries leg reached the same verdict.
- REFUTED (lenses memo) "n.threads has 3 failures on length-2 input". PROBE `n.threads=c(1,2)` on all five
  entries: four name the argument; ONE offender (`xbart`, M10).
- REFUTED (lenses memo F5) "SamplerBase declares 36 pure virtuals and has one implementation", "ResidT
  serves one non-default instantiation", "three leaf concepts constrain no template". PROBE reading-engine
  sec 1: 61 virtual declarations, ONE source-level impl instantiated FIVE ways; `float` is reachable from R
  via `control@storage == "single"`; all nine leaf concepts are live. Do not run a deletion pass off the
  memo's numbers.
- QUALIFIED (matrix F3's last paragraph) "with `sample` omitted, `extract(type='trees')` silently ignores
  the filter". PROBE the plain call returns the full 75x5 table: nothing is forwarded and nothing is
  ignored. The defect is the SUPPLIED-argument one only, as E1/G2 state.
- QUALIFIED (G1) "the data-frame branch is correct". PROBE a data.frame TRAINING input fails identically;
  only a data.frame `newdata` works. Broader than reported, not narrower.
- QUALIFIED (matrix F5) "`offset.category` is unreachable as a keyword on every fitting entry". PROBE
  `bart2(family="multinomial", offset=<n x K>)` is ACCEPTED; it is `dbarts()`'s matrix interface that
  reaches neither spelling (E8 corrects this).
- QUALIFIED (memo, and both matrix legs) the (class x generic) empty-cell count: the memo says 18 of 63, the
  executed grid says 21 of 63 by `getS3method`. Use the executed number.
- NOT RE-PLANTED: mutation-B's equivalent-mutant claim R1 (`moves.hpp:280`, the dead birth-rejection
  restore) - this pass's pattern did not match, so the leg's own reasoning stands unchecked. It is a
  code-health item, not a defect.
- FIXED-PRE-REVIEW, not open (matrix F6-F8): the heteroscedastic silently-wrong loglik/ppd/summary the 08-24
  value scan reported was repaired at 221ec7af ("Score heteroscedastic fits at their own variance surface",
  2026-08-24 13:57, an ancestor of b102e17c; it added `inst/tinytest/test-heteroscedastic-channels.R`, 214
  lines). VERIFIED NUMERICALLY here against an independently coded oracle: on a `bart2(variance = ~x1)` fit,
  `extract(type = "loglik")` equals `dnorm(y, f(x), s(x), log = TRUE)` draw by draw, max abs diff 0, summed
  -751.3719; the old scalar-sigma rule would give -1798.922, the value scan's exact symptom shape. An
  earlier draft of this report listed it as open; it is not.

## (g) Addendum - mutation legs C and D

Legs C and D cover the 102 tinytest files leg A did not table: C files 1-51 (test-aft.R through
test-ppd-sigma-pairing.R, 132 mutations planted, 129 scored, 1105 expect_ sites), D files 52-102
(test-predict-sparse.R through test-zero-weights.R, 66 scored before a budget stop, 684 sites). D scored
against its 51-file half only and could not afford the full-suite replants that separate a MAJOR from a
BLOCKER; C scored per target file and escalated 6 survivors wide over its own 51. THIS PASS closed both
gaps: 31 mutations were re-planted against the FULL 167-file suite (pristine control 0 failing files / 6946
assertions), and D's key engine mutation was additionally planted against tests/cpp. Verdicts below are
therefore whole-suite, not half-suite. Driver and per-mutation logs: `addlog/` under the scratch dir named
in (e); the two legs' own catalogues (`mutation-C-evidence/mutations.jsonl`,
`mutation-D-evidence/mutations.jsonl`) recording not retained.

METHOD CAVEAT, recorded: the earlier header-mutation lane in section (b) was staged from a copy taken while
the R lane was mid-flight, so it carried leg A's `mA09` (the aft-weights refusal deleted) as a background
mutation. That lane's own pristine control was 0 failing files, `mA09` alone is a proven 0-failure mutation,
and it removes a refusal no test in the suite triggers - so no assertion's execution path changes and the
B4/B5/B6 R-side verdicts stand. Noted rather than hidden.

### (g.1) Full-suite replant verdicts

31 mutations re-planted here against all 167 tinytest files (pristine control 0 failing files / 6946
assertions): 16 SURVIVED and 15 were KILLED. A verdict of KILLED on a leg's zero-killer means the leg's
half-suite scope hid a real gate, not that the tree is untested; a verdict of SURVIVED promotes the leg's
finding to whole-suite.

- D's 14 zero-killer mutations, now scored over all 167 files:
  - `r03` rbart in-core: keepTrainingFits = FALSE still stores training fits -> **SURVIVED**
  - `s02` sliceSample: the upper boundary no longer clamps the stepping-out interval -> **SURVIVED**
  - `p02` sparse test column recode: category codes off by one -> **KILLED** by a run-aborting error,
    categorical test predictors must hold existing category codes
  - `p03` sparse test column recode: the reference level is left on the training coding -> **KILLED** by
    test-sparse-factor.R (3 assertions)
  - `p04` sparse test column: unseen-level refusal dropped -> **KILLED** by test-sparse-factor.R (13
    assertions)
  - `pp3` samplePriorPredictive: sigma is held at sigest rather than drawn from its prior -> **SURVIVED**
  - `o01` dbartsControl: n.threads silently forced to 1 -> **KILLED** by test-capi.R, test-control-errors.R,
    test-control-valuesAreUsed.R (5 assertions)
  - `ro2` rbart_vi: DART split probabilities are dropped from the result -> **SURVIVED**
  - `sp2` uniform split probabilities are no longer canonicalized to the empty vector -> **SURVIVED**
  - `e01` empty-leaf veto priced at a finite -1e7 instead of -HUGE_VAL -> **SURVIVED** (and KILLED by
    tests/cpp: "a vetoed current against an admissible proposal loses outright")
  - `rb2` student(): non-positive df is accepted -> **SURVIVED**
  - `sd1` setData: the replacement data is quantized onto a single cut per column -> **SURVIVED**
  - `jo2` joint update: a numeric column resolves to the first column's name -> **SURVIVED**
  - `rr1` rbart_vi in-core path: a supplied seed no longer fixes the draw stream -> **SURVIVED**
- D's unproven-GAP files, their aimed mutation run here for the first time:
  - `pr1` constant-leaf prior draw uses half the prior scale -> **KILLED** by test-bartcore.R,
    test-calibration-midchain.R, test-calibration-prior-draws.R, test-sampler-prior.R (15 assertions)
  - `sr1` getSumsOfSquaredResiduals reports forest 0's own total, not the combined location -> **KILLED** by
    test-sampler-residuals.R (2 assertions)
  - `po1` per-observation setPredictor installs every value without the empty-leaf check -> **KILLED** by a
    run-aborting error, bartcore updatePredictorPerObservation produced a tree with an empty l
  - `po2` per-observation setPredictor scans in index order, not randomized -> **KILLED** by
    test-composition-sequences.R, test-sampler-setPredictorPerObservation.R (2 assertions)
  - `sm1` NEON addVectorsInPlace prefix loop diverges from the scalar kernel -> **KILLED** by
    test-bartcore.R, test-calibration-creation.R, test-capi.R, test-composition-sequences.R and 8 more (28
    assertions)
  - `st1` bridge sum-to-one tolerance back at the old 1e-10 -> **KILLED** by test-sum-to-one-tolerance.R (2
    assertions)
  - `st2` dbartsModel validity: proposal-probability tolerance back at 1e-10 -> **KILLED** by a run-aborting
    error, invalid class "dbartsModel" object: rule proposal probabilities must s
  - `wl1` logistic weights: the positive-integer count check is dropped -> **KILLED** by test-bcf-family.R,
    test-capi.R, test-logistic-weight-swap.R (8 assertions)
  - `wl2` R-side logistic count-weight validation dropped -> **SURVIVED**
- C's 8 wide-escalated survivors, now scored over all 167 files:
  - `gm1` predict.bart passes a nonsense n.threads down -> **SURVIVED**
  - `gm2` sampler$predict default n.threads pinned to 1 -> **SURVIVED**
  - `pd4` plot.pdbart renders only the first requested predictor -> **SURVIVED**
  - `bg3` bart2 grows one sweep regardless of n.grow.sweeps -> **SURVIVED**
  - `fr3` 3+-level factor refusal dropped from the single-forest entries -> **SURVIVED**
  - `gs1` setResponse(updateScale=TRUE) no longer refused under a grouped decorator -> **KILLED** by
    test-bcf-r5-surface.R, test-forest-basis-r5.R (2 assertions)
  - `se1` sampleTreesFromPrior never reaches the engine -> **KILLED** by a run-aborting error, 'names'
    attribute [10] must be the same length as the vector [0]
  - `dcm2` x/y interface silently discards the user offset -> **KILLED** by a run-aborting error, cannot
    replicate NULL to a non-zero length

### (g.2) New BLOCKER

- **B8. `fitted()` and `residuals()` on an `rbart` fit SEGFAULT the R session when the random-effect array
  has no dimnames.** `R/generics.R:2040-2046` (`fitted.rbart`, `type = "ev"`) reads `ranefNames <-
  dimnames(object$ranef)[[length(dimnames(object$ranef))]]` and passes `match(object$group.by, ranefNames)`
  straight into `.Call(C_rbart_fitted, ...)`, which indexes on `NA_INTEGER`. REPRO (this pass, no source
  mutation at all - one assignment on an ordinary fit): `f <- rbart_vi(y ~ x1 + x2 + x3, data = d, group.by
  = g, ...)`; `dimnames(f$ranef) <- NULL`; `fitted(f)` -> `*** caught segfault *** address 0x..., cause
  'invalid permissions'`, traceback `1: fitted.rbart(f)`, R aborts, exit 139. `residuals(f)` crashes
  identically (it routes through `fitted.rbart`). Mechanism verified: after the strip, `match(f$group.by,
  NULL)` is 80 NAs. BLAST RADIUS measured, one generic per subprocess: `fitted` SEGFAULT, `residuals`
  SEGFAULT, `predict` OK, `summary` OK, `extract` raises a clean R error (`no 'dimnames' attribute for
  array`). STRONGER THAN THE LEG REPORTED: leg D reached it only by mutating `R/bart.R:40` and called it
  "not reachable from the public surface today". It is reachable from the public R surface with a single
  assignment, and the flat C API invites hand-built rbart-shaped objects besides. Against the repo's own
  "safe over fast in R" rule a crash is not an acceptable answer to bad input. LEG mutation-D "Real defects
  exposed incidentally". VERDICT CONFIRMED and UPGRADED (leg filed it as a robustness note; it is a crash on
  the exported surface). BIN AGENT-FIX - refuse an all-NA / any-NA group match R-side before the `.Call`,
  and bounds-check the group index in `src/R_interface_rbart.cpp`'s `rbart_fitted` besides. COST Sonnet
  R-side (~4 lines) + Opus for the C-side bound (~6 lines), plus one tinytest. Does NOT move draws.

### (g.3) New VD-JUDGEMENT item

- **VD-I. `predict`'s `n.threads`: wire it, or drop the formal and say so on `bart.Rd` too?** DECIDES
  `predict.bart`'s `n.threads` formal, the `dbartsSampler$predict` formal, and the fate of
  `inst/tinytest/test-generics-multithreaded.R`. EVIDENCE, all verified here rather than inferred:
  `n.threads` occurs exactly ONCE in the R5 `predict` method (`R/dbarts.R:1079`) - in its own formal list -
  and the body's only `.Call` is `.Call(C_dbarts_bartcore_predict, ptr, x.test, offset.test)` against
  `DEF_FUNC("dbarts_bartcore_predict", bartcore_predict, 3)` (`src/R_interface.cpp:224`). `predict.bart`
  coerces `as.integer(n.threads)[1L]` (`R/generics.R:294`) and hands it to that method, which discards it.
  BEHAVIORAL PROOF without any mutation: on a keepTrees fit, `n.threads = 2L, 8L, 0L, -99L, NA_integer_` and
  even the string `"zzz"` are ALL accepted and all return a result `identical()` to `n.threads = 1L`; an
  invalid thread count that reached the engine could not be silent. QUALIFICATION OF THE LEG'S FRAMING, and
  it changes the fix: leg C says "do not leave a documented argument that the engine cannot see", but
  `man/dbartsSampler-class.Rd:152-157` ALREADY discloses it - "Currently has no effect: `run` and `predict`
  both execute serially regardless of the value passed here ... reserved for a future per-call override."
  That Rd is honest and AGREES with the code. `man/bart.Rd` is the one that does not: `:49` puts `n.threads`
  in `predict.bart`'s `\usage` while its only `\arguments` item, the shared `nthread, n.threads` at
  `:149-151`, describes a live thread count ("Integer specifying how many threads to use") written for the
  fitting entries. ALTERNATIVES (a) wire it - a 4-argument `bartcore_predict`, or a `setNumThreads()` around
  the call - which changes the ONE shipped header and so touches `LinkingTo: dbarts` consumers, and then
  needs a real discriminator in the test file; (b) keep the reserved formal, copy
  `dbartsSampler-class.Rd:152-157`'s disclosure into `bart.Rd`'s item, and make the coercion refuse a
  non-positive or NA value by name instead of swallowing it; (c) drop the formal from both signatures and
  delete the test file. RECOMMEND (b): the sampler Rd already settled the semantics, the argument is a
  declared future extension point, and (a) spends the shipped-header budget on a per-call override nothing
  has asked for. Under EVERY branch `test-generics-multithreaded.R` is a defect: verified by diff, its two
  blocks are `test-generics-correctValues.R`'s first two, differing only in the header comment, two integer
  literal suffixes, and `n.threads = 2L` for `1L`. It claims to "test that predict gives same result when
  single or multi-threaded" and it never runs multi-threaded, because nothing can.

### (g.4) What the replants changed, leg finding by leg finding

- **CONFIRMED whole-suite, MAJOR stands (the planted defect survives all 167 files).**
  - **M14. `plot.pdbart` can render only the FIRST requested predictor.** `R/plot.R:219`'s `for (i in xind)`
    narrowed to `xind[1L]` -> SURVIVED. `test-pdbart.R:88-90` states the claim in a comment ("the plot
    method renders each requested predictor into a null device") and then asserts
    `expect_silent(plot(pdb1))`, which sees neither how many panels were drawn nor which (verified: the file
    has 12 `expect_`, 2 of them `expect_silent`, 0 content pins on the plot). Same shape as leg A's
    `print.bart` gap. [C 2] FIX count the panels, or add `expect_error(plot(pdb1, xind = 99L))` plus a
    per-panel `expect_silent(plot(pdb1, xind = 2L))`.
  - **M15. The 3+-level factor refusal in the single-forest entries can be deleted.** `R/data.R:536-548` ->
    SURVIVED. `test-factor-response.R:85-103` is the file that exists to hold it - five `expect_error`s
    naming "multinomial" - and every one still passes, because the SECOND `stop()` at `R/data.R:549` (the
    `family == "auto"` arm) fires instead and carries the same word. The guard for the EXPLICIT-family
    callers, `xbart` and `rbart_vi`, is the one left unheld. [C 4] FIX anchor the non-auto refusals on their
    own text.
  - **M16. `test-sampler-state-emptyLeafVeto.R` is not the regression its header says it is.** Line 10 reads
    "Regression for the -1e7 -> -HUGE_VAL fix."; reverting `moves.hpp:120-123`'s `-HUGE_VAL` to the pre-fix
    finite `-1.0e7` -> SURVIVED all 167 files, because its fixture (n = 2000, `resid.prior = fixed(1e-6)`,
    50 trees) no longer drives a valid branch's score past -1e7. The invariant IS owned - by
    `tests/cpp/test_moves.cpp:1212-1218` - and this pass proved it by planting the same mutation there:
    KILLED, `FAIL: a vetoed current against an admissible proposal loses outright`. MAJOR, not BLOCKER, and
    the defect is the R file's claim. Distinct from BLOCKERS B4/B5, which are the veto WEIGHT channels in
    `combiner.hpp`; this is the veto RANK pricing in `moves.hpp`. [D 1] FIX rescale the fixture until its
    four assertions (`:43`, `:48`, `:52`, `:56`) separate, or drop the claim and point at tests/cpp.
  - **M17. `samplePriorPredictive`'s sigma prior draw can be replaced by the point estimate.**
    `R/dbarts.R:836`, `sqrt(df * sigest^2 * rawScale / rchisq(n, df))` -> `rep_len(sigest, n)` -> SURVIVED.
    `test-prior-predictive.R:104-111`'s (d) block asserts only `var(ppd.pinned) >= var(pinned1)`, which
    holds for any positive noise scale - while carrying its own prior uncertainty is the whole point of a
    prior predictive. [D 4] FIX at the pinned seed, assert a nonzero per-draw sigma spread.
  - **M18. `test-sampler-setData.R`'s ONE assertion is self-cancelling.** Quantizing the replacement data
    onto a single cut per column (`R/bartcore.R:441`) -> SURVIVED: it compares `sd(rowMeans(samples1$train) -
    y)` against `sd(rowMeans(samples2$train) - y)` at tol 1e-2 and both fits degrade together (verified:
    the file has exactly 1 `expect_`). [D 8] FIX add an absolute band on the post-setData fit beside the
    paired comparison.
  - **M19. `rbart_vi(keepTrainingFits = FALSE)` still storing training fits, `rbart_vi`'s in-core seed no
    longer fixing the draws, the DART split-probability report, the uniform split-probability
    canonicalization, and the joint per-observation update's numeric-column name resolution all survive the
    whole suite** (`r03`, `rr1`, `ro2`, `sp2`, `jo2`). Five more one-assertion gaps of the same shape, each
    in the file that names the feature. [D 9, 10 and the zero-killer set] FIX one assertion each; `rr1` also
    leaves `test-rbart-reproducibility.R` (7 `expect_`, 6 value comparisons) the one genuinely unproven GAP
    file of leg D's ten.
  - **M20. The R-side logistic count-weight validation can be deleted with all 167 files green** (`wl2`,
    `R/spec.R`) - while the C-side twin (`wl1`) is caught three ways (`test-bcf-family.R`, `test-capi.R`,
    `test-logistic-weight-swap.R`). The R layer is redundant defense that no test holds; that is a decision
    (keep it and pin it, or drop it), not a bug. FIX pin the R message, or delete it.
  - **M21. `bart2`'s forwarding of `n.grow.sweeps` is unpinned** (`bg3`: `growFromRoot(1L, ...)` regardless)
    -> SURVIVED. This is BLOCKER B3 one layer up - leg A pinned the sampler METHOD, this is the user's
    NUMBER on the way in - so it folds into B3's fix rather than standing alone. [C 3]
- **NARROWED by the replant: the leg's file-level finding stands, the tree-level one does not.**
  - `gs1` (`setResponse(updateScale = TRUE)`'s R-layer refusal deleted) is KILLED whole-suite by
    `test-bcf-r5-surface.R` and `test-forest-basis-r5.R`. Leg C's claim that the R guard is "unheld" is
    REFUTED at tree scope; what survives is that `test-grouped-swap.R:64-74`, written to hold exactly it,
    does not - its three `expect_error`s pass on the C bridge's message. [C 5 narrowed to a misplaced-gate
    finding]
  - `se1` (`$sampleTreesFromPrior()` made a no-op) and `dcm2` (`validateXYOffset` returning NULL) are both
    KILLED whole-suite. Leg C said as much for `se1`; for `dcm2` the tree catches it downstream (`cannot
    replicate NULL to a non-zero length`). What stands in both cases is the FILE:
    `test-generics-sequentialExecution.R` is 2 assertions, one of them `expect_inherits`, and its `:42`
    comparison is common-mode blind because BOTH sides call the mutated function;
    `test-data-compatibility.R` is 8 `expect_` with 7 `expect_inherits` and `test-data-formula.R` 13 with 11
    (both verified), so neither can see whether the offset/weights/subset they ingest arrived. [C 6, C 7
    narrowed]
  - `p02`/`p03`/`p04` (the sparse test-column level recode) are all KILLED whole-suite, by
    `test-sparse-factor.R` and by a run-aborting engine refusal. Leg D's "NO live coverage" is REFUTED at
    tree scope; `test-predict-sparse.R`'s own reach is the finding. [D 5 narrowed]
  - `o01` (`dbartsControl`'s `n.threads` forced to 1) is KILLED by `test-capi.R`, `test-control-errors.R`
    and `test-control-valuesAreUsed.R`. `test-rbart-multithreaded.R` is still a 1-assertion
    `expect_inherits` file that cannot fail, but the capability is gated. [D 6 narrowed]
  - `x05`/`x03`/`x04` were never zero-killer: leg D's own run shows `test-xbart-oracle.R` catching each. The
    findings are that `test-xbart-method.R`'s section headed "test that k-fold subdivides data correctly
    when data do not divide evenly by k" (verified at `:114`) asserts only `expect_inherits(xval, "array")`,
    and that `test-xbart-loss.R`'s 8 assertions never read a loss value. MISPLACED GATES, not coverage
    holes. [D 2, D 3]
- **PROVEN NOT GAPS - six of leg D's ten "unproven GAP" files.** Their aimed mutations, run here for the
  first time, are all KILLED: `pr1` -> `test-bartcore.R`, `test-calibration-midchain.R`,
  `test-calibration-prior-draws.R` + 1; `sr1` -> `test-sampler-residuals.R`; `po1` -> an engine refusal
  aborts the run; `po2` -> `test-composition-sequences.R`, `test-sampler-setPredictorPerObservation.R`;
  `sm1` -> 12 files; `st1` -> `test-sum-to-one-tolerance.R`; `st2` -> S4 validity aborts the run; `wl1` -> 3
  files. Leg D's GAP column was a scope artifact, not a coverage verdict. Still unprobed:
  `test-reproducibility-binaryResponse.R` and `-continuousResponse-singleThreaded.R`, for which no aimed
  mutation was ever written.

### (g.5) MINORS added by legs C and D

- `R/plotTree.R:9-35` decodeCategoricalSplits' PADDING BRANCH IS DEAD CODE: skipping the decode and padding
  with "R" instead of "L" both leave the three categorical files green, including
  `test-data-categorical-declared.R:53`'s `nchar(directions) == 4L` on a factor whose 4th level is never
  observed - the C side already emits one character per DECLARED level, so `padding > 0` never holds on any
  path they reach. Folded into VD-H (delete the padding loop and the factorLevels lookup, or find the path
  that needs it and test that path). [C 8]
- `R/model.R:1554` `resolvePriorScale`'s `prior.sd * k` conversion is unreached: every sampler in
  `test-calibration-prior-draws.R` spells `scale =`, never `sd =`. AGENT-FIX. [C 9]
- `R/generics.R:125` pointwiseLogLikelihood's zero-weight NaN flag can report -Inf instead; the file's
  zero-weight arm asserts only that the value is not finite, which is exactly the distinction the comment
  says the channel exists to preserve. AGENT-FIX. [C 10]
- `R/bart.R:2709` the burn side of `n.burn %/% n.thin` is unpinned on both `bart()` and `bart2()`;
  `test-control-valuesAreUsed.R`'s keepevery arm pins only the SAMPLE side. AGENT-FIX. [C 11]
- `R/model.R:1614-1627` `student()`'s own `df <= 0.0` refusal is redundant with `dbartsStudentDist` S4
  validity (`rb2` SURVIVED the whole suite; the S4 twin `st2` is caught). Keep and pin, or delete - same
  shape as M20's R-side logistic check. DEFER. [D 7]
- Equivalent mutants, recorded so a later leg does not re-file them: `sliceSample.R:185`'s upper clamp
  (`s02` SURVIVED - the beta target is 0 outside [0, 1], so every out-of-range proposal is shrunk away),
  `stateFormatVersion` 3 -> 4 (writer and reader move together, which is what
  `minReadableStateFormatVersion` is for), and leg A's `formForestVetoWeights` snap seen a third time (leg
  C's `bz2` again moved nothing, and `test-bcf-zero-multiplier.R` - the one untouched file that names the
  snap - does not see it either, which re-confirms BLOCKER B4 from a third direction). DEFER. [C 13, D 11]
- Reach gaps whose sole killer is another file, so no new assertion is owed - only a note not to trim the
  killer: leg C's six (`dcm3`, `mp3`, `me2`/`ml1`/`gl1`, `fo1`) and leg D's `ro2`-class items. DEFER. [C 12,
  D 10]
- Two more structurally unfalsifiable files, the leg A `test-bart-formula.R` shape:
  `test-rbart-multithreaded.R` (1 assertion, `expect_inherits`) and `test-mutate-then-serialize.R`
  (verified: ZERO `expect_` calls of its own - all 111 of its assertions live in
  `inst/common/stateContinuation.R`, which no leg has mutated, and which is the obvious next mutation
  TARGET). Folded into VD-E. [D 6, C static notes]

## (h) Calibration lane - ensemble-scale SBC

`calibration-sbc.md` landed as this addendum closed: benchmarks/R/sbc.R UNCHANGED, re-parameterised only
through `nTrees`, run at n.trees = 200 (`bart()`'s shipped default; `dbartsControl` defaults to 75, so 200
brackets both) against the 50 every recorded SBC had used. 3 h 43 m wall, 9.45 h CPU at 3 concurrent.
Admission band Bonferroni'd over the matrix's 83 functionals, 0.05/83 = 6.02e-4.

- RESULT 11 families carry an ensemble-scale verdict; 75 of 83 functionals PASS at the 5% band and 3 FLAG at
  the Bonferroni band. Covered and clean at m = 200: gaussian, student-t, weighted gaussian, probit,
  logistic, ordinal, nbinom's identified mean, multinomial, BCF's identified effect and prognostic
  functions, and grouped random intercepts under both bases.
- THE THREE FLAGGED FUNCTIONALS COME FROM TWO PRE-ADJUDICATED CAUSES. (1) nbinom `r`/`agg.psi` (ecdf 0.2000
  / 0.1785) reproduces the recorded flag at four times the forest size and is the demonstrated r-vs-psi
  mixing ridge - integrated autocorrelation 2199-6359 sweeps, two chains holding disjoint r basins for
  100000 sweeps with identical `avg.mu`, the n = 20 control passing everything. Forest-size-invariant, as a
  ridge the trees never touch would be. (2) grouped-gaussian `sigma` (0.2916) is a HARNESS artifact
  confirmed three ways in this run - the arm's own tau/b1/b2 and all five f* pass, grouped-PROBIT (no
  response scaling) is 9/9, and plain gaussian's sigma is 0.0513. Not a sampler finding.
- ONE RECORDED OPEN ITEM CLOSES, ONE IMPROVES. multinomial's three raw per-forest `f_ik` cells - the finding
  `sbc-family-tiers.md` could not discharge at m = 50 (persistent U, chisq p 0.000, not shrinking under 3x
  spacing) - are CLEAN at m = 200 (0.0779 / 0.0701 / 0.0439 against a 0.0924 band, chisq p 0.326-0.958, no
  U). It does not clear the suspected centering approximation, but it bounds its reach: at the shipped
  default the raw levels calibrate. And BCF's `sigma`, which read as a calibration DEFECT at m = 50, is
  clean at m = 200 at both spacings.
- NOT COVERED, and this is the honest answer to "which families are trusted from a single-tree derivation
  only": aft and heteroscedastic. Both own real sampling code and neither reduces to a covered family.
  hazard and hurdle are calibration-INHERITED - each adds no engine code and each has a draw-for-draw
  bitwise reduction gate to a covered arm (hazard-reduction.R at 40 trees, hurdle-reduction.R at 30) -
  though hurdle's combine/retransform step still owes its analytic oracle.
- CORRECTION TO (b)'s M13, from the lane's own census: the memo's "14 of 16 exact gates are single-tree" is
  CORRECT but GENEROUS. The two exceptions exceed one tree only where the predictor is CONSTANT
  (`heteroscedastic-exact` partA is `x <- matrix(0.5, n, 1L)`; `multinomial-exact`'s multi-tree arms are
  intercept-only), so every tree is root-only and tree-structure MCMC plus the backfit residual bookkeeping
  are gated at m = 1 in 16 of 16, not 14 of 16. Also: `logistic-reference` part 2 (ntree = 50/200 against
  BART::lbart/pbart) is the ONLY true ensemble-scale comparison against an independent implementation on the
  push path - and it does not gate (`anyFailure` is fixed before it runs) and is skipped unless BART is
  installed. That is a concrete AGENT-FIX for slice S4.
