# Matrix review 1: fitting entries, constructors, family resolvers, and their Rd

Branch `bartcore`, tip `b102e17c`. Read-only. Built from `git archive HEAD` into a
private lib (`R CMD INSTALL --preclean`); every call below ran against that lib.

## 0. Reproduction and cells added

- `matrix.R` + `matrix-sampler-census.R` re-run against the fresh lib: the emitted
  `matrix-grid.csv` is BYTE-IDENTICAL to the checked-in one (877 cells, `diff` 0
  lines). The prior grid reproduces.
- New and re-runnable: `matrix-entries.R` here (`R_LIBS=<lib> Rscript
  matrix-entries.R <outDir>` -> `matrix-entries-grid.csv`). **655 cells added**,
  1532 across both grids, closing the prior grid's two stated coarsenings:
  - Part D: all 13 union family tokens x 6 entries (81 cells) - each entry probed
    with tokens OUTSIDE its own vocabulary, where the alias/redirect gaps live
    (accept/refuse/raw: dbarts 11/2/0, bart2 13/0/0, bart 3/4/6, rbart_vi 2/1/10,
    xbart 4/4/5, dbartsSpec 6/5/2).
  - Part E: 12 conditioning columns x EVERY accepted (entry, family) pair, not one
    representative (490 cells); plus sparse predictors (6), bart2's `forest()` term
    route x 13 families (13), 62 resolver/ladder/alias probes.
- Rd sweep: **314 testable claims** extracted with line numbers from the eight
  topic pages. 79 have a bearing grid cell: **71 AGREE, 8 CONTRADICT**; 235 are
  SILENT against it (\value component presence, numeric identities, bitwise and
  timing claims - beyond an accept/refuse grid's reach).

## 1. Findings

### BLOCKER

- **E1. `extract(fit, type="trees", sample=<numeric>)` silently returns a FILTERED
  tree table; `sample="train"/"test"` dies on an NA.** `R/generics.R:363-382`
  (`extract.bart`) rewrites the call onto `object$fit$getTrees` and strips `object`
  (:379) and `type` (:380) but NOT `sample`; `getTrees`'s formals are
  `(treeNums, chainNums, sampleNums, ...)`, so R partial-matches `sample` ->
  `sampleNums`. `R/generics.R:1851-1866` (`extract.rbart`) has the same rewrite.

      extract(f, type="trees")            -> data.frame 106x5
      extract(f, type="trees", sample=1)  -> data.frame  20x5   # SILENTLY filtered
      extract(f, type="trees", sample="train")
        -> "missing value where TRUE/FALSE needed" (warning: NAs introduced by coercion)

  CONFIRMS and SHARPENS F3: the prior report had the numeric case as "silently
  ignores the filter". It does not ignore it - it applies a filter nobody asked for,
  so the tree table is wrong output on a call `man/bart.Rd:258` describes ("extract
  will accept ... chainNums, sampleNums, and treeNums"). 33 cells, one cause. Fix:
  **agent-fix** - `sample` belongs in the strip list beside `object` and `type`.

### MAJOR

- **E2. `man/bart.Rd:174` claims all ten non-`bart` family tokens are refused BY
  NAME; six of the ten get R's bare `match.arg` message.** Rd: "The other ten
  `bart2` tokens are refused BY NAME rather than falling through to `match.arg`'s
  generic 'should be one of' message, which names neither the token nor `bart2`".
  `bart(x, y, family=F)` for F in `gaussian`, `probit`, `hazard`, `hazard.probit`,
  `hazard.logistic`, `twopart` -> `'arg' should be one of "auto", "logistic",
  "aft"`. Only `multinomial`, `ordinal`, `nbinom`, `hurdle.lognormal` are named
  (`bart() does not fit family = "F"; use bart2(x.train, y.train, family = "F")`) -
  the literal contents of `bartOwnClassFamilies` (`R/bart.R:2587-2592`).
  `docs/design/feature-matrix.md` [f1] states the OPPOSITE of the Rd ("The
  remaining six tokens ... fall to `match.arg`'s generic message if typed"), and
  the code follows the design doc. Fix: **agent-fix** - the design doc says what is
  right; bart.Rd:174's first sentence is the wrong one.
- **E3. `"twopart"` is the tenth token and is missing from BOTH the Rd's
  enumeration and `bartOwnClassFamilies`.** bart.Rd:174 lists nine tokens across
  its three reason-groups; `twopart` appears in none. `R/bart.R:2587-2592` is
  checked against the RAW token before `bart2()`'s alias fold runs, so
  `bart(family="hurdle.lognormal")` gets a named redirect to bart2() while
  `bart(family="twopart")` gets bare `match.arg` - for the identical request, the
  two spellings being documented as the same family (`man/dbarts.Rd:111`,
  `man/bart2.Rd:237`). CONFIRMS F9. Fix: **agent-fix** (add `"twopart"` to
  `bartOwnClassFamilies`; the alias equivalence is already documented).
- **E4. `dbartsData()` mis-reports EVERY multi-column response with a factually
  false row-count message; four surfaces inherit it.** `R/data.R:1179-1181` and
  `R/data.R:1238-1241` compare `NROW(formula)` against `NROW(codeResponse(data)$y)`,
  and `codeResponse` flattens an n x 2 response (a `Surv` included) to length 2n.
  With `NROW(x) == NROW(sv) == 40` and `dim(sv) == 40x2`, all of
  `dbartsData(x, sv)`, `dbartsData(x, cbind(yG, yG))`, `dbartsData(x, counts)`
  (positional), `rbart_vi(x, sv, family="aft")`, `xbart(x, sv)` and
  `dbartsSpec(dbartsData(x, sv), ...)` (aft, hazard, hazard.* alike) report
  `'x' must have the same number of observations as 'y'`.
  GENERALIZES F2, which read this as an rbart_vi/aft accounting bug: the defect is
  in the shared ingest and rbart_vi is one of four callers. `dbarts(x, sv,
  family="aft")` succeeds because it extracts log-time+status BEFORE calling
  `dbartsData()`. Fix: **agent-fix** for the message (it asserts something false);
  **VD-judgement** (G-B) for whether the ingest should name a survival/matrix
  response - `man/rbart.Rd:95` already says survival "enter[s] only through the
  formula interface", the shape of refusal the Rd implies.
- **E5. `man/bart.Rd:153` states the `combinechains = TRUE` result shape
  backwards.** Rd: "if TRUE, samples will be returned in arrays of dimensions equal
  to nchain x ndpost x number of observations". Behavior:
  `bart(nchain=2, combinechains=TRUE)$yhat.train` -> 12x40 (collapsed matrix);
  `FALSE` -> 2x6x40 (the array); bart2's `combineChains` is identical. Its own
  \value at `bart.Rd:277` and `man/bart2.Rd:164` state it correctly. Fix: **agent-fix**.
- **E6. `man/bart.Rd:165` transposes the `proposalprobs` defaults.** Rd lists
  "birth_death, change, and swap ... Defaults are 0.5, 0.1, 0.4, and 0.5
  respectively", i.e. change = 0.1, swap = 0.4. The live default on both `dbarts`
  and `bart2` is `c(birth_death=0.5, swap=0.1, change=0.4, birth=0.5)`, as
  `man/bart2.Rd:200` states. Fix: **agent-fix**.
- **E7. `man/dbartsSpec.Rd:40`'s "a family can never resolve two ways" is falsified,
  and `:48` contradicts the function's own signature.** CONFIRMS F1.
  `dbartsSpec(dbartsData(x, yMultinom), ctl, family="multinomial")` -> `family
  "multinomial" cannot fit a 3-level factor response; ...` while
  `dbarts(x, yMultinom, family="multinomial", ...)` is accepted, and
  `dbartsSpec(dbartsData(x, counts=cnt), ctl, family="multinomial")` is accepted.
  `:48` calls `"multinomial"` unavailable because it describes "more than one
  sampler" - false (one K-forest sampler), and the token IS in dbartsSpec's
  vocabulary. `dbartsSpec.Rd:25` documents the counts route correctly.
  Fix: **agent-fix** - :25 is right; :40 and :48 are wrong.
- **E8. The per-category offset has three spellings and `dbarts()`'s own matrix
  interface reaches none.** REFINES F5 (which read it as keyword-unreachable
  everywhere).

      bart2(x, yM, family="multinomial", offset=<n x K>)           -> accepted
      bart2(..., offset.category=<n x K>)                          -> unknown argument 'offset.category'
      dbarts(x, yM, family="multinomial", offset.category=<n x K>) -> unused argument (raw)
      dbarts(x, yM, family="multinomial", offset=<n x K>)          -> 'offset' must have the same length as 'y'
      dbartsData(x, counts=, offset.category=<n x K>)              -> accepted
      dbarts(<that dbartsData>, family="multinomial")              -> accepted

  `man/dbarts.Rd:103` says "the per-category shift is `offset.category`" without
  saying it is not a `dbarts()` argument; the `offset=<n x K>` attempt gets a length
  message rather than a pointer at the channel that serves the caller.
  Fix: **VD-judgement** (G-B).
- **E9. Entries disagree on whether family or response is validated first.** With a
  bad family AND a `Surv` response, `dbarts`/`bart2`/`bart`/`rbart_vi` all report
  `'arg' should be one of ...` (family first); `xbart(x, <Surv>, family="zzz")`
  reports `'x' must have the same number of observations as 'y'`.
  Fix: **agent-fix** (resolve family first, as the other four do).
- **E10. `n.threads = c(1, 2)`: four entries name the argument, `xbart` raises a raw
  R coercion error.** `dbarts`/`bart2`/`rbart_vi`/`bart` -> `'n.threads' must be of
  length 1`; `xbart` -> `'length = 2' in coercion to 'logical(1)'`. Fix: **agent-fix**.
- **E11. `n.samples = 0` gets three different answers across five entries.**
  `dbartsControl`/`dbarts` accept; `bart2`/`xbart` -> `'n.samples' must be a positive
  integer`; `rbart_vi` -> `no posterior draws will be taken after thinning`.
  `man/dbartsControl.Rd:33` (a per-`run()` return count) and `bart2.Rd:152` (a sweep
  budget) are both silent on zero. Fix: **VD-judgement** (G-E).
- **E12. `monotone` accepts three undocumented spellings and is case-folded.**
  `R/model.R:548-570` (`parseMonotoneSign`) switches on `tolower(value)` and accepts
  `"inc"`, `"dec"`, `"0"`; its own refusal at `R/model.R:559-562` enumerates only the
  documented set. `dbarts(monotone=list(x1=S))` is accepted for S in
  `"inc"`, `"dec"`, `"0"`, `"INC"`, `"Increasing"`, and refused for `"up"` with
  `'direction' must be one of '+'/'-', 'increasing'/'decreasing', or +1/-1`.
  `man/dbarts.Rd:72` and `man/bart2.Rd:203` document only `"+"`/`"increasing"`/1 and
  `"-"`/`"decreasing"`/-1. CONFIRMS the options-audit item. Fix: **VD-judgement** (G-F).
- **E13. `dbartsDrawLatents()` refuses its own formal default when written out.**
  `R/augmentation.R:65` declares `sigma = 1`; `:79-81` guards with
  `if (!missing(sigma))`, so `dbartsDrawLatents("probit", fit=f, y=y, sigma=1)` ->
  `'sigma' applies only to family "aft" and "student", not "probit"`, while omitting
  it returns a numeric len40. Fix: **agent-fix** (guard on the value, or default NULL).

### MINOR

- **E14. `makeind(x, all=TRUE)` is a live formal that does nothing.**
  `R/bart.R:2833-2836` binds `all` to `ignored` and calls
  `makeModelMatrixFromDataFrame(x, TRUE)` unconditionally; both values give
  byte-identical output. PARTIAL REFUTATION of the options-audit framing:
  `man/makeind.Rd:26` DOES document it, as "Not currently implemented" - inert but
  documented. Fix: **VD-judgement** (implement vs delete from \usage; narrow).
- **E15. F4's raw-`unused argument` wall extends to two more entries.** F4 named
  `bart()` and `xbart()`; `dbarts()` has no `...` either, nor does `dbartsSpec()`.
  New-grid raw counts: dbartsSpec 42, xbart 32, dbarts 22, bart 12. e.g.
  `dbarts(..., samplerOnly=TRUE)` -> `unused argument (samplerOnly = TRUE)` on all
  11 accepted families; `bart(..., offset=o)` -> `unused argument` (bart spells it
  `binaryOffset`). Contrast `bart2`/`rbart_vi`: `unknown argument 'X'`. Fix:
  **VD-judgement** (G-C).
- **E16. A refusal never echoes the token typed when an alias folded.**
  `dbarts(family="twopart")` -> `family "hurdle.lognormal" fits two component
  samplers ...`; likewise every bart2 `"twopart"` refusal. Fix: **VD-judgement** (G-D).
- **E17. A refusal names the RESOLVED family, not the requested one.**
  `dbarts(family="hazard", resid.dist=student())` -> `... family "probit" has its
  own fixed error scale`; `bart2(family="hurdle.lognormal", resid.dist=student())`
  -> the same, naming only the occupancy component. Separately,
  `dbarts(family="multinomial", variance=~x1)` gives the generic gaussian-only
  message where `bart2` gives a multinomial-specific one. Fix: **VD-judgement** (G-D).
- **E18. `bart(keepevery=-1)` -> `'n.thin' must be a positive integer`** names an
  argument `bart()` does not have. Fix: **VD-judgement** (G-D).
- **E19. Two sibling R helpers 42 lines apart disagree on unknown-family policy.**
  `R/model.R:400-414` `defaultNodeScale` has NO default arm, so
  `dbarts:::defaultNodeScale("hazard")` and `("student")` return **NULL** silently,
  while `R/model.R:442-453` `defaultAmplitudePriorScale` `stop()`s by name. The C++
  mirror at `R_interface_bartcore.cpp:2291-2306` claims in comment to mirror the R
  one but has `default: return 0.5` and no multinomial arm (its enum has six).
  Fix: **agent-fix** (give `defaultNodeScale` the same `stop()`).
- **E20. The four C++ `default:` arms on a family switch - what each absorbs, and
  none is R-reachable.** `chain.hpp:756` (K-forest ctor) -> `GaussianResponse`,
  absorbing aft/ordinal/nbinom and any 7th enumerator; `chain.hpp:5026`
  (`latentScaleAnchor`) -> `scaledResponseSd()`, deliberate for gaussian (comment
  says why), also absorbing aft/ordinal/nbinom; `R_interface_bartcore.cpp:2291`
  (`defaultNodeScale`) -> 0.5 (gaussian's), absorbing gaussian and aft;
  `R_interface_bartcore.cpp:2812` (`validateResponseSupport`) -> `default: break`,
  deliberate. All three K-forest arms sit behind `refusedAmplitudeFamilyReason`
  (`R_interface_bartcore.cpp:2266`), called from both creation routes (:2313, :3150)
  and the factory (`facade.hpp:809-826`). The grid confirms by execution:
  `dbarts(forests=list(forest(), forest(basis=~z)), family=)` is accepted for
  gaussian/probit/logistic and refused BY NAME for aft, ordinal, nbinom,
  multinomial. No R call reaches any of the four. Fix: **VD-judgement** (G-G).
- **E21. `resolveFamily` vs `augmentationFamily` DO have disjoint vocabularies, but
  the gap is unreachable from R.** `R_interface_bartcore.cpp:1582` takes
  `""`/gaussian/probit/logistic/ordinal/nbinom/aft; `:6151` takes
  probit/logistic/ordinal/aft/nbinom/student and refuses gaussian. CONFIRMS the
  memo. But `R/augmentation.R:7`'s `augFamilies` is byte-identical to
  `augmentationFamily`'s set and gates first: `dbartsDrawLatents(family="gaussian")`
  -> `'arg' should be one of "probit","logistic","ordinal","aft","nbinom","student"`.
  Only the flat `dbarts_drawLatents` (`src/C_interface.cpp:1090`) reaches the C++
  refusal. Fix: **defer** (each table is total over its own callers).
- **E22. `defaultAmplitudePriorScale` has no multinomial arm - confirmed, inert.**
  It errors by name, but both call sites (`R/bartcore.R:777`, `R/model.R:1146`) sit
  behind the K-forest gate, which refuses multinomial first. Fix: **defer**.

## 2. Error-without-reason cells, grouped by root cause

243 across both grids (66 prior + 177 new). Four causes, one fix each:

- **G1, bare `match.arg` on a token outside the entry's vocabulary** - 70 cells
  (16 prior + 54 new; rbart_vi 10, bart 6, xbart 5, dbartsSpec 2, 31 in the new
  resolver block). ONE fix: a shared by-name family redirect ahead of `match.arg`
  on each entry, the shape `refuseBartOwnClassFamily` already has. `bart.Rd:174`
  already asserts this fix as if it were done (E2/E3).
- **G2, R's raw `unused argument` from an entry with no `...`** - 137 cells
  (16 + 121; dbartsSpec 42, xbart 32, dbarts 22, bart 12, plus resolver probes).
  ONE fix: `...` plus the shared `rejectUnknownDotsArgs` that `bart2`/`rbart_vi`
  run, or an explicit decision that these four keep R's wall (E15, G-C).
- **G3, `sample` partial-matching into `getTrees(sampleNums)`** - 33 cells; fix is
  E1's strip-list line. **G4, a raw coercion error** - 1 cell,
  `xbart(n.threads=c(1,2))`; fix is E10. Residue: 2 cells are this pass's own probe
  artifacts (`monotone` takes a named vector, not a formula), 1 is the prior grid's
  fixture note.

## 3. Prior findings: confirmed / refuted

- **F1 CONFIRMED** (E7), with the counts-built route shown to work.
- **F2 CONFIRMED and GENERALIZED** (E4) - not rbart_vi-specific; four surfaces.
- **F3 CONFIRMED and SHARPENED** (E1) - a numeric `sample` applies a filter rather
  than ignoring one.
- **F4 CONFIRMED and EXTENDED** (E15) - four entries lack `...`, not two.
- **F5 PARTLY REFUTED** (E8) - `bart2` DOES reach the per-category shift as a
  keyword, spelled `offset = <n x K>`; `dbarts()`'s matrix interface reaches
  neither spelling.
- **F6 CONFIRMED** - `rbart_vi(subset=, group.by=<unsubsetted>)` -> `'group.by' not
  of length equal to that of data`, on gaussian and auto alike; no fix named.
- **F7 / F8 carried, not re-probed** (plot-gap existence is not an entry/family
  question; F8 is value-level, and this grid checks no numbers).
- **F9 CONFIRMED and SUBSUMED by E2** - `twopart` is one of six tokens, not one.
- **Options audit `makeind(all=)`: CONFIRMED inert, framing PARTLY REFUTED** -
  `man/makeind.Rd:26` documents it as "Not currently implemented." **Monotone
  `"inc"`/`"dec"`: CONFIRMED and EXTENDED** - `"0"` and case-folding too (E12).
- **Options audit `family="probit"` unspellable on bart()/rbart_vi(): CONFIRMED** -
  both give bare `match.arg`; feature-matrix [f1] and `man/rbart.Rd:95` make the
  narrow vocabularies intended, so the defect is the message, not the set (E2).
- **Options audit `defaultAmplitudePriorScale` no multinomial case: CONFIRMED,
  unreachable** (E22). **bart()'s redirect misses twopart and the hazard tokens:
  CONFIRMED** (E2, E3).
- **Memo "rbart_vi tests == 0L not <= 0L, so a NEGATIVE thinned count passes":
  REFUTED.** All five refuse identically (`bart` via `keepevery`) with `'n.thin'
  must be a positive integer`.
- **Memo "n.threads ... 3 failures on length-2 input": REFUTED as stated; one live
  offender** (`xbart`, E10).
- **Memo "resolveFamily accepts gaussian where augmentationFamily rejects it":
  CONFIRMED in C++, shown unreachable from R** (E21).

## 4. Rd claims that AGREE (spot list of the 71)

- `bart2.Rd:306` "gaussian, probit, and logistic are the only families a term can
  join" - `bart2(y ~ x1 + x2 + z:forest(x2), family=F)` accepted for
  auto/gaussian/probit/logistic, refused by name for the other nine; and "no
  `forests =` formal of its own" - `unknown argument 'forests'` on all 13.
- `rbart.Rd:95` survival "do not support `weights` or `subset` ... and enter only
  through the formula interface" - formula route accepted, weights and subset
  refused by name, bare `cbind(t, s)` + explicit `family="aft"` accepted.
- `dbartsSpec.Rd:25` counts data resolves to multinomial from `"auto"`, refused
  under any other family; `:31` `survival=` is aft-only AND aft requires it (both by
  name); `:43` the empty-string C-API dispatch set matches `resolveFamily`'s shape
  gates exactly; `:56` `family` never `"auto"`. `dbartsData.Rd:24`/`:27` counts
  needs >= 2 categories and offset.category is refused without counts - by name.
- `bart.Rd:114`/`dbarts.Rd:42` weight policy (probit refuses, logistic requires
  positive integer counts, unit counts accepted) holds on all four entries taking
  `weights`, across all 13 tokens. `bart2.Rd:188` `predict` requires `keepTrees` - by
  name on all seven families probed; `dbarts.Rd:103` flat `offset` refused for
  multinomial - by name on `dbarts` and `bart2` alike.
- Documented family-specific conditioning refusals all fire BY NAME across the full
  family cross: `test` on hazard/hurdle (6 cells), `subset` on survival/multinomial/
  hurdle (12), `variance` off gaussian (18), `resid.dist` off gaussian (24),
  `forests` off gaussian/probit/logistic (4).

## 5. VD-judgement groups

- **G-A. May every fitting entry SPELL every family its siblings spell?** 30 cells
  (bart 6 unnamed tokens, rbart_vi 10, xbart 9, dbartsSpec 5). Decides whether E2's
  fix is "write the redirects" or "rewrite bart.Rd".
- **G-B. May every fitting entry spell every `dbartsData` channel?**
  `offset.category` on dbarts/bart2/rbart_vi/xbart, `weights`/`offset`/`test`/
  `subset` on dbartsSpec, and whether the ingest names a survival response (E4, E8).
- **G-C. Does an entry with no `...` owe a named unknown-argument diagnostic?** 137
  raw `unused argument` cells on dbarts, bart, xbart, dbartsSpec (E15).
- **G-D. Does a refusal echo the token TYPED or the one it RESOLVED to?**
  twopart -> hurdle.lognormal, hazard -> probit, keepevery -> n.thin (E16-E18).
  **G-E. Is a zero-sized run legal, and on which entries?** (E11)
- **G-F. Do the undocumented short and case-folded `monotone` spellings stay (then
  document) or go (delete)?** (E12)  **G-G. Should the C++ family switches lose their
  `default:` arms, so a 7th enumerator is a compile error, not a silent gaussian?** (E20)
