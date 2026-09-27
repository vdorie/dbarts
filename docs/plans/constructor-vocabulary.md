# constructor-vocabulary

agent: sonnet (R, tests, Rd, sister edits); opus (diff review)
rng: neutral. No draw changes: only name resolution moves. The same objects reach
  [`resolveSamplerSpec`](../../R/spec.R).
window: before the merge (dec-A67, TODO "constructor-vocabulary"); lockstep with stan4bart and bartCause.
budget: dbarts ~450 changed lines (R ~130 including the ~60-line resolver, Rd/pkgdown/NEWS ~90,
  tinytest ~210 including one new file of ~120, benchmarks ~10, docs ~10); stan4bart ~15;
  bartCause ~40.

## Goal

`interactions`, `blocks`, `forest` and `varianceForest` are no longer exported. Inside the arguments
that take them they resolve by bare name, whatever the caller has attached, under the rule in
[Decision: resolving a name the caller also binds (ruled)](#decision-resolving-a-name-the-caller-also-binds-ruled). Outside a call they are
reached through one exported list, `dbartsForests`. `?interactions` and its siblings still work, and
a `forest()` formula term fits as before.

## Context

Decision: dec-A67 in [docs/decisions.md](../decisions.md) ("Yeah, treat them like the priors"), and
dec-A56 (blocks stays beside interactions). The no-pollution rule is in
[3a. Prior specification](../design/public-surface.md#3a-prior-specification). None of the four
reached a release, so there is no tombstone.

How the vocabularies work today.
- [`vocabularyEnv`](../../R/family.R) builds a child of the caller's environment that binds each
  vocabulary name.
  [`evalInVocabulary`](../../R/family.R)`(expr, vocabulary, evalEnv, resolve)` evaluates the
  argument's unevaluated expression there and passes the value through `resolve`. For a `..N`
  forwarded through a wrapper's dots, it retries with the expression
  [`recoverForwardedArgument`](../../R/family.R) recovers.
  [`resolvedAs`](../../R/family.R) is the usual `resolve`: a bare constructor name means its
  defaults, and a wrong class is refused with "see ?topic".
- Prior sites: [`parsePriors`](../../R/model.R) (`tree.prior`, `node.prior`, `resid.prior`, over
  [`dbartsPriors`](../../R/model.R) plus `num.vars`/`numvars`). [`xbart`](../../R/xbart.R) reads
  `node.prior` and `tree.prior` over subsets of `dbartsPriors`, plus `k`.
  [`bart`](../../R/bart.R)'s multinomial branch reads `tree.prior`, and
  [`resolveConsolidatedArgs`](../../R/tombstones.R) handles the retired `resid.prior`.
- Family site: [`resolveFamily`](../../R/family.R) evaluates over
  `c(dbartsFamilies, dbartsPriors)`. It is called by [`dbarts`](../../R/dbarts.R),
  [`dbartsSpec`](../../R/spec.R), `bart` and `xbart`.
- Precedent test: ["viaDots <- function"](../../inst/tinytest/test-family-objects.R) checks the list
  shape, that no member is exported, the bare-name form, and forwarding through dots.

Where the four are taken today. Each one is forced as an ordinary promise, so today it resolves on
the search path.
- [`dbarts`](../../R/dbarts.R): `interactions`, `blocks`, `variance` (a shorthand or a
  `varianceForest()`) and `forests` (a list of `forest()`, each of which can nest `interactions()`
  and `blocks()`). All four are passed to [`resolveSamplerSpec`](../../R/spec.R). `forests` is
  forced first, by the formula-term collision check and by
  [`forestBasisDeclarations`](../../R/model.R).
- [`dbartsSpec`](../../R/spec.R): the same four, evaluated against `parentEnv` rather than
  `parent.frame()`.
- [`bart`](../../R/bart.R): `interactions`, `blocks` and `variance`, with no `forests` formal. These
  reach `dbarts` as the caller's own expressions:
  [`buildHostSamplerCall`](../../R/bart.R) uses [`redirectCall`](../../R/utility.R) and evaluates
  the call in `callingEnv`. The hurdle composition in [`bart2Hurdle`](../../R/bart.R) redirects the
  same way. The one place bart forces anything itself is the multinomial refusal
  `"'variance'" = !is.null(variance)`.
- Formula `forest()` terms. [`isForestCall`](../../R/formulaTerms.R) detects the call
  SYMBOLICALLY, as the name `forest` or `dbarts::forest`, so unexporting costs nothing there.
  However, [`processHit`](../../R/formulaTerms.R) evaluates the knob arguments with a plain
  `eval(..., envir = env)`, so a term such as `forest(x1, interactions = interactions(max.order = 1))`
  needs the vocabulary. Also, [`finalizeTermForests`](../../R/formulaTerms.R) calls
  `dbarts::forest()` and `do.call(dbarts::forest, args)`, and both fail once the name is unexported.
  That would break every formula-term fit.
- Internal only: [`bartcoreBCFSampler`](../../R/bartcore.R) takes `mu.interactions`,
  `tau.interactions`, `mu.blocks` and `tau.blocks` as values, and tests reach it through `dbarts:::`.
  It stays value-typed.
- No sites: `bart2` and `rbart_vi` are tombstones in [`rbart_vi`](../../R/tombstones.R). `xbart`,
  `bartBT` and [`dbartsControl`](../../R/dbarts.R) take none of the four.
- Consumers downstream check by S3 class (`dbartsInteractions`, `dbartsBlocks`, `dbartsForest`,
  `dbartsVarianceForest`) in [`resolveInteractions`](../../R/model.R),
  [`resolveBlocks`](../../R/model.R), [`resolveForests`](../../R/model.R) and
  [`resolveSamplerSpec`](../../R/spec.R). The class strings are unaffected. Function identity is
  checked only by `isForestCall` and in `finalizeTermForests`.

Sites outside an argument. I classified every call token by parse tree (script run over
inst/tinytest, inst/common, benchmarks/R, the Rd examples extracted by `tools::Rd2ex`, and R/).
"Outside" means the value is built somewhere other than inside a door argument: assigned first,
passed to a local helper's named formal, spliced by `do.call`, or called at top level.
- inst/tinytest, 71 sites in 9 files. Each needs `dbartsForests$...` or a helper reshaped to forward
  through `...`:

  | File | Sites | Where |
  | --- | --- | --- |
  | test-blocks | 22 | `doFitBlocks` 14, `do.call` lists 5 (one an `interactions()`), `bartcoreBCFSampler` 2, `expect_error` 1 |
  | test-interactions | 18 | `doFitInteractions` 14, `do.call` 1, `bartcoreBCFSampler` 2, `expect_error` 1 |
  | test-bcf-family | 10 | `transportParams` 6, assigned lists 4 |
  | test-predict-blend | 6 | `fitFromPredictBlend` 4, assigned 2 |
  | test-argument-surface | 5 | `varianceAttr` 3, `print` 1, `format` 1, all `dbarts::varianceForest` |
  | test-bcf-forest-channel | 4 | `packageFrom`, assigned |
  | test-bcf-creation | 2 | |
  | test-predict-forest | 2 | |
  | test-proposal-probs | 2 | `dbarts::forest` |

  Two further sites are forwarded through dots, ["fitOf <- function(...)"](../../inst/tinytest/test-fit-descriptors.R)
  and ["anchorSampler <- function(response, ...)"](../../inst/tinytest/test-calibration-prior-draws.R).
  Both should resolve unchanged through dots recovery and stay as live pins.
- inst/tinytest, 25 more sites across 12 files. These sit inside a door argument but are spelled
  `dbarts::varianceForest(`, `dbarts::forest(` or `dbarts::interactions(`, and `::` fails on an
  unexported name. The fix is mechanical: drop the prefix.
- benchmarks/R, 9 sites. `composition-matrix.R` has 5 (`switch` arms and `buildBase` argument
  lists). `equivalence.R` has 4 (`samplerArgs` lists spliced by `do.call`).
- Rd examples, 1 site: the last line of the varianceForest.Rd example prints a bare
  `varianceForest(...)`. vignettes: 0. The two mentions in gibbs_sampler_mixture_model.Rmd are
  inline prose. inst/common: 0. R/: the `finalizeTermForests` pair above.
- Total: 82 outside sites (71 + 9 + 1 + 1), plus 25 prefix drops. The remaining 271 door-argument
  sites and 52 formula-term sites need no change.

Sister packages (swept with `git grep` on each branch).
- stan4bart (bartcore), 3 sites plus one guard.
  - R/stan4bart.R binds each `dbartsPriors` name that appears in call position in `bart_args`
    (its `called_names`) to a quoting function, so the call reaches `dbartsSpec` unevaluated. That
    call is then evaluated with `parentEnv` set to `stan4bart_fit`'s frame, so a quoted
    `interactions(groups = g)` written inside a user function would miss `g` or pick up a local.
    Fix: in `defn_env`, bind called `interactions`/`blocks` directly to
    `dbarts::dbartsForests$interactions`/`$blocks`, with no quoting. The value is then built in
    the user's own frame and reaches `dbartsSpec` as an object. The priors have the same `g`
    hazard, which is out of scope and worth a stan4bart TODO line.
  - `bart_args$forests` today rides the generic forwarding loop (stan4bart TODO
    "bart-args-forests-guard"). Reserve `forests` beside `variance`, rather than adding `forest`
    to the bound set.
  - man/stan4bart.Rd (the `bart_args` item): today it says `interactions` and `blocks` are "ordinary
    values". Rewrite it: they resolve in call position, and `dbartsForests$...` builds one
    beforehand.
  - inst/tinytest/test-09-bartArgs.R: `dbarts::interactions(max.order = 0L)` inside
    `fitWith(list(...))` becomes `dbarts::dbartsForests$interactions(...)`. Add one case with a
    local `g` in a user function.
- bartCause (dbarts-1.0), 10 sites plus bcf's argument path.
  - R/bcf.R `fitBCF` builds `samplerEnv[["forests"]]` from 2 `dbarts::forest(` calls. Instead,
    build the `forests` EXPRESSION, `list(forest(vars = muVars, ..., interactions = <expr>), ...)`,
    into `samplerCall` and evaluate it in `samplerEnv`, whose chain reaches `callingEnv`. dbarts
    then applies the resolution rule below, including recovery of a nested `..N`. That reuses
    [`evalInVocabulary`](../../R/family.R) without a `:::` call and without a bare `eval`.
  - `bcf()` passes `matchedCall$mu.interactions`, `$tau.interactions` and `$tau.blocks` (NULL when
    absent) in `fitArgs`, in place of the forced values. `do.call(..., quote = TRUE)` already keeps
    language intact.
  - `bartc(method.rsp = "bcf", ...)` forces its dots through `list(...)`, and so does `bcf()` for
    its own `...` extras (`extraArgs`). man/bcf.Rd and man/bartc.Rd say that an argument passed
    through either `...` is evaluated before dbarts sees it, and so is spelled
    `dbarts::dbartsForests$...`.
  - Tests: in test-03-responseFit.R and test-14-bcf.R, 4 `dbarts::forest(` inside
    `dbarts::dbarts(forests = )` drop the prefix, and 2 `dbarts::blocks(` become
    `dbarts::dbartsForests$blocks(`. Add a case for `bcf(tau.interactions = interactions(max.order = 1))`
    and one called through a wrapper's dots.
  - man/bcf.Rd: the `\link[dbarts]{interactions}`/`{blocks}` links survive through the kept
    aliases. The argument text gives the spelling.
- treatSens (dbarts-1.0 worktree): 0 sites. bairrtt (main): 0 sites.
- Timing: once dbarts unexports, bartCause's `dbarts::forest(` breaks every `bcf` fit, and no
  spelling works on both sides of the change (`dbartsForests` does not exist before it). So the
  sister commits push alongside the dbarts unexport commit. A staggered push would need
  `export(dbartsForests)` and its Rd moved into the first commit group.

## Decision: resolving a name the caller also binds (ruled)

VD, 2026-09-26: "Yes, it should default to the user's likely intent. We do this in other call
redirection functions, either here or in the sister packages, or maybe in blme. It can be tricky."

Orchestrator adjudication, 2026-09-26. This was agent-made after an independent critique probed the
first draft of the rule; it is not a VD ruling.
- Dropped the stack-walk recovery of a failed wrapper formal. Through `lapply`/`Map` it silently
  picked the wrong environment, since `substitute` reached the original expression while the
  caller frame was `lapply`'s. A defaulted formal's default evaluates in the function's own frame.
  The prior vocabulary offers no such recovery either.
- Rule 2 looks up only names in value position, and stops at `topenv(P)`.
- "Missing" means an empty default in `formals()`.
- A failed promise is re-forced, not cached, so the restart warning is muffled.

Prior art.
- stan4bart's `called_names` binding in R/stan4bart.R binds a vocabulary name only where it is
  CALLED. Its comment: "In call position a prior name has no competing reading".
- blme's R/priorEval.R `evaluateCovPriors` ("check to see if it refers to a variable in the calling
  environment") and `evaluateFixefPrior` look a bare symbol up in the caller first.
- dbarts' own [`recoverForwardedArgument`](../../R/family.R) recovers a `..N`.

The rule is implemented as a new `evalInForestVocabulary`, used only at the four constructors'
sites. The priors and families keep [`evalInVocabulary`](../../R/family.R) unchanged; there a bare
`cgm` still silently means `cgm()` even when the caller binds `cgm`, which is out of scope.

Rule, for an argument expression E, site vocabulary V and evaluation environment P. P is the
door's caller frame, `parentEnv` for `dbartsSpec`, or the formula's environment for a term's knobs.
1. Call position is always the constructor. Every call in E whose head is a bare symbol in V gets
   the constructor object inlined as its head. The walk skips `~`, `quote` and `bquote`. This also
   gives the mask protection.
2. Value position is the caller's when the caller binds it. Collect the names of V that appear as
   bare symbols in value position in E, skipping `~`, `quote` and `bquote`; no other name is
   looked at, so an unrelated caller formal is never forced. For each collected name n, walk P's
   enclosing chain (`parent.env`) up to and including `topenv(P)`:
   - A package wrapper therefore never sees a user's global, and the search path never counts.
   - The dbarts namespace and the base namespace never count.
   - The first binding found in a function frame is tested for a formal with an EMPTY default
     (`formals()` there, not `missing()`, which is TRUE for a defaulted formal too). When that
     formal was not supplied, it means NULL, the door default.
   - Any other binding is forced, and a defaulted formal evaluates normally in its own frame.
   - A function value is not the caller's value, since no site takes a function, so n stays the
     constructor.
   - Anything else, NULL included, is the caller's value.

   Names the caller does not claim are bound to the constructor in a child environment of P. E,
   after step 1, is evaluated there. This holds at any depth, so
   `forests = list(forest(interactions = interactions))` sees the caller's `interactions`.
3. A top-level value identical to the site's constructor (a bare `variance = varianceForest`) is
   called for its defaults, as [`resolvedAs`](../../R/family.R) does. A bare `interactions` or
   `blocks` then fails with the constructor's own "needs"/"requires" message.
4. Dots. A `..N` anywhere in E, top level or nested, is resolved before E is evaluated, keeping
   [`evalInVocabulary`](../../R/family.R)'s order: evaluated as it stands and, on failure,
   recovered and resolved under this rule in its recovered environment. The value is inlined into
   E. A forwarded argument that evaluates without error keeps its value, so a forwarded argument
   gets no mask protection. A probe shows `..N` recovery is correct through `lapply`.
5. No other recovery. A constructor written somewhere that is not its argument fails with R's own
   `could not find function "interactions"`.
   - Where dbarts can append a hint cheaply: at the three places it forces the caller's code. That
     is rule 2's forcing of a binding, rule 4's `..N` (after recovery also fails), and the
     evaluation of E itself. There, one handler matches the message against
     `gettextf("could not find function \"%s\"", n, domain = "R")` for each n in
     `names(dbartsForests)` and re-signals it with "; outside the argument that takes it, write
     dbartsForests$interactions(...)". This covers a wrapper formal passed on
     (`function(ints) dbarts(..., interactions = ints)` called with `interactions(...)`).
   - Where it cannot: a constructor evaluated in user code before any dbarts frame exists, such as
     `x <- interactions(...)` at top level, or a helper that forces through `list(...)` (the test
     `doFit*` helpers, and bartCause's `bartc`). R's plain message stands there, and
     dbartsForests.Rd and each constructor's page say so.
   - Every forcing point muffles "restarting interrupted promise evaluation". A failed promise is
     not cached, and a later site may force it again, side effects included.

Test matrix (new test file). "Obj" is `dbartsForests$interactions(max.order = 1)`. "Hint" is R's
could-not-find message with the dbartsForests hint. Each Obj row compares against the fit given Obj
directly.

| # | caller writes | caller binds | expect |
| --- | --- | --- | --- |
| 1 | `interactions = interactions(max.order = 1)`, direct | nothing | Obj |
| 2 | same, with an attached env exporting `interactions`, `blocks`, `forest` | attached only | Obj |
| 3 | `interactions = interactions` at top level | global `interactions <- Obj` | Obj |
| 4 | wrapper `function(interactions) dbarts(..., interactions = interactions)` | called with Obj | Obj |
| 5 | same wrapper | called with `interactions(max.order = 1)` | Hint |
| 6 | wrapper with formal default `NULL`, not supplied | - | same as no constraint |
| 7 | wrapper with formal, empty default, not supplied | - | same as no constraint |
| 8 | wrapper `function(k, x = interactions(max.order = k)) dbarts(..., interactions = x)` | not supplied | Hint (the default evaluates in the wrapper's frame) |
| 9 | `lapply(1, function(i, x) dbarts(..., interactions = x), x = interactions(max.order = 1))` | - | Hint |
| 10 | `lapply(1, function(i, ...) dbarts(..., interactions = ...), interactions(max.order = 1))` | - | Obj (`..N` recovery) |
| 11 | `interactions = interactions` | `interactions <- NULL` | no constraint |
| 12 | `blocks = blocks` | an unrelated function `blocks` | "requires 'groups'" (rule 3) |
| 13 | `variance = varianceForest` | an unrelated function | same as `varianceForest()` |
| 14 | `forests = list(forest(), forest(interactions = interactions))` | global Obj | forest 2 carries Obj |
| 15 | row 14 | nothing | refused, "see ?dbartsForests" |
| 16 | `forests = forest` | a list of forest objects named `forest` | that list |
| 17 | `viaDots(interactions = interactions(max.order = 1))`, and `viaNested` | - | Obj |
| 18 | a nested `..N`: `f <- function(...) dbarts(..., forests = list(forest(), forest(interactions = ..1)))` | - | Obj |
| 19 | `do.call(viaDots, list(interactions = quote(interactions(max.order = 1))))` | - | Obj |
| 20 | row 19 with `envir = new.env()` | - | Hint |
| 21 | formula term `forest(x1, interactions = interactions)` | global Obj | Obj |
| 22 | wrapper enclosed by a namespace (a child of `asNamespace("stats")`) writing `interactions = interactions` | global Obj | the global is NOT seen, so "needs at least one" (rule 3) |
| 23 | a site evaluated from a frame enclosed by the dbarts namespace | the namespace binds the constructor | constructor, not a value |
| 24 | wrapper `function(interactions, blocks) dbarts(..., blocks = blocks)`, `interactions` supplied as `stop("forced")` | - | fits; the formal is never forced (rule 2) |
| 25 | a `..N` whose expression errors, read by two sites | - | one error, no "restarting interrupted promise" warning |

## Constraints

- Neutral: every existing door-argument spelling fits bitwise as before.
- Per-site vocabularies, following xbart's subsetting of dbartsPriors:
  - `interactions` = {interactions}
  - `blocks` = {blocks}
  - `variance` = {varianceForest}
  - `forests` = {forest, interactions, blocks}
  - a formula term's knob arguments = {interactions, blocks}

  No priors or families are added. The existing class refusals downstream validate, and their
  messages gain "see ?dbartsForests".
- `bartcoreBCFSampler` and the other internal creators stay value-typed.
- Out of scope:
  - renaming, and a tombstone;
  - the priors' and families' bare-name shadowing;
  - any change to the forest-term grammar beyond the knob evaluation and dropping the
    `dbarts::forest` spelling (step 5).

## Steps

Commit grouping:
- Steps 1-4 may land first. They are neutral while the four are still exported, since dbarts
  sits on the search path, which the rule ignores.
- Steps 6-9 (NAMESPACE, the new Rd and pkgdown entry, the Rd edits, tests, benchmarks, NEWS) land
  together as one unexport commit.
- Step 5 lands with its test-formula-terms line.
- Step 11 (sisters) pushes alongside the unexport commit.

1. [R/family.R](../../R/family.R): add `evalInForestVocabulary` implementing the rule, about 60
   lines. It reuses `recoverForwardedArgument` for `..N` and adds the rule-5 hint handler.
2. [R/model.R](../../R/model.R):
   - Add the exported list `dbartsForests <- list(interactions = , blocks = , forest = , varianceForest = )`
     beside `dbartsPriors`.
   - Add `resolveForestArguments(matchedCall, evalEnv)`, which applies the rule to each of the
     four names in `matchedCall`. A missing name gives NULL, which is every door's default.
   - Rewrite the constructors' "Exported ..." comments.
   - Add "see ?dbartsForests" to the refusals in `resolveInteractions`, `resolveBlocks` and
     `resolveForests`.
3. Wire the sites:
   - [`dbarts`](../../R/dbarts.R): call the helper right after `resolveFamily`, before
     `ingestFormulaTerms`, and rebind the four locals.
   - [`dbartsSpec`](../../R/spec.R): the same, against `parentEnv`, before
     `forestBasisDeclarations`.
   - [`bart`](../../R/bart.R): the multinomial `unsupported` check reads `variance` through the
     helper.
4. [R/formulaTerms.R](../../R/formulaTerms.R): `finalizeTermForests` calls the internal `forest`
   twice, and `processHit` evaluates its knob arguments under the rule over {interactions, blocks}
   in the formula's environment.
5. `isForestCall` drops `quote(dbarts::forest)`, since the term grammar names `forest` only. The
   `dbarts::forest(x1 + x2)` term inside the first
   ["expectSameForest"](../../inst/tinytest/test-formula-terms.R) block becomes a plain `forest`.
   The alternative, keeping it, costs a spelling that works only inside a formula.
6. NAMESPACE: remove the four `export()` lines and add `export(dbartsForests)`.
7. Rd:
   - New man/dbartsForests.Rd, `\docType{data}`, in the shape of dbartsPriors.Rd. It lists the four
     under `\format`, each `\link`ed, and states the rule in two sentences.
   - interactions.Rd, blocks.Rd, forest.Rd and varianceForest.Rd each keep their `\alias` and
     `\usage`. `?interactions` and bartCause's links keep resolving, and `tools::codoc` checks
     namespace objects as well as exports. Each `\description` gains a sentence.
   - The varianceForest.Rd example's last line becomes `dbartsForests$varianceForest(...)`.
   - dbarts.Rd, bart.Rd and dbartsSpec.Rd argument items say the four resolve in the vocabulary.
   - _pkgdown.yml adds `dbartsForests`.
   - inst/NEWS.Rd's unreleased constraints entry gains one sentence.
8. Tests:
   - Rewrite the 71 outside sites. For `doFitInteractions` and `doFitBlocks`, the alternative is
     to forward `...` into `bart`.
   - Drop the 25 `dbarts::` prefixes.
   - New inst/tinytest/test-constructor-vocabulary.R, modelled on test-family-objects.R, covering:
     - the list's shape;
     - none of the four in `getNamespaceExports("dbarts")`;
     - bare equals list spelling at each door;
     - the matrix above. The mask rows use DIRECT calls, since forwarded dots carry no mask
       protection.
9. benchmarks/R: switch the 9 sites to the list spelling.
10. Docs: public-surface.md gets a one-sentence note in section 8 (dec-A67), plus a
    docs/plans/INDEX.md row. The TODO entry is removed at landing.
11. Sisters: the stan4bart and bartCause edits in Context.

## Verification

Neutral gates, against the slice's own library (`R_LIBS=<lib>`):
- `cd tests/cpp && make && ./test_bartcore` passes.
- `Rscript -e 'tinytest::test_package("dbarts")'` has zero failures.
- Bitwise equivalence (equivalence.R is edited, and CI's cpp-tests job runs the same command), on a
  `--preclean` reference build: `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-d49e2103.rds --bitwise --strict-coverage`
  reports every scenario identical.
- The exact gates reached through `variance = varianceForest(...)`, as
  .github/workflows/exact-gates.yaml invokes them, each exiting 0:
  - `Rscript benchmarks/R/aft-exact.R quick`
  - `Rscript benchmarks/R/aft-hetero-pit.R quick`
  - `Rscript benchmarks/R/heteroscedastic-exact.R quick`
- A smoke run of `benchmarks/R/composition-matrix.R` builds every cell.

Reviewer checklist (a NAMESPACE edit, R/ and man/ touched, a new Rd topic): each of these exits 0.
- `Rscript -e 'lintr::lint_package()'`
- `air format --check .`
- `Rscript -e 'pkgdown::check_pkgdown(".")'`
- a non-NULL `tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")`
- `Rscript tools/check-rc-codoc.R .`
- `Rscript tools/check-win-drift.R .`
- `Rscript tools/check-doc-freshness.R .`
- `R CMD check --as-cran` on a tarball built outside the tree, with no codoc or
  "Undocumented code objects" findings.

Probes:
- `?interactions` opens its page with dbarts attached.
- `dbarts::interactions` errors "not an exported object".
- `library(metafor)` then a direct `dbarts(..., forests = list(forest(), forest(basis = ~ z)))`
  fits identically to the same call without metafor.

Sister suites, against the slice's library:
- stan4bart: `Rscript -e 'tinytest::test_package("stan4bart")'`
- bartCause: `Rscript -e 'testthat::test_local("~/Repositories/bartCause")'`
- treatSens and bairrtt: one install and test run each, to confirm the zero-site finding.

## Landing

LANDED (pending hash), 2026-09-26, three commits: "Resolve the forest constructors by bare name
inside the arguments that take them" (steps 1-4), "Recognize a formula's forest() term by its bare
name only" (step 5), and "Unexport interactions, blocks, forest and varianceForest behind
dbartsForests" (steps 6-10). The sister edits (step 11) are a separate pass against this build.

Deviations:
- Rule 1 inlines the constructor as a call head only for a name the caller also claims in value
  position. Elsewhere a called name is bound to the constructor in the child environment, which
  gives the same mask protection and keeps a constructor's own error naming `interactions(...)`
  rather than a deparsed function.
- Rule 4 leaves a `..N` inside a function literal in E alone: it names that function's own dots,
  not the caller's.
- `dbarts` resolves the four just before `ingestFormulaTerms`, after the hurdle refusal, so that
  refusal still comes first.
- test-formula-terms case (16) became a refusal of a `dbarts::forest` head. The plain-`forest`
  rewrite would compare a formula with itself.
- Matrix row 10 is written with `..1`, not `interactions = ...`.
- The composition-matrix smoke run was not possible: the harness fails at its own table parse (it
  looks for a `bart2()` header that feature-matrix.md no longer has), at the base commit too. Its
  five edited sites are list spellings only.
- The metafor probe was not run (metafor is not installed). Matrix row 2's attached mask of
  `interactions`, `blocks` and `forest` covers the same case.

Gates, on the slice's library: full tinytest 8927/0; lintr, air; check-rc-codoc.R,
check-win-drift.R, check-doc-freshness.R; pkgdown::check_pkgdown; the NEWS parse gate; on a
`--preclean` reference build, equivalence.R 53 of 53 "identical draws (same RNG stream)" with no
"max |z|" line, bcf-equivalence 15/15 and multinomial-equivalence 11/11 bitwise; aft-exact.R,
aft-hetero-pit.R and heteroscedastic-exact.R in quick mode; `R CMD check --as-cran` on a tarball
built outside the tree.
