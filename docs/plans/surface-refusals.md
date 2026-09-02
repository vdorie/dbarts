# R-surface refusals, shapes and argument order (pre-RC)

Status: LANDED 2026-08-26 at d48aef8a (design record 37ab6ea9; see the landing note); designed at 9d0ee10f, revised the same day against an independent blind critique (findings 1
and 11 blocking; the adjudications are folded into sections 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15 and 17
rather than left as a change list). Every anchor below was read live at that tip.

Spec: the R-surface half of "breaking-if-deferred lands before the RC is declared" (VD). It extends
docs/plans/prerc-surface-freeze.md past its D1-D9 - the sweep that found these ran AFTER D1/D9/D5/D2 landed
(docs/plans/predict-surface.md) and D6/D8 landed (docs/plans/composition-refusals.md), and every finding is a
gap those slices did not cover rather than a regression they introduced. Four commits, one slice. No sampling
code moves, no RNG consumption changes on any path an existing test walks; the equivalence trio is expected
bitwise identical and zero baselines re-record. Budget: ~185 R + ~25 Rd + ~5 NEWS + ~175 test lines.

Wording is governed by docs/design/error-style.md (ADOPTED 2026-08-17, "every NEW message from this point on
follows this rule"). Every string this slice introduces is shown verbatim below and is born conformant: R1
single-quoted names, R2 lowercase initial, R3 no terminal period, R5 comma-concatenated interpolation, R6 one
main clause plus at most ONE extra (`:` for an explanation XOR `;` for a remedy, never both), R13's
`"<subject> does not support <feature>[: <reason>]"` for a composition refusal.

## 1. What this slice is, and what stays put

Three defect classes, all of them "the caller wrote something and the package answered as if they had not":

- SILENT NO-OP. A name lands in `...` (or in a formal only one arm reads) and is discarded. Fourteen sites.
- OUTPUT SHAPE. `type = "class"` is unreachable whenever `ci.level` is supplied, on four methods.
- POSITIONAL SPLIT. The same positional slot means two different things across sibling methods. Two sites.

Plus one input-validation asymmetry: a fractional `n.threads` is refused at `predict` and silently truncated
at fit construction.

Out of scope, stated so it is not rediscovered: the reference class's own `$predict`/`$predictForests`/
`$setTestPredictorAndOffset` vocabulary (`offset.test`, out of scope since D1 section 1); `makeind(all = )`,
kept as documented BayesTree signature compatibility (VD); `as_draws_array`/`as_draws_df.bartMultinomial`'s
`vars`, documented as inert rather than refused (VD, section 17); `run(n.threads)`'s documented-reserved
status, unchanged (VD); the variance-forest scale-leaf staleness and the `growFromRoot` all-zero alignment,
whose fix branch is chosen but whose code is deferred post-1.0 (VD) - recorded in section 17 only. And
man/dbarts.Rd, which this slice does not touch at all (section 15).

## 2. Census: where a name goes today

Read off the code at 9d0ee10f. `refuseUnusedGenericArgs` ([[R/generics.R:1908-1923@9d0ee10f]]) is the whole mechanism:
`intersect(names(reasons), names(dots))`, stop on `supplied[1L]` with
`"'<x>' is not used by <generic> on a <class> fit: <reason>"`. Under R6 that main clause plus the reason is
already the full budget, so a reason clause may contain no further `:` or `;`. Coverage is hand-listed per
class (`multinomialUnusedArgs` [[R/generics.R:1006@9d0ee10f]], `ordinalUnusedArgs` [[R/generics.R:1337@9d0ee10f]], `negbinUnusedArgs` [[R/generics.R:1614@9d0ee10f]], `hurdleUnusedArgs`
:2045) composed with three per-generic lists (`predictOffsetUnusedArgs` [[R/generics.R:254@9d0ee10f]], `predictWeightsUnusedArgs`
:281, `predictNoOffsetUnusedArgs` [[R/generics.R:291@9d0ee10f]]).

Methods with NO refusal call at all: `extract.bart` ([[R/generics.R:463-575@9d0ee10f]]; only the `trees` arm is guarded, by
`refuseTreesArguments` at [[R/generics.R:482@9d0ee10f]]), `extract.rbart` ([[R/generics.R:2426@9d0ee10f]]), `extract.dbartsSampler` ([[R/generics.R:2578@9d0ee10f]]), `fitted.bart`
([[R/generics.R:883@9d0ee10f]]), `fitted.rbart` ([[R/generics.R:2586@9d0ee10f]]), `residuals.bart` ([[R/generics.R:923@9d0ee10f]], `sample` only), `residuals.rbart` ([[R/generics.R:2643@9d0ee10f]], `sample`
only), `survivalProbabilities.bart` ([[R/bart.R:2547@9d0ee10f]]), `survivalProbabilities.rbart` ([[R/bart.R:2598@9d0ee10f]]),
`plotTree.dbartsSampler` ([[R/bart.R:2654@9d0ee10f]]).

The individual holes, each verified by reading the control flow:

| # | site | what is discarded |
|---|---|---|
| N1 | `predict.bart` [[R/bart.R:340-356@9d0ee10f]] | `ci.level`: `type = "forest"` returns at [[R/bart.R:348@9d0ee10f]], above the `ci.level` block at [[R/bart.R:427@9d0ee10f]] |
| N2 | `predict.bart` [[R/bart.R:361-373@9d0ee10f]], [[R/bart.R:386@9d0ee10f]] | `forest`: read only inside the `"forest"` arm; `predictBlend` ([[R/bart.R:808-818@9d0ee10f]]) and the plain path take it nowhere |
| N3 | `extract.bart` [[R/bart.R:496-497@9d0ee10f]] | `forest`, `contribution`: read only by `extractForest`; on `"ev"/"ppd"/"bart"/"loglik"` both are dropped |
| N4 | five `predict` methods [[R/bart.R:1216@9d0ee10f]], [[R/bart.R:1502@9d0ee10f]], [[R/bart.R:1754@9d0ee10f]], [[R/bart.R:2147@9d0ee10f]], [[R/bart.R:2194@9d0ee10f]] | `bases`: a formal only on `predict.bart` ([[R/bart.R:305@9d0ee10f]]) |
| N5 | all six `predict` methods | `sample`: never a predict formal, never refused |
| N6 | all six `predict` methods | `group.by`: a formal only on `predict.rbart` ([[R/bart.R:2204@9d0ee10f]]) |
| N7 | `residuals` on multinomial/ordinal/negbin ([[R/bart.R:1174@9d0ee10f]], [[R/bart.R:1478@9d0ee10f]], [[R/bart.R:1734@9d0ee10f]]) | `sample`: `refuseResidualsSample` ([[R/bart.R:913@9d0ee10f]]) is called only at [[R/bart.R:926@9d0ee10f]], [[R/bart.R:2133@9d0ee10f]], [[R/bart.R:2646@9d0ee10f]] |
| N8 | `residuals`, all six | `ci.level`: LIVE on bart/rbart/hurdle (forwarded through `...` into `fitted`, so `y` minus a 3-column band recycles into a 3-column object), dropped on the other three - the one live-and-working behavior this slice removes, argued as a fork in section 7 |
| N9 | all six `fitted` methods | `combineChains`: no `fitted.*` declares it; multinomial/ordinal/negbin/hurdle discard it outright, bart/rbart forward it into `extract` where it is provably value-neutral (both terminal reductions - `apply(result, length(dim(result)), mean)` [[R/bart.R:902@9d0ee10f]] and `posteriorInterval` [[R/bart.R:176-214@9d0ee10f]] - pool every margin but the last) |
| N10 | `survivalProbabilities.bart` [[R/bart.R:2547-2553@9d0ee10f]] | `group.by`, and every other fit-surface name; on the `.rbart` sibling `group.by` is a mandatory named formal ([[R/bart.R:2604@9d0ee10f]], enforced [[R/bart.R:2623-2625@9d0ee10f]]) |
| N11 | `plotTree.dbartsSampler` [[R/bart.R:2654-2656@9d0ee10f]] | `sample`/`chain` PARTIAL-MATCH onto the R5 method's `sampleNum`/`chainNum` ([[R/dbarts.R:2097-2102@9d0ee10f]]), where `plotTree.bart`/`.rbart` refuse them via `refusePlotTreeArgs` ([[R/dbarts.R:2664-2678@9d0ee10f]], called [[R/dbarts.R:2681@9d0ee10f]], [[R/dbarts.R:2702@9d0ee10f]]) |
| N12 | `extract.dbartsSampler` [[R/dbarts.R:2578-2584@9d0ee10f]] | everything but `type` |
| N13 | `fitted.bartHurdle` [[R/dbarts.R:2113-2122@9d0ee10f]] | `sample` is a bare `"train"` string (every sibling uses a choice vector) passed unvalidated into `extract`, which then validates it and answers in `extract`'s own voice |
| N14 | `fitted`/`extract`/`residuals`, all classes | `offset`, `weights`, `n.threads`, `newdata`: never formals there, never refused |

Shape defect S1: `predict.bartMultinomial` [[R/dbarts.R:1277-1287@9d0ee10f]], `predict.bartOrdinal` [[R/dbarts.R:1572-1583@9d0ee10f]],
`fitted.bartMultinomial` [[R/dbarts.R:1149-1156@9d0ee10f]] and `fitted.bartOrdinal` [[R/dbarts.R:1458-1465@9d0ee10f]] each take the `ci.level` return
BEFORE the class reduction, so `type = "class", ci.level = 0.9` silently returns the `type = "ev"` band. This
is currently DOCUMENTED ([[man/bart2.Rd:387@9d0ee10f]], [[man/bart2.Rd:391@9d0ee10f]] - the only two Rd promises of it) and PINNED
([[inst/tinytest/test-multinomial-generics.R:239-242@9d0ee10f]], [[test-ordinal.R:109-112@9d0ee10f]]), which is why section 5 argues it
rather than asserting it.

Order defects. `plot.bartMultinomial` is `(x, cols, plquants, ...)` ([[R/plot.R:226-231@9d0ee10f]]). The six sibling plot
methods that put `plquants` second are `plot.bart` ([[R/plot.R:53@9d0ee10f]]), `plot.rbart` ([[R/plot.R:123@9d0ee10f]]), `plot.bartOrdinal` ([[R/plot.R:308@9d0ee10f]]),
`plot.bartNegbin` ([[R/plot.R:368@9d0ee10f]]), `plot.bartHurdle` ([[R/plot.R:417@9d0ee10f]]) and `plot.pd2bart` ([[R/plot.R:505@9d0ee10f]], which carries no `cols` at all).
`plot.pdbart` is NOT one of them - it is `(x, xind, plquants, cols, ...)` ([[R/plot.R:476-482@9d0ee10f]]), a different
shape for a different object, and this slice does not touch it. And `fitted`'s positional slot 3 is `sample`
on bart ([[R/plot.R:883@9d0ee10f]]), hurdle ([[R/plot.R:2113@9d0ee10f]]) and rbart ([[R/plot.R:2586@9d0ee10f]]) but `ci.level` on multinomial ([[R/plot.R:1134@9d0ee10f]]), ordinal ([[R/plot.R:1437@9d0ee10f]]) and
negbin ([[R/plot.R:1700@9d0ee10f]]).

Input validation: `validatePredictThreads` ([[R/plot.R:234-248@9d0ee10f]]) refuses `n.threads != round(n.threads)`, while
`dbartsControl` coerces through `coerceOrError(n.threads, "integer")` ([[R/dbarts.R:244@9d0ee10f]]), whose `as.integer`
([[R/utility.R:167-178@9d0ee10f]]) truncates 2.7 to 2 without a warning; the S4 validity checks only length, NA and
positivity ([[R/A_class.R:295-296@9d0ee10f]], [[R/A_class.R:349-350@9d0ee10f]]). The legacy path pre-coerces the same way ([[R/bart.R:2780-2788@9d0ee10f]]).
Every other integer count in `dbartsControl` truncates identically - this is not an `n.threads` fact.

## 3. Mechanism: derived coverage, not longer hand-lists

The fork the sweep forces: extend the hand-written per-method lists to the methods that have none, or change
how coverage is computed.

RECOMMENDED: derive each method's coverage from the signatures. Define, once, the SURFACE VOCABULARY - every
name that is a formal of some `predict`/`extract`/`fitted`/`residuals`/`survivalProbabilities` method - with
one reason per (generic, name); a method's refusal list is that generic's table MINUS the method's own
formals, composed AFTER its class-specific list so the specific fact still wins.
`refuseUnusedGenericArgs` reports `names(reasons)`'s first hit, so composition order IS priority order,
exactly as the four existing call sites already rely on (`c(multinomialUnusedArgs, predictOffsetUnusedArgs,
predictWeightsUnusedArgs)` at [[R/bart.R:1232@9d0ee10f]]).

    # A name that is a formal on one method of this surface and not on another
    # is a caller mistake wherever it is foreign, not an argument to discard.
    # Deriving each method's list from its own formals - rather than listing
    # the foreign names by hand per class - is what keeps a name added to one
    # signature refused on every sibling that does not take it.
    foreignArgsFor <- function(reasons, own) {
      reasons[setdiff(names(reasons), own)]
    }

Each method then composes one line, e.g.

    refuseUnusedGenericArgs(
      list(...), "fitted", "bart",
      c(bartUnusedArgs, foreignArgsFor(fittedForeignReasons, names(formals(fitted.bart))))
    )

Cost: five reason tables (17 entries), one 3-line helper, 27 one-line call sites. Rejected alternatives:

- Hand-list the foreign names per method. ~27 literal lists, 2-10 names each. Same behavior, four times the
  text, and the sweep found fourteen holes in the four hand-lists that exist - the evidence that hand-lists
  do not stay complete is that they did not.
- A catch-all: refuse ANY unmatched name in `...`. Strictly stronger (it catches typos too), but it closes
  a DOCUMENTED forwarding path: `extract(type = "trees")` rewrites its own call onto the sampler's `getTrees`
  ([[R/bart.R:481-491@9d0ee10f]]), whose formals include `newdata` ([[R/dbarts.R:1999-2004@9d0ee10f]]); `refuseTreesArguments`'s own message
  advertises it verbatim ([[R/dbarts.R:456@9d0ee10f]], "'newdata' instead (see 'Extracting Trees' in ?bart)"), it is pinned at
  [[inst/tinytest/test-sampler-trees.R:60-70@9d0ee10f]] and documented at [[man/bart.Rd:221-223@9d0ee10f]]. It also closes
  `fitted` -> `extract` ([[man/bart.Rd:893@9d0ee10f]]) and `residuals` -> `fitted` ([[man/bart.Rd:927@9d0ee10f]]), so a catch-all at `extract` would fire
  with `extract` in the message on a call the user made to `residuals`. Since this slice does close every
  link and delete those three forwards (section 4), a terminal catch-all becomes a cheap ADDITIVE follow-on
  and is recorded as residue rather than taken here, where it would add a second message shape to a surface
  with 40 pinned messages of the first shape. The one place it IS taken is `extract.dbartsSampler`, whose
  `...` has no forwarding purpose and no delegating caller (section 4).

The `newdata`/`getTrees` collision is the reason section 4 states PLACEMENT, not only membership: on
`extract.bart` and `extract.rbart` the refusal call sits BELOW the `type == "trees"` branch.

## 4. The refusal lists, method by method

Wording family, unchanged: `refuseUnusedGenericArgs`'s `"'<x>' is not used by <generic> on a <class> fit:
<reason>"`. Under R6 that is main clause plus ONE explanatory clause, so every NEW reason below is a single
clause carrying no `:` and no `;`. `refuseResidualsSample` ([[man/bart.Rd:913-921@9d0ee10f]]), `refusePlotTreeArgs` ([[man/bart.Rd:2664-2678@9d0ee10f]]) and
`singleForestReason` ([[man/bart.Rd:997-1000@9d0ee10f]]) keep their strings VERBATIM - they are existing corpus, pinned by seven
tests between them, and error-style.md's slice L rewrites the corpus in one sweep so nothing is reworded
twice.

### 4a. Class lists, composed first

Existing and unchanged: `multinomialUnusedArgs` ([[man/bart.Rd:1006@9d0ee10f]]), `ordinalUnusedArgs` ([[man/bart.Rd:1337@9d0ee10f]]), `negbinUnusedArgs`
([[man/bart.Rd:1614@9d0ee10f]]), `hurdleUnusedArgs` ([[man/bart.Rd:2045@9d0ee10f]]), each holding `forest`/`contribution`.

Existing and unchanged, and the reason they cannot fold into a name-keyed table: the two offset lists are
CLASS-dependent, not generic-dependent. `predictOffsetUnusedArgs` ([[man/bart.Rd:254-256@9d0ee10f]], "this fit's out-of-sample offset
argument is named 'offset'") is right where `offset` IS a formal (predict.bart [[man/bart.Rd:300@9d0ee10f]], .bartMultinomial [[man/bart.Rd:1220@9d0ee10f]],
.bartNegbin [[man/bart.Rd:1758@9d0ee10f]]) and only `offset.test` is foreign; `predictNoOffsetUnusedArgs` ([[man/bart.Rd:291-294@9d0ee10f]]) is right where
neither is (predict.bartOrdinal, predict.bartHurdle). One name-keyed entry cannot hold both, so both lists
SURVIVE as composed-first class lists at their five existing call sites, and `offset`/`offset.test` appear in
no derived table.

NEW, because five cells in the additions table below have no class list to reach them:

    # Every own-class family has its own list; bart and rbart did not, so the
    # two names a fit-reduction never selects among had nowhere to be refused.
    bartUnusedArgs <- list(
      forest = "the reduction is over the combined location, in which every forest is already included",
      contribution = "the reduction is over the combined location, in which every forest is already included"
    )
    rbartUnusedArgs <- list(
      forest = singleForestReason,
      contribution = singleForestReason
    )

`bartUnusedArgs` carries new text because a `bart` fit MAY be amplitude-coupled and multi-forest, so
`singleForestReason` would be false there; it is composed on `fitted.bart` and `residuals.bart` only, the two
methods where both names are foreign (`forest` is a formal on `predict.bart` and `extract.bart`, where the
conditional refusal in 4c applies instead). `rbartUnusedArgs` reuses the existing clause verbatim - an rbart
fit is always single-forest - and is composed on all four rbart methods.

### 4b. The derived tables (verbatim strings)

`predictForeignReasons`
- `sample`: `"the fit's stored train and test channels are extract's 'sample'"`
- `weights`: `"this family's posterior-predictive draw takes no per-observation weight"` (the existing
  `predictWeightsUnusedArgs` clause, moved into the table verbatim; `foreignArgsFor` removes it on bart and
  rbart, where it is a formal)
- `bases`: `"only an amplitude-coupled multi-forest fit takes bases at the predicted rows"`
- `group.by`: `"'group.by' is the grouped (rbart_vi) fit's own predict argument"`
- `contribution`: `"the per-observation contribution decomposition belongs to extract(type = \"forest\")"`
  (reaches `predict.bart` alone - every other predict method's class list holds `contribution` and is
  composed first)

`extractForeignReasons`
- `ci.level`: `"extract returns the draws that fitted() and predict() take a band over"`
- `newdata`: `"predict(object, newdata) is the read at new rows"`
- `offset`, `weights`, `n.threads`, `bases`: `"extract reads stored channels and replays nothing"`
- `group.by`: `"the stored channels already carry the fit's own grouping"`

`fittedForeignReasons`
- `combineChains`: `"the per-chain draws are extract(object, combineChains = FALSE)"`
- `sample`: `"fitted values are always the fit's training rows"` - fires on `bartHurdle` alone, since
  `foreignArgsFor` removes it on bart and rbart (a formal there) and the three own-class
  `<class>FittedSampleReason` lists ([[man/bart.Rd:1121@9d0ee10f]], [[man/bart.Rd:1425@9d0ee10f]], [[man/bart.Rd:1688@9d0ee10f]]) are composed first
- `newdata`, `offset`, `weights`, `n.threads`, `bases`: `"fitted summarizes stored channels and replays
  nothing"`
- `group.by`: `"the stored channels already carry the fit's own grouping"`

`residualsForeignReasons`
- `ci.level`: `"residuals are the observed response minus the posterior-mean fit"` (section 7)
- `combineChains`: `"the per-chain draws are extract(object, combineChains = FALSE)"`
- `newdata`, `offset`, `weights`, `n.threads`, `bases`: `"residuals summarize stored channels and replay
  nothing"`
- `sample` is NOT in this table: `refuseResidualsSample` is called first on all six and keeps its own string.

`survivalProbabilitiesForeignReasons`
- `group.by`: `"'group.by' is the grouped (rbart_vi) fit's own argument"`
- `type`, `sample`, `ci.level`: `"survivalProbabilities returns the draws of S(t | x) at 'times'"`
- `offset`, `weights`, `n.threads`, `forest`, `contribution`, `bases`: `"survivalProbabilities takes 'times'
  and 'newdata' alone"`

### 4c. Additions per method, and where each call sits

Names newly refused; each method's existing coverage is not repeated. Unless a placement is named, the call
sits at the method head beside the existing validation.

| method | added | placement |
|---|---|---|
| `predict.bart` [[man/bart.Rd:296@9d0ee10f]] | `sample`, `contribution`, `group.by` | head ([[man/bart.Rd:313@9d0ee10f]], extend the existing call) |
| `predict.rbart` [[man/bart.Rd:2194@9d0ee10f]] | `sample`, `bases`, + `forest`/`contribution` via `rbartUnusedArgs` | head ([[man/bart.Rd:2219@9d0ee10f]]) |
| `predict.bartMultinomial` [[man/bart.Rd:1216@9d0ee10f]] | `sample`, `bases`, `group.by` | head ([[man/bart.Rd:1228@9d0ee10f]]) |
| `predict.bartOrdinal` [[man/bart.Rd:1502@9d0ee10f]] | `sample`, `bases`, `group.by` | head ([[man/bart.Rd:1512@9d0ee10f]]) |
| `predict.bartNegbin` [[man/bart.Rd:1754@9d0ee10f]] | `sample`, `bases`, `group.by` | head ([[man/bart.Rd:1765@9d0ee10f]]) |
| `predict.bartHurdle` [[man/bart.Rd:2147@9d0ee10f]] | `sample`, `bases`, `group.by` | head ([[man/bart.Rd:2157@9d0ee10f]]) |
| `extract.bart` [[man/bart.Rd:463@9d0ee10f]] | NEW call: `ci.level`, `newdata`, `offset`, `weights`, `n.threads`, `bases`, `group.by` | BELOW the `type == "trees"` branch (after [[man/bart.Rd:492@9d0ee10f]]), so `newdata` still forwards to `getTrees` |
| `extract.rbart` [[man/bart.Rd:2426@9d0ee10f]] | NEW call: the same seven, + `forest`/`contribution` via `rbartUnusedArgs` | BELOW the trees branch (after [[man/bart.Rd:2485@9d0ee10f]]), same reason |
| `extract.bartMultinomial` [[man/bart.Rd:1017@9d0ee10f]] | `ci.level`, `newdata`, `offset`, `weights`, `n.threads`, `bases`, `group.by` | head ([[man/bart.Rd:1030@9d0ee10f]]) - no trees arm |
| `extract.bartOrdinal` [[man/bart.Rd:1342@9d0ee10f]] | the same seven | head ([[man/bart.Rd:1351@9d0ee10f]]) |
| `extract.bartNegbin` [[man/bart.Rd:1619@9d0ee10f]] | the same seven | head ([[man/bart.Rd:1628@9d0ee10f]]) |
| `extract.bartHurdle` [[man/bart.Rd:2050@9d0ee10f]] | the same seven | head ([[man/bart.Rd:2059@9d0ee10f]]) |
| `extract.dbartsSampler` [[man/bart.Rd:2578@9d0ee10f]] | NEW call, its own catch-all (below) | head |
| `fitted.bart` [[man/bart.Rd:883@9d0ee10f]] | NEW call: `combineChains`, `newdata`, `offset`, `weights`, `n.threads`, `bases`, `group.by`, + `forest`/`contribution` via `bartUnusedArgs` | head |
| `fitted.rbart` [[man/bart.Rd:2586@9d0ee10f]] | NEW call: the same seven, + `forest`/`contribution` via `rbartUnusedArgs` | head |
| `fitted.bartMultinomial` [[man/bart.Rd:1134@9d0ee10f]] | `combineChains`, `newdata`, `offset`, `weights`, `n.threads`, `bases`, `group.by` | head ([[man/bart.Rd:1142@9d0ee10f]]) |
| `fitted.bartOrdinal` [[man/bart.Rd:1437@9d0ee10f]] | the same seven | head ([[man/bart.Rd:1444@9d0ee10f]]) |
| `fitted.bartNegbin` [[man/bart.Rd:1700@9d0ee10f]] | the same seven | head ([[man/bart.Rd:1707@9d0ee10f]]) |
| `fitted.bartHurdle` [[man/bart.Rd:2113@9d0ee10f]] | the same seven, plus `sample` (its formal is dropped - section 8) | head ([[man/bart.Rd:2121@9d0ee10f]]) |
| `residuals.bart` [[man/bart.Rd:923@9d0ee10f]] | `ci.level`, `combineChains`, `newdata`, `offset`, `weights`, `n.threads`, `bases`, + `forest`/`contribution` via `bartUnusedArgs` | after `refuseResidualsSample` ([[man/bart.Rd:926@9d0ee10f]]) |
| `residuals.rbart` [[man/bart.Rd:2643@9d0ee10f]] | the same, `bartUnusedArgs` replaced by `rbartUnusedArgs` | after [[man/bart.Rd:2646@9d0ee10f]] |
| `residuals.bartHurdle` [[man/bart.Rd:2129@9d0ee10f]] | `ci.level`, `combineChains`, `newdata`, `offset`, `weights`, `n.threads`, `bases` | extend the existing call ([[man/bart.Rd:2134@9d0ee10f]]) |
| `residuals.bartMultinomial` [[man/bart.Rd:1174@9d0ee10f]] | `sample` (NEW `refuseResidualsSample` call), plus the seven above | `refuseResidualsSample` first, then extend [[man/bart.Rd:1175@9d0ee10f]] |
| `residuals.bartOrdinal` [[man/bart.Rd:1478@9d0ee10f]] | the same eight | same shape ([[man/bart.Rd:1479@9d0ee10f]]) |
| `residuals.bartNegbin` [[man/bart.Rd:1734@9d0ee10f]] | the same eight | same shape ([[man/bart.Rd:1735@9d0ee10f]]) |
| `survivalProbabilities.bart` [[R/bart.R:2547@9d0ee10f]] | NEW call: `group.by`, `type`, `sample`, `ci.level`, `offset`, `weights`, `n.threads`, `forest`, `contribution`, `bases` | head, above the family check |
| `survivalProbabilities.rbart` [[R/bart.R:2598@9d0ee10f]] | NEW call: the same minus `group.by` | head |
| `plotTree.dbartsSampler` [[R/bart.R:2654@9d0ee10f]] | `refusePlotTreeArgs(sys.call())`, one line, the string the fit methods already raise | head |

`extract.dbartsSampler` gets the one catch-all this slice takes, because `refuseUnusedGenericArgs`'s "on a
`<class>` fit" noun is wrong for a sampler, its `...` forwards nowhere, and nothing delegates through it:

    refuseSamplerExtractArgs <- function(dots) {
      if (length(names(dots)) > 0L) {
        stop(
          "'", names(dots)[1L], "' is not used by extract on a dbartsSampler: ",
          "this method returns the sampler's coded predictor matrix"
        )
      }
      invisible(NULL)
    }

Three refusals are CONDITIONAL on `type`, because the name is a real formal on one arm. R13 shape, one
explanatory clause:

    # 'forest' selects among the per-forest channels only the "forest" arm
    # reports; every other arm has already recombined them into the reported
    # location, so a selection there would silently choose nothing.
    if (type != "forest" && !is.null(forest)) {
      stop(
        "type = \"", type, "\" does not support 'forest': every forest is ",
        "already recombined into the location it reports"
      )
    }

in `predict.bart` (above the `"forest"` arm at [[R/bart.R:340@9d0ee10f]]) and `extract.bart` (above [[R/bart.R:496@9d0ee10f]]), and in `extract.bart`
for `contribution` under `type != "forest" && isTRUE(contribution)`:

    stop(
      "type = \"", type, "\" does not support 'contribution': the ",
      "per-observation decomposition applies to the per-forest channel alone"
    )

Because every foreign name is refused at each link, three internal `...` forwards become dead and are
DELETED, which also removes the wrong-generic-name hazard the catch-all alternative fails on:
`fitted.bart` [[R/bart.R:893@9d0ee10f]] becomes `extract(object, type, sample)`, `fitted.rbart` [[R/bart.R:2605@9d0ee10f]] and [[R/bart.R:2631@9d0ee10f]] the same, and
`residuals.bart` [[R/bart.R:927@9d0ee10f]] / `residuals.rbart` [[R/bart.R:2647@9d0ee10f]] / `residuals.bartHurdle` [[R/bart.R:2140@9d0ee10f]] call `fitted.*` without `...`.

## 5. Fork: `type = "class"` with `ci.level`

Today the two are silently mutually exclusive - `type` is inert whenever `ci.level` is given.

Options.

(a) REFUSE the combination on all four methods. Cost ~16 R lines, 2 Rd sentences, 4 tests rewritten, 4 added.

(b) IMPLEMENT a class band. There is no interval on a nominal label; the nearest defined object is the
posterior distribution of the argmax (a per-observation probability vector over levels), which is a NEW
channel, not a band, and would need its own name, its own shape and its own Rd. Additive, and post-1.0 by
rule.

(c) KEEP and document. It is already documented ([[man/bart2.Rd:387@9d0ee10f]], [[man/bart2.Rd:391@9d0ee10f]]) and pinned
([[test-multinomial-generics.R:239-242@9d0ee10f]], [[test-ordinal.R:109-112@9d0ee10f]]), so this is the zero-cost option.

RECOMMEND (a). The rule the package already enforces everywhere else is that a supplied argument is honored
or refused, never inert; here `type = "class"` is inert, and the documentation of an inert argument is what
the sweep is closing on `fitted(combineChains = )` and on `bases`. (c)'s own consistency argument cuts the
other way: `fitted(type = "class")`'s class reduction ALREADY fires only at `ci.level = NULL` ([[test-ordinal.R:1152-1156@9d0ee10f]],
:1461-1465), so (a) states at `predict` the rule `fitted` follows in fact. And the direction is the
reversible one - refusing now leaves (b) additive later, while shipping (c) makes the ev-band-under-a-class-
request an interface promise that (b) would have to break.

Message, R13 verb, one explanatory clause:

    stop(
      "type = \"class\" does not support 'ci.level': a class prediction is ",
      "a label rather than a quantity with a credible band"
    )

This overturns nothing settled: the D1 landing pinned the ARM ORDER (the class reduction below the `ci.level`
block) as an implementation invariant, not the combination's meaning as an adjudicated fork.

## 6. Fork: `ci.level` on `type = "forest"`

`predict.bart(type = "forest")` returns at [[test-ordinal.R:348@9d0ee10f]] with a trailing forest margin (`predictForest` [[test-ordinal.R:690-696@9d0ee10f]]), so
a band is computable in two lines - `posteriorInterval(result, ci.level, trailing = 2L)` gives (obs x forest
x 3), the same widening the multinomial K margin already gets ([[test-ordinal.R:1150@9d0ee10f]], [[test-ordinal.R:1459@9d0ee10f]]).

Options: (a) REFUSE by name, in the shape of the `bases` refusal that already sits on this arm ([[test-ordinal.R:341-347@9d0ee10f]]);
(b) IMPLEMENT the (obs x forest x 3) band.

RECOMMEND (a). The arm's contract is stated in its own comment ([[test-ordinal.R:337-339@9d0ee10f]]) and in [[man/bart.Rd:204@9d0ee10f]] - it
"reports the raw per-forest total BEFORE any basis, which is what leaves the recombination to the caller".
Every other `ci.level` band summarizes the quantity the model reports; this one would summarize an
intermediate whose scale the amplitudes still carry, and a per-forest band read as a component's uncertainty
is exactly the misreading the arm's own `bases` refusal exists to prevent. Cost is a wash (6 R lines either
way), so the tie-breaker is direction: (a) leaves (b) additive, (b) does not leave (a) available.

Message, in R13's verb and in the exact shape of the sibling refusal ten lines above it in the same file
(`"type = \"forest\" does not support sample = \"test\": ..."`, [[man/bart.Rd:606-611@9d0ee10f]]):

    stop(
      "type = \"forest\" does not support 'ci.level': that arm reports each ",
      "forest's own total before any basis"
    )

## 7. Fork: `ci.level` on `residuals`

The one live, working behavior this slice removes, so it is argued rather than assumed. Today
`residuals(fit, ci.level = 0.9)` on bart, rbart and hurdle forwards through `...` into `fitted` ([[man/bart.Rd:893@9d0ee10f]], [[man/bart.Rd:2605@9d0ee10f]],
:2140) and returns `y` minus a 3-column `est`/`ci.lower`/`ci.upper` matrix, recycled down the columns; on
multinomial, ordinal and negbin the same call silently drops it (N8). No test and no Rd pins it anywhere
(verified zero hits).

Options: (a) REFUSE on all six; (b) plumb it into the three that drop it, making it live everywhere.

RECOMMEND (a), and the reason is what the object IS rather than what it costs. `y - band` is not a residual
band: subtracting an increasing triple from a constant reverses it, so the returned `ci.lower` column holds
the UPPER end of the residual's interval and vice versa - a labelled object whose labels are wrong on two of
three columns. Making that live on three more classes propagates a mislabelling; refusing it costs a caller
one line (`y - fitted(object, ci.level = )`, where the same reversal is at least visible at the call site).
Direction decides the rest: refusing leaves (b) additive, and a later (b) would want the columns re-labelled,
which shipping (a) does not prevent and shipping the status quo does.

Reason clause (composed into the standard main clause): `"residuals are the observed response minus the
posterior-mean fit"`.

## 8. Fork: `fitted()`'s positional slot 3

Target: ONE meaning for slot 3 on all six methods.

(a) `sample` third everywhere. Requires giving multinomial, ordinal and negbin a `sample` formal whose only
job is to be refused - which contradicts prerc-surface-freeze D1's settled sub-choice 1 ("refuse by name and
gain NO dummy formal"). Rejected.

(b) `ci.level` third everywhere; `sample` moves to slot 4 on the two methods that keep it. RECOMMENDED.

    fitted.bart(object, type = c("ev", "ppd", "bart"), ci.level = NULL,
                sample = c("train", "test"), ...)
    fitted.rbart(object, type = c("ev", "ppd", "bart", "ranef"), ci.level = NULL,
                 sample = c("train", "test"), ...)
    fitted.bartMultinomial(object, type = c("ev", "class", "bart"), ci.level = NULL, ...)
    fitted.bartOrdinal(object, type = c("ev", "class", "bart"), ci.level = NULL, ...)
    fitted.bartNegbin(object, type = c("ev", "ppd", "bart"), ci.level = NULL, ...)
    fitted.bartHurdle(object, type = c("ev", "ppd", "prob", "bart"), ci.level = NULL, ...)

(c) `ci.level` third and `sample` named-only after `...`, the `group.by` device D1 used. Rejected: D1 reached
for that device because `group.by` sat in a slot that meant something ELSE on five sibling classes; `sample`
means the same thing on both methods that carry it, so slot 4 is honest and a formal after `...` would be
machinery without a collision to justify it.

(d) Keep the split, document it. Rejected: it is the defect.

`fitted.bartHurdle` loses its `sample` formal outright (N13). A hurdle fit has no separate test channel at
all - `extract.bartHurdle` [[man/bart.Rd:2060-2065@9d0ee10f]] refuses `sample = "test"` unconditionally - so the argument has exactly
one legal value, and the three own-class siblings already refuse `sample` by name for the same reason. Its
body's `extract(object, type = type, sample = sample, ...)` ([[man/bart.Rd:2122@9d0ee10f]]) becomes `sample = "train"` literally,
`residuals.bartHurdle` [[man/bart.Rd:2140@9d0ee10f]] drops its `sample = "train"` pass, and `fittedForeignReasons$sample` (4b) is
what refuses a supplied one afterwards.

What breaks: `fitted(fit, "ev", "train")` now binds `"train"` to `ci.level` and raises `posteriorInterval`'s
`"'ci.level' must be a single number in (0, 1)"` ([[man/bart.Rd:184@9d0ee10f]]) - loud, not silent. Positional slot-3 `fitted` call
sites in the repo: ZERO (section 13's parse). Every in-repo `fitted(..., sample = )` names it:
[[test-generics-correctValues.R:53@9d0ee10f]], [[test-generics-posteriorPredictiveDistribution.R:125@9d0ee10f]],
[[test-rbart-generics.R:137@9d0ee10f]], [[test-rbart-groupby.R:184@9d0ee10f]], [[test-rbart-groupby.R:192@9d0ee10f]], [[test-rbart-groupby.R:242@9d0ee10f]], [[test-rbart-groupby.R:250@9d0ee10f]], [[test-rbart-groupby.R:292@9d0ee10f]], and the three refusal probes at
[[test-multinomial-generics.R:445@9d0ee10f]], [[test-ordinal.R:488@9d0ee10f]], [[test-nbinom.R:469@9d0ee10f]].

## 9. `plot.bartMultinomial`'s argument order

`(x, cols, plquants, ...)` -> `(x, plquants, cols, ...)`, matching the six sibling methods enumerated in
section 2 (`plot.pdbart` is a different shape and is not touched). Body unchanged. Zero in-repo positional
callers: test-multinomial-generics.R calls `plot(fitCombined)` bare at [[test-nbinom.R:132@9d0ee10f]], [[test-nbinom.R:133@9d0ee10f]], [[test-nbinom.R:481@9d0ee10f]], [[test-nbinom.R:483@9d0ee10f]], and the man
examples call `plot(bartFit)` bare ([[man/bart.Rd:397@9d0ee10f]]). One `\usage` line, [[man/bart2.Rd:147@9d0ee10f]]. This is the whole
item.

## 10. Fractional counts at construction

The asymmetry is one argument name with two behaviors one level apart, so the fix belongs where the
truncation is, not at `n.threads` alone: every `coerceOrError(x, "integer")` truncates.

RECOMMENDED: refuse a non-integral double inside `coerceOrError`'s integer branch ([[R/utility.R:160-179@9d0ee10f]]) -
one edit, no call-site changes, and `coerceOrError(x, "numeric")` (xbart's proportional `n.test`,
[[R/xbart.R:385@9d0ee10f]]) is untouched by construction.

The message's subject cannot be `mc[[2L]]` unguarded: at [[R/spec.R:544-548@9d0ee10f]] the argument is
`if (is.null(varianceNTrees)) 40L else varianceNTrees`, and deparsing a call into a quoted name renders
garbage. Guard it:

    func <- switch(type, logical = as.logical, integer = as.integer, numeric = as.numeric)
    # as.integer TRUNCATES a fractional double silently, so a mistyped count
    # would take a different value than the caller wrote; the per-call thread
    # argument already refuses one, and the two must not disagree. The subject
    # is the caller's own spelling only when the argument arrived as a bare
    # name - an expression has no name to quote.
    if (type == "integer" && is.double(x) && any(is.finite(x) & x != trunc(x))) {
      subject <- if (is.symbol(mc[[2L]])) paste0("'", mc[[2L]], "'") else "value"
      stop(subject, " must be a whole number, not ", deparse(x)[1L])
    }
    result <- tryCatch(func(x), warning = function(e) e)

Renderings: `dbartsControl(n.threads = 2.7)` gives `"'n.threads' must be a whole number, not 2.7"`;
`dbartsSpec(variance = varianceForest(n.trees = 40.5))` gives `"value must be a whole number, not 40.5"`,
which names no wrong argument. The predicate falls through correctly at the edges: `NA`, `Inf` and `2^31` are
not `is.finite() & x != trunc(x)` or are caught by the existing coercion error; `200` passes; `200.5`
refuses; `is.double` excludes integer and logical inputs.

Coverage acquired at a stroke - 32 integer call sites: `dbartsControl`'s ten counts ([[R/dbarts.R:239-248@9d0ee10f]]),
`dbarts()`'s own `seed` mirror ([[R/dbarts.R:611@9d0ee10f]]), `bart()`'s nine legacy spellings ([[R/bart.R:2780-2788@9d0ee10f]], which
pre-coerce precisely so the message names the caller's own argument), `rbart_vi`'s `n.chains`/`n.threads`
([[R/rbart.R:74@9d0ee10f]], [[R/rbart.R:79@9d0ee10f]]), `xbart`'s SEVEN ([[R/xbart.R:43@9d0ee10f]] `n.cuts`, [[R/xbart.R:51@9d0ee10f]] `n.thin`, [[R/xbart.R:197@9d0ee10f]] `n.trees`, [[R/xbart.R:370@9d0ee10f]] `n.test`,
:401 `n.reps`, [[R/xbart.R:405@9d0ee10f]] `n.burn`, [[R/xbart.R:409@9d0ee10f]] `n.threads`), `dbartsSpec`'s `n.trees` ([[R/spec.R:545@9d0ee10f]]) and `seed` ([[R/spec.R:841@9d0ee10f]]),
and `max.leaf.size` ([[R/model.R:1475@9d0ee10f]]).

Alternative rejected: an `n.threads`-only check in `dbartsControl` reusing `validatePredictThreads`'s exact
sentence. Cheaper by nothing, and it leaves `n.trees = 200.7` truncating next to `n.threads = 2.7`
refusing - recreating one level down the very asymmetry being fixed.

Wording: `validatePredictThreads` keeps its richer `"'n.threads' must be a single positive integer, not
2.7"` (it also owns length and positivity). Same R11/R12 family - subject named, constraint stated, offending
value echoed after a `, not` - and the S4 validity ([[R/A_class.R:295-296@9d0ee10f]], [[R/A_class.R:349-350@9d0ee10f]]) already owns length, NA
and positivity at construction, so the constructor's message must not claim them.

## 11. The two R5 readers with a vestigial `result`

VD ruling: refuse `result` by name on `$getSigmas` and `$getSumsOfSquaredResiduals` ([[R/dbarts.R:1675@9d0ee10f]], [[R/dbarts.R:1686@9d0ee10f]];
[[man/dbartsSampler-class.Rd:290-291@9d0ee10f]] "Accepted but not used"). It holds positional slot 1, and honoring it
later would write into a caller's buffer.

Implementation: keep the formal, refuse when supplied.

    getSigmas = function(result) {
      "Return current residual error term on original, standard deviation scale."
      # the formal is held so it cannot be quietly repurposed: this reader
      # allocates its own vector, and filling a caller's buffer in place is
      # getLatents' contract alone
      if (!missing(result)) {
        stop("'result' is not used by getSigmas: this reader allocates its own vector")
      }
      ...

and the same with `getSumsOfSquaredResiduals` in place of `getSigmas`. Main clause plus one explanatory
clause; the cure (`getLatents`) is one clause too many under R6 and lives in the Rd instead.

Alternative rejected: making both nullary (`function()`), matching `getDispersion` ([[man/dbartsSampler-class.Rd:1681@9d0ee10f]]). Shorter, but R's
own "unused argument" wall carries no reason and no stable string to pin, and the Rd usage lines
([[man/dbartsSampler-class.Rd:89@9d0ee10f]], [[man/dbartsSampler-class.Rd:92@9d0ee10f]]) would have to lose an argument the `\item{result}` block still
describes for `getLatents`.

In-repo callers passing an argument to either: none. 46 call sites across R/, inst/tinytest/, benchmarks/ and
vignettes/, all nullary.

## 12. Rd plan

- [[man/bart.Rd:45-53@9d0ee10f]] (`predict.bart` usage) unchanged; [[man/bart.Rd:55-70@9d0ee10f]] (`extract`/`fitted` usage) - `fitted.bart`'s
  `sample`/`ci.level` swap, section 8. [[man/bart.Rd:206@9d0ee10f]] `\item{sample}` gains a sentence: it is `extract`'s and
  `fitted`'s only, and is refused by name on `predict`. [[man/bart.Rd:209@9d0ee10f]] `\item{forest}` gains the conditional refusal (a
  selection outside `type = "forest"`). [[man/bart.Rd:218@9d0ee10f]] `\item{ci.level}` gains the `type = "forest"` refusal (section
  6) and the `residuals` refusal (section 7). [[man/bart.Rd:221-223@9d0ee10f]] `\item{\dots}` currently says "Not used in
  `predict`" - restate as: not used by `predict`, `fitted` or `residuals`, and a name belonging to a sibling
  method is refused rather than ignored; the `type = "trees"` forwarding to `getTrees` (which is where
  `newdata` remains legal) is the one exception. [[man/bart.Rd:247@9d0ee10f]] `\item{residuals}` prose gains the `ci.level`/`sample`
  refusal. One sentence on the legacy count items (`ntree`, `nskip`, `ndpost`, `numcut`, `keepevery`,
  `printevery`, `nthread`, `nchain`, `printcutoffs`): a fractional value is refused, not truncated.
- [[man/bart2.Rd:72-145@9d0ee10f]] - the six own-class `\usage` blocks: `fitted.bartHurdle` loses `sample`, all six
  `fitted` blocks carry `ci.level` third. [[man/bart2.Rd:147@9d0ee10f]] `plot.bartMultinomial` reordered. [[man/bart2.Rd:344@9d0ee10f]] `\item{sample}` and
  [[man/bart2.Rd:359@9d0ee10f]] `\item{ci.level}` gain the refusals. The four "Generics for a X fit" paragraphs ([[man/bart2.Rd:387@9d0ee10f]], [[man/bart2.Rd:391@9d0ee10f]], [[man/bart2.Rd:395@9d0ee10f]],
  [[man/bart2.Rd:397@9d0ee10f]]) each lose the sentence promising an ev band under `type = "class"` (only [[man/bart2.Rd:387@9d0ee10f]] and [[man/bart2.Rd:391@9d0ee10f]] carry it) and
  gain the `combineChains`-on-`fitted` and `ci.level`-on-`residuals` refusals in the sentence that already
  enumerates that class's refused names.
- [[man/rbart.Rd:53-62@9d0ee10f]] - `fitted.rbart`'s usage swap; [[man/rbart.Rd:107@9d0ee10f]] `\item{sample}`, [[man/rbart.Rd:110@9d0ee10f]] `\item{ci.level}` as bart's.
- [[man/survivalProbabilities.Rd:83-85@9d0ee10f]] `\item{...}` currently reads "Not used; for compatibility with the
  generic" - restate as refused by name; [[man/survivalProbabilities.Rd:71@9d0ee10f]] `\item{group.by}` gains "refused by name on the `bart` method,
  which has no grouping".
- [[man/plotTree.Rd:63@9d0ee10f]] `\item{\dots}` gains the `sample`/`chain` refusal, now that it holds on the sampler
  method too. Its `\usage` at [[man/plotTree.Rd:33@9d0ee10f]] (`plotTree.dbartsSampler(object, \dots)`) is unchanged.
- [[man/dbartsSampler-class.Rd:284-291@9d0ee10f]] `\item{result}` - "Accepted but not used by `getSigmas` and
  `getSumsOfSquaredResiduals`" becomes "Refused by name by ...", naming `getLatents` as the reader that fills
  a buffer in place (the clause R6 kept out of the message).
- man/dbartsControl.Rd - one sentence on the count items, the SHARED home for the whole-number rule; the
  fitting pages cross-reference it rather than repeating it, except man/bart.Rd's legacy spellings above,
  which `dbartsControl` does not name.
- man/xbart.Rd, man/rbart.Rd, man/dbartsSpec.Rd - no count sentence: each already defers to
  `\link{dbartsControl}`/`\link{dbarts}` for the shared counts. man/dbarts.Rd is NOT touched (section 15).
- inst/NEWS.Rd `\subsection{UPGRADING}` (:5) - four `\item`s, one per commit: the refusal coverage (naming
  the classes of name now refused, and that `residuals(ci.level = )` is among them), the `type = "class"` and
  `type = "forest"` `ci.level` refusals, the `fitted` and `plot.bartMultinomial` argument-order changes, and
  the whole-number count rule plus the two sampler readers. The file holds 396 `\item` entries at this tip
  and should hold 400 after - a landing-note observation, not a gate.

## 13. Test plan

House battery per commit: `R CMD INSTALL .`, `tinytest::test_package("dbarts")`, `tools/check-rc-codoc.R`
(commits 1 and 3 especially - commit 3 moves formals), `Rscript tools/check-doc-freshness.R .`,
`air format --check .`, `lintr::lint_package()`, `R CMD check --as-cran`. No C/C++ moves at all, so
`tests/cpp` must stay at its current count on every commit and cannot change.

Verification scope for every "no call site does X" claim below: R/, inst/tinytest/, tests/, benchmarks/,
tools/, AND vignettes/ and the `\examples` in man/ - the last two because `R CMD check` executes them. Swept
clean: the only affected-surface calls there are `plot(bartFit)` ([[man/bart.Rd:397@9d0ee10f]]), `plot(pdb1, ylim = )`
([[man/pdbart.Rd:159@9d0ee10f]]), `plotTree(fit, treeNum = , sampleNum = )` ([[man/plotTree.Rd:97@9d0ee10f]]), a nullary `$getSigmas()`
(vignettes/dbarts-as-a-component.Rmd) and a named-only `$plotTree` (vignettes/working_with_saved_trees.Rmd).

Existing tests that MUST change:

- [[inst/tinytest/test-multinomial-generics.R:239-242@9d0ee10f]] and [[inst/tinytest/test-ordinal.R:109-112@9d0ee10f]] - the two
  blocks asserting `predict(fit, x, type = "class", ci.level = 0.9)` EQUALS the `type = "ev"` band. Under
  section 5 they become `expect_error(..., "does not support 'ci.level'")`, and the comment above each
  ([[inst/tinytest/test-ordinal.R:229-230@9d0ee10f]], [[inst/tinytest/test-ordinal.R:103@9d0ee10f]]) is rewritten to state the rule rather than the arm order.
- [[inst/tinytest/test-hurdle.R:80-82@9d0ee10f]] - `residuals(h, sample = "train")` already expects the
  `refuseResidualsSample` message; unchanged, but verify it still fires FIRST now that
  `residualsForeignReasons` is composed behind it.
- No other existing assertion changes. Verified by reading: every `fitted(..., sample = )` names the
  argument (section 8's list); no test passes `combineChains` or `ci.level` to `fitted` or `residuals`; no
  test calls `plot` on a multinomial fit with a positional second argument; no test passes an argument to
  `$getSigmas` or `$getSumsOfSquaredResiduals` (46 call sites, all nullary); and the ONLY fractional count
  anywhere in the repo is [[inst/tinytest/test-generics-multithreaded.R:313@9d0ee10f]]'s
  `predict(gaussianFit, friedman$x, n.threads = 2.7)`, whose message section 10 leaves unchanged. The
  implementer must still run the full suite - this is a reading, and several files hardcode values that
  depend on their full execution history.

New tests, all `expect_error(..., fixed = TRUE)` on the full message where it is short and stable, `pattern =`
on a short stem otherwise (`"is not used by"`, `"does not support"`, `"must be a whole number"`) - the
discipline the four existing refusal files already use.

- inst/tinytest/test-generics-errors.R, extending its per-class blocks ([[inst/tinytest/test-generics-multithreaded.R:110-300@9d0ee10f]]). Count-exact: 6 methods x
  `sample` = 6; 5 methods x `bases` = 5; 6 methods x `group.by` (5 refused, 1 accepted on rbart) = 6; the
  `predict.bart` `contribution` and `predict.rbart` `forest`/`contribution` = 3. **20 new assertions.**
- inst/tinytest/test-generics-errors.R, extract arm: `ci.level`, `newdata`, `n.threads` refused on
  `extract.bart` and `extract.rbart` = 6; one representative name on each of the four own-class extracts = 4;
  `extract(sampler, sample = "train")` and `extract(sampler, combineChains = FALSE)` = 2; and the
  ACCEPTANCE half that finding 1 turns on - `extract(fit, type = "trees", newdata = x)` still returns a data
  frame whose `n` column counts the supplied rows, proving the refusal sits below the trees branch = 1.
  **13 assertions.**
- inst/tinytest/test-generics-errors.R, fitted/residuals arm: `combineChains` refused on all six `fitted` =
  6; `ci.level` refused on all six `residuals` (section 7) = 6; `sample` refused on the three `residuals`
  methods that did not refuse it = 3; `n.threads` on one `fitted` and one `residuals` = 2. **17 assertions.**
- inst/tinytest/test-rbart-aft.R (or test-generics-errors.R): `survivalProbabilities(bartAftFit, times,
  newdata = x, group.by = g)` refused naming `group.by`; the same call on the rbart fit still works; one
  foreign name (`type`) refused on each of the two methods. **4 assertions.**
- inst/tinytest/test-plot-generics.R, beside [[inst/tinytest/test-generics-multithreaded.R:127-142@9d0ee10f]]: `plotTree(sampler, sample = 5L)` and
  `plotTree(sampler, chain = 2L)` raise the same two strings the fit methods raise, and
  `plotTree(sampler, sampleNum = 5L)` still draws. **3 assertions.**
- inst/tinytest/test-multinomial-generics.R and test-ordinal.R: the two rewritten blocks (above) plus, on
  `fitted`, `type = "class", ci.level = 0.9` refused. **4 assertions** (2 rewritten, 2 new).
- inst/tinytest/test-predict-blend.R, at the END of the file (it is seeded and several blocks depend on the
  file's execution history): `predict(fit, xNew, type = "forest", ci.level = 0.9)` refused;
  `predict(fit, xNew, type = "ev", forest = 1L)` refused; `extract(fit, type = "ev", forest = 1L)` refused;
  `extract(fit, type = "ev", contribution = TRUE)` refused; and, the acceptance half,
  `predict(fit, xNew, type = "forest", forest = 1L)` still returns the one-forest slice unchanged.
  **5 assertions.**
- inst/tinytest/test-generics-intervals.R: `fitted(fit, "ev", 0.9)` is `identical` to
  `fitted(fit, type = "ev", ci.level = 0.9)` on all six classes, plus `fitted(fit, "ev", "train")` raising
  the `ci.level` validity message. **7 assertions**, and their power is ASYMMETRIC, stated so nobody reads
  them as six equal guards: on bart, rbart and hurdle they genuinely discriminate (before the swap the
  positional `0.9` reached `validateSample` [[inst/tinytest/test-generics-multithreaded.R:1893-1898@9d0ee10f]], which refuses a non-character), while on
  bartMultinomial, bartOrdinal and bartNegbin `ci.level` is ALREADY third ([[inst/tinytest/test-generics-multithreaded.R:1136@9d0ee10f]], [[inst/tinytest/test-generics-multithreaded.R:1439@9d0ee10f]], [[inst/tinytest/test-generics-multithreaded.R:1702@9d0ee10f]]), so those
  three assert nothing this slice could break and are kept as plain regression guards.
- inst/tinytest/test-plot-generics.R: `plot(multinomialFit, c(0.1, 0.9))` runs (a `plquants` in slot 2), one
  `expect_silent`. **1 assertion.**
- inst/tinytest/test-control-errors.R: `dbartsControl(n.threads = 2.7)`, `dbarts(..., n.threads = 2.7)`,
  `bart2(..., n.trees = 50.5)`, `bart(..., nthread = 1.5)`, `rbart_vi(..., n.chains = 2.5)`,
  `xbart(..., n.trees = 50.5)` each refused naming the caller's own spelling; and, the acceptance half,
  `n.threads = 2` as a double still constructs (an integral double is not fractional). **7 assertions.**
- inst/tinytest/test-sampler-residuals.R and test-bartcore.R: `$getSigmas(numeric(1))` and
  `$getSumsOfSquaredResiduals(numeric(1))` refused; both nullary reads unchanged. **4 assertions.**

Total: **85 new or rewritten assertions**, split 57 / 9 / 8 / 11 across the four commits. No new warning is
emitted anywhere, so no `expect_warning` count moves. No new assertion draws: every one is a refusal or a
shape check on channels the file already built, so none may sit above an RNG-locked block - each goes at the
END of its file or carries its own `set.seed`.

Regression risk of the two reorders, stated plainly: an argument-order change is invisible to a suite whose
calls are all named, which this one is - so the suite CANNOT catch a mistake here, and the three
discriminating identity assertions above are the only thing that does.

## 14. Consumer sweep (read-only)

Four repos, swept for every affected surface INCLUDING `do.call` and args-list construction (the site class a
previous sweep missed - docs/plans/predict-surface.md's landing note).

- bartCause, /Users/vdorie/Repositories/bartCause branch `dbarts-1.0` (d825cfc). R API only, no `src/`.
  `predict` on a dbarts fit at [[R/generics.R:162@9d0ee10f]], [[R/generics.R:169@9d0ee10f]], [[R/generics.R:176@9d0ee10f]]: object and `newdata` positional, everything
  after named (`group.by =`, `combineChains =`). The `do.call` path is [[R/generics.R:193@9d0ee10f]], [[R/generics.R:199@9d0ee10f]], [[R/generics.R:202@9d0ee10f]], which
  build `list(object$fit.rsp, x.new, group.by = , combineChains = FALSE, ...)`, mutate positional slot 2 at
  [[R/generics.R:222@9d0ee10f]] and [[R/generics.R:229@9d0ee10f]], and splat at [[R/generics.R:212@9d0ee10f]], [[R/generics.R:223@9d0ee10f]], [[R/generics.R:230@9d0ee10f]] - again positional only through slot 2, which this slice does
  not move. `extract` at [[R/responseFit.R:244@9d0ee10f]], [[R/responseFit.R:245@9d0ee10f]] and [[R/treatmentFit.R:189@9d0ee10f]] (and
  [[tests/testthat/test-05-generics.R:28@9d0ee10f]], [[tests/testthat/test-05-generics.R:29@9d0ee10f]]) SKIPS `type` and names `sample`/`combineChains` - `extract`'s
  signature is unchanged here, so those hold. No `fitted`, `residuals`, `plot`, `plotTree` or
  `survivalProbabilities` on a dbarts fit at all; no argument to `$getSigmas`/`$getSumsOfSquaredResiduals`.
  MIGRATION LINES: ZERO.
- stan4bart, /Users/vdorie/Repositories/stan4bart branch `bartcore` (f9bca65). R API plus `LinkingTo`.
  Exactly two hits on the whole affected surface, both `fitted` on a dbarts fit and both NAMED with `type`
  skipped: [[inst/tinytest/test-02-binary.R:60@9d0ee10f]] `fitted(base_bart_fit, sample = "test")` and [[inst/tinytest/test-02-binary.R:68@9d0ee10f]]
  `fitted(rbart_fit, sample = "test")`. Under section 8 `sample` stays a formal (slot 4) and stays named, so
  both hold. Nothing else: no `predict`/`extract`/`residuals`/`plot`/`plotTree`/`survivalProbabilities` on a
  dbarts fit, no argument to either reader, no `do.call` into any dbarts entry.
  `makeTestModelMatrix(object$bartData, mf.bart)` ([[R/new_data.R:135@9d0ee10f]]) is positional but that function's
  signature is untouched. MIGRATION LINES: ZERO.
- treatSens. The `master` branch (d1da1dd) still targets the DELETED pre-1.0 C++ ABI
  (`R_GetCCallable` on `dbarts::BARTFit`, [[src/bartTreatmentModel.cpp:38-52@9d0ee10f]]), so it is not a 1.0 consumer and
  is out of scope; the `dbarts-1.0` branch (1db3d89) is the one that matters and is on the flat C API
  ([[src/R_interface.cpp:31@9d0ee10f]]). Zero hits on every affected R surface: no `predict`/`fitted`/`extract`/
  `residuals`/`plot`/`plotTree`/`survivalProbabilities` on a dbarts object, no argument to either reader, no
  `do.call` into any dbarts entry. MIGRATION LINES: ZERO.
- bairrtt, /Users/vdorie/Repositories/bairrtt branch `main` (6167423). R API only. Its six predict calls
  ([[R/irt_causal_bart.R:568@9d0ee10f]], [[R/irt_causal_bart.R:571@9d0ee10f]], [[R/irt_causal_bart.R:636@9d0ee10f]], [[R/irt_causal_bart.R:637@9d0ee10f]], [[R/irt_causal_bart.R:709@9d0ee10f]], [[R/irt_causal_bart.R:713@9d0ee10f]]) are `model$predict(frame)` on the reference class,
  which this slice does not touch, and it passes no other argument to any of them. It re-exports
  `dbarts::extract` ([[R/diagnostics.R:69-71@9d0ee10f]], NAMESPACE:18) but every `extract` call targets its own
  `irt_causal_fit` method. No argument to either reader. MIGRATION LINES: ZERO.

Total lockstep migration cost: zero lines in any consumer.

Section 10's whole-number rule is a different matter from migration LINES: it converts silent truncation into
an error on every unlaundered pass-through, and three consumers have them. Full exposure list, all four
verified by reading:

- bartCause carries its OWN verbatim copy of the pre-change `coerceOrError` (bartCause [[R/utility.R:1-12@9d0ee10f]]), so
  the counts it launders itself (`n.samples`, `n.burn`, `n.chains`, `n.threads`, `n.trees`, `seed`,
  [[R/bcf.R:102-105@9d0ee10f]], [[R/bcf.R:210-212@9d0ee10f]]) keep TRUNCATING after this slice - its behavior drifts from dbarts's, and that
  is its maintainer's fix, not ours. What is NOT laundered is `extraControl`
  ([[R/bcf.R:210-212@9d0ee10f]]) and the `samplerCall` fill loop ([[R/bcf.R:232-238@9d0ee10f]]), both splatted into
  `dbarts::dbartsControl`/`dbarts::dbarts`, so `bcf(n.cuts = 100.5)` moves from truncate to error.
- stan4bart fills `bart_args` names into the control call ([[R/stan4bart_fit.R:517-520@9d0ee10f]]) and the spec call
  ([[R/stan4bart_fit.R:554-560@9d0ee10f]]) without coercion, so `stan4bart(bart_args = list(n.trees = 75.5))` moves to error; the same for
  `mvbart(n.trees = 75.5)` ([[R/mvbart.R:120@9d0ee10f]], `n.trees` uncoerced).
- treatSens `dbarts-1.0` passes `nburn` straight to `pdbart`'s `nskip` ([[R/treatSensBART.R:155@9d0ee10f]], [[R/treatSensBART.R:159@9d0ee10f]]), so
  `treatSens.BART(nburn = 200.5)` moves to error. Its own `nthreads` validation ([[R/treatSensBART.R:90-97@9d0ee10f]])
  warns and truncates before dbarts sees it, so that one argument is unaffected.
- bairrtt truncates `n_threads` with its own `as.integer` ([[R/irt_causal_bart.R:433@9d0ee10f]]) before `dbartsControl`,
  so it is unaffected.

None of those is a mistyped-count regression a suite would hit by accident (every in-repo consumer test
passes integer literals), but "invisible to them" would be false, and the NEWS UPGRADING item should say
"refused, not truncated" plainly for exactly this reason.

Gate: run bartCause's and stan4bart's suites after commit 3 (the reorders) and again after commit 4 (the
count rule); treatSens `dbarts-1.0` and bairrtt need only the branch's habitual rebuild.

Two observations, neither a dbarts change:

- bartCause's own generics carry the same positional-split defect this slice fixes in dbarts:
  `extract.bartcFit`'s slot 3 is `sample`, matched positionally at [[R/summary.R:232@9d0ee10f]], [[R/summary.R:235@9d0ee10f]], [[R/summary.R:236@9d0ee10f]], [[R/summary.R:238@9d0ee10f]], [[R/summary.R:239@9d0ee10f]]
  and [[tests/testthat/test-05-generics.R:287@9d0ee10f]], [[tests/testthat/test-05-generics.R:320@9d0ee10f]], while `extract.bartBCF` ([[R/bcf.R:500@9d0ee10f]]) puts `forest` in
  that slot. Its maintainer's call, exactly as `predict.bartcFit`'s shape was at D1.
- treatSens `dbarts-1.0` reaches two UNEXPORTED dbarts internals through `asNamespace("dbarts")`:
  `estimateSigmaFromLinearModel` ([[R/cibart.R:64@9d0ee10f]]) and `parsePriors` ([[R/cibart.R:76-80@9d0ee10f]], called with two
  positional arguments). Neither is touched here, but neither is a contract; recorded so the next slice that
  moves either knows it has a caller.

## 15. Relation to the two parallel slices

- capi-shape.md is effectively DISJOINT from this one: its only R-side citation is [[R/dbarts.R:1968@9d0ee10f]],
  which this slice does not move. The two share exactly two artifacts, and both are integration-time rather
  than authoring-time: the docs/design anchor re-alignment (run LAST, ONCE, over the stacked tree rather than
  per slice, since two independent line maps cannot be applied twice to the same anchors) and the
  freshness/codoc gates. inst/NEWS.Rd is a shared file but an append-only one; the two slices' `\item`
  additions are resolved at integration stacking, not by either author.
- rd-records.md OVERLAPS substantially - R/generics.R, R/A_class.R, R/spec.R, R/xbart.R, R/rbart.R,
  man/bart.Rd, man/bart2.Rd, man/rbart.Rd, man/dbartsControl.Rd, man/dbartsSampler-class.Rd, man/xbart.Rd -
  and the resolution is ORDER, not partition: this slice lands FIRST, and rd-records rebases onto its landed
  tip. That is the cheaper direction because this slice moves formals (section 8) and rd-records rewrites the
  prose around them; the reverse order would rewrite the same prose twice.
- One boundary is drawn by content rather than by order: rd-records' item 2 rewrites man/dbarts.Rd's
  `n.samples`, and THIS SLICE DOES NOT TOUCH man/dbarts.Rd at all. Section 12's whole-number sentence lands
  in man/dbartsControl.Rd (the shared home for the control counts) plus man/bart.Rd's legacy spellings, and
  nowhere else, so the two slices never edit the same `\item`.

## 16. Commit plan and gates

Four commits.

1. `surface-refusal-coverage`. R/generics.R (the `foreignArgsFor` helper, five reason tables, the two new
   class lists, 27 call sites with the two BELOW-the-trees-branch placements, the three conditional
   `forest`/`contribution` refusals, the three dead `...` forwards deleted, `plotTree.dbartsSampler`'s one
   line, `refuseSamplerExtractArgs`); [[R/bart.R:2547@9d0ee10f]], [[R/bart.R:2598@9d0ee10f]] (`survivalProbabilities`); man/bart.Rd,
   man/bart2.Rd, man/rbart.Rd, man/survivalProbabilities.Rd, man/plotTree.Rd; one NEWS UPGRADING item; the 57
   new assertions from section 13's first five bullets.
2. `surface-shape-refusals`. R/generics.R (sections 5 and 6: four `class`/`ci.level` checks, one
   `forest`/`ci.level` check); [[man/bart.Rd:218@9d0ee10f]], [[man/bart2.Rd:359@9d0ee10f]], [[man/bart2.Rd:387@9d0ee10f]], [[man/bart2.Rd:391@9d0ee10f]]; one NEWS item;
   test-multinomial-generics.R and test-ordinal.R rewritten blocks, test-predict-blend.R's new block.
3. `surface-argument-order`. R/generics.R (six `fitted` signatures, `fitted.bartHurdle`'s `sample` formal
   dropped and its two body uses, `residuals.bartHurdle`'s pass); [[R/plot.R:226-231@9d0ee10f]]; man/bart.Rd,
   man/bart2.Rd, man/rbart.Rd `\usage`; one NEWS item; test-generics-intervals.R and test-plot-generics.R
   additions. This is the codoc-sensitive commit: a formal and its `\usage` line move together.
4. `construction-whole-numbers`. [[R/utility.R:160-179@9d0ee10f]]; [[R/dbarts.R:1675@9d0ee10f]], [[R/dbarts.R:1686@9d0ee10f]]; man/dbartsControl.Rd and
   man/bart.Rd's legacy count items; [[man/dbartsSampler-class.Rd:284-291@9d0ee10f]]; one NEWS item;
   test-control-errors.R, test-sampler-residuals.R, test-bartcore.R additions.

Doc-freshness re-anchoring: per section 15 the anchor pass runs LAST and ONCE over the stacked tree, not per
commit - re-align every strict miss from the stacked `git diff -U0` line map by editing the docs/design
anchors in place so each file's line count is invariant. Anchor exposure from this slice alone: 21 anchors
into R/generics.R (15 past [[man/dbartsSampler-class.Rd:230@9d0ee10f]], so effectively all of them move on commit 1), 41 into R/dbarts.R (6 past
:1675), 34 into R/bart.R (2 past [[man/dbartsSampler-class.Rd:2547@9d0ee10f]]), 6 into R/A_class.R, 1 into R/utility.R, 0 into R/plot.R, plus 41 into
the man pages this slice edits (bart.Rd 6, bart2.Rd 9, rbart.Rd 3, survivalProbabilities.Rd 1, plotTree.Rd 1,
dbartsSampler-class.Rd 21). Take the live `Rscript tools/check-doc-freshness.R .` count as the baseline
before commit 1 rather than quoting an earlier slice's.

Equivalence, every commit: `benchmarks/R/equivalence.R compare` against
`benchmarks/baselines/equivalence-736bfb05.rds`, plus `benchmarks/R/bcf-equivalence.R` against
`bcf-equivalence-6e3b9fb8.rds` and `benchmarks/R/multinomial-equivalence.R` against
`multinomial-equivalence-4d9a3337.rds` (43 / 12 / 11 channels). ALL BITWISE IDENTICAL, zero re-records: no
file under src/ is touched, no R change reaches the sampler, and the one behavior on a draw path
(`predict(type = "ppd")`) is unchanged in every arm. If any channel moves, something in this slice reached
the sampler and the commit is wrong. The harness calls `predict` with named arguments only, so it needs no
edit for commit 3 - verify by reading before running.

`bench-sampler` is not required: no sampling path is touched and the only added work is a `setdiff` over at
most twelve names on generic calls that are not in any inner loop. Run
`benchmarks/R/bench-sampler.R compare` against the current baseline once at the end of the slice on a quiet
machine if the RC tip wants a clean sheet.

## 17. Residue

Observed, out of scope, recorded so it is not rediscovered.

- The terminal catch-all (section 3): landed 74e2e050 (the pre-review defect
  slice: "positional foreign arguments"), not post-1.0 residue - the
  integration tip this doc's landing was branched from.
- `coerceOrError`'s two EXISTING messages ([[R/utility.R:164@74e2e050]], [[R/utility.R:175@74e2e050]]) interpolate `mc[[2L]]` unguarded and have
  the same latent defect section 10 fixes for the new one - at [[R/spec.R:544-548@74e2e050]] both would render an
  expression as a quoted name. Left alone here: they are existing corpus, and error-style.md's slice L
  rewrites the corpus in one sweep.
- `extract(type = "trees")` forwards its whole call to the sampler's `getTrees`, which declares no `...`
  ([[R/dbarts.R:1999-2004@74e2e050]]), so an unknown name there raises R's own "unused argument" rather than a package
  message. Left: the two mechanisms disagree in voice, not in outcome, and closing it means rewriting the
  call-rewrite at [[R/dbarts.R:481-491@74e2e050]].
- `fitted(type = "ppd")` on bart is a Monte Carlo mean over a fresh draw, so it consumes RNG where
  `type = "ev"` does not. Unchanged here; noted because section 13's new `fitted` assertions must not sit
  above an RNG-locked block for that reason as well.
- `as_draws_array`/`as_draws_df.bartMultinomial`'s `vars` ([[R/diagnostics.R:226@74e2e050]], [[R/diagnostics.R:229@74e2e050]]) stays a formal and
  stays inert, documented rather than refused (VD). It is the one place on the S3 surface where an inert
  argument survives this slice, and it survives because the `posterior` generics fix the signature.
- `makeind(all = )` keeps its BayesTree signature (VD). `run(n.threads)` keeps its documented-reserved
  status (VD).
- Variance-forest scale-leaf staleness and the `growFromRoot` all-zero alignment: fix branch chosen, code
  deferred post-1.0 (VD). Neither is an R-surface fact and neither is touched here.
- `fitted.bartHurdle`'s `type` vocabulary includes `"prob"`, which no other `fitted` carries; that is a real
  channel, not drift, and stays.
- `plot.pdbart` ([[R/plot.R:476-482@74e2e050]]) is `(x, xind, plquants, cols, ...)` - a fourth shape on the plot surface,
  correct for an object whose panels are indexed. Not aligned, not touched, recorded so section 9's "six
  siblings" is not later read as "seven".
- No `plot` method gains a refusal: every one of them forwards `...` to `graphics::plot` deliberately.

## Landing note (2026-08-26)

LANDED at d48aef8a (design record 37ab6ea9), cherry-picked from its
gated worktree build; both gate batteries green on its own base
9d0ee10f: tinytest 7458/0 (85 new or rewritten assertions, split
57/9/8/11), equivalence trio bitwise 43/12/11, lintr::lint_package()
0 tree-wide, R CMD check --as-cran 1 NOTE, NEWS db 345, rc-codoc
green. Behavioral probes confirmed at gate: extract(type = "trees",
newdata =) still forwards; fitted(fit, "ev", 0.9) takes ci.level in
slot 3; the class + ci.level and residuals(ci.level =) refusals fire
with the designed wording. One probe correction for the record:
bart() carries BayesTree argument names and no dots, so the
whole-number refusal for a thread count fires at
dbartsControl()/dbarts() and at the legacy interface's own coerced
sites, not at a bart(n.threads =) spelling (which R itself refuses as
unused). Deviation from section 13: the design's dbarts(n.threads =
2.7) probe named a formal dbarts() does not have; the landed test
probes dbarts(seed = 2.7), the direct call site the section 10 census
names. Two review touch-ups landed with the slice: the
refuseSamplerExtractArgs helper sits above the comment block that
describes extract.dbartsSampler, and the NEWS entry claims the
plotTree catch-all for the dbartsSampler method alone (the fit-class
refusals predate this slice).
