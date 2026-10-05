# pdbart-on-bart: pdbart and pd2bart fit through bart, with type, newdata and the survival families

Status: PLANNED 2026-10-04. Ruled: dec-B203 to dec-B229 in
[decisions.md](../decisions.md); one agents' call (9, the default time's
fallback) pending the maintainer before slice 3.

agent: one implementer per slice (R only); an Opus reviewer per slice.
rng: posterior-changing at the R layer only, with the gate exception stated
under Constraints. Every pdbart and pd2bart call fits a different model or
reads a different fit path than before, so its draws change; there is no
engine change and `bart`'s own draws are untouched.
window: pre-release. Three slices in order; each leaves the package
releasable, and all three land before 1.0-0.
budget: ~3,300 lines (R ~1,300, tinytest ~1,550, man ~320, NEWS ~40,
records ~80, treatSens ~10). Plans have run 1.5-2x low.

## Goal

pdbart and pd2bart fit their model through `bart`, under `bart`'s names and
defaults, with BayesTree's spellings translated for one release. They average
any prediction a fit makes on the scale `type` chooses, over the fit's own rows,
a subsample of them, or `newdata`, optionally weighted; in a formula fit they
vary variables of the data. They serve every family whose prediction is one
value per row per draw, and on aft and hazard fits plot survival, the event
probability or the cumulative hazard at chosen times. The result says which
scale it holds, and the plots label it.

## Context

- Today [`pdbart`](../../R/partialDependence.R) and
  [`pd2bart`](../../R/partialDependence.R) redirect their call to
  [`bartBT`](../../R/bart.R) in
  [`pdbart.prologue`](../../R/partialDependence.R), dropping any name `bartBT`
  does not take, take a sampler from `samplerOnly`, and run it in
  [`pdbart.getAndInitializeSampler`](../../R/partialDependence.R). A fit
  passed in is predicted from its saved trees, or refit from its stored call
  through the same `samplerOnly` route. `xind` names model-matrix columns
  ([`pdbart.resolveXind`](../../R/partialDependence.R)); the default grid is
  [`pdbart.defaultLevs`](../../R/partialDependence.R), which fails on a
  missing value; the result is built by
  [`pdbart.buildResult`](../../R/partialDependence.R); a sampler's draws are
  reduced by [`pdbart.drawMeans`](../../R/partialDependence.R).
- [`bart`](../../R/bart.R) sets up the initial forest (`warm.start`,
  `n.grow.sweeps`, or a draw from the prior) after its `samplerOnly` return,
  and runs burn-in through [`runWithBurnIn`](../../R/bart.R) with the callback
  and its warning handling. A sampler taken from `samplerOnly` and run by
  pdbart therefore starts differently from `bart`'s own fit at the same seed.
- `bart`'s handling of 0.9-x calls lives in [`forwardToLegacyDoor`](../../R/tombstones.R),
  [`bartBTOnlyFormals`](../../R/tombstones.R),
  [`consolidatedArgsFor`](../../R/tombstones.R) (`power`, `base`, `sigdf`,
  `sigquant` and the rest, which `bart` accepts with its own warning) and
  [`noteFrontDoorDefaults`](../../R/tombstones.R) with its key
  [`frontDoorDefaultsKey`](../../R/tombstones.R) in `onceWarnState`. The
  registry every transition notice joins is
  [`dbartsTombstones`](../../R/tombstones.R), checked by
  [test-tombstones.R](../../inst/tinytest/test-tombstones.R): an
  `"argument"` row's successor must be a formal of its owner, so pdbart's
  translations are `"behaviour"` rows, as `bart`'s own forwarding is.
- Prediction on new rows: [`predict.bart`](../../R/generics.R),
  [`predict.bartNegbin`](../../R/generics.R),
  [`predict.bartHurdle`](../../R/generics.R), with
  [`validateType`](../../R/generics.R) (which folds `"link"` into
  `"bart"`), [`predictTermOffset`](../../R/generics.R) and
  [`preparePredictRows`](../../R/data.R) coding rows and offsets as
  `predict.lm` does. Coded rows replay through
  [`predictCodedTest`](../../R/dbarts.R).
- Survival: [`survivalProbabilities.bart`](../../R/bart.R) and
  [`hazardSurvivalProbabilities`](../../R/bart.R). The latter expands every
  subject to every period and holds draws x subjects x periods at once, which
  does not run at ordinary sizes on a fine period grid; pdbart replays in
  chunks instead (dec-B209, dec-B220). An aft fit stores log times in `fit$y`
  and the event indicator in `fit$status`; a hazard fit stores its period
  grid in `fit$periods` and no original times.
- Response classification before a fit: [`autoRawResponse`](../../R/bart.R),
  [`detectAutoCounts`](../../R/bart.R),
  [`detectAutoMultinomial`](../../R/bart.R) and
  [`detectAutoOrdinal`](../../R/bart.R). A `Surv` response under `family =
  "auto"` resolves to aft, which none of the detectors reports.
- Plots: [`plot.pdbart`](../../R/plot.R) and [`plot.pd2bart`](../../R/plot.R).
  `plot.pdbart` passes `type = "n"`, `xlab` and `ylab` to its first `plot`
  call beside `...`, so a user's `type`, `xlab` or `ylab` collides with them
  (0.9-34 too); `plot.pd2bart` hands `...` to `image`.
- Tests today: [test-pdbart.R](../../inst/tinytest/test-pdbart.R) (every call
  in BayesTree names), [test-pdbart-keeptrees.R](../../inst/tinytest/test-pdbart-keeptrees.R),
  the burn-in sigma pin in [test-input-guards.R](../../inst/tinytest/test-input-guards.R),
  and `pdbart.drawMeans` in
  [test-packaging-copies.R](../../inst/tinytest/test-packaging-copies.R),
  which stays: the sampler route still uses it.

## Decision

Ruled, and not restated here:

- Fit through `bart`: dec-B203. Data arguments named `formula` and `data`:
  dec-B204. BayesTree spellings translated for one release: dec-B205.
  `bart`'s defaults and the chain margin under `combineChains = FALSE`:
  dec-B206.
- `type` and `newdata`: dec-B207 (amended for the link default and the hurdle
  default). Variables in a formula fit: dec-B208. Families served, survival on
  aft and hazard, multinomial and ordinal refused: dec-B209. Transform each
  row, then average: dec-B210. The offset included: dec-B211. Arguments
  pdbart sets itself, and `...` until 1.1-0: dec-B212. The once-per-session
  message: dec-B213. Averaging weights: dec-B214. The scale recorded and
  labelled: dec-B215. The grid from `newdata`: dec-B216. Offsets on
  `newdata`: dec-B217. `n.average.rows`: dec-B218. A formula fit's rows from
  its stored call: dec-B219. The hazard size check and `n.max.predictions`:
  dec-B220. The default time: dec-B221. The survival result's shape:
  dec-B222. The survival scales: dec-B223. The plot views and `plot.type`:
  dec-B224. Subject rows on hazard: dec-B225. Aft's survival default:
  dec-B226. `type = "auto"`: dec-B227. Weights under a subsample: dec-B228.
  The translation warning for package callers: dec-B229.

Calls made by the agents in drafting. Those marked (*) change what a user
sees and go on the register as agent-made decisions (step 13). Those marked
(pending the maintainer) depart from, or cut against, a ruling's wording:
each is put to the maintainer, one per message, before the slice that would
implement it, and that behaviour is not built until ruled.

1. A data call fits `bart(..., keepTrees = TRUE, keepSampler = TRUE)` and
   predicts each grid value from the saved trees, replacing the `samplerOnly`
   route; a fit kept without its sampler is refit by re-evaluating its stored
   call through the function that made it, with trees kept. So `pdbart(x, y,
   seed = s)`, `pdbart(bart(x, y, seed = s, keepTrees = TRUE))` and the refit
   of `bart(x, y, seed = s)` agree exactly; a fit seeded only by `set.seed`
   is not reproduced by a refit. Measured at 1,000 rows, five predictors and
   55 grid values, about a tenth of the peak memory of stacking test rows.
2. The `bart` call is pdbart's own matched call rewritten - the function
   replaced, names translated, pdbart's own arguments removed - and
   evaluated in the caller's frame, never forwarded through `...`. So an
   argument the caller wrote unevaluated, such as treatSens' `k = chi(df,
   scale)` beside a local `chi` stub, reaches `bart` as written and resolves
   in `bart`'s prior vocabulary.
3. (*) The translation table covers every name `bartBT` takes and `bart`
   does not - `bartBTOnlyFormals` together with `power`, `base`, `sigdf` and
   `sigquant` - less `x.test` and `sampleronly`, which are refused; a test
   pins it to that set. pdbart warns once per name itself and holds back
   `bart`'s own consolidated-name warning, leaving `bart`'s once-per-session
   keys as they were. `power`, `base` and `splitprobs` become one
   `tree.prior = cgm(...)`, refused beside a caller's `tree.prior`;
   `proposalprobs` is set on a copy of the caller's `control` when one is
   given, otherwise `control = dbartsControl(proposal.probs = )`; `sigdf` and
   `sigquant` reach `bart` under its consolidated names, since the residual
   prior rides a family pdbart cannot build before the response is known.
   The translation warning fires for package callers too, as `bart`'s
   forwarding warning does (dec-B229); treatSens changes in lockstep.
4. In a formula fit a column number in `xind` is refused, the grid is taken
   on the variable's own scale, and an offset expression built from the
   varied variable is rebuilt with it.
5. (*) A `dbartsSampler` passed in takes only `type = "bart"`; it carries no
   fit object to transform through. A hazard sampler, recognised by
   `attr(sampler$control, "bartcore.hazard.periods")` (its model's family
   reads `"probit"`), stays refused after slice 3, since its rows are
   person-period rows. An aft sampler is served only with an explicit `type
   = "bart"` and refused under `"auto"`, which on an aft fit means survival.
6. (*) pd2bart's two-predictor shortcut (each grid point predicted as one
   row) applies only when every averaged row is identical once the two
   variables are set - a fit with exactly two predictor variables and no
   offset or a constant one - and `type` is not `"ppd"`, whose rows each
   draw their own noise, so their mean has a narrower spread than one row's
   draw. Then `newdata`, `n.average.rows` and
   `average.weights` have no effect and a warning names whichever was given;
   otherwise the general route runs.
7. (*) `n.average.rows` is sampled once, without replacement, before the grid
   loop, from the rows with a positive averaging weight (all rows when no
   `average.weights`), so every grid value averages the same rows; on
   survival fits it counts subjects. `average.weights` must be finite,
   non-negative and not all zero. Its length is the number of the fit's rows
   before the 0-weight exclusion (one per `newdata` row with `newdata`); an
   entry for a row the fit gives weight 0 is ignored.
8. The hazard size check counts in double precision, subjects after
   `n.average.rows` and the 0-weight exclusion. The replay runs in chunks of
   whole subjects, sized so a chunk's rows x draws stay under about 5e6
   doubles (40 MB per array); the bound is an argument of an unexported
   helper, which the tests call with a small value, not a user setting.
9. The Kaplan-Meier default time is computed inside the package (`survival`
   is only suggested): on an aft fit from `exp(fit$y)` and `fit$status`; on a
   hazard fit from its person-period rows, on its period grid, so with
   `hazard(breaks = )` it is the Kaplan-Meier median of the coarsened times.
   (pending the maintainer, before slice 3) The condition for falling back
   to the median follow-up: dec-B221 says "when fewer than half the subjects
   have the event"; the plan's draft was "when the Kaplan-Meier curve does
   not reach 0.5", the condition under which no median exists. Under
   censoring the two differ in both directions. (*) The median follow-up is
   the median of all subjects' observed times, events and censored alike.
10. (*) `plot.type = "curves"` is refused when fewer than three times were
    computed, and on a result with no times margin.
11. (*) `plot.pdbart` takes `type`, `xlab` and `ylab` out of `...` for its
    own lines and axes (a bug fix). `plot.pd2bart` passes `...` to `image`
    like any other argument, with no handling of `type` of its own
    (dec-B224 as corrected).
12. (*) `times` is refused by name on a family without survival, and on aft with
    a non-survival `type`. On a hazard fit the `period` column is left out of
    the default `xind`, and `xind = "period"` is refused by name.
13. (*) On a hurdle fit the averaged rows are all observations (the zero part's
    rows); the result carries the two component samplers as `fit`, dropped
    under `keepSampler = FALSE`, and no copied draws.

## Constraints

- No engine or bridge change. The compiled hazard routine is TODO
  pdbart-hazard-compiled; individual curves pdbart-ice; per-category partial
  dependence pdbart-per-category; replacing `...` with formals
  pdbart-closed-arguments. All after 1.0-0.
- Each slice is releasable: a family a later slice serves is refused by name,
  before anything is fitted, until that slice lands - from an explicit
  `family`, from the response under `"auto"` (a `Surv` response is aft), or
  from a fit or sampler passed in. Samplers follow agents' call 5 after
  slice 3. Slice 1 refuses negative binomial, hurdle,
  aft and hazard; slice 2 lifts negative binomial and hurdle; slice 3 lifts
  aft and hazard. Multinomial and ordinal stay refused.
- The default rows are the fit's own (its `subset`, less rows dropped for a
  missing response or under `na.action`, less rows it gives a 0 weight),
  found by its row names; never every row of `data`.
- Each averaged row's offset is applied exactly once: the fit's offset
  expression or `offset()` term evaluated on the (varied) row, or, for a
  plain-vector training offset, that row's stored value.
- Every value is computed per row (per subject on survival fits) and averaged
  last.
- `fd` stays a draws x levels matrix (pd2bart: draws x grid points) under the
  default `combineChains = TRUE`; treatSens and glossa read it unchanged.
- The `samplerOnly` route stays only for a `dbartsSampler` passed without
  `keepTrees`.
- Exception to the RNG-class gates. The plan README's posterior-changing
  class asks for regenerated snapshots, a re-recorded equivalence baseline,
  the statistical compare and the exact-posterior gates. None applies here:
  no seed-locked snapshot, equivalence scenario, benchmark or exact gate
  calls pdbart or pd2bart (checked by search of `inst/tinytest`,
  `benchmarks` and `.github/workflows`), and the plan changes nothing those
  gates exercise - `bart`, `bartBT` and the engine are untouched. The design
  note is kept. A slice that turns out to touch `bart`, `bartBT`, the bridge
  or the engine loses this exception and runs the full class.
- Out of scope: `predict`'s own handling of a caller's `offset` beside the
  fit's.

## Steps

Slice 1 - fit through `bart`.

1. Formals `formula, data, xind = NULL, levs = NULL, levquants = c(0.05,
   seq(0.1, 0.9, 0.1), 0.95), pl = TRUE, plquants = c(0.05, 0.95), ...` on
   both functions. In the prologue, per agents' calls 2 and 3: translate
   BayesTree names in the matched call; refuse a setting given in both
   spellings, naming both; a named `x.train` holding a fit or sampler warns
   to pass it first, unnamed. Refuse by name `keepTrees = FALSE`,
   `samplerOnly`, `test`, `offset.test` and their spellings. Refuse before
   fitting multinomial and ordinal (pointing at `predict`) and, until their
   slices, negative binomial, hurdle, aft and hazard. Fit through `bart` with
   trees and sampler kept; take the fit path. Refit a fit without its
   sampler through its own function. `pdbart.defaultLevs` drops missing
   values before taking quantiles or counting unique values.
2. Migrate [test-pdbart.R](../../inst/tinytest/test-pdbart.R) to `bart`'s
   names and add: equality at one seed of a data call, a fit with trees and
   the refit of a fit without (with `n.grow.sweeps` too); a sampler passed in
   checked against the row mean of its own `$predict` at a grid point; each
   translated name asserting the setting itself (the control's
   `n.threads`, `printEvery`, `n.cuts`, `useQuantiles`, `keepTrainingFits`,
   the result's `bartcall`, the presence or absence of `fit` and
   `yhat.train`), not only equal draws; the table pinned to `formals(bartBT)`
   less `formals(bart)`, less the two refused names; `power` with
   `tree.prior` refused, `proposalprobs` merged into a given `control`;
   warnings counted with `withCallingHandlers`; a call with a local `chi`
   stub and `k = chi(1.25, Inf)` whose fit's leaf prior records that
   hyperprior; `n.trees` and `family` honoured; a misspelled name refused;
   every family refusal before any fit (an explicit family, a three-level
   factor, an ordered factor, a count matrix, a `Surv` response under
   `"auto"`, and a fit and a sampler of each refused family passed in);
   the translation warning firing for a call from package code (dec-B229); a training column
   with a missing value giving a grid.
3. Result: `n.chains` added; `bartcall` the fit's call; under
   `combineChains = FALSE` `fd` gains a leading chain margin (pdbart each list
   entry, pd2bart its matrix), merged otherwise; `keepSampler = FALSE` drops
   `fit`; each row's offset added before averaging; rows with a 0 fit weight
   left out. pd2bart's shortcut limited per agents' call 6. Plot methods
   merge a chain margin before taking quantiles, and handle `type`, `xlab`
   and `ylab` per agents' call 11. Tinytest: `n.chains`; shapes under both
   `combineChains`, the split array equal to the merged matrix after merging;
   both plot methods on both shapes into a null device; `plot(pd, type =
   "l", xlab = "a", ylab = "b")` without error; `plot.pd2bart` handing a
   `type` in `...` to `image` unchanged (dec-B224); `keepSampler = FALSE`; an offset fit and a `binaryOffset` call
   shifted by their offset; a probit fit with a 0/1 row mask averaging its
   active rows only; a two-variable fit with a varying offset taking the
   general route and equal to the row mean of `predict`; a two-variable fit
   with `type = "ppd"` taking the general route (its spread across draws
   matching the general route's, not one row's).
4. pdbart's own once-per-session message for a data call with no BayesTree
   name, exempt for package callers; `bart`'s message held back inside pdbart
   with its key saved and restored. `"behaviour"` registry rows: a
   BayesTree-spelled pdbart call, the same for pd2bart, and the defaults
   message, each owned by its function, expiring 1.1-0. Tinytest, clearing
   and restoring both keys in `onceWarnState` as test-tombstones does: the
   message once; `bart(x, y)` afterwards still shows its own.
5. Retarget [test-pdbart-keeptrees.R](../../inst/tinytest/test-pdbart-keeptrees.R)
   to a translated `keeptrees = TRUE, nskip = 0` call; the burn-in sigma pin
   compares against a `bart` fit. Help page, NEWS for this slice, saying
   which families are refused until later; `docs/design/pdbart-on-bart.md`
   started with the fit route, the averaged rows and the translation, Status
   PLANNED, and grown by each later slice.

Slice 2 - `type`, `newdata`, variables, subsamples and weights.

6. `type`, formal default `"auto"`, resolved once the fit is known; the
   result records the name `validateType` returns (`"bart"` for `"link"`),
   which the plot labels key on. Each grid value's rows are predicted
   through the fit's own `predict` method with that type and averaged per
   draw; `"forest"` refused; `"sigma"` only on a heteroscedastic fit; values
   a family does not take refused by name. Negative binomial and hurdle
   lifted (agents' call 13 for hurdle). Tinytest: `"auto"` per family
   (`"bart"` on gaussian, student, probit, logistic, negative binomial;
   `"ev"` on hurdle); each value against the row mean of `predict(fit,
   newdata, type = )` at a grid point, binary `"ev"` included (which pins
   transform-then-average); `"prob"` on hurdle; `"sigma"` on a
   heteroscedastic fit and refused elsewhere.
7. Variables in a formula fit. The raw rows, for a data call and a fit
   passed in alike, are the formula's variables as `get_all_vars` collects
   them from the call's `data` (re-evaluated from the stored call for a fit
   passed in), cut to the fit's rows by name; refused with a request for
   `newdata` when that fails. For each grid value the variable is set in the
   raw rows and the rows are coded through the stored formula. Tinytest:
   `xind = "a"` on `y ~ log(a) + poly(b, 2) + a:b` equal to the row mean of
   `predict` with `a` set; the grid on `a`'s scale; a factor variable by
   level; a column number refused; an offset expression in `a` rebuilt; a fit
   passed in reading its stored data, and refused with `keepCall = FALSE`; a
   response with 5 of 50 values missing averaging 45 rows.
8. `newdata` coded as `predict` codes it, the default grid from its rows
   (missing values dropped), a factor showing every level the fit knows;
   offsets as `predict` evaluates them, a plain-vector offset refused on
   other rows. `n.average.rows` refused with `newdata`. `average.weights`
   per agents' call 7, normalized; under `n.average.rows` given for all the
   fit's rows and renormalized over the sample (dec-B228); each row's offset
   applied once (Constraints), checked against `predictTermOffset`'s rule
   for a plain vector so that a subsample's stored offsets are neither
   refused nor added twice. Tinytest: a subgroup equal to the training run
   restricted to it; an unseen level refused; missing predictors kept; the
   subsample reproducible under `set.seed` and equal to `newdata` of the
   same rows; weights against a hand-weighted mean, a rescaled vector giving
   the same `fd`, each refusal by name, an entry for a 0-weight row ignored;
   weights with a subsample equal to the renormalized weighted mean over the
   sampled rows; an offset fit (an expression and a plain vector) under a
   subsample equal to the hand computation.
9. Plot labels by scale and family; a user's `ylab` or `main` still wins.
   Tinytest: the label per type. Help page, NEWS and the design note for
   this slice.

Slice 3 - survival on aft and hazard fits.

10. `times`, `NULL` meaning the default time of agents' call 9; refused per
    agents' call 12. `type` on aft and hazard: `"survival"` (the `"auto"`
    default), `"event"`, `"cumhaz"`, each per subject then averaged; aft also
    takes `"bart"` and its other `predict` values, which give no times
    margin. A survival result's `fd` is draws x times x levels (pd2bart
    draws x times x grid points), chains leading when split, the times margin
    kept at length 1, its dimnames the times. The hazard `period` column per
    agents' call 12. Aft and hazard lifted.
11. Hazard: each grid value set in every subject row (training subjects, or
    `newdata` and `n.average.rows` as subject rows), latents replayed
    through `predictCodedTest` for periods up to the largest time only, in
    chunks per agents' call 8, each chunk turned into per-subject survival by
    the cumulative product, transformed by `type`, and summed per draw.
    Before any replay, the size check: subjects x periods replayed x draws x
    grid values, in double precision, against `n.max.predictions` (default
    5e9), the refusal giving the count, the limit and each way out. Aft:
    `survivalProbabilities` per grid value in chunks of subjects. Tinytest:
    hazard and aft `fd` at a grid value equal to the subject mean of
    `survivalProbabilities` with the variable set, across chunk boundaries
    (the helper called with a small bound); `"event"` and `"cumhaz"` against
    the subject means of 1 - S and -log S, the latter differing from -log of
    the mean; the default time against `survival::survfit` when installed -
    on aft's times, on a hazard fit's default grid, and under `hazard(breaks
    = )` on the coarsened times - and its fallback on a fit whose curve stays
    above 0.5; the size check refusing with its levers, not firing on a
    20-period fit, and counting a product above 2^31 without overflow;
    `n.max.predictions` raising it and refused when not a positive number;
    `xind = "period"` refused and left out of the default; shapes under both
    `combineChains`; a hazard sampler refused after this slice, and an aft
    sampler served with `type = "bart"` and refused under `"auto"`; the
    fallback condition as the maintainer rules on agents' call 9.
12. Plots: `plot.type = c("dependence", "curves")` on both methods; the
    default puts the variable on the axis with one line per time, the time
    named; `"curves"` puts time on the axis with one curve per grid value;
    pd2bart one image per time. Tinytest: both views into a null device, the
    refusals of agents' call 10, the labels. Help page (the type table; the
    pictures named: a partial dependence plot, adjusted curves or direct
    adjustment; `type` chooses what is computed and `plot.type` how it is
    drawn), NEWS and the design note for this slice.

Records, with slice 1 and at the end.

13. Before slice 3 is built: agents' call 9 put to the maintainer, its
    ruling recorded in the register and this plan updated. Before slice 1
    lands: the agents' calls marked (*) - 3's translation table, 5, 6, 7,
    9's median follow-up, 10, 11's `plot.pdbart` half, 12 and 13 - entered
    on the register as agent-made decisions with a blank "Marked:" line. At each
    slice's landing: this plan's Landing note, the INDEX row's status, and
    the design note's Status line. After slice 3: "Marked:" filled on dec-B203
    to dec-B228 as the maintainer marks them, and the TODO item pdbart-on-bart
    removed.

Lockstep, with slice 1: treatSens, on its `dbarts-1.0` branch, in
`treatSensBART.R`, changes its two `pdbart` calls from `ntree = , nskip = ,
ndpost = ` to `n.trees = , n.burn = , n.samples = ` and adds `n.chains = 1L`,
which keeps its run length and the shapes of `fd` and `yhat.train` it reads.
Its `k`, a number or `chi()` through a local stub, reaches `bart` by agents'
call 2, which the slice-1 test pins. About ten lines and treatSens' own check.

Note for the maintainer, no code: glossa (CRAN) calls pdbart on a probit fit
and applies the normal CDF to `fd` itself. The link-scale default keeps that
correct. If its users pass `binaryOffset` through its `...`, its plots move
by the offset under dec-B211; worth a line to its maintainer before the
release.

## Help pages and NEWS

- `man/pdbart.Rd` rewritten over the three slices: usage with the full
  signature; arguments (`formula`, `data`, `...` as `bart`'s, `type` with the
  per-family table, `times`, `newdata`, `n.average.rows`, `average.weights`,
  `n.max.predictions`, `plot.type`); value (`fd`'s shapes, `n.chains`,
  `type`, the times, components from the fit); details (the averaged rows,
  offsets, transform-then-average, thinning is `n.thin`, what each slice
  still refuses); examples in `bart`'s names with `n.chains` and `n.threads`
  small and the survival example inside `\donttest`. `man/dbarts-deprecated.Rd`:
  the behaviour rows.
- `inst/NEWS.Rd`, against 0.9-34 only, pruned:
  - UPGRADING: pdbart and pd2bart fit through `bart`, with its names and
    defaults (75 trees, four chains merged, so `fd` has 2000 rows by
    default; a factor is one predictor); the data arguments are `formula`
    and `data`, positional calls unchanged; a BayesTree name (`ntree`,
    `keepevery`, `x.train`, ...) is translated with a warning until 1.1-0,
    and an argument neither takes is an error where it was ignored. The
    0.9-34 model is `pdbart(bartBT(x, y, keeptrees = TRUE))`, its results
    still subject to the changes below. In a formula fit `xind` names
    variables of the data, where it named model-matrix columns.
  - USER-VISIBLE CHANGES: `fd` includes the fit's offset, as BayesTree's
    included `binaryOffset`; rows a fit gives a 0 weight are left out of the
    average; the result carries `n.chains` and its `type`, and
    `combineChains = FALSE` gives `fd` a chain margin.
  - NEW FEATURES: `type` (`"auto"`: the link scale, 0.9-34's, except the
    mean response on hurdle fits and survival at the median survival time on
    aft and hazard fits, with `"event"`, `"cumhaz"` and `times`); `newdata`,
    `n.average.rows`, `average.weights`, `n.max.predictions`; `plot.type`.
  - BUG FIXES: `plot` on a pdbart result accepts `type`, `xlab` and `ylab`,
    which collided with its own; a training column with a missing value no
    longer stops the default grid.
- The existing factor bullet stands.

## Verification

Per slice, against the slice's own library:

```sh
R CMD INSTALL -l "$LIB" .
R_LIBS="$LIB" Rscript -e 'tinytest::run_test_file("inst/tinytest/test-pdbart.R")'
R_LIBS="$LIB" Rscript -e 'tinytest::run_test_file("inst/tinytest/test-pdbart-keeptrees.R")'
R_LIBS="$LIB" Rscript -e 'tinytest::run_test_file("inst/tinytest/test-tombstones.R")'
R_LIBS="$LIB" Rscript -e 'tinytest::test_package("dbarts")'
R_LIBS="$LIB" Rscript -e 'lintr::lint_package()' && air format --check . && \
  R_LIBS="$LIB" Rscript tools/check-rc-codoc.R . && \
  R_LIBS="$LIB" Rscript tools/check-win-drift.R . && \
  Rscript tools/check-doc-freshness.R .
Rscript -e 'stopifnot(!is.null(tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")))'
```

and `R CMD check --as-cran` on a tarball built from a clean copy outside the
tree. Expected: every file passes; NEWS parses with its entry count up by
the new bullets; each slice's examples run in seconds. Per the gate
exception under Constraints, no snapshot, baseline or exact gate is
re-recorded or rerun; the reviewer confirms by search that the slice's diff
touches none of `bart`, `bartBT`, the bridge or the engine. treatSens: its
own test suite and `R CMD check` against the slice-1 library.
