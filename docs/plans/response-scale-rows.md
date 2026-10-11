# response-scale-rows: what the response gives a prior is taken from the rows in the likelihood at creation

Status: PLANNED. Revised from its critique and to dec-B405 on 2026-10-09; rechecked on 2026-10-11
against bartcore a2ef56a9, where every slice queued ahead of it has landed
([Recheck against bartcore](#recheck-against-bartcore)). Implements dec-B296, dec-B302,
dec-B309, dec-B314 (sigest held), dec-B315 (a hazard fit's periods), dec-B362 (several forests keep
refusing), dec-B363 (slice B's half: no documented by-hand call yet), dec-B364 (for `setResponse` and
`setOffset`), dec-B405 (the keep is a status: the engine returns it, the R interface warns, the flat
entries return it), dec-B424 (the hash alone moves, the C API's version stays 1/0), dec-B381 (the
planning calls approved) and dec-B386 (a constant response's window where every row is in); rests on
dec-B357 (a survival fit's range from observed times), dec-B418 (every install goes in as stored),
dec-B429, dec-B433 and dec-B439 (the starting sigma as landed), dec-B437 (stan4bart fixes its own range)
and, for slice C, dec-B327; revises dec-A140. No open call. This file plans slice B, creation. Slice A is
[aft-reanchor-observed-times.md](aft-reanchor-observed-times.md), landed. Slices C and D are stated at
the end and have no steps yet. dec-B303 stays in [forest-multiplier.md](forest-multiplier.md).

agent: one push. Opus implementer for the engine, the keep's status through the facade, the bridge and
the flat entries, dbarts.h, and their tests; sonnet for the R code (the hazard grid included), the
tinytest file, the gate scripts, the help once the engine is fixed, and stan4bart's part (a few lines in
its range install and their test, on its bartcore branch); opus reviewer told to refute, and to read
the diff three times: once for bits that must not move where no row is out, once for which rows each
number reads, once for the ABI (the header, the token, the C API test). Tier: changes draws
([Process by risk](README.md#process-by-risk)). Its blind critique ran at 8ad1fa58 and a second, scoped
to the keep's status, on 2026-10-09; both are answered in the rule, the constraints, the keep's status
and the calls below. The recheck changed no rule.
rng: by call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines
them. A row is IN when its weight is positive and no mask switches it off, OUT otherwise; a hazard
subject is in when its weight is positive.
- NEUTRAL, bit for bit: every sampler created with every row in, whatever it does afterwards (a mask or
  zero weights installed after creation included); every probit, logistic, ordinal and multinomial
  sampler; a hazard fit with every subject in, with stated boundaries, or on the default grid whose
  subjects at weight 0 hold only times some subject in holds; a single-forest gaussian or Student-t
  sampler created with rows out whose smallest and largest response less offset are both on rows in;
  a stan4bart fit, which installs its own range over creation's (dec-B437; item 9 of
  [The keep's status (dec-B405)](#the-keeps-status-dec-b405)). The keep (a re-derivation in rule 4's
  case) is the `updateScale = FALSE` call, literally: the sampler makes that call on every chain. Its
  status, its warning and the flat entries' new return value move no draw.
- POSTERIOR-CHANGING, a default changes, for a sampler created with rows out: a gaussian or Student-t fit
  with an extreme on a row out (the range); a gaussian model of several forests with any row out (the
  unit, so every forest's default sd, and what a stated `sd` means until
  [forest-sd-unit.md](forest-sd-unit.md) lands); a count fit under 0/1 weights (the center); `xbart` with
  a weight of 0 (each fold's range); any of them where the rows in leave sigest's regression no residual
  degrees of freedom (the sd that stands in) or take its floor; a hazard fit with a subject at weight 0
  whose time no subject in holds, or under `breaks = K` with any subject at weight 0 (the periods, and
  the rows of the subjects in). And `setResponse` or `setOffset` with `updateScale = TRUE` while rows
  are out, in every family that has a range or a center, through the R methods and the two flat C
  entries: the rows in, or, where they hold fewer than two distinct values, the transform kept in force
  where the tip reads every row. A deprecated `rbart_vi` fit with a zero weight is in this class: it
  calls `setOffset(..., isWarmup)` every warm-up sweep.
- Moves one recorded baseline scenario, `zeroweights` of the main corpus; nothing else recorded moves
  ([Context](#context)).
window: before the merge to main: it changes draws against 0.9-34. Everything dec-A191 queued ahead of it
has landed, and so has sparse-starting-sigma. After stan4bart-creation-mapping has landed on stan4bart's
bartcore branch (item 9). Before [forest-sd-unit.md](forest-sd-unit.md), so the unit's tests are written
once against the rows that stay. Serial with any other work in model.hpp, the response setters of
[`Chain`](../../src/bartcore/chain.hpp) and [`Sampler`](../../src/bartcore/sampler.hpp), the K-forest
constructor of [`Chain`](../../src/bartcore/chain.hpp), [`resolveSamplerSpec`](../../R/spec.R),
[`estimateSigmaFromLinearModel`](../../R/utility.R), [`expandDiscreteTimeHazard`](../../R/dbarts.R) or
dbarts.h.
budget: ~2065 lines planned: engine ~175; facade, bridge, flat entries, dbarts.h and the token ~75; R
~170; tests/cpp ~335; tinytest ~580; the C API test ~75; benchmarks and the workflow ~265; help, header
comments and NEWS ~125; design note, indexes, amended notes and TODO ~265. Forecast 3700 to 4750; stop at
~4750, or at ~630 of package code ([Budget and stops](#budget-and-stops)). stan4bart's part, ~25 lines,
is outside it.

## Goal

A fit created with rows out of the likelihood has the leaf prior, the prior center, the default forest
sds and, for a hazard fit, the periods of the fit on the remaining rows: the response of a row out
reaches no number and no draw at the rows in. Those numbers are then held, as they are today. A fit
created with every row in is the fit it is today, bit for bit. A re-derivation that finds too little to
read keeps the scale in force, and says so: as a typed status from the engine, as one warning per
sampler from the R interface, and as a distinct return value from the flat C entries.

## Context

Marks: (ran) rerun on 2026-10-11 on a shipped build of bartcore a2ef56a9, "the tip" below; (stands)
measured on the design's prototype library, on 0.9-34 or by a critic at 8ad1fa58, and not rerun; (read)
in the source at a2ef56a9. One fixture unless a line says otherwise: 200 rows, two predictors, rows 101
to 200 at weight 0 with their response 50 higher.

What reads which rows today.
- (ran) The range ([`GaussianResponse::rescale`](../../src/bartcore/model.hpp)): every row. Recorded
  [-2.165, 56.371] where the rows in span [-2.165, 6.029] and `subset` gives the latter; `k.scale` 29.27
  against 4.10 on the prototype (stands). Student-t and aft hold a Gaussian response and inherit it
  (read).
- (ran) The count center ([`NBResponse::computeShift`](../../src/bartcore/model.hpp)): every row. Half
  the rows structural zeros under 0/1 weights: recorded 0.2964, the rows in give 0.9895, the prototype
  0.9895 (stands). The 0/1 weights become the mask and the weights slot is empty.
- (ran) The unit of a forest's sd in a model of several forests
  ([`scaledResponseSd`](../../src/bartcore/chain.hpp)): every row. 0.4313 of the range, 25.25 in y,
  with or without the weights; (stands) 1.466 in y under `subset` and on the prototype, and held across
  a later `setWeights` on both builds.
- (ran) sigest is already the weighted regression over the rows of positive weight, by each routine:
  [`residualStandardError`](../../R/utility.R) for a matrix or a data frame and
  [`sparseResidualStandardError`](../../R/utility.R) for a caller's sparse x both give 1.0741, `lm`'s.
  What stands in for it is not. Where the rows in leave no residual degrees of freedom (three rows in,
  two predictors) sigest is 25.245, the sd of every row, where the rows in give 0.890
  ([`floorSigmaEstimate`](../../R/utility.R)); 40 rows, 20 at weight 0 with their response 100 higher,
  25 predictors: 50.84 for 1.12. The floor is not either: four rows in at one value give 8.4e-07, from
  the largest residual of every row ([`floorMarginalSigma`](../../R/utility.R)). No warning says so
  since sparse-starting-sigma (dec-B429): a `dbartsSigmaFallback` condition is signaled, `bart` records
  `sigest.fallback`, and a verbose fit prints a line that already counts the rows of positive weight
  ("in 3 rows"). No estimate is made under a fixed residual sd (dec-B439) or beside a supplied number.
- (ran) Held already: `setWeights` and `setActiveRows` move nothing. `setResponse` and `setOffset` with
  `updateScale = TRUE` under zero weights or a mask read every row ([-2.165, 56.371], with three rows in
  at one value as well). On several forests `updateScale = TRUE` is refused ("every forest keeps its
  leaf calibration stated against the scale fixed at creation").
- (ran) No row in. Every weight 0: the sampler is created in silence, the range every row's, sigest
  25.245 (the sd of every row; `lm` without the weights gives 25.006), the verbose line reading "in 0
  rows". A count fit under a mask of zeros is created in silence.
- (ran) One row in. The tip reads every row, in silence. (stands) The prototype's window of width 1
  gave, once every row was back in, rmse 2.41 and sigma 2.66 where the fallback gives 0.46 and 0.90.
- (ran) Fractional weights with a zero among them raise the data object's validity warning ("rows with
  'weights' of 0 are left out of the likelihood", [`dbartsData`](../../R/A_class.R)) at creation, and
  nothing else: one row in at weight 0.5, and half the rows out, each give that one warning. 0/1 weights
  raise none.
- (ran) `xbart` with a weight of 0: the response of the rows out reaches the result (1.38, 1.23 against
  1.15, 1.20 with those responses ordinary); (stands) on the prototype the two are identical, the engine
  edit alone closing it. `xbart` never calls `resolveSamplerSpec`: it calls
  [`enforceWeightPolicy`](../../R/spec.R) and [`estimateStartingSigma`](../../R/spec.R) itself and fixes
  the routine every fold runs (`sigmaPerFold`, [`xbart`](../../R/xbart.R): "dense", "sparse", "marginal"
  where all rows leave no residual degrees of freedom, or "function") by reading the
  `dbartsSigmaFallback` condition; each fold's estimate, its stand-in and its floor are made in
  [`xbartRunUnits`](../../R/xbart.R) (read). With every weight 0, and with one row in, it stops in its
  loss ("loss function must return non-missing numeric values"), with no warning; with 60 rows in at one
  value it runs, silent.
- (ran) A hazard fit's periods: 40 subjects on times 1 to 5, two at weight 0 holding 3.5 and 20: 7
  periods and 167 person-period rows, where the fit without the two has 5 and 134 and the two moved onto
  the grid (3 and 5) give 5 and 142. `breaks = 3`: periods 3, 4, 20 against 3, 4, 5 without the two.
  Every weight 0: created in silence on every subject's 7 periods. `hazard.logistic` refuses a zero
  weight, so only the probit link reaches the case (dec-B315 says the same). A masked probit row's
  response reaches no draw at the rows in: flipping it leaves the rows in `identical()` (3.3e-15 at the
  rows out). The periods are taken in [`resolveHazardGrid`](../../R/dbarts.R) from `gridTime`, which
  already narrows to the subjects the na.action keeps, a later time placed in the last period (read).
- (ran) `rbart_vi` is deprecated ("is deprecated and is removed in dbarts 1.1-0"), so its `rel.scale`
  is not in this slice.

What the setters return today (read).
- [`SamplerBase`](../../src/bartcore/facade.hpp)'s `setResponse` and `setOffset` return nothing, and so
  do [`Sampler::setResponse`](../../src/bartcore/sampler.hpp),
  [`Chain::setResponse`](../../src/bartcore/chain.hpp) and the response models' virtuals. A chain's
  setter always runs its tail: [`reanchorVarianceForest`](../../src/bartcore/chain.hpp) under a variance
  forest and [`restateSavedDraws`](../../src/bartcore/chain.hpp), and, at `updateScale` true only,
  [`translateNamedSds`](../../src/bartcore/chain.hpp), which restates a named leaf sd against the moved
  scale; the sampler then runs [`Sampler::recordAnchor`](../../src/bartcore/sampler.hpp). The engine's
  one typed outcome of a setter is [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp) (accepted,
  rolled back, invalid cut points), which the bridge converts.
- Every part of that tail is a no-op when the transform has not moved: `reanchorVarianceForest` returns
  at an unchanged sigma scale and `restateSavedDraws` at a ratio of 1 and a shift of 0. So the keep
  needs no branch in the chain: it is the false call.
- [`bartcore_setResponse`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setOffset`](../../src/R_interface_bartcore.cpp) return NULL. The holder
  [`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp) carries one warned-once flag,
  `slowCountWarned`, last by the positional-aggregate rule its comment states.
- [`dbarts_sampler_setResponse`](../../src/C_interface.cpp) and
  [`dbarts_sampler_setOffset`](../../src/C_interface.cpp) return a capability status: 1 on a swap, 0
  touching nothing where the coupling admits no conduit; refusals raise. dbarts.h's token
  ([`dbarts_apiToken`](../../src/C_interface.cpp)) folds the signatures, three enumerator lists and the
  struct layouts, so a new return value of an `int` entry moves nothing unless it is named in a folded
  list. [`DBARTS_C_API_MINOR`](../../inst/include/dbarts/dbarts.h) is 0, and the header says the version
  pair does not move before the first release (dec-B60); the C API test pins it
  (["expect_equal(versions, c(1L, 0L))"](../../inst/tinytest/test-capi.R)). A consumer's handshake is
  the version pair unless it defines `DBARTS_REQUIRE_EXACT_ABI`, so a consumer built against another
  hash still loads.
- Membership and served precision differ for Student-t: `TResponse` serves w x lambda x a to its
  Gaussian part, lambda drawn per chain, and counts the rows in by membership (`w * a > 0.0`). A product
  that underflows would make served precision and membership disagree, and chains with it.
- Consumers (read on 2026-10-11, each on the branch named). stan4bart, as its creation-mapping slice
  leaves it ([stan4bart-creation-mapping.md](stan4bart-creation-mapping.md), in build): it resolves its
  spec by `dbartsSpec` in the calling process and creates each chain's sampler by
  `new("dbartsSampler", ...)` in that chain's worker; it then installs one response range for every
  chain through the R method `setResponse(updateScale = TRUE)` on a made-up response (the range's lower
  end on the first row of positive weight, its upper end on every other row), gives the sampler back its
  own response with the range held, and reads the range back, stopping if it did not land; its restore
  does the same per chain. Both of its `dbarts_sampler_setOffset` calls pass the re-derive flag false;
  a continuous fit fixes the residual variance in its family, so no starting sigma is estimated; and it
  does not define `DBARTS_REQUIRE_EXACT_ABI`. treatSens (dbarts-1.0) calls
  `dbarts_sampler_setOffset` on a probit sampler, and `dbarts_sampler_setResponse(fit, ..., 1)` every
  sweep on its gaussian outcome sampler, created by `dbartsData(x, y)` with no weights; it ignores both
  returns. bartCause and bairrtt make no flat call. The deprecated `rbart_vi` calls
  `$setOffset(..., isWarmup)` every warm-up sweep, under the user's weights.

What a correct change must not move.
- (stands) The prototype's rewrite of the unit's loop changed its last bit on unweighted data: at
  n = 64 its base gave 0x1.c828f537ceda5p-3 and the prototype ...da4p-3. On a reference build that
  failed 3 of the 15 two-forest scenarios bitwise at max abs z 0.00; the suite did not see it.
- (ran) Of the 55 scenarios of the main corpus's current file, `equivalence-2b48939b.rds`, one is
  created with a row out: `zeroweights`, 80 of 500 at weight 0 with the smallest response on one of
  them, so [1.3742, 25.2117] becomes [1.6243, 25.2117]. The two that carry a mask (`maskprobit`,
  `maskordinal`) install it after creation, on families that take nothing from the response; the other
  five that pass weights pass positive ones; `hazard` passes none.
- (stands) On a reference build of the prototype, re-ported onto 8ad1fa58, against the main corpus's
  file of that day (re-recorded twice since, for other scenarios): 54 of 55 identical and `zeroweights`
  at max abs z 2.85; two-forest 15 of 15; multinomial 11 of 11; the snapshot files and the exact gates
  pass. The hazard grid was not in the prototype.
- (stands) The suite on the re-ported prototype: 2 failures of 19506 results, both the two-forest pin
  of test-active-rows-pins.R (["Bayesian causal forest, bitwise"](../../inst/tinytest/test-active-rows-pins.R)),
  which compares a sampler created with weights w and then masked against one created with w * a.
- (ran) The suite as it stands, the 93 installed files that name weights or `updateScale`, with the R
  entry points traced: ten create a gaussian or count sampler with a row out (test-active-rows-pins.R,
  test-empty-leaf-veto-weights.R, test-gp-leaves.R, test-nbinom.R, test-pdbart-averaging.R,
  test-pdbart.R, test-pointwise-loglik.R, test-predict-training-rows.R,
  test-prior-init-empty-leaves.R, test-zero-weights.R); all but test-gp-leaves.R are as the prototype's
  suite ran them, and its additions compare a copy, a `setState` and a reload with their own source.
  None reaches rule 4's case at creation, re-derives with a row out on a family that has a range, or
  takes sigest's stand-in or floor under a zero weight; no hazard or `xbart` test passes a weight of 0.

Gates that can fail.
- (ran) [`zeroWeightArm`](../../benchmarks/R/bd-balance.R) creates with every row and so proves nothing
  here. With the weights given at creation and the zeroed cells' response 3 higher, about 1 second in
  `quick`: the tip matches the oracle whose range is over every row (largest abs z 1.4) and misses the
  rows-in oracle by up to 106.6, at the (none) partition; (stands) the prototype matches the rows-in
  oracle (2.1) and misses the other by 40.1. bd-balance.R holds four arms, `zeroWeightArm`,
  `narrowArm`, `constMissingArm` and `weightedArm`, behind one `stopifnot` (read).
- (ran) A count twin, from the fixed-shape arm of negbin-exact.R: 4 rows in of 25 a cell, the rows out 30
  higher, left out by 0/1 weights at creation, the oracle's likelihood over the rows in and its center
  over every row (3.334) or the rows in (1.056). 2 seconds for two seeds. Tip: gap 0.012 to the every-row
  oracle, 0.493 to the rows-in one, tolerance 0.12. (stands) Prototype: 0.484 and 0.003.

Released behavior. (read) 0.9-34 took the range over every row, at creation and at
`setOffset(updateScale = TRUE)`, and had no `updateScale` on `setResponse`; (stands) on its build two
seeded fits that differ only in the responses of the zero-weight rows differ by up to 3.65 in a drawn
fitted value at the rows in. Counts, several forests, masks, Student-t and hazard fits reached no
release, and neither did dbarts.h: 0.9-34's C interface was the R_C_interface.hpp ABI, deleted for 1.0-0.

## The rule

1. Rows. The range (gaussian, Student-t, aft, and under a variance forest or any leaf model), the count
   center and the unit of a forest's sd are computed over the rows in when the sampler is created,
   unweighted: a row of positive weight counts in full.

       dbarts(y ~ x, d, weights = w)    # range of (y - offset)[w > 0]; unit sd((y - offset)[w > 0])

2. Held. No later call moves them unless it carries `updateScale = TRUE` (as today), and none of
   `setWeights`, `setActiveRows`, `setState` in either form, a copy, a reload or a warm start does. No
   argument is added in this slice.
3. Re-derived. `updateScale = TRUE` on `setResponse` and `setOffset` reads the rows in at that call.
   `setData` replaces the data and derives as creation does, rules 1 and 4, its warning included.
   Several forests still refuse `updateScale = TRUE` (dec-B362).
4. Too little to read. Where a row is out and the rows in hold fewer than two distinct values (none in;
   one in; several in at one value):
   - At creation and at `setData`, each number is taken over every row, as if no weight or mask had
     been given, sigest included where one is estimated (dec-B302, dec-B309). One warning, of its own
     class; the starting sigma records no stand-in for that cause (no `dbartsSigmaFallback`, no
     `sigest.fallback`, no verbose line of it), its estimate being the regression on every row:

         fewer than two distinct values of the response are in the likelihood (<m> of <n> observations are in it), so the response's scale and the prior's center are taken from all <n>

   - At a re-derivation through `setResponse` or `setOffset`, the transform and the center in force are
     kept: the call completes as it does with `updateScale = FALSE`, bit for bit (dec-B364). The engine
     returns that outcome as a status and raises nothing; the R methods warn once per sampler (a copy or
     a reload is a new sampler; `setData` is not), with the same class; the flat entries return a value
     of their own for it and warn nothing (dec-B405):

         fewer than two distinct values of the response are in the likelihood, so the response's scale and the prior's center are kept as they were (said once per sampler)

   The value tested is the response less its offset for a range and the count for a center. With every
   row in, a constant response is not this case: it takes dec-B386's window c - 0.5 to c + 0.5 at
   creation and at a re-derivation, with no warning of this class and no kept status, so dec-B364 and
   dec-B386 never apply to the same call.
5. Nothing out, nothing changed. With every row in, every number above has the bits it has today.
6. sigest stays the weighted regression over the rows in and is re-estimated by no call (dec-B314).
   Where that regression leaves no residual degrees of freedom, the sd that stands in for it (dec-B433)
   reads the rows in, unweighted, and so does the floor, by the dense and the sparse routine alike; that
   stand-in is announced as it is today (dec-B429). Nothing here applies under a fixed residual sd
   (dec-B439) or beside a supplied number, where no estimate is made.
7. `xbart`: each fold's range and sigma follow 1, 4 and 6 over the fold's training rows. A function
   given as `sigest` is handed the fold's weights as today and nothing here touches it.
8. Hazard periods. A hazard fit's default periods, and under `breaks = K` its quantiles, are read from
   the times of the subjects in the likelihood when it is created (weight positive, after `subset` and
   the na.action); stated boundaries are as they were, still checked against every subject's time. A
   subject at weight 0 is expanded on that grid, its time in the interval it falls in; one past the last
   period is at risk in every period with no event (dec-B315). With no subject in, the periods are every
   subject's, with one warning of rule 4's class:

       no subject is in the likelihood (0 of <n> subjects are in it), so the hazard periods are taken from all <n>

## Constraints

- No new argument; no bridge entry changes arity; no new facade virtual. Two virtuals of
  [`SamplerBase`](../../src/bartcore/facade.hpp), `setResponse` and `setOffset`, change their return type
  from void to the status (dec-B405): `R CMD INSTALL --preclean` is mandatory on that commit, a stale
  object on the old vtable bus-erroring, as it is on the engine commit for the edited headers.
- The engine never warns (dec-B405): the keep is a returned status, not an exception and not a message.
  Errors stay C++ exceptions the bridge converts.
- The keep is the `updateScale = FALSE` call, literally: the sampler decides it once and calls every
  chain with `updateScale` false, so no chain grows a branch for it; the sampler skips `recordAnchor`.
- Rule 5 is a constraint on the code, not only a test: the unit's loop and the center's loop run as they
  are written today whenever no row is out, and the rows-out case is a separate path. A reviewer checks
  that neither present loop gained a branch.
- One decision of which rows are read at a re-derivation, on membership (weight positive, not masked),
  never on served precision. The sampler makes it from chain 0's response model, by the same read routine
  the model's own `rescale` uses; the unit follows the response model's answer, and the bridge and R read
  the status, none of them testing the rows again. At creation and `setData` R also tests rule 4, for
  sigest and the creation warning; step 8's literals hold the two in agreement.
- dbarts.h keeps every signature. It gains one enumerator list, folded into the token, so
  `DBARTS_C_API_HASH` moves; `DBARTS_C_API_MINOR` stays 0 (dec-B424); the comments of the two
  re-anchoring entries and the contract paragraph on capability statuses change (item 6). An ABI event:
  the lockstep consumers' CI rebuilds against it.
- No state format change and no new record: the range is the model's `response.range`, the unit the
  forests' `anchor`, sigest the data object's, the periods the fit's `$periods`, as now. A re-creation
  restates the record and derives nothing.
- Cut points, a hazard fit's stated boundaries, the number of categories, leaf covariate
  standardization and a multiplier's scale are not touched.
- [`rbart_vi`](../../R/rbart.R) is not touched; it reaches rule 3 and the keep's warning through
  `$setOffset(..., isWarmup)` under the user's weights, as any caller of the R method does.

## The keep's status (dec-B405)

How the keep is decided and how it travels to each caller. Steps 1 and 2 build the read; step 4 builds the
rest.

1. The read. Each response model with a transform or a center (`GaussianResponse`, `NBResponse`, and
   through their Gaussian parts `TResponse` and `AFTResponse`) gains one const read: given the incoming
   vector (a new response, or a new offset against the response in force), whether a row is out and the
   rows in hold fewer than two distinct values of what a re-derivation reads (the response less the
   offset, the observed times for aft, the count for nbinom). Rows in are by membership: user weight
   positive and the mask not zero, each model reading its own (`TResponse` its `w * a`, never its
   served `w * lambda * a`). It is the routine `rescale` and `computeShift` call at creation and
   `setData` for rule 4's every-row fallback, so creation and a re-derivation decide alike. Families
   with neither answer no.
2. Sampler, the one decision. `ResponseUpdateResult { accepted, scaleKept }`, an `enum class` beside
   [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp).
   [`Sampler::setResponse`](../../src/bartcore/sampler.hpp) and
   [`Sampler::setOffset`](../../src/bartcore/sampler.hpp), at `updateScale` true, ask chain 0's response
   model the read against the incoming vector, once (through a one-line forward on `Chain`); on a keep
   they call every chain with `updateScale` false and skip `recordAnchor`, and return `scaleKept`;
   otherwise they run as today and return `accepted`. Chains cannot disagree, since one answer drives
   them all; no chain, response model or tail gains a branch, and on a keep no named sd is restated,
   the arm that does it not being entered. Under per-chain responses
   ([5. Recommendation](../design/per-chain-responses.md#5-recommendation)) the pooled answer is that
   design's to state. [`SamplerBase`](../../src/bartcore/facade.hpp)'s two virtuals return the status and
   [`SamplerFacade`](../../src/bartcore/facade.hpp) forwards it.
3. Bridge. `bartcore_setResponse` and `bartcore_setOffset` return a length-one logical, "warn now": TRUE
   when the engine returned `scaleKept` and the holder's new `scaleKeptWarned` was clear, which the call
   then sets; FALSE otherwise. So the R methods return nothing new, and tinytest sees a later keep only
   through the readers and the twin comparison. `scaleKeptWarned` goes after `slowCountWarned`, last, by
   the holder's positional-aggregate rule. It lives as long as the engine handle: a copy, a reload or any
   other re-creation is a new sampler and warns again on its first keep; `setData` keeps the flag, as
   dec-B364's "once per sampler" reads, and the help says both.
4. R warning. The `setResponse` and `setOffset` methods of [dbarts.R](../../R/dbarts.R) carry the
   bridge's answer out of `bartcoreSamplerSetResponse` and `bartcoreSamplerSetOffset` and raise rule 4's
   keep text through [`warnClassed`](../../R/utility.R), class `dbartsScaleFallbackWarning`, as their
   last act: after the data mirror, the test offset, `recordAnchor` and `storeState()`, so a warning
   turned into an error leaves the R object and its stored state agreeing with the engine, as the flat
   run raises its slow-count warning last.
5. Flat entries. [`dbarts_sampler_setResponse`](../../src/C_interface.cpp) and
   [`dbarts_sampler_setOffset`](../../src/C_interface.cpp) return `DBARTS_SWAP_SCALE_KEPT` (2) when the
   engine returned `scaleKept`, `DBARTS_SWAP_DONE` (1) on any other swap and `DBARTS_SWAP_NONE` (0)
   where the coupling admits no conduit, as today. They never read or set `scaleKeptWarned` and raise no
   warning: a flat keep does not spend the R interface's once. Several forests, and gp leaves holding
   saved draws, refuse `updateScale` true before the engine decides, as today.
6. Contract. dbarts.h gains `DBARTS_SWAP_STATUS_LIST` and `typedef enum { ... } dbarts_swap_status`,
   generated by the list macro as the other three are and placed before `#undef DBARTS_ENUMERATOR`, so
   the token folds it. Both entries keep their `int` return. Their comments: 0 and 1 as before; 2 a swap
   that asked for `updateScale`, on a sampler with a row out of the likelihood (weight 0 or masked),
   whose rows in hold fewer than two distinct values of the response less the offset (of the count, for
   nbinom), so the transform and the count center in force were kept and the swap done as with
   `updateScale` false; a constant response with every row in answers 1 (dec-B386's window); a caller
   that tests whether the swap happened tests for nonzero, not for 1; `updateScale` reads the rows in.
   The contract paragraph on returns ("Contracts common to all entry points": a capability 0 is a fixed
   property, "test it once at setup rather than every sweep") gains that the two swap entries' 2 is a
   per-call outcome and not a capability. "All three ABI enums" becomes four in the header and in the
   fold's comment ([`dbarts_apiToken`](../../src/C_interface.cpp)); `DBARTS_C_API_HASH` is re-baked.
7. The version (dec-B424): "Hash only, we're still pre-release." The pair stays 1/0 and test-capi.R's
   version pins and the header's version paragraph stay. A consumer built against the old header loads
   and runs against the new library; one built against the new header and run against an older
   development build never sees a 2.
8. The C API test. [consumer.c](../../inst/tinytest/capi/consumer.c) exports the three constants;
   [test-capi.R](../../inst/tinytest/test-capi.R): the constants are 0, 1, 2; the hash its new literal;
   the version pins unchanged (item 7). A gaussian sampler with three rows in at one value, the rest at
   weight 0: each entry with `updateScale` true returns 2 and raises no condition (counted with
   `withCallingHandlers`), the transform and five sweeps bit for bit an `updateScale`-false twin's; then
   `$setOffset(updateScale = TRUE)` on the same sampler warns once (the flat call left the flag), and a
   second call none. Two rows in at two values: 1 and the rows-in transform. `updateScale` false: 1.
   Every row in at a constant response: 1. nbinom under a mask keeping one count: 2. probit: 1. The
   stale-minor and wrong-hash handshake arms, and the arm that drives the flat `setOffset` under a named
   sd, as they are.
9. stan4bart's part. stan4bart-creation-mapping (dec-B437) lands on stan4bart's bartcore branch first,
   as its [Against response-scale-rows](stan4bart-creation-mapping.md#against-response-scale-rows)
   recommends; this part is one commit on top of it, pushed with this slice's header commit. After that
   slice stan4bart never asks dbarts to re-derive during a run and installs its own range (Context,
   Consumers). So:
   - Nothing of stan4bart's reads the new value: no flat call of its own can return 2. The read sites,
     the counts and the per-fit warning this item once planned are gone, and dec-B405's "stan4bart's
     branch reads it" has nothing left to read under dec-B437. The value serves treatSens's kind of
     caller, which re-derives every sweep.
   - The install under rule 3. With two or more rows of positive weight the made-up response holds both
     ends of the range on rows in, so the install is exact under this slice as on the tip ((ran) on the
     tip: it lands with the first 100 rows at weight 0, with two rows of positive weight and with one).
     (ran, emulated) Twenty samplers built by hand under zero weights and given one range this way drew
     15 sweeps bit for bit the same whatever range their creation had derived, which is why item 10
     expects no stan4bart draw to move, zero weights or not.
   - One row of positive weight. The rows in then hold one made-up value beside rows out: rule 4's keep,
     the R method's warning inside a worker, and a range that did not land, which stan4bart's read-back
     turns into "internal error: the BART response range was not installed". That plan expects this of
     slice C; it arrives here, where dec-B364 is built for `setResponse`. Only a supplied `bart_range`
     reaches the install with one such row (an estimated range has no width there and nothing is
     installed: read, its plan), and under dec-B445 a supplied pair is the range. The change, in stan4bart's
     install: where fewer than two of the sampler's weights are positive it sets the weights to 1 for
     its two `setResponse` calls and then puts them back, rule 2 holding the range across both
     `setWeights`. (ran) On the tip, ten samplers with one row of positive weight: the range lands, the
     weights are back, and 10 sweeps are the plain install's bit for bit. So nothing a stan4bart user
     sees changes.
   - The test, in stan4bart, on its install: a sampler with one row of positive weight takes the range
     exactly, holds its own weights afterwards, and the install raises no warning. R ~8 lines, the test
     ~15; this plan's records amend that plan's sentence on slice C.
   - Had this slice come first, the stan4bart before that slice would move its zero-weight fits here
     and again there, and its restore, by a response alternating the two ends, would hold one value
     among the rows in wherever the weights are zero on every second row: a keep, and dbarts's warning
     the only sign. Hence the stop in [Budget and stops](#budget-and-stops).
   - treatSens and bartCause change no source: treatSens's probit `setOffset` and its gaussian outcome
     sampler's `setResponse`, created with no weights and no mask, always answer 1, and bartCause makes
     no flat call. Their lockstep CI rebuilds against the new hash.
10. Gates of this part: tests/cpp (step 4); the C API test; step 8's keep arm through the R methods;
    stan4bart's suite and its new test against the slice's library; one seeded continuous stan4bart
    fit with zero weights and one without, each on the base and the slice libraries with one stan4bart
    build, both pairs identical; and once, reported, a stan4bart built against the base header loaded
    and fitted against the slice's library.

## Steps

"Fails today" is what the tip does where the test expects otherwise.

1. Engine, the range. [`GaussianResponse::rescale`](../../src/bartcore/model.hpp) reads the rows in by
   membership (user weight positive, mask not zero), from the response or from the vector slice A has it
   take. The read (the range over the rows in, and whether they hold two distinct values) is one routine,
   separate from the write: item 1 of [The keep's status (dec-B405)](#the-keeps-status-dec-b405) is that
   routine asked of an incoming vector, and a later pooled re-anchor can call it per chain. In rule 4's
   case the constructor and `setData` read every row; at a re-derivation the sampler has already ruled
   the keep out (step 4). Student-t hands its Gaussian part its membership (`w * a`), not the served `w *
   lambda * a`; aft its mask. The response model answers which rows it read, for step 3. tests/cpp,
   beside [`testActiveRowsGaussianDf`](../../tests/cpp/test_model.cpp), two chains where a chain is
   involved:
   - gaussian, Student-t and aft (aft by a mask then `setOffset(true)`), the weights of the rows in
     fractional: the scale and shift equal the unweighted literal over the rows in (fails today: every
     row);
   - no weight at all, and weights all positive: scale, shift and the working response are bit for bit
     those of the present expression, written out in the test;
   - created with none in, one in, four in at one value: the every-row literal; `setData` to those rows:
     the same;
   - the read: with a row out and none in, one in, four in at one value, it answers yes for an incoming
     response and for an incoming offset; with every row in at a constant response, no; Student-t with a
     row in whose served precision underflows (a subnormal weight) still counts it in;
   - held: after `setWeights` and `setActiveRows` the transform is unchanged and 20 sweeps equal an
     untouched twin's under the same weights;
   - `setResponse(true)` and `setOffset(true)` under zero weights and under a mask read the rows in;
     `setData` reads the rows in under the new weights.
2. Engine, the center. [`NBResponse::computeShift`](../../src/bartcore/model.hpp) takes the mask where
   one is installed; the unmasked call is the present loop. Its two re-deriving callers pass the mask,
   and its read (item 1 of the keep's status) answers yes where the rows the mask keeps hold fewer than
   two distinct counts. tests/cpp, beside [`testNBLogMeanAnchor`](../../tests/cpp/test_model.cpp): the
   center after a mask and `setOffset(true)` equals the literal over the rows the mask keeps, with and
   without an exposure offset; unmasked, bit for bit the present value; the read on a mask keeping one
   row, or rows of one count: yes; held across `setActiveRows`. The comments of
   [`NBResponse::setActiveRows`](../../src/bartcore/model.hpp) and
   [`AFTResponse::setActiveRows`](../../src/bartcore/model.hpp) stop saying "full-data".
3. Engine, the unit. [`scaledResponseSd`](../../src/bartcore/chain.hpp) is left as it is and is what
   runs when the response model read every row; otherwise a second routine takes the sample sd over the
   rows it read, n - 1 of those rows in the divisor. tests/cpp, beside
   [`testForestCalibration`](../../tests/cpp/test_sampler.cpp): with rows out the map's anchor equals the
   literal sd to 1e-14; with none out it is bit for bit the present loop's value at n = 64, 100, 137 and
   250 (the size at which the prototype moved a bit among them); a finite carried anchor still wins.
4. The keep's status, engine to C API, as [The keep's status (dec-B405)](#the-keeps-status-dec-b405)
   items 2 to 8: the sampler's decision and its status, the facade, the bridge's logical and holder flag,
   the R warning, the flat entries' values, the header's list and comments and the token's fold, and the
   version left at 1/0 (item 7). `--preclean`. tests/cpp, two chains: re-derived with a row out and none
   in, one in, four in at one value, through `setResponse` and `setOffset`: `scaleKept`, and the
   transform, the working response and 20 sweeps of both chains bit for bit an `updateScale`-false
   twin's, `getAnchor` unchanged (fails today: every row); the same under a variance forest, with kept
   draws of constant leaves (gp leaves holding saved draws refuse the call), and with a named sd, which
   the keep leaves where it was; rows in at two values: `accepted` and the rows-in transform; every row
   in at a constant response: dec-B386's window and `accepted`; the facade forwards the status
   ([`SpySampler`](../../tests/cpp/test_facade.cpp) takes the new return type). The C API test as item 8.
5. R, the rows and the warnings. One function in [spec.R](../../R/spec.R) answers rule 4 for a data
   object and a mask and raises the creation warning, for the families that take something from the
   response (gaussian with any residual distribution, nbinom) and no other. Its callers:
   [`resolveSamplerSpec`](../../R/spec.R), after [`enforceWeightPolicy`](../../R/spec.R) (so `dbartsSpec`
   and every door that resolves a spec, stan4bart's included);
   [`bartcoreSamplerSetData`](../../R/bartcore.R), for `setData`; and [`xbart`](../../R/xbart.R), once,
   after its own `enforceWeightPolicy` call and before the all-rows estimate that fixes `sigmaPerFold`.
   It is not called at a re-derivation: there the engine's status answers (step 4).
   - sigest. In rule 4's case [`estimateStartingSigma`](../../R/spec.R) is handed a copy of the data
     without its weights, so the regression, the row count its verbose lines print, the stand-in and
     the floor all read every row, and a `dbartsSigmaFallback` is signaled only if that regression too
     leaves no residual degrees of freedom. Otherwise the regression keeps the weights, which drops the
     rows out as today, and [`estimateSigmaFromLinearModel`](../../R/utility.R) hands
     [`floorSigmaEstimate`](../../R/utility.R) and [`floorMarginalSigma`](../../R/utility.R) the
     residual of the rows of positive weight, by either routine. With no weight the vector they are
     handed is the one they are handed today (rule 5).
   - `xbart`. In rule 4's case the all-rows estimate is the regression without the weights, so
     `sigmaPerFold` is that regression's routine unless every row too leaves it nothing. Each fold tests
     rule 4 over its training rows with no warning ([`xbartRunUnits`](../../R/xbart.R)): in its case the
     fold's regression is unweighted; otherwise weighted, its stand-in and its floor over the training
     rows in. So `xbart` warns once and a fold falls back in silence; the verbose count of units that
     started from the response's sd counts as it does today.
6. R, creation under a mask. Of the families whose 0/1 weights become the mask
   ([`isMaskedWeightFamily`](../../R/spec.R)) only nbinom takes anything from the response. For it,
   where rule 4's case does not hold, [`dbarts`](../../R/dbarts.R) and
   [`bartcoreSamplerSetData`](../../R/bartcore.R) follow their install of the mask with the sampler's own
   `setOffset(offset in force, updateScale = TRUE, updateState = FALSE)`, before the sampler is handed
   back: it draws nothing ((ran) nbinom twins, one given the call with every row in, draw alike, under a
   named sd too) and records the range; a named sd is restated by the engine at that call (read). In
   rule 4's case the follow-up is skipped and the every-row center of creation stands, under step 5's
   warning, so no keep can arise there. A sampler built by hand from [`dbartsSpec`](../../R/spec.R) and
   `new("dbartsSampler", ...)` that installs `spec$active` itself is created with every row in and, for
   a count fit, keeps the every-row center; in this slice its help says so and gives no call. The
   documented call, `setActiveRows(spec$active, updateScale = TRUE)` with its loop form, comes with
   slice C, which builds it (dec-B363); the `setOffset` route is not documented.
7. R, the hazard periods. [`expandDiscreteTimeHazard`](../../R/dbarts.R)'s two callers in
   [`dbarts`](../../R/dbarts.R), the direct response and the formula's, pass as `gridTime` the times of
   the subjects kept by the na.action whose weight is positive, all of them where no weight is given or
   none is positive (rule 8's fallback, raising its warning there). [`resolveHazardGrid`](../../R/dbarts.R)
   takes the default grid and the `breaks = K` quantiles from it unchanged; stated boundaries keep their
   check over every subject. A subject whose time is past the last period is expanded to every period
   with a response of 0. Nothing in the engine.
8. tinytest, a new file test-response-scale-rows.R, fixtures with the rows out extreme in the response
   and no missing value in any predictor (an install then draws no direction, so states compare bit for
   bit):
   - Inert. Two fits that differ only in the response of the rows out at creation are `identical()` at
     the rows in and within 1e-12 at the rows out, in draws, sigma and every reader: gaussian, a fixed
     sigma, Student-t, a variance forest, linear and gp leaves (covariates equal across the pair), a
     named sd, two forests on a factor and on a number, an offset, nbinom, and `bart`, `bartBT` and
     `xbart` once each. Fails today on every arm.
   - Literals, with 0/1 weights. `response.shift`, `response.scale`, `k.scale`, `prior.mean`, the
     forests' `anchor` and sigest equal values written from the rows in, and equal the `subset` fit's.
     `k.scale` is the data's scale over the rows in under the k and the sd spelling alike, and an
     sd-spelled fit's k is `k.scale` over the named sd. With fractional weights in, sigest is the
     weighted fit and the `subset` fit's is not, so that arm asserts the range, the center and the unit
     only.
   - Held. The readers are unchanged across `setWeights`, `setActiveRows`, a copy and a reload, read from
     the engine (`getLeafPrior`) and not from the record alone; five sweeps after a copy agree with the
     original's to 1e-10.
   - `updateScale = TRUE` on `setResponse` and `setOffset` under zero weights and under a mask gives the
     readers of a sampler created under those rows, and no warning; several forests still refuse.
   - The keep (dec-B364, dec-B405). Gaussian under zero weights and nbinom under a mask, re-derived
     through each of `setResponse` and `setOffset` with none in, one in, and several in at one value: the
     readers and five sweeps equal an `updateScale = FALSE` twin's; exactly one warning of the new class
     on the first such call, its text, counted with `withCallingHandlers`, and none on a second through
     either setter, nor after a `setData` (the flag survives it); a copy and a reload warn again on their
     own first keep; a constant response with every row in: dec-B386's window and no warning. A warning
     turned into an error (`options(warn = 2)`) leaves `$state` and the test offset agreeing with the
     engine. The flat entries' side is in the C API test (step 4).
   - Rule 4 at creation. Every weight 0; a mask of zeros (nbinom); one row in; rows in at one value:
     created, exactly one warning of the new class, counted with `withCallingHandlers`, its text as above
     with m and n; the readers and sigest equal those of the fit with no weights; no
     `dbartsSigmaFallback` is signaled, a `bart` fit has no `sigest.fallback`, and a verbose fit's
     estimating line counts every row; no warning where every row is in and the response is constant;
     none for probit under a mask of zeros. Under fractional weights with a zero among them the data
     object's validity warning is raised as well, as today, and pinned by its own pattern exactly where
     it fires. `setData` to such rows: every row, the same one warning. `dbartsSpec` with one row in:
     the one warning, from the spec.
   - sigest. Three rows in and two predictors: the sd of the rows in (fails today: every row's), no
     warning, one `dbartsSigmaFallback`, `sigest.fallback` TRUE on a `bart` fit, the verbose line naming
     3 rows; a response the predictors fit exactly among the rows in: the floor from the largest
     response less offset among them. The same under `sigest = "dense"` and `"sparse"`, for a caller's
     `dgCMatrix` x and for an indicator expansion (`factors = "indicators"`). With rows out and residual
     degrees of freedom left, both routines give the dense design's sigest over the rows in, as today (a
     pin). A fixed residual sd and a supplied number: the slot holds them, with rows out and in rule 4's
     case.
   - `xbart` with 60 rows in at one value: one warning of the new class, none per fold, and each fold's
     sigma the unweighted fit over its training rows.
   - Creation routes for a count fit: `dbarts` with 0/1 weights and `setData` with 0/1 weights give the
     rows-in center; `dbartsSpec` then `new("dbartsSampler", ...)` and `setActiveRows(spec$active)` gives
     every row's, pinned as the documented difference of this slice.
   - Hazard. Two subjects at weight 0 holding times no subject in holds: `$periods` equals the fit's
     without them, and the subjects in expand to that fit's rows (covariates, period, response) bit for
     bit; the subject past the last period has every row at 0. Two fits that differ only in one such
     subject's time within one period of that grid (3.5 against 3.7 on periods 1 to 5) are
     `identical()` (fails today: the periods differ). `breaks = 3` gives the quantiles of the subjects
     in; stated boundaries are unchanged. Every weight 0: every subject's periods and exactly one
     warning of the new class, its text. `hazard.logistic` still refuses a zero weight.
   - A mask over w against weights w * a: equal when both samplers are created under w, two fits when
     not (the readers differ), stated as such.
   - A state across row sets. A single-forest, a gp and a two-forest state from a sampler created under
     other zero weights, whose range therefore differs, installs as stored (dec-B418): `setState(state)`
     returns TRUE and `setState(state, forceUpdate = TRUE)` NULL, with no warning; the recipient's range
     and `getLeafPrior` are unmoved; its stored forests are the donor's bit for bit, a gp forest's
     compared by `tree.params` and structure. (ran) On the tip, across two ranges, all of this holds; a
     gp forest's `tree.values` differ after any install, a sampler's own state included.
9. Restated pins. In test-active-rows-pins.R the two-forest pin creates both samplers under w and masks
   or reweights after; the gaussian and Student-t pins above it get the same repair, since they pass only
   while the fixture's extremes are on rows in. The comments
   ["The log-mean shift is the full-data one by design"](../../inst/tinytest/test-active-rows-pins.R) and
   ["the response transform is the FULL-data one by design"](../../inst/tinytest/test-active-rows-pins.R)
   say "creation's, every row then being in". Another existing test that fails is first looked for
   among the ten files Context names.
10. Exact gates. bd-balance.R gains a fifth arm beside [`zeroWeightArm`](../../benchmarks/R/bd-balance.R),
    `narrowArm`, `constMissingArm` and `weightedArm`: the weights at creation, the zeroed cells' response
    3 higher, the oracle's range over the rows in; its `stopifnot` over the arms takes the fifth.
    negbin-exact.R gains a third arm beside its fixed-shape one, as measured in Context, with the
    oracle's center ([`shiftOf`](../../benchmarks/R/negbin-exact.R)) over the rows in.
    [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) runs the bd-balance arm beside its
    `zeroweight` one; the negbin arm runs inside its script. About 3 seconds more in `quick`.
    benchmarks/README.md names both.
11. The baseline, in its own commit after review, as d820097c did it. `zeroweights` is recorded fresh on
    a `--preclean` reference install at the slice's newest package-code commit, which names the file,
    and merged into a copy of `equivalence-2b48939b.rds` in its scenario order; the other 54 are
    carried, being bitwise identical in the partition and identical objects in the two files; a
    `--bitwise --strict-coverage` compare from a second `--preclean` reference install reproduces all
    55. The [MANIFEST](../../benchmarks/baselines/MANIFEST) row is in its predecessor's form: the class
    (posterior-changing: the range of a fit created with a zero weight) and its rulings; the partition
    against 2b48939b with its counts (54 of 55 identical, 55 compared, 0 skipped) and `zeroweights`'s
    number of summaries and max abs z, which size the move; the two-forest and multinomial files 15 of
    15 and 11 of 11 and keeping their recordings; the snapshot files; the oracle (P17), step 12's row 3
    run on the scenario's own data, where the slice's draws are bitwise the base build's given the 80
    out rows' responses moved inside [1.6243, 25.2117]; the exact gates' count; the build mode and where
    it was recorded; 2b48939b demoted to historical. The file's name is repinned where d820097c
    repinned it: cpp-tests.yaml, equivalence.yaml, mutation-battery.R and feature-matrix.md.
12. The pair script, benchmarks/R, run on the base and the slice builds with the same seeds, 20 sweeps:
    1. no row out, each family and composition (hazard included), weights all positive, a mask after
       creation, `updateScale = TRUE` with no row out, `xbart`: `identical()`;
    2. rows out, no extreme among them, one forest; a hazard fit whose subjects at weight 0 hold only
       times of subjects in: `identical()`;
    3. rows out with an extreme, gaussian and Student-t: the slice's fit equals, at the rows in, the base
       build's fit with the out rows' responses moved inside the range of the rows in: `identical()`;
    4. a count fit whose masked rows are structural zeros, as many out as in, no offset: the slice's fit
       equals the base build's with the masked rows given the counts of the rows in: `identical()`;
    5. several forests with rows out, a hazard fit with an out subject off the grid, and each mover of
       the `rng:` block: must differ, and by how much.
13. Mutations (Verification).
14. Records. Help: the `weights` item of [bart.Rd](../../man/bart.Rd) and [dbarts.Rd](../../man/dbarts.Rd)
    (the sentences holding ["still counts toward"](../../man/bart.Rd) and
    ["the response scale held equal"](../../man/bart.Rd)) say that a row of weight 0 takes no part in the
    response's scale, the count center or a forest's default sd, that it stays in the design and the fit
    is close to the fit on the remaining rows without being it, that the step is at 0 and not at a small
    weight, and that a mask or weights set after creation leave the numbers creation set, so the three
    ways to leave a row out differ by when it left; the `sigest` item of the pages that take it says
    the response's sd that stands in for the regression is over the observations of positive weight;
    the nbinom paragraph names the rows of c; the hazard paragraphs of both and of
    [dbartsFamilies.Rd](../../man/dbartsFamilies.Rd)
    (["sorted distinct observed times"](../../man/bart.Rd)) say the subjects in the likelihood;
    bart.Rd's paragraph of warning classes (["Warning classes."](../../man/bart.Rd)) gains
    `dbartsScaleFallbackWarning`; [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the
    response item (["which reads every row"](../../man/dbartsSampler-class.Rd) becomes the rows in),
    the weights item, the `active` item
    (["stay the FULL-data calibration"](../../man/dbartsSampler-class.Rd) becomes creation's),
    `updateScale` with the keep and its once-per-sampler warning, and that `setData` and a sampler made
    by hand from another's parts derive again; [dbartsSpec.Rd](../../man/dbartsSpec.Rd), the `active`
    value: a count sampler made by hand keeps the every-row center in this slice; xbart.Rd, `weights`.
    dbarts.h: the two entries' comments, the new list and the contract paragraph (step 4).
    NEWS, user-visible changes, one item: "An observation of weight zero no longer takes part in the
    response's scale, which the leaf prior is stated against, so a fit with zero weights has the leaf
    prior of the fit on the remaining rows; 0.9-34 took the scale over every row. The observation is
    still in the design." The keep, its status, the flat entries' value, the stand-in sd and the hazard
    periods reached no release and take no item. A design note, docs/design/response-scale-rows.md,
    with its index row: the rule, the inventory, the measurements, the named-sd case (a prior centered
    at 52 where the rows in sit near 2), why membership and not weight, the keep and how its status
    travels (engine, R, flat entries), the hazard grid, the four slices.
    [active-rows-mask.md](../design/active-rows-mask.md),
    [negative-binomial.md](../design/negative-binomial.md), [nbinom-log-mean.md](nbinom-log-mean.md) and
    [survival.md](../design/survival.md) take a dated amendment where they say "full-data", "all rows" or
    the hazard grid's source; [error-style.md](../design/error-style.md)'s list of warning classes takes
    the new one; [stan4bart-creation-mapping.md](stan4bart-creation-mapping.md)'s sentence on slice C
    names this slice and the lifted weights. TODO: the `response-scale-rows` entry is cut to slices C and D;
    `hazard-periods-rows-in` is removed; `starting-sigma-fallback-zero-weights` is cut to its open half,
    whether the stand-in sd should be weighted. The ledger entry for the calls made. The landing note
    here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, clean under ASan
  and UBSan; the full tinytest suite, expected clean after step 9 (the two-forest pin's 2 failures
  without it); the new file and test-capi.R under ASan on the R-loaded path; tests/cpp plain and under
  the address, undefined-behavior and leak sanitizers and tinytest once on x86-64 Linux, as the recent
  engine landings ran them.
- On a reference build, the three compares, counted per scenario with no `max |z|` line but the one
  named:
  - main corpus against `equivalence-2b48939b.rds`: 54 identical; `zeroweights` differs and is the only
    one, being the one scenario created with a row out and its smallest response on such a row
    (Context). Its statistical verdict against that file is reported, not gated: the posterior changed
    (the re-ported prototype's was max abs z 2.85). Then all 55 identical against the new file.
  - two-forest corpus against `bcf-equivalence-1b7d730c.rds`: 15 of 15 identical; its weighted scenarios
    pass positive weights and its mask goes in after creation (read). Fewer is rule 5 broken in the
    unit's loop.
  - multinomial corpus against `multinomial-equivalence-80b1c8d4.rds`: 11 of 11 identical; the family
    takes nothing from the response and no weight.
  - the four `test-reproducibility-*.R` files pass unchanged, 27 results; none is regenerated. Their
    weights are 1 and 2, and they install no mask and re-derive nothing (read).
- Every gate of [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) in `quick`, unchanged (41
  runs at the current baseline's recording; each that leaves a row out does so after creation, read),
  and the two new arms. Each new arm is also run once on the base build, where it must FAIL (measured:
  by 106.6 standard errors and by 0.49 against 0.12), and that run is reported.
- The pair script: rows 1 to 4 `identical()`, row 5 differing.
- The C API test passes on the slice's library; it fails on the base library (the constants, the hash,
  and 1 where 2 is expected).
- stan4bart: item 10 of [The keep's status (dec-B405)](#the-keeps-status-dec-b405).
- Mutations, each to fail the named check:
  - the range over every row: step 8 "Inert" and "Literals"; the bd-balance arm;
  - the center ignoring the mask: step 2's literal; step 8's nbinom arm; the negbin arm;
  - the unit over every row: step 3's literal; step 8's two-forest arms;
  - the unit's present loop replaced by the rows-in loop when no row is out: step 3's bits; the
    two-forest compare;
  - the center's present loop given the mask test: step 2's bits (and the main corpus's count scenario
    if a bit moves);
  - weighted, not membership (a row counted by its weight): step 1's literal with fractional weights;
  - rule 4 off (one row in gives the window of width 1): step 1's "one in"; step 8 "Rule 4 at creation";
  - rule 4 on with every row in (a constant response warns): step 8's "no warning";
  - the keep dropped (the sampler re-derives in rule 4's case, reading every row): step 4's twin tests;
    step 8 "The keep"; the C API keep arms;
  - the read on served precision, not membership: step 1's Student-t underflow arm;
  - the status dropped (the keep returns `accepted`): step 4's status checks; step 8's warning count;
    the C API's 2;
  - a named sd restated on a keep: step 4's named-sd twin;
  - the R warning raised before `storeState()`: step 8's `options(warn = 2)` arm;
  - the flat entries warning (reading or setting `scaleKeptWarned`): the C API arm counting no condition,
    and the R call after it warning once;
  - the keep's warning raised on every such call: step 8's "none on a second";
  - the creation fallback at `setData` replaced by the keep: step 1's `setData` test; step 8's `setData`
    arm;
  - creation's re-derivation under the mask dropped: step 8 "Creation routes"; the negbin arm; and run
    for probit too: the suite's seeded probit-under-mask pins, if the call draws;
  - sigest's stand-in or its floor left on every row, in the dense and in the sparse routine: step 8
    "sigest"; the stand-in weighted: its literal under fractional weights;
  - the stand-in recorded in rule 4's case (the estimate left weighted there): step 8's absent
    `sigest.fallback` and its sigest literal;
  - `xbart`'s folds warning each: step 8's `xbart` count of one;
  - the hazard grid over every subject: step 8 "Hazard", the periods and the `identical()` pair;
  - `updateScale = TRUE` left reading every row: step 1's last bullet; step 8's fourth;
  - in stan4bart, the weights left in place at the install: its test's range.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; inst/NEWS.Rd parses with a non-NULL result; dbarts.h compiles as C (`gcc -std=c99 -pedantic`,
  as its own comment requires); `R CMD check --as-cran` on a tarball from a clean copy. No new Rd topic.
- Cost: creation and the two response setters at `updateScale` true gain a comparison per row where a
  weight vector or a mask is installed, and nothing where neither is (no row can be out); a hazard fit
  one subset of its times; no sweep path. treatSens's per-sweep `setResponse(..., 1)` on a sampler with
  no weights and no mask gains nothing, and stan4bart makes no re-deriving call in a sweep. No bench
  compare.

## Budget and stops

Planned ~2065 lines, by step:

| step | package code | tests | the rest |
|---|---|---|---|
| 1 the range | engine 60 | tests/cpp 130 | |
| 2 the center | engine 30 | tests/cpp 60 | |
| 3 the unit | engine 45 | tests/cpp 50 | |
| 4 the keep's status | engine 40; facade, bridge, flat entries, dbarts.h, token 75; R 15 | tests/cpp 95; C API test 75 | |
| 5 rows, warnings, sigest, `xbart` | R 105 | | |
| 6 creation under a mask | R 20 | | |
| 7 hazard periods | R 30 | | |
| 8 the new tinytest file | | 540 | |
| 9 restated pins | | 40 | |
| 10 exact arms | | | benchmarks, workflow, README 180 |
| 11 the baseline | | | MANIFEST row and repins 15 |
| 12 the pair script | | | 70 |
| 14 records | | | help, header comments, NEWS 125; design note, indexes, amendments, TODO, ledger, landing 265 |

Package code ~420, tests ~990, the rest ~655. The three most recent slices ran 1.8, 2.2 and 2.3 times
their plans (setState's two forms, the sparse starting sigma, the state's record of missable columns),
with 68, 52 and 83 percent of what they added in tests. Assumed here: the tests run to about
2.5 times their figure and everything else to about 1.3, which is ~3900 lines, inside a forecast of 3700
to 4750 (1.8 to 2.3 times). An overrun in tests is expected and is not a reason to stop; one in package
code means the design moved. stan4bart's part, ~25 lines, is outside these counts.

Stop and report when:
- the diff passes ~4750 lines, or the package code (engine, facade, bridge, flat entries, header, R)
  passes ~630, 1.5 times its figure;
- a recorded scenario other than `zeroweights` moves, the two-forest or the multinomial corpus is not
  bitwise, or a snapshot moves;
- the keep is not the `updateScale = FALSE` twin bit for bit on any arm (two chains, a variance forest,
  kept draws, a named sd);
- the keep needs a branch in a chain or a response model's setter, or a decision anywhere but the
  sampler;
- the status needs a facade virtual beyond the two return types; or dbarts.h needs more than the
  enumerator list, the hash, the comments and the contract paragraph (any signature change is a stop);
- an existing test fails that step 9 does not name and whose file is not among Context's ten;
- stan4bart-creation-mapping has not landed on stan4bart's bartcore branch when this slice's header
  commit is ready to push (item 9 says what landing first would do); stan4bart's part needs more than
  the lifted weights in its install and their test; a stan4bart draw moves; or a stan4bart built
  against the base header fails to load or to fit;
- a hazard fit with every subject in changes.

Expected verdicts: tests/cpp, tinytest and the C API test clean; the main corpus 54 identical and
`zeroweights` moved, then 55 identical against its new file; two-forest 15 and multinomial 11 identical;
the snapshots unchanged; every exact gate passing and the two new arms passing on the slice and failing
on the base; the pair script's rows 1 to 4 identical and row 5 differing; stan4bart's suite clean and
its two fit pairs identical.

## What a consumer sees

Read in each package's source (weights, `subset`, masks, `updateScale`, the flat setters' return); no
consumer suite was run.

- stan4bart, once its creation-mapping slice has landed (item 9 of
  [The keep's status (dec-B405)](#the-keeps-status-dec-b405)). Each chain's sampler is created by hand
  under the fit's weights, zeros allowed, and stan4bart at once installs its own range over creation's;
  nothing re-derives afterwards, so no draw of a continuous fit moves. One thing a user can meet: where
  the rows of positive weight share one response value `dbartsSpec` raises rule 4's creation warning,
  once, in the calling process, which is true of the sampler as created and not of a range stan4bart
  then installs from the response less its offset and parametric part. Its binary fits install a mask after
  creation on a probit sampler, which takes nothing from the response. It rebuilds against the new
  header in lockstep; a build made before loads and runs unchanged.
- bartCause passes the user's weights and `subset` at creation, to `bcf` too, and changes neither
  afterwards. A fit with a zero weight moves as the `rng:` block says, under `bcf` by the unit as well.
  No test of its suite passes a zero weight (searched). No flat call.
- treatSens fits a probit treatment model and re-anchors it through the flat `setOffset`, and swaps its
  gaussian outcome sampler's response through the flat `setResponse` with `updateScale` true every
  sweep, ignoring both returns. The probit sampler takes nothing from the response, and the gaussian one
  is created with no weights and no mask, so no row is out: both always answer 1 and its draws do not
  move. bairrtt uses no weights, no mask and no `updateScale`.
- `rbart_vi`, deprecated, calls `$setOffset(..., isWarmup)` every warm-up sweep: a fit with a zero
  weight re-derives over the rows in, and one whose rows in hold fewer than two values keeps the scale
  with one warning.

## Interactions

Landed and built against ([Recheck against bartcore](#recheck-against-bartcore)): the per-gap
default rule, the k-internal slice, setState's two forms, the state's record of missable columns and the
sparse starting sigma. Still ahead:

- stan4bart-creation-mapping ([stan4bart-creation-mapping.md](stan4bart-creation-mapping.md)): first, in
  stan4bart; no dbarts code. Item 9 of [The keep's status (dec-B405)](#the-keeps-status-dec-b405).
- [forest-sd-unit.md](forest-sd-unit.md): after this slice, whose unit it restates.
- Per-chain responses ([per-chain-responses.md](../design/per-chain-responses.md)): dec-B437 took the
  range fixed at creation for stan4bart; if the per-chain data is built it follows this slice, whose
  range readers it reuses. This slice keeps the read (the range over the rows in and whether they hold
  two values) one routine apart from the write, so a pooled re-anchor can call it per chain; the status
  stays one answer per sampler.
- The record slices of TODO `state-frame-prior`: slice D waits on them.

## Out of scope, and where it goes

- Slices C and D, below.
- `rbart_vi`'s `rel.scale`: left as it is; the function goes in 1.1-0.
- A sampler made by hand from another sampler's control, model and data derives the range again under
  the weights in force and keeps the unit from the control: an existing seam, opened further by this
  rule; named in the help here, closed by the record slices of TODO `state-frame-prior`.
- A setter for the residual prior, which is how a caller moves what sigest set (dec-B314): TODO
  `sigma-prior-setter`.
- Whether the sd that stands in for an undefined starting-sigma regression should be weighted: TODO
  `starting-sigma-fallback-zero-weights`, cut to that question. This slice's rulings say which rows.
- A documented method that sets a sampler's response range, which stan4bart's made-up response stands
  in for: deferred in that plan.

## Slices C and D, not planned here

C. `updateScale` on `setWeights` and `setActiveRows`. Waits on: this slice. Settled, to carry into its
   plan:
   - What a re-derivation does to the live chain on a call that changes no response: as
     `setResponse(updateScale = TRUE)` does, the live leaves keeping their internal numbers and the
     kept draws converted (dec-B327). The critic measured the live fit at the rows in moving by up to
     15.07 with no data change until the next sweep; the ruling accepts it.
   - Position. The new argument goes after `updateState`, so 0.9-34's positional `setWeights(w, TRUE)`
     keeps its meaning (dec-B327, dec-B381).
   - A sampler of several forests keeps refusing `updateScale = TRUE` in 1.0-0, the help saying the
     response's scale and a re-derived basis scale (dec-B303) may span different rows (dec-B362).
   - A re-derivation whose rows in hold fewer than two distinct values keeps the prior in force and warns
     once per sampler (dec-B364), as this slice builds for the two response setters; slice C extends it
     to the two new arguments, the sampler deciding from the same read and returning
     `ResponseUpdateResult` through the same facade path and holder flag (dec-B405). No flat entry sets
     weights or a mask, so the C API gains nothing, and stan4bart's range install, which goes through
     `setResponse`, is settled here.
   - The by-hand count sampler's call, `setActiveRows(spec$active, updateScale = TRUE)`, and the loop
     form `setActiveRows(mask, updateScale = i <= n.burn)`, written into dbartsSpec's help (dec-B363).
   - A help sentence, and a test, that `setState` does not undo a re-derivation.
D. dec-B304: a linear or gp leaf's covariate center and scale, and the default lengthscale, over the
   rows in at creation, held, re-derived by `updateStandardization`; `xbart`'s folds taking their own.
   Waits on: the record slices of TODO `state-frame-prior` (a state install today replaces the
   recipient's center and scale and reports a clean install, and a re-creation without a state derives
   them afresh, so there is nowhere to hold them); `updateStandardization` itself (dec-B233, not built).
   Open: creation under a mask must derive the standardization under it too, which this slice's
   `setOffset` step does not name, so either creation takes the mask or the install re-derives two
   things; the leaf is built from the column store, which holds no weights.

## Open calls

None. The dependency held here on 2026-10-09, whether stan4bart stops re-deriving during warm-up, is
ruled (dec-B437) and item 9 is written to it.

## Calls made in planning

Agent-made, each for the maintainer's later mark. Those ruled or approved since are in the ledger and not
repeated here: dec-B309 (the every-row fallback, not a window of width 1), dec-B327 and dec-B381 (slice
C's argument and its live chain; the count tested for a count fit; only families that take something
from the response warn; in rule 4's case sigest is the fit without the weights, under the one warning;
`xbart` warns once; `rbart_vi` cut) and dec-B424 (the hash alone). Still the agents', from the design's
critique (2026-10-07), the planner, the revisions (2026-10-09) and the recheck (2026-10-11):

- Four slices: A the aft re-anchor fix; B creation; C the two new arguments; D dec-B304 with the record
  slices. Not taken: the design's one slice, estimated at 1380 lines before dec-B304 was ruled.
- B is creation only, with `updateScale = TRUE` on the two existing response setters reading the rows
  in. No new argument or bridge arity.
- With no row out every number keeps its bits, as a constraint on the loops and a gate; two exact arms
  that fail on the base build.
- Creation under a mask is R's install followed by the sampler's own `setOffset`, for nbinom only.
  Not taken: the engine's creation taking the mask, which changes the creation entry's arity. The
  by-hand route's documented call is dec-B363's, with slice C.
- The unit follows the response model's decision of which rows were read; it does not test again.
- The warning's class is `dbartsScaleFallbackWarning` (with `dbartsWarning`), for the creation text, the
  keep's and the hazard grid's; it is listed in bart.Rd's paragraph of warning classes and in
  error-style.md, as its neighbors are.
- The negbin arm runs inside negbin-exact.R; the bd-balance arm is a named arm the workflow runs beside
  `zeroweight`, not in the loop of the two cut-grid arms.
- NEWS names the response's scale only, with what 0.9-34 did.
- dec-B364 is built here for `setResponse` and `setOffset` (the orchestrator), replacing the plan's
  every-row reading in silence at a re-derivation; slice C extends it.
- Rule 4's "where a row is out" stays, so dec-B364 (fewer than two distinct values) and dec-B386 (a
  constant response's window) do not overlap: a constant response with every row in takes the window and
  no warning, at creation and at a re-derivation (the orchestrator).
- `setData` in rule 4's case is treated as creation: every row, with dec-B302's warning. dec-B364 names
  four setters and not `setData`, which replaces the data (the orchestrator). The alternative was the
  keep, which would hold a scale from data no longer there.
- The hazard periods (dec-B315) come into this slice: TODO `hazard-periods-rows-in` schedules them "R
  only; with the creation slice of response-scale-rows, whose rule it follows". `breaks = K` takes its
  quantiles from the subjects in as the default grid does, dec-B315 leaving only stated boundaries as
  they were; the alternative was every subject's quantiles, which a far-out subject at weight 0 moves
  (20 against 5). A hazard fit's only fallback is no subject in: one subject in, or several at one time,
  is a grid of one period, a defined fit.
- To dec-B405. The sampler decides the keep once, from chain 0's response model by the read the model's
  own `rescale` uses, against the incoming vector, on membership; it calls every chain with
  `updateScale && !keep`, skips `recordAnchor` on a keep, and returns
  `ResponseUpdateResult { accepted, scaleKept }` through the facade (the coordinator, adopting the
  scoped critique's design). Not taken: each response model deciding and returning the status up through
  the chain, where chains could disagree under Student-t's served precision; and a flag read after the
  call through a new facade read, which dec-B405's "a status the engine's setResponse and setOffset
  return" rules out.
- To dec-B405. The once-per-sampler flag is the holder's (`scaleKeptWarned`), spent only by the R
  entries, so a flat keep leaves the R interface's warning for the R methods' first keep. The
  alternative, a flag on the R object, needs a new field on the reference class.
- To dec-B405. The flat value is named by a fourth enumerator list, `dbarts_swap_status`
  (`DBARTS_SWAP_NONE`, `DBARTS_SWAP_DONE`, `DBARTS_SWAP_SCALE_KEPT`), folded into the token, the entries'
  return type left `int`: a bare new value would not move the hash dec-B405 names. The alternative, the
  return type changed to the enum, changes two signatures for no reader's gain.
- The R warning is the dbarts.R methods' last act, after `storeState()`; raised in the bartcore.R
  helpers it left `$state` and the test offset stale when the warning became an error.
- Recheck. stan4bart reads no flat return: under dec-B437 it makes no re-deriving flat call, so
  dec-B405's clause on its branch has nothing to apply to. The alternative, read sites and counts for a
  value that cannot occur, was this plan's until dec-B437.
- Recheck. stan4bart-creation-mapping lands first and this slice stops if it has not (item 9), and
  stan4bart's install lifts the weights where fewer than two are positive, so a supplied range is still
  used there (dec-B445) and no user sees a change. Not taken: leaving its read-back's "internal error"
  for a legal call; and a refusal in plain words, about ten lines fewer, which no ruling covers.
- Recheck. `dbartsSpec`'s creation warning passes through stan4bart unmuffled. Not taken: stan4bart
  muffling the class where it will install a range, for the one case (the rows in share a response
  value and not their parametric part) in which the scale the warning speaks of is then replaced.
- Recheck. The sd that stands in for an undefined regression is unweighted over the rows in, as rule 1
  counts a row and as dec-B433 words it ("sd(y) otherwise"). Not taken here: the weighted residual sd
  of the intercept alone (1.28 for 1.12 on one run under weights between 0.5 and 2), which is on the
  regression's scale where weights are far from 1 and which no ruling of this plan covers; it stays in
  TODO.
- Recheck. In rule 4's case the starting sigma signals no stand-in: the estimate is handed the data
  without its weights, so its verbose lines count every row. The alternative, the weighted call
  followed by a repair of its value, would announce a stand-in the fit does not use.
- Recheck. The state arm asserts the landed values of `setState` (TRUE unforced, NULL forced) and
  compares a gp forest by its leaf parameters; before the landing it asserted only a clean install.

## Recheck against bartcore

On 2026-10-11, against a2ef56a9: the notes of what landed since this plan was last revised were read
(setState's two forms and the state's record of missable columns,
[install-surface.md](../design/install-surface.md) and
[`SamplerStateData`](../../src/bartcore/sampler.hpp); the k-internal slice; the sparse starting sigma, its
[Landing note](sparse-starting-sigma.md#landing-note); stan4bart's
[stan4bart-creation-mapping.md](stan4bart-creation-mapping.md), with dec-A197), then the code, and the
runs below were made on a shipped build. The text above is written to that code throughout.

Runs: the planner's and the critics' behavior probes again (the range, the center, the unit, sigest, the
setters, no row in, one row in, fractional weights, the sparse sources, the by-hand route, `xbart`, the
hazard grid), every value as before and every starting-sigma warning gone; the two exact arms against
both oracles; the corpus scan and the traced suite of Context; a state across two ranges and two weight
sets, in both forms, for three leaf models; the stand-in under 0/1 and fractional weights; the nbinom
follow-up call; stan4bart's install, plain and with the weights lifted, emulated on samplers built by
hand.

Not verified: nothing of the slice is built, so no mutant and no new test was run. The numbers marked
(stands) were not rerun, and the prototype's suite run predates about 3,800 lines of tests/cpp, which
were not run against the rule. stan4bart was read on its creation-mapping branch, in build and not yet
landed, and not built or run: its install was emulated on dbarts samplers, and whether a fit on one
weighted row otherwise runs is not known. treatSens was read, not built; bartCause and bairrtt were
not read again. No reference build, sanitizer, x86 run or `R CMD check` was made.
