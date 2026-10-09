# response-scale-rows: what the response gives a prior is taken from the rows in the likelihood at creation

Status: PLANNED, revised from its critique 2026-10-09, then the same day to dec-B405 and onto the
bartcore of dec-A191's order. Implements dec-B296, dec-B302, dec-B309, dec-B314 (sigest held), dec-B315
(a hazard fit's periods), dec-B362 (several forests keep refusing), dec-B363 (slice B's half: no
documented by-hand call yet), dec-B364 (for `setResponse` and `setOffset`), dec-B405 (the keep is a
status: the engine returns it, the R interface warns, the flat entries return it), dec-B381 (the planning
calls approved) and dec-B386 (a constant response's window where every row is in); rests on dec-B357 (a
survival fit's range from observed times) and, for slice C, dec-B327; revises dec-A140. Nothing is held
for the maintainer. One open dependency and one premise for the orchestrator are recorded in
[Open items](#open-items). This file plans slice B, creation. Slice A is
[aft-reanchor-observed-times.md](aft-reanchor-observed-times.md), landed. Slices C and D are stated at
the end and have no steps yet. dec-B303 stays in [forest-multiplier.md](forest-multiplier.md).

agent: one push. Opus implementer for the engine, the keep's status through the facade, the bridge and
the flat entries, dbarts.h, and their tests; sonnet for the R code (the hazard grid included), the
tinytest file, the gate scripts and the help once the engine is fixed; opus reviewer told to refute, and
to read the diff three times: once for bits that must not move where no row is out, once for which rows
each number reads, once for the ABI (the header, the token, the C API test). stan4bart's part: a sonnet
implementer on its bartcore branch, reviewed by the same opus reviewer, landing in lockstep. Tier:
changes draws ([Process by risk](README.md#process-by-risk)). Its blind critique ran at 8ad1fa58; the
sections that dec-B405 and the rebase add (the keep's status, Interactions, stan4bart's part) take a
blind critique scoped to them before the build.
rng: by call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines
them. A row is IN when its weight is positive and no mask switches it off, OUT otherwise; a hazard
subject is in when its weight is positive.
- NEUTRAL, bit for bit: every sampler created with every row in, whatever it does afterwards (a mask or
  zero weights installed after creation included); every probit, logistic, ordinal and multinomial
  sampler; a hazard fit with every subject in, with stated boundaries, or on the default grid whose
  subjects at weight 0 hold only times some subject in holds; a single-forest gaussian or Student-t
  sampler created with rows out whose smallest and largest response less offset are both on rows in.
  The keep (a re-derivation in rule 4's case) is the `updateScale = FALSE` call, bit for bit; its status,
  its warning and the flat entries' new return value move no draw.
- POSTERIOR-CHANGING, a default changes, for a sampler created with rows out: a gaussian or Student-t fit
  with an extreme on a row out (the range); a gaussian model of several forests with any row out (the
  unit, so every forest's default sd, and what a stated `sd` means until
  [forest-sd-unit.md](forest-sd-unit.md) lands); a count fit under 0/1 weights (the center); `xbart` with
  a weight of 0 (each fold's range); any of them where sigest takes its fallback; a hazard fit with a
  subject at weight 0 whose time no subject in holds, or under `breaks = K` with any subject at weight 0
  (the periods, and the rows of the subjects in). And `setResponse` or `setOffset` with
  `updateScale = TRUE` while rows are out, in every family that has a range or a center, through the R
  methods and the two flat C entries: the rows in, or, where they hold fewer than two distinct values,
  the transform kept in force where the tip reads every row. A stan4bart gaussian fit with a zero weight
  is in this class.
- Moves one recorded baseline scenario, `zeroweights` of the main corpus; nothing else recorded moves
  (the main corpus's `hazard` scenario carries no weight, read).
window: before the merge to main: it changes draws against 0.9-34. dec-A191's engine order puts it after
the k-internal slice (k-internal-parameter.md, branch wt/k-internal-plan) and setstate-force-update Parts A
and B, with default-rule-per-gap landing ahead of both; and after sparse-starting-sigma, which is R only
and not queued behind the engine ([Interactions](#interactions)). Before
[forest-sd-unit.md](forest-sd-unit.md), so the unit's tests are written once against the rows that stay.
Serial with any other work in model.hpp, the response setters of
[`Chain`](../../src/bartcore/chain.hpp) and [`Sampler`](../../src/bartcore/sampler.hpp), the K-forest
constructor of [`Chain`](../../src/bartcore/chain.hpp), [`resolveSamplerSpec`](../../R/spec.R),
[`estimateSigmaFromLinearModel`](../../R/utility.R), [`expandDiscreteTimeHazard`](../../R/dbarts.R) or
dbarts.h. stan4bart's part lands with this slice's header commit, on its bartcore branch.
budget: ~2050 lines (engine ~180, of which the status ~45; facade, bridge, flat entries, dbarts.h and the
token ~75; tests/cpp ~340; R ~165, of which the hazard grid ~30; tinytest ~580, of which ~540 a new file;
the C API test ~80; benchmarks and the workflow ~250; help, header comments and NEWS ~110; design note,
indexes, amended notes and TODO ~260), upper figure 3700 ([Budget and stops](#budget-and-stops)).
stan4bart's part, ~60 lines, is outside it.

## Goal

A fit created with rows out of the likelihood has the leaf prior, the prior center, the default forest
sds and, for a hazard fit, the periods of the fit on the remaining rows: the response of a row out
reaches no number and no draw at the rows in. Those numbers are then held, as they are today. A fit
created with every row in is the fit it is today, bit for bit. A re-derivation that finds too little to
read keeps the scale in force, and says so: as a typed status from the engine, as one warning per
sampler from the R interface, and as a distinct return value from the flat C entries.

## Context

Marks: (ran) by the planner or the reviser on the tip's build, shipped mode, and on the design's
prototype library where a line says so; (stands) measured by the design, its critic or the plan's
critic and not rerun; (read) in the source. One fixture unless a line says otherwise: 200 rows, two
predictors, rows 101 to 200 at weight 0 with their response 50 higher. The measurements were taken at
8ad1fa58; the rebase onto a8884ef4 moved none of the code they read (read).

What reads which rows today.
- (ran) The range ([`GaussianResponse::rescale`](../../src/bartcore/model.hpp)): every row. Recorded
  [-2.165, 56.371] where the rows in span [-2.165, 6.029] and `subset` gives the latter; `k.scale` 29.27
  against 4.10 on the prototype. Student-t and aft hold a Gaussian response and inherit it (read).
- (ran) The count center ([`NBResponse::computeShift`](../../src/bartcore/model.hpp)): every row. Half
  the rows structural zeros under 0/1 weights: recorded 0.2964, the rows in give 0.9895, the prototype
  0.9895. The 0/1 weights become the mask and the weights slot is empty.
- (ran) The unit of a forest's sd in a model of several forests
  ([`scaledResponseSd`](../../src/bartcore/chain.hpp)): every row. 0.4313 of the range, 25.25 in y,
  with or without the weights; 1.466 in y under `subset` and on the prototype. Held across a later
  `setWeights` on both builds.
- (ran) sigest ([`estimateSigmaFromLinearModel`](../../R/utility.R)) is already the weighted fit over
  the rows in (1.0741, `lm`'s). Its fallback and its floor are not: with three rows in and two predictors
  it is 25.245, the sd of every row, where the rows in give 0.890
  ([`floorSigmaEstimate`](../../R/utility.R)); with four rows in at one value it is 8.4e-07, the floor
  from the largest residual of every row ([`floorMarginalSigma`](../../R/utility.R)).
- (ran) sigest's sparse paths (dec-B370), as on the tip. A caller's sparse x (a `dgCMatrix`) takes the
  marginal fallback with its own warning: 25.245, the sd of every row, where the rows in give 1.466. An
  indicator expansion R built sparse goes through [`sparseResidualStandardError`](../../R/utility.R),
  already over positive weights (read). (read) sparse-starting-sigma (branch wt/sparse-sigma-plan,
  d1be8904) replaces both: one routine for every sparse source, which drops the rows of weight 0 as
  `lm.wfit` does, the caller's-sparse fallback and its warning class retired. After it lands, the sparse
  paths read the rows in already, and only rule 4's case and the floor are this slice's.
- (ran) Held already: `setWeights` and `setActiveRows` move nothing. `setResponse` and `setOffset` with
  `updateScale = TRUE` under zero weights or a mask read every row ([-2.165, 56.371], with three rows in
  at one value as well); on the prototype the rows in. On several forests `updateScale = TRUE` is
  refused ("every forest keeps its leaf calibration stated against the scale fixed at creation").
- (stands) A re-derivation on the prototype, which reads the rows in and falls back to every row: a
  component created over its rows [-2.06, 2.02], left with one row and re-derived, jumps silently to
  every row's [-2.06, 22.29]; left with 25 rows, [-1.91, 1.51]. That jump is what dec-B364 turned down.
- (ran) No row in. Every weight 0: the sampler is created, with one warning that names another cause
  ("starting sigma estimate falls back to the marginal response sd"), the range every row's, sigest 25.245
  (the marginal sd; `lm` without the weights gives 25.006). A count fit under a mask of zeros is created
  in silence.
- (ran) One row in. The tip reads every row. The prototype gives the range [1.666, 1.666] and `k.scale`
  0.5, and a count fit with one row in at zero the center log(1/2). (stands) Then with every row brought
  back in and 600 sweeps the fit's rmse is 2.41 and sigma 2.66, where the fallback gives 0.46 and 0.90.
- (ran) Fractional weights with a zero among them raise the data object's validity warning ("rows with
  'weights' of 0 are left out of the likelihood", [`dbartsData`](../../R/A_class.R)) at creation: one
  row in at weight 0.5 gives two warnings today, that one and the sigma fallback's; half the rows out
  at fractional weights gives that one alone. 0/1 weights raise none.
- (ran) `xbart` with a weight of 0: the response of the rows out reaches the result (1.38, 1.23 against
  1.15, 1.20 with those responses ordinary); on the prototype the two are identical, the engine edit
  alone closing it. `xbart` never calls `resolveSamplerSpec`: it calls
  [`enforceWeightPolicy`](../../R/spec.R) and [`estimateStartingSigma`](../../R/spec.R) itself and
  chooses each fold's sigma (`sigmaPerFold`, [`xbart`](../../R/xbart.R)) by catching
  `dbartsSigmaFallbackWarning` (read). With every weight 0, and with one row in, it stops in its loss
  ("loss function must return non-missing numeric values") after the sigma fallback's warning; with 60
  rows in at one value it runs, silent today. Only several rows in at one value reach a result.
- (ran) A hazard fit's periods: 40 subjects on times 1 to 5, two at weight 0 holding 3.5 and 20: 7
  periods and 167 person-period rows, where the fit without the two has 5 and 134 and the two moved onto
  the grid (3 and 5) give 5 and 142. `breaks = 3`: periods 3, 4, 20 against 3, 4, 5 without the two.
  Every weight 0: created in silence on every subject's 7 periods. `hazard.logistic` refuses a zero
  weight, so only the probit link reaches the case (dec-B315 says the same). A masked probit row's
  response reaches no draw at the rows in: flipping it leaves the rows in `identical()` (3.3e-15 at the
  rows out). The periods are taken in [`resolveHazardGrid`](../../R/dbarts.R) from `gridTime`, which
  already narrows to the subjects the na.action keeps, a later time placed in the last period (read).
- (ran) `rbart_vi` is deprecated on the tip ("is deprecated and is removed in dbarts 1.1-0"), so its
  `rel.scale` is not in this slice.

What the setters return today (read).
- [`SamplerBase`](../../src/bartcore/facade.hpp)'s `setResponse` and `setOffset` return nothing, and so
  do [`Sampler::setResponse`](../../src/bartcore/sampler.hpp),
  [`Chain::setResponse`](../../src/bartcore/chain.hpp) and the response models' virtuals. A chain's
  `updateScale` arm always runs its tail: [`reanchorVarianceForest`](../../src/bartcore/chain.hpp) under
  a variance forest and [`restateSavedDraws`](../../src/bartcore/chain.hpp); the sampler then runs
  [`Sampler::recordAnchor`](../../src/bartcore/sampler.hpp). The engine's one typed outcome of a setter
  is [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp) (accepted, rolled back, invalid cut
  points), which the bridge converts.
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
  (["expect_equal(versions, c(1L, 0L))"](../../inst/tinytest/test-capi.R)).
- Consumers (read). stan4bart (bartcore 963956b) builds its sampler through `dbartsSpec` and
  `new("dbartsSampler", ...)` in its parent R process, then in each chain's worker (PSOCK or FORK
  processes) calls `dbarts_sampler_setOffset` with `updateScale` true once after creation and on a
  halving schedule through warm-up, and false after; it ignores the return. treatSens (dbarts-1.0) calls
  it on a probit sampler and ignores the return. bartCause and bairrtt make no flat call.

What a correct change must not move.
- (ran) The prototype's rewrite of the unit's loop changed its last bit on unweighted data: at n = 64 the
  tip gives 0x1.c828f537ceda5p-3 and the prototype ...da4p-3. (stands) On a reference build that failed 3
  of the 15 two-forest scenarios bitwise at max abs z 0.00; the suite did not see it.
- (stands; the plan's critic re-ported the prototype onto 8ad1fa58) On a reference build of the
  prototype: main corpus against equivalence-51201107 54 of 55 identical and `zeroweights` at max abs z
  2.85, 0 skipped; two-forest 15 of 15; multinomial 11 of 11; the four snapshot files pass; the exact
  gates and the zero-weight arm pass in `quick`. `zeroweights` is the one recorded scenario that holds a
  zero weight when its sampler is created (80 of 500, the smallest response on one of them:
  [1.3742, 25.2117] becomes [1.6243, 25.2117]). The hazard grid was not in the prototype; the `hazard`
  scenario carries no weight (read). The main corpus's current file has moved since (equivalence-734441f1
  on bartcore; default-rule-per-gap and sparse-starting-sigma each re-record other scenarios, neither
  `zeroweights`, read): the compare is against whichever file is current at this slice's base.
- (stands) The suite on the re-ported prototype: 2 failures of 19506 results, both the two-forest pin
  of test-active-rows-pins.R (["Bayesian causal forest, bitwise"](../../inst/tinytest/test-active-rows-pins.R)),
  which compares a sampler created with weights w and then masked against one created with w * a.
- (read) No test hazard fit passes a weight of 0 (the hazard files pass weights of 1, refused
  fractional ones, or none).

Gates that can fail.
- (ran) [`zeroWeightArm`](../../benchmarks/R/bd-balance.R) creates with every row and so proves nothing
  here. With the weights given at creation and the zeroed cells' response 3 higher, about 1 second in
  `quick`: the tip matches the oracle whose range is over every row (largest abs z 1.4) and misses the
  rows-in oracle by up to 106.6, at the (none) partition; (stands) the prototype matches the rows-in
  oracle (2.1) and misses the other by 40.1. bd-balance.R holds three arms already, `zeroWeightArm`,
  `narrowArm` and `constMissingArm`, behind one `stopifnot` (read).
- (ran) A count twin, from the fixed-shape arm of negbin-exact.R: 4 rows in of 25 a cell, the rows out 30
  higher, left out by 0/1 weights at creation, the oracle's likelihood over the rows in and its center
  over every row (3.334) or the rows in (1.056). 2 seconds for two seeds. Tip: gap 0.012 to the every-row
  oracle, 0.493 to the rows-in one, tolerance 0.12. Prototype: 0.484 and 0.003.

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
   `setWeights`, `setActiveRows`, `setState`, a copy, a reload or a warm start does. No argument is added
   in this slice.
3. Re-derived. `updateScale = TRUE` on `setResponse` and `setOffset` reads the rows in at that call.
   `setData` replaces the data and derives as creation does, rules 1 and 4, its warning included.
   Several forests still refuse `updateScale = TRUE` (dec-B362).
4. Too little to read. Where a row is out and the rows in hold fewer than two distinct values (none in;
   one in; several in at one value):
   - At creation and at `setData`, each number is taken over every row, as if no weight or mask had
     been given, sigest included where one is estimated (dec-B302, dec-B309). One warning, of its own
     class, and the sigma fallback's warning is not raised for that cause:

         fewer than two distinct values of the response are in the likelihood (<m> of <n> observations are in it), so the response's scale and the prior's center are taken from all <n>

   - At a re-derivation through `setResponse` or `setOffset`, the transform and the center in force are
     kept: the call completes as it does with `updateScale = FALSE`, bit for bit (dec-B364). The engine
     returns that outcome as a status and raises nothing; the R methods warn once per sampler, with the
     same class; the flat entries return a value of their own for it and warn nothing (dec-B405):

         fewer than two distinct values of the response are in the likelihood, so the response's scale and the prior's center are kept as they were (said once per sampler)

   The value tested is the response less its offset for a range and the count for a center. With every
   row in, a constant response is not this case: it takes dec-B386's window c - 0.5 to c + 0.5 at
   creation and at a re-derivation, with no warning of this class and no kept status, so dec-B364 and
   dec-B386 never apply to the same call.
5. Nothing out, nothing changed. With every row in, every number above has the bits it has today.
6. sigest stays `lm`'s weighted estimate over the rows in and is re-estimated by no call (dec-B314). Its
   fallback sd and its floor read the rows in, on each of its paths: the dense linear fit and the sparse
   routine sparse-starting-sigma leaves. Under a fixed residual prior no estimate is made (that plan's
   Change 4), and nothing here applies.
7. `xbart`: each fold's range and sigma follow 1, 4 and 6 over the fold's training rows.
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
- The keep is the `updateScale = FALSE` call, bit for bit: on it no variance-forest re-anchor, no
  restatement of kept draws, no named-sd translation (the k-internal slice's routine, which runs only
  where the mapping moves) and no `recordAnchor`.
- Rule 5 is a constraint on the code, not only a test: the unit's loop and the center's loop run as they
  are written today whenever no row is out, and the rows-out case is a separate path. A reviewer checks
  that neither present loop gained a branch.
- One decision of which rows are read at a re-derivation. The response model makes it and returns it;
  the unit follows the response model's answer, and the bridge and R read the status, none of them
  testing the rows again. At creation and `setData` R also tests rule 4, for sigest and the creation
  warning, as before; step 8's literals hold the two in agreement.
- dbarts.h keeps every signature. It gains one enumerator list, folded into the token, so
  `DBARTS_C_API_HASH` moves; `DBARTS_C_API_MINOR` goes from 0 to 1 (dec-B405; see
  [Open items](#open-items) for dec-B60); the comments of the two re-anchoring entries gain "over the rows
  in the likelihood", the keep and the return values. An ABI event: every LinkingTo consumer is rebuilt
  in lockstep.
- No state format change and no new record: the range is the model's `response.range`, the unit the
  forests' `anchor`, sigest the data object's, the periods the fit's `$periods`, as now. A re-creation
  restates the record and derives nothing.
- Cut points, a hazard fit's stated boundaries, the number of categories, leaf covariate
  standardization and a multiplier's scale are not touched.
- [`rbart_vi`](../../R/rbart.R) is not touched.

## The keep's status (dec-B405)

How the keep travels from the response model to each caller. Step 4 builds it; steps 1 and 2 produce it.

1. Engine type. `ResponseUpdateResult { accepted, scaleKept }`, an `enum class` beside the
   `ResponseModel` interface in model.hpp, on the pattern of
   [`PredictorUpdateResult`](../../src/bartcore/sampler.hpp). `accepted` is every call that is not a
   keep, `updateScale` false included. The response models' `setResponse` and `setOffset` return it:
   probit, ordinal, logistic and multinomial always `accepted`; `GaussianResponse` and `NBResponse`
   `scaleKept` in rule 4's case at `updateScale` true, having run their false branch; `TResponse` and
   `AFTResponse` their Gaussian part's answer, one line each. `setData` returns nothing (creation's
   fallback, not the keep).
2. Chain and sampler. [`Chain::setResponse`](../../src/bartcore/chain.hpp) and
   [`Chain::setOffset`](../../src/bartcore/chain.hpp) return the model's answer and skip their
   re-anchor tail on `scaleKept`. [`Sampler::setResponse`](../../src/bartcore/sampler.hpp) and
   [`Sampler::setOffset`](../../src/bartcore/sampler.hpp) return it, and skip `recordAnchor` on it. The
   chains share the response, the offset, the weights and the mask, so they agree; the sampler returns
   chain 0's answer and a debug assertion checks the rest. Under per-chain responses
   ([5. Recommendation](../design/per-chain-responses.md#5-recommendation)) the pooled answer is that
   design's to state. [`SamplerFacade`](../../src/bartcore/facade.hpp) forwards.
3. Bridge. `bartcore_setResponse` and `bartcore_setOffset` return a length-one logical: TRUE when the
   engine returned `scaleKept` and the holder's new `scaleKeptWarned` was clear, which the call then
   sets; FALSE otherwise. `scaleKeptWarned` goes after `slowCountWarned`, last, by the holder's
   positional-aggregate rule. It lives as long as the engine handle, so a copy, a reload or any other
   re-creation is a new sampler and warns again on its first keep.
4. R warning. `bartcoreSamplerSetResponse` and `bartcoreSamplerSetOffset`
   ([bartcore.R](../../R/bartcore.R)) raise rule 4's keep text through [`warnClassed`](../../R/utility.R)
   with class `dbartsScaleFallbackWarning` when the bridge returns TRUE, after the data mirror is
   updated, so a warning turned into an error leaves the R object agreeing with the engine.
5. Flat entries. [`dbarts_sampler_setResponse`](../../src/C_interface.cpp) and
   [`dbarts_sampler_setOffset`](../../src/C_interface.cpp) return `DBARTS_SWAP_SCALE_KEPT` (2) when the
   engine returned `scaleKept`, `DBARTS_SWAP_DONE` (1) on any other swap and `DBARTS_SWAP_NONE` (0)
   where the coupling admits no conduit, as today. They never read or set `scaleKeptWarned` and raise no
   warning: a flat keep does not spend the R interface's once.
6. Contract. dbarts.h gains `DBARTS_SWAP_STATUS_LIST` and `typedef enum { ... } dbarts_swap_status`,
   generated by the list macro as the other three are, so the token folds it. The return type of both
   entries stays `int`. Its comment: 0 and 1 as before; 2 a swap that asked for `updateScale` and kept
   the transform and the count center in force because the rows in the likelihood (positive weight,
   not masked) hold fewer than two distinct values of the response less the offset, or of the count
   for nbinom; the swap itself done, as with `updateScale` false; a caller that tests whether the swap
   happened tests for nonzero, not for 1. The two entries' comments: `updateScale` reads the rows in.
   "All three ABI enums" becomes four in the header and in the fold's comment
   ([`dbarts_apiToken`](../../src/C_interface.cpp)); `DBARTS_C_API_HASH` is re-baked;
   `DBARTS_C_API_MINOR` is 1 and the version paragraph names this bump.
7. The C API test. [consumer.c](../../inst/tinytest/capi/consumer.c) exports the three constants;
   [test-capi.R](../../inst/tinytest/test-capi.R): the constants are 0, 1, 2; the versions `c(1L, 1L)`,
   their comment rewritten; the hash its new literal. A gaussian sampler with one row in at positive
   weight: each entry with `updateScale` true returns 2 and raises no condition (counted with
   `withCallingHandlers`), the transform and five sweeps bit for bit an `updateScale`-false twin's; then
   `$setOffset(updateScale = TRUE)` on the same sampler warns once (the flat call left the flag), and a
   second call none. Two rows in at two values: 1 and the rows-in transform. `updateScale` false: 1.
   nbinom under a mask keeping one count: 2. probit: 1. The stale-minor and wrong-hash handshake arms as
   they are; one new arm compiles the consumer with `-DDBARTS_C_API_MINOR=0` and without
   `DBARTS_REQUIRE_EXACT_ABI`, a consumer built before the bump, and shows it passes the handshake and
   reads 2 as nonzero (the append-compatibility claim).
8. stan4bart's part, in lockstep on its bartcore branch. What it does with the value:
   - At each call with `updateScale` true (src/init.cpp: the call after creation and the warm-up loop's
     scheduled calls), compare the return with `DBARTS_SWAP_SCALE_KEPT` and count it on its sampler
     struct. The chain's results carry the count out of the worker, since a warning raised inside a
     PSOCK worker does not reach the user reliably; the parent, after the chains return, raises one
     warning per fit when any chain kept: "the BART component's response scale was kept at <k> of <r>
     re-derivations in <c> of <C> chains: fewer than two distinct values of the response less the offset
     are in the likelihood". Whether it warns or only records is stan4bart's call (dec-B405, "stan4bart
     deciding what its users see"); the plan recommends the warning, since dbarts's R methods warn for
     the same event.
   - When it fires (read). stan4bart's weights are fixed for a fit and its gaussian fits take no mask,
     so with two rows in at distinct responses a keep needs every response less offset among the rows in
     to be equal, which a continuous offset does not give. In practice it fires with one row in, where
     `dbartsSpec` has already raised dbarts's creation warning once in the parent. The pair is the usual
     case, and the test pins both.
   - Size: C++ ~15 (two read sites, the count, its export), R ~15, one test ~30 (one row at positive
     weight, warm-up 10, two chains on one core: dbarts's creation warning and stan4bart's one warning,
     each pinned by class or pattern; two rows in at distinct values: neither), ~60 lines, and the
     rebuild against the new header.
   - Either way of the open dependency ([Open items](#open-items)). If stan4bart stops re-deriving
     during warm-up, the loop's read site goes and the call after creation is the only one: the same
     count and warning, worded for that call, ~45 lines. If it also stops calling with `updateScale`
     true after creation (its sampler created at the initial fit less its intercept), nothing asks to
     re-derive: its part is the rebuild against the new header, with no read, and dbarts's creation
     warning is its only message. Nothing on dbarts's side depends on how often stan4bart asks.
   - treatSens and bartCause change no source: treatSens's probit sampler always answers 1, and
     bartCause makes no flat call. Their lockstep CI rebuilds against the new hash.
9. Gates of this part: tests/cpp at each level (step 4); the C API test; step 8's keep arm through the
   R methods; stan4bart's suite and its new test against the slice's library; and once, reported, a
   stan4bart binary built against the base header run against the slice's library, which must pass the
   handshake and fit (it ignores the return).

## Steps

"Fails today" is what the tip does where the test expects otherwise.

1. Engine, the range. [`GaussianResponse::rescale`](../../src/bartcore/model.hpp) reads the rows whose
   served precision is positive, from the response or from the vector slice A has it take. The read
   (the range over the rows in, and whether they hold two distinct values) is one routine, separate from
   the write, so a later pooled re-anchor can call it per chain. In rule 4's case the constructor and
   `setData` read every row, and `setResponse(true)` and `setOffset(true)` run their false branch and
   return `scaleKept`. Student-t (its composite precision is zero exactly where the weight or the mask
   is) and aft (the mask only) follow with no edit but passing the status on. The response model
   answers which rows it read, for step 3. tests/cpp, beside
   [`testActiveRowsGaussianDf`](../../tests/cpp/test_model.cpp), two chains where a chain is involved:
   - gaussian, Student-t and aft (aft by a mask then `setOffset(true)`), the weights of the rows in
     fractional: the scale and shift equal the unweighted literal over the rows in (fails today: every
     row);
   - no weight at all, and weights all positive: scale, shift and the working response are bit for bit
     those of the present expression, written out in the test;
   - created with none in, one in, four in at one value: the every-row literal; `setData` to those rows:
     the same;
   - re-derived with a row out and none in, one in, four in at one value: the transform, the working
     response and 20 sweeps bit for bit those of a twin given the same call with `updateScale` false, and
     `scaleKept` returned (fails today: every row); the same under a variance forest, with kept draws,
     and with a named sd, each bit for bit its twin's; with every row in and a constant response:
     dec-B386's window and `accepted`;
   - held: after `setWeights` and `setActiveRows` the transform is unchanged and 20 sweeps equal an
     untouched twin's under the same weights;
   - `setResponse(true)` and `setOffset(true)` under zero weights and under a mask read the rows in and
     return `accepted`; `setData` reads the rows in under the new weights.
2. Engine, the center. [`NBResponse::computeShift`](../../src/bartcore/model.hpp) takes the mask where
   one is installed; the unmasked call is the present loop. Its two re-deriving callers pass the mask,
   and where the rows the mask keeps hold fewer than two distinct counts they keep the center in force
   and return `scaleKept`. tests/cpp, beside [`testNBLogMeanAnchor`](../../tests/cpp/test_model.cpp): the
   center after a mask and `setOffset(true)` equals the literal over the rows the mask keeps, with and
   without an exposure offset; unmasked, bit for bit the present value; a mask keeping one row, or rows of
   one count, then `setOffset(true)`: the center in force, bit for bit the `updateScale`-false twin's, and
   `scaleKept` (fails today: every row); held across `setActiveRows`. The comments of
   [`NBResponse::setActiveRows`](../../src/bartcore/model.hpp) and
   [`AFTResponse::setActiveRows`](../../src/bartcore/model.hpp) stop saying "full-data".
3. Engine, the unit. [`scaledResponseSd`](../../src/bartcore/chain.hpp) is left as it is and is what
   runs when the response model read every row; otherwise a second routine takes the sample sd over the
   rows it read, n - 1 of those rows in the divisor. tests/cpp, beside
   [`testForestCalibration`](../../tests/cpp/test_sampler.cpp): with rows out the map's anchor equals the
   literal sd to 1e-14; with none out it is bit for bit the present loop's value at n = 64, 100, 137 and
   250 (the size at which the prototype moved a bit among them); a finite carried anchor still wins.
4. The keep's status, engine to C API, as [The keep's status (dec-B405)](#the-keeps-status-dec-b405)
   items 1 to 7: the chain and sampler returns and their skipped tails, the facade, the bridge's logical
   and holder flag, the R warning, the flat entries' values, the header's list, comments, minor and hash,
   and the token's fold. `--preclean`. tests/cpp: the status at the chain and the sampler on two chains,
   which agree; the facade forwards it ([`SpySampler`](../../tests/cpp/test_facade.cpp) takes the new
   return type); after a keep `getAnchor` is unchanged. The C API test as item 7.
5. R, the rows and the warnings. One function in [spec.R](../../R/spec.R) answers rule 4 for a data
   object and a mask and raises the creation warning, for the families that take something from the
   response (gaussian with any residual distribution, nbinom) and no other. Its callers:
   [`resolveSamplerSpec`](../../R/spec.R), after [`enforceWeightPolicy`](../../R/spec.R) (so `dbartsSpec`
   and every door that resolves a spec, stan4bart's included);
   [`bartcoreSamplerSetData`](../../R/bartcore.R), for `setData`; and [`xbart`](../../R/xbart.R), once,
   after its own `enforceWeightPolicy` call and before the starting-sigma estimate whose warning chooses
   `sigmaPerFold`. It is not called at a re-derivation: there the engine's status answers (step 4).
   - sigest, written against sparse-starting-sigma's routine.
     [`estimateSigmaFromLinearModel`](../../R/utility.R), in rule 4's case, fits without the weights on
     the dense path and through the sparse routine; otherwise both keep the weights, which drops the rows
     out, and [`floorSigmaEstimate`](../../R/utility.R)'s fallback and
     [`floorMarginalSigma`](../../R/utility.R) read the rows in. The fallback's warning is not raised for
     rule 4's cause.
   - `xbart`. In rule 4's case the all-rows estimate is the unweighted fit, which raises no fallback
     warning, so `sigmaPerFold` is "linear". Each fold tests rule 4 over its training rows with no
     warning: in its case the fold's fit is unweighted; otherwise the weighted fit, and the fallback and
     the floor over the training rows in. So `xbart` warns once and a fold falls back in silence.
6. R, creation under a mask. Of the families whose 0/1 weights become the mask
   ([`isMaskedWeightFamily`](../../R/spec.R)) only nbinom takes anything from the response. For it,
   where rule 4's case does not hold, [`dbarts`](../../R/dbarts.R) and
   [`bartcoreSamplerSetData`](../../R/bartcore.R) follow their install of the mask with the sampler's own
   `setOffset(offset in force, updateScale = TRUE, updateState = FALSE)`, before the sampler is handed
   back: it draws nothing and records the range; a named sd is retranslated by the engine at that call
   (the k-internal slice; no R write-back remains). In rule 4's case the follow-up is skipped and the
   every-row center of creation stands, under step 5's warning, so no keep can arise there. A sampler
   built by hand from [`dbartsSpec`](../../R/spec.R) and `new("dbartsSampler", ...)` that installs
   `spec$active` itself is created with every row in and, for a count fit, keeps the every-row center;
   in this slice its help says so and gives no call. The documented call,
   `setActiveRows(spec$active, updateScale = TRUE)` with its loop form, comes with slice C, which builds
   it (dec-B363); the `setOffset` route is not documented.
7. R, the hazard periods. [`expandDiscreteTimeHazard`](../../R/dbarts.R)'s two callers in
   [`dbarts`](../../R/dbarts.R), the direct response and the formula's, pass as `gridTime` the times of
   the subjects kept by the na.action whose weight is positive, all of them where no weight is given or
   none is positive (rule 8's fallback, raising its warning there). [`resolveHazardGrid`](../../R/dbarts.R)
   takes the default grid and the `breaks = K` quantiles from it unchanged; stated boundaries keep their
   check over every subject. A subject whose time is past the last period is expanded to every period
   with a response of 0. Nothing in the engine.
8. tinytest, a new file test-response-scale-rows.R, fixtures with the rows out extreme in the response:
   - Inert. Two fits that differ only in the response of the rows out at creation are `identical()` at
     the rows in and within 1e-12 at the rows out, in draws, sigma and every reader: gaussian, a fixed
     sigma, Student-t, a variance forest, linear and gp leaves (covariates equal across the pair), a
     named sd, two forests on a factor and on a number, an offset, nbinom, and `bart`, `bartBT` and
     `xbart` once each. Fails today on every arm.
   - Literals, with 0/1 weights. `response.shift`, `response.scale`, `k.scale`, `prior.mean`, the
     forests' `anchor` and sigest equal values written from the rows in, and equal the `subset` fit's.
     `k.scale` is the data's scale over the rows in under the k and the sd spelling alike (the k-internal
     slice), and an sd-spelled fit's k is `k.scale` over the named sd. With fractional weights in, sigest
     is the weighted fit and the `subset` fit's is not, so that arm asserts the range, the center and the
     unit only.
   - Held. The readers are unchanged across `setWeights`, `setActiveRows`, a copy and a reload, read from
     the engine (`getLeafPrior`) and not from the record alone; five sweeps after a copy agree with the
     original's to 1e-10.
   - `updateScale = TRUE` on `setResponse` and `setOffset` under zero weights and under a mask gives the
     readers of a sampler created under those rows, and no warning; several forests still refuse.
   - The keep (dec-B364, dec-B405). Gaussian under zero weights and nbinom under a mask, re-derived
     through each of `setResponse` and `setOffset` with none in, one in, and several in at one value: the
     readers and five sweeps equal an `updateScale = FALSE` twin's; exactly one warning of the new class
     on the first such call, its text, counted with `withCallingHandlers`, and none on a second through
     either setter; a copy warns again on its own first keep; a constant response with every row in:
     dec-B386's window and no warning. The flat entries' side is in the C API test (step 4).
   - Rule 4 at creation. Every weight 0; a mask of zeros (nbinom); one row in; rows in at one value:
     created, exactly one warning of the new class, counted with `withCallingHandlers`, its text as above
     with m and n; the readers and sigest equal those of the fit with no weights; no warning where every
     row is in and the response is constant; none for probit under a mask of zeros. Under fractional
     weights with a zero among them the data object's validity warning is raised as well, as today, and
     pinned by its own pattern exactly where it fires. `setData` to such rows: every row, the same one
     warning. `dbartsSpec` with one row in: the one warning, from the spec.
   - sigest: three rows in and two predictors gives the sd of the rows in with the existing fallback
     warning; the floor reads the rows in. Sparse: a caller's `dgCMatrix` x and an indicator expansion
     (`factors = "indicators"`) under the fixture give the dense design's sigest over the rows in; with
     three rows in, the fallback over the rows in and its one warning; in rule 4's case both give the
     unweighted fit over every row and the one new warning.
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
   - A state across row sets: a single-forest, a gp and a two-forest state from a sampler created under
     other zero weights install as stored on the internal scale (dec-B418), with no conversion and no
     refusal; the recipient's range is unmoved and its internal leaves are the donor's bit for bit. Its
     report is whatever setstate-force-update landed; a cross-mapping install is clean (read: the
     k-internal plan).
9. Restated pins. In test-active-rows-pins.R the two-forest pin creates both samplers under w and masks
   or reweights after; the gaussian and Student-t pins above it get the same repair, since they pass only
   while the fixture's extremes are on rows in. The comments
   ["The log-mean shift is the full-data one by design"](../../inst/tinytest/test-active-rows-pins.R) and
   ["the response transform is the FULL-data one by design"](../../inst/tinytest/test-active-rows-pins.R)
   say "creation's, every row then being in".
10. Exact gates. bd-balance.R gains a fourth arm beside [`zeroWeightArm`](../../benchmarks/R/bd-balance.R),
    `narrowArm` and `constMissingArm`: the weights at creation, the zeroed cells' response 3 higher, the
    oracle's range over the rows in; its `stopifnot` over the arms takes the fourth. negbin-exact.R gains
    a third arm beside its fixed-shape one, as measured in Context, with the oracle's center
    ([`shiftOf`](../../benchmarks/R/negbin-exact.R)) over the rows in.
    [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) runs the bd-balance arm in its
    `for arm in narrow constmissing` loop, whose comment widens; the negbin arm runs inside its script.
    About 3 seconds more in `quick`. benchmarks/README.md names both.
11. The baseline. `zeroweights` is recorded fresh on a reference build at the slice's newest code commit
    and the other 54 scenarios are carried from the main corpus's current file at the slice's base after
    they compare bitwise; the [MANIFEST](../../benchmarks/baselines/MANIFEST) row says so in its
    predecessors' form and names the oracle (step 12's pair script, rows 3 and 4). The workflows that
    pin the file's name follow.
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
    ways to leave a row out differ by when it left; the nbinom paragraph names the rows of c; the hazard
    paragraphs of both and of [dbartsFamilies.Rd](../../man/dbartsFamilies.Rd)
    (["sorted distinct observed times"](../../man/bart.Rd)) say the subjects in the likelihood;
    [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the response item
    (["which reads every row"](../../man/dbartsSampler-class.Rd) becomes the rows in), the weights item,
    the `active` item (["stay the FULL-data calibration"](../../man/dbartsSampler-class.Rd) becomes
    creation's), `updateScale` with the keep and its once-per-sampler warning, and that `setData` and a
    sampler made by hand from another's parts derive again; [dbartsSpec.Rd](../../man/dbartsSpec.Rd), the
    `active` value: a count sampler made by hand keeps the every-row center in this slice; xbart.Rd,
    `weights`. dbarts.h: the two entries' comments, the new list and the version paragraph (step 4).
    NEWS, user-visible changes, one item: "An observation of weight zero no longer takes part in the
    response's scale, which the leaf prior is stated against, so a fit with zero weights has the leaf
    prior of the fit on the remaining rows; 0.9-34 took the scale over every row. The observation is
    still in the design." The keep, its status, the flat entries' value and the hazard periods reached no
    release and take no item. A design note, docs/design/response-scale-rows.md, with its index row: the
    rule, the inventory, the measurements, the named-sd case (a prior centered at 52 where the rows in sit
    near 2), why membership and not weight, the keep and how its status travels (engine, R, flat
    entries), the hazard grid, the four slices. [active-rows-mask.md](../design/active-rows-mask.md),
    [negative-binomial.md](../design/negative-binomial.md), [nbinom-log-mean.md](nbinom-log-mean.md) and
    [survival.md](../design/survival.md) take a dated amendment where they say "full-data", "all rows" or
    the hazard grid's source. TODO: the `response-scale-rows` entry is cut to slices C and D;
    `hazard-periods-rows-in` is removed. The ledger entry for the calls made. The landing note here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, clean under ASan
  and UBSan; the full tinytest suite, expected clean after step 9 (2 failures without it); the new file
  and test-capi.R under ASan on the R-loaded path.
- On a reference build, the three compares, counted per scenario with no `max |z|` line but the one
  named:
  - main corpus against its current file at the slice's base: 54 identical; `zeroweights` differs and is
    the only one. Its statistical verdict against that file is reported, not gated: the posterior
    changed (the re-ported prototype's was max abs z 2.85). Then all 55 identical against the new file.
  - two-forest corpus against its current file (`bcf-equivalence-1b7d730c.rds` on bartcore; the
    k-internal slice expects no move): 15 of 15 identical. Fewer is rule 5 broken in the unit's loop.
  - multinomial corpus against `multinomial-equivalence-80b1c8d4.rds` or its successor: 11 of 11
    identical.
  - the four `test-reproducibility-*.R` files pass unchanged; none is regenerated.
- Every gate of [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) in `quick`, unchanged, and
  the two new arms. Each new arm is also run once on the base build, where it must FAIL (measured: by
  106.6 standard errors and by 0.49 against 0.12), and that run is reported.
- The pair script: rows 1 to 4 `identical()`, row 5 differing.
- The C API test passes on the slice's library; it fails on the base library (the constants list, the
  version and the hash, and 1 where 2 is expected).
- stan4bart, on the slice's library: its suite unchanged and its new test passing; one seeded gaussian
  fit with zero weights on the base and slice libraries, differing, beside one without, identical; the
  stale binary of item 9, passing.
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
  - the keep dropped (a re-derivation in rule 4's case reads every row): step 1's and step 2's
    re-derivation tests; step 8 "The keep"; the C API keep arms;
  - the keep running the chain's tail (the variance-forest re-anchor, the restatement or the named-sd
    translation): step 1's twin arms under a variance forest, with kept draws and with a named sd;
  - the status dropped (the keep returns `accepted`): step 1's and step 2's status checks; step 4's
    sampler check; step 8's warning count; the C API's 2;
  - the flat entries warning (reading or setting `scaleKeptWarned`): the C API arm counting no condition,
    and the R call after it warning once;
  - the keep's warning raised on every such call: step 8's "none on a second";
  - the creation fallback at `setData` replaced by the keep: step 1's `setData` test; step 8's `setData`
    arm;
  - creation's re-derivation under the mask dropped: step 8 "Creation routes"; the negbin arm; and run
    for probit too: the suite's seeded probit-under-mask pins, if the call draws;
  - sigest's fallback left on every row: step 8 "sigest", the sparse arm included; both warnings raised
    with no row in: step 8's count of one;
  - `xbart`'s folds warning each: step 8's `xbart` count of one;
  - the hazard grid over every subject: step 8 "Hazard", the periods and the `identical()` pair;
  - `updateScale = TRUE` left reading every row: step 1's last bullet; step 8's fourth.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; inst/NEWS.Rd parses with a non-NULL result; dbarts.h compiles as C (`gcc -std=c99 -pedantic`,
  as its own comment requires); `R CMD check --as-cran` on a tarball from a clean copy. No new Rd topic.
- Cost: creation and the two response setters gain a comparison per row; a hazard fit one subset of its
  times; no sweep path. stan4bart's warm-up calls gain that comparison at each scheduled re-derivation,
  which already passes over every row. No bench compare.

## Budget and stops

Planned ~2050 lines: engine ~180 (range ~60, center ~30, unit ~45, status ~45); facade, bridge, flat
entries, dbarts.h and the token ~75; tests/cpp ~340; R ~165; tinytest ~580; the C API test ~80
(consumer.c ~15, test-capi.R ~65); benchmarks and the workflow ~250; help, header comments and NEWS ~110;
docs ~260. Against the revision of faf8ca7d: +~45 engine and +~40 tests/cpp for the status, +~75 for the
ABI and +~80 for the C API test, where answer A had held ~40 of bridge; -~20 R and -~20 tinytest, the
caller's-sparse fallback gone with sparse-starting-sigma and the named-sd write-back with the k-internal
slice. Forecast: recent engine slices ran 1.8 to 2.3 times their plans (leaf-conversions 2.3, almost all
tests), so 3700 is the upper figure. stan4bart's part, ~60 lines, is outside it.

Stop and report when:
- the diff passes ~3700 lines;
- a recorded scenario other than `zeroweights` moves, the two-forest or the multinomial corpus is not
  bitwise, or a snapshot moves;
- the keep is not the `updateScale = FALSE` twin bit for bit on any arm (two chains, a variance forest,
  kept draws, a named sd);
- the chains of one sampler disagree on the status (they share all the decision reads; a disagreement is
  a defect, not a case);
- the status needs a facade virtual beyond the two return types, or dbarts.h needs more than the
  enumerator list, the minor, the hash and the comments (any signature change is a stop);
- stan4bart's part needs more than the read sites, the count's way out and one warning, or moves a
  stan4bart draw; or the stale stan4bart binary fails the handshake or a fit;
- a hazard fit with every subject in changes;
- sparse-starting-sigma has not landed when step 5 starts (step 5 is not written against the
  caller's-sparse fallback), or the k-internal slice or setstate-force-update Part B has not landed when
  the slice starts (the window).

Expected verdicts: tests/cpp, tinytest and the C API test clean; the main corpus 54 identical and
`zeroweights` moved, then 55 identical against its new file; two-forest 15 and multinomial 11 identical;
the snapshots unchanged; every exact gate passing and the two new arms passing on the slice and failing
on the base; the pair script's rows 1 to 4 identical and row 5 differing; stan4bart's suite clean.

## What a consumer sees

Read in each package's source (weights, `subset`, masks, `updateScale`, the flat setters' return); no
consumer suite was run.

- stan4bart. A gaussian fit hands its weights to the sampler at creation, zeros allowed (it refuses only
  all zero), through `dbartsSpec` and `new("dbartsSampler", ...)`, and then calls the flat `setOffset`
  entry with `updateScale` true once and through warm-up. Today both read every row. After this slice
  creation reads the rows in, and so must every one of those calls, or the first would undo creation's
  scale: that is why rule 3 is in this slice. A gaussian stan4bart fit with a zero weight therefore
  changes its leaf prior and draws; one without is bit for bit unchanged. Where the rows in hold fewer
  than two distinct values a call keeps the scale in force and returns 2, which stan4bart reads
  ([The keep's status (dec-B405)](#the-keeps-status-dec-b405), item 8); `dbartsSpec` has by then raised the
  creation warning once in stan4bart's parent process. Its binary fits install a mask after creation on
  a probit sampler, which takes nothing from the response and always answers 1. No test of its suite
  shows the change: test-15-weights_offset_k.R draws its weights from (0.5, 2) and from {4, 0.25}, and
  test-22-binary-weight-mask.R is probit. It rebuilds against the new header in lockstep.
- bartCause passes the user's weights and `subset` at creation, to `bcf` too, and changes neither
  afterwards. A fit with a zero weight moves as the `rng:` block says, under `bcf` by the unit as well.
  No test of its suite passes a zero weight (searched). No flat call.
- treatSens fits a probit treatment model and re-anchors it through the flat entry, ignoring the return;
  a probit sampler takes nothing from the response and answers 1. bairrtt uses no weights, no mask and
  no `updateScale`.

## Interactions

- The k-internal slice (k-internal-parameter.md, branch wt/k-internal-plan, bfe38871; lands before this
  under dec-A191). The engine holds a named sd and retranslates it wherever a chain's mapping moves,
  among them the `updateScale` arms of `Chain::setResponse` and `Chain::setOffset`, and
  `reissueNamedLeafSd` and its R callers go. Here: step 6's follow-up `setOffset` restates the named sd
  in the engine; the keep moves no mapping, so the translation does not run on it (that plan: "a kept
  scale moves no mapping"); step 8's `k.scale` literals read the data's scale over the rows in under
  both spellings. Installs convert nothing (dec-B418), so step 8's state arm expects no conversion and no
  gp or two-forest refusal, and Out of scope loses that bullet. This slice's edits sit beside that
  slice's translation call in the same two chain arms: a rebase, read for order.
- setstate-force-update (branch wt/install-surface-plan; Parts A and B land before this). Its plan
  records this slice as following it, meeting only through `setData`'s re-quantize: a rebase. Step 8
  asserts no `setState` report beyond a clean cross-mapping install.
- default-rule-per-gap (landing now). It re-records eight main-corpus scenarios, not `zeroweights`, and
  changes no response setter: the base's baseline file name is the only effect.
- sparse-starting-sigma (wt/sparse-sigma-plan, d1be8904; R only, may run beside an engine slice, before
  the merge). It lands first: this slice is fourth in dec-A191's engine queue, and it is unqueued.
  After it, every sparse source goes through one routine that already drops the rows of weight 0, the
  caller's-sparse marginal fallback and its warning class are retired, and a fixed residual prior makes
  no estimate. So this slice's sigest work narrows to rule 4's case (the routine called without the
  weights) and to the floor and the non-finite fallback, which that plan leaves at sd(y) of every row
  and this slice moves to the rows in; the call of faf8ca7d that a caller's sparse x in rule 4's case
  gives two warnings is gone. That plan's tests of the fallback through `dbarts` ("sd(y) with one
  dbartsSigmaFallbackWarning") are reread here: with no weight they hold, with a zero weight they read
  the rows in. It re-records eight scenarios, not `zeroweights`. Had this slice landed first, that plan
  would have had to carry the rows-in fallback into its routine's tests instead; it does not.
- Per-chain responses ([per-chain-responses.md](../design/per-chain-responses.md), for decision under
  dec-B419): if built, it follows this slice, whose range readers it reuses. This slice keeps the read
  (the range over the rows in and whether they hold two values) one routine apart from the write, so a
  pooled re-anchor can call it per chain; the status stays one answer per sampler, the pooled meaning
  left to that design.

## Out of scope, and where it goes

- Slices C and D, below.
- `rbart_vi`'s `rel.scale`: left as it is; the function goes in 1.1-0.
- A sampler made by hand from another sampler's control, model and data derives the range again under
  the weights in force and keeps the unit from the control: an existing seam, opened further by this
  rule; named in the help here, closed by the record slices of TODO `state-frame-prior`.
- A setter for the residual prior, which is how a caller moves what sigest set (dec-B314): TODO
  `sigma-prior-setter`.
- What stan4bart does about warm-up re-derivation itself ([Open items](#open-items)).

## Slices C and D, not planned here

C. `updateScale` on `setWeights` and `setActiveRows`. Waits on: this slice; leaf-conversions, for the
   conversion of the kept draws. Settled, to carry into its plan:
   - What a re-derivation does to the live chain on a call that changes no response: as
     `setResponse(updateScale = TRUE)` does, the live leaves keeping their internal numbers and the
     kept draws converted (dec-B327). The critic measured the live fit at the rows in moving by up to
     15.07 with no data change until the next sweep; the ruling accepts it.
   - Position. The new argument goes after `updateState`, so 0.9-34's positional `setWeights(w, TRUE)`
     keeps its meaning (dec-B327, dec-B381).
   - A sampler of several forests keeps refusing `updateScale = TRUE` in 1.0-0, the help saying the
     response's scale and a re-derived basis scale (dec-B303) may span different rows (dec-B362).
   - A re-derivation whose rows in hold fewer than two distinct values keeps the prior in force and
     warns once per sampler (dec-B364), as this slice builds for the two response setters; slice C
     extends it to the two new arguments, with `ResponseUpdateResult` returned through the same chain,
     sampler and facade path and the same holder flag (dec-B405). No flat entry sets weights or a mask,
     so the C API gains nothing.
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

## Open items

- Open dependency, recorded and not resolved here: the maintainer is deciding whether stan4bart stops
  re-deriving the response scale during warm-up
  ([5. Recommendation](../design/per-chain-responses.md#5-recommendation), rule (c)). stan4bart's part is
  written to hold either way (item 8 of [The keep's status (dec-B405)](#the-keeps-status-dec-b405)): if
  it stops, the warm-up read site goes (~60 to ~45 lines); if no call with `updateScale` true remains,
  the part is the rebuild alone. dbarts's side is unchanged either way; the C API's value still serves
  any host that re-derives.
- A premise for the orchestrator: dec-B405 names "a minor version bump"; dec-B60 holds the version pair
  through the pre-release ("No need to increment versions"), dbarts.h says the constants do not move
  before the first release, and dec-A139 held them at 1/0 for the last ABI event. The plan builds the
  bump as dec-B405 says, the newer and more specific ruling, and amends the header's paragraph and the C
  API test's pin. If dec-B405's words were not meant to reverse dec-B60, the alternative is the hash
  alone, moved by the new enumerator list: 1/0 stays, the test's pin and the paragraph stay, and a
  rebuilt stan4bart runs against an older dbarts 1.0 build and reads its 1 as "no keep". Worth one line
  to the maintainer before the build.

## Calls made in planning

Agent-made, each for the maintainer's later mark. The first seven are the coordinator's rulings on the
critique of the design (2026-10-07); then the planner's; then the revisions' (2026-10-09).

- Four slices: A the aft re-anchor fix; B creation; C the two new arguments; D dec-B304 with the record
  slices. Not taken: the design's one slice, estimated at 1380 lines before dec-B304 was ruled.
- B is creation only, with `updateScale = TRUE` on the two existing response setters reading the rows
  in. No new argument or bridge arity.
- Fewer than two distinct values among the rows in takes dec-B302's fallback and warning, not the window
  of width 1. Ruled by the maintainer since, as dec-B309; the measurement is in Context.
- With no row out every number keeps its bits, as a constraint on the loops and a gate.
- Two exact arms that fail on the base build.
- For slice C: the argument after `updateState`; a re-derivation through it treats the live chain as
  `setResponse(updateScale = TRUE)` does. Both ruled since (dec-B327, dec-B381).
- `rbart_vi` is cut: deprecated on the tip.
- The value tested for "fewer than two distinct" is the count for a count fit. So rows in that all hold
  one positive count take the fallback although their center is defined; the alternative is a second
  test for one family. Approved (dec-B381).
- The warning is raised only by families that take something from the response. A probit or ordinal
  sampler under a mask of zeros is created in silence, as today; dec-B302's sentence would be false of
  it. The alternative is a warning in every family with other words. Approved (dec-B381).
- In rule 4's case sigest is the fit without the weights (25.006 on the fixture) where today it is the
  marginal sd (25.245), and the sigma fallback's warning is not raised: one warning, as dec-B302 says.
  Approved (dec-B381).
- Creation under a mask is R's install followed by the sampler's own `setOffset`, for nbinom only.
  Not taken: the engine's creation taking the mask, which changes the creation entry's arity. The
  by-hand route's documented call is dec-B363's, with slice C.
- `xbart` warns once, from the estimate on all rows, and a fold in rule 4's case falls back in silence.
  Approved (dec-B381).
- The unit follows the response model's decision of which rows were read; it does not test again.
- The warning's class is `dbartsScaleFallbackWarning` (with `dbartsWarning`), for the creation text, the
  keep's and the hazard grid's.
- The negbin arm runs inside negbin-exact.R; the bd-balance arm is a named arm run by the workflow.
- The test of a state across row sets does not assert `setState`'s return value beyond a clean install.
- NEWS names the response's scale only, with what 0.9-34 did.
- Revision. dec-B364 is built here for `setResponse` and `setOffset` (the orchestrator), replacing the
  plan's every-row reading in silence at a re-derivation; slice C extends it.
- Revision. Rule 4's "where a row is out" stays, so dec-B364 (fewer than two distinct values) and
  dec-B386 (a constant response's window) do not overlap: a constant response with every row in takes
  the window and no warning, at creation and at a re-derivation (the orchestrator, the critique's M1).
- Revision. `setData` in rule 4's case is treated as creation: every row, with dec-B302's warning. dec-B364
  names four setters and not `setData`, which replaces the data (the orchestrator, the critique's M3).
  The alternative was the keep, which would hold a scale from data no longer there.
- Revision. The hazard periods (dec-B315) come into this slice: TODO `hazard-periods-rows-in` schedules
  them "R only; with the creation slice of response-scale-rows, whose rule it follows". The
  reviser's reading of that entry and of dec-B315, which puts the periods in dec-B296's class.
- Revision. `breaks = K` takes its quantiles from the subjects in as the default grid does: the
  positions are read from the response, and dec-B315 leaves only stated boundaries as they were. The
  alternative was every subject's quantiles, which a far-out subject at weight 0 moves (20 against 5).
- Revision. A hazard fit's only fallback is no subject in. One subject in, or several at one time, is a
  grid of one period, which is a defined fit; dec-B309's window problem has no counterpart.
- Against the design and the critique. Neither said what sigest is in the no-row case, where the tip
  already warns once for another cause, so "one warning" needed the call above. The aft behavior
  both report as found is recorded as a property in
  [aft-status-setter.md](../design/aft-status-setter.md). Every measured claim rerun here held: the
  per-chain aft transforms (1.872, 2.098), the unit's last bit at n = 64, the two pin failures, one row
  in giving `k.scale` 0.5 on the prototype.
- Revision to dec-B405. The status type is `ResponseUpdateResult { accepted, scaleKept }`, returned by
  the response model, the chain, the sampler and the facade; the alternative was a flag on the response
  model read after the call through a new facade read, which dec-B405's "a status the engine's
  setResponse and setOffset return" rules out.
- Revision to dec-B405. The once-per-sampler flag is the holder's (`scaleKeptWarned`), spent only by the
  R entries, so a flat keep leaves the R interface's warning for the R methods' first keep. The
  alternative, a flag on the R object, needs a new field on the reference class; the holder already
  carries the flat run's equivalent and has the sampler's lifetime.
- Revision to dec-B405. The flat value is named by a fourth enumerator list, `dbarts_swap_status`
  (`DBARTS_SWAP_NONE`, `DBARTS_SWAP_DONE`, `DBARTS_SWAP_SCALE_KEPT`), folded into the token, the entries'
  return type left `int`. The token folds signatures, lists and layouts only, so a bare new value would
  not move the hash dec-B405 names; a named list gives the value a name a consumer compares against.
  The alternative, the return type changed to the enum, changes two signatures for no reader's gain.
- Revision to dec-B405. stan4bart counts keeps per chain and warns once per fit from its parent process,
  recommended to it; the alternative, recording the count in the fit silently, is its call.
- Revision. sparse-starting-sigma lands first and this slice is written against its routine (the
  coordinator's reading of dec-A191's queue, which it is not in).

## Revision from critique

The plan's critique (2026-10-09, at 8ad1fa58) said build after corrections; its nine edits were made at
faf8ca7d. Step numbers in this list are that revision's; the dec-B405 revision inserted step 4, so its
steps 4 to 13 are now 5 to 14.

1. dec-B364 is built for the two response setters: rule 3, rule 4's keep, steps 1, 2 and 7, two
   mutations; the warning's route was held then and is ruled since (dec-B405).
2. The baseline is the main corpus's current file; the (none) partition misses by 106.6, rerun.
3. Step 7's creation routes drop the by-hand call, which slice C builds (dec-B363); step 13 says what
   dbartsSpec.Rd carries meanwhile; step 5's sentence is tidied.
4. The ruled questions (dec-B314, dec-B315, dec-B362) leave Out of scope and slice C; slice C's first
   bullet is dec-B327.
5. The hazard periods come in: rule 8, step 6, the `rng:` block, a fallback and tests.
6. Step 4 names `xbart`'s own call and its `sigmaPerFold` choice, and the sparse paths; step 7 has a
   sparse arm; the order against `sparse-starting-sigma` is stated.
7. Rule 4's single warning is stated for 0/1 weights, the validity warning pinned where it fires, the
   Literals fixture 0/1.
8. bd-balance's new arm is the fourth, with its `stopifnot` and the workflow's loop.
9. Status cites the rulings; M1 and M3 are recorded as calls.

Rerun by the reviser on a shipped build of 8ad1fa58: the bd-balance arm (106.6; 1.4 against every
row), the hazard grid, the sparse, fractional-weight and `xbart` cases, a masked probit row's response.
Not rerun: the prototype's compares and suite, which are the critic's.

## Revision to dec-B405 and the engine order (2026-10-09)

On bartcore a8884ef4. The held question is replaced by the ruled design
([The keep's status (dec-B405)](#the-keeps-status-dec-b405)) and step 4; the constraints on the facade
and the header follow it; stan4bart's part is planned and sized. Read for the rebase: the k-internal plan (bfe38871),
setstate-force-update (wt/install-surface-plan), default-rule-per-gap's baseline row,
sparse-starting-sigma (d1be8904) and per-chain-responses.md; [Interactions](#interactions) says what each
changes here. Landed since faf8ca7d's window: leaf-conversions and slice A. The baselines are named as
the current files at the slice's base. Nothing was rerun; the new lines are marked (read).
