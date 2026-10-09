# response-scale-rows: what the response gives a prior is taken from the rows in the likelihood at creation

Status: PLANNED, revised from its critique 2026-10-09. Implements dec-B296, dec-B302, dec-B309, dec-B314
(sigest held), dec-B315 (a hazard fit's periods), dec-B362 (several forests keep refusing), dec-B363
(slice B's half: no documented by-hand call yet), dec-B364 (for `setResponse` and `setOffset`),
dec-B381 (the planning calls approved) and dec-B386 (a constant response's window where every row is
in); rests on dec-B357 (a survival fit's range from observed times) and, for slice C, dec-B327; revises
dec-A140. One question is held:
[Held for the maintainer: the flat entries' warning](#held-for-the-maintainer-the-flat-entries-warning).
This file plans slice B, creation. Slice A is
[aft-reanchor-observed-times.md](aft-reanchor-observed-times.md). Slices C and D are stated at the end
and have no steps yet. dec-B303 stays in [forest-multiplier.md](forest-multiplier.md).

agent: one push. Opus implementer for the engine and its tests; sonnet for the R code (the hazard grid
included), the tinytest file, the gate scripts and the help once the engine is fixed; opus reviewer
told to refute, and to read the diff twice: once for bits that must not move where no row is out, once
for which rows each number reads.
rng: by call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines
them. A row is IN when its weight is positive and no mask switches it off, OUT otherwise; a hazard
subject is in when its weight is positive.
- NEUTRAL, bit for bit: every sampler created with every row in, whatever it does afterwards (a mask or
  zero weights installed after creation included); every probit, logistic, ordinal and multinomial
  sampler; a hazard fit with every subject in, with stated boundaries, or on the default grid whose
  subjects at weight 0 hold only times some subject in holds; a single-forest gaussian or Student-t
  sampler created with rows out whose smallest and largest response less offset are both on rows in.
- POSTERIOR-CHANGING, a default changes, for a sampler created with rows out: a gaussian or Student-t fit
  with an extreme on a row out (the range); a gaussian model of several forests with any row out (the
  unit, so every forest's default sd, and what a stated `sd` means until
  [forest-sd-unit.md](forest-sd-unit.md) lands); a count fit under 0/1 weights (the centre); `xbart` with
  a weight of 0 (each fold's range); any of them where sigest takes its fallback, a caller's sparse x
  included; a hazard fit with a subject at weight 0 whose time no subject in holds, or under
  `breaks = K` with any subject at weight 0 (the periods, and the rows of the subjects in). And
  `setResponse` or `setOffset` with `updateScale = TRUE` while rows are out, in every family that has a
  range or a centre, through the R methods and the two flat C entries: the rows in, or, where they hold
  fewer than two distinct values, the transform kept in force where the tip reads every row.
- Moves one recorded baseline scenario, `zeroweights` of the main corpus; nothing else recorded moves
  (the main corpus's `hazard` scenario carries no weight, read).
window: after [leaf-conversions.md](leaf-conversions.md) and slice A, which work in the same routines;
before [forest-sd-unit.md](forest-sd-unit.md), so the unit's tests are written once against the rows
that stay. Before the merge to main: it changes draws against 0.9-34. Serial with any other work in
model.hpp, the K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp),
[`resolveSamplerSpec`](../../R/spec.R), [`estimateSigmaFromLinearModel`](../../R/utility.R) or
[`expandDiscreteTimeHazard`](../../R/dbarts.R). TODO `sparse-starting-sigma` may land either side; step 4
says how the sparse arm reads under each outcome.
budget: ~1850 lines (engine ~135; the bridge ~40 under answer A of the held question, ~0 under B;
tests/cpp ~300; R ~170, of which the hazard grid ~30; tinytest ~600, of which ~560 a new file;
benchmarks and the workflow ~250; help, header comments and NEWS ~110; design note, indexes, amended
notes and TODO ~245), upper figure 3300. What it rests on: the prototype's four edits were 60 lines and
proved the engine half; the R half (sigest, the warnings, creation under a mask, the hazard grid) was
not prototyped; leaf-conversions, in the same files, came to 2.3 times its plan, almost all of it tests.

## Goal

A fit created with rows out of the likelihood has the leaf prior, the prior centre, the default forest
sds and, for a hazard fit, the periods of the fit on the remaining rows: the response of a row out
reaches no number and no draw at the rows in. Those numbers are then held, as they are today. A fit
created with every row in is the fit it is today, bit for bit.

## Context

Marks: (ran) by the planner or the reviser on the tip's build, shipped mode, and on the design's
prototype library where a line says so; (stands) measured by the design, its critic or the plan's
critic and not rerun; (read) in the source. One fixture unless a line says otherwise: 200 rows, two
predictors, rows 101 to 200 at weight 0 with their response 50 higher.

What reads which rows today.
- (ran) The range ([`GaussianResponse::rescale`](../../src/bartcore/model.hpp)): every row. Recorded
  [-2.165, 56.371] where the rows in span [-2.165, 6.029] and `subset` gives the latter; `k.scale` 29.27
  against 4.10 on the prototype. Student-t and aft hold a Gaussian response and inherit it (read).
- (ran) The count centre ([`NBResponse::computeShift`](../../src/bartcore/model.hpp)): every row. Half
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
- (ran) sigest's sparse paths (dec-B370). A caller's sparse x (a `dgCMatrix`) takes the marginal
  fallback with its own warning: 25.245, the sd of every row, where the rows in give 1.466; the same
  25.245 and one warning with every weight 0 or one row in. An indicator expansion R built sparse goes
  through [`sparseResidualStandardError`](../../R/utility.R), already over positive weights (read);
  with every weight 0 it returns a non-finite value and `floorSigmaEstimate` falls back to 25.245, one
  warning. Since the stray-warnings slice no sparse fit in the suite reaches the fallback (stands).
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
  0.5, and a count fit with one row in at zero the centre log(1/2). (stands) Then with every row brought
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
  scenario carries no weight (read).
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
  higher, left out by 0/1 weights at creation, the oracle's likelihood over the rows in and its centre
  over every row (3.334) or the rows in (1.056). 2 seconds for two seeds. Tip: gap 0.012 to the every-row
  oracle, 0.493 to the rows-in one, tolerance 0.12. Prototype: 0.484 and 0.003.

Released behaviour. (read) 0.9-34 took the range over every row, at creation and at
`setOffset(updateScale = TRUE)`, and had no `updateScale` on `setResponse`; (stands) on its build two
seeded fits that differ only in the responses of the zero-weight rows differ by up to 3.65 in a drawn
fitted value at the rows in. Counts, several forests, masks, Student-t and hazard fits reached no
release.

## The rule

1. Rows. The range (gaussian, Student-t, aft, and under a variance forest or any leaf model), the count
   centre and the unit of a forest's sd are computed over the rows in when the sampler is created,
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
     been given, sigest included (dec-B302, dec-B309). One warning, of its own class, and the sigma
     fallback's warning is not raised for that cause:

         fewer than two distinct values of the response are in the likelihood (<m> of <n> observations are in it), so the response's scale and the prior's centre are taken from all <n>

   - At a re-derivation through `setResponse` or `setOffset`, the transform and the centre in force are
     kept: the call completes as it does with `updateScale = FALSE`, bit for bit, and a warning of the
     same class is raised once per sampler (dec-B364); where it is raised from is the held question:

         fewer than two distinct values of the response are in the likelihood, so the response's scale and the prior's centre are kept as they were (said once per sampler)

   The value tested is the response less its offset for a range and the count for a centre. With every
   row in, a constant response is not this case: it takes dec-B386's window c - 0.5 to c + 0.5 at
   creation and at a re-derivation, with no warning of this class, so dec-B364 and dec-B386 never
   apply to the same call.
5. Nothing out, nothing changed. With every row in, every number above has the bits it has today.
6. sigest stays `lm`'s weighted estimate over the rows in and is re-estimated by no call (dec-B314). Its
   fallback sd and its floor read the rows in, on each of its three paths: the dense linear fit, an
   indicator expansion's sparse QR, and a caller's sparse x's marginal sd.
7. `xbart`: each fold's range and sigma follow 1, 4 and 6 over the fold's training rows.
8. Hazard periods. A hazard fit's default periods, and under `breaks = K` its quantiles, are read from
   the times of the subjects in the likelihood when it is created (weight positive, after `subset` and
   the na.action); stated boundaries are as they were, still checked against every subject's time. A
   subject at weight 0 is expanded on that grid, its time in the interval it falls in; one past the last
   period is at risk in every period with no event (dec-B315). With no subject in, the periods are every
   subject's, with one warning of rule 4's class:

       no subject is in the likelihood (0 of <n> subjects are in it), so the hazard periods are taken from all <n>

## Constraints

- No new argument; no bridge entry changes arity; no facade virtual, unless the held question is
  answered A, which adds exactly one read. `--preclean` on the engine commit, for the edited headers.
- Rule 5 is a constraint on the code, not only a test: the unit's loop and the centre's loop run as they
  are written today whenever no row is out, and the rows-out case is a separate path. A reviewer checks
  that neither present loop gained a branch.
- One decision of which rows are read. The response model makes it; the unit follows the response
  model's answer and does not test the rows again. Rule 4's keep at a re-derivation is engine-local:
  the range read in `GaussianResponse` and `NBResponse::computeShift` at its two re-deriving callers
  leave the transform where they find fewer than two values (stands, the plan's critic).
- The flat C header keeps every signature; the comments of its two re-anchoring entries gain "over the
  rows in the likelihood" and the keep. The API hash does not move.
- No state format change and no new record: the range is the model's `response.range`, the unit the
  forests' `anchor`, sigest the data object's, the periods the fit's `$periods`, as now. A re-creation
  restates the record and derives nothing.
- Cut points, a hazard fit's stated boundaries, the number of categories, leaf covariate
  standardization and a multiplier's scale are not touched.
- [`rbart_vi`](../../R/rbart.R) is not touched.

## Steps

"Fails today" is what the tip does where the test expects otherwise.

1. Engine, the range. [`GaussianResponse::rescale`](../../src/bartcore/model.hpp) reads the rows whose
   served precision is positive from the vector slice A has it take. In rule 4's case the constructor
   and `setData` read every row, and `setResponse(true)` and `setOffset(true)` keep the transform in
   force and record that they kept it, for the warning. Student-t (its composite precision is zero
   exactly where the weight or the mask is) and aft (the mask only) follow with no edit of their own.
   The response model answers which rows it read, for step 3. tests/cpp, beside
   [`testActiveRowsGaussianDf`](../../tests/cpp/test_model.cpp), two chains where a chain is involved:
   - gaussian, Student-t and aft (aft by a mask then `setOffset(true)`), the weights of the rows in
     fractional: the scale and shift equal the unweighted literal over the rows in (fails today: every
     row);
   - no weight at all, and weights all positive: scale, shift and the working response are bit for bit
     those of the present expression, written out in the test;
   - created with none in, one in, four in at one value: the every-row literal; `setData` to those rows:
     the same;
   - re-derived with a row out and none in, one in, four in at one value: the transform, the working
     response and 20 sweeps bit for bit those of a twin given the same call with `updateScale` false,
     and the kept flag set (fails today: every row); with every row in and a constant response:
     dec-B386's window and the flag clear;
   - held: after `setWeights` and `setActiveRows` the transform is unchanged and 20 sweeps equal an
     untouched twin's under the same weights;
   - `setResponse(true)` and `setOffset(true)` under zero weights and under a mask read the rows in;
     `setData` reads the rows in under the new weights.
2. Engine, the centre. [`NBResponse::computeShift`](../../src/bartcore/model.hpp) takes the mask where
   one is installed; the unmasked call is the present loop. Its two re-deriving callers pass the mask,
   and where the rows the mask keeps hold fewer than two distinct counts they keep the centre in force
   and record it. tests/cpp, beside [`testNBLogMeanAnchor`](../../tests/cpp/test_model.cpp): the centre
   after a mask and `setOffset(true)` equals the literal over the rows the mask keeps, with and without
   an exposure offset; unmasked, bit for bit the present value; a mask keeping one row, or rows of one
   count, then `setOffset(true)`: the centre in force, bit for bit the `updateScale`-false twin's (fails
   today: every row); held across `setActiveRows`. The comments of
   [`NBResponse::setActiveRows`](../../src/bartcore/model.hpp) and
   [`AFTResponse::setActiveRows`](../../src/bartcore/model.hpp) stop saying "full-data".
3. Engine, the unit. [`scaledResponseSd`](../../src/bartcore/chain.hpp) is left as it is and is what
   runs when the response model read every row; otherwise a second routine takes the sample sd over the
   rows it read, n - 1 of those rows in the divisor. tests/cpp, beside
   [`testForestCalibration`](../../tests/cpp/test_sampler.cpp): with rows out the map's anchor equals the
   literal sd to 1e-14; with none out it is bit for bit the present loop's value at n = 64, 100, 137 and
   250 (the size at which the prototype moved a bit among them); a finite carried anchor still wins.
4. R, the rows and the warnings. One function in [spec.R](../../R/spec.R) answers rule 4 for a data
   object and a mask and raises the creation warning, for the families that take something from the
   response (gaussian with any residual distribution, nbinom) and no other. Its callers:
   [`resolveSamplerSpec`](../../R/spec.R), after [`enforceWeightPolicy`](../../R/spec.R);
   [`bartcoreSamplerSetData`](../../R/bartcore.R), for `setData`; and [`xbart`](../../R/xbart.R), once,
   after its own `enforceWeightPolicy` call and before the starting-sigma estimate whose warning chooses
   `sigmaPerFold`.
   - sigest. [`estimateSigmaFromLinearModel`](../../R/utility.R), in rule 4's case, fits without the
     weights on the dense path and on the sparse QR path
     ([`sparseResidualStandardError`](../../R/utility.R) given no weights), and a caller's sparse x
     takes the sd of every row; otherwise the dense fit and the QR keep the weights, and the caller's
     sparse x's marginal sd, [`floorSigmaEstimate`](../../R/utility.R)'s fallback and
     [`floorMarginalSigma`](../../R/utility.R) read the rows in. The non-finite fallback's warning is
     not raised for rule 4's cause; the sparse fallback's is, its cause being the storage, so a caller's
     sparse x in rule 4's case gives two warnings.
   - `xbart`. In rule 4's case the all-rows estimate is the unweighted fit, which raises no fallback
     warning, so `sigmaPerFold` is "linear" (or "marginal" for a caller's sparse x, as today). Each fold
     tests rule 4 over its training rows with no warning: in its case the fold's fit is unweighted
     ("linear") or its sd is over every training row ("marginal"); otherwise the weighted fit, and the
     sd and the floor over the training rows in. So `xbart` warns once and a fold falls back in silence.
   - TODO `sparse-starting-sigma` is open with the maintainer. If it keeps the caller's-sparse fallback,
     the sparse arm of step 7 tests that fallback reading the rows in and the indicator expansion's QR.
     If it removes the fallback, every sparse x takes the QR path, which already reads the rows in,
     and the arm tests the QR alone for both sources. Whichever lands second rereads the other's sparse
     tests; neither changes this slice's rule.
   - The re-derivation's warning. Answer A of the held question: the bridge raises it from the engine's
     kept flag, for the R methods and the flat entries alike. Answer B: the R methods `setResponse` and
     `setOffset` run this step's function before the bridge call and warn once from a flag on the
     sampler object; the flat entries keep in silence.
5. R, creation under a mask. Of the families whose 0/1 weights become the mask
   ([`isMaskedWeightFamily`](../../R/spec.R)) only nbinom takes anything from the response. For it,
   where rule 4's case does not hold, [`dbarts`](../../R/dbarts.R) and
   [`bartcoreSamplerSetData`](../../R/bartcore.R) follow their install of the mask with the sampler's own
   `setOffset(offset in force, updateScale = TRUE, updateState = FALSE)`, before the sampler is handed
   back: it draws nothing, records the range and restates a named sd. In rule 4's case the follow-up is
   skipped and the every-row centre of creation stands, under step 4's warning. A sampler built by hand
   from [`dbartsSpec`](../../R/spec.R) and `new("dbartsSampler", ...)` that installs `spec$active`
   itself is created with every row in and, for a count fit, keeps the every-row centre; in this slice
   its help says so and gives no call. The documented call, `setActiveRows(spec$active, updateScale =
   TRUE)` with its loop form, comes with slice C, which builds it (dec-B363); the `setOffset` route is
   not documented.
6. R, the hazard periods. [`expandDiscreteTimeHazard`](../../R/dbarts.R)'s two callers in
   [`dbarts`](../../R/dbarts.R), the direct response and the formula's, pass as `gridTime` the times of
   the subjects kept by the na.action whose weight is positive, all of them where no weight is given or
   none is positive (rule 8's fallback, raising its warning there). [`resolveHazardGrid`](../../R/dbarts.R)
   takes the default grid and the `breaks = K` quantiles from it unchanged; stated boundaries keep their
   check over every subject. A subject whose time is past the last period is expanded to every period
   with a response of 0. Nothing in the engine.
7. tinytest, a new file test-response-scale-rows.R, fixtures with the rows out extreme in the response:
   - Inert. Two fits that differ only in the response of the rows out at creation are `identical()` at
     the rows in and within 1e-12 at the rows out, in draws, sigma and every reader: gaussian, a fixed
     sigma, Student-t, a variance forest, linear and gp leaves (covariates equal across the pair), a
     named sd, two forests on a factor and on a number, an offset, nbinom, and `bart`, `bartBT` and
     `xbart` once each. Fails today on every arm.
   - Literals, with 0/1 weights. `response.shift`, `response.scale`, `k.scale`, `prior.mean`, the
     forests' `anchor` and sigest equal values written from the rows in, and equal the `subset` fit's.
     With fractional weights in, sigest is the weighted fit and the `subset` fit's is not, so that arm
     asserts the range, the centre and the unit only.
   - Held. The readers are unchanged across `setWeights`, `setActiveRows`, a copy and a reload, read from
     the engine (`getLeafPrior`) and not from the record alone; five sweeps after a copy agree with the
     original's to 1e-10.
   - `updateScale = TRUE` on `setResponse` and `setOffset` under zero weights and under a mask gives the
     readers of a sampler created under those rows; through the flat entries once
     ([`dbarts_sampler_setOffset`](../../src/C_interface.cpp)); several forests still refuse.
   - The keep (dec-B364). Gaussian under zero weights and nbinom under a mask, re-derived through each
     of `setResponse` and `setOffset` with none in, one in, and several in at one value: the readers and
     five sweeps equal an `updateScale = FALSE` twin's; exactly one warning of the new class on the first
     such call, counted with `withCallingHandlers`, and none on a second; a constant response with every
     row in: dec-B386's window and no warning. Through `dbarts_sampler_setOffset` once: under answer A
     one warning, and none from a following R call (one flag); under answer B none.
   - Rule 4 at creation. Every weight 0; a mask of zeros (nbinom); one row in; rows in at one value:
     created, exactly one warning of the new class, counted with `withCallingHandlers`, its text as above
     with m and n; the readers and sigest equal those of the fit with no weights; no warning where every
     row is in and the response is constant; none for probit under a mask of zeros. Under fractional
     weights with a zero among them the data object's validity warning is raised as well, as today, and
     pinned by its own pattern exactly where it fires. `setData` to such rows: every row, the same one
     warning.
   - sigest: three rows in and two predictors gives the sd of the rows in with the existing fallback
     warning; the floor reads the rows in. Sparse: a caller's `dgCMatrix` x under the fixture gives the
     sd of the rows in (1.466, where the tip gives 25.245) beside the sparse fallback's warning; an
     indicator expansion (`factors = "indicators"`) with three rows in reaches the QR's fallback and
     reads the rows in; in rule 4's case both read every row, the caller's sparse x with its sparse
     warning beside the new one.
   - `xbart` with 60 rows in at one value: one warning of the new class, none per fold, and each fold's
     sigma the unweighted fit over its training rows.
   - Creation routes for a count fit: `dbarts` with 0/1 weights and `setData` with 0/1 weights give the
     rows-in centre; `dbartsSpec` then `new("dbartsSampler", ...)` and `setActiveRows(spec$active)` gives
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
   - A state across row sets: a single-forest state from a sampler created under other zero weights
     installs converted, the recipient's range unmoved and the installed fit equal to the donor's to
     1e-12; a gp and a two-forest state are refused by the existing message. The return value of
     `setState` is not asserted (dec-B305 revises it).
8. Restated pins. In test-active-rows-pins.R the two-forest pin creates both samplers under w and masks
   or reweights after; the gaussian and Student-t pins above it get the same repair, since they pass only
   while the fixture's extremes are on rows in. The comments
   ["The log-mean shift is the full-data one by design"](../../inst/tinytest/test-active-rows-pins.R) and
   ["the response transform is the FULL-data one by design"](../../inst/tinytest/test-active-rows-pins.R)
   say "creation's, every row then being in".
9. Exact gates. bd-balance.R gains a fourth arm beside [`zeroWeightArm`](../../benchmarks/R/bd-balance.R),
   `narrowArm` and `constMissingArm`: the weights at creation, the zeroed cells' response 3 higher, the
   oracle's range over the rows in; its `stopifnot` over the arms takes the fourth. negbin-exact.R gains
   a third arm beside its fixed-shape one, as measured in Context, with the oracle's centre
   ([`shiftOf`](../../benchmarks/R/negbin-exact.R)) over the rows in.
   [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) runs the bd-balance arm in its
   `for arm in narrow constmissing` loop, whose comment widens; the negbin arm runs inside its script.
   About 3 seconds more in `quick`. benchmarks/README.md names both.
10. The baseline. `zeroweights` is recorded fresh on a reference build at the slice's newest code commit
    and the other 54 scenarios are carried from equivalence-51201107 after they compare bitwise; the
    [MANIFEST](../../benchmarks/baselines/MANIFEST) row says so in its predecessors' form and names the
    oracle (step 11's pair script, rows 3 and 4). The workflows that pin the file's name follow.
11. The pair script, benchmarks/R, run on the base and the slice builds with the same seeds, 20 sweeps:
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
12. Mutations (Verification).
13. Records. Help: the `weights` item of [bart.Rd](../../man/bart.Rd) and [dbarts.Rd](../../man/dbarts.Rd)
    (the sentences holding ["still counts toward"](../../man/bart.Rd) and
    ["the response scale held equal"](../../man/bart.Rd)) say that a row of weight 0 takes no part in the
    response's scale, the count centre or a forest's default sd, that it stays in the design and the fit
    is close to the fit on the remaining rows without being it, that the step is at 0 and not at a small
    weight, and that a mask or weights set after creation leave the numbers creation set, so the three
    ways to leave a row out differ by when it left; the nbinom paragraph names the rows of c; the hazard
    paragraphs of both and of [dbartsFamilies.Rd](../../man/dbartsFamilies.Rd)
    (["sorted distinct observed times"](../../man/bart.Rd)) say the subjects in the likelihood;
    [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the response item
    (["which reads every row"](../../man/dbartsSampler-class.Rd) becomes the rows in), the weights item,
    the `active` item (["stay the FULL-data calibration"](../../man/dbartsSampler-class.Rd) becomes
    creation's), `updateScale` with the keep and its warning, and that `setData` and a sampler made by
    hand from another's parts derive again; [dbartsSpec.Rd](../../man/dbartsSpec.Rd), the `active`
    value: a count sampler made by hand keeps the every-row centre in this slice; xbart.Rd, `weights`.
    NEWS, user-visible changes, one item: "An observation of weight zero no longer takes part in the
    response's scale, which the leaf prior is stated against, so a fit with zero weights has the leaf
    prior of the fit on the remaining rows; 0.9-34 took the scale over every row. The observation is
    still in the design." The keep and the hazard periods reached no release and take no item. A design
    note, docs/design/response-scale-rows.md, with its index row: the rule, the inventory, the
    measurements, the named-sd case (a prior centred at 52 where the rows in sit near 2), why
    membership and not weight, the keep, the hazard grid, the four slices.
    [active-rows-mask.md](../design/active-rows-mask.md),
    [negative-binomial.md](../design/negative-binomial.md), [nbinom-log-mean.md](nbinom-log-mean.md) and
    [survival.md](../design/survival.md) take a dated amendment where they say "full-data", "all rows" or
    the hazard grid's source. dbarts.h's two comments. TODO: the `response-scale-rows` entry is cut to
    slices C and D; `hazard-periods-rows-in` is removed. The landing note here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, clean under ASan
  and UBSan; the full tinytest suite, expected clean after step 8 (2 failures without it); the new file
  under ASan on the R-loaded path.
- On a reference build, the three compares, counted per scenario with no `max |z|` line but the one
  named:
  - main corpus against `equivalence-51201107.rds`: 54 identical; `zeroweights` differs and is the only
    one. Its statistical verdict against the old baseline is reported, not gated: the posterior changed
    (the re-ported prototype's was max abs z 2.85). Then all 55 identical against the new file.
  - two-forest corpus against `bcf-equivalence-1b7d730c.rds`: 15 of 15 identical. Fewer is rule 5 broken
    in the unit's loop.
  - multinomial corpus against `multinomial-equivalence-80b1c8d4.rds`: 11 of 11 identical.
  - the four `test-reproducibility-*.R` files pass unchanged; none is regenerated.
- Every gate of [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) in `quick`, unchanged, and
  the two new arms. Each new arm is also run once on the base build, where it must FAIL (measured: by
  106.6 standard errors and by 0.49 against 0.12), and that run is reported.
- The pair script: rows 1 to 4 `identical()`, row 5 differing.
- Mutations, each to fail the named check:
  - the range over every row: step 7 "Inert" and "Literals"; the bd-balance arm;
  - the centre ignoring the mask: step 2's literal; step 7's nbinom arm; the negbin arm;
  - the unit over every row: step 3's literal; step 7's two-forest arms;
  - the unit's present loop replaced by the rows-in loop when no row is out: step 3's bits; the
    two-forest compare;
  - the centre's present loop given the mask test: step 2's bits (and the main corpus's count scenario
    if a bit moves);
  - weighted, not membership (a row counted by its weight): step 1's literal with fractional weights;
  - rule 4 off (one row in gives the window of width 1): step 1's "one in"; step 7 "Rule 4 at creation";
  - rule 4 on with every row in (a constant response warns): step 7's "no warning";
  - the keep dropped (a re-derivation in rule 4's case reads every row): step 1's and step 2's
    re-derivation tests; step 7 "The keep";
  - the keep's warning raised on every such call: step 7's "none on a second";
  - the creation fallback at `setData` replaced by the keep: step 1's `setData` test; step 7's `setData`
    arm;
  - creation's re-derivation under the mask dropped: step 7 "Creation routes"; the negbin arm; and run
    for probit too: the suite's seeded probit-under-mask pins, if the call draws;
  - sigest's fallback left on every row: step 7 "sigest", the sparse arm included; both warnings raised
    with no row in: step 7's count of one;
  - `xbart`'s folds warning each: step 7's `xbart` count of one;
  - the hazard grid over every subject: step 7 "Hazard", the periods and the `identical()` pair;
  - `updateScale = TRUE` left reading every row: step 1's last bullet; step 7's fourth.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; inst/NEWS.Rd parses with a non-NULL result; `R CMD check --as-cran` on a tarball from a clean
  copy. No new Rd topic.
- Cost: creation and the two response setters gain a comparison per row; a hazard fit one subset of its
  times; no sweep path. No bench compare.

## What a consumer sees

Read in each package's source (weights, `subset`, masks, `updateScale`); no consumer suite was run.

- stan4bart. A gaussian fit hands its weights to the sampler at creation, zeros allowed (it refuses only
  all zero), and then calls the flat `setOffset` entry with `updateScale` true at the start and through
  warm-up. Today both read every row. After this slice creation reads the rows in, and so must every one
  of those calls, or the first would undo creation's scale: that is why rule 3 is in this slice. A
  gaussian stan4bart fit with a zero weight therefore changes its leaf prior and draws; one without is
  bit for bit unchanged. Where the rows in hold fewer than two distinct values a call keeps the scale
  in force, with a warning or without it by the held question. Its binary fits install a mask after
  creation on a probit sampler, which takes nothing from the response. No test of its suite shows the
  change: test-15-weights_offset_k.R draws its weights from (0.5, 2) and from {4, 0.25}, and
  test-22-binary-weight-mask.R is probit. To verify: its suite on the slice's library, unchanged; and
  one seeded gaussian fit with zero weights on the base and slice libraries, differing, beside one
  without, identical.
- bartCause passes the user's weights and `subset` at creation, to `bcf` too, and changes neither
  afterwards. A fit with a zero weight moves as the `rng:` block says, under `bcf` by the unit as well.
  No test of its suite passes a zero weight (searched).
- treatSens fits a probit treatment model and re-anchors it through the flat entry; a probit sampler
  takes nothing from the response. bairrtt uses no weights, no mask and no `updateScale`.

## Out of scope, and where it goes

- Slices C and D, below.
- `rbart_vi`'s `rel.scale`: left as it is; the function goes in 1.1-0.
- A sampler made by hand from another sampler's control, model and data derives the range again under
  the weights in force and keeps the unit from the control: an existing seam, opened further by this
  rule; named in the help here, closed by the record slices of TODO `state-frame-prior`.
- A gp or several-forest state between two samplers on the same data created under different zero
  weights is refused where the tip installs it: kept, tested, and said in the help.
- A setter for the residual prior, which is how a caller moves what sigest set (dec-B314): TODO
  `sigma-prior-setter`.

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
     extends it to the two new arguments, through the same warning route.
   - The by-hand count sampler's call, `setActiveRows(spec$active, updateScale = TRUE)`, and the loop
     form `setActiveRows(mask, updateScale = i <= n.burn)`, written into dbartsSpec's help (dec-B363).
   - A help sentence, and a test, that `setState` does not undo a re-derivation.
D. dec-B304: a linear or gp leaf's covariate centre and scale, and the default lengthscale, over the
   rows in at creation, held, re-derived by `updateStandardization`; `xbart`'s folds taking their own.
   Waits on: the record slices of TODO `state-frame-prior` (a state install today replaces the
   recipient's centre and scale and reports a clean install, and a re-creation without a state derives
   them afresh, so there is nowhere to hold them); `updateStandardization` itself (dec-B233, not built).
   Open: creation under a mask must derive the standardization under it too, which this slice's
   `setOffset` step does not name, so either creation takes the mask or the install re-derives two
   things; the leaf is built from the column store, which holds no weights.

## Held for the maintainer: the flat entries' warning

dec-B364 says a re-derivation that keeps the scale warns once per sampler. stan4bart calls
`dbarts_sampler_setOffset` with `updateScale` true at the start and at every warm-up step. The keep
itself is engine-local and built either way; the question is only whether the flat entries say so. The
bridge cannot tell a kept transform from one re-derived to the same values through today's facade
(stands, the plan's critic; read: no `SamplerBase` member reports it).

- A. Lift "no facade virtual" for one read. [`SamplerBase`](../../src/bartcore/facade.hpp) gains one
  pure virtual const read, true when the last `setResponse` or `setOffset` with `updateScale` kept the
  transform, and [`SamplerFacade`](../../src/bartcore/facade.hpp) forwards it to the response model's
  flag of steps 1 and 2. The holder [`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp)
  gains one `bool` after `slowCountWarned` (the positional-aggregate rule its comment states), and one
  helper beside [`raiseSlowCountWarning`](../../src/C_interface.cpp) evaluates the warning from the
  dbarts namespace. [`dbarts_sampler_setResponse`](../../src/C_interface.cpp),
  [`dbarts_sampler_setOffset`](../../src/C_interface.cpp) and the R bridge's
  [`bartcore_setResponse`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setOffset`](../../src/R_interface_bartcore.cpp) call it after the swap, so the R methods
  and the flat entries share one flag and one decision. About 40 lines of engine and bridge, 20 of
  tests. ABI: the facade is internal; no signature in dbarts.h changes, `DBARTS_C_API_HASH` does not
  move and no consumer rebuilds; `--preclean` is mandatory on the commit, a stale object on the old
  vtable bus-erroring. Cost: a stan4bart fit sees an R warning raised from inside its warm-up loop the
  first time a component's rows in fall below two values; under `options(warn = 2)` it is an error
  unwinding through the caller's frames, as the flat run's slow-count warning already can be.
- B. A silent keep on the flat entries. No facade change: the R methods `setResponse` and `setOffset`
  run step 4's function before the bridge call and warn once from a flag on the sampler object; the
  flat entries keep in force without a word, and dbarts.h says so. About 20 lines of R. Cost: a
  stan4bart user whose component loses its rows mid warm-up is not told; the R methods' test is a
  second reading of which rows are in beside the engine's, held in agreement by step 7's tests against
  the constraint "one decision of which rows are read"; and the flat entries depart from dec-B364's
  "once per sampler", which needs the maintainer's leave.

Not offered: the bridge testing the rows from its own copies, a third reading, and it holds no copy of
the mask. The plan builds either way: steps 1 and 2 are the same; step 4's last bullet and step 7's
flat-entry assertion take the answer.

## Calls made in planning

Agent-made, each for the maintainer's later mark. The first seven are the coordinator's rulings on the
critique of the design (2026-10-07); then the planner's; then the revision's (2026-10-09).

- Four slices: A the aft re-anchor fix; B creation; C the two new arguments; D dec-B304 with the record
  slices. Not taken: the design's one slice, estimated at 1380 lines before dec-B304 was ruled.
- B is creation only, with `updateScale = TRUE` on the two existing response setters reading the rows
  in. No new argument or bridge arity; no facade virtual unless the held question is answered A.
- Fewer than two distinct values among the rows in takes dec-B302's fallback and warning, not the window
  of width 1. Ruled by the maintainer since, as dec-B309; the measurement is in Context.
- With no row out every number keeps its bits, as a constraint on the loops and a gate.
- Two exact arms that fail on the base build.
- For slice C: the argument after `updateState`; a re-derivation through it treats the live chain as
  `setResponse(updateScale = TRUE)` does. Both ruled since (dec-B327, dec-B381).
- `rbart_vi` is cut: deprecated on the tip.
- The value tested for "fewer than two distinct" is the count for a count fit. So rows in that all hold
  one positive count take the fallback although their centre is defined; the alternative is a second
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
- The test of a state across row sets does not assert `setState`'s return value, which dec-B305 changes.
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
- Revision. A caller's sparse x in rule 4's case raises the sparse fallback's warning beside the new
  one, the storage being another cause.
- Against the design and the critique. Neither said what sigest is in the no-row case, where the tip
  already warns once for another cause, so "one warning" needed the call above. The aft behaviour
  both report as found is recorded as a property in
  [aft-status-setter.md](../design/aft-status-setter.md). Every measured claim rerun here held: the
  per-chain aft transforms (1.872, 2.098), the unit's last bit at n = 64, the two pin failures, one row
  in giving `k.scale` 0.5 on the prototype.

## Revision from critique

The plan's critique (2026-10-09, at 8ad1fa58) said build after corrections; its nine edits are made.

1. dec-B364 is built for the two response setters: rule 3, rule 4's keep, steps 1, 2 and 7, two
   mutations; the warning's route is held (above).
2. The baseline is equivalence-51201107; the (none) partition misses by 106.6, rerun.
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
