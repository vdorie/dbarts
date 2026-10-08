# response-scale-rows: what the response gives a prior is taken from the rows in the likelihood at creation

Status: PLANNED (dec-B296, dec-B302; revises dec-A140). This file plans slice B, creation. Slice A is
[aft-reanchor-observed-times.md](aft-reanchor-observed-times.md). Slices C and D are stated at the end
and have no steps yet. dec-B303 stays in [forest-multiplier.md](forest-multiplier.md).

agent: one push. Opus implementer for the engine and its tests; sonnet for the R code, the tinytest
file, the gate scripts and the help once the engine is fixed; opus reviewer told to refute, and to read
the diff twice: once for bits that must not move where no row is out, once for which rows each number
reads.
rng: by call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines
them. A row is IN when its weight is positive and no mask switches it off, OUT otherwise.
- NEUTRAL, bit for bit: every sampler created with every row in, whatever it does afterwards (a mask or
  zero weights installed after creation included); every probit, logistic, ordinal, hazard and
  multinomial sampler; a single-forest gaussian or Student-t sampler created with rows out whose smallest
  and largest response less offset are both on rows in.
- POSTERIOR-CHANGING, a default changes, for a sampler created with rows out: a gaussian or Student-t fit
  with an extreme on a row out (the range); a gaussian model of several forests with any row out (the
  unit, so every forest's default sd, and what a stated `sd` means until
  [forest-sd-unit.md](forest-sd-unit.md) lands); a count fit under 0/1 weights (the centre); `xbart` with
  a weight of 0 (each fold's range); any of them where sigest takes its fallback. And for
  `setResponse` or `setOffset` with `updateScale = TRUE` while rows are out, in every family that has a
  range or a centre, through the R methods and the two flat C entries.
- Moves one recorded baseline scenario, `zeroweights` of the main corpus; nothing else recorded moves.
window: after [leaf-conversions.md](leaf-conversions.md) and slice A, which work in the same routines;
before [forest-sd-unit.md](forest-sd-unit.md), so the unit's tests are written once against the rows
that stay. Before the merge to main: it changes draws against 0.9-34. Serial with any other work in
model.hpp, the K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp),
[`resolveSamplerSpec`](../../R/spec.R) or [`estimateSigmaFromLinearModel`](../../R/utility.R).
budget: ~1500 lines (engine ~115; tests/cpp ~260; R ~110; tinytest ~460, of which ~420 a new file;
benchmarks and the workflow ~245; help, header comments and NEWS ~90; design note, indexes, amended
notes and TODO ~220), upper figure 2700. What it rests on: the prototype's four edits were 60 lines and
proved the engine half; the R half (sigest, the warning, creation under a mask as planned here) was not
prototyped; leaf-conversions, in the same files, came to 2.3 times its plan, almost all of it tests.

## Goal

A fit created with rows out of the likelihood has the leaf prior, the prior centre and the default
forest sds of the fit on the remaining rows: the response of a row out reaches no number and no draw at
the rows in. Those numbers are then held, as they are today. A fit created with every row in is the fit
it is today, bit for bit.

## Context

Marks: (ran) by the planner on the tip's build, shipped mode, and on the design's prototype library
where a line says so; (stands) measured by the design or its critic and not rerun; (read) in the source.
One fixture unless a line says otherwise: 200 rows, two predictors, rows 101 to 200 at weight 0 with
their response 50 higher.

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
- (ran) Held already: `setWeights` and `setActiveRows` move nothing. `setResponse` and `setOffset` with
  `updateScale = TRUE` under zero weights or a mask read every row ([-2.165, 56.371]); on the prototype
  the rows in. On several forests `updateScale = TRUE` is refused ("every forest keeps its leaf
  calibration stated against the scale fixed at creation").
- (ran) No row in. Every weight 0: the sampler is created, with one warning that names another cause
  ("starting sigma estimate falls back to the marginal response sd"), the range every row's, sigest 25.245
  (the marginal sd; `lm` without the weights gives 25.006). A count fit under a mask of zeros is created
  in silence.
- (ran) One row in. The tip reads every row. The prototype gives the range [1.666, 1.666] and `k.scale`
  0.5, and a count fit with one row in at zero the centre log(1/2). (stands) Then with every row brought
  back in and 600 sweeps the fit's rmse is 2.41 and sigma 2.66, where the fallback gives 0.46 and 0.90.
- (ran) `xbart` with a weight of 0: the response of the rows out reaches the result (1.38, 1.23 against
  1.15, 1.20 with those responses ordinary); on the prototype the two are identical, the engine edit
  alone closing it.
- (ran) `rbart_vi` is deprecated on the tip ("is deprecated and is removed in dbarts 1.1-0"), so its
  `rel.scale` is not in this slice.

What a correct change must not move.
- (ran) The prototype's rewrite of the unit's loop changed its last bit on unweighted data: at n = 64 the
  tip gives 0x1.c828f537ceda5p-3 and the prototype ...da4p-3. (stands) On a reference build that failed 3
  of the 15 two-forest scenarios bitwise at max abs z 0.00; the suite did not see it.
- (stands) On a reference build of the prototype: main corpus 54 of 55 identical and `zeroweights` at
  max abs z 2.85; multinomial 11 of 11; the four snapshot files pass; the 27 exact gates and the
  zero-weight arm pass in `quick`. `zeroweights` is the one recorded scenario that holds a zero weight
  when its sampler is created (80 of 500, the smallest response on one of them: [1.3742, 25.2117] becomes
  [1.6243, 25.2117]).
- (stands; one file ran) The suite on the prototype: 2 failures of 17757 results, both the two-forest pin
  of test-active-rows-pins.R (["Bayesian causal forest, bitwise"](../../inst/tinytest/test-active-rows-pins.R)),
  which compares a sampler created with weights w and then masked against one created with w * a. Rerun:
  54 results in that file, those 2 failing.

Gates that can fail.
- (stands) [`zeroWeightArm`](../../benchmarks/R/bd-balance.R) creates with every row and so proves
  nothing here. With the weights given at creation and the zeroed cells' response 3 higher, about 2
  seconds in `quick`: the tip matches the oracle whose range is over every row (largest abs z 1.4) and
  misses the rows-in oracle by up to 44.5; the prototype matches the rows-in oracle (2.1) and misses the
  other by 40.1.
- (ran) A count twin, from the fixed-shape arm of negbin-exact.R: 4 rows in of 25 a cell, the rows out 30
  higher, left out by 0/1 weights at creation, the oracle's likelihood over the rows in and its centre
  over every row (3.334) or the rows in (1.056). 2 seconds for two seeds. Tip: gap 0.012 to the every-row
  oracle, 0.493 to the rows-in one, tolerance 0.12. Prototype: 0.484 and 0.003.

Released behaviour. (read) 0.9-34 took the range over every row, at creation and at
`setOffset(updateScale = TRUE)`, and had no `updateScale` on `setResponse`; (stands) on its build two
seeded fits that differ only in the responses of the zero-weight rows differ by up to 3.65 in a drawn
fitted value at the rows in. Counts, several forests, masks and Student-t reached no release.

## The rule

1. Rows. The range (gaussian, Student-t, aft, and under a variance forest or any leaf model), the count
   centre and the unit of a forest's sd are computed over the rows in when the sampler is created,
   unweighted: a row of positive weight counts in full.

       dbarts(y ~ x, d, weights = w)    # range of (y - offset)[w > 0]; unit sd((y - offset)[w > 0])

2. Held. No later call moves them unless it carries `updateScale = TRUE` (as today), and none of
   `setWeights`, `setActiveRows`, `setState`, a copy, a reload or a warm start does. No argument is added
   in this slice.
3. Re-derived. `updateScale = TRUE` on `setResponse` and `setOffset`, and `setData`, read the rows in at
   that call. Several forests still refuse `updateScale = TRUE`.
4. Too little to read. Where a row is out and the rows in hold fewer than two distinct values (none in;
   one in; several in at one value), each number is taken over every row, as if no weight or mask had
   been given, sigest included. At creation there is one warning, of its own class, and the sigma
   fallback's warning is not raised for that cause:

       fewer than two distinct values of the response are in the likelihood (<m> of <n> observations are in it), so the response's scale and the prior's centre are taken from all <n>

   The value tested is the response less its offset for a range and the count for a centre. At a
   re-derivation the same fallback applies with no message, which is what such a call does today.
5. Nothing out, nothing changed. With every row in, every number above has the bits it has today.
6. sigest stays `lm`'s weighted estimate over the rows in and is re-estimated by no call. Its fallback sd
   and its floor read the rows in.
7. `xbart`: each fold's range and sigma follow 1, 4 and 6 over the fold's training rows.

## Constraints

- No new argument; no bridge entry changes arity; no facade virtual. `--preclean` on the engine commit,
  for the edited headers.
- Rule 5 is a constraint on the code, not only a test: the unit's loop and the centre's loop run as they
  are written today whenever no row is out, and the rows-out case is a separate path. A reviewer checks
  that neither present loop gained a branch.
- One decision of which rows are read. The response model makes it; the unit follows the response
  model's answer and does not test the rows again.
- The flat C header keeps every signature; the comments of its two re-anchoring entries gain "over the
  rows in the likelihood". The API hash does not move.
- No state format change and no new record: the range is the model's `response.range`, the unit the
  forests' `anchor`, sigest the data object's, as now. A re-creation restates the record and derives
  nothing.
- Cut points, a hazard fit's periods, the number of categories, leaf covariate standardization and a
  multiplier's scale are not touched.
- [`rbart_vi`](../../R/rbart.R) is not touched.

## Steps

"Fails today" is what the tip does where the test expects otherwise.

1. Engine, the range. [`GaussianResponse::rescale`](../../src/bartcore/model.hpp) reads the rows whose
   served precision is positive from the vector slice A has it take; with rule 4's case, every row. The
   constructor, `setResponse`, `setOffset` and `setData` follow with no edit of their own, and so do
   Student-t (its composite precision is zero exactly where the weight or the mask is) and aft (the mask
   only). The response model answers which rows it read, for step 3. tests/cpp, beside
   [`testActiveRowsGaussianDf`](../../tests/cpp/test_model.cpp), two chains where a chain is involved:
   - gaussian, Student-t and aft (aft by a mask then `setOffset(true)`), the weights of the rows in
     fractional: the scale and shift equal the unweighted literal over the rows in (fails today: every
     row);
   - no weight at all, and weights all positive: scale, shift and the working response are bit for bit
     those of the present expression, written out in the test;
   - none in, one in, four in at one value: the every-row literal;
   - held: after `setWeights` and `setActiveRows` the transform is unchanged and 20 sweeps equal an
     untouched twin's under the same weights;
   - `setResponse(true)` and `setOffset(true)` under zero weights and under a mask read the rows in;
     `setData` reads the rows in under the new weights.
2. Engine, the centre. [`NBResponse::computeShift`](../../src/bartcore/model.hpp) takes the mask where
   one is installed; the unmasked call is the present loop. Its two re-deriving callers pass the mask.
   tests/cpp, beside [`testNBLogMeanAnchor`](../../tests/cpp/test_model.cpp): the centre after a mask and
   `setOffset(true)` equals the literal over the rows the mask keeps, with and without an exposure
   offset; unmasked, bit for bit the present value; a mask keeping one row, or rows of one count: the
   every-row value; held across `setActiveRows`. The comments of
   [`NBResponse::setActiveRows`](../../src/bartcore/model.hpp) and
   [`AFTResponse::setActiveRows`](../../src/bartcore/model.hpp) stop saying "full-data".
3. Engine, the unit. [`scaledResponseSd`](../../src/bartcore/chain.hpp) is left as it is and is what
   runs when the response model read every row; otherwise a second routine takes the sample sd over the
   rows it read, n - 1 of those rows in the divisor. tests/cpp, beside
   [`testForestCalibration`](../../tests/cpp/test_sampler.cpp): with rows out the map's anchor equals the
   literal sd to 1e-14; with none out it is bit for bit the present loop's value at n = 64, 100, 137 and
   250 (the size at which the prototype moved a bit among them); a finite carried anchor still wins.
4. R, the rows and the warning. One function in [spec.R](../../R/spec.R), called by
   [`resolveSamplerSpec`](../../R/spec.R) after [`enforceWeightPolicy`](../../R/spec.R), answers rule 4
   for a data object and a mask and raises its warning, for the families that take something from the
   response (gaussian with any residual law, nbinom) and no other.
   [`estimateSigmaFromLinearModel`](../../R/utility.R) fits without the weights in rule 4's case, and
   otherwise hands the rows in to the sparse fallback, to
   [`floorSigmaEstimate`](../../R/utility.R) and to [`floorMarginalSigma`](../../R/utility.R).
   [`foldData`](../../R/xbart.R) does the same within a fold, with no warning per fold.
5. R, creation under a mask. Of the families whose 0/1 weights become the mask
   ([`isMaskedWeightFamily`](../../R/spec.R)) only nbinom takes anything from the response. For it,
   where rule 4's case does not hold, [`dbarts`](../../R/dbarts.R) and
   [`bartcoreSamplerSetData`](../../R/bartcore.R) follow their install of the mask with the sampler's own
   `setOffset(offset in force, updateScale = TRUE, updateState = FALSE)`, before the sampler is handed
   back: it draws nothing, records the range and restates a named sd. A sampler built by hand from
   [`dbartsSpec`](../../R/spec.R) and `new("dbartsSampler", ...)` that installs `spec$active` itself was
   created with every row in; the help of `dbartsSpec` gives the call that makes it `dbarts()`'s. Revised by dec-B363: that call is `setActiveRows(spec$active, updateScale = TRUE)`, written into the help when slice C builds it, with the loop form; the `setOffset` route is not documented.
6. tinytest, a new file test-response-scale-rows.R, fixtures with the rows out extreme in the response:
   - Inert. Two fits that differ only in the response of the rows out at creation are `identical()` at
     the rows in and within 1e-12 at the rows out, in draws, sigma and every reader: gaussian, a fixed
     sigma, Student-t, a variance forest, linear and gp leaves (covariates equal across the pair), a
     named sd, two forests on a factor and on a number, an offset, nbinom, and `bart`, `bartBT` and
     `xbart` once each. Fails today on every arm.
   - Literals. `response.shift`, `response.scale`, `k.scale`, `prior.mean`, the forests' `anchor` and
     sigest equal values written from the rows in, and equal the `subset` fit's.
   - Held. The readers are unchanged across `setWeights`, `setActiveRows`, a copy and a reload, read from
     the engine (`getLeafPrior`) and not from the record alone; five sweeps after a copy agree with the
     original's to 1e-10.
   - `updateScale = TRUE` on `setResponse` and `setOffset` under zero weights and under a mask gives the
     readers of a sampler created under those rows; through the flat entries once
     ([`dbarts_sampler_setOffset`](../../src/C_interface.cpp)); several forests still refuse.
   - Rule 4. Every weight 0; a mask of zeros (nbinom); one row in; rows in at one value: created, exactly
     one warning, counted with `withCallingHandlers`, of the new class, its text as above with m and n;
     the readers and sigest equal those of the fit with no weights; no warning where every row is in and
     the response is constant; none for probit under a mask of zeros.
   - sigest: three rows in and two predictors gives the sd of the rows in with the existing fallback
     warning; the floor reads the rows in.
   - Creation routes for a count fit: `dbarts` with 0/1 weights, `setData` with 0/1 weights, and
     `dbartsSpec` then the documented call all give the rows-in centre; the by-hand route without the
     call gives every row's, pinned as the documented difference.
   - A mask over w against weights w * a: equal when both samplers are created under w, two fits when
     not (the readers differ), stated as such.
   - A state across row sets: a single-forest state from a sampler created under other zero weights
     installs converted, the recipient's range unmoved and the installed fit equal to the donor's to
     1e-12; a gp and a two-forest state are refused by the existing message. The return value of
     `setState` is not asserted (dec-B305 revises it).
7. Restated pins. In test-active-rows-pins.R the two-forest pin creates both samplers under w and masks
   or reweights after; the gaussian and Student-t pins above it get the same repair, since they pass only
   while the fixture's extremes are on rows in. The comments
   ["The log-mean shift is the full-data one by design"](../../inst/tinytest/test-active-rows-pins.R) and
   ["the response transform is the FULL-data one by design"](../../inst/tinytest/test-active-rows-pins.R)
   say "creation's, every row then being in".
8. Exact gates. bd-balance.R gains a third arm beside [`zeroWeightArm`](../../benchmarks/R/bd-balance.R):
   the weights at creation, the zeroed cells' response 3 higher, the oracle's range over the rows in.
   negbin-exact.R gains a third arm beside its fixed-shape one, as measured in Context, with the oracle's
   centre ([`shiftOf`](../../benchmarks/R/negbin-exact.R)) over the rows in.
   [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) runs the bd-balance arm as it runs the
   zero-weight one; the negbin arm runs inside its script. About 4 seconds more in `quick`.
   benchmarks/README.md names both.
9. The baseline. `zeroweights` is recorded fresh on a reference build at the slice's newest code commit
   and the other 54 scenarios are carried from the current file after they compare bitwise; the
   [MANIFEST](../../benchmarks/baselines/MANIFEST) row says so in its predecessors' form and names the
   oracle (step 10's pair script, rows 3 and 4). The workflows that pin the file's name follow.
10. The pair script, benchmarks/R, run on the base and the slice builds with the same seeds, 20 sweeps:
    1. no row out, each family and composition, weights all positive, a mask after creation,
       `updateScale = TRUE` with no row out, `xbart`: `identical()`;
    2. rows out, no extreme among them, one forest: `identical()`;
    3. rows out with an extreme, gaussian and Student-t: the slice's fit equals, at the rows in, the base
       build's fit with the out rows' responses moved inside the range of the rows in: `identical()`;
    4. a count fit whose masked rows are structural zeros, as many out as in, no offset: the slice's fit
       equals the base build's with the masked rows given the counts of the rows in: `identical()`;
    5. several forests with rows out, and each mover of the `rng:` block: must differ, and by how much.
11. Mutations (Verification).
12. Records. Help: the `weights` item of [bart.Rd](../../man/bart.Rd) and [dbarts.Rd](../../man/dbarts.Rd)
    (the sentences holding ["still counts toward"](../../man/bart.Rd) and
    ["the response scale held equal"](../../man/bart.Rd)) say that a row of weight 0 takes no part in the
    response's scale, the count centre or a forest's default sd, that it stays in the design and the fit
    is close to the fit on the remaining rows without being it, that the step is at 0 and not at a small
    weight, and that a mask or weights set after creation leave the numbers creation set, so the three
    ways to leave a row out differ by when it left; the nbinom paragraph names the rows of c;
    [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): the response item
    (["which reads every row"](../../man/dbartsSampler-class.Rd) becomes the rows in), the weights item,
    the `active` item (["stay the FULL-data calibration"](../../man/dbartsSampler-class.Rd) becomes
    creation's), `updateScale`, and that `setData` and a sampler made by hand from another's parts derive
    again; [dbartsSpec.Rd](../../man/dbartsSpec.Rd), the `active` value; xbart.Rd, `weights`. NEWS,
    user-visible changes, one item: "An observation of weight zero no longer takes part in the
    response's scale, which the leaf prior is stated against, so a fit with zero weights has the leaf
    prior of the fit on the remaining rows; 0.9-34 took the scale over every row. The observation is
    still in the design." A design note, docs/design/response-scale-rows.md, with its index row: the
    rule, the inventory, the measurements, the named-sd case (a prior centred at 52 where the rows in sit
    near 2), why membership and not weight, the four slices.
    [active-rows-mask.md](../design/active-rows-mask.md),
    [negative-binomial.md](../design/negative-binomial.md) and [nbinom-log-mean.md](nbinom-log-mean.md)
    take a dated amendment where they say "full-data" or "all rows". dbarts.h's two comments. TODO: the
    entry is cut to slices C and D. The landing note here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, clean under ASan
  and UBSan; the full tinytest suite, expected clean after step 7 (2 failures without it); the new file
  under ASan on the R-loaded path.
- On a reference build, the three compares, counted per scenario with no `max |z|` line but the one
  named:
  - main corpus against `equivalence-1b7d730c.rds`: 54 identical; `zeroweights` differs and is the only
    one. Its statistical verdict against the old baseline is reported, not gated: the posterior changed
    (the prototype's was max abs z 2.85). Then all 55 identical against the new file.
  - two-forest corpus against `bcf-equivalence-1b7d730c.rds`: 15 of 15 identical. Fewer is rule 5 broken
    in the unit's loop.
  - multinomial corpus against `multinomial-equivalence-80b1c8d4.rds`: 11 of 11 identical.
  - the four `test-reproducibility-*.R` files pass unchanged; none is regenerated.
- Every gate of [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) in `quick`, unchanged, and
  the two new arms. Each new arm is also run once on the base build, where it must FAIL (measured: by
  44.5 standard errors and by 0.49 against 0.12), and that run is reported.
- The pair script: rows 1 to 4 `identical()`, row 5 differing.
- Mutations, each to fail the named check:
  - the range over every row: step 6 "Inert" and "Literals"; the bd-balance arm;
  - the centre ignoring the mask: step 2's literal; step 6's nbinom arm; the negbin arm;
  - the unit over every row: step 3's literal; step 6's two-forest arms;
  - the unit's present loop replaced by the rows-in loop when no row is out: step 3's bits; the
    two-forest compare;
  - the centre's present loop given the mask test: step 2's bits (and the main corpus's count scenario
    if a bit moves);
  - weighted, not membership (a row counted by its weight): step 1's literal with fractional weights;
  - rule 4 off (one row in gives the window of width 1): step 1's "one in"; step 6 "Rule 4";
  - rule 4 on with every row in (a constant response warns): step 6's "no warning";
  - creation's re-derivation under the mask dropped: step 6 "Creation routes"; the negbin arm; and run
    for probit too: the suite's seeded probit-under-mask pins, if the call draws;
  - sigest's fallback left on every row: step 6 "sigest"; both warnings raised with no row in: step 6's
    count of one;
  - `updateScale = TRUE` left reading every row: step 1's last bullet; step 6's fourth.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; inst/NEWS.Rd parses with a non-NULL result; `R CMD check --as-cran` on a tarball from a clean
  copy. No new Rd topic.
- Cost: creation and the two response setters gain a comparison per row; no sweep path. No bench compare.

## What a consumer sees

Read in each package's source (weights, `subset`, masks, `updateScale`); no consumer suite was run.

- stan4bart. A gaussian fit hands its weights to the sampler at creation, zeros allowed (it refuses only
  all zero), and then calls the flat `setOffset` entry with `updateScale` true at the start and through
  warm-up. Today both read every row. After this slice creation reads the rows in, and so must every one
  of those calls, or the first would undo creation's scale: that is why rule 3 is in this slice. A
  gaussian stan4bart fit with a zero weight therefore changes its leaf prior and draws; one without is
  bit for bit unchanged. Its binary fits install a mask after creation on a probit sampler, which takes
  nothing from the response. No test of its suite shows the change: test-15-weights_offset_k.R draws its
  weights from (0.5, 2) and from {4, 0.25}, and test-22-binary-weight-mask.R is probit. To verify: its
  suite on the slice's library, unchanged; and one seeded gaussian fit with zero weights on the base and
  slice libraries, differing, beside one without, identical.
- bartCause passes the user's weights and `subset` at creation, to `bcf` too, and changes neither
  afterwards. A fit with a zero weight moves as the `rng:` block says, under `bcf` by the unit as well.
  No test of its suite passes a zero weight (searched).
- treatSens fits a probit treatment model and re-anchors it through the flat entry; a probit sampler
  takes nothing from the response. bairrtt uses no weights, no mask and no `updateScale`.

## Out of scope, and where it goes

- Slices C and D, below.
- `rbart_vi`'s `rel.scale`: left as it is; the function goes in 1.1-0.
- A hazard fit's periods, which are taken from every subject's time, weight 0 included: for the
  maintainer (below).
- A sampler made by hand from another sampler's control, model and data derives the range again under
  the weights in force and keeps the unit from the control: an existing seam, opened further by this
  rule; named in the help here, closed by the record slices of TODO `state-frame-prior`.
- A gp or several-forest state between two samplers on the same data created under different zero
  weights is refused where the tip installs it: kept, tested, and said in the help.

## Slices C and D, not planned here

C. `updateScale` on `setWeights` and `setActiveRows`. Waits on: this slice; leaf-conversions, for the
   conversion of leaf values; answers to the points below. Open, to settle in its plan:
   - What a re-derivation does to the fit on a call that changes no response. Measured by the critic: by
     the response setters' path, which keeps every leaf's internal number, the live fit at the rows in
     moved by up to 15.07 with no data change and `predict` on the same kept draws by up to 14.50; under
     a mask redrawn each sweep the fit was rescaled whenever one extreme row crossed. The call recorded
     here: it should keep the fitted function and restate the prior, converting the live leaves through
     the conversion leaf-conversions lands, not follow the response setters.
   - Position. The new argument goes after `updateState`, so 0.9-34's positional `setWeights(w, TRUE)`
     keeps its meaning.
   - A sampler of several forests: its unit cannot be re-derived today, while dec-B303 lets
     `setForestBasis(updateBasisScale = TRUE)` re-derive a multiplier's scale over the rows in at the
     call; the two would then span different rows.
   - A help sentence, and a test, that `setState` does not undo a re-derivation.
   - dec-B302's refusal of a re-derivation with no row in (rule 4 of this slice keeps today's fallback
     until then), on all four setters at once.
   For the maintainer, plainly:
   - A hazard fit takes its periods from every subject's time, subjects at weight 0 included; a far-out
     subject at weight 0 adds a period and moves the draws at the subjects in. Cut points, or response?
   - sigest can never be re-derived: a sampler created with every row and masked afterwards keeps 50.07
     where the rows in give 0.949, and `updateScale = TRUE` moves the range and leaves it.
   - The several-forests collision above.
D. dec-B304: a linear or gp leaf's covariate centre and scale, and the default lengthscale, over the
   rows in at creation, held, re-derived by `updateStandardization`; `xbart`'s folds taking their own.
   Waits on: the record slices of TODO `state-frame-prior` (a state install today replaces the
   recipient's centre and scale and reports a clean install, and a re-creation without a state derives
   them afresh, so there is nowhere to hold them); `updateStandardization` itself (dec-B233, not built).
   Open: creation under a mask must derive the standardization under it too, which this slice's
   `setOffset` step does not name, so either creation takes the mask or the install re-derives two
   things; the leaf is built from the column store, which holds no weights.

## Calls made in planning

Agent-made, each for the maintainer's later mark. The first seven are the coordinator's rulings on the
critique of the design (2026-10-07); the rest are the planner's.

- Four slices: A the aft re-anchor fix; B creation; C the two new arguments; D dec-B304 with the record
  slices. Not taken: the design's one slice, estimated at 1380 lines before dec-B304 was ruled.
- B is creation only, with `updateScale = TRUE` on the two existing response setters reading the rows
  in. No new argument, bridge arity or facade virtual.
- Fewer than two distinct values among the rows in takes dec-B302's fallback and warning, not the window
  of width 1. Ruled by the maintainer since, as dec-B309; the measurement is in Context.
- With no row out every number keeps its bits, as a constraint on the loops and a gate.
- Two exact arms that fail on the base build.
- For slice C: the argument after `updateState`. A re-derivation through it treats the live chain as
  `setResponse(updateScale = TRUE)` does, ruled by the maintainer as dec-B327; keeping the function
  was the orchestrator's open call.
- `rbart_vi` is cut: deprecated on the tip.
- The value tested for "fewer than two distinct" is the count for a count fit. So rows in that all hold
  one positive count take the fallback although their centre is defined; the alternative is a second
  test for one family.
- The warning is raised only by families that take something from the response. A probit or ordinal
  sampler under a mask of zeros is created in silence, as today; dec-B302's sentence would be false of
  it. The alternative is a warning in every family with other words.
- In rule 4's case sigest is the fit without the weights (25.006 on the fixture) where today it is the
  marginal sd (25.245), and the sigma fallback's warning is not raised: one warning, as dec-B302 says.
- At a re-derivation through the two existing setters rule 4's case reads every row with no message and
  no refusal, as today. dec-B302 says such a re-derivation is refused; a refusal before anything is
  touched needs the bridge to ask the engine which rows are in, a facade read this slice is barred from,
  and stan4bart makes the call every warm-up iteration. Left to slice C. The coordinator decides.
- Creation under a mask is R's install followed by the sampler's own `setOffset`, for nbinom only; a
  sampler made by hand keeps the every-row centre unless its caller makes that call, and the help says
  so. Not taken: the engine's creation taking the mask, which changes the creation entry's arity.
- `xbart` warns once, from the estimate on all rows, and a fold in rule 4's case falls back in silence.
- The unit follows the response model's decision of which rows were read; it does not test again.
- The warning's class is `dbartsScaleFallbackWarning` (with `dbartsWarning`), and its text is rule 4's.
- The negbin arm runs inside negbin-exact.R; the bd-balance arm is a named arm run by the workflow.
- The test of a state across row sets does not assert `setState`'s return value, which dec-B305 changes.
- NEWS names the response's scale only, with what 0.9-34 did.
- Against the design and the critique. Neither said what sigest is in the no-row case, where the tip
  already warns once for another cause, so "one warning" needed the third call above. The aft behaviour
  both report as found is recorded as a property in
  [aft-status-setter.md](../design/aft-status-setter.md). Every measured claim rerun here held: the
  per-chain aft transforms (1.872, 2.098), the unit's last bit at n = 64, the two pin failures, one row
  in giving `k.scale` 0.5 on the prototype.
