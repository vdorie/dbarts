# response-scale-rows, slice D: a leaf's covariate center, scale and lengthscale over the rows in

Status: PLANNED 2026-10-17 against 9e79f4a3 (bartcore; state-frame-prior designed and planned, not built;
forest-multiplier push 2 built on wt/forest-mult-p2, not landed); revised from its blind critique
(scratch/rsrd/crit/crit-report.md, sound after changes; [Revision from critique](#revision-from-critique))
and then to the maintainer's rulings dec-B464 and dec-B465 ([Rulings](#rulings)). No call open.
Implements dec-B304 as dec-B464 confirms it and dec-B465 extends it. Rests on dec-B231, dec-B233,
dec-B237, dec-B296, dec-B302, dec-B303, dec-B309, dec-B329, dec-B330, dec-B364, dec-B457 to dec-B461 and
dec-A217. Builds on state-frame-prior slices 1 to 3 ([state-frame-prior.md](state-frame-prior.md)); the
standing statement is the design note [response-scale-rows.md](../design/response-scale-rows.md), which
this slice extends; rechecked 2026-10-23 against 030b39f0, and against slice 3's fix round at 9740f9c2
([Recheck against state-frame-prior slice 3's fix round](#recheck-against-state-frame-prior-slice-3s-fix-round)).

agent: two pushes, serial. Per push an opus implementer for the engine, the facade, the bridge, tests/cpp
and (push 1) the exact arm; sonnet for the R lines, the tinytest files, the help and the records; one opus
reviewer told to refute, who runs that push's mutants ([Mutants](#mutants)) and reads the diff twice: once
for bits that must not move where every row is in, once for which rows each door reads, under which mask,
and what it refuses or keeps. Tier: changes draws ([Process by risk](README.md#process-by-risk)) for both.
rng, against equivalence-32418b6e, bcf-equivalence-20eed7ed (current since forest-multiplier push 2,
dec-A219) and multinomial-equivalence-80b1c8d4
([RNG class and expected verdicts](#rng-class-and-expected-verdicts)):
- push 1, the rows. POSTERIOR-CHANGING for a sampler with linear or gp leaves at a door that derives the
  standardization (creation, `setData` deriving, `setPredictor(updateStandardization = TRUE)`) while a row
  is out of the likelihood, and for every `xbart` fit with a linear or gp leaf, zero weights or not, each
  fold reading its own training rows. REFUSED where accepted today: a derivation whose rows in hold fewer
  than two distinct observed values of a leaf covariate, a covariate constant over every row included,
  except a gp leaf whose lengthscale is stated. NEUTRAL, bit for bit, everywhere else.
- push 2, the burn-in arguments. NEUTRAL for every existing call (`updateStandardization` is FALSE by
  default on both setters). NEW, defined by an identity: `setWeights(w, updateStandardization = TRUE)` and
  `setActiveRows(a, updateStandardization = TRUE)` leave the sampler where the plain call followed by a
  re-derivation of every leaf covariate over the rows then in leaves it.
No recorded scenario reaches a changed sequence in either push: expected 55, 19 and 11 identical, no
re-record.
window: directly after state-frame-prior slice 3, on top of its records; push 2 after push 1. Serial with any
engine slice in the leaves' derivation, [`Chain::applyNewData`](../../src/bartcore/chain.hpp), the
setData, predictor, weight and mask paths of [`Sampler`](../../src/bartcore/sampler.hpp) or
`ColumnStore::buildFromParent`. Before 1.0-0, and before the merge to main if it can be (TODO
`response-scale-rows`).
budget: lines added, as rechecked 2026-10-23 against 9740f9c2: raw itemization ~2575 (push 1 ~2075, push 2
~500); forecast ~4430 at landing, a fix round included (~3590 and ~840), package code ~720 (~585 and ~135);
stops at 1.5 times the forecast ([Budget and stops](#budget-and-stops)). No ABI event.

## Goal

A linear or gp leaf's covariate center and scale, and the default lengthscale derived with them, are read
from the rows in the likelihood wherever the sampler derives them, as the response's scale already is, and
held everywhere else through state-frame-prior's record. A caller whose rows in change, an outer sampler's
mixture or imputation loop, re-derives them during burn-in through the call that changes the rows, and then
stops, so every kept draw comes from one fixed prior (dec-B465). Where a covariate's rows in hold one value,
a derivation is refused and a re-derivation keeps that covariate's numbers. A fit with every row in and no
constant leaf covariate is bit for bit what it was, `xbart` aside. Cut points still span every row.

## Rulings

Both 2026-10-17, put during planning; the ledger holds the alternatives and their costs.

- dec-B464 confirms dec-B304: the numbers stay taken over the rows in, at creation. Shown that clustered
  rows in give a tiny scale and an arbitrarily wide prior at the rows out (posterior sd near 720 with the
  rows in spread over plus or minus 0.005), against the measured gain where the rows out lie far off, the
  maintainer: "So I think we soldier on. It's a coherent model, and if people try and fit it on data that
  don't support it they can't expect good things."
- dec-B465 adds adaptation in burn-in. Asked "What about treating it like other scales that can be
  'learned' or 'adapted' as the sampler runs, with a fixed value after burn-in? Like we do elsewhere.", and
  shown that rule, the maintainer: "Let's do it." The rule: `setWeights` and `setActiveRows` gain
  `updateStandardization`, as `setData` and `setPredictor` have it (dec-B460), so the code driving the
  sampler re-derives over the current rows in during burn-in and then stops; a re-derivation that finds the
  rows in holding one value of a covariate keeps that covariate's numbers with dec-B364's warning once per
  sampler, the others re-deriving; creation, a `setData` that derives and an `xbart` fold whose rows in hold
  one value of a leaf covariate are refused, there being no number in force to keep. A plain `bart()` fit,
  whose mask never changes, is unchanged by the rule.

dec-B465 supersedes the first versions' open question (every row, warned or not, or the placeholder, for a
covariate the rows in hold at one value) and their call that `setWeights` and `setActiveRows` take no new
argument.

## Context

Marks: (ran) on 2026-10-17 on a shipped build of 9e79f4a3 in scratch/libs/rsrd-plan, scripts
scratch/rsrd/p1.R, p1b.R, p2.R, p3.R and p4.R, seconds each; (critic) by the blind critique on the same
commit, scratch/rsrd/crit/q1.R to q5.R, q4.R rerun here; (0.9-34) on the released build; (read) in the
source. "Emulated" means the rows-in numbers were put into a tip sampler through `setState`, which at the
tip installs a state's standardization; the draws are not slice D's, the prior is.

Today (read), at 9e79f4a3; what the landings since and state-frame-prior slices 2 and 3 as built change is
in [Recheck against bartcore 030b39f0 and state-frame-prior slices 2 and 3 as built](#recheck-against-bartcore-030b39f0-and-state-frame-prior-slices-2-and-3-as-built).
- The derivation is [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp), the mean and sample
  sd of a column's observed values over every row, a constant column centered at its value exactly with
  scale 1 (dec-B358), the placeholder. [`LinearGaussianLeaf`](../../src/bartcore/model.hpp) and
  [`GPGaussianLeaf`](../../src/bartcore/model.hpp) call it from `initialize` and `reinitialize`; the gp
  leaf then takes a lengthscale the model does not supply from
  [`GPGaussianLeaf::medianPairwiseDistance`](../../src/bartcore/model.hpp), the median pairwise distance of
  the standardized values, subsampled to at most 512 rows by a fixed stride, a row missing the covariate
  counted at its standardized value 0, and 1 where the distances are all 0.
- The leaf is initialized in the [`Chain`](../../src/bartcore/chain.hpp) constructor from the column store
  alone, before the response model exists. Weights are there; a latent family's mask is not, being
  installed after creation by [`installWeightsMask`](../../R/bartcore.R). The engine holds the mask in each
  response model and drops it at a data swap, whose row count may differ.
- [`Chain::applyNewData`](../../src/bartcore/chain.hpp) (setData) calls `reinitialize`: every row of the
  new data. The predictor setters regather under the constants in force. The weight and mask setters touch
  no leaf constant.
- A fold of `xbart` is a view, `ColumnStore::buildFromParent` in [data.hpp](../../src/bartcore/data.hpp),
  which copies the parent's every-row constants and serves them through
  [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp): "every fold runs the prior a
  full-data fit would", which dec-B304 overturns.
- Families (ran, p3.R; critic, q4.R rerun): linear and gp leaves are taken under gaussian, Student-t, aft,
  probit, logistic, ordinal, nbinom and hazard.probit; a multinomial sampler and every model of several
  forests refuse them, as does a variance forest. At creation a row can be out only by a weight of 0
  (gaussian, Student-t) or by 0/1 weights that become the mask (probit, ordinal, nbinom, and hazard.probit,
  whose subjects at weight 0 mask their person-period rows); logistic and hazard.logistic weights are
  positive counts and aft and multinomial take none ([`enforceWeightPolicy`](../../R/spec.R)). `xbart`
  refuses 0/1 weights under its latent families (critic, q5.R), so a fold has rows out only by gaussian
  weights.
- The creation warning ([`scaleRowsFallBack`](../../R/spec.R)) fires only under gaussian, Student-t and
  nbinom, the families that read the response; a probit sampler with every weight 0 and a linear leaf is
  created with no warning (critic, q1.R). Forest-multiplier push 2, landed (dec-A219), raises the warning
  under probit and logistic for a multiplier's column from an engine status, in
  [`recordBasisScales`](../../R/dbarts.R) called by the sampler's `initialize`, and names the column in the
  gaussian text, one joint warning through [`dbarts()`](../../R/dbarts.R)'s held warning and its check on
  [`basisScaleFallbacks`](../../R/dbarts.R).
- Under a mask family the sampler's own data object carries no weights: `s$data@weights` is NULL, and
  `s$setData(s$data)` drops the mask (critic, q2.R).
- `setModel` accepts a changed gp lengthscale and the bridge has no path for it (critic, q1.R).
- A constant designated column is created today with the placeholder, under linear and gp leaves alike,
  and pinned so in [test-boundary-inputs.R](../../inst/tinytest/test-boundary-inputs.R) ("a constant
  designated column under linear and gp leaves") and in sections 1 and 6 of
  [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R), which create samplers on constant
  columns to reach a covariate without spread.
- The burn-in loop works today for the response (ran, p4.R): a gp sampler keeping trees, the weights
  redrawn each pass, `setWeights(w, updateScale = i <= 3)` with `run(1L, 0L)` for three passes and
  `run(0L, 1L)` after, saves its draws; a later `updateScale = TRUE` is refused (dec-B330).
- (0.9-34) There are no linear or gp leaves: `exists("linear", asNamespace("dbarts"))` is FALSE. Nothing
  this slice touches existed in the release.

What moves on the tip (ran). Fixture: 200 rows, x1 uniform on rows 1 to 100 and on [5, 15] on rows 101 to
200, which carry weight 0; y = 2 sin(3 x1) + noise; 20 trees, one chain.

| case | tip | rows in (literal, or emulated) |
|---|---|---|
| gaussian, `linear("x1")`: center, scale | 5.346, 5.210 (every row) | 0.518, 0.268; the slope prior per unit of x1 19.5 times wider |
| the same, rmse of the posterior mean to the truth at the rows in | 0.176 | 0.135 emulated; a fit by `subset` on the rows in 0.128 |
| gaussian, `gp("x1")`: kernel width in x1 units (scale times lengthscale) | 5.66 | 0.27 |
| the same, rmse at the rows in | 0.185 | 0.101 emulated; by `subset` 0.112 |
| probit, ordinal, nbinom with 0/1 weights (the mask) | 5.346, 5.210 (ordinal on its own fixture 5.54, 4.99) | the rows in |
| every weight 0 | every row; warned under gaussian and nbinom only | every row, warned under every family |
| the rows in at one value of x1 (p2.R) | every row's numbers, created | refused, unless gp with a stated lengthscale |
| x1 constant over every row | the placeholder, created | refused, unless gp with a stated lengthscale |
| `setWeights(w)`, `setActiveRows(a)`, with or without `updateScale = TRUE`; `setPredictor(3 * x1)` | held | held; re-derived with push 2's argument |
| `setData` of a fresh or edited data object | derived again over every row (x1 + 100: 105.35, 5.21) | per dec-B460 and this slice's rows |

Why the one-value case cannot be read (ran, p2.R, emulated): the placeholder reads the rows out at their raw
distance from the one value, so their fitted values depend on the column's units (the largest draw at a row
out 40.0 with x1 as given, 3.24 with x1 divided by 1000), and every row's numbers put the rows in at a
standardized value of about -0.91, the slope no row in informs entering their fit. dec-B465 refuses the
derivation instead.

The nearest precedent, the numeric multiplier's scale (forest-multiplier push 2, landed as dec-A219, its
"Push 2 as built"): the scale is the column's sd over the rows the response's range was read from, under
probit and logistic the rows the creation mask keeps, every row where none is in, with the fallback reported
by the engine and warned by R under every family; a column constant over those rows refused where a default
would need its spread (dec-B303); taken once, recorded on the data object as `basis.scale`, handed to every
later construction; the creation mask reaches the engine through
`new("dbartsSampler", control, model, data, active = )` and a fifth argument of `bartcore_create`, at the
first creation alone, every re-creation installing it after the install (dec-A222). Slice D
takes the mask route, the warning's route, the refusal's route (an engine error naming the column, worded by
the bridge) and the record pattern, and not the rows (the first call in
[Calls made in planning](#calls-made-in-planning)).

Earlier rulings on the subject (the ledger, grepped for standardization, lengthscale, leaf covariate and the
rows in), and what stands:
- dec-B231 (derived data are data; setData's rule) and dec-B233 (`updateStandardization` on setData and
  setPredictor, three arguments for three parts; the record names): stand; dec-B465 adds the argument to
  the two row setters.
- dec-B304, as dec-B464 confirms it and dec-B465 extends it: the rule built here. Its "until
  updateStandardization is built" is superseded by state-frame-prior slice 3; its "with no row in at
  creation they are taken over every row under dec-B302's warning" stands, a covariate then constant over
  every row being refused by dec-B465.
- dec-B303 (a multiplier's scale over the rows in, held; a column constant over the rows in refused where a
  default needs its spread): stands. Its "no argument is added to setWeights or setActiveRows for it" is the
  multiplier's and stays so; for a leaf, dec-B465 adds one.
- dec-B296 (membership, not weight), dec-B302 (no row in: every row, warned) and dec-B309 (fewer than two
  distinct response values in: every row, warned): stand for the response; for a leaf covariate dec-B465
  replaces the every-row fallback with a refusal, except where no row at all is in.
- dec-B364 (a re-derivation with too little to read keeps the numbers in force and warns once per sampler):
  stands; dec-B465 applies it per covariate at every leaf re-derivation.
- dec-B237 and dec-B330 (a gp sampler holding kept draws refuses a change of standardization or scale):
  stand; every deriving door, the two new arguments among them, keeps the refusal.
- dec-B329 and dec-B328 (slopes converted by the formula at setData and onto or off a single-valued
  column): stand; a derivation converts as today, and a column moved to one value by `setPredictor` under
  held numbers is reached as before.
- dec-B358 (a single-valued covariate centered at its value exactly): stands where a placeholder is still
  made: a gp leaf with a stated lengthscale, and a record carried with no spread. dec-A177's filed floor is
  TODO `leaf-covariate-single-value-off-centre`, left out here.
- dec-B344, dec-A146, dec-A148 (a supplied lengthscale is model): stand.
- dec-B457 to dec-B461 and dec-A217: stand; this slice changes no install and no record check.
- dec-B200: its home for the range is revised by dec-B232 and its conversion on install reversed by
  dec-B418; neither is touched here. dec-B354 (a held coefficient is 1, the basis carries the coding):
  stands; it is the multiplier's, whose column is never centered or divided, and no sampler holds both.
- dec-A106 (sd is stated on the leaf's own parameter, a linear leaf's coefficients on standardized
  covariates): stands; what "standardized" is taken over changes.

The division: center, scale and the heuristic lengthscale are data, derived and then held (dec-B231); a
supplied lengthscale is model. They change only by a call that names them (`updateStandardization` on the
four setters, or `setData` of a data object without the record) and no state installs them (dec-B457,
slice 2).

## What this slice extends in state-frame-prior

As built: slice 1 landed (dec-A222); slices 2 and 3 read on wt/sfp-3 at 90af5380, slice 2 as its
"Slice 2 as built" records it, slice 3 built and not yet reviewed. Names that exist only there are written
without a cite until they land. The builder checks them again on the landed slices first
([Assumptions the builder checks first](#assumptions-the-builder-checks-first)).

- The record: slice 2's `standardization` slot on `dbartsData`, one entry per forest, `list(center = ,
  scale = )` named by the leaf covariates, a scale NA for a covariate without spread, and `lengthscale`
  beside them on a gp leaf whose model supplies none; NULL where no forest's leaves take a covariate. Read
  through `dataStandardization` and held to the leaves by `checkStandardizationRecord` inside
  [`checkDerivedRecord`](../../R/data.R) at creation and `setData`; written by
  [`recordDerived`](../../R/data.R) from [`bartcore_derivedRecords`](../../src/R_interface_bartcore.cpp),
  which returns the engine's [`DerivedRecords`](../../src/bartcore/chain.hpp) (`leafCenters`,
  `leafScales`, `leafLengthscales`). This slice adds no slot and changes no shape; the numbers its derive
  path writes there are the rows-in ones, and push 2 writes the record after each re-derivation. The record
  is why this slice could not land earlier: a copy, a reload and a re-creation take it and derive nothing,
  so numbers read under a mask survive a re-creation that has no mask at construction. Any record lost on a
  re-creation path becomes a silent change of numbers, or now a refusal, so the tests take two levels of
  copy and reload.
- The doors. Construction takes a supplied standardization through the standard factory only (the
  multinomial factory's leaves take no covariate and a record there is refused), as a borrowed field of
  [`SamplerOptions`](../../src/bartcore/chain.hpp) (`SuppliedStandardization`) consumed in the leaf's
  derivation by `suppliedMoments`, not through
  [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp), which still serves only a view's
  parent constants; absent derives. This slice edits the derive arm only, and a carried record is never
  refused for its values. Slice 3's [`Sampler::setData`](../../src/bartcore/sampler.hpp) takes a
  `DataReplacement`, a `DerivedValueChoice` (derive, keep or supply) for the grid, the standardization and
  the range, and the facade virtual takes it as a defaulted last argument;
  [`bartcore_setData`](../../src/R_interface_bartcore.cpp) takes five arguments, the choices as words.
  This slice edits the derive arm, hands it the mask on
  `DataReplacement` and refuses there, before the store is replaced. Slice 3's `updateStandardization` on
  `setPredictor` is the bridge calling `Sampler::rederiveStandardization`, a facade virtual, after an
  accepted update of [`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp) or
  [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp); it re-derives every leaf covariate over
  every row and converts live and kept linear slopes (`Chain::rederiveStandardization`), and a gp sampler
  holding kept draws is refused before anything moves (`standardizationUpdateFromExpression`). Slice 3's
  fix round narrowed it to the replaced columns, through `Chain::heldStandardization`'s supply of the
  numbers in force with not-a-number at the covariates to derive; this slice hands it the rows and gives it
  the keep; push 2 calls the same virtual with null columns, for every leaf covariate. Slice 1's data handle ignores every record, and slice 2 drops a standardization
  for it twice, in the bridge and in the over-store constructor of [`Sampler`](../../src/bartcore/sampler.hpp);
  this slice changes what a fold then derives, its own rows, and refuses there.
- The tests: slice 2's "copy holds it" (test-standardization-record.R); slice 3's "fresh object neutral",
  "arguments win", "new row count", "linear converts", "gp narrowed" and "fallback where derived"
  (test-setdata-records.R), whose fixtures have every row in and one leaf covariate; and slice 1's
  "xbart ignores". Slice 1's "fresh object neutral" of 16 fixtures was built as a two-library compare in
  scratch (18 fixtures), and slice 2's base-build literals likewise (6 fixtures), so no tinytest literal of
  state-frame-prior's records a base-build number for a fixture this slice moves. Each is named below where
  it is extended or where a pin of it legitimately moves.

## The rule at every door

Rows in: a row of positive user weight that no mask switches off, unweighted (dec-B296): a weight of 0.001
counts as a weight of 3 does. "Derived over the rows in" means, per leaf covariate: the center and scale
over the observed values of the rows in, and for a gp leaf whose lengthscale the model does not supply, the
heuristic over those rows' standardized values, a row in that misses the covariate counted at 0 as today.
"Too little" for a covariate means its rows in hold fewer than two distinct observed values: one value,
or none observed.

- A derivation (creation, `setData` deriving, an `xbart` fold): where no row at all is in, every row is
  read, under the creation warning under every family, its text naming the leaf covariates (dec-B304,
  dec-B302). Then, where any leaf covariate has too little over the rows read, the call is refused and
  nothing is created or changed, unless the leaf is gp with a stated lengthscale, which takes the
  placeholder for that covariate. In an `xbart` fold "every row" is the fold's training rows (rule 7).
- A re-derivation (`setPredictor`, and push 2's `setWeights` and `setActiveRows`, with
  `updateStandardization = TRUE`): each covariate with too little keeps its numbers in force, with
  dec-B364's warning once per sampler, the flag shared with the response's keeps; every other covariate
  re-derives.
- `setData(..., updateStandardization = TRUE)` replaces the data, so there is no number in force for the
  new rows: it is a derivation, and refuses (dec-B465 names "a setData that derives").

With every row in and no constant covariate, today's numbers, bit for bit.

| door | rule after this slice | what moves against slice 3's tip |
|---|---|---|
| creation: `dbarts()`, `bart()`, `bart2()`, `dbartsSpec()` and `new()` from a data object without the record | derived over the rows in at creation, the mask a latent family's 0/1 weights resolve to included, handed to construction as forest-multiplier push 2 hands it; too little refused | posterior where a row is out; refused where too little; a warning where none is in under probit, ordinal and hazard.probit |
| creation from a data object carrying the record (`s$data`, `samplePriorPredictive`) | the record, as carried (slice 2); no row is read, nothing is warned or refused | nothing |
| `new("dbartsSampler", control, model, data)` by hand without `active` | every row under a mask family, the mask going on later; `new(..., active = spec$active)` is `dbarts()`'s | nothing; the help gains a sentence |
| `setData(d)`, `updateStandardization = NULL`, `d` without the record | derived over `d`'s rows in: its positive weights, or under a mask family the mask its 0/1 weights resolve to; too little refused, the sampler as it was | posterior where a row is out; refused where too little |
| `setData(d)`, NULL, `d` carrying the record | the record (dec-B458) | nothing |
| `setData(d, updateStandardization = TRUE)` | derived over `d`'s rows in, slopes converted (dec-B329); too little refused; a gp sampler holding kept draws refused (dec-B237) | as two rows up |
| `setData(d, updateStandardization = FALSE)` | the sampler's own (dec-B460) | nothing |
| `setData(s$data, ...)` under a mask family: the sampler's own data object, or an edit of it | `s$data` carries no weights, so the call drops the mask, as today, and a derivation reads every row; the help says to hand back the 0/1 weights on the data object or to call `setActiveRows` after | nothing new for the mask |
| `setPredictor(x, column, ..., updateStandardization = FALSE)`, a matrix or columns; the per-observation sessions (`forceUpdate = "partial"`, [`updatePredictorPerObservationJointly`](../../R/updatePredictorPerObservationJointly.R)), which refuse `TRUE` (slice 3) | held, regathered under the numbers in force | nothing |
| `setPredictor(..., updateStandardization = TRUE)` | the replaced leaf covariates only, re-derived over their new values at the rows in at that call (the weights and the mask in force); every other leaf covariate held (slice 3's fix round); too little keeps, warned once per sampler; gp with kept draws refused | posterior where a row is out; a replaced covariate with too little kept, where slice 3 makes the placeholder |
| `setWeights(weights, updateState = NULL, updateScale = FALSE, updateStandardization = FALSE)` (push 2) | FALSE: held. TRUE: after the install, every leaf covariate re-derived over the rows then in (the new weights, the mask in force); too little keeps, warned once per sampler; gp with kept draws refused before the install; a mask family forwards both arguments to `setActiveRows` | nothing at FALSE; new at TRUE |
| `setActiveRows(active, updateState = NULL, updateScale = FALSE, updateStandardization = FALSE)` (push 2) | the same, over the weights in force and the new mask | nothing at FALSE; new at TRUE |
| `setResponse`, `setOffset`, `setCutPoints`, `setModel` | held | nothing |
| live `setState` | the sampler's own; slopes as stored; gp kept draws under another standardization refused (slice 2) | nothing |
| warm start | the recipient's own; the donor's slopes as stored (dec-B457) | nothing |
| `copy()` of a live sampler; a reload; a dead-pointer `setState`; a copy of a copy and a reload of a copy | re-created from the record (slices 1 and 2), the live chain or the stored state installed (dec-B459), the mask put back through the bridge's mask entry | nothing |
| a re-creation of a data object without the record (saved by a development build, or the slot cleared by hand) | derived afresh from `(control, model, data)`, the weights on the data object read and no mask (dec-B461); too little refused | nothing owed |
| `xbart`, a fold | derived over the fold's training rows of positive weight, by the creation rule; the view supplies no constants; a fold with too little refuses the call | posterior, every fit with a linear or gp leaf; refused where a fold has too little |
| flat C API | no entry creates a sampler, replaces data, sets a predictor, weights or a mask | nothing |

Cut points, the response's range and center, sigest and a supplied lengthscale are not touched by any of
these. A plain `bart()` fit, whose mask never changes, meets only the creation row.

## The refusal's message

Raised by R from the engine's error, which names the forest and the covariate's column (forest-multiplier
push 2's [`BasisScaleError`](../../src/R_interface_bartcore.cpp) route, worded as
[`basisScaleMessage`](../../src/R_interface_bartcore.cpp) words it); the bridge supplies the column's name
and the counts. Three texts:

- linear, or gp with no stated lengthscale, at creation or `setData`:
  `the leaf covariate "x1" holds one value over the rows in the likelihood (100 of 200 observations are in
  it), so its center and scale cannot be taken from them; drop it from linear(), or change the rows in
  (the weights, the 0/1 weights' mask, or subset)` and, under gp, `...; state its lengthscale in gp(), drop
  it from gp(), or change the rows in (...)`. "holds no observed value" where every row in misses it. At
  `setData` the text adds `, or keep the sampler's with updateStandardization = FALSE`, slice 3's keep,
  which has a number in force and is served (slice 3's gp refusal there names FALSE the same way).
- an `xbart` fold: `in a fold of xbart, the leaf covariate "x1" holds one value over the fold's training
  rows in the likelihood, so its center and scale cannot be taken from them; drop it from linear(), or
  change the rows (the weights, subset or n.test)`, with the gp clause under gp.
- every row in, the covariate constant: the first text, its count "200 of 200".

The keep's warning, class `dbartsScaleFallbackWarning`, raised last by
[`warnScaleKept`](../../R/bartcore.R) on the sampler's first keep of either kind, the holder's one flag
(`scaleKeptWarned`, answered through [`firstScaleKept`](../../src/R_interface_bartcore.cpp)). Its text today
names the response only, so a leaf keep takes its own words (added at the recheck): where leaf covariates
were kept, `fewer than two distinct values of the leaf covariate "x1" are in the likelihood, so
its center and scale are kept as they were (said once per sampler)`, several named together with "their
centers and scales"; where one call keeps the response's scale and a leaf covariate (push 2, both
arguments TRUE), one text naming both; a response keep alone keeps today's text. The bridge's answer
carries which covariates were kept, for R to name them.

## The burn-in idiom, as the help states it

dbartsSampler-class.Rd, beside the `updateScale` item, and dbarts-embedding.Rd's "Whose prior is in force":
"A loop whose rows in the likelihood change, such as a mixture component or an imputation step, can let the
response's scale and the leaf covariates' center and scale follow the rows in during burn-in and then hold
them, so every kept draw comes from one prior: pass `updateScale` and `updateStandardization` as `TRUE` on
the call that changes the rows while burning in, and `FALSE` after. A gp leaf that has kept a draw refuses
both, so burn in with sweeps that keep nothing."

    sampler <- dbarts(y ~ x1 + x2, data, weights = w, leaf.prior = linear("x1"),
                      control = dbartsControl(n.chains = 1L, keepTrees = TRUE))
    for (i in seq_len(n.burn + n.keep)) {
      w <- drawMembership(sampler)            # the outer sampler's step
      burning <- i <= n.burn
      sampler$setWeights(w, updateScale = burning, updateStandardization = burning)
      draw <- if (burning) sampler$run(1L, 0L) else sampler$run(0L, 1L)
    }

The help adds: a covariate the rows in hold at one value at such a call keeps its numbers, with one warning
per sampler; a sampler is refused at creation when its rows in hold one value of a leaf covariate.

## RNG class and expected verdicts

Push 1 POSTERIOR-CHANGING and, for the one-value case, REFUSED; push 2 NEUTRAL and NEW, both as the `rng:`
line states.

Expected after each push, on a `--preclean` reference install, counted per scenario with no "max |z|" line:
- equivalence-32418b6e, 55 of 55 identical. The scenarios with leaf covariates are linear, gp,
  leaffactor and leaffactormixed, none weighted and none on a constant covariate, and wtgp, whose weights
  lie in [0.5, 2]; the xbart scenarios fit constant leaves. No scenario derives under a row out, folds a
  leaf covariate or passes the new argument.
- bcf-equivalence-20eed7ed, current since forest-multiplier push 2 (dec-A219), 19 of 19 identical.
  Several forests refuse leaf covariates.
- multinomial-equivalence-80b1c8d4, 11 of 11 identical. A multinomial sampler refuses them (ran).
- The four `test-reproducibility-*.R` files unchanged; none fits a linear or gp leaf.
- exact-gates.yaml at quick: every recorded command as before (43 at slice 1's landing), plus the new arm of
  [linear-exact.R](../../benchmarks/R/linear-exact.R) (push 1) passing. The SBC leaf arms use positive
  weights and are not run.
Any other count is a stop. No baseline is re-recorded. No hot path: the derivation runs at creation,
`setData`, a deriving setter and a fold's build, one pass per leaf covariate; no bench compare.

## Engine and bridge

The engine decides which rows are in, from the weights and the mask it is handed, in one place, reports
where it read every row and which covariates it kept, refuses by a typed error naming the column, and never
warns or names an R type. The bridge hands it the mask, words the refusal and carries the status; R raises
the warnings. Names that exist only in state-frame-prior slices 2 and 3 as built are written without a cite
([What this slice extends in state-frame-prior](#what-this-slice-extends-in-state-frame-prior)).

Engine.
- [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp) takes an optional row set (null for
  every row, today's loop unchanged) and answers whether the rows read held fewer than two distinct
  observed values. [`GPGaussianLeaf::medianPairwiseDistance`](../../src/bartcore/model.hpp) takes the same
  rows, its stride counted over the rows read.
- The leaves' `initialize` and `reinitialize` (slice 2 gave both a `SuppliedStandardization`, which still
  wins over any rows and is never refused) take the row set and, for a derivation, throw an error naming
  the covariate where it has too little, unless the gp leaf's lengthscales are supplied; a re-derivation
  takes which covariates it derives (the replaced columns for `setPredictor`, all of them for the row
  setters), holds the others as they are, keeps a covariate with too little and reports which it kept. Each
  leaf reports, per covariate, whether its numbers came from every row for no row in, beside
  forest-multiplier's [`basisScaleFromAllRows`](../../src/bartcore/chain.hpp), as a field of
  [`DerivedRecords`](../../src/bartcore/chain.hpp) that R reads and does not write to the data object.
- A derivation at `setData` is tested before the swap, so a refusal leaves the sampler as it was: in slice
  3's [`Sampler::setData`](../../src/bartcore/sampler.hpp), under a derived standardization only, the
  incoming values of the leaf columns are read over the incoming rows in before
  [`ColumnStore::setData`](../../src/bartcore/data.hpp) replaces the store; the derivation itself stays in
  [`Chain::applyNewData`](../../src/bartcore/chain.hpp)'s derive arm, handed the row set. A kept or a
  supplied standardization reads no row and refuses nothing.
- [`Sampler`](../../src/bartcore/sampler.hpp) forms the row set once, from the weights and a mask, null
  where no row is out or where none is in, and hands the same bytes to every chain, as slice C decides the
  keep once from chain 0; at construction in `initializeChains`, which both constructors reach (the one
  over a pre-built store with the fold's weights). The mask comes from:
  - at creation, a new borrowed, consumed field of [`SamplerOptions`](../../src/bartcore/chain.hpp), as
    slice 2's supplied standardization rides it;
  - at `setData`, a new borrowed field of `DataReplacement`, which replaces the sampler's copy (null clears
    it) before any chain derives, as the swap drops each response model's;
  - at a re-derivation, the sampler's copy, set by `Sampler::setActiveRows` only once every chain has
    accepted the mask, normalized as the chains normalize it (null and all ones are no mask), and cleared
    by a call that clears it. A re-creation starts with no copy and gains the mirrored mask when R puts it
    back through the bridge's mask entry. The weights in force are chain 0's, read through a new accessor.
- Slice 3's `Sampler::rederiveStandardization`, void and over every row, takes the replaced columns since
  its fix round (null for every one). Push 1 has it form the row set from chain 0's weights and the
  sampler's copy of the mask, hand it to `Chain::rederiveStandardization`, which leaves a replaced
  covariate with too little out of `Chain::heldStandardization`'s not-a-number marks so that it is held,
  and return which covariates it kept; the linear slopes, live and kept, convert as they do there
  (dec-B329). Push 2 calls it with null and adds no engine method. `Sampler::setWeights` is not edited;
  `Sampler::setActiveRows` keeps the copy of the mask (below).
- [`ColumnStore::buildFromParent`](../../src/bartcore/data.hpp) stops gathering the parent's constants
  into a view, so a fold's leaf derives over the fold's rows by the creation rule. Slice 2 serves a carried
  record through `SamplerOptions`, so [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp)
  then serves nothing and goes, with `gatheredMeans` and `gatheredSds`; `isView` stays, read by
  `hasRequantizeSource`.
- Facade: push 1 adds the mask to `DataReplacement` (no signature changes, the struct's layout does) and
  changes the `rederiveStandardization` virtual's return, to the covariates kept, its argument being slice
  3's; `SamplerOptions` and `DerivedRecords` gain a field each. Push 2 changes no facade. `--preclean` on
  both.

Bridge.
- `createHolder` in [R_interface_bartcore.cpp](../../src/R_interface_bartcore.cpp) checks the creation
  mask once, outside the arm for several forests where forest-multiplier push 2 checks it today, and hands
  it to `SamplerOptions` for one forest as it hands it to `AmplitudeSpec` for several; it catches the typed
  error around the single-forest factory as it catches
  [`BasisScaleError`](../../src/R_interface_bartcore.cpp) around the amplitude factory.
- [`bartcore_setData`](../../src/R_interface_bartcore.cpp), five arguments since slice 3, takes the mask as
  a sixth, set on the `DataReplacement` it builds; its registration goes to 6. It copies the new response,
  weights, offset and test offset into buffers of its own and moves them into the holder only once the
  engine has taken the replacement: today it resizes and overwrites the holder's buffers, which the
  sampler borrows, before the engine call, so a refusal raised after that would leave the sampler reading
  the new values or freed memory. It catches the engine's typed error and raises the refusal after.
- The refusal: the engine's error, worded with the column's name and the counts, raised as an R error; the
  fold's build in [`bartcore_createFromHandle`](../../src/R_interface_bartcore.cpp) raises the fold text.
- The leaf's fallback flags reach R through [`bartcore_derivedRecords`](../../src/R_interface_bartcore.cpp).
  The predictor entries ([`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp), the matrix, and
  [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp), the columns), which since slice 3 call
  `rederiveStandardization` after an accepted update behind slice 3's gp refusal, the columns entry passing
  its replaced columns since slice 3's fix round, carry beside the accepted logical they return (as an attribute R reads and strips, so the
  value a caller gets is unchanged) "warn now" on a sampler's first keep and which covariates were kept,
  from the holder's one `scaleKeptWarned` flag, which the four response and row setters already share
  (response-scale-rows-c rule 5), through [`firstScaleKept`](../../src/R_interface_bartcore.cpp)
  generalized to a keep of either kind.
- Push 2: [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp), three arguments each, take a fourth. With
  it true the gp kept-draws refusal joins
  [`refuseRowsRederivation`](../../src/R_interface_bartcore.cpp) before the install, in the words slice 3's
  `standardizationUpdateFromExpression` gives it (dec-B237; the pattern of
  [`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp)), and never the refusal of several
  forests there, which have no leaf covariate and take the argument as nothing to do; after the install and
  [`rederiveAfterInstall`](../../src/R_interface_bartcore.cpp)'s scale where asked they call
  `rederiveStandardization` with null columns and fold its keep, with what was kept, into the answer they
  return.
- No other entry changes arity; no flat entry, [dbarts.h](../../inst/include/dbarts/dbarts.h) or token
  moves.

R.
- [`bartcoreSamplerSetData`](../../R/bartcore.R) passes the mask it resolves from the new weights to the
  bridge as the sixth argument, before the engine swaps, and still installs it afterwards through
  [`installWeightsMask`](../../R/bartcore.R) as slice 3 does.
- The creation warning for no row in: where the engine reports a leaf covariate read from every row, the
  sampler's `initialize` raises it from [`createEngine`](../../R/dbarts.R)'s records, class
  `dbartsScaleFallbackWarning`, under every family, in forest-multiplier push 2's form: under gaussian,
  Student-t and nbinom one joint text naming the response, the prior's center and the leaf covariates,
  through [`dbarts()`](../../R/dbarts.R)'s held warning, whose check on
  [`basisScaleFallbacks`](../../R/dbarts.R) extends to the leaf's flag (`bart()` and `bart2()` reach it
  through `dbarts()`; `dbartsSpec()` and then `new()` by hand raise the spec's and the leaf's apart, as for
  a multiplier); under probit, ordinal and hazard.probit "no observation is in the likelihood (0 of n have
  a positive weight), so the center and scale of x1 are taken from all n". Several forests refuse leaf
  covariates, so a multiplier's column never joins the text. `setData` raises the same, last, where the
  standardization was derived, as slice 3 raises the range's only where the range was derived.
- The setters raise [`warnScaleKept`](../../R/bartcore.R)'s warning, with the leaf's words
  ([The refusal's message](#the-refusals-message)), last on a first keep, after the mirror, the record and
  the store; [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) strips the bridge's attribute.
- Push 2: `setWeights` and `setActiveRows` in [dbarts.R](../../R/dbarts.R) take `updateStandardization`
  after `updateScale`, checked by one helper with slice 3's check of `setPredictor`'s argument ("must be
  TRUE or FALSE", as [`checkUpdateScale`](../../R/bartcore.R) checks `updateScale`), and write the
  standardization after a re-derivation through the records reader without the grids, as
  [`recordResponseRange`](../../R/data.R) writes the range (dec-A222's cost for a call made every
  iteration), not through [`recordDerived`](../../R/data.R); `setWeights` forwards both to
  `setActiveRows` on a mask family and keeps its rollback on a refused call. `reapplyActiveRows` passes
  FALSE for both.
- Creation needs nothing else past forest-multiplier push 2: `dbarts()` already hands `active` to `new()`,
  which hands it to the bridge.

## Constraints

- With every row in and no leaf covariate with too little, every derivation runs today's code on today's
  inputs: no such fixture moves a bit, state-frame-prior's two-library compares' fixtures included (none
  of its tinytest pins records a base-build number); `xbart` with a leaf covariate is the one exception,
  its folds reading their own rows.
- The row set is membership, never the weights' sizes.
- A refused call creates nothing or changes nothing: the data object, the mirrors, the record and the
  stored state are as they were.
- No new slot on `dbartsData`, no state block, no state format change; the engine's new fields are the
  mask on `SamplerOptions` and `DataReplacement` and the no-row flag on `DerivedRecords`. The new arguments
  are push 2's two, and `bartcore_setData`'s sixth; the new texts are the refusal's three with `setData`'s
  clause, the leaf's part of the creation warning and the leaf's words of the keep's warning, in the
  existing class; the keep's warning reaches `setPredictor` and, for a leaf, the two row setters.
- Cut points, the response's transform, sigest, a supplied lengthscale and every install path are not
  touched.
- No subprocess in tinytest; base-build comparisons by literals written into the test or a two-library
  compare in scratch recorded in the landing note.
- Out of scope: TODO `leaf-covariate-single-value-off-centre` (the floor); weights on the grid; kept gp
  draws that store their kernel (dec-B237's later option); the mask riding the sampler's own data object
  (the TODO entry at the end); `setModel` accepting a changed gp lengthscale it does not honor (dec-B254
  calls it a parameter; if it is ever honored, a change to the heuristic becomes a door that derives);
  stan4bart, whose binary mask goes on after creation, so a linear leaf passed through its `bart_args`
  reads every row, as today, and is refused only on a covariate constant over every row.

## Steps

Estimates are the raw itemization of lines added, as rechecked 2026-10-23; [Budget and stops](#budget-and-stops)
gives the forecast.

### Push 1: the rows, the refusal and the keep

1.1 Engine, the derivation (~85). The row set and the too-little answer in
    `standardizationMomentsForColumn`; the gp heuristic over the rows read; the leaves' `initialize` and
    `reinitialize` take the row set and which covariates they derive, the others held, beside slice 2's
    supplied standardization; they throw the typed error at a derivation (the stated-lengthscale gp exempt),
    keep at a re-derivation, and record the no-row flag in `DerivedRecords`.
1.2 Engine, the doors (~100). `SamplerOptions` gains the creation mask; `Sampler` forms the row set in
    `initializeChains`; `DataReplacement` gains the mask, and slice 3's `Sampler::setData` replaces the
    sampler's copy, tests a derived standardization before `ColumnStore::setData` replaces the store and
    hands the row set to [`Chain::applyNewData`](../../src/bartcore/chain.hpp)'s derive arm; the sampler's
    copy of the mask by its lifecycle above, and chain 0's weights through an accessor; slice 3's
    `rederiveStandardization` (sampler, chain and facade virtual), which takes the replaced columns since
    its fix round, reads the row set, keeps and reports; the view's parent constants go with `ColumnStore::suppliedStandardization`,
    `gatheredMeans` and `gatheredSds` (~30 removed); the fold's build raises the refusal. `--preclean`.
1.3 Bridge (~85). In `createHolder` the creation mask checked once and handed to `SamplerOptions`, and the
    typed error caught around the single-forest factory; `bartcore_setData`'s sixth argument and its
    registration, its copies moved into the holder only after the engine accepts, and the typed error
    caught; the refusal's wording at creation, `setData` (with its FALSE clause) and the fold; the no-row
    flags in `bartcore_derivedRecords`'s list; the predictor entries carry "warn now" and the kept
    covariates as an attribute (their columns reach `rederiveStandardization` since slice 3's fix round);
    `firstScaleKept` generalized.
1.4 R (~110). `bartcoreSamplerSetData` passes the mask; the creation warning for no row in raised by
    `initialize`, joint through `dbarts()`'s held warning, and `setData`'s where it derived, under every
    family, folded into forest-multiplier's text; `setPredictor`'s keep warning in the leaf's words, the
    attribute stripped; the docstrings of `setData` (the mask `s$data` does not carry; the refusal and
    FALSE as its remedy) and of `setPredictor`, which since slice 3's fix round derives for the replaced
    columns over their new values, gaining the rows in and the keep.
1.5 tests/cpp (~315). `testLeafStandardizationRows` beside the leaf tests: the moments over a row set
    equal literals; null and all-ones row sets equal today's bits; one row in, a column constant over the
    rows in, a column constant over every row and NA-only rows in each answer too little; no row in reads
    every row and sets the flag; a derivation with too little throws naming the column, a stated-lengthscale
    gp takes the placeholder, a re-derivation keeps and reports; the heuristic over more than 512 rows in,
    interleaved with rows out and with an NA among them, equals the literal of the stride over the rows
    read. A chain created with weight-0 rows and one created with the same numbers supplied through slice
    2's path draw identically for 20 sweeps; the same with the creation mask under probit and nbinom.
    `setData` with a mask, and with a mask after `setActiveRows` at another row count: the numbers over
    the rows the new mask keeps; a refused `setData`: the sampler bit for bit an untouched twin. A predictor
    update after `setActiveRows`: the mask in force, not the creation one. `rederiveStandardization` given
    column j alone with a row out: j over the rows in, k bit for bit as held (slice 3's fix round pins the
    hold with every row in). The view pins of test_data.cpp ("subset view gathers a leaf covariate through the
    column map", and `gatheredRawAgrees`'s "the constants are the parent's full-data ones") and
    `testSuppliedStandardization`'s view case, whose comment says a view takes the parent's constants, turn
    into: a view supplies no constants and a fold's leaf derives over its own rows.
1.6 tinytest (~780 new, ~280 edited). A new file test-leaf-standardization-rows.R ([Tests](#tests)). Edits
    where a pin creates a sampler on a constant leaf covariate, now refused: test-boundary-inputs.R's
    constant designated column, refused under linear and gp and created under gp with a stated
    lengthscale; sections 1 and 6 of test-leaf-conversions.R; in test-standardization-record.R (slice 2)
    `twoFit` on `flat` and the warm start of `flatDonor` into `flatRecipient`; test-warm-start-forms.R's
    single-value block (`withX2`, `slopedOn`, `linear23`); test-calibration-prior-draws.R's linear sampler on
    the constant `x4`; test-setdata-records.R's `flat` gp sampler on a constant `b` (slice 3's fix round,
    "an NA scale compared as the leaf divides"). These six files are every file in which the suite derives
    a standardization on a column constant over every row at 9740f9c2, by a run of the suite on a build
    that logs each such derivation (scratch/rsrd1/probe-src; creations, `setData` derivations and the
    view's gather, none from `setPredictor`), and no creation in it leaves a leaf covariate at one value
    over its rows of positive weight (a mask's rows were not logged). Each reaches a covariate without spread through a carried record whose scale is NA,
    or a stated-lengthscale gp, instead, keeping what it pins. Edits where a pin reads a number this slice
    moves: the `checkDataHandleViews` callers of [leafPriorChecks.R](../../inst/common/leafPriorChecks.R)
    and the xbart arms of test-linear-leaves.R and test-gp-leaves.R where they pin one. test-setdata-records.R's
    unforced `setPredictor(rep(0.5, n), "b", ..., updateStandardization = TRUE)`, where accepted, now keeps
    and warns once. No tinytest literal of state-frame-prior moves (its base-build comparisons are
    two-library compares in scratch); slice D runs its own, scratch/sfp2/compare.R's six fixtures (`flat`
    now refused) and fixtures with a row out, against its base, recorded in the landing note.
    test-gp-leaves.R's zero-weight sampler keeps rows in, so it warns as before; its pins are finiteness
    and stand.
1.7 The exact arm (~80). [linear-exact.R](../../benchmarks/R/linear-exact.R) gains an arm with rows at
    weight 0 on the existing four values in another mix (so the tree space does not grow), sized so that
    the oracle standardized over every row fails it; the oracle standardizes over the rows in. Its header
    says what it does not see: the gp heuristic, the masks, the refusal and the keep, which the literals of
    1.5 and 1.6 hold.
1.8 Help and NEWS (~110). man/dbarts.Rd's `leaf.prior` item ("standardized internally" over the rows in at
    creation, held; refused where the rows in hold one value of a covariate, a constant column included,
    unless a gp lengthscale is stated); dbartsSampler-class.Rd's `setData` (the mask `s$data` does not carry
    under a mask family; the refusal and FALSE), `setPredictor` (the replaced columns only, where slice 3's
    item says the values the update leaves; the keep and its warning), the `updateStandardization` items
    and the manual's table row for the standardization; man/xbart.Rd (each fold its own rows, and the
    refusal); man/dbartsSpec.Rd. NEWS: a linear or gp leaf reads its covariates' center, scale and
    lengthscale from the rows in the likelihood; a constant leaf covariate is refused.
1.9 Records (~150). The design note [response-scale-rows.md](../design/response-scale-rows.md): rule 9 for
    the leaf, section 2's table, section 4 for the leaf's no-row flag and keep, section 8's slice D;
    [state-frame-prior.md](../design/state-frame-prior.md)'s "What response-scale-rows slice D gets" closed
    by a dated line; TODO; the ledger entry for the calls made; the indexes; this plan's landing note.

### Push 2: the burn-in arguments

2.1 Engine (~5) and tests/cpp (~25). No engine method or virtual: the row setters reach push 1's
    `rederiveStandardization` with null columns. tests/cpp `testRederiveStandardization`: after
    `setWeights` or `setActiveRows`, the call leaves the sampler bit for bit where a twin given a predictor
    update of every leaf column in force followed by `rederiveStandardization` of those columns leaves it;
    a keep reports the covariate and moves none of its numbers.
2.2 Bridge (~35). The fourth argument of `bartcore_setWeights` and `bartcore_setActiveRows`, their
    registration, the gp refusal in `refuseRowsRederivation` before the install and outside the refusal of
    several forests, the call after `rederiveAfterInstall`, the folded answer with what was kept.
2.3 R (~45). The two methods' argument, its check (one helper with `setPredictor`'s) and position, the
    forwarding on a mask family, the standardization written through the reader without the grids beside
    `recordResponseRange`, the store, the warning last; `setWeights`' rollback; `reapplyActiveRows` passes
    FALSE for both.
2.4 tinytest (~290). In the new file:
    - the identity: for gaussian, Student-t, probit, ordinal, nbinom and hazard.probit, linear and gp with no
      kept draw, the call with `updateStandardization = TRUE` against the plain call and then
      `setPredictor` of every leaf column in force with the argument: `identical()` in the record, the stored
      state and 5 sweeps of 2 chains; with kept linear draws, `predict` unchanged to 1e-12;
    - both arguments TRUE: the range and the standardization each the literal over the rows in; a keep of
      both, one warning naming both;
    - the argument: absent and FALSE move nothing; positional `setWeights(w, TRUE)` still names
      `updateState`; `updateStandardization = NA` and `= 1` refused with nothing installed; on a sampler of
      several forests and on one with no leaf covariate, TRUE accepted and nothing moved;
    - the keep: rows in holding one value of x1 and spread of x2: x1 kept, x2 re-derived, one warning in the
      leaf's words naming x1, none on a second keep or after a keep through `setOffset`;
    - gp with a kept draw: refused before the install, the weights and the mask as they were;
    - the burn-in loop of the help, by `setWeights` and by `setActiveRows`: during burn-in the record is the
      literal over each pass's rows in, after it the last burn-in pass's, and every kept draw replays under
      it; the loop equals a twin written with the identity's calls; a copy and a reload after the loop hold
      the record.
2.5 Help and NEWS (~50). The two usage lines and the `updateStandardization` item on all four setters; the
    burn-in idiom and its example ([The burn-in idiom, as the help states it](#the-burn-in-idiom-as-the-help-states-it));
    dbarts-embedding.Rd's paragraph; NEWS: `setWeights` and `setActiveRows` take `updateStandardization`.
2.6 Records (~50). The design note's rule; TODO drops `response-scale-rows`; the ledger entry; the indexes;
    the landing note.

## Tests

Push 1's new file, fixtures as slice B's (rows out extreme in the covariate, no missing values but where
named), two chains. Every literal is written in the test from R over the rows named (`mean(x[kept])`,
`sd(x[kept])`, the heuristic's stride computed in R); "fails today" on every arm where a row is out or a
covariate has too little. Warnings are counted with `withCallingHandlers`, by class and text; refusals by
text, then the readers and 3 sweeps against an untouched twin.

- At creation, each family that takes leaf covariates: gaussian and Student-t by weight 0, probit, ordinal,
  nbinom and hazard.probit by 0/1 weights, linear and gp: the record and the engine's reader equal the
  literals over the rows in to 1e-12; logistic and aft, every row.
- Fractional weights 0, 0.001 and 3 under gaussian: the unweighted literals over the rows of positive
  weight.
- The subset twin: a sampler created with rows at weight 0 has the record a sampler created by `subset` on
  the rows in has, `identical()`.
- Neutral: every row in and every covariate spread, each leaf kind and family: the record equal to the
  base build's, written as literals to 1e-12; the draws by the two-library compare in scratch recorded in
  the landing note, as state-frame-prior's slices found draws recorded on arm64 do not hold on the other CI
  platforms.
- No row in, each family that can have one (gaussian, Student-t, probit, ordinal, nbinom, hazard.probit):
  every row, exactly one warning of the class, its text naming the leaf covariates; under gaussian and
  nbinom the one joint text, not two, through `dbarts()` and `bart()` (by hand through `dbartsSpec()` and
  `new()` the spec's and the leaf's apart, as for a multiplier, pinned so); then a covariate constant over
  every row refused.
- Refused at creation, with the text: rows in at one value of x1 under linear and under gp; x1 constant
  over every row; x1 missing on every row in; one row in; at `dbarts()`, `bart()`, `dbartsSpec()` with
  `new()`, and `new()` with `active`; nothing created. Created: the same under gp with a stated lengthscale,
  its record carrying the value and NA, and the next 5 sweeps finite.
- NA among the rows in: the observed rows in only for the moments; the heuristic counting them at 0.
- By hand: `new(control, model, data, active = spec$active)` is `dbarts()`'s under probit and nbinom with a
  linear leaf, record and 5 sweeps; without `active`, every row (the help's sentence).
- Held: `setWeights`, `setActiveRows`, each with and without `updateScale = TRUE`, `setResponse` and
  `setOffset` with `updateScale = TRUE`, `setPredictor` without the argument on a matrix and on columns,
  and a per-observation session: the record `identical()` and 5 sweeps equal a twin's.
- `setPredictor(x[, j], j, updateStandardization = TRUE)`:
  - after `setWeights(w2)`: the literals over the rows `w2` keeps for j, live and kept linear predictions
    unchanged to 1e-12 (dec-B329, slice 3's "linear converts");
  - a second leaf covariate k held, `identical()`, while j re-derives over the rows in (slice 3's fix round
    pins the hold with every row in);
  - after `setActiveRows(a)` under gaussian and under probit: the rows `a` keeps;
  - after `setActiveRows(a)` then `setData` to a data object of another row count: the new data's rows in,
    and no read past them;
  - on a `copy()` of a sampler holding a mask: the mask put back on the copy;
  - no row in, and j at one value over the rows in: kept, exactly one warning on the sampler's first keep,
    in the leaf's words naming j, none on a second, none after a keep through `setOffset` (the shared
    flag); a response keep's text as today; FALSE holds.
- `setData`: a fresh data object with weight-0 rows, and one with 0/1 weights under nbinom: its rows in; a
  carried record: held, also one whose scale is NA; TRUE: rows in; FALSE: the sampler's (slice 3's
  "arguments win" over a row out); a fresh object whose rows in hold one value of x1, at NULL and at TRUE:
  refused, the sampler as it was, the text naming `updateStandardization = FALSE`, and at FALSE served;
  under probit, `setData(s$data, updateStandardization = TRUE)`: the mask
  dropped and every row read, as the help says.
- The record: after creation with rows out under a mask, linear and gp: `copy()` of a live sampler, a copy
  of that copy, `saveRDS` and `readRDS` with a first use, a reload of a copy, and `new()` from `s$data`: the
  record, the heuristic lengthscale among it, and 5 sweeps equal the source's (slice 2's "copy holds it",
  with numbers no derivation at re-creation could give).
- gp with kept draws: `setPredictor(updateStandardization = TRUE)` and `setData(updateStandardization =
  TRUE)` refused, nothing moved (slice 3's "gp narrowed" stands).
- xbart: through the data handle, a fold with weight-0 rows: the fold's numbers are the literals over its
  training rows of positive weight; a fold whose training rows in hold one value of the covariate: the
  call refused with the fold text; a full-rows view without weights still matches the raw path bitwise
  (the `checkDataHandleViews` pin); `xbart(s$data)` equals `xbart` of a fresh object (slice 1's "xbart
  ignores").

Push 2's arms are its step 2.4.

## Mutants

Each applied alone to a clean export, installed with `--preclean`, `touch` after revert.

Push 1:

| mutant | killed by |
|---|---|
| the row set ignored (every row) | the creation literals; the exact arm; the subset twin |
| the weights' sizes read (a weighted mean) | the fractional-weight arm |
| the creation mask not handed to construction | the probit, ordinal, nbinom and hazard.probit literals; by hand with `active` |
| the mask read and the weights ignored | the gaussian and Student-t literals |
| the predictor derivation reading the creation mask, not the mask in force | setPredictor after `setActiveRows` |
| the sampler's copy of the mask not replaced at setData | setPredictor after `setData` to another row count; tests/cpp under ASAN |
| the copy set before every chain accepts, or all ones kept as a mask | tests/cpp's refused mask; the all-ones arm |
| the placeholder made at creation instead of the refusal | the creation refusals; test-boundary-inputs.R |
| the stated-lengthscale gp refused too | its creation arm |
| a creation refused where only rows out are constant | the creation literals (rows out at one value, rows in spread) |
| setData refusing after the swap | the refused-setData twin |
| one row in refused where the response falls back to every row (the call before the 2026-10-23 ruling) | the one-row-in creation arms (fix round) |
| no fallback where no row is in | tests/cpp's flag or the no-row literal |
| the leaf's no-row fallback not warned under probit | the no-row arm under probit and ordinal |
| two warnings where response and leaf both fall back | the joint-text count |
| a replaced covariate with too little re-derived (the placeholder) instead of held | its keep arm; tests/cpp's keep report |
| `bartcore_setData` adopting its copies before the engine call, as today | the refused-setData twin under ASAN, the sampler's response read after the refusal |
| setPredictor's keep re-deriving over every row, or warned on every keep | its keep arm; the second-keep count |
| the keep's flag not shared with the response setters | the `setOffset` then `setPredictor` count |
| a leaf keep raising the response's text | the keep arm's text |
| `setData`'s refusal without its FALSE clause, or FALSE refused too | the `setData` refusal arm |
| the heuristic over every row while the moments read the rows in | the gp literal |
| the heuristic's stride counted over every row's index | tests/cpp's 512-row case |
| the heuristic dropping NA rows in | tests/cpp's NA case |
| the view still supplying the parent's constants | the xbart fold literal; tests/cpp view pin |
| a fold with too little run instead of refused | the xbart refusal |
| setData deriving with no mask, or with the old one | the nbinom setData literal |
| a re-creation deriving instead of taking the record, at any level | the record arms |
| a carried record with an NA scale refused | the carried-record arm |
| NA rows counted among the rows in | the NA arm |

Push 2:

| mutant | killed by |
|---|---|
| the argument dropped in the bridge | the identity; the burn-in loop |
| the re-derivation before the install | the identity under weights that move rows |
| the keep re-deriving the kept covariate, or keeping all | the keep arm's x1 and x2 |
| the gp refusal after the install | the refused arm's mirrors |
| the argument placed before `updateScale`, or `NA` let through | the argument arm |
| the argument refused on several forests, with the scale's refusal | the several-forest arm |
| a joint keep raising two warnings, or the response's text alone | the both-arguments arm |
| a mask family's `setWeights` not forwarding it | the probit and nbinom identity |
| the record not written after a re-derivation | the burn-in loop's copy and reload |
| `reapplyActiveRows` re-deriving | the record arms after a reload |

Expected equivalent, said in the landing note if alive: the row-set loop taken with an all-ones set where
none is out (same sums in the same order).

## Verification

Per push, against its own library:

- tests/cpp plain and under ASAN and UBSAN; `R CMD INSTALL --preclean`; the full suite; the R-loaded address
  sanitizer over the new file, test-gp-leaves.R, test-linear-leaves.R, test-leaf-conversions.R,
  test-boundary-inputs.R, test-data-handle.R, test-standardization-record.R, test-setdata-records.R,
  test-warm-start-forms.R and test-calibration-prior-draws.R; the x86-64 leg (tests/cpp under the
  sanitizers and the suite) before landing.
- The two-library compare in scratch against the push's base (1.6), every fixture with every row in and
  every covariate spread identical, recorded in the landing note.
- Reference build: the four snapshot files and the three bitwise compares, counted, as expected above.
- exact-gates.yaml's commands at quick, push 1's new arm among them; the arm's every-row mutant run once to
  show it fails.
- Lint, air, rc-codoc, win-drift and doc-freshness, each on its own exit status; inst/NEWS.Rd parses with a
  non-NULL result; `R CMD check --as-cran` from a clean tarball.
- Consumers read, not run: bartCause fits several forests (refused leaf covariates), stan4bart's mask goes
  on after creation, treatSens and bairrtt pass no linear or gp leaf (searched).

## Budget and stops

Lines added, removals reported beside them, the baseline files aside; as rechecked 2026-10-23.

| part | push 1 raw | push 1 forecast | push 2 raw | push 2 forecast |
|---|---|---|---|---|
| engine and facade | ~170 | ~270 | ~5 | ~10 |
| bridge | ~85 | ~135 | ~35 | ~55 |
| R | ~110 | ~175 | ~45 | ~70 |
| tests/cpp | ~315 | ~630 | ~25 | ~50 |
| tinytest | ~1060 | ~1905 | ~290 | ~520 |
| exact arm | ~80 | ~130 | - | - |
| help and NEWS | ~105 | ~140 | ~50 | ~65 |
| records | ~150 | ~195 | ~50 | ~65 |
| total | ~2075 | ~3590 | ~500 | ~840 |

Package code: push 1 ~365 raw, ~585 forecast; push 2 ~85 raw, ~135 forecast (rechecked against 9740f9c2;
the figures before it were ~2110, ~3650, ~370 and ~590). dec-B465 estimated the
arguments at 50 to 100 lines beside slice D; push 2's package code is inside that and its tests beyond it,
the identity and the loop being where its defects would show.

From the landed history: slice B planned 2065 lines and landed 3830 (1.85 times); slice C planned 1120 and
landed 2537 (2.27), its package code 165 against 280 and its tests 710 against 1530; state-grid-install
1475 against 2092 (1.42). Since the plan: forest-multiplier push 2 planned ~2950 and landed 3526 added and
1110 removed with fix rounds of 65 and 84 (1.25 added, 1.57 added and removed), push 1 ~2300 against 2016
and 871, push 3, R surface only, ~520 against about 1240 (2.4); warm-start-salvage slice D ~950 against 1003
and a fix round of 233 (1.3); state-frame-prior slice 1 ~1860 against 2026 added and 623 removed at its
landing, the fix round of about 350 included (1.3 to 1.6 by the coordinator's count), slice 2 ~1100 against
1480 added and 468 removed before its review (1.35), and slice 3 ~1650 against 1656 added and 234 removed
before its record and review (its engine 306 against 230, 1.33; its tests 932 against 1000). In lines
added those ran 0.9 to 2.4 times their plans, most between 1.25 and 1.85, the low ones before review.
Tests are forecast at about 1.8 to 2 times, package code at 1.6, the rest at 1.3, about 1.73 overall, the
upper half of that, and kept: the coordinator's reading of the same history is 1.5 to 2. The forecast is at
landing, a fix round included; the range is 3200 to 4200 for push 1 and 750 to 1000 for push 2 (1.5 to 2
times the raw figure).
Against the first version (~1460 raw): the critique's warnings and mask lifecycle add about 370, dec-B465's
refusal and the rewritten constant-column tests about 200, push 2 about 500. The recheck adds about 80 to
push 1 (the keep's words and the joint creation warning in R, three more test files on constant
covariates, the narrowing pins, the larger removal in the view) and takes about 20 from push 2's engine,
the virtual being slice 3's.

Stop and report, building nothing further, when:
- push 1's diff passes ~5390 lines added or its package code ~880 (~5450 and ~890 before the recheck
  against 9740f9c2); push 2's passes ~1250 or ~200;
- a bitwise compare prints "max |z|", a snapshot moves, or an exact gate other than the new arm moves;
- a fixture with every row in and every covariate spread moves a bit outside `xbart`, or a pin of
  state-frame-prior's tests moves other than one named in 1.6 or one with a linear or gp leaf and a row
  out;
- a test outside those named in 1.6 fails because it created a sampler on a constant leaf covariate (list
  it, do not edit it);
- the derivation needs a line in a response model, a state block, a flat entry or dbarts.h, or a facade
  change beyond push 1's mask on `DataReplacement`, its fields on `SamplerOptions` and `DerivedRecords`
  and the `rederiveStandardization` signature, or any facade change in push 2;
- state-frame-prior's record or doors differ from
  [Assumptions the builder checks first](#assumptions-the-builder-checks-first) in a way that moves a step
  by more than a quarter of its estimate.

## Calls made in planning

Agent-made, each for the maintainer's later mark.

- The rows are the rows in by membership, every row where none is in. Taken from dec-B304's own words, and
  one rule at every door, those that read no response among them. Not taken: forest-multiplier push 2's
  rows, those the response's range was read from, which differ only where the rows in hold several rows but
  one response value, and which at a setter would need membership anyway.
- A gp leaf whose lengthscale is stated takes the placeholder for a covariate with too little, and is not
  refused: its kernel width is then the stated lengthscale in the covariate's own units, fixed by the model
  and not read from the rows, so there is a number in force (dec-B465's reason for the refusal does not
  reach it), and the refusal's message names stating it as a remedy. Not taken: refusing it too, which
  makes that remedy false and leaves no gp fit on such a column.
- One row in, or any rows in holding fewer than two distinct values of the response, under a family that reads
  it: the leaf's numbers fall back to every row with the response's, under dec-B309's one warning (the
  maintainer's ruling of 2026-10-23, which replaced this call; [Push 1 fix round](#push-1-fix-round)). The call
  first made, refusing one row in as "one value" while the response read every row, is push 1 as built and its
  fix round's mutant. Under a family that reads no response one row in is one value and refused (dec-B465); no
  row in at all keeps dec-B304's every row.
- A covariate missing on every row in has "no observed value" and is refused with the one-value case. Not
  taken: the placeholder, which is the case dec-B465 refuses.
- A carried record is never refused for its values, an NA scale included: dec-B231 uses a supplied value,
  and the refusal is of a derivation. The no-spread paths of dec-B328 and dec-B329 stay reachable that way
  and through a stated-lengthscale gp, and test-leaf-conversions.R reaches them so. Not taken: refusing a
  record whose scale is NA, which leaves those paths reachable only by `setPredictor` onto one value.
- `setData(updateStandardization = TRUE)` refuses, not keeps: the data are replaced, so no number in force
  belongs to them (dec-B465 lists "a setData that derives"). Not taken: keeping the sampler's numbers for
  that covariate, which is `updateStandardization = FALSE`'s meaning.
- The creation mask reaches the leaf by forest-multiplier push 2's route, read only for these rows; the
  mask is still installed after creation as today, so no other draw moves. Not taken: initializing the leaf
  after the mask's install, which reorders creation for every family.
- The mask in force at a re-derivation is one copy on the sampler with the lifecycle of
  [Engine and bridge](#engine-and-bridge). Not taken: a virtual on the response model read from chain 0, one
  source of the mask but an edit in every family's model.
- `setPredictor(TRUE)` re-derives the replaced columns only; the row setters re-derive every leaf
  covariate. Not taken: a `setPredictor(TRUE)` that re-derives every covariate, which would move columns the
  caller did not name.
- `updateStandardization` on the row setters defaults to FALSE and follows `updateScale`, so 0.9-34's
  positional `setWeights(w, TRUE)` still names `updateState`. Not taken: NULL as on `setData`, which has a
  carried record to defer to and these setters do not.
- On a sampler with no leaf covariate the new argument is accepted and does nothing, as `updateScale` is on
  a probit sampler. Not taken: a refusal, which would make a generic outer loop branch on the leaf model.
- The keep at every leaf re-derivation shares the holder's one warned-once flag with the four response and
  row setters, as slice C's rule 5 counts samplers and not setters; with both arguments TRUE one call raises
  at most one warning. Not taken: a flag of its own.
- Two pushes: the rows and the refusal (engine, draw-changing) apart from the two arguments, which is an
  identity with push 1's `setPredictor` path and reviews against it. Not taken: one push of about 4300
  lines.
- A re-creation of a record-less object reads the data object's weights and no mask (dec-B461), and is
  refused where that leaves too little. Not taken: handing the mirrored mask to re-creation, for objects no
  release wrote.
- The creation warning for no row in is forest-multiplier push 2's: one per creation, joint where the
  response and a leaf both fall back.
- The `setData(s$data)` case under a mask family is named in the table and the help, and the mask riding the
  data object is filed, not built: it is the state of every mask today, and its fix moves `setData` for
  every masked sampler.
- TODO `leaf-covariate-single-value-off-centre` stays out: a column constant to the eye and not to the bit
  holds two distinct values and is not refused; its floor is a change to the spread mark
  ([`standardizedColumnHasSpread`](../../src/bartcore/model.hpp)) and the conversions, not to the rows.

## Assumptions the builder checks first

The recheck of 2026-10-23 read each of the first eight on 030b39f0 and on wt/sfp-3 at 90af5380
([Recheck against bartcore 030b39f0 and state-frame-prior slices 2 and 3 as built](#recheck-against-bartcore-030b39f0-and-state-frame-prior-slices-2-and-3-as-built));
the answers are folded into the sections above. On the tip after state-frame-prior slice 3 lands, before
writing a line, the builder confirms them against the landed code, slice 3's review and fix round being
still to come:

1. Slice 3's path for `updateStandardization = TRUE` on `setPredictor`: as built, the predictor entries
   call `Sampler::rederiveStandardization` after an accepted update, and it re-derives every leaf
   covariate. Step 1.2 narrows it to the replaced columns. (At 9740f9c2 slice 3's fix round has narrowed
   it and pins the narrowing; step 1.2 adds the rows and the keep.) Slice 3's tests as built pin `setPredictor`'s
   argument on one leaf covariate only, so none pins every one; any that does after its fix round is
   listed, not edited, and stops the build for the coordinator. Push 2's identity is written against the
   narrowed path.
2. The record's name and shape: as step 1.6 and the sections above assume (`dataStandardization`,
   `checkStandardizationRecord` inside [`checkDerivedRecord`](../../R/data.R), NULL without leaf
   covariates, an NA scale accepted, a carried heuristic lengthscale taken as the record's and not as the
   model's).
3. A carried standardization enters as a construction argument (`SamplerOptions`), so
   [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp) serves only a view and goes with
   the view's parent constants.
4. Slice 3's `Sampler::setData` takes `DataReplacement` and no mask; its level-code refusals come before
   the store is replaced and its gp refusal is the bridge's, before the call. If its fix round changes
   that, step 1.2 follows it.
5. Forest-multiplier push 2 has landed (dec-A219): nothing is added for it.
6. The baselines: equivalence-32418b6e, bcf-equivalence-20eed7ed and multinomial-equivalence-80b1c8d4, as
   the MANIFEST marks them current. State-frame-prior's base-build comparisons are scratch compares; the one
   fixture of scratch/sfp2/compare.R on a constant covariate, `flat`, is refused here, and none has a row
   out.
7. [`ColumnStore::buildFromParent`](../../src/bartcore/data.hpp) still gathers the parent's constants.
8. The suite's creations on a constant leaf covariate: the recheck found three files beyond the plan's
   (step 1.6) by reading the 37 test files that name `linear(` or `gp(`, not by running them; the builder
   greps again on the landed slice 3 and sizes step 1.6 from the count.
9. Slices 2 and 3 land rebased onto 030b39f0 or later, so the per-iteration setters write through
   [`recordResponseRange`](../../R/data.R) and the records reader without the grids: that read returns the
   standardization, or push 2 adds it there rather than reading the grids.

## Not verified

Nothing is built. The rows-in numbers were emulated on the tip by a `setState` install of the
standardization, which state-frame-prior slice 2 removes; the draws are therefore not slice D's. The xbart
fold's constants were read, not measured: no R reader reaches a fold. The heuristic's subsample over rows in
was read, not run. Which tests create a sampler on a constant leaf covariate was read, at the recheck over
the 37 files that name a linear or gp leaf, not run. The burn-in loop was run for `updateScale` only. The
baselines' verdicts are read from the scenario sources and the MANIFEST, not run; no reference build,
sanitizer, x86 run, `R CMD check` or consumer suite was run. The line estimates are an itemization scaled by
the landed history's ratios. The recheck read state-frame-prior slices 2 and 3 on wt/sfp-3 at 90af5380 and
ran nothing; slice 3 there is built and not reviewed.

## Revision from critique

Item by item, the critique's labels:

- M1, the no-row warning under probit, ordinal and hazard.probit. Taken: the engine reports a per-covariate
  fallback, R raises the class's warning under every family in forest-multiplier push 2's form, joint where
  the response falls back too; tests per family.
- M2, the sampler's copy of the mask. Taken: its lifecycle (set once every chain accepts, normalized,
  replaced at `setData`, put back on a re-creation by the bridge's mask entry); arms after `setActiveRows`
  under gaussian and probit, after `setData` to another row count and on a copy; the two mutants.
- S1, rule 4 without its warning. Taken: dec-B309 cited; every leaf keep warns once per sampler
  (dec-B364). The per-column case is now settled by dec-B465.
- S2, the sampler's own data object carries no mask. Taken: named in the door table, the help and a test;
  filed as a TODO entry.
- S3, a lost record becomes a silent change. Taken: copy of a copy and reload of a copy for linear and gp
  under a mask; assumption 2 checks the heuristic lengthscale's write-back.
- S4, which covariates `setPredictor(TRUE)` re-derives. Taken: the replaced columns only, a second covariate
  held and pinned; assumption 1 checks slice 3's path first.
- S5, test gaps. Taken: a fractional-weight arm (0, 0.001, 3); hazard.probit in the families and tests; an
  xbart fold's "every row" is its training rows, and a fold with too little is now refused (dec-B465).
- NITs. Taken: the `rng:` line names every `xbart` fit with a leaf covariate; null alone means every row;
  `updatePredictor` and the per-observation sessions in the table; the flat C API row names no mask entry;
  the heuristic counts rows in that miss the covariate at 0, as today, pinned; dec-A217 cited; `setModel`'s
  unhonored lengthscale out of scope.
- The question. Put to the maintainer and answered by dec-B464 and dec-B465; options A, B and C and the
  every-row fallback for one value are gone.

Nothing rejected.

## Revision to the rulings

dec-B464 confirms the rows; dec-B465 changes three things here. The one-value case is refused at a
derivation (creation, `setData`, a fold) and kept at a re-derivation, where the plan had every row at
creation and a keep at `setPredictor` only; a covariate constant over every row, created today with the
placeholder, is refused with it, and the tests that created samplers on constant columns are rewritten
(step 1.6). `setWeights` and `setActiveRows` take `updateStandardization`, push 2, with the burn-in idiom in
the help, where the plan had added no argument. The plan becomes two pushes; the estimate grows from about
3100 to about 4300.

## Recheck against bartcore 030b39f0 and state-frame-prior slices 2 and 3 as built

Read 2026-10-23, nothing run beyond reading. Landed since 46e1a644: warm-start-salvage slice D (dec-A220),
forest-multiplier pushes 2 and 3 (dec-A219, dec-A221) and state-frame-prior slice 1 with its fix round
(dec-A222). State-frame-prior slices 2 (ac9f7acf to 1294c756) and 3 (4d15a88d and 90af5380) read on
wt/sfp-3 at 90af5380; slice 3 there is built and not yet reviewed, and has no "as built" record. No ruling
is reopened; the plan's calls stand.

### The plan's cites and reasoning, checked

Holds:
- [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp): every row, shared with the view's
  gather; [`standardizedColumnHasSpread`](../../src/bartcore/model.hpp).
- [`LinearGaussianLeaf`](../../src/bartcore/model.hpp) and [`GPGaussianLeaf`](../../src/bartcore/model.hpp)
  derive in `reinitialize`, which `initialize` calls; slice 2 gives both a supplied standardization.
  [`GPGaussianLeaf::medianPairwiseDistance`](../../src/bartcore/model.hpp): static, a stride over n, 1 where
  degenerate.
- The leaf is initialized in the [`Chain`](../../src/bartcore/chain.hpp) constructor before the response
  model; the response models hold the mask and drop it at a data swap.
- [`SamplerOptions`](../../src/bartcore/chain.hpp) takes borrowed, consumed fields (slice 2's supplied
  standardization is one).
- `Sampler::setActiveRows` fans out and refuses all or none; `Sampler::rederiveScale`; `Sampler::setWeights`.
- `createHolder` and `AmplitudeSpec`: the creation mask is read today only inside the arm for several
  forests, as the plan's "checked once" step assumes.
- `bartcore_create`'s fifth argument and `new(..., active = )`, passed by `dbarts()` through
  [`createEngine`](../../R/dbarts.R), at the first creation alone (dec-A222).
- [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp): three arguments each, so a fourth.
- The holder's `scaleKeptWarned` flag, shared by the four response and row setters through
  [`firstScaleKept`](../../src/R_interface_bartcore.cpp).
- [`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp), called by
  [`refuseRowsRederivation`](../../src/R_interface_bartcore.cpp), the nearer home for push 2's refusal.
- [`installWeightsMask`](../../R/bartcore.R) (slice 3 gives it a `rederive` argument),
  [`bartcoreSamplerSetData`](../../R/bartcore.R) (resolves the 0/1 weights' mask before the bridge call),
  [`warnScaleKept`](../../R/bartcore.R), [`checkUpdateScale`](../../R/bartcore.R), `reapplyActiveRows`
  (passes FALSE to the bridge), [`enforceWeightPolicy`](../../R/spec.R), [`scaleRowsFallBack`](../../R/spec.R)
  (slice 3 reads it at `setData` only where the range is derived).
- `s$data` under a mask family carries no weights and `s$setData(s$data)` drops the mask; with slice 3's
  NULL defaults that call now takes the record, so it derives nothing unless asked.
- [`ColumnStore::buildFromParent`](../../src/bartcore/data.hpp) still gathers the parent's every-row
  constants, served by [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp); the quoted
  "every fold runs the prior a full-data fit would" is `LinearGaussianLeaf::reinitialize`'s comment in
  model.hpp.
- [test-boundary-inputs.R](../../inst/tinytest/test-boundary-inputs.R)'s constant designated column and
  sections 1 and 6 of [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R);
  `checkDataHandleViews` in [leafPriorChecks.R](../../inst/common/leafPriorChecks.R); test_data.cpp's two view
  pins; [linear-exact.R](../../benchmarks/R/linear-exact.R) (four values, 15 trees); the four snapshot files
  fit no linear or gp leaf; [dbarts.h](../../inst/include/dbarts/dbarts.h) untouched.
- Test names: slice 2's "copy holds it"; slice 3's "arguments win", "new row count", "linear converts", "gp
  narrowed" and "fallback where derived"; slice 1's "xbart ignores".
- The design notes' sections: response-scale-rows.md sections 2, 4 and 8 (rule 9 is new, the note having 1
  to 8); state-frame-prior.md's "What response-scale-rows slice D gets". The help files named in 1.8 and
  2.5, dbarts-embedding.Rd's "Whose prior is in force" among them. TODO `response-scale-rows`,
  `mask-on-own-data-object` and `leaf-covariate-single-value-off-centre`.

Moved:
- bcf-equivalence-47838362 (15) to bcf-equivalence-20eed7ed (19), current since dec-A219.
- Forest-multiplier push 2, "built on wt/forest-mult-p2", landed as dec-A219: its warning is raised by
  [`recordBasisScales`](../../R/dbarts.R) from the engine's status, joint through `dbarts()`'s held warning
  and [`basisScaleFallbacks`](../../R/dbarts.R).
- [`Chain::applyNewData`](../../src/bartcore/chain.hpp) "calls `reinitialize` over every row": slice 3
  derives, keeps or takes the standardization per `DataReplacement` (`replacementStandardization`); the
  derive arm still reads every row.
- "`Sampler::setData` and the facade virtual take the mask": slice 3's
  [`Sampler::setData`](../../src/bartcore/sampler.hpp) takes `DataReplacement` as a defaulted last argument;
  the mask becomes a field of it.
- "[`bartcore_setData`](../../src/R_interface_bartcore.cpp) takes the mask as a third argument": slice 3
  gives it five (the data object and three choices as words), so a sixth.
- "the predictor update's derive arm": slice 3 has none; the bridge calls `Sampler::rederiveStandardization`
  after an accepted update of [`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp) or
  [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp).
- "the reader slice 2 adds": `Sampler::derivedRecords` filling
  [`DerivedRecords`](../../src/bartcore/chain.hpp) through
  [`bartcore_derivedRecords`](../../src/R_interface_bartcore.cpp), with slice 1's fix round's read without
  the grids for the per-iteration setters ([`recordResponseRange`](../../R/data.R)).
- `updatePredictor` as a door: no R method has the name; it is `setPredictor` on columns, the bridge's
  `bartcore_updatePredictor`; the per-observation sessions refuse the argument (slice 3).
- "slice 2's 'fresh object neutral' (whose 16 fixtures include weights and a mask)": the 16 fixtures are
  slice 1's row, built as a two-library compare in scratch (18 fixtures); slice 3's row of that name is a
  twin compare with every row in. No tinytest literal of state-frame-prior records a base-build number.

No longer true:
- "Push 2: `Sampler::rederiveStandardization()` ... a facade virtual forwards it" and "push 2 adds one
  virtual": slice 3 added both (every covariate, every row, void). Push 1 changes the virtual's signature
  for the columns, the rows and the keep; push 2 adds no engine method and no facade change.
- "Whatever slice 2 serves through `ColumnStore::suppliedStandardization` for a carried record stays":
  slice 2 serves a carried record through `SamplerOptions` and `suppliedMoments` instead, so that channel
  serves only a view and goes with the view's parent constants.
- "test-gp-leaves.R's zero-weight sampler, which warns now where it did not": its weights keep four rows in
  five in (as at 9e79f4a3), so no fallback is warned; it pins finiteness and stands.
- "the creation warning ... (and a multiplier's column, after that push)": several forests refuse leaf
  covariates, so the multiplier's fallback and the leaf's never meet.
- The suite's constant-covariate creations "beyond test-boundary-inputs.R and sections 1 and 6 of
  test-leaf-conversions.R" were assumed few: test-standardization-record.R (slice 2's `twoFit` on `flat`,
  `flatDonor` into `flatRecipient`), test-warm-start-forms.R's single-value block and
  test-calibration-prior-draws.R's `x4` create samplers on constant leaf covariates, refused here.

### Where state-frame-prior slices 2 and 3 as built differ from what the plan assumed

- The standardization slot's shape is as assumed, NULL where no forest's leaves take a covariate. It is
  read by `dataStandardization`, checked by `checkStandardizationRecord` inside
  [`checkDerivedRecord`](../../R/data.R) (which now takes the model), parsed by the bridge's
  `readStandardizationRecord`, taken by the standard factory only, and dropped for a data handle twice.
  Nothing to change.
- The carried record enters as a construction argument, not as a store channel (assumption 3's first
  branch): the plan's view removal stands and takes the store channel with it.
- `setData` takes per-value choices on `DataReplacement`, the facade virtual's last argument; refusals of
  level codes come before the store is replaced and the gp refusal is the bridge's. The plan now puts the
  mask on that struct (no virtual signature change), tests a derived standardization before
  `ColumnStore::setData` and reads no row under keep or supply.
- `setPredictor(updateStandardization = TRUE)` is a facade virtual, `rederiveStandardization`, called by
  the bridge after an accepted update, re-deriving every leaf covariate over every row: the case assumption
  1 foresaw. The plan narrows it in push 1 and reuses it in push 2. Slice 3's tests use one leaf covariate,
  so none pins every one. Slice 3's help says it derives "over the values the update leaves"; 1.4 and 1.8
  rewrite that. Narrowing in slice 3's own fix round, before it lands, would save that rewrite and the
  every-covariate behavior's brief life on bartcore; that is the coordinator's call, not this plan's.
- `setWeights` and `setActiveRows` are untouched by slices 2 and 3; on wt/sfp-3 they still write the record
  through `recordDerived` after `updateScale`, which slice 1's fix round replaced on bartcore by the read
  without the grids. Push 2 writes the standardization by that read.
- The warned-once flag is unchanged; the predictor entries return only the accepted logical, so the plan
  adds an attribute for "warn now" and the kept covariates. The keep's text named the response only and
  the plan gave no leaf words; [The refusal's message](#the-refusals-message) now gives them.
- The gp refusals: slice 3 narrows `setData`'s to calls that change the standardization or the range and
  refuses `setPredictor`'s TRUE before anything moves (`standardizationUpdateFromExpression`); push 2 reuses
  its words.
- Slice 3's keep at `setData` (FALSE) has a number in force, so the `setData` refusal now names it as a
  remedy, as slice 3's gp refusal does.
- Base-build literals: slices 1 and 2 replaced them by two-library compares in scratch, draws recorded on
  arm64 not holding on the other CI platforms; this plan's neutral arm does the same for its draws.

### Budget

Raw itemization from ~2530 to ~2610 (push 1 from ~2030 to ~2110, push 2 ~500 as before), forecast from
~4330 to ~4490 (~3500 to ~3650, ~830 to ~840), package code from ~665 to ~725 (~525 to ~590, ~140 to ~135),
stops at push 1 ~5450 lines or ~890 of package code and push 2 ~1250 or ~200
([Budget and stops](#budget-and-stops)). The measure is now named: lines added. Push 1 grows by the keep's
words and the joint creation warning, three more test files, the narrowing pins and a larger view removal;
push 2's engine and facade work is slice 3's already. The ratios are kept, at 1.73 overall with a fix round
counted: the slices landed and built since ran 0.9 to 2.4 times their plans in lines added, the coordinator's
reading 1.5 to 2.

### Open questions

None: nothing the recheck found changes what a user sees outside what dec-B464, dec-B465 and the plan's
calls settle. A consequence of a standing call, for the reviewer: under the one warned-once flag (the call
read from dec-B465's "once per sampler"), a sampler that has warned of one kind of keep says nothing of a
later keep of the other kind. (The second consequence named here at the recheck, narrowing slice 3's
every-covariate re-derivation, is slice 3's own since its fix round.)

## Recheck against state-frame-prior slice 3's fix round

Read 2026-10-23 at 9740f9c2 (wt/sfp-3, slice 3's fix round, not yet landed), before building push 1, with a
run of the suite on a build of 9740f9c2 that logs every derivation on a column constant over every row
and every creation whose rows of positive weight hold one value of a leaf covariate (scratch/rsrd1/probe-src,
scratch/rsrd1/probe-logs). The sections above are amended in place; this records what moved.

- The narrowing is slice 3's. Its fix round has `setPredictor(updateStandardization = TRUE)` derive the
  covariates among the replaced columns alone, over every row, the others held:
  [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp) marks its columns and calls
  `rederiveStandardization` with them, [`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp)
  replaces every column and passes null, and `Chain::rederiveStandardization` hands each forest's leaf the
  numbers in force with not-a-number at the marked covariates (`Chain::heldStandardization`). The facade
  virtual already takes the marks, so push 1 changes its return only, to the covariates kept, and the row
  set is the sampler's own read, not an argument. The keep is that same supply with a too-little covariate
  left unmarked. Slice 3's tests pin the narrowing on linear and gp leaves (test-setdata-records.R, "TRUE
  derives the replaced columns' alone") and tests/cpp pins a mark on another column deriving nothing, so
  push 1's arm for a second covariate held runs with a row out, its tests/cpp pin likewise, and the mutant
  of every covariate re-derived is replaced by one of a too-little covariate re-derived. The help and the
  docstring of `setPredictor` already say "the covariates among the columns the update replaces"; push 1
  adds the rows and the keep. Assumption 1 holds in its narrowed form: no test pins every covariate.
- `bartcore_setData`'s copies. The bridge resizes and overwrites the holder's response, weights, offset and
  test-offset buffers, which the sampler borrows, before it calls `Sampler::setData`, so the plan's refusal
  before the store is replaced would come after the sampler's borrowed buffers had moved. Push 1 has the
  bridge copy into buffers of its own and move them into the holder once the engine accepts (a move keeps
  each buffer, as `createHolder` relies on), and catch the typed error there; bridge ~75 to ~85, one
  mutant added.
- The suite's creations on constant leaf covariates: by the logging run, six files derive a
  standardization on a column constant over every row, the five step 1.6 names and test-setdata-records.R's
  `flat` gp sampler (slice 3's fix round, which pins an NA scale compared as the leaf divides); it reaches
  its covariate without spread through a stated lengthscale instead. No creation in the suite has its rows
  of positive weight at one value of a leaf covariate (the run logged weights, not masks). Assumption 8 is
  answered by that run; the stop on a seventh file stands. The run also shows test-mutate-then-serialize.R
  failing at 9740f9c2 itself: its last section's `setState` of a state without the leaf blocks is refused
  as holding gp draws under another standardization. It is the base's, not this slice's; push 1's suite
  gate reports it beside its own count. test-setdata-records.R's unforced `setPredictor(rep(0.5, n), "b", ...,
  updateStandardization = TRUE)` is pinned there as rolled back (FALSE), so it derives nothing and keeps
  nothing, and stands unedited.
- Assumptions 2 to 7 and 9 hold as the recheck of 030b39f0 read them: the record is `dataStandardization`
  and `checkStandardizationRecord` inside [`checkDerivedRecord`](../../R/data.R), its keys since slice 2's
  third fix round a name only where it names exactly one column; a carried standardization enters through
  `SamplerOptions`; `Sampler::setData` takes `DataReplacement` with no mask, its level-code refusals before
  the store and the gp refusal in the bridge; [`ColumnStore::buildFromParent`](../../src/bartcore/data.hpp)
  still gathers the parent's constants; the MANIFEST's current files are the three named; the read without
  the grids returns the standardization.
- Slice 3's fix round adds `setData`'s refusal of an ordered factor's other cut points and the count
  family's equal range refusal; neither touches a leaf covariate or a door of this slice.
- Budget: push 1 raw ~2110 to ~2075 (engine ~185 to ~170, the narrowing done; bridge ~75 to ~85, the
  copies; tests/cpp ~330 to ~315 and tinytest ~1070 to ~1060, the narrowing's pins slice 3's, one more
  constant-covariate edit; help ~110 to ~105), forecast ~3650 to ~3590, package code ~590 to ~585 forecast;
  the stops at 1.5 times, ~5390 lines and ~880 of package code. Push 2 unchanged.

No question changes what a user sees outside dec-B464, dec-B465 and the plan's calls; nothing stops the
build.

## Push 1 as built

Built on 9740f9c2 (wt/sfp-3, state-frame-prior slice 3 in its fix round, not landed) with this plan's two
rechecks on top, the steps rechecked against it first
([Recheck against state-frame-prior slice 3's fix round](#recheck-against-state-frame-prior-slice-3s-fix-round),
e5cd37e8). Commits: a69c8bca (the push), 6440e06c (the help, NEWS and notes) and this record.

What it does. A linear or gp leaf derives its covariates' center and scale, and a gp leaf its default
lengthscale, over the rows in the likelihood, rows of positive weight that no mask leaves out, by membership
([`Sampler::leafRowsFrom`](../../src/bartcore/sampler.hpp)): null where every row is in, so those runs take
today's loops bit for bit, and every row where none is in, the leaf marking the covariates so read
(`DerivedRecords::leafFromAllRows`). [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp) takes the
row set and answers too little, fewer than two distinct observed values;
[`GPGaussianLeaf::medianPairwiseDistance`](../../src/bartcore/model.hpp) strides over the rows read. At creation
the rows are the weights' and the creation mask's, a new field of [`SamplerOptions`](../../src/bartcore/chain.hpp),
`activeRows`, which `createHolder` now reads once for either arm; at `setData` the new data's weights and a new
field of [`DataReplacement`](../../src/bartcore/data.hpp), set from `bartcore_setData`'s sixth argument; in an
`xbart` fold the fold's weights over its gathered rows, `ColumnStore::buildFromParent` gathering the values alone
(`suppliedStandardization`, `gatheredMeans` and `gatheredSds` are gone). A derivation with too little throws
[`LeafStandardizationError`](../../src/bartcore/model.hpp) from [`deriveLeafMoments`](../../src/bartcore/model.hpp),
naming the covariate, its column and the counts, unless a gp leaf's model states its lengthscales, which takes
the placeholder; `Sampler::setData` tests every covariate it will derive before anything moves
([`Sampler::refuseLeafDerivation`](../../src/bartcore/sampler.hpp)), and the bridge words the error at creation,
at `setData` with the FALSE remedy and in a fold (`leafStandardizationMessage`). Slice 3's narrowed
[`Sampler::rederiveStandardization`](../../src/bartcore/sampler.hpp) reads the rows the weights and the sampler's
copy of the mask leave in ([`Sampler::keepActiveRows`](../../src/bartcore/sampler.hpp): the mask every chain
holds, cleared by a replacement), holds a replaced covariate with too little by leaving it out of
`Chain::heldStandardization`'s marks, and returns the covariates kept; the facade virtual returns them and the
predictor entries carry them, with "warn now" from the holder's one flag, on their answer
([`predictorUpdateAnswer`](../../src/R_interface_bartcore.cpp), `firstKeep`), which R strips and warns of last in
the leaf's words ([`warnScaleKept`](../../R/bartcore.R)). R raises the creation warning where a leaf read every
row, one text with the response's under gaussian, Student-t and nbinom, `dbarts()` dropping its held spec warning
then, and `setData` the same where it derived ([`warnRowsFallBack`](../../R/bartcore.R),
[`leafFromAllRows`](../../R/bartcore.R)). `bartcore_setData` copies the new response, weights, offset and test
offset into buffers of its own and moves them into the holder once the engine has taken the replacement.

Departures from the steps and the test table:

- The sampler's copy of the mask is the mask every chain holds: a replacement clears it, and the mask on
  `DataReplacement` is read for the derivation alone, the host installing it after as before. The plan had the
  replacement's mask replace the copy, which between the swap and the host's install would claim a mask no chain
  holds.
- `Sampler::setData` tests a supplied record's unmatched covariates (a not-a-number center) as well as a derived
  standardization's, so a partial record cannot refuse after the swap.
- hazard.probit with no subject in warns twice: the hazard periods' text, counted in subjects and raised before
  the sampler exists, and the leaf's, counted in person-period rows. The test table had one warning.
- test-response-scale-rows.R's gp state arm, not named in step 1.6, leaves the standardization blocks out of its
  forest compare: donor and recipient are created under other rows in, so each holds its own (a gp leaf with rows
  out, which the stops allow to move). The direct `.Call` sites of `bartcore_setData` in test-data-handle.R,
  test-multinomial-counts-mutation.R and test-sampler-boundary-refusals.R gain the sixth argument.
- The held arm pins the record and the engine's reader before and after each call, not 5 sweeps against a twin;
  the record arm compares a route's next sweeps with the source's to 1e-10, the live fits being rebuilt from the
  trees, and its record and reader bit for bit.
- The xbart arms of test-linear-leaves.R and test-gp-leaves.R and leafPriorChecks.R's `checkDataHandleViews` pass
  unedited; test-leaf-conversions.R's sections 1 and 6 reach a constant through a carried record (`onRecord`), its
  seeded `bart()` call from a data object.
- linear-exact.R's arm puts its 200 rows out at the first value with the mean response; its draws are the first
  arm's, bit for bit.
- Step 1.9's TODO entry, ledger entry, indexes and landing note are the landing's; the design notes are written.

Calls made in building, each with its alternative:

- The refusal where no row is in and a covariate is constant over every row reads "holds one value over all n
  observations, which are read because none is in the likelihood". Alternative: the plan's "(0 of n observations
  are in it)", which reads as though no row was read.
- An unnamed column is "the leaf covariate in column j" in the refusal and the keep. Alternative: `x[, j]`.
- Several kept covariates are named "each of the leaf covariates x1 and x2 ... their centers and scales".
  Alternative: one warning per covariate, against dec-B465's one per sampler.
- A stated-lengthscale gp's placeholder records its scale as NA only where no training row reads off the value,
  and 1 where rows out do, the spread mark reading the values it holds. Alternative: a placeholder mark of its own.
- The no-row flags are a fifth entry of `bartcore_derivedRecords`'s list present only where one is set, so the
  reader equals the record a data object carries, which slices 2 and 3's tests compare with `identical()`.
  Alternatives: always present, which breaks those compares, or a new entry point.
- `bartcoreSamplerSetPredictor` returns `list(value, visible, kept)` to the method. Alternative: an attribute on
  the value, which a forced update's `NULL` cannot carry.
- The bridge answer is a vector of its own: `Rf_ScalarLogical` hands out R's shared TRUE and FALSE, and an
  attribute set on one leaked onto every later TRUE (found by the new file, before commit).
- Every mutant installed with `--preclean`, R-only ones included, the export holding the previous build's objects.

Tests. tinytest test-leaf-standardization-rows.R (225 expectations): at creation for gaussian, Student-t,
probit, ordinal, nbinom, hazard.probit, logistic and aft, linear and gp, the record and the reader against
literals over the rows in (every row for logistic and aft); fractional weights; the subset twin; every row in for
each family; no row in, each family's warning and the joint text through `dbarts()` and `bart()`, apart by hand;
the refusals at `dbarts()`, `bart()` and `new()` with `active` for rows in at one value, a constant column, none
observed and one row in, linear and gp, and the stated-lengthscale gp created; rows out at one value created;
missing values among the rows in; `new()` with and without `active`; nine calls that hold; `setPredictor`'s
re-derivation after `setWeights`, after `setActiveRows` under gaussian and probit, after a `setData` of more rows,
on a copy, a second covariate held, the keep and its one warning, the flag shared with `setOffset`'s keep, several
kept, FALSE; `setData` by weight and by nbinom mask, a carried record, TRUE, FALSE, the refusal at NULL and TRUE
leaving data, state, response and draws as an untouched twin's, served at FALSE, a probit sampler's own data
object; the record through a copy, a copy of a copy, a reload, a reload of a copy and `new()` from `s$data`; gp
kept-draw refusals; an `xbart` fold's record and its refusal. tests/cpp `testLeafStandardizationRows`: the moments
over a row set, the all-ones set's bits, too little; at creation by weight and by mask under gaussian, probit and
nbinom, linear and gp, the literals with the heuristic over 960 rows in with a missing row among those its stride
takes, and 20 sweeps identical to a record's twin; no row in; the refusal and the stated-lengthscale placeholder;
the re-derivation over the rows in, the mask in force, a refused mask, the keep reported, no row in; a replacement
deriving over its own mask's rows, dropping the mask, and a refused replacement untouched.

Existing tests changed: test-boundary-inputs.R (a constant designated column refused under linear and gp, the
stated-lengthscale gp created); test-leaf-conversions.R sections 1 and 6, test-standardization-record.R's `twoFit`
and `flatDonor`/`flatRecipient`, test-warm-start-forms.R's single-value block, test-calibration-prior-draws.R's
`x4` and test-setdata-records.R's `flat` (a carried record or a stated lengthscale); test-response-scale-rows.R's gp
state arm; the three direct `.Call` sites; tests/cpp's view pins in test_data.cpp and test_model.cpp, which now pin
a view's leaf deriving over the view's rows and a constant view refused, testNoSpreadStandardization,
testLinearLeafOneSidedConversion and testNoSpreadLiveConversion (constants through a supplied record), and the
facade spy.

Mutants, each alone on a clean export, installed with `--preclean`; failures in test-leaf-standardization-rows.R
unless named, "stops" where the file stops at a refusal or an error: the row set ignored, 76, tests/cpp 16, and
linear-exact.R's arm (largest |z| 39.9); membership read off the weight's size (0.001 out), 1; the creation mask
not handed, 21; the weights ignored, the mask read, 51; the copy keeping the first mask, 1 (it survived a first run,
the probit arm's sampler having had no creation mask; the arm now creates it under one); the copy not cleared at
`setData`, 1, tests/cpp 1; the copy set before the chains accept, tests/cpp 2; the placeholder for the refusal, 34,
test-boundary-inputs.R 2; the stated-lengthscale gp refused, stops in the new file, test-boundary-inputs.R and
test-setdata-records.R; refused where only rows out are constant, stops; `setData` refusing after the swap, 4,
tests/cpp 1; one row in read as none, 6; no fallback where none is in, stops, tests/cpp aborts; the probit no-row
warning dropped, 4; two warnings for a joint fallback, 4; the re-derivation over every row, 7; every keep warned, 2;
the response's keeps not spending the flag, 1; the response's text for a leaf keep, 3; FALSE refused, stops; the
heuristic over every row, 10; its stride over every row's index, 1, tests/cpp 3 (tests/cpp's first two-in-three
pattern took the same rows under either stride; now four in five); the heuristic dropping a missing row in, 1,
tests/cpp 3 (the missing row moved into the subsample); a view supplying its parent's constants, 2; a fold with too
little run, 1; `setData` handed no mask, 1; a re-creation deriving, 1 then stops, and test-standardization-record.R
stops; an NA scale refused, stops in the new file and test-leaf-conversions.R; missing rows counted among the rows
in, 1; `bartcore_setData` adopting its copies before the engine call, 2 (it survived a first run: the response
reads the raw values again only at a re-derivation, which the refusal arm now makes). Alive as expected: the row
set built where every row is in (the same sums in the same order), and an all-ones mask kept as a copy, which
`leafRowsFrom` reads as every row by construction (not run).

Size, against e5cd37e8 before this record: engine 376 added and 120 removed, bridge 251 and 67, R 184 and 23,
tests/cpp 472 and 46, tinytest 948 and 50, the exact arm 108 and 60, help and NEWS 25 and 13, notes 56 and 19;
2420 added and 398 removed in all, against ~3590 forecast at landing and ~2075 raw (stop ~5390); package code 811
against ~585 forecast (stop ~880), its engine and bridge over their itemization (the typed error and its wording,
the copy's lifecycle, the buffer staging, the keep's answer), so a fix round has about 70 lines before that stop.

Gates, on the push's libraries, from a `--preclean` install of 6440e06c (this record touches the plan alone; the
reference, exact, check and R-loaded sanitizer runs were made on the push's earlier spellings, whose package
code is 6440e06c's, the later changes touching only test-leaf-standardization-rows.R and test_state.cpp): tests/cpp plain and under ASAN and
UBSAN, 411 ok, no failure and no diagnostic; the suite in one process file by file as `test_package` runs them,
271 files, 28982 expectations, none failed, test-mutate-then-serialize.R stopping on an error that was 9740f9c2's own
(restated after the rebase onto f58d1482: a7189447, in slice 3's landing, fixed it, and the file passes on the
rebased push; [Push 1 fix round](#push-1-fix-round)); lint (0 lints), air, rc-codoc (48 methods), win-drift,
doc-freshness and verify-anchors clean; on a `--preclean` reference install the four snapshot files (14, 3, 7 and 3
expectations, none failed) and the three `--bitwise` compares, 55 of 55 identical against equivalence-32418b6e
(`--strict-coverage`, 55 compared, 0 skipped), 19 of 19 against bcf-equivalence-20eed7ed and 11 of 11 against
multinomial-equivalence-80b1c8d4, no "max |z|" line, as expected; the 43 commands of exact-gates.yaml at `quick`
on the shipped build, all passing, linear-exact.R's new arm within |z| 0.8, and the two cross-host compares on the
reference passing; the R-loaded run over 13 touched and nearby files with the package built under ASAN and UBSAN
through `CXX20FLAGS` (the bridge's compile line carrying `-fsanitize=address,undefined`), 2332 expectations, none
failed, no diagnostic, and the final test-leaf-standardization-rows.R, 225, none failed; `R CMD check --as-cran
--no-manual` from a clean tarball, the Date field's NOTE, its tests then stopping on the same file (restated: the
Date NOTE alone after the rebase); scratch/rsrd1/compare.R against 9740f9c2, 12 fixtures with every row in and every covariate spread
(linear, gp and two covariates under gaussian, probit, Student-t, nbinom, logistic, aft, ordinal, hazard.probit and
constant leaves, with `setData` of a fresh object and `setPredictor` deriving), 96 of 96 parts identical, and the 6
that move all moved (24 parts: rows out by weight and by mask, centers 7.18 and scale 9.50 over every row against
0.491 and 0.323 over the rows in, a gp lengthscale 0.080 against 1.03, `xbart` with linear and gp leaves, a constant
column refused); stan4bart's suite (its bartcore branch as it then stood, built against this push), 1545
expectations, none failed (restated: its bartcore branch is now at 3716b51, 1555 expectations, none failed). Not
run: the x86-64 leg. These gates ran before the rebase; the fix round's are in its record.

Left open. Package code stands at 811 of the 880 stop. The x86-64 leg waits for landing. (Restated: the
test-mutate-then-serialize.R failure this paragraph named is gone with the rebase.)

### Push 1 fix round

The orchestrator rebased the push onto bartcore f58d1482, state-frame-prior slice 3's landing among it (its fix
round's a7189447 fixes the test-mutate-then-serialize.R failure "Push 1 as built" named): the push is ab2d350a,
e5cd37e8, a69c8bca (the code), 6440e06c (the help and notes) and 8a887960 (the record). Its review, FIX FIRST: eight
mutants surviving (R12, R13, R26, R07, R08, R05b, R21, R22) and R14 killed by tests/cpp alone; the no-row warning
garbling an unnamed column; the help naming "the hazard families" among the mask families; hazard.probit with no
subject in warning twice; nits (a fixed refusal buffer, a part-way failure leaving pointers into freed staging
buffers, two spellings); stale hashes. Then the maintainer's ruling of 2026-10-23 below. Commits: eff97e92 (the
review's fixes and pins), 4b618c6c (the ruling) and this record, which also restates the hashes and the claims
above.

The maintainer's ruling of 2026-10-23, on dec-A216 and this plan's planning call (one row in refused as "one
value"): where the response falls back to every row because fewer than two distinct values of it are in the
likelihood at creation (dec-B309's case, one row in say), every quantity the prior reads from the rows in falls
back with it, a linear or gp leaf covariate's center, scale and default lengthscale among them, under dec-B309's
one warning; dec-B465's refusal stays for its own case, the response spread over the rows in and a leaf covariate
at one value there; re-derivations keep dec-B364 and dec-B465's keep. Built in
[`Sampler::withResponseFallback`](../../src/bartcore/sampler.hpp): under a family that reads the response
(gaussian with either residual distribution, nbinom), with a row out and the response less its offset (the counts
under nbinom) at one value over the rows in, the leaf reads every row and marks its covariates, at creation, a
`setData` that derives the range, or an `xbart` fold (whose sampler is built through the same constructor tail);
the creation text counts the rows in. The planning call, the mutant table,
the help and rule 9 are restated. Read as the ruling's doors: a family that reads no response (probit, ordinal,
hazard.probit) has no fallback to follow, and still refuses one row in as one value. A multiplier forest takes no
linear leaf ("a treatment forest does not support a linear leaf prior"), refused before any warning, so a sampler
with both cannot warn of both; that refusal is pinned.

What changed besides:
- Pins in test-leaf-standardization-rows.R, each on the reviewer's probe: a partial carried record refused before
  the swap, the reader and 3 sweeps an untouched twin's (R12); a stated-lengthscale gp at a deriving `setData`, at
  NULL and TRUE, the placeholder taken (R13); a refused `setData` leaving the mask in force for a later derivation
  (R26); after `setWeights`, a kept covariate leaving one not replaced held (R07); a covariate at one value that the
  call does not replace neither reported nor kept (R08); a gp leaf's no-row warning under gaussian and probit
  (R05b); `setData`'s no-row warning under gaussian and probit (R21); the unnamed no-row and keep texts (R22); a
  refused covariate whose position is not its column, the reader unmoved (R14); a 600-character column name named
  whole, the remedy after it.
- The no-row text names an unnamed column "the leaf covariate in column j" ([`leafCentersAndScales`](../../R/bartcore.R),
  [`warnRowsFallBack`](../../R/bartcore.R) composed once).
- hazard.probit with no subject in: `dbarts()` holds the periods' warning (`holdPeriods`) and hands it to the
  creation, which folds the leaves' numbers into its text, one warning; rule 9 says so. By hand, `dbartsSpec()` then
  `new()`, the two stay apart, as the response's do.
- The refusal is a `std::string` the closure owns at the three sites, so a long name is not cut and the jump leaks
  nothing.
- A `Sampler::setData` failure other than the refusal (an allocation) can leave chains part-swapped onto the
  closure's buffers; the bridge now releases the sampler (`holder.sampler.reset()`) and raises the failure, so
  nothing points into freed buffers, and the R object re-creates the sampler from its stored state at its next use,
  as after a reload. No test reaches it (it needs an allocation failure part-way).
- The help says hazard.probit; data.hpp's comment and test-leaf-conversions.R's say centered; the refusal over every
  row says "read because too few are in the likelihood", true of no row in and of the response's fallback.

Mutants, each alone on a clean export of 4b618c6c, `--preclean`; failures in test-leaf-standardization-rows.R: the
reviewer's R12 2, R13 stops (its setData refused), R26 1, R07 2, R08 1, R05b 4, R21 3, R22 2, R14 1; the periods not
folded 1; the unnamed label unfixed 1; the refusal cut at 511 bytes 1; one row in refused where the response falls
back (the call before the ruling) stops at its creation; the setData fallback blind to the range choice 1.

Size, against 8a887960 before this record: engine 27 added and 3 removed, bridge 36 and 28, R 83 and 55, tinytest
238 and 18, help 2 and 2, rule 9 4 and 2; 390 added and 108 removed. The push's package code stands at 877 lines
added against f58d1482, under the 880 stop, the ruling's engine lines paid for by composing the fallback text
once.

Gates, from a `--preclean` install of 4b618c6c (this record touches the plan alone): tests/cpp plain and under ASAN
and UBSAN, 411 ok, no failure and no diagnostic; `tinytest::test_package` in one process, 266 files with
expectations, 29662, none failed, its 12 warnings those of the reviewer's run at 8a887960; lint (0 lints), air,
rc-codoc (48 methods), win-drift, doc-freshness and verify-anchors clean, NEWS parsing with 227 entries; on a
`--preclean` reference install the four snapshot files (14, 3, 7 and 3, none failed) and the bitwise compares, 55 of
55 identical against equivalence-32418b6e (`--strict-coverage`, 55 compared, 0 skipped), 19 of 19 against
bcf-equivalence-20eed7ed and 11 of 11 against multinomial-equivalence-80b1c8d4, no "max |z|" line, neutral as
expected (no recorded scenario has one row in with a leaf covariate), and both cross-host compares passing; the 43
commands of exact-gates.yaml at `quick`, all passing, linear-exact.R's arm within |z| 0.8; the R-loaded run with
the package built under ASAN and UBSAN through `CXX20FLAGS` (the bridge's compile line carrying
`-fsanitize=address,undefined`) over 17 files, the 13 of "Push 1 as built", test-mutate-then-serialize.R and the
three hazard files, 2692 expectations, none failed, no diagnostic; `R CMD check --as-cran --no-manual` from a clean
tarball, the Date field's NOTE alone; stan4bart's suite, its bartcore branch at 3716b51 built against this push, 25
files and 1555 expectations, none failed. Not run: the x86-64 leg.

Left open: the ruling's reading for families that read no response, above; package code is 3 lines from the stop.
