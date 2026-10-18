# response-scale-rows, slice D: a leaf's covariate center, scale and lengthscale over the rows in

Status: PLANNED 2026-10-17 against 9e79f4a3 (bartcore; state-frame-prior designed and planned, not built;
forest-multiplier push 2 built on wt/forest-mult-p2, not landed); revised from its blind critique
(scratch/rsrd/crit/crit-report.md, sound after changes; [Revision from critique](#revision-from-critique))
and then to the maintainer's rulings dec-B464 and dec-B465 ([Rulings](#rulings)). No call open.
Implements dec-B304 as dec-B464 confirms it and dec-B465 extends it. Rests on dec-B231, dec-B233,
dec-B237, dec-B296, dec-B302, dec-B303, dec-B309, dec-B329, dec-B330, dec-B364, dec-B457 to dec-B461 and
dec-A217. Builds on state-frame-prior slices 1 to 3 ([state-frame-prior.md](state-frame-prior.md)); the
standing statement is the design note [response-scale-rows.md](../design/response-scale-rows.md), which
this slice extends.

agent: two pushes, serial. Per push an opus implementer for the engine, the facade, the bridge, tests/cpp
and (push 1) the exact arm; sonnet for the R lines, the tinytest files, the help and the records; one opus
reviewer told to refute, who runs that push's mutants ([Mutants](#mutants)) and reads the diff twice: once
for bits that must not move where every row is in, once for which rows each door reads, under which mask,
and what it refuses or keeps. Tier: changes draws ([Process by risk](README.md#process-by-risk)) for both.
rng, against equivalence-32418b6e, bcf-equivalence-47838362 (or push 2's re-record) and
multinomial-equivalence-80b1c8d4 ([RNG class and expected verdicts](#rng-class-and-expected-verdicts)):
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
No recorded scenario reaches a changed sequence in either push: expected 55, 15 and 11 identical, no
re-record.
window: directly after state-frame-prior slice 3, on top of its records; push 2 after push 1. Serial with any
engine slice in the leaves' derivation, [`Chain::applyNewData`](../../src/bartcore/chain.hpp), the
setData, predictor, weight and mask paths of [`Sampler`](../../src/bartcore/sampler.hpp) or
`ColumnStore::buildFromParent`. Before 1.0-0, and before the merge to main if it can be (TODO
`response-scale-rows`).
budget: raw itemization ~2530 lines (push 1 ~2030, push 2 ~500); forecast ~4300 from the landed history
(~3500 and ~830), package code ~665 (~525 and ~140); stops at 1.5 times the forecast
([Budget and stops](#budget-and-stops)). No ABI event.

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

Today (read).
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
  created with no warning (critic, q1.R). Push 2 of forest-multiplier adds the warning under probit and
  logistic for a multiplier's column, raised by the sampler's `initialize` from an engine status, and names
  the column in the gaussian text.
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

The nearest precedent, the numeric multiplier's scale (forest-multiplier push 2, built on wt/forest-mult-p2,
its "Push 2 as built"): the scale is the column's sd over the rows the response's range was read from, under
probit and logistic the rows the creation mask keeps, every row where none is in, with the fallback reported
by the engine and warned by R under every family; a column constant over those rows refused where a default
would need its spread (dec-B303); taken once, recorded on the data object as `basis.scale`, handed to every
later construction; the creation mask reaches the engine through
`new("dbartsSampler", control, model, data, active = )` and a fifth argument of `bartcore_create`. Slice D
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

Assumed as that plan states it; each is checked first
([Assumptions the builder checks first](#assumptions-the-builder-checks-first)).

- The record: slice 2's `standardization` slot on `dbartsData`, per forest `center` and `scale` (NA for a
  column without spread) named by covariate, and `lengthscale` only where the model supplies none, written
  by the sampler from the engine's reader. This slice adds no slot and changes no shape; the numbers its
  derive path writes there are the rows-in ones, and push 2 writes the record after each re-derivation. The
  record is why this slice could not land earlier: a copy, a reload and a re-creation take it and derive
  nothing, so numbers read under a mask survive a re-creation that has no mask at construction. Any record
  lost on a re-creation path becomes a silent change of numbers, or now a refusal, so the tests take two
  levels of copy and reload.
- The doors: slice 2's construction, which takes a supplied standardization or derives one ("absent
  derives"), through both creating factories; this slice edits the derive arm only, and a carried record is
  never refused for its values. Slice 3's `Sampler::setData`, which takes derive, keep or a supplied value
  per record; this slice edits the derive arm, hands it the mask and refuses there. Slice 3's
  `updateStandardization` on `setPredictor` and the predictor entries; this slice edits its derivation,
  pins that it re-derives only the replaced columns and gives it the keep. Slice 1's data handle, which
  ignores every record; this slice changes what a fold then derives, its own rows, and refuses there.
- The tests: slice 2's "copy holds it" and "fresh object neutral" (whose 16 fixtures include weights and a
  mask), slice 3's "arguments win", "new row count", "linear converts", "gp narrowed" and "fallback where
  derived", and slice 1's "xbart ignores", each named below where it is extended or where a literal of it
  legitimately moves.

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
| `setPredictor(x, column, ..., updateStandardization = FALSE)`; `updatePredictor`; the per-observation sessions | held, regathered under the numbers in force | nothing |
| `setPredictor(..., updateStandardization = TRUE)` | the replaced leaf covariates only, re-derived over their new values at the rows in at that call (the weights and the mask in force); every other leaf covariate held; too little keeps, warned once per sampler; gp with kept draws refused | posterior where a row is out |
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

Raised by R from the engine's error, which names the forest and the covariate's column (push 2's
`BasisScaleError` route); the bridge supplies the column's name and the counts. Three texts:

- linear, or gp with no stated lengthscale, at creation or `setData`:
  `the leaf covariate "x1" holds one value over the rows in the likelihood (100 of 200 observations are in
  it), so its center and scale cannot be taken from them; drop it from linear(), or change the rows in
  (the weights, the 0/1 weights' mask, or subset)` and, under gp, `...; state its lengthscale in gp(), drop
  it from gp(), or change the rows in (...)`. "holds no observed value" where every row in misses it.
- an `xbart` fold: `in a fold of xbart, the leaf covariate "x1" holds one value over the fold's training
  rows in the likelihood, so its center and scale cannot be taken from them; drop it from linear(), or
  change the rows (the weights, subset or n.test)`, with the gp clause under gp.
- every row in, the covariate constant: the first text, its count "200 of 200".

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
- bcf-equivalence-47838362, 15 of 15 identical; if forest-multiplier push 2 has landed first, its
  re-record (on the branch bcf-equivalence-85bc9ffb, 19 scenarios), 19 of 19. Several forests refuse leaf
  covariates.
- multinomial-equivalence-80b1c8d4, 11 of 11 identical. A multinomial sampler refuses them (ran).
- The four `test-reproducibility-*.R` files unchanged; none fits a linear or gp leaf.
- exact-gates.yaml at quick: every recorded command as before, plus the new arm of
  [linear-exact.R](../../benchmarks/R/linear-exact.R) (push 1) passing. The SBC leaf arms use positive
  weights and are not run.
Any other count is a stop. No baseline is re-recorded. No hot path: the derivation runs at creation,
`setData`, a deriving setter and a fold's build, one pass per leaf covariate; no bench compare.

## Engine and bridge

The engine decides which rows are in, from the weights and the mask it is handed, in one place, reports
where it read every row and which covariates it kept, refuses by a typed error naming the column, and never
warns or names an R type. The bridge hands it the mask, words the refusal and carries the status; R raises
the warnings.

Engine.
- [`standardizationMomentsForColumn`](../../src/bartcore/data.hpp) takes an optional row set (null for
  every row, today's loop unchanged) and answers whether the rows read held fewer than two distinct
  observed values. `GPGaussianLeaf::medianPairwiseDistance` takes the same rows, its stride counted over
  the rows read.
- The leaves' `initialize` and `reinitialize` take the row set and, for a derivation, throw an error naming
  the covariate where it has too little, unless the gp leaf's lengthscales are supplied; a re-derivation
  takes the columns it replaces (all of them for the row setters), keeps a covariate with too little and
  reports which it kept. Each leaf reports, per covariate, whether its numbers came from every row for no
  row in, beside forest-multiplier's `basisScaleFromAllRows`, through the reader slice 2 adds.
- A derivation at `setData` is tested before the swap, so a refusal leaves the sampler as it was.
- [`Sampler`](../../src/bartcore/sampler.hpp) forms the row set once, from the weights and a mask, null
  where no row is out or where none is in, and hands the same bytes to every chain, as slice C decides the
  keep once from chain 0. The mask comes from:
  - at creation, a new borrowed, consumed field of [`SamplerOptions`](../../src/bartcore/chain.hpp);
  - at `setData`, its new argument, which replaces the sampler's copy (null clears it) before any chain
    derives, as the swap drops each response model's;
  - at a re-derivation, the sampler's copy, set by `Sampler::setActiveRows` only once every chain has
    accepted the mask, normalized as the chains normalize it (null and all ones are no mask), and cleared
    by a call that clears it. A re-creation starts with no copy and gains the mirrored mask when R puts it
    back through the bridge's mask entry.
- Push 2: `Sampler::rederiveStandardization()`, after slice C's `Sampler::rederiveScale`, re-derives every
  leaf covariate over the rows in, converting live and kept linear slopes (dec-B329) and returning which it
  kept; a facade virtual forwards it. `Sampler::setWeights` and `Sampler::setActiveRows` are not edited.
- `ColumnStore::buildFromParent` stops copying the parent's constants into a view, so a fold's leaf
  derives over the fold's rows by the creation rule. Whatever slice 2 serves through
  `ColumnStore::suppliedStandardization` for a carried record stays.
- Facade: push 1 changes the setData virtual's signature (the mask); push 2 adds one virtual.
  `--preclean` on both.

Bridge.
- `createHolder` in [R_interface_bartcore.cpp](../../src/R_interface_bartcore.cpp) hands forest-multiplier
  push 2's creation mask to `SamplerOptions` for one forest as it hands it to `AmplitudeSpec` for several,
  checked once.
- [`bartcore_setData`](../../src/R_interface_bartcore.cpp) takes the mask as a third argument.
- The refusal: the engine's error, worded with the column's name and the counts, raised as an R error; the
  data handle's fold build raises the fold text.
- The leaf's fallback flags reach R through the reader slice 2 adds; the predictor entries that take
  `updateStandardization` answer "warn now" on a sampler's first keep, from the holder's one
  `scaleKeptWarned` flag, which the four response and row setters already share (response-scale-rows-c
  rule 5).
- Push 2: [`bartcore_setWeights`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setActiveRows`](../../src/R_interface_bartcore.cpp) take a fourth argument. With it true, after
  their refusals and the gp kept-draws refusal ([`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp)'s
  pattern, dec-B237's text), the install, and the scale's re-derivation where asked, they call
  `rederiveStandardization` and fold its keep into the logical "warn now" they return.
- No other entry changes arity; no flat entry, [dbarts.h](../../inst/include/dbarts/dbarts.h) or token
  moves.

R.
- [`bartcoreSamplerSetData`](../../R/bartcore.R) passes the mask it resolves from the new weights to the
  bridge before the engine swaps, and still installs it afterwards as today.
- The creation warning for no row in: where the engine reports a leaf covariate read from every row, the
  sampler's `initialize` raises one warning of class `dbartsScaleFallbackWarning` under every family, in
  forest-multiplier push 2's form: under gaussian, Student-t and nbinom one joint text naming the response,
  the prior's center and the leaf covariates (and a multiplier's column, after that push); under probit,
  ordinal and hazard.probit "no observation is in the likelihood (0 of n have a positive weight), so the
  center and scale of x1 are taken from all n". `setData` raises the same, last, as it raises the range's.
- The setters raise [`warnScaleKept`](../../R/bartcore.R)'s warning last on a first keep, after the
  mirror, the record and the store.
- Push 2: `setWeights` and `setActiveRows` in [dbarts.R](../../R/dbarts.R) take `updateStandardization`
  after `updateScale`, checked as `updateScale` is ([`checkUpdateScale`](../../R/bartcore.R)), and write the
  record after a re-derivation; `setWeights` forwards both to `setActiveRows` on a mask family and keeps its
  rollback on a refused call. `reapplyActiveRows` passes FALSE.
- Creation needs nothing else past forest-multiplier push 2: `dbarts()` already hands `active` to `new()`,
  which hands it to the bridge.

## Constraints

- With every row in and no leaf covariate with too little, every derivation runs today's code on today's
  inputs: no such fixture moves a bit, the base-build literals of state-frame-prior's tests included;
  `xbart` with a leaf covariate is the one exception, its folds reading their own rows.
- The row set is membership, never the weights' sizes.
- A refused call creates nothing or changes nothing: the data object, the mirrors, the record and the
  stored state are as they were.
- No new slot or field; no state format change. The new arguments are push 2's two; the new texts are the
  refusal's three and the leaf's part of the creation warning, in the existing class; the keep's warning
  reaches three more setters.
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

Estimates are the raw itemization; [Budget and stops](#budget-and-stops) gives the forecast.

### Push 1: the rows, the refusal and the keep

1.1 Engine, the derivation (~75). The row set and the too-little answer in
    `standardizationMomentsForColumn`; the gp heuristic over the rows read; the leaves' `initialize` and
    `reinitialize` take the row set and the replaced columns, throw the typed error at a derivation (the
    stated-lengthscale gp exempt), keep at a re-derivation, and record the no-row flag for the reader.
1.2 Engine, the doors (~95). `SamplerOptions` gains the creation mask; `Sampler` forms the row set at
    construction; `Sampler::setData` and the facade virtual take the mask, replace the sampler's copy, test
    the derivation before the swap and hand the row set to
    [`Chain::applyNewData`](../../src/bartcore/chain.hpp)'s derive arm; the sampler's copy of the mask by
    its lifecycle above; the predictor update's derive arm reads that copy, re-derives the replaced columns
    only and keeps; the view's parent constants go (~15 removed); the fold's build raises the refusal.
    `--preclean`.
1.3 Bridge (~70). The creation mask to `SamplerOptions`; `bartcore_setData`'s third argument and its
    registration; the refusal's wording at creation, setData and the fold; the fallback flags in slice 2's
    reader; the predictor entries' "warn now".
1.4 R (~90). `bartcoreSamplerSetData` passes the mask; the creation and setData warnings for no row in
    under every family, folded into forest-multiplier's text; `setPredictor`'s keep warning; the docstrings
    of `setData` and `setPredictor`, the first naming the mask `s$data` does not carry.
1.5 tests/cpp (~330). `testLeafStandardizationRows` beside the leaf tests: the moments over a row set
    equal literals; null and all-ones row sets equal today's bits; one row in, a column constant over the
    rows in, a column constant over every row and NA-only rows in each answer too little; no row in reads
    every row and sets the flag; a derivation with too little throws naming the column, a stated-lengthscale
    gp takes the placeholder, a re-derivation keeps and reports; the heuristic over more than 512 rows in,
    interleaved with rows out and with an NA among them, equals the literal of the stride over the rows
    read. A chain created with weight-0 rows and one created with the same numbers supplied through slice
    2's path draw identically for 20 sweeps; the same with the creation mask under probit and nbinom.
    `setData` with a mask, and with a mask after `setActiveRows` at another row count: the numbers over
    the rows the new mask keeps; a refused `setData`: the sampler bit for bit an untouched twin. A predictor
    update after `setActiveRows`: the mask in force, not the creation one. The view pins of test_data.cpp
    ("subset view gathers a leaf covariate through the column map", and `gatheredRawAgrees`'s "the
    constants are the parent's full-data ones") turn into: a view supplies no constants and a fold's leaf
    derives over its own rows.
1.6 tinytest (~800 new, ~240 edited). A new file test-leaf-standardization-rows.R ([Tests](#tests)). Edits
    where a pin expects every row or a constant covariate created: slice 2's "fresh object neutral"
    literals for a leaf fixture with a row out, each re-recorded and named in the landing note;
    test-gp-leaves.R's zero-weight sampler, which warns now where it did not; the `checkDataHandleViews`
    callers of [leafPriorChecks.R](../../inst/common/leafPriorChecks.R) and the xbart arms of
    test-linear-leaves.R and test-gp-leaves.R where they pin a number; test-boundary-inputs.R's constant
    designated column, now refused under linear and gp and created under gp with a stated lengthscale;
    sections 1 and 6 of test-leaf-conversions.R, whose samplers on constant columns reach a covariate
    without spread through a carried record whose scale is NA, or a stated-lengthscale gp, instead.
1.7 The exact arm (~80). [linear-exact.R](../../benchmarks/R/linear-exact.R) gains an arm with rows at
    weight 0 on the existing four values in another mix (so the tree space does not grow), sized so that
    the oracle standardized over every row fails it; the oracle standardizes over the rows in. Its header
    says what it does not see: the gp heuristic, the masks, the refusal and the keep, which the literals of
    1.5 and 1.6 hold.
1.8 Help and NEWS (~100). man/dbarts.Rd's `leaf.prior` item ("standardized internally" over the rows in at
    creation, held; refused where the rows in hold one value of a covariate, a constant column included,
    unless a gp lengthscale is stated); dbartsSampler-class.Rd's `setData` (the mask `s$data` does not carry
    under a mask family), `setPredictor` (the replaced columns only; the keep and its warning), the
    `updateStandardization` items and the manual's table row for the standardization; man/xbart.Rd (each
    fold its own rows, and the refusal); man/dbartsSpec.Rd. NEWS: a linear or gp leaf reads its covariates'
    center, scale and lengthscale from the rows in the likelihood; a constant leaf covariate is refused.
1.9 Records (~150). The design note [response-scale-rows.md](../design/response-scale-rows.md): rule 9 for
    the leaf, section 2's table, section 4 for the leaf's no-row flag and keep, section 8's slice D;
    [state-frame-prior.md](../design/state-frame-prior.md)'s "What response-scale-rows slice D gets" closed
    by a dated line; TODO; the ledger entry for the calls made; the indexes; this plan's landing note.

### Push 2: the burn-in arguments

2.1 Engine (~45). `Sampler::rederiveStandardization` and its facade virtual, after `rederiveScale`;
    `--preclean`. tests/cpp `testRederiveStandardization`: after `setWeights` or `setActiveRows`, the call
    leaves the sampler bit for bit where a twin given `setPredictor` of every leaf column in force with
    `updateStandardization` true leaves it; a keep reports the covariate and moves none of its numbers.
2.2 Bridge (~35). The fourth argument of `bartcore_setWeights` and `bartcore_setActiveRows`, their
    registration, the gp refusal before the install, the folded "warn now".
2.3 R (~60). The two methods' argument, its check and position, the forwarding on a mask family, the record,
    the store, the warning last; `reapplyActiveRows` passes FALSE.
2.4 tinytest (~290). In the new file:
    - the identity: for gaussian, Student-t, probit, ordinal, nbinom and hazard.probit, linear and gp with no
      kept draw, the call with `updateStandardization = TRUE` against the plain call and then
      `setPredictor` of every leaf column in force with the argument: `identical()` in the record, the stored
      state and 5 sweeps of 2 chains; with kept linear draws, `predict` unchanged to 1e-12;
    - both arguments TRUE: the range and the standardization each the literal over the rows in;
    - the argument: absent and FALSE move nothing; positional `setWeights(w, TRUE)` still names
      `updateState`; `updateStandardization = NA` and `= 1` refused with nothing installed;
    - the keep: rows in holding one value of x1 and spread of x2: x1 kept, x2 re-derived, one warning, none
      on a second keep or after a keep through `setOffset`;
    - gp with a kept draw: refused before the install, the weights and the mask as they were;
    - the burn-in loop of the help, by `setWeights` and by `setActiveRows`: during burn-in the record is the
      literal over each pass's rows in, after it the last burn-in pass's, and every kept draw replays under
      it; the loop equals a twin written with the identity's calls.
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
- Neutral: every row in and every covariate spread, each leaf kind and family: the record and 20 sweeps
  equal literals recorded on the base build.
- No row in, each family that can have one (gaussian, Student-t, probit, ordinal, nbinom, hazard.probit):
  every row, exactly one warning of the class, its text naming the leaf covariates; under gaussian and
  nbinom the one joint text, not two; then a covariate constant over every row refused.
- Refused at creation, with the text: rows in at one value of x1 under linear and under gp; x1 constant
  over every row; x1 missing on every row in; one row in; at `dbarts()`, `bart()`, `dbartsSpec()` with
  `new()`, and `new()` with `active`; nothing created. Created: the same under gp with a stated lengthscale,
  its record carrying the value and NA, and the next 5 sweeps finite.
- NA among the rows in: the observed rows in only for the moments; the heuristic counting them at 0.
- By hand: `new(control, model, data, active = spec$active)` is `dbarts()`'s under probit and nbinom with a
  linear leaf, record and 5 sweeps; without `active`, every row (the help's sentence).
- Held: `setWeights`, `setActiveRows`, each with and without `updateScale = TRUE`, `setResponse` and
  `setOffset` with `updateScale = TRUE`, `updatePredictor` and a per-observation session: the record
  `identical()` and 5 sweeps equal a twin's.
- `setPredictor(x[, j], j, updateStandardization = TRUE)`:
  - after `setWeights(w2)`: the literals over the rows `w2` keeps for j, live and kept linear predictions
    unchanged to 1e-12 (dec-B329, slice 3's "linear converts");
  - a second leaf covariate k held, `identical()`;
  - after `setActiveRows(a)` under gaussian and under probit: the rows `a` keeps;
  - after `setActiveRows(a)` then `setData` to a data object of another row count: the new data's rows in,
    and no read past them;
  - on a `copy()` of a sampler holding a mask: the mask put back on the copy;
  - no row in, and j at one value over the rows in: kept, exactly one warning on the sampler's first keep,
    none on a second, none after a keep through `setOffset` (the shared flag); FALSE holds.
- `setData`: a fresh data object with weight-0 rows, and one with 0/1 weights under nbinom: its rows in; a
  carried record: held, also one whose scale is NA; TRUE: rows in; FALSE: the sampler's (slice 3's
  "arguments win" over a row out); a fresh object whose rows in hold one value of x1, at NULL and at TRUE:
  refused, the sampler as it was; under probit, `setData(s$data, updateStandardization = TRUE)`: the mask
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
| every row read where one row in (dec-B309's fallback) instead of the refusal | the one-row-in refusal |
| no fallback where no row is in | tests/cpp's flag or the no-row literal |
| the leaf's no-row fallback not warned under probit | the no-row arm under probit and ordinal |
| two warnings where response and leaf both fall back | the joint-text count |
| setPredictor re-deriving every leaf covariate | the second covariate k held |
| setPredictor's keep re-deriving over every row, or warned on every keep | its keep arm; the second-keep count |
| the keep's flag not shared with the response setters | the `setOffset` then `setPredictor` count |
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
| a mask family's `setWeights` not forwarding it | the probit and nbinom identity |
| the record not written after a re-derivation | the burn-in loop's copy and reload |
| `reapplyActiveRows` re-deriving | the record arms after a reload |

Expected equivalent, said in the landing note if alive: the row-set loop taken with an all-ones set where
none is out (same sums in the same order).

## Verification

Per push, against its own library:

- tests/cpp plain and under ASAN and UBSAN; `R CMD INSTALL --preclean`; the full suite; the R-loaded address
  sanitizer over the new file, test-gp-leaves.R, test-linear-leaves.R, test-leaf-conversions.R,
  test-boundary-inputs.R and test-data-handle.R; the x86-64 leg (tests/cpp under the sanitizers and the
  suite) before landing.
- Reference build: the four snapshot files and the three bitwise compares, counted, as expected above.
- exact-gates.yaml's commands at quick, push 1's new arm among them; the arm's every-row mutant run once to
  show it fails.
- Lint, air, rc-codoc, win-drift and doc-freshness, each on its own exit status; inst/NEWS.Rd parses with a
  non-NULL result; `R CMD check --as-cran` from a clean tarball.
- Consumers read, not run: bartCause fits several forests (refused leaf covariates), stan4bart's mask goes
  on after creation, treatSens and bairrtt pass no linear or gp leaf (searched).

## Budget and stops

| part | push 1 raw | push 1 forecast | push 2 raw | push 2 forecast |
|---|---|---|---|---|
| engine and facade | ~170 | ~270 | ~25 | ~40 |
| bridge | ~70 | ~110 | ~35 | ~55 |
| R | ~90 | ~145 | ~30 | ~45 |
| tests/cpp | ~330 | ~660 | ~20 | ~40 |
| tinytest | ~1040 | ~1870 | ~290 | ~520 |
| exact arm | ~80 | ~130 | - | - |
| help and NEWS | ~100 | ~130 | ~50 | ~65 |
| records | ~150 | ~195 | ~50 | ~65 |
| total | ~2030 | ~3500 | ~500 | ~830 |

Package code: push 1 ~330 raw, ~525 forecast; push 2 ~90 raw, ~140 forecast. dec-B465 estimated the
arguments at 50 to 100 lines beside slice D; push 2's package code is at the top of that and its tests
beyond it, the identity and the loop being where its defects would show.

From the landed history: slice B planned 2065 lines and landed 3830 (1.85 times); slice C planned 1120 and
landed 2537 (2.27), its package code 165 against 280 and its tests 710 against 1530; state-grid-install
1475 against 2092 (1.42). Tests are forecast at about 1.8 to 2 times, package code at 1.6, the rest at 1.3.
Against the first version (~1460 raw): the critique's warnings and mask lifecycle add about 370, dec-B465's
refusal and the rewritten constant-column tests about 200, push 2 about 500.

Stop and report, building nothing further, when:
- push 1's diff passes ~5250 lines or its package code ~800; push 2's passes ~1250 or ~210;
- a bitwise compare prints "max |z|", a snapshot moves, or an exact gate other than the new arm moves;
- a fixture with every row in and every covariate spread moves a bit outside `xbart`, or a literal of
  state-frame-prior's tests moves other than one with a linear or gp leaf and a row out or a constant
  covariate;
- a test outside those named in 1.6 fails because it created a sampler on a constant leaf covariate (list
  it, do not edit it);
- the derivation needs a line in a response model, a state block, a flat entry or dbarts.h, or a facade
  change beyond push 1's setData signature and push 2's one virtual;
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
- One row in is "one value" and refused, though the response reads every row there (dec-B309). Not taken:
  dec-B309's every-row fallback for the leaf, which dec-B465 replaces for one value; no row in at all keeps
  dec-B304's every row.
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

On the tip after state-frame-prior slice 3, before writing a line:

1. Slice 3's path for `updateStandardization = TRUE` on `setPredictor` and the predictor entries: which
   function it calls, and whether it re-derives the replaced columns only or every leaf covariate (through
   `reinitialize`, it would be every one). This slice pins the replaced columns only; if slice 3 re-derives
   every one, step 1.2 narrows it and slice 3's tests that pin otherwise are listed, not edited, and stop the
   build for the coordinator. Push 2's identity is written against whatever that path is after 1.2.
2. The record's name and shape: `data@standardization`, per forest, `center` and `scale` named by
   covariate with NA for a column without spread, `lengthscale` only where the model supplies none; its
   accessor and check in R/data.R; that a re-creation writes back a carried heuristic lengthscale (not
   taking it as a supplied one and dropping it); and that its check accepts an NA scale.
3. Where a carried standardization enters the engine: slice 2 plans
   `ColumnStore::suppliedStandardization` on a top-level store. If it is a construction argument instead,
   step 1.2's removal of the view's parent constants is unchanged; if it reuses the view's
   gathered-constants storage, the removal keeps that storage and drops only the parent derivation.
4. Slice 3's `Sampler::setData`: the per-value choice, the facade virtual's signature and where its
   refusals sit against the swap. If slice 3 already passes the mask, step 1.2's signature change goes.
5. Forest-multiplier push 2 has landed: `bartcore_create`'s fifth argument, `new(..., active = )`,
   `dbarts()` passing it, the fallback warning raised by `initialize` from the engine's status, and
   `BasisScaleError`'s route for a refusal naming a column. If not, push 1 adds them for one forest (~80 more
   lines) and that push rebases onto them.
6. The baselines named in benchmarks/baselines/MANIFEST, and slice 2's 16 fixtures: which carry a linear or
   gp leaf with a row out or a constant covariate, the only literals expected to move.
7. Whether `ColumnStore::buildFromParent` still copies parent constants after slices 1 and 2.
8. The suite's creations on a constant leaf covariate beyond test-boundary-inputs.R and sections 1 and 6 of
   test-leaf-conversions.R: grep for `linear(` and `gp(` beside a constant column before building, and size
   step 1.6 from the count.

## Not verified

Nothing is built. The rows-in numbers were emulated on the tip by a `setState` install of the
standardization, which state-frame-prior slice 2 removes; the draws are therefore not slice D's. The xbart
fold's constants were read, not measured: no R reader reaches a fold. The heuristic's subsample over rows in
was read, not run. Which tests create a sampler on a constant leaf covariate was read for three files, not
counted over the suite. The burn-in loop was run for `updateScale` only. The baselines' verdicts are read
from the scenario sources and the MANIFEST, not run; no reference build, sanitizer, x86 run, `R CMD check`
or consumer suite was run. The line estimates are an itemization scaled by three landings' ratios.

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
