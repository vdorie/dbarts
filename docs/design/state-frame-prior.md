# Derived data on the data object: the cut grid, the leaf standardization and the response range

Status: PROPOSED 2026-10-17, no code. Plan: [state-frame-prior.md](../plans/state-frame-prior.md). Rulings:
dec-B231, dec-B232 and dec-B233, with dec-B200 as dec-B418 revises it, dec-B230, dec-B237, dec-B304, dec-B449 and
dec-B450 in docs/decisions.md; dec-A159 (the forest calibration anchor goes to the data object beside the
range) was confirmed under dec-B378. Calls no ruling covers are [Open questions](#open-questions), each with a
recommendation; nothing below settles them.

A sampler holds a few numbers it derived from the data when it was made and then keeps while the data move:
the cut grid and its weights, a linear or Gaussian-process (gp) leaf's covariate center and scale with the
default lengthscale derived with them, the response range, and, in a model of several forests, the
calibration anchor and the unit a forest's sd is stated in. dec-B231 classes them as data, "albeit one that
is initially derived from the user input". So they belong on the data object, change only through a call that
names them, and are never installed by a state. Today the engine holds them, and the only other copy is
scattered: the grid in the stored state, the range on the model, the anchor on the control, the
standardization in the state again. This note gives each one record on `dbartsData` and says what every door
does with it.

## Today

On 6fe5237a. "Held" means the call keeps the value in force.

| value | derived | `setData` | `setPredictor` | live `setState` | copy, reload | recorded |
|---|---|---|---|---|---|---|
| cut grid and its weights | at creation, every row | derived again | held unless `updateCutPoints` | never written (dec-B449) | taken from the state | the state's `cutPoints` and `cutMass` |
| leaf center, scale; heuristic lengthscale | at creation, every row | derived again; linear slopes converted | held, no argument | installed from the state | derived, then the state's installed | the state, per forest |
| response range | at creation, rows in | derived again | held | never moved (dec-B418) | the model's record applied | the model's `response.range` attribute |
| forest anchor and unit | at creation | refused | held | not in a state | the control's record | the control's `bartcore.forests` attribute |

Two values already live on `dbartsData` and show the pattern this note extends: `sigma`, the starting estimate
`setData` carries over, and `missing.seen`, written from the engine after each predictor change and read when
a sampler is created, so a copy and a reload see what their source saw ([`dataMissingSeen`](../../R/data.R)).

A gp sampler that holds kept draws refuses `setData` outright
([`bartcore_setData`](../../src/R_interface_bartcore.cpp)), since the call always re-derives the
standardization and the range and a kept gp draw replays only under the kernel and shift it was drawn with
(dec-B237; [leaf-conversions.md](leaf-conversions.md)).

## The records

Four slots on [`dbartsData`](../../R/A_class.R), each `NULL` where absent. Each is read through an accessor
that answers `NULL` for an object saved before the slot existed, as [`dataMissingSeen`](../../R/data.R) is,
and checked by a plain function at the doors that use it, as [`checkBasisRecord`](../../R/data.R) is, never by
the class's validity. Names: the first three are dec-B233's.

| slot | holds |
|---|---|
| `cutPoints` | a list, one entry per predictor column in the object's column order: the strictly increasing points; `NULL` for an unordered factor. A column whose points carry weights has them as the entry's `mass` attribute, the increasing sequence a state's `cutMass` holds for it |
| `standardization` | a list, one entry per forest, `NULL` where the forest's leaves take no covariate: `center`, `scale` (`NA` for a column without spread) and `lengthscale` (`NA` where the model supplies one, or on a linear leaf), each named by the leaf covariates |
| `response.range` | `c(min, max)`, what [`bartcore_anchor`](../../src/R_interface_bartcore.cpp) reports today; a constant response `c(c, c)` |
| `forest.calibration` | `c(anchor = , unit = )` on a model of several forests, what the control's attribute holds today |

The record lives in three places, each holding what its layer can:

- R. The slots on `sampler$data` are the record a re-creation reads. The sampler writes them from the engine
  after creation and after every call that can change one. The object a caller handed in is never written:
  two fits made from one `dbartsData()` object each derive their own.
- Bridge. [`parseData`](../../src/R_interface_bartcore.cpp) reads the slots into plain vectors and hands them
  to construction and to `setData`; one reader returns all four from the engine for R to store. Every R type
  stops here.
- Engine. It keeps the values where it keeps them now: the grid and its weights on the column store, the
  center, scale and lengthscale on the leaf ([`LinearGaussianLeaf`](../../src/bartcore/model.hpp),
  [`GPGaussianLeaf`](../../src/bartcore/model.hpp)), the range at
  [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp), the anchor and unit in the calibration map. It gains
  construction from supplied values (each optional, absent derives), a `setData` that derives, keeps or takes
  each value, and the one reader. Nothing in it names R; [dbarts.h](../../inst/include/dbarts/dbarts.h) does not
  change, no flat entry replacing data or installing a state.

A state keeps its `cutPoints`, `cutMass` and standardization blocks: they record what its thresholds and
slopes were stored under, which an install compares. [`stateFormatVersion`](../../src/R_interface_bartcore.cpp)
does not move. The model loses its `response.range` attribute and the control its `anchor` and `unit`.

## The rule at every door

dec-B231's rule, in the maintainer's words: "re-derive if no values are present (to recreate old behavior)
and if values are supplied, use them. If explicit arguments are given, we could then follow them directly."

- Creation: `dbarts()`, `bart()`, `bart2()`, `dbartsSpec()` and `new("dbartsSampler", control, model, data)`.
  Each record the data object carries is used; an absent one is derived. The fitting functions build a
  fresh data object, so they derive as they do now. A record that does not fit the data (another column
  count or kind, a grid not strictly increasing, an ordered factor's grid other than its levels', a
  standardization for columns the leaf does not take) is refused, naming the slot and saying that setting it
  to `NULL` derives it. dec-B232 names the one cost: a range derived under one family and offset is used
  under another.
- `setData(newData, updateState = NULL, updateCutPoints = NULL, updateScale = NULL, updateStandardization =
  NULL)`. Each record separately: an argument given is followed (derive, or keep what the sampler holds); with
  none, a record `newData` carries is used, and one it lacks is derived. Where `setData` puts the columns in
  the sampler's order (dec-A215) the grid's entries move with their columns. A linear leaf's slopes, live and
  kept, are converted wherever the standardization changes, and kept draws are rewritten wherever the range
  changes, as today ([`Chain::applyNewData`](../../src/bartcore/chain.hpp)). The argument values are
  question 1 of [Open questions](#open-questions).
- `setPredictor(x, column, forceUpdate, updateCutPoints = "none", updateState = NULL, updateStandardization =
  FALSE)`. `TRUE` derives the center, scale and heuristic lengthscale of the leaf covariates over the new
  values and converts linear slopes, live and kept; a gp sampler holding kept draws refuses `TRUE` (dec-B237).
  dec-B233 has it default to keeping them.
- `setCutPoints`, and `updateScale = TRUE` on `setResponse`, `setOffset`, `setWeights` and `setActiveRows`:
  unchanged, except that the record they write is the data object's.
- A live `setState`. Writes no record. The grid: dec-B449 as landed. The standardization: compared forest by
  forest. Equal, the install is today's, bit for bit. Different, the sampler keeps its own; a linear leaf's
  live and kept coefficients are converted from the state's standardization into it
  ([`convertLinearCoefficients`](../../src/bartcore/model.hpp),
  [`convertFlatLinearTree`](../../src/bartcore/chain.hpp)), so every leaf is the function it was stored as, as
  a split keeps its value; a gp leaf's live fits go in as stored, and kept gp draws made under another center,
  scale or heuristic lengthscale are refused, as kept draws under another supplied lengthscale are
  ([`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp)). The range: as stored (dec-B418).
  Whether the conversion is clean is question 3 of [Open questions](#open-questions).
- `copy()`, a reload's first use, and `setState` on a sampler whose engine is gone. The engine is re-created
  from `(control, model, data)` by creation's rule, so it takes every record from the data object and none
  from the state, and the stored state then goes in by the live rule, forced as these always are (dec-B305).
  The grid adopted from the state ([`installStateOnto`](../../R/dbarts.R)'s flag, `adoptGrid` in
  [`Sampler::setState`](../../src/bartcore/sampler.hpp)) remains only for an object saved without the record
  (question 6 of [Open questions](#open-questions)). The field refresh after an install onto another grid
  ([`heldAfterAdoption`](../../R/dbarts.R)) stays.
- A warm start: unchanged. The recipient keeps its own values; a donor's slopes are converted from the
  standardization its state records ([`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp)).
- `xbart`: unchanged. Its folds are views of one store inside the engine and take the full data's
  standardization ([`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp)); dec-B304 moves them
  to each fold's own rows in response-scale-rows slice D.
- Saved objects. A 0.9-34 sampler is refused at its first use ([`refuseLegacyState`](../../R/tombstones.R)), as
  today; a 0.9-34 data object has no slots and derives, as 0.9-34 did. An object saved by an earlier bartcore
  build is question 6 of [Open questions](#open-questions).
- The flat C API has no record. A host that re-anchors through it changes the engine and leaves the R record
  stale, as it leaves the model's record today ([state-not-model.md](state-not-model.md)).

dec-B230 is untouched: a data object, like a state, is not checked against the data it came from.

A gp sampler holding kept draws is then refused by `setData` only where the call would change the
standardization or the range in force, and served where it keeps both; the grid never refuses, kept draws
holding their thresholds by value. Which calls count as a change is question 4 of [Open questions](#open-questions).

## What changes for users

Runs in scratch/sfd (p1 to p3) on the tip's shipped build and on 0.9-34; 200 rows, two uniform columns, 20
trees, one chain. "Held" was emulated on the tip with the calls that exist: `setCutPoints` by value after
`setData`, and the engine's anchor entry with the old range.

Against 0.9-34:

- `setData` handed the sampler's own data object, edited (`d <- s$data; d@x[, 1] <- ...; s$setData(d)`), now
  holds what 0.9-34 derived again, the object carrying the sampler's record. Column 1 doubled: 0.9-34 and the
  tip move its grid from [0.023, 0.983] to [0.046, 1.966]; held, it stays on [0.023, 0.983], and the posterior
  mean fit after 300 draws is 0.177 from the truth against 0.081 (fits up to 0.52 apart). Response replaced by
  3 y + 10: the range moves from [0.038, 2.294] to [10.11, 16.88]; held, the fit is 0.324 from 3 f + 10
  against 0.244 (up to 0.72 apart). A fresh `dbartsData()` carries no record and derives, as 0.9-34 did; that
  is the form the equivalence harness uses. Of the suite's 126 `setData` calls, about 20 hand back a sampler's
  data object or an edit of it, most to pin a refusal that comes first; the 3 in test-response-scale-rows.R
  pin the derivation's fallback warning and change.
  Question 2 of [Open questions](#open-questions).
- `new("dbartsSampler", control, model, s$data)` and `dbarts(s$data)` take the grid and range the sampler
  held; 0.9-34 derived them.
- `copy()` and a reload whose stored state predates a grid change: 0.9-34 took the state's grid (and over-read
  `numCuts`, state-grid-install.md); now the sampler's.

Against the tip, besides the three above (the tip behaves as 0.9-34 in the first two):

- A copy or reload after `setCutPoints(c(0.25, 0.5, 0.75), 1)` with `updateState = FALSE`: the tip's copy has
  the state's 100 points and its 17 splits on the column, where the sampler has 3 points and 16; under the
  record the copy has the 3 points and the stored trees moved onto them by value, 16 splits. With
  `updateState = TRUE` both give 3 points and 16 splits, unchanged.
- A forced `setState` across a change of standardization (linear leaf, a column replaced by 3 a + 5 through
  `setData`): the tip installs the state's center and scale (0.518, 0.269) over data whose own are (6.55,
  0.807), which makes the slope's prior three times as wide in the column's units; under the record the
  sampler keeps (6.55, 0.807) and the slopes are converted. Kept draws predict the same either way (the tip 0
  apart, the conversion 1e-14 as at `setData`).
- A gp sampler holding kept draws, `setData(s$data)`: refused on the tip, served.
- `attr(s$model, "response.range")` and the control's `anchor` and `unit` are gone; 20 lines in 13 test files
  and 2 lines read them.

Nothing moves for a fit that never calls `setData` on a carried record, `setState` across a standardization,
or copies a sampler whose stored state predates a grid change: every creation from a fresh data object
derives what it derives today, and a re-creation from a current record gives the values today's re-creation
reaches by adopting from the state.

## The manual's table

The help's Saving section gains one table, which says what puts each value back before `setState`. A data
object held in a variable is the general answer: `setData(d)` puts back the observations and, because `d`
carries the record, every derived value with them, weights included.

| value | changed by | put back by |
|---|---|---|
| cut points | `setCutPoints`; `setPredictor(updateCutPoints = )`; `setData` deriving | `setCutPoints(d@cutPoints)`, the whole set; or `setData(d)` |
| their weights | the same | `setPredictor(x[, j], j, updateCutPoints = "value")` while the column holds the values they came from; or `setData(d)` |
| leaf center, scale, lengthscale | `setPredictor(updateStandardization = TRUE)`; `setData` deriving | `setPredictor` of the old values with `updateStandardization = TRUE`; or `setData(d)` |
| response range | `updateScale = TRUE` on the four setters; `setData` deriving | the same setter on the old response, offset, weights or mask with `updateScale = TRUE`; or `setData(d)` |
| forest anchor and unit | nothing | - |

`setData(d)` also gives dec-B449's weighted case its by-hand answer without weights on `setCutPoints`, which
dec-B450 defers past the merge.

## The engine and the bridge

- Construction takes optional supplied values beside the options it takes now: per column a grid and its
  weights, per forest a center, scale and lengthscale, a range, an anchor and unit. A grid goes in before the
  columns are quantized, so a carried grid costs no second pass; a carried standardization goes in where
  [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp) now serves a view; a range through
  the path `applyAnchor` uses. Absent derives, through today's code.
- [`Sampler::setData`](../../src/bartcore/sampler.hpp) takes, per value, derive, keep or a supplied value.
  Keep re-quantizes the new rows on the grid in force, regathers the leaf covariates under the constants in
  force for the new row count, and installs the response without moving the transform. Derive is today's
  path, unchanged.
- [`Chain::setState`](../../src/bartcore/chain.hpp) stops calling `restoreCalibration` on a live install
  and converts instead; `adoptGrid` in [`Sampler::setState`](../../src/bartcore/sampler.hpp) stays for the
  legacy path only.
- One reader of the four values, read by the bridge after every call that can change one.
- The bridge resolves the R arguments and the carried slots into the per-value choice; the engine does not
  know which door asked.

## What response-scale-rows slice D gets

dec-B304 has a leaf's center and scale taken over the rows in at creation, held, and re-derived by
`updateStandardization`; slice D waits because "a state install today replaces the recipient's center and
scale and ... a re-creation without a state derives them afresh, so there is nowhere to hold them". After
this item a copy, a reload and a state install hold them through the record, and `updateStandardization`
exists on both setters. What D still owns: the rows (creation's leaf derivation reads the weights or mask; the
column store holds none), the fold views in `xbart`, and `updateStandardization = TRUE` reading the rows in at
that call. D lands after the slices that give the record and the argument (the plan's order).

## Alternatives

- Records on the model, where dec-B200 put the range: rejected by dec-B232, a model handed to other data
  carrying numbers derived from one response.
- Records in the state, as the grid is today: a restore would then move a prior (dec-B231, dec-B449).
- A hidden field on the sampler: every re-creation path would need its own carrying, and a host could not hand
  a sampler a range or grid with its data. The data object is already what creation, copy and reload read.
- Holding at every plain `setData`, or requiring the arguments on every call: dec-B231's two rejected options.
- `setState` installing the standardization, as the tip does: no tree changes, but restoring a state moves the
  slope prior and kernel, which is the model; and a gp sampler takes kept draws whose kernel then reads under
  another frame, which dec-B237 refuses elsewhere.

## Open questions

Each for the maintainer; the plan is written to the recommendation and names what a different answer
changes.

1. setData's argument values. dec-B233 names them, not their values. (a) `NULL` by default, the rule;
   `updateCutPoints` takes `setPredictor`'s words, `"none"` keeping the sampler's grid and `"position"` or
   `"value"` deriving and placing the splits, today's `setData` being `"value"`; `updateScale` and
   `updateStandardization` take `TRUE` (derive) or `FALSE` (keep the sampler's). An explicit keep is the
   sampler's value, not the carried one, which `NULL` alone uses. The three go after `updateState`, so 0.9-34's
   `setData(d, TRUE)` keeps its meaning. (b) Logical throughout, `updateCutPoints = TRUE` placing by value.
   (c) A fourth value naming the carried record. Recommended: (a), one vocabulary for the grid across the three
   setters and `NULL` as base R's absent; (b) loses `"position"`, (c) adds a word for what `NULL` already does.
2. The edited `s$data` idiom of [What changes for users](#what-changes-for-users). The rule holds a
   carried record, which is what a host that moves predictors each sweep wants and what dec-B231 rules; the cost is a caller who edits the sampler's own
   object expecting 0.9-34's re-derivation. (a) Hold, with a NEWS line and the help naming `updateCutPoints =
   "value"`, `updateScale = TRUE`, `updateStandardization = TRUE` or a `NULL` slot to derive. (b) As (a) and
   a warning when a held grid leaves a column's new values beyond its last point, or a held range leaves the
   new response outside it. (c) Read a carried record identical to the sampler's own as absent: the idiom
   derives as in 0.9-34, but a host handing back `s$data` with a new response then re-derives each sweep,
   which dec-B296 rules out. Recommended: (a); (b) fires in a host's loop where the hold is intended, and
   costs about 40 lines and a warning class.
3. A state under another standardization. That a state does not install the standardization is the
   register's conclusion from dec-B231, not the maintainer's words; dec-B449 ruled it for the grid alone. (a)
   The sampler keeps its own, linear coefficients live and kept are converted, the install is clean (`TRUE`
   unforced): no tree changes and every leaf is the function it was, as on a finer grid (state-grid-install,
   "other weights on equal points are clean"). (b) The same, but not clean: `FALSE` unforced, a conversion
   only when forced. (c) Install it, as the tip does. Recommended: (a). (b) leaves no unforced way to restore
   a state stored before `updateStandardization = TRUE` on a linear sampler; (c) moves the model on a restore. Kept draws converted, not as
   stored: dec-B418's "never convert" is about the response's internal scale; a slope's covariate units are
   the analog of a split's value, which moves by value.
4. The gp `setData` refusal. (a) Refused where the call derives the standardization or the range (an argument
   `TRUE`, or `NULL` with no record carried) or brings a carried one other than the sampler's; served
   otherwise. (b) Refused only where the values in force would differ, compared after deriving. Recommended:
   (a), as the `updateScale` refusal reads no values (leaf-conversions.md): a derivation that happens to land
   on the same numbers would otherwise pass on one sweep and stop on the next.
5. Copy and reload when the stored state predates a grid change. The record makes the copy's grid the
   sampler's and moves the stored trees onto it by value ([What changes for users](#what-changes-for-users)).
   (a) That, with the help saying to `storeState()` before copying or saving after a grid change. (b) As (a), and `copy()` reads the live
   state where the source's engine is alive, so a copy is the sampler as it stands; a sampler kept with
   `updateState = FALSE` would no longer copy its stored state. Recommended: (a); (b) changes what `copy()`
   is for every `updateState = FALSE` user, for a case the help can prevent.
6. Objects saved by an earlier bartcore build, whose data object has no record. (a) Migrate at the first
   re-creation: the grid and weights adopted from the stored state (today's path), the standardization
   installed from it, the range from the model's attribute and the anchor and unit from the control's, then
   written to the data object; the reload is today's bit for bit. (b) Derive afresh: a grid set by
   `setCutPoints` is lost and the stored trees are moved. (c) Refuse, with a refit message. Recommended: (a),
   about 50 lines of R and the adopting path kept rather than deleted; no released version wrote such
   objects, so (b) is defensible and saves those lines.
7. Smaller calls, each recommended as written in the records and the rule above: the weights as the entry's
   `mass` attribute rather than a fifth slot; the standardization keyed by forest and covariate; the slot name
   `forest.calibration`; a misfit record refused, not dropped; a record carried to a sampler of another
   family used as it is (dec-B232's cost); and no `dbartsData()` argument for any record in 1.0-0, the slots
   documented and settable, an argument being an addition later. A host such as stan4bart could then hand a
   range as `d@response.range <- pair` instead of its two `setResponse` calls.

## Not here

Weights on `setCutPoints` (dec-B450, TODO `cut-point-weights`); the rows a leaf's center and scale are read
from, and `xbart`'s folds (response-scale-rows slice D, dec-B304); kept gp draws that store their own kernel
(dec-B237's later option); a multiplier's basis scale, a fifth record of the same kind, which forest-multiplier
push 2 adds to the data object itself and whose read-back folds into this one's whichever lands second; the
column's level tables (factor-levels-across-setdata.md).
