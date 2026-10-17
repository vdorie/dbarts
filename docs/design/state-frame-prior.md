# Derived data on the data object: the cut grid, the leaf standardization and the response range

Status: PROPOSED 2026-10-17, no code; revised the same day from its blind critique. Plan:
[state-frame-prior.md](../plans/state-frame-prior.md). Rulings: dec-B231, dec-B232 and dec-B233, with
dec-B200 as dec-B418 revises it, dec-B230, dec-B237, dec-B304, dec-B305, dec-B310, dec-B329, dec-B330, dec-B449
and dec-B450 in docs/decisions.md; dec-A159 (the forest calibration anchor goes to the data object beside the
range) was confirmed under dec-B378. Calls no ruling covers are either [Open questions](#open-questions), each
self-contained for the maintainer, or [Smaller calls](#smaller-calls), each with its alternative. Nothing
below settles an open question; where the rule depends on one, it says so.

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

"Held" means the call keeps the value in force.

| value | derived | `setData` | `setPredictor` | live `setState` | copy, reload | recorded |
|---|---|---|---|---|---|---|
| cut grid and its weights | at creation, every row | derived again | held unless `updateCutPoints` | never written (dec-B449) | taken from the state | the state's `cutPoints` and `cutMass` |
| leaf center, scale; heuristic lengthscale | at creation, every row | derived again; linear slopes converted (dec-B329) | held, no argument | installed from the state | derived, then the state's installed | the state, per forest |
| response range | at creation, rows in | derived again | held | never moved (dec-B418) | the model's record applied | the model's `response.range` attribute |
| forest anchor and unit | at creation | refused | held | not in a state | the control's record | the control's `bartcore.forests` attribute |

Two values already live on `dbartsData` and show the pattern this note extends: `sigma`, the starting estimate
`setData` carries over, and `missing.seen`, written from the engine after each predictor change and read when
a sampler is created, so a copy and a reload see what their source saw ([`dataMissingSeen`](../../R/data.R)).

A gp sampler that holds kept draws refuses `setData` outright
([`bartcore_setData`](../../src/R_interface_bartcore.cpp)), since the call always re-derives the
standardization and the range and a kept gp draw replays only under the kernel and shift it was drawn with
(dec-B237, dec-B330; [leaf-conversions.md](leaf-conversions.md)).

## The records

Four slots on [`dbartsData`](../../R/A_class.R), each `NULL` where absent. Each is read through an accessor
that answers `NULL` for an object saved before the slot existed, as [`dataMissingSeen`](../../R/data.R) is,
and checked by a plain function at the doors that use it, as [`checkBasisRecord`](../../R/data.R) is, never by
the class's validity. Names: the first three are dec-B233's.

| slot | holds |
|---|---|
| `cutPoints` | a list, one entry per predictor column in the object's column order: the strictly increasing points; `NULL` for an unordered factor. A column whose points carry weights has them as the entry's `mass` attribute, the increasing sequence a state's `cutMass` holds for it |
| `standardization` | a list, one entry per forest, `NULL` where the forest's leaves take no covariate: `center` and `scale` (`NA` for a column without spread), each named by the leaf covariates, and `lengthscale` only on a gp leaf whose lengthscale the model does not supply |
| `response.range` | `c(min, max)`, what [`bartcore_anchor`](../../src/R_interface_bartcore.cpp) reports today, under a family that reads one from the response: gaussian, Student-t, aft, and nbinom, whose window of width 1 about its count center is that center's record. `NULL` under probit, logistic, ordinal and multinomial, which today report `c(0, 0)`. A constant response is `c(c, c)` |
| `forest.calibration` | `c(anchor = , unit = )` on a model of several forests with a calibration map, what the control's attribute holds today; `NULL` elsewhere, a multinomial sampler's K forests included |

The record lives in three places, each holding what its layer can:

- R. The slots on `sampler$data` are the record a re-creation reads. The sampler writes them from the engine
  after creation and after every call that can change one. The object a caller handed in is never written:
  two fits made from one `dbartsData()` object each derive their own.
- Bridge. [`parseData`](../../src/R_interface_bartcore.cpp) reads the slots into plain vectors and hands them
  to every factory that builds an engine, the multinomial one and the data handle included; one reader
  returns all four from the engine for R to store. Every R type stops here.
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

- Creation: `dbarts()`, `bart()`, `bart2()`, `dbartsSpec()`, `new("dbartsSampler", control, model, data)`,
  and the sampler `samplePriorPredictive` makes from `sampler$data`. Each record the data object carries is
  used; an absent one is derived. The fitting functions take a `dbartsData` object as it is (`dbartsData()`
  returns one unchanged), so `bart(s$data)` and `dbarts(s$data)` carry the sampler's record; a formula and a
  frame build a fresh object, which derives as now. A record that does not fit the data is refused, naming
  the slot and saying that setting it to `NULL` derives it ([Smaller calls](#smaller-calls)).
- `xbart`. Its folds are views of one store inside the engine, built through the data handle. The handle
  ignores every carried record and derives as today: dec-B304 has each fold take its scale from its own
  training rows, which a whole-data record would override. `xbart(s$data)` is therefore today's `xbart`.
- `setData(newData, updateState = NULL, updateCutPoints = NULL, updateScale = NULL, updateStandardization =
  NULL)`. Each record separately: an argument given is followed (derive, or keep what the sampler holds); with
  none, a record `newData` carries is used, and one it lacks is derived. A record the sampler wrote into its
  own data object counts as carried; whether it should is open question 2. Where `setData` puts the columns
  in the sampler's order (dec-A215) the grid's entries move with their columns. A linear leaf's slopes, live
  and kept, are converted wherever the standardization changes (dec-B329), and kept draws are rewritten
  wherever the range changes, as today ([`Chain::applyNewData`](../../src/bartcore/chain.hpp)). The argument
  values are open question 4.
- `setPredictor(x, column, forceUpdate, updateCutPoints = "none", updateState = NULL, updateStandardization =
  FALSE)`. `TRUE` derives the center, scale and heuristic lengthscale of the leaf covariates over the new
  values and converts linear slopes, live and kept; a gp sampler holding kept draws refuses `TRUE` (dec-B237).
  dec-B233 has it default to keeping them.
- `setCutPoints`, `setModel`, and `updateScale = TRUE` on `setResponse`, `setOffset`, `setWeights` and
  `setActiveRows`: unchanged, except that the record they write or carry over is the data object's.
- A live `setState`. Writes no record. The grid: dec-B449 as landed. The standardization: compared forest by
  forest. Equal, the install is today's, bit for bit. Different, the sampler keeps its own; what happens to
  the slopes is open question 1. A gp leaf's live fits go in as stored, and kept gp draws made under another
  center, scale or heuristic lengthscale are refused (dec-B237), beside kept draws under another supplied
  lengthscale ([`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp)). The range: as stored
  (dec-B418).
- `copy()`, a reload's first use, and `setState` on a sampler whose engine is gone. The engine is re-created
  from `(control, model, data)` by creation's rule, so it takes every record from the data object and none
  from the state, and the stored state then goes in by the live rule, forced as these always are (dec-B305).
  A dead-pointer `setState` with no saved state takes the record's grid too, where today it adopts the handed
  state's. Open question 3 is the copy whose stored state is older than a grid change. The field refresh
  after an install onto another grid ([`heldAfterAdoption`](../../R/dbarts.R)) stays.
- A warm start: unchanged. The recipient keeps its own values; a donor's slopes are converted from the
  standardization its state records ([`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp),
  dec-B329).
- Saved objects. A 0.9-34 sampler is refused at its first use ([`refuseLegacyState`](../../R/tombstones.R)), as
  today; a 0.9-34 data object has no slots and derives, as 0.9-34 did. An object saved by an earlier bartcore
  build is open question 5.
- The flat C API has no record. A host that re-anchors through it changes the engine and leaves the R record
  stale, as it leaves the model's record today ([state-not-model.md](state-not-model.md)).

dec-B230 is untouched: a data object, like a state, is not checked against the data it came from.

A gp sampler holding kept draws is then refused by `setData` only where the call would change the
standardization or the range in force, and served where it keeps both; the grid never refuses, kept draws
holding their thresholds by value. Where the call derives either, dec-B237 and dec-B330 refuse it already;
the one open case is a carried record other than the sampler's, open question 6.

## What changes for users

Measured on the tip and on 0.9-34 with 200 rows, two uniform columns, 20 trees and one chain (the plan names
the runs). "Held" was emulated on the tip with the calls that exist: `setCutPoints` by value after
`setData`, and the engine's anchor entry with the old range.

| door | 0.9-34 | tip | proposed |
|---|---|---|---|
| `setData` of the sampler's own object, edited (`d <- s$data; d@x[, 1] <- 2 * x1; s$setData(d)`) | grid derived again: [0.023, 0.983] to [0.046, 1.966] | the same | held at [0.023, 0.983]; posterior mean fit 0.177 from the truth against 0.081 derived, fits up to 0.52 apart |
| the same with `d@y <- 3 * y + 10` | range derived again | [0.038, 2.294] to [10.11, 16.88] | held; 0.324 from 3 f + 10 against 0.244, up to 0.72 apart |
| `setData` of a fresh `dbartsData()` | derived | derived | derived, unchanged |
| creation from `s$data`: `new()`, `dbarts()`, `bart()`, `dbartsSpec()`, `samplePriorPredictive` | grid and range derived from the data | the same; the range from the model where the door re-creates | the sampler's grid and range; `samplePriorPredictive` on the sampler's grid where it derived one from the data |
| a control taken from a sampler of several forests, handed other data | no anchor | the control's anchor (0.190 where the data give 0.0796) | the data's own, derived |
| `xbart(s$data)` | derived | derived | derived: the handle ignores records |
| copy or reload whose stored state is older than `setCutPoints(c(0.25, 0.5, 0.75), 1)` (`updateState = FALSE`) | the state's grid | the state's 100 points and its 17 splits, where the sampler has 3 and 16; the same for a multinomial sampler | open question 3; as recommended, the 3 points with the stored splits moved onto them, 16 |
| the same with `updateState = TRUE` | 3 points | 3 points, 16 splits | the same |
| dead-pointer `setState(st)` with no saved state | the state's grid | the handed state's grid | the record's grid |
| forced `setState` across a standardization (linear leaf, a column replaced by 3 a + 5 through `setData`) | no leaf models | the state's (0.518, 0.269) installed over data whose own are (6.55, 0.807): the slope's prior three times as wide in the column's units | the sampler's (6.55, 0.807) kept; the slopes by open question 1 |
| gp sampler holding kept draws, `setData(s$data)` | no leaf models | refused | served |
| `attr(s$model, "response.range")`, the control's `anchor` and `unit` | absent | read by tests | gone; the data object's slots |

Nothing moves for a fit that never reaches one of these doors: every creation from a fresh data object
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
  the path `applyAnchor` uses. Absent derives, through today's code. The data handle is not handed any.
- [`Sampler::setData`](../../src/bartcore/sampler.hpp) takes, per value, derive, keep or a supplied value.
  Keep re-quantizes the new rows on the grid in force, regathers the leaf covariates under the constants in
  force for the new row count, and installs the response without moving the transform. Derive is today's
  path, unchanged.
- [`Chain::setState`](../../src/bartcore/chain.hpp) stops calling `restoreCalibration` on an install; the
  slopes as open question 1 has them. `adoptGrid` in [`Sampler::setState`](../../src/bartcore/sampler.hpp)
  goes, unless open question 5 keeps it.
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
that call.

## Alternatives

- Records on the model, where dec-B200 put the range: rejected by dec-B232, a model handed to other data
  carrying numbers derived from one response.
- Records in the state, as the grid is today: a restore would then move a prior (dec-B231, dec-B449).
- A hidden field on the sampler: every re-creation path would need its own carrying, and a host could not hand
  a sampler a range or grid with its data. The data object is already what creation, copy and reload read.
- Holding at every plain `setData`, or requiring the arguments on every call: dec-B231's two rejected options.
- `setState` installing the standardization, as the tip does: open question 1, option (c).

## Open questions

Each is written to be put to the maintainer alone, costliest first. The plan is written to each
recommendation and says what another answer changes.

### 1. A restored state whose slopes were drawn under another standardization

Background. A linear leaf fits a line in each leaf against its covariate after centering and scaling it, so a
stored slope is a number per scaled unit. The center and scale are data and stay with the sampler. When a
state stored under one standardization is installed into a sampler holding another (after `setData` with
new values, or `setPredictor(updateStandardization = TRUE)`), its slopes can either be read as they are, which
changes the line each leaf draws, or rewritten so that each leaf draws the same line. Two rulings speak.
dec-B329, 2026-10-08, on `setData` and the warm start: "Convert by the formula here too, in both
directions." dec-B418, 2026-10-09, on every install: a warm start, `setState`, `copy()` and a reload "install
the chain's parameters (tree structure, leaves, k, sigma and the rest) as stored on the internal scale ...
and never convert them", the maintainer adding "The mental model of the sampler is a linear one through
time, with each sample being a snapshot." The two meet on this case; the code today follows dec-B329 at a
warm start and installs the state's standardization at `setState`. 0.9-34 had no linear leaves.

Options:

- (a) The sampler keeps its own; live and kept slopes are converted so each leaf draws its stored line; the
  install reports clean. About 60 engine lines. A gp leaf's live fits cannot be converted and go in as
  stored, so the "same function" promise holds for linear leaves only.
- (b) As (a), but not clean: unforced it declines, forced it converts. No unforced way to restore a state
  stored before `updateStandardization = TRUE`.
- (c) Install the state's standardization, as the tip does: restoring a state then changes the slope prior,
  which is model, and a state carries data.
- (d) The sampler keeps its own and the slopes go in as stored: bit for bit, no conversion code, clean. The
  live leaves and the kept draws then draw other lines until the next sweeps move the live ones; `predict`
  reads the kept draws under the sampler's standardization, as it reads kept draws across a change of
  response range under dec-B418.

Recommended: (d). dec-B418 is the later ruling and the one about installs, and its intent, a chain as a
sequence of snapshots under one model and data, with a restore across a change legal and of no special
meaning, is what (d) does; it is also how the package already treats the response range, where a re-anchor
rewrites kept draws but an install does not. dec-B329 governs the calls that move the units themselves,
`setData` and the warm start, and stays. The cost of (d): after a restore across a change of standardization,
kept linear draws predict other values until replaced, which the help says, as it says for the range.

### 2. The sampler's own data object handed back to setData

Background. Under the new records, `s$data` carries what the sampler derived. A caller who edits it and hands
it back (`d <- s$data; d@x[, 1] <- newValues; s$setData(d)`) then gets the old cut grid, range and
standardization held, where 0.9-34 derived them again from the new values: on the measured fixture the fit
was twice as far from the truth. dec-B231 gives the default: "re-derive if no values are present ... and if
values are supplied, use them". It does not say whether a record the sampler itself wrote counts as
supplied. Its other words cut both ways: "You don't want to update derived data in an ambiguous case. It
needs to be explicit, because a caller might not understand that setData is changing a different part of
the sampler." Holding is not changing; but the 0.9-34 caller expects the change. A host that moves predictors
or the response every sweep wants the hold.

Options:

- (a) A record the sampler wrote counts as supplied: held. A NEWS line, and the help names the three
  arguments and a `NULL` slot as the ways to derive.
- (b) As (a), and a warning where a held grid leaves a column's new values past its last point or a held range
  leaves the new response outside it. About 40 lines and a warning class; it fires inside a host's loop where
  the hold is meant.
- (c) A record identical to the sampler's own counts as absent: the edited idiom derives as in 0.9-34, but a
  host handing back `s$data` with a new response re-derives every sweep, which dec-B296 rules out, and the
  record is "supplied" only when the caller changed it.

Recommended: (a). It is dec-B231's rule read plainly, it serves the host case the ruling was made for, and the
change from 0.9-34 is one idiom the help can name.

### 3. A copy or reload whose stored state is older than a grid change

Background. A copy and a reload are rebuilt from the sampler's control, model and data, then given its stored
state. With `updateState = FALSE` the stored state can predate a `setCutPoints` call. 0.9-34 and the tip both
rebuild on the stored state's grid. Under the records the rebuilt engine has the data object's grid, which is
the live sampler's.

Options:

- (a) The record's grid; the stored trees are moved onto it by value, as a forced `setState` moves them. The
  copy has the sampler's grid and the stored chain's trees moved: neither the live sampler nor the stored
  state exactly (measured: 17 stored splits become 16 on 3 points). The help says to `storeState()` before
  copying or saving after a grid change.
- (b) As (a), and `copy()` reads the live state where the source's engine is alive, so a copy is the sampler as
  it stands. A sampler kept with `updateState = FALSE` no longer copies its stored state; a reload has no
  live engine and stays (a). About 20 lines.
- (c) Today's: the stored state's grid. The copy is the stored state exactly, on a grid the sampler no longer
  holds; a state installs a grid on these two paths, against dec-B231, and the adopting code stays.

Recommended: (a). It keeps "a state never installs a grid" on every path; the case is reached only with
`updateState = FALSE` and a grid change since the last store, which the help can prevent.

### 4. setData's three argument values

Background. dec-B233 names `updateCutPoints`, `updateScale` and `updateStandardization` on `setData`, not their
values. `setPredictor`'s `updateCutPoints` already takes "none", "position" and "value"; the other setters'
`updateScale` is `TRUE` or `FALSE`. 0.9-34's `setData` took `(newData, updateState)`.

Options:

- (a) Each `NULL` by default, meaning the rule; `updateCutPoints` takes `setPredictor`'s words, "none" keeping
  the sampler's grid and "position" or "value" deriving and placing the splits (today's `setData` is
  "value"); the other two take `TRUE` (derive) or `FALSE` (keep the sampler's). An explicit keep is the
  sampler's value, never the carried one. They follow `updateState`, so `setData(d, TRUE)` keeps its 0.9-34
  meaning.
- (b) Logical throughout, `updateCutPoints = TRUE` placing by value: loses "position".
- (c) A fourth value naming the carried record: a word for what `NULL` already does.

Recommended: (a): one vocabulary for the grid across the three setters, and `NULL` as base R's absent.

### 5. Samplers saved by a bartcore build before this one

Background. A sampler saved before the records existed has a data object without them. No release wrote such
objects: 1.0-0 is not released and 0.9-34 samplers are refused already. Of the consumers, stan4bart saves its
own control, model and data with the states and restores by making a new sampler from them, a data object
that never carried a record; bartCause and treatSens save no sampler (read with `git grep`).

Options:

- (a) Migrate at the first rebuild: grid and weights adopted from the stored state, standardization installed
  from it, range from the model's old attribute, anchor from the control's; the reload is today's bit for
  bit. About 50 lines of R, and the adopting paths stay in 1.0-0 for development objects only.
- (b) Derive afresh, as a data object without records does anywhere: no code. Such a reload keeps every value
  that was never changed after creation; one with a grid set by `setCutPoints` or a re-anchored range comes
  back on derived values, its trees moved.
- (c) Refuse, with a refit message: about 10 lines.

Recommended: (b). The objects exist only on development machines, no consumer needs them, and (b) needs
nothing beyond the rule every other record-less object follows; (c) would also refuse the many such objects
(b) reloads exactly.

### 6. A gp sampler with kept draws handed a record other than its own

Background. A kept gp draw replays only under the standardization and response range it was drawn with.
dec-B237 refuses a change of standardization on such a sampler and dec-B330 refuses a re-derived range
"outright, whatever the new values", so a `setData` that derives either is refused, and one that keeps both
is served. The case no ruling names: a data object carrying a standardization or range other than the
sampler's, used under the default rule.

Options:

- (a) Refused when any number of the carried record differs from the sampler's, with the remedy (`FALSE` for
  the argument, or a new sampler).
- (b) The carried record ignored and the sampler's kept, silently.

Recommended: (a): it is the same change dec-B237 and dec-B330 refuse, reached without an argument; (b) drops
what the caller supplied without a word, against dec-B231's "use them".

### 7. A range carried to a sampler of another family

Background. A range is read from the response under the family that made it: an aft fit's from the observed
times, a count fit's window about its center. A data object made for one family can be handed to a sampler of
another (the test suite hands a gaussian sampler's data to an aft one). dec-B232 named the cost that a range
"derived under the settings of the sampler that made it" is used under others.

Options:

- (a) Used as it is, the cost dec-B232 named.
- (b) The record carries its family; under another family it counts as absent and is derived. A string on the
  record, about 15 lines.

Recommended: (b). The cost dec-B232 named was a carried value a fresh derivation would not give; across
families the value does not mean the same thing, and (b) costs little.

## Smaller calls

Agent-made, each with its alternative, for the maintainer's later mark:

- The weights are the entry's `mass` attribute. Alternative: a fifth slot parallel to `cutPoints`.
- The standardization is keyed by forest and by covariate name. Alternative: by column index, which a
  reordered `setData` would have to remap.
- The anchor and unit's slot is `forest.calibration`. Alternative: inside `response.range`.
- A record that does not fit the data (another column count or kind, an ordered factor's grid other than its
  levels', a standardization for columns the leaf does not take) is refused, naming the slot. Alternative:
  dropped and derived, which hides a mistake.
- A carried grid out of order is sorted and one that repeats a point is refused, as dec-A170 as the
  maintainer revised it has `setCutPoints` do. Alternative: refuse any unsorted grid.
- A carried grid is used whatever `n.cuts` and the quantile rule say; those two say how to derive. Alternative:
  refuse a grid whose length differs from `n.cuts`.
- Where `setData` uses a carried grid other than the sampler's, the splits go to the nearest point by value,
  as a derivation at `setData` places them today. Alternative: by position.
- No `dbartsData()` argument for a record in 1.0-0; the slots are documented and settable, an argument being
  an addition later. A host such as stan4bart could then hand a range as `d@response.range <- pair` instead
  of its two `setResponse` calls. Alternative: arguments now.

## Not here

Weights on `setCutPoints` (dec-B450, TODO `cut-point-weights`); the rows a leaf's center and scale are read
from, and `xbart`'s folds (response-scale-rows slice D, dec-B304); kept gp draws that store their own kernel
(dec-B237's later option); the conversion of a warm start's slopes, which dec-B329 and dec-B418 also both
reach and which this item does not change; a multiplier's basis scale, a fifth record of the same kind, which
forest-multiplier push 2 adds to the data object itself and whose write-back folds into this one's whichever
lands second; the columns' level tables (factor-levels-across-setdata.md).
