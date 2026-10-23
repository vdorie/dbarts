# Derived data on the data object: the cut grid, the leaf standardization and the response range

Status: DECIDED 2026-10-17; slice 1 of the plan LANDED 2026-10-23 (954235cb; calls dec-A222); revised from its blind critique and then to the maintainer's rulings.
Plan: [state-frame-prior.md](../plans/state-frame-prior.md). Rulings: dec-B231, dec-B232 and dec-B233, with
dec-B200 as dec-B418 revises it, dec-B230, dec-B237, dec-B304, dec-B305, dec-B310, dec-B329, dec-B330, dec-B449,
dec-B450, and for this design dec-B457 to dec-B461 ([Rulings](#rulings)) in docs/decisions.md; dec-A159 (the
forest calibration anchor goes to the data object beside the range) was confirmed under dec-B378. The calls
no ruling names are [Smaller calls](#smaller-calls), each with its alternative.

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
| `standardization` | `NULL` as a whole where no forest's leaves take a covariate, else a list, one entry per forest, `NULL` where the forest's leaves take none: `center` and `scale` (`NA` for a column without spread), each named by the leaf covariates (by their positions in `columns` where one of those names is missing, empty or another column's too), and `lengthscale` only on a gp leaf whose lengthscale the model does not supply |
| `response.range` | `c(min, max)`, what [`bartcore_anchor`](../../src/R_interface_bartcore.cpp) reports today, under a family that reads one from the response: gaussian, Student-t, aft, and nbinom, whose window of width 1 about its count center is that center's record. `NULL` under probit, logistic, ordinal and multinomial, which today report `c(0, 0)` and read none. A constant response is `c(c, c)`. It carries no family: one carried to a sampler of another family is used as it is ([Smaller calls](#smaller-calls)) |
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
  none, a record `newData` carries is used, and one it lacks is derived. Presence alone decides: a record
  the sampler wrote into its own data object counts as carried, and nothing is compared (dec-B458). Where `setData` puts the columns
  in the sampler's order (dec-A215) the grid's entries move with their columns. A linear leaf's slopes, live
  and kept, are converted wherever the standardization changes (dec-B329), and kept draws are rewritten
  wherever the range changes, as today ([`Chain::applyNewData`](../../src/bartcore/chain.hpp)). Each argument
  defaults to `NULL`, the rule; `updateCutPoints` takes "none", "position" or "value", the other two `TRUE`
  (derive) or `FALSE` (keep the sampler's own, never the carried one), all after `updateState` (dec-B460).
- `setPredictor(x, column, forceUpdate, updateCutPoints = "none", updateState = NULL, updateStandardization =
  FALSE)`. `TRUE` derives the center, scale and heuristic lengthscale of the leaf covariates over the new
  values and converts linear slopes, live and kept; a gp sampler holding kept draws refuses `TRUE` (dec-B237).
  dec-B233 has it default to keeping them.
- `setCutPoints`, `setModel`, and `updateScale = TRUE` on `setResponse`, `setOffset`, `setWeights` and
  `setActiveRows`: unchanged, except that the record they write or carry over is the data object's.
- A live `setState`. Writes no record. The grid: dec-B449 as landed. The standardization: compared forest by
  forest. Equal, the install is today's, bit for bit. Different, the sampler keeps its own and a linear
  leaf's slopes, live and kept, go in as stored (dec-B457); the next draws move the live ones, and kept
  linear draws predict under the sampler's standardization until replaced. A gp leaf's live fits go in as stored, and kept gp draws made under another
  center, scale or heuristic lengthscale are refused (dec-B237), beside kept draws under another supplied
  lengthscale ([`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp)). The range: as stored
  (dec-B418).
- `copy()`, a reload's first use, and `setState` on a sampler whose engine is gone. The engine is re-created
  from `(control, model, data)` by creation's rule, so it takes every record from the data object and none
  from the state. `copy()` of a live sampler, one whose engine exists in this R session
  (`bartcore_isValidPointer`), then installs its engine's current chain; a sampler that is not live installs
  its stored state by the live rule, forced as these always are (dec-B305), the stored trees moved onto the
  data object's grid by value, and is refused without one, as today (dec-B459). A dead-pointer `setState`
  with no saved state takes the record's grid too, where today it adopts the handed state's. The field
  refresh after an install onto another grid (retired: [`heldAfterAdoption`](../../R/dbarts.R), gone with
  adoption; the field is then the engine's state) stays.
- A warm start. The recipient keeps its own values, and a linear donor's slopes go in as stored (dec-B457),
  where today they are converted from the standardization the donor's state records
  (retired: [`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp), gone in slice 2). The seeded fit then
  starts away from the donor's function where the two standardizations differ, and the next draws move it.
- Saved objects. A 0.9-34 sampler is refused at its first use ([`refuseLegacyState`](../../R/tombstones.R)), as
  today; a 0.9-34 data object has no slots and derives, as 0.9-34 did. An object saved by an earlier bartcore
  build has no record and derives afresh, as any record-less data object does; no migration (dec-B461:
  compatibility is owed to shipped versions only).
- The flat C API has no record. A host that re-anchors through it changes the engine and leaves the R record
  stale, as it leaves the model's record today ([state-not-model.md](state-not-model.md)).

dec-B230 is untouched: a data object, like a state, is not checked against the data it came from.

A gp sampler holding kept draws is then refused by `setData` only where the call would change the
standardization or the range in force, and served where it keeps both; the grid never refuses, kept draws
holding their thresholds by value. Where the call derives either, dec-B237 and dec-B330 refuse it already;
a carried record other than the sampler's is refused too ([Smaller calls](#smaller-calls)).

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
| `copy()` of a live sampler whose stored state is older than `setCutPoints(c(0.25, 0.5, 0.75), 1)` (`updateState = FALSE`) | the state's grid and chain | the state's 100 points and its 17 splits, where the sampler has 3 and 16; the same for a multinomial sampler | the live chain on the 3 points: the sampler as it stands |
| a reload of the same object | the state's grid | the state's 100 points, 17 splits | the 3 points, the stored splits moved onto them, 16 |
| the same with `updateState = TRUE` | 3 points | 3 points, 16 splits | the same |
| dead-pointer `setState(st)` with no saved state | the state's grid | the handed state's grid | the record's grid |
| forced `setState` across a standardization (linear leaf, a column replaced by 3 a + 5 through `setData`) | no leaf models | the state's (0.518, 0.269) installed over data whose own are (6.55, 0.807): the slope's prior three times as wide in the column's units; kept draws predict as stored | the sampler's (6.55, 0.807) kept, slopes as stored; kept linear draws predict under it until replaced |
| warm start from a linear donor under another standardization | no leaf models | slopes converted: the seeded fit is the donor's | slopes as stored: the seeded fit is not the donor's until the draws move it (leaf-conversions.md measured up to 6.30 between the two) |
| a sampler saved by an earlier development build | refused (0.9-34's own) | its state's grid adopted | derived afresh |
| gp sampler holding kept draws, `setData(s$data)` | no leaf models | refused | served |
| `attr(s$model, "response.range")`, the control's `anchor` and `unit` | absent | read by tests | gone; the data object's slots |
| `setData` of a gp sampler holding kept draws, the data carrying another standardization or range | no leaf models | refused | refused, with dec-B237's remedy |

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
  slopes go in as stored, and the warm start's conversion,
  retired: [`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp), goes with it (gone in slice 2). `adoptGrid` in
  [`Sampler::setState`](../../src/bartcore/sampler.hpp) goes (dec-B461).
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
- `setState` installing the standardization, as the tip does, and converting slopes at every install: not
  taken (dec-B457).

## Rulings

Put to the maintainer one at a time on 2026-10-17, each with its background, 0.9-34's behavior and its
options; the entries hold the alternatives and their costs.

- A restored state's slopes (dec-B457). A linear leaf's stored slope is a number per unit of the scaled
  covariate. dec-B329 converts slopes at `setData`, which moves the data under a live chain; dec-B418 has
  every install put parameters in "as stored ... and never convert them", and the code converted at a warm
  start. The maintainer: "install as stored", at every install, the warm start included. The sampler keeps
  its standardization; kept linear draws predict under it until replaced, which the help says.
- The sampler's own data handed back (dec-B458). Presence of a record alone decides; nothing is compared.
  The edited `s$data` idiom holds where 0.9-34 derived; NEWS and the help name it. "Hold is fine."
- A copy whose stored state is older than a grid change (dec-B459). `copy()` of a live sampler copies its
  engine's current chain; a sampler that is not live is rebuilt from its stored state on the data object's
  grid, its trees moved by value; the help says to `storeState()` before saving after a grid change. Against
  this design's first recommendation, which installed the stored state in both cases.
- setData's argument values (dec-B460). `NULL` by default, then the setters' existing vocabularies.
- Objects saved by a development build (dec-B461). No migration; such an object derives afresh. A standing
  rule: compatibility is owed to shipped versions only.

## Smaller calls

Two are the orchestrator's readings of rulings already made; the rest are agent-made, each with its
alternative, for the maintainer's later mark.

- A gp sampler holding kept draws, handed by `setData` a carried standardization or range other than its
  own, is refused with dec-B237's remedy (`updateStandardization = FALSE` and `updateScale = FALSE`, or a new
  sampler): it is the change dec-B237 and dec-B330 refuse, reached without an argument. Read as settled by
  dec-B237. The alternative, ignoring the carried record silently, drops what the caller supplied.
- A range carried to a sampler of another family is used as it is: dec-B232 accepted exactly that cost, a
  range "derived under the settings of the sampler that made it" used under others. No family is recorded.
  A family that reads no range (probit, logistic, ordinal, multinomial) does not read it, and its sampler
  writes `NULL`. The alternative, a family tag under which another family derives, is not taken.
- The copy's own `state` field after `copy()` of a live sampler holds the chain it was made from, so a
  reload of the copy rebuilds what the copy was. Alternative: the source's field, as today, which can be older
  than the copy's chain.
- The weights are the entry's `mass` attribute. Alternative: a fifth slot parallel to `cutPoints`.
- The standardization is keyed by forest and by covariate name, and by column position, in a `columns`
  element, where a leaf covariate's name is missing, empty or another column's too, so the two kinds of key
  never meet and a name key matches only the one column it names; a leaf whose name repeats elsewhere in the
  data derives under a named record. Alternative: by column index
  throughout, which a reordered `setData` would have to remap.
- The anchor and unit's slot is `forest.calibration`. Alternative: inside `response.range`.
- A record that does not fit the data (another column count or kind, an ordered factor's grid other than its
  levels') is refused, naming the slot. Alternative: dropped and derived, which hides a mistake.
- A standardization describes the covariates a leaf model takes, a model choice and not the data's, so it is
  matched to them by name: an entry for a covariate the leaves do not take is ignored and a covariate they take
  without an entry derived; only a malformed one is refused. Alternative: refused as a misfit, where a change
  of leaf model hides no mistake.
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
(dec-B237's later option); a multiplier's basis scale, a fifth record of the same kind, which
forest-multiplier push 2 adds to the data object itself and whose write-back folds into this one's whichever
lands second; the columns' level tables (factor-levels-across-setdata.md).
