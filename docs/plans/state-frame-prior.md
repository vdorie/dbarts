# state-frame-prior: the derived data's records on dbartsData, and every door reading them

Status: PLANNED 2026-10-17; revised from its blind critique (scratch/sfdc-crit.md, sound after revisions;
[Revision from critique](#revision-from-critique)) and then to the maintainer's rulings dec-B457 to dec-B461
([Revision to the rulings](#revision-to-the-rulings)). No call open. Design:
[state-frame-prior.md](../design/state-frame-prior.md). Rests on dec-B231 to dec-B233, dec-B237, dec-B304,
dec-B305, dec-B310, dec-B329, dec-B330, dec-B418, dec-B449, dec-B450, dec-B457 to dec-B461 and dec-A159
(confirmed under dec-B378).

agent: per slice, one opus implementer and one opus reviewer who runs that slice's mutants. Tier "Changes
draws" ([Process by risk](README.md#process-by-risk)) for all three: each edits the engine's creation,
`setData` or install path.
rng, per slice, against equivalence-32418b6e, bcf-equivalence-47838362 (current since forest-multiplier
push 1 landed; bcf-equivalence-1b7d730c before it) and multinomial-equivalence-80b1c8d4. No recorded scenario reaches a changed
sequence: the `setdata` and `hetsetdata` scenarios hand a fresh `dbartsData()`, the `setPredictor` scenarios
name no new argument, none copies, reloads, calls `setState` or `samplePriorPredictive`, or hands one
sampler's control or data to another, and none makes a warm start. Expected 55, 15 and 11 identical in every slice, no re-record; a moved
stream is a stop.
- 1, grid, range and anchor. NEUTRAL, bit for bit, for every creation from a data object that carries no
  record, every live call, `xbart` (its handle ignores records), and every copy, reload and dead-pointer
  `setState` whose stored state records the grid and chain in force. POSTERIOR-CHANGING on seven sequences,
  each a changed tree prior or leaf prior or a changed chain: `copy()` of a live sampler whose stored state is
  older than its chain (the engine's current chain and grid, dec-B459); a reload of such an object (the
  record's grid, the stored trees moved by value); a dead-pointer `setState` with no saved state (the record's grid, not
  the handed state's); a creation from a data object carrying a record other than a fresh derivation's
  (`new()`, `dbarts()`, `bart()`, `bart2()` or `dbartsSpec()` handed `s$data`); `samplePriorPredictive` (the
  sampler's grid, where it derived one from `data@x`); a control of several forests handed other data (the
  data's anchor, not the control's); a reload of an object saved before this slice whose grid or range was
  changed after creation (derived, dec-B461).
- 2, standardization. NEUTRAL for every install, warm start and re-creation under the standardization in
  force. POSTERIOR-CHANGING for a `setState`, copy or reload across a standardization: the sampler keeps its
  own where the tip installs the state's, and the slopes go in as stored (dec-B457). SHIFTING for a warm start
  from a linear donor under another standardization: the slopes go in as stored where the tip converts
  them, so the seeded fit and the draws after it differ and the recipient's posterior does not; such a donor
  that the tip declined unforced for a column without spread (`keptAcrossValues`) is clean. Refused where the
  tip installs: kept gp draws under another standardization or heuristic lengthscale.
- 3, setData and setPredictor. NEUTRAL for `setData` of a data object carrying no record and `setPredictor`
  without `updateStandardization`. POSTERIOR-CHANGING for `setData` of a carried record (held where the tip
  derives) and for the new arguments. Accepted where refused: a gp sampler's `setData` that keeps the
  standardization and range. Still refused: a carried standardization or range other than the gp sampler's.
window: before the merge to main. Engine slices stay serial (Order and interactions).
budget: 1 about 1860 lines, 2 about 1100, 3 about 1650; [Budget and stops](#budget-and-stops).

## Goal

Each derived value has one record, on the data object: `cutPoints` (with its weights), `standardization`,
`response.range` and `forest.calibration`. Creation and `setData` use a carried record, derive an absent one
and follow an explicit argument; `setPredictor` takes `updateStandardization`; every install, the warm start included, keeps the
sampler's standardization and puts slopes in as stored; a copy of a live sampler copies its current chain,
and every re-creation takes every value from the record and none from the state; the
gp `setData` refusal covers only calls that would change a kept draw; the help carries the table of what puts
each value back. Closes TODO `state-frame-prior` and unblocks response-scale-rows slice D.

## Context

- The rule and its doors: [The rule at every door](../design/state-frame-prior.md#the-rule-at-every-door).
  What moves: [What changes for users](../design/state-frame-prior.md#what-changes-for-users), measured by
  scratch/sfd/p1.R (the edited data object), p2.R (a state across a standardization, the gp refusal) and
  p3.R (a copy with a stale state), on the tip's shipped build and on 0.9-34; the critique's q1.R to q9.R in
  scratch/sfdc.
- Today's homes: the grid adopted from a state at re-creation ([`installStateOnto`](../../R/dbarts.R),
  `adoptGrid` in [`Sampler::setState`](../../src/bartcore/sampler.hpp)); the standardization installed from a
  state ([`Chain::setState`](../../src/bartcore/chain.hpp), `restoreCalibration` on
  [`LinearGaussianLeaf`](../../src/bartcore/model.hpp) and [`GPGaussianLeaf`](../../src/bartcore/model.hpp))
  and derived at `setData` ([`Chain::applyNewData`](../../src/bartcore/chain.hpp)); the range on the model
  (retired: [`recordAnchor`](../../R/dbarts.R), gone in slice 1, and [`applyAnchor`](../../R/dbarts.R), copied by `setModel`); the anchor
  and unit on the control (the sampler's `initialize` in [R/dbarts.R](../../R/dbarts.R), read by
  [`applyForestAttributes`](../../src/R_interface_bartcore.cpp)); [`samplePriorPredictive`](../../R/dbarts.R)
  creates from `sampler$data` and applies the model's range.
- Factories that build an engine: [`bartcore_create`](../../src/R_interface_bartcore.cpp),
  [`createMultinomialDataHolder`](../../src/R_interface_bartcore.cpp), and the
  data handle, [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_createFromHandle`](../../src/R_interface_bartcore.cpp), which `xbart` uses.
- A record pattern to follow: `missing.seen`, read through [`dataMissingSeen`](../../R/data.R), written in
  [`bartcoreSamplerSetData`](../../R/bartcore.R) from [`bartcore_getMissingSeen`](../../src/R_interface_bartcore.cpp);
  checks by a plain function, [`checkBasisRecord`](../../R/data.R).
- The gp refusals: [`holdsSavedGPDraws`](../../src/R_interface_bartcore.cpp),
  [`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp), the outright one in
  [`bartcore_setData`](../../src/R_interface_bartcore.cpp), and
  [`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp).
- Tests that read what moves (counted on 6fe5237a): `attr(model, "response.range")` on 7 lines in 5 files and
  through 24 calls of the `recordedRange` helpers of
  [test-response-scale-rows.R](../../inst/tinytest/test-response-scale-rows.R) and
  [test-update-scale-rows.R](../../inst/tinytest/test-update-scale-rows.R); the control's unit through 14
  calls of `unitOf` in [test-forest-sd-unit.R](../../inst/tinytest/test-forest-sd-unit.R). Of the suite's
  126 `setData` calls in 37 files about 20 hand back a sampler's data object or an edit of it, most to pin a
  refusal that comes first; five in test-response-scale-rows.R expect the derivation and its fallback
  warning, and ["donor$setData(donor$data)"](../../inst/tinytest/test-sampler-boundary-refusals.R), after
  a doubled column, will hold the grid and may change the `installTrees` refusal it leads to.

## Decision

Ruled, 2026-10-17: install as stored at every install, the warm start included (dec-B457); a carried record
holds at `setData` by presence alone (dec-B458); `copy()` of a live sampler copies its current chain, one not
live rebuilds from its stored state on the data object's grid (dec-B459); `NULL` defaults and the existing
vocabularies (dec-B460); no migration for development objects (dec-B461). The orchestrator's readings: a gp
sampler's carried record other than its own is refused with dec-B237's remedy; a range carried across
families is used as it is (dec-B232). The design's [Smaller calls](../design/state-frame-prior.md#smaller-calls)
are agent-made and the steps follow them.

## Constraints

- A creation from a data object without a record runs the code that runs today; so does `setData` of one.
- The range goes into a re-created engine by the call `applyAnchor` makes today, after construction, so a
  copy and a reload stay bit for bit; construction takes supplied values for the grid and the
  standardization only.
- The data handle is handed no record (dec-B304's folds).
- The engine names no R type and no door; [dbarts.h](../../inst/include/dbarts/dbarts.h) and the C API hash do
  not move. [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) stays at 1; a state keeps every block.
- No message or warning beyond the refusals named here; `setState` keeps dec-B305's two forms and dec-B310's
  return values.
- Tests that compare to the base build do it with literals recorded from the base build and written into
  the test, or with a two-library compare in scratch recorded in the landing note, as response-scale-rows-c
  recorded its count case; no subprocess in tinytest.
- Out of scope: the design's [Not here](../design/state-frame-prior.md#not-here).

## Slice 1: the grid, range and anchor records, and every re-creation reading them

### Steps

1.1 Engine (about 170 lines, [data.hpp](../../src/bartcore/data.hpp), [sampler.hpp](../../src/bartcore/sampler.hpp),
    [facade.hpp](../../src/bartcore/facade.hpp)). Optional supplied grids and weights at construction,
    installed before the columns are quantized, through the standard and the multinomial factories. One
    reader of the grids, weights, range, anchor and unit, through a facade virtual. `adoptGrid` and its path
    in [`Sampler::setState`](../../src/bartcore/sampler.hpp) go (dec-B461). `--preclean`.
1.2 Bridge (about 160 lines). [`parseData`](../../src/R_interface_bartcore.cpp) reads `cutPoints`,
    `response.range` and `forest.calibration`, each guarded for an object without it; the two creating
    factories take them; [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp) skips them. The
    anchor and unit come from `data@forest.calibration`. The reader entry. A range record carries no family
    and is used under any family that reads one (dec-B232).
1.3 R (about 280 lines). The three slots on [`dbartsData`](../../R/A_class.R), accessors in
    [R/data.R](../../R/data.R), `checkDerivedRecord` at the creation doors and `setData` (a grid out of order
    sorted, a repeated point refused, a misfit refused naming the slot). The sampler writes them after
    creation and after `setCutPoints`, `setPredictor`, `setData` and each `updateScale = TRUE`;
    retired: [`recordAnchor`](../../R/dbarts.R) (gone, as built) writes `data@response.range` (`NULL` under probit, logistic, ordinal and
    multinomial) and `setModel` no longer carries the attribute; the anchor and unit leave the control and
    its `initialize` block. `copy`, [`getPointer`](../../R/dbarts.R), [`recreatePointer`](../../R/dbarts.R),
    the dead-pointer `setState` and [`samplePriorPredictive`](../../R/dbarts.R) re-create from the record.
    `copy()` of a live sampler (`C_dbarts_bartcore_isValidPointer`) installs its engine's current chain, read
    without touching the source's field, and the copy's field holds that chain (about 20 lines, dec-B459); a
    sampler that is not live installs its stored state, forced, onto the record's grid.
    Docstrings of `copy`, `setState` and `storeState`.
1.4 Existing tests: the range and unit reads of Context; the adoption pins of
    [test-state-grid-install.R](../../inst/tinytest/test-state-grid-install.R) ("re-creation adopts", "reload,
    first call") rewritten to the record.
1.5 Tests (about 1100 lines with 1.4): tinytest `test-derived-records.R`; tests/cpp in `test_state.cpp`.
1.6 Help, NEWS, notes (about 150 lines): man/dbartsData.Rd, the sampler's Saving section and `copy`; NEWS
    (creation from `s$data`); one sentence each in [state-not-model.md](../design/state-not-model.md),
    [install-surface.md](../design/install-surface.md), [cut-grid.md](../design/cut-grid.md) and
    [architecture.md](../architecture.md).

### Tests

| Test | Pins | Mutant |
|---|---|---|
| written at creation | each door: the records equal the engine's reader; the caller's object carries none; probit, logistic, ordinal, multinomial: range `NULL`; multinomial: calibration `NULL` | the caller's object written; `c(0, 0)` recorded |
| fresh object neutral | 16 fixtures (each leaf kind, weights, a mask, a factor, a sparse store, several forests, multinomial): readers and the next 20 draws equal literals recorded from the base build | the supplied path taken with derived values |
| copy and reload bitwise | after `storeState()`, the same fixtures: state, readers and next 20 draws equal the recorded literals | the range supplied at construction |
| copy is the live sampler | `updateState = FALSE`, 3 points set, 20 more sweeps: `copy()` has 3 points, the source's current trees, generator and readers, and its next 20 draws equal the source's; its field holds that chain; the source's field untouched; multinomial too | the stored state installed |
| reload moves the stored trees | the same object saved and read: 3 points, the stored trees moved by value (16 splits) | adopting from the state |
| copy of a live sampler with a current field | `updateState = TRUE`: the copy's state and next 20 draws equal the recorded literals | the field and the live read differing |
| dead pointer, no saved state | `setState(st)`: the record's grid | the handed state's grid adopted |
| carried at creation | `new()`, `dbarts()`, `bart()` and `dbartsSpec()` handed `s$data` after `setCutPoints`, `setPredictor(updateCutPoints = "value")` and `setResponse(updateScale = TRUE)`: the sampler's grid, weights and range; `samplePriorPredictive` on the sampler's grid | the record ignored at one door |
| xbart ignores | `xbart(s$data)` after `setResponse(updateScale = TRUE)`: each fold's range its own, equal to `xbart` of a fresh object | the handle handed the record |
| control anchor | a control of several forests handed other data: the data's anchor (q7's 0.0796, not 0.190) | the control's read |
| family as is | a gaussian sampler's data handed to an aft sampler: the range used as carried; to a probit sampler: not read, `NULL` written | derived across families; `c(0, 0)` read under probit |
| weights ride | a 30-value integer column: copy and reload keep its weights | weights dropped |
| misfit and order | wrong column count, an ordered factor's other grid, a range of length 1: refused, nothing created; an unsorted grid sorted; a repeated point refused | the check skipped at one door |
| record-less object | the slots removed (as test-missingness-first-seen.R removes one): derived, as a fresh object | migration |
| constant response | several forests on a constant response (forest-multiplier push 1's case): the anchor recorded and a copy created | the record's anchor of 0 |

## Slice 2: the standardization record, and every install keeping it

### Steps

2.1 Engine (about 120 lines, [model.hpp](../../src/bartcore/model.hpp), [chain.hpp](../../src/bartcore/chain.hpp),
    [data.hpp](../../src/bartcore/data.hpp)). Supplied center, scale and heuristic lengthscale per forest at
    construction, through both creating factories, served as
    [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp) serves a view. The reader returns
    them. [`Chain::setState`](../../src/bartcore/chain.hpp) no longer calls `restoreCalibration`: the
    sampler's standardization stays and the slopes go in as stored (dec-B457); a gp forest whose kept
    draws were made under another center, scale or heuristic lengthscale is refused beside
    [`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp). The warm start's conversion goes:
    [`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp), its call in `Sampler::installForests`
    and its `keptAcrossValues` reason for an install that is not clean (about 70 lines removed).
2.2 Bridge and R (about 110 lines): the slot, its accessor and check, its write-back; the refusal's text
    beside the lengthscale one; [`readWarmStartState`](../../src/R_interface_bartcore.cpp) no longer reads a
    donor's standardization blocks; docstrings of `setState`, `installTrees` and `bart`'s `warm.start`.
2.3 Existing tests: the pins of a restore that expects the state's standardization, and of the warm start's
    conversion: section 4 of [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R) ("a warm
    start between standardizations") and its donor without spread, the conversion pins of
    [test-warm-start-forms.R](../../inst/tinytest/test-warm-start-forms.R), and tests/cpp
    [`testWarmStartStandardization`](../../tests/cpp/test_state.cpp), each rewritten to slopes as stored.
    warm-start-salvage slice D's pins ("linear, other covariates", "linear into constant") already expect no
    restatement and stand.
2.4 Tests (about 750 lines with 2.3).
2.5 Help, NEWS, notes (about 120 lines): the `setState`, `installTrees` and `warm.start` sentences (kept
    linear draws predict under the sampler's standardization after a restore across a change; a seeded fit
    starts away from the donor's where the standardizations differ); NEWS; the standardization row of
    [state-not-model.md](../design/state-not-model.md); [leaf-conversions.md](../design/leaf-conversions.md)'s
    warm-start rule and its table row; [warm-start-salvage.md](../design/warm-start-salvage.md).

### Tests

| Test | Pins | Mutant |
|---|---|---|
| same frame exact | every restore of slice 1's fixtures under the standardization in force: equal to the recorded literals | a branch on a standardization that does not differ |
| keeps its own | p2's sequence: after a forced `setState(st)` the sampler reports (6.55, 0.807); slopes as stored, by the state block; NULL forced, TRUE unforced where the trees fit | the state's installed |
| copy holds it | a copy and a reload after `setPredictor(updateStandardization = TRUE)` emulated through `setData`: the record's standardization | derived at re-creation |
| gp kept refused | kept gp draws under another standardization and under another heuristic lengthscale: error, sampler untouched; live fits alone install | live fits refused too |
| supplied lengthscale | a gp leaf whose lengthscale the model supplies: no `lengthscale` in the record, the model's in force | the record's read over the model's |
| warm start as stored | a linear donor under another standardization, live and from a kept draw: the recipient's slopes are the donor's stored numbers, its standardization its own; TRUE unforced, a column without spread on both sides included | the conversion left in; FALSE kept for `keptAcrossValues` |
| warm start, same frame | a donor under the recipient's standardization: equal to the recorded literals | a branch on equal standardizations |

## Slice 3: setData and setPredictor under the rule, and the manual's table

### Steps

3.1 Engine (about 230 lines). [`Sampler::setData`](../../src/bartcore/sampler.hpp) takes, for the grid, the
    standardization and the range each, derive, keep or a supplied value, and the placement for a derived
    grid. Keep: the new rows quantized on the grid in force; the leaf covariates regathered under the
    constants in force for the new row count; the response installed without moving the transform, the kept
    draws then not rewritten. A supplied grid other than the sampler's places the splits by value. Derive:
    today's path. `setPredictor` and `updatePredictor` take `updateStandardization`: the leaf's derivation over
    the new values, then the conversion of live and kept linear coefficients (dec-B329).
3.2 Bridge (about 90 lines). [`bartcore_setData`](../../src/R_interface_bartcore.cpp) takes the three choices;
    its outright gp refusal becomes: refused where the call derives the standardization or the range (dec-B237,
    dec-B330) or brings a carried one whose numbers differ from the sampler's (the orchestrator's reading of
    dec-B237), with dec-B237's remedy.
    The predictor entries take the flag; a gp sampler holding kept draws refuses `TRUE`.
3.3 R (about 160 lines). `setData(newData, updateState = NULL, updateCutPoints = NULL, updateScale = NULL,
    updateStandardization = NULL)` resolves each choice (dec-B460); a record the sampler wrote counts as
    carried, by presence alone (dec-B458); the carried records are aligned with the columns before the call
    ([`alignDataFactorLevels`](../../R/bartcore.R)'s order); the fallback warning
    ([`scaleRowsFallBack`](../../R/spec.R)) is tested only where the range is derived.
    `setPredictor(..., updateStandardization = FALSE)`.
3.4 Existing tests: the five setData sites of test-response-scale-rows.R (the sequences at its comments "the
    flag survives setData", "setData to such rows", "an nbinom sampler's setData" and the two after) state
    `updateScale = TRUE`; the boundary-refusals donor states `updateCutPoints = "value"`; the gp refusal pins
    of [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R) are split into refused and
    served.
3.5 Tests (about 1000 lines with 3.4).
3.6 Help, NEWS, notes (about 170 lines): the arguments; the design's table,
    [The manual's table](../design/state-frame-prior.md#the-manuals-table), in the sampler's Saving section;
    NEWS: `setData` uses the record a data object carries, with the ways to derive; the design note's Status;
    TODO drops `state-frame-prior` and says response-scale-rows slice D is unblocked.

### Tests

| Test | Pins | Mutant |
|---|---|---|
| fresh object neutral | `setData(dbartsData(...))` on slice 1's fixtures: equal to the recorded literals | the keep path on an absent record |
| carried held | p1's two sequences: grid and range as before the call, by reader; kept draws not rewritten | deriving a carried record |
| arguments win | each argument at each value over a carried and an absent record: 2 x 3 x 3 cases, by reader | `NULL` read as keep |
| keep is the sampler's | `"none"` and `FALSE` with a carried record other than the sampler's: the sampler's values | the carried one used |
| placement | `"position"` and `"value"` on a derived grid; a carried other grid by value | one placement for all |
| new row count | keep with 40 rows added: codes on the old grid, covariates under the old constants, fits finite | covariates left at the old length |
| linear converts | `updateStandardization = TRUE` at both setters: live and kept predictions unchanged to 1e-12 | conversion skipped for kept draws |
| gp narrowed | kept gp draws: `setData(s$data)` served; a fresh object refused; `FALSE, FALSE` served; a carried record differing in one number refused; `setPredictor(updateStandardization = TRUE)` refused; a refusal touches nothing | the refusal by argument alone |
| record written | after each call the slots equal the reader | a stale slot |
| fallback where derived | the five response-scale-rows sequences: warned with `updateScale = TRUE`, silent held | warning on a held range |
| table sequences | each row of the manual's table end to end: the value back, `setState` TRUE, next 5 draws equal a twin that only restored | (pins the help) |

## Verification

Per slice, against the slice's own library ([Landing](README.md#landing)):

- `cd tests/cpp && make && ./test_bartcore`, plain and under ASAN and UBSAN; `R CMD INSTALL --preclean`; the
  full tinytest suite; the R-loaded sanitizer run over the new and changed test files.
- Reference build: the four `test-reproducibility-*.R` files unedited, and the three `--bitwise` compares
  (`--strict-coverage` for the main harness), 55, 15 and 11 "identical draws" lines and no "max |z|" line.
- The exact-gate commands of exact-gates.yaml at `quick`; expected unchanged.
- The base-build literals: recorded once per fixture on the slice's base, written into the test, and the
  two-library compare run in scratch and recorded in the landing note.
- Lint, air, rc-codoc, win-drift and doc-freshness; `R CMD check --as-cran` from a clean tarball; the x86-64
  leg; stan4bart's suite (its restores are first creations from its own data object, which carries no
  record, and forced same-grid installs).

## Budget and stops

| | 1 | 2 | 3 |
|---|---|---|---|
| engine | about 170 | about 120 added, 70 removed | about 230 |
| bridge and R | about 430 | about 110 | about 250 |
| tests | about 1100 | about 750 | about 1000 |
| help, NEWS, notes | about 150 | about 120 | about 170 |
| planned | about 1860 | about 1100 | about 1650 |
| stop, engine lines added | 320 | 220 | 380 |
| stop, all | 2800 | 1700 | 2500 |

Sized from the landed history: state-grid-install planned 250 engine lines and landed 275, and 1475 in all
against 2092 at the implementer's tip (1.4 times), tests 1549 against 1000. Slice 1 names three factories and
a facade reader where that slice had one install path. Stop and report, building nothing further: a bitwise compare prints "max |z|" or a
snapshot moves; a fixture's recorded literal differs where the slice is to be neutral; a test outside those
named in steps 1.4, 2.3 and 3.4 fails because it expected a value derived where the record now holds it
(list it, do not edit it); a slice's stop in the table.

## Order and interactions

The serial engine order: warm-start-salvage slice C (wt/warm-start-c), forest-multiplier push 1
(wt/forest-mult-p1), warm-start-salvage slice D, then slices 1, 2 and 3 here, then response-scale-rows slice D.

- warm-start-salvage slice C: the install paths in sampler.hpp and chain.hpp; slice 1 rebases on it.
- forest-multiplier push 1: [`parseData`](../../src/R_interface_bartcore.cpp) and the anchor on a constant
  response; slice 1 moves that anchor and carries its value ("constant response").
- warm-start-salvage slice D: leaf kinds, [`readWarmStartState`](../../src/R_interface_bartcore.cpp) and
  [`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp). It lands first, as planned, and needs
  no change: where it drops a donor's slopes for other covariates it already makes no restatement, and its
  pins expect none. Slice 2 here then removes the conversion and its flag for every warm start, and rewrites
  the pins slice A of that plan and leaf-conversions left on it (step 2.3).
- forest-multiplier push 2: adds `basis.scale` to the data object and writes it where the anchor is written.
  After slice 1 it writes beside `forest.calibration`; before it, whichever lands second folds the two
  write-backs into one (that plan's "With the state-frame arc's data record"). Either order.
- forest-multiplier push 3: R surface only, the printed block; any time.
- response-scale-rows slice D: after slice 3. It needs slice 2 (the record holds the center and scale through
  a copy, a reload and a state install) and slice 3 (`updateStandardization`).
- The control-migration arc (dec-B239) moves the forests' record to the model; the anchor and unit leave the
  control in slice 1, so that arc moves the rest.

## Not verified

No code is built. The design's "held" numbers are the tip's behavior and emulations of holding with existing
calls (p1 to p3), not runs of the rule. Not read: bartCause and treatSens beyond `git grep` (neither calls
`setData`, `setPredictor`, `setState`, saves a sampler or reads a record). The line estimates are from the code
read and the last state slice's ratios, not an itemization.

## Revision from critique

Item by item, the critique's numbering; the question numbers are the design's before the rulings:

1. Question 1 lacked the install as stored and narrowed dec-B418. Taken: option (d), dec-B329's and dec-B418's
   words and their conflict, recommended (d) from dec-B418's intent; slice 2 written to it.
2. Doors missing from the draw classes. Taken: `bart()`, `samplePriorPredictive`, a control's anchor handed
   other data, `xbart` (its handle ignores records), the dead-pointer `setState` without a saved state, in the
   design's table and the rng block.
3. Question 4 (now 6) mostly ruled. Taken: shrunk to the carried record, citing dec-B237 and dec-B330.
4. A range under families that have none. Taken: `NULL` there; nbinom's window kept as its center's record;
   the family tag is question 7.
5. Factories. Taken: the multinomial factory and the handle in steps 1.1 and 1.2, a multinomial fixture.
6. Today's copy. Taken: option (c) of question 3, with its cost.
7. Migration. Taken: argued from release scope and the consumers; recommended (b), derive.
8. Uncovered calls. Taken: question 2 states that a record the sampler wrote counts as supplied, with both of
   dec-B231's quotes; the placement, `n.cuts`, and the sorted grid (following dec-A170 as revised) are
   smaller calls with alternatives.
9. Order. Taken: warm-start-salvage slice D and forest-multiplier pushes 2 and 3 placed.
10. Base-build compares. Taken: recorded literals or a two-library scratch compare (Constraints).
11. Inventory. Taken: recounted (Context); five response-scale-rows sites and the boundary-refusals donor in
    step 3.4; `setModel` in step 1.3.
12. Slice 1 too large. Taken: split into grid, range and anchor (1) and standardization (2); re-estimated.
13. to 16. Taken: lengthscale absent rather than `NA`; multinomial's calibration `NULL`; dec-B305 and dec-B310
    in place of dec-B234; the probes, counts and commit moved here from the design note.

## Revision to the rulings

dec-B457 to dec-B461, 2026-10-17, settle the design's five open questions; dec-B237 and dec-B232 settle the
other two as the orchestrator reads them. What changed here:

- dec-B457: slice 2 installs slopes as stored and also removes the warm start's conversion
  (`convertDonorStandardization`, the donor's blocks in `readWarmStartState`, `keptAcrossValues`); its draw
  class gains the warm start (shifting); its tests and help cover the warm start.
- dec-B458: step 3.3 counts presence alone.
- dec-B459: step 1.3 copies a live sampler's current chain (about 20 lines, against the design's first
  recommendation); slice 1's draw class and tests split the live copy from the reload.
- dec-B460: step 3.3's vocabulary.
- dec-B461: no migration; `adoptGrid` goes in slice 1 and the install's `restoreCalibration` in slice 2.
- dec-B232 read: no family tag, slice 1's test "family as is". dec-B237 read: step 3.2's refusal of a carried
  record other than a gp sampler's own.



## Slice 1 as built

Built on bartcore 46e1a644, with warm-start-salvage slice C landed; the plan's steps held against that tip,
warm-start-salvage slice D and forest-multiplier push 2 not yet landed. Commits, as rebased onto c1d1473f
(below): 48e4428f (the slice), 7b5dd467 (setData hands no carried record on; below), feb54aa3 (cites of the
retired helpers), 37d9e8cb (the forest-prior exact gate reads the unit from the data object) and this record.

What it does. Three slots on `dbartsData`, `cutPoints` (the weights as each entry's `mass` attribute),
`response.range` and `forest.calibration`, read through accessors in R/data.R that answer `NULL` for an object
saved before them. One engine reader, [`Sampler::derivedRecords`](../../src/bartcore/sampler.hpp) behind a facade
virtual, returns the grids, weights, range (absent where the family carries no units) and a calibration map's
anchor and unit; [`recordDerived`](../../R/data.R) writes the three slots from it after creation, `setCutPoints`,
`setPredictor` with a grid update, `setData` and every `updateScale = TRUE`. Every engine is made by
[`createEngine`](../../R/dbarts.R): the grids go in at construction
([`SuppliedCutGrid`](../../src/bartcore/data.hpp), before the columns are quantized, through the
single-forest, combining and multinomial constructors), the range by the anchor move after construction, the
anchor and unit through the data object in place of the control's attribute, which no longer carries them.
The model no longer carries `response.range`; a model saved with it has it stripped at creation.
[`checkDerivedRecord`](../../R/data.R) holds a record to the data at creation and at `setData`.
`copy()` of a live sampler installs its engine's current chain and keeps it as its field (dec-B459); a reload, a
copy of a sampler whose engine is gone and a dead-pointer `setState` rebuild on the record and move the stored
splits onto it by value. `adoptGrid` is gone from `Sampler::setState`, the facade, the bridge entry (one argument
fewer) and R, and `heldAfterAdoption` with it (dec-B461); the field after a re-creation that moved a state is
the engine's (status 2), which replaces it.

Departures from the steps and the test table:

- parseData reads `cutPoints` and `forest.calibration`, and not `response.range`: the range goes in by the
  anchor move after construction (Constraints), which R makes from the slot, so the bridge never needs it.
- The plan's base-build literals are a two-library compare in scratch, not literals in tinytest: an install
  is not a continuation of its source's stream (a twin handed a sampler's live state draws other numbers than
  the sampler continuing, on the tip as here), and draws recorded on arm64 do not hold on the other CI
  platforms. scratch/sfp1/compare.R records, per fixture, the state, anchor and leaf prior at creation, the
  next 20 draws, a copy's and a reload's state and next 20 draws after `storeState()`, a twin's draws after
  `setState` and `samplePriorPredictive`, on 46e1a644 and on this slice. tinytest pins the doors against each
  other instead: a copy, a reload and a twin handed the stored state hold one state and draw alike.
- "copy is the live sampler" pins the copy's state, generator included, to the source's current one, and its
  draws to a second copy's, not to the source's next draws, for the reason above.
- An ordered factor's carried grid is held to its levels' by reading the engine's grid after creation
  ([`refuseCarriedOrderedGrid`](../../R/data.R)): the level count the grid follows is the engine's.
- Two mutants of the table cannot be written: "the caller's object written" (R's value semantics give the
  sampler its own data object; there is nothing to mutate), and "the range supplied at construction" and
  "the supplied path taken with derived values" (no such path exists; the neutrality they guard is the
  two-library compare's). "migration" has no code to mutate (dec-B461).

Calls made in building, each with its alternative:

- A record the engine refuses to hold for an ordered factor whose data object declares other levels than
  the engine was made with (assigned by hand, saved by an earlier build, or reached with no hand edit, below)
  now makes `copy()` and a reload
  refuse, the state's grid on that column being another column's ("kinds" in
  test-state-grid-install.R). The tip adopted the state and kept the field current through
  `heldAfterAdoption`; with adoption gone (dec-B461) and the live install refusing another grid on an ordered
  factor, there is nothing left to install it with. Alternative: move an ordered factor's splits by value as a
  numeric column's are (engine work and a ruling); not taken, the case first read as reachable only by hand.
  It is not (corrected in the fix round): a sampler made from a matrix of codes declaring the column ordered
  through `attr(x, "varTypes")` and handed fewer codes by `setData` keeps the record of the levels it was made
  with, its own grid, so `copy()` and a reload refuse in R before any install, "the data's 'cutPoints' record
  gives predictor 2, an ordered factor, other cut points than those of its levels; set it to NULL to derive it
  from the data", and with the record set to NULL as that says both refuse "state is not consistent with this
  sampler"; a000d239 made the copy and ran the reload (TODO code-matrix-category-count). The move above does
  not reach that case; a carried level grid setting the level count at creation would. For the maintainer.
- The copy of a sampler whose engine is gone and that has no stored state is made fresh, as on the tip,
  not refused (refused since the fix round, below): dec-B459's "refused, as today" read as the reload's
  refusal, which stands. Alternative: refuse the copy too.
- A NULL entry in `cutPoints` for a numeric column is refused at the R surface (a kind mismatch) and derives
  that column at the bridge, which is only the backstop. Alternative: let it derive at both.
- The record is written after `setPredictor` only where it updated a grid, so a predictor loop pays no read.
  Alternative: after every call.
- A carried range is checked in R (two finite numbers, the second not below the first) and the count
  family's strict order by the anchor entry, whose message now names the data's record.
- `copy()` of a live sampler whose install moved its chain (status 2) keeps the engine's state as its field;
  never reached from a live chain, which is on the record's grid.

- `setData` derives every value in this slice, so it hands the call no carried record: the column alignment
  (dec-A215) can leave a carried grid in another column order than the bridge reads it in, which refused a
  `setData` of a sampler's data object with a factor in another place (found while building, pinned in
  test-derived-records.R). Slice 3 aligns the records with the columns (its step 3.3).

Tests. tinytest test-derived-records.R (177 expectations): the records at every door equal the reader and
the handed object carries none; no range under probit, logistic, ordinal and multinomial; copy, reload and a
twin agree on six fixtures; the live copy; the reload's moved splits; a current field; a dead pointer with no
saved state; the record carried by `new()`, `dbarts()`, `bart()` and `dbartsSpec()`, and
`samplePriorPredictive` on the sampler's grid; `xbart` deriving; the control's anchor not read, and the
record holding the anchor and unit across a response swap; a range carried to nbinom and not read by probit;
weights through a copy and a reload; misfits refused at creation and at `setData`, a grid sorted, an ordered
factor's other grid and an unordered factor's grid refused; a record-less object deriving; a constant response
of several forests copied. tests/cpp `testSuppliedCutGrid` (the three constructors) and a facade row for the
reader.

Existing tests changed beyond step 1.4. Each expected the adoption dec-B461 removes or the stored state that
dec-B459 has a live copy no longer take, none a value derived where the record now holds it:
test-copy-nsamples-edit.R (the copy's field is the live chain), test-cut-grid-distinct.R (a live copy of a
sampler holding a refused field is made; a reload of a sampler whose field records no weights keeps the
record's), test-data-mixed-mutation.R and test-mutate-sparse-valued.R (a reload after `setCutPoints` is on the
recorded points, equal to a live forced install), test-setstate-force-update.R (the repairs pinned on a reload
and on a copy of a dead sampler), test-state-empty-leaf-merge.R (the cold restore is a reload),
test-state-family.R (a live copy takes the chain), test-multiforest-leaf-prior-writer.R (the anchor read from
the data object), test-data-handle.R and test-sampler-boundary-refusals.R (the setState entry's argument
fewer), test-update-scale-rows.R (the internal writer's name). tests/cpp: the adoption pins of
`testWeightedGridState`, `testRestoreStatus` and `testLiveGridInstall` and the facade's setState row
rewritten to the install that keeps the store's grid, the ordered-factor adoption pin dropped.

Mutants, each alone on an export (scratch/sfp1/mut-src), test-derived-records.R unless named: the record of
a family that reads no range written c(0, 0), 13 failures; the range record not applied at creation, 7; a live
copy installing the stored state, 3; a reload adopting the state's grid, 3; a dead-pointer setState with no
saved state adopting the handed state's grid, 1; the copy's field the source's, 1; the creation check skipped,
the file stops at the bridge's backstop; the control's anchor written and read again, 3; the setData check
skipped, 1; the anchor recorded as 0, the copy refused; the weights dropped by the bridge, 7; the data handle
handed the record, 1 (xbart); the build ignoring a supplied grid, 18, and tests/cpp `testSuppliedCutGrid` 2 of
its checks; and the setData record handed on, the file stops (run on the build before 7b5dd467). Survived:
`Sampler::setAnchor` without its family guard, an equivalent mutant, a family without units having no
transform for the move to change. Not written: the caller's object written, the range supplied at
construction, the supplied path taken with derived values, and migration (above).

Size, against 46e1a644 before this record: engine 106 added and 86 removed (stop 320), bridge 163 and 48, R
319 and 178, tests and the gate harness 887 and 215, help, NEWS and notes 64 and 50; 1539 added and 577
removed in all (stop 2800), against about 1860 planned.

Gates, on the slice's libraries: tests/cpp plain and under ASAN and UBSAN, all passed, no diagnostic; the full
tinytest suite in one process, 27855 expectations, none failed; lint, air, rc-codoc, win-drift and
doc-freshness clean; on the reference build the four snapshot files unedited (27 expectations) and the three
`--bitwise` compares, 55 (`--strict-coverage`, 55 compared, 0 skipped), 15 and 11 identical, no "max |z|"
line; the two-library compare (scratch/sfp1/compare.R), 18 fixtures and 162 parts, every one identical but
`samplePriorPredictive` after `setCutPoints`, which this slice moves onto the sampler's grid; the exact gates
of exact-gates.yaml at `quick` on the shipped build, every gate and arm passing, forest-prior-exact.R once its
oracle read the unit from the data object (it read the control's, which no longer carries one); the R-loaded
ASAN run over the new and the twelve most changed test files, no diagnostic; `R CMD check --as-cran
--no-manual` from a clean tarball, the Date field's NOTE alone; stan4bart's suite (its bartcore branch at
ac80108, built against this slice), 1545 expectations, none failed. Not run: the x86-64 leg.

Left open. The plan's "fresh object neutral" and "copy and reload bitwise" literals, held instead by the
scratch compare above, which the reviewer reruns against 46e1a644. The ordered-factor refusal of a copy and a
reload of a record declaring other levels than its engine holds (first call above) is for the reviewer's and the
maintainer's mark.

Rebased onto a000d239 (2026-10-22), forest-multiplier push 2 having landed, which `git rebase` stopped on in
R/dbarts.R and man/dbartsData.Rd; the fold is in the slice's commit and the notes' and the harness's, and this
paragraph and the repointed cites are their own commit. Rebased again onto a563d753 (2026-10-22),
warm-start-salvage slice D having landed, with no stop, and in the fix round onto c1d1473f, forest-multiplier
push 3 having landed (below); the slice is now on c1d1473f, and the cites above are its hashes there. The
calls, each with its alternative:

- The first creation hands [`createEngine`](../../R/dbarts.R) push 2's mask, `active`, the rows a numeric
  column's scale is read over; every re-creation (`copy`, a reload, a dead-pointer `setState`,
  `samplePriorPredictive`) hands none and installs its mask after the install, as push 2's re-creations did.
  Alternative: the first creation calling the bridge itself, beside `createEngine`, two ways to make an
  engine.
- `data@basis.scale` keeps its own writer, [`recordBasisScales`](../../R/dbarts.R), called after
  [`recordDerived`](../../R/data.R) where the control declares forests, and is not folded into it. The derived
  records are rewritten from one engine reader after every call that can change one; the scale is written
  once, at the first creation and only where the object has no entry, is read off `getLeafPrior`, and raises a
  fallback warning that is the first creation's alone (Constraints: no warning beyond those named). Folding it
  in would add a second reader to `recordDerived` or a sixth value to `Sampler::derivedRecords`, an engine
  change, and would need a guard to keep the warning from firing again at every later write. Alternative: fold
  it in as a step that fills absent entries, the warning split out. Every engine `createEngine` makes is
  handed the recorded scale, as push 2's were: the bridge reads the slot at every parse and
  `checkDerivedRecord` leaves it alone. The Constraints' "construction takes supplied values for the grid and
  the standardization only" predates push 2, whose scale also goes in at construction; the fold adds nothing
  to it.
- Push 2's `recordAnchor` and its write of the anchor and unit onto the control's forests' record are
  dropped. Push 2 added no reader of the model's range; the unit and anchor its tests and its gate read off
  the control now come from `data@forest.calibration` in test-forest-units-twin.R, test-bcf-family.R, test-forest-leaf-prior-reader.R,
  test-forest-scale-record.R, test-forest-held.R, test-forest-sd-unit.R and forest-prior-exact.R's numeric
  oracle.
- man/dbartsData.Rd keeps both paragraphs, the second naming `basis.scale` as the record that keeps its first
  creation's value. The notes forest-multiplier.md, state-not-model.md, multiplier-combiner.md and
  forest-sd-unit.md put the anchor and the unit on the data object; the last two said the control before push
  2 and were missed here.

The fold is about 15 lines of R beyond the slice's, 7 test lines, one help paragraph, 9 lines of notes and one
of the harness; no engine or bridge line. Its two mutants, each on a clean export: `createEngine` not handing
the mask, test-forest-scale-record.R 6 failures; `recordBasisScales` not called at creation, 22.

Gates on the rebased tip, on its own libraries: tests/cpp 407 ok, plain and under ASAN and UBSAN, no
diagnostic; the full suite in one process, 28228 expectations, none failed; lint (0 lints), air, rc-codoc,
win-drift, doc-freshness and verify-anchors clean; on a `--preclean` reference install the four snapshot files
(14, 3, 7 and 3 expectations) and the three `--bitwise` compares, 55 of 55 identical against
equivalence-32418b6e (`--strict-coverage`, 55 compared, 0 skipped), 19 of 19 against bcf-equivalence-20eed7ed
and 11 of 11 against multinomial-equivalence-80b1c8d4, no "max |z|" line; the 43 commands of exact-gates.yaml at
`quick`, all passing; the R-loaded ASAN run over the 23 test files the slice and the fold touch, 5765
expectations, none failed, no diagnostic, with the engine and the bridge instrumented through `CXX20FLAGS`
(the run above instrumented the C files alone, R compiling the package's C++20 under `CXX20FLAGS` and not
`CXXFLAGS`); `R CMD check --as-cran --no-manual` from a clean tarball, the Date field's NOTE alone;
scratch/sfp1/compare.R against a000d239, 18 fixtures and 162 parts, every one identical but
`samplePriorPredictive` after `setCutPoints`, as against 46e1a644; scratch/sfp1/rb-compare-fold.R, six
fixtures of numeric forests (a default, a held one, probit with a 0/1 mask at creation, logistic, three forests
under weights, and a swapped basis with a response swap at `updateScale = FALSE`), 54 parts from creation
through copy, reload and a twin, every one identical to a000d239's; stan4bart's suite (its bartcore branch at
ac80108, built against the rebased tip), 1545 expectations, none failed. Not run: the x86-64 leg.

### Slice 1 fix round

On the slice's review (opus, 2026-10-22, on the tip before the second rebase above), rebased with the slice onto
c1d1473f (below): ae835866 (the fixes), 451de695 (the pins) and this record. Every pin is in
test-derived-records.R.

- A reload's first call handing back its own field (the review's B1): a reload of an `updateState = FALSE`
  sampler saved after `setCutPoints` with no store since takes `setState(r$state)`, TRUE, and its field is then
  its engine's.
- A dead copy's field (B2): a copy of that reload, made before its engine is rebuilt, holds its engine's state
  as its field and takes it back, TRUE.
- The prior predictive's grid (B3): its draws differ from those of a sampler made on the same data object with
  the grid record set to NULL, the same range on derived points.
- `setPredictor(updateCutPoints = "value")` in "carried at creation" is handed twice the creation's values; the
  record equals the reader right after the call and differs from the creation's grid.
- A carried grid out of order is sorted with its weights: `checkDerivedRecord` checks the `mass` attribute
  first, then reorders its successive differences with the points and accumulates them again from its first
  value. (0.75, 0.5, 0.25) weighing 1, 2 and 3 records (0.25, 0.5, 0.75) with mass c(0, 3, 5, 6), and the engine
  holds that.
- A hand-edited grid out of order at a reload is sorted by `recreatePointer`'s check, the stored splits going in
  unmoved.
- The record's cost. A response, offset, weights or active-row change at `updateScale = TRUE` writes the range
  alone, through [`recordResponseRange`](../../R/data.R), and only where it moved; the read takes a second
  argument of the reader's bridge entry, FALSE leaving the grids unread (`DerivedRecords::withGrids` in the
  engine; no virtual changes). A calibration map's anchor and unit are set at construction only, so none of
  those calls moves them. `setCutPoints`, `setPredictor` with a grid update and `setData` write every record from
  the full read. tests/cpp's facade row for the reader adds the read without grids. Each call alone in a fresh
  process on this machine, n = 1000, 50 trees, the best of ten runs of 400 calls, ms per call; base a563d753,
  before the rebased slice, after this round:

  | call | p = 500: base | before | after | p = 50: base | before | after |
  |---|---|---|---|---|---|---|
  | `setResponse(updateScale = TRUE)` | 0.065 | 0.170 | 0.087 | 0.065 | 0.097 | 0.088 |
  | `setResponse(updateScale = FALSE)` | 0.048 | 0.047 | 0.047 | 0.048 | 0.048 | 0.048 |
  | `setPredictor(updateCutPoints = "value")` | 0.580 | 0.450 | 0.450 | 0.220 | 0.267 | 0.267 |
  | `setPredictor()` | 0.572 | 0.578 | 0.585 | 0.215 | 0.215 | 0.227 |
  | `setOffset(updateScale = TRUE)` | 0.075 | 0.177 | 0.095 | 0.075 | 0.108 | 0.095 |
  | `setOffset(updateScale = FALSE)` | 0.057 | 0.057 | 0.058 | 0.057 | 0.057 | 0.058 |
  | `setWeights(updateScale = TRUE)` | 0.055 | 0.163 | 0.053 | 0.057 | 0.092 | 0.052 |
  | `setWeights(updateScale = FALSE)` | 0.038 | 0.037 | 0.037 | 0.038 | 0.038 | 0.037 |

  A sweep is 0.21 ms at p = 500 and 0.17 at p = 50 on all three. What is left over the base at
  `updateScale = TRUE`, about 0.02 ms, is the write of the data field where the range moved, most of it the
  reference class's check of the field's class; these weights leave the range where it was, and nothing is
  written. `setPredictor` with a grid update still reads every grid; its differences from the base are the
  slice's.
- A copy continues its source: its next 20 draws equal the source's to a tolerance of 1e-12, beside the bitwise
  pin against a second copy.
- Call (1) is corrected above: the code matrix reaches the ordered-factor refusal with no hand edit. The
  behavior is unchanged, for the maintainer.
- Call (2) is built as the orchestrator's reading of dec-B459, whose "refused, as today, without a stored state"
  covers the copy: `copy()` of a sampler whose engine is gone and that has no stored state is refused with the
  reload's message, from one helper, before anything is created, so the two agree. Alternative: the fresh copy
  the slice kept, the ruling's "as today" having been true of the reload alone. Pinned with
  `updateState = FALSE`, never run and run: the field is NULL and `copy()` and `run()` refuse alike. A live
  sampler's copy and a dead one's with a stored state are unchanged. The `copy` docstring, the sampler's help
  (`shallow` and Saving), bart's help, which still said a copy installs the stored state and that one without
  starts over, and NEWS say so.
- TODO: code-matrix-category-count adds the ordered case; legacy-object-message records the review's finding on
  0.9-34 objects, the same before the slice.
- Forest-multiplier push 3 (c1d1473f, a printed block of each forest's prior) landed during the round, and the
  slice and the round were rebased onto it. `git rebase` stopped on the copy method in the slice's commit and in
  the fixes, push 3 wrapping the creation in `withoutForestBlock` beside the slice's docstring and this round's
  refusal; both are kept. Push 3 read the unit off the control's forests record, which the slice empties:
  [`packageBartResults`](../../R/bart.R) now sets a fit's `forest.unit` from `dataForestCalibration(fit$data)`
  and [`samplerForestLines`](../../R/forestPrint.R) reads it off `sampler$data`, the control's record kept for
  the labels and basis terms; forest-multiplier.md says a sampler holds the unit on its data object. The fold
  is in the slice's commit, two reads and a comment in R and one line of the note. A verbose creation prints
  the block after the record is written. No other read of the control's anchor or unit is in R/ or
  inst/tinytest. test-forest-print.R, 138 results, passes, its `forest.unit` pins included; each of the two
  reads put back on the control fails it, 10 and 14 failures. The timings above were taken on a563d753, push 3
  touching no setter.

Mutants, each alone on an export with a plain install (scratch/sfp1fix/mut.py, on the review's definitions),
test-derived-records.R (195 expectations) unless named: M19, a dead copy's field the installed state, 2
failures; M22, no own-field refresh in a dead-pointer `setState`, 2; M29, `samplePriorPredictive` on a NULL
grid record, 1; M32, no check in `recreatePointer`, the file stops at the bridge's refusal; M4, `setPredictor`
writing no record, 2, this file now killing it alone. The range writes, each removed: `setResponse`'s, 1;
`setOffset`'s, test-aft.R 4; `setWeights`'s, test-update-scale-rows.R 7; `setActiveRows`'s, the same file 5. The
sort reordering no weight, 2. Not written, both equivalent but for cost: the grids read under FALSE, and the
range written where it did not move.

Size: package code 81 added and 41 removed (engine 7 and 2, bridge 12 and 8), tests 55 and 5, help, NEWS and
notes 17 and 14, beside this record and TODO.

Gates on 451de695, the folded tip, on its own libraries: tests/cpp 408 ok, plain and under ASAN and UBSAN, no
diagnostic; the full suite in one process, 28439 expectations, none failed; lint (0 lints), air, rc-codoc,
win-drift, doc-freshness and verify-anchors clean, NEWS.Rd parsing to 226 entries; on a `--preclean` reference
install the four snapshot files (14, 3, 7 and 3 expectations) and the three `--bitwise` compares, 55 of 55
identical against equivalence-32418b6e (`--strict-coverage`, 55 compared, 0 skipped), 19 of 19 against
bcf-equivalence-20eed7ed and 11 of 11 against multinomial-equivalence-80b1c8d4, no "max |z|" line; the 43
commands of exact-gates.yaml at `quick` on the shipped build, all passing; against c1d1473f,
scratch/sfp1/compare.R, 18 fixtures and 162 parts, every one identical but `samplePriorPredictive` after
`setCutPoints`, scratch/sfp1/rb-compare-fold.R, 6 fixtures and 54 parts, all identical, and a fit's
`forest.unit` with its print, summary and a sampler's show, identical; the R-loaded ASAN run over the review's
23 files and test-forest-print.R, 5921 expectations, none failed, no diagnostic, the bridge built under
`CXX20FLAGS`; `R CMD check --as-cran` from a clean tarball, with the manual, the Date field's NOTE alone;
stan4bart's suite (its bartcore branch at ac80108, built against this tip), 1545 expectations, none failed. Not
run: the x86-64 leg.
