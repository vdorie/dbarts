# state-frame-prior: the derived data's records on dbartsData, and every door reading them

Status: PLANNED 2026-10-17; revised the same day from its blind critique (scratch/sfdc-crit.md, sound after
revisions; [Revision from critique](#revision-from-critique)). Seven questions open with the maintainer
([Open questions](../design/state-frame-prior.md#open-questions)), the steps written to each recommendation.
Design: [state-frame-prior.md](../design/state-frame-prior.md). Rests on dec-B231 to dec-B233, dec-B237,
dec-B304, dec-B305, dec-B310, dec-B329, dec-B330, dec-B418, dec-B449, dec-B450 and dec-A159 (confirmed under
dec-B378).

agent: per slice, one opus implementer and one opus reviewer who runs that slice's mutants. Tier "Changes
draws" ([Process by risk](README.md#process-by-risk)) for all three: each edits the engine's creation,
`setData` or install path.
rng, per slice, against equivalence-32418b6e, bcf-equivalence-1b7d730c and multinomial-equivalence-80b1c8d4
(bcf-equivalence-0fdeb37e if forest-multiplier push 1 lands first). No recorded scenario reaches a changed
sequence: the `setdata` and `hetsetdata` scenarios hand a fresh `dbartsData()`, the `setPredictor` scenarios
name no new argument, none copies, reloads, calls `setState` or `samplePriorPredictive`, or hands one
sampler's control or data to another. Expected 55, 15 and 11 identical in every slice, no re-record; a moved
stream is a stop.
- 1, grid, range and anchor. NEUTRAL, bit for bit, for every creation from a data object that carries no
  record, every live call, `xbart` (its handle ignores records), and every copy, reload and dead-pointer
  `setState` whose stored state records the grid in force. POSTERIOR-CHANGING on six sequences, each a
  changed tree prior or leaf prior: a copy or reload whose stored state is older than a grid change (the
  record's grid, design question 3); a dead-pointer `setState` with no saved state (the record's grid, not
  the handed state's); a creation from a data object carrying a record other than a fresh derivation's
  (`new()`, `dbarts()`, `bart()`, `bart2()` or `dbartsSpec()` handed `s$data`); `samplePriorPredictive` (the
  sampler's grid, where it derived one from `data@x`); a control of several forests handed other data (the
  data's anchor, not the control's); a reload of an object saved before this slice whose grid or range was
  changed after creation (derived, design question 5).
- 2, standardization. NEUTRAL for every install and re-creation under the standardization in force.
  POSTERIOR-CHANGING for an install across a standardization: the sampler keeps its own where the tip
  installs the state's, and the slopes go in as stored (design question 1, option (d)). Refused where the tip
  installs: kept gp draws under another standardization or heuristic lengthscale.
- 3, setData and setPredictor. NEUTRAL for `setData` of a data object carrying no record and `setPredictor`
  without `updateStandardization`. POSTERIOR-CHANGING for `setData` of a carried record (held where the tip
  derives) and for the new arguments. Accepted where refused: a gp sampler's `setData` that keeps the
  standardization and range.
window: before the merge to main. Engine slices stay serial (Order and interactions).
budget: 1 about 1840 lines, 2 about 980, 3 about 1650; [Budget and stops](#budget-and-stops).

## Goal

Each derived value has one record, on the data object: `cutPoints` (with its weights), `standardization`,
`response.range` and `forest.calibration`. Creation and `setData` use a carried record, derive an absent one
and follow an explicit argument; `setPredictor` takes `updateStandardization`; a state install keeps the
sampler's standardization; a copy and a reload take every value from the record and none from the state; the
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
  ([`recordAnchor`](../../R/dbarts.R), [`applyAnchor`](../../R/dbarts.R), copied by `setModel`); the anchor
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

The design's seven open questions. What another answer changes:

1. (a) or (b): slice 2 adds the conversion of live and kept linear coefficients (about 60 engine lines, 150
   of tests), and under (b) a FALSE verdict. (c): slice 2 keeps `restoreCalibration` and shrinks to the gp
   refusal.
2. (b): about 40 lines and a warning class in slice 3. (c): slice 3 compares a carried record with the
   sampler's before using it.
3. (b): about 20 lines in slice 1 and a sentence of `copy`'s help. (c): slice 1 keeps `adoptGrid` on the
   copy and reload paths.
4. (b): slice 3 loses `"position"`. (c): one more branch in step 3.3.
5. (a): about 50 lines of R in slice 1, and slice 1 keeps `adoptGrid` and slice 2 `restoreCalibration` for
   record-less objects. (c): one refusal in slice 1.
6. (b): slice 3's refusal ignores a carried record instead.
7. (a): slice 1's record carries no family, about 15 lines fewer.

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
    in [`Sampler::setState`](../../src/bartcore/sampler.hpp) go (design question 5, (b)). `--preclean`.
1.2 Bridge (about 160 lines). [`parseData`](../../src/R_interface_bartcore.cpp) reads `cutPoints`,
    `response.range` and `forest.calibration`, each guarded for an object without it; the two creating
    factories take them; [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp) skips them. The
    anchor and unit come from `data@forest.calibration`. The reader entry. A range record carries its family
    (question 7) and is skipped under another.
1.3 R (about 260 lines). The three slots on [`dbartsData`](../../R/A_class.R), accessors in
    [R/data.R](../../R/data.R), `checkDerivedRecord` at the creation doors and `setData` (a grid out of order
    sorted, a repeated point refused, a misfit refused naming the slot). The sampler writes them after
    creation and after `setCutPoints`, `setPredictor`, `setData` and each `updateScale = TRUE`;
    [`recordAnchor`](../../R/dbarts.R) writes `data@response.range` (`NULL` under probit, logistic, ordinal and
    multinomial) and `setModel` no longer carries the attribute; the anchor and unit leave the control and
    its `initialize` block. `copy`, [`getPointer`](../../R/dbarts.R), [`recreatePointer`](../../R/dbarts.R),
    the dead-pointer `setState` and [`samplePriorPredictive`](../../R/dbarts.R) re-create from the record.
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
| copy holds the sampler's grid | `updateState = FALSE`, 3 points set: copy and reload have 3 points and the stored trees moved by value (16 splits), multinomial too | adopting from the state |
| dead pointer, no saved state | `setState(st)`: the record's grid | the handed state's grid adopted |
| carried at creation | `new()`, `dbarts()`, `bart()` and `dbartsSpec()` handed `s$data` after `setCutPoints`, `setPredictor(updateCutPoints = "value")` and `setResponse(updateScale = TRUE)`: the sampler's grid, weights and range; `samplePriorPredictive` on the sampler's grid | the record ignored at one door |
| xbart ignores | `xbart(s$data)` after `setResponse(updateScale = TRUE)`: each fold's range its own, equal to `xbart` of a fresh object | the handle handed the record |
| control anchor | a control of several forests handed other data: the data's anchor (q7's 0.0796, not 0.190) | the control's read |
| family tag | a gaussian sampler's data handed to an aft and a probit sampler: derived | the range used across families |
| weights ride | a 30-value integer column: copy and reload keep its weights | weights dropped |
| misfit and order | wrong column count, an ordered factor's other grid, a range of length 1: refused, nothing created; an unsorted grid sorted; a repeated point refused | the check skipped at one door |
| record-less object | the slots removed (as test-missingness-first-seen.R removes one): derived, as a fresh object | migration |
| constant response | several forests on a constant response (forest-multiplier push 1's case): the anchor recorded and a copy created | the record's anchor of 0 |

## Slice 2: the standardization record, and a state install that keeps it

### Steps

2.1 Engine (about 120 lines, [model.hpp](../../src/bartcore/model.hpp), [chain.hpp](../../src/bartcore/chain.hpp),
    [data.hpp](../../src/bartcore/data.hpp)). Supplied center, scale and heuristic lengthscale per forest at
    construction, through both creating factories, served as
    [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp) serves a view. The reader returns
    them. [`Chain::setState`](../../src/bartcore/chain.hpp) no longer calls `restoreCalibration`: the
    sampler's standardization stays and the slopes go in as stored (question 1, (d)); a gp forest whose kept
    draws were made under another center, scale or heuristic lengthscale is refused beside
    [`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp).
2.2 Bridge and R (about 110 lines): the slot, its accessor and check, its write-back; the refusal's text
    beside the lengthscale one; docstrings.
2.3 Existing tests: the pins of a restore that expects the state's standardization
    ([test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R), tests/cpp
    [`testWarmStartStandardization`](../../tests/cpp/test_state.cpp) is the warm start's and stays).
2.4 Tests (about 650 lines with 2.3).
2.5 Help, NEWS, notes (about 100 lines): the `setState` sentences (kept linear draws predict under the
    sampler's standardization after a restore across a change); the standardization row of
    [state-not-model.md](../design/state-not-model.md); [leaf-conversions.md](../design/leaf-conversions.md).

### Tests

| Test | Pins | Mutant |
|---|---|---|
| same frame exact | every restore of slice 1's fixtures under the standardization in force: equal to the recorded literals | a branch on a standardization that does not differ |
| keeps its own | p2's sequence: after a forced `setState(st)` the sampler reports (6.55, 0.807); slopes as stored, by the state block; NULL forced, TRUE unforced where the trees fit | the state's installed |
| copy holds it | a copy and a reload after `setPredictor(updateStandardization = TRUE)` emulated through `setData`: the record's standardization | derived at re-creation |
| gp kept refused | kept gp draws under another standardization and under another heuristic lengthscale: error, sampler untouched; live fits alone install | live fits refused too |
| supplied lengthscale | a gp leaf whose lengthscale the model supplies: no `lengthscale` in the record, the model's in force | the record's read over the model's |

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
    dec-B330) or brings a carried one whose numbers differ from the sampler's (question 6), with the remedy.
    The predictor entries take the flag; a gp sampler holding kept draws refuses `TRUE`.
3.3 R (about 160 lines). `setData(newData, updateState = NULL, updateCutPoints = NULL, updateScale = NULL,
    updateStandardization = NULL)` resolves each choice by question 4; a record the sampler wrote counts as
    carried (question 2); the carried records are aligned with the columns before the call
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
| engine | about 170 | about 120 | about 230 |
| bridge and R | about 420 | about 110 | about 250 |
| tests | about 1100 | about 650 | about 1000 |
| help, NEWS, notes | about 150 | about 100 | about 170 |
| planned | about 1840 | about 980 | about 1650 |
| stop, engine lines added | 320 | 220 | 380 |
| stop, all | 2800 | 1500 | 2500 |

Sized from the landed history: state-grid-install planned 250 engine lines and landed 275, and 1475 in all
against 2092 at the implementer's tip (1.4 times), tests 1549 against 1000. Slice 1 names three factories and
a facade reader where that slice had one install path. Under question 1 (a) or (b) slice 2 grows by about 60
engine and 150 test lines. Stop and report, building nothing further: a bitwise compare prints "max |z|" or a
snapshot moves; a fixture's recorded literal differs where the slice is to be neutral; a test outside those
named in steps 1.4, 2.3 and 3.4 fails because it expected a value derived where the record now holds it
(list it, do not edit it); a slice's stop in the table.

## Order and interactions

The serial engine order: warm-start-salvage slice C (wt/warm-start-c), forest-multiplier push 1
(wt/forest-mult-p1), warm-start-salvage slice D, then slices 1, 2 and 3 here, then response-scale-rows slice D.

- warm-start-salvage slice C: the install paths in sampler.hpp and chain.hpp; slice 1 rebases on it.
- forest-multiplier push 1: [`parseData`](../../src/R_interface_bartcore.cpp) and the anchor on a constant
  response; slice 1 moves that anchor and carries its value ("constant response").
- warm-start-salvage slice D: leaf kinds, `readWarmStartState` and
  [`Chain::convertDonorStandardization`](../../src/bartcore/chain.hpp), beside slice 2's state reader; it lands
  first so slice 2 rebases on one reader of the standardization block.
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

Item by item, the critique's numbering:

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
