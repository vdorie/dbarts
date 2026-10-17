# state-frame-prior: the derived data's records on dbartsData, and every door reading them

Status: PLANNED 2026-10-17, not critiqued; seven calls open with the maintainer
([Open questions](../design/state-frame-prior.md#open-questions)), the steps written to each recommendation.
Design: [state-frame-prior.md](../design/state-frame-prior.md). Rests on dec-B231 to dec-B233, dec-B237,
dec-B304, dec-B418, dec-B449, dec-B450 and dec-A159 (confirmed under dec-B378).

agent: per slice, one opus implementer and one opus reviewer who runs that slice's mutants; a blind critique
of this plan before slice 1. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)) for all
three: each edits the engine's creation, `setData` or install path.
rng, per slice, against equivalence-32418b6e, bcf-equivalence-1b7d730c and multinomial-equivalence-80b1c8d4
(bcf-equivalence-0fdeb37e if forest-multiplier push 1 lands first):
- 1 is NEUTRAL, bit for bit, for every creation from a data object that carries no record (every fitting
  function, every fresh `dbartsData()`), every live call, and every copy, reload and dead-pointer `setState`
  whose stored state records the grid and standardization in force. POSTERIOR-CHANGING on two sequences: a
  copy or reload whose stored state predates a grid change (the grid, part of the tree prior, is the
  sampler's and no longer the state's); a creation from a data object that carries a record other than a
  fresh derivation's (`new()`, `dbarts()` or `dbartsSpec()` handed `s$data`).
- 2 is NEUTRAL for `setData` of a data object carrying no record and for `setPredictor` without
  `updateStandardization`. POSTERIOR-CHANGING for `setData` of a carried record (held where the tip derives)
  and for the new arguments. Accepted where refused: a gp sampler's `setData` that keeps the standardization
  and range.
- 3 is NEUTRAL for every install under the standardization in force. POSTERIOR-CHANGING for an install across
  a standardization: a linear sampler keeps its own and the slopes are converted where the tip installs the
  state's; refused where the tip installs: kept gp draws under another standardization.
- No recorded scenario reaches a changed sequence: the `setdata` and `hetsetdata` scenarios hand a fresh
  `dbartsData()`, the `setPredictor` scenarios name no new argument, and none copies, reloads or calls
  `setState`. Expected 55, 15 and 11 identical, no re-record; a moved stream is a stop.
window: before the merge to main. Engine slices stay serial: 1 after warm-start-salvage slice C and
forest-multiplier push 1 land, then 2, then 3, then response-scale-rows slice D.
budget: 1 about 1660 lines, 2 about 1480, 3 about 830; [Budget and stops](#budget-and-stops).

## Goal

Each derived value has one record, on the data object: `cutPoints` (with its weights), `standardization`,
`response.range` and `forest.calibration`. Creation and `setData` use a carried record, derive an absent one
and follow an explicit argument; `setPredictor` takes `updateStandardization`; a state install compares the
standardization and converts linear slopes instead of installing it; a copy and a reload take every value
from the record and none from the state; the gp `setData` refusal covers only calls that would change a
kept draw; the help carries the table of what puts each value back. Closes TODO `state-frame-prior` and
unblocks response-scale-rows slice D.

## Context

- The rule and its doors: [The rule at every door](../design/state-frame-prior.md#the-rule-at-every-door).
  What moves, measured: [What changes for users](../design/state-frame-prior.md#what-changes-for-users)
  (probes p1 to p3 in scratch/sfd, run on the tip's shipped build and on 0.9-34).
- Today's homes: the grid adopted from a state at re-creation ([`installStateOnto`](../../R/dbarts.R),
  `adoptGrid` in [`Sampler::setState`](../../src/bartcore/sampler.hpp)); the standardization installed from a
  state ([`Chain::setState`](../../src/bartcore/chain.hpp), `restoreCalibration` on
  [`LinearGaussianLeaf`](../../src/bartcore/model.hpp) and [`GPGaussianLeaf`](../../src/bartcore/model.hpp))
  and derived at `setData` ([`Chain::applyNewData`](../../src/bartcore/chain.hpp)); the range on the model
  ([`recordAnchor`](../../R/dbarts.R), [`applyAnchor`](../../R/dbarts.R)); the anchor and unit on the control
  (the sampler's `initialize` in [R/dbarts.R](../../R/dbarts.R), read by
  [`applyForestAttributes`](../../src/R_interface_bartcore.cpp)).
- A record pattern to follow: `missing.seen`, read through [`dataMissingSeen`](../../R/data.R), written in
  [`bartcoreSamplerSetData`](../../R/bartcore.R) from [`bartcore_getMissingSeen`](../../src/R_interface_bartcore.cpp);
  checks by a plain function, [`checkBasisRecord`](../../R/data.R).
- Conversions that exist: [`convertLinearCoefficients`](../../src/bartcore/model.hpp),
  [`convertFlatLinearTree`](../../src/bartcore/chain.hpp), the kept-draw rewrite at a re-anchor
  ([`Chain::restateSavedDraws`](../../src/bartcore/chain.hpp)).
- The gp refusals: [`holdsSavedGPDraws`](../../src/R_interface_bartcore.cpp),
  [`refuseSavedGPDrawReanchor`](../../src/R_interface_bartcore.cpp), the outright one in
  [`bartcore_setData`](../../src/R_interface_bartcore.cpp), and
  [`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp).

## Decision

The seven open questions of the design note. What a different answer changes:

1. Logical `updateCutPoints`: step 2.3 loses `"position"`. A fourth value: one more branch in step 2.3.
2. A warning on a held grid or range that does not span the new data: about 40 lines and a class in slice 2.
3. Not clean: slice 3's verdict returns FALSE unforced where a slope converts. Installing as the tip does:
   slice 3 shrinks to the gp refusal and the table, and the design's sentence on the frame changes.
4. By value: slice 2's refusal compares after deriving, in the engine.
5. A copy reading the live state: about 20 lines in slice 1 and a sentence of `copy`'s help.
6. No migration: slice 1 deletes the adopting path (about 50 lines fewer, two tests of reloading a
   record-less object rewritten to expect a derived grid). Refusal: one message in its place.
7. A fifth slot for the weights, or a `dbartsData()` argument: slice 1's record and its check.

## Constraints

- A creation from a data object without a record runs the code that runs today; so does `setData` of one.
- The range goes into a re-created engine by the call `applyAnchor` makes today, after construction, so a
  copy and a reload stay bit for bit; construction takes supplied values for the grid and the
  standardization only.
- The engine names no R type and no door; [dbarts.h](../../inst/include/dbarts/dbarts.h) and the C API hash do
  not move. [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) stays at 1; a state keeps every block.
- No message or warning beyond the refusals named here (dec-B234 for `setState`).
- Out of scope: the design's [Not here](../design/state-frame-prior.md#not-here).

## Slice 1: the records, and every re-creation reading them

### Steps

1.1 Engine (about 110 lines, [data.hpp](../../src/bartcore/data.hpp), [model.hpp](../../src/bartcore/model.hpp),
    [sampler.hpp](../../src/bartcore/sampler.hpp), [facade.hpp](../../src/bartcore/facade.hpp)). Optional
    supplied values at construction: per column a grid and its weights, installed before the columns are
    quantized; per forest a center, scale and lengthscale, served to the leaf's derivation as
    [`ColumnStore::suppliedStandardization`](../../src/bartcore/data.hpp) serves a view's. Absent derives
    through today's code. One reader of the grids, weights, standardizations, range, anchor and unit.
    `adoptGrid` stays. A new virtual: `--preclean`.
1.2 Bridge (about 140 lines). [`parseData`](../../src/R_interface_bartcore.cpp) reads the four slots, each
    guarded for an object without it, and hands the supplied values to
    [`bartcore_create`](../../src/R_interface_bartcore.cpp); the anchor and unit come from
    `data@forest.calibration`, from the control's attribute only for an object saved without the slot. A
    reader entry returns the four records.
1.3 R (about 260 lines). The four slots and prototypes on [`dbartsData`](../../R/A_class.R); accessors in
    [R/data.R](../../R/data.R); `checkDerivedRecord`, a plain function, at the creation doors and `setData`
    (refusal text names the slot and says `NULL` derives it). The sampler writes the records into its own
    data object after creation and after `setCutPoints`, `setPredictor`, `setData` and each `updateScale =
    TRUE`; [`recordAnchor`](../../R/dbarts.R) writes `data@response.range` and the model's attribute goes;
    the anchor and unit leave the control. `copy`, [`getPointer`](../../R/dbarts.R) and
    [`recreatePointer`](../../R/dbarts.R) re-create from the record and install with `adoptGrid` false, except
    the migration (question 6): a data object without the record adopts from the state as today, takes the
    range from the model's attribute and the anchor and unit from the control's, and is written the record.
    Docstrings of `copy`, `setState` and `storeState`.
1.4 Existing tests. The 20 lines in 13 files reading `attr(model, "response.range")` and the 2 reading the
    control's anchor read the data object; the copy and reload cases of
    [test-state-grid-install.R](../../inst/tinytest/test-state-grid-install.R) that pin adoption are split:
    a record-less object keeps the pin, a recorded one pins the record.
1.5 Tests (about 1000 lines with 1.4): tinytest `test-derived-records.R`, tests/cpp in `test_state.cpp`.
1.6 Help, NEWS, notes (about 150 lines): man/dbartsData.Rd (the slots), the sampler's Saving section and
    `copy`; NEWS (the edited-object creation the design measures); one sentence each in
    [state-not-model.md](../design/state-not-model.md) (the frame bullet and the first-creation bullet),
    [install-surface.md](../design/install-surface.md), [cut-grid.md](../design/cut-grid.md) and
    [architecture.md](../architecture.md).

### Tests

| Test | Pins | Mutant |
|---|---|---|
| written at creation | every door: the four records equal the engine's reader; the caller's object carries none | the caller's object written |
| fresh object neutral | 16 fixtures (each leaf kind, weights, a mask, a factor, a sparse store, several forests): creation from a fresh object `identical()` in every reader and the next 20 draws to the base build | the supplied path taken with derived values |
| copy and reload bitwise | after `storeState()`: state, readers, next 20 draws equal the base build's copy and reload, on the same fixtures | the range supplied at construction |
| copy holds the sampler's grid | `updateState = FALSE`, 3 points set: the copy and the reload have 3 points and the stored trees moved by value (p3: 16 splits) | adopting from the state |
| carried at creation | `new()` and `dbarts()` handed `s$data` after `setCutPoints`, after `setPredictor(updateCutPoints = "value")` and after `setResponse(updateScale = TRUE)`: the sampler's grid, weights and range | the record ignored |
| weights ride | a 30-value integer column: copy and reload keep its weights from the record | weights dropped |
| misfit refused | wrong column count, a grid not increasing, an ordered factor's other grid, a standardization for other columns, a range of length 1: refused with the text, nothing created | the check skipped on one door |
| migration | an object with the slots removed (as test-missingness-first-seen.R removes one), a grid set and a state stored: the reload is the base build's bit for bit and the record is written | deriving afresh |
| anchor and unit | several forests: the record is the creation's after `setForestBasis`, `copy` and a reload, on a constant response (forest-multiplier push 1's case) too | the control's attribute read |
| column order | `setData` of the same columns in another order (dec-A215): the record's entries follow their names | the record kept by position |

## Slice 2: setData and setPredictor follow the rule

### Steps

2.1 Engine (about 230 lines). [`Sampler::setData`](../../src/bartcore/sampler.hpp) takes, for the grid, the
    standardization and the range each, derive, keep or a supplied value, and the placement for a derived
    grid. Keep: the new rows quantized on the grid in force; the leaf covariates regathered under the
    constants in force for the new row count; the response installed without moving the transform, the kept
    draws then not rewritten. A supplied value goes in as a carried one at creation does. Derive: today's path.
    `setPredictor` and `updatePredictor` take `updateStandardization`: the leaf's derivation over the new
    values, then the conversion of live and kept linear coefficients, as at `setData`.
2.2 Bridge (about 90 lines). [`bartcore_setData`](../../src/R_interface_bartcore.cpp) takes the three
    choices; its outright gp refusal becomes question 4's: refused where the call derives the
    standardization or the range or brings one other than the sampler's, in today's words with the remedy
    (`updateStandardization = FALSE, updateScale = FALSE`). The predictor entries take the flag; a gp sampler
    holding kept draws refuses `TRUE` (dec-B237).
2.3 R (about 140 lines). `setData(newData, updateState = NULL, updateCutPoints = NULL, updateScale = NULL,
    updateStandardization = NULL)` resolves each choice by question 1; the carried records are aligned with
    the columns before the call ([`alignDataFactorLevels`](../../R/bartcore.R)'s order); the fallback
    warning ([`scaleRowsFallBack`](../../R/spec.R)) is tested only where the range is derived.
    `setPredictor(..., updateStandardization = FALSE)`.
2.4 Existing tests: the three setData pins of
    [test-response-scale-rows.R](../../inst/tinytest/test-response-scale-rows.R) that hand back the
    sampler's object and expect the fallback state `updateScale = TRUE`; the gp refusal pins of
    [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R) are split into refused and served.
2.5 Tests (about 900 lines with 2.4).
2.6 Help, NEWS, notes (about 120 lines): the arguments; NEWS: `setData` uses the record a data object
    carries, with the three ways to derive; TODO.

### Tests

| Test | Pins | Mutant |
|---|---|---|
| fresh object neutral | `setData(dbartsData(...))` on the slice 1 fixtures: `identical()` to the base build in readers and 20 draws | the keep path on an absent record |
| carried held | p1's two sequences: grid and range as before the call, by reader; kept draws not rewritten | deriving a carried record |
| arguments win | each argument at each value over a carried and an absent record: 2 x 3 x 3 cases, by reader | `NULL` read as keep |
| keep is the sampler's | `"none"` and `FALSE` with a carried record other than the sampler's: the sampler's values | the carried one used |
| placement | `"position"` and `"value"` on a derived grid: splits as `setPredictor` places them | one placement for both |
| new row count | keep with 40 rows added: codes on the old grid, covariates under the old constants, fits finite | covariates left at the old length |
| linear converts | `updateStandardization = TRUE` at both setters: live and kept predictions unchanged to 1e-12 | conversion skipped for kept draws |
| gp narrowed | kept gp draws: `setData(s$data)` served; a fresh object refused; `FALSE, FALSE` served; a carried record of another standardization refused; `setPredictor(updateStandardization = TRUE)` refused; nothing touched by a refusal | the refusal by argument alone |
| record written | after each call the slots equal the reader | a stale slot |
| fallback where derived | the three response-scale-rows sequences: warned with `updateScale = TRUE`, silent held | warning on a held range |

## Slice 3: a state install compares the standardization, and the manual's table

### Steps

3.1 Engine (about 110 lines). [`Chain::setState`](../../src/bartcore/chain.hpp) on a live install: per forest,
    an equal standardization runs today's code; another one leaves the leaf's own and converts the live
    and kept linear coefficients from the state's into it; a gp forest whose kept draws were made under
    another center, scale or heuristic lengthscale is refused beside
    [`Chain::lengthscaleStateFeasible`](../../src/bartcore/chain.hpp). The verdict by question 3 (clean). The
    adopting install (the migration) keeps `restoreCalibration`.
3.2 Bridge and R (about 30 lines): the refusal's text, beside the lengthscale one; docstrings.
3.3 Tests (about 550 lines).
3.4 Help, NEWS, notes (about 140 lines): the design's table, [The manual's table](../design/state-frame-prior.md#the-manuals-table),
    in the sampler's Saving section, and the `setState` sentences; NEWS; [leaf-conversions.md](../design/leaf-conversions.md)
    and [state-not-model.md](../design/state-not-model.md) (the standardization row of its table); the design
    note's Status; TODO drops `state-frame-prior` and says slice D of response-scale-rows is unblocked.

### Tests

| Test | Pins | Mutant |
|---|---|---|
| same frame exact | every restore of the slice 1 fixtures under the standardization in force: byte-identical to the base build | conversion run on equal standardizations |
| keeps its own | p2's sequence: after `setState(st)` the sampler reports its own center and scale; live fit and kept predictions equal the stored ones to 1e-12; TRUE unforced | the state's installed |
| without spread | a covariate without spread on one side, then on both: the conversion of leaf-conversions.md's rules | slope set to zero |
| gp kept refused | kept gp draws under another standardization, and under another heuristic lengthscale: error, sampler untouched; live fits alone install | live fits refused too |
| migration installs | a record-less reload installs the state's standardization as today | the live rule on adoption |
| table sequences | each row of the table, end to end: the value back, `setState` TRUE, next 5 draws equal a twin that only restored | (pins the help) |

## Verification

Per slice, against the slice's own library ([Landing](README.md#landing)):

- `cd tests/cpp && make && ./test_bartcore`, plain and under ASAN and UBSAN; `R CMD INSTALL --preclean`; the
  full tinytest suite; the R-loaded sanitizer run over the new and changed test files.
- Reference build: the four `test-reproducibility-*.R` files unedited, and the three `--bitwise` compares
  (`--strict-coverage` for the main harness), 55, 15 and 11 "identical draws" lines and no "max |z|" line.
- The 42 exact-gate commands at `quick`; expected unchanged, none creating from a sampler's data, copying,
  reloading or restoring across a standardization.
- Lint, air, rc-codoc, win-drift and doc-freshness; `R CMD check --as-cran` from a clean tarball; the x86-64
  leg; stan4bart's suite (its restores are first creations from its own data object, which carries no
  record, and forced same-grid installs).

## Budget and stops

| | 1 | 2 | 3 |
|---|---|---|---|
| engine | about 110 | about 230 | about 110 |
| bridge and R | about 400 | about 230 | about 30 |
| tests | about 1000 | about 900 | about 550 |
| help, NEWS, notes | about 150 | about 120 | about 140 |
| planned | about 1660 | about 1480 | about 830 |
| stop, engine lines added | 220 | 380 | 180 |
| stop, all | 2600 | 2300 | 1300 |

The last state slices ran 1.5 to 2 times over in tests, and the stops follow. Stop and report, building
nothing further: a bitwise compare prints "max |z|" or a snapshot moves; a copy or a reload of a fixture with
a current state differs from the base build in any byte or draw; a test outside those named in steps 1.4 and
2.4 fails because it expected a value derived where the record now holds it (list it, do not edit it); a
slice's stop in the table.

## Order and interactions

- warm-start-salvage slice C (built on wt/warm-start-c, the install paths in sampler.hpp and chain.hpp) and
  forest-multiplier push 1 (built on wt/forest-mult-p1, [`parseData`](../../src/R_interface_bartcore.cpp) and
  the anchor on a constant response): both land before slice 1, which rebases on them. Slice 1 moves the
  anchor push 1 fixes and carries its value; its test "anchor and unit" runs push 1's constant-response case.
- forest-multiplier push 2 adds `basis.scale` to the data object and writes it where the anchor is written;
  whichever lands second folds the two write-backs into one (that plan's "With the state-frame arc's data
  record").
- response-scale-rows slice D follows slice 3. It needs slice 1 (a copy and a reload hold the center and
  scale through the record), slice 2 (`updateStandardization`) and slice 3 (a state install no longer
  replaces them). If D is wanted sooner it can follow slice 2, at the cost that until slice 3 a state stored
  before an `updateStandardization` call restores the older numbers.
- The control-migration arc (dec-B239) moves the forests' record to the model; the anchor and unit leave the
  control here first, so that arc moves the rest.

## Not verified

No code is built. The design's measured numbers are the tip's behavior and emulations of holding with existing calls
(probes p1 to p3), not runs of the rule; the conversion's predictions under slice 3 are argued from
leaf-conversions.md's measured 1e-14, not run. Not read: bartCause and treatSens beyond `git grep` (neither
calls `setData`, `setPredictor`, `setState` or reads a record), and the line estimates are from the code read
and the last state slices' ratios, not an itemization.
