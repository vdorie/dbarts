# control-migration: a control holds run settings only

Status: PLANNED 2026-10-30 on c55a1477, against main as it will be after tree-prior-flat (5b050cdf) lands. Not
critiqued yet. One call held for the maintainer ([Open calls](#open-calls)). Rulings: dec-B239, dec-B240,
dec-B241 (the move), dec-B250 (when the mirrored slots go), dec-B334, dec-B335, dec-B385 (a control copied from a
fit), dec-B247 (a live sampler's count), dec-B249 (the forest column), dec-B503 and dec-B505 (counts among
forests with no basis), dec-B518 (the live forest door), with dec-A02 (relocated arguments), dec-B116 (a flat
argument beats the control's), dec-B231 (split points are data), dec-B254 (structure is refused by setModel),
dec-B274 (defaults by kind) and dec-B512 (one rule at two levels). Earlier work: a design of 2026-10-05 and two
critiques of it, outside the tree; this plan re-derives every row against main.

agent: slices 1 and 2 one opus implementer each (effort high) and one opus review (effort high) over both, which
runs the mutants below and checks the fix round: a slip in either fits another model with no message (a bridge
read falling back to its default, a door taking the wrong count). Slices 3 and 4 and the two consumer branches
one sonnet implementer and one sonnet review (effort medium). Tier: slice 1 edits the bridge's creation path, so
it takes the engine tier's review and gates (sanitizers on the R-loaded path) although no engine file moves;
slices 2 to 4 are "R surface, no draw change" ([Process by risk](README.md#process-by-risk)).
rng: NEUTRAL, bit for bit, for every call main accepts with the same values, except the calls listed under
[What changes for a user](#what-changes-for-a-user), each by a ruling or by the rule a ruling states: a control
copied from a fit no longer sets the tree count (dec-B334), a sampler assembled by `new()` takes the model's count
(dec-B335), a data object's own cut count is used by `dbarts()`, and a control's stated cut count beside a data
object that carries a grid re-derives it, as 0.9-34 did. No engine file changes; nothing is re-recorded.
window: for 1.0-0; after tree-prior-flat; slice 1 before forest-door (dec-B518) is built
([Order of landing](#order-of-landing)).
budget: ~5600 lines in the tree. Slice 1 ~1900 (R ~650, bridge ~250, tinytest ~900, docs ~100). Slice 2 ~1250
(R ~520, tinytest ~730). Slice 3 ~2100 changed lines, mechanical (tinytest, inst/common, benchmarks, man
examples, vignettes). Slice 4 ~350 (help ~200, NEWS ~40, design docs and records ~110). Consumers: stan4bart ~60,
bartCause ~20, treatSens and bairrtt none. Out of the tree ~500 lines of scratch: the A/B script and the
creation-trace driver. Agent time, estimated from comparable landings, not measured: slice 1 ~14 agent-hours,
slice 2 ~9, review of both with mutants and the fix round ~6, slice 3 ~5 with its gate, slice 4 ~3, consumers
~2. Plans here run 1.5 to 2 times low: budget 60 to 75.

## Summary for the maintainer

You ruled that the control holds run settings only: the facts a fit's control carries today on hidden
attributes (the variance forest, the forests and their settings, the ordinal category count, the
negative-binomial shape, the hazard periods, the survival status) move to the model and data objects; the tree
count becomes each forest's, with a plain `n.trees` on `dbarts()` and `dbartsSpec()`; the cut count and the
quantile rule become `n.cuts` and `useQuantiles` on `dbarts()` and `dbartsData()`; the binary flag is the
family's. The control's `n.trees`, `useQuantiles` and `binary` stay as copies the sampler writes, for code that
reads them, until no CRAN package does; its `n.cuts` goes. A control taken from a fit sets no tree count for the
next fit, and a warning says when its count differs and is not used; the same at `new("dbartsSampler", ...)`,
where the model's count wins, and at `$setControl`, which keeps the sampler's count. This plan builds all of it.

What a user sees. `dbarts(y ~ ., d, n.trees = 50, n.cuts = 20)` works, as `bart()` always has, and
`forests = forest(n.trees = 50)` works without `list()`. A 0.9-34 script that writes
`dbartsControl(n.trees = 50)` fits 50 trees and gets one warning a session saying where the setting now lives;
the same for `n.cuts` and `useQuantiles`. Code that reads `fit$control@n.trees` or `@binary` reads the right value
with no warning (R cannot warn on a slot read). Passing one fit's control to another fit, to `bart()`, `xbart()`
or `dbarts()` alike, carries its run settings and nothing else; today `bart()` refuses such a control from a
multi-forest fit and `dbarts()` quietly takes the old fit's tree count. A model object now describes its whole
model, so `$setModel` refuses one of another structure by name where today the structure silently survives on
the control.

Nothing a user fits with the same values draws differently, except the copied-control cases you ruled on and two
cut-count cases where `dbarts()` today drops what the caller stated: a data object built with its own cut count,
and a control stating one beside a data object taken from a fit (0.9-34 honored the second; main does not).

One call is yours: what a control whose tree count the caller edited by hand does at `new("dbartsSampler",
...)`. 0.9-34 used the edited count; your ruling for a copied count is that the model's wins with a warning. I
recommend the same for an edited one.

The consumers: stan4bart routes its `bart_args` count and quantile rule to their new homes (about 60 lines);
bartCause's `bcf()` stops sending the cut settings through the control (about 20); treatSens and bairrtt already
use the new spellings and need nothing.

## Inventory

Counts on main c55a1477, from grep and a balanced-call scan (`dbartsControl(` calls whose arguments name the
setting); a comment that names a quantity counts. Probes ran on a private build of main and on CRAN's 0.9-34.

The control's hidden attributes, written by [`resolveSamplerSpec`](../../R/spec.R) and read by the bridge or R:
`bartcore.variance`, `bartcore.forests` (params, sd, vars, labels, heldAmplitude, interactions, blocks,
basisTerms, and `control.n.trees`, the carried count of [`resolveSamplerSpec`](../../R/spec.R)),
`bartcore.forestsDeclared`, `bartcore.n.categories`, `bartcore.shape`, `bartcore.hazard.periods`,
`bartcore.survival`, `bartcore.moves` (per-forest move settings, added by per-forest-move-settings after the
2026-10-05 design). Lines naming one: R 44, src 17, tinytest 80. Bridge readers: [`parseControl`](../../src/R_interface_bartcore.cpp)
(categories, shape, the moves), [`applySurvivalAttribute`](../../src/R_interface_bartcore.cpp),
[`applyVarianceAttributes`](../../src/R_interface_bartcore.cpp), [`applyForestAttributes`](../../src/R_interface_bartcore.cpp)
and the forest-count read of [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp)'s callers. The
calibration anchor already lives on the data object (`forest.calibration`).

The control's slots: `@n.trees` R 21, src 2 (the bridge's `numTrees`), tinytest 28, man 1; `@binary` R 22
(the bridge's `responseIsBinary`; xbart, rbart, partial dependence), tinytest 8; `@useQuantiles` R 3 (the bridge's
`useQuantiles`, the xbart data handle), tinytest 5; `@n.cuts` on a control R 4 (copied to `data@n.cuts` by
`dbarts()`, `dbartsSpec()` and `xbart()`), never the bridge. Provenance records: `dbarts.supplied` (the names a
`dbartsControl()` call carried, read by [`controlSuppliedSlots`](../../R/dbarts.R)) and `dbarts.stated` (what
the caller stated of the plain forest's arguments, [`plainForestStated`](../../R/model.R)).

Constructor calls: of 749 `dbartsControl(` calls in R/, inst/, benchmarks/, man/ and vignettes/, those naming
`n.trees`, `n.cuts`, `useQuantiles` number 400, 19, 10 in tinytest; 59, 24, 7 in benchmarks; 15, 1, 1 in man;
7, 0, 0 in vignettes; 2 in inst/common; one call in R/ (`bartBT`, which names all three), in 228 files. `bart()`
builds its control through a redirected call that names all three on every fit.

Doors on main (probes p03, p05):

| call | main | 0.9-34 |
|---|---|---|
| `dbarts(..., n.trees = 30)`, `dbartsSpec(..., n.trees = 30)` | unused argument | unused argument |
| `dbarts(..., control = dbartsControl(n.trees = 30))` | 30 | 30 |
| `forests = forest(n.trees = 30)` without `list()` | refused | - |
| `y ~ forest(x1 + x2, n.trees = 30)` with the control's 40 | 30, silent | - |
| `dbarts(..., control = fit50$control)` | 50, the old fit's | error |
| `dbarts(...)` given a `bart()` fit's control, every forest with a basis | refused as naming `n.trees` | - |
| `new("dbartsSampler", fit50$control, fit75$model, fit75$data)` | 50 | the control's |
| `fit75$setControl(fit50$control)`; `fit200$setControl(dbartsControl())` | refused | resized |
| `dbarts(..., n.cuts = 10)`; `dbartsData(..., n.cuts = 10)` | unused argument | unused argument |
| `dbarts(dd)`, `dd@n.cuts` set to 10 by hand | 100, the control's default | 100 |
| `dbartsSpec(dd, dbartsControl(n.cuts = 20))` | 10, the control's dropped | - |
| `dbarts(fit$data, control = dbartsControl(n.cuts = 20))` | the fit's grid, the 20 dropped | 20 |
| a probit fit's control into a gaussian `dbarts()` | gaussian | error |
| one forest declared by `forests =` or a lone forest term: `getTrees()` | a `forest` column | - |

The last row is a defect against dec-B249 ("a single-forest fit returns 0.9-34's columns"), kept alive by
`bartcore.forestsDeclared`; it goes with the attribute.

Consumers, read at their main branches (stan4bart 566a1d9, bartCause f7edf5b, treatSens cd59d0f, bairrtt
6575ba0, all clean):

- stan4bart R/stan4bart_fit.R forwards every `bart_args` name that is a `dbartsControl` formal into the control,
  `n.trees` and `useQuantiles` included, then every `dbartsSpec` formal into `dbartsSpec`; `n.cuts` is already
  written on the data object. R/mvbart.R builds its control from `bart_args` by `dbartsControl`'s formals and
  refuses `n.trees` there. R/generics.R rebuilds a sampler from a fit's stored triple through `new()`.
- bartCause R/bcf.R routes `...` names that are `dbartsControl` formals into the control (`n.cuts` and
  `useQuantiles` reach the sampler that way) and states each forest's count on `forest()`. Its tests read
  `fit$control@n.trees` of a `bart()` fit (test-21). R/responseFit.R and R/treatmentFit.R filter by the formals of
  `bart` and `dbartsControl`.
- treatSens R/cibart.R builds the triple through `dbartsSpec(forests = list(forest(n.trees = )))` and writes
  `data@n.cuts`; R/treatSensBART.R's control names none of the three.
- bairrtt passes `forests = list(forest(n.trees = ))` and a control naming none.
- CRAN (dec-B250's scan of 27 reverse dependencies): bartMan and tidypredict read `fit$control@binary`;
  countSTAR, mcmcsae and nlfh pass `dbartsControl(n.trees = )` to `dbarts`; no package reads `@n.trees` or
  `@useQuantiles`.

## Design

### Homes

| quantity | home after | stated by | fixed at creation |
|---|---|---|---|
| tree count, per mean forest | the model's forest record | `forest(n.trees = )`, the door's `n.trees` for the forests with no basis, the defaults by kind | yes (dec-B247) |
| the rest of a forest's declaration: vars, sd, base, power, label, held amplitude, interactions, blocks, basis term, move settings | the same record | `forest()`, the door's level | as today: `sd` and the priors through the live doors, the rest structure |
| the variance forest | the model's variance record (count, base, power, columns, moves) | `varianceForest()`, `variance =` | yes |
| ordinal category count, negative-binomial shape setting, hazard periods | the model | the family | yes; the held shape a parameter (setmodel-parameters) |
| survival status | the data object, one value per row | the response, `setResponse(status = )` | no |
| cut count | the data object (`n.cuts`, exists) | `dbarts(n.cuts = )`, `dbartsData(n.cuts = )` | the grid is data (dec-B231) |
| quantile rule | the data object (`useQuantiles`, new) | `dbarts(useQuantiles = )`, `dbartsData(useQuantiles = )` | as the cut count |
| binary | none: the model's family | the family | yes |

New homes are typed S4 slots with validity, set after `new()` where `initialize` is closed. Every model built by a
door carries one forest record per mean forest, a single-forest model one and a multinomial model one that every
category forest takes. The record class has an unstated count (`NULL`) as its prototype, so a model built by
`new("dbartsModel")` by hand states none ([Every door](#every-door), the direct door). No reader guards the new
slots against an object saved by a development build: only shipped versions are owed compatibility, and a fit
saved by 0.9-34 fails as today (TODO legacy-object-message, dec-B461).

The control keeps `binary`, `n.trees` and `useQuantiles` as mirrored slots and loses the `n.cuts` slot. A
`dbartsControl(n.cuts = )` value, the tombstoned spelling, rides an attribute named `n.cuts`, where 0.9-34 kept
it: `ctl@n.cuts` reads back what the caller wrote, as in 0.9-34, and a sampler's stored control carries none, as
in 0.9-34.

### Stated, carried and silent

A control's value for `n.trees`, `useQuantiles` or `n.cuts` is one of three:

- written: the caller's `dbartsControl()` call named it (`dbarts.supplied`), or its slot differs from the value
  the sampler last mirrored into it, or, on a control with no mirror record, from a fresh control's;
- carried: the control holds a mirror record (`dbarts.mirror`, written with the slots) and the slot equals it;
- silent: anything else, a fresh control's default.

A written value is the relocated argument of dec-A02: it warns once a session (`dbartsDeprecatedWarning`, one key
per name shared by `dbartsControl()` and the doors, since an edited slot or a control read from a file reaches a
door with no constructor call) and is read exactly as the door's flat argument would be. A carried value is never
a statement (dec-A158, dec-B334). A silent value is never read. A written value equal to the mirror cannot be told
from no edit and is carried; the warning below then says so whenever it differs from what is taken, which closes
control-carried-count-edges' last case loudly. `binary` is never read as a statement at a front door: 0.9-34 had
no argument for it and overwrote the slot from the response.

The not-used warning (`dbartsWarning`, every call, not once a session): where a written or carried value on a
control differs from what the call takes, "the control's n.trees (50) is not used: this fit has 75 trees in its
forest with no basis; state the count as dbarts(n.trees = ) or forest(n.trees = )", worded per setting and per
door. That is dec-B334's warning; dec-B335 and dec-B385 name the same one.

### The rule at each door

At a door with flat arguments the levels are, for the forests with no basis: a forest's own `n.trees`; the door's
`n.trees`, or a written control value standing for it; the defaults by kind (dec-B274). For the cut settings: the
door's `n.cuts` and `useQuantiles`, or a written control value standing for them; what the data object carries
(a grid, else its own count and rule); the defaults (100, uniform).

1. A flat argument and a written control value in one call: the flat one is used (dec-B116), the control's draws
   the not-used warning where it differs.
2. The door's count and a forest's own: refused as today where every forest with no basis states one (dec-B241,
   dec-B503), and the same for a written control count, the control's being the same level ("bart's n.trees or
   the control's", man/forest.Rd). Today the forest's governs the control's in silence.
3. A count, flat or written, where every forest has a basis: refused by name (dec-B274). A carried one is not a
   statement and is never refused, which retires control-carried-count-edges' first case.
4. An explicit cut setting, flat or written, beside a data object that carries a grid: the grid is re-derived at
   the stated setting (dec-B231: "If explicit arguments are given, we could then follow them directly"), as
   0.9-34 re-cut there. A data object stating its own count or rule with no grid: used, at `dbarts()` as at
   `dbartsSpec()` (dec-B241: the cut settings describe the predictors).
5. The mirrors are written last, from the model and data the sampler takes.

#### Every door

| door | tree count | cut count, quantile rule |
|---|---|---|
| `dbarts()` | `n.trees = NULL` after `forests`; the rule above | `n.cuts = NULL`, `useQuantiles = NULL` after `factors`; the rule above |
| `dbartsSpec()` | `n.trees = NULL` after `forests` (dec-B241's call); the rule | no flat arguments: a written control value is read as `dbarts()` would read it beside the same data object, which `?dbartsSpec` promises; today the control's is dropped |
| `dbartsData()` | - | `n.cuts = NULL`, `useQuantiles = NULL`, last; NA refused, `NULL` states none |
| `bart()`, `pdbart()`, `rbart_vi()` | their `n.trees`, passed flat to `dbarts()`; the control they build names none of the three | their `n.cuts`, `useQuantiles`, passed flat |
| `bartBT()` | `ntree`, passed flat where its caller named it, as dec-A179 records `bartBT` doing for the tree prior | `numcut`, `usequants`, passed flat |
| `xbart()` | its grid; a written control count is a grid of one; each cell's count goes on the cell's model record | flat; the shared data handle takes the data object's rule |
| a formula's `forest()` terms, `forests = ` | per forest; a bare `forest()` is one forest (dec-B241's call) | - |
| a data object as `data` | the rule; a data object carries no count | the rule, row 4 |
| `new("dbartsSampler", control, model, data)` | the model's count (dec-B335); a differing control count draws the not-used warning; a model stating none takes a written control count with the deprecation warning, else the default (open call 1 covers a written count beside a model that states one) | the data object's; a differing carried or written control value draws the not-used warning |
| re-creation (`getPointer`, `setState`, `copy`, a reload, stan4bart's restore) | the stored model's; the mirror agrees, so no warning | the stored data's |
| `$setControl(new)` | a carried count that differs: taken, the sampler's kept, the not-used warning (dec-B385); a written one that differs: refused, "changing 'n.trees' is not available on an existing sampler" (dec-B247) | the same two cases for each, the sampler's data holding them; the `n.cuts` attribute is stripped after the check |
| `$setModel(model)` | a model whose forest records, variance record, category count or hazard periods differ is refused by name before anything is stored (dec-B254); the moved facts join setmodel-parameters' structure check | - |
| `$setData(data)` | - | `useQuantiles` joins `n.cuts` and `sigma` among what the sampler keeps from its own data, as `setData` documents for `n.cuts` |
| `$setForest()`, `$getForest()` (dec-B518, unbuilt) | `n.trees` is structure, refused by name at `$setForest`; `$getForest` reads it from the record, so giving it back changes nothing | neither: the cut settings are data, changed through `setCutPoints` and `setPredictor(updateCutPoints = )` |

`binary`: the model's family decides at every door. At the direct door a control whose `binary` differs from the
model's family draws the not-used warning, where 0.9-34 took the control's (dec-B335's rule, applied to the flag).
At the front doors a control's `binary` is never read and never warned about.

`bart()` and `xbart()` lose [`refuseFitStateControl`](../../R/dbarts.R)'s attribute refusal: no control carries a
fit's model any more, so every entry point takes any control (dec-B239). [`mergeFrontDoorControl`](../../R/dbarts.R)
keeps its rule for run settings and stops covering the three.

#### Checked against scripts

| script | 0.9-34 | main | after |
|---|---|---|---|
| `ctl <- dbartsControl(n.trees = 50); dbarts(..., control = ctl)`, twice | 50, 50 | 50, 50 | 50, 50; one warning a session |
| `ctl@n.trees <- 200L` on a fresh control, then `dbarts()` | 200 | 200 | 200; the deprecation warning |
| `dbarts(..., control = fit50$control)` | error | 50 | 75; the not-used warning |
| `cc <- fit$control; cc@n.trees <- 30L; dbarts(..., control = cc)` | error | 30 | 30; the deprecation warning |
| `bart(..., n.trees = 12, control = dbartsControl(n.trees = 9))` | - | 12 | 12; both warnings |
| `y ~ forest(x, n.trees = 30)` with `dbartsControl(n.trees = 50)` | - | 30 | refused, rule 2 |
| a `bart()` fit's control into a `dbarts()` whose every forest has a basis | - | refused | taken; 50 and 50 by kind |
| `new("dbartsSampler", fit200$control, fit75$model, fit75$data)` | 200 | 200 | 75; the not-used warning |
| `fit75$setControl(fit200$control)` | resized to 200 | refused | taken, 75 kept; the not-used warning |
| `fit75$setControl(dbartsControl(n.trees = 20))` | resized to 20 | refused | refused |
| `fit200$setControl(dbartsControl(n.samples = 100))` | resized to 75 | refused | taken, 200 kept, silent |
| `dbarts(dbartsData(y ~ ., d, n.cuts = 10))` | - | - | 10 |
| `dbarts(fit$data, control = dbartsControl(n.cuts = 20))` | 20 | the fit's grid | 20; the deprecation warning |
| `dbarts(fit$data, n.cuts = 20)` | - | - | 20 |
| `dbarts(fit$data, control = dbartsControl(useQuantiles = TRUE))` | quantiles | the fit's grid | quantiles; the deprecation warning |
| `dbartsSpec(dd, dbartsControl(n.cuts = 20))`, `dd` stating 10 | - | 10 | 20; the deprecation warning |
| a variance fit's control into a plain `dbarts()` or `bart()` | - | plain; `bart()` refuses | plain at both |
| `fit$control@n.trees`, `@binary`, `@useQuantiles` after any fit | the fit's | the fit's (first forest's count) | the same values |

### Mirror mechanics

One way, from the model and data to the control, in one helper called at the end of
[`resolveSamplerSpec`](../../R/spec.R) (so `dbartsSpec()`'s returned control reads right), when a sampler stores
its control, and in `$setControl`: `n.trees` the first forest's count (as man/forest.Rd states today),
`useQuantiles` the data object's rule, `binary` whether the family is binary (probit and logistic, a hazard fit
being one of them); the mirror record written beside them, the three names taken out of `dbarts.supplied`, the
`n.cuts` attribute stripped. The bridge reads none of the three: after slice 1 it reads the model and the data
alone for every model fact.

Expiry. `dbartsControl()`'s three arguments are tombstones that expire with the others at 1.1-0 (dec-A02). The
slots stay until no CRAN package reads them, checked at the 1.1-0 freeze (dec-B250): on today's landscape
`n.trees` and `useQuantiles` go at 1.1-0 and `binary` waits on bartMan and tidypredict. Registry rows in
[`dbartsTombstones`](../../R/tombstones.R): the three arguments; "a count, cut count or quantile rule read from a
control"; the three slots, each with dec-B250's condition. Deleted with the arguments: the written-value reading,
the deprecation key, the `n.cuts` attribute. Deleted with each slot: its mirror write, and the not-used warning's
case for it.

### The bridge

[`parseControl`](../../src/R_interface_bartcore.cpp) stops reading `binary`, `n.trees`, `useQuantiles`, the
category count, the shape and the moves; [`parseModel`](../../src/R_interface_bartcore.cpp) reads the forest
records (the first forest's count among them), the variance record, the category count, the shape and each
forest's moves; [`parseData`](../../src/R_interface_bartcore.cpp) reads `useQuantiles` and the status. The binary
flag the bridge uses (the summary, `sigmaIsFixed`, [`resolveFamily`](../../src/R_interface_bartcore.cpp)) comes
from the parsed family. [`applyVarianceAttributes`](../../src/R_interface_bartcore.cpp),
[`applySurvivalAttribute`](../../src/R_interface_bartcore.cpp) and [`applyForestAttributes`](../../src/R_interface_bartcore.cpp)
read the model and data and are renamed for what they read. [`bartcore_setControl`](../../src/R_interface_bartcore.cpp)
loses its tree-count refusal (the control has none to compare); `bartcore_setModel` gains it as a backstop for
callers that skip R's check (xbart's cell install). The data handle takes `useQuantiles` from the data object and
`sparseDensityThreshold` from the control as now. A slot holding NULL reaches the bridge as R's S4 NULL marker;
every new NULL-able slot is read through the normalizing helper `parseData` uses for `bases`. Unchanged:
src/bartcore, the facade, `SamplerOptions` and what goes into it, the order of construction and generator
creation, inst/include/dbarts/dbarts.h.

### What changes for a user

Refused where main runs: a written control count beside a forest that states its own (rule 2); a written
`useQuantiles` or `n.cuts` on `$setControl` that differs from the sampler's (main refuses `useQuantiles` already;
`n.cuts` is stored there and never reaches the engine); `$setModel` with a model of another structure in a fact
that today rode the control.

Taken where main refuses: `forests = forest()`; `dbarts(n.trees = )`, `dbartsSpec(n.trees = )`,
`dbarts(n.cuts = , useQuantiles = )`, `dbartsData(n.cuts = , useQuantiles = )`; any control at `bart()` and
`xbart()`; a carried control count at `$setControl` and beside forests that all have a basis.

Different with no error: a copied control's count (dec-B334) and the direct door (dec-B335), each with the
not-used warning; `dbarts()` given a data object stating its own cut count or rule; a written cut setting beside a
data object carrying a grid (0.9-34's result restored); a single-forest fit's `getTrees()` loses the `forest`
column however the forest was written (dec-B249); `dbartsControl()`'s three arguments warn.

## Slices

### Slice 1: the records (R and bridge)

The model's forest records and variance record, its category count, shape and hazard periods; the data's status
and `useQuantiles`; the mirror helper and record; the bridge reads of [The bridge](#the-bridge);
[`resolveSamplerSpec`](../../R/spec.R) writes the records in place of the attributes, from what main resolves;
the readers move: [`fitTreeCounts`](../../R/generics.R), [`forestPriorLines`](../../R/forestPrint.R), the fit
packagers in R/bart.R, the hazard periods of partial dependence, `getSigmas`, `samplePriorPredictive`,
`getTrees` (its column by dec-B249), [`resolveForestSpreads`](../../R/dbarts.R) and
[`writeLeafPrior`](../../R/dbarts.R), `setResponse(status = )`, xbart's cell models and handle; `$setControl`'s
attribute loops and carried-count branch go, with dec-B385's case; `$setModel`'s structure check gains the moved
facts; [`refuseFitStateControl`](../../R/dbarts.R)'s attribute test goes. The carried count dies with
`bartcore.forests`: dec-B334 and dec-B335 land here, the doors otherwise reading the count where main reads it
(the slot, the plain forest's statement). ~1900 lines: R ~650, bridge ~250, tinytest ~900 (the 80 attribute
lines follow the objects; the read-backs and the door rows of dec-B334, dec-B335 and dec-B385; a new
test-control-migration.R), docs ~100.

### Slice 2: the doors (R)

`n.trees` on `dbarts()` and `dbartsSpec()`, `n.cuts` and `useQuantiles` on `dbarts()` and `dbartsData()`, the
bare `forest()`; written, carried and silent with both warnings at every door of
[Every door](#every-door); `bart()`, `bartBT()`, `rbart_vi()`, `xbart()` and `pdbart()` pass the three flat and
build a control naming none, so control-carried-count-edges' first case cannot arise; `dbartsControl()`'s three
arguments default to `NULL` and warn; the `n.cuts` attribute; the deprecation key in
[`warnOnce`](../../R/utility.R)'s scheme; the registry rows. ~1250 lines: R ~520, tinytest ~730 (every row of
[Checked against scripts](#checked-against-scripts) at every door that has it; each way a value is written,
carried or silent; a refusal leaving the sampler untouched; the mirrors after every family's fit; the warnings
counted with `withCallingHandlers`).

### Slice 3: respell (tests, benchmarks, examples)

The ~520 constructor calls naming the three, their later uses and the wrapper helpers that forward to
`dbartsControl()`, moved to the flat arguments or `forest(n.trees = )`; man examples (8 pages) and vignettes (3).
No R/ or src/ change. Tests that pin a 0.9-x spelling on purpose keep it and reset the warning key. ~2100 changed
lines. Its gate is mechanical ([Tests and gates](#tests-and-gates)).

### Slice 4: help, NEWS, docs; consumers

Help: dbarts.Rd, dbartsSpec.Rd, dbartsData.Rd, dbartsControl.Rd (the three arguments as tombstones, the slots as
copies), forest.Rd (the first-forest paragraph: the control's count is the fitting function's level and is
refused beside a forest's), dbartsForests.Rd, dbartsSampler-class.Rd (`setControl`, `setModel`, the direct door),
bart.Rd, xbart.Rd. NEWS against 0.9-34: the new arguments; the three `dbartsControl()` arguments warn and expire at
1.1-0; the direct door takes the model's count; `$setControl` keeps the sampler's count (dec-B247's item reworded:
a control naming another count is refused, one carrying another is taken with a warning); a control taken from a
fit can be passed to any door. Design: docs/design/feature-matrix.md, the forest-defaults-by-kind design's "Where
the numbers live", docs/architecture.md's control paragraph; records. ~350 lines.

Consumers, each on its own branch, green against this branch's tip before the merge:

- stan4bart: `bart_args$n.trees` to `dbartsSpec(n.trees = )` and `useQuantiles` to the data object, both kept out
  of the control loop and named in its known arguments; R/mvbart.R sends `n.cuts` and `useQuantiles` to
  `dbarts()`. Its restore through `new()` reads agreeing mirrors and stays silent. ~60 lines with tests.
- bartCause: `bcf()` keeps `n.cuts` and `useQuantiles` out of the control and lets its `dbarts()` loop pass them;
  test-21's `@n.trees` read stays right. ~20 lines.
- treatSens, bairrtt: nothing; both suites run against the tip.

## Tests and gates

Slice 1, per moved quantity at a non-default value, read back from the engine, each with a mutant that replaces
the read with its default:

| quantity | read back by |
|---|---|
| the variance forest's count, columns, base, power | the stored state's variance trees; seeded draws against the default's; equal after `copy()` and a reload |
| a fixed negative-binomial shape; the ordinal category count | `getShape()`; the thresholds' length |
| per-forest count, vars, interactions, blocks, moves, three forests | `bartcoreForestTreeCount`, `getTrees(forest = )` split variables, the move install |
| per-forest sd, held amplitude, labels | `getLeafPrior()`, `getForestAmplitudes()`, a label door |
| the status after `setResponse(status = )`; the hazard periods | a reloaded sampler's seeded draws against its unsaved twin; the fit's `$periods` and `pdbart` after a reload |
| the family alone deciding binary | each family's run channels, a hand-built model included |
| `useQuantiles` at each door and under `xbart()` | the state's cut points, quantile against uniform |

Slice 2: the script table and the door table as tests; the deprecation warning once per key per session and the
not-used warning on every call, counted; none for a silent or agreeing value. Slice 3: a creation-trace driver logs
each sampler creation of the full suite and the three equivalence harnesses (per file, in order: each forest's
engine tree count, `data@n.cuts`, the rule), before and after; the tables must be equal, the expectation count
unchanged.

RNG gates, on every slice's tip and on the merged tree: tests/cpp; the full tinytest suite in one process at the
default ncpu; the equivalence trio bit for bit on a reference build (55, 19 and 11 "identical draws" lines, no
`max |z|`); the four snapshot files unchanged; a scratch A/B digesting seeded draws for a matrix of calls (every
door, family, forest layout, cut setting, a copy, a reload, `xbart`) on a c55a1477-plus-tree-prior-flat build and
on the slice, identical except the listed rows, which must match their expected new values. The fit object
carries new slots, so `.github/workflows/exact-gates.yaml`'s gates run with `quick`. Slice 1: sanitizers on
tests/cpp and on the R-loaded path. Every slice: `lintr::lint_package()` on the slice's library, `air format
--check .`, tools/check-rc-codoc.R, tools/check-win-drift.R, tools/check-doc-freshness.R, `R CMD check --as-cran`
from a clean tarball, and the mutation battery's `verify-anchors` (slice 1 touches src/).

## Mutants

The review runs each and requires a failing test:

1. The bridge reads a fixed 75 for the first forest's count.
2. The bridge reads the variance forest's count as 40 whatever the record says.
3. `binary` taken as FALSE for logistic (the parsed family's flag dropped).
4. `useQuantiles` read as FALSE from the data object.
5. The status not read on re-creation (the reload twin differs).
6. A forest's moves dropped from the record (the move install reads the control's).
7. `getTrees()`' column shown for one forest.
8. A carried count treated as written (the copied control's 50 used).
9. A written count treated as carried (the edited 30 dropped).
10. The not-used warning skipped at the direct door.
11. `$setControl` refusing a carried differing count.
12. `$setControl` taking a written differing count.
13. The flat `n.trees` losing to a written control count.
14. A written control count beside a forest's taken in silence.
15. `dbarts()` overwriting a data object's own count with the default.
16. A written cut count beside a carried grid ignored.
17. The mirror written before the forest's statement resolves (the control reads the door's count).
18. `$setModel` accepting a model with another variance forest.
19. The deprecation key not shared (two warnings for one name).
20. xbart's cell count left on the control.

## Calls made

1. A per-forest move setting lives in the forest record: a forest is declared once, and a control passed on must
   need no handling (dec-B239). `$setModel` refuses a record whose moves differ, naming `$setControl` and
   `$setForest`. Alternative: keep `bartcore.moves` on the control.
2. `dbartsControl()` keeps the three formals, defaulting to `NULL`, until 1.1-0, so callers that route by its
   formals (bartCause, stan4bart) keep reaching them and the name resolves (dec-A02). Alternative: move them to
   `...`, which makes such callers drop them.
3. `dbartsControl(n.cuts = )` rides an attribute named `n.cuts`, as in 0.9-34, so `ctl@n.cuts` reads what the
   caller wrote. Alternative: a private name, under which a 0.9-x read errors.
4. A control count written beside a forest's own is refused, the control's being the fitting function's level
   (dec-B241's refusal, dec-B512's one rule). Today the forest's governs; no release had `forest()`.
5. A written value beaten by a flat argument draws the not-used warning where it differs, beside the deprecation
   warning. Alternative: the deprecation warning alone.
6. `useQuantiles` and `n.cuts` follow dec-B334 and dec-B385's rule for the count (a copy never speaks; it warns
   where it differs), as dec-A158 treated the three alike.
7. `binary` is never a statement at a front door; at the direct door the model's family wins with the not-used
   warning (dec-B335's rule).
8. The mirrored `n.trees` is the first forest's count, as man/forest.Rd says today. Alternative: the forest with
   no basis's, which a model with none lacks.
9. At `dbartsSpec()` a written control cut setting is read as at `dbarts()` beside the same data object, which
   `?dbartsSpec` promises; today it is dropped. Alternative: refuse it.
10. A model with no stated count, only buildable by `new("dbartsModel")`, takes a written control count at the
    direct door with the deprecation warning, else the default.
11. No guarded accessors for development-build objects (compatibility is owed to shipped versions only).
12. dec-B334 and dec-B335 land in slice 1, where the carried count's attribute dies, rather than keeping it alive
    a slice longer.
13. The flat arguments default to `NULL`, meaning not stated, the help giving the values that then apply (75 by
    kind, 100, `FALSE`).
14. Slice order: records before doors, so the doors are written once against the records.

## Open calls

1. A tree count edited by hand on a control at `new("dbartsSampler", control, model, data)`. Background: a
   sampler can be assembled from a control, a model and a data object. In 0.9-34 the control's count decided, and
   a script could edit `ctl@n.trees <- 20L` and build a 20-tree sampler that way. You ruled (dec-B335) that a
   control copied from a fit loses to the model's count there, with a warning; a count the caller edited is not
   named, and every model now states its count. Options: (a) the model's count, with the not-used warning naming
   where to set it (`dbartsSpec(n.trees = )`); one rule at the door, a 0.9-x script of this shape gets the
   model's count and is told; (b) the control's count, with the deprecation warning, as 0.9-34 and as at the
   front doors; the direct door then has two rules, by provenance, and a model that states 75 is overruled by an
   object that is not the model; (c) refused by name; loud, but a script that worked in 0.9-34 stops. I
   recommend (a).

## Order of landing

1. tree-prior-flat lands first. This plan is written against it: its two-level rule for `power` and `base` is
   the rule this plan gives the count (rule 2), and its edits to [`resolveSamplerSpec`](../../R/spec.R),
   [`plainForestStated`](../../R/model.R), [`forestParams`](../../R/model.R), R/bart.R, R/xbart.R, R/rbart.R,
   R/tombstones.R, test-front-door-control.R, test-bcf-creation.R and man/forest.Rd are textual collisions;
   this branch branches from main after it.
2. forests-shared-total (slices B and 1 built, slice 2 unbuilt): its slice 2's [`forestParams`](../../R/model.R)
   split of the count among sharers writes into the forest records if slice 1 here lands first, and into the
   params vector otherwise, one site either way. Its R item 2's last sentence, that the forests' record keeps the
   count the next fit inherits where forest 1 is not plain, is the carried count dec-B334 retires: that sentence
   goes from its plan whichever lands first. Its bridge edits (`createHolder`, `printInitialSummary`,
   `bartcore_getLeafPrior`) meet slice 1's in src/R_interface_bartcore.cpp: whichever is second rebases. Engine
   features merge serially; slice 1 edits no engine file.
3. setmodel-parameters (slices 1 to 3 built, slice 4 after forests-shared-total): its structure check
   (`modelStructureRefusal`) reads model attributes, and slice 1 here adds the moved facts to it; its `$setModel`
   writes the negative-binomial shape onto `bartcore.shape` for a copy to read, which slice 1 moves to the
   model's shape slot. Whichever lands second carries both edits.
4. forest-after-merge (slice 2 built): `setForestBasis` names a forest from `bartcore.forests`' labels, which
   slice 1 moves to the record. Whichever is second moves the read.
5. forest-door (dec-B518, unbuilt): built after slice 1, so `$getForest` reads the record and `$setForest` refuses
   `n.trees` from it.
6. Consumers: stan4bart's and bartCause's branches green against this branch's tip, merged in the same sitting as
   dbarts; between the merges stan4bart's `bart_args$useQuantiles` warns once (it still reaches the control).

## Not verified

- Nothing is prototyped: the slots, the records, the mirror record, the warnings, the bridge reads.
- The slice sizes and agent-hours; the respell count is from a balanced-call scan of `dbartsControl(` and misses
  calls built by `do.call` or aliases.
- That the creation-trace driver can trace the reference class's `initialize` from outside the package (the
  2026-10-05 critique ran one on an older tree).
- That xbart's cell install, the one caller of `bartcore_setModel` that skips R's check, needs only the backstop.
- The consumers and the CRAN packages were read, not run; dec-B250's CRAN scan was not repeated.
- A control saved by 0.9-34 fails `validObject` on main (slots added since); this plan does not change that, and
  TODO legacy-object-message may want to name it.
- `rbart_vi`'s path was not probed (its fit object did not expose a sampler in the probe).
