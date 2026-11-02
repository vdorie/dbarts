# control-migration: a control holds run settings only

Status: PLANNED 2026-10-30 on c55a1477, against main as it stood once coefficient-prior-opt-in (dec-A245) and
tree-prior-flat (dec-A246) landed; amended the same day after its blind critique (scratch/cm-crit/report.md, build
after amendments) and the coordinator's calls on it. Slices 1a, 1b and 2 built on the branch and reviewed, the review's fix batch in; slices 3 and 4 and the consumer branches built, not reviewed. Both calls
held for the maintainer are ruled: dec-B526, dec-B527 ([Open calls](#open-calls)). Rulings: dec-B239, dec-B240,
dec-B241 (the move), dec-B250 (when the mirrored slots go), dec-B334, dec-B335, dec-B385 (a control copied from a
fit), dec-B247 (a live sampler's count), dec-B249 (the forest column), dec-B503 and dec-B505 (counts among
forests with no basis), dec-B518 (the live forest door), with dec-A02 (relocated arguments), dec-B116 (a flat
argument beats the control's), dec-B231 (split points are data), dec-B254 (structure is refused by setModel),
dec-B274 (defaults by kind) and dec-B512 (one rule at two levels). Earlier work: a design of 2026-10-05 and two
critiques of it, outside the tree; this plan re-derives every row against main.

agent: slices 1a, 1b and 2 one opus implementer each (effort high) and one opus review (effort high) over the
three, which runs the mutants below and checks the fix round: a slip in any fits another model with no message
(a bridge read falling back to its default, a door taking the wrong count). Slices 3 and 4 and the two consumer
branches one sonnet implementer and one sonnet review (effort medium). Tier: slice 1a edits the bridge's creation
path, so it takes the engine tier's review and gates (sanitizers on the R-loaded path) although no engine file
moves; slices 1b to 4 are "R surface, no draw change" ([Process by risk](README.md#process-by-risk)).
rng: slice 1a NEUTRAL, bit for bit, for every call main accepts, except a control from one fit beside another
fit's model and data at the direct door, which follows the model rather than the control's hidden facts
([What changes for a user](#what-changes-for-a-user)): it moves the records and changes no other reading. Slices 1b
and 2 NEUTRAL, bit for bit, for every call main accepts with the same
values, except the calls listed under [What changes for a user](#what-changes-for-a-user), each by a ruling, by
the rule a ruling states: a control copied from a fit no longer sets the tree count or the cut settings
(dec-B334, call 5), the direct door's count and quantile rule (dec-B335, dec-B526), a data object's own cut count
used by `dbarts()`, a control's stated cut setting beside a data object that carries a grid re-deriving it, as
0.9-34 did, and at the direct door what a `dbartsSpec()` data object or a fit's model holds beside a fresh control.
No engine file changes;
nothing is re-recorded.
window: for 1.0-0; after coefficient-prior-opt-in and tree-prior-flat; forest-door (dec-B518) stacks on slice
1a's tip and lands after this branch ([Order of landing](#order-of-landing)).
budget: ~5750 lines in the tree. Slice 1a ~1450 (R ~500, bridge ~250, tinytest ~620, docs ~80). Slice 1b ~600
(R ~200, tinytest ~380, docs ~20). Slice 2 ~1350 (R ~560, tinytest ~790). Slice 3 ~2000 changed lines,
mechanical (tinytest, inst/common, benchmarks, man examples, vignettes). Slice 4 ~350 (help ~200, NEWS ~40,
design docs and records ~110). Consumers: stan4bart ~60, bartCause ~20, treatSens and bairrtt none. Out of the
tree ~500 lines of scratch: the A/B script and the creation-trace driver. Agent time, estimated from comparable
landings, not measured: slice 1a ~12 agent-hours and ~3 more for the rebases onto the branches that land first
(coefficient-prior-opt-in, tree-prior-flat, likely forests-shared-total's slice 2 and setmodel-parameters),
slice 1b ~5, slice 2 ~12, review of the three with mutants and the fix round ~7, slice 3 ~5 with its gate,
slice 4 ~3, consumers ~2: ~49. Plans here run 1.5 to 2 times low: budget 75 to 100.

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

Nothing a user fits with the same values draws differently, except:
- the copied-control cases you ruled on (a fit's control passed to a new fit, to `new()` or to `$setControl`);
- two cut-count cases where `dbarts()` today drops what the caller stated: a data object built with its own cut
  count, and a control stating one beside a data object taken from a fit (0.9-34 honored the second; main does
  not);
- at `new("dbartsSampler", ...)`, a sampler assembled from a `dbartsSpec()` result whose data object asked for
  quantile cuts, which main builds with uniform cuts when the control is fresh;
- at `new("dbartsSampler", ...)`, a tree count or quantile rule written on the control, as dec-B526 rules: it fills
  what the model and data object leave unstated, with the one-release warning, and is refused where they state
  another;
- at `new("dbartsSampler", ...)`, a control from one fit beside another fit's model and data, which follows the
  model.

You ruled the two calls held here. A tree count or quantile rule written on a control at
`new("dbartsSampler", control, model, data)` is used where the model or data object states none and refused by name
where it states another (dec-B526). A count written on a control beside a forest's own `n.trees` is refused by name,
as the plain `n.trees` beside a forest's is (dec-B527).

The consumers: stan4bart routes its `bart_args` count and quantile rule to their new homes (about 60 lines);
bartCause's `bcf()` stops sending the cut settings through the control (about 20); treatSens and bairrtt already
use the new spellings and need nothing.

## Inventory

Counts on slice 1a's tip 60b960ed, from grep and a balanced-call scan (`dbartsControl(` calls whose arguments name
the setting); a comment that names a quantity counts. Probes ran on a private build of main and on CRAN's 0.9-34.

The control's hidden attributes, written by [`resolveSamplerSpec`](../../R/spec.R) and read by the bridge or R:
`bartcore.variance`, `bartcore.forests` (params, sd, vars, labels, heldAmplitude, interactions, blocks,
basisTerms, and `control.n.trees`, the carried count of [`resolveSamplerSpec`](../../R/spec.R)),
`bartcore.forestsDeclared`, `bartcore.n.categories`, `bartcore.shape`, `bartcore.hazard.periods`,
`bartcore.survival`, `bartcore.moves` (per-forest move settings, added by per-forest-move-settings after the
2026-10-05 design); coefficient-prior-opt-in (dec-A245) added `amplitude.df` to `bartcore.forests`, read by the
bridge. Slice 1a moved every one; 144 lines named them on c55a1477. Bridge readers: [`parseControl`](../../src/R_interface_bartcore.cpp)
(categories, shape, the moves), `applySurvivalAttribute`,
`applyVarianceAttributes`, `applyForestAttributes` (since slice 1a [`applySurvivalStatus`](../../src/R_interface_bartcore.cpp), [`applyVarianceRecord`](../../src/R_interface_bartcore.cpp), [`applyForestRecords`](../../src/R_interface_bartcore.cpp))
and the forest-count read of [`bartcore_createDataHandle`](../../src/R_interface_bartcore.cpp)'s callers. The
calibration anchor already lives on the data object (`forest.calibration`).

The control's slots, lines naming the slot on any object: `@n.trees` R 39, tinytest 35, man 1; `@binary` R 26
(xbart, rbart, partial dependence), tinytest 11; `@useQuantiles` R 15, tinytest 13; `@n.cuts` R 23 (a control's is
copied to `data@n.cuts` by `dbarts()`, `dbartsSpec()` and `xbart()`). Since slice 1a the bridge reads none on a
control. Provenance records: `dbarts.supplied` (the names a
`dbartsControl()` call carried, read by [`controlSuppliedSlots`](../../R/dbarts.R)) and `dbarts.stated` (what
the caller stated of the plain forest's arguments, [`plainForestStated`](../../R/model.R)).

Constructor calls: of 753 `dbartsControl(` calls in R/, inst/, benchmarks/, man/ and vignettes/, those naming
`n.trees`, `n.cuts`, `useQuantiles` number 401, 19, 10 in tinytest; 59, 25, 7 in benchmarks; 15, 1, 1 in man;
7, 0, 0 in vignettes; 2 in inst/common; one call in R/ (`bartBT`, which names all three), in 229 files. `bart()`
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
| `dbartsSpec(dd, dbartsControl(n.cuts = 20))` | 10, the control's dropped | - |
| `dbarts(fit$data, control = dbartsControl(n.cuts = 20))` | the fit's grid, the 20 dropped | 20 |
| a probit fit's control into a gaussian `dbarts()` | gaussian | error |
| one forest with no basis declared by `forests =`: `getTrees()` | a `forest` column | - |
| `new("dbartsSampler", dbartsControl(n.trees = 20), fit75$model, fit75$data)` | 20 | 20 (the model has no count) |
| `new("dbartsSampler", dbartsControl(useQuantiles = TRUE), fit$model, fit$data)` | uniform, the fit's grid | quantiles |
| the same with `dbartsSpec()`'s model and data (no grid yet) | quantiles | - |
| `new("dbartsSampler", dbartsControl(n.cuts = 20), model, data)` | the data object's 100 | the data object's 100 |
| `dbartsSpec(..., control = dbartsControl(useQuantiles = TRUE))`, then `new()` with a fresh control | uniform | - |
| `bartBT()` with no `ntree` | 200 | 200 |

The `getTrees()` row is a defect against dec-B249 ("a single-forest fit returns 0.9-34's columns"; a single
forest with a basis keeps the column, its forest being combined through an amplitude), kept alive by
`bartcore.forestsDeclared`; it goes with the attribute. The direct-door rows are the critique's runs (p02, p04,
p05 in scratch/cm-crit/probes): 0.9-34's model carried no tree count and its data object no quantile rule, so a
control was the only place to state either at that door, while its cut count was always the data object's.
`rbart_vi()` reaches the sampler through `dbarts()`, so it follows `dbarts()`'s rule.

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
| tree count, per mean forest | the model's forest record: the count it runs and whether its caller stated it | `forest(n.trees = )`, the door's `n.trees` for the forests with no basis, the defaults by kind | yes (dec-B247) |
| the rest of a forest's declaration: vars, sd, base, power, label, held amplitude and its `amplitude.df`, interactions, blocks, basis term, move settings | the same record | `forest()`, the door's level | structure: count, vars, interactions, blocks, basis, label; parameters: sd, base, power, the amplitude's hold and prior (dec-B254); move settings are run settings |
| the variance forest | the model's variance record (count, base, power, columns, moves) | `varianceForest()`, `variance =` | yes |
| ordinal category count, negative-binomial shape setting, hazard periods | the model | the family | yes; the held shape a parameter (setmodel-parameters) |
| survival status | the data object, one value per row | the response, `setResponse(status = )` | no |
| cut count | the data object (`n.cuts`, exists; dec-B545) | `dbarts(n.cuts = )`, `dbartsData(n.cuts = )` | the grid is data (dec-B231) |
| quantile rule | the data object (`useQuantiles`, new, dec-B545; plain `TRUE` or `FALSE` on anything a sampler holds, `NA` only on a fresh object stating none, dec-B548) | `dbarts(useQuantiles = )`, `dbartsData(useQuantiles = )` | as the cut count |
| binary | none: the model's family | the family | yes |

New homes are typed S4 slots with validity, set after `new()` where `initialize` is closed. Every model built by a
door carries one forest record per mean forest, a single-forest model one and a multinomial model one that every
category forest takes. The record class has an unstated count (`NULL`) as its prototype, so a model built by
`new("dbartsModel")` by hand states none ([Every door](#every-door), the direct door). No reader guards the new
slots against an object saved by a development build: only shipped versions are owed compatibility, and a fit
saved by 0.9-34 fails as today (TODO legacy-object-message, dec-B461).

The record keeps whether a count was stated beside the count it runs, as setmodel-parameters' records keep a
stated base and power: a door-built model states a count only where its caller wrote one (the door's `n.trees`,
`forest(n.trees = )`, or a written control count); the default, by kind, is run and not stated. dec-B526 reads it;
the count is in neither `$getForest`'s list nor its statement.

The control keeps `binary`, `n.trees` and `useQuantiles` as mirrored slots and loses the `n.cuts` slot. A fresh
`dbartsControl()` holds the values a default fit mirrors (`n.trees` 75, `useQuantiles` `FALSE`, `binary`
`FALSE`), so `dbartsControl()@n.trees` reads as on main. A `dbartsControl(n.cuts = )` value, the tombstoned
spelling, rides an attribute named `n.cuts`, where 0.9-34 kept it: `ctl@n.cuts` reads back what the caller
wrote, as in 0.9-34, and a sampler's stored control carries none, as in 0.9-34. `dbarts.stated`, `bart()`'s
record of what its own caller named, stays as it is: written by the wrappers, read and removed by
[`resolveSamplerSpec`](../../R/spec.R); it now also names the count and the cut settings, which the wrappers
pass flat on every call (rule 6).

### Stated, carried and silent

A control's value for `n.trees`, `useQuantiles` or `n.cuts` is one of three:

- written: the caller's `dbartsControl()` call named it (`dbarts.supplied`), or its slot differs from the value
  the sampler last mirrored into it, or, on a control with no mirror record, from a fresh control's;
- carried: the control holds a mirror record (`dbarts.mirror`, written with the slots) and the slot equals it;
- silent: anything else, a fresh control's default.

From slice 2 a cut count has no slot: written, it rides the `n.cuts` attribute; carried, it is a field of the mirror
record.

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
door. That is dec-B334's warning; dec-B335 and dec-B385 name the same one. The remedy names what the call can do:
a wrapper given a data object has no forest term to write, so `bart(dd, control = )` says "state the count as
bart(n.trees = )", and "fit with dbarts() to state the count" where every forest of that data object has a basis.

### The rule at each door

At a door with flat arguments the levels are, for the forests with no basis: a forest's own `n.trees`; the door's
`n.trees`, or a written control value standing for it; the defaults by kind (dec-B274). For the cut settings: the
door's `n.cuts` and `useQuantiles`, or a written control value standing for them; what the data object carries
(a grid, else its own count and rule); the defaults (100, uniform).

1. A flat argument and a written control value in one call: the flat one is used (dec-B116), the control's draws
   the not-used warning where it differs.
2. The door's count and a forest's own: refused as today where every forest with no basis states one (dec-B241,
   dec-B503). A written control count there is refused the same way (dec-B527); today the forest's governs the
   control's in silence.
3. A count, flat or written, where every forest has a basis: refused by name (dec-B274). A carried one is not a
   statement and is never refused, which retires control-carried-count-edges' first case.
4. An explicit cut setting, flat or written, beside a data object that carries a grid: the grid is re-derived at
   the stated setting (dec-B231: "If explicit arguments are given, we could then follow them directly"), as
   0.9-34 re-cut there. A data object stating its own count or rule with no grid: used, at `dbarts()` as at
   `dbartsSpec()` (dec-B241: the cut settings describe the predictors). A data object stating neither, with no
   grid: the defaults, 100 and uniform.
5. The mirrors are written last, from the model and data the sampler takes, and whenever a sampler stores a
   control.
6. A wrapper that fits through `dbarts()` (`bart()`, `bartBT()`, `pdbart()`, `rbart_vi()`) passes its resolved
   count and cut settings flat on every call, `bartBT()`'s 200 trees and `bart()`'s 75 included, and names in
   `dbarts.stated` which of them its own caller wrote, as `bart()` does today for the tree prior. `dbarts()`
   reads statedness from that record, so a wrapper's unstated default is never refused beside a forest's count or
   beside forests that all have a basis, and is not stated in the record.

#### Every door

| door | tree count | cut count, quantile rule |
|---|---|---|
| `dbarts()` | `n.trees = NULL` after `forests`; the rule above | `n.cuts = NULL`, `useQuantiles = NULL` after `factors`; the rule above |
| `dbartsSpec()` | `n.trees = NULL` after `forests` (dec-B241's call); the rule | no flat arguments: a written control value is read as `dbarts()` would read it beside the same data object, which `?dbartsSpec` promises; today the control's is dropped |
| `dbartsData()` | - | `n.cuts = NULL`, `useQuantiles = NULL`, last; NA refused, `NULL` states none |
| `bart()`, `pdbart()`, `rbart_vi()` | their `n.trees`, passed flat to `dbarts()` by rule 6; the control they build names none of the three | their `n.cuts`, `useQuantiles`, passed flat by rule 6 |
| `bartBT()` | `ntree`, 200 unless its caller names another, passed flat by rule 6 | `numcut`, `usequants`, passed flat by rule 6 |
| `xbart()` | its grid; a written control count is a grid of one; each cell's count goes on the cell's model record | flat; the shared data handle takes the data object's rule and a grid it carries (dec-B549); a data object's own count and rule stand as at `dbarts()`, row 4 |
| a formula's `forest()` terms, `forests = ` | per forest; a bare `forest()` is one forest (dec-B241's call) | - |
| a data object as `data` | the rule; a data object carries no count | the rule, row 4 |
| `new("dbartsSampler", control, model, data)` | a copied control: the model's count, with the not-used warning where it differs (dec-B335); a fresh control: the model's, a silent value never being read; a written control count: used where the model states none, with the one-release warning, refused by name where it states another (dec-B526) or where every forest has a basis (rule 3) | a copied control: the data object's rule and count, with the not-used warning where they differ; a fresh control: the data object's; a written rule: dec-B526's reading, the grid re-derived; a written count: the data object's, with the not-used warning, as 0.9-34 and main use the data object's |
| re-creation (`getPointer`, `setState`, `copy`, a reload, stan4bart's restore) | the stored model's; the mirror agrees, so no warning | the stored data's |
| `$setControl(new)` | a carried count that differs: taken, the sampler's kept, the not-used warning (dec-B385); a written one that differs: refused, "changing 'n.trees' is not available on an existing sampler" (dec-B247) | the same two cases for each, the sampler's data holding them; the `n.cuts` attribute is stripped after the check |
| `$setModel(model)` | refused by name before anything is stored where a structure field differs (dec-B254): a forest's count, vars, interactions, blocks, basis or label, the variance forest's count and columns, the category count, the hazard periods, a count compared only where the model's record states it, one it does not state leaving the sampler's (dec-B551); the moved facts join setmodel-parameters' structure check. A parameter field (a forest's sd, base, power, the amplitude's hold and prior) goes where setmodel-parameters sends it (main refuses a sampler of several forests wholesale until its slice 4). Per-forest move settings are installed, as `$setControl` installs the control's (dec-B254: the moves "are not model"; dec-B518's `$setForest` changes them) | - |
| `$setData(data)` | - | a count or rule the data object states wins, one carrying a grid brings what made it, one stating neither keeps the sampler's, and the sampler's record follows; `updateCutPoints = "none"` keeps the sampler's (dec-B546) |
| `$setForest()`, `$getForest()` (dec-B518, unbuilt) | `n.trees` is structure, refused by name at `$setForest`; it is not in `$getForest`'s list | neither: the cut settings are data, changed through `setCutPoints` and `setPredictor(updateCutPoints = )` |

`binary`: the model's family decides at every door. At the direct door a control whose `binary` differs from the
model's family draws the not-used warning, where 0.9-34 took the control's (dec-B335's rule, applied to the flag).
At the front doors a control's `binary` is never read and never warned about.

`bart()` and `xbart()` lose the attribute refusal, keeping [`requireDbartsControl`](../../R/dbarts.R)'s class check: no control carries a
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
| `y ~ forest(x, n.trees = 30)` with `dbartsControl(n.trees = 50)` | - | 30 | refused (dec-B527) |
| `bartBT(x, y)`, no `ntree`; `pdbart()` on a BayesTree call | 200 | 200 | 200 |
| a `bart()` fit's control into a `dbarts()` whose every forest has a basis | - | refused | taken; 50 and 50 by kind |
| `new("dbartsSampler", fit200$control, fit75$model, fit75$data)` | 200 | 200 | 75; the not-used warning |
| `new("dbartsSampler", dbartsControl(n.trees = 20), fit75$model, fit75$data)` | 20 | 20 | 20 with the deprecation warning, the model stating no count (dec-B526) |
| `new("dbartsSampler", dbartsControl(useQuantiles = TRUE), fit$model, fit$data)` | quantiles | the fit's grid | quantiles, the grid re-derived, the data object stating no rule (dec-B526) |
| `sp <- dbartsSpec(..., control = dbartsControl(useQuantiles = TRUE))`, then `new("dbartsSampler", dbartsControl(), sp$model, sp$data)` | - | uniform | quantiles: the data object states the rule |
| `fit75$setControl(fit200$control)` | resized to 200 | refused | taken, 75 kept; the not-used warning |
| `fit75$setControl(dbartsControl(n.trees = 20))` | resized to 20 | refused | refused |
| `fit200$setControl(dbartsControl(n.samples = 100))` | resized to 75 | refused | taken, 200 kept, silent |
| `dbarts(dbartsData(y ~ ., d, n.cuts = 10))` | - | - | 10 |
| `dbarts(fit$data, control = dbartsControl(n.cuts = 20))` | 20 | the fit's grid | 20; the deprecation warning |
| `dbarts(fit$data, n.cuts = 20)` | - | - | 20 |
| `dbarts(fit$data, control = dbartsControl(useQuantiles = TRUE))` | quantiles | the fit's grid | quantiles; the deprecation warning |
| `dbartsSpec(dd, dbartsControl(n.cuts = 20))`, `dd` stating 10 | - | 10 | 20; the deprecation warning |
| a variance fit's control into a plain `dbarts()` or `bart()` | - | plain; `bart()` refuses | plain at both |
| `fit$control@n.trees`, `@binary`, `@useQuantiles` after any fit | the fit's | the fit's (first forest's count) | the same values on every model of one forest with no basis; on others the fitting function's level ([Mirror mechanics](#mirror-mechanics)) |

### Mirror mechanics

One way, from the model and data to the control, in one helper called at the end of
[`resolveSamplerSpec`](../../R/spec.R) (so `dbartsSpec()`'s returned control reads right), when a sampler stores
its control (every copy, reload and stan4bart restore), and in `$setControl`: `n.trees` the count at the fitting
function's level, the summed counts of the forests with no basis (the one forest's on a model of one; the shared
total under dec-B503, where three sharers of 25 read 75), or the first forest's where every forest has a basis;
`useQuantiles` the data object's rule, `binary` whether the family is binary (probit and logistic, a hazard fit
being one of them); the mirror record written beside them, the three names taken out of `dbarts.supplied`, the
`n.cuts` attribute stripped. The not-used warning compares at the same level, so a sharers' control
given back to the same formula draws none. man/forest.Rd's "A sampler's control@n.trees ... hold its first
forest's" changes with it. The bridge reads none of the three: after slice 1a it reads the model and the data
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
from the parsed family. `applyVarianceAttributes`,
`applySurvivalAttribute` and `applyForestAttributes`
read the model and data and are renamed for what they read: [`applyVarianceRecord`](../../src/R_interface_bartcore.cpp), [`applySurvivalStatus`](../../src/R_interface_bartcore.cpp), [`applyForestRecords`](../../src/R_interface_bartcore.cpp). [`bartcore_setControl`](../../src/R_interface_bartcore.cpp)
loses its tree-count refusal (the control has none to compare); `bartcore_setModel` gains it as a backstop for
callers that skip R's check (xbart's cell install). The data handle takes `useQuantiles` and a grid it records
(dec-B549) from the data object and
`sparseDensityThreshold` from the control as now. A slot holding NULL reaches the bridge as R's S4 NULL marker;
every new NULL-able slot is read through the normalizing helper `parseData` uses for `bases`. `bartcore_setData` hands
the engine the count and rule R chose, on `DataReplacement` (dec-B546), the one change under src/bartcore. A survival fit's data
object under another family runs: [`takeControlReadings`](../../R/records.R), which every door passes through, leaves
the status behind unless the model's family is `aft`, so no sampler stores it. The bridge's refusal stays, as the
backstop `$setResponse(status = )` reaches on a family that is not `aft`. Otherwise unchanged:
src/bartcore, the facade, `SamplerOptions` and what goes into it, the order of construction and generator
creation, inst/include/dbarts/dbarts.h.

### What changes for a user

Refused where main runs: a written control count beside a forest that states its own (dec-B527); a written
control count or rule at the direct door beside a model or data object that states another (dec-B526), and a
written count there where every forest has a basis; a written
`useQuantiles` or `n.cuts` on `$setControl` that differs from the sampler's (main refuses `useQuantiles` already;
`n.cuts` is stored there and never reaches the engine); `$setModel` with a model of another structure in a fact
that today rode the control, a tree count among them only where the model's record states it (dec-B551; a spec
stating none is taken and keeps the sampler's, as on main; A/B rows `b551.*`); at the direct door, a basis forest's model beside a data object with no bases,
`new("dbartsSampler", run(), bcfFit$model, f75$data)`, which the bridge refuses ("a basis forest was configured but
the data carry no forest bases") where main followed the control and fit gaussian with 75 trees; an `n.cuts` naming
no column or term, leaving a column with no value, or unnamed and neither one value nor one per column (dec-B547), so
`c(x2 = 10)` and `c(25, 10)` on three columns, which main and 0.9-34 take (a hazard fit's period column is read by its
own name only, so that its refusals count the caller's predictors, dec-B553); `dbartsControl()@n.cuts` and a sampler's `control@n.cuts`, which read
100 on main and error here, "no slot of name" (dec-B552; a written `dbartsControl(n.cuts = 20L)@n.cuts` reads 20 on
both; A/B rows `b552.*`).

Taken where main refuses: `forests = forest()`; `dbarts(n.trees = )`, `dbartsSpec(n.trees = )`,
`dbarts(n.cuts = , useQuantiles = )`, `dbartsData(n.cuts = , useQuantiles = )`; any control at `bart()` and
`xbart()`; a carried or fresh control's differing count or rule at `$setControl`; a carried control count beside
forests that all have a basis; at the direct door, a control from one fit beside another fit's model and data where
the control's hidden facts refused the pair (100 of the 156 crossed triples main was run on, slice 1a).

Different with no error: a copied control's count, rule and cut count at `dbarts()`, `dbartsSpec()`, `bart()` and
`xbart()` (dec-B334, call 5) and its count and rule at the direct door (dec-B335), each with the not-used warning;
at the direct door, from slice 1a, a control from one fit beside another fit's model and data follows the model and
not the control's hidden facts (21 of the 156 crossed triples differ among those main ran, 35 identical, 12 of
those an unused-argument error on both builds, `binary` not being a `dbartsControl()` argument): a plain
fit's control beside a variance fit's model runs the variance forest, and a binary control beside a multinomial
model runs the multinomial model, with its draws, and reads `binary` FALSE; a fresh control beside a fit's model or
a `dbartsSpec()` data object takes their count and rule (main: the control's 75 and uniform); `dbarts()` and `xbart()` given a
data object stating its own cut count or rule; a written cut setting beside a data object carrying a grid (0.9-34's result restored); at the direct door, a written
rule beside a fit's data object re-deriving its grid (dec-B526; 0.9-34's result restored); a
fit of one forest with no basis loses `getTrees()`' `forest` column however the forest was written (dec-B249);
on a model whose forest with no basis is not first or which has several, `control@n.trees` reads the fitting
function's level; `$setModel` installs the move settings of a model's records, where main left the sampler's;
`dbartsControl()`'s three arguments warn. `$setData` with a data object stating a cut count or rule, or carrying a
grid (dec-B546): what it states wins and the sampler's record follows, one carrying a grid bringing the count and rule
that made it, an unstated rule reading as even, where main keeps the sampler's 25 even and after another fit's object
gives 25 even points beside 10 quantile on the first refresh; one stating neither keeps the sampler's, and
`updateCutPoints = "none"` keeps the sampler's whatever is stated (A/B rows `b546.*`). `fit$data@useQuantiles` reads
plain `TRUE` or `FALSE` on anything a sampler holds (dec-B548), where main has no slot to read (A/B row `b548.read`).
A per-column `n.cuts` at every door that takes it reads as `split.probs` reads (dec-B547): `c(x3 = 5, x1 = 25, x2 = 10)`
gives each column its own, where main and 0.9-34 read it by position (A/B rows `b547.*`). On a hazard fit the appended `period` column takes `n.cuts` only by its own name and is otherwise
at the default 100 (dec-B553, revising dec-B547 in part): `c(5, 7, 9)` on three predictors gives 5, 7, 9, 100, `5`
gives 5, 5, 5, 100, `c(x1 = 5, .default = 9)` gives 5, 9, 9, 100 and `c(period = 5, .default = 5)` gives 5, 5, 5, 5; four
unnamed values stay by position, `c(x2 = 7)` and `c(5, 7)` are refused counting the caller's predictors, and a count a
data object carries for its predictors leaves the period at 100. `dbartsData()` takes a count a hazard sampler could read, one value short or by name without `period`, as written, and a sampler of any other kind refuses it where it reads it; no column's name switches a check, and it refuses at once only a count no sampler could read. Main recycled `5` over the period and took
`c(5, 7, 9)`, reading the period at 5 and 100; a hazard fit at `n.cuts = 5` no longer matches the binary fit on the
expanded design at `5` unless it names `period`, and draws move only on a hazard fit stating a count (A/B rows
`b553.*`). `split.probs` keeps its refusal naming `period`. A data object's slots are read-only copies (dec-B554): a
count is stated through `dbartsData(n.cuts = )` or `dbarts(n.cuts = )`, and a count written on the slot by hand is no
way in, so this plan has no row for it. `xbart()` given a data object carrying a grid
takes it for every fold (dec-B549), a stated `n.cuts` or `useQuantiles` deriving it again: a fit's own object with its
predictors as they were gives the number the formula gives, to the bit, and a hand grid or edited predictors stay on
the grid carried, where main and 0.9-34 derive (A/B rows `b549.*`).

## Slices

### Slice 1a: the records (R and bridge), bitwise neutral

The model's forest records (count and statedness among them) and variance record, its category count, shape and
hazard periods; the data's status and `useQuantiles`; coefficient-prior-opt-in's `amplitude.df` with the rest of
the forest record; the mirror helper and record; the bridge reads of [The bridge](#the-bridge);
[`resolveSamplerSpec`](../../R/spec.R) writes the records in place of the attributes, from what main resolves;
the readers move: [`fitTreeCounts`](../../R/generics.R), [`forestPriorLines`](../../R/forestPrint.R), the fit
packagers in R/bart.R, the hazard periods of partial dependence, `getSigmas`, `samplePriorPredictive`,
`getTrees` (its column still as main), [`resolveForestSpreads`](../../R/dbarts.R) and
[`writeLeafPrior`](../../R/dbarts.R), `setResponse(status = )`, xbart's cell models and handle, forest-after-merge's
label read if it has landed; `$setControl`'s attribute loops go. No reading changes: the doors read the count and
the cut settings from the control exactly where main does, `control.n.trees`, the carried count of a model whose
first forest has a basis, riding the mirror record until slice 1b, and `$setControl` keeps main's results. The
A/B is identical on every call but a control from one fit beside another fit's model and data at the direct door,
which follows the model (accepted by the orchestrator, the maintainer told). ~1450 lines: R ~500, bridge ~250,
tinytest ~620 (the 80 attribute lines follow the objects; the read-backs of [Tests and gates](#tests-and-gates);
a new test-control-migration.R), docs ~80.

### Slice 1b: the copied control and the structure check (R)

A copied control stops speaking at every door: dec-B334's not-used warning; the direct door's copied count and
rule (dec-B335); `$setControl` taking a carried differing value and refusing a written one (dec-B385, dec-B247);
the mirror at the fitting function's level and the comparison there; `control.n.trees` and the doors' carried
branches go; `$setModel`'s structure check gains the moved structure fields, and per-forest move settings are
installed there; `getTrees`' column by dec-B249; the attribute refusal of a control at `bart()` and `xbart()` goes
([`requireDbartsControl`](../../R/dbarts.R) is the class check left), so both take any control; a fresh control at
the direct door reads as silent. A written value keeps main's reading at every door but `$setControl` until slice
2. ~600 lines: R ~200, tinytest ~380 (the dec-B334, dec-B335 and dec-B385 rows at every door; `$setModel` per structure field, sampler and stored model unchanged after a
refusal; the column per model kind), docs ~20. The A/B's differing rows are exactly these.

### Slice 2: the doors (R)

`n.trees` on `dbarts()` and `dbartsSpec()`, `n.cuts` and `useQuantiles` on `dbarts()` and `dbartsData()`, the
bare `forest()`; written values with both warnings at every door of [Every door](#every-door), the direct door by
dec-B526 and a written count beside a forest's refused by dec-B527; rule 6's wrappers (`bart()`,
`bartBT()`, `rbart_vi()`, `pdbart()`) and `xbart()` pass the three flat and build a control naming none, so
control-carried-count-edges' first case cannot arise; `dbartsControl()`'s three arguments default to `NULL` and
warn; the `n.cuts` attribute; the deprecation key in [`warnOnce`](../../R/utility.R)'s scheme; the registry rows.
~1350 lines: R ~560, tinytest ~790 (every row of [Checked against scripts](#checked-against-scripts) at every
door that has it; each way a value is written, carried or silent; a refusal leaving the sampler untouched; the
mirrors after every family's fit; the warnings counted with `withCallingHandlers`).

### Slice 3: respell (tests, benchmarks, examples)

The ~520 constructor calls naming the three, their later uses and the wrapper helpers that forward to
`dbartsControl()`, moved to the flat arguments or `forest(n.trees = )`; man examples (8 pages) and vignettes (3).
No R/ or src/ change. Tests that pin a 0.9-x spelling on purpose keep it and reset the warning key. The tinytest
calls move now and not at 1.1-0: the suite runs in one process, so the first stray deprecation warning spends the
session's key and a later file that asserts it, or counts warnings, then depends on file order. ~2000 changed
lines. Its gate is mechanical ([Tests and gates](#tests-and-gates)).

### Slice 4: help, NEWS, docs; consumers

Help: dbarts.Rd, dbartsSpec.Rd, dbartsData.Rd, dbartsControl.Rd (the three arguments as tombstones, the slots as
copies), forest.Rd (the first-forest paragraph: the control's count is the fitting function's level and is
refused beside a forest's), dbartsForests.Rd, dbartsSampler-class.Rd (`setControl`, `setModel`, the direct door),
bart.Rd, xbart.Rd. NEWS against 0.9-34: the new arguments; the three `dbartsControl()` arguments warn and expire at
1.1-0; the direct door takes the model's count and the data object's rule where they state them (dec-B526);
`$setControl` keeps the sampler's count (dec-B247's item reworded:
a control naming another count is refused, one carrying another is taken with a warning); a control taken from a
fit can be passed to any door. Design: docs/design/feature-matrix.md, the forest-defaults-by-kind design's "Where
the numbers live", docs/architecture.md's control paragraph; records. ~350 lines.

Consumers, each on its own branch, green against this branch's tip before the merge:

- stan4bart: `bart_args$n.trees` to `dbartsSpec(n.trees = )` and `useQuantiles` to the data object, both kept out
  of the control loop and named in its known arguments; R/mvbart.R sends `n.cuts` and `useQuantiles` to
  `dbarts()`. Its restore through `new()` reads agreeing mirrors and stays silent. ~60 lines with tests.
- bartCause: `bcf()` keeps `n.cuts` and `useQuantiles` out of the control and lets its `dbarts()` loop pass them;
  test-21's `@n.trees` read stays right. ~20 lines.
- treatSens: its `data.bart@n.cuts` write at `makeBartSpecs` moves to `dbartsData(n.cuts = )` (a data object's
  slot is no way in, dec-B554); the rest of the slot writes in the consumers wait on their own entry in the TODO
  (consumer-slot-writes). bairrtt: nothing; the suite runs against the tip.

## Tests and gates

Slice 1a, per moved quantity at a non-default value, read back from the engine on the live sampler, after `copy()`
and after a save and reload (each record field, not only the variance forest and the status), each with a mutant
that replaces the read with its default:

| quantity | read back by |
|---|---|
| the variance forest's count, columns, base, power | the stored state's variance trees; seeded draws against the default's; equal after `copy()` and a reload |
| a fixed negative-binomial shape; the ordinal category count | `getShape()`; the thresholds' length |
| per-forest count, vars, interactions, blocks, moves, three forests | `bartcoreForestTreeCount`, `getTrees(forest = )` split variables, the move install |
| per-forest sd, held amplitude, `amplitude.df`, labels, basis terms | `getLeafPrior()`, `getForestAmplitudes()`, a label door, `predict()` on new rows |
| a multinomial model's one record | each category forest's tree count |
| the mirror after creation, `copy()`, a reload, stan4bart's restore | `control@n.trees`, `@binary`, `@useQuantiles` against the model and data, and no not-used warning on any of them |
| the status after `setResponse(status = )`; the hazard periods | a reloaded sampler's seeded draws against its unsaved twin; the fit's `$periods` and `pdbart` after a reload |
| the family alone deciding binary | each family's run channels, a hand-built model included |
| `useQuantiles` at each door and under `xbart()` | the state's cut points, quantile against uniform |

Slice 1b: the dec-B334, dec-B335 and dec-B385 rows and the `$setModel` fields as tests. Slice 2: the script table
and the door table as tests; the deprecation warning once per key per session and the
not-used warning on every call, counted; none for a silent or agreeing value. Slice 3: a creation-trace driver logs
each sampler creation of the full suite and the three equivalence harnesses (per file, in order: each forest's
engine tree count, `data@n.cuts`, the rule), before and after; the tables must be equal, the expectation count
unchanged.

RNG gates, on every slice's tip and on the merged tree: tests/cpp; the full tinytest suite in one process at the
default ncpu; the equivalence trio bit for bit on a reference build (55, 19 and 11 "identical draws" lines, no
`max |z|`); the four snapshot files unchanged; a scratch A/B digesting seeded draws for a matrix of calls (every
door, family, forest layout, cut setting, a copy, a reload, `xbart`) on a c55a1477-plus-tree-prior-flat build and
on the slice, identical except the listed rows, which must match their expected new values. The matrix and the door
tests hold a survival fit's data object under another family at each door, and a hand model with an empty forest
list at the direct door, `$setControl`, `$copy` and a reload. The fit object
carries new slots, so `.github/workflows/exact-gates.yaml`'s gates run with `quick`. On slice 1a's tip the A/B
is identical but for a control crossed with another fit's model and data at the direct door; on 1b's and 2's it differs exactly on
[What changes for a user](#what-changes-for-a-user)'s rows. Slice 1a: sanitizers on
tests/cpp and on the R-loaded path. Every slice: `lintr::lint_package()` on the slice's library, `air format
--check .`, tools/check-rc-codoc.R, tools/check-win-drift.R, tools/check-doc-freshness.R, `R CMD check --as-cran`
from a clean tarball, and the mutation battery's `verify-anchors` (slice 1a touches src/).

## Mutants

The review runs each and requires a failing test:

1. The bridge reads a fixed 75 for the first forest's count.
2. The bridge reads the variance forest's count as 40 whatever the record says.
3. `binary` taken as FALSE for logistic (the parsed family's flag dropped).
4. `useQuantiles` read as FALSE from the data object.
5. The status not read on re-creation (the reload twin differs).
6. A forest's moves dropped from the record (the move install reads the control's).
7. `getTrees()`' column shown for a fit of one forest with no basis, or dropped for one forest with a basis.
8. A carried count treated as written (the copied control's 50 used).
9. A written count treated as carried (the edited 30 dropped).
10. The not-used warning skipped at the direct door.
11. `$setControl` refusing a carried differing count.
12. `$setControl` taking a written differing count.
13. The flat `n.trees` losing to a written control count.
14. A written control count beside a forest's taken in silence (dec-B527).
15. `dbarts()` overwriting a data object's own count with the default.
16. A written cut count beside a carried grid ignored.
17. The mirror written before the forest's statement resolves (the control reads the door's count).
18. `$setModel` accepting a model with another variance forest.
19. The deprecation key not shared (two warnings for one name).
20. xbart's cell count left on the control.
21. `bartBT()` passing no count where its caller named none (75 trees in place of 200).
22. A wrapper's unstated default marked stated (an all-basis model refused, a forest's count refused beside it).
23. The mirrored `binary` taken from the control's own slot rather than the family (bartMan and tidypredict read
    it); the mirrored `useQuantiles` taken from the control rather than the data object.
24. The mirror not written when a sampler stores its control (every copy, reload and restore then warns).
25. A multinomial model's category forests taking the default count in place of the record's.
26. `amplitude.df` dropped from the record (a student() coefficient read as the normal's).
27. The mirrored `n.trees` taken from the first forest on a model whose forest with no basis is second.
28. `$setModel` refusing a model that differs only in a forest's move settings, or in a parameter field.

## Calls made

1. A per-forest move setting lives in the forest record: a forest is declared once, and a control passed on must
   need no handling (dec-B239). `$setModel` installs a record's move settings, as `$setControl` installs the
   control's (dec-B254: the moves "are not model"; dec-B518 lets `$setForest` change them), and refuses only the
   structure fields. Alternative: keep `bartcore.moves` on the control.
2. `dbartsControl()` keeps the three formals, defaulting to `NULL`, until 1.1-0, so callers that route by its
   formals (bartCause, stan4bart) keep reaching them and the name resolves (dec-A02). Alternative: move them to
   `...`, which makes such callers drop them.
3. `dbartsControl(n.cuts = )` rides an attribute named `n.cuts`, as in 0.9-34, so `ctl@n.cuts` reads what the
   caller wrote. Alternative: a private name, under which a 0.9-x read errors.
4. A written value beaten by a flat argument draws the not-used warning where it differs, beside the deprecation
   warning. Alternative: the deprecation warning alone.
5. `useQuantiles` and `n.cuts` follow dec-B334 and dec-B385's rule for the count (a copy never speaks; it warns
   where it differs), as dec-A158 treated the three alike.
6. `binary` is never a statement at a front door; at the direct door the model's family wins with the not-used
   warning (dec-B335's rule).
7. The mirrored `n.trees` is the count at the fitting function's level, the summed counts of the forests with no
   basis (the shared total under dec-B503), or the first forest's where every forest has a basis, and the
   not-used warning compares there, so a sharers' control given back to its own formula draws no warning. On a
   model of one forest with no basis, every 0.9-34 shape, it is that forest's count, as on main. Alternative: the
   first forest's, as on main, which reads 25 on three sharers of 75 and warns on its own reuse.
8. At `dbartsSpec()` a written control cut setting is read as at `dbarts()` beside the same data object, which
   `?dbartsSpec` promises; today it is dropped. Alternative: refuse it.
9. A wrapper passes its resolved count and cut settings flat on every call and records in `dbarts.stated` which
   its own caller wrote (rule 6), so `bartBT()` keeps BayesTree's 200 trees and an unstated default is never
   refused. Alternative: pass only what the caller named, which gives `bartBT()` 75.
10. At the direct door a written control cut count is not read: the data object's count is used, with the
    not-used warning where it differs, as 0.9-34 and main both read the data object's there. The written rule is
    dec-B526's.
11. A data object that states no rule and carries no grid takes the default, uniform; `NA` in a fresh one's slot means not
    stated.
12. The forest record keeps whether its count was stated, for dec-B526.
13. No guarded accessors for development-build objects (compatibility is owed to shipped versions only).
14. Slice 1 is split: 1a moves the records and changes no reading but the crossed direct door; 1b carries dec-B334, dec-B335, dec-B385, the mirror's level, the `$setModel` fields and dec-B249's column.
15. The flat arguments default to `NULL`, meaning not stated, the help giving the values that then apply (75 by
    kind, 100, `FALSE`).
16. Slice order: records before doors, so the doors are written once against the records.
17. A fresh control at the direct door is silent, so the model's count and the data object's rule stand, from slice
    1b: a silent value is never read. Alternative: keep main's 75 and uniform until slice 2, overruling a model that
    states a count with an object that states nothing. dec-B335; NEWS says it.
18. Slice 1b's control still has its `n.cuts` slot, so the mirror records it and writes it from the data object, and
    a carried one stops speaking with the count and the rule (call 5). Alternative: leave `n.cuts` speaking until
    slice 2 removes the slot, one more copied value carried a slice longer.
19. A carried value set aside at a door rides the attribute `dbarts.carried` to the end of
    [`resolveSamplerSpec`](../../R/spec.R), where the warning compares against the model and data the call took; a
    data object handed in keeps its own rule and count beside a carried one. Alternative: warn at the door, before
    the count the call takes is known.
20. The not-used warning names the flat arguments slice 2 adds (`dbarts(n.trees = )`), since one feature merges once
    (dec-B494). Alternative: name them from slice 2, rewording every message and pin twice.
21. `$setModel`'s check names a model with no forest record, one made by hand, as stating no change, and by dec-B551
    (ruled) a record not stating its count, so it keeps the
    sampler's structure as main does, and compares only the fields this plan names: `stated` and what forest-door
    adds to a record do not refuse. Alternative: compare whole records, which would refuse forest-door's additions.
22. A parameter of the moved records (a variance forest's base and power, the shape setting) given to `$setModel`
    stays the sampler's as on main until setmodel-parameters sends it elsewhere. Alternative: refuse it now, against
    dec-B254.
23. A written count at the direct door where every forest has a basis is refused by name, rule 3 at this door. This
    narrows dec-B526(d), whose "used where they state none" would otherwise take it; rule 3 and dec-B527 refuse it, and
    no release ran the shape. Alternative: the first forest takes it, as main, a count on a forest of another kind.
24. A carried cut count is a field of the mirror record, since the slot goes, so a copied control's differing count
    draws the not-used warning at every door. `@n.cuts` on a fit's control then errors, as a fresh control's does, on
    this branch only: both read 100 on main, and on 0.9-34 a fresh control's reads 100 and a fit's errors (dec-B552,
    which leaves the error).
    Alternative: carry nothing, a copied control's cut count changing the grid with no message.
25. An edit of a fit's value back to a fresh control's reads as written, so it warns and is read. Alternative: treat
    it as silent, dropping what the caller wrote.
26. RULED, dec-B544. `xbart()` given a data object keeps the data object's count and rule under rule 4, as `dbarts()`
    does, and dec-B549 adds the grid it carries. Alternative: main's overwrite, one fitting function reading a data
    object two ways. NEWS has it already: a data object saved by 0.9-34 is rebuilt with `dbartsData()`.
27. A written value agreeing with the model or data object at the direct door or `$setControl` is taken with the
    deprecation warning alone. Alternative: refuse any written value there, an agreeing script stopping for nothing.
28. `dbartsData(dd, n.cuts = )` on a data object joins that function's existing "are ignored" warning, a `NULL` not
    counting. Alternative: re-derive the grid, which `dbarts(dd, n.cuts = )` does already.
29. A written `NULL` at `bart()`'s, `xbart()`'s and `rbart_vi()`'s three and `bartBT()`'s `ntree`, `numcut` and
    `usequants` is refused (dec-B519): their defaults are not `NULL`.
30. A hurdle fit's positive part muffles the not-used warning, so the two parts give one between them. Alternative:
    the same message twice from one call.
31. The three slot rows of the registry carry dec-B250's `condition`, the form the registry's test requires of a row
    with no fixed expiry.
32. `resolveSamplerSpec` loses its `written` argument: the door stamps `dbarts.stated` on the control, so a wrapper
    and a door are read one way.
33. REPLACED by dec-B548. A fit whose call stated no rule reads `FALSE` in its data object, with a private marker
    that no caller stated it, so a written rule beside that data object at the direct door is still used and not
    refused. The marker goes with `dbartsControl(useQuantiles = )` at 1.1-0.
34. At the direct door a data object from `dbartsData()` stating no count still stops at the bridge ("maximum number
    of cuts cannot be NA"), as on main. Alternative: apply the default there, a change no row of this plan asks for.
35. RULED, dec-B553. A hazard fit's period column takes `n.cuts` only by its own name and is otherwise at the
    default 100, so one value or one per predictor leaves it there. Alternatives: refuse a count that leaves it out
    (as built before), or let a single value cap it.
36. `dbartsData()` tests no column's name: it refuses at once only a count no sampler could read, too long, naming a
    column that is not there, or missing a column other than the last, and keeps the rest for the sampler that
    reads it. A count one value short is then refused later, by any sampler but a hazard one, and on a plain matrix
    of three columns `c(7, 8)` could be two predictors and a period, so it waits too. Alternative: a test of the last
    column's name, under which a person-period matrix with that column named `t` was refused at once.
37. A hazard door given a count with names on a matrix whose last column is not named `period` refuses with
    "unrecognized variable name 'period'", a name the caller did not write. Left as is: it takes that matrix and
    that count, and the sampler is unchanged. Alternative: a message of its own.
38. Beside several forests with no basis that all state a count, a count on the control is refused by name, as dec-B527
    has it beside one and dec-B503 for the fitting function's own; one forest leaving it out reads it. Main ran the
    call and ignored the control's count. No release had `forests`. Alternative: main's silence.

## Open calls

Both are ruled; the text is what the maintainer was shown.

1. RULED, dec-B526, option (d). A tree count or quantile rule written on a control at `new("dbartsSampler", control, model, data)`.
   Background: a sampler can be assembled by hand from a control, a model and a data object. In 0.9-34 the model
   held no tree count and the data object no quantile rule, so the control was the only place to state either
   at that door: `new("dbartsSampler", dbartsControl(n.trees = 20), fit$model, fit$data)` built 20 trees, and
   `dbartsControl(useQuantiles = TRUE)` there gave quantile cuts; the cut count was always the data object's.
   Main still takes the written count (20) and takes the written rule where the data object carries no grid
   yet. Now every model a door builds holds its count and every data object can hold its rule. Your ruling
   dec-B335 settles a control copied from a fit: the model's count wins, with a warning ("a sampler assembled by
   hand takes the model's tree count"); a value the caller wrote is not named there. The agents' earlier call
   dec-A158, in the part dec-B335 and dec-B385 left standing, says a written value is "used, with the
   once-per-session warning, or refused by name where the model, the data or the sampler states another".
   Options:
   (a) The model's count and the data object's rule always win, with the not-used warning. One rule at the
       door; every 0.9-34 script of that shape changes its count or its cuts and is told, and it reads dec-B335's
       rule text as covering written values; against dec-A158's words.
   (b) The written value always wins, with the one-release warning, as at the front doors. 0.9-34's result for
       one release; a model that states 75 trees is overruled by an object that is not the model, and the door
       then has two rules by where the value came from.
   (c) A written value that differs from what the model or data object holds is refused by name. Loud, never
       silent, dec-A158's words; but every model a door builds holds a count, so nearly every 0.9-34 script of
       that shape stops.
   (d) The model and data object record whether their caller stated the count and the rule (the record keeps
       this anyway, call 12): where they state one, a differing written value is refused by name, as (c); where
       they state none (a default fit's model, a fit's data object that asked for no rule), the written value is
       used with the one-release warning, the grid re-derived under a written rule. Keeps every 0.9-34 script of
       that shape working for one release, refuses only a real conflict, follows dec-A158's words and reads
       dec-B335 as covering the copied value alone. Costs about 40 lines beyond (a) and one more branch to test.
   I recommend (d).
2. RULED, dec-B527, option (a). A count written on a control beside a forest's own `n.trees`, as `dbarts(y ~ forest(x, n.trees = 30), d,
   control = dbartsControl(n.trees = 50))`. Background: today the forest's 30 is used and the control's 50
   dropped with no message, and man/forest.Rd says a forest's count "governs it over the control's n.trees",
   while `bart()`'s and now `dbarts()`'s own `n.trees` beside the forest's is refused as the same count given
   twice (dec-B241, extended by dec-B503 to the case where every forest with no basis states one). A
   surface-audit question held for you, whether a setting stated at two levels is refused everywhere or the
   forest's wins everywhere, was settled for the tree prior by dec-B512, which refuses a name stated at both
   levels; this is what is left of it. Under dec-A02 a relocated argument is the argument: the control's count is
   the fitting function's count spelled the 0.9-x way. No release had `forest()`, so no released script is of
   this shape. Options:
   (a) Refused by name, as the plain `n.trees` beside it is: one rule for one level, in either spelling; main's
       pin and the forest.Rd sentence change, and a script that mixes a 0.9-x control with a new forest term
       stops with a message naming the remedy.
   (b) The forest's count wins, with the not-used warning: no call stops, but the same count is refused in one
       spelling and taken in the other, for one release.
   (c) The forest's count wins in silence, as main: no change, and a stated value is dropped with no message.
   I recommend (a).

## Order of landing

1. coefficient-prior-opt-in landed (dec-A245): slice 1a moved its `amplitude.df` with the rest of the forest
   record.
2. tree-prior-flat landed (dec-A246). This plan is written against it: its two-level rule for `power` and `base` is
   the rule dec-B527 gives the count. This branch is rebased onto main after both.
3. forests-shared-total (slices B and 1 built, slice 2 unbuilt): its slice 2's [`forestParams`](../../R/model.R)
   split of the count among sharers writes into the forest records if slice 1a lands first, and into the params
   vector otherwise, one site either way. Its R item 2's last sentence, that the forests' record keeps the count
   the next fit inherits where forest 1 is not plain, is the carried count dec-B334 retires: that sentence goes
   from its plan whichever lands first. Its bridge edits (`createHolder`, `printInitialSummary`,
   `bartcore_getLeafPrior`) meet slice 1a's in src/R_interface_bartcore.cpp. Engine features merge serially;
   slice 1a edits no engine file.
   forest-amplitude-normal lands after forests-shared-total, and nothing forces it before slice 1a: only
   forest-door's slice 2 and setmodel-parameters' slice 4 need it first.
4. setmodel-parameters has not landed, its structure check (`modelStructureRefusal`) and its `$setModel` shape
   write are on its own branch. Slice 1b checks the moved fields in
   [`movedStructureRefusal`](../../R/records.R), which the second to land joins to `modelStructureRefusal`; whichever
   is second also sends the shape write to the model's shape slot and stops keeping a shape the first keeps.
5. forest-after-merge landed (dec-A244): slice 1a's records hold the labels `setForestBasis` names a forest from.
6. Every "whichever is second" above is expected to fall to this branch, which lands after them: the budget
   carries ~3 agent-hours of rebasing in slice 1a, each rebase re-running the slice's A/B and suite.
7. forest-door (dec-B518) is built in its own worktree, stacked on slice 1a's tip, so `$setForest` refuses `n.trees`
   whatever the record says; it rebases as this branch moves and lands after this
   branch merges (one merge per feature, dec-B494).
8. Consumers: stan4bart's and bartCause's branches green against this branch's tip, merged in the same sitting as
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
- Verified since the first draft: `rbart_vi()` fits through `dbarts()` (the critique read its two calls), so it
  follows `dbarts()`'s rule.
