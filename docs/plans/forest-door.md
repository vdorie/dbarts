# forest-door: a live sampler changes one forest through $setForest, by label, in forest()'s names

Status: PLANNED 2026-10-30 on c55a1477 (dec-B518); amended the same day after its blind critique (REVISE FIRST,
scratch/fd-crit/report.md, probes c1.R to c9.R there) and the coordinator's calls on it, and again after its
re-check (BUILD AFTER AMENDMENTS, scratch/fd-crit/recheck.md, probe c10.R). Six calls are held for
the maintainer ([Open calls](#open-calls)). Rulings: dec-B518 (the door), with dec-B254 (parameters, never
structure), dec-B276 (a forest by label or position), dec-B477 and dec-A239 (per-forest move settings), dec-B511
(the coefficient's prior live), dec-B512 (plain power and base), dec-B514 and dec-B515 (a setter means what
creation means; dec-B515 revised in part by dec-B518), dec-B516 (revised in part by dec-B518), dec-B129 (leaf
where it means a leaf), dec-B247 (a live count does not change), dec-B459 (copy of a live chain), dec-B500 (no
`treeShift` on the variance forest), dec-B505 and dec-B509 (forests that share a total), dec-B513.

agent: one opus implementer (effort high) for slices 1 and 2, one opus review over both that runs the mutants
below, the same reviewer checking the fix round unless the fixes change engine behavior; slice 3, if open call 6
takes it, one opus implementer and one opus review. Tier "Changes draws"
([Process by risk](README.md#process-by-risk)): slice 1 adds an engine install and slice 2 a door that changes
priors between sweeps; slice 3 edits the shared-total group in the engine.
rng: NEUTRAL, bit for bit, for every call main accepts and keeps. Three main calls, none shipped, stop being
taken: `$setLeafPrior(forests = )`, removed and respelled `$setForest(forest, sd = )`, which draws bitwise as the
old spelling did (an A/B against a c55a1477 build); and, under open call 5's (a), `$setLeafPrior(normal(k = 3))`
on a `forest(sd = 2)` sampler and `$setModel` with a changed power on a `forest(x, power = )` sampler, each
refused by name. Every new call is held to the twin oracle: before any run its next draws are
those of a sampler created with the same statement, one named handoff excepted ([RNG class](#rng-class)). No
call draws a random number.
window: for 1.0-0, after the merge to main (dec-B518 names none; nothing per-forest shipped in 0.9-34).
budget: ~2900 lines in the tree (~2500 if open call 6 refuses). Slice 1 ~250 (engine ~60, bridge ~60, tests/cpp
~130). Slice 2 ~2150 (R ~650 added and ~130 removed, the levels' records among them; help ~170, design notes ~25,
NEWS ~10, tinytest ~1300, test-forest-door.R ~1000 of it over the kinds-by-names grid). Slice 3 ~460 (engine
~100, bridge ~30, R ~40, help ~10, tests/cpp ~120, tinytest ~160), or ~40 for the refusal. Records ~50. Agent
time, estimated from comparable landings and not measured: slices 1 and 2 ~22, their review with mutants and the
fix round ~7, slice 3 ~9 (or ~1), gates ~4, rebases onto control-migration ~4: ~46 (~38). Plans here run 1.5 to 2
times low, so 55 to 90; the critique put it at 3200 to 4200 lines and 50 to 70 hours.

## Summary for the maintainer

You ruled that a live sampler changes one forest at a time, naming it by its label, in the names `forest()`
takes: `s$setForest("tau", sd = 0.5, power = 3)`. Its reader, `s$getForest("tau")`, gives back what the forest
states, so that handing it back changes nothing. The variance forest gets its own pair. The sampler-wide doors
keep only the fitting function's level.

What a call does. It changes only what it names, from the next sweep, and leaves the trees, leaf values and
coefficients already drawn as they are. Each name means what it means when the sampler is created: `sd` on a
held forest is its size, on a drawn one its prior median; `power` and `base` are the forest's tree depth; the
move settings stop following the control once a forest states its own. A forest's columns, tree count, basis,
interaction limits and blocks are refused by name, as you ruled: they change what the model is, and the route is
a new sampler. Anything creation would refuse, the door refuses too, with the same reason. Before the sampler has
run, an accepted call gives exactly the draws a sampler created with the same statement gives. A call that states
only what is already in force touches nothing at all.

One note you may veto: the door also takes `treeShift`, the third per-forest move setting `forest()` takes,
though your ruling's list leaves it out.

What the reader gives. The forest's own statement, one entry per name, empty where the forest states nothing of
its own and takes the level above (the fitting function's settings for a forest with no basis, its kind's
defaults for one with a basis, the control's for the move settings). Beside it, for reading, every value in
force, including the tree count, columns and depth. Today no reader shows a forest's depth.

What this needs that main lacks. Today the sampler forgets, once created, which level stated a setting:
`dbarts(y ~ x, power = 1.5)` and `dbarts(y ~ forest(x, power = 1.5))` leave identical records. The door cannot
return a forest to "the fitting function's setting" if that setting was not kept. This item keeps it, and keeps
for every forest whether it stated its own size and depth.

What goes. `$setLeafPrior(forests = )`, built on main and never released, is removed; `$setLeafPrior` and the
planned `$setTreePrior` take only the fitting function's level, and the per-forest lists planned for
`$setTreePrior` are never built. Nothing changes for bartCause, stan4bart, treatSens or bairrtt.

Changing a coefficient between held and drawn, or between the normal and the t, needs the engine work you ruled
for in dec-B511; that item stays its own and adds it to this door afterwards.

Six questions for you, below:
1. The door's name: I recommend `setForest`.
2. On a multinomial sampler, whose category forests creation can only state together: I recommend the door take
   them together with the forest left out, `s$setForest(power = 3)`, and that a name on the one forest entry
   creation takes be refused, since it names nothing. Labeling those forests by the response's levels in the
   readers is a separate, older gap, listed apart.
3. How to return a forest to the level above: I recommend `NULL`, as `options()` and `modifyList()` use it, a
   name left out meaning no change.
4. Whether a variance forest that states no depth follows a live change of the depth it took at creation: I
   recommend it follows, which is how your ruling's words read.
5. A forest that only the fitting function's setting reaches (a one-forest sampler's forest, for one), where the
   fitting function stated a value: I recommend refusing a forest-level statement there, as creation does, and
   naming the fitting function's door; to move the setting to the forest, the caller first clears the fitting
   function's (`$setLeafPrior(NULL)`, `$setTreePrior(power = NULL)` under question 3's recommendation).
6. A forest that shares the total with other forests with no basis, stating its own size live: I recommend
   refusing it for 1.0-0, as the shared-total plan does, because it would move the other forests' sizes too;
   taking it is about 420 lines of engine work and can be added later without breaking a call.

## Goal

A live sampler has `$setForest(forest, ...)`, `$getForest(forest)`, `$setVarianceForest(...)` and
`$getVarianceForest()`. `$setForest` changes on the forest named by label or position only the names it is given,
among `sd`, `power`, `base`, `amplitude`, `proposal.probs`, `n.perturb.cuts` and `treeShift`, each meaning what
it means in `forest()` on that forest at creation; refuses structure by name; refuses a statement where creation
refuses the same; and changes nothing where nothing stated differs from the statement in force. `$getForest`
returns that statement, so `do.call(s$setForest, c(list(f), s$getForest(f)))` changes nothing, draws included,
and a statement saved earlier restores. The sampler records, from creation on, the fitting function's level apart
from every forest's own, so each level can be read, returned to and checked as creation checks it. Every change
is recorded where re-creation reads it, so a copy and a reload run what the live sampler runs.
`$setLeafPrior(forests = )` is gone.

## Constraints

- Out of scope: a changed coefficient hold or prior (coefficient-prior-live adds it to this door);
  `$setTreePrior` and an unstated forest following it (setmodel-parameters' slice 4); `$setLeafPrior`'s
  fitting-function level on several forests and `getLeafPrior`'s record fields (leaf-prior-writer-by-level,
  forests-shared-total); per-category multinomial statements; labels for a multinomial sampler's category
  forests in the readers
  ([A separate item: multinomial reader labels](#a-separate-item-multinomial-reader-labels)); per-forest handles.
- No change to inst/include/dbarts/dbarts.h, the state format or any consumer.
- The engine sees numbers, never a record or an R name: the bridge converts.

## What exists today

Run on a c55a1477 build (scratch/fd-plan/p1.R to p4.R; the critique's c1.R to c9.R):

- The per-forest writers are [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd)'s `forests` list,
  which restates `sd` alone (anything else "is fixed at creation"), refuses itself on one forest, and is used at
  45 call sites in 9 tinytest files; and `setForestBasis` and `setForestWeights`, by label or position through
  [`samplerForestIndex`](../../R/bartcore.R). The forest-indexed readers are `getLeafPrior`, `getK`,
  `getForestFits`, `getForestAmplitudes`, `getForestVariableCounts`, `getTrees(forest = )` and `predictForests`.
- The forests' record, `attr(control, "bartcore.forests")`, holds per forest the resolved tree count, base,
  power and drawn flag ([`forestParams`](../../R/model.R)), the stated `sd` or NA
  ([`forestStatements`](../../R/model.R)), vars, labels, the hold, interactions, blocks and basis terms. A
  one-forest sampler has no record. The move settings ride `attr(control, "bartcore.moves")`
  ([`forestMoveRecord`](../../R/model.R)), NULL where a forest states none; the variance forest's base and power
  ride `attr(control, "bartcore.variance")`, resolved.
- The levels are not kept apart (c8, c3): `tree.prior = cgm(power = 1.5)` and `forest(power = 1.5)` on one
  forest, the fitting function's power and forest 1's on two forests, and `leaf.prior = normal(sd = 2)` and
  `forest(sd = 2)` on one forest each give identical models and records; only a `bartcore.forestsDeclared` flag
  differs, naming no setting. The bridge reads forest 1's base and power from the model's tree prior on every
  sampler ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)), the others' from the record.
- The engine keeps a tree prior on every forest and the variance forest, caching nothing derived from base or
  power ([`CGMTreePrior`](../../src/bartcore/model.hpp)); nothing installs one live except `$setModel`, which
  writes the first mean forest's. Live sd and move installs exist:
  [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setMoveSettings`](../../src/R_interface_bartcore.cpp).
- Creation: a lone forest's `sd` beside any named `leaf.prior` is refused, "states its size twice", a bare
  `normal()` and a named `k = NULL` included (c1, c5); `amplitude = fixed()` on a lone forest is refused; a
  variance forest beside several forests is refused; an unstated variance forest takes the mean forest's
  resolved depth, a lone forest's own included (c9; tree-prior-flat's D5 keeps it). Main takes
  `$setLeafPrior(normal(k = 3))` on a `forest(sd = 2)` sampler (c3).
- A multinomial sampler (matrix interface only) takes one `forests` entry for every category forest; a name on
  that entry is taken and ignored (c6); the category forests have no labels (`getLeafPrior("q")` is refused),
  the engine sharing one [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp) across them.
- `getLeafPrior(2)` on a heteroscedastic sampler is out of range; no reader shows a forest's power or base.

## Design

### The two levels, recorded

The sampler keeps, from creation on:

- the fitting function's level, apart from any forest's: `power`, `base` and the leaf prior, each its resolved
  value always, with a flag beside it saying whether its caller stated it. The value is kept whether stated or
  not, so a sampler built by `new("dbartsSampler", control, model, data)` whose model carries power 3, and its
  reload, return a forest to 3 and not to the default (c10). A named NULL is recorded as not stated, as NULL
  means a name not given (dec-B140); whether creation should stop counting a named NULL as stated (c5) is held
  apart, open call 3's background;
- for every forest, a lone forest and every multinomial model's one entry included, whether it stated its own
  `sd`, `power` and `base`, beside the values it runs under, a named NULL counted as not stated; the move
  settings' record already holds NULL where
  unstated, the hold and the coefficient's df where stated;
- for the variance forest, whether it stated `power` and `base`.

What each forest runs under stays where the bridge reads it, so no creation read moves: forest 1's base and
power in the model's tree prior, the others' in the record. A write to forest 1's depth therefore writes the
model's tree prior too, on every sampler, while the bridge reads forest 1 from it, and the forest-1 cases of the
tests cover a forest 1 of two as well as a lone one.

Where it lives. The fields are slots of the model's forest and variance records that control-migration's slice 1a
creates (its Homes table), written in [`resolveSamplerSpec`](../../R/spec.R), the one place where which level
stated what is still known, and they are carried by this plan's slice 2. Not slice 1a: that slice is a move, held
bitwise neutral by changing no reading, and fields nothing reads could only be pinned there by reading them back;
the door is their first reader, so they are built and tested with it. Not slice 1b: its structure check compares
models, and these fields are statedness of parameters, which it must not compare. A row for control-migration's
writer is in [Interaction with the other pending plans, and order](#interaction-with-the-other-pending-plans-and-order).

### The arguments, by kind of forest

"Unstated" is the level above: the fitting function's for a forest with no basis (dec-B274, dec-B509), the
kind's defaults for a forest with a basis, the control's for a move setting (dec-B477). Every value is checked by
the validator creation uses for it.

| name | no basis, lone forest | no basis beside forests with a basis | basis: factor | basis: numeric | sharer (several with no basis) | multinomial category forests | variance forest |
|---|---|---|---|---|---|---|---|
| `sd` | the forest's size, held; open call 5 where the fitting function stated a leaf prior | its size, held after forests-shared-total's flip; before it, its coefficient's prior median | the level difference's size, held, or its prior median where drawn | the size per unit of the column, held or prior median; unstated is the default from the column's recorded scale | open call 6 | refused, as at creation (the softmax map sets the scale) | not a name here |
| `power`, `base` | the forest's depth; open call 5 where the fitting function stated it; an unstated variance forest follows (open call 4) | the forest's depth; open call 5 likewise | its depth; unstated 3 and 0.25 | as factor | its depth; unstated the fitting function's | one statement for all K, open call 2 | its depth; unstated the mean forest's (open call 4) |
| `amplitude` | `fixed()` taken as the default restated after forests-shared-total's flip (dec-B513), refused before it, as creation; `student()` refused | `fixed()`, `student(df)`; a change waits on coefficient-prior-live | `fixed()`, `student(df)`; as left | as left | `fixed()` taken as the default restated; `student()` refused among two or more, as creation | not a name | not a name |
| `proposal.probs`, `n.perturb.cuts` | the forest's own; unstated, the control's | as left | as left | as left | as left | one statement for all K | the variance forest's own |
| `treeShift` | as above | as above | as above | as above | as above | one statement for all K | not a name (dec-B500) |

Held against drawn follows the coefficient's hold, as at creation: `sd` is the exact size where held and the
prior median where drawn, under the normal or under `student(df)` (dec-B510). On a lone forest a stated `sd`
replaces the fitting function's leaf prior for that forest, a drawn k included, and creation's twin is
`y ~ forest(x, sd = v)`.

### Forests one statement reaches from two levels

A lone forest, a multinomial sampler's category forests, and any forest with no basis that is the last one a
stated fitting-function setting reaches take each such setting from one of two levels, and creation refuses it
stated at both ("states its size twice"; tree-prior-flat's D5 for `power` and `base`). Live, under open call 5's
recommendation, the same refusal holds in both directions:

- `$setForest(f, sd = v)` (or `power`, `base`) is refused where the fitting function states that setting and no
  other forest would take it, naming `$setLeafPrior(normal(sd = v))` or `$setTreePrior(power = )` for the change
  and, for a move to the forest's own level, the spelling that clears the fitting function's statement: under
  open call 3's (a) `$setTreePrior(power = NULL)` (slice 4) and `$setLeafPrior(NULL)`, each returning that
  level to its default recorded as not stated; a bare `normal()` stays a statement of the default, as at
  creation (c1);
- `$setLeafPrior(...)`, `$setTreePrior(...)` and `$setModel`'s tree prior are refused where every forest they
  would reach states its own, naming `$setForest(f, sd = NULL)`; so the forest record is never left stale against
  the engine (c3's case today). On a one-forest sampler `$setModel(s$model)` restates what is in force and is
  taken; `$setModel` with a changed `power` where the forest states its own is refused, which main takes for a
  `forest(x, power = )` term (never released).

The refusal goes into the writers this item touches (`$setLeafPrior`, `$setModel`) and into slice 4's
`$setTreePrior`, whose collision rule (tree-prior-flat's) already reads that way.

### Parameter and structure

Parameters, taken: `sd`, `power`, `base`, the coefficient's hold and prior (dec-B254, dec-B511), and the move
settings, which are run settings and not model (dec-B477). Structure, refused by name whatever its value, before
anything is written (dec-B254, dec-B518): `vars`, `n.trees` (dec-B247), `interactions`, `blocks`: "'n.trees' is
fixed at creation: a forest's tree count is part of the model; make a new sampler". `basis` is refused naming
`$setForestBasis`, the basis being data. On the variance forest, `vars` and `n.trees` are refused alike. A
`$setVarianceForest` on a sampler with no variance forest is refused: whether there is one is structure.

### The reader and its round trip

`$getForest(forest)` returns a plain named list in `$setForest`'s argument order (`sd`, `power`, `base`,
`amplitude`, `proposal.probs`, `n.perturb.cuts`, `treeShift`): the forest's statement in force, each entry what
it states and NULL where it states none, as `forest()` objects hold NULL for a name not given. An attribute
`in.force` carries every value the forest runs under, in `forest()`'s names, with its label, what it multiplies,
its tree count and its columns: the depth a forest takes from the fitting function is read there. Giving the
statement back changes nothing, statedness included, because a stated entry restates its own value and a NULL
entry returns the forest to the level it is already at (open call 3). Because the fitting function's level is
recorded apart, a saved statement restores after later changes, as `options(op)` does: after
`op <- s$getForest(1)` on `dbarts(y ~ x)`, then `$setForest(1, power = 4)`, `do.call(s$setForest, c(list(1),
op))` returns forest 1 to the fitting function's power, which the record says was never stated (and under open
call 5's (b), where a forest-level write may cover a stated fitting-function value, returns it to that value).
`$getVarianceForest()` is the same for `power`, `base`,
`proposal.probs` and `n.perturb.cuts`. The setters return NULL invisibly, as every other sampler setter does.

`getLeafPrior` keeps its fields (dec-B515's reader is leaf-prior-writer-by-level's); its help stops saying its
`forest(sd = v)` goes back into `$setLeafPrior(forests = )` and says its `sd` is the one `$setForest` states.

### The forest argument

A label or a position (dec-B276), read by [`samplerForestIndex`](../../R/bartcore.R); a one-forest sampler's is
1 or `"forest1"`, or left out. Under open call 2's recommendation a multinomial sampler's category forests are
reached together with `forest` left out. Elsewhere it must be given. `"variance"` is an ordinary label for a mean
forest, so the variance forest has its own methods.

### How a call is carried out

1. Resolve the forest. Read the statement in force and the fitting function's level from the records.
2. Merge: a name given replaces its entry; a name left out keeps it (`missing()`, the formals carrying no
   defaults); under open call 3's recommendation an explicit NULL clears it, the forest then taking the level
   above as it stands now (the fitting function's current value, not creation's).
3. Refuse structure by name. Check the merged statement, not the names given, as creation checks it: the per-name
   validators ([`validateForestSd`](../../R/model.R), [`validateForestAmplitude`](../../R/model.R), the base,
   power and move checks) and the cross-forest and cross-level rules creation applies (two levels on one forest,
   held and drawn on one forest, `student()` among sharers, the multinomial refusals of
   [`refuseMultinomialLeafPrior`](../../R/dbarts.R)), factored out of [`resolveSamplerSpec`](../../R/spec.R)
   where they are inline so that one function serves both doors.
4. If the merged statement equals the one in force, return: no engine call, no record write.
5. Install each changed part, every check done first: `sd` through
   [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp), or on a lone forest through the existing
   single-forest leaf-prior install; `power` and `base` through slice 1's entry; move settings by rewriting the
   move record and calling [`bartcore_setMoveSettings`](../../src/R_interface_bartcore.cpp), as `$setControl`'s
   move install does; under open call 4's recommendation, an unstated variance forest's depth with a lone forest's.
   A failure after a write rolls back to the statement read in step 1, as `$setControl` rolls its control back.
6. Record: the forest record as creation with the merged statement would write it, field for field, and the
   model's tree prior for forest 1 on every sampler and its leaf prior on a lone forest, so `identical()` holds
   between the live sampler's records and a creation twin's.

### What goes and what is never built

- `$setLeafPrior(forests = )` goes with [`resolveForestSpreads`](../../R/dbarts.R) and
  [`writeForestSpreads`](../../R/dbarts.R); `forests = ` becomes R's unused-argument error, as
  `forest(leaf.prior = )` did under dec-B508. It never shipped, so it takes no tombstone (dec-A02 keeps those
  for shipped names). `$setLeafPrior` keeps `leaf.prior` (dec-B515's first level, leaf-prior-writer-by-level's
  item), with the refusal above. Its 45 call sites are respelled.
- `$setTreePrior(forests = , variance = )` is never built: setmodel-parameters' slice 4 builds
  `$setTreePrior(power, base, split.probs)`, the fitting function's level only, as tree-prior-flat's amendment
  already has it, reading and writing the fitting function's level this item records.
- Per-forest handles stay possible later as sugar over the door (dec-B518).

### RNG class

- Inert: a call whose merged statement equals the one in force writes nothing and moves no draw, a round trip
  and a NULL on an unstated name included.
- The twin oracle, before any run, on two chains: a taken call's next draws are bitwise those of a sampler
  created on the same data and seed with the merged statement. Exempt by name: a handoff from held to drawn,
  where the drawn spread starts at the one in force rather than at creation's start (c4: `sd = NULL` on a lone
  probit forest that stated its `sd`, its k drawn from the spread the held size gave; the coefficient's in
  coefficient-prior-live), as setmodel-parameters' call 8 exempts its handoffs. From drawn to held is not
  exempt: it already draws as creation does (c4).
- After a run: the change acts from the next sweep; `copy()` (which copies the live chain, dec-B459) and a
  `saveRDS` and `readRDS` after `storeState` continue `identical()` to the live sampler.
- No call draws: every install restates a scale, a depth or a mixture.

### Copy, reload and state

Every change rides the records and the model's slots; re-creation reads them. A state carries no prior (dec-B195)
and installs across statements as today. No state block, no `stateFormatVersion` move, no change to
inst/include/dbarts/dbarts.h.

## Slices

### Slice 1: a forest's tree prior installed live (engine and bridge)

1. A `Sampler` method, behind a facade virtual, that checks a forest index (a mean forest or the variance forest)
   and a base in (0, 1) and positive power before writing any chain, then writes that forest's tree prior on
   every chain; nothing else moves. [`CGMTreePrior`](../../src/bartcore/model.hpp) caches nothing derived from
   either (the critique read it).
2. A bridge entry, `bartcore_setForestTreePrior(ptr, forest, base, power)`, forest -1 the variance forest; R
   calls it in slice 2. It reads explicit numbers, not records, so it does not depend on where the records live.
3. tests/cpp: the install on two chains and three forests, forest 2 alone moved; the variance forest's; the
   refusals writing nothing; a sampler installed before any sweep drawing as one constructed with the value
   (c4's case D already holds through `$setModel` for forest 1).

Built off main before control-migration's slice 1a lands, and taken from this branch by setmodel-parameters'
slice 4, whose item 1 needed the same install. Neutral: nothing in R reaches it until slice 2.

### Slice 2: the door and the levels' records (R), on control-migration's slice 1a

0. Before it starts, re-read slice 1a as built: where the forest and variance records live, which slots exist,
   whether the bridge still reads forest 1's depth from the model's tree prior; amend step 6 and item 1 to it.
1. The records of [The two levels, recorded](#the-two-levels-recorded), written at creation.
2. The four methods, their docstrings and help, by [How a call is carried out](#how-a-call-is-carried-out).
3. The writers' side of [Forests one statement reaches from two levels](#forests-one-statement-reaches-from-two-levels)
   in `$setLeafPrior` and `$setModel`.
4. `$setLeafPrior(forests = )` removed; the 45 call sites respelled to `$setForest`.
5. The cross-forest checks of step 3 factored out of [`resolveSamplerSpec`](../../R/spec.R) where they are
   inline.
6. Under open call 2's recommendation, a name on a multinomial model's one forest entry refused at creation.

Slice 2 builds on slice 1a's tip (control-migration's Order of landing, item 7) because the door reads and writes
the records slice 1a moves onto the model, a one-forest model's included, which main does not have.

### Slice 3: a sharer leaving or rejoining the share (under open call 6)

Under (a): `$setForest(f, sd = v)` on a sharer takes it out of the share with size v, and `sd = NULL` on a forest
with no basis that states one returns it, so the K that divide the total change (dec-B505: each sharer an equal
part, 1/K). The engine restates every sharer's leaf scale against the new K, the group's k pool and the probit
rescaling step's membership; the leaver's k is pinned at 1 as any forest stating sd. `sd = v, amplitude =
student(1)` in one call is taken, as `forest(sd = v, amplitude = student(1))` is. Built by whichever of this item
and forests-shared-total lands second, on the other's tip. Under (b), the recommendation: both are refused by
name, "this forest shares the total with forests ...; a forest's own size is stated at creation: make a new
sampler", ~40 lines with its pins, in slice 2 or by forests-shared-total, whichever lands second.

## A separate item: multinomial reader labels

Labeling a multinomial sampler's category forests by the response's levels in every forest-indexed reader
(`getLeafPrior("q")`, `getTrees(forest = "q")`, `getForestFits`), `forest<i>` still position i, is dec-B276's
debt ("named by label everywhere"), not dec-B518's, and is not built here. About 50 lines and tests (~120 in
all); recorded in TODO as multinomial-forest-labels. The maintainer may fold it into this item.

## Help and NEWS

- man/dbartsSampler-class.Rd: usage and value for the four methods, the argument table above in prose, the
  refusals, the round trip with `do.call`, that no call touches the trees; `setLeafPrior` loses `forests` and
  gains the two-level refusal; `setModel` the same; `getLeafPrior`'s sentence on `forest(sd = v)` changes; the
  "variance forest ... is not addressable" sentence points to `setVarianceForest`.
- man/forest.Rd and man/dbarts.Rd: their `$setLeafPrior(forests = )` sentences become `$setForest`.
- inst/NEWS.Rd: the 1.0-0 multi-forest item's "which `$setLeafPrior(forests = )` restates" becomes `$setForest`,
  and one item adds the four methods. Nothing about `forests = ` itself, which never reached a release.
- docs/design: forest-defaults-by-kind.md, forest-multiplier.md, multiplier-combiner.md and forest-sd-unit.md
  name `$setLeafPrior(forests = )`; each is respelled or marked as the history it is.
- tools/check-rc-codoc.R covers the four methods with no change to the script; their usage carries no defaults
  for the stated names, `updateState = NULL` last, as `setLeafPrior` has.

## Tests and gates

tinytest, a new test-forest-door.R (~1000 lines), every R-side twin on two chains:

- Each argument on each kind of the table: the twin oracle before a run (draws and records `identical()` to a
  creation twin); after a run, the change in force (the reader, `getK`, the next draws against a twin given the
  same state), `copy()` and a reload `identical()` to the live sampler; forest 1 of two as well as a lone forest.
- Inert: the round trip on every kind, a subset restated, a NULL on an unstated name: draws bitwise equal to no
  call, and a following `$setControl` or the fitting function's door still reaching the forest.
- Return to the level above: a stated `power` cleared, then reading the fitting function's value in force, set at
  creation and changed since, and on a sampler built by `new()` from a model carrying power 3, and its reload,
  landing on 3; `$setLeafPrior(NULL)` and `$setTreePrior(power = NULL)` clearing the fitting function's level so
  a forest may then state its own; a move setting cleared then following `$setControl`; a numeric forest's `sd`
  cleared to the recorded default, and to one derived over the rows in where creation stated it.
- A saved statement restored after two later changes.
- The two levels: both orders on a lone forest, `$setForest(1, sd = )` after `leaf.prior = normal(k = 3)` and
  `$setLeafPrior` after `forest(sd = 2)`, and the same for `power` with `$setModel`; the records agreeing with the
  engine after each.
- Refusals: a table of statements run at creation and at the door, the same refused or taken at both with the
  same reason; structure by name; a refused call with one good and one bad name leaving the reader, records and
  draws as they were.
- The variance forest: its own door; under open call 4's (a), an unstated one moving with the lone forest's
  `$setForest(1, power = )` and a stated one not, name by name: one stating only `power` keeps it and still
  follows the lone forest's `base` (c9).
- Labels: a forest selected by label where the formula and a list put it at different positions; the multinomial
  door with `forest` left out.
- The old spelling: on a c55a1477 build `$setLeafPrior(forests = )` and on the branch `$setForest(label, sd = )`,
  bitwise equal draws over the bcf, numeric and three-forest fixtures (scratch A/B, ~150 lines).

Slice 1: tests/cpp plain and under ASan/UBSan, `--preclean`, the mutation battery's `verify-anchors`; the
equivalence trio bitwise (55, 19, 11 identical against benchmarks/baselines/MANIFEST's baselines), the four
snapshot files unchanged on the reference build, the full suite. Slice 2: the above, the R-loaded sanitizer run
over test-forest-door.R (slice 1's entry is reachable only by `.Call`), the full suite once in one process, lint
against the slice's library, air, check-rc-codoc, check-win-drift, check-doc-freshness, `R CMD check --as-cran`
from a clean tarball, and the exact gates of .github/workflows/exact-gates.yaml with `quick`, the records
changing what a fit carries. No new exact gate: the twin oracle ties each new statement's posterior to a creation
the exact gates cover. Slice 3 under (a): tests/cpp sanitized; the twin oracle against creation with the forest
stating `sd`, at sharers of unequal counts (38 and 37); forest-prior-exact.R's sharer arms quick; the A/B neutral.

Mutants, each to fail a named pin:

1. The no-change test removed, so every call writes: the inert pins on statedness (a round trip that states the
   values in force stops the fitting function's door reaching the forest).
2. A name left out read as NULL: a call stating `power` alone clears the forest's `sd`.
3. NULL read as no change: the return-to-level pins.
4. The record written and the engine not: the twin oracle. The engine written and the record not: copy and reload.
5. The tree-prior install writing forest 1 whatever the index, or base and power swapped at the bridge: forest
   2's twin.
6. The move record rewritten without the move install: the twin.
7. A check moved after the first write: the atomic refusal pin.
8. Structure equal to the value in force taken: the refusal table.
9. Labels read as positions: the label pin.
10. `$setVarianceForest` writing mean forest 1: the variance twin.
11. The model's tree prior not written for forest 1, on a lone forest or on forest 1 of two: copy after
    `$setForest(1, power = )`.
12. A creation check left out of the shared function (student() among sharers, the lone forest's two sizes): the
    refusal table.
13. Checks applied to the names given rather than the merged statement: a call stating `amplitude = student(1)`
    on a forest whose stated `sd` it needs, and one clearing `sd` beside a held amplitude.
14. A stale record after `$setLeafPrior` or `$setModel` (the refusal removed): the two-level pins' records against
    the engine.
15. A NULL return reading creation's fitting-function value instead of the current one: the return pin after a
    change at the fitting function's level.
16. The fitting function's record keeping NULL where unstated instead of the resolved value: a direct-door
    sampler built with power 3, `$setForest(1, power = 4)`, then `power = NULL` landing on 2; and a named NULL
    recorded as stated.
17. An R-side twin run on one chain: the pins assert two chains, and a per-chain write to chain 0 alone fails them.
18. Slice 3 under (a): the other sharers' leaf scales left at 1/K after one leaves, or the leaver kept in the k
    pool: the unequal-count twin.

## Interaction with the other pending plans, and order

1. coefficient-prior-opt-in (acee2dbf) lands first. Its `student(df)` and the record's `amplitude.df` are what
   `amplitude` reads and gives back; a `student(3)` forest's round trip and `sd` restated keep its df (its pin
   "a restated sd keeps a forest's student() prior" moves to `$setForest`). Shared:
   [`validateForestAmplitude`](../../R/model.R), [`reportLeafPrior`](../../R/dbarts.R).
2. tree-prior-flat (fix round in progress) lands next. Its plain `power` and `base` are this door's names; its
   collision rule's message names `$setForest`; its D5 (an unstated variance forest takes the mean forest's
   resolved depth, a lone forest's own included) is the creation rule open call 4's (a) carries live. Its table
   row "live door, one forest: `$setForest(1, power = 2)` changes that forest only" disagrees with open call 5's
   (a) where the fitting function stated `power`, and with open call 4's (a) where an unstated variance forest is
   present. Its `$setTreePrior(power = NULL, ...)` reads NULL as no change, where open call 3's (a) reads it as a
   return to the default. Records step for its writer, after the maintainer rules.
3. forest-after-merge (f251f037, built) lands whenever it is green. Its derivation of a numeric forest's scale
   over the rows in is what `sd = NULL` uses on a numeric forest that stated `sd` at creation (Calls made 8).
4. control-migration (93f95ac0): slice 2 is built on its slice 1a's tip and this branch lands after it merges.
   Its row "`$getForest` reads the stated count from the record, so giving it back changes nothing" holds: the
   count is in the reader's `in.force` part and `n.trees` is never a statement name. Collision row for its
   writer: this item's slice 2 adds slots to 1a's forest and variance records (each forest's `sd`, `power` and
   `base` stated or not) and a record of the fitting function's level as stated, written in
   `resolveSamplerSpec`; 1a's record classes and validity take them without change of meaning; 1b's structure
   check must not compare them, being parameters' statedness; and the bridge keeps reading forest 1's depth from
   the model's tree prior unless 1a changes that, in which case slice 2's step 6 follows 1a. Slice 1 is built
   before slice 1a, off main.
5. forests-shared-total (f93810e2): its D10 refusal of a sharer's restated `sd` stands under open call 6's (b),
   moved from `$setLeafPrior(forests = )` (its R item 6) to `$setForest`; under (a) it is superseded and slice 3
   is built by whichever lands second. Both edit [`resolveForestSpreads`](../../R/dbarts.R), which this item
   deletes: if forests-shared-total lands first, its edit there moves to `$setForest`. Records step for its writer
   after the maintainer rules.
6. setmodel-parameters (abdcea31): slice 4's open call 1 is answered by dec-B518; its item 1 per-forest install is
   this item's slice 1, which it takes (this item lands first; if slice 4 is ready first it builds the install as
   planned and slice 1 drops); its item 4 door is the fitting function's level only, reading and writing the
   level this item records, with open call 5's refusal; its item 5 refusal "naming `$setLeafPrior(forests = )`"
   names `$setForest`; its open call 2 is open call 4 here. Records step for its writer.
7. coefficient-prior-live (dec-B511) follows this item and coefficient-prior-opt-in, adding a changed `amplitude`
   to `$setForest` (Calls made 9).
8. leaf-prior-writer-by-level (dec-B515's remaining half) and forest-list-by-label (predict's half of dec-B516):
   independent; textual meetings in R/dbarts.R only.

Order of landing: coefficient-prior-opt-in; tree-prior-flat; forest-after-merge when green; control-migration;
this item (slice 1 built early, slice 2 on slice 1a, merged once after control-migration); forests-shared-total
and slice 3 by whichever is second; setmodel-parameters' slice 4; coefficient-prior-live. Engine features merge
serially.

## Cost

- Four methods, a reader attribute and records of the fitting function's level; `$setLeafPrior` loses an
  argument main has.
- An explicit NULL and a name left out mean different things, under open call 3's recommendation.
- A refusal list that must track creation's; the shared check function and the refusal table carry it.
- Under open call 5's recommendation, two writers main lets overwrite each other refuse instead.
- Until coefficient-prior-live lands, a changed `amplitude` is refused by name on main (never in a release).

## Calls made

1. The reader returns the statement, NULL where unstated, and the values in force as an attribute, plain list, no
   class. Alternatives: values in force alone, which state everything when given back; a `dbartsForest` object,
   which `forests = ` would take with its basis dropped.
2. Structure names are refused whatever their value; the reader carries none in the statement, so the round trip
   needs no equality rule. Alternative: take a structure name equal to the value in force.
3. All checks before any write, on the merged statement; a failure after one rolls back. Alternative: write per
   name as checked.
4. A call that changes nothing writes nothing. Alternative: write and rely on bitwise-inert installs, which still
   rewrite statedness.
5. One check function serves creation and the door, pinned by a two-door refusal table. Alternative: the door's
   own checks, which drift.
6. `treeShift` is taken, the third move setting of dec-B477 and `forest()`'s name; not on the variance forest
   (dec-B500). In the summary as a note the maintainer may veto. Alternative: leave it to `$setControl` alone.
7. `updateState = NULL` last, as `setLeafPrior`; the setters return NULL invisibly, as every sampler setter.
   Alternative: the previous statement, as `options()`, which `$getForest` already gives.
8. `sd = NULL` on a numeric forest returns the default from the column's recorded scale, else from the column in
   force over the rows in at the call, recorded, as creation would derive it (dec-B514's principle, dec-B303's
   rows). Alternative: refuse where no scale is recorded.
9. coefficient-prior-live stays its own item, adding to this door; until then a changed `amplitude` is refused and
   one equal to the statement in force is taken. Alternative: fold its ~400 lines in as a slice, delaying this
   item behind it.
10. Slice 1's install is built here and taken by slice 4. Alternative: wait for slice 4, which waits on
    forests-shared-total.
11. The levels' records are slots of slice 1a's records built by this item's slice 2; forest 1's depth stays read
    from the model's tree prior, which every write to it updates. Alternatives: in slice 1a, where nothing reads
    them; the bridge reading forest 1 from the record, a creation read moved.
12. Removal of `forests = ` is R's unused-argument error, as dec-B508's removal was. Alternative: a message naming
    `$setForest`, for a spelling no release had.

## Open calls

1. The door's name. Background: the sampler's methods are `set*` and `get*`, and `setForestBasis` and
   `setForestWeights` already change a forest's data; in a simulated reading one reader feared `setForest` resets
   the trees.
   - (a) Recommended: `setForest` and `getForest`, `setVarianceForest` and `getVarianceForest`. Cost: the help
     must say the trees are left alone.
   - (b) `updateForest`, after `update()`'s "change what you state". Cost: in a sampler "update" reads as a draw
     (the package's own `updateState`), and no sampler method is named so.
   - (c) `setForestPrior`. Cost: it carries move settings, which are not prior.
2. The door on a multinomial sampler. Background: creation states the category forests together, one `forests`
   entry for all K, and silently ignores a name on that entry (`forests = list(q = forest(power = 3))` gives power
   3 to every category); the engine shares one specification, so no per-category change can be made.
   - (a) Recommended: the door takes them together with `forest` left out, `s$setForest(power = 3)` and
     `s$getForest()`, a category position refused naming that form; a name on the one creation entry names
     nothing and is refused at creation by name ("a multinomial model's one forest entry states every category
     forest and takes no name"). Cost: a creation call main takes is refused (never released); ~30 lines.
   - (b) Any category position reaches all K. Cost: a call naming one category changes every one.
   - (c) No door on a multinomial sampler: its forests change through the fitting function's doors only. Cost: a
     setting stated on its forest entry at creation can then never change live.
   - (d) Per-category statements at creation and live. Cost: a creation feature and engine specifications per
     category; after 1.0-0.
3. The spelling that returns a forest to the level above (it covers dec-B511's open spelling, a coefficient back
   to the normal or a hold released). Background: a name left out changes nothing, so another spelling is needed.
   At creation a NULL named at the fitting function counts as stated (`bart(k = NULL)` beside a lone forest's
   `sd` is refused, c5), against dec-B140's NULL for absent; tree-prior-flat's planned `$setTreePrior(power =
   NULL)` reads NULL as no change. Whether creation keeps counting a named NULL as stated is a separate decision
   held for the maintainer (the coordinator's decisions-pending #24, which also covers tree-prior-flat's
   `power = NULL`); this call depends on it and does not settle it. The live records count a named NULL as not
   stated either way.
   - (a) Recommended: an explicit NULL, `s$setForest("tau", sd = NULL)`, as `options(x = NULL)` and
     `modifyList(x, list(a = NULL))` clear an entry and as NULL means a name not stated in `forest()`; the reader
     writes NULL for an unstated name, so a saved statement restores. Cost: NULL and a name left out differ, which
     few R functions do; slice 4's `$setTreePrior` reads NULL the same way, a return to the default as dec-B514
     reads a bare `normal()`, and `$setLeafPrior(NULL)` clears the fitting function's leaf prior (open call 5);
     and it needs the fitting function's level recorded apart from the forests' (built
     in slice 2, ~150 lines of it), since a NULL returns to that level as it stands.
   - (b) A constructor, `sd = inherit()`. Cost: a new exported name, the reader returning it; the same records.
   - (c) An argument naming what to clear, `unstate = "sd"`. Cost: the reader cannot carry it, so a saved
     statement cannot be restored.
4. Held decision 6 re-read under dec-B518: does a variance forest that states no depth follow a live change of
   the depth it took at creation? Background: at creation it takes the mean forest's resolved depth, the fitting
   function's or a lone forest's own (dec-B512; tree-prior-flat's D5; c9); `$setModel` today leaves it at
   creation's value. dec-B518's text, the sampler's doors "each reaching the forests that state none of their
   own", reads as (a); this plan keeps it open because held decision 6 is.
   - (a) Recommended: it follows, through the fitting function's doors (`$setTreePrior`, `$setModel`, slice 4) and
     through `$setForest(1, power = , base = )` on a lone forest whose own depth it took at creation (this item),
     so the twin oracle holds. Cost: `$setForest(1, ...)` moves a second forest, the unstated variance one, which
     the help says; one call main honors changes, a heteroscedastic `$setModel` with a new tree prior,
     posterior-changing for that call, with a new arm in heteroscedastic-exact.R.
   - (b) It follows the fitting function's doors only. Cost: after `$setForest(1, power = )` the sampler differs
     from creation's twin unless creation stops passing a lone forest's own depth on, which moves a main call.
   - (c) It keeps creation's value and changes only through `$setVarianceForest`. Cost: unlike every other
     unstated setting, it follows at creation and not after.
5. A forest one setting reaches from two levels (a lone forest; a multinomial sampler's forests; the last forest
   with no basis a stated fitting-function setting reaches). Background: creation refuses the setting stated at
   both levels; main's `$setLeafPrior(normal(k = 3))` on a `forest(sd = 2)` sampler is taken and, after this
   item, would leave the forest's stated `sd` in the record while the engine runs k = 3 (c3).
   - (a) Recommended: refuse as creation does, in both directions: `$setForest` naming the fitting function's door
     where the fitting function states the setting, and `$setLeafPrior`, `$setTreePrior` and `$setModel` naming
     `$setForest(f, ... = NULL)` where the forest states its own. A setting moves to the forest's level by first
     clearing the fitting function's: `$setTreePrior(power = NULL)` and `$setLeafPrior(NULL)` under open call 3's
     (a), `$setLeafPrior` taking NULL for it (~10 lines; main refuses it, "'leaf.prior' must be a leaf prior
     specification", c10); under open call 3's (b) or (c) that call's spelling. Cost: two calls to move a
     setting between levels;
     `$setModel` with a changed power refused on a `forest(x, power = )` sampler, which main takes; needs the
     levels recorded (open call 3's cost, shared).
   - (b) The later write wins: a forest-level statement clears the fitting function's for that forest, and a
     fitting-function statement clears the forest's own. Cost: a live sampler reaches states creation refuses,
     and a NULL afterwards lands on the default rather than the value written over.
6. A forest that shares the total, stating its own `sd` live. Background: at creation a forest with no basis that
   states `sd` leaves the share, and the forests left share the total among fewer (dec-B505); dec-B518 says the
   door refuses where creation refuses, and also that it changes "only what it states", while leaving the share
   moves every other sharer's size. forests-shared-total planned a refusal.
   - (a) Take it, and `sd = NULL` to rejoin (slice 3). Cost: ~420 lines more than (b) and ~8 agent-hours, engine
     tier with its own review; a call that changes other forests' sizes.
   - (b) Recommended: refuse both for 1.0-0, naming a new sampler, ~40 lines. Cost: one statement creation takes
     is refused live; a caller who wants a forest's own size states it at creation. Taking it later changes no
     call.

## Not verified

- control-migration's slice 1a is unbuilt: where its records live, its slots and whether the bridge keeps reading
  forest 1's depth from the model's tree prior are from its plan; slice 2's step 0 re-reads them.
- The twin oracle's exemptions beyond the held-to-drawn handoff were checked on the cases in c4 only; creation may
  consume random numbers by setting in a way an install cannot reproduce for a kind not run there.
- forests-shared-total's engine for slice 3 is read from its plan and slice 1, not prototyped; its size is an
  estimate.
- The design pass's timing (about 100 microseconds a call by label, a sweep about 350) is for `setForestBasis`;
  `$setForest`'s checks were not timed.
