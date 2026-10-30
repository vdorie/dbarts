# forest-door: a live sampler changes one forest through $setForest, by label, in forest()'s names

Status: PLANNED 2026-10-30 on c55a1477 (dec-B518). Four calls are held for the maintainer
([Open calls](#open-calls)); the blind critique has not run. Rulings: dec-B518 (the door), with dec-B254
(parameters, never structure), dec-B276 (a forest by label or position), dec-B477 and dec-A239 (per-forest move
settings), dec-B511 (the coefficient's prior live), dec-B512 (plain power and base), dec-B514 and dec-B515 (a
setter means what creation means; dec-B515 revised in part by dec-B518), dec-B516 (revised in part by dec-B518),
dec-B129 (leaf where it means a leaf), dec-B247 (a live count does not change), dec-B459 (copy of a live chain),
dec-B500 (no `treeShift` on the variance forest), dec-B505 and dec-B509 (forests that share a total), dec-B513.

agent: one opus implementer (effort high) for slices 1 and 2, one opus review over both that runs the mutants
below, the same reviewer checking the fix round unless the fixes change engine behavior; slice 3 one opus
implementer and one opus review. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)): slice 1
adds an engine install and slice 2 a door that changes priors between sweeps; slice 3 edits the shared-total
group in the engine.
rng: NEUTRAL, bit for bit, for every call main accepts. The one main spelling this item removes,
`$setLeafPrior(forests = )`, is respelled `$setForest(forest, sd = )` and draws bitwise as the old spelling did
(an A/B against a c55a1477 build). Every new call is held to the twin oracle: before any run its next draws are
those of a sampler created with the same statement ([RNG class](#rng-class)). No call draws a random number.
window: for 1.0-0, after the merge to main (dec-B518 names none; nothing per-forest shipped in 0.9-34).
budget: ~2100 lines in the tree. Slice 1 ~250 (engine ~60, bridge ~60, tests/cpp ~130). Slice 2 ~1350 (R ~430
added and ~130 removed, help ~130, NEWS ~10, tinytest ~780). Slice 3 ~460 (engine ~100, bridge ~30, R ~40, help
~10, tests/cpp ~120, tinytest ~160). Records ~40. Agent time, estimated from comparable landings and not
measured: slices 1 and 2 ~14, their review and fix round ~5, slice 3 ~6 and its review ~3, gates ~3, rebases
onto control-migration ~3: ~34; plans here run 1.5 to 2 times low.

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
run, any accepted call gives exactly the draws a sampler created with the same statement gives. A call that
states only what is already in force touches nothing at all.

What the reader gives. The forest's own statement, one entry per name, empty where the forest states nothing of
its own and takes the level above (the fitting function's settings for a forest with no basis, its kind's
defaults for one with a basis, the control's for the move settings). Beside it, for reading, every value in
force, including the tree count, columns and depth. Today no reader shows a forest's depth.

What goes. `$setLeafPrior(forests = )`, built on main and never released, is removed; `$setLeafPrior` and the
planned `$setTreePrior` take only the fitting function's level, and the per-forest lists planned for
`$setTreePrior` are never built. Nothing changes for bartCause, stan4bart, treatSens or bairrtt: none calls a
per-forest writer, and bartCause's read of the first forest's leaf prior stays.

Two cases need more than this item. Changing a coefficient between held and drawn, or between the normal and the
t, needs the engine work you ruled for in dec-B511; that item stays its own and adds it to this door afterwards.
A forest that shares the total with other forests with no basis (the planned several-forest models) can state
its own `sd` at creation, which takes it out of the share; under your rule the door must take that too, so the
others' shares rebalance. That is a third slice, built by whichever of this item and the shared-total item lands
second. The shared-total plan had refused it; your rule overrides that.

Four questions for you, below: the door's name (I recommend `setForest`), labels for a multinomial sampler's
category forests (I recommend labeling them by the response's levels in every reader, and having the door take
the category forests together with the forest left out, since creation can only state them together), the
spelling that returns a forest to the level above (I recommend `NULL`, as `options()` and `modifyList()` use it,
with a name left out meaning no change), and whether a variance forest that states no depth follows a live
change of the fitting function's (I recommend it follows, as it does at creation).

## Goal

A live sampler has `$setForest(forest, ...)`, `$getForest(forest)`, `$setVarianceForest(...)` and
`$getVarianceForest()`. `$setForest` changes on the forest named by label or position only the names it is given,
among `sd`, `power`, `base`, `amplitude`, `proposal.probs`, `n.perturb.cuts` and `treeShift`, each meaning what
it means in `forest()` on that forest at creation; refuses structure by name; refuses a statement where creation
refuses the same; and changes nothing where nothing stated differs from the statement in force. `$getForest`
returns that statement, so `do.call(s$setForest, c(list(f), s$getForest(f)))` changes nothing, draws included.
Every change is recorded where re-creation reads it, so a copy and a reload run what the live sampler runs.
`$setLeafPrior(forests = )` is gone.

## Constraints

- Out of scope: a changed coefficient hold or prior (coefficient-prior-live adds it to this door);
  `$setTreePrior` and an unstated variance forest following it (setmodel-parameters' slice 4); `$setLeafPrior`'s
  fitting-function level on several forests and `getLeafPrior`'s record fields (leaf-prior-writer-by-level,
  forests-shared-total); per-category multinomial statements; per-forest handles.
- No change to inst/include/dbarts/dbarts.h, the state format or any consumer.
- The engine sees numbers, never a record or an R name: the bridge converts.

## What exists today

Run on a c55a1477 build (scratch/fd-plan/p1.R to p4.R):

- The per-forest writers are [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd)'s `forests` list,
  which restates `sd` alone (anything else "is fixed at creation"), refuses itself on one forest ("only a sampler
  whose forests carry amplitudes has" one), and is used at 45 call sites in 9 tinytest files; and
  `setForestBasis` and `setForestWeights`, by label or position through
  [`samplerForestIndex`](../../R/bartcore.R). The forest-indexed readers are `getLeafPrior`, `getK`,
  `getForestFits`, `getForestAmplitudes`, `getForestVariableCounts`, `getTrees(forest = )` and `predictForests`.
- The forests' record, `attr(control, "bartcore.forests")`, holds per forest the resolved tree count, base,
  power and drawn flag ([`forestParams`](../../R/model.R)), the stated `sd` or NA
  ([`forestStatements`](../../R/model.R)), vars, labels, the hold, interactions, blocks and basis terms. Whether
  `base` and `power` were stated is not recorded. A one-forest sampler has no record. The move settings ride
  `attr(control, "bartcore.moves")` ([`forestMoveRecord`](../../R/model.R)), NULL where a forest states none;
  the variance forest's base and power ride `attr(control, "bartcore.variance")`, resolved.
- The engine keeps a tree prior on every forest and the variance forest; nothing installs one live except
  `$setModel`, which writes the first mean forest's. Live sd and move installs exist:
  [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setMoveSettings`](../../src/R_interface_bartcore.cpp).
- Creation: a lone forest's `sd` beside `leaf.prior` is refused, "states its size twice"; `amplitude = fixed()`
  on a lone forest is refused; a variance forest beside several forests is refused. A multinomial sampler
  (matrix interface only) takes one `forests` entry for every category forest, whatever its name, and its
  category forests have no labels (`getLeafPrior("q")` is refused), the engine sharing one
  [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp) across them.
- `getLeafPrior(2)` on a heteroscedastic sampler is out of range; no reader shows a forest's power or base.

## Design

### The arguments, by kind of forest

"Unstated" is the level above: the fitting function's for a forest with no basis (dec-B274, dec-B509), the
kind's defaults for a forest with a basis, the control's for a move setting (dec-B477). Every value is checked by
the validator creation uses for it.

| name | no basis, lone forest | no basis beside forests with a basis | basis: factor | basis: numeric | sharer (several with no basis) | multinomial category forests | variance forest |
|---|---|---|---|---|---|---|---|
| `sd` | the forest's size, held; the fitting function's `leaf.prior` then reaches no forest | its size, held after forests-shared-total's flip; before it, its coefficient's prior median | the level difference's size, held, or its prior median where drawn | the size per unit of the column, held or prior median; unstated is the default from the column's recorded scale | leaves the share with its own size (slice 3) | refused, as at creation (the softmax map sets the scale) | not a name here |
| `power`, `base` | the forest's depth | the forest's depth | its depth; unstated 3 and 0.25 | as factor | its depth; unstated the fitting function's | one statement for all K | its depth; unstated the fitting function's (open call 4) |
| `amplitude` | `fixed()` taken as the default restated after forests-shared-total's flip (dec-B513), refused before it, as creation; `student()` refused | `fixed()`, `student(df)`; a change waits on coefficient-prior-live | `fixed()`, `student(df)`; as left | as left | `fixed()` taken as the default restated; `student()` refused among two or more, as creation | not a name | not a name |
| `proposal.probs`, `n.perturb.cuts` | the forest's own; unstated, the control's | as left | as left | as left | as left | one statement for all K | the variance forest's own |
| `treeShift` | as above | as above | as above | as above | as above | one statement for all K | not a name (dec-B500) |

Held against drawn follows the coefficient's hold, as at creation: `sd` is the exact size where held and the
prior median where drawn, under the normal or under `student(df)` (dec-B510). On a lone forest a stated `sd`
replaces the fitting function's leaf prior for that forest, a drawn k included, and creation's twin is
`y ~ forest(x, sd = v)`.

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
its tree count and its columns: the depth a forest takes from the fitting function is read there, which no
reader shows today. Giving the statement back changes nothing, statedness included, because a stated entry
restates its own value and a NULL entry returns the forest to the level it is already at (open call 3). It
restores a saved statement after later changes, as `options(op)` does. `$getVarianceForest()` is the same for
`power`, `base`, `proposal.probs` and `n.perturb.cuts`. `$setForest` and `$setVarianceForest` return the
statement from before the call, invisibly, as `options()` and `par()` do.

`getLeafPrior` keeps its fields (dec-B515's reader is leaf-prior-writer-by-level's); its help stops saying its
`forest(sd = v)` goes back into `$setLeafPrior(forests = )` and says its `sd` is the one `$setForest` states.

### The forest argument

A label or a position (dec-B276), read by [`samplerForestIndex`](../../R/bartcore.R); a one-forest sampler's is
1 or `"forest1"`. Under open call 2's recommendation `forest` may be left out where the sampler has one forest
statement, a one-forest sampler or a multinomial one, and must be given elsewhere. `"variance"` is an ordinary
label for a mean forest, so the variance forest has its own methods.

### How a call is carried out

1. Resolve the forest. Read the statement in force from the record.
2. Merge: a name given replaces its entry; a name left out keeps it (`missing()`, the formals carrying no
   defaults); under open call 3's recommendation an explicit NULL clears it.
3. Refuse structure by name. Check the merged statement as creation checks it: the per-name validators
   ([`validateForestSd`](../../R/model.R), [`validateForestAmplitude`](../../R/model.R), the base, power and move
   checks) and the cross-forest rules creation applies (a lone forest's `sd` beside a stated `leaf.prior`, held
   and drawn on one forest, `student()` among sharers, the multinomial refusals of
   [`refuseMultinomialLeafPrior`](../../R/dbarts.R)), factored out of [`resolveSamplerSpec`](../../R/spec.R)
   where they are inline so that one function serves both doors.
4. If the merged statement equals the one in force, return: no engine call, no record write.
5. Install each changed part, every check done first: `sd` through
   [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp), or on a lone forest through the existing
   single-forest leaf-prior install; `power` and `base` through slice 1's entry; move settings by rewriting the
   move record and calling [`bartcore_setMoveSettings`](../../src/R_interface_bartcore.cpp), as `$setControl`'s
   move install does. A failure after a write rolls back to the statement read in step 1, as `$setControl` rolls
   its control back.
6. Record: the forest record is written as creation with the merged statement would write it, field for field,
   and on a one-forest sampler the model's leaf prior and tree prior slots too, which there are that forest's
   (tree-prior-flat's D1), so `identical()` holds between the live sampler's records and a creation twin's.

### What goes and what is never built

- `$setLeafPrior(forests = )` goes with [`resolveForestSpreads`](../../R/dbarts.R) and
  [`writeForestSpreads`](../../R/dbarts.R); `forests = ` becomes R's unused-argument error, as
  `forest(leaf.prior = )` did under dec-B508. It never shipped, so it takes no tombstone (dec-A02 keeps those
  for shipped names). `$setLeafPrior` keeps `leaf.prior` (dec-B515's first level, leaf-prior-writer-by-level's
  item). Its 45 call sites are respelled.
- `$setTreePrior(forests = , variance = )` is never built: setmodel-parameters' slice 4 builds
  `$setTreePrior(power, base, split.probs)`, the fitting function's level only, as tree-prior-flat's amendment
  already has it. Where every forest a fitting-function statement could reach states its own, that statement is
  refused naming `$setForest(f, power = NULL)` (tree-prior-flat's collision rule, respelled).
- Per-forest handles stay possible later as sugar over the door (dec-B518).

### RNG class

- Inert: a call whose merged statement equals the one in force writes nothing and moves no draw, a round trip
  and a NULL on an unstated name included.
- The twin oracle, before any run: a taken call's next draws are bitwise those of a sampler created on the same
  data and seed with the merged statement. Exempt by name: a handoff between drawn and held, where the value in
  force carries over rather than creation's start (on a lone forest a drawn k replaced by a stated `sd`; the
  coefficient's handoffs in coefficient-prior-live), as setmodel-parameters' call 8 exempts its handoffs.
- After a run: the change acts from the next sweep; `copy()` (which copies the live chain, dec-B459) and a
  `saveRDS` and `readRDS` after `storeState` continue `identical()` to the live sampler.
- No call draws: every install restates a scale, a depth or a mixture.

### Copy, reload and state

Every change rides the forest record and, on one forest, the model; re-creation reads them. A state carries no
prior (dec-B195) and installs across statements as today. No state block, no `stateFormatVersion` move, no change
to inst/include/dbarts/dbarts.h.

## Slices

### Slice 1: a forest's tree prior installed live (engine and bridge)

1. A `Sampler` method, behind a facade virtual, that checks a forest index (a mean forest or the variance forest)
   and a base in (0, 1) and positive power before writing any chain, then writes that forest's tree prior on
   every chain; nothing else moves. Any value derived from base and power that the tree prior caches is
   refreshed.
2. A bridge entry, `bartcore_setForestTreePrior(ptr, forest, base, power)`, forest -1 the variance forest; R
   calls it in slice 2. It reads explicit numbers, not records, so it does not depend on where the records live.
3. tests/cpp: the install on two chains and three forests, forest 2 alone moved; the variance forest's; the
   refusals writing nothing; a sampler installed before any sweep drawing as one constructed with the value.

Built off main before control-migration's slice 1a lands, and taken from this branch by setmodel-parameters'
slice 4, whose item 1 needed the same install. Neutral: nothing in R reaches it until slice 2.

### Slice 2: the door (R), on control-migration's slice 1a

1. The four methods, their docstrings and help, by [How a call is carried out](#how-a-call-is-carried-out).
2. The record gains whether each forest's `base` and `power` and the variance forest's were stated, beside the
   resolved values the bridge reads. (slice 1a's record is assumed to carry the stated count and vars; if it
   carries statedness for base and power already, this item reads it.)
3. Under open call 2's recommendation, a multinomial sampler's category forests take labels from the response's
   levels in every forest-indexed reader, and the door takes them together.
4. `$setLeafPrior(forests = )` removed; the 45 call sites respelled to `$setForest`; man/dbarts.Rd and
   man/forest.Rd lose their mentions of it; NEWS, 1.0-0, one item for the door.
5. The cross-forest checks of step 3 factored out of [`resolveSamplerSpec`](../../R/spec.R) where they are
   inline.

Slice 2 builds on slice 1a's tip (control-migration's Order of landing, item 7) because the door reads and writes
the forest record that slice 1a moves onto the model, a one-forest model's included, which main does not have.
Built before slice 1a it would write the control's attributes and be rewritten at every site by the rebase.

### Slice 3: a sharer leaving or rejoining the share (engine, bridge, R)

`$setForest(f, sd = v)` on a sharer takes it out of the share with size v, and `sd = NULL` on a forest with no
basis that states one returns it, so the K that divide the total change (dec-B505: each sharer an equal part,
1/K). The engine restates every sharer's leaf scale against the new K, the group's k pool and the probit
rescaling step's membership; the leaver's k is pinned at 1 as any forest stating sd. A sharer stating
`amplitude = student()` while it shares stays refused, as at creation; `sd = v, amplitude = student(1)` in one
call is taken, as `forest(sd = v, amplitude = student(1))` is. Built by whichever of this item and
forests-shared-total lands second, on the other's tip.

## Help and NEWS

- man/dbartsSampler-class.Rd: usage and value for the four methods, the argument table above in prose, the
  refusals, the round trip with `do.call`, that no call touches the trees; `setLeafPrior` loses `forests`;
  `getLeafPrior`'s sentence on `forest(sd = v)` changes; the "variance forest ... is not addressable" sentence
  points to `setVarianceForest`.
- man/forest.Rd and man/dbarts.Rd: their `$setLeafPrior(forests = )` sentences become `$setForest`.
- inst/NEWS.Rd, 1.0-0: one item, `$setForest`, `$getForest`, `$setVarianceForest`, `$getVarianceForest`. Nothing
  about `forests = `, which never reached a release.
- tools/check-rc-codoc.R covers the four methods with no change to the script; their usage carries no defaults
  for the stated names, `updateState = NULL` last, as `setLeafPrior` has.

## Tests and gates

tinytest, a new test-forest-door.R (~550 lines):

- Each argument on each kind of the table: the twin oracle before a run (draws and records `identical()` to a
  creation twin); after a run, the change in force (the reader, `getK`, the next draws against a twin given the
  same state), `copy()` and a reload `identical()` to the live sampler.
- Inert: the round trip on every kind, a subset restated, a NULL on an unstated name: draws bitwise equal to no
  call, and a following `$setTreePrior` or `$setControl` still reaching the forest.
- Return to the level above: a stated `power` cleared then following `$setModel`'s tree prior on a one-forest
  sampler; a move setting cleared then following `$setControl`; a numeric forest's `sd` cleared to the recorded
  default, and to one derived over the rows in where creation stated it.
- Saved statement restored after two later changes, `options(op)` style; the invisible return.
- Refusals: a table of statements run at creation and at the door, the same refused or taken at both with the
  same reason; structure by name; a refused call with one good and one bad name leaving the reader, records and
  draws as they were.
- Labels: a forest selected by label where the formula and a list put it at different positions; multinomial
  labels in every reader (open call 2).
- The old spelling: on a c55a1477 build `$setLeafPrior(forests = )` and on the branch `$setForest(label, sd = )`,
  bitwise equal draws over the bcf, numeric and three-forest fixtures (scratch A/B, ~150 lines).

Slice 1: tests/cpp plain and under ASan/UBSan, `--preclean`, the mutation battery's `verify-anchors`; the
equivalence trio bitwise (55, 19, 11 identical against benchmarks/baselines/MANIFEST's baselines), the four
snapshot files unchanged on the reference build, the full suite. Slice 2: the above, the R-loaded sanitizer run
over test-forest-door.R (slice 1's entry is reachable only by `.Call`), the full suite once in one process, lint
against the slice's library, air, check-rc-codoc, check-win-drift, check-doc-freshness, `R CMD check --as-cran`
from a clean tarball, and the exact gates of .github/workflows/exact-gates.yaml with `quick`, the forest record
changing what a fit carries. No new exact gate: the twin oracle ties each new statement's posterior to a
creation the exact gates cover. Slice 3: tests/cpp sanitized; the twin oracle against creation with the forest
stating `sd`, at sharers of unequal counts (38 and 37); forest-prior-exact.R's sharer arms quick; the A/B neutral.

Mutants, each to fail a named pin:

1. The no-change test removed, so every call writes: the inert pins on statedness (a round trip that states the
   values in force stops a later `$setTreePrior` reaching the forest).
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
11. On one forest, the model's tree prior slot not mirrored: copy after `$setForest(1, power = )`.
12. A creation check left out of the shared function (student() among sharers, the lone forest's two sizes): the
    refusal table.
13. Slice 3: the other sharers' leaf scales left at 1/K after one leaves, or the leaver kept in the k pool: the
    unequal-count twin.

## Interaction with the other pending plans, and order

1. coefficient-prior-opt-in (acee2dbf) lands first. Its `student(df)` and the record's `amplitude.df` are what
   `amplitude` reads and gives back; a `student(3)` forest's round trip and `sd` restated keep its df (its pin
   "a restated sd keeps a forest's student() prior" moves to `$setForest`). Shared:
   [`validateForestAmplitude`](../../R/model.R), [`reportLeafPrior`](../../R/dbarts.R).
2. tree-prior-flat (d5ad5213, a fix round amending its live-door text now) lands next. Its plain `power` and
   `base` are this door's names; its collision rule's message names `$setForest`. Collision: its `$setTreePrior(power
   = NULL, base = NULL, split.probs = NULL)` reads NULL as no change, where open call 3's recommendation reads it
   as a return to the level above; under that recommendation slice 4's door uses `missing()` and NULL returns
   the fitting function's value to its default, as dec-B514 reads a bare `normal()`. Records step for its
   writer, before slice 4 is built.
3. forest-after-merge (f251f037, built) lands whenever it is green. Its derivation of a numeric forest's scale
   over the rows in is what `sd = NULL` uses on a numeric forest that stated `sd` at creation (Calls made 9).
4. control-migration (93f95ac0): slice 2 is built on its slice 1a's tip and this branch lands after it merges.
   Its row "`$getForest` reads the stated count from the record, so giving it back changes nothing" holds: the
   count is in the reader's `in.force` part and `n.trees` is never a statement name, refused by name. Slice 1 is
   built before slice 1a, off main.
5. forests-shared-total (f93810e2): its D10 refusal of a sharer's restated `sd`, and its R item 6
   (`$setLeafPrior(forests = )` refusing it), are superseded by dec-B518's rule: creation takes a sharer stating
   `sd`, so the door does (Calls made 12). Records step for its writer: D10's bullet and item 6 point here. Slice 3
   is built by whichever lands second. Both edit [`resolveForestSpreads`](../../R/dbarts.R), which this item
   deletes: if forests-shared-total lands first, its edit there moves to `$setForest`.
6. setmodel-parameters (abdcea31): slice 4's open call 1 is answered by dec-B518; its item 1 per-forest install is
   this item's slice 1, which it takes (this item lands first; if slice 4 is ready first it builds the install as
   planned and slice 1 drops); its item 4 door is the fitting function's level only; its item 5 refusal "naming
   `$setLeafPrior(forests = )`" names `$setForest`; its open call 2 is open call 4 here, built by slice 4 (this
   item records the variance forest's statedness). Records step for its writer.
7. coefficient-prior-live (dec-B511) follows this item and coefficient-prior-opt-in, adding a changed `amplitude`
   to `$setForest` (Calls made 10).
8. leaf-prior-writer-by-level (dec-B515's remaining half) and forest-list-by-label (predict's half of dec-B516):
   independent; textual meetings in R/dbarts.R only.

Order of landing: coefficient-prior-opt-in; tree-prior-flat; forest-after-merge when green; control-migration;
this item (slice 1 built early, slice 2 on slice 1a, merged once after control-migration); forests-shared-total
and slice 3 by whichever is second; setmodel-parameters' slice 4; coefficient-prior-live. Engine features merge
serially.

## Cost

- Four methods and a reader attribute; `$setLeafPrior` loses an argument main has.
- An explicit NULL and a name left out mean different things, under open call 3's recommendation.
- A refusal list that must track creation's; the shared check function and the refusal table carry it.
- Slice 3 adds engine work forests-shared-total had avoided by refusing.
- Until coefficient-prior-live lands, a changed `amplitude` is refused by name on main (never in a release).

## Calls made

1. The reader returns the statement, NULL where unstated, and the values in force as an attribute, plain list, no
   class. Alternatives: values in force alone, which state everything when given back; a `dbartsForest` object,
   which `forests = ` would take with its basis dropped.
2. Structure names are refused whatever their value; the reader carries none in the statement, so the round trip
   needs no equality rule. Alternative: take a structure name equal to the value in force.
3. All checks before any write; a failure after one rolls back. Alternative: write per name as checked.
4. A call that changes nothing writes nothing. Alternative: write and rely on bitwise-inert installs, which still
   rewrite statedness.
5. One check function serves creation and the door, pinned by a two-door refusal table. Alternative: the door's
   own checks, which drift.
6. `treeShift` is taken, the third move setting of dec-B477, `forest()`'s name; the ruling's list read as its
   parameter names. Not on the variance forest (dec-B500). Alternative: leave it to `$setControl` alone.
7. `updateState = NULL` last, as `setLeafPrior`; the setters return the prior statement invisibly, as `options()`
   and `par()`. Alternative: NULL invisibly, as the other setters.
8. Slice 1's install is built here and taken by slice 4. Alternative: wait for slice 4, which waits on
   forests-shared-total.
9. `sd = NULL` on a numeric forest returns the default from the column's recorded scale, else from the column in
   force over the rows in at the call, recorded, as creation would derive it (dec-B514's principle, dec-B303's
   rows). Alternative: refuse where no scale is recorded.
10. coefficient-prior-live stays its own item, adding to this door; until then a changed `amplitude` is refused and
    one equal to the statement in force is taken. Alternative: fold its ~400 lines in as slice 4, delaying this
    item behind it.
11. The record stores statedness for base and power beside the resolved values the bridge reads, so re-creation is
    unchanged and slice 4 resolves the unstated ones on a change. Alternative: NULL in the record, resolved at the
    bridge, which moves a creation read.
12. A sharer stating `sd` live is taken (slice 3), as creation takes it; the share's membership is a parameter, the
    trees staying valid (dec-B254). Read as settled by dec-B518's rule; forests-shared-total's refusal was an agent
    call. Alternative: refuse, which contradicts the rule.
13. Removal of `forests = ` is R's unused-argument error, as dec-B508's removal was. Alternative: a message naming
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
2. Labels for a multinomial sampler's category forests. Background: creation states them together (one `forests`
   entry for all K) and the engine shares one specification, so no per-category change can be made; the readers
   index them by position only.
   - (a) Recommended: label them by the response's levels in every forest-indexed reader (`getLeafPrior("q")`,
     `getTrees(forest = "q")`), `forest<i>` still position i; the door takes them together with `forest` left
     out, `s$setForest(power = 3)`, a category label there refused naming that form. Cost: ~50 lines; a
     per-category change stays impossible, as at creation.
   - (b) Labels by level, and any category's label reaches all. Cost: a call naming one category changes K.
   - (c) No labels, positions as today, the door as (a). Cost: dec-B276's "by label everywhere" unmet there.
   - (d) Per-category statements at creation and live. Cost: a creation feature and engine specifications per
     category; after 1.0-0.
3. The spelling that returns a forest to the level above (it covers dec-B511's open spelling, a coefficient back
   to the normal or a hold released). Background: a name left out changes nothing, so another spelling is needed.
   - (a) Recommended: an explicit NULL, `s$setForest("tau", sd = NULL)`, as `options(x = NULL)` and
     `modifyList(x, list(a = NULL))` clear an entry and as NULL means a name not stated in `forest()`; the reader
     writes NULL for an unstated name, so a saved statement restores. Cost: NULL and a name left out differ, which
     few R functions do, and slice 4's `$setTreePrior` must read NULL the same way (a return to the default).
   - (b) A constructor, `sd = inherit()`. Cost: a new exported name, the reader returning it.
   - (c) An argument naming what to clear, `unstate = "sd"`. Cost: the reader cannot carry it, so a saved
     statement cannot be restored.
4. Held decision 6 re-read under dec-B518: does a variance forest that states no depth follow a live change of the
   fitting function's? Background: at creation it takes the fitting function's (dec-B512); `$setModel` today
   leaves it at creation's value; dec-B518 has the sampler's `$setTreePrior` reach "the forests that state none of
   their own", and `$setVarianceForest` now lets a caller state its own.
   - (a) Recommended: it follows, as at creation and as an unstated move setting follows `$setControl`. Cost: one
     call main honors changes, a heteroscedastic `$setModel` with a new tree prior, posterior-changing for that
     call, with a new arm in heteroscedastic-exact.R; built in setmodel-parameters' slice 4.
   - (b) It keeps creation's value and changes only through `$setVarianceForest`. Cost: unlike every other
     unstated setting, it follows at creation and not after.

## Not verified

- control-migration's slice 1a is unbuilt: that its record carries a one-forest model's record and the stated
  count, and that the bridge then reads every forest's base and power from the record rather than forest 1's
  from the model's tree prior, are from its plan.
- Whether the engine's tree prior caches anything derived from base and power was not read.
- The twin oracle's exemptions beyond the drawn-to-held handoffs were not enumerated by running every kind;
  creation may consume random numbers by setting in a way an install cannot reproduce.
- forests-shared-total's engine for slice 3 is read from its plan and slice 1, not prototyped; its size is an
  estimate.
- The design pass's timing (about 100 microseconds a call by label, a sweep about 350) is for `setForestBasis`;
  `$setForest`'s checks were not timed.
