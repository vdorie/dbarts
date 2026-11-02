# forest-door: a live sampler changes one forest through $setForest, by label, in forest()'s names

Status: PLANNED 2026-10-30 on c55a1477 (dec-B518); amended after its blind critique (scratch/fd-crit/report.md,
probes c1.R to c10.R there) and its re-check (scratch/fd-crit/recheck.md); rewritten the same day for the
maintainer's rulings dec-B519 to dec-B525. Slice 1 is built and reviewed LAND (branch forest-door-s1, da875712).
Slice 2 is built on branch forest-door over 341a1c1ed, its step 0 amending this plan; a fresh review and its mutants follow.
No call is held; the rulings' own open questions are ruled too, in [Open calls](#open-calls). Rulings:
dec-B518 (the door), dec-B519 (a written NULL), dec-B520 (the names), dec-B521 (multinomial), dec-B522 (no return
spelling), dec-B523 (later write wins), dec-B524 (the variance forest), dec-B525 (a sharer stating its size),
dec-B535 (the multinomial reader), dec-B538 (the value in force handed back), dec-B539 (a drawn normal coefficient),
dec-B540 (`$setLeafPrior` where the coefficient is drawn), dec-B550 (the trees carry a forest's `sd`, held or drawn), with
dec-B254 (parameters, never structure), dec-B276 (a forest by label or position), dec-B477 and dec-A239 (per-forest
move settings), dec-B511 (the coefficient's prior live), dec-B512 (plain power and base), dec-B514 and dec-B515 (a
setter means what creation means), dec-B516, dec-B129 (leaf where it means a leaf), dec-B247 (a live count does not
change), dec-B459 (copy of a live chain), dec-B500 (no `treeShift` on the variance forest), dec-B505 and dec-B509
(forests that share a total), dec-B513.

agent: slice 1 done. One opus implementer (effort high) for slice 2, one opus review that runs the mutants below,
the same reviewer checking the fix round unless the fixes change engine behavior. Slice 3: one opus implementer and
its own opus review. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)): slice 2 is a door that
changes priors between sweeps, slice 3 edits the shared-total group in the engine.
rng: NEUTRAL, bit for bit, for every call main accepts and keeps. Two main calls, neither shipped, stop being taken:
`$setLeafPrior(forests = )`, removed and respelled `$setForest(forest, sd = )`, which draws bitwise as the old
spelling did (an A/B against a c55a1477 build); and a name on a multinomial model's one creation entry, refused.
Every new call is held to the twin oracle: before any run its next draws are those of a sampler created with the
same statement, the exemptions named in [RNG class](#rng-class). Slice 3's state has no creation call and takes the
oracles named there. No call draws a random number.
window: for 1.0-0, after the merge to main. Slice 3 is required for 1.0-0 (dec-B525) and lands after
setmodel-parameters' slice 4; until then its call is refused on main, in no release.
budget: ~3230 lines in the tree. Slice 1, built: 527 (engine 44, bridge 37, tests/cpp 314, tinytest 132). Slice 2
~1930 (R ~470 added and ~130 removed, help ~150, design notes ~20, NEWS ~10, tinytest ~1150, test-forest-door.R ~880
of it over the kinds-by-names grid). Slice 3 ~720 (engine ~95, bridge ~45, R ~70, help and design note ~30, tests/cpp
~180, tinytest ~200, exact arms ~100). Records ~50. Agent time, estimated from comparable landings and not measured:
slice 2 ~16, its review with mutants and the fix round ~6, slice 3 ~10, its review and fix round ~4, gates ~4,
rebases ~2: ~42. Plans here run 1.5 to 2 times low (slice 1 ran 2.1 times its lines, nearly all tests), so 60 to 85.

## Summary for the maintainer

You ruled that a live sampler changes one forest at a time, naming it by its label, in the names `forest()` takes:
`s$setForest("tau", sd = 0.5, power = 3)`. Its reader, `s$getForest("tau")`, gives every value that forest runs
under, as a plain list, and handing the list back changes nothing. The variance forest has its own pair.

What a call does. It changes only what it names, from the next sweep, and leaves the trees, leaf values and
coefficients already drawn as they are. Each name means what it means when the sampler is created. A forest's
columns, tree count, basis, interaction limits and blocks are refused by name: they change what the model is, and
the route is a new sampler. Anything creation would refuse, the door refuses too, with the same reason, except that
two doors reaching one setting are not a refusal: the later write wins.

There is no way to send a forest back to following the fitting function. A setting a forest states is changed by
stating another value. A NULL, like a name left out, changes nothing. The variance forest takes the mean forest's depth
when created and is on its own after that.

A forest sharing the prior total with other forests that states its own size leaves the group; the others stay where
they are, and a later change of the total scales them by the shares they were created with. It has no way back in.

What goes. `$setLeafPrior(forests = )`, built on main and never released, is removed. Nothing changes for
bartCause, stan4bart, treatSens or bairrtt.

Changing a coefficient between held and drawn, or between the normal and the t, needs the engine work you ruled for
in dec-B511; that item stays its own and adds it to this door afterwards. Until then the door takes the prior in
force and refuses a different one.

One note you may veto: the door also takes `treeShift`, the third per-forest move setting `forest()` takes, though
your ruling's list leaves it out.

## Goal

A live sampler has `$setForest(forest, ...)`, `$getForest(forest)`, `$setVarianceForest(...)` and
`$getVarianceForest()`. `$setForest` changes on the forest named by label or position only the names it is given,
among `sd`, `power`, `base`, `amplitude`, `proposal.probs`, `n.perturb.cuts` and `treeShift`, each meaning what it
means in `forest()` on that forest at creation; refuses structure by name; refuses a statement where creation
refuses the same; and changes nothing where the value given is the one in force. `$getForest` returns every value
in force in those names, so `do.call(s$setForest, c(list(f), s$getForest(f)))` changes nothing, statedness and draws
included. Every change is recorded where re-creation reads it, so a copy and a reload run what the live sampler
runs. `$setLeafPrior(forests = )` is gone.

## Constraints

- Out of scope: a changed coefficient hold or prior (coefficient-prior-live adds it to this door);
  `$setTreePrior`, and `$setModel` and `$setLeafPrior` on samplers of several forests (setmodel-parameters' slice
  4); per-category multinomial statements, not offered (dec-B521); per-forest handles.
- The fitting function's level is not recorded apart from the forests' (dec-B522). No field holds it. A forest's
  value in force and whether it states it are all a level ever needs.
- No change to inst/include/dbarts/dbarts.h, the state format or any consumer.
- The engine sees numbers, never a record or an R name: the bridge converts.

## What exists today

Run on a c55a1477 build (scratch/fd-plan/p1.R to p4.R; the critique's c1.R to c10.R), read again against
control-migration's slice 1a (60b960ed) and forests-shared-total (c8ad05e2):

- The per-forest writers are [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd)'s `forests` list,
  which restates `sd` alone, refuses itself on one forest, and is used at 45 call sites in 9 tinytest files; and
  `setForestBasis` and `setForestWeights`, by label or position through
  [`samplerForestIndex`](../../R/bartcore.R). The forest-indexed readers are `getLeafPrior`, `getK`,
  `getForestFits`, `getForestAmplitudes`, `getForestVariableCounts` and `getTrees(forest = )`. `predictForests` takes
  no forest and refuses a multinomial model.
- After slice 1a the model holds one `dbartsForestRecord` per forest of a model of several (`n.trees`, `stated`,
  `base`, `power`, `held.amplitude`, `sd`, `amplitude.df`, `vars`, `interactions`, `blocks`, `label`, `basis.term`,
  `moves`) and a `dbartsVarianceRecord`. `stated` is a character slot, today only the count's. A record's `sd` is NA
  where none is stated, its `moves` NULL where none is. On a model of one forest the forest's depth stays in
  `model@tree.prior`, and its size in `model@leaf.prior`. The bridge reads forest 1's base and power from the model's
  tree prior on every sampler ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)), the others' from the
  records.
- The engine keeps a tree prior on every forest and the variance forest, caching nothing derived from base or power
  ([`CGMTreePrior`](../../src/bartcore/model.hpp)). Live sd and move installs exist:
  [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp) and
  [`bartcore_setMoveSettings`](../../src/R_interface_bartcore.cpp). Slice 1 adds the tree prior's.
- Creation: a lone forest's `sd` beside any stated `leaf.prior` is refused, "states its size twice"; `student()` on
  a lone forest and on a sharer is refused; a variance forest beside several forests is refused; an unstated
  variance forest takes the mean forest's resolved depth, a lone forest's own included (c9).
- Main takes `$setLeafPrior(normal(k = 3))` on a `forest(sd = 2)` sampler, `$setLeafPrior(normal())` on a probit
  `forest(sd = 2)` sampler, k drawn again, and `$setModel` with a changed power on a `forest(x, power = 3)` sampler.
  dec-B523 keeps all three.
- A multinomial sampler (matrix interface only) takes one `forests` entry for every category forest; a name on that
  entry is taken and ignored (c6); the category forests have no labels. The levels are the column names of the count
  matrix ([`dataCounts`](../../R/bartcore.R)), there being no `data@response.levels` and no levels block in print.
  Slice 1's engine install takes a per-category value and moves draws, so a
  per-category change could be made; it is simply not offered (dec-B521).
- `forest(amplitude = normal())` is taken (forest-amplitude-normal, built over forests-shared-total); `fixed()` and
  `student(df)` too. On a forest with a basis it records and draws what leaving it out does.
- `getLeafPrior(2)` on a heteroscedastic sampler is out of range; no reader shows a forest's power or base.

## Design

### What a forest's record says

Each forest's values in force and whether it states them, nothing else:

- `sd`: the record's `sd`, NA where none is stated (as slice 1a has it).
- `power`, `base`: the record's values; `stated` gains the names `"power"` and `"base"`, written by
  `resolveSamplerSpec` at creation, a real value named counting as stated, the default's own included, a NULL not
  (dec-B519). Forest 1's depth is also in the model's tree prior, which the bridge reads; every write to it writes
  both, on every sampler.
- The move settings: `moves`, as slice 1a has it.
- The coefficient: `held.amplitude` and `amplitude.df`, as slice 1a has it.
- A model of one forest, and a multinomial model's one entry, need no statedness: every door reaches that one
  setting ([Forests one setting reaches from two doors](#forests-one-setting-reaches-from-two-doors)). Their values
  stay in the model's tree and leaf priors, and in the one record's `moves`.
- The variance forest needs none either: nothing follows anything (dec-B524).
- Slice 3 adds the record's `share` ([Slice 3: a sharer leaves the group](#slice-3-a-sharer-leaves-the-group)).

These slots are carried by this plan's slice 2, not slice 1a: slice 1a is a move, held neutral, and the door is the
first reader. control-migration's slice 1b must not compare `stated`'s new names or `share` in its structure check,
these being parameters' statedness.

### The arguments, by kind of forest

Every value is checked by the validator creation uses for it. A forest with no basis stating no `sd` is a sharer
wherever forests-shared-total names it (slice F included: a lone forest with no basis beside forests with a basis
is one).

| name | lone forest | sharer | no basis, stating `sd` | basis: factor | basis: numeric | multinomial, every category | variance forest |
|---|---|---|---|---|---|---|---|
| `sd` | its size, written as `$setLeafPrior(normal(sd = v))` writes it | a group of two or more: leaves the group, slice 3, refused until then; a lone forest with no basis beside forests with a basis: taken, leaving the leaf prior with no forest as creation's `sd` does | the trees' scale: its size, held, or under `student(df)` or `normal()` its size at a coefficient of 1, the coefficient's prior median (dec-B550) | the level difference's size, the same way | the size per unit of the column, the same way | refused, as at creation | not a name |
| `power`, `base` | its depth | its depth | its depth | its depth | as factor | its depth | its depth |
| `amplitude` | `fixed()` | `fixed()`; `student()` refused, as at creation | `fixed()`, `student(df)`, `normal()` | `fixed()`, `normal()`, `student(df)` | as factor | not a name | not a name |
| `proposal.probs`, `n.perturb.cuts` | its own | its own | its own | its own | its own | its own | its own |
| `treeShift` | its own | its own | its own | its own | its own | its own | not a name (dec-B500) |

`amplitude` takes the prior in force, as no change; a different one is refused until coefficient-prior-live.
`forest(amplitude = normal())` becomes taken at creation and at the door on a forest with a basis, drawing and
recording exactly as leaving it out does (~10 lines in [`validateForestAmplitude`](../../R/model.R)); only a bare
`normal()`. On a forest with no basis it is a drawn normal coefficient, taken where `student(df)` is (dec-B539; built by
TODO forest-amplitude-normal after slice F, not here), and at this door the value in force only until
coefficient-prior-live. A coefficient has no way back to following a default: a hold is released by writing a shape
(dec-B522).

Held against drawn follows the coefficient's hold, as at creation: a forest's size is its coefficient times `sd`, the coefficient
being exactly 1 where held and a multiplier with prior median 1 where drawn, under the normal or under `student(df)` (dec-B510, in
dec-B550's words). On every kind and every option the write sets the trees' scale and nothing else: no leaf, no coefficient and no
auxiliary moves, before the first sweep too (dec-B550; forest-amplitude-normal builds that and lands before this door's slice 2).
The help says a write changes the model, in those words at most (dec-B550). One exemption from the twin oracle follows: a
coefficient's start is set once at creation, and on a drawn forest with no basis it is the unit over the `sd` then in force, so an
`sd` written there before any run is the same model and not creation's draws; every other forest's start does not read `sd`, so
its twin holds. A numeric forest's default `sd` re-derives on
`setForestBasis(updateBasisScale = TRUE)` only while unstated, and reaches the same place.

### Forests one setting reaches from two doors

The later write wins (dec-B523). Where a door of the fitting function and `$setForest` reach one setting of one
forest, each writes the one value in force. Creation still refuses a setting stated at both levels of one call; a
live sampler has no such refusal, because no level is recorded for it to contradict.

- `$setForest` is taken on any forest, whatever the fitting function stated.
- A fitting-function door (`$setLeafPrior`, `$setTreePrior`, `$setModel`) has as its candidates the forests with no
  basis. Where there is exactly one, a lone forest or a lone forest with no basis beside forests with a basis or a
  multinomial model's one entry, it is written, whether it states its own or not. Where there are several, the door
  reaches those stating none of their own, and where all of them state their own it is refused, naming `$setForest`,
  since the call would do nothing. Where there is none it is refused per name (tree-prior-flat's call 13). A forest
  stating `student(df)` or `normal()` has no k: there the door writes the forest's sd, the trees' scale holding it, by
  dec-B540 and dec-B550, in setmodel-parameters' slice 4.
- A value at a fitting-function door equal to the one in force is no change, which is how `$setModel(s$model)` and a
  bare `normal()` beside forests that all state `sd` stay taken.
- `$setLeafPrior(normal())` on a probit `forest(sd = 2)` sampler returns k to drawn, from the spread in force (the
  exemption in [RNG class](#rng-class)). `$setLeafPrior(NULL)` stays refused, as on main. `$setTreePrior(power =
  NULL)` reads as left out (dec-B519).

This plan's writers are the one-forest ones: `$setLeafPrior` and `$setModel` on one forest already behave so on main,
and slice 2 only pins it. Samplers of several forests get the rule in setmodel-parameters' slice 4, which owns those
routes.

### Parameter and structure

Parameters, taken: `sd`, `power`, `base`, the coefficient's hold and prior (dec-B254, dec-B511), and the move
settings, which are run settings and not model (dec-B477). Structure, refused by name whatever its value, before
anything is written (dec-B254, dec-B518): `vars`, `n.trees` (dec-B247), `interactions`, `blocks`: "'n.trees' is
fixed at creation: a forest's tree count is part of the model; make a new sampler". `basis` is refused naming
`$setForestBasis`, the basis being data. On the variance forest, `vars` and `n.trees` are refused alike. A
`$setVarianceForest` on a sampler with no variance forest is refused: whether there is one is structure. These names
are not formals; they are read from `...`, which ends the signature.

### The reader and its round trip

`$getForest(forest)` returns a plain named list, no class and no attribute, of every value in force in `$setForest`'s
argument order, those names the kind takes in the table above and no others: the sd in force, the depth, the
coefficient's prior as `fixed()`, `normal()` or `student(df)`, the move settings. A value the forest takes from above
reads like any other. There is one exception, told to the maintainer: where k is drawn, on a lone forest or a sharer,
there is no single sd in force, so `sd` is NULL there, which by dec-B519 changes nothing when given back.

Giving the list back changes nothing, statedness included: a value equal to the one in force is no change, on every
name, so a forest following the fitting function goes on following it, a sharer stays in its group, a numeric
forest's default goes on re-deriving. The equality is exact, the reader's numbers being the ones the door compares.

`$getVarianceForest()` is the same for `power`, `base`, `proposal.probs` and `n.perturb.cuts`. The setters return
NULL invisibly, as every other sampler setter does. Tree counts and columns are not in the list, being refused on
the way back; print and `getTrees` show them.

`getLeafPrior` keeps its fields (leaf-prior-writer-by-level's); its help stops saying its `forest(sd = v)` goes back
into `$setLeafPrior(forests = )` and says its `sd` is the one `$setForest` states.

### The variance forest

It takes the mean forest's depth at creation; after that `$setTreePrior`, `$setModel` and `$setForest(1, )` leave it
as it is and `$setVarianceForest` changes it (dec-B524). The help says so, on `$setForest`, `$setModel` and
`$setVarianceForest`, and on `$setTreePrior` with setmodel-parameters' slice 4, which builds it; ~15 lines with its pin. No heteroscedastic-exact arm: main's call does not
move. So the twin of `$setForest(1, power = 4)` on a heteroscedastic sampler created at 2 is a creation stating
forest power 4 and variance forest power 2.

### Multinomial models

The category forests take one statement (dec-B521):

- `s$setForest(power = 4)` and `s$getForest()` take them together with `forest` left out.
- `s$getForest(2)` and `s$getForest("b")` are taken like every other reader and return the one statement, a plain list
  like any other; one help sentence says it covers every category (dec-B535).
- `s$setForest(2, power = 4)` and `s$setForest("b", power = 4)` are refused, position 1 too, a position or a level alike: "a
  multinomial model's forests take one statement, so this would change every category; write
  s$setForest(power = 4)".
- A name on the one creation entry is refused at creation ("a multinomial model's one forest entry states every
  category forest and takes no name"), ~10 lines, a call main takes and no release had.
- The category forests take the response's levels as labels on every forest-indexed reader (`getLeafPrior("b")`,
  `getK("b")`, `getForestFits`, `getForestAmplitudes`, `getForestVariableCounts`, `getTrees(forest = "b")`,
  and `getForest`), `forest<i>` still position i, a level that collides with another position's
  `forest<i>` refused by `selectForest`'s existing conflict message. The NULL-reads that name their list by labels
  name it by levels. This pays dec-B276's debt; ~120 lines, ~50 of code. TODO multinomial-forest-labels closes.

### The forest argument

A label or a position (dec-B276), read by [`samplerForestIndex`](../../R/bartcore.R); a one-forest sampler's is 1 or
`"forest1"`, or left out; a multinomial sampler's as above. Elsewhere it must be given. `"variance"` is an ordinary
label for a mean forest, so the variance forest has its own methods.

### How a call is carried out

The signature is `setForest = function(forest = NULL, sd = NULL, power = NULL, base = NULL, amplitude = NULL,
proposal.probs = NULL, n.perturb.cuts = NULL, treeShift = NULL, updateState = NULL, ...)`, every setting's default
NULL, so a written NULL is the same as leaving the name out (dec-B519).

1. Resolve the forest, the multinomial refusal among it. Refuse structure by name.
2. Drop the names NULL. Drop each name whose validated value equals the one in force, which the reader's function
   gives. If none is left, return: no engine call, no record write.
3. Check the forest's statement with the names left merged in, not the names alone, as creation checks it: the
   per-name validators ([`validateForestSd`](../../R/model.R), [`validateForestAmplitude`](../../R/model.R), the base,
   power and move checks) and the cross-forest rules creation applies (held and drawn on one forest, `student()` among
   sharers, the multinomial refusals of [`refuseMultinomialLeafPrior`](../../R/dbarts.R)), factored out of
   [`resolveSamplerSpec`](../../R/spec.R) where they are inline so that one function serves both doors. Creation's
   rule against one setting at two levels stays creation's.
4. Install each changed part, every check done first, in the order that leaves only the `sd` able to refuse after a
   write: move settings by rewriting the move record and calling
   [`bartcore_setMoveSettings`](../../src/R_interface_bartcore.cpp), as `$setControl`'s move install does; `power` and
   `base` through slice 1's entry, every category on a multinomial model; last `sd`, through
   [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp), or on a model of one forest through the leaf-prior
   install. What the `sd` refuses undoes the two before it: on a model of one forest, an `sd` whose k overflows, which
   the leaf-prior install refuses in `$setLeafPrior`'s words; on a model of several nothing a sampler can reach now, the
   bridge's refusal being a backstop.
5. Record, once the engine holds the values and as creation with that statement would write it, so `identical()`
   holds between the live sampler's records and a creation twin's. Forest 1's depth in the model's tree prior on every
   sampler, the bridge still reading it there, and on a model of several every forest's in its record too, its
   `stated` gaining the name in creation's order; a model of one has no depth or statedness in its record (NA). A
   model of one's `sd` in the leaf prior, by what `$setLeafPrior` writes in that leaf model; on a model of several in
   the record, and where it was the last forest with no basis stating none, the fitting function's leaf prior is left
   with no forest, as creation leaves it.

### What goes and what is never built

- `$setLeafPrior(forests = )` goes with `resolveForestSpreads` and
  `writeForestSpreads`; `forests = ` becomes R's unused-argument error, as `forest(leaf.prior = )`
  did under dec-B508. It never shipped, so it takes no tombstone. `$setLeafPrior` keeps `leaf.prior`. Its 45 call
  sites are respelled.
- `$setTreePrior(forests = , variance = )` is never built: setmodel-parameters' slice 4 builds `$setTreePrior(power,
  base, split.probs)`.
- Per-forest handles stay possible later as sugar over the door (dec-B518).

### RNG class

- Inert: a call that drops every name at step 2 writes nothing and moves no draw: the round trip, a subset restated,
  a NULL.
- The twin oracle, before any run, on two chains: a taken call's next draws are bitwise those of a sampler created
  on the same data and seed with the statement it leaves. Exempt by name: a handoff from held to drawn, where the
  drawn spread starts at the one in force rather than at creation's start (c4: `$setLeafPrior(normal())` on a lone
  probit forest that stated its `sd`; the coefficient's in coefficient-prior-live), as setmodel-parameters' call 8
  exempts its handoffs; and `sd` on a drawn forest with no basis, whose coefficient keeps creation's start (dec-B550). From drawn to held is not exempt: `$setForest(1, sd = v)` on a drawn k draws as
  `forest(x, sd = v)` does (c4).
- A value equal to the one in force has the draws and the record of no call at all, so the twin of a round trip is the
  sampler itself, not a creation stating those values, whose draws agree and whose statedness does not.
- After a run: the change acts from the next sweep; `copy()` (dec-B459) and a `saveRDS` and `readRDS` after
  `storeState` continue `identical()` to the live sampler.
- No call draws: every install restates a scale, a depth or a mixture.

### Copy, reload and state

Every change rides the records and the model's slots; re-creation reads them. A state carries no prior (dec-B195) and
installs across statements as today. No state block, no `stateFormatVersion` move, no change to
inst/include/dbarts/dbarts.h.

## Slices

### Slice 1: a forest's tree prior installed live (engine and bridge)

Built and reviewed LAND, branch forest-door-s1, da875712, on main 0766f069: `Sampler::setForestTreePrior` behind a
facade virtual, `bartcore_setForestTreePrior(ptr, forest, base, power)`, forest -1 the variance forest, its tests/cpp
pins and test-forest-tree-prior-install.R. It merges with slice 2. Neutral: nothing in R reaches it until slice 2.

### Slice 2: the door (R), on control-migration's slice 1a

0. Before it starts, re-read slice 1a as built and forests-shared-total with slice F: the records' slots, whether the
   bridge still reads forest 1's depth from the model's tree prior, whether `forestSharers` names a lone sharer, where
   multinomial levels are; amend step 5 of the call to them.
1. `stated`'s `"power"` and `"base"` written at creation; the reader's function.
2. The four methods, their docstrings and help, by [How a call is carried out](#how-a-call-is-carried-out).
3. `$setLeafPrior(forests = )` removed; the 45 call sites respelled to `$setForest`.
4. The cross-forest checks factored out of [`resolveSamplerSpec`](../../R/spec.R) where they are inline.
5. [Multinomial models](#multinomial-models), the labels included.
6. `forest(amplitude = normal())`: built already by forest-amplitude-normal; the door only reads it.
7. The variance forest's help sentences and pin.
8. A sharer's `sd`, where `samplerSharers` names two or more, refused, a value equal to the one in force excepted:
   "this forest shares one prior total with forests ...; a forest's own size cannot be changed live yet: make a new
   sampler". On main only; slice 3 replaces it. It moves here from forests-shared-total's `resolveForestSpreads`. A
   lone forest with no basis beside forests with a basis is named by `samplerSharers` too, at a share of 1; its `sd` is
   taken, as forests-shared-total's own route took it, dec-B525's end state with no group to leave, and its twin is a
   creation stating that `sd`.

Slice 2 builds on slice 1a's tip because the door reads and writes the records slice 1a moves onto the model.

### Slice 3: a sharer leaves the group

The rule (dec-B525). `$setForest(f, sd = v)` on a sharer, v not the value in force, takes it out of the group with
size v, its k pinned as any forest stating sd's is. The others keep their leaf scales, shares unchanged: a share is
set once, at creation. A later statement of the total on the fitting function's door scales the forests still sharing
by their creation shares and does not re-divide. The group's k draw and the rescale step read the forests still in
it. The share is recorded on the forest, so a copy and a reload continue. There is no way back in but a new sampler.
`sd = v, amplitude = student(df)` in one call is refused with every changed `amplitude`, until coefficient-prior-live.

A statement of the total needs no new engine: `translateNamedSd`'s sharer branch and setmodel slice 4's group install
write k, and a sharer's spread is `options_.nodeScale * sqrt(share) / sqrt(m) / k`, so each stays at its creation
share of whatever total is stated, and the sizes stop adding up to it. Help says so.

Before it starts, re-read slice F as built: whether a lone sharer held at 2 on dec-B290's entry reads as a sharer in the
engine, which of `$setForest`'s two routes it takes and what its record's `share` then says; and setmodel slice 4's
group install.

Engine, on forests-shared-total as built (c8ad05e2):

1. `ForestSpec` gains `share`, NaN where none is handed. `ForestMapEntry` has one, read by every path below, but
   `ForestSpec` has none, and `forestSharers` derives the group from what the forests state and clears a group of
   fewer than two, so it cannot reach this state from a record. `forestSharers`: a finite `share` on any forest names
   the group, whatever its size, one included; otherwise it derives as built. `Chain`'s K-forest constructor takes
   each sharer's share from the spec where handed, otherwise 1 / K as built, and hands it to
   `forestMapEntry` per forest. Slice F's held-2 exception for a lone sharer must key on a share of 1, not on one
   sharer: one left of two at a held k of 2 keeps its share, and dec-B290's entry would give it the whole size.
2. `Chain::leaveSharedTotal(f, sd)`: false, writing nothing, off a sharer, on an sd not positive and finite, or where
   the unit is not (as `setForestMapSd`). Erases f from `sharers_`, with no re-clearing below two; states the sd and
   re-derives the forest's `forestMapEntry` with no share, a held entry whose coefficient prior is a sharer's, so no
   combiner write; `leaf.scale = mapLeafScale(f)`; `updateK = false`, `k = 1`, `namedSd` NaN. Touches no other forest
   and no leaf value. `Sampler` method behind a facade virtual, checked on every chain before writing any, as slice
   1's is.
3. State. `fs.k` is absent for a forest not drawing k (the capture in `Chain`), and both installs write k per forest
   where `updateK`. A state captured after a leave and installed on a sampler where that forest still shares would leave
   the group's k unequal, which `drawSharedK` and the rescale step assume it is not. Both installs write the first k
   the state holds among the sharers to every sharer. Bitwise on every state where the sharers' k agree, which is
   every state today.

Every other `sharers_` path needs no edit and takes a pin:

| path | one left | none left |
|---|---|---|
| sweep's k draw, both sites (`!sharers_.empty() && isSharer(f)`, `f == sharers_.front()`) | draws at it; a leaver that was front moves the draw to the next | else-branch, `updateK` false, no draw |
| `drawSharedK` | one forest's leaves standardized, the group's rule | never called |
| `forestRescaleApplies` | the sharer branch, one forest scaled, every other read as offset, the leaver included | the empty branch, whose `forests_.size() != 1 \|\| combiner_` declines with no draw, as a creation twin does |
| `drawForestRescale` | as above | never reached |
| `draw.k`, `results.k` | the front's | forest 0's k of 1, and the channel is off (bridge 2) |
| `getK`, `forestCalibration`, `mapLeafScale`, `translateNamedSd` | per forest by `entry.share` | as any forest stating sd |

Bridge:

1. `bartcore_leaveSharedTotal(ptr, forest, sd)`, every check first; the unit's message as `bartcore_setForestSd`.
   `applyForestStatements` reads the record's `share` and hands it on; a forest with no basis stating no sd, no df and
   no share beside a forest with a share is refused.
2. `run()`'s `kIsDrawn` reads `forestCalibration(0, 0)` alone, so forest 1 leaving would drop the group's k column.
   It reads any forest's `kHasHyperprior`.
3. `printInitialSummary`'s "each 1/K" counts the sharers; it reads the share.

R:

1. The record's `share`, numeric, NA where none, positive and at most 1, and NA where `sd` is stated; written at
   creation by `resolveSamplerSpec` on every sharer. `$setModel` with a differing share is refused: a share is set
   once.
2. `$setForest` routes `sd` on a sharer to the entry, writing the record's `sd` and clearing its `share`. Slice 2's
   refusal goes.
3. The last leave. `refuseUnsupportedAmplitudeComposition` refuses a k hyperprior, a k other than 2 and a named sd
   where no forest shares, so a reload with the group's leaf prior still on the model would fail. When the last sharer
   leaves the model's leaf prior and hyperprior are written as creation with every forest stating sd writes them (the
   fixed default), so the records are `identical()` to that twin's and the reload runs.
4. `forestPriorLines` and `sharerLeafPrior` read the share rather than counting. `samplerSharers` keeps reading the
   engine, so one left is named.

Help: `$setForest` says a sharer's `sd` leaves the group, that the others keep their sizes, and that there is no way
back; `$setLeafPrior` that the total reaches the forests still sharing, each at its creation share. docs/design's
forests-shared-total note: its live-sampler sentence and its "each reading what the forests state; no flag passes
between them" both change. NEWS: within the `$setForest` item.

The oracle. No creation call reaches a group with one forest gone: creation stating that forest's `sd` divides the
total among fewer. So:

- The recorded twin, bitwise: a leave before any run against a sampler created from the live sampler's records
  (`new("dbartsSampler", control, model, data)`, the path a reload takes), two chains. It ties the install to the
  carried share, and says nothing of whether the carried share is right.
- The stated-sizes twin, for that, at a held k: creation with every remaining sharer stating `forest(sd = )` at its
  part in force reaches the same sizes by another path. Leaf scales within 1e-12 and draws within 1e-8 over a short
  seeded run, not bitwise, the map's operands being others. Independent of the carried share; fails a re-division.
- The creation twin proper where one exists: none left against every forest stating sd at the counts the split gave;
  a lone sharer leaving against `forest(x, sd = v)` beside its forest with a basis.
- For a drawn k, which no stated size reaches, exact arms: forest-prior-exact.R, three forests, one given `sd` by the
  door before the run, the two left at T / sqrt(3) and the oracle `sharedOracle`'s; probit-k-scale-exact.R, one sharer
  of two left under a drawn k beside the leaver's stump read as offset, integrated in the outer quadrature, and two of
  three left at unequal counts. Each without the door's share fails.
- After a run: `copy()`, a reload and a state saved before and after the leave continue `identical()` to the live
  sampler.

Built by a new branch off main after setmodel-parameters' slice 4 lands, so its "later statement of the total" pins run
through the real door. Merged on its own: forest-door's feature merges twice, slices 1 and 2 first, because slice 4
reads them.

## Help and NEWS

- man/dbartsSampler-class.Rd: usage and value for the four methods, the argument table above in prose, the refusals,
  the round trip with `do.call`, that no call touches the trees, the variance forest's sentence; `setLeafPrior` loses
  `forests`; `getLeafPrior`'s sentence on `forest(sd = v)` changes; the "variance forest ... is not addressable"
  sentence points to `setVarianceForest`; the readers say the levels are labels on a multinomial model.
- man/forest.Rd and man/dbarts.Rd: their `$setLeafPrior(forests = )` sentences become `$setForest`; `amplitude` lists
  `normal()`.
- inst/NEWS.Rd: the 1.0-0 multi-forest item's "which `$setLeafPrior(forests = )` restates" becomes `$setForest`, and
  one item adds the four methods. Nothing about `forests = ` itself, which never reached a release.
- docs/design: forest-defaults-by-kind.md, forest-multiplier.md, multiplier-combiner.md and forest-sd-unit.md name
  `$setLeafPrior(forests = )`; each is respelled or marked as the history it is.
- tools/check-rc-codoc.R covers the four methods with no change to the script. Their usage carries the NULL defaults
  of the signature.

## Tests and gates

tinytest, a new test-forest-door.R (~880 lines), every R-side twin on two chains:

- Each argument on each kind of the table: the twin oracle before a run (draws and records `identical()` to a
  creation twin); after a run, the change in force (the reader, `getK`, the next draws against a twin given the same
  state), `copy()` and a reload `identical()` to the live sampler; forest 1 of two as well as a lone forest.
- Inert: the round trip on every kind, a subset restated, a NULL on any name: draws bitwise equal to no call, and,
  statedness holding, a following `$setControl` or fitting-function door still reaching a forest that took its value
  from above; a numeric forest's default still re-deriving after its round trip; a drawn k's `sd` NULL in the list and
  the round trip inert.
- The reader: a plain list, no attribute and no class, every name the kind takes, the depth a forest takes from the
  fitting function read there, on a sampler built by `new()` from a model carrying power 3 and on its reload.
- Later write wins, both orders on a lone forest: `$setForest(1, sd = )` after `leaf.prior = normal(k = 3)`,
  `$setLeafPrior` after `forest(sd = 2)`, `$setLeafPrior(normal())` on a probit `forest(sd = 2)` returning k to drawn,
  the same for `power` with `$setModel`; the records agreeing with the engine after each. `$setLeafPrior(NULL)`
  refused.
- Refusals: a table of statements run at creation and at the door, the same refused or taken at both with the same
  reason, the two-level one excepted; structure by name; a refused call with one good and one bad name leaving the
  reader, records and draws as they were.
- The variance forest: its own door; `$setForest(1, power = )` and `$setModel` with a changed depth leaving it, its
  draws those of a twin stating its old depth; a sampler with none refusing the door.
- Labels and multinomial: a forest selected by label where the formula and a list put it at different positions; the
  multinomial door with `forest` left out; `getForest(2)` and `getForest("b")` the one list; both setter spellings
  refused with the reason; a name on the creation entry refused; a level read by every forest-indexed reader, the same
  values as its position; `forest2` still position 2; a level spelled `forest2` at position 1 refused by the conflict
  message.
- `amplitude = normal()`: at creation on a factor and a numeric forest, draws and records identical to left out; at
  the door the value in force; on a held forest refused until coefficient-prior-live.
- The old spelling: on a c55a1477 build `$setLeafPrior(forests = )` and on the branch `$setForest(label, sd = )`,
  bitwise equal draws over the bcf, numeric and three-forest fixtures (scratch A/B, ~150 lines).

Slice 2's gates: the above, the R-loaded sanitizer run over test-forest-door.R, the equivalence trio bitwise (55, 19,
11 identical against benchmarks/baselines/MANIFEST's baselines), the four snapshot files unchanged on the reference
build, the full suite once in one process, lint against the slice's library, air, check-rc-codoc, check-win-drift,
check-doc-freshness, `R CMD check --as-cran` from a clean tarball, and the exact gates of
.github/workflows/exact-gates.yaml with `quick`, the records changing what a fit carries. No new exact gate: the twin
oracle ties each new statement's posterior to a creation the exact gates cover.

Slice 3's gates: tests/cpp plain and sanitized, `--preclean`, the mutation battery's `verify-anchors`; tests/cpp pins
on two chains at unequal counts (38, 37, 37), one leaving, two, all; the carried-share constructor; the k draw and the
rescale step over one and over none; the state both ways; the `sharers_` table above, row by row. Tinytest: the oracles
above; forest 1 leaving with the k column kept and `fit$k` stored; `getK` and `getLeafPrior` on the leaver and those
left; a statement of the total after a leave, by `$setLeafPrior` and `$setModel`, leaving each at its creation share;
all left, then a differing total refused naming `$setForest` and the default taken; `$setModel` with a model saved
before the leave refused; print. The trio bitwise, shares being written at creation; forest-prior-exact.R and
probit-k-scale-exact.R with their new arms; the full suite; `R CMD check`.

Mutants, each to fail a named pin. Slice 2:

1. The no-change test removed, so every call writes: the inert pins on statedness.
2. A value equal to the one in force marked stated: a round trip stops the fitting function's door reaching the forest.
3. A NULL read as a clearing: `s$setForest(1, sd = NULL)` moves the forest.
4. The record written and the engine not: the twin oracle. The engine written and the record not: copy and reload.
5. Base and power swapped at the door, or forest 1 written whatever the index: forest 2's twin.
6. The move record rewritten without the move install: the twin.
7. A check moved after the first write: the atomic refusal pin.
8. Structure equal to the value in force taken: the refusal table.
9. Labels read as positions, a level mapped off by one: the label pins.
10. `$setVarianceForest` writing mean forest 1: the variance twin.
11. `$setForest(1, power = )` moving the variance forest: the variance pin.
12. The model's tree prior not written for forest 1, on a lone forest or on forest 1 of two: copy after
    `$setForest(1, power = )`.
13. A creation check left out of the shared function: the refusal table.
14. Checks applied to the names given rather than the merged statement: `amplitude = student(1)` on a forest whose
    stated `sd` it needs.
15. A multinomial position writing one category: the refused pin. `$getForest(2)` reading one category's value after a
    hand write: the one-statement pin.
16. A drawn k's `sd` read as a number or NA: the round trip's NULL pin.
17. An R-side twin run on one chain: the pins assert two chains, and a per-chain write to chain 0 alone fails them.

Slice 3:

18. The others' leaf scales re-divided after a leave: the stated-sizes twin. A constructor deriving the share from who
    is left: the recorded twin, copy and reload.
19. The leaver kept in the k pool, or scaled by the rescale step: unequal-count pins and the exact arm.
20. `forestSharers` clearing a group of one: the one-left row. One left at a held 2 on dec-B290's entry: the
    stated-sizes twin at k 2.
21. `kIsDrawn` left reading forest 0: forest 1 leaving keeps the k column.
22. The installs left writing k per forest: the state-after-the-leave pin.
23. The last leave leaving the model's leaf prior: the reload of all left.
24. A statement of the total re-dividing: the total-after-a-leave pin.

## Interaction with the other pending plans, and order

Order of landing: forests-shared-total with slice F; control-migration through slice 1a; this item's slices 1 and 2;
setmodel-parameters' slice 4; this item's slice 3; coefficient-prior-live. Engine features merge serially.

1. Landed: coefficient-prior-opt-in (dec-A245), tree-prior-flat (dec-A246), forest-after-merge (dec-A244). A
   `student(3)` forest's round trip keeps its df. tree-prior-flat's live-door row "`$setForest(1, power = 2)` changes
   that forest only" now holds on a heteroscedastic sampler too, by dec-B524. Its `$setTreePrior(power = NULL, ...)` as
   no change holds. Its call 13, a statement reaching no forest refused per name, is revised for a forest that states
   its own and is the only candidate, written by dec-B523; the writer of records notes it.
2. control-migration: slice 2 is built on slice 1a's tip, 60b960ed; `stated` takes the new names (its report agrees).
   Its row "`$getForest` reads the stated count from the record, so giving it back changes nothing", its call 12's "for
   `$getForest`" and item 7's "so `$getForest` reads the record" are stale: `n.trees` is in neither the list nor the
   statement; the count's flag still serves its open call 1 and dec-B526. Slice 1b's structure check must not compare
   `stated`'s new names or `share`. dec-B526, dec-B527: nothing here reads them; `$setForest(n.trees = )` is refused
   whatever they say. Records step for its writer.
3. forests-shared-total: its D10 "a restated sd" and R item 6 are replaced by dec-B525; slice 2 moves its refusal to
   `$setForest` and slice 3 removes it. Both edit `resolveForestSpreads`, which this item deletes.
   dec-B528 touches nothing here. Records step for its writer.
4. setmodel-parameters' slice 4: takes slice 1's install and reads `stated`; owns the fitting-function doors on several
   forests under dec-B522 and dec-B523; closes its open call 2 by dec-B524. Its R reads the group from the engine, not
   `sharingForests`, so slice 3 needs no change there.
5. surface-sweep (dec-B519): creation's reading of a written NULL lands there; `stated`'s new names read a NULL as not
   written whichever lands first.
6. coefficient-prior-live (dec-B511) follows, adding a changed `amplitude` to `$setForest`, and with it
   `sd = v, amplitude = student(df)` on a sharer.
7. leaf-prior-writer-by-level and forest-list-by-label: independent; textual meetings in R/dbarts.R only.

## Cost

- Four methods. `$setLeafPrior` loses an argument main has.
- A forest that states a value can never follow the fitting function again; only a new sampler does that (dec-B522).
- A value equal to the one in force cannot pin a forest: a forest following the fitting function goes on following it
  until it states something different.
- Two doors can overwrite each other in silence (dec-B523).
- A refusal list that must track creation's; the shared check function and the refusal table carry it.
- After a leave a sampler's forests run sizes that no longer add up to the total it reads, and it is a state no fitting
  function can write (dec-B525).
- Forest-door's feature merges twice, and main refuses a sharer's `sd` between the two.
- Until coefficient-prior-live lands, a changed `amplitude` is refused by name.

## Calls made

1. The reader's list holds the names the door takes for that kind, tree count and columns left to print and `getTrees`.
   Alternative: every value in force in `forest()`'s names, which the door would refuse handed back.
2. Structure names are refused whatever their value. Alternative: take a structure name equal to the value in force.
3. All checks before any write; a failure after one rolls back. Alternative: write per name as checked.
4. A call that changes nothing writes nothing. Alternative: write and rely on bitwise-inert installs, which rewrite
   statedness.
5. One check function serves creation and the door, pinned by a two-door refusal table. Alternative: the door's own
   checks, which drift.
6. `treeShift` is taken, the third move setting of dec-B477 and `forest()`'s name; not on the variance forest
   (dec-B500). In the summary as a note the maintainer may veto. Alternative: leave it to `$setControl` alone.
7. `updateState = NULL` after the settings, `...` last; the setters return NULL invisibly. Alternative: the previous
   list, as `options()`, which `$getForest` already gives.
8. coefficient-prior-live stays its own item; until then a changed `amplitude` is refused and one equal to the value in
   force taken. Alternative: fold its ~400 lines in as a slice, delaying this item behind it.
9. `stated`'s new names are slots of slice 1a's records built by slice 2; forest 1's depth stays read from the model's
   tree prior, which every write to it updates. Alternative: the bridge reading forest 1 from the record, a creation
   read moved.
10. Removal of `forests = ` is R's unused-argument error, as dec-B508's removal was. Alternative: a message naming
    `$setForest`, for a spelling no release had.
11. A sharer's share is written at creation on every sharer, not only after a leave, so one slot says it in every
    state; the trio pins that this moves no draw. Alternative: NA until the first leave, which leaves a state to read
    two ways.
12. The last leave writes the model's leaf prior to the fixed default creation writes there. Alternative: relax the
    bridge's refusal for a model with no sharers on a flag saying they all left, a record of history a creation twin
    would not have.
13. Slice 3 after setmodel-parameters' slice 4, on the coordinator's order, its total pins through the real door.
    Alternative: slices 1 to 3 merged once before slice 4, slice 4 carrying those pins, which delays slice 4 by slice
    3's ~14 agent-hours and removes main's interim refusal.
14. `share` on `ForestSpec` and not `ForestStatement`: a statement is what a forest's caller wrote, and no caller wrote
    this. Alternative: on the statement, where `setForestMapSd`'s copy would carry it past a leave.

## Calls made at slice 2

Each the builder's, none in a ruling; the alternative after it.

1. The reader's `sd` on a numeric forest carries no column name: the value `forest(sd = )` takes, whose named form is
   refused. Alternative: the column's name, as `getLeafPrior` reads it.
2. The reader's `sd` is NULL wherever no one positive finite number is in force, beyond the drawn k the plan names:
   `normal(k = Inf)`, a unit that is not positive. Alternative: the 0 or Inf, which the door refuses handed back.
3. `$getVarianceForest()` reads NULL on a sampler with none. Alternative: refused, as its setter is.
4. An inert call with `updateState = TRUE` stores nothing, writing nothing being the rule. Alternative: store, as
   `$setLeafPrior`'s inert path does.
5. A multinomial `s$setForest(1, power = )` is refused as position 2 is. Alternative: position 1 taken, which reads as
   one category.
6. `$getForest` on a lone forest and a sharer reads `amplitude = fixed()`, the plan's table. A lone forest refuses every
   other `amplitude` in creation's words.
7. What a forest has no name for in its door is "unused argument 'x': $setForest takes forest, sd, ...", not R's bare
   form, because `...` catches it. Alternative: R's own message, which cannot be had past `...`.
8. A multinomial `amplitude` at the door is refused "a multinomial model's forests have no coefficient: its forests take
   no 'amplitude'", and at creation `forest(amplitude = fixed())` on that entry is still taken and ignored, as on main,
   while `student()` and `normal()` are refused for having one forest. The two doors differ on `fixed()` only. Recommended
   for a later sweep: refuse it at creation by the door's words. Alternative: take the value in force at the door, of a
   name the kind does not have.

## Open calls

1. RULED, dec-B520 (2026-10-30). The door's name. Shown `setForest` (recommended), `updateForest` and `setForestPrior`,
   the maintainer: "Use (a)." `setForest`, `getForest`, `setVarianceForest`, `getVarianceForest`; the help says the
   trees are left alone.
2. RULED, dec-B521 (2026-10-30). The door on a multinomial sampler. The earlier background, that the engine shares one
   specification so no per-category change can be made, was wrong: slice 1's install takes a per-category value and
   moves draws. Per-category settings are simply not offered. Shown the door taking the forests together, the
   maintainer: "IF we use your recommendation, will the interface behave similarly to other responses? Will a user know
   what they're getting?", then, with the additions in [Multinomial models](#multinomial-models): "OK, use that."
3. RULED, dec-B522 (2026-10-30). The spelling that returns a forest to the level above: none. The maintainer, of NULL:
   "That seems confusing. Will users intuitively understand that? Is there a simpler form that captures the same
   intent?", then, shown no return spelling: "Use your revised recommendation." With dec-B519 (a written NULL is taken
   where the signature's default is NULL, reading as not written; "Yeah, use your revised recommendation.").
4. RULED, dec-B524 (2026-10-30). A variance forest stating no depth and a live change: it keeps its value. The
   maintainer: "Will that confuse users? Personally, and maybe this is just because I tend to know how these things are
   coded, I wouldn't assume that the objects were tied together just because they shared something at creation. Feel
   free to push back and tell me if I'm wrong.", then "Yes, use your revised recommendation."
5. RULED, dec-B523 (2026-10-30). One setting of one forest from two doors. The maintainer: "I think the later write
   wins? No one would expect otherwise, would they?"
6. RULED, dec-B525 (2026-10-30). A sharer stating its own size live: taken, the others left where they are. The
   maintainer: "I think you could allow it, and leave the others where they are. Option (c). It would match my
   understanding: we size you at creation and after that it's your object to play with." and "It would behave similarly
   to how we set things from the data when construction priors."

7. RULED, dec-B535 (2026-10-30). How `$getForest(2)` on a multinomial sampler says it covers every category: a help
   sentence alone. The maintainer: "A help sentence only is fine, but if we're only going to accept factor levels as
   arguments isn't a little redundant? What am I missing here?"
8. RULED, dec-B538 (2026-10-30). A value handed back that a forest took from above changes nothing, statedness
   included; dec-B522's entry corrected. The maintainer: "Use (a)."
9. RULED, dec-B539 (2026-10-30). `amplitude = normal()` on a forest with no basis. The maintainer: "Offer a drawn normal
   coefficient here - I like the consistency."
10. RULED, dec-B540 (2026-10-30), its "median" being the forest's sd, which the trees carry by dec-B550. A
    fitting-function door on a forest with no basis that draws a coefficient and is the only candidate: it writes that forest's
    sd, reading the spelling from creation and not its calibration. The
    maintainer: "I would point out that we do have some differences between calls creation and after. Some things (e.g.
    data-drive defaults, diving the total sd) only make sense at creation. As long your recommendation is sensitive to
    that distinction, we can use it."

## Not verified

- Slice 3's engine is from reading forests-shared-total as built, not a prototype: every row of its table is a reading.
  Slice F is built and slice 2 ran its lone sharer's paths.
- The stated-sizes twin's 1e-8 is a judgment: rounding could flip an accept within a short run on some seed, and the
  builder fixes the seed and says so.
- control-migration's slice 1a is built and awaits a ruling on its direct-door rows, with no effect on what this plan
  reads.
- The twin oracle's exemptions beyond the held-to-drawn handoff were checked on the cases in c4 only; creation may
  consume random numbers by setting in a way an install cannot reproduce for a kind not run there.
- The design pass's timing (about 100 microseconds a call by label, a sweep about 350) is for `setForestBasis`;
  `$setForest`'s checks were not timed.
