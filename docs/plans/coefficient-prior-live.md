# coefficient-prior-live: a live sampler changes a forest's coefficient between held and drawn, and between normal() and student(df)

Status: LANDED 2026-11-05, slices 1 and 2 on one branch, reviewed LAND with a fix round (dec-A261, [Landing
note](#landing-note)). PLANNED 2026-11-03 on forest-door-s3 (64c8b0816), for dec-B511; amended the same day after its blind critique
(scratch/cpl-crit/report.md, probes c1.R to c5.R there), verdict BUILD AFTER AMENDMENTS, with the coordinator's calls on its
amendments 2, 3, 4 and 6 ([Amendments after the critique](#amendments-after-the-critique)). Rulings, quoted where
they settle a call: dec-B511 (the coefficient's prior live), dec-B518 (the door, `$setForest(amplitude = )`, refusing only where
creation refuses the same), dec-B550 (the trees carry `sd`; a drawn coefficient is a multiplier with prior median 1, `fixed()`
exactly 1; a write of `sd` moves no leaf, coefficient or auxiliary), dec-B539 (`normal()` on a forest with no basis), dec-B510
(`student()`), dec-B513 (`fixed()` on a lone forest), dec-B522 (no return spelling), dec-B524 (the variance forest), dec-B525 (a
share set once), dec-B538 (the value in force handed back), dec-B540 (`$setLeafPrior` on the one forest that draws), dec-B254
(parameters, never structure), with dec-B247 (a tree count is structure), dec-B354 (a held coefficient is 1, or 0 and 1 on a
factor), dec-B509 (a forest with no basis held by default), dec-B519 (a written NULL), dec-B523 (later write wins), dec-B515 (a
setter means what creation means), dec-B418 (a state installs as stored). The runs behind it: scratch/cpl-plan/runs, p1.R to
p11.R, on a library of 64c8b0816 and on a prototype of the engine write (scratch/cpl-plan/prototype.diff, about 80 lines, not
committed), and the critique's.

agent: slice 1 one opus implementer (effort high); slice 2 the same or a fresh one; one opus review over both slices that runs
the mutants below, the same reviewer checking the fix round unless the fixes change engine behavior. Tier "Changes draws"
([Process by risk](README.md#process-by-risk)): an engine write between sweeps.
rng: neutral, bit for bit, for every call taken today: the equivalence trio and the four snapshot files do not move, a
`$setForest(amplitude = )` naming the prior in force writes nothing, as now, and every call forest-door's slice 2 takes compares
its names as it does now. Every new call changes the model it is asked to, from the next sweep, and draws no random number. Each
is held to a creation twin, bitwise before any run where the start is creation's, and otherwise to the state twin of
[RNG class](#rng-class). The `$setModel` refusals of [Refusals](#refusals) turn calls taken today, and ignored, into errors.
window: for 1.0-0 (dec-B511). After forest-door's slice 3 lands; one branch, two slices, merged once.
budget: ~1,530 lines in the tree, ~300 of them package code: engine ~115 (chain ~80, combiner ~20, sampler and facade ~15),
bridge ~25, R ~160 (the door ~95, records ~25, `$setModel`'s refusals ~30, the `$setLeafPrior` words ~10), help ~40, NEWS ~5,
design notes ~20, tests/cpp ~380, tinytest ~640 (a new file ~590, 5 refusal pins respelled, the `$setModel` and `$setLeafPrior`
pins ~50), exact arms ~110, records ~40. Agent time, estimated and not measured: slice 1 ~7.5, slice 2 ~11, review with mutants
and a fix round ~6, gates ~3, rebase onto main after forest-door ~1: ~28.5. Plans here run 1.5 to 2 times low, so 43 to 57.

## Summary for the maintainer

You ruled that on a live sampler a forest's coefficient can be switched between held and drawn, between the normal prior and the
Student-t, and a t's degrees of freedom changed, each from the next sweep. You also ruled the spelling: the forest door's
`amplitude`. This plan builds it.

What a call does. `fixed()` holds the coefficient at its held value, 1, or 0 and 1 on a factor's two levels, at once.
`normal()` or `student(df)` on a held forest releases it, and the chain draws it from the held value. Between `normal()` and
`student(df)`, and between two degrees of freedom, the coefficient keeps the value it has. No call moves a leaf, and no call
moves the trees' scale: the forest's size stays what it was. The t's hidden variance starts at the normal's, as you ruled; going
back to the normal resets it; a new df keeps it.

The spelling that returns a coefficient to the normal or releases a hold needs no design: your ruling that the door has no
return spelling already says "a hold is released by writing a shape". A coefficient has no level above it, since no fitting
function states one, so every default is a value a forest can be given: `normal()` on a forest with a basis, `fixed()` on one
with none. Leaving `amplitude` out, or writing NULL, changes nothing.

Where it is refused. Only where creation refuses the same statement, with creation's reason: a model of one forest and a
multinomial model take only `fixed()`, which changes nothing; two or more forests still sharing one prior total refuse a drawn
coefficient unless the same call states the forest's `sd`, which takes it out of the group; the variance forest has no
coefficient.

Forests whose size comes from the fitting function's leaf prior. A forest with no basis that holds the whole of that prior,
alone of its kind or beside forests stating their own size, takes the size creation would give it when released; where k was
drawn that is the size at the prior's median k, so its size can jump, and the run then reports no k. Held again, it goes back to
following the leaf prior. The last forest left of a group, holding only its part of the total, keeps that part when released:
the part is written down as its own size, so that its statement as the reader gives it creates the same model.

As you already accepted: holding a factor changes the model, rows at the first level losing the forest; a hold at a large `sd`
jumps; a release from a held chain at a large `sd` takes a while to settle (about 200 sweeps at an `sd` of 7 in our runs, none
visible at 0.7).

Defects found and fixed here. `$setModel` handed a model whose forest records a different coefficient, or a different tree
count, takes it today and changes nothing, silently (the sampler keeps its trees). The coefficient will be refused naming
`$setForest`, the tree count as part of the model, naming a new sampler. And `$setLeafPrior`'s refusal on a sampler where every
forest with no basis states its own size or a coefficient's prior says only "its own sd"; its words are corrected.

There are no open calls. Nothing changes for bartCause, stan4bart, treatSens or bairrtt.

## Goal

`$setForest(f, amplitude = )` takes `fixed()`, `normal()` and `student(df)` wherever creation takes them on that forest, changes
the coefficient's hold and prior from the next sweep as [Every transition](#every-transition) says, records the change where a
copy and a reload read it, and leaves every leaf and the trees' scale. `$getForest(f)` reads it back, handing the list back
changes nothing, and the list creates the same model. `$setModel` refuses every forest-record difference it now drops in silence.

## Context

- The door as built: [`writeForestDoor`](../../R/forestDoor.R), whose [`checkForestStatement`](../../R/forestDoor.R) ends in
  today's refusal, "a forest's coefficient prior cannot be changed live yet"; [`changedForestSettings`](../../R/forestDoor.R),
  which drops a name equal to the value in force; the reader [`forestCoefficientInForce`](../../R/forestDoor.R). forest-door's
  [Slice 3: a sharer leaves the group](forest-door.md#slice-3-a-sharer-leaves-the-group) and its calls 8 and 14.
- The engine as forest-amplitude-normal left it: [`forestMapEntry`](../../src/bartcore/combiner.hpp) resolves a statement
  into the node factor and the coefficient's `ForestAmplitudePrior`; the node factor is the multiple m on every row but a
  sharer's, held or drawn, so a hold and a release cross nothing between the two
  ([Every live path](forest-amplitude-normal.md#every-live-path), [The start](forest-amplitude-normal.md#the-start)).
  [`holdAtKindValue`](../../src/bartcore/combiner.hpp) sets a block to its kind's value;
  [`drawAmplitudes`](../../src/bartcore/combiner.hpp) refreshes a mixture's variance v right after its block;
  [`restoreGlue`](../../src/bartcore/combiner.hpp) installs a block where it is drawn and v where it mixes.
- A forest with no basis on the leaf prior: [`forestSharers`](../../src/bartcore/combiner.hpp) and the K-forest constructor
  decide a group, a share and the held entry at k 2; [`Chain::installFittingLeafPrior`](../../src/bartcore/chain.hpp) moves the
  one candidate both ways between the k path ([`Chain::resolveLoneSharer`](../../src/bartcore/chain.hpp)) and a coefficient's df
  (its size from [`Chain::leafPriorSize`](../../src/bartcore/chain.hpp)).
- `$setModel`'s records: [`keepModelRecords`](../../R/records.R) keeps the sampler's forest and variance records and takes only
  the given `moves`, so a field no comparison reads is dropped.
- docs/design/forest-multiplier.md, docs/design/multiplier-combiner.md, docs/plans/coefficient-prior-opt-in.md (open call 3).

## Constraints

- No change to inst/include/dbarts/dbarts.h, the state format or `stateFormatVersion`. One facade virtual: `--preclean`.
- The engine sees numbers: a held flag, a df (NaN none, +inf the normal's internal code, finite a t) and an sd. The bridge
  converts; R decides what is refused.
- Out of scope: a start by kind (TODO drawn-coefficient-start-by-kind); a held coefficient in a state (TODO
  held-coefficients-in-state); a value other than 1 in `fixed()` (dec-B354); several columns or three levels (dec-B281,
  dec-B282); a forest rejoining a group (dec-B525); `$setLeafPrior` reaching a forest that states a df beside others stating
  their sd (it refuses there today, as creation's level has no forest to reach; only its words change).

## What exists today

Run on a 64c8b0816 build (p1.R, p2.R, p7.R, p8.R, p11.R; the critique's c1.R to c3.R):

- Creation takes `fixed()`, `normal()`, `student(1)` and `student(3)` on a factor, a numeric column, a forest with no basis
  stating `sd`, and a forest with no basis stating no `sd` where it is the only such forest: alone of its kind beside forests
  with a basis, or beside forests with no basis that state their `sd`, with or without a basis forest. It refuses `normal()` and
  `student()` on a model of one forest ("this model has one forest, which has no coefficient to give a prior") and where two or
  more forests with no basis state no `sd` ("these forests share one total: state each forest's 'sd' to give it its own"), and
  `student()` with no `sd` where the leaf prior gives no size (`k = chi(1.25, Inf)`: "which here gives no size").
- That forest under gaussian's default (k 2, held) has `sd` 1.415 held and 1.415 under `student(1)`; under probit's default (a
  drawn k) `sd` NA held and 1.574 under `student(1)`, the size at the chi prior's median k. Beside a forest stating its `sd` under
  probit it is a sharer of share 1 with a drawn k, not the held entry (c4.R).
- The door takes the prior in force and refuses any other: "a forest's coefficient prior cannot be changed live yet: forest 2
  keeps normal(); make a new sampler".
- `$setModel` takes a model whose forest record differs in `held.amplitude`, `amplitude.df` or `n.trees` (a forest's, on every
  sampler, and the variance forest's) and changes nothing: the record and the engine keep the old values, the sampler its trees
  (p8.R, p11.R, c2.R). It refuses `vars`, `label`, `basis.term`, `interactions`, `blocks`, `stated`, `sd`, `base`, `power` and the
  variance forest's `columns`, and a differing `share` except on the one candidate (forest-door's call 18).
- `$setLeafPrior` on `yb ~ forest(x1, sd = 1) + forest(x2, amplitude = student(3))` is refused saying "every forest with no
  basis stating its own sd", which forest 2 does not (c3.R).
- The model reorders forests with no basis first: `y ~ forest(x2, basis = ~ z) + forest(x1)` has the factor at position 2.

## Every transition

Kinds, each read off the forest's own state and not off how many forests have no basis:

- a factor; a numeric column;
- a forest with no basis stating `sd`, a forest that left its group among them;
- a whole-total forest: a forest with no basis stating no `sd` and no df that holds the whole of the fitting function's leaf
  prior, a group of one at share 1 or the held entry at k 2, whatever other forests the model has (`y ~ forest(x1) + forest(x2,
  basis = ~ z)`, `y ~ forest(x1, sd = 1) + forest(x2)`, with or without a basis forest); and, once released, the same forest
  stating a df and no `sd`;
- a sharer of a group of two or more;
- the last left of a group: a group of one at a share below 1, the others having stated their `sd`;
- a model of one forest; a multinomial model; the variance forest.

| kind | hold, `fixed()` | release, `normal()` or `student(df)` | `normal()` and `student(df)`, a df changed |
|---|---|---|---|
| factor | taken | taken | taken |
| numeric | taken | taken | taken |
| no basis, stating `sd` | taken | taken | taken |
| whole-total | taken: back on the leaf prior | taken: its size the leaf prior's, as creation gives it; refused where that gives none | taken |
| a sharer of two or more | the value in force | refused, as at creation; with `sd` in the call, taken: it leaves the group and is drawn | not reachable |
| the last left of a group | the value in force | taken: its part in force written as its `sd`, at the chi prior's median k where k is drawn | not reachable |
| a model of one forest | the value in force | refused, as at creation | not reachable |
| multinomial | the value in force | refused, as at creation | not reachable |
| variance forest | not a name | not a name | not a name |

What each moves, at the call, on every chain:

| transition | coefficient block | t's mixing variance v | leaves | trees' scale |
|---|---|---|---|---|
| hold | set to its kind's value: 1, or (0, 1) on a factor | none: the held prior's | unmoved | unmoved |
| release to `normal()` | stays at the held value, drawn from there next sweep | none: the normal's fixed variance | unmoved | unmoved |
| release to `student(df)` | the same | starts at the normal's variance: kappa / 0.674^2 with a basis, 1 / z^2 with none | unmoved | unmoved |
| `normal()` to `student(df)` | unmoved | starts at the normal's variance, which is the variance in force | unmoved | unmoved |
| `student(df)` to `normal()` | unmoved | reset: the normal's fixed variance | unmoved | unmoved |
| df to another df | unmoved | kept | unmoved | unmoved |

z is the normal's upper quartile and kappa 1/2 on a level, 1 on a numeric column ([`forestMapEntry`](../../src/bartcore/combiner.hpp)'s
table in forest-amplitude-normal). "Unmoved" for the trees' scale holds bitwise: the node factor is m on held and drawn rows alike
(`k.scale` identical to a creation twin's on seven cases, p5.R). The exceptions are the forests whose size the leaf prior gives,
each by creation's meaning of the statement (dec-B515) or, where creation would re-divide, by dec-B525's share set once:

- A whole-total forest released leaves the k path: its k pinned at 1, its trees' scale the size the fitting function's leaf
  prior gives at the scale creation set, as creation and `$setLeafPrior` read it (dec-B540): the same number at a held k, the size
  at the chi prior's median k where k is drawn. A leaf prior giving no size is refused with creation's words, naming `sd`. Where
  k was drawn and no other forest draws it, the run's k column goes off, as after a group's last leave (forest-door slice 3).
- Held again, it follows the leaf prior again, as creation of the merged statement places it: a group of one at share 1, or the
  held entry at k 2 only where creation gives that entry, a drawn k starting at the k that gives the spread in force, as every
  handoff into a drawn k does ([`Chain::resolveLoneSharer`](../../src/bartcore/chain.hpp)'s `sdWas`).
- The last left of a group released states its part in force as its `sd`: sqrt(share) times the leaf prior's size, the same
  number at a held k and at the chi prior's median k where k is drawn. Creation takes the merged statement but would give the
  forest the whole total; the written `sd` keeps the part (dec-B525), and `$getForest`'s list then creates the same model. Held
  again, it states that `sd` and stays at it.

A release takes effect on the coefficient's next draw; a hold takes effect at the call, so `getForestAmplitudes()` reads the held
value at once. A hold on a factor changes the model: rows at the first level lose the forest's part (dec-B550's cost). A hold at a
large `sd` jumps and a release at a large `sd` from a held chain settles slowly (dec-B550's costs; measured below).

### Names in one call

A name is dropped as no change where it equals the value in force after the call's other names apply, not before
([`changedForestSettings`](../../R/forestDoor.R) compares before today). Only `sd` reads another name: where `amplitude`
changes, the `sd` given is compared with the `sd` the forest would have under the new coefficient. Under probit, after a release
of a whole-total forest, `$getForest(1)$sd` is 1.574; `$setForest(1, sd = 1.574, amplitude = fixed())` then holds the forest at
`sd` 1.574 with k pinned at 1, as `forest(x1, sd = 1.574)` creates it, where today's comparison would drop the `sd` and return it
to a drawn k. Depth, moves and an `sd` beside an unchanged `amplitude` read no other name, and forest-door's slice 2 never took a
changed `amplitude`, so no call it takes compares differently; `$setForest(2, sd = <in force>, power = 3)` still drops the `sd`
and writes `power` alone, a pin. The round trip stays inert: every name equals the value in force after the others apply.

### Measured: where v starts, and the release transient

The TODO asked whether v's start shows. It does not: v is refreshed straight after its block's first draw, so its start reaches
one draw. A release to `student(3)` on a forest with no basis, after 100 held sweeps, 16 seeds, the fit's RMSE against a long
drawn chain's posterior mean over that chain's late RMSE, sweeps 1-5, 6-20, 21-80 (p6.R, p6b.R):

| | v at 1 / z^2 (ruled) | v at 1 / m^2 (creation's) | a chain created drawn |
|---|---|---|---|
| gaussian, `sd` 0.7 | 0.96 0.99 1.03 | 0.96 0.99 1.03 | 1.07 1.05 1.03 |
| gaussian, `sd` 7 | 1.35 1.30 1.20 | 1.34 1.30 1.22 | 1.02 1.00 0.99 |
| probit, `sd` 0.7 | 1.03 1.00 1.07 | 1.03 1.00 1.06 | 1.01 1.13 1.08 |
| probit, `sd` 7 | 2.16 2.13 1.92 | 2.05 2.05 1.80 | 0.99 0.98 0.99 |

At `sd` 7, sweeps 81-200 and 201-400 read 1.04 and 0.97 gaussian, 1.18 and 1.01 probit. The transient is the held chain sitting
in another posterior, not the start of v. The help says a write changes the model and no more (dec-B550).

## The spelling

dec-B522: "forest(amplitude = normal()) is taken, at creation and at the door ... so the normal coefficient prior has a value to
write; a hold is released by writing a shape." So `normal()` returns a coefficient to the normal and releases a hold, `student(df)`
releases one, `fixed()` holds. A coefficient has no level above it: no fitting function takes an `amplitude`, so no door of the
fitting function reaches it and no statedness is recorded for it. Every default is a value of the forest's kind, `normal()` with a
basis and `fixed()` with none, so a forest is returned to its default by writing it, and the help says so: `normal()` writes a
basis forest back to its default. An `amplitude` left out or NULL changes nothing (dec-B519); the value in force changes nothing,
statedness included (dec-B538), so `$getForest`'s list goes back inert.

## Refusals

Every check before any write, on the merged statement, by the validators creation uses. dec-B518: the door refuses where creation
refuses the same statement.

1. A model of one forest: `normal()` and `student()` in creation's words; `fixed()` taken, changing nothing (dec-B513).
2. Multinomial: as today, `fixed()` alone.
3. A sharer of a group of two or more, with no `sd` in the call: creation's words, "forest 1: amplitude = student() draws one
   forest's coefficient; these forests share one total: state each forest's 'sd' to give it its own".
4. A whole-total forest released with no `sd` where the leaf prior gives no size: creation's words. The last left of a group
   under that prior is refused the same way, its part having no size either.
5. A df below 1 or not finite, `student(Inf)` (naming `normal()`), `fixed(v)` with v not 1, a bare `student`: the validator's
   words, as today.
6. `$setModel`, every forest-record field [`keepModelRecords`](../../R/records.R) drops in silence:
   - a forest's `held.amplitude` or `amplitude.df`: through [`refuseOwnDifference`](../../R/fittingLevel.R), naming
     `$setForest(f, amplitude = ...)` in the value the model states; a model saved before a call and handed back after it is
     refused the same way;
   - a forest's `n.trees`, on every sampler a model of one forest included, and the variance forest's: "$setModel cannot change
     forest 2's 'n.trees': a forest's tree count is part of the model; make a new sampler instead", in the words of the other
     structure refusals (dec-B247, dec-B254);
   - the one candidate's `share` stays unread, as forest-door's call 18 has it: the fitting function's leaf prior reaches that
     forest whatever it states, so a model returning it to the leaf prior carries another share. `moves` stays taken.
7. `$setLeafPrior` where no forest takes the fitting function's leaf prior: the refusal's words become "every forest with no basis
   stating its own sd or a coefficient's prior" where any states a df, so they are true of c3.R's sampler. Behavior unchanged.

The `$setModel` refusals belong to setmodel-parameters' scope (dec-B254's refusal of structure by name, and that plan's
"a forest's own difference is refused naming setForest"); they are built here because this item makes the coefficient a door's
and the same comparison serves the tree count, ~30 lines.

## RNG class

- Inert: a call whose `amplitude` is the prior in force writes nothing and moves no draw, as now; so does the round trip.
- Twin, bitwise before any run, two chains: the next draws are a creation's with the merged statement. Holds by the prototype on
  every kind, a hold under probit, a factor standing second, `normal()` and `student(df)` both ways and a df changed on a basis
  forest, and `student(df)` to `normal()` and df to df on a forest with no basis stating `sd` (p4.R, p9.R, p10.R). A hold of a
  whole-total forest beside a forest stating `sd` is held to its creation twin under gaussian and probit; the prototype's route 1
  fails the probit one (c4.R), which routes 2 and 3 exist for.
- Exempt by name, where the start cannot be creation's: a release, whose block starts at the held value where creation starts
  at its start scale, 1 / 0.674 on a basis forest and 1 / m with none, a factor first by position; and `normal()` to
  `student(df)` on a forest with no basis, whose v starts at 1 / z^2 where creation's starts at 1 / m^2. Both are ruled (dec-B550:
  "a release draws from 1"; dec-B511). A release before the first sweep therefore is not creation's chain either; it is the same
  model, and a start at creation's value would be a second rule keyed on whether the sampler has run.
- The state twin, for every call, exempt or not: after the call, store the live sampler's state and install it into the live
  sampler and into a creation with the merged statement; their next draws are bitwise equal. It ties the model the call leaves to
  creation's and says nothing of the block and v the call set, which value pins read from the stored state's `glue`. Bitwise on
  seven cases after five sweeps, holds, releases, a whole-total forest released, a df changed and `student(df)` to `normal()`
  (p5.R).
- The reader-to-creation twin, for the last left of a group released: a creation in which every forest states what `$getForest`
  reads of it, installed with the live sampler's state as above, draws bitwise as the live sampler.
- With `sd` in the same call on a sharer of two or more: the recorded twin (`new("dbartsSampler", control, model, data)` from the
  live sampler's records) and the stated-sizes twin of forest-door slice 3.
- Later writes against a twin: after each transition, a `$setForest(f, sd = )` and a `$setLeafPrior(normal(k = 3))` (or its
  refusal) on the live sampler and on its creation twin give the same records, refusals and draws, so a statement left stale in
  the engine resolves the next write from the old df and fails it.
- After a run: `copy()` and a `saveRDS` and `readRDS` after `storeState` continue `identical()` to the live sampler.
- No call draws: every write sets a prior, a block or a variance.

## Copy, reload and state

A call writes the forest's record as creation with the merged statement writes it (`held.amplitude`, `amplitude.df` by
[`forestStatements`](../../R/model.R)'s rule, NA for `normal()` on a basis forest and Inf with none; `share` cleared where a
whole-total forest or the last left of a group is released, 1 where a whole-total forest is held again; `sd` where the last left
is released), so `identical()` holds between the live records and a creation twin's. Re-creation reads the records and installs
the stored state; a held block is the destination's, a drawn one and a mixing v the state's
([`restoreGlue`](../../src/bartcore/combiner.hpp)), which is what the state twin runs.

A state carries no prior (dec-B195) and goes in as stored (dec-B418): one stored while held and installed after a release gives
the block its held value and v the held prior's variance, a valid start that the next sweep redraws. One stored while drawn and
installed after a hold leaves the block held (TODO held-coefficients-in-state). One stored while a whole-total forest was
released stores no k for it, its k being pinned, so installed after it is held again under a drawn k it leaves that k where route 3
started it, a state installing a k only where it holds one (the critique's re-check, d1.R: 2.622 before and after). No migration: only shipped versions are owed.

## Design

### Engine (slice 1)

1. [`ForestCombiner`](../../src/bartcore/combiner.hpp) gains `setForestAmplitudePrior(f, prior, kind, normalVariance)`, false
   by default; `AmplitudeForestCombiner` overrides it. It replaces the block's prior whole (dec-B550's "one combiner write").
   Held: [`holdAtKindValue`](../../src/bartcore/combiner.hpp). Drawn and mixing: v is the variance in force where the old prior
   was drawn and mixing, and `normalVariance` otherwise. Drawn and not mixing: the new prior's variance. The block is not touched
   on a release.
2. `Chain::setForestCoefficientPrior(f, held, df, sd)` and its predicate `takesForestCoefficientPrior`, beside
   [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp). Each route keys on forest f's own state (its statement, its entry's
   share, the held-entry flag), never on the number of forests with no basis. A finite `sd` is stated first through
   `setForestMapSd`, which takes a sharer out of its group. The statement's `updateAmplitude` and `amplitudeDf` are written (df
   NaN where held) on `MapForest::statement` before any route resolves, so every later write (`setForestMapSd`,
   `installFittingLeafPrior`, `restateForestBasisScale`) reads the new df. Then one of four routes:
   - route 1, a forest that states an `sd`, or has a basis: the entry resolved again by
     [`forestMapEntry`](../../src/bartcore/combiner.hpp), its node factor and multiple unchanged bit for bit (a pin, not a copy:
     the prototype copied them and the pin is what checks the map agrees);
   - route 2, a whole-total forest released (a sharer at share 1, or the held entry at k 2):
     [`Chain::installFittingLeafPrior`](../../src/bartcore/chain.hpp)'s df branch factored out and called: the entry resolved with
     the df, f out of `sharers_`, k 1 and held, `namedSd` NaN, the held-entry flag cleared, the leaf scale from `mapLeafScale`;
   - route 3, a forest with no basis stating a df and no `sd`, held: placed as the K-forest constructor places a lone sharer of
     that model, as a group of one at share 1 or the held entry at k 2 exactly where creation gives that entry
     ([`Chain::resolveLoneSharer`](../../src/bartcore/chain.hpp) decides the entry from the options alone, so the decision is
     factored out of the constructor and both call it), a drawn k from the spread in force;
   - route 4, the last left of a group released: its part in force, sqrt(share) times `leafPriorSize()` in the response's units,
     computed by the engine from its own entry and stated through `setForestMapSd`, then route 1. R hands no number: it reads
     the `sd` back from the engine's row for the record, so no R double is compared with a C++ one.
   Then the combiner write, `normalVariance` being the variance a +inf df's entry carries. The predicate is false where f is no map
   forest, where a df is neither NaN, +inf nor at least 1, where a release names a sharer of a group of two or more with no `sd`,
   and where routes 2 and 4 find `leafPriorSize()` not positive and finite.
3. `Sampler::setForestCoefficientPrior`: every chain's predicate first, as [`Sampler::setForestMapSd`](../../src/bartcore/sampler.hpp)
   does, then every chain. A facade virtual and its forwarding.
4. Comments: the combiner constructor's "A held block is set here to its kind's value", which a hold now does live too, and
   the block above `forestMapEntry`, whose "a live write resolves the entry again and neither reaches the combiner" holds for the
   start scale only.

### Bridge (slice 1)

`bartcore_setForestCoefficient(ptr, forest, held, df, sd)`, registered beside
[`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp): sd NA where none; the unit check `bartcore_setForestSd` makes where
an sd is given; the engine's false as a backstop message. Nothing in R reaches it until slice 2, so slice 1 is neutral.

### R (slice 2)

1. [`checkForestStatement`](../../R/forestDoor.R): today's live refusal goes; refusals 3 and 4 come in, reading each forest's
   share, kind and statement from the engine's row and the records, as [`samplerSharers`](../../R/dbarts.R) does, never a count of
   forests with no basis; creation's words.
2. [`changedForestSettings`](../../R/forestDoor.R): where `amplitude` is among the changed names, `sd` is compared with the `sd`
   in force under the new coefficient ([Names in one call](#names-in-one-call)); otherwise as now.
3. [`writeForestDoor`](../../R/forestDoor.R): a changed `amplitude` goes to `bartcore_setForestCoefficient`, with a changed `sd`
   folded into the same entry so that one engine call checks both on every chain before writing either, and on the last left of a group with no
   `sd` of R's (route 4, the engine computing the part); an `sd` alone goes where it goes now. The engine order is moves, depth, then this, the undo of the
   two before it as now.
4. `modelWithCoefficient`, beside [`modelWithForestSd`](../../R/forestDoor.R): the record as creation writes it, by
   [`forestStatements`](../../R/model.R)'s rule and creation's share, and the last left's `sd` as the engine's row reads it after the call; where a whole-total forest is
   released, the leaf prior stays the model's, as creation with a df there leaves it (step 0 checks that against `identical()`).
5. `$setModel`'s refusals of [Refusals](#refusals) item 6, beside the `sd` and depth comparisons of
   [`levelDepths`](../../R/fittingLevel.R) and [`levelLeafPrior`](../../R/fittingLevel.R), through
   [`setForestCall`](../../R/fittingLevel.R), which learns to write an `amplitude`, and on a model of one forest in `$setModel`'s
   own path. Refusal 7's words in [`levelLeafPrior`](../../R/fittingLevel.R).
6. No reader changes: [`forestCoefficientInForce`](../../R/forestDoor.R) reads the engine's row, and `$getLeafPrior`'s `amplitude`
   and `prior.sd.of` follow. `print` reads the records.

## Slices

### Slice 1: the engine write (engine, bridge, tests/cpp)

Step 0: re-read forest-door slice 3 as landed: `setForestMapSd` on a sharer, the held-entry flag, the constructor's placement of a
lone sharer, `installFittingLeafPrior`'s df branch; amend the routes above to it. Then items 1 to 4, the bridge entry, and
tests/cpp. Neutral: no R path reaches it.

### Slice 2: the door (R, help, NEWS, tinytest, exact arms)

Step 0: run creation's records for a whole-total forest stating a df, alone of its kind and beside a forest stating `sd`, under
gaussian and probit defaults, against what item 4 writes. Then items 1 to 5, help, NEWS, design notes,
test-coefficient-prior-live.R, the five pins respelled, the `$setModel` and `$setLeafPrior` pins, the exact arms.

## Tests and gates

tests/cpp, on two chains at unequal tree counts, against a chain built from the merged statement where one exists:

- each transition of the table on each kind: the block, v, the entry's node factor and multiple bitwise, the leaf scale bitwise,
  the leaves untouched; a hold's twin; a release's block at the held value and v at the normal's variance, 2.1981 with no basis,
  1.1007 on a factor and 2.2013 on a numeric column; `student(df)` to `normal()` resetting v; df to df keeping it after a sweep;
- a whole-total forest, alone of its kind beside a basis forest and beside a forest stating `sd`, with and without a basis
  forest, released from the held entry at k 2 and from a drawn-k share of 1: out of `sharers_`, k 1, the leaf scale at
  `leafPriorSize()`; held again: placed as the constructor places it, the held entry only where creation gives it, a drawn k
  starting at the spread in force;
- the last left of a group of three released at a held and at a drawn k: its leaf scale its part in force (at the median k where
  drawn), its statement's `sd` that part, the other two untouched;
- the rescaling step after routes 2, 3 and 4: `forestRescaleApplies` declining with no draw where no forest is left on the k path,
  applying to the forest put back on it by route 3, its k and leaves scaled and every other forest read as offset;
- the statement written on `MapForest`: after each route, a `setForestMapSd` and an `installFittingLeafPrior` resolve as on a chain
  built with the new statement;
- the predicate's refusals, every chain checked before any write (one chain refusing leaves both as they were);
- `sd` and a release in one call on a sharer of three: the leave and the release, the others' shares and leaf scales unchanged;
- a store, a call, a restore: the block and v the state's, the prior the call's; a state stored while released installed after a
  hold-again, the drawn k left at route 3's start.

tinytest, test-coefficient-prior-live.R (~590 lines), every R twin on two chains:

- the twin before any run for every non-exempt transition, draws and records `identical()` to creation's;
- the state twin for every transition, the exempt ones included, with the block and v read off the stored `glue`;
- the reader-to-creation twin for the last left of a group released, at a held and a drawn k;
- after a run: the change in force through `$getForest`, `$getLeafPrior`, `getForestAmplitudes()`, `getK`; `copy()` and a reload
  `identical()` to the live sampler; the recorded twin;
- later writes after each transition, against the creation twin: `$setForest(f, sd = )` and `$setLeafPrior(normal(k = 3))`, same
  records, refusals and draws (c4b.R's case, `y ~ forest(x1, sd = 0.5) + forest(x2) + forest(x2, basis = ~ z)` with forest 2 held
  again, among them);
- inert: the round trip on every kind, `amplitude = NULL`, and the value in force, bitwise equal to no call;
- names in one call: under probit, after a release, `$setForest(1, sd = 1.574..., amplitude = fixed())` at the reader's `sd`
  holding at that `sd` with k pinned, against `forest(x1, sd = )`'s creation twin; `$setForest(2, sd = <in force>, power = 3)`
  writing `power` alone, records and draws as on forest-door's slice 2;
- whole-total forests under gaussian's and probit's defaults, alone of their kind and beside a forest stating `sd`, released and
  held again: `sd`, k, the run's k column, and a later `$setLeafPrior` reaching them as dec-B540 says;
- each refusal by its words, the sampler, records and draws unchanged after it, and the same statement at creation refused or
  taken alike;
- `$setModel`: a model with a differing coefficient refused naming `$setForest`; a differing `n.trees` refused as structure on a
  model of several, a model of one and the variance forest; the sampler's own model taken; `moves` still taken;
- `$setLeafPrior`'s corrected words on c3.R's sampler;
- the five pins of today's refusal respelled: "a forest's coefficient prior cannot be changed live yet" in test-forest-door.R
  and its siblings in test-forest-arguments.R, test-forest-amplitude-student.R and test-forest-amplitude-normal.R.

Exact arms (~110 lines). In benchmarks/R/forest-prior-exact.R, the pattern of its write arms (`writeFrom`, `writeMid`) with an
`amplitudeFrom`: a numeric forest created `fixed()` and released to `normal()` before the first sweep; a factor created drawn and
held after 50 sweeps; a numeric forest created `normal()` and moved to `student(3)` after 50 sweeps; a forest with no basis stating
`sd` created `student(3)` and moved to `normal()` after 50 sweeps. In benchmarks/R/probit-k-scale-exact.R, in place of the fifth
arm first planned (a factor released to `student(3)`, which the factor and numeric arms bracket): a variant of its `lone-factor`
arm, the frozen forest with no basis taking the whole total under the drawn k beside the held factor stump, released to `normal()`
and held again by `$setForest` after 50 sweeps, read against `lone-factor`'s own oracle; it runs routes 2 and 3 and the k's restart at the spread in force, the rescaling
step included. Each reads the oracle of the statement in force at the end. The state twin ties each call's model to a creation the
existing arms cover; these arms check that a chain after a live write targets it, v's and k's handling included.

Gates. Slice 1: tests/cpp plain and sanitized, `--preclean`, `Rscript benchmarks/R/mutation-battery.R verify-anchors`. Slice 2,
over both: the R-loaded sanitizer run over the new file; the equivalence trio bitwise against benchmarks/baselines/MANIFEST's
baselines; the four snapshot files unchanged on the reference build; the full suite once in one process; lint against the slice's
library, air, check-rc-codoc, check-win-drift, check-doc-freshness; `R CMD check --as-cran` from a clean tarball; the exact gates of
.github/workflows/exact-gates.yaml in `quick` (the records change what a sampler carries), and forest-prior-exact.R and
probit-k-scale-exact.R in full for the new arms; the x86 leak-sanitizer leg before landing.

## Mutants

Each to fail a named pin.

1. A hold leaving the block where it was: the hold twin and the block pin.
2. A release setting the block to creation's start: the block pin (the state twin does not see it).
3. A release to `student(df)` starting v at the held variance, or at 1 / m^2: the v pin.
4. `student(df)` to `normal()` keeping v: the twin before a run and the exact arm.
5. A df changed resetting v, or ignored (the scale written, `mixtureDf` not): the v-kept pin; the factor's df twin.
6. The entry's node factor re-resolved off the in-force one: the bitwise leaf-scale and `k.scale` pins.
7. One chain written: the two-chain twins.
8. The record written and the engine not: the twin. The engine and not the record: copy and reload.
9. `amplitude.df` Inf for `normal()` on a basis forest: the records against the creation twin.
10. A whole-total forest released and left in `sharers_`: `getK`, the k column, the twin.
11. Routes 2 and 3 keyed on the number of forests with no basis: the probit hold-again beside a forest stating `sd` (c4.R) and the
    later `$setLeafPrior` (c4b.R).
12. Route 3 placing the held entry where creation gives a group of one at share 1: a hold-again under a held k other than 2
    (`leaf.prior = normal(k = 3)`) and under a named sd, beside a forest stating `sd`, each bitwise against its creation twin. The
    two placements agree to rounding there, so only a bitwise twin catches it. (First aimed at gaussian's default beside a forest
    stating `sd`, where creation gives the held entry too: the constructor's test counts no forests with no basis, d4.R.)
13. The combiner written and `MapForest::statement` not: the later-writes pins, and tests/cpp's `setForestMapSd` after a route.
14. `sd` and `amplitude` on a sharer written in two engine calls with the check between: the atomic refusal pin.
15. Route 4 resolving the last left at the whole total: the part's value pins, the sd the forest had under a held k and the
    prior's median k under a drawn one (the reader-to-creation twin states whatever part the reader gives, so it cannot see this).
    Its `sd` not recorded: the reader-to-creation twin.
16. `sd` compared before `amplitude` applies: the probit `sd = 1.574, amplitude = fixed()` pin. The sd under the coefficient
    taken as no sd under a held k too: the gaussian pin, the sd dropped under each held leaf prior.
17. `$setModel` ignoring a coefficient or a tree-count difference: its refusal pins.
18. A check on the names given rather than the merged statement: `sd` and `student(1)` on a sharer taken, `student(1)` alone
    refused.
19. The rescaling step reading a forest route 2 took off the k path: tests/cpp's rescale pins and the probit exact arm.

## Help and NEWS

- man/dbartsSampler-class.Rd, `setForest`: "any other change is refused for now" becomes what [Every transition](#every-transition)
  says in a few sentences: a hold sets the coefficient to its held value, a release draws it from there, the t's variance starts at
  the normal's, no leaf or trees' scale moves, a write changes the model; that `normal()` writes a basis forest back to its default
  and `fixed()` one with no basis; a forest sized by the leaf prior takes creation's size when released, the last left of a group
  its part as its `sd`; the refusals. `getForest` unchanged. `setModel`: a forest's own coefficient is refused naming `setForest`,
  a tree count as part of the model. `setLeafPrior`: nothing beyond the message.
- man/forest.Rd, `amplitude`: that a live sampler changes it through `$setForest(forest, amplitude = )`.
- inst/NEWS.Rd: the `$setForest` item's names already list `amplitude`; one clause that it switches a coefficient between held and
  drawn and between `normal()` and `student(df)`.
- docs/design/forest-multiplier.md section 1's "The prior is fixed at creation; a live change is its own item" becomes the rule;
  docs/design/multiplier-combiner.md gains the combiner write beside the basis swap.

## Consumers

Grepped read-only on each repository's checked-out branch: no consumer calls `$setForest` or `$setModel` with a coefficient or a
tree count. bartCause's `bcf()` states `amplitude = student(1)` or `fixed()` at creation; stan4bart reaches one forest through the
flat C API, which has no coefficient door; treatSens and bairrtt touch neither. No consumer changes.

## Interaction with the other pending plans, and order

1. forest-door slice 3 (branch forest-door-s3, in review): this branch builds on its tip and rebases onto main after it lands.
   Both edit [`writeForestDoor`](../../R/forestDoor.R), [`checkForestStatement`](../../R/forestDoor.R) and
   [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp)'s neighborhood. Serial. Its interaction row 6 and its call 8
   ("until coefficient-prior-live") are records to update at this landing.
2. setmodel-parameters (slice 4 landed in this base): the `$setModel` refusals extend its rule; its plan's records step notes them.
3. collapsed-move-main (engine, chain.hpp's k move): textual meeting at most; engine features merge serially.
4. drawn-coefficient-start-by-kind, not planned: would move creation's start, not a release's, which is the held value.
5. held-coefficients-in-state, not planned: a hold makes the case it names more common; nothing here changes it.
6. leaf-prior-writer-by-level and sbc-bcf-probit-burn: independent.

## Cost

- A release is never creation's chain, even before the first sweep: the same model from another start.
- A whole-total forest released under a drawn k takes the size at the prior's median k, a jump, and the run's k column can go; held
  again it draws k from the spread in force.
- The last left of a group, released, states its part as its `sd`, so it is reached by no fitting-function door afterwards that
  would not reach a forest stating `sd`.
- Costs dec-B550 ruled: a hold at a large `sd` jumps, a release at a large `sd` from a held chain is slow (about 200 sweeps at 7),
  holding a factor changes the model.
- Two `$setModel` calls taken today, and ignored, become errors.
- A facade virtual and a combiner virtual.

## Calls made

1. A release keeps the block at the held value at every point in a chain, before the first sweep included. Alternative: creation's
   start before any sweep, so the twin holds there, at the price of two rules keyed on whether the sampler has run.
2. A df changed keeps v: v is chain state valid under any df, as a parameter change leaves the state (dec-B254). Alternative:
   reset to the normal's variance, as a release starts it.
3. A hold takes effect at the call, the block set at once. Alternative: at the next sweep, which no reader could tell apart but
   `getForestAmplitudes()` read between the call and the sweep.
4. A whole-total forest, released, takes creation's size: the leaf prior's at the chi median where k is drawn (dec-B515,
   dec-B540). Alternative: refuse under a drawn k, naming `sd`, which refuses what creation takes.
5. The same forest, held again, follows the leaf prior again, placed as creation places it (dec-B513), its drawn k from the spread
   in force. Alternative: keep it at its size, stating it as its `sd`, which a creation twin does not do.
6. The last left of a group, released, states its part in force as its `sd`, at the chi median k where k is drawn, so that its
   statement as `$getForest` gives it creates the same model (dec-B518 takes what creation takes; dec-B525 sets a share once; the
   coordinator's call). Alternative: refuse it naming `sd`, as a sharer, the plan's first form, which refuses a statement creation
   takes.
7. `sd` and `amplitude` in one call go to one engine entry, every chain checked for both before either is written. Alternative:
   two entries, `sd` first, with R's checks the only guard between them.
8. A name is compared with the value in force after the call's other names apply (the coordinator's call). Alternative: before,
   as today, which drops an `sd` that the new coefficient gives a different meaning.
9. `$setModel` refuses every forest-record field it drops: the coefficient naming `$setForest`, the tree count as structure, the
   one candidate's share left unread as forest-door's call 18 left it. Built here, though setmodel-parameters' scope.
   Alternative: the coefficient alone, the tree count left to a TODO item.
10. Exact arms are kept though the state twin ties every call to a creation the gates cover, one moved to probit for routes 2 and
    3 under a drawn k. Alternative: none, as forest-door argued for its door.
11. Two slices, one review over both. Alternative: one slice, as forest-amplitude-normal was; the engine slice is neutral alone and
    reviews more easily against tests/cpp.

## Open calls

None. The TODO's one item to design, the spelling that returns a coefficient to the normal or releases a hold, is ruled by dec-B522
([The spelling](#the-spelling)). Every refusal is one creation makes of the same statement (dec-B518); the last left of a group is
taken by dec-B518 and sized by dec-B525.

## Amendments after the critique

The blind critique (scratch/cpl-crit/report.md) found the design sound and asked for eight amendments; the coordinator ruled 2, 3,
4 and 6. All are folded into the text above.

1. A whole-total forest is any forest with no basis holding the whole leaf prior, not only the only forest with no basis: routes 2
   and 3, the checks and the tests key on its own state (c4.R, c4b.R).
2. The last left of a group is released at its part in force, written as its `sd` (Calls made 6), in place of the refusal; "creation
   cannot reach this state" was wrong: creation takes the statement and would re-divide.
3. `$setModel` refuses every field `keepModelRecords` drops (c2.R, c2b.R, p11.R).
4. Names are compared after the call's other names apply, with the probit pin and a slice-2 call pinned unchanged.
5. The statement is written on `MapForest` before a route resolves; a mutant and the later-writes pins.
6. The fifth forest-prior-exact arm moves to probit-k-scale-exact.R; tests/cpp checks the rescaling step after the routes.
7. A state stored while released, installed after a hold-again: said in [Copy, reload and state](#copy-reload-and-state). The
   critique's "starts k at 1" was corrected by its re-check (d1.R): the state holds no k for a pinned forest, so k stays where it
   was.
8. The help names `normal()` as what writes a basis forest back to its default. And `$setLeafPrior`'s refusal words (c3.R),
   carried in slice 2.

## Not verified

- The prototype covered route 1 and the held entry at k 2; routes 2 to 4 (a whole-total forest under a drawn k, its hold back onto
  the k path, the last left released) and `sd` with `amplitude` on a sharer are from reading `installFittingLeafPrior`,
  `resolveLoneSharer`, `forestSharers` and `setForestMapSd`, not run. The probit hold-again failing on route 1 is the critique's
  run (c4.R).
- The state twin is bitwise after five sweeps on seven cases at one seed; a `storeState` and `setState` on one sampler alone moves a
  draw by 5.6e-14 after a run, so the twin must put both samplers through the install, as p5.R does.
- Records `identical()` to a creation twin's after a call were not run: the prototype wrote no record.
- The transient table is one design, 16 seeds; the exact arms' runtime in `quick` mode is estimated, about a minute an arm, and
  whether probit-k-scale-exact.R's `lone-factor` arm takes a `$setForest` coefficient write without an oracle change is a reading.

## Landing note

Landed 2026-11-05: slice 1 (the engine write, its bridge entry and tests/cpp) and slice 2 (the door, help, NEWS, tinytest and
the exact arms) on one branch, rebased onto main after the host-neutral engine's slice S1 with no conflict, merged once.

Review: one Opus reviewer over both slices ran the plan's mutants and two of its own: LAND, with three test gaps and one help
sentence. Every named mutant was killed but 14, which is equivalent (once the sd write is taken the coefficient write cannot
refuse). The fix round pinned the three gaps (an sd beside a hold under gaussian, a hold-again under a named sd, the probit part
at the median k), named mutant 15a's real killers, and narrowed setForest's sentence on what nothing returns to following the
fitting function to a stated sd or depth. A second Opus reviewer checked the fix round in place of the first: LAND, the three
mutants the gaps let through (R1, 12b, 15a) each killed by the new pins.

The builder dropped one planned refusal: `$setModel` refusing any differing tree count refused a count the model does not
state, which the sampler keeps (Calls made 9); a stated count is refused at both levels, as dec-B551 rules.

Gates on the final tip: tests/cpp 427 ok, plain and under ASan and UBSan with no reports; the full suite in one process, 40,995
passed, 0 failed; equivalence 55 of 55, bcf 19 of 19 and multinomial 11 of 11 bitwise on the reference build, none skipped; the
four snapshot files; the exact gates in quick, 45 of 45, and in full the new arms of forest-prior-exact.R (largest |z| 2.48) and
probit-k-scale-exact.R's lone-factor-released (largest |z| 1.9, 4e6 sweeps); R-loaded ASan and UBSan over the seven door files
before the fix round, no reports; R CMD check --as-cran with one NOTE, the date; lint, air, codoc, Windows drift, doc freshness and
the mutation anchors. The x86 leg (Linux, x86-64): install clean; tests/cpp 427 ok, plain and with ASan, UBSan and
LeakSanitizer; the full suite 41,123 passed, 0 failed; equivalence in statistical mode, bcf and multinomial --cross-host, all
within bounds. CI, all nine workflows, on the branch tip.

The calls are in dec-A261.
