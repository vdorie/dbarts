# forest-amplitude-normal: the trees carry a forest's `sd`, a drawn coefficient is a multiplier with median 1, and a drawn normal on a forest with no basis

Status: LANDED 2026-11-03 (dec-B562), one slice: built, reviewed LAND with its fix round. sbc's bcf arms under the rule dec-B557
set tripped on one functional of 70, bcf-probit's `abs.a`, which a diagnosis read as shared with the base; the maintainer landed it
on that reading; a longer-burn rerun after landing read the U as the Cauchy tail's mixing, shared by both libraries
([Landing note](#landing-note)). Planned
2026-10-30 (dec-B539), rewritten for dec-B550 on 555906d37: the trees carry a forest's `sd` on every forest, held or
drawn; a drawn coefficient is a multiplier with median 1, `fixed()` exactly 1; a write of `sd` sets the trees' scale and moves no
leaf, coefficient or auxiliary, before the first sweep too. Earlier blind critiques read other forms; a fresh one
read this form and this text takes its findings, each reread against the code at 7f53d1c78. Rulings: dec-B562 (landing), dec-B557
(the sbc rule), dec-B550, dec-B539 (the coefficient), with dec-B510 (`student()`, its words revised by dec-B550),
dec-B511 (the coefficient's prior live, its own item), dec-B513 (`fixed()` on a lone forest), dec-B522 (`normal()` has a value to
write, no return spelling), dec-B528 (a constant basis column refused), dec-B540 (`$setLeafPrior` on the one forest that draws).
The runs behind this form: a port of the ruled form onto slice F in scratch/fan-plan3 (patch-f.py, patch-f-r.py, fpe-f.R, runs/),
the ruling's own in scratch/fan-sw/out and scratch/fan-crit, and the critique's in scratch/fan-crit3/runs.

agent: one opus implementer (effort high), one opus review that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)).
rng: no draw is added, removed or reordered on any path. Every sampler with a drawn forest moves by rounding, so a seed can fork
where rounding flips a decision. Bitwise for every call with no drawn forest: held forests of every kind, sharers, a lone forest,
multinomial, the snapshot files. Posterior-defining for the one new model, `normal()` on a forest with no basis, gated by the exact
arms below. Moves by design: a write of `sd` on a drawn forest with no basis, before its first sweep or in mid-chain: the same model
at the written `sd`, gated by the exact arms that follow a write, and not the chain the base runs, which puts that write in a different
place. After 30 sweeps a mid-chain write there reads 2 to 6 fit sds from the base (scratch/fan-crit3/runs/attack1.R). One corner: a caller's own function named `normal`, called inside `forest()` or `forests =`, stops being reached
([Calls made](#calls-made) 6).
window: for 1.0-0, after forests-shared-total's slice F (dec-B550). Merged once, one slice.
budget: ~1,200 lines in the tree, ~165 of them package code: engine ~100 touched (59 measured on the port, then the comments,
`nodeScaleDivisor`'s removal and the finite start), R ~65 (39 measured, then the refusals, predicate and decoder), help and NEWS ~70,
design notes ~60, tests/cpp ~150 new and ~130 re-pinned (44 messages), tinytest ~260 new and ~95 re-pinned (42 sites and
`prior.sd.of`'s words), exact arms ~150 (138 measured, ~10 for the four that follow a write) and the other gates' restated reads ~60,
records ~25. Agent-hours, built up from the counts: engine and R 2, tests/cpp 5 (44 messages at ~5 minutes, and the new tests),
tinytest 5.5 (42 sites at ~5 minutes, and the new file), gates' oracles 1.5, help and notes 2, gate runs and re-records 2, sbc's bcf
arms 2 (restating the generator, running base and tree), the fork read 0.5, review of 39 mutants and a fix round 5: ~26, range 22 to 31.
The rates are guesses; the counts are not. Counts are on slice F's tip 7f53d1c78.

## Summary for the maintainer

You ruled that the trees carry `sd` on every forest (dec-B550). This item builds that and adds `forest(amplitude = normal())` on a
forest with no basis. After it, doubling `sd` doubles the trees' scale on every kind under every coefficient option; a drawn
coefficient has prior median 1; a write of `sd` changes the trees' scale and moves nothing already drawn.

Nothing changes in the model for what runs today. The prior on the coefficient times the trees is the same, and the coefficient's
start is set once at creation where the chain stays. So the chain is the same until a rounding flips a decision. Under gaussian and logistic: 24 runs of 4,000
sweeps, none past 1e-6, and 24 of 40,000, none past 1.1e-7; but a run of 20,000 sweeps can part there too, on sbc's first replication in 4 of 6 runs of
`bcf-weak`, 2 of 6 of `bcf-logistic`, 5 of 6 of `bcf-logistic-weak` and 6 of 6 of `bcf-probit`. Under probit it is the same for about 1,000 sweeps, and then a run can fork where rounding flips a decision:
a half-Cauchy forest with no basis forked in 4 of 4 runs, at sweeps 811 to 2,266; a numeric forest in 1 of 4, at 3,723; a factor in 0 of 4. A
forked chain is another run of the same chain: after the fork the fits differ draw by draw by 1.5 to 3.3 and not on average, over 288 seeds
z -0.39 on the grand mean. No exact arm crosses a fork, so that is a comparison of chains and not an exact reading: [Forks](#forks). The engine has about 100 lines touched and two things go: the mixture-scale write and the 0.674 node
divisor.

`sd` is the trees' prior scale, so it is the forest's size at a coefficient of 1, and not the prior median of the forest's contribution,
which is 0.54 `sd` under `normal()` since the contribution is a product of two random numbers (run: 0.541). The reader's `prior.sd.of` said
"(prior median)" on a drawn forest and now says "(at a coefficient of 1)": [Calls made](#calls-made) 11.

What it costs, each measured on a port to slice F:

- Raw reads change: a basis forest's coefficient 1.48 times larger and its trees that much smaller; a forest with no basis's
  coefficient divided by `sd` over the unit and its trees multiplied by it. Fits, contributions and bartCause's effects do not move.
- Seeded draws move by rounding wherever a forest is drawn, and a write of `sd` on a drawn forest with no basis leaves the base's
  chain by design. bcf-equivalence and one equivalence cell re-record, judged by the cross-host check, which passes 19 of 19 against
  both baselines.
- 88 of 32,947 tinytest results and 121 tests/cpp checks re-pin. All but 8 are units. None is a
  regression. Against main the counts were 121 and 123: slice F holds a forest with no basis that states no prior, so most writes
  before a first sweep in the tests stop reaching a drawn coefficient.
- bcf-exact-weak reads the coefficient alone and must have its oracle restated. benchmarks/R/sbc.R's generator draws the coefficient in
  the old units and must be restated from the rule; it is judged against the base's verdicts, since CI skips those arms.
- About 26 agent-hours, range 22 to 31.

Smaller things you would want to know:

- Inside the engine `normal()` is stated as a df of +inf, an internal code. No t is drawn, no auxiliary kept, no quantile taken, and
  users never see it.
- A function of your own named `normal`, called inside `forest()`, is no longer reached, as is already so for `fixed`, `student`
  and `invchi`. No released call is affected.
- dec-B528's refusal of a constant column gains a clause naming `forest(x2, sd = , amplitude = normal())`, words on a ruled refusal.
- forest.Rd's `indexSd` example breaks silently under this form: it would report 1.06 for a true 1.57, with no error. This item
  fixes it and adds the pin no test had, which reads the installed help.
- At an extreme `sd` the auxiliary's start would overflow, and the base runs there: the start is kept finite in the engine, with no
  refusal, since a start is not a model. [Calls made](#calls-made) 17.

## Goal

`forest(amplitude = normal())` on a forest with no basis draws its coefficient under a normal prior with median 1, taken where
`student(df)` is taken and refused where it is refused, with no `sd` taking the size `student(df)` takes there. On every forest,
held or drawn, the trees have prior scale `sd`. A forest's size is its coefficient times `sd`: `fixed()` is exactly 1, and
`normal()` and `student(df)` draw it with median 1, exactly on a forest with no basis and 1.0007 on a forest with a basis, which
keeps bcf's 0.674. A write of `sd` sets the trees' scale and nothing else. At the live door `amplitude` is the value in force only,
until coefficient-prior-live.

## What 0.9-34 did

It had no `forest()`, no `amplitude`, no coefficient and no model of several forests. Nothing is owed. `normal(k)` existed there
only as a leaf prior, which it still is at `leaf.prior`. A one-forest write moves no leaf there or on main, which is the shape this
item gives every forest.

## What exists on the base

Read on slice F's tip, 7f53d1c78, which this item builds on; run on a port of the ruled form onto it, against F's library.

- [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) has one conditional, `priorPrecision = 1 / prior.variance`, and does not
  read a forest's kind. [`drawAmplitudes`](../../src/bartcore/combiner.hpp) adds the variance's inverse-gamma draw only where
  `mixtureScale > 0`.
- What `sd` sets, by doubling it: the trees under `fixed()` on every kind, under the normal and `student(df)` on a basis forest
  (through 0.674), and on a sharer; the coefficient only under `student(df)` on a forest with no basis.
- [`forestMapEntry`](../../src/bartcore/combiner.hpp) is resolved at three places: creation,
  [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) and [`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp).
- [`restoreGlue`](../../src/bartcore/combiner.hpp) restores a block where it updates and a variance where it mixes. A fixed variance
  is model, not state, and stays the destination's.
- The product: [`forestMultiplier`](../../src/bartcore/combiner.hpp) contracts the block by the row; the fit and every reader of a
  contribution read the product. The row bound's `zeroMultiplierTolerance` is relative to the forest's largest, so it does not see
  a scale.
- On main and on F `forest(amplitude = normal())` is "could not find function". `student(Inf)` is refused.

### What slice F changes in the base

Slice F is on branch forests-shared-total, tip 7f53d1c78, built and in review, not landed. This item builds on it.

- A forest with no basis that states no df is held, so [`forestMapEntry`](../../src/bartcore/combiner.hpp)'s kind-none drawn rows
  are reached only by `student(df)` here and `normal()` after this item. Its hold keys on `isnan(df)`, which a +inf df passes.
- The map loses `family` and `oneForestSize` for `leafPriorSize`. No sd on a drawn forest with no basis takes the leaf prior's size.
- A sharer's entry is its own, `share` finite, node factor sqrt(share), and takes no part here.
- Start: where two or more forests have no basis a drawn one starts at 1 rather than by position.
- Baselines: bcf-equivalence-e346b018 and equivalence-e346b018 are the current ones.
- Pins: forests with no basis stating `sd` and no prior are held, so their write twins pass unchanged here. F's tests added 3 results
  and 6 messages that this item re-pins, and its `hold: an infinite df is stated, and is not held` flips: +inf draws after this.
- Its D8 sentence, written for dec-B542, is stale: under dec-B550 "reports k 1 and `k.scale` that sd" holds on every forest stating
  an sd, drawn or held, so no scope change is owed. Its plan's two dec-B542 mentions read dec-B550. This item's records step edits
  them, once F has landed.

## Did we build the wrong thing, and is there a simpler form?

No to the first: dec-B539's option was "a new prior in the engine with its own exact gate", and there is no new prior. The ruling's
form is itself the simplification of three earlier ones, and it removes more than it adds: no write rule, no write into a block, no
new state, no new virtual. The earlier prototype's `writeAmplitudePrior`, `kindStart_` and `startVariance_` are not built.

Two simpler forms fail, as the ruling records. A start of 1 on every kind with no scale moves every chain for no gain, and at a
large `sd` under a link a forest with no basis takes up to 80 sweeps to settle (probit, `normal()`, `sd` 6 read 1.61, 2.30, 2.40, 1.72
at sweeps 1-5, 6-20, 21-40, 41-80, and read 0.88, 0.90, 0.98, 0.97 with the scaled start). 0.674 replaced by z on a basis forest
moves every treatment forest's draws for a 0.07 percent reading.

### Same prior, same chain, other units

(a, F) on slice F is (a / g, g F) under this plan, g being 0.674 on a basis forest and the multiple with none: [The start](#the-start). The scale between a
coefficient and its trees is a gauge. Each block's conditional is equivariant in it, a row's zero bound is relative, a mixture's
variance scales by the square and its scale by the factor. So if the block starts at the rescaled value the two chains are one kernel,
and one chain until a rounding flips a decision ([Forks](#forks)).

## Design

m is the forest's multiple, `sd` over the unit, a numeric forest's default being 1 / basisScale and a default with no basis the leaf
prior's size; s is the anchor. c = 0.674, bcf's rounding; z = 0.6744897501960817, the normal's exact upper quartile. kappa is the
kind's variance, 1/2 on levels and 1 otherwise. t is t_nu,.75. df is the forest's `amplitudeDf`: NaN where none is stated, +inf the
normal's code, finite for `student(df)`.

### forestMapEntry

Every row is [`forestMapEntry`](../../src/bartcore/combiner.hpp) as it reads after this item. The node factor is m on every row
but a sharer's, which is the rule; the node divisor is gone.

| kind | coefficient | df | node factor | variance (v's start where it mixes) | mixture scale | mixture df | start scale | against slice F |
|---|---|---|---|---|---|---|---|---|
| any | held | any | m | kappa | 0 | NaN | 1 | unchanged, bit for bit |
| none | sharer | NaN | sqrt(share) | 1 | 0 | NaN | 1 | unchanged |
| none | held | NaN, no prior stated | m | 1 | 0 | NaN | 1 | unchanged |
| levels | drawn | NaN or +inf | m | kappa / c^2 | 0 | NaN | 1 / c | node was m / c, variance kappa |
| levels | drawn | nu | m | kappa / c^2 | sqrt(kappa) z / (t c) | nu | 1 / c | node was m / c, v's start kappa, scale without c |
| numeric | drawn | NaN or +inf | m | 1 / c^2 | 0 | NaN | 1 / c | node was m / c, variance 1 |
| numeric | drawn | nu | m | 1 / c^2 | z / (t c) | nu | 1 / c | node was m / c, v's start 1, scale without c |
| none | drawn | +inf | m | 1 / z^2 | 0 | NaN | 1 / m | new |
| none | drawn | 1 | m | 1 / m^2 | 1 | NaN, the half-Cauchy arm | 1 / m | node was 1, v's start 1, scale m |
| none | drawn | nu > 1 | m | 1 / m^2 | 1 / t | nu | 1 / m | node was 1, v's start 1, scale m / t |

- A coefficient's prior does not read `sd` on any row: the variance of a fixed-variance row and every mixture scale are constants of
  the kind and df. Only v's start and the start scale read m, and both are read once.
- The medians. Kind none: exactly 1 under `normal()`, `student(df)` and the half-Cauchy (run: `coefficient@0.5` +0.50029 against
  0.5 on a normal's prior arm). A basis forest: z / c = 1.0007 under each, for a numeric coefficient's size and a factor's
  difference, whose variance is 2 kappa / c^2.
- A basis forest's +inf resolves as its NaN does: forest-door's record does not use it and a fixture may. Below 1 and -inf stay
  refused.
- `sd` is one thing on every row: the trees' prior scale. The forest's size is the coefficient times `sd`. The median of abs(a F)
  under a drawn coefficient is smaller than `sd`, 0.54 `sd` under the normal, since a F is a product of two random numbers, and
  sd(a F) under the normal is 1 / z = 1.48 `sd` against `sd` held; the help says so once.
- The entry carries `startScale`, beside `amplitudePrior` and not in it (dec-B550: a field separate from the prior).

### The start

Each block starts at what it starts at on slice F, times its start scale, once, in
[`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor, after
[`holdAtKindValue`](../../src/bartcore/combiner.hpp) and the position rule.

| block | slice F | after |
|---|---|---|
| numeric, drawn | 1 | 1 / c |
| levels, drawn, standing second | (0, 1) | (0, 1 / c) |
| levels, drawn, standing first or third | by position, (1, 1) | (1 / c, 1 / c) |
| none, drawn, standing first | 1 | 1 / m |
| none, drawn, a basis forest standing first, the only forest with no basis | 0 | 0 |
| none, drawn, two or more forests with no basis | 1 | 1 / m |
| held, any | its kind's value | the same |
| a mixture's v | kappa, or 1 on kind none | kappa / c^2, or 1 / m^2 |

On kind none both reads of m are held: the start scale to 1e-100 through 1e100 and v's start to the squares of those, and the chain runs from there, as the base runs. Nothing else
reads them. A finite v alone does not do: a sweep squares the coefficient, so a start past about 1e154 overflows there. The bound is on the multiple m, `sd` over the unit, and not on `sd`: at a unit of 0.55, `sd` 1e99 and 1e-100 hold (worst difference 3.7e-12) and 1e100 and 1e-101 part (by 1.2 to 1.3 from sweep 1), and at m = 1e100 and 1e-100 exactly the chain holds (worst 1e-8). Where it holds the chain is the base's until a rounding flips a decision. [Calls made](#calls-made) 17.

Why this keeps the chain: with g = c on a basis forest and g = m with none, (a, F) on the base is (a / g, g F) here, the same state in
other units, so the first conditional is the same draw and so is each one after it, v scaling by 1 / g^2. Run on the port against
F's library, 19 calls at one seed, each kind and option under gaussian, probit and logistic: `train` and `sigma` within 3.1e-11, varcount identical, and the held call bitwise (runs/chain.R).

Set once. The chain's K-forest constructor copies each entry's `startScale` into a local vector and hands it to the combiner's
constructor, which multiplies and keeps nothing. A default empty vector leaves the neutral (1, 0, 1), so a fixture building a combiner
directly starts as today and 3 tests/cpp checks stay unpinned that the prototype moved. No later path reads the field: a write
re-resolves the entry and its start scale is ignored; a copy and a reload run the constructor at the `sd` in force and
[`restoreGlue`](../../src/bartcore/combiner.hpp) then overwrites the block and v; a sampler made again from its control is a
creation at the written `sd`, as it is today. No state, record or format moves.

A forest with no basis keeps its start by position under slice F's rule: [Calls made](#calls-made) 3.

### Every live path

Run on the port over ten kind-and-option cases: drawn under the normal and `student(3)` on none, levels and numeric, the half-Cauchy on
none, and held on each kind (runs/live.R). No case moved a leaf or a coefficient on a write.

| path | what moves, every kind and option | what does not | why it holds |
|---|---|---|---|
| `$setForest(f, sd = v)`, `$setLeafPrior(forests = )` | the trees' scale, from the next sweep; `$sd`, `k.scale`, the record | every leaf, coefficient, v | [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) re-resolves the entry and calls [`Chain::mapLeafScale`](../../src/bartcore/chain.hpp), and never reaches the combiner |
| the same before the first sweep | the same | the block keeps creation's start, as a one-forest write keeps its leaves | the same model, not creation's draws, on a drawn forest with no basis only: twin `identical()` on levels, numeric and every held case |
| a write of the `sd` in force | nothing | everything | the same multiple; draws `identical()` to no write |
| `$setLeafPrior(normal(sd = v))` on the one candidate | the same, by setmodel-parameters' slice 4 under dec-B540 | | no branch |
| `$setForestBasis(updateBasisScale = TRUE)`, numeric, sd unstated | the trees' scale, 1 / the new scale | coefficient, leaves, v | [`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp) unchanged: its comment's "no amplitude prior moves" is now true of every kind. Run: `k.scale` 0.5763 to 0.2882, reads otherwise unmoved |
| `$setForestBasis` otherwise | nothing | everything | no entry is re-resolved |
| a copy | rebuilt at the `sd` in force, then restores block and v | | a written unswept sampler's copy holds the live block, run |
| a reload | as a copy, from the stored state | | run; the record carries the written `sd` |
| a state stored, a write, a restore | block and v go back to the state; the prior and trees' scale stay the written ones | | a state carries no prior. Store, write, sweep, restore, write back returns the block `identical()`: run on all ten cases, the ratio the critic's probe read 10 on the form it replaced |
| a donor warm start | refused on several forests, as today | | test-multiforest-warmstart-refusal.R |
| `growFromRoot`, `sampleTreesFromPrior`, `sampleLeafParametersFromPrior` | run from the block's start, the trees at `sd` | | product and conditional as a sweep's |
| the live coefficient write of dec-B511, not built | below | | |
| a saved state from a development build | loads as it stands, the old units' same product; the next draws settle. No migration, `stateFormatVersion` does not move | | only shipped versions are owed compatibility |

The live coefficient write (coefficient-prior-live, dec-B511, TODO). Its plan does not exist yet; TODO's entry takes these words when
this lands, and its ~400 lines should fall since no leaf scale moves:

- A hold sets the block to its kind's value through [`holdAtKindValue`](../../src/bartcore/combiner.hpp) and leaves the leaves and
  the trees' scale. A release leaves the block at that value and the trees' scale, and draws from there: 1, or (0, 1) on levels. The
  entry re-resolves with `updateAmplitude` flipped, whose node factor is m both ways, so nothing crosses between the two.
- That item needs one combiner write, of `ForestAmplitudePrior` whole, which is what the probe's hook did, and a v: its TODO's "a t's
  mixing variance starts at the normal's" reads 1 / z^2 with no basis and kappa / c^2 with one, never the entry's 1 / m^2, which is
  creation's chain only. The probe started no-basis at 1 and read 1.75, 1.57, 1.13 gaussian and 4.24, 3.68, 1.64 probit at `sd` 7;
  its plan measures whether the difference shows.
- Costs it inherits, as ruled: a hold at a very large `sd` jumps, numeric gaussian `sd` 7 reading 10.49, 2.29, 1.20 (main 14.54,
  2.56, 1.16); a release with no basis at `sd` 7 is slower, as above against 2.02, 1.01, 0.98 on main. Holding a factor changes the
  model, level-0 rows losing their part, under any form.

### What goes

- `ForestCombiner::setAmplitudePriorScale`, its override and its call in
  [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp): nothing writes a mixture scale once no prior reads `sd`. The
  [`setForestMapSd`](../../src/bartcore/chain.hpp) branch's `return false` goes with it.
- `ForestMapEntry::nodeScaleDivisor`, its use in `ForestMapEntry::nodeScale`, and its assignment: always 1. x / 1 is x, so a held
  forest keeps its bits. 8 tests/cpp sites read it, among the re-pins.
- Comments stating the old split: the block above [`forestMapEntry`](../../src/bartcore/combiner.hpp), "Association and operand
  order are the K = 2 instance's", which no longer keeps a drawn instance bitwise, and
  [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp)'s.
- In the prototype, not built: `forceAmplitudePrior`, `holdForestAtKind`, `divideCoefficient`, `coefficientForAbsorb`,
  `setCoefficientAux`, `writeAmplitudePrior`, `kindStart_`, `startVariance_`, `scaleForestLeaves` and the two probe environment
  variables.

Stays: `aPriorScale` and `sdModerate` as a fixture's multiples, though the first no longer reaches a mixture ([Calls
made](#calls-made) 10); the half-Cauchy arm of [`drawAmplitudes`](../../src/bartcore/combiner.hpp); the `amplitude.prior.scale`
column.

### Engine, ~100 lines touched

1. [`ForestStatement::amplitudeDf`](../../src/bartcore/combiner.hpp): a +inf df states the normal coefficient, and on a forest with a
   basis means what no df means (3 lines of comment).
2. [`forestMapEntry`](../../src/bartcore/combiner.hpp), as the table above: `normalStated` kept apart from a finite df's `dfStated`,
   with the hold still on `isnan(df)`, so no local excludes the normal; one node factor; the two constants; `startScale`.
3. `ForestMapEntry::startScale`, and v's start, each held to 1e100 and its square (Calls made 17); `startScales` on [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor; the
   K-forest [`Chain`](../../src/bartcore/chain.hpp) constructor's vector.
4. What goes, above. No facade virtual, so [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move; `--preclean`
   regardless.

Not [`drawForestAmplitude`](../../src/bartcore/combiner.hpp), [`drawAmplitudes`](../../src/bartcore/combiner.hpp),
[`restoreGlue`](../../src/bartcore/combiner.hpp), [`glueIsValid`](../../src/bartcore/combiner.hpp),
[`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp), the k draw, the rescaling step. A new slot from a width change stays at
1, the multiplier's median. `stateFormatVersion` does not move.

### Bridge

No code. [`applyForestStatements`](../../src/R_interface_bartcore.cpp) passes a +inf through, and its comment says what it means. The
calibration reader's columns keep their names; `prior.scale` reads `sd`, and `amplitude.prior.scale` a constant of the kind and df.

### One spelling inside the engine

The statement is `ForestStatement::amplitudeDf = +inf`, not a new field. It is an internal code for normal: no t is drawn, no auxiliary
kept, no quantile taken, and the engine comment says so. R's `student(Inf)` stays refused. What reads "does this forest state a
coefficient prior" as `!isnan(amplitudeDf)` then needs no engine line: slice F's `forestSharers` and hold,
`ForestCalibration::amplitudeDf`, [`applyForestStatements`](../../src/R_interface_bartcore.cpp), the `amplitude.df` column,
control-migration's `dbartsForestRecord` slot, whose validity checks only length.

What R does need: the df is a three-valued code (NaN, finite, +inf), and one prior reads NA or Inf by kind. One decoder and one
predicate hold that in one place (R items 4 and 5). A third prior that is not a t would force a real field. [Calls made](#calls-made)
2 has the alternative.

### R, ~65 lines

1. `normal` in [`forestConstructors`](../../R/model.R), [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R),
   [`processHit`](../../R/formulaTerms.R)'s vocabulary and [`forceCallerCode`](../../R/family.R)'s hint. A caller's own variable or
   column named `normal` stays theirs in `basis =` and in a forest's predictors, as `student` does.
2. [`validateForestAmplitude`](../../R/model.R): a bare `normal` is the constructor; `normal()` kept as the object; `k` or `prior.sd`
   refused, by the words in [What the user sees](#what-the-user-sees); the refusal's list.
3. [`forestStatements`](../../R/model.R) gains the kinds and records +inf for `normal()` on a forest with no basis and NA on one with a
   basis. Its caller in [`resolveSamplerSpec`](../../R/spec.R) passes `hasBasis`. A default fit's record stays as it was.
4. One predicate, `statesCoefficientPrior(amplitude)`, true of a `student()` or `normal()` object, in place of every `is(x,
   "dbartsStudentPrior")` that asks whether a forest states its coefficient's prior: on F's tip, the lone branch and sharer refusal in
   model.R, the flat-prior refusal in spec.R and forestBasis.R's suggested spelling. The refusals' words take the spelling written.
   [`validateForestAmplitude`](../../R/model.R) and [`forestStatements`](../../R/model.R) keep the class test, which is about which
   class.
5. One decoder, `statedCoefficientPrior(df, held)`: `fixed()` where held, `normal()` at Inf, `studentPrior(df)` at a finite df, NULL
   where none is stated. [`reportLeafPrior`](../../R/dbarts.R) uses it, and forest-door's `$getForest` and its "value in force"
   comparison must, turning NULL into the kind's default. [`printedForest`](../../R/forestPrint.R): "; coefficient normal".
6. [`studentPriorRefusal`](../../R/model.R)'s words. Where forests-shared-total's slice 3 has landed, dec-B528's refusal with `sd`
   stated gains a clause: [Calls made](#calls-made) 5.
7. No reader changes code. `$getLeafPrior`'s `k.scale` reads the leaf scale and `sd` the map's `map.sd`, so both are right. The size
   reader stays TODO forest-size-reader.

## What the user sees

Ten taken-and-refused rows run on the port against `student(3)`, `normal()` matching it in every one (runs/taken.R), the rest on an
earlier prototype; the port's refusals still say `student()`. Step 4 re-runs them all on the tree.

| call | `normal()` |
|---|---|
| beside a factor; a numeric; a factor and a numeric, gaussian | taken, no sd: the leaf prior's size, as `student(df)` |
| the same with `sd = 0.7` | taken, trees at 0.7 |
| probit, logistic, no sd and with sd | taken; no sd: 1.574 and 2.85, slice F's values |
| a weighted fit | taken, exact arm with rows of weight zero |
| beside two sharers with `sd`, plain, and also beside a factor | taken, outside the share |
| with no `sd` on any forest of two or more unsized | refused, below |
| a lone forest, with or without `sd`; the multinomial entry | refused, below |
| a variance forest | refused, whole model, today's words |
| `bart()`, `dbarts()`, a `forests` list, `predict(bases = )` | taken |
| `xbart()`, `rbart_vi()` | refused, no `forest()` term, today's words |

The refusals, each as `student()`'s at the same place, tested by exact pattern:

- A lone forest and a multinomial entry: "this model has one forest, which has no coefficient to give a prior; a lone forest draws its
  size through the fitting function's leaf prior, as leaf.prior = normal(k = chi(1.5, 2)) or normal(sd = invchi(1.5, s)) with s a size
  in the response's units", one branch.
- A sharer: "forest 1: amplitude = normal() draws one forest's coefficient; these forests share one total: state each forest's 'sd' to
  give it its own".
- `normal(k = 2)`, `normal(sd = 2)`, `normal(k = chi(1.5, 2))`, on every kind: "'amplitude = normal()' takes no 'k' or 'sd': a drawn
  coefficient's size is the forest's own 'sd'". Only a bare `normal()`, as forest-door has it.
- An `amplitude` of another kind, the string `"normal"` among them: "a forest's 'amplitude' must be fixed(), which holds its
  coefficient, normal() or student(df), which draw it under that prior, or left out; ...".
- `student(Inf)`: "a coefficient's student(df) takes one finite df of 1 or more; student() is student(1), the half-Cauchy, and normal()
  is the normal".
- A flat prior on k with no `sd`: F's words, "amplitude = normal() with no 'sd' takes ...".
- A live change of `amplitude`: today's "$setLeafPrior restates only a forest's 'sd': 'amplitude' is fixed at creation", and at
  `$setForest` forest-door's words for every changed `amplitude`, until coefficient-prior-live.

### Readers

Run on the port at `sd` 0.7, a forest's unit 0.576: every `k.scale` read 0.7, and a drawn block with no basis started at 0.823, the unit over `sd`.

| reader | a basis forest, drawn | no basis, drawn | held |
|---|---|---|---|
| `$getLeafPrior(f)$k.scale`, `fit$leaf.prior`'s | `sd`; was `sd` / 0.674 | `sd`; was the unit | `sd`, as today |
| `$getLeafPrior(f)$sd`, `$sd.stated` | unchanged | unchanged | unchanged |
| `$getLeafPrior(f)$prior.sd.of` | "... (at a coefficient of 1)"; was "... (prior median)" | "forest (at a coefficient of 1)"; was "forest (prior median)" | unchanged |
| `$getLeafPrior(f)$amplitude` | unchanged, NULL at the default and under a written `normal()` | `student(df)` or `normal()` | `fixed()` |
| `getK(f)` | 1 | 1 | as today |
| the map's `amplitude.prior.scale`, read by tests only | NaN under the normal; z / (t c) sqrt(kappa) under a t; was without c | 1 at the half-Cauchy, 1 / t under a t, NaN under the normal; was m / t, and no longer moves with `sd` | NaN |
| the map's `amplitude.prior.variance` | kappa / c^2; was kappa | 1 / z^2 under the normal; NaN under a mixture | kappa |
| `$getForestAmplitudes(f)`, `fit$glue`, raw extract's coefficients | 1.48 times today's; prior median size 1.0007 | today's over m; prior median size 1 | its kind's value |
| raw `$getForestFits`, `$predictForests`, `extract(fit, "forest")`, `predict(type = "forest")`, `fit$forestFits` | today's times 0.674 | today's times m | as today |
| `extract(contribution = TRUE)`, the fit, `predict(type = "ppd")` | unchanged until a rounding flips a decision | unchanged until a rounding flips a decision; new under the normal | unchanged |
| `print(sampler)`, `print(fit)` | unchanged: `$sd` | "1 forest1: no basis; sd 0.7" | as today |
| `summary(fit)`, the verbose creation summary | unchanged | "...sd 0.7 (stated); coefficient normal" | as today |

- `extract(type = "leaf.prior.sd")` on several forests reads `$sd`, unchanged, and `k.scale` only on a sharer under a drawn k, whose
  entry does not move. No other R reads `k.scale` on a forest of several.
- `forest(sd = 2, amplitude = normal())` prints `amplitude = normal()`. `fit$fixed$amplitude`: no entry, as for anything drawn.
- `$getForest(f)`: `amplitude = normal()`, through this slice's decoder, forest-door's "value in force" comparison likewise.
- bartCause, no code. Its `bcf()` states `amplitude = student(1), sd = sd.control` and `sd = sd.moderate`, and reads only `response.scale`
  and `response.shift` from `$getLeafPrior`. `icate` is `response.scale * forests$tau * bcfGlueDifference(glue)` and the counterfactuals
  use (b.1 - b.0) times tau: products. Run on the ruling's prototype: effects match draw for draw on gaussian to 1e-10; on a binary
  response, which forks through rounding, over 72 seeds, largest |z| 2.26 against 4.39 for main compared with itself. Its test-14-bcf.R,
  which reconstructs through those products, passes 184 of 184 on the port. What moves: `fit$glue[, , "b.0"]`, `"b.1"` 1.48 times larger and `fit$forests$tau`
  smaller by that factor; `fit$glue[, , "a"]` over m, 0.5 at the default gaussian `sd.control`, and `fit$forests$mu` times m. On its own
  branch when this lands: man/bcf.Rd's `sd.control`, "A drawn amplitude takes a half-Cauchy prior whose median is that size", becomes
  "takes a half-Cauchy prior with median 1, the forest's trees at that size", and `glue` and `forests` say a coefficient has prior median 1.
- stan4bart: not touched. Its mvbart makes one lone forest per equation, with no basis or `amplitude`, and reads only `response.scale`
  and `response.shift` from `$getLeafPrior()`. treatSens and bairrtt: not touched. No shipped header changes.

### Help

- [`forest`](../../man/forest.Rd), the main change:
  - `sd`: "Where it is drawn the size is itself uncertain, and `sd` is its prior median" becomes "Where it is drawn it is a multiplier
    with prior median 1: the forest's size is the coefficient times `sd`, uncertain, and a coefficient of 1 means the size stated".
    "The size of what this forest contributes" and the kinds above it stand.
  - `amplitude`: `normal()` beside `student(df)`. "With a basis the t is scaled so that its median size is the normal's, so the forest's
    `sd` keeps its meaning, the prior median of what the forest contributes, and only the tails move" becomes "On every kind a drawn
    coefficient's prior median size is 1, under `normal()` and `student(df)`, so only the tails move". "With no basis the coefficient's
    size is the forest's size, its median the `sd` or, where none is stated, the size ..." becomes "With no basis `sd` is the size of
    the trees or, where none is stated, ...". "Bayesian causal forests' prognostic prior is written `amplitude = student(1), sd = 2 *
    sd(y)`" stands: it is the same prior on the product.
  - One sentence: at one `sd`, a drawn normal's contribution has about 1.48 times the standard deviation of a `fixed()` forest's, since
    the coefficient and the trees are both random, and thinner tails than `student(df)`'s.
  - The size section's 0.29, 1.049 sqrt(K) and 0.699 read products, and hold: reread, not rewritten.
  - `prior.sd.of`'s "(at a coefficient of 1)" where the reader's list says it, once.
  - The example: `indexSd <- sqrt(calibration$k.scale^2 * 0.5 * rowSums(basis^2))` becomes `sqrt((calibration$k.scale / 0.674)^2 * 0.5 *
    rowSums(basis^2))`, its comment "its coefficients' prior variance is 1/2 over 0.674^2". As it stands it gives 1.06 at `sd` 1.5 for a
    true 1.57, with no error: `index-sd` below reads this line from the installed help and is the pin.
  - A caller's own function named `normal` is not reached inside `forest()`.
- [`dbartsSampler-class`](../../man/dbartsSampler-class.Rd) and the docstring in [`R/dbarts.R`](../../R/dbarts.R), for codoc:
  - `getLeafPrior`: "`k.scale` is the map's leaf scale" gets "the forest's `sd`, on every forest but a sharer"; "A drawn forest with no
    basis carries its `sd` as the median of its coefficient's prior, which does not move with the mixture's own variance auxiliary, which
    is drawn; a level's coefficient is drawn under prior variance 1/2 and a numeric column's under 1" becomes "A drawn coefficient is a
    multiplier with prior median 1, a level's difference and a numeric column's through bcf's 0.674".
  - `forests` argument: "Each `sd` lands in the channel creation gave it - the median of a drawn forest's coefficient where it was
    created without a basis, the forest's own leaf scale otherwise" becomes "Each `sd` sets the forest's leaf scale".
  - `setLeafPrior`'s entry: "restating each named forest's spread ... in the channel its creation gave it" loses the channel, and gains
    "A write changes the model." Those words and no more: the maintainer, "I wouldn't over-emphasize that, you can just say that it
    changes the model, if anything." forest-door's `$setForest` entry says the same, once.
  - `getForestAmplitudes`: a drawn coefficient has prior median 1.
- [`dbartsPriors`](../../man/dbartsPriors.Rd)'s `normal` says the coefficient use in two lines. [`bartBT`](../../man/bartBT.Rd)'s `glue`
  and `forestFits`: the coefficient has prior median 1 and the raw total prior scale `sd`.
- NEWS: one clause in 1.0-0's forest item. Nothing released moves.

## Slices

One slice: engine, R, help, tests, arms, design note, then one review and the gate list. Steps, each verifiable:

1. Rebase on slice F's tip, 7f53d1c78; the orchestrator rebases onto main at landing. Read F as built, its hold and no-sd median each keyed on `isnan(df)`. Grep R as built for
   `dbartsStudentPrior` and `amplitude.df`: each site goes to the predicate or the decoder. Re-count the failures below before
   re-pinning anything.
2. Engine items, then the chain A/B and `bcf-equivalence.R compare --cross-host` before any test is touched: tier 1 on 19 of 19.
3. tests/cpp, re-pinned and new. Run plain and under ASan and UBSan.
4. R items, `--preclean`, tinytest re-pins and file. Run the A/B against the base.
5. The exact arms, quick and full; the other gates' reads restated; sbc's bcf arms against the base's verdicts; the fork read.
6. Re-record bcf-equivalence and equivalence, each with a MANIFEST row.
7. Help, NEWS, design notes, lint, air, codoc, freshness.
8. Review with mutants, fix round, gates, records, CI on the branch, merge.

## Tests and gates

### What the failures pin

Every failure read, on two runs. On the ruling's prototype over main's engine: 121 of 32,380 tinytest results and 123 tests/cpp checks.
On the port to slice F at 7f53d1c78, the base this builds on: 88 of 32,947 and 121, the same 41 sites (at home, on eb1ed0a2b, 112 of 32,994). A rebase onto a moved F re-counts. The
port hands the start to the combiner's constructor, so 3 prototype failures on a directly built combiner do not arise.

| class | what it pins | tinytest, main / F | tests/cpp, main / F |
|---|---|---|---|
| units | a reader or an entry in the old channel: `k.scale` x 0.674, `amplitude.prior.scale` x unit, a literal variance, a start | 116 / 87 | 115 / 114 |
| a write before the first sweep | a sampler written before its first sweep draws as one created there, on a drawn forest with no basis; dec-B550 gives this up | 4 / 0 | 5 / 3 |
| dec-B539's | `normal()` or +inf was refused, and is taken | 1 / 1 | 3 / 4 |
| a regression | | 0 / 0 | 0 / 0 |

No failure is a regression. Each was read to the check and its value.

**tinytest on F, 88 results at 41 sites (24 more on one site's at-home arm).**

| file | results, sites | pins | class |
|---|---|---|---|
| test-bcf-family.R | 25, 3 | `k.scale` x 0.674 is `sd`; `priorScales` against c(s, s) / 0.674; the index sd from hard-coded variances | units |
| test-forest-sd-unit.R | 23, 18 | `expectStated`'s `inForce` reads `k.scale` x 0.674 at every door: term, list, spec, data object, `bart()`, rows in the likelihood, `setWeights`, probit, logistic, writer, copy, reload, `new()`, `setState`, `setResponse`, other data; one helper | units |
| test-forest-leaf-prior-reader.R | 15, 1 | `expectEngine`'s `k.scale` x 0.674 is `$sd`, in a loop | units |
| test-forest-basis-r5.R | 9, 9 | `priorScale(.)` at s / 0.674 on a swap, an install, a restore | units |
| test-forest-scale-record.R | 4, 3 | `expectRule` reads the column rule through the leaf scale | units |
| test-forest-units-twin.R | 3, 1 (27, 2 at home) | `sizeCheck` reads the trees' spread x z, which is now `sd` exactly | units |
| test-forest-lone-no-basis.R | 3, 1 | a drawn forest's `k.scale` is the unit, F's own pin | units |
| test-forest-amplitude-student.R | 2, 1 | `coefficientSizes`: the coefficient's median 0.674 at a numeric default; 1.0007 after | units |
| test-bcf-creation.R | 1, 1 | `k.scale` is `$sd` / 0.674 | units |
| test-calibration-midchain.R | 1, 1 | the same | units |
| test-prior-init-composed-distribution.R | 1, 1 | the start c(1, 0, 1) on the K = 2 sampler | units, a start |
| test-forest-arguments.R | 1, 1 | its refusal list has `normal()`, taken now | dec-B539's, rewritten with R item 2 |

**tinytest on main's engine, 121 results, on the ruling's prototype.**

| file | results | class |
|---|---|---|
| test-forest-sd-unit.R | 44, and `namedWrite`'s 1 | units 44; write 1, a draw against a twin created at the written `sd` |
| test-bcf-family.R | 27 | units |
| test-forest-leaf-prior-reader.R | 20 | units |
| test-forest-basis-r5.R | 10 | units |
| test-forest-scale-record.R | 4 | units |
| test-multiforest-leaf-prior-writer.R | 1, and `twinCases`' 3 | units 1, a write on forest 1 leaving its `k.scale`; write 3, a gaussian scale-mixture median, probit with both at once, K = 3 |
| test-forest-units-twin.R | 3 | units, `sizeCheck` |
| test-forest-amplitude-student.R, test-calibration-midchain.R, test-bcf-creation.R | 2 each | units |
| test-prior-init-composed-distribution.R | 1 | units, the start c(1, 0, 1) |
| test-forest-arguments.R | 1 | dec-B539's |

On F the 4 write results and the unit in test-multiforest-leaf-prior-writer.R stop failing, each a forest with no basis stating `sd` and no prior, which F holds. The state-before-write pin in
that file passes on both: a state holds no prior and the block is not moved.

**tests/cpp on F, 121 checks in 44 messages.**

| messages, by group | checks | pins | class |
|---|---|---|---|
| `map statements`: its leaf scale is the expression's bits; a factor takes one unit over 0.674 | 36 | the entry's node factor and divisor | units |
| `calibration map`: the reported scale, fixed-variance report, scale-mixture report, a zero column | 18 | the reader's product | units |
| `sd unit`: nine messages | 24 | quotient on the anchor; carried in the leaf scale or the median; a write of the `sd` in force keeping its multiple | units |
| `numeric rule`: factor beside it; one unit per sd; prior precision 1; after a write; after a swap; standing second | 13 | quotient, precision 1, start 1 | units, starts |
| `student`: scale at df 1 and 3; reports; v starts; median matched; no df is the kind's entry; donors' v; no basis keeps the unit, divides by t; df 1 draws | 14 | the t's scale without c, v's start kappa | units |
| `held coefficients: the drawn numeric forest beside them`, `and the forests' defaults`, `median: at a held k` 2, `shared total: each drawn coefficient starts at 1`, `bcf shape`, `drawGlue`, `a pinned b block`, `the q-variate sweep` | 9 | a leaf scale over 0.674, a mixture's m, or a start read in the old unit | units, starts |
| `sd unit: draw for draw` | 3 | a write before the first sweep | write |
| `student: a df below 1 or not finite throws` 3; `hold: an infinite df is stated, and is not held` 1 | 4 | the bad list has +inf, the normal's code now | dec-B539's |

**tests/cpp on main's engine, 123 checks, on the ruling's prototype.**

| group | checks | class |
|---|---|---|
| `map statements`: the two above, 18 each, and its two no-basis forms, 3 each | 42 | units |
| `calibration map` | 18 | units |
| `sd unit`, nine messages | 20 | units |
| `numeric rule` | 13 | units, starts |
| `student` | 14 | units |
| `held coefficients`, `the forests' defaults`, `relative row bound`, `a unit treatment multiplier`, `a constructed BCF combiner holds the neutral (1, 0, 1) glue`, `drawGlue`, `a pinned b block`, `the q-variate sweep`, 1 each | 8 | units, starts |
| `sd unit: draw for draw` 3; `map writers: a written twin draws bitwise as the constructed one` 2 | 5 | write |
| `student: a df below 1 or not finite throws` | 3 | dec-B539's |

Between the two runs: F's held forests retire the 2 `map writers` checks; a start handed to the constructor, neutral by default, keeps `relative row bound`, `a unit treatment multiplier` and the
neutral glue passing on a fixture; F's rules add the `shared total` start, the `median` pair, and its `hold` pin, which flips.

### What is given up, and what replaces its pin

A write before the first sweep on a drawn forest with no basis is the same model and not creation's draws, settling by about sweep 80 rather than 20 after a tenfold write
(dec-B550's run): the block starts at the unit over the `sd` in force at creation. On every other forest the start does not read `sd` and the twin holds, `identical()` (run). The 3
`sd unit: draw for draw` checks re-pin on a forest whose start does not read `sd`, and the new `write` checks pin that this block stays. No twin and no bcf-equivalence scenario then reads
this path, so what replaces the twin as evidence of the distribution is the exact arms that follow a write ([The exact gate](#the-exact-gate)).

### tests/cpp, new

`testForestAmplitudeNormal` beside [`testForestAmplitudeStudent`](../../tests/cpp/test_sampler.cpp), ~150 lines:

- `map: every row`: each row of [forestMapEntry](#forestmapentry), to the last place, by a sample-free read, kind by kind, at no sd and at
  a stated sd through the unit; no field of a drawn entry's prior moves with `sd`, only v's start and `startScale`.
- `normal: a basis keeps its own entry`: +inf and NaN `sameEntry` on levels and numeric, 0.674 and not z. `normal: a held forest ignores
  it`. `normal: what is refused`: -inf, 0, 0.5 and -1 throw.
- `start: each block`: every row of [The start](#the-start), with and without a mixture, standing first, second and third; through the
  chain, and neutral on a combiner built directly.
- `start: the chain stays`: entries built by hand on main's map (node m / c, variance kappa, kind's start; node 1, scale m, v 1) against
  the map's: the fit's draws within 1e-9, on none, a factor, a numeric and each under `student(3)`.
- `write: nothing drawn moves`, normal and student, on every kind, before the first sweep, after a prior draw of the trees and after a
  sweep: block, v and leaves `==`, the leaf scale the new one's bits, the prior's entries unmoved.
- `write: the sd in force`: nothing moves. `write: a sampler stated at creation and written back`: the same bits.
- `write: restore commutes`: store, write, sweep, restore, write back is the stored block and v.
- `write: restate`: `restateForestBasisScale` writes the leaf scale, block unmoved, a sampler made at that scale drawing equal on a forest
  whose start does not read it.
- `normal: the same chain as a column of ones`: kind none under +inf at multiple m against kind numeric on a ones basis at m c / z, compared
  at m = z, the one m where the two starts agree (1 / z against 1 / c, with the trees at z and c), under gaussian, probit and logistic: same
  counts, draws equal to rounding. Run: 7.8e-15 gaussian, 1.2e-13 probit, 1.9e-13 logistic. It is `normal()`'s only oracle under a link. At m z / c it differs
  by 5.5e-4 and at another m it forks. It survives dec-B528, which refuses at R.
- `extreme sd`: `student(3)`, the half-Cauchy and `normal()` on a forest with no basis at multiples 1e-160 and 1e300: v's start and the start scale finite and
  positive; a sweep's fit and sigma finite.
- `normal: no sharer`: `forestSharers` excludes it, with sd and without. `normal: the no-sd size`: equal to df 3's multiple at a held k, a
  drawn k and a named sd.
- `normal: the calibration reader`: `amplitudePriorVariance` the entry's, `amplitudePriorScale` NaN, `amplitudeDf` +inf.

### tinytest, new

A new test-forest-amplitude-normal.R, small n, ~260 lines:

- `sd-carried`: the sibling sweep dec-B550 rests on. On one data set, doubling `sd` under `fixed()`, `normal()` and `student(df)`, on every
  kind where each is taken: `k.scale` and `$sd` doubling, and each equal to the other; the coefficient's prior columns not moving. So no
  exception can come back.
- `live`: `$setLeafPrior(forests = )` in mid-chain on each kind and option: `getForestFits` and `getForestAmplitudes` `identical()`, v
  through the state `identical()`; before the first sweep and after `sampleTreesFromPrior`, the block `identical()`; through
  `$setForestBasis(updateBasisScale = TRUE)` likewise; a copy and a reload continue as the live sampler does. runs/live.R is the start.
- `commute`: store, write, sweep, restore, write back on each kind and option: block `identical()`.
- `index-sd`: forest.Rd's `indexSd` line as installed, read through `tools::Rd_db` and evaluated, against the product's spread at prior, by simulation, within 2
  percent. So the Rd is what is read and not a copy in the test: nothing else reads the example's number.
- `extreme-sd`: a drawn `student()` and `normal()` forest with no basis at `sd` 1e-160 and 1e300: no refusal, `train` and `sigma` finite.
- `normal-median`: `$sd` and `$prior.sd.of` equal to `student(3)`'s, the latter "forest (at a coefficient of 1)", `$k.scale` the `$sd`, no sd and sd stated, under gaussian, probit,
  logistic.
- `normal-reader`: `$amplitude` the class; print and summary lines; `forest()`'s print; `getForestAmplitudes` moves; no
  `fit$fixed$amplitude`.
- `normal-record`: `amplitude.df` Inf on the forest stating it, NA beside it, absent on a fit stating none, on a basis forest's `normal()`
  too.
- `normal-doors`: formula, `forests` list, a data object's `bases`, `bart()`, a bare `normal`, draws `identical()` across the spellings.
- `normal-siblings`: beside one basis, beside two; weights; at a stated `sd` beside sharers, sharers' parts and shares unmoved.
- `normal-refusals`: each pattern in [What the user sees](#what-the-user-sees), counted with `withCallingHandlers`; a variance forest; the
  multinomial entry; `normal(sd = 2)` and `normal(k = 2)` on a basis forest too.
- `normal-basis`: on a factor and a numeric forest, draws and record `identical()` to left out.
- `normal-shadow`: with a caller's `normal <- function(v) 2 * v`, `forest(x1, sd = normal(.5))` refused, "forest 'sd' must be a number, not
  a dbartsNormalPrior", in a formula and in `forests =`; with a caller's `normal <- function() dbartsPriors$fixed()`, `amplitude = normal()`
  is the package's and draws; a caller's variable or column named `normal` stays theirs in `basis =` and the predictors.
- `normal-column`, once forests-shared-total's slice 3 has landed: `y ~ a + forest(b, basis = one, sd = 2)` refused with the clause; its
  advice as written meets the plain-terms refusal, which names `forest(a) + forest(b, sd = 2, amplitude = normal())`, and that is taken.
- `normal-live`: `$setLeafPrior(forests = list(forest(amplitude = normal()), forest()))` refused as today.

### Mutants

39, each with a name in `mutation-battery.R`'s form where one fits, else a scratch patch.

| name | change | killed by |
|---|---|---|
| start-basis | a basis block starts at its kind's value, no 1 / c | `start: the chain stays`; tier 1 |
| start-none | kind none starts at today's value, no 1 / m | `start: each block`; tier 1 |
| start-held | the start scale reaches a held block | `start: each block`; the trio bitwise |
| start-fixture | a combiner built directly takes a scale | `start: each block`, neutral |
| v-basis | v starts at kappa, not kappa / c^2 | `start: each block`; tier 1 |
| v-none | v starts at 1, not 1 / m^2 | `start: each block`; tier 1 |
| mixture-basis | no c in a basis forest's mixture scale | `map: every row`; arms under `student(df)` on a basis |
| variance-basis | no c^2 in a basis forest's variance | `map: every row`; every arm of a drawn basis |
| node-divisor | a basis arm keeps a node of m / c | `map: every row`; `sd-carried`; every arm |
| node-none | kind none's node stays 1 | `map: every row`; `sd-carried` |
| prior-reads-sd | kind none's mixture keeps m | `map: every row`; `sd-carried`; arms at 0.7 and 3 of the unit |
| normal-reads-sd | kind none's normal variance (m / z)^2 | `map: every row`; arms at 0.7 and 3 of the unit |
| basis-z | z on a basis | `map: every row`; `start: the chain stays` |
| none-674 | 0.674 on kind none's normal | `map: every row`; the arms cannot see 0.07 percent, and say so |
| write-start | a write moves an unswept block | `write: nothing drawn moves`; `live` |
| write-v | a write resets v | `write: nothing drawn moves`; `commute` |
| write-leaves | a write scales the leaves and divides the block, the form first shown | `write: nothing drawn moves`; `commute` |
| write-none-scale | a write on kind none leaves the leaf scale | `write: nothing drawn moves`; `sd-carried` |
| mixes | a t path at +inf | `map: every row`; the suite on NaN |
| held | drawn switched off on +inf, slice F's hold keying on a local | `normal-reader`; prior arm |
| not-recorded | Inf not recorded on a forest with no basis, which holds silently | `normal-record`; `normal-reader` |
| basis-inf | +inf on a basis reaching the t arm | `normal: a basis keeps` |
| sharer | `forestSharers` taking it | `normal: no sharer` |
| r-sharer | the sharer refusal without normal | `normal-siblings`, its print |
| no-sd-size | a no-sd normal not taking the leaf prior's size | `normal: the no-sd size`; `normal-median` |
| reader-student | the decoder builds `studentPrior(Inf)` | `normal-reader` |
| door-reader | `$getForest` or the "value in force" comparison reads Inf as a t: for the door's writer, whichever lands second | forest-door's no-basis round trip |
| predicate | a site left testing `dbartsStudentPrior` alone | `normal-refusals`; `normal-siblings` |
| shadow | `normal` left out of a vocabulary list | `normal-doors`; `normal-shadow` |
| lone-taken | the lone branch student only | `normal-refusals` |
| basis-recorded | Inf on a basis forest's `normal()` | `normal-record`; `normal-basis` |
| spare-args | `normal(sd = 2)` taken | `normal-refusals` |
| example | forest.Rd's example keeps `k.scale` undivided, in the Rd and not the test | `index-sd`, which reads the installed help |
| dead-write | `setAmplitudePriorScale` kept, writing a mixture scale that reads `sd` | `prior-reads-sd`'s killers; a grep at review |
| restart-after-restore | the start scale re-applied after a restore, in the copy or reload path | `write: restore commutes`; `live`'s copy and reload |
| v-not-restored | a mixture row's v not restored | `write: restore commutes`; `commute`, through the state |
| inforce-resolves | a write of the `sd` in force re-resolving the entry | `write: the sd in force`, by the bits |
| kappa-dropped | a factor's kappa dropped from its variance or mixture scale | `map: every row`; the factor arms of the exact gate; bcf-exact |
| start-overflow | v's start or the start scale not kept finite | `extreme sd`; `extreme-sd` |

### The exact gate

[`forest-prior-exact.R`](../../benchmarks/R/forest-prior-exact.R) restates every drawn arm and gains the new ones. It has every kind already:
[`numericPosterior`](../../benchmarks/R/forest-prior-exact.R) and [`factorPosterior`](../../benchmarks/R/forest-prior-exact.R) integrate a drawn
coefficient on a normal or t grid, and [`priorFit`](../../benchmarks/R/forest-prior-exact.R) reads a size against its quantiles. The oracle is
written from the rule and nothing else: every forest's trees have sd `size`, held or drawn; a coefficient's scale is 1 / 0.674 with a basis and
1 / z with none (`spread`), a factor's sqrt(1/2) of that, and a t's scale is that spread's z / t_nu,.75. A prior read of a kind-none coefficient
has median 1, and a basis's z / 0.674. scratch/fan-plan3/fpe-f.R is that gate on slice F, 138 lines against F's.

| arms | how many | reads |
|---|---|---|
| held on each kind, alone and beside others, and sharers: unchanged | 18 | terms or fits, splits |
| a 0/1 column and four values, drawn, at the default and stated; rows of weight zero; three forests of three trees: restated | 6 | fits, splits, within 0.5, 1, 2 |
| four values under `student(1)` and `student(3)`, a factor at the default and under `student(3)`: restated | 4 | the same |
| four values and a factor under `normal()` written out, default and at 3 of the unit; four values under `student(3)` at 3 and `student(1)` at 0.6 of the unit; a 0/1 column under `student(3)`; a factor under `student(3)` at 0.6 of the unit and under `student(1)`: new | 9 | the same |
| no basis under each of `normal()`, `student(1)`, `student(3)`: default, 0.7 and 3 of the unit, beside a numeric held, three forests of three trees, rows of weight zero: new, or restated where `student(df)` had one | 18 | the same |
| prior: a number and a factor at the default and under `student(1)` and `student(3)`: restated | 3 | the size within the prior's quantiles |
| prior: no basis under `student(1)` at its default, `student(3)` at 0.7 of the unit, and `normal()` at the default, 0.7 and 3 of the unit | 5 | the same |
| following a write, no basis: `normal()` 0.7 to 3 of the unit before the sweep and in mid-chain; `student(3)` 3 to 0.7 before the sweep; `student(1)` 0.7 to 3 in mid-chain: new, the oracle the written `sd`'s, written from the rule | 4 | fits, splits, within 0.5, 1, 2 |

Run on the port, quick, 32 seeds of 25,000 kept: all 63 arms pass, largest |z| 3.45, one `within2`; the 4 following a write pass, largest |z| 1.85 (scratch/fan-crit3/runs/fpe-write.R). Full is owed at step 5.

The other gates on exact-gates.yaml's list that read a forest, quick on the port with no edit:

| gate | on the port, quick | owes |
|---|---|---|
| bcf-exact (also run full), bcf-exact-restricted, bcf-latent-exact, backfit-exact, blocks-exact, probit-k-scale-exact, forests-blocks-equivalence | pass, unedited: each reads a fit or a product. bcf-exact in full prints the base's numbers to four places in every mode | nothing |
| bcf-exact-weak | fails, `P(\|a\|<=t)` 0.2197 0.4423 0.7343 against exact 0.4335 0.7381 0.9277 | its read goes in units: `abs(aVec) <= q / aPriorScale`, one line. Run so: 0.4423 0.7343 0.9190, gap 0.0088, passes |

### sbc's bcf arms

[`sbcBCFGlueDraw`](../../benchmarks/R/sbc.R) is a generator, not a read: it draws a ~ Cauchy(sd.control) and b0, b1 ~ N(0, 1/2), installs them, and draws the trees at the prior, so
under this form its coefficient is in the old units and its trees in the new, and every bcf arm's prior is wrong. Its self-check still reads an IQR of 4, since that reads what the generator says.

- The arms, none in sbc.yaml's matrix: `bcf`, `bcf-weak`, `bcf-probit`, `bcf-logistic`, `bcf-probit-weak`, `bcf-logistic-weak`. The ladder: `burn-bcf-probit`, `burn-bcf-logistic`, and their prescribed |a| strata, which run at
  an `sd.control` of 1 and so stand. It runs on one arm of each link, as a cost check, and does not gate.
- The generator, restated from the rule: a ~ Cauchy(1) and b0, b1 ~ N(0, 0.5 / 0.674^2), the trees at `sd` as the engine draws them. The moment check [`sbcCheckBCFGlue`](../../benchmarks/R/sbc.R) reads those. The `glue-sd` poison's a-scale is
  its 2 over `sd.control`, so it stays a mismatch of the same size. The fixed-glue path's (1, 0, 1) stands.
- The four read sites, [`sbcBCFSample`](../../benchmarks/R/sbc.R), [`sbcMakeBCF`](../../benchmarks/R/sbc.R)'s map, `sbcCheckBCFLatent` and the family spec's `sample`, read the product or a rank in matching units, so they stand.
- The pass rule, since CI skips these arms and the base already flags some: the same arm at the same R, L, thin, burn and seed on the base's library and the tree's. Each functional's verdict on the tree is the base's verdict, or the difference is explained
  by the fork read. A functional that the tree flags and the base does not fails, and stops the build. A run of the old generator against the new engine is the discrimination check: a small bcf-weak run read abs.a's ecdf gap 0.31 against the base's 0.18.
- Sizes: chosen so both libraries' runs take about 2 hours in all; the sizes and numbers go in this plan's Landing note. Agent-hours 2.

### Forks

A fork is a decision that rounding flips, after which the two chains are two runs of one transition kernel, so a forked chain has nothing to bias it. That is the argument. What was run, on the critique's port at 7f53d1c78 against the base's library:

- No exact arm reads a fork. bcf-latent-exact in quick, both links, a drawn coefficient: its chains captured on both libraries, 62,500 to 150,000 sweeps a seed over 9 modes, no chain past 1e-6, every z the same to the digits printed. Its forests are one tree over two cells, with too few decisions to flip. So the exact arms pass on chains that did not fork, and say nothing of one.
- Forked chains against each other: the half-Cauchy forest with no basis under probit, the call that forks most, 200 rows, 20 trees, 4,000 sweeps. On 24 seeds every chain forks, at sweeps 666 to 2,123, median 1,426. Over 288 further seeds on each library, reading sweeps 2,500 to 4,000: the grand mean of the fit differs by -0.0002, standard error 0.0005, z -0.39; per row, median |t| 0.66 against 0.68, 2.0 percent past the 97.5 percent point, largest |t| 2.62 against a base-against-base split's typical 2.86. After a fork the rms difference of a seed's mean is 0.037, against 0.042 between two seeds of the base. (Smaller sets read further out, 24 seeds a grand-mean z of 2.8 and 72 seeds z 0.19 with 10.5 percent of rows past; the 288 do not repeat either.)
- So by this reading a fork is what a new seed is. No gate catches a bias in one: the arms that would fail on a wrong kernel are the exact ones, which pass.

Probit and logistic are not gated by a new arm, as `student(df)` is not (dec-A245): the update is the one conditional, the latent refresh reads
the combined fit, and the prior enters as one precision. [Calls made](#calls-made) 7. The chain's agreement with F's on probit and logistic is
tier 1's.

### Baselines

| baseline | status | run on the port |
|---|---|---|
| bcf-equivalence-e346b018, F's | re-recorded, with a MANIFEST row | `compare --cross-host`: tier 1 on 19 of 19, worst ratio 4.2e-4, varcount identical. Against 20eed7ed the same, the check the ruling ran on its prototype. `numeric_held` and any scenario with no drawn forest bitwise |
| equivalence-e346b018, F's | re-recorded with a MANIFEST row, for one cell, `bart2twoforest`, whose factor forest is drawn | `compare --bitwise --strict-coverage` on a reference build: 54 of 55 bitwise; `bart2twoforest` not, 209 summaries, largest \|z\| 1.48 |
| multinomial-80b1c8d4 | no re-record | `compare --bitwise` on a reference build: 11 of 11 |
| the four snapshot files | unchanged | run whole on a reference build: 7, 3, 14 and 3 results, none failing |
| bench-sampler | not re-recorded, no compare owed | the sweeps are one conditional's: the same up to rounding |

Tier 1 stands in for a bitwise compare in bcf-equivalence, where only rounding moves: `sigma` and `train` within tier 1 and varcount identical.
Not z mode: at 40 draws it fails an unchanged sampler, main with engine seed + 1 failing 17 of 19. The snapshot channels `mu`, `tau` and `glue` are
compared by neither, being in new units. `numeric_swap_write_copy` writes `sd` on a running sampler: under this form that write is the leaf scale
as on main, so it passes with the rest, 3.2e-7 on the ruling's prototype, where dec-B542's form could not.

### The transient on a write

No gate: a pin says what a write moves, and no pin says how long the fit takes. The ruling's numbers, one data set, 8 seeds, rms over settled at
sweeps 1-5, 6-20, 21-80: a forest with no basis, 0.7 to 7, 1.32 1.26 1.03 gaussian and 1.76 1.83 1.03 probit, main near 1.0; a numeric forest
under probit, 0.7 to 7, 1.07 1.19 1.09, as main's. Burn-in at a large `sd` with no basis: probit, `normal()`, `sd` 6, 0.88 0.90 0.98. A larger
sampler's `sd` mixes slowly where it draws `sd` from the leaves, 34 to 41 per 5,000, the fit mixing the same. The help says only that a write changes
the model.

### The rest of the gate list

By the tier: an A/B against the base's library over the channels coefficient-prior-opt-in's A/B read, on the shipped and reference builds, none
differing past rounding and only with a drawn forest, but for a write of `sd` on a drawn forest with no basis, before the first sweep or in mid-chain, which differs by design and has the exact arms; the 19 calls of runs/chain.R; multinomial bitwise; the snapshot files unchanged; tests/cpp plain
and under ASan and UBSan; the R-loaded ASan run; exact-gates.yaml's list quick; [sbc's bcf arms](#sbcs-bcf-arms) against the base's verdicts; the full suite once in one process at default
ncpu; `R CMD check --as-cran`; `lintr::lint_package()` against the slice's library, air, check-rc-codoc, check-win-drift, check-doc-freshness and
`mutation-battery.R verify-anchors`; the x86 leak-sanitizer leg before landing; bartCause's bcf tests against this library.

### What stays bitwise

Every call with no drawn forest: held forests of every kind, a sharer, a lone forest, multinomial 11 of 11, the four snapshot files, equivalence's 54 of
55. What moves by rounding: every call with a drawn forest, whatever stands beside it, student(df) on a forest with no basis included, bcf() among them,
a seed able to fork. What moves by design: a write of `sd` on a drawn forest with no basis, before the first sweep or in mid-chain.

## Landing order

forests-shared-total with slice F comes first, and only forest-door's slice 2 and setmodel-parameters' slice 4 need this first, for the write rule. Nothing forces this before control-migration's slice 1a: 1a edits no engine file, Inf passes its record's validity and its
bridge, and the R edits here apply to its R but for F's own sharer refusal (the critique's trial). Readiness decides between the two, and control-migration's "Order of landing" takes a line. coefficient-prior-live and forest-door slice 3 come after. Engine features merge serially.
[Calls made](#calls-made) 8 has the alternative on forest-door.

- forests-shared-total: read, not edited. It lands first. Its D8 line and its two dec-B542 mentions are edited by this item's records step, once it has
  landed. Its R sharer refusal and flat-prior refusal take `normal()` here.
- forest-door: reworded under dec-B550 on its branch. Its table's "`fixed()`, `student(df)`, `normal()`" row holds; its reader does not as forest-door-s1
  has it, and takes this slice's decoder, its "value in force" comparison too. Records step for its writer: slice 2's item 6 at creation has moved here,
  its door half stays; its tests' no-basis round trip is new and is the door-reader mutant's killer. Its `validateForestAmplitude` hunk drops out.
- setmodel-parameters' slice 4: reworded under dec-B550 on its branch. Its record of the scale creation set is still needed, since `k.scale` reads the
  `sd` in force. Records step for its writer: its pin of the one candidate gets a `normal()` row, and a pin of a bare `$setLeafPrior(normal())` on a forest
  stating `amplitude = normal()`, which dec-B540 does not name and which reads like a restatement of the coefficient's prior.
- control-migration: `dbartsForestRecord`'s `amplitude.df` holds Inf, length 1; textual meeting in `forestStatements`.
- coefficient-prior-live (dec-B511): [Every live path](#every-live-path). Its TODO entry is rewritten at this item's records step.
- bartCause: no code. man/bcf.Rd's lines on its own branch when this lands; its tests run against this library as a gate above. Slice F's test fix,
  branch bcf-tests-state-prior, comes first. stan4bart, treatSens, bairrtt: not touched. No shipped header changes.

## Design notes and records

- [1. The rule](../design/forest-multiplier.md#1-the-rule) of forest-multiplier.md: the table's rows and columns: the trees have `sd` on every row, the coefficient's
  prior by kind and option as above with median 1, the start, a row for `normal()`, the two constants. Its pushes keep the record of how the map got here.
  [7. What the engine keeps general, and what is left for later](../design/forest-multiplier.md#7-what-the-engine-keeps-general-and-what-is-left-for-later)
  keeps its start item.
- multiplier-combiner.md: "SCALE-MIXTURE amplitude carries it in `amplitudePriorScale` ... FIXED-VARIANCE amplitude carries it in the node scale" is now one
  channel, rewritten, with the half-Cauchy lines around it.
- bcf.md: the mapping paragraph, which gives the node scales and a's median, and a note that dec-B18's parity holds for the prior and not for any
  coefficient's units: a, b0, b1 and the raw forests do not read on the bcf package's scale. bcf 2.0.2 itself puts `sd_control` on the trees, with a
  coefficient starting at 1.
- forest-sd-unit.md: "The calibration map is not touched" is history; the note says where a forest's `sd` goes now.
  [The rule](../design/forest-defaults-by-kind.md#the-rule) of forest-defaults-by-kind.md if it names the opt-ins. No note of its own.
- [`forest`](../../man/forest.Rd) and the rest of Help as above. TODO forest-size-reader stays. TODO coefficient-prior-live takes the words above. TODO
  written-sd-fresh-gap stays.
- A ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, TODO forest-amplitude-normal removed.

## Calls made

1. The start is a local handed to the combiner's constructor, from a field of the entry, and not a field the combiner keeps. Nothing after creation reads
   it, and a fixture's direct combiner stays neutral, 3 tests/cpp checks fewer. Alternative: `ForestAmplitudePrior::startScale`, as the prototype had it:
   one line fewer, in the struct a live write copies, the wrong side of dec-B550's "a field separate from the prior".
2. The statement is a +inf df, not a new field, because it is the smallest way to state it: nothing is added at the bridge, the record, the reader's column,
   `forestSharers` or control-migration. The maintainer treats normal against t as structure, so it is an internal code only, its comment says so, no t arm
   runs on it, and the summary tells him. Alternative: `ForestStatement::amplitudeNormal`, a bridge vector, a record slot and a reader column, ~40 more lines
   in five places and three other plans. Cost taken: `amplitudeDf` is not only a df, three-valued, with NA or Inf for one prior by kind, which the decoder and
   predicate hold in one place.
3. A forest with no basis keeps slice F's start by position times 1 / m: 0 where it is the only forest with no basis and a basis forest stands first. Moving
   it would shift `student(df)`'s draws there, against the ruling's chain. A start of 0 is not absorbing. TODO drawn-coefficient-start-by-kind stays open.
4. It lands after slice F, as dec-B550 says, and before forest-door's slice 2, see call 8. Against control-migration 1a nothing forces an order.
5. dec-B528 refuses a constant column naming `forest(x2)`. Its stated cost was that a sign-symmetric coefficient on a constant can no longer be written; this
   gives it a spelling, so the refusal with `sd` stated says "...; for a drawn coefficient, forest(x2, sd = , amplitude = normal())". Words only, ~3 lines, and
   the summary tells him since the words are a ruled refusal's. Where the call has plain terms the clause's advice meets a second refusal, which names the
   right spelling (`normal-column`). Alternative: leave the words as ruled.
6. `normal` joins the forest vocabulary, so a caller's own function of that name called inside `forest()` or `forests =` is no longer reached. With `normal <-
   function(v) 2 * v`, `forest(x1, sd = normal(.5))` runs on main and is refused after; with `normal <- function() dbartsPriors$fixed()`, `amplitude = normal()`
   goes from held to drawn with no message. `fixed`, `student` and `invchi` already shadow this way, `forest()` never shipped, and a caller's variable or column
   named `normal` stays theirs. Pinned (`normal-shadow`), said in the help and the summary. Alternative: resolve `normal` only at `amplitude =`, ~15 lines
   against `forest()`'s one vocabulary, with `normal` then differing from the other three.
7. No latent exact mode. Alternative: ~60 lines and a few minutes in CI; the check it adds is one precision in code the arms and bcf-latent-exact's other modes
   already run.
8. It lands before forest-door's slice 2, taking that slice's item 6 at creation, ~10 lines and a pin. One validator serves both kinds. Alternative: after it, no
   lines moved, waiting on control-migration 1a and ~16 agent-hours.
9. `student(Inf)` stays refused, its message naming `normal()`, as dec-B510 ruled the rule. Alternative: take it as the normal, 2 lines.
10. `aPriorScale` keeps its name on the K = 2 fixture spelling though it is now forest 0's multiple, not a mixture's scale: a rename touches ~40 test lines and
    nothing a user reads. Its comment says so. Alternative: rename to the multiple it is.
11. `prior.sd.of` on a drawn forest says "(at a coefficient of 1)" and not "(prior median)". `sd` is not the prior median of the forest's contribution, which is 0.54 `sd` under `normal()` (dec-B550's own correction; rerun 0.541), and the
    words said it was. "level difference (at a coefficient of 1)", "forest, per unit of basis (at a coefficient of 1)", "forest (at a coefficient of 1)". ~8 pins. Alternative: keep "(prior median)", which is false, or drop the qualifier, losing the mark of a drawn forest.
12. On a basis forest `normal()` records NA, as forest-door has it, so its print and `$getLeafPrior(f)$amplitude` show nothing there and show `normal()` on a forest
    with no basis, where it is not the default. `$getForest` shows it on both. Alternative: Inf on both, a default spelled out changing the record, no gain.
13. `k.scale` keeps being the leaf scale read back, now equal to `$sd` on every forest but a sharer, and no reader of the coefficient's prior scale is added: it is a
    constant. The size reader stays TODO forest-size-reader.
14. The `amplitude.prior.scale` and `amplitude.prior.variance` columns stay: tests read them and they say which row a forest is on. Alternative: drop, ~10 pins.
15. A new slot from a width change stays 1, the multiplier's median, against the earlier form's scaled slot. No R path widens a forest today.
16. No `larger-sampler` pin: a short run has no power, at an effective sample size of 34 to 41 per 5,000. The full run stays a scratch script (scratch/fan-sw/outer2.R, sd's quartiles 0.869, 1.045, 1.286 against 0.875, 1.059, 1.318), not a gate. Alternative: a benchmarks/R gate, ~80 lines.
17. At an extreme `sd` a drawn forest with no basis keeps its start finite in the engine and is not refused. The base runs at `sd` 1e-160 and 1e300 on a drawn `student()` forest with no basis; with the start as written, 1 / m^2 overflows v's start or underflows it, and `train` and `sigma` are NaN with no message. A start is
    not a model, and a refusal would be a new restriction on a value the prior takes. The bound is 1e100 on the start scale and 1e200 on v's start: a finite v alone still gave NaN, since a sweep squares the coefficient, which overflows past about 1e154. At the multiples 1e-100 and 1e100 the fit reads the base's (rmse .289 and .284, .278 and .28, gaussian); past them the chain runs and the forest may be off for a while (.39, .324), with no NaN at `sd` 1e-320 to 1e300, gaussian and probit, df 1 and 3. Pinned at both values (`extreme sd`, `extreme-sd`). Alternative: refuse past some multiple, a limit with no reason in the model.

## Open calls

None. The earlier open call, a live switch between held and drawn, is answered by dec-B550: a hold sets the coefficient to its held value and leaves the leaves, a
release draws from 1.

## Not verified

- Every count of re-pins is on slice F's 7f53d1c78. F landed with only comments, words, one pin and renamed baselines past that, and this branch was rebased onto it.
- The port has the map, the start, what goes but `nodeScaleDivisor` (kept at 1 so tests/cpp compile) and the R that takes `normal()`. Not the predicate and decoder
  as one function each, the refusals' words, or the column clause.
- The door, `$getForest` and slice 4 are not built, so their rows are from those plans' text. coefficient-prior-live is not built: its rows are from reading and one
  hook on a prototype.
- sbc's bcf arms: run, see the [Landing note](#landing-note). They have no power on a fork beyond what a new seed has.
- The transient on a write is one data set, 8 seeds. Coverage of f or of a treatment effect at equal time is not measured: the kernel is the same up to rounding, and the
  start stays.
- How often a seed forks is counted on nine calls of 4 seeds at 4,000 sweeps; past that and on other data it is not. The ruling ran a binary bcf() over 72 seeds, a probit at `sd` 0.05 forked, and bcf-equivalence's 19 did not.
- No exact arm of this gate reaches probit or logistic (`normal: the same chain as a column of ones` and bcf-latent-exact stand nearest), a drawn sigma (every arm holds it), weights other than 0 and 1, or a
  forest beside sharers. No arm sees a median to better than about 0.5 percent, the literal pin's job. No new arm is planned for these.
- Test predictors with `normal()` were not run: the refusal is of every model of several forests.
- x86 was not run by the build. The leg at landing ran on ed24323eb, whose code is the landed code: tests/cpp 420 of 420 plain and under ASan, UBSan and leak checks with no report, the full suite 37043 of 37043, equivalence 55 of 55 within the statistical bound, bcf and multinomial within the cross-host tier 1 bound. New tests/cpp and tinytest line counts are estimates, raised for the record of plans here running low.

## Landing note

sbc's bcf arms, R 200 and L 150 on every arm: `bcf` thin 120 burn 72000, `bcf-weak` thin 120 burn 18000, the four latent arms thin 50 burn 12000. Both libraries' 12 runs took about 2.5 hours, the replications below 1.5 more.

The pass rule as written trips. On the driver's own seed the tree flags, past the band of 0.0924, what the base does not: `prog4` 0.1034 and `prog5` 0.0955 on `bcf-weak`, `prog2` 0.1237 and `p2` 0.1041 on `bcf-probit`, `p3` 0.1221 on `bcf-probit-weak`, and the ill-posed raw `a` 0.0994 on `bcf-weak`. The base flags what the tree does not on the same runs: `prog2` 0.1086 on `bcf-weak`, `p3` 0.1029 on `bcf-probit`, `a` and `abs.a` on `bcf-probit-weak`, `prog2` and `prog3` on `bcf-logistic-weak`, `b1.minus.b0`, `eff1` and `prog3` on `bcf`. `bcf-logistic` flags nothing on either.

The rule assumed one seed gives the two libraries one comparison. It does not: one sampler serves all 200 replications, so once the two chains part every later replication has another stream, and the two runs are two draws of a noisy verdict. Read on `bcf-weak`: the first replication's ranks stand near each other, `sigma` 57 and 65, `abs.a` 60 and 63, so its data agree and the chains parted inside its 36,000 sweeps; the second's are apart, `sigma` 32 and 145. So a gaussian chain too is the same to rounding only until it parts, a longer run than the summary's counts read. So those three arms were run on two fresh seeds on both libraries:

| arm | seed | base flags | tree flags |
|---|---|---|---|
| `bcf-weak` | 20261101 | none | none |
| `bcf-weak` | 20261102 | `sigma` 0.0969, `eff2` 0.1025 | `b1.minus.b0` 0.0958, `eff4` 0.1055 |
| `bcf-probit` | 20261101 | `b1.minus.b0` 0.0925 | `eff1` 0.1029 |
| `bcf-probit` | 20261102 | none | `p1` 0.1121 |
| `bcf-probit-weak` | 20261101 | none | none |
| `bcf-probit-weak` | 20261102 | `a` 0.0965, `b1.minus.b0` 0.0944, `p3` 0.0966 | `a` 0.1019, `b1.minus.b0` 0.1046, `prog1` 0.0942 |

No functional the tree flagged on the first seed flags again, and no identified functional flags twice on the tree in any arm. Over the 12 runs, leaving out the two ill-posed functionals, the base flags 10 of 140 and the tree 9 of 140. The mean gap on `bcf-probit` is 0.0609 and 0.0633 over the three seeds, on `bcf-probit-weak` 0.0645 and 0.0641. So the read is a seed's noise, which the base shows as often, and not a bias. That reading replaces the rule's verdict and is a call for the review: the rule as stated would stop the build.

The discrimination check, the old generator against the new engine on `bcf-weak` at the same size: `abs.diff` 0.1740 flags against 0.0671, its chi-square p 0.000. `abs.a` does not separate them at this size, 0.0576 against 0.0697, and the small run's gap of 0.31 against 0.18 does not repeat. The generator's self-check reads the rule's IQR of 2 and b's sd of 1.049.

The burn ladders as a cost check, not gating, 40,000 sweeps over 28 data sets: `burn-bcf-probit` 181 us a sweep, `burn-bcf-logistic` 206. Read as cost only: no burn was adjudicated.

### The rule run and its trip

dec-B557 replaced the pass rule above: a replay of all six arms, four fresh seeds pooled, and a flag count. The replay first
read as failing on every arm because its length counted thinned blocks, so it ran 20 times longer than stated, past where
the two chains part; it counts sweeps now.

- The replay, one block of 20 sweeps over 40 replications: all six arms agree, worst 5.5e-10.
- Seeds 20261113 to 20261116 on both libraries, 48 runs, 8.9 core-hours, pooled to R 800: one functional of 70 fails,
  `bcf-probit`'s `abs.a`, the branch's gap 0.0773 against the band of 0.0693 at .05/70 (z +2.94), the base's 0.0464 (z
  +0.06). `bcf`'s `sigma` is past the z bar on both (-3.75 and -4.11), so it does not fail. The five functionals flagged on
  the driver's seed are all under their .01 band of 0.0561. The flag count at 0.0924: 25 on the branch, 23 on the base, of
  280, sign-flip p 0.931.

The diagnosis found no difference between the generator and the engine on either library and no inexact move. `abs.a`'s
ranks are U-shaped on both libraries, equally, on every bcf arm, most on `bcf-probit` and `bcf`. On eight further seeds
named before any fit, 20261123 to 20261130 (R 1600, band 0.0477), the split reversed: the base past the band at 0.0538,
the branch not at 0.0441. Over all 22 seeds (R 4400, band 0.0295) both are past it, the branch at 0.0469 (z +2.59) and the
base at 0.0493 (z +1.61). The rule's "past the band and the base is not" turns a deficit both share into a coin flip over
which library crosses first.

The reading: the harness starts each chain at the truth in the glue with fresh trees, and the burn is too short for the
trees to take up the signal at large |a|, so the drawn |a| runs small (z positive) and `bcf`'s sigma large. Two stale
starts turned up on the way, a's mixing variance and, under gaussian, sigma, each carried over from the last chain; each
BCF replication now installs them with the truth, a's variance drawn given the installed a. The runs above predate those
fixes, but both are redrawn every sweep from their exact conditionals, so neither can move ranks at these burns. dec-B562 landed the item on that reading. The rerun it called for, after landing, read the U otherwise: three times the
burn (36000, 4 seeds on each library) left the end-decile share at .275 and .259, and a start at the full truth, trees
included, gave .243 with no inverted U. The U is the Cauchy tail's mixing, which both libraries share: past |a0| of 3,
about a fifth of the prior, |a| has an effective sample of 5 to 23 in 150 draws. TODO bcf-latent-mixing carries it.
