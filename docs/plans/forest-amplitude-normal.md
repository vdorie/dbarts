# forest-amplitude-normal: `sd` on every drawn coefficient's prior, and a drawn normal on a forest with no basis

Status: PLANNED 2026-10-30 (dec-B539), rewritten for dec-B542 on 9f86fc6ff: wherever a coefficient is drawn, `sd` sets that
coefficient's prior scale, on every kind of forest and under `normal()` and `student(df)`, the trees staying at the response's
unit; where it is held, `sd` sets the trees. A fresh blind critique reads this form: the earlier ones read other designs. One call
is open, and it does not hold this item. Rulings: dec-B542, dec-B539 (the coefficient), with dec-B510 (`sd` the prior median,
`student()`), dec-B511 (the coefficient's prior live, its own item), dec-B513 (`fixed()` on a lone forest), dec-B522 (`normal()`
has a value to write, no return spelling), dec-B528 (a constant basis column refused), dec-B540 (`$setLeafPrior` on the one forest
that draws). The runs behind this form are in scratch/fan-cons and scratch/fan-plan/runs4.

agent: one opus implementer (effort high), one opus review that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)).
rng: no draw is added, removed or reordered on any path. Every sampler with a drawn basis forest moves by rounding, so a seed can
fork where rounding flips a decision: one of 16 calls forks near draw 100, its fits matching over 8 seeds. Bitwise for every call
with no drawn basis forest, `student(df)` on a forest with no basis included, with one corner: a caller's own function named
`normal`, called inside `forest()` or `forests =`, stops being reached ([Calls made](#calls-made) 8). Posterior-defining for the
one new model, `normal()` on a forest with no basis, gated by the exact arms below.
window: for 1.0-0, after forests-shared-total's slice F (dec-B542). Merged once, one slice.
budget: ~1,100 lines in the tree, ~150 of them package code: engine ~95 touched (69 measured on a patched copy, 63 in combiner.hpp
and 6 in chain.hpp, then the layout slot and comments, less `setAmplitudePriorScale`), R ~55 (29 measured, then the refusals,
predicate and decoder), help and NEWS ~60, design notes ~60, tests/cpp ~200 new and ~100 re-pinned (35 messages), tinytest ~280
new and ~90 re-pinned (43 sites), exact arms and readers ~130 (125 measured), records ~25. Agent-hours, built up from the measured
counts: engine and R 2, tests/cpp 5.5 (35 messages at ~5 minutes, and the new file), tinytest 6.5 (43 sites at ~5 minutes, and the
new file), arms 1, help and notes 2, gates 2, review of 31 mutants and a fix round 5: ~24, range 20 to 28, above the ledger's 15
to 22. The rates are guesses; the counts are not. Slice F may lower the counts.

## Summary for the maintainer

You ruled that wherever a coefficient is drawn, `sd` sets its prior scale, and that where it is held `sd` sets the trees
(dec-B542). So this item does two things. It adds `forest(amplitude = normal())` on a forest with no basis. And it moves a basis
forest's drawn coefficient, under `normal()` and `student(df)`, onto the rule `student(df)` on a forest with no basis already
follows. After it, doubling `sd` doubles the drawn coefficient's prior scale on every kind, the trees stay at the response's unit,
and a live write moves no leaf and does not jump the fit.

Nothing about the model changes for what runs today. A basis forest's product of coefficient and trees has the same prior and the
same chain to rounding, because its coefficient starts where the chain stays. `student(df)` on a forest with no basis does not
move by a bit. The engine gets about 95 lines touched and no new prior: the normal on a forest with no basis is one fixed variance
on the one coefficient update there is.

What it costs:

- A basis forest's coefficients and raw trees read in new units. At the default on a factor, b0 and b1 are 1.48 times larger and
  the trees smaller by that factor. Fits, contributions and bartCause's effects do not move.
- Seeded draws move by rounding wherever a basis forest is drawn. bcf-equivalence and equivalence, where one of 55 cells moves,
  are re-recorded, 86 tinytest and 100 tests/cpp checks re-pinned on units, and it comes to about 24 agent-hours, above the 15 to
  22 you were shown.
- The record understated one thing: a write before the first sweep has to move the start, and that reached 13 checks, not three.
  They all pass on a patched copy once the write moves the start. One pin of a pre-write state changes meaning: see [What a write
  before the first sweep pins](#what-a-write-before-the-first-sweep-pins).
- A live write on a basis forest does better than main in every cell written up, and worse in one written down: a numeric basis
  under probit, `sd` 7 to 0.7, reads 1.20 1.11 1.07 over sweeps 1-5, 6-20, 21-80, where main reads 1.09 0.99 1.00.
- forest.Rd's own example, `indexSd`, reads `k.scale` as the trees' scale and would go wrong with no error. No test catches it;
  this item fixes it and adds one.

Three smaller things you would want to know:

- Inside the engine `normal()` is stated as a df of +inf, an internal code. No t is drawn, no auxiliary kept, no quantile taken,
  and users never see it: you have said one does not sample a normal through a t, and this does not.
- Adding `normal` to what `forest()` reads means a function of your own named `normal`, called inside `forest()`, is no longer
  reached, as is already so for `fixed`, `student` and `invchi`. No released call is affected.
- dec-B528's refusal of a constant column gains a clause naming `forest(x2, sd = , amplitude = normal())`, words on a ruled
  refusal.

One call needs you, and it does not hold this item: dec-B542 says `fixed()` is chosen at creation so that nothing moves on a
running sampler, but dec-B254 and dec-B511 rule a live switch between held and drawn. Under dec-B542 that switch moves the trees'
scale by about `sd` over the unit. See [Open calls](#open-calls).

## Goal

`forest(amplitude = normal())` on a forest with no basis draws its coefficient under a normal prior, taken where `student(df)` is
taken and refused where it is refused, with no `sd` taking the median `student(df)` takes there. `sd` is the coefficient's prior
median in the response's units on every drawn forest: exactly on a forest with no basis, and through bcf's 0.674 on a forest with
a basis. The trees stay at the response's unit wherever a coefficient is drawn. Where it is held, `sd` is the trees' size. At the
live door `amplitude` is the value in force only, until coefficient-prior-live.

## What 0.9-34 did

It had no `forest()`, no `amplitude`, no coefficient and no model of several forests. Nothing is owed. `normal(k)` existed there
only as a leaf prior, which it still is at `leaf.prior`.

## What exists today

Run on main's engine at 9f86fc6ff, on private libraries built from scratch copies.

- [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) has one conditional: `priorPrecision = 1 / prior.variance`, a
  row-contracted design, one standard normal per coordinate. It does not read a forest's kind. A forest of kind none has an
  all-ones basis column ([`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor), so it is the block of width
  1 a numeric forest has. [`drawAmplitudes`](../../src/bartcore/combiner.hpp) adds the variance's inverse-gamma draw only where
  `mixtureScale > 0`, and `mixtureScale == 0` is the fixed-variance normal. No code decides that a forest of kind none must mix
  except the last arm of [`forestMapEntry`](../../src/bartcore/combiner.hpp).
- What `sd` sets on main, by doubling it: the trees under `fixed()` on every kind, under the normal and `student(df)` on a basis
  forest, and on a sharer; the coefficient only under `student(df)` on a forest with no basis, and main's left-out default, which
  slice F removes.
- [`forestMapEntry`](../../src/bartcore/combiner.hpp) is resolved at three places only: creation,
  [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) and [`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp).
  `setResponse` does not re-resolve.
- Starts, by position, from the (1, 0, 1) slots [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp) inherits: a forest of
  kind none starts at 1 where it stands first and at 0 where a factor or numeric forest stands first, a factor's levels at (0, 1)
  standing second; a numeric forest's is set to 1 in the constructor (TODO drawn-coefficient-start-by-kind).
- A write today. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) re-resolves the entry, writes a mixture's scale through
  [`setAmplitudePriorScale`](../../src/bartcore/combiner.hpp), which is false at `mixtureScale == 0`, and always re-derives the
  leaf scale. [`restoreGlue`](../../src/bartcore/combiner.hpp) restores a prior variance only where it mixes: a fixed variance is
  model, not state, and stays the destination's.
- Multiplier bookkeeping. [`forestMultiplier`](../../src/bartcore/combiner.hpp) contracts the block by the row; the fit and every
  reader read the product. The row bound's `zeroMultiplierTolerance` is relative to the forest's largest, so it does not see a
  scale. Nothing reads a forest's kind after creation but [`holdAtKindValue`](../../src/bartcore/combiner.hpp) and a numeric
  column's scale readers.
- On main `forest(amplitude = normal())` is "could not find function". `student(Inf)` is refused.

## Did we build the wrong thing, and is there a simpler form?

No to the first: dec-B539's option was "a new prior in the engine with its own exact gate", and there is no new prior. The normal
path exists and takes any forest, so a forest with no basis needs one fixed variance and a statement.

The one thing the rule costs in code is that a basis forest's start has to keep its chain. Simpler forms fail:

- A basis coefficient starting at its kind's value with no scale: the chain moves for every basis forest, at every default, for no
  gain. Starting at the scaled value costs one field and a write rule.
- No write rule at all: a sampler whose `sd` is written before its first sweep stops matching one created at that `sd`, which 11
  checks pin.
- A flag saying whether a sweep has run: it needs a new field in state and format. The rule below reads the block itself, and a
  copy and a reload of an unswept sampler stay at the start with nothing added.

### Same prior, same chain, other units

(a, F) today is (a' / q, q F') for (a', F') under this plan on a basis forest, with q = m / 0.674. The conditional of each block
is equivariant in that rescaling, a row's zero bound is relative, a mixture's variance scales by q^2 and its scale by q. So if a
starts at q times its old start, the two chains are one chain up to rounding.

Run on a patched copy against main at one seed, 16 calls: a factor and a numeric basis, a column in other units, `student()`,
three forests, held forests, probit and logistic. 15 read main's chain to rounding with varcount identical. The 16th, a numeric
basis under probit, forks through rounding near draw 100, its fits matching main's over 8 seeds, largest |z| 2.56 over 200 rows.
`bcf-equivalence.R compare --cross-host` passes tier 1 on 18 of 19 with varcount identical, worst ratio 4.3e-4. The 19th,
`numeric_swap_write_copy`, writes `sd` on a running sampler and moves by design ([The transient on a
write](#the-transient-on-a-write)).

| read, numeric basis at the default | main | this plan |
|---|---|---|
| `k.scale` | unit / (0.674 sd(column)) | unit |
| `$sd`, and `print`'s "sd" | unit / sd(column) | the same |
| coefficient's prior | N(0, 1), starting at 1 | N(0, q^2), starting at q |
| `getForestAmplitudes` | a | q a |
| raw `getForestFits`, `extract(fit, "forest")` | F | F / q |
| their product, `extract(contribution = TRUE)`, the fit | | the same, to rounding |
| `$setLeafPrior` or `$setForest(sd = v)` writes | the leaf scale | the coefficient's prior |

### The two constants

z = 0.6744897501960817 on a forest with no basis, 0.674 on a basis (dec-B542). 0.674 is bcf's rounding; moving a basis to z would
shift every treatment forest's draws for no gain. Kind none under `student(df)` is today's arm, which reads z through t_nu,.75,
and the normal divides by z. The difference is one token, in the map.

### One spelling inside the engine

The statement is `ForestStatement::amplitudeDf = +inf`, not a new field. It is an internal code for normal: a forest of kind none
with a +inf df is drawn, with no mixture; no t is drawn, no auxiliary kept, no quantile taken, and the engine comment says so. R's
`student(Inf)` stays refused. What reads "does this forest state a coefficient prior" as `!isnan(amplitudeDf)` then needs no
engine line: forests-shared-total's `forestSharers`, `ForestCalibration::amplitudeDf`,
[`applyForestStatements`](../../src/R_interface_bartcore.cpp), the `amplitude.df` column, control-migration's `dbartsForestRecord`
slot, whose validity checks only length. Slice F's items 1 and 3 need no line only if the hold keys on `isnan(df)`, on the
statement and not on a local: forests-shared-total has been told.

What R does need: the df is a three-valued code (NaN, finite, +inf), and one prior reads NA or Inf by kind. One decoder and one
predicate hold that in one place (R items 4 and 5). `R/model.R` on the forests-shared-total tip already tests `dbartsStudentPrior`
at four sites beyond the validator, and forest-door-s1's reader still has main's `else if (!is.nan(df)) studentPrior(df)`, which
throws on Inf. A third prior that is not a t would force a real field. [Calls made](#calls-made) 2 has the alternative.

## Design

m is `sd` over the unit s, a numeric forest's default being 1 / basisScale; q = m / 0.674; z as above; kappa is the kind's
variance, 1/2 on levels and 1 otherwise.

### forestMapEntry

Every row is [`forestMapEntry`](../../src/bartcore/combiner.hpp) as it reads after this item. `df` is the forest's `amplitudeDf`:
NaN when none is stated, +inf for the normal's code, finite for `student(df)`.

| kind | coefficient | df | node scale | prior variance, or v's start | mixture scale | mixture df | start scale | against main |
|---|---|---|---|---|---|---|---|---|
| any | held | any | s m | kappa | 0 | NaN | 1 | unchanged |
| none | sharer, slice F | NaN | s m x share | 1 | 0 | NaN | 1 | unchanged by this item |
| levels | drawn | NaN or +inf | s | kappa q^2 | 0 | NaN | q | node was s q, variance kappa |
| levels | drawn | nu | s | kappa q^2 | q sqrt(kappa) z / t_nu,.75 | nu | q | node was s q, v and scale without q |
| numeric | drawn | NaN or +inf | s | q^2 | 0 | NaN | q | node was s q, variance 1, start 1 |
| numeric | drawn | nu | s | q^2 | q z / t_nu,.75 | nu | q | node was s q, v and scale without q |
| none | drawn | +inf | s | (m / z)^2 | 0 | NaN | 1 | new |
| none | drawn | nu | s | 1 | m / t_nu,.75 (m at nu = 1, the half-Cauchy arm, NaN df) | nu | 1 | unchanged, bit for bit |

A basis forest's +inf resolves as its NaN does: forest-door's record does not use it and a fixture may. Below 1 and -inf stay
refused. `ForestAmplitudePrior` gains `startScale`, and that is the only new field in the entry.

So `sd` is one thing on every drawn row: the prior median of the coefficient's size, in the response's units. It is not the median
of what the forest contributes, though dec-B542's first sentence reads that way: a x F is a product, and over 2e6 draws the median
of abs(a F) is 0.542 `sd` under the normal, 0.556 under `student(3)` and 0.587 under `student(1)`, and sd(a F) under the normal is
1 / z = 1.48 `sd`, against `sd` held. On a basis, the 0.674 makes it 1.0007 `sd` for the coefficient's size; with no basis it is
exact.

### The start

Each block starts at its kind's value, as today's position rule gives it, times `startScale`:

| block | today | after |
|---|---|---|
| numeric, drawn | 1 | q |
| levels, drawn, standing second | (0, 1) | (0, q) |
| levels, drawn, standing elsewhere | by position | that, times q |
| none, `student(df)` or `normal()` | 1 standing first, 0 where a basis forest stands first | the same |
| held, any | its kind's value | the same |

Why this keeps the chain: (a, F) and (q a, F / q) are one state in other units, so the first conditional is the same draw and so
is each one after it. A mixture's v starts at the entry's variance, kappa q^2, for the same reason.

A forest with no basis keeps its start by position: [Calls made](#calls-made) 4.

The rule for a write before the first sweep. `writeAmplitudePrior`, new, takes the re-resolved entry's prior, and is the only way
a write reaches the block:

1. If the block holds exactly, to the bit, its kind's start times the old `startScale`, it is still at its start: every coordinate
   takes its kind's start times the new `startScale`. Under a mixture, v takes the new variance too if v is still its old start.
2. Otherwise it keeps its coefficients and its v: a write moves no drawn coefficient.
3. Either way it takes the new `startScale`, `mixtureScale` and, where it does not mix, the new fixed variance.

No sweep flag and no state: the block says it itself, a chain lands on that exact value with probability zero, and a copy and a
reload of an unswept sampler rebuild the same start and restore it. A prior draw of the trees leaves the block at its start, so a
write after one moves the start as a write before it does: run, `identical()` to a twin on a factor and a numeric. Grow-from-root
draws the block, so a write after it is a write in mid-chain (run). A forest with no basis has `startScale` 1, so rule 1 changes
nothing there.

The combiner holds the kind's start beside the amplitudes (`kindStart_`) and v's start per forest (`startVariance_`), set in its
constructor. A width change through [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp) gives a new slot the forest's
`startScale` and refreshes `kindStart_`, since a new slot of 1 would otherwise read in the old units.

### Every live path

| path | drawn | held or sharer | why it holds |
|---|---|---|---|
| `$setForest(f, sd = v)` | writes the coefficient's prior through `writeAmplitudePrior`. No leaf, no drawn coefficient, no fit moves; before the first sweep the start moves too | writes the leaf scale, as today | [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) calls it in place of `setAmplitudePriorScale`, and still re-derives the leaf scale, which on a drawn forest does not change |
| `$setLeafPrior(forests = )` | the same, through the same call | the same | one writer |
| `$setLeafPrior(normal(sd = v))` on the one candidate | the same, by setmodel-parameters' slice 4 under dec-B540 | | no branch |
| `$setForestBasis(updateBasisScale = TRUE)`, numeric, sd unstated | re-resolves at the new scale, so q moves and the coefficient's prior takes it through `writeAmplitudePrior`; the leaf scale stays | leaf scale, as today | [`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp) calls it. The doc comment there, "no amplitude prior moves", is reversed |
| `$setForestBasis` otherwise | nothing moves | nothing moves | no entry is re-resolved |
| a write of the `sd` in force | nothing moves | nothing moves | the same entry; a start scale equal to itself |
| a copy | rebuilds the entry at the sd in force and at `basisScale`, then restores block and v | as today | the fixed variance is the destination's, and the destination's is the written one |
| a reload | as a copy | as today | run: test-forest-basis-scale-update.R's copy and reload against live, which a patch without the restate write failed |
| a state restored, taken before a write | the block and v go back to what the state holds: the old start, if unswept. The prior stays the written one | as today | a state carries no prior. A block restored this way is not at its new start, so it is one more legal start, not the twin's. The pin moves: [What a write before the first sweep pins](#what-a-write-before-the-first-sweep-pins) |
| a donor warm start | refused on several forests, as today | | test-multiforest-warmstart-refusal.R |
| `growFromRoot`, `sampleTreesFromPrior` | run from the block's start, F at the unit: main's chain to rounding, within 5e-14 with varcount identical, on a factor and a numeric, at the default and under `student(3)` | as today | product and conditional as a sweep's |
| a saved state from another build of this development version | loads as it stands: its block and trees are the old units' same product, and the next draws settle into the new ones. No migration, and `stateFormatVersion` does not move | unchanged | only shipped versions are owed compatibility, and 0.9-34 had no coefficient |

### Engine, ~95 lines touched

1. [`ForestStatement::amplitudeDf`](../../src/bartcore/combiner.hpp): a +inf df states the normal coefficient, and on a forest
   with a basis means what no df means (3 lines of comment).
2. [`forestMapEntry`](../../src/bartcore/combiner.hpp), as the table above: `normalStated` and `dfStated` kept apart, so no local
   that reads "states a df" excludes the normal; the basis arm writes q into the variance, the mixture scale and `startScale`, and
   no longer sets a node factor or divisor; the +inf arm on kind none.
3. `ForestAmplitudePrior::startScale`, `AmplitudeForestCombiner`'s constructor applying it after
   [`holdAtKindValue`](../../src/bartcore/combiner.hpp), and `kindStart_` and `startVariance_`.
4. `writeAmplitudePrior`, new, on [`ForestCombiner`](../../src/bartcore/combiner.hpp) and
   [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp), replacing `setAmplitudePriorScale`. No facade virtual, so
   [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move.
5. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) and [`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp)
   call it.
6. [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp): a new slot's start.
7. Comments that state the old rule: the block above [`forestMapEntry`](../../src/bartcore/combiner.hpp), and
   [`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp)'s.

Not `drawForestAmplitude`, `drawAmplitudes`, `restoreGlue`, `glueIsValid`, the k draw, the rescaling step. The half-Cauchy arm
stays as it is. `aPriorScale` keeps its name: on the K = 2 fixture it is still the mu coefficient's half-Cauchy median.
`stateFormatVersion` does not move.

### Bridge

No code. [`applyForestStatements`](../../src/R_interface_bartcore.cpp) passes a +inf through, and its comment says what it means.

### R, ~55 lines

1. `normal` in [`forestConstructors`](../../R/model.R), [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R),
   [`processHit`](../../R/formulaTerms.R)'s vocabulary and [`forceCallerCode`](../../R/family.R)'s hint. Run: a caller's own
   variable or column named `normal` stays theirs in `basis =` and in a forest's predictors, as `student` does.
2. [`validateForestAmplitude`](../../R/model.R): a bare `normal` is the constructor; `normal()` kept as the object; `k` or
   `prior.sd` refused, as the adviser's second item has it, by the words in [What the user sees](#what-the-user-sees); the
   refusal's list.
3. [`forestStatements`](../../R/model.R) gains the kinds and records +inf for `normal()` on a forest with no basis and NA on one
   with a basis. Its caller in [`resolveSamplerSpec`](../../R/spec.R) passes `hasBasis`. A default fit's record stays as it was.
4. One predicate, `statesCoefficientPrior(amplitude)`, true of a `student()` or `normal()` object, in place of every `is(x,
   "dbartsStudentPrior")` that asks whether a forest states its coefficient's prior: [`resolveForests`](../../R/model.R)'s lone
   branch, forests-shared-total's sharer refusal and `sharingForests`, and whatever step 1's grep finds. The refusal's words take
   the spelling written. [`validateForestAmplitude`](../../R/model.R) and [`forestStatements`](../../R/model.R) keep the class
   test, which is about which class.
5. One decoder, `statedCoefficientPrior(df, held)`: `fixed()` where held, `normal()` at Inf, `studentPrior(df)` at a finite df,
   NULL where none is stated. [`reportLeafPrior`](../../R/dbarts.R) uses it, and forest-door's `$getForest` and its "value in
   force" comparison must, turning NULL into the kind's default. [`printedForest`](../../R/forestPrint.R): "; coefficient normal".
6. [`studentPriorRefusal`](../../R/model.R)'s words. Where forests-shared-total's slice 3 has landed, dec-B528's refusal with `sd`
   stated gains a clause: [Calls made](#calls-made) 5.
7. No reader changes code. [`$getLeafPrior`](../../R/dbarts.R)'s `k.scale` reads the map's `prior.scale`, the unit on a drawn
   forest, and `sd` reads the map's `map.sd`, so both are already right. The size reader stays TODO forest-size-reader.

## What the user sees

Every taken-or-refused row below ran on a patched copy against `student(3)` and `fixed()` beside it, and `normal()` matched
`student(3)` in every one.

| call | `normal()` | run |
|---|---|---|
| beside a factor; a numeric; a factor and a numeric, gaussian | taken, no sd | median in force 2.5994 on main's map, slice F's value after |
| the same with `sd = 0.7` | taken, median 0.7 | |
| probit, logistic, no sd and with sd | taken | sd 1 and 1.8138 on main's map; 1.574 and 2.85 after slice F, its item 3 |
| a weighted fit | taken | exact arm with rows of weight zero; unit read over rows in the likelihood |
| beside two sharers with `sd`, plain, and also beside a factor | taken, outside the share | sharers 1.1993 each |
| with no `sd` beside a forest stating `sd` | taken | |
| with no `sd` on any forest of two or more unsized | refused | below |
| a lone forest, with or without `sd`; the multinomial entry | refused | below |
| a variance forest | refused | whole model, today's words |
| `bart()`, `dbarts()`, a `forests` list, `predict(bases = )` | taken | |
| `xbart()`, `rbart_vi()` | refused | no `forest()` term, today's words |

The refusals, each as `student()`'s at the same place, tested by exact pattern:

- A lone forest and a multinomial entry: "this model has one forest, which has no coefficient to give a prior; a lone forest draws
  its size through the fitting function's leaf prior, as leaf.prior = normal(k = chi(1.5, 2)) or normal(sd = invchi(1.5, s)) with
  s a size in the response's units", or forests-shared-total D10's rewording, the one branch either way.
- A sharer: "forest 1: amplitude = normal() draws one forest's coefficient; these forests share one total: state each forest's
  'sd' to give it its own".
- `normal(k = 2)`, `normal(sd = 2)`, `normal(k = chi(1.5, 2))`, on every kind: "'amplitude = normal()' takes no 'k' or 'sd': a
  drawn coefficient's size is the forest's own 'sd'". Only a bare `normal()`, as forest-door has it.
- An `amplitude` of another kind, the string `"normal"` among them: "a forest's 'amplitude' must be fixed(), which holds its
  coefficient, normal() or student(df), which draw it under that prior, or left out; ...".
- `student(Inf)`: "a coefficient's student(df) takes one finite df of 1 or more; student() is student(1), the half-Cauchy, and
  normal() is the normal".
- A live change of `amplitude`: today's "$setLeafPrior restates only a forest's 'sd': 'amplitude' is fixed at creation", and at
  `$setForest` forest-door's words for every changed `amplitude`, until coefficient-prior-live.

### Readers

| reader | a basis forest, drawn | no basis, `normal()` | no basis, `student(df)` | held |
|---|---|---|---|---|
| `$getLeafPrior(f)$k.scale` | the unit; was sd / 0.674 | the unit | the unit, as today | sd, as today |
| `$getLeafPrior(f)$sd`, `$sd.stated`, `$prior.sd.of` | unchanged | as `student(df)` | as today | as today |
| `$getLeafPrior(f)$amplitude` | unchanged, NULL at the default | a `dbartsNormalPrior` | as today | `fixed()` |
| `getK(f)` | 1 | as `student(df)`: 1 on main's map, slice F decides | as today | as today |
| `$getForestAmplitudes(f)`, `fit$glue` | q times today's; prior median z q | prior median sd / unit | as today | held at kind's value |
| raw `$getForestFits`, `extract(fit, "forest")`, `predict(type = "forest")`, `fit$forests` | today's over q | at the unit | as today | as today |
| `extract(contribution = TRUE)`, the fit, `predict(type = "ppd")` | unchanged to rounding | new | unchanged | unchanged |
| `print(sampler)`, `print(fit)` | unchanged: `$sd` | "1 forest1: no basis; sd 0.7" | as today | as today |
| `summary(fit)`, the verbose creation summary | unchanged | "...sd 0.7 (stated); coefficient normal" | as today | as today |

- `forest(sd = 2, amplitude = normal())` prints `amplitude = normal()`. `fit$fixed$amplitude`: no entry, as for anything drawn.
- `$getForest(f)`: `amplitude = normal()`, through this slice's decoder, forest-door's "value in force" comparison likewise.
- bartCause, no code. `icate` is `response.scale * forests$tau * bcfGlueDifference(glue)` and the counterfactuals use (b.1 - b.0)
  times tau: products, so they do not move, and its test-14-bcf.R, which reconstructs through those products and compares with the
  samples, passes on a patched copy, 184 expectations. `fit$glue[, c("b.0", "b.1")]` and `fit$forests$tau` read in new units, 1.48
  times larger and smaller by that factor at the default. `fit$glue[, "a"]` and `fit$forests$mu` do not move: kind none under
  `student(1)`. Its `bcf` help's words on `sd.control` stay true; a line on b and tau goes in on its own branch.
- stan4bart: not touched. Its mvbart forests couple only across equations, and it uses no `forest(basis = )` or `amplitude`.
  treatSens and bairrtt: not touched. No shipped header changes.

### Help

- [`forest`](../../man/forest.Rd), the main change:
  - `amplitude`: `normal()` beside `student(df)`, its median and its reach; "Where it is drawn ... `sd` is its prior median" gets
    "of the coefficient's size, in the response's units, the trees staying at the unit".
  - The half-Cauchy paragraph: "the prior median of what the forest contributes" is wrong for a product (above) and goes, and
    "With no basis the coefficient's size is the forest's size" is rewritten, a basis forest under the same sentence.
  - The adviser's sentence, once: at one `sd`, a drawn normal's contribution has about 1.48 times the standard deviation of a
    `fixed()` forest's, since the coefficient and the trees are both random, and thinner tails than `student(df)`'s.
  - The size section's 1.049 sqrt(K) budget, reread: it reads products, so it should hold.
  - The example: `indexSd <- sqrt(calibration$k.scale^2 * 0.5 * rowSums(basis^2))` becomes `calibration$sd / 0.674`, with its
    comment. It reads the unit otherwise and is wrong with no error.
  - A live write of `sd` on a drawn forest changes the coefficient's prior, not the leaves, and the fit settles as the coefficient
    does. A caller's own function named `normal` is not reached inside `forest()`.
- [`dbartsPriors`](../../man/dbartsPriors.Rd)'s `normal` says the same in two lines.
- [`dbartsSampler-class`](../../man/dbartsSampler-class.Rd) and the docstring in [`R/dbarts.R`](../../R/dbarts.R), for codoc:
  "`k.scale` is the map's leaf scale" gets "the response's unit where the coefficient is drawn"; "A drawn forest with no basis
  carries its sd as the median..." and "a level's coefficient is drawn under prior variance 1/2 and a numeric column's under 1"
  become one sentence on every kind; the spread line reread.
- [`bartBT`](../../man/bartBT.Rd): `forestFits` and raw glue say they read in these units.
- NEWS: one clause in 1.0-0's forest item. Nothing released moves.

## Slices

One slice: engine, R, help, tests, arm, design note, then one review and the gate list. Steps, each verifiable:

1. Rebase on a main carrying slice F. Re-read F as built: its hold predicate and its no-sd median, each to key on `isnan(df)`.
   Grep R as built for `dbartsStudentPrior` and `amplitude.df`: slice F's R item 4 may have added sites, and each goes to the
   predicate or the decoder. Re-count the failures below on F's tip before re-pinning anything.
2. Engine items, then the 16 calls and `bcf-equivalence.R compare --cross-host` before any test is touched: tier 1 on 18 of 19.
3. tests/cpp, re-pinned and new. Run plain and under ASan and UBSan.
4. R items, `--preclean`, tinytest re-pins and file. Run the A/B against the base.
5. The exact arms, quick and full; every exact gate's read restated; sbc's bcf arm by hand.
6. Re-record bcf-equivalence and equivalence, each with a MANIFEST row.
7. Help, NEWS, design notes, lint, air, codoc, freshness.
8. Review with mutants, fix round, gates, records, CI on the branch, merge.

## Tests and gates

### What the failures pin

Counted on a probe of the map and start with no write rule, no `normal()` on a forest with no basis and no reader change: 93 of
32,380 tinytest results in 12 files, 106 tests/cpp checks. Then on a fuller patched copy with the write rule, the restate write
and the normal, which is this design's: every count below is read from those two runs, and all were on main's map before slice F,
which re-pins some. Step 1 re-counts.

**tinytest: 86 pin units and move by design; 7 pin something else, and all 7 pass on the fuller copy.**

| file | checks | sites | pins | class |
|---|---|---|---|---|
| test-bcf-family.R | 25 | 3 | `k.scale` x 0.674 is sd; `priorScales(fit)`; the index sd from hard-coded variances | units |
| test-forest-sd-unit.R | 23 | 18 | `inForce` and `expectStated` read `prior.scale` x divisor at every door: term, list, spec, data object, `bart()`, rows in the likelihood, `setWeights`, probit, logistic, writer, copy, reload, `new()`, `setState`, `setResponse`, other data | units |
| test-forest-leaf-prior-reader.R | 15 | 1 | `k.scale` x 0.674 is `$sd`, in a loop | units |
| test-forest-basis-r5.R | 10 | 10 | `priorScale(.)` on a swap; line 476's donor against recipient must read the coefficient's prior to keep its meaning | units |
| test-forest-scale-record.R | 4 | 3 | `expectRule` reads the column rule through the leaf scale | units |
| test-forest-units-twin.R | 3 | 1 | `sizeCheck` reads the trees' spread x z, and must read the coefficient's scale | units |
| test-forest-amplitude-student.R | 2 | 1 | the coefficient's median z at a numeric default; 0.40 after | units |
| test-multiforest-leaf-prior-writer.R | 1 | 1 | line 180, `k.scale` follows a write | units |
| test-bcf-creation.R | 1 | 1 | `k.scale` is `$sd` / 0.674 | units |
| test-calibration-midchain.R | 1 | 1 | `k.scale` is `$sd` / 0.674 | units |
| test-prior-init-composed-distribution.R | 1 | 1 | the start c(1, 0, 1) | start value |
| test-forest-basis-scale-update.R | 1 | 1 | line 76: a restate before the first sweep draws as its twin | write before the first sweep |
| test-forest-sd-unit.R | 1 | 1 | line 790: a write before the first sweep draws as its twin | write before the first sweep |
| test-multiforest-leaf-prior-writer.R | 3 | 1 | line 133: the same, gaussian with a fixed variance's sd, probit with both written, K = 3 | write before the first sweep |
| test-forest-basis-scale-update.R | 2 | 2 | lines 174 and 175: copy and reload against live after a restate before any sweep | a real regression, in the probe only |

The two at lines 174 and 175 are the one real regression the runs found in the engine: the probe's
[`Chain::restateForestBasisScale`](../../src/bartcore/chain.hpp) wrote no variance, so a live sampler and its copy disagreed.
Engine item 5 cures it and they pass. No other check names a regression.

The fuller copy adds two to re-pin: 111 of 32,475, on a run at home, whose 23 further failures are `sizeCheck`'s at-home arm, the
same reader.

- test-forest-arguments.R line 249: its refusal list has `normal()`, taken now. Rewritten with R item 2, as dec-B539 needs.
- test-multiforest-leaf-prior-writer.R line 251: [What a write before the first sweep
  pins](#what-a-write-before-the-first-sweep-pins).

**tests/cpp: 100 pin units or start values and move by design; 6 pin a write before the first sweep, and all 6 pass on the fuller
copy.** The fuller copy shows 103, the 100 and 3 for the normal, below.

| messages, by group | checks | pins | class |
|---|---|---|---|
| `map statements`: its leaf scale; a factor takes one unit over 0.674 | 36 | the entry's node factor and divisor | units |
| `calibration map`: the reported scale, fixed-variance report, a zero column | 14 | the reader's product | units |
| `sd unit`: six messages at 3 each; `a write of the sd in force` | 20 | quotient on the anchor, carried in the leaf scale | units |
| `numeric rule`: factor beside it; one unit per sd; prior precision 1; starts at 1; made again; after a write; after a swap; standing second | 14 | quotient, precision 1, start 1 | units and start values |
| `student`: scale at df 1 and 3; reports; v starts; median matched; no df is today's entry; donors' v | 11 | the t's scale without q, v's start kappa | units |
| `held coefficients`, `and the forests' defaults`, `relative row bound` (fails on `start[0] == 1`, not on the bound), `a unit treatment multiplier`, `BCF combiner neutral (1, 0, 1)` | 5 | units, or a start read before the test's point | units and start values |
| `map writers: a written twin` 2; `sd unit: draw for draw` 3; `basis rescale: the rescaled sampler draws as its twin` 1 | 6 | a write before the first sweep | write before the first sweep |
| `student: a df below 1 or not finite throws` | 3, fuller copy only | the `bad` list has +inf, the normal's code now | dec-B539's |

So the ledger's "three of them pin" a write before the first sweep is 11 in all: 5 tinytest and 6 tests/cpp, and 13 with the copy
and reload checks that follow a restate.

### What a write before the first sweep pins

Twin pins hold: a sampler whose `sd` is written before its first sweep matches one created at that `sd`, `identical()` on a factor
and a numeric, on a chain, after a prior draw of the trees, and through a copy and a reload (run).

One pin moves. test-multiforest-leaf-prior-writer.R line 251 stores a state, writes, restores, and expects a twin that wrote,
stored and restored. The state holds the old start, and the twin's holds the new one, so the draws differ; reader and prior still
match. The pin was that a state holds no prior, so restoring undoes no write, and it still holds: it is re-pinned on a twin that
stored before writing too. A state restored at the old start is a legal start, not the twin's: [Calls made](#calls-made) 10.

### tests/cpp, new

`testForestAmplitudeNormal` beside [`testForestAmplitudeStudent`](../../tests/cpp/test_sampler.cpp), ~200 lines:

- `map: every row`: each row of [forestMapEntry](#forestmapentry), to the last place, by a sample-free read, kind by kind, family
  by family at no sd, and a stated sd through the unit; the node factor and divisor 1 on every drawn entry.
- `normal: a basis keeps its own entry`: +inf and NaN `sameEntry` on levels and numeric, 0.674 and not z. `normal: a held forest
  ignores it`. `normal: what is refused`: -inf, 0, 0.5 and -1 throw.
- `student: no basis does not move`: the entry `sameEntry` as main's values, nu 1 and 3.
- `start: each block`: every row of [The start](#the-start), with and without a mixture, standing first, second and third.
- `start: the chain stays`: entries built by hand at node scale s q, variance kappa and the kind's start, against the map's at
  node scale s, variance kappa q^2 and start times q: the fit's draws within 1e-9, on a factor, a numeric and each under
  `student(3)`.
- `write: before the first sweep`, normal and student, block and v: moved, the twin's bits. `write: after a prior draw of the
  trees`: the same. `write: after a sweep`: no coefficient, no v, no leaf moves; the prior is the new one; the leaf scale's bits
  stay. `write: the sd in force`: nothing moves. `write: a held forest`: the leaf scale, as today, block unmoved.
- `write: restate`: `restateForestBasisScale` writes the prior, the leaf scale stays, then a sampler made at that scale draws
  equal.
- `layout: a new slot`: widened through the combiner's `setForestBasis`, the new slot starts at `startScale`, and a write sees the
  block at its start.
- `normal: the same chain as a column of ones`: kind none under +inf at a multiple of z against kind numeric on a ones basis at a
  multiple of 0.674, same counts: each reads node scale s, variance 1 and start 1, and the draws are equal. It survives dec-B528,
  which refuses at R.
- `normal: no sharer`: `forestSharers` excludes it, with sd and without. `normal: the no-sd median`: equal to df 3's multiple at a
  held k, a drawn k and a named sd (slice F's three cases).
- `normal: the calibration reader`: `amplitudePriorVariance` the entry's, `amplitudePriorScale` NaN, `amplitudeDf` +inf.

### tinytest, new

A new test-forest-amplitude-normal.R, small n, ~280 lines:

- `sd-carried`: the sibling sweep dec-B542 rests on. On one data set, doubling `sd` under `fixed()`, `normal()` and `student(df)`,
  on every kind where each is taken: `k.scale` doubling where held and unmoved where drawn, the coefficient's prior doubling where
  drawn, `$sd` doubling everywhere. So no exception can come back.
- `drawn-live`: `$setLeafPrior(forests = )` in mid-chain on each kind and option: `getForestFits` and `getForestAmplitudes`
  unmoved, `k.scale` unmoved, `$sd` moved; before the first sweep and after `sampleTreesFromPrior`, `identical()` to a twin;
  through `$setForestBasis(updateBasisScale = TRUE)` likewise; then a copy and a reload continue.
- `index-sd`: the example's `calibration$sd / 0.674` against the product's spread at prior, by simulation, within 2 percent: no
  check covers it now, and it is how forest.Rd's example could go wrong silently again.
- `normal-median`: `$sd` and `$prior.sd.of` equal to `student(3)`'s, `$k.scale` the unit, no sd and sd stated, under gaussian,
  probit, logistic.
- `normal-reader`: `$amplitude` the class; print and summary lines; `forest()`'s print; `getForestAmplitudes` moves; no
  `fit$fixed$amplitude`.
- `normal-record`: `amplitude.df` Inf on the forest stating it, NA beside it, absent on a fit stating none, on a basis forest's
  `normal()` too.
- `normal-doors`: formula, `forests` list, a data object's `bases`, `bart()`, a bare `normal`, draws `identical()` across the
  spellings.
- `normal-siblings`: beside one basis, beside two; weights; at a stated `sd` beside sharers, sharers' parts and shares unmoved.
- `normal-refusals`: each pattern in [What the user sees](#what-the-user-sees), counted with `withCallingHandlers`; a variance
  forest; the multinomial entry; `normal(sd = 2)` and `normal(k = 2)` on a basis forest too.
- `normal-basis`: on a factor and a numeric forest, draws and record `identical()` to left out.
- `student-unmoved`: `student(df)` on a forest with no basis, beside a held forest, `identical()` to the base's recorded draws.
- `normal-shadow`: with a caller's `normal <- function(v) 2 * v`, `forest(x1, sd = normal(.5))` refused, "forest 'sd' must be a
  number, not a dbartsNormalPrior", in a formula and in `forests =`; with a caller's `normal <- function() dbartsPriors$fixed()`,
  `amplitude = normal()` is the package's and draws; a caller's variable or column named `normal` stays theirs in `basis =` and
  the predictors.
- `normal-column`, once forests-shared-total's slice 3 has landed: `y ~ a + forest(b, basis = one, sd = 2)` refused with the
  clause; its advice as written meets the plain-terms refusal, which names `forest(a) + forest(b, sd = 2, amplitude = normal())`,
  and that is taken.
- `normal-live`: `$setLeafPrior(forests = list(forest(amplitude = normal()), forest()))` refused as today.

### Mutants

| name | change | killed by |
|---|---|---|
| start-unscaled | a basis block starts at its kind's value | `start: the chain stays`; tier 1 |
| v-unscaled | v starts at kappa, not kappa q^2 | `start: each block`; tier 1 |
| mixture-unscaled | no q in a basis forest's mixture scale | `map: every row`; arms under student(df) on a basis |
| variance-unscaled | no q^2 in a basis forest's variance | `map: every row`; every arm of a drawn basis |
| trees-carry | a basis arm keeps its node factor | `map: every row`; `sd-carried` |
| basis-z | z on a basis | `map: every row`; `start: the chain stays` |
| none-674 | 0.674 on kind none's normal | `map: every row`; the arms cannot see 0.07 percent, and say so |
| none-q | `student(df)` on kind none moved by q | `student: no basis does not move`; `student-unmoved` |
| held-q | a held forest takes q | the trio bitwise |
| no-write-start | a write leaves the start | `write: before the first sweep`; `drawn-live` twin |
| always-start | a write moves the block whatever it holds | `write: after a sweep`; `drawn-live` |
| write-leaf | a drawn write moves the leaf scale | `write: after a sweep`; `sd-carried` |
| restate-no-write | the restate writes no prior | `write: restate`; test-forest-basis-scale-update.R lines 174 and 175 |
| restate-leaf | the restate moves the leaf scale | `write: restate` |
| v-kept | a write before the first sweep leaves v | `write: before the first sweep`, student |
| layout-slot | a new slot stays 1 | `layout: a new slot` |
| mixes | a t path at +inf | `map: every row`; the suite on NaN |
| held | drawn switched off on +inf (slice F's hold keying on a local) | `normal-reader`; prior arm |
| not-recorded | Inf not recorded on a forest with no basis, which holds silently after slice F | `normal-record`; `normal-reader` |
| basis-inf | +inf on a basis reaching the t arm | `normal: a basis keeps` |
| sharer | `forestSharers` taking it | `normal: no sharer` |
| r-sharer | `sharingForests` without normal | `normal-siblings`, its print |
| no-sd-median | the default's 2 U kept | `normal: the no-sd median`; `normal-median` |
| reader-student | the decoder builds `studentPrior(Inf)` | `normal-reader` |
| door-reader | `$getForest` or the "value in force" comparison reads Inf as a t, not through the decoder: for the door's writer, whichever lands second | forest-door's no-basis round trip |
| predicate | a site left testing `dbartsStudentPrior` alone | `normal-refusals`; `normal-siblings` |
| shadow | `normal` left out of a vocabulary list | `normal-doors`; `normal-shadow` |
| lone-taken | the lone branch student only | `normal-refusals` |
| basis-recorded | Inf on a basis forest's `normal()` | `normal-record`; `normal-basis` |
| spare-args | `normal(sd = 2)` taken | `normal-refusals` |
| example | forest.Rd's `k.scale` kept in `index-sd`'s read | `index-sd` |

### The exact gate

[`forest-prior-exact.R`](../../benchmarks/R/forest-prior-exact.R) restates every drawn arm and gains the new ones. It has every
kind already: [`numericPosterior`](../../benchmarks/R/forest-prior-exact.R) and
[`factorPosterior`](../../benchmarks/R/forest-prior-exact.R) integrate a drawn coefficient on a normal or t grid, and
[`priorFit`](../../benchmarks/R/forest-prior-exact.R) reads a size against its quantiles. The oracle's `sd` becomes the unit on a
drawn forest, the grid gains the coefficient's scale (size / unit / 0.674 on a basis, size / unit / z with no basis), and a prior
read's median follows. A patched copy has it at 125 added lines.

| arms | how many | reads |
|---|---|---|
| held on every kind, alone and beside others; unchanged | 13 | fits, splits |
| a 0/1 column and four values, drawn, at the default and stated; rows of weight zero; three forests of three trees, restated | 6 | fits, splits, within 0.5, 1, 2 |
| a basis under `student(1)` and `student(3)`, at the default and stated: four values, a 0/1 column; restated and new | 5 | the same |
| a basis under `normal()` written out, default and stated, new | 2 | the same |
| a factor at the default, `normal()`, `student(1)`, `student(3)`, stated; restated and new | 6 | the same |
| no basis under `normal()`: default, 0.7 and 3 of the unit, beside a numeric held, three forests of three trees, rows of weight zero, new | 6 | the same |
| no basis under `student(1)` and `student(3)`: the same six, new, the oracle main's | 12 | the same |
| prior: a number and a factor, at the default and under `student(1)` and `student(3)`, restated | 3 | size within the prior's quantiles |
| prior: no basis under `student(1)`, `student(3)`, and `normal()` at the default, 0.7 and 3 of the unit | 5 | the same |
| whole gate | 58 | full, 64 seeds of 50,000 kept: largest \|z\| 2.70, no failure; quick, 32 seeds of 25,000 kept: 3.49, on a `within0.5` that reads -2.14 at full |

The oracle discriminates, and the two readings of the product agree. Run on the unrestated harness against the fuller copy, 48
reads on 13 of 28 arms fail, every one a read of a coefficient alone, and every fit and split read passes, largest |z| 3.16: the
same posterior on the product, read from both sides. An oracle whose coefficient scale is 5 percent off fails 82 reads on 28 arms,
12 of them fits, |z| up to 23.

The other exact gates, quick on the fuller copy, each on exact-gates.yaml's list: bcf-exact, bcf-exact-weak, bcf-exact-restricted
and bcf-latent-exact all pass unedited, since they read a forest with no basis, a held forest or a product; a comment or two in
them that describes a drawn treatment forest's scale is restated. benchmarks/R/sbc.R's bcf arm calls `getForestAmplitudes` at four
sites and has not been run under this: it is owed before landing, by hand, with any read of a coefficient alone restated.

Probit and logistic are not gated by a new arm, as `student(df)` is not (dec-A245): the update is the one conditional, the latent
refresh reads the combined fit, and the prior enters as one precision. A mode in bcf-latent-exact.R is ~60 lines and its own
tangent measure, which the plan does not take: [Calls made](#calls-made) 6. The basis chain's agreement with main's on probit and
logistic is tier 1's.

### Baselines

| baseline | status | run on the fuller copy |
|---|---|---|
| bcf-equivalence-20eed7ed, or what slice F re-records it to | re-recorded, with a MANIFEST row and the run below | `compare --cross-host`: tier 1 on 18 of 19, worst ratio 4.3e-4, varcount identical |
| equivalence-32418b6e, or slice F's | re-recorded with a MANIFEST row, for one cell, `bart2twoforest`, whose factor forest is drawn. equivalence.R has no cross-host mode, so that cell is judged by its statistical compare, largest \|z\| 2.71 over 209 summaries, and the 16 calls' reading to rounding | `--bitwise`: 54 of 55 |
| multinomial-80b1c8d4 | no re-record | `--bitwise`: 11 of 11 |
| the four snapshot files | unchanged | no `forest(` term in them, by grep |
| bench-sampler | not re-recorded, no compare owed | the sweeps are one conditional's: the same up to rounding |

Tier 1 stands in for a bitwise compare in bcf-equivalence, where only rounding moves: `sigma` and `train` within tier 1 on
`compare --cross-host` and varcount identical. Not z mode: bcf-equivalence's statistical mode at 40 draws fails an unchanged
sampler, main with engine seed + 1 reaching |z| 2.2 to 15.3 and failing 17 of 19. The snapshot channels `mu`, `tau` and `glue` are
compared by neither, being in new units.

`numeric_swap_write_copy` is the one scenario tier 1 does not cover. It writes `sd` on a running sampler, which now changes the
coefficient's prior and not the trees', so its chain moves by design: tier 1 ratio 8e6, varcount differing in 241 of 320 cells,
its weak z bar 5.90. Over 2,000 more sweeps the two builds' fits differ by rms 0.011 in the mean over 16 seeds, each further from
a sampler created fresh at that prior, 0.025 on main and 0.031 here. That gap is on main too and unexplained: TODO
written-sd-fresh-gap. It is read before 1.0-0, by its own entry, and does not hold this item.

### The transient on a write

No gate: a pin says what a write moves, and no pin says how long the fit takes. Measured, rms about the settled mean over its
settled value, sweeps after the write, 1-5, 6-20, 21-80, 81-320, 8 seeds, the last rows of each kind 6:

| cell | main | this plan |
|---|---|---|
| gaussian, factor, 0.7 to 7 | 0.98 1.02 1.07 1.01 | 0.98 0.99 0.99 0.99 |
| gaussian, numeric, 0.7 to 7 | 1.06 1.03 1.05 1.02 | 1.00 1.01 1.00 0.99 |
| probit, factor, 0.7 to 7 | 1.00 1.05 0.99 0.98 | 0.91 0.93 0.98 1.00 |
| probit, numeric, 0.7 to 7 | 1.09 1.20 1.08 1.00 | 0.96 0.96 1.02 1.01 |
| probit, numeric, 0.7 to 1.4 | 0.98 1.01 1.00 1.01 | 0.95 0.99 0.97 1.00 |
| logistic, factor, 0.7 to 7 | 1.13 1.10 0.99 0.98 | 0.89 0.92 1.01 0.99 |
| logistic, numeric, 0.7 to 7 | 1.09 1.17 1.00 0.99 | 0.97 0.98 0.96 1.00 |
| gaussian, numeric, 7 to 0.7 | 1.20 1.02 1.02 1.02 | 0.96 0.98 1.00 1.00 |
| probit, numeric, 7 to 0.7 | 1.09 0.99 1.00 0.99 | 1.20 1.11 1.07 1.02 |
| no basis, `normal()`, gaussian, 0.7 to 7 | | 0.96 1.01 0.99 0.99 |
| no basis, `normal()`, probit, 0.7 to 7 | | 1.09 1.06 1.01 0.97 |
| no basis, `normal()`, logistic, 0.7 to 7 | | 0.94 1.00 0.94 0.98 |
| no basis, `normal()`, probit, 7 to 0.7 | | 0.95 0.94 1.03 1.02 |

A write moves no leaf, so the fit does not jump. A write up improves on main in every cell. The one cell that is worse is a
numeric basis under probit written down, 1.20 through 5 sweeps and settled by about 320, where main reads 1.09 and then 1: leaves
and a coefficient drawn at the larger size stay until the coefficient comes down. Kind none under `student(1)` writes as it does
on main, 1.00 0.99 1.00 1.01 at 0.7 to 1.4. Burn-in at a large `sd` on a forest with no basis: probit, `normal()`, sd 6, sweeps
1-5, 6-20, 21-40, 41-80, 0.90 0.98 1.06 1.02. The help says a live write settles over some sweeps.

### The rest of the gate list

By the tier: an A/B against the base's library over the channels coefficient-prior-opt-in's A/B read, on the shipped and reference
builds, none differing past rounding and only with a drawn basis forest; multinomial bitwise; the snapshot files unchanged;
tests/cpp plain and under ASan and UBSan; the R-loaded ASan run; exact-gates.yaml's list quick; sbc's bcf arm by hand, with any
read of a coefficient alone restated; the full suite once in one process at default ncpu; `R CMD check --as-cran`;
`lintr::lint_package()` against the slice's library, air, check-rc-codoc, check-win-drift, check-doc-freshness and
`mutation-battery.R verify-anchors`; the x86 leak-sanitizer leg before landing; bartCause's bcf tests against this library.

### What stays bitwise

Every call with no drawn basis forest: held forests of every kind, `student(df)` on a forest with no basis beside held ones,
multinomial 11 of 11, the four snapshot files, a lone forest. equivalence's 54 of 55. bcf-equivalence's `glue_toggle` and
`numeric_held`, no channel differing. What moves by rounding: every call with a drawn basis forest, whatever stands beside it.
What moves by design: a write of `sd` on a running drawn basis forest.

## Landing order

forests-shared-total with slice F, then this, then control-migration 1a and forest-door slice 2, setmodel-parameters' slice 4,
forest-door slice 3, coefficient-prior-live. Engine features merge serially. [Calls made](#calls-made) 3 has the alternative on
forest-door.

- forests-shared-total: read, not edited. It lands first and keeps its code and numbers. Its hold keys on `isnan(df)`, the
  statement's, not on a local of "states a df", so that a +inf df draws. It replaces its z-mode gate on bcf-equivalence with tier
  1 and varcount. Its new pins avoid `k.scale` x 0.674, a basis start of (1, 0, 1) and literal variances on a drawn basis forest,
  which this item would re-pin. Its D8 line on `k.scale` holds until this item, which edits it. Its R item 5 and `sharingForests`
  take `normal()` here, not there. Its `forestSharers` needs no change. dec-B528's refusal gains a clause if its slice 3 has
  landed; else slice 3's writer adds it.
- forest-door: reworded under dec-B542 on its branch. Its table's "`fixed()`, `student(df)`, `normal()`" row holds; its reader
  does not as forest-door-s1 has it, and takes this slice's decoder, its "value in force" comparison too. Records step for its
  writer: slice 2's item 6 at creation has moved here, its door half stays; its tests' no-basis round trip is new and is the
  door-reader mutant's killer; its table has no column for a forest stating a df and no sd, a gap `student()` shares. Its
  `validateForestAmplitude` hunk drops out.
- setmodel-parameters' slice 4: reworded under dec-B542 on its branch. `writeAmplitudePrior` takes either prior, no branch. Its
  record of the scale is still needed: `k.scale` reads the unit, not that scale. Records step for its writer: its pin of the one
  candidate gets a `normal()` row, and a pin of a bare `$setLeafPrior(normal())` on a forest stating `amplitude = normal()`, which
  dec-B540 does not name and which reads like a restatement of the coefficient's prior.
- control-migration: `dbartsForestRecord`'s `amplitude.df` holds Inf, length 1; textual meeting in `forestStatements`.
- coefficient-prior-live (dec-B511): [Open calls](#open-calls) 1. Between the normal and a t on any kind the switch writes
  `mixtureScale` and `mixtureDf`, and a t to normal resets v, with no leaf moving. Not built, so not verified.
- bartCause: no code. Its `bcf` help gets a line on what `fit$glue[, c("b.0", "b.1")]` and `fit$forests$tau` read, on its own
  branch when this lands. Its tests run against this library as a gate above. stan4bart, treatSens, bairrtt: not touched. No
  shipped header changes.

## Design notes and records

- [1. The rule](../design/forest-multiplier.md#1-the-rule) of forest-multiplier.md: the table's rows and columns for nothing, two
  levels and numeric, on what `sd` sets: the coefficient's prior where drawn, the trees where held; the start; a row for
  `normal()`; the two constants. [7. What the engine keeps general, and what is left for
  later](../design/forest-multiplier.md#7-what-the-engine-keeps-general-and-what-is-left-for-later) keeps its start item.
- multiplier-combiner.md: "Two branches, which are the two ways a magnitude can be carried" is now one per held and drawn, not per
  kind; rewritten, and the half-Cauchy lines around it, and its "a basis-free forest's own contribution is Cauchy with no sd".
- bcf.md: the mapping paragraph, and a note that dec-B18's parity holds for the prior and not for any coefficient's units: b0, b1
  and tau no longer read on the bcf package's scale.
- forest-sd-unit.md: where a drawn forest's `sd` goes. [The rule](../design/forest-defaults-by-kind.md#the-rule) of
  forest-defaults-by-kind.md if it names the opt-ins. No note of its own.
- [`forest`](../../man/forest.Rd) and the rest of Help as above. TODO forest-size-reader stays. TODO coefficient-prior-live takes
  the open call's answer. TODO written-sd-fresh-gap stays.
- A ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, TODO forest-amplitude-normal removed.

## Calls made

1. The write rule reads the block itself: exact equality with the kind's start times the old scale. A chain reaches that with
   probability zero, and nothing is added to state. Alternative: a flag of whether a sweep has run, ~15 lines and a state-format
   field, with copies and reloads to carry it. Alternative: no rule, so `$setLeafPrior` before `$run` stops matching creation, 11
   pins changing meaning.
2. The statement is a +inf df, not a new field, because it is the smallest way to state it: nothing is added at the bridge, the
   record, the reader's column, `forestSharers` or control-migration. The maintainer treats normal against t as structure, so it
   is an internal code only, its comment says so, no t arm runs on it, and the summary tells him. Alternative:
   `ForestStatement::amplitudeNormal`, a bridge vector, a record slot and a reader column, ~40 more lines in five places and three
   other plans; it reads plainer at the statement. Cost taken: `amplitudeDf` is not only a df, three-valued, with NA or Inf for
   one prior by kind, which the decoder and predicate hold in one place; a third prior that is not a t would force the field.
3. It lands after slice F, as dec-B542 says, and before forest-door's slice 2, taking that slice's item 6 at creation, ~10 lines
   and a pin. One validator serves both kinds, and writing a refusal on the basis kind only to remove it is waste. Alternative:
   land after forest-door's slice 2, no lines moved, waiting on control-migration 1a and ~16 agent-hours, slice 4's and the door's
   rows then written from reading.
4. A forest with no basis keeps its start by position: 1 standing first, as a formula puts it, and 0 where a basis forest stands
   first in a `forests` list. dec-B542 says "starts at 1", which holds only standing first. Moving it to 1 everywhere would shift
   `student(df)`'s draws there, and dec-B542 shows that shape as moving no draw. A start of 0 is not absorbing: 0.26, 0.42, 0.47
   over three sweeps. Alternative: by kind, one word in the constructor; TODO drawn-coefficient-start-by-kind stays open either
   way.
5. dec-B528 refuses a constant column naming `forest(x2)`. Its stated cost was that a sign-symmetric coefficient on a constant can
   no longer be written; this gives it a spelling, so the refusal with `sd` stated says "...; for a drawn coefficient, forest(x2,
   sd = , amplitude = normal())". Words only, ~3 lines, and the summary tells him since the words are a ruled refusal's. Where the
   call has plain terms the clause's advice meets a second refusal, which names the right spelling (`normal-column`). Alternative:
   leave the words as ruled, and a user who wanted that model does not learn it exists.
6. No latent exact mode: above. Alternative: ~60 lines and a few minutes in CI; the check it adds is one precision in code the
   arms and bcf-latent-exact's other modes already run.
7. `student(Inf)` stays refused, its message naming `normal()`, as dec-B510 ruled the rule. Alternative: take it as the normal, 2
   lines, which a reader of the engine would find natural and dec-B510 did not.
8. `normal` joins the forest vocabulary, so a caller's own function of that name called inside `forest()` or `forests =` is no
   longer reached. With `normal <- function(v) 2 * v`, `forest(x1, sd = normal(.5))` runs on main and is refused after, "forest
   'sd' must be a number, not a dbartsNormalPrior"; with `normal <- function() dbartsPriors$fixed()`, `amplitude = normal()` goes
   from held to drawn with no message. `fixed`, `student` and `invchi` already shadow this way on main, `forest()` never shipped,
   and a caller's variable or column named `normal` stays theirs. The behavior stays, is pinned (`normal-shadow`) and is said in
   the help and the summary. Alternative: resolve `normal` only at `amplitude =`, ~15 lines against `forest()`'s one vocabulary,
   with `normal` then differing from the other three.
9. On a basis forest `normal()` records NA, as forest-door has it, so its print and `$getLeafPrior(f)$amplitude` show nothing
   there and show `normal()` on a forest with no basis, where it is not the default. `$getForest` shows it on both. Alternative:
   Inf on both, a default spelled out changing the record, no gain.
10. A state taken before a write, restored after it, brings back the old start under the new prior, a legal start, so its draws
    differ from a twin's that wrote first. The state is data, and restoring it undoes no write of the prior. Alternative: have a
    restore move the block if it reads as an old start, which needs the old scale in the state; no reason to carry it.
11. `writeAmplitudePrior` replaces `setAmplitudePriorScale`, which has no other caller, rather than sitting beside it.
    Alternative: keep both, the old one dead.
12. `k.scale` keeps being the leaf scale read back, the unit on a drawn forest, and no reader of the coefficient's prior scale is
    added. `$sd` says what `sd` sets; forest.Rd's example computes its spread from `$sd / 0.674`, and the size reader stays TODO
    forest-size-reader, as the adviser's and the ledger's item. Alternative: a `$coefficient.scale` field now, ~5 lines plus a
    codoc row, when no caller has asked.
13. A new slot from a width change starts at `startScale`. No R path widens a forest today: tests/cpp only. Alternative: leave it
    1, one wrong unit nobody can reach yet.
14. `aPriorScale` keeps its name, and the half-Cauchy arm of `drawAmplitudes` stays: `student(df)` on a forest with no basis does
    not move, so both still say what they do. The earlier design's rename and removal are dropped.

## Open calls

1. A live switch between held and drawn. dec-B542 says "fixed() being chosen at creation so that nothing moves on a running
   sampler". But dec-B254 and dec-B511 rule that `$setForest(amplitude = )` changes a forest between held and drawn on a live
   sampler, and TODO coefficient-prior-live, ~400 lines and queued after forest-door, is that item. dec-B542 did not list it as
   replaced, and the adviser said such a planned switch is what would flip its recommendation. It does not hold this item, whose
   map and write rule are the same either way; it decides that item's shape.

   Under dec-B542 a switch moves the trees' scale by about `sd` over the unit, since held puts `sd` on the trees and drawn on the
   coefficient. No release has it.

   - (a) Keep the live switch. A release leaves the coefficient at its held value and writes the leaf scale, as a held forest's
     `sd` write does on main today, whose transient runs 1.02 1.07 1.06 1.05; the fit stays where it is and the trees' scale
     changes under leaves already drawn. On main that step is only the factor of 1 / 0.674; here it is about `sd` over the unit,
     so that item's plan measures it and its 400 lines may rise. Moving between `normal()` and `student(df)` moves no leaf in any
     case.
   - (b) Narrow it to normal against t, retiring held against drawn as the ledger's sentence reads. It saves part of that item,
     not sized, and a caller inside a larger sampler that holds for a burn-in and then releases has to create a new sampler.

   Recommended: (a). dec-B254 and dec-B511 were ruled for a reason nothing here takes away, the cost is a write of a kind that
   already exists, and only the ledger's one clause needs correcting. Not run: there is no switch to run yet.

## Not verified

- Slice F is not built. The no-sd median is read from its item 3; the prototype's no-sd medians are main's, 2 U, U and U under
  logistic (1.8138). Its hold predicate and whether it keys on a local are from reading. Every count of re-pins is on main's map
  before slice F.
- Sharers ran on the forests-shared-total tip without slice F, and on an earlier form of this item.
- The door, `$getForest` and slice 4 are not built, so their rows are from those plans' text.
- coefficient-prior-live is not built, so its rows are from reading.
- sbc under this is owed, not run. bartCause's test-14-bcf.R passes on this form, 184 expectations; its other files were not run.
- Item 6's layout slot and the column-of-ones pin are not built. No R path reaches the layout slot.
- The transient on a write is one data set, 6 or 8 seeds, the cells above. A write down on a factor was not run, nor a write on
  `student(df)` on a basis forest.
- Coverage of f or of a treatment effect at equal time is not measured: the kernel is the same up to rounding, and the start
  stays.
- The fork at a numeric basis under probit is one call at one seed, its fits over 8 seeds. How often a seed forks is not counted.
- Neutrality for the rest is the structure of the map; the A/B and the trio are the gate.
- No exact arm reaches probit or logistic, a drawn sigma (every arm holds it), weights other than 0 and 1, a restated sd (the twin
  only) or a forest beside sharers. No arm sees a median to better than about 0.5 percent, the literal pin's job. No new arm is
  planned for these.
- Test predictors with `normal()` were not run: the refusal is of every model of several forests.
- x86 was not run. New tests/cpp and tinytest line counts are estimates, already raised for the record of plans here running low.
