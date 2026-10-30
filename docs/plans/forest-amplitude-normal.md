# forest-amplitude-normal: a drawn normal coefficient on a forest with no basis

Status: PLANNED 2026-10-30 on ecee8bdb0 (dec-B539), not yet critiqued. One call held for the maintainer, not blocking
([Open calls](#open-calls)). Rulings: dec-B539 (the coefficient), with dec-B510 (`sd` the prior median, `student()`), dec-B511
(the coefficient's prior live, its own item), dec-B513 (`fixed()` on a lone forest), dec-B522 (`normal()` has a value to
write, no return spelling), dec-B528 (a constant basis column refused), dec-B540 (`$setLeafPrior` on the one forest that
draws).

agent: one opus implementer (effort high), one opus review that runs the mutants below, the same reviewer checking the fix
round unless the fixes change engine behavior. A blind critique of this plan first. Tier "Changes draws"
([Process by risk](README.md#process-by-risk)): a new prior reaches users, though no new sampler code does.
rng: NEUTRAL, bit for bit, for every call that runs today: the only new branch tests a stated df of +inf, which no accepted
call can state (run: nine fits, main against the prototype, `identical()`). Posterior-defining for the one new model, gated
by the exact arms below and a design note. No RNG draw is added, removed or reordered on any path.
window: for 1.0-0, after forests-shared-total's slice F. Merged once, one slice.
budget: ~490 lines in the tree, ~80 of them package code: engine ~35, R ~45, help and NEWS ~20, design notes ~15, tests/cpp
~120, tinytest ~190, exact arms ~50, records ~25. The prototype's own diff on main, with both channels below, was 24 engine
lines and 29 R lines. Agent time, estimated from comparable landings: build 4, review and fix round 3, gates 2: ~9; plans
here run 1.5 to 2 times low, so 12 to 18.

## Summary for the maintainer

You asked for a drawn normal coefficient on a forest with no basis, and the guess was a few hundred lines for a new prior.
It is not a new prior. The sampler's one coefficient update already draws a normal coefficient on any forest; a t is that
same update plus a variance redrawn after it. What is missing is a way to say so. About 35 lines of the engine, most of
them comments, and about 45 of R. The rest of the work is tests.

The model has been running all along. Today, `forest(x, basis = one, sd = v)` on a column of ones with the same tree count
and depth gives a normal coefficient on a forest that multiplies nothing, and a prototype of `amplitude = normal()` draws
identically to it. dec-B528 refuses that spelling; this item gives it a proper one.

`sd` means what it means under `student(df)`: the prior median of the coefficient's size in the response's units. With no
`sd` it takes the median `student(df)` takes there. It mixes as `student(3)` and `student(1)` do at the same median, and in
the same time. A posterior for this prior matches its exact answer on nine arms.

One choice needs you, and it does not hold the work up. The same prior can be written two ways that give the same fits and
differ in what you read back. Written the way `student(df)` is written, the coefficient carries the size. Written the way a
basis forest's normal is, the trees carry it. The plan follows `student(df)`, since you ruled `sd` means what it means there
and a later switch between the normal and the t then moves nothing but the coefficient's prior.

Nothing that runs today changes. bartCause, stan4bart, treatSens and bairrtt are not touched.

## Goal

`forest(amplitude = normal())` on a forest with no basis draws its coefficient under a normal prior, taken where
`student(df)` is taken and refused where it is refused, with `sd` its prior median and, with no `sd`, the median
`student(df)` takes there. At the live door it is the value in force only, until coefficient-prior-live. No call that runs
today moves.

## What 0.9-34 did

It had no `forest()`, no `amplitude`, no coefficient and no model of several forests. Nothing is owed. `normal(k)` existed
there only as a leaf prior, which it still is at `leaf.prior`.

## What exists today

Run on ecee8bdb0 unless marked, on private libraries built from scratch copies (scratch/fan-plan: `fan-main`, main;
`fan-proto`, main with both channels; `fan-fst`, the forests-shared-total tip c8ad05e22 with the variance channel).

- [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) has one conditional: `priorPrecision = 1 / prior.variance`, a
  row-contracted design, one standard normal per coordinate. It does not read a forest's kind. A forest of kind none has an
  all-ones basis column ([`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor), so it is the block of
  width 1 a numeric forest has. [`drawAmplitudes`](../../src/bartcore/combiner.hpp) adds the variance's inverse-gamma draw
  only where `mixtureScale > 0`, and `mixtureScale == 0` is the fixed-variance normal. No code decides that a forest of kind
  none must mix except the last arm of [`forestMapEntry`](../../src/bartcore/combiner.hpp).
- What sd means in each arm of the map, m being sd in units of s and z = 0.6744897501960817:

  | | node scale | coefficient prior | median of abs(a) | median of abs(a) x node scale |
  |---|---|---|---|---|
  | numeric, drawn | s m / 0.674 | N(0, 1) | z | 1.0007 m s |
  | none, `student(nu)` | s | t_nu, scale m / t_nu,.75 | m | m s |
  | none, this plan | s | N(0, (m / z)^2) | m | m s |

  So sd is one thing in all of them, the prior median of what the forest contributes, and the two kinds carry it in
  different places. The rows of kind none are exact; a basis forest's 0.674 is bcf's rounding, kept for its bits.
- Starts. A forest of kind none starts at 1 where it stands first and 0 where a factor or numeric forest stands first in a
  `forests` list ([`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp), TODO drawn-coefficient-start-by-kind), under
  the half-Cauchy and `student(df)` alike. A formula puts it first. Run: sigma and the fit's spread are stationary within about
  5 sweeps from a start of 1 at sd 0.05, 0.7, 5 and 30, under either channel. A start of 0 is not absorbing: 0.26, 0.42, 0.47
  over the first three sweeps.
- `k`, leaf scale. A map forest's k is 1 and its leaf scale is the node scale over the root of its tree count
  ([`mapLeafScale`](../../src/bartcore/chain.hpp)). Run: under `student(3)`, sd 0.7 on a unit of 1.3132, `getLeafPrior`
  reports `k.scale` 1.3132, the unit, and `prior.sd.of` "forest (prior median)". So dec-B540's clause "this forest's reported
  k.scale is its median" does not hold on main; see [Landing order](#landing-order).
- Restatement. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) re-resolves the entry, writes a mixture's scale through
  [`setAmplitudePriorScale`](../../src/bartcore/combiner.hpp), which is false at `mixtureScale == 0`, and re-derives the leaf
  scale. A fixed variance has no writer. [`restoreGlue`](../../src/bartcore/combiner.hpp) restores a prior variance only where
  it mixes, so a fixed one stays the destination's.
- Multiplier bookkeeping. [`forestMultiplier`](../../src/bartcore/combiner.hpp) contracts a by the row; the fit and every reader
  read a x F. Nothing else reads a forest's kind after creation but [`holdAtKindValue`](../../src/bartcore/combiner.hpp) and a
  numeric column's scale readers.
- On main `forest(amplitude = normal())` is "could not find function". `student(Inf)` is refused.

## Did we build the wrong thing, and is there a simpler form?

Yes to the second. The ruling's option was "a new prior in the engine with its own exact gate"; there is no new prior. The
normal path exists and takes any forest, so the item is a statement and one branch in the map. Three runs say so:

1. A prototype on main with `forest(x, sd = v', amplitude = normal(), n.trees = 40, power = 2, base = 0.95)` beside a factor,
   at v' = v z / 0.674, against `forest(x, basis = one, sd = v, n.trees = 40, power = 2, base = 0.95)` beside the same factor:
   `identical()` train and coefficient draws over 350 sweeps.
2. The two channels below at sd = z U, where both are node scale s and variance 1: `identical()`.
3. The exact arms pass on both channels (below).

### Two channels for one prior

- Channel V, the variance. Node scale s, as under `student(df)`; `variance = (m / z)^2`, `mixtureScale = 0`. It is the limit
  of `student(nu)`, in the same parameterization.
- Channel N, the node scale. Node scale s m / z, `variance = 1`, as a numeric basis forest. No new writer: `setForestMapSd`
  works as it stands.

They are one prior on a x F. (a, F) under V is (m a' / z, z F' / m) for (a', F') under N, the conditional of each block is
equivariant in that rescaling, and the chain differs only by where a's start of 1 falls. Run: identical at sd = z U, and each
passes the same oracle at every other sd. What differs is what is read back, run at sd 0.7 on a unit of 1.3132:

| read | V | N | `student(3)` |
|---|---|---|---|
| `getLeafPrior(1)$k.scale` | 1.3132 | 1.0378 | 1.3132 |
| size of `getForestAmplitudes`, prior median | sd / U = 0.533 | z = 0.674 | sd / U = 0.533 |
| `extract(fit, "forest")`, prior scale | U | sd / z | U |
| `$setForest(1, sd = v)` writes | the coefficient's variance | the leaf scale | the coefficient's scale |

This plan takes V: [Open calls](#open-calls) 1.

### One spelling inside the engine

The statement is `ForestStatement::amplitudeDf = +inf`, not a new field. A forest of kind none with a +inf df takes the fixed
variance; no t is drawn, no auxiliary kept, no quantile taken. The maintainer's "One doesn't sample normal by using infinite
degrees of freedom on a t" (dec-B511) is kept: this samples nothing that way, and the engine comment says so. R's
`student(Inf)` stays refused. Everything already reading "does this forest state a coefficient prior" then needs no line:
forests-shared-total's `forestSharers` and R record, slice F's items 1 and 3 (a forest of kind none draws only where
`!isnan(amplitudeDf)`, and its no-sd median is keyed the same way), `ForestCalibration::amplitudeDf`,
[`applyForestStatements`](../../src/R_interface_bartcore.cpp), the `amplitude.df` record and its column, control-migration's
`dbartsForestRecord` slot, whose validity checks only length. [Calls made](#calls-made) 2 has the alternative.

## What the user sees

Every row below ran on `fan-proto` or `fan-fst` against `student(3)` and `fixed()` beside it, and `normal()` matched
`student(3)`'s taken or refused in every one.

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

- A lone forest and a multinomial entry: "this model has one forest, which has no coefficient to give a prior; a lone forest
  draws its size through the fitting function's leaf prior, as leaf.prior = normal(k = chi(1.5, 2)) or normal(sd = invchi(1.5,
  s)) with s a size in the response's units", or forests-shared-total D10's rewording, the one branch either way.
- A sharer: "forest 1: amplitude = normal() draws one forest's coefficient; these forests share one total: state each forest's
  'sd' to give it its own".
- `normal(k = 2)`, `normal(sd = 2)`, `normal(k = chi(1.5, 2))`: "'amplitude = normal()' takes no 'k' or 'sd': a drawn
  coefficient's size is the forest's own 'sd'". Only a bare `normal()`, as forest-door has it.
- An `amplitude` of another kind, the string `"normal"` among them: "a forest's 'amplitude' must be fixed(), which holds its
  coefficient, normal() or student(df), which draw it under that prior, or left out; ...".
- `student(Inf)`: "a coefficient's student(df) takes one finite df of 1 or more; student() is student(1), the half-Cauchy, and
  normal() is the normal".
- A live change of `amplitude`: today's "$setLeafPrior restates only a forest's 'sd': 'amplitude' is fixed at creation", and at
  `$setForest` forest-door's words for every changed `amplitude`, until coefficient-prior-live.

The readers:

- `forest(sd = 2, amplitude = normal())` prints `amplitude = normal()`.
- `print(sampler)`, `print(fit)`: the short line unchanged, "1 forest1: no basis; sd 0.7".
- `summary(fit)` and the verbose creation summary: "1 forest1: no basis; sd 0.7 (stated); coefficient
  normal", beside `student(3)`'s "coefficient student(3)".
- `$getLeafPrior(f)`: `$amplitude` a `dbartsNormalPrior`; `$sd`, `$sd.stated`, `$prior.sd.of` "forest (prior median)", `$k.scale`
  all as `student(df)`.
- `getK(f)`: as `student(df)`. 1 on main's map; slice F decides what a forest with no k reports, and this follows it.
- `getForestAmplitudes(f)`, `fit$glue`: the draws, median of the size sd / U under the prior. `fit$fixed$amplitude`: no entry, as it
  has none for anything drawn.
- `extract(fit, "forest")` and `predict(type = "forest")`: each forest's own total before its coefficient, as for any forest
  (run: train less a F1 plus b F2 is the shift); `extract` has no type for coefficients.
- `$getForest(f)`: `amplitude = normal()`, as `fixed()` and `student(df)` read; the list handed back changes nothing. `$setForest(f,
  sd = v)`: writes the median, twin-identical. `$setForest(f, amplitude = normal())`: no change; `fixed()` or `student(df)` refused
  until coefficient-prior-live. `$setLeafPrior(normal(sd = v))` on the one candidate: setmodel-parameters' slice 4, under dec-B540,
  through `setForestMapSd`, which then needs no branch.
- A copy and a reload continue the live sampler within tolerance, as on main (run: 2e-14 and 1e-14); a copy against a reload exactly is a
  pin below, not run here.
- Help: [`forest`](../../man/forest.Rd)'s `amplitude` gets normal() beside `student(df)`, its median and its reach; a sentence that
  the normal is `student(df)`'s thin-tailed end with the same `sd`, and how it differs from `fixed()`. [`dbartsPriors`](../../man/dbartsPriors.Rd)'s `normal` says so too. NEWS, in 1.0-0's forest item, one clause; nothing released moves.

## Design

### Engine, ~35 lines

1. [`ForestStatement::amplitudeDf`](../../src/bartcore/combiner.hpp): a +inf df states the normal coefficient on any kind (3 lines
   of comment).
2. [`forestMapEntry`](../../src/bartcore/combiner.hpp): `normalStated = df == +inf`, `dfStated = !isnan(df) && !normalStated`.
   On kind none, drawn: `variance = (multiple / z)^2`, `mixtureScale = 0`, node factor and divisor 1, z the literal
   `normalUpperQuartile` already there. On a basis, the kind's own entry: NaN and +inf resolve alike, which forest-door's record
   does not use and a fixture may. Held: as today. -inf, below 1 and NaN as today. ~6 lines and ~8 of comment.
3. `setAmplitudePriorVariance` on [`ForestCombiner`](../../src/bartcore/combiner.hpp) and
   [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp), beside `setAmplitudePriorScale`: false where the forest mixes.
   No facade virtual; [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move. ~10 lines.
4. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp): where the entry does not mix and its variance moves, write it. A basis
   forest's never moves, so no other path reads this. ~5 lines.

Nothing else. Not `drawForestAmplitude`, `drawAmplitudes`, `restoreGlue`, `glueIsValid`, the layout, the state format, the k draw, the
rescaling step. Starts stay by position, as `student(df)`'s; `stateFormatVersion` does not move.

### Bridge

No code. [`applyForestStatements`](../../src/R_interface_bartcore.cpp) passes a +inf through, and its comment says what it means.

### R, ~45 lines

1. `normal` in [`forestConstructors`](../../R/model.R), [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R), [`processHit`](../../R/formulaTerms.R)'s
   vocabulary and [`forceCallerCode`](../../R/family.R)'s hint. Run: a caller's own variable or column named `normal` stays theirs
   in `basis =` and in a forest's predictors, as `student` does.
2. [`validateForestAmplitude`](../../R/model.R): a bare `normal` is the constructor; `normal()` kept as the object; `k` or `prior.sd`
   refused; the refusal's list.
3. [`forestStatements`](../../R/model.R) gains the kinds and records +inf for `normal()` on a forest with no basis and NA on one
   with a basis. Its caller in [`resolveSamplerSpec`](../../R/spec.R) passes `hasBasis`. A default fit's record stays as it was.
4. [`resolveForests`](../../R/model.R): the lone branch tests either class. forests-shared-total's sharer refusal and
   `sharingForests` test either, the refusal's words taking the spelling written.
5. [`reportLeafPrior`](../../R/dbarts.R): an infinite df reads `normal()`. [`printedForest`](../../R/forestPrint.R): "; coefficient normal".
6. [`studentPriorRefusal`](../../R/model.R)'s words. Where forests-shared-total's slice 3 has landed, dec-B528's refusal with `sd`
   stated gains a clause: [Calls made](#calls-made) 7.

## Slices

One slice: engine, R, help, tests, arm, design note, then one review and the gate list. Steps, each verifiable:

1. Rebase on a main carrying slice F. Re-read F as built: its hold predicate, its no-sd median, and whether either keys on a local
   that +inf would miss.
2. Engine items 1 to 4 and their tests/cpp. Run tests/cpp plain and under ASan and UBSan.
3. R items 1 to 6, `--preclean`, tinytest file. Run the A/B against the base.
4. The exact arms, quick and full.
5. Help, NEWS, design notes, lint, air, codoc, freshness.
6. Review with mutants, fix round, gates, records, CI on the branch, merge.

## Tests and gates

tests/cpp, `testForestAmplitudeNormal` beside [`testForestAmplitudeStudent`](../../tests/cpp/test_sampler.cpp), ~120 lines:

- `normal: the entry`: kind none, +inf: `sdMultiple` m, node factor and divisor 1, `variance` (m / 0.6744897501960817)^2 to the last
  place, `mixtureScale` 0, `update` true. Family by family at no sd, and a stated sd through the unit.
- `normal: a basis keeps its own entry`: +inf and NaN `sameEntry` on levels and numeric. `normal: a held forest ignores it`.
- `normal: what is refused`: -inf, 0, 0.5 and -1 throw; [`testForestAmplitudeStudent`](../../tests/cpp/test_sampler.cpp)'s `bad` list
  loses +inf.
- `normal: the same chain as a column of ones`: kind none with sd multiple z against kind numeric on a ones basis with sd multiple
  0.674, same counts: node scale s and variance 1 on both, bitwise equal draws. It survives dec-B528, which refuses at R.
- `normal: the calibration reader`: `amplitudePriorVariance` finite, `amplitudePriorScale` NaN, `amplitudeDf` +inf.
- `normal: a restated sd`: `setForestMapSd` writes the combiner's variance, the leaf scale and every leaf unmoved, then a chain made
  at that sd draws equal; the sd in force writes nothing.
- `normal: restore keeps the destination's variance`: a state of one median into a chain of another.
- `normal: no sharer`: `forestSharers` excludes it, with sd and without. `normal: the no-sd median`: equal to df 3's multiple at a held
  k, a drawn k and a named sd (slice F's three cases).

tinytest, a new test-forest-amplitude-normal.R, small n, ~190 lines:

- `normal-median`: `$sd`, `$k.scale`, `$prior.sd.of` equal to `student(3)`'s, no sd and sd stated, under gaussian, probit, logistic.
- `normal-reader`: `$amplitude` the class; print and summary lines; `forest()`'s print; `getForestAmplitudes` moves; no
  `fit$fixed$amplitude`.
- `normal-record`: `amplitude.df` Inf on the forest stating it, NA beside it, absent on a fit stating none, on a basis forest's
  `normal()` too.
- `normal-doors`: formula, `forests` list, a data object's `bases`, `bart()`, a bare `normal`, draws `identical()` across the spellings.
- `normal-siblings`: beside one basis, beside two; weights; at a stated `sd` beside sharers, sharers' parts and shares unmoved.
- `normal-refusals`: each pattern in [What the user sees](#what-the-user-sees), counted with `withCallingHandlers`; a variance
  forest; the multinomial entry.
- `normal-basis`: on a factor and a numeric forest, draws and record `identical()` to left out; or forest-door's pin if it lands first.
- `normal-twin`: a restated sd against a new sampler, seeded before creation and before the run, `identical()`.
- `normal-copy`: a copy exactly equal to a reload, the live sampler within 1e-12; a restated sd then reloaded.
- `normal-live`: `$setLeafPrior(forests = list(forest(amplitude = normal()), forest()))` refused as today.
- Three existing pins of the amplitude refusal's words reread: test-forest-amplitude-student.R's, the two prefixes in
  test-forest-arguments.R and test-bcf-creation.R stay.

Mutants, each with the gate expected to kill it:

| name | change | killed by |
|---|---|---|
| median-is-sd | `variance = m^2` | prior arm (run: 54 functionals out, z up to 658); `normal: the entry` |
| median-674 | z as 0.674 | `normal: the entry` only; the arms cannot see 0.07 percent, and say so |
| mixes | a t path at +inf | `normal: the entry`; the suite on NaN |
| held | drawn switched off on +inf (slice F's hold keying on the local) | `normal-reader`; prior arm |
| setter-dropped | `setForestMapSd` leaves the variance | `normal-twin`; `normal: a restated sd` |
| restore-writes | `restoreGlue` writes a fixed variance | `normal: restore keeps` |
| basis-inf | +inf on a basis reaching the t arm | `normal: a basis keeps` |
| sharer | `forestSharers` taking it | `normal: no sharer` |
| r-sharer | `sharingForests` without normal | `normal-siblings`, its print |
| no-sd-median | the default's 2 U kept | `normal: the no-sd median`; `normal-median` |
| reader-student | the reader builds `studentPrior(Inf)` | `normal-reader` |
| lone-taken | the lone branch student only | `normal-refusals` |
| basis-recorded | Inf on a basis forest's `normal()` | `normal-record`; `normal-basis` |
| spare-args | `normal(sd = 2)` taken | `normal-refusals` |

### The exact gate

[`forest-prior-exact.R`](../../benchmarks/R/forest-prior-exact.R), new arms, not bcf-exact.R. Its
[`numericPosterior`](../../benchmarks/R/forest-prior-exact.R) already integrates a drawn coefficient on a normal grid; the arms need
a spread on [`coefficientGrid`](../../benchmarks/R/forest-prior-exact.R), the normal in
[`numericForest`](../../benchmarks/R/forest-prior-exact.R) and the channel's F sd in
[`numericOracle`](../../benchmarks/R/forest-prior-exact.R), ~50 lines. It is on exact-gates.yaml's list already, so no workflow edit.
bcf-exact.R's modes each carry a half-Cauchy or a fixed value and a sigma quadrature, K = 2 only: a mode there is ~80 lines for less
coverage, and dec-A245 made the same choice.

Run on the prototype, two cores, whole new block with its oracles, the checks on the wider grid included:

| arm | reads | quick | full |
|---|---|---|---|
| prior: no basis under normal() at its default, at 0.7 and at 3 of the unit | size within the prior's 0.25, 0.5, 0.9 quantiles | 0.7 to 1.0 s each | 3.7 to 4.1 s each |
| no basis under normal() at its default, at 0.7 and at 3 of the unit, beside a factor | fits, splits, within 0.5, 1, 2 | 0.9 to 1.6 s each | 4.4 to 5.2 s each |
| beside a numeric column held | the same | 0.9 s | 4.5 s |
| three forests of three trees | the same | 2.5 s | 8.2 s |
| rows of weight zero | the same | 0.9 s | 3.1 s |
| whole block | | 40 s, no failure, channel V and channel N | 75 s, largest z 2.43, channel V |

The sampler's seconds exclude the oracles, about 30 s of each total. After slice F the no-sd arms read the held size, a quarter of the
range, by the oracle's existing branch. Discriminates: sd read as an sd rather than a median fails 54 functionals, z up to 658.

Probit and logistic are not gated exactly, as `student(df)` is not (dec-A245): the update is the one conditional, the latent refresh
reads the combined fit, and the prior enters as one precision. A mode in bcf-latent-exact.R is ~60 lines and its own tangent measure,
which the plan does not take: [Calls made](#calls-made) 6.

### The rest of the gate list

By the tier and the neutral class: an A/B against the base's library over the channels coefficient-prior-opt-in's A/B read, on the
shipped and reference builds, none differing; the trio bitwise; the four snapshot files unchanged; tests/cpp plain and under ASan and
UBSan; the R-loaded ASan run; exact-gates.yaml's list quick; the full suite once in one process at default ncpu; `R CMD check --as-cran`;
`lintr::lint_package()` against the slice's library, air, check-rc-codoc, check-win-drift, check-doc-freshness and
`mutation-battery.R verify-anchors`; the x86 leak-sanitizer leg before landing.

### What stays bitwise

Every call that runs today, on every baseline in benchmarks/baselines/MANIFEST as slice F leaves it: equivalence 55 of 55,
bcf-equivalence, multinomial-equivalence, the four snapshot files. Nothing is re-recorded. bench-sampler is not re-recorded and not
owed a compare: the one new branch is at creation and at a restatement, off the hot path; run: normal() and `student(3)` take
the same seconds per sweep.

## Landing order

forests-shared-total with slice F, then this, then control-migration 1a and forest-door slice 2, setmodel-parameters' slice 4,
forest-door slice 3, coefficient-prior-live. Engine features merge serially. [Calls made](#calls-made) 4 has the alternative.

- forests-shared-total: read, not edited. It lands first. Its R item 5 and `sharingForests` take `normal()` here, not there. Its
  `forestSharers` needs no change. dec-B528's refusal gains a clause if its slice 3 has landed; else slice 3's writer adds it.
- forest-door: its table's "`fixed()`, `student(df)`, `normal()`" row and reader hold as written, the door reading `amplitude.df`.
  Records step for its writer: slice 2's item 6 at creation has moved here, its door half stays; its tests' no-basis row is new. Its
  `validateForestAmplitude` hunk drops out.
- setmodel-parameters' slice 4: `setForestMapSd` writes either, no branch. Records step for its writer: its pin of the one candidate
  gets a `normal()` row. Its clause, from dec-B540, "this forest's reported k.scale being its median", does not hold on main, where
  `k.scale` reads the unit (run above) and slice F does not move it: the scale it carries on the record is needed all the more, and its
  reader's words want rereading. Told to the orchestrator; dec-B540's rule is not touched.
- control-migration: `dbartsForestRecord`'s `amplitude.df` holds Inf, length 1; textual meeting in `forestStatements`.
- coefficient-prior-live (dec-B511): channel V makes a switch on kind none a change of `variance` and `mixtureScale` alone, as dec-B511
  describes, with no leaf scale moving. Its estimate is unchanged and the setter here is its writer.
- bartCause, stan4bart, treatSens, bairrtt: not touched. No shipped header changes.

## Design notes and records

- [1. The rule](../design/forest-multiplier.md#1-the-rule) of forest-multiplier.md: its table gains the normal on kind none, and its list the one sentence on
  channel and median. [The rule](../design/forest-defaults-by-kind.md#the-rule) of forest-defaults-by-kind.md if it names the
  opt-ins. No note of its own; this adds a cell.
- A ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, TODO forest-amplitude-normal removed.

## Calls made

1. Channel V: held as [Open calls](#open-calls) 1, the plan written to it.
2. The statement is a +inf df, not a new field. It is what the prior is, the t's limit, and nothing is added at the bridge, the record, the
   reader's column, `forestSharers`, slice F or control-migration. Alternative: `ForestStatement::amplitudeNormal`, a bridge vector, a
   record slot and a reader column, ~40 more lines in five places and three other plans, each a place to forget; it reads plainer at
   the statement. Cost taken: `amplitudeDf` is not only a df, which its comment states, and the reader must turn Inf into `normal()`.
3. `normal()` is bare: `k`, `sd` and a chi prior refused, as forest-door has it. Alternative: `normal(sd = 2)` as a second spelling of
   `sd`, which dec-B508 removed from `forest()` for leaf priors.
4. It lands after slice F and before forest-door's slice 2, taking that slice's item 6 at creation, ~10 lines and a pin. It is ready as soon
   as F is, one validator serves both kinds, and writing a refusal on the basis kind only to remove it is waste. Alternative: land after
   forest-door's slice 2, no lines moved, waiting on control-migration 1a and ~16 agent-hours; slice 4's and the door's rows are then
   written from reading.
5. Starts stay by position, as `student(df)`'s. A start at 1 everywhere is one condition and moves no call that runs; it would differ
   from `student(df)` in a `forests` list that puts the plain forest second, for no gain the run could see, and TODO
   drawn-coefficient-start-by-kind holds the cleaner rule for all of them.
6. No latent exact mode: above. Alternative: ~60 lines and a few minutes in CI; the check it adds is one precision in code the arms and
   bcf-latent-exact's other modes already run.
7. dec-B528 refuses a constant column naming `forest(x2)`. Its stated cost was that a sign-symmetric coefficient on a constant can no
   longer be written; this gives it a spelling, so the refusal with `sd` stated says "...; for a drawn coefficient, forest(x2, sd = ,
   amplitude = normal())". Words only, ~3 lines. Alternative: leave the words as ruled, and a user who wanted that model does not learn it
   exists.
8. `student(Inf)` stays refused, its message naming `normal()`, as dec-B510 ruled the rule. Alternative: take it as the normal, 2 lines,
   which a reader of the engine would find natural and dec-B510 did not.
9. On a basis forest `normal()` records NA, as forest-door has it, so its print and `$getLeafPrior(f)$amplitude` show nothing there and show
   `normal()` on a forest with no basis, where it is not the default. `$getForest` shows it on both. Alternative: Inf on both, a default
   spelled out changing the record, no gain.

## Open calls

1. How the coefficient's prior is written, which does not change any fit and does change what is read back. 0.9-34 had no `forest()` and
   no coefficient. Run on the prototype at sd 0.7, unit 1.3132: draws identical at sd = z U; both pass the exact arms; with V,
   `k.scale` 1.3132, the coefficient's prior median 0.533, as `student(3)`; with N, `k.scale` 1.0378, the coefficient's prior median 0.674, as a
   basis forest's normal.
   - V, recommended and planned. The coefficient carries the size, as under `student(df)`: `sd` reads the same across the two values
     of `amplitude` on one forest, `student(1000)` converges to it with no jump in any reader, and `$setForest(sd = )` and
     coefficient-prior-live change a prior and never the scale of leaves already drawn. Cost: ~15 lines more than N (the writer), and
     a coefficient that reads on a different scale from a basis forest's normal, as `student(df)` already does.
   - N. The trees carry the size, as on a basis forest: `amplitude = normal()` means one thing on every kind. Cost: `k.scale` and
     `extract(fit, "forest")` jump between `student(df)` and `normal()` on one forest, `sd` restated rescales the leaf prior under
     leaves already drawn, and coefficient-prior-live's switch on kind none must move a leaf scale, which dec-B511's estimate did
     not count.

   The plan needs a yes or a change of letter; moving to N is a smaller diff.

## Not verified

- Slice F is not built. The no-sd median is read from its item 3; the prototype's no-sd medians are main's, 2 U, U and U under logistic
  (1.8138), `student(df)`'s there today. Its hold predicate and whether it keys on a local are from reading.
- Sharers ran on the forests-shared-total tip without slice F: the refusals and shares above are its slices 1 and 2.
- The door, `$getForest` and slice 4 are not built, so their rows are from those plans' text.
- Coverage of f or of a treatment effect at equal time is not measured. Mixing is one run: n 300, 4 seeds, 4000 kept, sd 1.2, sigma ESS
  1123 normal, 1098 `student(3)`, 1120 `student(1)`, 2225 `fixed()`; the fit's ESS 1670, 1634, 1680, 1969; 0.66 to 0.71 s. It says the
  normal is no worse than the t here and costs what the t costs against a hold, nothing more.
- Neutrality is nine fits and the structure of one branch; the A/B and the trio are the gate.
- Probit and logistic are not exact-gated, and no exact arm sees a median to better than about 0.1 percent, the literal pin's job.
- Test predictors with `normal()` were not run: the refusal is of every model of several forests. A restated sd then a reload was not run.
- x86 was not run. tests/cpp, tinytest and arm counts are estimates; plans here run low.
