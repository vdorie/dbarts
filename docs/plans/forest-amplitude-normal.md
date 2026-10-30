# forest-amplitude-normal: `sd` on the trees' scale for a drawn normal and for `student(df)` on a forest with no basis

Status: PLANNED 2026-10-30 (dec-B539), rewritten for dec-B541 on bdc9b8f34: `sd` sets the trees' scale for `normal()` and for
`student(df)` on a forest with no basis. Its blind critique (scratch/fan-crit/report.md, probes p1 to p12) read the earlier form,
and the evidence for this one is scratch/fan-plan/c-report.md with the runs in scratch/fan-plan/runs3. No call is held. Rulings:
dec-B541 (what `sd` sets), dec-B539 (the coefficient), with dec-B510 (`sd` the prior median, `student()`), dec-B511 (the
coefficient's prior live, its own item), dec-B513 (`fixed()` on a lone forest), dec-B522 (`normal()` has a value to write, no
return spelling), dec-B528 (a constant basis column refused), dec-B540 (`$setLeafPrior` on the one forest that draws).

agent: one opus implementer (effort high), one opus review that runs the mutants below, the same reviewer checking the fix round
unless the fixes change engine behavior. A blind critique of the rewrite first, since the earlier one read another design. Tier
"Changes draws" ([Process by risk](README.md#process-by-risk)).
rng: DRAW-SHIFTING, not posterior-changing, for `student(df)` on a forest with no basis, bcf()'s stated `student(1)` and the K = 2
fixtures included: the prior on a x F stays exactly, the chain moves because a starts at 1 on another scale. Bitwise for every
other call that runs today, with one corner: a caller's own function named `normal`, called inside `forest()` or `forests =`,
stops being reached ([Calls made](#calls-made) 10). Posterior-defining for the one new model, `normal()`, gated by the exact arms
below. No RNG draw is added, removed or reordered on any path.
window: for 1.0-0, after forests-shared-total's slice F. Merged once, one slice.
budget: ~850 lines in the tree, ~90 of them package code: engine ~35 touched, net about -10 (`forestMapEntry` +7 -13, measured on
a patched copy; the half-Cauchy arm in `drawAmplitudes` -6; the start +1; the rest comments and a rename), R ~55, help and NEWS
~40, design notes ~50, tests/cpp ~180 new and ~60 re-pinned, tinytest ~250 new and ~50 re-pinned, exact arms and readers ~100,
records ~25. coefficient-prior-opt-in landed 1.7 times its tests/cpp estimate and 1.4 times its tinytest estimate, which these
sizes carry. Agent time, estimated from comparable landings and the re-pins: 15 to 22.

## Summary for the maintainer

You ruled that `sd` sets the trees' scale on every forest and every coefficient option (dec-B541). So this item now does two
things: it adds `forest(amplitude = normal())` on a forest with no basis, and it moves `student(df)` there from the coefficient to
the trees. After it, doubling `sd` doubles the trees' scale everywhere, and the coefficient under `normal()` or `student(df)` has
a fixed prior.

It is not a new prior. The sampler's one coefficient update already draws a normal coefficient on any forest, and a t is that same
update plus a variance redrawn after it. The engine gets smaller: its special branch for a forest with no basis goes. Inside the
engine `normal()` is stated as a df of +inf, an internal code; no t is drawn, and users never see it. The size is in the tests and
re-pins.

The prior on what the forest adds does not change for `student(df)`, and a forest with no basis passes the exact arms on 14 t arms
and 9 normal arms. With today's start rescaled, the chain matches today's to rounding on all 19 bcf scenarios, and mixing and time
per sweep match. `normal()` at one `sd` adds about 1.48 times the standard deviation of a held forest, which the help says.

What it costs, each as you were shown or found since:

- A coefficient and a raw forest no longer read in bcf's units: `$getForestAmplitudes`, raw `extract`, bartCause's `fit$glue[,
  "a"]` and `fit$forests$mu`. The prior on their product stays bcf's, and bartCause's effects do not move.
- An oversized `sd` takes longer to burn in, most under a link, since the coefficient starts at 1 far from where it settles. A
  small version shows under gaussian at bcf's default `sd`.
- After a live change of `sd` the coefficient takes a few hundred sweeps to settle, and the fit's spread is up to a third wider
  while it does. This is one finding beyond what you were shown: the entry says the fit is unaffected, which came from a
  first-sweep check.
- Seeded draws move for `student(df)` here and for bcf(). One baseline is re-recorded and about 90 test checks re-pinned.

A forest with no basis divides by the normal's exact upper quartile and a forest with a basis keeps bcf's 0.674. That one token is
the difference left between the kinds, so the drawn rule is not fully one rule.

Adding `normal` to what `forest()` reads means a function of your own named `normal`, called inside `forest()`, is no longer
reached, as is already so for `fixed`, `student` and `invchi`. No released call is affected.

## Goal

`forest(amplitude = normal())` on a forest with no basis draws its coefficient under a normal prior, taken where `student(df)` is
taken and refused where it is refused, with no `sd` taking the median `student(df)` takes there. `sd` sets the trees' scale for it
and for `student(df)` on that forest, so `sd` is the prior median of the coefficient's size times the trees' scale on a forest
with no basis, exactly, and the coefficient starts at 1. A basis forest does not move. At the live door `amplitude` is the value
in force only, until coefficient-prior-live.

## What 0.9-34 did

It had no `forest()`, no `amplitude`, no coefficient and no model of several forests. Nothing is owed. `normal(k)` existed there
only as a leaf prior, which it still is at `leaf.prior`.

## What exists today

Run on ecee8bdb0 unless marked, on private libraries built from scratch copies (scratch/fan-plan: `fan-main`, main; `fan-proto3`,
main with this design behind switches; `fan-fst`, the forests-shared-total tip c8ad05e22).

- [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) has one conditional: `priorPrecision = 1 / prior.variance`, a
  row-contracted design, one standard normal per coordinate. It does not read a forest's kind. A forest of kind none has an
  all-ones basis column ([`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor), so it is the block of width
  1 a numeric forest has. [`drawAmplitudes`](../../src/bartcore/combiner.hpp) adds the variance's inverse-gamma draw only where
  `mixtureScale > 0`, and `mixtureScale == 0` is the fixed-variance normal. No code decides that a forest of kind none must mix
  except the last arm of [`forestMapEntry`](../../src/bartcore/combiner.hpp).
- What sets what, by doubling `sd` (scratch/decisions-evidence/fan-carries.R on main): the trees under `fixed()` on every kind,
  under the normal and `student(df)` on a basis forest, and on a sharer; the coefficient only under `student(df)` on a forest with
  no basis, and main's left-out default, which slice F removes. In each arm, m being sd in units of s and z = 0.6744897501960817:

  | | node scale | coefficient prior | median of abs(a) | median of abs(a) x node scale |
  |---|---|---|---|---|
  | numeric, drawn | s m / 0.674 | N(0, 1) | z | 1.0007 m s |
  | none, `student(nu)`, today | s | t_nu, scale m / t_nu,.75 | m | m s |
  | none, `student(nu)`, this plan | s m / z | t_nu, scale z / t_nu,.75 | z | m s |
  | none, `normal()`, this plan | s m / z | N(0, 1) | z | m s |

  So sd is one thing in all of them, the prior median of abs(a) times the node scale. It is not the median of what the forest
  contributes: a x F is a product, and p9 over 2e6 draws has the median of abs(a F) at 0.542 m s under the normal, 0.556 under
  `student(3)` and 0.587 under `student(1)`, and sd(a F) at 1 / z = 1.48 m s under the normal, against m s held. The rows of kind
  none are exact; a basis forest's 0.674 is bcf's rounding, kept for its bits.
- Starts. A forest of kind none starts at 1 where it stands first and 0 where a factor or numeric forest stands first in a
  `forests` list, its slot inheriting (1, 0, 1) from [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp); a numeric
  forest's is set to 1 in the constructor (TODO drawn-coefficient-start-by-kind). A formula puts the plain forest first. Under
  gaussian, sigma and the fit's spread are stationary within about 5 sweeps from a start of 1 at sd 0.05, 0.7, 5 and 30, at either
  place `sd` could sit. That does not hold under a link or for the node-scale form in general: see [Costs
  measured](#costs-measured). A start of 0 is not absorbing: 0.26, 0.42, 0.47 over the first three sweeps.
- `k`, leaf scale. A map forest's k is 1 and its leaf scale is the node scale over the root of its tree count
  ([`mapLeafScale`](../../src/bartcore/chain.hpp)). Today, under `student(3)`, `getLeafPrior` reports `k.scale` the unit, whatever
  sd is, and `prior.sd.of` "forest (prior median)". After this, `k.scale` is sd / z.
- Restatement. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) re-resolves the entry, writes a mixture's scale through
  [`setAmplitudePriorScale`](../../src/bartcore/combiner.hpp), which is false at `mixtureScale == 0`, and always re-derives the
  leaf scale. So it needs no change: a restated sd reaches the leaf scale and the mixture scale it writes is the one in force.
  [`restoreGlue`](../../src/bartcore/combiner.hpp) restores a prior variance only where it mixes.
- Multiplier bookkeeping. [`forestMultiplier`](../../src/bartcore/combiner.hpp) contracts a by the row; the fit and every reader
  read a x F. Nothing else reads a forest's kind after creation but [`holdAtKindValue`](../../src/bartcore/combiner.hpp) and a
  numeric column's scale readers.
- On main `forest(amplitude = normal())` is "could not find function". `student(Inf)` is refused.

## Did we build the wrong thing, and is there a simpler form?

Yes to the second. dec-B539's option was "a new prior in the engine with its own exact gate"; there is no new prior. The normal
path exists and takes any forest. Moving `student(df)` to the trees removes a branch, so the normal is the drawn rule at no df and
needs only a statement.

The prior is the one a column of ones gives, and the chain matches that column (p3, 350 sweeps, n 120, seeded before creation and
before the run): the node-scale form against `forest(x, basis = one, sd = v, n.trees = 40, power = 2, base = 0.95)` beside a
factor, at sd = v z / 0.674, `identical()` at v 1.3 and 0.5494, 2e-12 apart at 0.5 and 0.9. dec-B528 refuses that spelling; this
item gives it a proper one, and that refusal gains a clause naming it ([Calls made](#calls-made) 7).

### Same prior on student(df), other units

(a, F) today is (m a' / z, z F' / m) for (a', F') under this plan, the conditional of each block is equivariant in that rescaling,
and the chain differs only by where a's start of 1 falls and by rounding. Run on `fan-proto3` at unit 1.083, `student(3)`, sd 0.7
and 1.4:

| read | today, 0.7 / 1.4 | this plan, 0.7 / 1.4 |
|---|---|---|
| `getLeafPrior(1)$k.scale` | 1.083 / 1.083 | 1.038 / 2.076 |
| `$sd`, and `print`'s "no basis; sd 0.7" | 0.7 / 1.4 | the same |
| coefficient's prior median | 0.646 / 1.293 | 0.674 / 0.674 |
| coefficient's posterior median, abs | 0.73 / 0.75 | 0.78 / 0.38 |
| `extract(fit, "forest")`, prior scale | U | sd / z |
| `extract(contribution = TRUE)`, the fit | the same, to rounding | the same |
| `$setForest(1, sd = v)` writes | the coefficient's scale | the leaf scale |

Runs that it is the same prior: the exact arms below pass one oracle on both (on today's only the `within` reads differ, in a's
unit, max z 1.98 on the rest); and with today's start rescaled (a = z / m, v = (z / m)^2) on all 19 bcf-equivalence scenarios tier
1 passes, worst ratio 9.4e-5, varcount identical. Posterior means over 12 cells (three families, four sds, 8 seeds, 10,000 kept)
agree: max row z 2.1 to 4.8 over 200 rows, sigma within 0.0006.

### The one token

z on kind none, 0.674 on a basis: `kind == ForestKind::none ? normalUpperQuartile : 0.674`. At 0.674 on kind none the median of
abs(a) x node scale is 1.00073 m s and the chain differs (bcf-equivalence's tier 1 fails 18 of 19, ratios 7 to 650). The exact
arms cannot see it (quick max z 2.83 on either oracle, oracle gap 0.00046). So the drawn rule is one branch of the map with one
token that differs. Moving a basis forest to z would shift every bcf treatment forest's draws for no gain.

### One spelling inside the engine

The statement is `ForestStatement::amplitudeDf = +inf`, not a new field. It is an internal code for normal: a forest of kind none
with a +inf df is drawn, with no mixture; no t is drawn, no auxiliary kept, no quantile taken, and the engine comment says so. R's
`student(Inf)` stays refused. What reads "does this forest state a coefficient prior" as `!isnan(amplitudeDf)` then needs no
engine line: forests-shared-total's `forestSharers`, `ForestCalibration::amplitudeDf`,
[`applyForestStatements`](../../src/R_interface_bartcore.cpp), the `amplitude.df` column, control-migration's `dbartsForestRecord`
slot, whose validity checks only length. Slice F's items 1 and 3 need no line only if the hold keys on `isnan(df)`, on the
statement and not on a local: forests-shared-total must be told, since it is rebasing now.

What R does need: the df is a three-valued code (NaN, finite, +inf), and one prior reads NA or Inf by kind. One decoder and one
predicate hold that in one place (R item 5). `R/model.R` on the forests-shared-total tip already tests `dbartsStudentPrior` at
four sites beyond the validator, and forest-door-s1's reader still has main's `else if (!is.nan(df)) studentPrior(df)`, which
throws on Inf. A third prior that is not a t would force a real field. [Calls made](#calls-made) 2 has the alternative.

## What the user sees

Every row below ran on `fan-proto` or `fan-fst`, the earlier form of this item, against `student(3)` and `fixed()` beside it, and
`normal()` matched `student(3)`'s taken or refused in every one. What is taken and refused does not read where `sd` sits, and the
medians in force are `$sd`, which `fan-proto3` reads alike.

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
- `summary(fit)` and the verbose creation summary: "1 forest1: no basis; sd 0.7 (stated); coefficient normal", beside
  `student(3)`'s "coefficient student(3)".
- `$getLeafPrior(f)`: `$amplitude` a `dbartsNormalPrior`; `$sd`, `$sd.stated`, `$prior.sd.of` "forest (prior median)", `$k.scale`
  all as `student(df)`, whose own `$k.scale` becomes sd / z.
- `getK(f)`: as `student(df)`. 1 on main's map; slice F decides what a forest with no k reports, and this follows it.
- `getForestAmplitudes(f)`, `fit$glue`: the draws, median of the size z under the prior, under `student(df)` too, whatever `sd`
  is. `fit$fixed$amplitude`: no entry, as it has none for anything drawn.
- `extract(fit, "forest")` and `predict(type = "forest")`: each forest's own total before its coefficient, as for any forest (run:
  train less a F1 plus b F2 is the shift); `extract` has no type for coefficients.
- `$getForest(f)`: `amplitude = normal()`, as `fixed()` and `student(df)` read, through this slice's decoder, the door's "value in
  force" comparison likewise; the list handed back changes nothing. `$setForest(f, sd = v)`: writes the trees' scale,
  twin-identical. `$setForest(f, amplitude = normal())`: no change; `fixed()` or `student(df)` refused until
  coefficient-prior-live. `$setLeafPrior(normal(sd = v))` on the one candidate: setmodel-parameters' slice 4, under dec-B540,
  through `setForestMapSd`, which then needs no branch.
- A copy and a reload continue the live sampler within tolerance, as on main (run: 2e-14 and 1e-14; p6, after a restated sd,
  4e-14); a copy against a reload exactly is a pin below.
- Help: [`forest`](../../man/forest.Rd)'s `amplitude` gets `normal()` beside `student(df)`, its median and its reach; that `sd`
  sets the trees' scale and is the prior median of the coefficient's size times that scale, not of what the forest contributes,
  the existing sentence "the prior median of what the forest contributes" corrected with it; "With no basis the coefficient's size
  is the forest's size" rewritten; that at one `sd` the normal's contribution has about 1.48 times the standard deviation of a
  held forest's, and thinner tails than `student(df)`'s; that an `sd` far above the default takes longer to burn in under a link,
  and that a live change of `sd` takes some sweeps to settle. [`dbartsPriors`](../../man/dbartsPriors.Rd)'s `normal` says so too,
  and that a caller's own function of that name is not reached inside `forest()`. NEWS, in 1.0-0's forest item, one clause;
  nothing released moves.

## Design

### Engine, ~35 lines touched, net about -10

1. [`ForestStatement::amplitudeDf`](../../src/bartcore/combiner.hpp): a +inf df states the normal coefficient on any kind (3 lines
   of comment).
2. [`forestMapEntry`](../../src/bartcore/combiner.hpp): `dfStated = !isnan(df)` stays, and `tStated = dfStated && df != +inf`
   takes its place in the t arm and the validity test, so no local that reads "states a df" excludes the normal. The drawn arm
   becomes one arm for every kind, and the last arm, kind none's, goes:
   - `nodeScaleFactor = multiple`, `nodeScaleDivisor = kind == ForestKind::none ? normalUpperQuartile : 0.674`, the one token;
   - `variance` the kind's (1, or 1/2 on levels), `mixtureScale` 0;
   - where `tStated`, `mixtureDf = df`, `mixtureScale = sqrt(variance) * normalUpperQuartile / t_df,.75`, df 1 no longer special.

   A forest of kind none at NaN never reaches it: slice F holds it. Held, and -inf, below 1 and NaN on a basis, as today. +7 -13
   on a patched copy of main, with the doc comment above it rewritten.
3. [`drawAmplitudes`](../../src/bartcore/combiner.hpp): the NaN `mixtureDf` arm goes, since no drawn entry has it. The general
   expressions at nu = 1 are the same numbers. -6, and the field's comment. `mixtureDf` stays NaN on an entry with no mixture.
4. [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s constructor: a drawn forest of kind none starts at 1 at any
   position, as a numeric forest does, `kind == ForestKind::numeric` becoming `kind != ForestKind::levels`. A factor stays by
   position. TODO drawn-coefficient-start-by-kind narrows to the factor. [Calls made](#calls-made) 5.
5. The K = 2 field: [`AmplitudeSpec::aPriorScale`](../../src/bartcore/combiner.hpp), "half-Cauchy median for the mu scalar a", is
   now an sd carried by the trees. Renamed `sdControl`, bcf's own word and bcf-exact.R's, with its
   [`expandForestSpecs`](../../src/bartcore/combiner.hpp) line and the fixtures' uses, in the same commit.
6. Comments that state the old rule: the block above [`forestMapEntry`](../../src/bartcore/combiner.hpp),
   [`defaultForestSd`](../../src/bartcore/combiner.hpp)'s, and the half-Cauchy mentions on
   [`expandForestSpecs`](../../src/bartcore/combiner.hpp) and [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp).

Nothing else. No variance setter and no change to [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp), which already writes
the leaf scale. No facade virtual; [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move. Not
`drawForestAmplitude`, `restoreGlue`, `glueIsValid`, the layout, the k draw, the rescaling step. `stateFormatVersion` does not
move: a state written by a development build reads in the old units, and only shipped versions are owed compatibility.

### Bridge

No code. [`applyForestStatements`](../../src/R_interface_bartcore.cpp) passes a +inf through, and its comment says what it means.

### R, ~55 lines

1. `normal` in [`forestConstructors`](../../R/model.R), [`FOREST_ARGUMENT_VOCABULARIES`](../../R/model.R),
   [`processHit`](../../R/formulaTerms.R)'s vocabulary and [`forceCallerCode`](../../R/family.R)'s hint. Run: a caller's own
   variable or column named `normal` stays theirs in `basis =` and in a forest's predictors, as `student` does.
2. [`validateForestAmplitude`](../../R/model.R): a bare `normal` is the constructor; `normal()` kept as the object; `k` or
   `prior.sd` refused; the refusal's list.
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
   stated gains a clause: [Calls made](#calls-made) 7.

## Slices

One slice: engine, R, help, tests, arm, design note, then one review and the gate list. Steps, each verifiable:

1. Rebase on a main carrying slice F. Re-read F as built: its hold predicate and its no-sd median, each to key on `isnan(df)`.
   Grep R as built for `dbartsStudentPrior` and `amplitude.df`: slice F's R item 4 may have added sites, and each goes to the
   predicate or the decoder.
2. Engine items 1 to 3 and 5 to 6 first, with a throwaway rescaled start (a = z / m, v = (z / m)^2): judge [The one-off
   build](#the-one-off-build), then delete it and take item 4.
3. tests/cpp, new and re-pinned. Run plain and under ASan and UBSan.
4. R items 1 to 6, `--preclean`, tinytest file and re-pins. Run the A/B against the base over everything but a forest stating a
   df.
5. The exact arms, quick and full; bcf-exact-weak's read; sbc's bcf arm by hand.
6. Re-record bcf-equivalence, with its MANIFEST row.
7. Help, NEWS, design notes, lint, air, codoc, freshness.
8. Review with mutants, fix round, gates, records, CI on the branch, merge.

## Tests and gates

tests/cpp, `testForestAmplitudeNormal` beside [`testForestAmplitudeStudent`](../../tests/cpp/test_sampler.cpp), ~200 lines:

- `normal: the entry`: kind none, +inf: `sdMultiple` m, node factor m, divisor 0.6744897501960817 to the last place, `variance` 1,
  `mixtureScale` 0, `update` true. Family by family at no sd, and a stated sd through the unit.
- `student: the entry on a forest with no basis`: node as above, `mixtureDf` nu, `mixtureScale` z / t_nu,.75, at nu 1 and 3, and v
  starting at 1.
- `normal: a basis keeps its own entry`: +inf and NaN `sameEntry` on levels and numeric, divisor 0.674. `normal: a held forest
  ignores it`.
- `normal: what is refused`: -inf, 0, 0.5 and -1 throw; [`testForestAmplitudeStudent`](../../tests/cpp/test_sampler.cpp)'s `bad`
  list loses +inf.
- `normal: the same chain as a column of ones`: kind none with sd multiple z against kind numeric on a ones basis with sd multiple
  0.674, same counts, `spec.anchor` handed as 1: node scale 1 and variance 1 on both, draws equal. Exact only where (0.674 s) /
  0.674 rounds back to s, which fails for 13 percent of s, so the pin first checks that of its fixture's anchor; run, 1 passes. It
  survives dec-B528, which refuses at R.
- `normal: the calibration reader`: `amplitudePriorVariance` 1, `amplitudePriorScale` NaN, `amplitudeDf` +inf.
- `a restated sd`, normal and student: `setForestMapSd` writes the leaf scale; the coefficient, its variance or t scale and every
  leaf unmoved; then a chain made at that sd draws equal; the sd in force writes nothing.
- `a coefficient of kind none starts at 1`, standing first and standing second, normal and student; a factor's start unmoved.
- `student: the median`: median of abs(a) x node scale equals sd to the last place on kind none, by a sample-free read of the
  entry.
- `normal: no sharer`: `forestSharers` excludes it, with sd and without. `normal: the no-sd median`: equal to df 3's multiple at a
  held k, a drawn k and a named sd (slice F's three cases).

Re-pinned in tests/cpp, 47 failing checks on `fan-proto3` before slice F, by message. Slice F will have re-spelled some; re-read
first: `map statements: a forest with no basis takes the half-Cauchy median 2 where sigma is drawn and 1 where it is pinned` and
`map statements: its leaf scale is the anchor over sqrt(m)`, 9 each; `sd unit: divided once, by the unit the chain reports`, `sd
unit: a forest with a basis carries it in its leaf scale and one with none in its median`, `sd unit: a forest stating none has the
base spelling's bits`, `sd unit: a default's number stated at creation is the default's multiple, bit for bit` and `sd unit: with
the anchor handed and the unit not, a stated sd is divided by the response's own unit`, 3 each; `sd unit: a write of the sd in
force on a chain spelled in multiples states the forest and keeps its multiple` and `calibration map: the reported scale is the
map's product, whatever the basis's units`, and `calibration map: a scale-mixture forest reports its half-Cauchy median and no
variance`, 2 each; one each of `sd unit: a multiple of the anchor is taken as it is`, `student: a forest stating no df resolves to
today's entry`, `student: a forest with no basis keeps its leaf scale at the unit`, `student: a forest with no basis divides its
median by t_nu,.75`, `calibration: the BCF prognostic scale carries its own tree count`, `drawGlue draws a, then aVariance, then
b0, then b1`, `the q-variate sweep reproduces the two-scalar conditional on an indicator shape, off the same rng stream in the
same order` and `BCF grow-from-root combined fit characteristic value`. The `bad` list above is the 48th to 50th, and is this
plan's whichever form. `test_facade.cpp`'s and common.hpp's half-Cauchy fixtures reread.

tinytest, a new test-forest-amplitude-normal.R, small n, ~250 lines:

- `normal-median`: `$sd` and `$prior.sd.of` equal to `student(3)`'s, `$k.scale` equal to it and to sd / z, no sd and sd stated,
  under gaussian, probit, logistic. The coefficient's size at prior read by the arm, not here.
- `normal-reader`: `$amplitude` the class; print and summary lines; `forest()`'s print; `getForestAmplitudes` moves; no
  `fit$fixed$amplitude`.
- `normal-record`: `amplitude.df` Inf on the forest stating it, NA beside it, absent on a fit stating none, on a basis forest's
  `normal()` too.
- `normal-doors`: formula, `forests` list, a data object's `bases`, `bart()`, a bare `normal`, draws `identical()` across the
  spellings.
- `normal-siblings`: beside one basis, beside two; weights; at a stated `sd` beside sharers, sharers' parts and shares unmoved.
- `normal-refusals`: each pattern in [What the user sees](#what-the-user-sees), counted with `withCallingHandlers`; a variance
  forest; the multinomial entry.
- `normal-basis`: on a factor and a numeric forest, draws and record `identical()` to left out; or forest-door's pin if it lands
  first.
- `normal-twin`: a restated sd against a new sampler, seeded before creation and before the run, `identical()`, on one chain and
  on two.
- `normal-copy`: a copy exactly equal to a reload, the live sampler within 1e-12; a restated sd then reloaded.
- `normal-shadow`: with a caller's `normal <- function(v) 2 * v`, `forest(x1, sd = normal(.5))` refused, "forest 'sd' must be a
  number, not a dbartsNormalPrior", in a formula and in `forests =`; with a caller's `normal <- function() dbartsPriors$fixed()`,
  `amplitude = normal()` is the package's and draws; a caller's variable or column named `normal` stays theirs in `basis =` and
  the predictors.
- `normal-column`, once forests-shared-total's slice 3 has landed: `y ~ a + forest(b, basis = one, sd = 2)` refused with the
  clause; its advice as written, `y ~ a + forest(b, sd = 2, amplitude = normal())`, meets the plain-terms refusal, which names
  `forest(a) + forest(b, sd = 2, amplitude = normal())`; that is taken (p11, the two steps, on the tip before slice 3).
- `normal-live`: `$setLeafPrior(forests = list(forest(amplitude = normal()), forest()))` refused as today.
- `sd-carried`: on one data set, `k.scale` doubling with `sd` under `fixed()`, `normal()` and `student(df)` on every kind where
  each is taken, the sibling sweep dec-B541 rests on, so the exception cannot come back.
- Three existing pins of the amplitude refusal's words reread: test-forest-amplitude-student.R's, the two prefixes in
  test-forest-arguments.R and test-bcf-creation.R stay.

Re-pinned in tinytest, 45 failing checks on `fan-proto3` before slice F, which re-pins some; re-read first. Every file not listed
passed under the new map, test-forest-amplitude-student.R among them:

- test-forest-sd-unit.R, 27: `expectStated` at creation, the writer, copy, reload, `new()`, `setState`, a write then a copy, a
  second write and other data; `inForce` in its two loops; `engineEntry(sampler, 1L, "amplitude.prior.scale")` read as the median,
  now `prior.scale`; `plainRows` on two chains.
- test-bcf-family.R, 8: `priorScales(fit)` against `c(s, s / 0.674)`; `engineSlot(.., 1L, "amplitude.prior.scale")`, in its
  gaussian, probit, logistic and `declaredShape` reads.
- test-forest-leaf-prior-reader.R, 5: `engineMedian(sampler, f) * unit` against `prior$sd`.
- test-forest-basis-r5.R, 3: `priorFree$k.scale`, `priorScale(stale, 1L)`, `priorScale(recipient, 1L)`.
- test-bcf-creation.R, 1: `plainSecondPriors[[2L]]$sd` against twice its `k.scale`.
- test-multiforest-leaf-prior-writer.R, 1: `k.scale` of the two samplers.

The 46th failure on the prototype, in test-forest-arguments.R, is the vocabulary's and fails at `FAN_MODE=main` too: its word pin
is rewritten with R item 2.

Mutants, each with the gate expected to kill it:

| name | change | killed by |
|---|---|---|
| median-is-sd | divisor 1 on kind none | posterior arms (run as an oracle at divisor 1: 108 functionals out, z up to 182); `normal: the entry` |
| median-674 | z as 0.674 on kind none | `normal: the entry`, `student: the median`, and the one-off build's tier 1; the arms cannot see 0.07 percent, and say so |
| basis-z | z on a basis | `normal: a basis keeps`; the trio bitwise |
| t-on-coefficient | `student(df)` keeps today's arm on kind none | `student: the entry`; `sd-carried` |
| mixes | a t path at +inf | `normal: the entry`; the suite on NaN |
| held | drawn switched off on +inf (slice F's hold keying on a local) | `normal-reader`; prior arm |
| not-recorded | Inf not recorded on a forest with no basis, which holds silently after slice F | `normal-record`; `normal-reader` |
| start-position | kind none standing second starting at 0 | `a coefficient of kind none starts at 1` |
| restate-moves | a restated sd rescales a or the leaves | `a restated sd` |
| leaf-unwritten | `setForestMapSd` leaves the leaf scale | `normal-twin`; `a restated sd` |
| half-cauchy-arm | the NaN arm kept and reached at nu = 1 with another rate | `student: the entry`; exact arms at df 1 |
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

### The exact gate

[`forest-prior-exact.R`](../../benchmarks/R/forest-prior-exact.R), new arms, not bcf-exact.R. Its
[`numericPosterior`](../../benchmarks/R/forest-prior-exact.R) already integrates a drawn coefficient on a normal or t grid; the
arms need the normal in [`numericForest`](../../benchmarks/R/forest-prior-exact.R), and in
[`numericOracle`](../../benchmarks/R/forest-prior-exact.R) a forest with no basis at sd / z with its no-sd size 2 U. The two
existing no-basis student prior arms read a median of z in place of sd / U. ~90 lines, as scratch/fan-plan/runs3/fpe-c.R added. It
is on exact-gates.yaml's list already, so no workflow edit. bcf-exact.R's modes each carry a half-Cauchy or a fixed value and a
sigma quadrature, K = 2 only: a mode there is ~80 lines for less coverage, and dec-A245 made the same choice.

Run on `fan-proto3`, two cores, whole new block with its oracles, the checks on the wider grid included:

| arm | reads | full |
|---|---|---|
| prior: no basis under student(1) at its default; student(3) at 0.7 of the unit | size within the prior's quantiles | 2.6 s each |
| prior: no basis under normal() at its default, at 0.7 and at 3 of the unit | the same | 2.4 to 2.6 s each |
| df 1 and df 3 each: student at its default, at 0.7 and at 3 of the unit, beside a factor | fits, splits, within 0.5, 1, 2 | 3 to 3.5 s each |
| df 1 and df 3 each: beside a numeric column held; three forests of three trees; rows of weight zero | the same | 2.9 to 7.3 s each |
| normal() at its default, at 0.7 and at 3 of the unit, beside a factor; beside a numeric column held; three forests of three trees; rows of weight zero | the same | 2.9 to 6.9 s each |
| whole block, 23 arms: the 12 new student arms and 2 prior arms, the nine normal | | largest z 2.44, no failure |

The same 14 student arms against the same oracle on today's map: largest z 1.98 on fit, split and prior reads, and only the
`within` reads differ, in a's unit. That is the same prior read from both sides. Quick against a divisor of 0.674: largest z 2.83
on either oracle, oracle gap 0.00046. Discriminates: an oracle at divisor 1, sd read as an sd, fails 108 functionals, z up to 182;
the prior arm no longer reads sd at all, since the coefficient's median is z, so the arms at 0.7 and 3 and the literal pin carry
it. After slice F the no-sd arms read the median slice F gives.

Probit and logistic are not gated by a new arm, as `student(df)` is not (dec-A245): the update is the one conditional, the latent
refresh reads the combined fit, and the prior enters as one precision. A mode in bcf-latent-exact.R is ~60 lines and its own
tangent measure, which the plan does not take: [Calls made](#calls-made) 6.

The other exact gates, quick under the new map on main's fixtures: bcf-exact, bcf-exact-restricted and bcf-latent-exact pass
unedited. bcf-exact-weak fails only its `P(|a|<=0.25,0.5,1)`, a read in the old unit (E[a mu] holds, gap 0.0008), and restates its
oracle's a at z / m, ~5 lines. benchmarks/R/sbc.R's bcf arm calls `getForestAmplitudes` at four sites and has not been run under
this: it is owed before landing, by hand, with any read of a alone restated.

### Baselines

| baseline | status | run |
|---|---|---|
| bcf-equivalence-20eed7ed, or what slice F re-records it to | re-recorded, DRAW-SHIFTING, MANIFEST row citing the one-off build | below |
| equivalence-32418b6e, or slice F's | no re-record: 54 of 55 identical on main's map, and the 55th, `bart2twoforest`, is the left-out default slice F holds | run, `--bitwise`; the worst offenders there were raw `forestfit` and `glue`, the units |
| multinomial-80b1c8d4 | no re-record | run, 11 of 11 identical |
| the four snapshot files | unchanged | no `forest(` term in them, by grep |
| bench-sampler | not re-recorded, no compare owed | equal seconds per fit, cell by cell, in the mixing run below |

### The one-off build

The re-record has to be judged before it is taken. Judge it by a build with today's start rescaled, a = z / m and v = (z / m)^2 on
a drawn forest of kind none, `bcf-equivalence.R compare <baseline> --cross-host`, `sigma` and `train` through tier 1 on 19 of 19
with varcount identical. Run on main's map: 19 of 19, worst ratio 9.4e-5. That build is then deleted, about 8 lines kept in
scratch, with the start of 1 taken, and the baseline re-recorded on the result. The run's output goes in the MANIFEST row.

Not z mode. bcf-equivalence's statistical mode at 40 draws fails an unchanged sampler: main with engine seed + 1 reaches |z| 2.2
to 15.3, failing 17 of 19 (scratch/fan-plan/runs3/bcf-eq-nullseed.R). The start of 1 lands in the same range, 4.2 to 13.6, failing
18 of 19, so the mode cannot tell them apart. The snapshot channels `mu`, `tau` and `glue` are compared by neither, being in new
units.

### Costs measured

Stated in the help and recorded in forest-multiplier.md, with the numbers below. Their gates are the pins of the start and of a
restated sd; there is no gate on how long a transient lasts.

Mixing. Stationary ESS, today / this plan, 8 seeds, 10,000 kept (scratch/fan-plan/runs3/mix3.log), the kernel being the same up to
rescaling. Time per fit matched within 0.04 s in every cell, so equal sweeps is equal time.

| sd | gaussian sigma | gaussian a | probit a | logistic a |
|---|---|---|---|---|
| bcf's default | 512 / 590 | 96 / 106 | 238 / 224 | 495 / 457 |
| 0.1 | 562 / 578 | 100 / 87 | 210 / 170 | 377 / 397 |
| 1.2 | 575 / 616 | 88 / 90 | 228 / 240 | 488 / 466 |
| 6 | 609 / 629 | 89 / 93 | 211 / 200 | 502 / 472 |

Burn-in. rms about the stationary mean over its stationary value, by sweep block; `dbarts()`'s default burn is 200 and `bart()`'s
500. 8 to 24 seeds.

| cell | sweeps | today | this plan, start 1 |
|---|---|---|---|
| gaussian, bcf's default | 6-20, 21-40, 41-80 | 1.26 1.11 1.07 | 1.38 1.19 1.10 |
| probit, its default | 6-20, 21-40, 41-80 | 0.94 1.00 1.02 | 0.94 0.99 1.04 |
| gaussian, sd 6 | 6-20, 21-40, 41-80 | 1.30 1.15 1.07 | 1.46 1.24 1.16 |
| probit, sd 6 | 21-40, 41-80, 81-160, 161-320 | 0.97 1.01 1.02 0.99 | 1.63 1.66 1.31 1.07 |
| probit, sd 30 | 41-80, 81-160, 161-320, 321-640 | 1.00 1.07 1.01 1.00 | 2.88 3.21 2.43 1.31 |
| logistic, sd 30 | 6-20, 21-40 | 1.06 1.02 | 2.00 1.50 |

`normal()` at this start pays the same at probit sd 30 (2.71 at sweeps 81-160). A factor forest pays about it today (1.61 and 1.64
at sweeps 81-320), a numeric one does not. The rescaled start removes it all and is not taken, by the ruling: it makes the start
depend on `sd`, and any write before the first sweep would have to reset it, which is the `map writers: a written twin` pin.

After a restated sd (scratch/fan-plan/runs3/settle.R, settle2.R, 12 seeds, gaussian, after 1000 sweeps; blocks 1-5, 6-20, 21-80,
81-320):

| | restated | abs(a) over settled | fit's rms over settled |
|---|---|---|---|
| `student(1)`, today | 0.7 to 1.4 | 0.98 1.00 0.99 1.05 | 1.00 0.99 1.00 1.01 |
| `student(1)`, this plan | 0.7 to 1.4 | 2.05 1.91 1.40 1.07 | 1.14 1.18 1.09 1.02 |
| `student(1)`, this plan | 0.7 to 7 | 9.35 7.55 3.03 1.05 | 1.20 1.26 1.17 1.00 |
| `normal()`, this plan | 0.7 to 7 | 9.52 8.25 4.02 1.08 | 1.28 1.35 1.24 1.04 |
| `fixed()`, today | 0.7 to 7 | held | 1.02 1.07 1.06 1.05 |
| factor forest, today | 0.7 to 7 | not read | 1.03 1.01 1.02 1.02 |

So what a restated sd writes matches `fixed()` and a basis forest, the node scale under leaves already drawn, and how fast the fit
settles does not: their fits show no such rise. Settled by about sweep 320 in every row. A rescale of a at the door would remove
it, a change the ruling did not take and this plan does not build: [Not verified](#not-verified).

### The rest of the gate list

By the tier: an A/B against the base's library over the channels coefficient-prior-opt-in's A/B read, on the shipped and reference
builds, none differing but a forest stating a df on no basis, and bcf's own; equivalence and multinomial bitwise; the four
snapshot files unchanged; tests/cpp plain and under ASan and UBSan; the R-loaded ASan run; exact-gates.yaml's list quick; sbc's
bcf arm by hand; the full suite once in one process at default ncpu; `R CMD check --as-cran`; `lintr::lint_package()` against the
slice's library, air, check-rc-codoc, check-win-drift, check-doc-freshness and `mutation-battery.R verify-anchors`; the x86
leak-sanitizer leg before landing.

### What stays bitwise

Every call that does not state a df on a forest with no basis: equivalence 55 of 55 and multinomial 11 of 11 as slice F leaves
them, the four snapshot files, every held forest, every basis forest beside a held plain one. What moves: `student(df)` on a
forest with no basis; bcf-equivalence's 18 drawn scenarios, `held_plain` identical; the K = 2 fixtures. A numeric or factor forest
beside a drawn no-basis forest moves with it, since it shares the chain.

## Landing order

forests-shared-total with slice F, then this, then control-migration 1a and forest-door slice 2, setmodel-parameters' slice 4,
forest-door slice 3, coefficient-prior-live. Engine features merge serially. [Calls made](#calls-made) 4 has the alternative.

- forests-shared-total: read, not edited. It lands first and keeps its code. Its hold keys on `isnan(df)`, the statement's, not on
  a local of "states a df", so that a +inf df draws. It replaces its z-mode gate on bcf-equivalence (above) with tier 1 and
  varcount. Its item 4's K = 2 re-spelling keeps its bits there and moves here. Its R item 5 and `sharingForests` take `normal()`
  here, not there. Its `forestSharers` needs no change. dec-B528's refusal gains a clause if its slice 3 has landed; else slice
  3's writer adds it.
- forest-door: rewritten under dec-B541 on its branch. Its table's "`fixed()`, `student(df)`, `normal()`" row holds; its reader
  does not as forest-door-s1 has it, and takes this slice's decoder, its "value in force" comparison too. Records step for its
  writer: slice 2's item 6 at creation has moved here, its door half stays; its tests' no-basis round trip is new and is the
  door-reader mutant's killer; its table has no column for a forest stating a df and no sd, a gap `student()` shares. Its
  `validateForestAmplitude` hunk drops out.
- setmodel-parameters' slice 4: rewritten under dec-B541 on its branch. `setForestMapSd` writes either, no branch. Its record of
  the scale is still needed: `k.scale` is sd / z, which moves with each restatement. Records step for its writer: its pin of the
  one candidate gets a `normal()` row, and a pin of a bare `$setLeafPrior(normal())` on a forest stating `amplitude = normal()`,
  which dec-B540 does not name and which reads like a restatement of the coefficient's prior.
- control-migration: `dbartsForestRecord`'s `amplitude.df` holds Inf, length 1; textual meeting in `forestStatements`.
- coefficient-prior-live (dec-B511): a switch between the normal and the t on kind none writes `mixtureScale` and `mixtureDf` and
  a t to normal resets v to 1; no leaf scale moves (read: `prior.scale` is 1.038 and 2.076 under both at sd 0.7 and 1.4, a's
  posterior median 0.39 and 0.38). Not built, so not verified; there is no writer for a variance, which that item adds, so its
  estimate may rise slightly. A hold released on kind none now moves the node scale by 1 / z, as on a basis forest, and not by 1 /
  m.
- bartCause: no code. Its `bcf` help's "amplitude takes a half-Cauchy prior whose median is that size", and anything saying what
  `fit$glue[, "a"]` and `fit$forests$mu` read, change on its own branch when this lands. Its effects use only b and tau; its tests
  read shapes, not run under this. stan4bart, treatSens, bairrtt: not touched. No shipped header changes.

## Design notes and records

- [1. The rule](../design/forest-multiplier.md#1-the-rule) of forest-multiplier.md: the table's "nothing" row, both drawn columns,
  F at sd / z and the coefficient's fixed prior, starting at 1; a new column or row for `normal()`; the list's sentence on what
  `sd` sets and the one token; the costs above, with numbers. [7. What the engine keeps general, and what is left for
  later](../design/forest-multiplier.md#7-what-the-engine-keeps-general-and-what-is-left-for-later)'s start item narrows to the
  factor.
- multiplier-combiner.md: "Two branches, which are the two ways a magnitude can be carried" is no longer two kinds; rewritten, and
  the half-Cauchy lines around it and its "a basis-free forest's own contribution is Cauchy with no sd". bcf.md: the mapping
  paragraph naming `aPriorScale` and its "The half-Cauchy scalar a has median", and a note that dec-B18's parity holds for the
  prior and not for a's units. Both renames of `aPriorScale` there.
- [The rule](../design/forest-defaults-by-kind.md#the-rule) of forest-defaults-by-kind.md if it names the opt-ins. No note of its
  own.
- [`forest`](../../man/forest.Rd), [`dbartsPriors`](../../man/dbartsPriors.Rd) as above, TODO drawn-coefficient-start-by-kind
  narrowed.
- A ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, TODO forest-amplitude-normal removed.

## Calls made

1. The re-record is judged by a one-off build with today's start rescaled and tier 1, then the build is deleted: the only run that
   can say the chain is the same, z mode having no power here. Alternative: keep the rescaled start in the package, no transient,
   ~8 lines, 5 more pins, the start depending on `sd`; dec-B541 chose a start of 1.
2. The statement is a +inf df, not a new field, because it is the smallest way to state it: nothing is added at the bridge, the
   record, the reader's column, `forestSharers` or control-migration. The maintainer treats normal against t as structure, so it
   is an internal code only, its comment says so, no t arm runs on it, and the summary tells him. Alternative:
   `ForestStatement::amplitudeNormal`, a bridge vector, a record slot and a reader column, ~40 more lines in five places and three
   other plans; it reads plainer at the statement. Cost taken: `amplitudeDf` is not only a df, three-valued, with NA or Inf for
   one prior by kind, which the decoder and predicate hold in one place; a third prior that is not a t would force the field.
3. `normal()` is bare: `k`, `sd` and a chi prior refused, as forest-door has it. Alternative: `normal(sd = 2)` as a second
   spelling of `sd`, which dec-B508 removed from `forest()` for leaf priors.
4. It lands after slice F and before forest-door's slice 2, taking that slice's item 6 at creation, ~10 lines and a pin. It is
   ready as soon as F is, one validator serves both kinds, and writing a refusal on the basis kind only to remove it is waste.
   Alternative: land after forest-door's slice 2, no lines moved, waiting on control-migration 1a and ~16 agent-hours; slice 4's
   and the door's rows are then written from reading.
5. A drawn coefficient of kind none starts at 1 wherever it stands, one word in the constructor, not prototyped. dec-B541 rules
   "starting the coefficient at 1"; the reason it stayed by position, not differing from `student(df)`, is gone since
   `student(df)` moves with it and every such chain moves anyway. A factor keeps its start by position: moving it would shift
   models this plan does not otherwise touch. Alternative: leave position alone, a plain forest second in a `forests` list
   starting at 0, as it does today.
6. No latent exact mode: above. Alternative: ~60 lines and a few minutes in CI; the check it adds is one precision in code the
   arms and bcf-latent-exact's other modes already run.
7. dec-B528 refuses a constant column naming `forest(x2)`. Its stated cost was that a sign-symmetric coefficient on a constant can
   no longer be written; this gives it a spelling, so the refusal with `sd` stated says "...; for a drawn coefficient, forest(x2,
   sd = , amplitude = normal())". Words only, ~3 lines, and the summary tells him since the words are a ruled refusal's. Where the
   call has plain terms the clause's advice meets a second refusal, which names the right spelling (p11, `normal-column`).
   Alternative: leave the words as ruled, and a user who wanted that model does not learn it exists.
8. `student(Inf)` stays refused, its message naming `normal()`, as dec-B510 ruled the rule. Alternative: take it as the normal, 2
   lines, which a reader of the engine would find natural and dec-B510 did not.
9. On a basis forest `normal()` records NA, as forest-door has it, so its print and `$getLeafPrior(f)$amplitude` show nothing
   there and show `normal()` on a forest with no basis, where it is not the default. `$getForest` shows it on both. Alternative:
   Inf on both, a default spelled out changing the record, no gain.
10. `normal` joins the forest vocabulary, so a caller's own function of that name called inside `forest()` or `forests =` is no
    longer reached. With `normal <- function(v) 2 * v`, `forest(x1, sd = normal(.5))` runs on main and is refused after, "forest
    'sd' must be a number, not a dbartsNormalPrior"; with `normal <- function() dbartsPriors$fixed()`, `amplitude = normal()` goes
    from held to drawn with no message (p8, p12). `fixed`, `student` and `invchi` already shadow this way on main (p8), `forest()`
    never shipped, and a caller's variable or column named `normal` stays theirs (p7). The behavior stays, is pinned
    (`normal-shadow`) and is said in the help and the summary. Alternative: resolve `normal` only at `amplitude =`, ~15 lines
    against `forest()`'s one vocabulary, with `normal` then differing from the other three.
11. `aPriorScale` is renamed `sdControl` in this commit, not left stale: it is a K = 2 fixture's field, no shipped surface, 8
    sites. Alternative: leave the name and say in a comment that it now sets the trees' scale.
12. The half-Cauchy arm of `drawAmplitudes` goes: dead after this, and the general expressions give the same numbers at nu = 1.
    dec-A245 kept it for its operand order and bits, which no longer matter since the chain moves. Alternative: leave it, 6 lines
    of dead code.

## Open calls

1. RULED, dec-B541 (2026-10-30). What `sd` sets, on the trees for `normal()` and for `student(df)` on a forest with no basis, the
   coefficient starting at 1. Shown the normal's `sd` set on the coefficient as `student(df)` has it, on the trees for the normal
   alone, and on the trees for both, recommended. The maintainer: "Use your revised recommendation."

## Not verified

- Slice F is not built. The no-sd median is read from its item 3; the prototype's no-sd medians are main's, 2 U, U and U under
  logistic (1.8138), `student(df)`'s there today. Its hold predicate and whether it keys on a local are from reading. Every count
  of re-pins is on main's map before slice F, which removes the left-out default many of them read.
- Sharers ran on the forests-shared-total tip without slice F: the refusals and shares above are its slices 1 and 2, and on an
  earlier form.
- The door, `$getForest` and slice 4 are not built, so their rows are from those plans' text.
- The live normal to t switch of coefficient-prior-live is not built and not verified: read from `prior.scale` and a's posterior
  median at the two priors.
- sbc under this is owed, not run. bartCause's tests under this were not run.
- Item 4's start at any position was not prototyped: the runs start a plain forest first, at 1, as a formula does.
- The settling after a restated sd is one gaussian data set, 12 seeds. No rescale at the door was tried. Whether a link makes it
  longer was not run.
- Coverage of f or of a treatment effect at equal time is not measured, and dec-B541 leaves it there: the kernel is the same up to
  rescaling, and only the start differs.
- Neutrality for the rest is the structure of the map; the A/B and the trio are the gate.
- No exact arm reaches probit or logistic, a drawn sigma (every arm holds it), weights other than 0 and 1, a restated sd (the twin
  only) or a forest beside sharers. No arm sees a median to better than about 0.5 percent, the literal pin's job. No new arm is
  planned for these.
- Test predictors with `normal()` were not run: the refusal is of every model of several forests.
- x86 was not run. New tests/cpp and tinytest line counts are estimates, already raised for the record of plans here running low.
- The tests/cpp column-of-ones pin has not been written; that an anchor of 1 rounds back is arithmetic run in R.
