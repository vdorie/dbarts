# forest-sd-unit: a forest's sd is stated and reported in the response's units, and the engine is handed the statement

Status: PLANNED (dec-B246 as revised, dec-B253, dec-B266; dec-B269 as dec-B280 and dec-B282 revise it).
Follows [forest-defaults-by-kind.md](forest-defaults-by-kind.md) and push 3 of
[written-surface.md](written-surface.md), neither of which has landed. Amended 2026-10-07 after the
critique of the multiplier law: this slice takes in what was that law's first push (the engine is handed
what was stated and holds the rule), gives the reader and `extract` the shape they keep, and takes a
named single sd (dec-B280). Amended again 2026-10-07 for dec-B296: the unit is taken over the rows in
the likelihood when the sampler is created, by slice B of [response-scale-rows.md](response-scale-rows.md),
which lands first; the Context bullet on L, rule 1, step 3's "L written out", row 13 of the pair
script and the mutation on the unit's rows are restated below.

agent: one push. Opus implementer for the engine, the bridge and the R code; sonnet for the respelled tests,
the benchmark scripts and the help once the code is fixed; opus reviewer told to refute, and to read the
diff twice: once for the change of shape, which must move no bit, and once for the division, which must
happen exactly once. The reason for opus: every slip is a prior off by a factor of sd(y) with no message.
The conversion is one division, and it can be made twice, or not at all, on any one of six paths
(creation at four doors, the writer, a re-creation); the review of written-surface's second push found
three such one-path faults in a slice that touched no arithmetic at all.
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- NEUTRAL for every fit that states no forest `sd`, and for every fit under probit, where the unit is 1.
  That covers the change of shape too: the engine resolves the law R resolved, from four numbers a forest
  where it was handed eight. Measured for the unit on a stand-in: 5 of 5 such pairs identical, and no
  seeded draw of the suite moved. Not measured for the change of shape, none of which is built: its
  gates are the three bitwise compares, the snapshot files, the pair script and the restated tests/cpp
  fixtures.
- POSTERIOR-CHANGING for one sequence: a model of several forests that states `sd` on a forest, at
  creation or through `$setLeafPrior(forests = )`, under a gaussian or a logistic response, and is left as
  written. The same number now states another prior: `sd = 0.7` was 0.7 standard deviations of the
  response, and is 0.7 units of it. The oracle is the same model on the base build with the number divided
  by the unit.
- SHIFTING for that sequence respelled into the response's units (`sd = 0.7 * sd(y)`): the prior is the
  same to rounding and the draws are not always the same bits. Measured on the stand-in, 18 gaussian and
  logistic pairs of 20 sweeps: 8 identical, 10 apart by 1e-14 to 1e-12.
- No scenario of the three equivalence baselines and no snapshot file states a forest `sd` (searched), so
  nothing is re-recorded. Four exact gates and `sbc.R` state one and are respelled.
- Accepted where refused: `sd` on the one forest of a single-forest model; a named single `sd` on a
  forest and in the writer, its name dropped.
- Not a released object: a sampler saved by an earlier build of this branch holds eight numbers a forest
  and is refused by the bridge's check of the record's length.
window: pre-release, after forest-defaults-by-kind and before the kind by class and the multiplier law, so
that every later test, pinned string and help page states an sd once, in the unit it keeps, and reads it
back in the shape it keeps. Serial with any other work in [`forestParams`](../../R/model.R), the
multi-forest block of [`resolveSamplerSpec`](../../R/spec.R), [`reportLeafPrior`](../../R/dbarts.R),
[`writeForestSpreads`](../../R/dbarts.R), [`extractParameter`](../../R/generics.R),
[`applyForestAttributes`](../../src/R_interface_bartcore.cpp),
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp), [`ForestSpec`](../../src/bartcore/combiner.hpp),
[`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) or the K-forest constructor of
[`Chain`](../../src/bartcore/chain.hpp); so serial with [leaf-conversions.md](leaf-conversions.md).
bartCause's edit lands the same day.
budget: ~2250 lines changed (engine ~300; bridge ~160; tests/cpp ~480; R ~300; tinytest ~620, of which
~450 in one new file and ~170 changed in about eight; benchmarks ~50; help ~130; design note,
architecture, public-surface, TODO and the two indexes ~210), upper figure 3900 (engine 500, bridge 270,
tests/cpp 800, R 520, tinytest 1100, benchmarks 90, help 250, records 370). As first planned the slice
was 1300 and the law's first push 1250; folded, less the kind's route (about 300, now the law's held
push's) and the refusals of development-build objects (about 120). The design estimated 800 and planned
for 1500; written-surface's two landed pushes ran at 1.8 and 2.7 times their plans.

## Goal

A forest's `sd` is a number of the response's units wherever it is written or read: `forest(sd = )` at
creation, `$setLeafPrior(forests = )`, what `$getLeafPrior()` returns, and
`extract(type = "leaf.prior.sd")`. The units are those of y for a continuous response and of the latent
index under probit and logistic, where the link's own error has standard deviation 1 and pi / sqrt(3). A
default does not move: it is the same multiple of the response's standard deviation as before, and is
reported as the number of response units that comes to. The engine is handed what the caller stated for
each forest (an sd or none, whether its coefficient is held, whether it has a basis) and holds the rule
and every default itself; R derives nothing. The reader returns the statement, with the number in force
and whether it was stated beside it. With one unit at every place, the one forest of a single-forest
model takes `sd` too, as its `leaf.prior = normal(sd = )`; and a single `sd` that carries a name is
taken, the name dropped.

## Context

Measured at the tip (4f2e79e6; written-surface pushes 1 and 2 landed, push 3 and forest-defaults-by-kind
not) on the shipped build, R 4.6.1, unless a line says otherwise. The spellings below are the tip's, with
a tilde on a basis written as code. L is the response's scale as the engine holds it; "size" is the prior
median of the absolute multiplier times the prior standard deviation of the forest's own fit, exact where
the coefficient is held.

- The suite. At the tip 231 files, of which 4 exit off a reference build and 3 ask for more than two
  threads; the other 224 give 16065 results, none failing, 102 seconds in one process. On push 3's build
  (4eeaf03d, traced 2026-10-07): 233 files, 230 without the three, 17121 results, none failing.
- L. Under gaussian it is the sample standard deviation of the response net of its offset, n - 1 in the
  divisor, unweighted, over the rows the sampler holds: 3.26736095623188 from the engine against R's
  3.26736095623189 for `sd(y)`, one unit in the last place apart; `sd(y - offset)` with an offset;
  `sd(y[subset])` under `subset`. Measured before response-scale-rows: unchanged by weights, rows of
  weight zero counting. Once that slice has landed it is over the rows of positive weight when the
  sampler is first created, still unweighted, and held when the weights or a mask change later; where
  those rows hold fewer than two distinct values it is over every row (dec-B296, dec-B302). Under probit it is
  1 and under logistic pi / sqrt(3), 1.8138. It is computed at a first creation
  ([`scaledResponseSd`](../../src/bartcore/chain.hpp), [`latentScaleAnchor`](../../src/bartcore/chain.hpp)),
  recorded by R as the forests' `anchor` and handed to every re-creation
  ([`AmplitudeSpec::anchor`](../../src/bartcore/combiner.hpp)), so it is the same number after
  `setResponse(10 * y + 5, updateScale = FALSE)` and on a copy; `updateScale = TRUE` is refused on such a
  sampler.
- What a stated or default number s is the size of today. Engine prior draws under a flat likelihood,
  5000 sweeps a row, Monte Carlo error about 0.02; each entry is the measured size over s L.

  | forest | what was measured | gaussian | probit | logistic |
  |---|---|---|---|---|
  | no basis, drawn, default (s = 2; 1 under a latent family) | a F | 0.996 | 0.977 | 1.002 |
  | no basis, drawn, `sd = 0.7` | a F | 0.967 | 0.987 | 0.955 |
  | no basis, held, default | F | 0.504 | 1.002 | 1.002 |
  | no basis, held, `sd = 0.7` | F | 1.423 | 1.415 | 1.430 |
  | factor, drawn, default (s = sqrt(2 / K): 1 at K = 2) | (a2 - a1) F | 1.014 | 1.012 | 1.023 |
  | factor, drawn, `sd = 0.7`; a logical | the same | 0.948; 1.018 | 1.006 | 0.990 |
  | factor, held, `sd = 0.7`; default | (a2 - a1) F | 1.483; 1.494 | 1.494; 1.498 | 1.474; 1.477 |
  | one numeric column w, drawn, default; `sd = 0.7` | a F, times the median nonzero abs(w) | 0.713; 0.691 | 0.701; 0.712 | 0.736; 0.693 |
  | the same column times 1000; times 0.001; an age (mean 50) | the same | 0.689; 0.665; 0.711 | | |
  | a 0/1 number at 50, 20, 5 percent treated | a F | 0.719; 0.700; 0.705 | | |
  | one numeric column, held | | refused (dec-A171) | refused | refused |

  So s is in units of L. It is the size of the forest's contribution where the forest has no basis or a
  factor and its coefficient is drawn. It is not where the coefficient is held (a held forest with no
  basis has size L whatever s is; a held factor 1.48 s L) or the basis is numeric (0.71 s L at a row whose
  basis has the median nonzero norm). Those are the multiplier law's to change; this slice changes the
  unit of s and nothing in that table but its last factor. (Rows for a factor of three levels and for
  two numeric columns were measured and agree; both shapes are refused before the release, dec-B281 and
  dec-B282, and are left out.)
- Who resolves the law today. R does: [`forestParams`](../../R/model.R) turns each forest's statement
  into eight numbers (tree count, base, power, a leaf-scale factor, its divisor, the coefficients'
  variance, a half-Cauchy median, the hold), the bridge passes them
  ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)), and
  [`ForestSpec`](../../src/bartcore/combiner.hpp) carries the last five as `nodeScaleFactor`,
  `nodeScaleDivisor`, `amplitudePriorVariance`, `amplitudePriorScale` and the hold flag. The K-forest
  constructor of [`Chain`](../../src/bartcore/chain.hpp) then divides by
  [`basisRowNorm`](../../src/bartcore/chain.hpp) ([`mapLeafScale`](../../src/bartcore/chain.hpp)). The
  defaults live in R ([`defaultAmplitudePriorScale`](../../R/model.R), and `sqrt(2 / K)` written in
  `forestParams`). Of a forest's kind the rule reads one thing, whether it has a basis.
- What reads and writes the number.

  | place | today | unit |
  |---|---|---|
  | `forest(sd = s)` at creation; `$setLeafPrior(forests = list(forest(sd = s)))` | s goes to the engine as given ([`forestParams`](../../R/model.R), [`writeForestSpreads`](../../R/dbarts.R)) | of L |
  | `$getLeafPrior(f)$leaf.prior` | `forest(sd = s)`, s the engine's own number, stated or not | of L |
  | `$getLeafPrior(f)$k.scale` | the forest's own prior standard deviation, before its coefficient: L for a forest with no basis, s L / 0.674 over the row norm for one with a basis | the response's |
  | `amplitude.prior.scale`, `leaf.scale.factor` entries | s again | of L |
  | `$getLeafPrior()` with no forest | a list with no names (read off push 3's build) | |
  | `extract(fit, type = "leaf.prior.sd")` | `k.scale` over k, a vector named `forest1`, `forest2`: 1.7408 and 2.5828 for a default two-forest fit whose sd(y) is 1.7408; with `forest = 2` the bare number ([`selectForests`](../../R/generics.R)) | the response's, of another quantity |
  | `print` of a fit, `show` of a sampler, the verbose summary | no forest's sd is printed | |
  | `leaf.prior = normal(sd = s)` on a single forest | the prior standard deviation of the forest's total, exactly | the response's |

  The reader and `extract` already disagree: for one default factor forest the reader says 1 and `extract`
  2.58. And the reader cannot say whether a number was stated.
- One forest. `y ~ forest(x1 + x2, sd = 2)` and `forests = list(forest(sd = 2))` are refused at `dbarts`
  and `bart` ("this model has one forest, so its size is the fitting function's leaf.prior = normal(sd = ),
  not forest(sd = )"). `leaf.prior = normal(sd = 2)` beside several forests is refused ("a multi-forest
  model does not support a named leaf-prior 'sd': the leaf prior's 'sd' is not a forest's 'sd'");
  `normal(k = 3)` beside them is refused and a bare `normal` accepted.
- A named single sd. `forest(basis = dose, sd = pars["s"])` is refused at creation and in
  `$setLeafPrior(forests = )` ("forest 'sd' must not be named"), on a forest with a basis and on one
  with none; `leaf.prior = normal(sd = pars["s"])` on a single forest is taken (read off push 3's build).
- A response with no spread. Two forests on a constant response are created with a warning; L is 0 and
  every forest's prior standard deviation is 0, with `sd = 0.7` stated or not.
- Where a forest's sd is stated, counted by running the suite on a build that logs it: of 776 models of
  several forests created in 53 test files, 43 forests state one, in 6 files (test-bcf-family.R 14,
  test-bcf-creation.R 8, test-multiforest-leaf-prior-writer.R 8, test-forest-basis-r5.R 7,
  test-forest-arguments.R 5, test-formula-terms.R 1): 20 under gaussian, 6 under logistic, 17 under
  probit. On push 3's build (1186 models): 51 forests state one, 20 gaussian, 10 logistic, 21 probit.
  `$setLeafPrior(forests = )` writes one 33 times in 3 files, 28 of them under gaussian. In
  benchmarks/R the four BCF exact gates and `sbc.R` state `sd = sdControl` and `sd = sdModerate`, 10
  lines; `bcf-equivalence.R`, `equivalence.R` and `multinomial-equivalence.R` state none. No vignette and
  no help example states one.
- What moves in the suite. For the unit, measured on a stand-in (the conversion done in R against
  `sd(y - offset)`): no file stops and 15 assertions fail in 7 files: pins of the stored numbers behind a
  stated sd (5), of the reader's `leaf.prior`, `leaf.scale.factor` and `amplitude.prior.scale` (5), of
  `k.scale` against a stated multiple (2) and of `extract` (3). For the change of shape, measured on the
  multiplier law's stand-in: 72 assertions pin the eight stored numbers by position (test-bcf-family.R
  46, test-forest-arguments.R 16, test-bcf-creation.R 8, test-multiforest-leaf-prior-writer.R 2). The two
  sets overlap in the first five. In tests/cpp the fixtures of
  [`testForestMapWriters`](../../tests/cpp/test_sampler.cpp),
  [`testBCFCalibrationMap`](../../tests/cpp/test_sampler.cpp) and
  [`testForestCalibration`](../../tests/cpp/test_sampler.cpp), and about 100 other lines, set the four
  derived numbers.
- The same model on both builds, the stand-in's side respelled as s L: with no sd stated, and under
  probit, identical; of 18 gaussian and logistic rows (each kind of forest, an offset, `subset`, weights,
  three forests, a formula term, the writer, a copy, a value read and written back, `bart`) 8 identical
  and 10 apart by at most 9.6e-13 after 20 sweeps. Left as written, `sd = 0.7` on both forests gives draws
  2.31 apart under gaussian and 3.63 under logistic.
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) passes `sd = sd.control`, `NULL` by default, and
  `sd = sd.moderate`, 1 by default, and two of its test files build the same sampler by hand with
  `sd = 1` on the treatment forest. At K = 2 the default of a forest with a basis is exactly 1, so
  `sd = 1` and no sd are one fit today. It reads `getLeafPrior(1L)$response.scale` and
  `$response.shift` and no other entry. stan4bart (bartcore, a9d081b), treatSens (dbarts-1.0, aecec71)
  and bairrtt (main, 3f57f61) state no forest sd and read no leaf prior of a model of several forests
  (searched again 2026-10-07).

## The rule

1. The unit. A forest's `sd` is in the response's units: the units of y under a gaussian response, and of
   the latent index under probit and logistic. One number of the sampler, the unit, converts it: L in the
   response's units, which under gaussian is `sd(y - offset)` over the rows in the likelihood when the
   sampler is created (the rows `subset` and the na.action keep, less those of weight zero; dec-B296),
   1 under probit and pi / sqrt(3) under logistic.
2. Stated. The engine is handed the number as written and divides it by the unit; what it then does with
   the quotient is what the tip does with a number in units of L. Nothing else of the law changes.
3. Not stated. The engine is told so and takes the default itself: a multiple of L (2 or 1 for a forest
   with no basis, sqrt(2 / K) for one with a basis), which never passes through the unit on the way in.
4. Who holds the law. The engine, in one function: a forest's statement (an sd or none, a hold, whether
   it has a basis, and a coefficient variance where `amplitude.prior.variance` states one) is resolved
   there, once for each construction, into the numbers the tip's eight carried. R derives none of them.
5. Read. `$getLeafPrior(f)` returns `leaf.prior`, the statement: `forest(sd = v)` where an sd was
   stated, v the number as written, and `forest()` where none was; `sd`, the number in force in the
   response's units, which is v itself, bit for bit, or the default's multiple times the unit; and
   `sd.stated`, `TRUE` or `FALSE`. With no forest it returns the list of those entries, named by the
   forests' labels.
6. Write. `$setLeafPrior(forests = )` with `forest(sd = v)` states v: the forest is a stated one from
   then on, also where v is the number in force, in which case no bit of the prior moves. An entry
   `forest()` leaves its forest as it is, so the reader's `leaf.prior` goes back as it came and changes
   nothing, to the bit; a caller who wants the default's number stated writes `forest(sd = p$sd)`.
7. `extract(type = "leaf.prior.sd")` on a fit of several forests returns the reader's `sd` entries: a
   list by forest, named by label, for several forests, and the entry itself for one, as `extract` gives
   every other quantity selected by `forest`.
8. The unit is derived once, by the engine, at a first creation, and recorded beside the anchor. Every
   re-creation is handed the record: a copy, a reload, a sampler made again from its own control, model
   and data. No later change of the response moves it, as none moves the anchor. A control or a `forests`
   list carried to other data states the same numbers of that data's units: the record is not carried,
   and the new sampler derives its own.
9. One forest. `forest(sd = s)` on the one forest of a single-forest model is that model's
   `leaf.prior = normal(sd = s)`, in every family, with that family's own verdict on a stated sd. Stated
   at both places it is refused by name.
10. A name on a single `sd` is dropped, at `forest()` and in the writer, on every forest (dec-B280 as
    dec-B282 leaves it). An `sd` of length greater than one stays refused, named or not.

Before and after, per door. "Several" is a model of several forests.

| door | before | after |
|---|---|---|
| a `forest()` term of a formula at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()`; a list of sds over a data object's `bases` | s times L | s response units; one model at every door, by one conversion in the engine |
| `$setLeafPrior(forests = )` | s times L | s response units; the forest is stated |
| `$getLeafPrior()`'s `leaf.prior` | `forest(sd = )` in units of L, stated or not | the statement: `forest(sd = )` in the response's units, or `forest()` |
| its `sd`, `sd.stated` | absent | the number in force, in the response's units; whether it was stated |
| its `k.scale`, `response.scale`, `response.shift` | the response's units | unchanged |
| its `amplitude.prior.scale`, `leaf.scale.factor`, `leaf.scale.divisor`, `amplitude.prior.variance`, `basis.row.norm`, `prior.sd.of` | the engine's channels, in units of L | unchanged, and said to be so; the multiplier law removes or renames them |
| `$getLeafPrior()` with no forest | a list with no names | named by label |
| `extract(type = "leaf.prior.sd")`, several | a vector `forest1`, `forest2` of each forest's own standard deviation before its coefficient | a list by label of each forest's `sd`; for one forest the number |
| the same on one forest; `extract(type = "k")` | | unchanged |
| `print`, `show`, the verbose summary | no forest's sd | unchanged; the printed block is the multiplier law's |
| `$setForestBasis` | the forest's own scale is derived again from the new block's row norm, under the number in force | unchanged, the number in force being the stated one or the default |
| `predict`, `fitted`, a state stored or installed | do not read a prior | unchanged |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | built from the eight numbers and the anchor | built from the statement, the anchor and the unit: the same prior, bit for bit |
| `sd` on the one forest of a single-forest model | refused | the fit with `leaf.prior = normal(sd = )` |
| the fitting function's `normal(sd = )` beside several forests | refused | refused, in words that say where the size is written |
| `forest(sd = pars["s"])`, a single number with a name | refused | taken as the number |

## Refused forms, with their texts

Base R's style, as in written-surface.

    this model has one forest and states its size twice: 'sd' on the forest and 'leaf.prior' given to the fitting function are one statement; give one            ('k', and the retired 'node.prior', named as written)
    a forest's 'sd' is in the response's units, and this response has no spread to state one against (its standard deviation is 0); leave 'sd' out, or give a response that varies
    a model of several forests states each forest's size on the forest: write the size of the forest with no basis as forest(x1 + x2, sd = 2), which is in the response's units as leaf.prior = normal(sd = ) is, and leave 'sd' out of 'leaf.prior'
    forest 'sd' must be a single number, not a vector of length 2: a forest states one sd; to size a column differently, rescale it in 'basis', as I(dose / 30)

Behind them, from the bridge: "the forests' unit must be a positive finite number", beside the anchor's
text; "the forests' 'sd' must hold one number, or NA, for each forest"; and the check of the record's
length, which reads four where it read eight and is what a record of another shape meets. Gone: "this
model has one forest, so its size is the fitting function's leaf.prior = normal(sd = ), not
forest(sd = )"; ["must not be named"](../../inst/tinytest/test-forest-arguments.R), with its hint to
`unname()`. Replaced by the third: "a multi-forest model does not support a named leaf-prior 'sd'", whose
reason, that the two sds are different things, stops being true. The fourth is today's text for a longer
`sd` without "for every column of its basis". Kept as they are: every other text of
[`validateForestSd`](../../R/model.R).

## Constraints

- A fit that states no forest sd draws what it draws on the base build, to the bit, in every family. No
  baseline is re-recorded and no snapshot file regenerated.
- The law is not touched: which channel a forest's number reaches, the 0.674, the coefficient's variance,
  the row norm, the held values and the starts by position, the defaults and the K in them. It moves
  from R to the engine with every expression in its order of operations. `amplitude.prior.variance`
  stays an argument.
- The law is resolved in one function of the engine and nowhere in R; a forest's law is resolved once
  for each construction, by the chain, and the combiner is handed the result.
- The engine learns nothing of a forest's kind here beyond whether it has a basis, which is all the
  tip's rule reads. The kind's route from the data object is the multiplier law's held push's.
- The conversion is made in one place, the engine, against one recorded number. R divides nothing and
  multiplies nothing: it hands over what was written and reports what the engine returns.
- A stated number is kept as written, so that the reader returns it bit for bit.
- The reader's `sd` is one unnamed number a forest in this slice. A numeric forest's takes its column's
  name with the multiplier law, when it becomes that many units per unit of the column.
- The forests' record holds four numbers a forest (tree count, base, power, the hold), `sd` and, where
  one is stated, the coefficient variance, each one number or NA a forest, the anchor and the unit.
- The flat C header, the stored state and every state block are untouched.
- A `forests` list written before the slice is still accepted, with every number read in the new unit.
- No test this slice adds fits a factor of three or more levels or a basis of several numeric columns
  (dec-B281, dec-B282: refused one slice on).
- Each commit leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the
reader's sake. Calls are in push 3's spelling of a basis; the fixture is 150 rows with x1, x2, a 0/1 z, a
factor zf, dose, an offset column and weights of which ten are 0.

1. The engine: the statement, the law and the unit. [`ForestSpec`](../../src/bartcore/combiner.hpp)
   loses its four derived numbers and states: `sd`, in the response's units, not a number where none is
   stated; the hold it has; and the coefficient variance, not a number where none is stated (the
   multiplier law's held push removes it). For fixtures and for
   [`expandForestSpecs`](../../src/bartcore/combiner.hpp), whose two-forest spelling states BCF's sizes
   as multiples of L, it may state a multiple of L in place of `sd`; nothing R reaches sets that.
   [`AmplitudeSpec`](../../src/bartcore/combiner.hpp) gains `unit`, not a number where the engine is to
   derive it. One function (`forestLaw`) turns a statement, whether the forest has a basis, K, the
   family, the unit and the row norm into the forest's leaf-scale factor and divisor, its coefficients'
   variance and its half-Cauchy median, and one (`defaultForestSd`) holds the defaults; both give the
   tip's numbers, the row norm and the starts by position included, with every expression of
   [`mapLeafScale`](../../src/bartcore/chain.hpp) in its order of operations. The K-forest constructor
   of [`Chain`](../../src/bartcore/chain.hpp) takes the unit it is handed, or derives it beside the
   anchor (under gaussian the anchor times the scale of the response transform at that moment, under
   probit and logistic the anchor), resolves each forest once and hands
   [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) the result. A forest that states an sd
   has `sd / unit` where its number in units of L stood: the half-Cauchy median of a forest with no
   basis, the leaf-scale factor of one with a basis. A stated sd with a unit that is not positive and
   finite, and an sd that is not, throw. [`ForestCalibration`](../../src/bartcore/chain.hpp) gains the
   unit, the sd in force in the response's units (the stated number itself, or the multiple of L times
   the unit) and whether it was stated. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) takes
   the response's units and marks the forest stated; where the number is the sd in force as reported it
   keeps the quotient it has, so no bit moves, and otherwise divides as creation does. Every struct
   here is read by objects that do not track headers: `--preclean`.
   Tests, tests/cpp:
   - The fixtures that set the four numbers are restated as statements, and every comparison they make
     holds unchanged, to the bit.
   - `testForestLawStatements`: for K = 2 and 3, a forest with no basis and one with a two-level block
     that state nothing have, under each family, the leaf scale, the coefficients' variance and the
     half-Cauchy median of the tip's expressions written out with literals, to the bit. These literals
     survive the multiplier law. A numeric column's law is held by the restated fixtures alone, which
     that law rewrites.
   - `testForestSdUnit`, beside [`testForestMapWriters`](../../tests/cpp/test_sampler.cpp): one fixture
     of three forests (none, a two-level block, a column of numbers) per family, with a response whose
     standard deviation is neither 1 nor its range. (a) The derived unit equals the literal (the
     fixture's `sd`, 1, pi / sqrt(3)) to 4 ulp. (b) Each forest stating v has the leaf scale, or the
     half-Cauchy median, that the base spelling gives at v over the literal unit, to 4 ulp; a forest
     stating none has the base spelling's bits. (c) The reader returns v itself and "stated", and a
     default as its multiple times the unit and "not stated". (d) A write of the value read leaves 20
     sweeps identical to an untouched twin's and the forest stated; a write of another value gives the
     leaf scale of a chain created with it. (e) A chain made again over the response times 10 plus 5
     with the first one's anchor and unit handed back has every leaf scale of the first, bit for bit.
     (f) The two throws.
2. The bridge. [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) reads four numbers a forest
   (tree count, base, power, the hold). [`applyForestAttributes`](../../src/R_interface_bartcore.cpp)
   reads three more entries of the forests' configuration: `sd` and `variance`, one number or NA for
   each forest, and `unit`, as it reads `anchor`, with the two texts for a malformed one.
   [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) reports the unit, the sd in force and
   the stated flag as three more columns.
   [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp) hands its number on unchanged. No
   reading of a record of eight numbers is kept and no refusal by name is written for one: the length
   check answers.
3. Creation. [`forestParams`](../../R/model.R) writes the four numbers and derives nothing;
   [`defaultAmplitudePriorScale`](../../R/model.R) goes, and the two cites of it and of
   `amplitudePriorScale` in docs/plans/bcf-latent-evidence.md and docs/design/feature-matrix.md are
   marked `retired:`. The multi-forest block of [`resolveSamplerSpec`](../../R/spec.R) writes the
   stated sds, as written, into the configuration's `sd`, and a stated `amplitude.prior.variance` into
   `variance`; the sampler's `initialize` records `unit` from the reader where it records `anchor`, at
   a first creation only. Tests, a new file test-forest-sd-unit.R:
   - Which model was fitted. `sd = 2` on a forest with no basis and `sd = 0.5` on a factor forest, at
     each of: a formula term at `dbarts` and at `bart`, a list on a formula, on a matrix and in
     `dbartsSpec`, a data object with `bases` and a list of sds. At every door the reader's `sd` is
     `identical()` to 2 and to 0.5 and `sd.stated` is `TRUE`, `k.scale * 0.674` of the factor forest
     equals 0.5 and the half-Cauchy median times the unit equals 2, both to 1e-12, and the train draws of
     the sampler doors are `identical()` (fails today: `k.scale` is 0.5 sd(y) / 0.674).
   - L written out. With an offset, `subset`, a dropped missing response and weights with zeros, a
     default forest's reader `sd` is `sqrt(2 / K) * sd((y - offset)[kept])` to 1e-12, where `kept` is
     written in the test and leaves out the rows of weight zero (dec-B296); the same reader is
     `identical()` after a later `$setWeights` that brings those rows back, the unit being held; and the
     same four with `sd = 0.5` stated give `k.scale * 0.674 == 0.5` to 1e-12, so a unit from the wrong
     rows cannot hide in a ratio.
   - Probit and logistic: `sd = 0.5` gives `k.scale * 0.674 == 0.5` under both (fails today under
     logistic: 0.907).
   - The size, from the engine's own prior draws under a flat likelihood. Always on, one a family, 1500
     sweeps and 10 percent: a forest with no basis and a factor forest, each drawn and stating `sd = v`,
     have size v (fails today under gaussian and logistic). `at_home`, 5000 sweeps a row and 6 percent:
     the same at three values of v. `tinytest::test_package()` runs with `at_home` false, in CI too, so
     the first is the one a push is gated on.
   - Not stated: a forest that states no sd has the engine's own multiple untouched: the reader's
     `leaf.scale.factor` is `identical()` to `sqrt(2 / K)` where the forest has a basis, and its
     `amplitude.prior.scale` to 2 (1 under a latent family) where it has none, at K = 2 and 3 under each
     family; `sd.stated` is `FALSE` and `sd` is that multiple times the unit. A default sent through the
     unit and back rounds, and fails the first two. The draws of unstated fits are the pair script's and
     the three bitwise compares' to hold.
   - A variance still stated: `amplitude.prior.variance = 2` on a factor forest gives the base build's
     draws (the pair script's row), and the reader's entry says 2.
4. The reader, the writer and `extract`. [`reportLeafPrior`](../../R/dbarts.R) builds `leaf.prior`,
   `sd` and `sd.stated` from what the bridge reports, by rule 5, and names the list it returns for
   every forest by the forests' labels; [`resolveForestSpreads`](../../R/dbarts.R) keeps the writer's
   form, one number a forest, an entry with no sd leaving its forest alone, and
   [`writeForestSpreads`](../../R/dbarts.R) writes the stated number to the engine and into the
   configuration's `sd`; [`extractParameter`](../../R/generics.R) returns, for a fit of several forests,
   the list of the stored reader's `sd` entries named by label, and
   [`selectForests`](../../R/generics.R) hands back the entry itself where one forest is selected, as it
   does for a vector. The three entries of the mutation battery that quote those lines move with them
   (["m45"](../../benchmarks/R/mutation-battery.R), m46, m47). Tests, same file. They pin `sd`,
   `sd.stated` and what a round trip does, and never the class or the text of `leaf.prior`: the reader
   and the writer are to trade `normal()` objects in a later slice (dec-B246, dec-B266), and a pin of
   `forest()` literals would be written to be changed.
   - One path at a time, each asserting the stated number from the reader's `sd`, `sd.stated` and
     `k.scale * 0.674` (or the half-Cauchy median times the unit) against it: after creation; after
     `$setLeafPrior(forests = )`; on `copy()`; after `saveRDS`, `readRDS` and a first use; on
     `new("dbartsSampler", control, model, data)` from the sampler's own three; after
     `setResponse(10 * y + 5, updateScale = FALSE)` and on a copy made then; after a write, a copy and a
     second write; after `$setForestBasis` to another factor of the same levels, which derives the
     forest's own scale again. A unit applied twice gives 0.5 / sd(y), and one not applied 0.5 sd(y), on
     that path alone.
   - Read then write. `s$setLeafPrior(forests = lapply(s$getLeafPrior(), function(p) p$leaf.prior))`
     after 5 sweeps leaves `sd`, `sd.stated` and the next 20 sweeps `identical()` to an untouched
     twin's, on a stated and on a default sampler. A default forest written `forest(sd = p$sd)` has
     `sd.stated` `TRUE`, the same `sd`, and the next 20 sweeps `identical()` to the twin's.
   - Other data: a `forests` list with `sd = 2`, and a control taken from a sampler, given to a response
     ten times as wide: the reader says 2 and `k.scale * 0.674` is 2 on the new sampler (a carried unit
     gives 20 or 0.2).
   - Names. `names(s$getLeafPrior())` and `names(extract(fit, type = "leaf.prior.sd"))` are the forests'
     labels, for a formula of three `forest()` terms and for a named `forests` list (fail today: no
     names, and `forest1`, `forest2`).
   - `extract(fit, type = "leaf.prior.sd")` on a `bart` fit of three forests is `identical()` to the list
     of the kept sampler's `sd` entries; with `forest = 2`, with the label and with `"forest2"` it is
     the one number; with two labels a list of two in the order asked (fails today: 1.74 and 2.58 where
     the reader says 2 and 1); on one forest it is what it is today.
5. One forest. [`resolveForests`](../../R/model.R) no longer refuses `sd` on the one forest; where the
   fitting function's own leaf prior is read, that forest's `sd` is `normal(sd = )`, and stated at both
   places it is refused with the first text. "Stated" is forest-defaults-by-kind's record of what the
   caller named. Tests: `y ~ forest(x1 + x2, sd = 2)` at `bart` and `dbarts`, and
   `forests = list(forest(sd = 2))` on a formula and on a matrix, are draw for draw
   `leaf.prior = normal(sd = 2)`, under gaussian and probit (fail today: refused); beside
   `leaf.prior = normal(sd = 3)`, `normal`, `normal(k = 2)`, `bart`'s `k` and a linear leaf, refused with
   the first text; `amplitude` on that forest keeps its refusal.
6. The two other texts. The engine's throw on a response with no spread is raised with the second text,
   naming nothing of the engine; [`resolveLeafPrior`](../../R/model.R) takes the third. Tests: two
   forests on a constant response are created with no sd (as today) and refused with `sd = 0.7` (fails
   today: created, every prior 0); the four pins of
   ["a named leaf-prior 'sd'"](../../inst/tinytest/test-calibration-creation.R) take the third text.
7. A named single sd (dec-B280). [`validateForestSd`](../../R/model.R) drops the name of a single
   number and refuses nothing for it; the text for a longer `sd` loses its words about columns. Tests:
   `sd = pars["s"]` on a forest with no basis, on a factor forest and on a numeric one, at a formula
   term, in a list and in `$setLeafPrior(forests = )`, and on the one forest of a single-forest model:
   each is created, the reader's `sd` is `identical()` to `unname(pars["s"])`, and the draws are those
   of the unnamed number (fail today: refused); `c(dose = 1, age = 2)` and `c(1, 2)` keep the length
   text at the three places; the pins of
   ["must not be named"](../../inst/tinytest/test-forest-arguments.R) are turned into these.
8. Respell and repair. Run the suite first and repair what it shows. Expected: the 72 pins of the
   stored numbers by position, each restated against the four numbers or against the reader; the 15
   assertions of the unit in 7 files (test-bcf-creation.R 3, test-bcf-family.R 3, test-fit-stores-k.R 3,
   test-calibration-midchain.R 2, test-multiforest-leaf-prior-writer.R 2, test-forest-arguments.R 1,
   test-forest-basis-r5.R 1), five of them among the 72; the pins of the reader's `leaf.prior` for a
   forest that states nothing, which now read `sd`; and the three pins of
   ["this model has one forest, so its size is"](../../inst/tinytest/test-bcf-creation.R), which become
   step 5's fits. Each gaussian and logistic forest that states an sd, and each of the 28 gaussian
   writes, is read: where the test means a multiple of the response's scale (the pins of the map in
   test-bcf-family.R, the writer's ratios, the prior-draw tests) it is written `r * sd(y)` or
   `r * pi / sqrt(3)` with the scale computed in the test; where it needs only some stated sd it is left,
   and now states that many units. Benchmarks: the four BCF exact gates and `sbc.R` write
   `sd = sdControl * s` and `sd = sdModerate * s`, s the scale their oracles already use, and their
   oracles are not edited.
9. Help and records. man/forest.Rd: the `sd` item, as push 3 leaves it, with the unit, what s is the size
   of for each kind as "Context" measured it, the defaults as multiples of the response's standard
   deviation and as the numbers they come to, that a name on the number is ignored, and one sentence
   that a single forest takes it as `leaf.prior = normal(sd = )`; the Details paragraph on the budget
   where it speaks of units.
   [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) and
   [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd) with their docstrings: `leaf.prior`
   is what was stated and goes back as it is; `sd` is the number in force, in the response's units, and
   `sd.stated` whether it was stated; the list is named by label; `amplitude.prior.scale` and
   `leaf.scale.factor` are multiples of the response's scale. The `extract` page's `leaf.prior.sd` item:
   the list by label, a number for one forest. man/bart.Rd and man/dbarts.Rd: `leaf.prior` beside
   several forests. docs/design/forest-sd-unit.md with its index row: the rule, who holds the law, the
   unit's record, the table of "Context", the changed sequence and its oracle.
   docs/design/multiplier-combiner.md, bcf.md, nameable-calibration.md and public-surface.md where they
   give a forest's sd in units of the latent scale or have R derive the law; docs/architecture.md where
   it lists the forests' configuration. TODO: `forest-prior-args` names this slice landed and
   `named-single-sd` closes.
10. bartCause, same day (its own commit on dbarts-1.0). `bcf()`'s `sd.moderate` defaults to `NULL` in both
    signatures of R/bcf.R, and man/bcf.Rd says of `sd.control` and `sd.moderate` that each is the forest's
    size in the response's units, `NULL` taking dbarts's default. The two hand-built samplers in
    tests/testthat/test-14-bcf.R and test-03-responseFit.R drop `sd = 1`. Both edits are the same fit to
    the bit on either side of this slice, a default being the same number as `sd = 1` at K = 2 before it.
    Nothing else of bartCause meets this slice: it reads `response.scale` and `response.shift`, which
    stay. stan4bart, treatSens and bairrtt are not edited.
11. Mutations (Verification): apply each, install with `--preclean` where it is the engine's, run the
    named test, record the failing count, revert, `touch` the file.

## Verification

Against the slice's own library (`R CMD INSTALL --preclean -l <lib> .`, `R_LIBS=<lib>` on every call;
check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series, at most two
cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`, and again built with
  `OPT="-O2 -g -fsanitize=address,undefined"` under `ASAN_OPTIONS=detect_container_overflow=0`; the
  R-loaded path under the address sanitizer for test-forest-sd-unit.R and
  test-multiforest-leaf-prior-writer.R, as [Gate hygiene](README.md#gate-hygiene) gives the commands.
- The full tinytest suite on the shipped build, in one process, counted per file: no failure, no file
  stopping, and at least the base build's count plus the new file's assertions; the landing note gives the
  figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded: no scenario states a forest sd, and a
  scenario that did, or a law that moved with its move to the engine, would show here as a `max |z|`
  line. These three are the gate of the change of shape.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick`: the class is posterior-changing and
  the four respelled BCF gates are the oracle for a stated sd in each channel. `sbc.R`'s respelled arm
  constructs and runs its smallest setting.
- `inst/include/dbarts/dbarts.h` has no diff and `tools/check-api-hash.sh` passes.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy.
- The consumers, each suite whole against a private install of the slice, none failing: bartCause on
  dbarts-1.0 with step 10's edit (1412 expectations at its last run), stan4bart on bartcore (582),
  treatSens on dbarts-1.0 (306), bairrtt on main (207), the last three unedited. bartCause unedited is run
  once too and passes, its `sd = 1` now one unit and its tests pinning no draw under it; if that run
  fails, the landing note says where.

The pair script. Each row is one model fitted on the base build as written there and on the slice's as
written here, the same seed, 20 sweeps: train draws, sigma and the coefficients. Rows 01 to 06 and 21 to
28 must be `identical()`; rows 07 to 20 are equal to 1e-10 and the count of identical ones is recorded;
rows 29 and 30 must differ.

| | base build | slice |
|---|---|---|
| 01 to 04 | no sd stated: two forests; three; under `subset` with weights and an offset; `bart` with a term | the same text |
| 05, 06 | logistic and probit, no sd stated | the same text |
| 07 to 10 | `sd = r` on the forest with no basis; on a factor forest; on both; on a numeric column | `sd = r * sd(y)` |
| 11 to 13 | r stated with an offset; under `subset`; with weights | `r * sd(y - offset)`; `r * sd(y[kept])`; `r * sd(y[w > 0])`, which is `r * sd(y)` only where no weight is zero (dec-B296: the base build already holds response-scale-rows, so its unit is over those rows too) |
| 14, 15 | a held factor; three forests each stating one | the same respelling |
| 16 to 18 | `$setLeafPrior(forests = )` after 5 sweeps; a copy after 5; a value read and written back | the same respelling; the same text |
| 19, 20 | logistic, r on both forests; `bart` with a term stating r | `r * pi / sqrt(3)`; `r * sd(y)` |
| 21 | probit, `sd = 0.7` on both forests | the same text |
| 22 to 25 | no sd stated: a held forest with no basis; a held two-level factor; a numeric column and a 0/1 number; no forest without a basis | the same text |
| 26 to 28 | `amplitude.prior.variance = 2` on a factor forest; a swap of a numeric column to three times it, then 20 sweeps; a copy, a reload and `new()` from the sampler's own three | the same text |
| 29, 30 | gaussian, `sd = 0.7` on both; logistic the same | the same text: another prior |

Mutations, each expected to fail the named test and no gate before it:

- the unit taken to be the anchor under gaussian, the transform's scale forgotten: step 1 (a), and step
  3's "which model was fitted";
- the logistic unit taken to be 1: step 1 (a), and step 3's logistic row;
- a stated sd divided by the unit in R as well as in the engine: step 3's "which model was fitted" at
  every door; in the writer alone: step 4's path "after `$setLeafPrior`";
- a stated sd not divided at all at the data-object door (the list of sds over `bases`): step 3's row
  for that door;
- the default sent through the unit and back: step 3's "not stated", and the three bitwise compares;
- the default of a forest with a basis computed for K - 1 forests: step 1's literals at K = 3;
- a forest with no basis given the multiplied forest's default: step 1's literals, and pair row 01;
- the hold read from the wrong one of the four numbers: pair rows 14, 22 and 23, and the suite's held
  pins;
- the stated variance dropped between R and the engine: step 3's variance test, and pair row 26;
- the row norm left out of the engine's law: the restated fixtures, and pair rows 10, 24 and 27;
- the reader returning the multiple of L: step 3's `identical()` to 2;
- the reader's `leaf.prior` built from the number in force, so that a default reads as stated: step 4's
  "read then write" on a default sampler (`sd.stated` turns `TRUE`);
- the unit derived again at a re-creation: step 1 (e), and step 4's path after `setResponse`;
- the scale derived after a swap from the stated number with no division: step 4's path after
  `$setForestBasis`;
- the unit carried on a control to other data: step 4's "other data";
- a write of the number in force skipped whole, the forest left unstated: step 4's `forest(sd = p$sd)`;
  the same write dividing again, so that a bit moves: the same test's twin;
- the unit computed with weights, over every row before `subset`, or over the rows of weight zero
  too: step 3's "L written out";
- `extract` left returning `k.scale` over k; left a vector; named by position; a list of one where one
  forest is asked for: step 4's `extract` and names tests;
- the stated number mirrored into the four numbers and not into `sd`: step 4's path "a write, a copy
  and a second write";
- the both-given refusal dropped: step 5;
- a name kept on the stored sd: step 7's `identical()`.

Not a hot-path change: the law is resolved when a forest is built or restated, never in a sweep.

## NEWS

No new item: `forest()` and its `sd` are new in 1.0-0 and nothing released changes.

## What this leaves for the multiplier law

- What an sd is the size of where the coefficient is held or the basis is numeric: the table of "Context"
  with its rows off 1 is still true after this slice, in the response's units. The help says so, in those
  words, until the law changes them.
- The engine after this slice is handed a statement and holds the tip's law in one function. The law's
  pushes change rows of that function and add one thing to the statement, the forest's kind; they keep
  `sd` not-a-number for "not stated", the recorded unit, the division, the reader's `sd` and
  `sd.stated`, and the four-number record.
- The reader's `amplitude.prior.variance`, `amplitude.prior.scale`, `leaf.scale.factor`,
  `leaf.scale.divisor`, `basis.row.norm` and the words of `prior.sd.of`; a numeric forest's `sd` named
  by its column; the entries for a forest's kind and hold; the printed block: none is built or renamed
  here.
- What the law may assume: every stated sd in the tree, the benchmarks and bartCause is in the response's
  units; the reader's and `extract`'s shapes are final but for the column's name; a fit of forests with
  no basis or a two-level factor, coefficients drawn, stated or not, is unchanged by the law to the bit,
  and can be gated so.

## What waits on what

- On push 3: the help's `sd` item is written over push 3's; every call in the tests and the help is in
  its spelling; the forests' labels, which name the reader's and `extract`'s lists.
- On forest-defaults-by-kind: step 5 reads its record of what the caller stated and its finder of the
  forest with no basis; step 3 edits [`forestParams`](../../R/model.R) after it does, and keeps the
  count, base and power it gives every forest as the first three of the four numbers; the third text
  stands beside its refusal of a leaf prior where no forest is without a basis, which comes first;
  selection by label, which step 4's `extract` by label uses; its interim refusal of held shapes, which
  this slice does not touch.
- To recheck once both have landed: the counts of "Context" (51 forests, 33 writes, 72 and 15
  assertions), by running the logging build again; that no baseline scenario has gained a stated sd;
  that [`resolveForests`](../../R/model.R) still holds the one-forest refusal step 5 removes; the names
  of the functions cited here.
- Order with the two slices after it: this one, then the kind by class, then the multiplier law. This
  slice and the kind by class share R files and are serial; neither needs the other. They are not one
  slice: this one changes arithmetic in the engine and a prior, that one is R only. It is not folded
  into the law's pushes that change a prior either: the law's cleanest gate is that forests with no
  basis or a factor are bit for bit under it, stated sds included, and that gate exists only if the unit
  moved first.

## Out of scope, and where it goes

- The fitting function's `leaf.prior = normal(sd = )` as the size of the forest with no basis in a model
  of several: refused here with the form to write. It belongs with the reader's `normal()` objects and
  `$setLeafPrior`'s single form, in TODO `forest-prior-args`.
- The kind by class, with the refusal of a factor of three or more levels and of several numeric
  columns; the multiplier law with everything listed above; `leaf.prior = normal(sd = )` on `forest()`
  and `updateBasisScale`, after the merge to main (dec-B276).
- A way to put a stated forest back to its default: none exists before `normal()` objects are traded.
- To TODO as a new entry: the verbose summary of a model of several forests prints one tree count and
  "k prior fixed to 2", the first forest's model slots, and nothing of the others.

## Calls made in planning

- The slice stands alone, before the kind by class and the law. Its draws do move, for a stated sd left
  as written, so it is not a renaming that could ride with either; the reasons against folding it are
  under "What waits on what".
- The multiplier law's first push as first planned is folded in here (the critique's finding 7, taken by
  the coordinator). That push handed the engine the statement and moved no draw; this slice as first
  planned added a stated sd beside four derived numbers it kept, a reader that returned a default as a
  number, a skip of an equal write, a reading of records saved before it and `extract` as a vector, and
  the law then removed or reversed each and repaired the same pins a second time. Folded: one gate
  battery fewer, the 72 pins repaired once, and nothing built to be removed. The cost is a review that
  holds two things, a change of shape that must move no bit and a division that must happen once; the
  bitwise compares gate the first and the unit's own tests the second, and the reviewer reads the diff
  once for each. Not folded: the kind's route to the engine, which the tip's law does not read.
- The engine converts, against a unit it derives and R records. The alternative, R dividing by
  `sd(y - offset)` before the bridge, is shorter and was built as the stand-in. It was not chosen: R's
  estimate differs from the engine's in the last place, so the reader could not return a stated number
  exactly; and a `forests` list or a control resolved against one data object and run on another would
  carry the first one's unit in silence, where the engine derives the unit of the data it is given.
- The record is one more number beside the anchor. It could be recomputed from the anchor and the
  recorded response range, which live in two places and would have to agree to the bit.
- The reader gives the statement, a default as `forest()`, with the number beside it (dec-B269's words:
  "returns what was stated, the default as the default"). The first plan returned
  `forest(sd = <the default>)`, which reads a default as a statement and, written to another sampler,
  states the first one's number there.
- A write of the number in force marks the forest stated and moves nothing. Skipped whole, as first
  planned, a caller could not state the number the reader gives.
- The tests pin `sd`, `sd.stated` and the round trip, not the class of `leaf.prior` (the coordinator, on
  the critique's question (e)): the reader and writer trade `normal()` objects in a slice not yet
  planned.
- `extract(type = "leaf.prior.sd")` takes its final shape here, with its new numbers: a list by label,
  and for one forest the entry itself (the coordinator, on the critique's question (c)). The law's plan
  had a list whatever was selected; `extract(type = "k", forest = 2)` returns a number, and one
  argument should not give a number for one quantity and a list for the other. At the tip `extract` is
  already in the response's units, of the forest's own standard deviation before its coefficient, which
  is neither the stated number nor what dec-B253 and dec-B275 say it returns.
- The reader's list and `extract`'s are named by label here and not with the law: dec-B275 names them
  so, forest-defaults-by-kind has landed selection by label by then, and a second renaming of the same
  list two slices on would repair its pins twice. A margin of draws keeps `forest<i>`.
- A named single sd is taken on every forest, with no condition (dec-B280 as dec-B282 leaves it). At
  this slice's tip a basis may still have several columns, where dec-B280 alone would refuse the name;
  that case is closed by the kind-by-class slice's second push, and a condition built for the tips
  between would be removed there.
- No refusal by name, and no reading, of a record written by an earlier build of this branch (the
  coordinator's cut): only a development build wrote one, and the bridge's check of the record's length
  stops it.
- One prior-size test a family runs always, at 1500 sweeps and 10 percent (the critique's finding 3):
  the first plan's were all `at_home`, which the package's test run never is.
- `amplitude.prior.scale` and `leaf.scale.factor` stay as the engine's multiples of L and are documented
  so. They are entries the law removes; converting them would change two pins to change them again.
- A stated sd on a response with no spread is refused. Today it is taken and every prior is 0. Dividing
  by a unit of 0 would hand the engine an infinite scale.
- One forest takes `sd` in this slice, as written-surface's table and the design say it does once the
  unit is one. The fitting function's `normal(sd = )` beside several forests is not served here: the
  design gives that to the slice that trades `normal()` objects, and the text now says where to write it.
- A test's stated sd is respelled only where the test means a multiple of the response's scale. Every
  other one is left to state units; the alternative, respelling all of them, keeps draws nothing pins.
- The pair script is run at landing and not tracked: its old side needs the base build. What stays in
  the suite is step 3's pins of unstated fits and step 4's identities between paths.
- Measured for this plan on stand-ins, not on the design: the 15 assertions, the 8 of 18 pairs and the
  two "left as written" distances come from a build that converts in R; the 72 pins from the law's
  stand-in, which keeps six numbers a forest. The engine's conversion divides by a number one unit in
  the last place away, so the counts of identical pairs may differ; the class does not. Nothing of the
  engine's statement is built.
