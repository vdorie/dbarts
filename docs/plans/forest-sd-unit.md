# forest-sd-unit: a forest's sd is stated and reported in the response's units, and the engine is handed the statement

Status: LANDED 2026-10-11 (e84c0dd0 to 37b59ea0; calls dec-A208; [Landing note](#landing-note)). Rests on
dec-B246 as revised, dec-B253, dec-B266, dec-B269 as dec-B280 and dec-B282 revise it, and dec-B296 and
dec-B302 for the rows the unit is read from; rechecked 2026-10-11 against d5b637ee, bartcore with the
creation slice of [response-scale-rows.md](response-scale-rows.md)
([Recheck against bartcore](#recheck-against-bartcore)). The one open call, on which slice builds dec-B291,
went to the multiplier slice ([Open calls](#open-calls)). Follows
[forest-defaults-by-kind.md](forest-defaults-by-kind.md) and [written-surface.md](written-surface.md),
landed, and the creation slice of response-scale-rows. It carries what was first planned as the multiplier
slice's neutral first push (the engine is handed what was stated and holds the rule), gives the reader and
`extract` the shape they keep, and takes a named single sd (dec-B280).

agent: one push. Opus implementer for the engine, the bridge and the R code; sonnet for the respelled tests,
the benchmark scripts and the help once the code is fixed; opus reviewer told to refute, and to read the
diff twice: once for the change of shape, which must move no bit, and once for the division, which must
happen exactly once. The reason for opus: every slip is a prior off by a factor of sd(y) with no message.
The conversion is one division, and it can be made twice, or not at all, on any one of six paths
(creation at four doors, the writer, a re-creation); the review of written-surface's second push found
three such one-path faults in a slice that touched no arithmetic at all. Tier: changes draws
([Process by risk](README.md#process-by-risk)).
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- NEUTRAL for every fit that states no forest `sd`, and for every fit under probit, where the unit is 1.
  That covers the change of shape too: the engine resolves the rule R resolved, from four numbers a
  forest where it was handed eight. Measured for the unit on a stand-in at d5b637ee: 6 of 6 such pairs
  identical, a pair created with ten rows of weight zero among them, and no seeded draw of the suite
  moved. Not measured for the change of shape, none of which is built: its gates are the three bitwise
  compares, the snapshot files, the pair script and the restated tests/cpp fixtures.
- POSTERIOR-CHANGING for one sequence: a model of several forests that states `sd` on a forest, at
  creation or through `$setLeafPrior(forests = )`, under a gaussian or a logistic response, and is left as
  written. The same number now states another prior: `sd = 0.7` was 0.7 standard deviations of the
  response, and is 0.7 units of it. The oracle is the same model on the base build with the number divided
  by the unit.
- SHIFTING for that sequence respelled into the response's units (`sd = 0.7 * sd(y)`): the prior is the
  same to rounding and the draws are not always the same bits. Measured on the stand-in, 18 gaussian and
  logistic pairs of 20 sweeps: 9 identical, 9 apart by 1e-14 to 3e-14.
- No scenario of the three equivalence baselines and no snapshot file states a forest `sd` (searched
  again at d5b637ee), so no recorded scenario moves and nothing is re-recorded. Four exact gates and
  `sbc.R` state one and are respelled.
- Accepted where refused: `sd` on the one forest of a single-forest model; a named single `sd` on a
  forest and in the writer, its name dropped.
- Not a released object: a sampler saved by an earlier build of this branch holds eight numbers a forest
  and is refused by the bridge's check of the record's length.
window: pre-release, after the creation slice of response-scale-rows, which fixed the rows the unit is
read from, and before the kind by class and the multiplier slice, so that every later test, pinned string
and help page states an sd once, in the unit it keeps, and reads it back in the shape it keeps. Serial
with any other work in [`forestParams`](../../R/model.R), the multi-forest block of
[`resolveSamplerSpec`](../../R/spec.R), [`reportLeafPrior`](../../R/dbarts.R),
[`writeForestSpreads`](../../R/dbarts.R), [`extractParameter`](../../R/generics.R),
[`applyForestAttributes`](../../src/R_interface_bartcore.cpp),
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp), [`ForestSpec`](../../src/bartcore/combiner.hpp),
[`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) or the K-forest constructor of
[`Chain`](../../src/bartcore/chain.hpp). bartCause's edit lands the same day.
budget: ~2340 lines planned: engine ~300; bridge ~160; R ~300; tests/cpp ~480; tinytest ~710, of which
~480 in one new file and ~230 changed in about a dozen; benchmarks ~50; help ~130; design note,
architecture, public-surface, TODO and the two indexes ~210. Forecast 4200 to 5400; stop at ~5400, or at
~1140 of package code ([Budget and stops](#budget-and-stops)). bartCause's part, ~10 lines, is outside it.

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

Measured at the tip (d5b637ee) on the shipped build, R 4.6.1, unless a line says otherwise. A basis
written as code carries a tilde. L is the response's scale as the engine holds it; "size" is the prior
median of the absolute multiplier times the prior standard deviation of the forest's own fit, exact where
the coefficient is held.

- The suite. 253 files, of which 4 exit off a reference build and 3 ask for more than two threads; the
  other 246 give 22978 results, none failing, 120 seconds in one process.
- L. Under gaussian it is the sample standard deviation of the response net of its offset, n - 1 in the
  divisor, unweighted, over the rows in the likelihood when the sampler is created: the rows `subset`
  and the na.action keep, less those of weight zero (dec-B296). 3.26736095623188 from the engine against
  R's 3.26736095623189 for `sd(y)`, one unit in the last place apart; `sd(y - offset)` with an offset;
  `sd(y[subset])` under `subset`; unchanged by positive weights; with 50 rows of weight zero
  3.32872882757977, which is `sd(y[weights > 0])`. It is held when the weights or a mask change later:
  the same number after `$setWeights` brings the rows back. Where a row is out and the rows in hold fewer
  than two distinct values it is over every row, under one warning of class
  `dbartsScaleFallbackWarning` (dec-B302, dec-B309). Under probit it is 1 and under logistic
  pi / sqrt(3), 1.8138. It is computed at a first creation
  ([`latentScaleAnchor`](../../src/bartcore/chain.hpp), which reads
  [`scaledResponseSdOver`](../../src/bartcore/chain.hpp) over the rows
  [`GaussianResponse::scaleRows`](../../src/bartcore/model.hpp) marks, and
  [`scaledResponseSd`](../../src/bartcore/chain.hpp) where the response model read every row), recorded
  by R as the forests' `anchor` and handed to every re-creation
  ([`AmplitudeSpec::anchor`](../../src/bartcore/combiner.hpp)), so it is the same number after
  `setResponse(10 * y + 5, updateScale = FALSE)` and on a copy; `updateScale = TRUE` is refused on such a
  sampler (dec-B362).
- What a stated or default number s is the size of today. Engine prior draws under a flat likelihood,
  5000 sweeps a row, Monte Carlo error about 0.02; each entry is the measured size over s L. Rerun at
  d5b637ee at 2000 sweeps a row, the two rescaled columns apart: every entry within three standard
  errors of the table.

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
  | one numeric column, held | | refused (dec-B353) | refused | refused |

  So s is in units of L. It is the size of the forest's contribution where the forest has no basis or a
  factor and its coefficient is drawn. It is not where the coefficient is held (a held forest with no
  basis has size L whatever s is; a held factor 1.48 s L) or the basis is numeric (0.71 s L at a row whose
  basis has the median nonzero norm). Those are the multiplier slice's to change; this slice changes the
  unit of s and nothing in that table but its last factor. (Rows for a factor of three levels and for
  two numeric columns were measured and agree; both shapes are refused before the release, dec-B281 and
  dec-B282, and are left out.)
- Who resolves the rule today. R does: [`forestParams`](../../R/model.R) turns each forest's statement
  into eight numbers (tree count, base and power, by the forest's kind where it states none; a leaf-scale
  factor, its divisor, the coefficients' variance, a half-Cauchy median, the hold), the bridge passes
  them ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)), and
  [`ForestSpec`](../../src/bartcore/combiner.hpp) carries the last five as `nodeScaleFactor`,
  `nodeScaleDivisor`, `amplitudePriorVariance`, `amplitudePriorScale` and the hold flag. The K-forest
  constructor of [`Chain`](../../src/bartcore/chain.hpp) then divides by
  [`basisRowNorm`](../../src/bartcore/chain.hpp) ([`mapLeafScale`](../../src/bartcore/chain.hpp)). The
  defaults live in R (retired: [`defaultAmplitudePriorScale`](../../R/model.R), an R function this slice removes, and `sqrt(2 / K)` written in
  `forestParams`). Of a forest's kind the rule reads one thing, whether it has a basis. The forests'
  record on the control holds `params` (the eight numbers), `vars`, `labels`, `heldAmplitude`,
  `interactions`, `blocks`, `basisTerms`, `control.n.trees` where the first forest has a basis, and
  `anchor`; [`dbartsSampler$setForestBasis`](../../R/dbarts.R) reads the hold at the eighth number.
- What reads and writes the number.

  | place | today | unit |
  |---|---|---|
  | `forest(sd = s)` at creation; `$setLeafPrior(forests = list(forest(sd = s)))` | s goes to the engine as given ([`forestParams`](../../R/model.R), [`writeForestSpreads`](../../R/dbarts.R)) | of L |
  | `$getLeafPrior(f)$leaf.prior` | `forest(sd = s)`, s the engine's own number, stated or not | of L |
  | `$getLeafPrior(f)$k.scale` | the forest's own prior standard deviation, before its coefficient: L for a forest with no basis, s L / 0.674 over the row norm for one with a basis | the response's |
  | `amplitude.prior.scale`, `leaf.scale.factor` entries | s again | of L |
  | `$getLeafPrior()` with no forest | a list with no names | |
  | `extract(fit, type = "leaf.prior.sd")` | `k.scale` over k, a vector named `forest1`, `forest2`: 1.7408 and 2.5828 for a default two-forest fit whose sd(y) is 1.7408; with `forest = 2`, with `"forest2"` or with the forest's label the bare number ([`selectForests`](../../R/generics.R)) | the response's, of another quantity |
  | `print` of a fit, `show` of a sampler, the verbose summary | no forest's sd is printed | |
  | `summary()` of a fit of several forests, the line under its table | `extract`'s vector, as `leaf.prior.sd[forest1] = 1.741, leaf.prior.sd[forest2] = 2.583` | the response's, of another quantity |
  | `leaf.prior = normal(sd = s)` on a single forest | the prior standard deviation of the forest's total, exactly; held by the engine as k = `k.scale` / s, `k.scale` the data's (dec-B416) | the response's |

  The reader and `extract` already disagree: for one default factor forest the reader says 1 and `extract`
  2.58. And the reader cannot say whether a number was stated.
- One forest. `y ~ forest(x1 + x2, sd = 2)` and `forests = list(forest(sd = 2))` are refused at `dbarts`
  and `bart`, on a formula and on a matrix ("this model has one forest, so its size is the fitting
  function's leaf.prior = normal(sd = ), not forest(sd = )"). `leaf.prior = normal(sd = 2)` beside several
  forests is refused ("a multi-forest model does not support a named leaf-prior 'sd': the leaf prior's
  'sd' is not a forest's 'sd'"); `normal(k = 3)` beside them is refused, and a bare `normal` and
  `normal(k = 2)` accepted.
- A named single sd. `forest(basis = dose, sd = pars["s"])` is refused at creation and in
  `$setLeafPrior(forests = )` ("forest 'sd' must not be named"), on a forest with a basis and on one
  with none; `leaf.prior = normal(sd = pars["s"])` on a single forest is taken.
- A response with no spread. Two forests on a constant response are created with a warning; L is 0 and
  every forest's prior standard deviation is 0, with `sd = 0.7` stated or not.
- Where a forest's sd is stated, counted by running the suite with the two R functions traced: of 1427
  models of several forests created in 66 test files, 59 forests state one, in 8 files
  (test-forest-selection.R 15, test-bcf-family.R 14, test-bcf-creation.R 8,
  test-multiforest-leaf-prior-writer.R 8, test-forest-basis-r5.R 7, test-forest-arguments.R 5,
  test-forest-defaults.R 1, test-formula-terms.R 1): 36 under gaussian, 6 under logistic, 17 under
  probit. `$setLeafPrior(forests = )` writes one 39 times, in 24 calls of 5 files, 34 of them under
  gaussian. 11 forests state `amplitude.prior.variance`, in 3 files. In benchmarks/R the four BCF exact
  gates and `sbc.R` state `sd = sdControl` and `sd = sdModerate`, 10 lines; `bcf-equivalence.R`,
  `equivalence.R` and `multinomial-equivalence.R` state none. No vignette and no help example states one.
- What moves in the suite. For the unit, measured on a stand-in (the conversion done in R against
  `sd(y - offset)` over the rows of positive weight): no file stops and 17 assertions fail in 9 files:
  pins of the stored numbers behind a stated sd (6), of the reader's `leaf.prior`, `leaf.scale.factor`
  and `amplitude.prior.scale` (5), of `k.scale` against a stated multiple (2) and of `extract` (4). For
  the change of shape, counted from the expectations the suite records at the tip, each read for
  whether it indexes the fourth to the eighth stored number or compares against a literal of eight:
  about 60 assertions in 7 files (test-bcf-family.R 33, test-forest-arguments.R 11,
  test-bcf-creation.R 10, test-calibration-midchain.R 2, test-forest-basis-r5.R 2,
  test-forest-capture.R 1, test-forest-defaults.R 1), and two records edited by hand in
  test-forest-basis-r5.R. The two sets share six. Pins of the first three numbers (55 in
  test-forest-defaults.R alone) and compares of one sampler's record with another's stand. In tests/cpp
  about 60 lines of test_sampler.cpp and test_facade.cpp set or read the four derived numbers, most of
  them in the fixtures of [`testForestMapWriters`](../../tests/cpp/test_sampler.cpp),
  [`testBCFCalibrationMap`](../../tests/cpp/test_sampler.cpp) and
  [`testForestCalibration`](../../tests/cpp/test_sampler.cpp).
- The same model on both builds, the stand-in's side respelled as s L: with no sd stated, and under
  probit, identical; of 18 gaussian and logistic rows (each kind of forest, an offset, `subset`, positive
  weights, ten weights of zero, three forests, a formula term, the writer, a copy, a value read and
  written back, `bart`) 9 identical and 9 apart by at most 2.8e-14 after 20 sweeps. Left as written,
  `sd = 0.7` on both forests gives draws 2.31 apart under gaussian and 3.63 under logistic.
- Consumers. bartCause's `bcf()` (dbarts-1.0, e833be7) passes `sd = sd.control`, `NULL` by default, and
  `sd = sd.moderate`, 1 by default, and two of its test files build the same sampler by hand with
  `sd = 1` on the treatment forest. At K = 2 the default of a forest with a basis is exactly 1, so
  `sd = 1` and no sd are one fit today. It reads `getLeafPrior(1L)$response.scale` and
  `$response.shift` and no other entry. stan4bart (bartcore, d6cdad7), treatSens (dbarts-1.0, babfaa6)
  and bairrtt (main, 3f57f61) fit one forest, state no forest sd and read no leaf prior of a model of
  several forests (searched again 2026-10-11).

## The rule

1. The unit. A forest's `sd` is in the response's units: the units of y under a gaussian response, and of
   the latent index under probit and logistic. One number of the sampler, the unit, converts it: L in the
   response's units, which under gaussian is `sd(y - offset)` over the rows in the likelihood when the
   sampler is created (the rows `subset` and the na.action keep, less those of weight zero; dec-B296),
   and over every row where those rows hold fewer than two distinct values, under creation's one warning
   (dec-B302, dec-B309); 1 under probit and pi / sqrt(3) under logistic.
2. Stated. The engine is handed the number as written and divides it by the unit; what it then does with
   the quotient is what the tip does with a number in units of L. Nothing else of the calibration map
   changes.
3. Not stated. The engine is told so and takes the default itself: the tip's multiple of L (2 or 1 for a
   forest with no basis, sqrt(2 / K) for one with a basis), which never passes through the unit on the
   way in. dec-B291 makes the second 1 whatever K; which slice builds that is the open call.
4. Who holds the map. The engine, in one function: a forest's statement (an sd or none, a hold, whether
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
   and data. No later change of the response, the weights or a mask moves it, as none moves the anchor;
   `updateScale = TRUE` stays refused on several forests (dec-B362). A control or a `forests` list
   carried to other data states the same numbers of that data's units: the record is not carried, and
   the new sampler derives its own.
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
| its `amplitude.prior.scale`, `leaf.scale.factor`, `leaf.scale.divisor`, `amplitude.prior.variance`, `basis.row.norm`, `prior.sd.of` | the engine's channels, in units of L | unchanged, and said to be so; the multiplier slice removes or renames them |
| `$getLeafPrior()` with no forest | a list with no names | named by label |
| `extract(type = "leaf.prior.sd")`, several | a vector `forest1`, `forest2` of each forest's own standard deviation before its coefficient | a list by label of each forest's `sd`; for one forest the number |
| the same on one forest; `extract(type = "k")` | | unchanged |
| `print`, `show`, the verbose summary | no forest's sd | unchanged; the printed block is the multiplier slice's |
| `summary()`'s line under the table, several forests | `leaf.prior.sd[forest<i>]`, `extract`'s vector | `leaf.prior.sd[<label>]`, the reader's `sd` entries as `extract` returns them |
| `$setForestBasis` | the forest's own scale is derived again from the new block's row norm, under the number in force | unchanged, the number in force being the stated one or the default |
| `$setWeights`, `$setActiveRows` | move no forest's scale | unchanged: the unit is creation's |
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
forest(sd = )"; retired: ["must not be named"](../../inst/tinytest/test-forest-arguments.R), with its hint to
`unname()`. Replaced by the third: "a multi-forest model does not support a named leaf-prior 'sd'", whose
reason, that the two sds are different things, stops being true. The fourth is today's text for a longer
`sd` without "for every column of its basis". Kept as they are: every other text of
[`validateForestSd`](../../R/model.R). The second is met only where every row of the response less its
offset holds one value: where the rows in the likelihood do and some row out differs, the unit is the
every-row one and creation's warning has already said so.

## Constraints

- A fit that states no forest sd draws what it draws on the base build, to the bit, in every family. No
  baseline is re-recorded and no snapshot file regenerated.
- The calibration map is not touched: which channel a forest's number reaches, the 0.674, the
  coefficient's variance, the row norm, the held values and the starts by position, the defaults and the
  K in them. It moves from R to the engine with every expression in its order of operations.
  `amplitude.prior.variance` stays an argument.
- The map is resolved in one function of the engine and nowhere in R; a forest's entry is resolved once
  for each construction, by the chain, and the combiner is handed the result.
- The engine learns nothing of a forest's kind here beyond whether it has a basis, which is all the
  tip's rule reads. The kind's route from the data object is the multiplier slice's held push's.
- The conversion is made in one place, the engine, against one recorded number. R divides nothing and
  multiplies nothing: it hands over what was written and reports what the engine returns.
- The rows the unit is read from are those the response model read its range from, asked of it once
  ([`GaussianResponse::scaleRows`](../../src/bartcore/model.hpp)); this slice adds no second test of
  which rows are in.
- A stated number is kept as written, so that the reader returns it bit for bit.
- The reader's `sd` is one unnamed number a forest in this slice. A numeric forest's takes its column's
  name with the multiplier slice, when it becomes that many units per unit of the column.
- The forests' record holds four numbers a forest (tree count, base, power, the hold), `sd` and, where
  one is stated, the coefficient variance, each one number or NA a forest, `hasBasis` (whether each forest
  was created with a basis, which `$setForestBasis` can change for a forest created with none while the
  four numbers do not record the channel its sd reached), the anchor and the unit; its other entries (`vars`, `labels`, `heldAmplitude`, `interactions`, `blocks`, `basisTerms`,
  `control.n.trees`) are as they are.
- The flat C header, the stored state and every state block are untouched.
- A `forests` list written before the slice is still accepted, with every number read in the new unit.
- No test this slice adds fits a factor of three or more levels or a basis of several numeric columns
  (dec-B281, dec-B282: refused one slice on).
- Each commit leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the
reader's sake. The fixture is 150 rows with x1, x2, a 0/1 z, a factor zf, dose, an offset column and
weights of which ten are 0.

1. The engine: the statement, the map and the unit. [`ForestSpec`](../../src/bartcore/combiner.hpp)
   loses its four derived numbers and states: `sd`, in the response's units, not a number where none is
   stated; the hold it has; and the coefficient variance, not a number where none is stated (the
   multiplier slice's held push removes it). For fixtures and for
   [`expandForestSpecs`](../../src/bartcore/combiner.hpp), whose two-forest spelling states BCF's sizes
   as multiples of L, it may state a multiple of L in place of `sd`; nothing R reaches sets that.
   [`AmplitudeSpec`](../../src/bartcore/combiner.hpp) gains `unit`, not a number where the engine is to
   derive it. One function (`forestMapEntry`) turns a statement, whether the forest has a basis, K, the
   family, the unit and the row norm into the forest's leaf-scale factor and divisor, its coefficients'
   variance and its half-Cauchy median, and one (`defaultForestSd`) holds the defaults; both give the
   tip's numbers, the row norm and the starts by position included, with every expression of
   [`mapLeafScale`](../../src/bartcore/chain.hpp) in its order of operations. The K-forest constructor
   of [`Chain`](../../src/bartcore/chain.hpp) takes the unit it is handed, or derives it beside the
   anchor (under gaussian the anchor [`latentScaleAnchor`](../../src/bartcore/chain.hpp) returns, which
   is already over the rows the response model read, times the scale of the response transform at that
   moment; under probit and logistic the anchor), resolves each forest once and hands
   [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) the result. A forest that states an sd
   has `sd / unit` where its number in units of L stood: the half-Cauchy median of a forest with no
   basis, the leaf-scale factor of one with a basis. A stated sd with a unit that is not positive and
   finite, and an sd that is not, throw. [`ForestCalibration`](../../src/bartcore/chain.hpp) gains the
   unit and the sd in force in the response's units (the stated number itself, or the multiple of L
   times the unit), and reports the stated number in its `namedSd`, which a map forest leaves not a
   number today, so that "stated" is one test for every forest.
   [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) takes the response's units and marks the
   forest stated; where the number is the sd in force as reported it keeps the quotient it has, so no
   bit moves, and otherwise divides as creation does. Every struct here is read by objects that do not
   track headers: `--preclean`.
   Tests, tests/cpp:
   - The fixtures that set the four numbers are restated as statements, and every comparison they make
     holds unchanged, to the bit; [`testForestUnitRowsIn`](../../tests/cpp/test_sampler.cpp) among them,
     whose anchor is unchanged.
   - `testForestMapStatements`: for K = 2 and 3, a forest with no basis and one with a two-level block
     that state nothing have, under each family, the leaf scale, the coefficients' variance and the
     half-Cauchy median of the tip's expressions written out with literals, to the bit. The literals
     for a forest with no basis, and those at K = 2, survive the multiplier slice; those for a block at
     K = 3 hold sqrt(2 / 3) and go with dec-B291. A numeric column's entry is held by the restated
     fixtures alone, which the multiplier slice rewrites.
   - `testForestSdUnit`, beside [`testForestMapWriters`](../../tests/cpp/test_sampler.cpp): one fixture
     of three forests (none, a two-level block, a column of numbers) per family, with a response whose
     standard deviation is neither 1 nor its range. (a) The derived unit equals the literal (the
     fixture's `sd`, 1, pi / sqrt(3)) to 4 ulp; with rows of weight zero that hold the response's
     extremes it is the literal over the rows in, far from the every-row one, and with the rows in at
     one value the every-row one. (b) Each forest stating v has the leaf scale, or the half-Cauchy
     median, that the base spelling gives at v over the literal unit, to 4 ulp; a forest stating none
     has the base spelling's bits. (c) The reader returns v itself and "stated", and a default as its
     multiple times the unit and "not stated". (d) A write of the value read leaves 20 sweeps identical
     to an untouched twin's and the forest stated; a write of another value gives the leaf scale of a
     chain created with it. (e) A chain made again over the response times 10 plus 5 with the first
     one's anchor and unit handed back has every leaf scale of the first, bit for bit. (f) The two
     throws.
2. The bridge. [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) reads four numbers a forest
   (tree count, base, power, the hold). [`applyForestAttributes`](../../src/R_interface_bartcore.cpp)
   reads three more entries of the forests' configuration: `sd` and `variance`, one number or NA for
   each forest, and `unit`, as it reads `anchor`, with the two texts for a malformed one.
   [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) reports the unit and the sd in force as
   two more columns, and a map forest's stated number in the `named.sd` column it has.
   [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp) hands its number on unchanged. No
   reading of a record of eight numbers is kept and no refusal by name is written for one: the length
   check answers.
3. Creation. [`forestParams`](../../R/model.R) writes the four numbers and derives nothing: the count,
   base and power it gives every forest by its kind stay its first three, and
   [`dbartsSampler$setForestBasis`](../../R/dbarts.R) reads the hold at the fourth.
   retired: [`defaultAmplitudePriorScale`](../../R/model.R) goes, and the cites of it in
   docs/plans/bcf-latent-evidence.md and of `ForestSpec::amplitudePriorScale` in
   docs/design/feature-matrix.md are marked `retired:`. The multi-forest block of
   [`resolveSamplerSpec`](../../R/spec.R) writes the stated sds, as written, into the configuration's
   `sd`, and a stated `amplitude.prior.variance` into `variance`; the sampler's `initialize` records
   `unit` from the reader where it records `anchor`, at a first creation only. Tests, a new file
   test-forest-sd-unit.R:
   - Which model was fitted. `sd = 2` on a forest with no basis and `sd = 0.5` on a factor forest, at
     each of: a formula term at `dbarts` and at `bart`, a list on a formula, on a matrix and in
     `dbartsSpec`, a data object with `bases` and a list of sds. At every door the reader's `sd` is
     `identical()` to 2 and to 0.5 and `sd.stated` is `TRUE`, `k.scale * 0.674` of the factor forest
     equals 0.5 and the half-Cauchy median times the unit equals 2, both to 1e-12, and the train draws of
     the sampler doors are `identical()` (fails today: `k.scale` is 0.5 sd(y) / 0.674).
   - L written out. With an offset, `subset`, a dropped missing response and weights with zeros, a
     default forest's reader `sd` is `sqrt(2 / K) * sd((y - offset)[kept])` to 1e-12, where `kept` is
     written in the test and leaves out the rows of weight zero (dec-B296); the same reader is
     `identical()` after a later `$setWeights` that brings those rows back and after a
     `$setActiveRows`, the unit being held; and the same four with `sd = 0.5` stated give
     `k.scale * 0.674 == 0.5` to 1e-12, so a unit from the wrong rows cannot hide in a ratio. With the
     rows of positive weight at one value the reader's `sd` is the multiple of the every-row
     `sd(y - offset)`, under exactly one warning of class `dbartsScaleFallbackWarning`, counted with
     `withCallingHandlers`.
   - Probit and logistic: `sd = 0.5` gives `k.scale * 0.674 == 0.5` under both (fails today under
     logistic: 0.907).
   - The size, from the engine's own prior draws under a flat likelihood. Always on, one a family, 1500
     sweeps and 10 percent: a forest with no basis and a factor forest, each drawn and stating `sd = v`,
     have size v (fails today under gaussian and logistic). `at_home`, 5000 sweeps a row and 6 percent:
     the same at three values of v. `tinytest::test_package()` runs with `at_home` false, in CI too, so
     the first is the one a push is gated on. It is seeded; at 2000 sweeps the rerun's entries ran from
     0.92 to 1.03 of v, so 10 percent at 1500 is about three standard errors, and the arm's sweep count
     is raised before its tolerance is widened.
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
   places it is refused with the first text. "Stated" is the record
   [`plainForestStated`](../../R/model.R) keeps of what the caller named. Tests:
   `y ~ forest(x1 + x2, sd = 2)` at `bart` and `dbarts`, and `forests = list(forest(sd = 2))` on a
   formula and on a matrix, are draw for draw `leaf.prior = normal(sd = 2)`, under gaussian and probit
   (fail today: refused); beside `leaf.prior = normal(sd = 3)`, `normal`, `normal(k = 2)`, `bart`'s `k`
   and a linear leaf, refused with the first text; `amplitude` on that forest keeps its refusal.
6. The two other texts. The engine's throw on a response with no spread is raised with the second text,
   naming nothing of the engine; [`resolveLeafPrior`](../../R/model.R) takes the third. Tests: two
   forests on a constant response are created with no sd (as today) and refused with `sd = 0.7` (fails
   today: created, every prior 0); the four pins of
   ["a named leaf-prior 'sd'"](../../inst/tinytest/test-calibration-creation.R), two there and two in
   test-bcf-family.R, take the third text.
7. A named single sd (dec-B280). [`validateForestSd`](../../R/model.R) drops the name of a single
   number and refuses nothing for it; the text for a longer `sd` loses its words about columns. Tests:
   `sd = pars["s"]` on a forest with no basis, on a factor forest and on a numeric one, at a formula
   term, in a list and in `$setLeafPrior(forests = )`, and on the one forest of a single-forest model:
   each is created, the reader's `sd` is `identical()` to `unname(pars["s"])`, and the draws are those
   of the unnamed number (fail today: refused); `c(dose = 1, age = 2)` and `c(1, 2)` keep the length
   text at the three places; the pins of
   retired: ["must not be named"](../../inst/tinytest/test-forest-arguments.R) are turned into these.
8. Respell and repair. Run the suite first and repair what it shows. Expected: the pins of the stored
   numbers by position, about 60 in the 7 files Context lists, each restated against the four numbers
   or against the reader, with the two hand-edited records of test-forest-basis-r5.R; the 17 assertions
   of the unit in 9 files (test-bcf-creation.R 3, test-bcf-family.R 3, test-fit-stores-k.R 3,
   test-calibration-midchain.R 2, test-multiforest-leaf-prior-writer.R 2, test-forest-arguments.R 1,
   test-forest-basis-r5.R 1, test-forest-defaults.R 1, test-forest-selection.R 1), six of them among
   the first; the pins of the reader's `leaf.prior` for a forest that states nothing, which now read
   `sd`; the pins of `extract`'s vector in test-fit-stores-k.R and test-forest-selection.R, which
   becomes a list; and the three pins of
   retired: ["this model has one forest, so its size is"](../../inst/tinytest/test-bcf-creation.R), one there,
   one in test-forest-arguments.R and one in test-formula-terms.R, which become step 5's fits. Each
   gaussian and logistic forest that states an sd, and each of the 34 gaussian writes, is read: where
   the test means a multiple of the response's scale (the pins of the map in test-bcf-family.R, the
   writer's ratios, the prior-draw tests) it is written `r * sd(y)` or `r * pi / sqrt(3)` with the
   scale computed in the test; where it needs only some stated sd it is left, and now states that many
   units. Benchmarks: the four BCF exact gates and `sbc.R` write `sd = sdControl * s` and
   `sd = sdModerate * s`, s the scale their oracles already use, and their oracles are not edited.
9. Help and records. man/forest.Rd: the `sd` item, with the unit and the rows it is read from, what s
   is the size of for each kind as "Context" measured it, the defaults as multiples of the response's
   standard deviation and as the numbers they come to, that a name on the number is ignored, and one
   sentence that a single forest takes it as `leaf.prior = normal(sd = )`; the Details paragraph on the
   budget where it speaks of units.
   [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) and
   [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd) with their docstrings: `leaf.prior`
   is what was stated and goes back as it is; `sd` is the number in force, in the response's units, and
   `sd.stated` whether it was stated; the list is named by label; `amplitude.prior.scale` and
   `leaf.scale.factor` are multiples of the response's scale. The `leaf.prior.sd` item of `extract` in
   man/bartBT.Rd: the list by label, a number for one forest. man/bart.Rd and man/dbarts.Rd:
   `leaf.prior` beside several forests. docs/design/forest-sd-unit.md with its index row: the rule, who
   holds the map, the unit's record, the table of "Context", the changed sequence and its oracle.
   docs/design/multiplier-combiner.md, bcf.md, nameable-calibration.md and public-surface.md where they
   give a forest's sd in units of the latent scale or have R derive the map;
   docs/design/response-scale-rows.md where it says what the unit sizes; docs/architecture.md where it
   lists the forests' configuration. TODO: `forest-prior-args` names this slice landed, and
   `forest-sd-unit` and `named-single-sd` close.
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
  stopping, and at least the base build's count (22978 over the 246 files that run in two threads) plus
  the new file's assertions; the landing note gives the figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-308c45be.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded: no scenario states a forest sd, and a
  scenario that did, or a map that moved with its move to the engine, would show here as a `max |z|`
  line. These three are the gate of the change of shape.
- Every command `.github/workflows/exact-gates.yaml` runs, 42 at the tip, in `quick`: the class is
  posterior-changing and the four respelled BCF gates are the oracle for a stated sd in each channel.
  `sbc.R`'s respelled arm constructs and runs its smallest setting.
- `inst/include/dbarts/dbarts.h` has no diff and `tools/check-api-hash.sh` passes.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy.
- The consumers, each suite whole against a private install of the slice, none failing: bartCause on
  dbarts-1.0 with step 10's edit, stan4bart on bartcore, treatSens on dbarts-1.0, bairrtt on main, the
  last three unedited; the landing note gives each count beside its count on the base build. bartCause
  unedited is run once too and passes, its `sd = 1` now one unit and its tests pinning no draw under it;
  if that run fails, the landing note says where.

The pair script. Each row is one model fitted on the base build as written there and on the slice's as
written here, the same seed, 20 sweeps: train draws, sigma and the coefficients. Rows 01 to 06 and 21 to
28 must be `identical()`; rows 07 to 20 are equal to 1e-10 and the count of identical ones is recorded;
rows 29 and 30 must differ.

| | base build | slice |
|---|---|---|
| 01 to 04 | no sd stated: two forests; three; under `subset` with weights, ten of them zero, and an offset; `bart` with a term | the same text |
| 05, 06 | logistic and probit, no sd stated | the same text |
| 07 to 10 | `sd = r` on the forest with no basis; on a factor forest; on both; on a numeric column | `sd = r * sd(y)` |
| 11 to 13 | r stated with an offset; under `subset`; with weights of which ten are zero | `r * sd(y - offset)`; `r * sd(y[kept])`; `r * sd(y[w > 0])`, the base build's unit being over those rows too (dec-B296) |
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
- the row norm left out of the engine's map: the restated fixtures, and pair rows 10, 24 and 27;
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
  too: step 3's "L written out", and step 1 (a);
- the unit taken over the rows in where they hold one value: step 3's fallback arm, and step 1 (a);
- the unit moved by `$setWeights`: step 3's "L written out";
- `extract` left returning `k.scale` over k; left a vector; named by position; a list of one where one
  forest is asked for: step 4's `extract` and names tests;
- the stated number mirrored into the four numbers and not into `sd`: step 4's path "a write, a copy
  and a second write";
- the both-given refusal dropped: step 5;
- a name kept on the stored sd: step 7's `identical()`.

Not a hot-path change: the map is resolved when a forest is built or restated, never in a sweep.

## Budget and stops

Planned ~2340 lines, by step:

| step | package code | tests | the rest |
|---|---|---|---|
| 1 the engine | engine 300 | tests/cpp 480 | |
| 2 the bridge | bridge 160 | | |
| 3 creation | R 110 | tinytest 260 | |
| 4 the reader, the writer, `extract` | R 120 | tinytest 190 | the mutation battery 10 |
| 5 one forest | R 40 | tinytest 30 | |
| 6 the two other texts | R 15 | tinytest 10 | |
| 7 a named single sd | R 15 | tinytest 20 | |
| 8 respell and repair | | tinytest 200 | benchmarks 40 |
| 9 help and records | | | help 130; design note, architecture, public-surface, TODO, indexes 210 |

Package code ~760, tests ~1190, the rest ~390. The four most recent slices ran 1.8, 2.2, 2.3 and 1.9
times their plans (setState's two forms, the sparse starting sigma, the state's record
of missable columns, the creation slice of response-scale-rows), with 68, 52, 83 and 64 percent of what
they changed in tests. Assumed here: the tests run to about 2.5 times their figure and everything else to
about 1.3, which is ~4470 lines, inside a forecast of 4200 to 5400 (1.8 to 2.3 times). An overrun in
tests is expected and is not a reason to stop; one in package code means the design moved. bartCause's
part, ~10 lines, is outside these counts.

Stop and report when:
- the diff passes ~5400 lines, or the package code (engine, bridge, R) passes ~1140, 1.5 times its
  figure;
- a recorded scenario moves in any of the three compares, or a snapshot file moves;
- a fit that states no forest sd is not bit for bit the base build's on any pair row (01 to 06, 21 to
  28), or a respelled row is apart by more than 1e-10;
- the conversion is needed anywhere but the one function of the engine, or R must divide or multiply by
  the unit on any path;
- the engine must be told more of a forest's kind than whether it has a basis;
- the unit needs a test of which rows are in beyond what the response model answers;
- dbarts.h needs any change, the ABI token moves, or a state block or the stored state changes;
- an existing test fails in a file that steps 5 to 8 do not name: the ten of Context (the seven and
  the nine, which share six), test-formula-terms.R and test-calibration-creation.R;
- bartCause unedited fails, or a consumer's count falls.

Expected verdicts: tests/cpp and tinytest clean; the three compares 55, 15 and 11 identical, no scenario
moved and nothing re-recorded; the four snapshot files unchanged; all 42 exact-gate commands passing, the
four respelled BCF gates among them; the pair script's rows 01 to 06 and 21 to 28 identical, 07 to 20
within 1e-10 and 29 and 30 differing; the consumers' suites clean.

## NEWS

No new item: `forest()` and its `sd` are new in 1.0-0 and nothing released changes.

## What this leaves for the multiplier slice

- What an sd is the size of where the coefficient is held or the basis is numeric: the table of "Context"
  with its rows off 1 is still true after this slice, in the response's units. The help says so, in those
  words, until the multiplier slice changes them.
- The engine after this slice is handed a statement and holds the tip's calibration map in one
  function. The multiplier slice's pushes change rows of that function and add one thing to the
  statement, the forest's kind; they keep `sd` not-a-number for "not stated", the recorded unit, the
  division, the reader's `sd` and `sd.stated`, and the four-number record.
- The defaults dec-B290 and dec-B291 revise (a held forest with no basis; a forest with a basis whatever
  K), unless the open call puts dec-B291 here: `defaultForestSd` is the one place either changes.
- The reader's `amplitude.prior.variance`, `amplitude.prior.scale`, `leaf.scale.factor`,
  `leaf.scale.divisor`, `basis.row.norm` and the words of `prior.sd.of`; a numeric forest's `sd` named
  by its column; the entries for a forest's kind and hold; the printed block: none is built or renamed
  here. A forest of several keeps k at 1 and a `k.scale` of its own, which moves with a stated sd
  (dec-B416 speaks of a forest with no map).
- What the multiplier slice may assume: every stated sd in the tree, the benchmarks and bartCause is in
  the response's units; the reader's and `extract`'s shapes are final but for the column's name; a fit
  of forests with no basis or a two-level factor, coefficients drawn, stated or not, is unchanged by it
  to the bit and can be gated so, at two forests if dec-B291 is built there and at any number if it is
  built here.

## What waits on what

- On forest-defaults-by-kind and written-surface, landed: step 5 reads the record of what the caller
  stated ([`plainForestStated`](../../R/model.R)) and the finder of the forest with no basis
  ([`plainForest`](../../R/model.R)); step 3 keeps the count, base and power `forestParams` gives every
  forest as the first three of the four numbers; the third text stands beside the refusal of a leaf
  prior where no forest is without a basis, which comes first; selection by label, which step 4's
  `extract` by label uses; the interim refusal of held shapes
  ([`refuseHeldShape`](../../R/model.R)), which this slice does not touch; the forests' labels, which
  name the reader's and `extract`'s lists.
- On the creation slice of response-scale-rows, which lands first: the rows the anchor is read from and
  the response model's answer to which rows it read. Its tests pin the anchor's rows; this slice's pin the
  unit's, which is the anchor times the transform's scale.
- Order with the two slices after it: this one, then the kind by class, then the multiplier slice. This
  slice and the kind by class share R files and are serial; neither needs the other. They are not one
  slice: this one changes arithmetic in the engine and a prior, that one is R only. It is not folded
  into the multiplier slice's pushes that change a prior either: that slice's cleanest gate is that
  forests with no basis or a factor are bit for bit under it, stated sds included, and that gate exists
  only if the unit moved first.
- TODO `multiforest-scale-rederive`, after 1.0-0: a re-derived scale on several forests restates each
  default on the new scale and keeps a stated `forest(sd = )` as given, which is what a number kept as
  written and a default kept as a multiple make possible. Nothing of it is built here.

## Out of scope, and where it goes

- The fitting function's `leaf.prior = normal(sd = )` as the size of the forest with no basis in a model
  of several: refused here with the form to write. It belongs with the reader's `normal()` objects and
  `$setLeafPrior`'s single form, in TODO `forest-prior-args`.
- The kind by class, with the refusal of a factor of three or more levels and of several numeric
  columns; the multiplier slice with everything listed above; `leaf.prior = normal(sd = )` on `forest()`
  and `updateBasisScale`, after the merge to main (dec-B276).
- `updateScale = TRUE` on a sampler of several forests: refused, as dec-B362 keeps it for 1.0-0.
- A sampler made by hand from another sampler's control, model and data after its weights changed
  derives the response's range again under the weights in force and keeps the anchor and the unit from
  the control, so its stated sds are off by the ratio of the two ranges: the seam response-scale-rows
  names, which this slice neither opens nor closes; closed by the record slices of TODO
  `state-frame-prior`.
- A way to put a stated forest back to its default: none exists before `normal()` objects are traded.
- To TODO as a new entry: the verbose summary of a model of several forests prints one tree count and
  "k prior fixed to 2", the first forest's model slots, and nothing of the others.

## Open calls

One, opened by the recheck. The steps above are written to option A.

1. Which slice makes the default sd of a forest with a basis 1 whatever the number of forests
   (dec-B291)? The maintainer's words there: "Drop K". The ruling is recorded against
   [forest-multiplier.md](forest-multiplier.md), which says the same holds for this plan wherever it
   states the default; this plan's rule 3 and its first constraint say a default does not move here.
   In a model of three forests on a response whose standard deviation is 10, a forest with a basis
   that states no sd has a default sd of 8.16 today (sqrt(2 / 3) of 10) and 10 under dec-B291; at two
   forests both are 10.
   - (A) The multiplier slice builds it; this slice moves sqrt(2 / K) into the engine as it is. Every
     fit that states no sd stays bit for bit through this slice, which is the gate of its change of
     shape. Costs: the help's sentence on the default, `defaultForestSd` and step 1's literals at K = 3
     are written here and rewritten one slice on; between the two slices the reader reports 8.16.
   - (B) This slice builds it. The default is written once, in the engine and in the help. Costs: rule
     3 and the first constraint change; the slice gains a second posterior-changing sequence, a model of
     three or more forests with a forest that has a basis and states no sd, so its bitwise gate for
     unstated fits holds at two forests only and pair row 02 is expected to differ; no recorded scenario
     moves (every one of several forests has two) and nothing is re-recorded. Measured on the tip with
     the default changed in R alone: 16 assertions fail, all in test-bcf-family.R, each a pin of
     sqrt(2 / K) or of the prior's size at three and four forests; no seeded draw of the suite moves,
     though it creates 220 such models in 13 files.

## Calls made in planning

- The slice stands alone, before the kind by class and the multiplier slice. Its draws do move, for a
  stated sd left as written, so it is not a renaming that could ride with either; the reasons against
  folding it are under "What waits on what".
- The multiplier slice's first push as first planned is folded in here (its critique's finding 7, taken
  by the coordinator). That push handed the engine the statement and moved no draw; this slice as first
  planned added a stated sd beside four derived numbers it kept, a reader that returned a default as a
  number, a skip of an equal write, a reading of records saved before it and `extract` as a vector, and
  the multiplier slice then removed or reversed each and repaired the same pins a second time. Folded:
  one gate battery fewer, the positional pins repaired once, and nothing built to be removed. The cost
  is a review that holds two things, a change of shape that must move no bit and a division that must
  happen once; the bitwise compares gate the first and the unit's own tests the second, and the reviewer
  reads the diff once for each. Not folded: the kind's route to the engine, which the tip's rule does
  not read.
- The engine converts, against a unit it derives and R records. The alternative, R dividing by
  `sd(y - offset)` before the bridge, is shorter and was built as the stand-in. It was not chosen: R's
  estimate differs from the engine's in the last place, so the reader could not return a stated number
  exactly; R would have to repeat the engine's choice of rows, the fallback included; and a `forests`
  list or a control resolved against one data object and run on another would carry the first one's
  unit in silence, where the engine derives the unit of the data it is given.
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
  and for one forest the entry itself (the coordinator, on the critique's question (c)). The multiplier
  slice's plan had a list whatever was selected; `extract(type = "k", forest = 2)` returns a number, and
  one argument should not give a number for one quantity and a list for the other. At the tip `extract`
  is already in the response's units, of the forest's own standard deviation before its coefficient,
  which is neither the stated number nor what dec-B253 and dec-B275 say it returns. dec-B376 approved
  `extract` as `k.scale` over k as built; that is the one forest's computation, which stays.
- The reader's list and `extract`'s are named by label here and not with the multiplier slice: dec-B275
  names them so, selection by label has landed, and a second renaming of the same list two slices on
  would repair its pins twice. A margin of draws keeps `forest<i>`.
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
  so. They are entries the multiplier slice removes; converting them would change two pins to change
  them again.
- A stated sd on a response with no spread is refused. Today it is taken and every prior is 0. Dividing
  by a unit of 0 would hand the engine an infinite scale.
- One forest takes `sd` in this slice, as written-surface's table and the design say it does once the
  unit is one. The fitting function's `normal(sd = )` beside several forests is not served here: the
  design gives that to the slice that trades `normal()` objects, and the text now says where to write it.
- A test's stated sd is respelled only where the test means a multiple of the response's scale. Every
  other one is left to state units; the alternative, respelling all of them, keeps draws nothing pins.
- The pair script is run at landing and not tracked: its old side needs the base build. What stays in
  the suite is step 3's pins of unstated fits and step 4's identities between paths.
- Measured for this plan on stand-ins, not on the design: the 17 assertions, the 9 of 18 pairs and the
  two "left as written" distances come from the tip with the conversion made in R; the positional pins
  are read off the suite's recorded expectations, not run. The engine's conversion divides by a number
  one unit in the last place away, so the counts of identical pairs may differ; the class does not.
  Nothing of the engine's statement is built.
- Recheck. The unit follows the response model's answer to which rows it read, as the anchor does since
  response-scale-rows; this slice tests no row itself. In the fallback a stated sd is divided by the
  every-row unit, under creation's warning and no second one.
- Recheck. A map forest's stated number is reported in the calibration's `namedSd` and the bridge's
  `named.sd` column, which the k-internal slice added for a forest with no map. Not taken: a third new
  field, which would give one thing two names.
- Recheck. The four numbers keep the hold though the record now holds `heldAmplitude` too: the bridge
  reads the numbers and a fit reads the other, and the two sibling plans count on four. Not taken: three
  numbers, the bridge reading `heldAmplitude`.
- Recheck. The prior-size arm that always runs keeps its 10 percent and takes more sweeps if its seed
  sits near the edge. Not taken: a wider tolerance, which lets more of a wrong unit through.

## Recheck against bartcore

On 2026-10-11, against d5b637ee: the notes of what landed since this plan was written were read
([forest-defaults-by-kind.md](forest-defaults-by-kind.md), [written-surface.md](written-surface.md),
[k-internal-parameter.md](k-internal-parameter.md), [default-rule-per-gap.md](default-rule-per-gap.md),
[setstate-force-update.md](setstate-force-update.md), [sparse-starting-sigma.md](sparse-starting-sigma.md),
[response-scale-rows.md](response-scale-rows.md) with its design note), then the code, the ledger for
every ruling cited here and those made since, and the runs below on a shipped build of that tip. The
text above is written to that code throughout.

What changed:

- What this plan waited on is in the tip, and its "to recheck once both have landed" list is done: the
  counts of Context are today's, no baseline scenario has gained a stated sd,
  [`resolveForests`](../../R/model.R) still holds the one-forest refusal step 5 removes, and every
  function cited exists under its name.
- The unit's rows. The anchor is read over the rows in the likelihood at creation, with the every-row
  fallback and its warning; Context's L, rule 1, the constraint on rows, step 1 (a) and (e), step 3's
  "L written out" with its fallback and `$setActiveRows` arms, pair rows 03 and 13 and three mutations
  say so. Rule 8 names the weights and a mask among what moves no unit.
- The suite: 253 files and 22978 results where the plan had 231 and 16065. 1427 models of several
  forests in 66 files where it had 776 in 53; 59 forests stating an sd in 8 files (51 in 6); 39 written
  in 5 files (33 in 3). The unit's stand-in fails 17 assertions in 9 files (15 in 7):
  test-forest-defaults.R and test-forest-selection.R are new. The positional pins are about 60 in 7
  files (72 in 4), counted another way.
- The pairs: 9 of 18 identical and the rest within 2.8e-14 (8 of 18, within 9.6e-13); the row of two
  numeric columns given as a matrix, which had the widest gap, is no longer a call the tip takes and was
  replaced by a row with ten weights of zero, which is identical.
- The record and its readers: `forestParams` takes the fitting function's count, base and power for the
  forest with no basis; the record holds `labels`, `heldAmplitude`, `basisTerms` and `control.n.trees`;
  `$setForestBasis` reads the hold by position (step 3). The reader's bridge matrix has `map.anchor`,
  `k.prior.scale` and `named.sd`, the last reused (steps 1 and 2).
- Names: the baseline is `equivalence-308c45be.rds`; the exact-gates workflow runs 42 commands; the
  `leaf.prior.sd` item of `extract` is in man/bartBT.Rd; the cite to retire in feature-matrix.md is
  `ForestSpec::amplitudePriorScale`; the refusal of a held numeric column is dec-B353's; the consumers'
  tips. The names this plan proposes for the engine's function and its test are `forestMapEntry` and
  `testForestMapStatements`.
- The budget: ~2340 planned for ~2250, a forecast and stops where there was an upper figure of 3900.
- dec-B291, ruled after the plan was amended, is the open call. dec-B290 and dec-B293 (a held forest
  with no basis) stay the multiplier slice's and touch no rule here.

What did not change: every rule but the words of 1, 3 and 8 named above; the refused texts; the
constraints but the one on rows; the steps' order and number; the gate expectation that no recorded
scenario moves. The rulings as they stand: no entry made since the plan was amended changes what
dec-B246 as dec-B266 revises it, dec-B253, dec-B269 as dec-B280 and dec-B282 revise it, and dec-B275
rule here (dec-B354 and dec-B290 refine dec-B253 for a held coefficient, the multiplier slice's);
dec-B296, dec-B302 and dec-B309 are built as the plan had restated them; dec-B362 and TODO
`multiforest-scale-rederive` read right beside rule 8; dec-B416 to dec-B418 concern a forest with no
map and a state, neither of which this slice touches; dec-B376's record of a held coefficient is the
`heldAmplitude` entry.

Runs: L with an offset, `subset`, positive weights, 50 and 10 rows of weight zero, the weights brought
back, one row in, both latent families, `setResponse` without and with `updateScale`; the prior-size
table at 2000 sweeps a row and one row three times at 5000; the reader, `extract` by position, by
`forest<i>`, by label and by two labels, the one-forest, named-sd, leaf-prior and constant-response
cases; the suite three times in one process, traced, under the unit's stand-in, and with dec-B291's
default; the 27 pairs. The stand-ins were made by replacing R functions in the loaded namespace of the
tip's build: nothing was compiled.

Not verified: nothing of the slice is built, so no mutant and no new test was run, and the change of
shape was run on no build: its count of pins is read, not measured. No reference build, bitwise compare,
exact gate, sanitizer, x86 run or `R CMD check` was made: that no recorded scenario moves rests on the
search of the three scripts and the four snapshot files. The consumers were searched at the tips named
and not built or run. The 5000-sweep table was not rerun whole.

## Landing note

Landed 2026-10-11 as e84c0dd0 to 37b59ea0, 21 commits on 5c3bcdc2, with bartCause's lockstep commit 70ad798
on its dbarts-1.0 branch (`bcf()`'s `sd.moderate` defaults to NULL, and its two hand-built samplers drop
`sd = 1`), pushed with it. An opus implementer built the engine, the bridge and the R code (41fc386c,
ba9f4a15, 9eda35d0, 5bb61ea7, dc34b927, 2de0eefe, 0787f8b1, 49ee1e44); a sonnet implementer the restated
tests, the BCF exact gates and sbc.R, the help and NEWS, the pin of summary's line and the design note
(a8e8684f, 4a7b7a77, 8ba72184, 612c61a8, 77562f66, e88b78eb); an opus review with 27 mutants; and a fix round
(f30f0920, 22f5acf9, 8b87950e, c7ed6d60, 37b59ea0).

Built as planned: the engine divides a stated number once, in `forestMapEntry`, by a unit it derives over the
rows in the likelihood at the first creation, and R divides nothing; the reader returns the statement with
`sd` and `sd.stated` beside it; `extract(type = "leaf.prior.sd")` is a list by label; the one forest of a
single-forest model takes `sd` as its `leaf.prior = normal(sd = )`, and the two stated together are refused;
a named single sd is taken; a stated sd on a response with no spread is refused. Departures, the calls of
dec-A208: the record keeps `hasBasis` and the engine's spec `createdWithoutBasis`; a stated sd equal to the
default's number takes the default's multiple; two engine branches no R call reaches are kept and tested in
tests/cpp; the three bcf-exact gates state `sd = sdControl * sd(y)`, the reviewer confirming that sd(y) and
not sd of the scaled response is the oracle's number to 2e-16.

The review found no defect in the shape or in the division (55 of 55 pairs identical to the build before,
every path checked per chain against a unit computed in R, no seventh path) and five findings, all closed:
the help called a forest's sd and `normal(sd = )` "not the same thing" where they are one statement; the
`setForestBasis` path for a stated sd hung on one expectation; two engine branches (the in-force shortcut in
`Chain::setForestMapSd`, the unit derived from the handed anchor) no test needed; each of the two bit-exact
rules was pinned in one tier only, the default-number rule by tests/cpp and the reader's return of a stated
number by tinytest; and, minor, a restated assertion that passed under two mutants, help that said
`extract(type = "k")` is named by label, and no test of the prior per chain after a write. The fix round
added a tests/cpp case and a tinytest arm for the basis swap, one tests/cpp case for each branch, a
tinytest that stretches the response until `sqrt(2 / 3) * u / u` differs from `sqrt(2 / 3)` and tests/cpp
values that the unit and back do not return, a pin of a stated variance with the forest with no basis
second, eight reader assertions and a two-chain arm. The eight mutants it targeted fail where they should
(E04, E05, E07, E08, E11, E14, E16 and B02); no mutant was run on the rebased tip.

Gates at the landed tree (tip 00764fc7 of the combined re-gate, whose tree equals 37b59ea0's but for the
docs of 5c3bcdc2): tests/cpp 397 ok plain and under ASAN+UBSAN; full tinytest 24124 results, 0 failed, 252
files; the three bitwise compares 55 of 55, 15 of 15, 11 of 11; the four snapshot files 27 / 0; the 42
exact-gate commands at `quick`, all exit 0; the composition matrix 269 checks and no disagreement; lintr,
air, rc-codoc, win-drift, doc-freshness and verify-anchors clean; `R CMD check --as-cran --no-manual` one
NOTE (the Date field). The interplay probe of this slice with `updateScale = TRUE` on `setWeights`,
`setActiveRows`, `setResponse` and `setOffset` ran 220 comparisons: a stated sd of 2 stays 2 while
`k.scale` is re-derived, several forests refuse and are left untouched, and nothing differed in the
combination; the four unexplained differences are a copy after `setResponse`, which this slice does not
cause (TODO copy-after-response-swap-ulp). The R-loaded ASAN run of the review's tip (14 files, 3224
results, 240 refusals, no report) was not repeated at the landed tip.

The x86-64 leg at 00764fc7: tests/cpp 397 ok plain and under ASAN+UBSAN with no sanitizer report; tinytest
24157 results, 0 failed; the main equivalence compare statistically indistinguishable on 55 scenarios.
bartCause's suite on the slice with its edit has the same 7 failures, with the same numbers, as on the
build before; unedited, `sd.moderate = 1` now means one response unit and not one sd(y), which is why its
commit was pushed with the landing.

Size against the plan: the plan estimated ~2340 lines and forecast 4200 to 5400, stopping at ~5400 or ~1140
of package code. The landed range is +3314 -1014 over 53 files, 4328 changed lines: package code 611 added
(src 455, R 156), tinytest 1161 added against ~710 planned, tests/cpp 763 against ~480, benchmarks 39, man
97, docs 608. Inside the forecast and the stops; the tests ran past the plan because the review added the
cases above.

Left: dec-B291's default of one unit of the response's scale for a forest with a basis, whatever the number
of forests, is built in the forest-multiplier slice, and until it lands the reader reports sqrt(2 / K) times
the unit at three forests or more. What `k` and `k.scale` report for a forest in a model of several forests
(today k = 1 for each forest and `k.scale` the forest's own prior sd) is decided with forests-shared-total.
`fit$leaf.prior` is still named `forest<i>` where the reader and `extract` use labels. The consumers' suites
other than bartCause's, pkgdown and the manual were not run on the landed tip.
