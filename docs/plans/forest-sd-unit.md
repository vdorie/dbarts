# forest-sd-unit: a forest's sd is stated and reported in the response's units

Status: PLANNED (dec-B246 as revised, dec-B253, dec-B266). Follows
[forest-defaults-by-kind.md](forest-defaults-by-kind.md) and push 3 of
[written-surface.md](written-surface.md), neither of which has landed.

agent: one push. Opus implementer for the engine, the bridge and the R code; sonnet for the respelled tests,
the benchmark scripts and the help once the code is fixed; opus reviewer told to refute. The reason for
opus: the slice is small and every slip in it is a prior off by a factor of sd(y) with no message. The
conversion is one division, and it can be made twice, or not at all, on any one of six paths (creation at
four doors, the writer, a re-creation); the review of written-surface's second push found three such
one-path faults in a slice that touched no arithmetic at all.
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- NEUTRAL for every fit that states no forest `sd`, and for every fit under probit, where the unit is 1.
  Measured on a stand-in for the slice: 5 of 5 such pairs identical, and no seeded draw of the suite moved.
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
- Accepted where refused: `sd` on the one forest of a single-forest model.
window: pre-release, after forest-defaults-by-kind and before the kind by class and the multiplier law, so
that every later test, pinned string and help page states an sd once, in the unit it keeps. Serial with
any other work in [`forestParams`](../../R/model.R), the multi-forest block of
[`resolveSamplerSpec`](../../R/spec.R), [`reportLeafPrior`](../../R/dbarts.R),
[`writeForestSpreads`](../../R/dbarts.R), [`applyForestAttributes`](../../src/R_interface_bartcore.cpp),
[`ForestSpec`](../../src/bartcore/combiner.hpp) or the K-forest constructor of
[`Chain`](../../src/bartcore/chain.hpp); so serial with [leaf-conversions.md](leaf-conversions.md) and
[cross-family-state-install.md](cross-family-state-install.md), in either order. bartCause's edit lands
the same day.
budget: ~1300 lines changed (engine ~90; bridge ~60; tests/cpp ~140; R ~170; tinytest ~470, of which
~350 in one new file and ~120 changed in seven; benchmarks ~45; help ~110; design note, architecture,
public-surface, TODO and the two indexes ~180), upper figure 2500 (engine 180, bridge 120, tests/cpp 250,
R 350, tinytest 900, benchmarks 80, help 220, records 300). The design estimated 800 and planned for
1500; written-surface's two landed pushes ran at 1.8 and 2.7 times their plans, and this slice has less
grammar and more arithmetic than either.

## Goal

A forest's `sd` is a number of the response's units wherever it is written or read: `forest(sd = )` at
creation, `$setLeafPrior(forests = )`, the `leaf.prior` entry `$getLeafPrior()` returns, and
`extract(type = "leaf.prior.sd")`. The units are those of y for a continuous response and of the latent
index under probit and logistic, where the link's own error has standard deviation 1 and pi / sqrt(3). A
default does not move: it is the same multiple of the response's standard deviation as before, and is
reported as the number of response units that comes to. With one unit at every place, the one forest of a
single-forest model takes `sd` too, as its `leaf.prior = normal(sd = )`.

## Context

Measured at the tip (4f2e79e6; written-surface pushes 1 and 2 landed, push 3 and forest-defaults-by-kind
not) on the shipped build, R 4.6.1. The spellings below are the tip's, with a tilde on a basis written as
code. L is the response's scale as the engine holds it; "size" is the prior median of the absolute
multiplier times the prior standard deviation of the forest's own fit, exact where the coefficient is
held.

- The suite. 231 files, of which 4 exit off a reference build and 3 ask for more than two threads; the
  other 224 give 16065 results, none failing, 102 seconds in one process.
- L. Under gaussian it is the sample standard deviation of the response net of its offset, n - 1 in the
  divisor, unweighted, over the rows the sampler holds: 3.26736095623188 from the engine against R's
  3.26736095623189 for `sd(y)`, one unit in the last place apart; `sd(y - offset)` with an offset;
  `sd(y[subset])` under `subset`; unchanged by weights, and rows of weight zero count. Under probit it is
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
  | factor, drawn, `sd = 0.7`; a logical; a character of 3 levels | the same | 0.948; 1.018; 0.994 | 1.006 | 0.990 |
  | factor, drawn, default, K = 3 (s = 0.8165) | the same | 1.011 | | |
  | factor, held, `sd = 0.7`; default | (a2 - a1) F | 1.483; 1.494 | 1.494; 1.498 | 1.474; 1.477 |
  | one numeric column w, drawn, default; `sd = 0.7` | a F, times the median nonzero abs(w) | 0.713; 0.691 | 0.701; 0.712 | 0.736; 0.693 |
  | the same column times 1000; times 0.001; an age (mean 50) | the same | 0.689; 0.665; 0.711 | | |
  | a 0/1 number at 50, 20, 5 percent treated | a F | 0.719; 0.700; 0.705 | | |
  | two numeric columns, default; `sd = 0.7` | each a_j F, times the median row norm | 0.69, 0.71; 0.69, 0.69 | | |
  | `cbind(1 - z, z)` as numbers, drawn, default | (a2 - a1) F | 0.975 | | |
  | one numeric column, held | | refused (dec-A171) | refused | refused |

  So s is in units of L. It is the size of the forest's contribution where the forest has no basis or a
  factor and its coefficient is drawn. It is not where the coefficient is held (a held forest with no
  basis has size L whatever s is; a held factor 1.48 s L) or the basis is numeric (0.71 s L at a row whose
  basis has the median nonzero norm). Those are the multiplier law's to change; this slice changes the
  unit of s and nothing in that table but its last factor.
- What reads and writes the number.

  | place | today | unit |
  |---|---|---|
  | `forest(sd = s)` at creation; `$setLeafPrior(forests = list(forest(sd = s)))` | s goes to the engine as given ([`forestParams`](../../R/model.R), [`writeForestSpreads`](../../R/dbarts.R)) | of L |
  | `$getLeafPrior(f)$leaf.prior` | `forest(sd = s)`, s the engine's own number | of L |
  | `$getLeafPrior(f)$k.scale` | the forest's own prior standard deviation, before its coefficient: L for a forest with no basis, s L / 0.674 over the row norm for one with a basis | the response's |
  | `amplitude.prior.scale`, `leaf.scale.factor` entries | s again | of L |
  | `extract(fit, type = "leaf.prior.sd")` | `k.scale` over k, a vector named `forest1`, `forest2`: 1.7408 and 2.5828 for a default two-forest fit whose sd(y) is 1.7408 | the response's, of another quantity |
  | `print` of a fit, `show` of a sampler, the verbose summary | no forest's sd is printed | |
  | `leaf.prior = normal(sd = s)` on a single forest | the prior standard deviation of the forest's total, exactly | the response's |

  The reader and `extract` already disagree: for one default factor forest the reader says 1 and `extract`
  2.58.
- One forest. `y ~ forest(x1 + x2, sd = 2)` and `forests = list(forest(sd = 2))` are refused at `dbarts`
  and `bart` ("this model has one forest, so its size is the fitting function's leaf.prior = normal(sd = ),
  not forest(sd = )"). `leaf.prior = normal(sd = 2)` beside several forests is refused ("a multi-forest
  model does not support a named leaf-prior 'sd': the leaf prior's 'sd' is not a forest's 'sd'");
  `normal(k = 3)` beside them is refused and a bare `normal` accepted.
- A response with no spread. Two forests on a constant response are created with a warning; L is 0 and
  every forest's prior standard deviation is 0, with `sd = 0.7` stated or not.
- Where a forest's sd is stated, counted by running the suite on a build that logs it: of 776 models of
  several forests created in 53 test files, 43 forests state one, in 6 files (test-bcf-family.R 14,
  test-bcf-creation.R 8, test-multiforest-leaf-prior-writer.R 8, test-forest-basis-r5.R 7,
  test-forest-arguments.R 5, test-formula-terms.R 1): 20 under gaussian, 6 under logistic, 17 under
  probit. `$setLeafPrior(forests = )` writes one 33 times in 3 files, 28 of them under gaussian. In
  benchmarks/R the four BCF exact gates and `sbc.R` state `sd = sdControl` and `sd = sdModerate`, 10
  lines; `bcf-equivalence.R`, `equivalence.R` and `multinomial-equivalence.R` state none. No vignette and
  no help example states one.
- What moves in the suite, measured on a stand-in (the conversion done in R against `sd(y - offset)`): no
  file stops and 15 assertions fail in 7 files: pins of the stored numbers behind a stated sd (5), of the
  reader's `leaf.prior`, `leaf.scale.factor` and `amplitude.prior.scale` (5), of `k.scale` against a
  stated multiple (2) and of `extract` (3).
- The same model on both builds, the stand-in's side respelled as s L: with no sd stated, and under
  probit, identical; of 18 gaussian and logistic rows (each kind of forest, an offset, `subset`, weights,
  three forests, a formula term, the writer, a copy, a value read and written back, `bart`) 8 identical
  and 10 apart by at most 9.6e-13 after 20 sweeps. Left as written, `sd = 0.7` on both forests gives draws
  2.31 apart under gaussian and 3.63 under logistic.
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) passes `sd = sd.control`, `NULL` by default, and
  `sd = sd.moderate`, 1 by default, and two of its test files build the same sampler by hand with
  `sd = 1` on the treatment forest. At K = 2 the default of a forest with a basis is exactly 1, so
  `sd = 1` and no sd are one fit today. stan4bart (bartcore, a9d081b), treatSens (dbarts-1.0, aecec71) and
  bairrtt (main, 3f57f61) state no forest sd.

## The rule

1. The unit. A forest's `sd` is in the response's units: the units of y under a gaussian response, and of
   the latent index under probit and logistic. One number of the sampler, the unit, converts it: L in the
   response's units, which is `sd(y - offset)` over the rows kept under gaussian, 1 under probit and
   pi / sqrt(3) under logistic.
2. Stated. The engine is handed the number as written and divides it by the unit; what it then does with
   the quotient is what it does today with a number in units of L. Nothing else of the law changes.
3. Not stated. The default is what it is today, a multiple of L (2 or 1 for a forest with no basis,
   sqrt(2 / K) for one with a basis), and never passes through the unit on the way in.
4. Read. `$getLeafPrior(f)$leaf.prior` is `forest(sd = v)` with v in the response's units: the stated
   number itself, bit for bit, or the default's multiple times the unit. It goes back into
   `$setLeafPrior(forests = )` as it is, and a write of the value read changes nothing, to the bit.
   `extract(type = "leaf.prior.sd")` on a fit of several forests returns the same numbers, in its shape
   of today.
5. The unit is derived once, by the engine, at a first creation, and recorded beside the anchor. Every
   re-creation is handed the record: a copy, a reload, a sampler made again from its own control, model
   and data. No later change of the response moves it, as none moves the anchor.
6. A control or a `forests` list carried to other data states the same numbers of that data's units: the
   record is not carried, and the new sampler derives its own.
7. One forest. `forest(sd = s)` on the one forest of a single-forest model is that model's
   `leaf.prior = normal(sd = s)`, in every family, with that family's own verdict on a stated sd. Stated
   at both places it is refused by name.

Before and after, per door. "Several" is a model of several forests.

| door | before | after |
|---|---|---|
| a `forest()` term of a formula at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()`; a list of sds over a data object's `bases` | s times L | s response units; one model at every door, by one conversion in the engine |
| `$setLeafPrior(forests = )` | s times L | s response units |
| `$getLeafPrior()`'s `leaf.prior` | `forest(sd = )` in units of L | in the response's units |
| its `k.scale`, `response.scale`, `response.shift` | the response's units | unchanged |
| its `amplitude.prior.scale`, `leaf.scale.factor` | units of L | unchanged, and said to be so; the multiplier law removes them |
| `extract(type = "leaf.prior.sd")`, several | each forest's own standard deviation before its coefficient | each forest's sd as the reader gives it |
| the same on one forest; `extract(type = "k")` | | unchanged |
| `print`, `show`, the verbose summary | no forest's sd | unchanged; the printed block is the multiplier law's |
| `$setForestBasis` | the forest's own scale is derived again from the new block's row norm, under the number in force | unchanged, the number in force being the stated one or the default |
| `predict`, `fitted`, a state stored or installed | do not read a prior | unchanged |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | built from the stored numbers and the anchor | built from the stored numbers, the anchor and the unit: the same prior, bit for bit |
| a sampler saved before the slice | | runs under the prior it had; its reader reports that prior in the response's units |
| `sd` on the one forest of a single-forest model | refused | the fit with `leaf.prior = normal(sd = )` |
| the fitting function's `normal(sd = )` beside several forests | refused | refused, in words that say where the size is written |

## Refused forms, with their texts

Base R's style, as in written-surface.

    this model has one forest and states its size twice: 'sd' on the forest and 'leaf.prior' given to the fitting function are one statement; give one            ('k', and the retired 'node.prior', named as written)
    a forest's 'sd' is in the response's units, and this response has no spread to state one against (its standard deviation is 0); leave 'sd' out, or give a response that varies
    a model of several forests states each forest's size on the forest: write the size of the forest with no basis as forest(x1 + x2, sd = 2), which is in the response's units as leaf.prior = normal(sd = ) is, and leave 'sd' out of 'leaf.prior'

Behind them, from the bridge: "the forests' unit must be a positive finite number", beside the anchor's
text, and "the forests' 'sd' must hold one number, or NA, for each forest". Gone: "this model has one
forest, so its size is the fitting function's leaf.prior = normal(sd = ), not forest(sd = )". Replaced
by the third: "a multi-forest model does not support a named leaf-prior 'sd'", whose reason, that the two
sds are different things, stops being true. Kept as they are: every text of
[`validateForestSd`](../../R/model.R), none of which names a unit; "must be a single number ... for every
column of its basis" stays without "per unit of each column" until a numeric forest's sd is that.

## Constraints

- A fit that states no forest sd draws what it draws on the base build, to the bit, in every family. No
  baseline is re-recorded and no snapshot file regenerated.
- The law is not touched: which channel a forest's number reaches, the 0.674, the coefficient's variance,
  the row norm, the held values, the defaults and the K in them. `amplitude.prior.variance` stays.
- The conversion is made in one place, the engine, against one recorded number. R divides nothing and
  multiplies nothing: it hands over what was written and reports what the engine returns.
- A stated number is kept as written, so that the reader returns it bit for bit.
- The shapes of the reader and of `extract` do not change; `sd.stated`, the per-column entries, the list
  by forest and the printed block are the multiplier law's.
- The eight numbers per forest keep their order. From this slice the fourth and seventh hold the default's
  multiple of L for a sampler created with it, and a stated sd rides beside them.
- The flat C header, the stored state and every state block are untouched.
- A `forests` list written before the slice is still accepted, with every number read in the new unit.
- Each commit leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the
reader's sake. Calls are in push 3's spelling of a basis; the fixture is 150 rows with x1, x2, a 0/1 z, a
factor zf, dose, an offset column and weights of which ten are 0.

1. The engine: a unit, and a stated sd. [`AmplitudeSpec`](../../src/bartcore/combiner.hpp) gains `unit`,
   not a number where the engine is to derive it; [`ForestSpec`](../../src/bartcore/combiner.hpp) gains
   `sd`, not a number where none is stated. The K-forest constructor of
   [`Chain`](../../src/bartcore/chain.hpp) takes the unit it is handed, or derives it beside the anchor:
   under gaussian the anchor times the scale of the response transform at that moment, under probit and
   logistic the anchor. A forest that states an sd has `sd / unit` where its number in units of L stands
   today: the half-Cauchy median of a forest with no basis, the leaf-scale factor of one with a basis.
   The expressions of [`mapLeafScale`](../../src/bartcore/chain.hpp) keep their order of operations, so
   a forest that states none keeps its bits. A stated sd with a unit that is not positive and finite, and
   an sd that is not, throw. [`ForestCalibration`](../../src/bartcore/chain.hpp) gains the unit and the sd
   in force in the response's units: the stated number itself, or the multiple of L times the unit.
   [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) takes the response's units, skips a write equal
   to the sd in force as reported, and otherwise divides as creation does. Every struct here is read by
   objects that do not track headers: `--preclean`.
   Tests, tests/cpp, a new `testForestSdUnit` beside
   [`testForestMapWriters`](../../tests/cpp/test_sampler.cpp): one fixture of three forests (none, a
   two-level block, a column of numbers) per family, with a response whose standard deviation is neither
   1 nor its range. (a) The derived unit equals the literal (the fixture's `sd`, 1, pi / sqrt(3)) to 4
   ulp. (b) Each forest stating v has the leaf scale, or the half-Cauchy median, that the base spelling
   gives at v over the literal unit, to 4 ulp; a forest stating none has the base spelling's bits. (c) The
   reader returns v itself, and a default as its multiple times the unit. (d) A write of the value read
   leaves 20 sweeps identical to an untouched twin's; a write of another value gives the leaf scale of a
   chain created with it. (e) A chain made again over the response times 10 plus 5 with the first one's
   anchor and unit handed back has every leaf scale of the first, bit for bit. (f) The two throws.
2. The bridge. [`applyForestAttributes`](../../src/R_interface_bartcore.cpp) reads two more entries of
   the forests' configuration: `sd`, one number or NA for each forest, and `unit`, as it reads `anchor`,
   with the two texts for a malformed one. [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)
   sets a stated forest's sd from the first. [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp)
   reports the unit and the sd in force as two more columns.
   [`bartcore_setForestSd`](../../src/R_interface_bartcore.cpp) hands its number on unchanged. A
   configuration with neither entry is a sampler saved before the slice: its fourth and seventh numbers
   are multiples of L whether stated or not, are read as today, and give the prior it had.
3. Creation. [`forestParams`](../../R/model.R) puts the default's multiple in the fourth and seventh
   numbers whatever the forest states; the multi-forest block of
   [`resolveSamplerSpec`](../../R/spec.R) writes the stated numbers, as written, into the configuration's
   `sd`; the sampler's `initialize` records `unit` from the reader where it records `anchor`, at a first
   creation only. Tests, a new file test-forest-sd-unit.R:
   - Which model was fitted. `sd = 2` on a forest with no basis and `sd = 0.5` on a factor forest, at
     each of: a formula term at `dbarts` and at `bart`, a list on a formula, on a matrix and in
     `dbartsSpec`, a data object with `bases` and a list of sds. At every door the reader's `leaf.prior`
     is `identical()` to `forest(sd = 2)` and `forest(sd = 0.5)`, `k.scale * 0.674` of the factor forest
     equals 0.5 and the half-Cauchy median times the unit equals 2, both to 1e-12, and the train draws of
     the sampler doors are `identical()` (fails today: the reader agrees and `k.scale` is 0.5 sd(y) /
     0.674).
   - L written out. With an offset, `subset`, a dropped missing response and weights with zeros, a
     default forest's reader sd is `sqrt(2 / K) * sd((y - offset)[kept])` to 1e-12, where `kept` is
     written in the test and includes the rows of weight zero; and the same four with `sd = 0.5` stated
     give `k.scale * 0.674 == 0.5` to 1e-12, so a unit from the wrong rows cannot hide in a ratio.
   - Probit and logistic: `sd = 0.5` gives `k.scale * 0.674 == 0.5` under both (fails today under
     logistic: 0.907).
   - The size, from the engine's own prior draws under a flat likelihood (`at_home`, 5000 sweeps a row):
     a forest with no basis and a factor forest, each drawn and stating `sd = v`, have size v within 6
     percent under each family (fails today under gaussian and logistic).
   - Not stated: a forest that states no sd has the engine's own multiple untouched, `identical()` to
     `sqrt(2 / K)` where it has a basis and to 2 (1 under a latent family) where it has none, at K = 2
     and 3 under each family: a default sent through the unit and back rounds, and fails this. The
     draws of unstated fits are the pair script's and the three bitwise compares' to hold.
4. The reader, the writer and `extract`. [`reportLeafPrior`](../../R/dbarts.R) builds `leaf.prior` from
   the sd in force as the bridge reports it; [`resolveForestSpreads`](../../R/dbarts.R) is unchanged and
   [`writeForestSpreads`](../../R/dbarts.R) writes the stated number to the engine and into the
   configuration's `sd`, leaving the eight numbers alone; [`extractParameter`](../../R/generics.R)
   returns, for a fit of several forests, each forest's sd from the stored reader's `leaf.prior`. The
   three entries of the mutation battery that quote those lines move with them
   (["m45"](../../benchmarks/R/mutation-battery.R), m46, m47). Tests, same file:
   - One path at a time, each asserting the stated number from the reader and `k.scale * 0.674` (or the
     half-Cauchy median times the unit) against it: after creation; after `$setLeafPrior(forests = )`;
     on `copy()`; after `saveRDS`, `readRDS` and a first use; on
     `new("dbartsSampler", control, model, data)` from the sampler's own three; after
     `setResponse(10 * y + 5, updateScale = FALSE)` and on a copy made then; after a write, a copy and a
     second write; after `$setForestBasis` to another factor of the same levels, which derives the
     forest's own scale again. A unit applied twice gives 0.5 / sd(y), and one not applied 0.5 sd(y), on
     that path alone.
   - Read then write: `s$setLeafPrior(forests = lapply(s$getLeafPrior(), function(p) p$leaf.prior))`
     after 5 sweeps leaves the next 20 `identical()` to an untouched twin's, on a stated and on a default
     sampler.
   - Other data: a `forests` list with `sd = 2`, and a control taken from a sampler, given to a response
     ten times as wide: the reader says 2 and `k.scale * 0.674` is 2 on the new sampler (a carried unit
     gives 20 or 0.2).
   - Saved before: a sampler whose configuration is rewritten by hand to the old shape (the stated
     multiples in the fourth and seventh numbers, no `sd`, no `unit`) runs draws `identical()` to the
     same model created fresh, and reads the same sd to 1e-12.
   - `extract(fit, type = "leaf.prior.sd")` on a `bart` fit of three forests equals the reader's three
     numbers, by position and by `forest = "forest2"` (fails today: 1.74 and 2.58 where the reader says 2
     and 1); on one forest it is what it is today.
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
7. Respell and repair. Run the suite first and repair what it shows; expected from the stand-in, 15
   assertions in 7 files (test-bcf-creation.R 3, test-bcf-family.R 3, test-fit-stores-k.R 3,
   test-calibration-midchain.R 2, test-multiforest-leaf-prior-writer.R 2, test-forest-arguments.R 1,
   test-forest-basis-r5.R 1), and the three pins of
   ["this model has one forest, so its size is"](../../inst/tinytest/test-bcf-creation.R), which become
   step 5's fits. Each of the 26 gaussian and logistic forests that state an sd, and the 28 gaussian
   writes, is read: where the test means a multiple of the response's scale (the pins of the map in
   test-bcf-family.R, the writer's ratios, the prior-draw tests) it is written `r * sd(y)` or
   `r * pi / sqrt(3)` with the scale computed in the test; where it needs only some stated sd it is left,
   and now states that many units. Benchmarks: the four BCF exact gates and `sbc.R` write
   `sd = sdControl * s` and `sd = sdModerate * s`, s the scale their oracles already use, and their
   oracles are not edited.
8. Help and records. man/forest.Rd: the `sd` item, as push 3 leaves it, with the unit, what s is the size
   of for each kind as "Context" measured it, the defaults as multiples of the response's standard
   deviation and as the numbers they come to, and one sentence that a single forest takes it as
   `leaf.prior = normal(sd = )`; the Details paragraph on the budget where it speaks of units.
   [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) and
   [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd) with their docstrings: the
   `leaf.prior` entry is in the response's units and goes back as it is; `amplitude.prior.scale` and
   `leaf.scale.factor` are multiples of the response's scale. The `extract` page's `leaf.prior.sd` item.
   man/bart.Rd and man/dbarts.Rd: `leaf.prior` beside several forests. docs/design/forest-sd-unit.md with
   its index row: the rule, the unit's record, the table of "Context", the changed sequence and its
   oracle, what a sampler saved before reads as. docs/design/multiplier-combiner.md, bcf.md,
   nameable-calibration.md and public-surface.md where they give a forest's sd in units of the latent
   scale; docs/architecture.md where it lists the forests' configuration. TODO: `forest-prior-args`
   names this slice landed.
9. bartCause, same day (its own commit on dbarts-1.0). `bcf()`'s `sd.moderate` defaults to `NULL` in both
   signatures of R/bcf.R, and man/bcf.Rd says of `sd.control` and `sd.moderate` that each is the forest's
   size in the response's units, `NULL` taking dbarts's default. The two hand-built samplers in
   tests/testthat/test-14-bcf.R and test-03-responseFit.R drop `sd = 1`. Both edits are the same fit to
   the bit on either side of this slice, a default being the same number as `sd = 1` at K = 2 before it.
10. Mutations (Verification): apply each, install with `--preclean` where it is the engine's, run the
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
  scenario that did would show here as a `max |z|` line.
- Every gate `.github/workflows/exact-gates.yaml` lists, in `quick`: the class is posterior-changing and
  the four respelled BCF gates are the oracle for a stated sd in each channel. `sbc.R`'s respelled arm
  constructs and runs its smallest setting.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy.
- The consumers, each suite whole against a private install of the slice, none failing: bartCause on
  dbarts-1.0 with step 9's edit (1412 expectations at its last run), stan4bart on bartcore (582),
  treatSens on dbarts-1.0 (306), bairrtt on main (207), the last three unedited. bartCause unedited is run
  once too and passes, its `sd = 1` now one unit and its tests pinning no draw under it; if that run
  fails, the landing note says where.

The pair script. Each row is one model fitted on the base build as written there and on the slice's as
written here, the same seed, 20 sweeps: train draws, sigma and the coefficients. Rows 01 to 06 and 21
must be `identical()`; rows 07 to 20 are equal to 1e-10 and the count of identical ones is recorded; rows
22 and 23 must differ.

| | base build | slice |
|---|---|---|
| 01 to 04 | no sd stated: two forests; three; under `subset` with weights and an offset; `bart` with a term | the same text |
| 05, 06 | logistic and probit, no sd stated | the same text |
| 07 to 10 | `sd = r` on the forest with no basis; on a factor forest; on both; on a numeric column and on two | `sd = r * sd(y)` |
| 11 to 13 | r stated with an offset; under `subset`; with weights | `r * sd(y - offset)`; `r * sd(y[kept])`; `r * sd(y)` |
| 14, 15 | a held factor; three forests each stating one | the same respelling |
| 16 to 18 | `$setLeafPrior(forests = )` after 5 sweeps; a copy after 5; a value read and written back | the same respelling; the same text |
| 19, 20 | logistic, r on both forests; `bart` with a term stating r | `r * pi / sqrt(3)`; `r * sd(y)` |
| 21 | probit, `sd = 0.7` on both forests | the same text |
| 22, 23 | gaussian, `sd = 0.7` on both; logistic the same | the same text: another prior |

Mutations, each expected to fail the named test and no gate before it:

- the unit taken to be the anchor under gaussian, the transform's scale forgotten: step 1 (a), and step
  3's "which model was fitted";
- the logistic unit taken to be 1: step 1 (a), and step 3's logistic row;
- a stated sd divided by the unit in R as well as in the engine: step 3's "which model was fitted" at
  every door; in the writer alone: step 4's path "after `$setLeafPrior`";
- a stated sd not divided at all at the data-object door (the list of sds over `bases`): step 3's row
  for that door;
- the default sent through the unit and back: step 3's "not stated", and the three bitwise compares;
- the reader returning the multiple of L: step 3's `identical()` to `forest(sd = 2)`;
- the unit derived again at a re-creation: step 1 (e), and step 4's path after `setResponse`;
- the scale derived after a swap from the stated number with no division: step 4's path after
  `$setForestBasis`;
- the unit carried on a control to other data: step 4's "other data";
- the skip of an equal write removed: step 4's "read then write";
- the unit computed with weights, or over every row before `subset`: step 3's "L written out";
- `extract` left returning `k.scale` over k: step 4's `extract` test;
- the stated number mirrored into the eight numbers and not into `sd`: step 4's path "a write, a copy
  and a second write";
- a configuration without `sd` read as one that states none: step 4's "saved before";
- the both-given refusal dropped: step 5.

Not a hot-path change: the division is made when a forest is built or restated, never in a sweep.

## NEWS

No new item: `forest()` and its `sd` are new in 1.0-0 and nothing released changes.

## What this leaves for the multiplier law

- What an sd is the size of where the coefficient is held or the basis is numeric: the table of "Context"
  with its rows off 1 is still true after this slice, in the response's units. The help says so, in those
  words, until the law changes them.
- The engine after this slice still takes the tip's four numbers per forest and one stated sd; the law
  replaces the four with the rule, and keeps `sd` not-a-number for "not stated", the recorded unit, the
  division and the reader's "sd in force". About 50 of this slice's engine lines are rewritten then.
- The reader's `amplitude.prior.scale` and `leaf.scale.factor`, `basis.row.norm`, the entry that says
  whether an sd was stated, `extract`'s list by forest and column, the printed block, and the per-column
  sds: none is built or renamed here.
- What the law may assume: every stated sd in the tree, the benchmarks and bartCause is in the response's
  units; a fit of forests with no basis or a factor, coefficients drawn, stated or not, is unchanged by
  the law to the bit, and can be gated so.

## What waits on what

- On push 3: the help's `sd` item is written over push 3's; every call in the tests and the help is in
  its spelling; a refusal that names a forest names it by its label where push 3 gives it one.
- On forest-defaults-by-kind: step 5 reads its record of what the caller stated and its finder of the
  forest with no basis; step 3 edits [`forestParams`](../../R/model.R) after it does; the third text
  stands beside its refusal of a leaf prior where no forest is without a basis, which comes first.
- To recheck once both have landed: the counts of "Context" (43 forests, 33 writes, 15 assertions), by
  running the logging build again; that no baseline scenario has gained a stated sd; that
  [`resolveForests`](../../R/model.R) still holds the one-forest refusal step 5 removes; the names of
  the functions cited here.
- Order with the two slices after it: this one, then the kind by class, then the multiplier law. This
  slice and the kind by class share R files and are serial; neither needs the other. They are not one
  slice: this one changes arithmetic in the engine and a prior, that one is R only and moves no draw but
  at one door. It is not folded into the law either: the law's cleanest gate is that forests with no
  basis or a factor are bit for bit under it, stated sds included, and that gate exists only if the unit
  moved first.

## Out of scope, and where it goes

- The fitting function's `leaf.prior = normal(sd = )` as the size of the forest with no basis in a model
  of several: refused here with the form to write. It belongs with the reader's `normal()` objects and
  `$setLeafPrior`'s single form, in TODO `forest-prior-args`.
- The kind by class; the multiplier law with everything listed above; `leaf.prior = normal(sd = )` on
  `forest()` and `updateBasisScale`, after the merge to main (dec-B276).
- A named single `sd` on a forest (dec-A171, open): refused as today.
- To TODO as a new entry: the verbose summary of a model of several forests prints one tree count and
  "k prior fixed to 2", the first forest's model slots, and nothing of the others.

## Calls made in planning

- The slice stands alone, before the kind by class and the law. Its draws do move, for a stated sd left
  as written, so it is not a renaming that could ride with either; the reasons against folding it are
  under "What waits on what".
- The engine converts, against a unit it derives and R records. The alternative, R dividing by
  `sd(y - offset)` before the bridge, is 40 lines shorter and was built as the stand-in. It was not
  chosen: R's estimate differs from the engine's in the last place, so the reader could not return a
  stated number exactly; and a `forests` list or a control resolved against one data object and run on
  another would carry the first one's unit in silence, where the engine derives the unit of the data it
  is given.
- The record is one more number beside the anchor. It could be recomputed from the anchor and the
  recorded response range, which live in two places and would have to agree to the bit.
- A stated number rides its own entry, `sd`, and the eight numbers keep the default. Putting the stated
  number of response units into the fourth or seventh number would make one slot mean two units by
  whether a flag is set, and a sampler saved before the slice would be misread.
- `extract(type = "leaf.prior.sd")` changes its numbers here and its shape with the law. The addendum to
  the design has every reader in units of L until this slice; at the tip `extract` is already in the
  response's units, of the forest's own standard deviation before its coefficient, which is neither the
  stated number nor what dec-B253 and dec-B275 say it returns. Leaving it would keep a reader and an
  `extract` that disagree, as they do today.
- `amplitude.prior.scale` and `leaf.scale.factor` stay as the engine's multiples of L and are documented
  so. They are entries the law removes; converting them would change two pins to change them again.
- A stated sd on a response with no spread is refused. Today it is taken and every prior is 0. Dividing
  by a unit of 0 would hand the engine an infinite scale.
- One forest takes `sd` in this slice, as written-surface's table and the design say it does once the
  unit is one. The fitting function's `normal(sd = )` beside several forests is not served here: the
  design gives that to the slice that trades `normal()` objects, and the text now says where to write it.
- A test's stated sd is respelled only where the test means a multiple of the response's scale. Every
  other one is left to state units; the alternative, respelling all 43, keeps draws nothing pins.
- The pair script is run at landing and not tracked: its old side needs the base build. What stays in
  the suite is step 3's pins of unstated fits and step 4's identities between paths.
- Measured for this plan on a stand-in, not on the design: the 15 assertions, the 8 of 18 pairs and the
  two "left as written" distances come from a build that converts in R. The engine's conversion divides
  by a number one unit in the last place away, so the counts of identical pairs may differ; the class
  does not.
